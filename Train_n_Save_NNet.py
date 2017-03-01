"""
IMPORTS
"""

import pickle
import numpy as np
from sentences import *
from DCS import *
from collections import defaultdict
import math

from romtoslp import *

from MatDB import *
import word_definite as WD
from heap_n_PrimMST import *
from nnet import *


"""
################################################################################################
###########################  LOAD SENTENCE AND DCS OBJECT FILES  ###############################
################################################################################################
"""
# loaded_SKT = pickle.load(open('../Simultaneous_CompatSKT_10K.p', 'rb'))
# loaded_DCS = pickle.load(open('../Simultaneous_DCS_10K.p', 'rb'))


"""
################################################################################################
######################  CREATE SEVERAL DATA STRUCTURES FROM SENTENCE/DCS  ######################
###########################  NODELIST, ADJACENCY LIST, GRAPH, HEAP #############################
"""
def GetTrainingKit(sentenceObj, dcsObj):
    nodelist = GetNodes(sentenceObj)
    
    # Nodelist with only the correct_nodes
    nodelist2 = GetNodes(sentenceObj)
    nodelist2_to_correct_mapping = {}
    nodelist_correct = []
    search_key = 0
    first_key = 0
    for chunk_id in range(len(dcsObj.lemmas)):
        while nodelist2[first_key].chunk_id != chunk_id:
            first_key += 1
        for j in range(len(dcsObj.lemmas[chunk_id])):
            search_key = first_key
            while (nodelist2[search_key].lemma != rom_slp(dcsObj.lemmas[chunk_id][j])) or (nodelist2[search_key].cng != dcsObj.cng[chunk_id][j]):
                search_key += 1
                if search_key >= len(nodelist2) or nodelist2[search_key].chunk_id > chunk_id:
                    break
    #         print((rom_slp(dcsObj.lemmas[chunk_id][j]), dcsObj.cng[chunk_id][j]))
    #         print(nodelist[search_key])
            nodelist2_to_correct_mapping[len(nodelist_correct)] = search_key
            nodelist_correct.append(nodelist2[search_key])
    return (nodelist, nodelist_correct, nodelist2_to_correct_mapping)
    

def GetGraph(nodelist, neuralnet):
    conflicts_Dict = Get_Conflicts(nodelist)

    featVMat = Get_Feat_Vec_Matrix(nodelist, conflicts_Dict)

    WScalarMat = Get_W_Scalar_Matrix_from_FeatVect_Matrix(featVMat, nodelist, conflicts_Dict, neuralnet)
    return (conflicts_Dict, featVMat, WScalarMat)

"""
################################################################################################
##############################  MAIN FUNCTION  #################################################
################################################################################################
"""
def main(loaded_SKT, loaded_DCS):
    # Train
    totalBatchToTrain = 1000
    filePerBatch = 30
    iterationPerBatch = 5
    for iterout in range(totalBatchToTrain):
        # Change current batch
        print('ITERATION OUT', iterout)
        perm = np.random.permutation(len(loaded_SKT))[0:filePerBatch]
        print('Permutation: ', perm)
        # trainer.Load('outputs/neuralnet_trained.p')
        
        # Run few times on same set of files
        for iterin in range(iterationPerBatch):
            print('ITERATION IN', iterin)        
            for fn in perm:
                sentenceObj = loaded_SKT[list(loaded_SKT.keys())[fn]]
                dcsObj = loaded_DCS[list(loaded_SKT.keys())[fn]]
                # trainer.Save('outputs/saved_trainer.p')
                trainer.Train(sentenceObj, dcsObj)     
                
def test(loaded_SKT, loaded_DCS, perm = None):
    if perm is None:
        perm = np.random.permutation(100)[0:10]
        print('Permutation: ', perm)
        for fn in perm:
            sentenceObj = loaded_SKT[list(loaded_SKT.keys())[fn]]
            dcsObj = loaded_DCS[list(loaded_SKT.keys())[fn]]   
            trainer.Test(sentenceObj, dcsObj)
    else:
        print('Permutation: ', perm)
        for fn in perm:
            sentenceObj = loaded_SKT[list(loaded_SKT.keys())[fn]]
            dcsObj = loaded_DCS[list(loaded_SKT.keys())[fn]]
            trainer.Test(sentenceObj, dcsObj)
                
class Trainer:
    def __init__(self):
        self._edge_vector_dim = WD._edge_vector_dim
        self._full_cnglist = list(WD.mat_cngCount_1D)
        self.neuralnet = NN(self._edge_vector_dim, 200)
        self.history = defaultdict(lambda: list())
        
    def Reset(self):
        self.neuralnet = NN(self._edge_vector_dim, 200)
        self.history = defaultdict(lambda: list())
        
    def Save(self, filename):
        pickle.dump({'nnet': self.neuralnet, 'history': dict(self.history)}, open(filename, 'wb'))
        
    
    def Load(self, filename):
        o = pickle.load(open(filename, 'rb'))
        self.neuralnet = o['nnet']
        self.history = defaultdict(lambda: list(), o['history'])
        
    def Train(self, sentenceObj, dcsObj):
        try:
            (nodelist, nodelist_correct, nodelist_to_correct_mapping) = GetTrainingKit(sentenceObj, dcsObj)
        except IndexError as e:
            # print('\x1b[31mError with {} \x1b[0m'.format(sentenceObj.sent_id))
            print(e)
            return
            
        
        """ FOR MST OF GRAPH WITH ONLY CORRECT SET OF NODES """
        (conflicts_Dict_correct, featVMat_correct, WScalarMat_correct) = GetGraph(nodelist_correct, self.neuralnet)
        source = 0
        (mst_nodes_correct, mst_adj_graph_correct_0) = MST(nodelist_correct, WScalarMat_correct, conflicts_Dict_correct, source)
        W_star = GetMSTWeight(mst_nodes_correct, WScalarMat_correct)
        
        """ FOR ALL POSSIBLE MST FROM THE COMPLETE GRAPH """
        (conflicts_Dict, featVMat, WScalarMat) = GetGraph(nodelist, self.neuralnet)        
        
        Total_Loss = 0

        # Convert correct spanning tree graph adj matrix to actual marix dimensions
        # Create full-size adjacency matrix for correct_mst
        nodelen = len(nodelist)
        mst_adj_graph_correct = np.ndarray((nodelen, nodelen), np.bool)*False
        for i in range(mst_adj_graph_correct_0.shape[0]):
            for j in range(mst_adj_graph_correct_0.shape[1]):
                mst_adj_graph_correct[nodelist_to_correct_mapping[i], nodelist_to_correct_mapping[j]] = \
                mst_adj_graph_correct_0[i, j]

        """ For each node - Find MST with that source"""
        dLdS = np.zeros(WScalarMat.shape)
        for source in range(len(nodelist)):
            (mst_nodes, mst_adj_graph) = MST(nodelist, WScalarMat, conflicts_Dict, source)
            # print('.', end = '')

            """ Gradient Descent"""
            
            Total_Loss += (W_star - GetMSTWeight(mst_nodes, WScalarMat))
            dLdS_inner = 1 / WScalarMat
            dLdS_inner[mst_adj_graph_correct == 1] *= -1
            dLdS_inner = dLdS_inner*(mst_adj_graph^mst_adj_graph_correct)
            dLdS += dLdS_inner        
        self.neuralnet.Back_Prop(dLdS/len(nodelist), len(nodelist), featVMat)

        Total_Loss /= len(nodelist)
        self.history[sentenceObj.sent_id].append(Total_Loss)
        print("\nFileKey: %s, Loss: %6.3f, Original MSTScore: %6.3f" % (sentenceObj.sent_id, Total_Loss, W_star))
        
    def Test(self, sentenceObj, dcsObj):
        neuralnet = self.neuralnet
        minScore = np.inf
        minMst = None
        try:
            (nodelist, nodelist_correct, _) = GetTrainingKit(sentenceObj, dcsObj)
            # nodelist = GetNodes(sentenceObj)
        except IndexError:
            print('\x1b[31mError with {} \x1b[0m'.format(sentenceObj.sent_id))
            return
            
        conflicts_Dict = Get_Conflicts(nodelist)
        conflicts_Dict_correct = Get_Conflicts(nodelist_correct)
        
        featVMat = Get_Feat_Vec_Matrix(nodelist, conflicts_Dict)
        featVMat_correct = Get_Feat_Vec_Matrix(nodelist_correct, conflicts_Dict_correct)
        
        WScalarMat = Get_W_Scalar_Matrix_from_FeatVect_Matrix(featVMat, nodelist, conflicts_Dict, neuralnet)
        # print(WScalarMat)
        # Get all MST
        for source in range(len(nodelist)):
            (mst_nodes, mst_adj_graph) = MST(nodelist, WScalarMat, conflicts_Dict, source)
            print('.', end = '')
            score = GetMSTWeight(mst_nodes, WScalarMat)
            if(score < minScore):
                print('.', end = '')
                minScore = score
                minMst = mst_nodes
        dcsLemmas = [[rom_slp(l) for l in arr]for arr in dcsObj.lemmas]
        full_match = 0
        partial_match = 0
        for chunk_id, wdSplit in mst_nodes.items():
            for wd in wdSplit:
                # Match lemma
                search_result = [i for i, j in enumerate(dcsLemmas[chunk_id]) if j == wd.lemma]
                if len(search_result) > 0:
                    partial_match += 1
                # Match CNG
                for i in search_result:
                    if(dcsObj.cng[chunk_id][i] == str(wd.cng)):
                        full_match += 1
                        # print(wd.lemma, wd.cng)
                        break
        dcsLemmas = [l for arr in dcsObj.lemmas for l in arr]
        print('\nFull Match: {}, Partial Match: {}, OutOf {}, NodeCount: {}, '.\
              format(full_match, partial_match, len(dcsLemmas), len(nodelist)))


trainer = None
def InitModule(_matDB):
    global WD, trainer
    WD.word_definite_extInit(_matDB)
    trainer = Trainer()

"""
################################################################################################
################################################################################################
################################################################################################
"""
if __name__ == '__main__':
    loaded_SKT = pickle.load(open('../Simultaneous_CompatSKT_10K.p', 'rb'))
    loaded_DCS = pickle.load(open('../Simultaneous_DCS_10K.p', 'rb'))
    main(loaded_SKT, loaded_DCS)
    # print ("Not Implemented")
