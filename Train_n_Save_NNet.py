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

import MatDB
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
################################################################################################
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

    (WScalarMat, SigmoidGateOutput) = Get_W_Scalar_Matrix_from_FeatVect_Matrix(featVMat, nodelist, conflicts_Dict, neuralnet)
    return (conflicts_Dict, featVMat, WScalarMat, SigmoidGateOutput)

# NEW LOSS FUNCTION
def GetLoss(_mst_adj_graph, _mask_de_correct_edges, _negLogLikelies):
    _negLogLikelies = _negLogLikelies.copy()
    _negLogLikelies[~_mst_adj_graph] = 0
    _negLogLikelies[~_mask_de_correct_edges] *= -1 # BAKA!!! Check before you try to fix this again
    return np.sum(_negLogLikelies)  

"""
################################################################################################
##############################  MAIN FUNCTION  #################################################
################################################################################################
"""
def main(loaded_SKT, loaded_DCS):
    # Train
    totalBatchToTrain = 50
    filePerBatch = 15
    iterationPerBatch = 5
    _skt_file_list = list(loaded_SKT.keys())
    for iterout in range(totalBatchToTrain):
        # Change current batch
        print('ITERATION OUT', iterout)
        perm = np.random.permutation(len(loaded_SKT))[0:filePerBatch]
        # print('Permutation: ', perm)
        # trainer.Load('outputs/neuralnet_trained.p')
        
        # Run few times on same set of files
        for iterin in range(iterationPerBatch):
            print('ITERATION IN', iterin)        
            for fn in perm:
                sentenceObj = loaded_SKT[_skt_file_list[fn]]
                dcsObj = loaded_DCS[_skt_file_list[fn]]
                # trainer.Save('outputs/saved_trainer.p')
                try:
                    trainer.Train(sentenceObj, dcsObj)
                except:
                    pass
                
def test(loaded_SKT, loaded_DCS, perm = None, n_testSet = 10):
    _skt_file_list = list(loaded_SKT.keys())
    total_lemma = 0;
    correct_lemma = 0;

    total_word = 0;
    correct_word = 0;
    file_counter = 0
    if perm is None:
        n_testSet = min(n_testSet, len(_skt_file_list))
        perm = np.random.permutation(len(_skt_file_list))[0:n_testSet]
        # print('Permutation: ', perm)
        for fn in perm:
            if file_counter % 50 == 0:
                print('Checkpoint: ', file_counter)
            file_counter += 1
            sentenceObj = loaded_SKT[_skt_file_list[fn]]
            dcsObj = loaded_DCS[_skt_file_list[fn]]
            try:
                (word_match, lemma_match, n_dcsWords) = trainer.Test(sentenceObj, dcsObj)
                total_word += n_dcsWords
                correct_lemma += lemma_match
                correct_word += word_match
            except:
                # print('Failed')
                pass
    else:
        # print('Permutation: ', perm)
        for fn in perm:
            if file_counter % 50 == 0:
                print('Checkpoint: ', file_counter)
            file_counter += 1
            sentenceObj = loaded_SKT[_skt_file_list[fn]]
            dcsObj = loaded_DCS[_skt_file_list[fn]]
            try:
                (word_match, lemma_match, n_dcsWords) = trainer.Test(sentenceObj, dcsObj)
                total_lemma += n_dcsWords
                total_word += n_dcsWords
                correct_lemma += lemma_match
            except:
                # print('Failed!')
                pass
            
            correct_word += word_match
    print()
    print('Micro Accuracy[lemma]: {}'.format(100*correct_lemma/total_word))
    print('Micro Accuracy[word]: {}'.format(100*correct_word/total_word))

"""
################################################################################################
#############################   TRAINER CLASS DEFINITION  ######################################
################################################################################################
"""
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
        
    def Test(self, sentenceObj, dcsObj):
        neuralnet = self.neuralnet
        minScore = np.inf
        minMst = None
        try:
            (nodelist, nodelist_correct, _) = GetTrainingKit(sentenceObj, dcsObj)
            # nodelist = GetNodes(sentenceObj)
        except IndexError:
            # print('\x1b[31mError with {} \x1b[0m'.format(sentenceObj.sent_id))
            return (0, 0, 0)
            
        conflicts_Dict = Get_Conflicts(nodelist)
        conflicts_Dict_correct = Get_Conflicts(nodelist_correct)
        
        featVMat = Get_Feat_Vec_Matrix(nodelist, conflicts_Dict)
        featVMat_correct = Get_Feat_Vec_Matrix(nodelist_correct, conflicts_Dict_correct)
        
        (WScalarMat, SigmoidGateOutput) = Get_W_Scalar_Matrix_from_FeatVect_Matrix(featVMat, nodelist, conflicts_Dict, neuralnet)
        # print(WScalarMat)
        # Get all MST
        for source in range(len(nodelist)):
            (mst_nodes, mst_adj_graph) = MST(nodelist, WScalarMat, conflicts_Dict, source)
            # print('.', end = '')
            score = GetMSTWeight(mst_nodes, WScalarMat)
            if(score < minScore):
                minScore = score
                minMst = mst_nodes
        dcsLemmas = [[rom_slp(l) for l in arr]for arr in dcsObj.lemmas]
        word_match = 0
        lemma_match = 0
        for chunk_id, wdSplit in mst_nodes.items():
            for wd in wdSplit:
                # Match lemma
                search_result = [i for i, j in enumerate(dcsLemmas[chunk_id]) if j == wd.lemma]
                if len(search_result) > 0:
                    lemma_match += 1
                # Match CNG
                for i in search_result:
                    if(dcsObj.cng[chunk_id][i] == str(wd.cng)):
                        word_match += 1
                        # print(wd.lemma, wd.cng)
                        break
        dcsLemmas = [l for arr in dcsObj.lemmas for l in arr]
        # print('\nFull Match: {}, Partial Match: {}, OutOf {}, NodeCount: {}, '.\
        #       format(word_match, lemma_match, len(dcsLemmas), len(nodelist)))
        return (word_match, lemma_match, len(dcsLemmas))
    
    def Train(self, sentenceObj, dcsObj):
        try:
            (nodelist, nodelist_correct, nodelist_to_correct_mapping) = GetTrainingKit(sentenceObj, dcsObj)
        except IndexError as e:
            # print('\x1b[31mError with {} \x1b[0m'.format(sentenceObj.sent_id))
            # print(e)
            return
        
        """ CREATE A MASK FOR ALL EDGES BETWEEN CORRECT NODE PAIRS"""
        mask_de_correct_edges = np.ndarray((len(nodelist), len(nodelist)), np.bool)*False
        for n1 in nodelist_to_correct_mapping.values():
            for n2 in nodelist_to_correct_mapping.values():
                if n1 != n2:
                    mask_de_correct_edges[n1, n2] = 1
                    mask_de_correct_edges[n2, n1] = 1
                
        
        """ FOR MST OF GRAPH WITH ONLY CORRECT SET OF NODES """
        (conflicts_Dict_correct, featVMat_correct, WScalarMat_correct, _) = GetGraph(nodelist_correct, self.neuralnet)
        source = 0
        (mst_nodes_correct, mst_adj_graph_correct_0) = MST(nodelist_correct, WScalarMat_correct, conflicts_Dict_correct, source)
        W_star = GetMSTWeight(mst_nodes_correct, WScalarMat_correct) # Total weight of correct MST
        
        """ FOR ALL POSSIBLE MST FROM THE COMPLETE GRAPH """
        (conflicts_Dict, featVMat, WScalarMat, SigmoidGateOutput) = GetGraph(nodelist, self.neuralnet)
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
        dLdOut = np.zeros(WScalarMat.shape)
        for source in range(len(nodelist)):
            (mst_nodes, mst_adj_graph) = MST(nodelist, WScalarMat, conflicts_Dict, source)
            # print('.', end = '')

            """ Gradient Descent"""
            Total_Loss += GetLoss(mst_adj_graph, mask_de_correct_edges, WScalarMat)
            
            # For new loss function
            dLdOut_inner = (1 - SigmoidGateOutput)
            dLdOut_inner[~mst_adj_graph] = 0 # Edges not in the current tree
            dLdOut_inner[mst_adj_graph^mask_de_correct_edges] = -dLdOut_inner[mst_adj_graph^mask_de_correct_edges]
            
            dLdOut += dLdOut_inner
        self.neuralnet.Back_Prop(dLdOut/len(nodelist), len(nodelist), featVMat)        
        Total_Loss /= len(nodelist)
        self.history[sentenceObj.sent_id].append(Total_Loss)
        # print("\nFileKey: %s, Loss: %6.3f, Original MSTScore: %6.3f" % (sentenceObj.sent_id, Total_Loss, W_star))


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
    matDB = MatDB.MatDB()
    InitModule(matDB)

    loaded_SKT = pickle.load(open('../Simultaneous_CompatSKT_10K.p', 'rb'))
    loaded_DCS = pickle.load(open('../Simultaneous_DCS_10K.p', 'rb'))

    print('PRE-TRAIN ACCURACIES:')
    test(loaded_SKT, loaded_DCS, n_testSet=1000)

    print('TRAINING ROUND 1::')
    main(loaded_SKT, loaded_DCS)

    print('POST-TRAIN ACCURACIES:')
    test(loaded_SKT, loaded_DCS, n_testSet=1000)

    print('TRAINING ROUND 2::')
    main(loaded_SKT, loaded_DCS)

    print('POST-TRAIN2 ACCURACIES:')
    test(loaded_SKT, loaded_DCS, n_testSet=1000)
    # print ("Not Implemented")
