"""
IMPORTS
"""

## bUILT-iN pACKAGES
import sys
import pickle
from collections import defaultdict
import json
import numpy as np
import math
import matplotlib.pyplot as plt
np.set_printoptions(suppress=True)
from IPython.display import display

## lAST sUMMER
from romtoslp import *
from sentences import *
from DCS import *

## tHIS sUMMER
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
    if not neuralnet.outer_relu:
        conflicts_Dict = Get_Conflicts(nodelist)

        featVMat = Get_Feat_Vec_Matrix(nodelist, conflicts_Dict)

        (WScalarMat, SigmoidGateOutput) = Get_W_Scalar_Matrix_from_FeatVect_Matrix(featVMat, nodelist, conflicts_Dict, neuralnet)
        return (conflicts_Dict, featVMat, WScalarMat, SigmoidGateOutput)
    else:
        conflicts_Dict = Get_Conflicts(nodelist)

        featVMat = Get_Feat_Vec_Matrix(nodelist, conflicts_Dict)

        WScalarMat = Get_W_Scalar_Matrix_from_FeatVect_Matrix(featVMat, nodelist, conflicts_Dict, neuralnet)
        return (conflicts_Dict, featVMat, WScalarMat)

# NEW LOSS FUNCTION
def GetLoss(_mst_adj_graph, _mask_de_correct_edges, _negLogLikelies):
    _negLogLikelies = _negLogLikelies.copy()
    _negLogLikelies[~_mst_adj_graph] = 0
    _negLogLikelies[~_mask_de_correct_edges] *= -1 # BAKA!!! Check before you try to fix this again
    return np.sum(_negLogLikelies)  

"""
################################################################################################
##############################  GET A FILENAME TO SAVE WEIGHTS  ################################
################################################################################################
"""
import time
st = str(int((time.time() * 1e6) % 1e13))
log_name = 'logs/train_nnet_t{}.out'.format(st)
p_name = 'outputs/train_nnet_t{}.p'.format(st)
print('nEURAL nET wILL bE sAVED hERE: ', p_name)

"""
################################################################################################
##############################  MAIN FUNCTION  #################################################
################################################################################################
"""
trainingStatus = defaultdict(lambda: bool(False))


"""
################################################################################################
##############################  TRAIN FUNCTION  ################################################
################################################################################################
"""

def train(loaded_SKT, loaded_DCS, n_trainset = -1, iterationPerBatch = 10, filePerBatch = 20, _debug = True):
    # Train
    if n_trainset == -1:
        totalBatchToTrain = 20
    else:
        totalBatchToTrain = math.ceil(n_trainset/filePerBatch)
    
    for iterout in range(totalBatchToTrain):
        # Change current batch
        trainer.Save(p_name)
        print('Batch: ', iterout)
        files_for_batch = TrainFiles[iterout*filePerBatch:(iterout + 1)*filePerBatch]
        print(files_for_batch)
        # trainer.Load('outputs/neuralnet_trained.p')
        
        # Run few times on same set of files
        for iterin in range(iterationPerBatch):
            print('ITERATION IN', iterin)        
            for fn in files_for_batch:
                sentenceObj = loaded_SKT[fn]
                dcsObj = loaded_DCS[fn]
                if trainingStatus[sentenceObj.sent_id]:
                    continue
                # trainer.Save('outputs/saved_trainer.p')
                try:
                    trainer.Train(sentenceObj, dcsObj, _debug)
                except (IndexError, KeyError) as e:
                    print('\x1b[31mFailed: {} \x1b[0m'.format(sentenceObj.sent_id))
            sys.stdout.flush() # Flush IO buffer 
    trainer.Save(p_name)
    
    sys.stdout.flush() # Flush IO buffer 
                
def test(loaded_SKT, loaded_DCS, n_testSet = -1, _testFiles = None, n_checkpt = 100):
    total_lemma = 0;
    correct_lemma = 0;

    total_word = 0;
    total_output_nodes = 0
    correct_word = 0;
    file_counter = 0
    if _testFiles is None:
        if n_testSet == -1:
            _testFiles = TestFiles
        else:
            _testFiles = TestFiles[0:n_testSet]
    else:
        if n_testSet == -1:
            _testFiles = _testFiles
        else:
            _testFiles = _testFiles[0:n_testSet]
            
    recalls = []
    recalls_of_word = []
    precisions = []
    precisions_of_words = []
    for fn in _testFiles:
        if file_counter % n_checkpt == 0:
            print(file_counter,' Checkpoint... ')
            sys.stdout.flush() # Flush IO buffer 
        file_counter += 1
        sentenceObj = loaded_SKT[fn]
        dcsObj = loaded_DCS[fn]        
        try:
            (word_match, lemma_match, n_dcsWords, n_output_nodes) = trainer.Test(sentenceObj, dcsObj)
            
            recalls.append(lemma_match/n_dcsWords)
            recalls_of_word.append(word_match/n_dcsWords)
            
            precisions.append(lemma_match/n_output_nodes)
            precisions_of_words.append(word_match/n_output_nodes)
            
            total_lemma += n_dcsWords
            total_word += n_dcsWords
            
            total_output_nodes += n_output_nodes            
            
            correct_lemma += lemma_match
            correct_word += word_match
        except (IndexError, KeyError) as e:
            print('Failed!')        

    print('Avg. Micro Recall of Lemmas: {}'.format(np.mean(np.array(recalls))))
    print('Avg. Micro Recall of Words: {}'.format(np.mean(np.array(recalls_of_word))))
    print('Avg. Micro Precision of Lemmas: {}'.format(np.mean(np.array(precisions))))
    print('Avg. Micro Precision of Words: {}'.format(np.mean(np.array(precisions_of_words))))
    
    return (recalls, recalls_of_word, precisions, precisions_of_words)

# NEW FUNCTION
def GetLoss(_mst_adj_graph, _mask_de_correct_edges, _WScalarMat):
    _WScalarMat = _WScalarMat.copy()
    _WScalarMat[_mst_adj_graph&(~_mask_de_correct_edges)] *= -1 # BAKA!!! Check before you try to fix this again
    _WScalarMat[~_mst_adj_graph] = 0
    return np.sum(_WScalarMat)

"""
################################################################################################
#############################   TRAINER CLASS DEFINITION  ######################################
################################################################################################
"""
class Trainer:
    def __init__(self):
        self.hidden_layer_size = 300
        self._edge_vector_dim = WD._edge_vector_dim
        # self._full_cnglist = list(WD.mat_cngCount_1D)
        self.neuralnet = NN(self._edge_vector_dim, self.hidden_layer_size, outer_relu=True)
        self.history = defaultdict(lambda: list())
        
    def Reset(self):
        self.neuralnet = NN(self._edge_vector_dim, self.hidden_layer_size)
        self.history = defaultdict(lambda: list())
        
    def Save(self, filename):
        #pickle.dump({'nnet': self.neuralnet, 'history': dict(self.history)}, open(filename, 'wb'))
        return
        
    
    def Load(self, filename):
        loader = pickle.load(open(filename, 'rb'))
        self.neuralnet.U = loader['U']
        self.neuralnet.W = loader['W']
        
    def Test(self, sentenceObj, dcsObj):
        neuralnet = self.neuralnet
        minScore = np.inf
        minMst = None
        
        # startT = time.time()
        try:
            (nodelist, nodelist_correct, _) = GetTrainingKit(sentenceObj, dcsObj)
            # nodelist = GetNodes(sentenceObj)
        except IndexError:
            print('\x1b[31mError with {} \x1b[0m'.format(sentenceObj.sent_id))
            return (0, 0, 0)
            
        conflicts_Dict = Get_Conflicts(nodelist)
        conflicts_Dict_correct = Get_Conflicts(nodelist_correct)
        
        # print('Load Nodelist + GetConflicts Time: ', time.time() - startT)
        # startT = time.time()
        
        featVMat = Get_Feat_Vec_Matrix(nodelist, conflicts_Dict)
        featVMat_correct = Get_Feat_Vec_Matrix(nodelist_correct, conflicts_Dict_correct)
        
        # print('Get Edge Features Time: ', time.time() - startT)
        # startT = time.time()
        
        if not self.neuralnet.outer_relu:
            (WScalarMat, SigmoidGateOutput) = Get_W_Scalar_Matrix_from_FeatVect_Matrix(featVMat, nodelist, conflicts_Dict, neuralnet)
        else:
            WScalarMat = Get_W_Scalar_Matrix_from_FeatVect_Matrix(featVMat, nodelist, conflicts_Dict, neuralnet)
        
        # print('NeuralNet Time: ', time.time() - startT)
        # startT = time.time()
        
        # Get all MST
        for source in range(len(nodelist)):
            (mst_nodes, mst_adj_graph, _) = MST(nodelist, WScalarMat, conflicts_Dict, source)
            # print('.', end = '')
            score = GetMSTWeight(mst_nodes, WScalarMat)
            if(score < minScore):
                minScore = score
                minMst = mst_nodes
        dcsLemmas = [[rom_slp(l) for l in arr]for arr in dcsObj.lemmas]
        word_match = 0
        lemma_match = 0
        n_output_nodes = 0
        for chunk_id, wdSplit in minMst.items():
            for wd in wdSplit:
                n_output_nodes += 1
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
        
        # print('All MST Time: ', time.time() - startT)
        # print('Node Count: ', len(nodelist))
#         print('\nFull Match: {}, Partial Match: {}, OutOf {}, NodeCount: {}, '.\
#               format(word_match, lemma_match, len(dcsLemmas), len(nodelist)))
        return (word_match, lemma_match, len(dcsLemmas), n_output_nodes)
    
    def Train(self, sentenceObj, dcsObj, _debug = True):
        # Hyperparameter for hinge loss: m
        M_hinge = 14
        
        """ Pre-Process DCS and SKT to get all Nodes etc. """
        try:
            (nodelist, nodelist_correct, nodelist_to_correct_mapping) = GetTrainingKit(sentenceObj, dcsObj)
        except IndexError as e:
            # print('\x1b[31mError with {} \x1b[0m'.format(sentenceObj.sent_id))
            # print(e)
            return
        
        """ FORM MAXIMUM(ENERGY) SPANNING TREE OF THE GOLDEN GRAPH : WORST GOLD STRUCTURE """
        if not self.neuralnet.outer_relu:
            raise Exception('Support for Non-outer-relu removed')
        else:
            (conflicts_Dict_correct, featVMat_correct, WScalarMat_correct) = GetGraph(nodelist_correct, self.neuralnet)
        source = 0
        
        """ Find the max spanning tree : negative Weight matrix passed """
        (max_st_gold_ndict, max_st_adj_gold_small, _) = MST(nodelist_correct, -WScalarMat_correct, conflicts_Dict_correct, source)
        energy_gold_max_ST = np.sum(WScalarMat_correct[max_st_adj_gold_small])
        
        """ Convert correct spanning tree graph adj matrix to full marix dimensions """
        """ Create full-size adjacency matrix for correct_mst_small """
        nodelen = len(nodelist)
        max_st_adj_gold = np.ndarray((nodelen, nodelen), np.bool)*False # T_STAR
        for i in range(max_st_adj_gold_small.shape[0]):
            for j in range(max_st_adj_gold_small.shape[1]):
                max_st_adj_gold[nodelist_to_correct_mapping[i], nodelist_to_correct_mapping[j]] = max_st_adj_gold_small[i, j]
        
        """ Delta(Margin) Function : MASK FOR WHICH NODES IN NODELIST BELONG TO DCS """
        gold_nodes_mask = np.array([False]*len(nodelist))
        for i in range(len(nodelist_correct)):
            gold_nodes_mask[nodelist_to_correct_mapping[nodelist_correct[i].id]] = True
        margin_f = lambda nodes_mask: np.sum(nodes_mask&gold_nodes_mask)**1.7
        
        """ FOR ALL POSSIBLE MST FROM THE COMPLETE GRAPH """
        if not self.neuralnet.outer_relu:
            raise Exception('Support for Non-outer-relu removed')
        else:
            (conflicts_Dict, featVMat, WScalarMat) = GetGraph(nodelist, self.neuralnet)

        """ For each node - Find MST with that source"""
        min_STx = None # Min Energy spanning tree with worst margin with gold_STx
        min_marginalized_energy = np.inf
        
        for source in range(len(nodelist)):
            (mst_nodes, mst_adj_graph, mst_nodes_bool) = MST(nodelist, WScalarMat, conflicts_Dict, source) # T_X
            # print('.', end = '')
           
            marginalized_dist = np.sum(WScalarMat[mst_adj_graph]) - margin_f(mst_nodes_bool)
            if marginalized_dist < min_marginalized_energy:
                min_marginalized_energy = marginalized_dist
                min_STx = mst_adj_graph
            # Energy diff should all be negative
#             print('Source: [{}], Del:{}, Energy_margin: {:.3f}, Energy: {:.3f}, GE:{:.3f}'.\
#                   format(source, margin_f(mst_nodes_bool), marginalized_dist,  np.sum(WScalarMat[mst_adj_graph]), energy_gold_max_ST))

        """ Gradient Descent """
        # FOR MOST OFFENdING Y
        doBpp = False
        
        Total_Loss = energy_gold_max_ST - min_marginalized_energy
#         print('Total Loss: ', Total_Loss)
        if Total_Loss > 0:
            doBpp = True
            dLdOut = np.zeros_like(WScalarMat)
            
            if not self.neuralnet.outer_relu:
                # For new loss function with sigmoid on neuralnet
                raise Exception('Support for relu has been removed')
            else:
                # For new loss function with ReLU on neuralnet
                dLdOut[max_st_adj_gold] = 1
                dLdOut[min_STx] = -1

        if doBpp:
            if _debug:
                print('{}. '.format(sentenceObj.sent_id), end = '')
            self.neuralnet.Back_Prop(dLdOut, len(nodelist), featVMat, _debug)
        else:
            trainingStatus[sentenceObj.sent_id] = True # Means tranining done on this file
            #print('Don\'t do BPP')
            pass
        
        Total_Loss /= len(nodelist)
        self.history[sentenceObj.sent_id].append(Total_Loss)
#         print("\nFileKey: %s, Loss: %6.3f, Loss2: %6.3f" % (sentenceObj.sent_id, Total_Loss, Loss_2))

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

    dataset_4k_1k = pickle.load(open('../SmallDataset_4K_1K.p', 'rb'))
    TrainFiles = dataset_4k_1k['TrainFiles']
    TestFiles = dataset_4k_1k['TestFiles']

    dataset_6k_3k = pickle.load(open('../SmallDataset_6K_3K.p', 'rb'))
    TrainFiles_2 = dataset_6k_3k['TrainFiles']
    TestFiles_2 = dataset_6k_3k['TestFiles']

    matDB = MatDB.MatDB()
    InitModule(matDB)

    print('PRE-TRAIN ACCURACIES:')
    _ = test(loaded_SKT, loaded_DCS, n_testSet = 1000, _testFiles = TestFiles, n_checkpt = 100)

    print('TRAINING ROUND 1::')
    train(loaded_SKT, loaded_DCS, n_trainset = 1000, filePerBatch = 20, iterationPerBatch = 8, _debug=False)

    print('POST-TRAIN ACCURACIES:')
    _ = test(loaded_SKT, loaded_DCS, n_testSet = 1000, _testFiles = TestFiles, n_checkpt = 100)

    # print('TRAINING ROUND 2::')
    # train(loaded_SKT, loaded_DCS, 1000, _debug = False)

    print('POST-TRAIN ACCURACIES DIFF SET:')
    _ = test(loaded_SKT, loaded_DCS, n_testSet = 3000, _testFiles = TestFiles_2, n_checkpt = 100)
    # print ("Not Implemented")
