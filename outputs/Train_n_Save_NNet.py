"""
IMPORTS
"""

## bUILT-iN pACKAGES
import sys, os, time, bz2, zlib, pickle, math, json, csv
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
np.set_printoptions(suppress=True)
from IPython.display import display

## lAST sUMMER
from romtoslp import *
from sentences import *
from DCS import *
import MatDB
from heap_n_PrimMST import *
from ECL_MST import *

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
###########################  OPENS AND EXTRACTS DCS AND SKT DATA STRUCTURES  ###################
######################################  FROM BZ2 FILES  ########################################
"""
def open_dsbz2(filename):
    with bz2.BZ2File(filename, 'r') as f:
        loader = pickle.load(f)
    
    conflicts_Dict_correct = loader['conflicts_Dict_correct']
    nodelist_to_correct_mapping = loader['nodelist_to_correct_mapping']
    nodelist_correct = loader['nodelist_correct']
    featVMat_correct = loader['featVMat_correct']
    featVMat = loader['featVMat']
    conflicts_Dict = loader['conflicts_Dict']
    nodelist = loader['nodelist']
    
    return (nodelist_correct, conflicts_Dict_correct, featVMat_correct, nodelist_to_correct_mapping,\
            nodelist, conflicts_Dict, featVMat)

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
##############################  MAIN FUNCTION  #################################################
################################################################################################
"""
trainingStatus = defaultdict(lambda: bool(False))


"""
################################################################################################
##############################  TRAIN FUNCTION  ################################################
################################################################################################
"""

def train_generator(loaded_SKT, loaded_DCS, bz2_input_folder, n_trainset = -1, iterationPerBatch = 10, filePerBatch = 20, _debug = True):
    # Train
    if n_trainset == -1:
        n_trainset = len(TrainFiles)
        totalBatchToTrain = math.ceil(n_trainset/filePerBatch)
    else:
        totalBatchToTrain = math.ceil(n_trainset/filePerBatch)
    
    register_nnet(trainer.neuralnet, bz2_input_folder)
    for iterout in range(totalBatchToTrain):
        # Add timer
        startT = time.time()

        # Change current batch
        if(iterout % 50 == 0):
            trainer.Save(p_name.replace('.p', '_i{}.p'.format(iterout)))
        else:
            trainer.Save(p_name)
        print('Batch: ', iterout)
        files_for_batch = TrainFiles[iterout*filePerBatch:(iterout + 1)*filePerBatch]
        print(files_for_batch)
        # trainer.Load('outputs/neuralnet_trained.p')
        try:
            # Run few times on same set of files
            for iterin in range(iterationPerBatch):
                print('ITERATION IN', iterin)        
                for fn in files_for_batch:
                    trainFileName = fn.replace('.ds.bz2', '.p2')
                    sentenceObj = loaded_SKT[trainFileName]
                    dcsObj = loaded_DCS[trainFileName]
                    if trainingStatus[sentenceObj.sent_id]:
                        continue
                    # trainer.Save('outputs/saved_trainer.p')
                    try:
                        trainer.Train(sentenceObj, dcsObj, bz2_input_folder, _debug)
                    except (IndexError, KeyError) as e:
                        print('\x1b[31mFailed: {} \x1b[0m'.format(sentenceObj.sent_id))
                sys.stdout.flush() # Flush IO buffer 
            finishT = time.time()
            print('Avg. time taken by 1 file(1 iteration): {:.3f}'.format((finishT - startT)/(iterationPerBatch*filePerBatch)))
        except KeyboardInterrupt:
            print('Training paused')
            trainer.Save(p_name)
            yield None
    trainer.Save(p_name)
                
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
            if file_counter > 0:
                print('Avg. Micro Recall of Lemmas: {}'.format(np.mean(np.array(recalls))))
                print('Avg. Micro Recall of Words: {}'.format(np.mean(np.array(recalls_of_word))))
                print('Avg. Micro Precision of Lemmas: {}'.format(np.mean(np.array(precisions))))
                print('Avg. Micro Precision of Words: {}'.format(np.mean(np.array(precisions_of_words))))
            sys.stdout.flush() # Flush IO buffer 
        
        file_counter += 1
        
        testFileName = fn.replace('.ds.bz2', '.p2')
        sentenceObj = loaded_SKT[testFileName]
        dcsObj = loaded_DCS[testFileName]    
        
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
    def __init__(self, modelFile = None):
        if modelFile is None:
            self.hidden_layer_size = 1200
            self._edge_vector_dim = 1500
            # self._edge_vector_dim = WD._edge_vector_dim
            # self._full_cnglist = list(WD.mat_cngCount_1D)
            
            self.neuralnet = NN(self._edge_vector_dim, self.hidden_layer_size, outer_relu=True)
            self.history = defaultdict(lambda: list())
        else:
            loader = pickle.load(open(filename, 'rb'))
            
            self.neuralnet.n = loader['n']
            self.neuralnet.d = loader['d']

            self.neuralnet = NN(self._edge_vector_dim, self.hidden_layer_size, outer_relu=True)

            self.neuralnet.U = loader['U']
            self.neuralnet.W = loader['W']
            self.neuralnet.B1 = loader['B1']
            self.neuralnet.B2 = loader['B2']
            
            self.history = defaultdict(lambda: list())
            
        # SET LEARNING RATES
        self.neuralnet.etaW = 3e-4
        self.neuralnet.etaB1 = 1e-4
        
        self.neuralnet.etaU = 1e-4
        self.neuralnet.etaB2 = 1e-4
            
    def Reset(self):
        self.neuralnet = NN(self._edge_vector_dim, self.hidden_layer_size)
        self.history = defaultdict(lambda: list())
        
    def Save(self, filename):
        print('Weights Saved: ', filename)
        pickle.dump({
                'U': self.neuralnet.U,
                'W': self.neuralnet.W,
                'n': self.neuralnet.n,
                'd': self.neuralnet.d,
                'B1': self.neuralnet.B1,
                'B2': self.neuralnet.B2
            }, open(filename, 'wb'))
        return
        
    
    def Load(self, filename):
        loader = pickle.load(open(filename, 'rb'))
        self.neuralnet.U = loader['U']
        self.neuralnet.W = loader['W']
        self.neuralnet.B1 = loader['B1']
        self.neuralnet.B2 = loader['B2']
        self.neuralnet.hidden_layer_size = loader['n']
        self.neuralnet._edge_vector_dim = loader['d']
        
    def Test(self, sentenceObj, dcsObj, dsbz2_name, ecl_mst = False):
        neuralnet = self.neuralnet
        minScore = np.inf
        minMst = None
        
        # dsbz2_name = sentenceObj.sent_id + '.ds.bz2'
        (nodelist_correct, conflicts_Dict_correct, featVMat_correct, nodelist_to_correct_mapping,\
            nodelist, conflicts_Dict, featVMat) = open_dsbz2(dsbz2_name)
        
        # if len(nodelist) > 50:
        #     return None

        if not self.neuralnet.outer_relu:
            (WScalarMat, SigmoidGateOutput) = Get_W_Scalar_Matrix_from_FeatVect_Matrix(featVMat, nodelist, conflicts_Dict, neuralnet)
        else:
            WScalarMat = Get_W_Scalar_Matrix_from_FeatVect_Matrix(featVMat, nodelist, conflicts_Dict, neuralnet)
        
        # print('NeuralNet Time: ', time.time() - startT)
        # startT = time.time()
        
        # Get all MST
        for source in range(len(nodelist)):
            if not ecl_mst:
                (mst_nodes, mst_adj_graph, _) = MST(nodelist, WScalarMat, conflicts_Dict, source)
            else:
                (mst_nodes_bool,mst_nodes,mst_adj_graph)=MST_ECL(nodelist,WScalarMat,conflicts_Dict,source)
            # print('.', end = '')
            score = GetMSTWeight(mst_adj_graph, WScalarMat)
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
    
    def Train(self, sentenceObj, dcsObj, bz2_input_folder, _debug = True):
        # Hyperparameter for hinge loss: m
        m_hinge_param = 14
        
        dsbz2_name = sentenceObj.sent_id + '.ds.bz2'
        (nodelist_correct, conflicts_Dict_correct, featVMat_correct, nodelist_to_correct_mapping,\
            nodelist, conflicts_Dict, featVMat) = open_dsbz2(bz2_input_folder + dsbz2_name)
        # Train for large graphs separately
#         if len(nodelist) < 40:
#             return
        
        """ FORM MAXIMUM(ENERGY) SPANNING TREE OF THE GOLDEN GRAPH : WORST GOLD STRUCTURE """
        WScalarMat_correct = Get_W_Scalar_Matrix_from_FeatVect_Matrix(featVMat_correct, nodelist_correct,\
                                                                      conflicts_Dict_correct, self.neuralnet)
        source = 0
        """ Find the max spanning tree : negative Weight matrix passed """
#         (max_st_gold_ndict, max_st_adj_gold_small, _) =\
#             MST(nodelist_correct, -WScalarMat_correct, conflicts_Dict_correct, source)
        (max_st_gold_ndict, max_st_adj_gold_small, _) =\
            MST(nodelist_correct, -WScalarMat_correct, conflicts_Dict_correct, source)
        energy_gold_max_ST = np.sum(WScalarMat_correct[max_st_adj_gold_small])
        
        """ Convert correct spanning tree graph adj matrix to full marix dimensions """
        """ Create full-size adjacency matrix for correct_mst_small """
        nodelen = len(nodelist)
        max_st_adj_gold = np.ndarray((nodelen, nodelen), np.bool)*False # T_STAR
        for i in range(max_st_adj_gold_small.shape[0]):
            for j in range(max_st_adj_gold_small.shape[1]):
                max_st_adj_gold[nodelist_to_correct_mapping[i], nodelist_to_correct_mapping[j]] =\
                    max_st_adj_gold_small[i, j]
        
        """ Delta(Margin) Function : MASK FOR WHICH NODES IN NODELIST BELONG TO DCS """
        gold_nodes_mask = np.array([False]*len(nodelist))
        gold_nodes_mask[list(nodelist_to_correct_mapping.values())] = True
        margin_f = lambda nodes_mask: np.sum(nodes_mask&(~gold_nodes_mask))**2
        
        """ FOR ALL POSSIBLE MST FROM THE COMPLETE GRAPH """
        WScalarMat = Get_W_Scalar_Matrix_from_FeatVect_Matrix(featVMat, nodelist, conflicts_Dict, self.neuralnet)

        """ For each node - Find MST with that source"""
        min_STx = None # Min Energy spanning tree with worst margin with gold_STx
        min_marginalized_energy = np.inf
        
        # Generate random set of nodes from which mSTs are to be considered
        n_nodes = len(nodelist)
        selection_prob = 0.4
        select_flag = np.random.rand(n_nodes) < selection_prob
        # Fix if all zeros
        if np.sum(select_flag) == 0:
            select_flag[np.random.randint(n_nodes)] = 1
        
        for source in range(len(nodelist)):
            (mst_nodes, mst_adj_graph, mst_nodes_bool) = MST(nodelist, WScalarMat, conflicts_Dict, source) # T_X
            # print('.', end = '')
           
            marginalized_en = np.sum(WScalarMat[mst_adj_graph]) - margin_f(mst_nodes_bool)
            # Minimum marginalized spanning tree : Randomization applied
            # if marginalized_en < min_marginalized_energy and select_flag[source]:
            if marginalized_en < min_marginalized_energy:
                min_marginalized_energy = marginalized_en
                min_STx = mst_adj_graph
            # Energy diff should all be negative
            if _debug:
                print('Source: [{}], Node_Diff:{}, Max_Gold_En: {:.3f}, Energy: {:.3f}'.\
                      format(source, np.sum((~gold_nodes_mask)&mst_nodes_bool), energy_gold_max_ST,  np.sum(WScalarMat[mst_adj_graph])))

        """ Gradient Descent """
        # FOR MOST OFFENdING Y
        doBpp = False
        
        Total_Loss = energy_gold_max_ST - min_marginalized_energy
        if Total_Loss > 0:
            dLdOut = np.zeros_like(WScalarMat)
            dLdOut[max_st_adj_gold&(~min_STx)] = 1
            dLdOut[(~max_st_adj_gold)&min_STx] = -1
            if _debug:
                print('{}. '.format(sentenceObj.sent_id), end = '')
            self.neuralnet.Back_Prop(dLdOut, len(nodelist), featVMat, _debug)
        else:
            trainingStatus[sentenceObj.sent_id] = True
        if _debug:
            print("\nFileKey: %s, Loss: %6.3f" % (sentenceObj.sent_id, Total_Loss))

TrainFiles = None
trainer = None
p_name = ''
odir = ''
def InitModule():
    global trainer
    trainer = Trainer()
    
def register_nnet(nnet, bz2_input_folder):
    if not os.path.isdir(odir):
        os.mkdir(odir)
    if not os.path.isfile('outputs/nnet_LOGS.csv'):
        with open('outputs/nnet_LOGS.csv', 'a') as fh:
            csv_r = csv.writer(fh)
            csv_r.writerow(['odir', 'p_name', 'hidden_layer_size', '_edge_vector_dim'])
    with open('outputs/nnet_LOGS.csv', 'a') as fh:
        csv_r = csv.writer(fh)
        csv_r.writerow([odir, p_name, nnet.n, nnet.d, bz2_input_folder])

"""
################################################################################################
################################################################################################
################################################################################################
"""
def main():
    global TrainFiles, p_name, odir
    """
    ################################################################################################
    ##############################  GET A FILENAME TO SAVE WEIGHTS  ################################
    ################################################################################################
    """
    st = str(int((time.time() * 1e6) % 1e13))
    log_name = 'logs/train_nnet_t{}.out'.format(st)
    odir = 'outputs/train_t{}'.format(st)
    p_name = 'outputs/train_t{}/nnet.p'.format(st)
    print('nEURAL nET wILL bE sAVED hERE: ', p_name)
    
    # Create Training File List
    excluded_files = []
    with open('inputs/Baseline4_advSample.csv', 'r') as f_handle:
        opener = csv.reader(f_handle)
        for line in opener:
            excluded_files.append(line[1].replace('.p', '.ds.bz2'))

    bz2_input_folder = '../NewData/skt_dcs_DS.bz2_4K_pmi_mir_10K/'
    # bz2_input_folder = '/home/rs/15CS91R05/vishnu/Data/skt_dcs_DS.bz2_compat_10k_check_again/'
    all_files = []
    skipped = 0
    for f in os.listdir(bz2_input_folder):
        if '.ds.bz2' in f:
            if f in excluded_files:
                skipped += 1
                continue
            all_files.append(f)

    print(skipped, 'files will not be used for training')
    print('Size of training set:', len(all_files))

    TrainFiles = all_files    
    
    # Load Simultaneous files
    print('Loading Large Files')
    loaded_SKT = pickle.load(open('../Simultaneous_CompatSKT_10K.p', 'rb'), encoding=u'utf-8')
    loaded_DCS = pickle.load(open('../Simultaneous_DCS_10K.p', 'rb'), encoding=u'utf-8')
    
    InitModule()
    trainingStatus = defaultdict(lambda: bool(False))
    
    train = train_generator(loaded_SKT, loaded_DCS, bz2_input_folder, n_trainset = -1, filePerBatch = 10, iterationPerBatch = 5, _debug=False)
    
    # Complete Training
    train.__next__()

    print('Training Complete')
    
if __name__ == '__main__':
    main()
