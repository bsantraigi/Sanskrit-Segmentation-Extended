"""
IMPORTS
"""

import pickle
import numpy as np
from sentences import *
from DCS import *
from collections import defaultdict
import math
from word_definite import *
from romtoslp import *
from nnet import *
from heap_n_PrimMST import *
from word_definite import *

loaded_SKT = pickle.load(open('../Simultaneous_CompatSKT_10K.p', 'rb'))
loaded_DCS = pickle.load(open('../Simultaneous_DCS_10K.p', 'rb'))

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

            while (nodelist2[search_key].lemma != rom_slp(dcsObj.lemmas[chunk_id][j])) or (nodelist2[search_key].cng != int(dcsObj.cng[chunk_id][j])):
                search_key += 1
    #         print((rom_slp(dcsObj.lemmas[chunk_id][j]), dcsObj.cng[chunk_id][j]))
    #         print(nodelist[search_key])
            nodelist2_to_correct_mapping[len(nodelist_correct)] = search_key
            nodelist_correct.append(nodelist2[search_key])
    return (nodelist, nodelist_correct, nodelist2_to_correct_mapping)
    
_edge_vector_dim = len(mat_cng2lem_1D)
_full_cnglist = list(mat_cng2lem_1D)

neuralnet = NN(_edge_vector_dim, 200)

# Train

for iterout in range(10):
    # Change batch
    print('ITERATION OUT', iterout)
    perm = np.random.permutation(len(loaded_SKT))[0:100]
    print('Permutation: ', perm)    
    pickle.dump(neuralnet, open('outputs/neuralnet_trained.p', 'wb'))
    
    # Run few times on same set of files
    for iterin in range(5):
        print('ITERATION IN', iterin)        
        for fn in perm:
            sentenceObj = loaded_SKT[list(loaded_SKT.keys())[fn]]
            dcsObj = loaded_DCS[list(loaded_SKT.keys())[fn]]
            # (chunkDict, lemmaList, wordList, revMap2Chunk, qu, cngList, verbs, tuplesMain, qc_pairs) = SentencePreprocess(sentenceObj)
            try:
                (nodelist, nodelist_correct, nodelist_to_correct_mapping) = GetTrainingKit(sentenceObj, dcsObj)
            except IndexError:
                # print('Error with ', fn)
                continue

            conflicts_Dict = Get_Conflicts(nodelist)
            conflicts_Dict_correct = Get_Conflicts(nodelist_correct)

            featVMat = Get_Feat_Vec_Matrix(nodelist, conflicts_Dict)
            featVMat_correct = Get_Feat_Vec_Matrix(nodelist_correct, conflicts_Dict_correct)

            WScalarMat = Get_W_Scalar_Matrix_from_FeatVect_Matrix(featVMat, nodelist, conflicts_Dict, neuralnet)
            WScalarMat_correct = Get_W_Scalar_Matrix_from_FeatVect_Matrix(featVMat_correct, nodelist_correct, conflicts_Dict_correct, neuralnet)
            Total_Loss = 0

            # Get MST for the correct nodelist
            source = 0
            (mst_nodes_correct, mst_adj_graph_correct_0) = MST(nodelist_correct, WScalarMat_correct, conflicts_Dict_correct, source)
        #     print('Correct MST Score: ', GetMSTWeight(mst_nodes_correct, WScalarMat_correct))
            W_star = GetMSTWeight(mst_nodes_correct, WScalarMat_correct)

            # Convert correct spanning tree graph adj matrix to actual marix dimensions
            nodelen = len(nodelist)
            mst_adj_graph_correct = np.ndarray((nodelen, nodelen), np.bool)*False
            for i in range(mst_adj_graph_correct_0.shape[0]):
                for j in range(mst_adj_graph_correct_0.shape[1]):
                    mst_adj_graph_correct[nodelist_to_correct_mapping[i], nodelist_to_correct_mapping[j]] = \
                    mst_adj_graph_correct_0[i, j]

            # Get all MST
            dLdS = np.zeros(WScalarMat.shape)
            for source in range(len(nodelist)):
            #     print('\nSOURCE: {}'.format(source))
                (mst_nodes, mst_adj_graph) = MST(nodelist, WScalarMat, conflicts_Dict, source)    
            #     print('Correct MST Score for Souce {}: {}'.format(source, GetMSTWeight(mst_nodes, WScalarMat)))
                print('.', end = '')
                Total_Loss += (W_star - GetMSTWeight(mst_nodes, WScalarMat))
                dLdS_inner = 1 / WScalarMat
                dLdS_inner[mst_adj_graph_correct == 1] *= -1
                dLdS_inner = dLdS_inner*(mst_adj_graph^mst_adj_graph_correct)
                dLdS += dLdS_inner        
            neuralnet.Back_Prop(dLdS/len(nodelist), len(nodelist), featVMat)

            Total_Loss /= len(nodelist)
            print("\nFileKey: %5d, Loss: %6.3f, Original MSTScore: %6.3f" % (fn, Total_Loss, W_star))
