from sentences import *
from DCS import *
from romtoslp import *
from collections import defaultdict
import numpy as np

class word_definite:
    def __init__(self,derived, lemma, cng, pos, chunk_id):
        self.lemma = lemma
        self.derived = derived
        self.cng = str(cng)
        self.tup = "{}_{}".format(self.lemma, self.cng)
        # self.form = form
        self.pos = pos
        self.chunk_id = chunk_id        
        # Fields Required for heap
        self.dist = np.inf
        self.src = -1
        self.id = -1
        self.isConflicted = False
    def __str__(self):
        return 'WD_Node[C: %d, P: %d, %s @(%s) => %s]' %(self.chunk_id, self.pos, self.lemma, self.cng, self.derived)
    def __repr__(self):
        return str(self)
        
def GetNodes(sentenceObj):
    def getCNGs(formsDict):
            if type(formsDict) == int or type(formsDict) == str:
                return [int(formsDict)]
            else:
                l = []
                for form, configs in formsDict.items():
                    for c in configs:
                        if(form == 'verbform'):
                            continue
                        else:
                            l.append(wtc_recursive(form, c))
                return list(set(l))
    nodeCount = 0
    for chunk in sentenceObj.chunk:
        for word_pos, words in chunk.chunk_words.items():
            for word in words:
                if len(word.forms) == len(word.lemmas):
                    for lemma_i in range(len(word.lemmas)):                    
                        cng_list = getCNGs(word.forms[lemma_i])
                        for cng_x in cng_list:
                            nodeCount += 1

    nodelist = [None]*nodeCount
    chunk_id = 0
    node_id = 0
    # chunk_dict = {}
    for chunk in sentenceObj.chunk:
        # chunk_dict[str(chunk_id)] = {}
        for word_pos, words in chunk.chunk_words.items():
            # chunk_dict[str(chunk_id)][str(word_pos)] = [];
            for word in words:
                if len(word.forms) == len(word.lemmas):
                    for lemma_i in range(len(word.lemmas)):
                        cng_list = getCNGs(word.forms[lemma_i])
                        for x in range(len(cng_list)):
                            cng_x = cng_list[x]
                            nodelist[node_id] = word_definite(
                                rom_slp(word.names), 
                                rom_slp(word.lemmas[lemma_i]),
                                cng_x, word_pos, chunk_id
                            )                        
                            # chunk_dict[str(chunk_id)][str(word_pos)].append(node_id)
                            node_id += 1

        chunk_id += 1
    
    # return nodelist, chunk_dict
    return nodelist
	
def Get_Conflicts(nodelist_new):
    nodesCount = len(nodelist_new)
    conflicts_Dict = defaultdict(lambda: [])
#     conflicts = np.ndarray((n, n), dtype = bool)
#     conflicts[:,:] = False
    for i in range(nodesCount):
        conflicts_Dict[i]; # Just to add a dummy entry
        chunk_id_1 = nodelist_new[i].chunk_id
        for j in range(i + 1, nodesCount):
            chunk_id_2 = nodelist_new[j].chunk_id
            # Ofcourse the two lemmas have to be in same chunk
            if(chunk_id_2 > chunk_id_1):
                # Not same chunk - get out of here
                break
            elif(chunk_id_2 == chunk_id_1):
                # Let's see then
                # Required: wd1 is on the left of wd2 in the chunk
                if(nodelist_new[i].pos < nodelist_new[j].pos):
                    wd1 = nodelist_new[i]
                    wd2 = nodelist_new[j]
                else:
                    wd2 = nodelist_new[i]
                    wd1 = nodelist_new[j]
                # Judge: Sandhi
                if not CanCoExist_sandhi(wd1.pos, wd2.pos, wd1.derived, wd2.derived):
                    # I am sorry guys - you two are not compatible
                    conflicts_Dict[i].append(j)
                    conflicts_Dict[j].append(i)
    return dict(conflicts_Dict)
   

"""
################################################################################################
##################### LOAD PROBABILITY MATRICES  ## LEMMA, CNG, CNG_GROUP # ####################
##############################  TUPLES ALSO  ###################################################
"""
mat_lem2lem_countonly = None
mat_lem2cng_countonly = None
mat_lem2tup_countonly = None

mat_cng2lem_countonly = None
mat_cng2cng_countonly = None
mat_cng2tup_countonly = None

mat_tup2lem_countonly = None
mat_tup2cng_countonly = None
mat_tup2tup_countonly = None

mat_lemCount_1D = None
mat_cngCount_1D = None
mat_tupCount_1D = None

_cg_count = None
_edge_vector_dim = None
_full_cnglist = None

def word_definite_extInit(matDB):
    global mat_lem2lem_countonly, mat_lem2cng_countonly, mat_lem2tup_countonly, \
     mat_cng2lem_countonly, mat_cng2cng_countonly, mat_cng2tup_countonly, \
     mat_tup2lem_countonly, mat_tup2cng_countonly, mat_tup2tup_countonly, \
     mat_lemCount_1D, mat_cngCount_1D, mat_tupCount_1D, \
     _cg_count, _full_cnglist, _edge_vector_dim

    mat_lem2lem_countonly = matDB.mat_lem2lem_countonly
    mat_lem2cng_countonly = matDB.mat_lem2cng_countonly
    mat_lem2tup_countonly = matDB.mat_lem2tup_countonly

    mat_cng2lem_countonly = matDB.mat_cng2lem_countonly
    mat_cng2cng_countonly = matDB.mat_cng2cng_countonly
    mat_cng2tup_countonly = matDB.mat_cng2tup_countonly
    
    mat_tup2lem_countonly = matDB.mat_tup2lem_countonly
    mat_tup2cng_countonly = matDB.mat_tup2cng_countonly
    mat_tup2tup_countonly = matDB.mat_tup2tup_countonly

    mat_lemCount_1D = matDB.mat_lemCount_1D
    mat_cngCount_1D = matDB.mat_cngCount_1D
    mat_tupCount_1D = matDB.mat_tupCount_1D

    # TODO: Change to actual value
    _cg_count = len(mat_cngCount_1D)
    # _edge_vector_dim = 9*_cg_count**2 + 9 * _cg_count + 8
    _edge_vector_dim = 4*_cg_count**2 + 9*_cg_count + 9
    _full_cnglist = list(mat_cngCount_1D)
    print(_cg_count)
"""
################################################################################################
######################  CREATING FEATURE MATRICES ##############################################
################################################################################################
"""
def tryProb_catchZero(mat2D, mat1D, key1, key2):
    try:
        v = float(mat2D[key1][key2])/mat1D[key1];
    except KeyError:
        v = 0
    return v
def Get_Features(node1, node2):
    feats = np.zeros((_edge_vector_dim, 1))
    # print('For ', node1, node2)
    

    # Path Constraint - Length 1 - # 8
    # tup2tup edge ignored
    fIndex = 0
    feats[fIndex] = tryProb_catchZero(mat_lem2lem_countonly, mat_lemCount_1D, node1.lemma, node2.lemma); fIndex += 1;
    feats[fIndex] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, node2.cng); fIndex += 1;
    feats[fIndex] = tryProb_catchZero(mat_lem2tup_countonly, mat_lemCount_1D, node1.lemma, node2.tup); fIndex += 1;

    feats[fIndex] = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, node1.cng, node2.lemma); fIndex += 1;
    feats[fIndex] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, node2.cng); fIndex += 1;
    feats[fIndex] = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, node1.cng, node2.tup); fIndex += 1;

    feats[fIndex] = tryProb_catchZero(mat_tup2lem_countonly, mat_tupCount_1D, node1.tup, node2.lemma); fIndex += 1;
    feats[fIndex] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, node2.cng); fIndex += 1;
    feats[fIndex] = tryProb_catchZero(mat_tup2tup_countonly, mat_tupCount_1D, node1.tup, node2.tup); fIndex += 1;

    # Path Constraint - Length 2 - # _cg_count

    # LEMMA->CNG->LEMMA
    for k in range(0, _cg_count):
        cng_k = _full_cnglist[k]
        # TODO: Some lemma's still missing
        pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, cng_k)
        pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, cng_k, node2.lemma)
        feats[fIndex + k] = pleft * pright
    fIndex += _cg_count

    # LEMMA->CNG->CNG
    for k in range(0, _cg_count):
        cng_k = _full_cnglist[k]
        pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, cng_k)
        pright = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, cng_k, node2.cng)
        feats[fIndex + k] = pleft * pright
    fIndex += _cg_count

    # LEMMA->CNG->TUP
    for k in range(0, _cg_count):
        cng_k = _full_cnglist[k]
        # TODO: Some lemma's still missing
        pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, cng_k)
        pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, cng_k, node2.tup)
        feats[fIndex + k] = pleft * pright
    fIndex += _cg_count

    # CNG->CNG->LEMMA
    for k in range(0, _cg_count):
        cng_k = _full_cnglist[k]
        pleft = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, cng_k)
        pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, cng_k, node2.lemma)
        feats[fIndex + k] = pleft * pright
    fIndex += _cg_count

    # CNG->CNG->CNG
    for k in range(0, _cg_count):
        cng_k = _full_cnglist[k]
        pleft = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, cng_k)
        pright = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, cng_k, node2.cng)
        feats[fIndex + k] = pleft * pright
    fIndex += _cg_count

    # CNG->CNG->TUP
    for k in range(0, _cg_count):
        cng_k = _full_cnglist[k]
        pleft = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, cng_k)
        pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, cng_k, node2.tup)
        feats[fIndex + k] = pleft * pright
    fIndex += _cg_count

    # TUP->CNG->LEMMA :: TOO MANY ZEROS
    for k in range(0, _cg_count):
        cng_k = _full_cnglist[k]
        # TODO: Some lemma's still missing
        pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, cng_k)
        pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, cng_k, node2.lemma)
        feats[fIndex + k] = pleft * pright
    fIndex += _cg_count

    # TUP->CNG->CNG :: TOO MANY ZEROS
    for k in range(0, _cg_count):
        cng_k = _full_cnglist[k]
        # TODO: Some lemma's still missing
        pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, cng_k)
        pright = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, cng_k, node2.cng)
        feats[fIndex + k] = pleft * pright
    fIndex += _cg_count

    # TUP->CNG->TUP :: TOO MANY ZEROS
    for k in range(0, _cg_count):
        cng_k = _full_cnglist[k]
        # TODO: Some lemma's still missing
        pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, cng_k)
        pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, cng_k, node2.tup)
        feats[fIndex + k] = pleft * pright
    fIndex += _cg_count

    # Path Constraint - Length 3 - # _cg_count^2

    # LEMMA->CGS->CGS->LEMMA
    for k1 in range(0, _cg_count):
        cng_k1 = _full_cnglist[k1]
        pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, cng_k1)
        for k2 in range(0, _cg_count): 
            cng_k2 = _full_cnglist[k2]
            # TODO: Some lemma's still missing
            pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, cng_k1, cng_k2)
            pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, cng_k2, node2.lemma)
            feats[fIndex + k1*_cg_count + k2] = pleft  * pmid * pright
    fIndex += _cg_count**2

    # # LEMMA->CGS->CGS->CNG
    # for k1 in range(0, _cg_count):
    #     cng_k1 = _full_cnglist[k1]
    #     pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, cng_k1)
    #     for k2 in range(0, _cg_count): 
    #         cng_k2 = _full_cnglist[k2]
    #         # TODO: Some lemma's still missing
    #         pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, cng_k1, cng_k2)
    #         pright = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, cng_k2, node2.cng)
    #         feats[fIndex + k1*_cg_count + k2] = pleft  * pmid * pright
    # fIndex += _cg_count**2

    # LEMMA->CGS->CGS->TUP
    for k1 in range(0, _cg_count):
        cng_k1 = _full_cnglist[k1]
        pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, cng_k1)
        for k2 in range(0, _cg_count): 
            cng_k2 = _full_cnglist[k2]
            # TODO: Some lemma's still missing
            pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, cng_k1, cng_k2)
            pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, cng_k2, node2.tup)
            feats[fIndex + k1*_cg_count + k2] = pleft  * pmid * pright
    fIndex += _cg_count**2

    # # CNG->CGS->CGS->LEMMA
    # for k1 in range(0, _cg_count):
    #     cng_k1 = _full_cnglist[k1]
    #     pleft = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, cng_k1)
    #     for k2 in range(0, _cg_count): 
    #         cng_k2 = _full_cnglist[k2]
    #         # TODO: Some lemma's still missing
    #         pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, cng_k1, cng_k2)
    #         pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, cng_k2, node2.lemma)
    #         feats[fIndex + k1*_cg_count + k2] = pleft  * pmid * pright
    # fIndex += _cg_count**2

    # # CNG->CGS->CGS->CNG
    # for k1 in range(0, _cg_count):
    #     cng_k1 = _full_cnglist[k1]
    #     pleft = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, cng_k1)
    #     for k2 in range(0, _cg_count): 
    #         cng_k2 = _full_cnglist[k2]
    #         # TODO: Some lemma's still missing
    #         pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, cng_k1, cng_k2)
    #         pright = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, cng_k2, node2.cng)
    #         feats[fIndex + k1*_cg_count + k2] = pleft  * pmid * pright
    # fIndex += _cg_count**2

    # # CNG->CGS->CGS->TUP
    # for k1 in range(0, _cg_count):
    #     cng_k1 = _full_cnglist[k1]
    #     pleft = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, cng_k1)
    #     for k2 in range(0, _cg_count): 
    #         cng_k2 = _full_cnglist[k2]
    #         # TODO: Some lemma's still missing
    #         pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, cng_k1, cng_k2)
    #         pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, cng_k2, node2.tup)
    #         feats[fIndex + k1*_cg_count + k2] = pleft  * pmid * pright
    # fIndex += _cg_count**2

    # TUP->CGS->CGS->LEM
    for k1 in range(0, _cg_count):
        cng_k1 = _full_cnglist[k1]
        pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, cng_k1)
        for k2 in range(0, _cg_count): 
            cng_k2 = _full_cnglist[k2]
            # TODO: Some lemma's still missing
            pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, cng_k1, cng_k2)
            pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, cng_k2, node2.lemma)
            feats[fIndex + k1*_cg_count + k2] = pleft  * pmid * pright
    fIndex += _cg_count**2

    # # TUP->CGS->CGS->CNG
    # for k1 in range(0, _cg_count):
    #     cng_k1 = _full_cnglist[k1]
    #     pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, cng_k1)
    #     for k2 in range(0, _cg_count): 
    #         cng_k2 = _full_cnglist[k2]
    #         # TODO: Some lemma's still missing
    #         pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, cng_k1, cng_k2)
    #         pright = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, cng_k2, node2.cng)
    #         feats[fIndex + k1*_cg_count + k2] = pleft  * pmid * pright
    # fIndex += _cg_count**2

    # TUP->CGS->CGS->TUP
    for k1 in range(0, _cg_count):
        cng_k1 = _full_cnglist[k1]
        pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, cng_k1)
        for k2 in range(0, _cg_count): 
            cng_k2 = _full_cnglist[k2]
            # TODO: Some lemma's still missing
            pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, cng_k1, cng_k2)
            pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, cng_k2, node2.tup)
            feats[fIndex + k1*_cg_count + k2] = pleft  * pmid * pright
    fIndex += _cg_count**2


    return feats
def Get_Feat_Vec_Matrix(nodelist_new, conflicts_Dict):
    nodesCount = len(nodelist_new)
    featVMat = [[None for _ in range(nodesCount)] for _ in range(nodesCount)]
    for i in range(nodesCount):
        for j in range(nodesCount):
            if j in conflicts_Dict[i] or i == j:                
                # featVMat[i][j][:] = 1e-35
                pass
            else:
                featVMat[i][j] = Get_Features(nodelist_new[i], nodelist_new[j])
    return featVMat

# This function will be some neural network or a linear function or something of sorts
def Get_W_Scalar_Matrix_from_FeatVect_Matrix(featVMat, nodelist_new, conflicts_Dict, _neuralnet):
    nodesCount = len(nodelist_new)
    WScalarMat = np.zeros((nodesCount, nodesCount))   
    sigmoid_25 = 1/(1+np.exp(-25))
    for i in range(nodesCount):
        for j in range(nodesCount):
            if featVMat[i][j] is None:
                pass
            else:
                # Since s is output of a sigmoid gate, it will always be greater than zero
                (_, _, s) = _neuralnet.Forward_Prop(featVMat[i][j])
                WScalarMat[i, j] = np.minimum(s, sigmoid_25)
    toinf = (WScalarMat == 0)
    WScalarMat[WScalarMat > 0] = -np.log2(WScalarMat[WScalarMat > 0])
    WScalarMat[toinf] = np.inf
    return WScalarMat
