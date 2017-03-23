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
    # _edge_vector_dim = 4*_cg_count**2 + 9*_cg_count + 9
    _edge_vector_dim = 1000
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
'''
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
'''
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
    if not _neuralnet.outer_relu:        
        sigmoid_25 = 1/(1+np.exp(-25)) # Upper limit for sigmoid gate to prevent overflow
        sigmoid_n25 = 1/(1+np.exp(25)) # Lower limit for sigmoid gate to prevent underflow
        for i in range(nodesCount):
            for j in range(nodesCount):
                if featVMat[i][j] is None:
                    pass
                else:
                    # Since s is output of a sigmoid gate, it will always be greater than zero
                    # We apply sigmoid as we require a probability value
                    (_, _, s) = _neuralnet.Forward_Prop(featVMat[i][j])
                    WScalarMat[i, j] = np.minimum(s, sigmoid_25) # Cap on the highest value of sigmoid game
        
        # Apply -log2(.) : We want to minimize the loss function rather that maximizing probability
        SigmoidGateOutput = WScalarMat.copy()
        toinf = (WScalarMat == 0)
        WScalarMat[WScalarMat > 0] = -np.log2(WScalarMat[WScalarMat > 0])
        WScalarMat[toinf] = np.inf
        return (WScalarMat, SigmoidGateOutput)
    else:
        for i in range(nodesCount):
            for j in range(nodesCount):
                if featVMat[i][j] is None:
                    pass
                else:
                    # Since s is output of a relu gate, it will always be greater than zero
                    (_, _, s) = _neuralnet.Forward_Prop(featVMat[i][j])
                    WScalarMat[i, j] = s   
        return WScalarMat



########################## ARTIFICIALLY GENERATED CODE FOR FEATURE GENERATION ###########################
############################# FOLLOWING FEATURES WERE SELECTED BASED ON A ###############################
############################# MUTUAL INFORMATION BASED SELECTION CRITERIA ###############################


def Get_Features(node1, node2):
    feats = np.zeros((_edge_vector_dim, 1))
    fIndex = 0
    
    # L->dat. du.->L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pright; fIndex += 1
            
    # L->-90->L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-90')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-90', node2.lemma)
    feats[fIndex] = pleft * pright; fIndex += 1
            
    # L->-57->L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-57')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-57', node2.lemma)
    feats[fIndex] = pleft * pright; fIndex += 1
            
    # L->8_fp->L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pright; fIndex += 1
            
    # L->169->L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '169')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '169', node2.lemma)
    feats[fIndex] = pleft * pright; fIndex += 1
            
    # L->54->L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pright; fIndex += 1
            
    # L->153->L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pright; fIndex += 1
            
    # L->-24->L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-24')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    feats[fIndex] = pleft * pright; fIndex += 1
            
    # L->137->L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pright; fIndex += 1
            
    # L->-112->L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pright; fIndex += 1
            
    # L->101->L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pright; fIndex += 1
            
    # L->91->L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '91')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    feats[fIndex] = pleft * pright; fIndex += 1
            
    # L->3_du->L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pright; fIndex += 1
            
    # L->91->T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '91')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '91', node2.tup)
    feats[fIndex] = pleft * pright; fIndex += 1
            
    # C->3_du->T
    pleft = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '3_du')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '3_du', node2.tup)
    feats[fIndex] = pleft * pright; fIndex += 1
            
    # T->8_fp->L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pright; fIndex += 1
            
    # T->54->L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pright; fIndex += 1
            
    # T->153->L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pright; fIndex += 1
            
    # T->137->L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pright; fIndex += 1
            
    # T->-112->L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pright; fIndex += 1
            
    # T->101->L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pright; fIndex += 1
            
    # T->91->L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '91')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    feats[fIndex] = pleft * pright; fIndex += 1
            
    # L -> 2_sp -> dat. du. -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_sp', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 2_sp -> 8_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_sp', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 2_sp -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_sp', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 2_sp -> pl -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_sp', 'pl')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 2_sp -> 28 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_sp', '28')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 2_sp -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_sp', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 2_sp -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_sp', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 2_sp -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_sp', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 5_du -> dat. du. -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_du', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 5_du -> -57 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_du', '-57')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-57', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 5_du -> 8_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_du', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 5_du -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_du', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 5_du -> pl -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_du', 'pl')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 5_du -> 153 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_du', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 5_du -> 28 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_du', '28')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 5_du -> -24 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_du', '-24')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 5_du -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_du', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 5_du -> -112 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_du', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 5_du -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_du', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 5_du -> 91 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_du', '91')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 5_du -> -35 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_du', '-35')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 5_du -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_du', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> gen. sg. -> dat. du. -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. sg.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. sg.', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> gen. sg. -> -90 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. sg.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. sg.', '-90')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-90', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> gen. sg. -> -57 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. sg.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. sg.', '-57')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-57', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> gen. sg. -> 8_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. sg.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. sg.', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> gen. sg. -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. sg.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. sg.', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> gen. sg. -> 153 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. sg.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. sg.', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> gen. sg. -> -24 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. sg.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. sg.', '-24')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> gen. sg. -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. sg.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. sg.', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> gen. sg. -> -112 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. sg.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. sg.', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> gen. sg. -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. sg.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. sg.', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> gen. sg. -> 91 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. sg.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. sg.', '91')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> gen. sg. -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. sg.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. sg.', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 15_sp -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_sp', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 5_sp -> dat. du. -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 5_sp -> -151 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', '-151')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-151', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 5_sp -> -57 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', '-57')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-57', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 5_sp -> 8_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 5_sp -> 169 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', '169')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '169', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 5_sp -> pl -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', 'pl')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 5_sp -> 153 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 5_sp -> -24 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', '-24')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 5_sp -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 5_sp -> -112 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 5_sp -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 5_sp -> 91 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', '91')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 5_sp -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> dat. du. -> 5_sp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. du.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. du.', '5_sp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> dat. du. -> dat. du. -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. du.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. du.', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> dat. du. -> -151 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. du.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. du.', '-151')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-151', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> dat. du. -> -90 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. du.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. du.', '-90')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-90', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> dat. du. -> 8_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. du.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. du.', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> dat. du. -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. du.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. du.', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> dat. du. -> pl -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. du.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. du.', 'pl')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> dat. du. -> 153 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. du.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. du.', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> dat. du. -> -24 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. du.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. du.', '-24')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> dat. du. -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. du.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. du.', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> dat. du. -> -112 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. du.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. du.', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> dat. du. -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. du.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. du.', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> dat. du. -> -35 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. du.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. du.', '-35')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> dat. du. -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. du.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. du.', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -37 -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-37')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-37', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -37 -> 153 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-37')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-37', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -37 -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-37')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-37', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -37 -> -112 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-37')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-37', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -37 -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-37')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-37', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 15_du -> dat. du. -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_du', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 15_du -> -151 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_du', '-151')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-151', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 15_du -> -90 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_du', '-90')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-90', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 15_du -> -57 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_du', '-57')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-57', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 15_du -> 8_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_du', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 15_du -> 169 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_du', '169')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '169', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 15_du -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_du', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 15_du -> pl -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_du', 'pl')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 15_du -> 153 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_du', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 15_du -> 28 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_du', '28')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 15_du -> -24 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_du', '-24')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 15_du -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_du', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 15_du -> -112 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_du', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 15_du -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_du', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 15_du -> 91 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_du', '91')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 15_du -> -35 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_du', '-35')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 15_du -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_du', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -151 -> dat. du. -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-151')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-151', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -151 -> -90 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-151')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-151', '-90')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-90', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -151 -> 8_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-151')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-151', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -151 -> 169 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-151')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-151', '169')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '169', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -151 -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-151')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-151', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -151 -> pl -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-151')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-151', 'pl')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -151 -> 153 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-151')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-151', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -151 -> 28 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-151')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-151', '28')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -151 -> -24 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-151')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-151', '-24')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -151 -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-151')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-151', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -151 -> -112 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-151')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-151', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -151 -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-151')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-151', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -151 -> 91 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-151')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-151', '91')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -151 -> -35 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-151')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-151', '-35')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -151 -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-151')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-151', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 39 -> dat. du. -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 39 -> 8_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 39 -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 39 -> 153 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 39 -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 39 -> -112 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 39 -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 39 -> 91 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '91')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 39 -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 12_fp -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_fp', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -306 -> 5_sp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '5_sp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -306 -> dat. du. -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -306 -> -151 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '-151')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-151', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -306 -> -90 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '-90')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-90', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -306 -> -57 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '-57')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-57', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -306 -> 8_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -306 -> 169 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '169')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '169', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -306 -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -306 -> pl -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', 'pl')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -306 -> 153 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -306 -> 28 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '28')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -306 -> -24 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '-24')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -306 -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -306 -> -112 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -306 -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -306 -> 91 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '91')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -306 -> -35 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '-35')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -306 -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -64 -> dat. du. -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-64')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-64', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -64 -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-64')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-64', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -64 -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-64')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-64', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -64 -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-64')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-64', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -64 -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-64')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-64', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -90 -> -151 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-90')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-90', '-151')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-151', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -90 -> -90 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-90')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-90', '-90')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-90', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -90 -> 8_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-90')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-90', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -90 -> 169 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-90')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-90', '169')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '169', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -90 -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-90')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-90', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -90 -> pl -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-90')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-90', 'pl')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -90 -> 153 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-90')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-90', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -90 -> 28 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-90')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-90', '28')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -90 -> -24 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-90')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-90', '-24')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -90 -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-90')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-90', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -90 -> -112 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-90')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-90', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -90 -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-90')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-90', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -90 -> 91 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-90')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-90', '91')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -90 -> -35 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-90')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-90', '-35')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -90 -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-90')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-90', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -210 -> 8_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -210 -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -210 -> 153 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -210 -> -112 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -210 -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -210 -> 91 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '91')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -210 -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -159 -> dat. du. -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-159')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-159', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -159 -> 8_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-159')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-159', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -159 -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-159')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-159', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -159 -> pl -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-159')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-159', 'pl')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -159 -> 153 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-159')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-159', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -159 -> -24 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-159')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-159', '-24')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -159 -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-159')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-159', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -159 -> -112 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-159')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-159', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -159 -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-159')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-159', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -159 -> 91 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-159')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-159', '91')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -159 -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-159')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-159', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -57 -> dat. du. -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-57')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-57', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -57 -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-57')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-57', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -57 -> pl -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-57')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-57', 'pl')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -57 -> 153 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-57')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-57', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -57 -> -24 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-57')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-57', '-24')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -57 -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-57')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-57', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -57 -> -112 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-57')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-57', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -57 -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-57')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-57', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -57 -> 91 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-57')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-57', '91')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -57 -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-57')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-57', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 10_sp -> 5_sp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sp', '5_sp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 10_sp -> dat. du. -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sp', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 10_sp -> 8_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sp', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 10_sp -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sp', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 10_sp -> pl -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sp', 'pl')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 10_sp -> 153 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sp', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 10_sp -> 28 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sp', '28')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 10_sp -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sp', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 10_sp -> -112 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sp', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 10_sp -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sp', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 10_sp -> -35 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sp', '-35')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 10_sp -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sp', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 6_sg -> dat. du. -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 6_sg -> 8_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 6_sg -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 6_sg -> pl -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', 'pl')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 6_sg -> 28 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', '28')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 6_sg -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 6_sg -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 6_sg -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 8_fp -> dat. du. -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_fp', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 8_fp -> -151 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_fp', '-151')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-151', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 8_fp -> -57 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_fp', '-57')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-57', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 8_fp -> 8_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_fp', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 8_fp -> 169 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_fp', '169')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '169', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 8_fp -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_fp', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 8_fp -> pl -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_fp', 'pl')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 8_fp -> 153 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_fp', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 8_fp -> 28 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_fp', '28')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 8_fp -> -24 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_fp', '-24')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 8_fp -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_fp', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 8_fp -> -112 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_fp', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 8_fp -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_fp', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 8_fp -> 91 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_fp', '91')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 8_fp -> -35 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_fp', '-35')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 8_fp -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_fp', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 4_pl -> dat. du. -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_pl')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_pl', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 4_pl -> 8_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_pl')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_pl', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 4_pl -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_pl')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_pl', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 4_pl -> pl -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_pl')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_pl', 'pl')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 4_pl -> 153 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_pl')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_pl', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 4_pl -> 28 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_pl')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_pl', '28')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 4_pl -> -24 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_pl')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_pl', '-24')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 4_pl -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_pl')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_pl', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 4_pl -> -112 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_pl')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_pl', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 4_pl -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_pl')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_pl', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 4_pl -> 91 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_pl')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_pl', '91')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 4_pl -> -35 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_pl')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_pl', '-35')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 4_pl -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_pl')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_pl', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 169 -> dat. du. -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '169')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '169', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 169 -> -151 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '169')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '169', '-151')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-151', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 169 -> -90 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '169')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '169', '-90')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-90', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 169 -> -57 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '169')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '169', '-57')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-57', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 169 -> 8_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '169')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '169', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 169 -> 169 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '169')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '169', '169')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '169', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 169 -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '169')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '169', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 169 -> 153 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '169')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '169', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 169 -> 28 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '169')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '169', '28')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 169 -> -24 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '169')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '169', '-24')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 169 -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '169')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '169', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 169 -> -112 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '169')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '169', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 169 -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '169')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '169', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 169 -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '169')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '169', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 54 -> 5_sp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '54')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '54', '5_sp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 54 -> dat. du. -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '54')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '54', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 54 -> -151 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '54')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '54', '-151')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-151', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 54 -> -90 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '54')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '54', '-90')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-90', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 54 -> -57 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '54')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '54', '-57')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-57', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 54 -> 8_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '54')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '54', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 54 -> 169 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '54')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '54', '169')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '169', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 54 -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '54')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '54', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 54 -> 114 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '54')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '54', '114')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '114', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 54 -> pl -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '54')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '54', 'pl')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 54 -> 153 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '54')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '54', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 54 -> -24 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '54')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '54', '-24')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 54 -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '54')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '54', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 54 -> -112 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '54')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '54', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 54 -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '54')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '54', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 54 -> 91 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '54')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '54', '91')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 54 -> -35 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '54')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '54', '-35')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 54 -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '54')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '54', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 114 -> dat. du. -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '114')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '114', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 114 -> -151 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '114')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '114', '-151')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-151', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 114 -> -90 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '114')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '114', '-90')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-90', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 114 -> 8_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '114')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '114', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 114 -> 169 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '114')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '114', '169')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '169', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 114 -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '114')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '114', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 114 -> pl -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '114')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '114', 'pl')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 114 -> 153 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '114')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '114', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 114 -> -24 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '114')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '114', '-24')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 114 -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '114')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '114', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 114 -> -112 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '114')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '114', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 114 -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '114')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '114', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 114 -> -35 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '114')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '114', '-35')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 114 -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '114')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '114', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> loc. sg. -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -139 -> 5_sp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '5_sp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -139 -> dat. du. -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -139 -> -151 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-151')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-151', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -139 -> -90 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-90')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-90', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -139 -> -57 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-57')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-57', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -139 -> 8_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -139 -> 169 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '169')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '169', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -139 -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -139 -> pl -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', 'pl')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -139 -> 153 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -139 -> 28 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '28')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -139 -> -24 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-24')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -139 -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -139 -> -112 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -139 -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -139 -> 91 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '91')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -139 -> -35 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-35')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -139 -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 8_du -> dat. du. -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_du', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 8_du -> 8_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_du', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 8_du -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_du', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 8_du -> pl -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_du', 'pl')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 8_du -> 28 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_du', '28')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 8_du -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_du', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 8_du -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_du', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 8_du -> -35 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_du', '-35')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 8_du -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_du', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -166 -> dat. du. -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-166')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-166', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -166 -> 8_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-166')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-166', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -166 -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-166')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-166', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -166 -> pl -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-166')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-166', 'pl')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -166 -> 153 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-166')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-166', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -166 -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-166')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-166', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -166 -> -112 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-166')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-166', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -166 -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-166')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-166', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -166 -> -35 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-166')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-166', '-35')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -166 -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-166')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-166', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 139 -> 5_sp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '139')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '139', '5_sp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 139 -> dat. du. -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '139')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '139', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 139 -> -90 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '139')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '139', '-90')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-90', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 139 -> -57 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '139')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '139', '-57')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-57', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 139 -> 8_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '139')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '139', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 139 -> 169 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '139')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '139', '169')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '169', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 139 -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '139')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '139', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 139 -> pl -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '139')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '139', 'pl')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 139 -> 153 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '139')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '139', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 139 -> 135 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '139')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '139', '135')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '135', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 139 -> 28 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '139')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '139', '28')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 139 -> -24 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '139')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '139', '-24')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 139 -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '139')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '139', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 139 -> -112 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '139')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '139', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 139 -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '139')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '139', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 139 -> 91 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '139')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '139', '91')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 139 -> -35 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '139')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '139', '-35')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 139 -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '139')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '139', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 92 -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '92')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '92', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> pl -> dat. du. -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> pl -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> pl -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> pl -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -91 -> dat. du. -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-91')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-91', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -91 -> -151 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-91')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-91', '-151')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-151', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -91 -> -90 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-91')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-91', '-90')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-90', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -91 -> 8_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-91')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-91', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -91 -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-91')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-91', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -91 -> pl -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-91')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-91', 'pl')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -91 -> 153 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-91')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-91', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -91 -> 28 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-91')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-91', '28')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -91 -> -24 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-91')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-91', '-24')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -91 -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-91')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-91', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -91 -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-91')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-91', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -91 -> 91 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-91')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-91', '91')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -91 -> -35 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-91')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-91', '-35')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -91 -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-91')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-91', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 28_sg -> 5_sp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_sg')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28_sg', '5_sp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 28_sg -> dat. du. -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_sg')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28_sg', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 28_sg -> 8_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_sg')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28_sg', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 28_sg -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_sg')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28_sg', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 28_sg -> pl -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_sg')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28_sg', 'pl')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 28_sg -> 153 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_sg')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28_sg', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 28_sg -> 28 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_sg')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28_sg', '28')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 28_sg -> -24 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_sg')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28_sg', '-24')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 28_sg -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_sg')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28_sg', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 28_sg -> -112 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_sg')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28_sg', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 28_sg -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_sg')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28_sg', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 28_sg -> 91 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_sg')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28_sg', '91')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 28_sg -> -35 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_sg')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28_sg', '-35')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 28_sg -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_sg')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28_sg', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 148 -> 5_sp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '148')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '148', '5_sp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 148 -> dat. du. -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '148')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '148', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 148 -> 8_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '148')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '148', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 148 -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '148')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '148', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 148 -> pl -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '148')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '148', 'pl')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 148 -> 153 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '148')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '148', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 148 -> 28 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '148')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '148', '28')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 148 -> -24 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '148')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '148', '-24')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 148 -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '148')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '148', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 148 -> -112 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '148')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '148', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 148 -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '148')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '148', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 148 -> 91 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '148')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '148', '91')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 148 -> -35 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '148')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '148', '-35')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 148 -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '148')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '148', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 5_sg -> -90 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sg')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sg', '-90')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-90', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 5_sg -> -57 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sg')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sg', '-57')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-57', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 5_sg -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sg')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sg', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 5_sg -> -112 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sg')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sg', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 5_sg -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sg')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sg', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 5_sg -> 91 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sg')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sg', '91')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 153 -> 5_sp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '153', '5_sp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 153 -> dat. du. -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '153', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 153 -> -151 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '153', '-151')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-151', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 153 -> -90 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '153', '-90')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-90', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 153 -> -57 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '153', '-57')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-57', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 153 -> 8_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '153', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 153 -> 169 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '153', '169')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '169', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 153 -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '153', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 153 -> pl -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '153', 'pl')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 153 -> 153 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '153', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 153 -> 135 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '153', '135')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '135', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 153 -> 28 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '153', '28')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 153 -> -24 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '153', '-24')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 153 -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '153', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 153 -> -112 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '153', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 153 -> 89 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '153', '89')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '89', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 153 -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '153', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 153 -> 91 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '153', '91')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 153 -> -35 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '153', '-35')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 153 -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '153', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -30 -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 11_sp -> dat. du. -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_sp', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 11_sp -> 8_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_sp', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 11_sp -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_sp', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 11_sp -> pl -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_sp', 'pl')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 11_sp -> 153 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_sp', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 11_sp -> 28 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_sp', '28')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 11_sp -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_sp', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 11_sp -> -112 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_sp', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 11_sp -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_sp', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 11_sp -> -35 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_sp', '-35')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 11_sp -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_sp', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 30_sp -> 8_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_sp', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 30_sp -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_sp', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 30_sp -> -24 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_sp', '-24')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 30_sp -> -112 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_sp', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 30_sp -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_sp', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 30_sp -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_sp', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 135 -> dat. du. -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '135')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 135 -> -151 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '135')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', '-151')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-151', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 135 -> -90 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '135')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', '-90')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-90', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 135 -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '135')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 135 -> -24 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '135')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', '-24')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 135 -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '135')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 135 -> -112 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '135')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 135 -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '135')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -86 -> dat. du. -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-86')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-86', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -86 -> -90 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-86')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-86', '-90')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-90', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -86 -> 8_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-86')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-86', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -86 -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-86')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-86', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -86 -> 153 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-86')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-86', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -86 -> -24 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-86')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-86', '-24')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -86 -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-86')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-86', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -86 -> -112 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-86')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-86', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -86 -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-86')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-86', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -86 -> 91 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-86')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-86', '91')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -86 -> -35 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-86')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-86', '-35')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -86 -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-86')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-86', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 7_fp -> 5_sp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_fp', '5_sp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 7_fp -> dat. du. -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_fp', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 7_fp -> -57 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_fp', '-57')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-57', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 7_fp -> 8_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_fp', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 7_fp -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_fp', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 7_fp -> pl -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_fp', 'pl')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 7_fp -> 153 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_fp', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 7_fp -> -24 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_fp', '-24')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 7_fp -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_fp', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 7_fp -> -112 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_fp', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 7_fp -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_fp', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 7_fp -> 91 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_fp', '91')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 7_fp -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_fp', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 111 -> dat. du. -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '111')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '111', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 111 -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '111')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '111', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 111 -> -24 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '111')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '111', '-24')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 111 -> -112 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '111')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '111', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 111 -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '111')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '111', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -92 -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -92 -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 177 -> dat. du. -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '177')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '177', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 177 -> -90 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '177')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '177', '-90')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-90', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 177 -> 8_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '177')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '177', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 177 -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '177')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '177', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 177 -> pl -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '177')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '177', 'pl')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 177 -> 153 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '177')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '177', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 177 -> 28 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '177')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '177', '28')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 177 -> -24 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '177')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '177', '-24')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 177 -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '177')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '177', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 177 -> -112 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '177')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '177', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 177 -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '177')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '177', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 177 -> 91 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '177')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '177', '91')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 177 -> -35 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '177')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '177', '-35')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 177 -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '177')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '177', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> neutr -> dat. du. -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'neutr')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'neutr', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> neutr -> -151 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'neutr')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'neutr', '-151')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-151', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> neutr -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'neutr')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'neutr', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> neutr -> pl -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'neutr')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'neutr', 'pl')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> neutr -> 153 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'neutr')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'neutr', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> neutr -> -24 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'neutr')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'neutr', '-24')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> neutr -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'neutr')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'neutr', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> neutr -> -112 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'neutr')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'neutr', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> neutr -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'neutr')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'neutr', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> neutr -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'neutr')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'neutr', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 59 -> dat. du. -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '59')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '59', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 59 -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '59')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '59', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 59 -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '59')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '59', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 59 -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '59')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '59', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 29 -> 5_sp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '29')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '29', '5_sp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 29 -> dat. du. -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '29')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '29', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 29 -> 8_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '29')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '29', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 29 -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '29')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '29', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 29 -> pl -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '29')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '29', 'pl')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 29 -> 153 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '29')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '29', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 29 -> 28 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '29')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '29', '28')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 29 -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '29')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '29', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 29 -> -112 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '29')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '29', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 29 -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '29')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '29', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 29 -> -35 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '29')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '29', '-35')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 29 -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '29')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '29', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 156 -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '156')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '156', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 156 -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '156')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '156', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -153 -> 5_sp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-153', '5_sp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -153 -> dat. du. -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-153', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -153 -> -151 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-153', '-151')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-151', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -153 -> -90 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-153', '-90')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-90', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -153 -> -57 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-153', '-57')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-57', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -153 -> 8_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-153', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -153 -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-153', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -153 -> pl -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-153', 'pl')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -153 -> 153 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-153', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -153 -> -24 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-153', '-24')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -153 -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-153', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -153 -> -112 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-153', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -153 -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-153', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -153 -> 91 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-153', '91')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -153 -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-153', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 8_tp -> 5_sp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_tp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_tp', '5_sp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 8_tp -> dat. du. -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_tp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_tp', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 8_tp -> 8_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_tp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_tp', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 8_tp -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_tp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_tp', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 8_tp -> pl -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_tp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_tp', 'pl')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 8_tp -> 153 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_tp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_tp', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 8_tp -> 28 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_tp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_tp', '28')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 8_tp -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_tp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_tp', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 8_tp -> -112 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_tp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_tp', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 8_tp -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_tp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_tp', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 8_tp -> -35 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_tp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_tp', '-35')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 8_tp -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_tp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_tp', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 4_tp -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_tp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_tp', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 28 -> dat. du. -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 28 -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 28 -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 28 -> -112 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 28 -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 28 -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -24 -> 5_sp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-24')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-24', '5_sp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -24 -> dat. du. -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-24')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-24', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -24 -> -151 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-24')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-24', '-151')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-151', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -24 -> -90 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-24')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-24', '-90')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-90', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -24 -> -57 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-24')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-24', '-57')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-57', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -24 -> 8_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-24')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-24', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -24 -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-24')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-24', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -24 -> pl -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-24')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-24', 'pl')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -24 -> 153 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-24')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-24', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -24 -> 28 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-24')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-24', '28')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -24 -> -24 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-24')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-24', '-24')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -24 -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-24')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-24', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -24 -> -112 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-24')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-24', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -24 -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-24')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-24', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -24 -> 91 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-24')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-24', '91')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -24 -> -35 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-24')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-24', '-35')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -24 -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-24')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-24', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 137 -> dat. du. -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '137')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '137', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 137 -> -151 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '137')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '137', '-151')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-151', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 137 -> -90 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '137')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '137', '-90')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-90', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 137 -> 8_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '137')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '137', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 137 -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '137')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '137', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 137 -> 153 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '137')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '137', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 137 -> -24 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '137')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '137', '-24')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 137 -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '137')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '137', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 137 -> -112 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '137')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '137', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 137 -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '137')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '137', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 137 -> 91 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '137')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '137', '91')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 137 -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '137')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '137', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -112 -> 5_sp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-112')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-112', '5_sp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -112 -> dat. du. -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-112')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-112', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -112 -> -151 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-112')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-112', '-151')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-151', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -112 -> -57 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-112')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-112', '-57')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-57', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -112 -> 8_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-112')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-112', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -112 -> 169 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-112')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-112', '169')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '169', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -112 -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-112')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-112', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -112 -> pl -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-112')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-112', 'pl')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -112 -> 153 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-112')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-112', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -112 -> 28 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-112')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-112', '28')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -112 -> -24 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-112')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-112', '-24')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -112 -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-112')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-112', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -112 -> -112 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-112')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-112', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -112 -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-112')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-112', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -112 -> 91 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-112')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-112', '91')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -112 -> -35 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-112')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-112', '-35')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -112 -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-112')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-112', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 89 -> -24 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '89')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '89', '-24')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 89 -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '89')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '89', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -161 -> 8_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-161')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-161', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -161 -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-161')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-161', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -161 -> pl -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-161')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-161', 'pl')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -161 -> 153 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-161')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-161', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -161 -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-161')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-161', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -161 -> -112 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-161')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-161', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -161 -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-161')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-161', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -161 -> -35 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-161')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-161', '-35')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -161 -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-161')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-161', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -96 -> dat. du. -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-96')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-96', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -96 -> 8_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-96')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-96', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -96 -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-96')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-96', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -96 -> pl -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-96')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-96', 'pl')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -96 -> 153 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-96')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-96', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -96 -> -24 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-96')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-96', '-24')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -96 -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-96')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-96', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -96 -> -112 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-96')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-96', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -96 -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-96')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-96', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -96 -> 91 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-96')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-96', '91')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -96 -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-96')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-96', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 174 -> 8_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '174')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '174', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 174 -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '174')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '174', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 174 -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '174')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '174', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 162 -> dat. du. -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '162')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '162', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 162 -> 8_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '162')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '162', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 162 -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '162')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '162', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 162 -> pl -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '162')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '162', 'pl')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 162 -> 28 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '162')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '162', '28')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 162 -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '162')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '162', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 162 -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '162')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '162', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 162 -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '162')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '162', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> gen. pl. -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. pl.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. pl.', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> gen. pl. -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. pl.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. pl.', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> pl_tp -> 5_sp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_tp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_tp', '5_sp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> pl_tp -> dat. du. -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_tp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_tp', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> pl_tp -> -57 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_tp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_tp', '-57')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-57', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> pl_tp -> 8_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_tp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_tp', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> pl_tp -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_tp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_tp', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> pl_tp -> pl -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_tp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_tp', 'pl')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> pl_tp -> 153 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_tp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_tp', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> pl_tp -> -24 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_tp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_tp', '-24')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> pl_tp -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_tp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_tp', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> pl_tp -> -112 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_tp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_tp', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> pl_tp -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_tp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_tp', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> pl_tp -> 91 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_tp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_tp', '91')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> pl_tp -> -35 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_tp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_tp', '-35')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> pl_tp -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_tp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_tp', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 101 -> 5_sp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '101')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '101', '5_sp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 101 -> dat. du. -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '101')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '101', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 101 -> -151 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '101')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '101', '-151')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-151', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 101 -> -90 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '101')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '101', '-90')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-90', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 101 -> -57 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '101')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '101', '-57')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-57', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 101 -> 8_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '101')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '101', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 101 -> 169 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '101')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '101', '169')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '169', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 101 -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '101')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '101', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 101 -> pl -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '101')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '101', 'pl')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 101 -> 153 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '101')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '101', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 101 -> 135 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '101')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '101', '135')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '135', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 101 -> 28 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '101')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '101', '28')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 101 -> -24 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '101')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '101', '-24')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 101 -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '101')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '101', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 101 -> -112 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '101')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '101', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 101 -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '101')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '101', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 101 -> 91 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '101')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '101', '91')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 101 -> -35 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '101')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '101', '-35')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 101 -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '101')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '101', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 91 -> 5_sp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '91')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '91', '5_sp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 91 -> dat. du. -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '91')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '91', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 91 -> -151 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '91')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '91', '-151')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-151', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 91 -> -90 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '91')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '91', '-90')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-90', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 91 -> -57 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '91')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '91', '-57')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-57', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 91 -> 8_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '91')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '91', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 91 -> 169 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '91')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '91', '169')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '169', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 91 -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '91')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '91', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 91 -> pl -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '91')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '91', 'pl')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 91 -> 153 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '91')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '91', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 91 -> 135 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '91')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '91', '135')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '135', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 91 -> 28 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '91')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '91', '28')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 91 -> -24 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '91')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '91', '-24')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 91 -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '91')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '91', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 91 -> -112 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '91')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '91', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 91 -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '91')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '91', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 91 -> 91 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '91')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '91', '91')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 91 -> -35 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '91')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '91', '-35')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 91 -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '91')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '91', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -35 -> -90 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-35')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-35', '-90')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-90', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -35 -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-35')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-35', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -35 -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-35')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-35', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 11_fp -> 5_sp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_fp', '5_sp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 11_fp -> dat. du. -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_fp', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 11_fp -> -57 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_fp', '-57')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-57', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 11_fp -> 8_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_fp', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 11_fp -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_fp', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 11_fp -> pl -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_fp', 'pl')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 11_fp -> 153 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_fp', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 11_fp -> 28 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_fp', '28')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 11_fp -> -24 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_fp', '-24')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 11_fp -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_fp', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 11_fp -> -112 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_fp', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 11_fp -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_fp', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 11_fp -> 91 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_fp', '91')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 11_fp -> -35 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_fp', '-35')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 11_fp -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_fp', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 3_du -> 5_sp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_du', '5_sp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 3_du -> dat. du. -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_du', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 3_du -> -151 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_du', '-151')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-151', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 3_du -> -90 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_du', '-90')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-90', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 3_du -> -57 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_du', '-57')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-57', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 3_du -> 8_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_du', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 3_du -> 169 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_du', '169')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '169', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 3_du -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_du', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 3_du -> 114 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_du', '114')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '114', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 3_du -> pl -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_du', 'pl')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 3_du -> 153 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_du', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 3_du -> 28 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_du', '28')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 3_du -> -24 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_du', '-24')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 3_du -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_du', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 3_du -> -112 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_du', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 3_du -> 89 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_du', '89')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '89', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 3_du -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_du', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 3_du -> 91 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_du', '91')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 3_du -> -35 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_du', '-35')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 3_du -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_du', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 8_sp -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_sp', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 8_sp -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_sp', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -296 -> dat. du. -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-296')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-296', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -296 -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-296')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-296', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -296 -> 28 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-296')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-296', '28')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -296 -> -24 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-296')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-296', '-24')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -296 -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-296')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-296', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -296 -> -112 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-296')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-296', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -296 -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-296')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-296', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -296 -> 91 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-296')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-296', '91')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -296 -> -35 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-296')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-296', '-35')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -296 -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-296')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-296', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 3_sp -> 5_sp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_sp', '5_sp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 3_sp -> dat. du. -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_sp', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 3_sp -> 8_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_sp', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 3_sp -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_sp', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 3_sp -> pl -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_sp', 'pl')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 3_sp -> 153 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_sp', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 3_sp -> 28 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_sp', '28')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 3_sp -> -24 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_sp', '-24')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 3_sp -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_sp', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 3_sp -> -112 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_sp', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 3_sp -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_sp', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 3_sp -> 91 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_sp', '91')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 3_sp -> -35 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_sp', '-35')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 3_sp -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_sp', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 120 -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '120')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '120', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 94 -> 5_sp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '94')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '94', '5_sp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 94 -> dat. du. -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '94')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '94', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 94 -> 8_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '94')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '94', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 94 -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '94')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '94', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 94 -> 153 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '94')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '94', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 94 -> -24 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '94')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '94', '-24')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 94 -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '94')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '94', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 94 -> -112 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '94')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '94', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 94 -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '94')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '94', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 94 -> 91 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '94')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '94', '91')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 94 -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '94')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '94', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -200 -> 5_sp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '5_sp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -200 -> dat. du. -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -200 -> 8_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -200 -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -200 -> pl -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', 'pl')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -200 -> 153 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -200 -> 28 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '28')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -200 -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -200 -> -112 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -200 -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -200 -> -35 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '-35')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -200 -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -81 -> 5_sp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', '5_sp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -81 -> dat. du. -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -81 -> 8_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -81 -> 54 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -81 -> pl -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', 'pl')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -81 -> 153 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -81 -> 28 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', '28')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -81 -> 137 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -81 -> -112 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -81 -> 101 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -81 -> -35 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', '-35')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -81 -> 3_du -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> tp -> dat. du. -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'dat. du.', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> tp -> 101 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '101')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '101', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 5_du -> 153 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_du', '153')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '153', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 157 -> 54 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '157')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '157', '54')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '54', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 157 -> 137 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '157')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '157', '137')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '137', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 157 -> 101 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '157')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '157', '101')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '101', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 157 -> 3_du -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '157')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '157', '3_du')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '3_du', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 5_sp -> -90 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', '-90')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-90', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 5_sp -> 91 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', '91')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '91', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -151 -> -24 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-151')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-151', '-24')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-24', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> acc. sg. -> 101 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. sg.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. sg.', '101')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '101', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -25 -> 137 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-25')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-25', '137')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '137', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -25 -> 101 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-25')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-25', '101')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '101', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -57 -> 91 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-57')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-57', '91')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '91', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 8_fp -> 91 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_fp', '91')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '91', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 54 -> 8_fp -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '54')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '54', '8_fp')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '8_fp', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 55 -> 101 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '55')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '55', '101')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '101', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -156 -> 101 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-156')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-156', '101')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '101', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -29 -> 54 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-29')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-29', '54')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '54', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -29 -> 137 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-29')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-29', '137')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '137', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -29 -> 101 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-29')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-29', '101')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '101', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 153 -> -24 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '153', '-24')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-24', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -24 -> pl -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-24')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-24', 'pl')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'pl', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> acc -> dat. du. -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'dat. du.', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> acc -> 101 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc', '101')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '101', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 13_pl -> 54 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '54')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '54', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 101 -> -151 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '101')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '101', '-151')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-151', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 101 -> 8_fp -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '101')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '101', '8_fp')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '8_fp', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 91 -> 91 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '91')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '91', '91')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '91', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 3_du -> 91 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_du', '91')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '91', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 35 -> 8_fp -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', '8_fp')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '8_fp', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 35 -> 54 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', '54')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '54', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 35 -> pl -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', 'pl')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'pl', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 35 -> 101 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', '101')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '101', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> gen. sg. -> 137 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'gen. sg.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. sg.', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> gen. sg. -> 101 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'gen. sg.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. sg.', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 5_sp -> 54 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '5_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 5_sp -> -24 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '5_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', '-24')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 5_sp -> 137 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '5_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 5_sp -> -112 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '5_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 5_sp -> 101 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '5_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 5_sp -> 91 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '5_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', '91')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 5_sp -> 3_du -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '5_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> dat. du. -> 54 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'dat. du.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. du.', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 15_du -> dat. du. -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '15_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_du', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 15_du -> 8_fp -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '15_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_du', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 15_du -> 54 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '15_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_du', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 15_du -> pl -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '15_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_du', 'pl')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 15_du -> 153 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '15_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_du', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 15_du -> 28 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '15_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_du', '28')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 15_du -> 137 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '15_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_du', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 15_du -> -112 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '15_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_du', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 15_du -> 101 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '15_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_du', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 15_du -> 3_du -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '15_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_du', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> -151 -> dat. du. -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-151')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-151', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> -151 -> 8_fp -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-151')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-151', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> -151 -> 54 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-151')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-151', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> -151 -> pl -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-151')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-151', 'pl')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> -151 -> 153 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-151')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-151', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> -151 -> 28 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-151')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-151', '28')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> -151 -> 137 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-151')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-151', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> -151 -> -112 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-151')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-151', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> -151 -> 101 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-151')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-151', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> -151 -> -35 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-151')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-151', '-35')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> -151 -> 3_du -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-151')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-151', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> -306 -> 137 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-306')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> -306 -> 101 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-306')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> -90 -> dat. du. -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-90')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-90', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> -90 -> 8_fp -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-90')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-90', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> -90 -> 54 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-90')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-90', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> -90 -> 153 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-90')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-90', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> -90 -> 137 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-90')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-90', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> -90 -> -112 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-90')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-90', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> -90 -> 101 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-90')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-90', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> -90 -> 3_du -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-90')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-90', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> -57 -> 8_fp -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-57')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-57', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> -57 -> 54 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-57')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-57', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> -57 -> pl -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-57')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-57', 'pl')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> -57 -> 153 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-57')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-57', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> -57 -> 137 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-57')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-57', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> -57 -> -112 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-57')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-57', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> -57 -> 101 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-57')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-57', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 8_fp -> dat. du. -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '8_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_fp', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 8_fp -> 8_fp -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '8_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_fp', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 8_fp -> 54 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '8_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_fp', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 8_fp -> pl -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '8_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_fp', 'pl')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 8_fp -> 153 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '8_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_fp', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 8_fp -> 28 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '8_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_fp', '28')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 8_fp -> -24 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '8_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_fp', '-24')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 8_fp -> 137 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '8_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_fp', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 8_fp -> -112 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '8_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_fp', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 8_fp -> 101 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '8_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_fp', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 8_fp -> 91 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '8_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_fp', '91')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 8_fp -> -35 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '8_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_fp', '-35')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 8_fp -> 3_du -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '8_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_fp', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 169 -> dat. du. -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '169')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '169', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 169 -> 8_fp -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '169')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '169', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 54 -> dat. du. -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '54')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '54', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 54 -> -151 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '54')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '54', '-151')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-151', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 54 -> -90 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '54')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '54', '-90')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-90', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 54 -> 8_fp -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '54')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '54', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 54 -> 54 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '54')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '54', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 54 -> pl -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '54')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '54', 'pl')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 54 -> 153 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '54')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '54', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 54 -> 135 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '54')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '54', '135')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '135', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 54 -> -24 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '54')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '54', '-24')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 54 -> 137 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '54')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '54', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 54 -> -112 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '54')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '54', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 54 -> 101 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '54')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '54', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 54 -> 91 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '54')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '54', '91')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 54 -> -35 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '54')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '54', '-35')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 54 -> 3_du -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '54')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '54', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 114 -> 54 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '114')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '114', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 114 -> 101 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '114')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '114', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> -139 -> 137 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-139')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> -139 -> 101 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-139')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 139 -> dat. du. -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '139')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '139', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 139 -> 8_fp -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '139')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '139', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 139 -> 54 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '139')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '139', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 139 -> pl -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '139')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '139', 'pl')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 139 -> 153 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '139')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '139', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 139 -> 28 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '139')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '139', '28')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 139 -> 137 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '139')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '139', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 139 -> -112 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '139')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '139', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 139 -> 101 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '139')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '139', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 139 -> -35 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '139')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '139', '-35')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 139 -> 3_du -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '139')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '139', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> pl -> -90 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'pl')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl', '-90')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-90', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> pl -> 8_fp -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'pl')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> pl -> 54 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'pl')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> pl -> -24 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'pl')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl', '-24')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> pl -> -112 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'pl')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> pl -> 101 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'pl')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 153 -> dat. du. -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '153', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 153 -> 8_fp -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '153', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 153 -> 54 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '153', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 153 -> pl -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '153', 'pl')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 153 -> 153 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '153', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 153 -> 28 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '153', '28')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 153 -> -24 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '153', '-24')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 153 -> 137 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '153', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 153 -> -112 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '153', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 153 -> 101 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '153', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 153 -> 91 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '153', '91')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 153 -> -35 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '153', '-35')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 153 -> 3_du -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '153', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 135 -> 54 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '135')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 135 -> pl -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '135')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', 'pl')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 135 -> 137 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '135')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 135 -> 3_du -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '135')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 7_fp -> 101 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '7_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_fp', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> neutr -> 137 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'neutr')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'neutr', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> neutr -> 101 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'neutr')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'neutr', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> -153 -> dat. du. -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-153', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> -153 -> 8_fp -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-153', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> -153 -> 54 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-153', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> -153 -> pl -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-153', 'pl')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> -153 -> 28 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-153', '28')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> -153 -> 137 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-153', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> -153 -> 101 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-153', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> -153 -> -35 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-153', '-35')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> -153 -> 3_du -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-153', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> -24 -> 8_fp -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-24')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-24', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> -24 -> 54 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-24')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-24', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> -24 -> pl -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-24')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-24', 'pl')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> -24 -> 137 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-24')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-24', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> -24 -> -112 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-24')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-24', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> -24 -> 101 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-24')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-24', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 137 -> dat. du. -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '137')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '137', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 137 -> -151 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '137')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '137', '-151')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-151', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 137 -> 54 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '137')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '137', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 137 -> 153 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '137')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '137', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 137 -> 28 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '137')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '137', '28')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 137 -> -24 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '137')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '137', '-24')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 137 -> 137 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '137')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '137', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 137 -> -112 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '137')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '137', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 137 -> 101 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '137')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '137', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 137 -> 3_du -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '137')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '137', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> -112 -> 153 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-112')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-112', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> -112 -> -112 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-112')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-112', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 101 -> dat. du. -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '101')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '101', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 101 -> -151 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '101')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '101', '-151')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-151', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 101 -> -90 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '101')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '101', '-90')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-90', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 101 -> -57 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '101')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '101', '-57')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-57', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 101 -> 8_fp -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '101')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '101', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 101 -> 169 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '101')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '101', '169')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '169', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 101 -> 54 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '101')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '101', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 101 -> 114 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '101')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '101', '114')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '114', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 101 -> pl -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '101')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '101', 'pl')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 101 -> 153 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '101')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '101', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 101 -> 28 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '101')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '101', '28')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 101 -> -24 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '101')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '101', '-24')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 101 -> 137 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '101')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '101', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 101 -> -112 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '101')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '101', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 101 -> 101 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '101')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '101', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 101 -> 91 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '101')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '101', '91')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 101 -> -35 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '101')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '101', '-35')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 101 -> 3_du -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '101')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '101', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 91 -> dat. du. -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '91')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '91', 'dat. du.')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 91 -> 54 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '91')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '91', '54')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 91 -> 153 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '91')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '91', '153')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 91 -> -24 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '91')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '91', '-24')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 91 -> 137 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '91')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '91', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 91 -> -112 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '91')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '91', '-112')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 91 -> 101 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '91')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '91', '101')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 91 -> 91 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '91')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '91', '91')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 91 -> -35 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '91')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '91', '-35')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 91 -> 3_du -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '91')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '91', '3_du')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 3_du -> 137 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '3_du')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_du', '137')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 153 -> 91 -> T
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '153')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '153', '91')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '91', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 28 -> 91 -> T
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '28')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28', '91')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '91', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> -296 -> 8_fp -> T
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-296')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-296', '8_fp')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '8_fp', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    feats[feats < 1e-25] = 1e-25
    return -np.log10(feats)
    