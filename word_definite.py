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
    _edge_vector_dim = 400
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
    
    # L->128->L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '128')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '128', node2.lemma)
    feats[fIndex] = pleft * pright; fIndex += 1
            
    # L->77->L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '77')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '77', node2.lemma)
    feats[fIndex] = pleft * pright; fIndex += 1
            
    # L->160->L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '160')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    feats[fIndex] = pleft * pright; fIndex += 1
            
    # L->50->L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '50')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '50', node2.lemma)
    feats[fIndex] = pleft * pright; fIndex += 1
            
    # L->128->T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '128')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '128', node2.tup)
    feats[fIndex] = pleft * pright; fIndex += 1
            
    # L->93->T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '93')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '93', node2.tup)
    feats[fIndex] = pleft * pright; fIndex += 1
            
    # L->50->T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '50')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '50', node2.tup)
    feats[fIndex] = pleft * pright; fIndex += 1
            
    # C->134->T
    pleft = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '134')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '134', node2.tup)
    feats[fIndex] = pleft * pright; fIndex += 1
            
    # C->128->T
    pleft = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '128')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '128', node2.tup)
    feats[fIndex] = pleft * pright; fIndex += 1
            
    # C->77->T
    pleft = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '77')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '77', node2.tup)
    feats[fIndex] = pleft * pright; fIndex += 1
            
    # C->160->T
    pleft = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '160')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '160', node2.tup)
    feats[fIndex] = pleft * pright; fIndex += 1
            
    # C->50->T
    pleft = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '50')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '50', node2.tup)
    feats[fIndex] = pleft * pright; fIndex += 1
            
    # T->128->L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '128')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '128', node2.lemma)
    feats[fIndex] = pleft * pright; fIndex += 1
            
    # T->77->L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '77')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '77', node2.lemma)
    feats[fIndex] = pleft * pright; fIndex += 1
            
    # T->160->L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '160')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    feats[fIndex] = pleft * pright; fIndex += 1
            
    # T->160->T
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '160')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '160', node2.tup)
    feats[fIndex] = pleft * pright; fIndex += 1
            
    # T->93->T
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '93')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '93', node2.tup)
    feats[fIndex] = pleft * pright; fIndex += 1
            
    # L -> 180 -> 134 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '180')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '180', '134')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '134', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 180 -> 160 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '180')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '180', '160')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 180 -> 50 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '180')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '180', '50')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '50', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 121 -> 160 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '121')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '121', '160')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 151 -> sg_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '151')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '151', 'sg_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 151 -> 160 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '151')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '151', '160')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 134 -> 128 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '134')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '134', '128')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '128', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 134 -> -46 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '134')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '134', '-46')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-46', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 134 -> 160 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '134')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '134', '160')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> gen. sg. -> sg_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. sg.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. sg.', 'sg_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> gen. sg. -> -46 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. sg.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. sg.', '-46')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-46', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> gen. sg. -> 77 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. sg.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. sg.', '77')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '77', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> gen. sg. -> 160 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. sg.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. sg.', '160')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 13_fp -> 128 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_fp', '128')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '128', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 13_fp -> 77 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_fp', '77')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '77', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 13_fp -> 160 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_fp', '160')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 13_fp -> 50 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_fp', '50')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '50', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> instr. masc. -> sg_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. masc.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. masc.', 'sg_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> instr. masc. -> 128 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. masc.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. masc.', '128')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '128', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> instr. masc. -> 77 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. masc.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. masc.', '77')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '77', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> instr. masc. -> 160 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. masc.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. masc.', '160')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> instr. masc. -> 50 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. masc.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. masc.', '50')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '50', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 74 -> 160 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '74')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '74', '160')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -79 -> 160 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-79')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-79', '160')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -27 -> 134 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-27')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-27', '134')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '134', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 28_tp -> 128 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_tp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28_tp', '128')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '128', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 28_tp -> 77 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_tp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28_tp', '77')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '77', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 28_tp -> 160 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_tp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28_tp', '160')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 28_tp -> 50 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_tp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28_tp', '50')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '50', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -59 -> 77 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-59')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-59', '77')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '77', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -59 -> 160 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-59')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-59', '160')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> sg_fp -> 77 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg_fp', '77')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '77', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> sg_fp -> 160 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg_fp', '160')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 6_sp -> 134 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sp', '134')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '134', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 6_sp -> 77 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sp', '77')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '77', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 6_sp -> 160 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sp', '160')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 132 -> 134 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '132')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '132', '134')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '134', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 132 -> 128 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '132')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '132', '128')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '128', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 132 -> 160 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '132')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '132', '160')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 128 -> 134 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '128')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '128', '134')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '134', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 128 -> 128 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '128')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '128', '128')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '128', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 128 -> -46 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '128')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '128', '-46')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-46', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 128 -> 77 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '128')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '128', '77')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '77', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 128 -> 130 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '128')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '128', '130')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '130', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 128 -> 160 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '128')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '128', '160')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 128 -> 50 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '128')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '128', '50')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '50', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> nom. sg. -> sg_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. sg.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. sg.', 'sg_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> nom. sg. -> 128 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. sg.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. sg.', '128')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '128', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> nom. sg. -> 77 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. sg.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. sg.', '77')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '77', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> nom. sg. -> 160 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. sg.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. sg.', '160')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> nom. sg. -> 50 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. sg.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. sg.', '50')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '50', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -159 -> 134 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-159')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-159', '134')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '134', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -159 -> 128 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-159')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-159', '128')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '128', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -159 -> 77 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-159')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-159', '77')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '77', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -159 -> 160 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-159')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-159', '160')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -159 -> 50 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-159')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-159', '50')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '50', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 3_pl -> 77 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_pl')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_pl', '77')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '77', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 3_pl -> 160 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_pl')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_pl', '160')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> loc. pl. -> 134 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. pl.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. pl.', '134')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '134', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> loc. pl. -> sg_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. pl.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. pl.', 'sg_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> loc. pl. -> 128 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. pl.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. pl.', '128')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '128', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> loc. pl. -> 5_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. pl.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. pl.', '5_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> loc. pl. -> 77 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. pl.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. pl.', '77')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '77', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> loc. pl. -> 160 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. pl.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. pl.', '160')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> loc. pl. -> 50 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. pl.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. pl.', '50')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '50', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -46 -> 128 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-46')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-46', '128')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '128', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 159 -> 128 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '128')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '128', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 159 -> 77 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '77')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '77', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 159 -> 160 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '160')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 2_sp -> 77 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_sp', '77')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '77', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 2_sp -> 160 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_sp', '160')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -111 -> 77 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-111')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-111', '77')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '77', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -111 -> 160 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-111')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-111', '160')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -158 -> sg_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-158')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-158', 'sg_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -158 -> 128 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-158')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-158', '128')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '128', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -158 -> 77 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-158')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-158', '77')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '77', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -158 -> 160 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-158')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-158', '160')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 42 -> 160 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '42')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '42', '160')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -263 -> sg_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-263')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-263', 'sg_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -263 -> 77 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-263')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-263', '77')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '77', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -263 -> 160 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-263')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-263', '160')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -123 -> 77 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-123')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-123', '77')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '77', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -123 -> 160 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-123')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-123', '160')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 2_pl -> 160 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_pl')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_pl', '160')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 172 -> 128 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '172')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '172', '128')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '128', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 172 -> 160 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '172')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '172', '160')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 16_tp -> 77 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', '77')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '77', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -157 -> 160 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '160')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 157 -> 128 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '157')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '157', '128')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '128', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 157 -> 77 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '157')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '157', '77')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '77', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 157 -> 160 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '157')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '157', '160')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 109 -> 160 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '109')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '109', '160')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 5_fp -> 128 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', '128')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '128', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 5_fp -> 130 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', '130')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '130', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 88 -> 160 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '88')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '88', '160')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 88 -> 50 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '88')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '88', '50')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '50', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 77 -> 128 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '77')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '77', '128')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '128', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 77 -> 77 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '77')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '77', '77')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '77', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 77 -> 160 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '77')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '77', '160')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 130 -> 128 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '130')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '130', '128')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '128', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 130 -> 77 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '130')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '130', '77')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '77', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 130 -> 160 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '130')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '130', '160')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 130 -> 50 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '130')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '130', '50')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '50', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 81 -> 160 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '81')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '81', '160')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 155 -> 134 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '155')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '155', '134')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '134', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 155 -> 77 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '155')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '155', '77')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '77', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 155 -> 160 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '155')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '155', '160')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 160 -> 134 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '160')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '160', '134')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '134', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 160 -> sg_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '160')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '160', 'sg_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 160 -> 128 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '160')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '160', '128')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '128', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 160 -> -46 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '160')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '160', '-46')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-46', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 160 -> 5_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '160')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '160', '5_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 160 -> 77 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '160')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '160', '77')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '77', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 160 -> 160 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '160')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '160', '160')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 160 -> 50 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '160')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '160', '50')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '50', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 93 -> 134 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '93')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '93', '134')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '134', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 93 -> sg_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '93')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '93', 'sg_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 93 -> 128 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '93')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '93', '128')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '128', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 93 -> 77 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '93')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '93', '77')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '77', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 93 -> 160 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '93')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '93', '160')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 93 -> 50 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '93')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '93', '50')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '50', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 8_fp -> 128 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_fp', '128')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '128', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -220 -> 134 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-220')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-220', '134')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '134', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -220 -> 77 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-220')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-220', '77')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '77', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -220 -> 160 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-220')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-220', '160')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -220 -> 50 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-220')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-220', '50')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '50', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 50 -> 134 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '50')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '50', '134')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '134', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 50 -> 128 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '50')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '50', '128')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '128', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 50 -> 5_fp -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '50')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '50', '5_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 50 -> 77 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '50')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '50', '77')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '77', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 50 -> 160 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '50')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '50', '160')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 50 -> 50 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '50')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '50', '50')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '50', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 14_tp -> 77 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_tp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_tp', '77')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '77', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 139 -> 134 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '139')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '139', '134')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '134', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 139 -> 77 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '139')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '139', '77')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '77', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 139 -> 160 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '139')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '139', '160')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 139 -> 50 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '139')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '139', '50')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '50', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -166 -> 160 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-166')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-166', '160')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -43 -> 160 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-43')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-43', '160')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -30 -> 77 -> L
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '77')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '77', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 108 -> sg_fp -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '108')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '108', 'sg_fp')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'sg_fp', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 108 -> 128 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '108')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '108', '128')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '128', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 108 -> -46 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '108')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '108', '-46')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-46', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 108 -> 77 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '108')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '108', '77')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '77', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 108 -> 160 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '108')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '108', '160')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '160', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 108 -> 50 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '108')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '108', '50')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '50', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 149 -> 77 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '77')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '77', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 149 -> 160 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '160')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '160', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 180 -> -158 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '180')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '180', '-158')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-158', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 180 -> 93 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '180')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '180', '93')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '93', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 121 -> -158 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '121')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '121', '-158')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-158', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 121 -> 93 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '121')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '121', '93')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '93', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> instr. pl. -> 160 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. pl.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. pl.', '160')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '160', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 11_fp -> 77 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_fp', '77')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '77', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 11_fp -> 160 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_fp', '160')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '160', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 151 -> sg_fp -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '151')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '151', 'sg_fp')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'sg_fp', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 151 -> 128 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '151')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '151', '128')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '128', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 151 -> 93 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '151')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '151', '93')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '93', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 151 -> 50 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '151')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '151', '50')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '50', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 134 -> sg_fp -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '134')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '134', 'sg_fp')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'sg_fp', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 134 -> 93 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '134')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '134', '93')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '93', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 134 -> 50 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '134')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '134', '50')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '50', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 2_tp -> sg_fp -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_tp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_tp', 'sg_fp')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'sg_fp', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 2_tp -> 93 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_tp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_tp', '93')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '93', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 13_pl -> 134 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '134')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '134', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 13_pl -> 77 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '77')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '77', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 13_pl -> 160 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '160')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '160', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> gen. sg. -> sg_fp -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. sg.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. sg.', 'sg_fp')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'sg_fp', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> gen. sg. -> 88 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. sg.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. sg.', '88')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '88', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> gen. sg. -> 160 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. sg.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. sg.', '160')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '160', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> gen. sg. -> 93 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. sg.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. sg.', '93')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '93', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 3_sp -> 160 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_sp', '160')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '160', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> pl_sp -> 134 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_sp', '134')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '134', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> pl_sp -> 77 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_sp', '77')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '77', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> pl_sp -> 160 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_sp', '160')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '160', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -139 -> 160 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '160')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '160', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 13_fp -> -46 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_fp', '-46')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-46', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 13_fp -> 88 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_fp', '88')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '88', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 13_fp -> 160 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_fp', '160')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '160', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -297 -> 134 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-297')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-297', '134')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '134', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -297 -> 77 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-297')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-297', '77')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '77', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -297 -> 160 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-297')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-297', '160')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '160', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 91 -> 160 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '91')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '91', '160')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '160', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> instr. masc. -> 128 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. masc.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. masc.', '128')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '128', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> instr. masc. -> 160 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. masc.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. masc.', '160')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '160', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> instr. masc. -> 93 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. masc.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. masc.', '93')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '93', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -79 -> -158 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-79')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-79', '-158')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-158', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -27 -> 93 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-27')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-27', '93')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '93', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -18 -> 77 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '77')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '77', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -18 -> 160 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '160')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '160', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 28_tp -> -158 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_tp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28_tp', '-158')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-158', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 28_tp -> 88 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_tp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28_tp', '88')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '88', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 28_tp -> 93 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_tp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28_tp', '93')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '93', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 28_tp -> 50 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_tp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28_tp', '50')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '50', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> sg_fp -> -158 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg_fp', '-158')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-158', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> sg_fp -> 160 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg_fp', '160')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '160', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> sg_fp -> 93 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg_fp', '93')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '93', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 132 -> -158 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '132')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '132', '-158')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-158', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 132 -> 8_fp -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '132')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '132', '8_fp')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '8_fp', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 15_tp -> 160 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_tp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_tp', '160')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '160', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 128 -> sg_fp -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '128')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '128', 'sg_fp')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'sg_fp', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 128 -> 128 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '128')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '128', '128')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '128', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 128 -> -158 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '128')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '128', '-158')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-158', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 128 -> 160 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '128')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '128', '160')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '160', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 128 -> 93 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '128')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '128', '93')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '93', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> nom. sg. -> 93 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. sg.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. sg.', '93')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '93', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 119 -> 93 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '119')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '119', '93')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '93', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -159 -> 128 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-159')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-159', '128')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '128', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -159 -> 160 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-159')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-159', '160')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '160', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -159 -> 50 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-159')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-159', '50')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '50', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -78 -> 134 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-78')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-78', '134')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '134', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -78 -> sg_fp -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-78')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-78', 'sg_fp')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'sg_fp', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -78 -> 128 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-78')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-78', '128')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '128', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -78 -> -46 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-78')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-78', '-46')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-46', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -78 -> 77 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-78')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-78', '77')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '77', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -78 -> 160 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-78')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-78', '160')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '160', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -78 -> 50 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-78')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-78', '50')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '50', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 3_pl -> 160 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_pl')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_pl', '160')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '160', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> nom. fem -> 134 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. fem')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. fem', '134')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '134', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> nom. fem -> 128 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. fem')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. fem', '128')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '128', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> nom. fem -> 5_fp -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. fem')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. fem', '5_fp')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '5_fp', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> nom. fem -> 77 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. fem')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. fem', '77')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '77', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> nom. fem -> 160 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. fem')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. fem', '160')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '160', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> nom. fem -> 50 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. fem')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. fem', '50')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '50', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> loc. pl. -> 128 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. pl.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. pl.', '128')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '128', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> loc. pl. -> 160 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. pl.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. pl.', '160')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '160', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> loc. pl. -> 93 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. pl.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. pl.', '93')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '93', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 54 -> -158 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '54')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '54', '-158')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-158', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 2_sp -> -158 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_sp', '-158')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-158', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -163 -> 134 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-163')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-163', '134')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '134', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -163 -> 77 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-163')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-163', '77')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '77', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -163 -> 160 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-163')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-163', '160')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '160', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -111 -> 128 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-111')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-111', '128')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '128', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -111 -> 93 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-111')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-111', '93')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '93', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -111 -> 50 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-111')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-111', '50')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '50', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -158 -> 88 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-158')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-158', '88')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '88', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -158 -> 93 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-158')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-158', '93')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '93', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -158 -> 50 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-158')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-158', '50')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '50', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 8_pl -> 128 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_pl')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_pl', '128')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '128', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 8_pl -> -158 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_pl')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_pl', '-158')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-158', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 8_pl -> 160 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_pl')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_pl', '160')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '160', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 177 -> sg_fp -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '177')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '177', 'sg_fp')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'sg_fp', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 9_sg -> 160 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sg')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sg', '160')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '160', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -123 -> -158 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-123')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-123', '-158')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-158', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 178 -> 134 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '178')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '178', '134')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '134', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 178 -> sg_fp -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '178')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '178', 'sg_fp')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'sg_fp', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 178 -> 128 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '178')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '178', '128')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '128', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 178 -> -46 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '178')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '178', '-46')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-46', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 178 -> 5_fp -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '178')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '178', '5_fp')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '5_fp', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 178 -> 77 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '178')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '178', '77')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '77', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 178 -> 160 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '178')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '178', '160')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '160', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 178 -> 8_fp -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '178')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '178', '8_fp')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '8_fp', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 178 -> 50 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '178')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '178', '50')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '50', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 172 -> sg_fp -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '172')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '172', 'sg_fp')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'sg_fp', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 172 -> 93 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '172')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '172', '93')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '93', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> nom. adj. -> sg_fp -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. adj.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. adj.', 'sg_fp')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'sg_fp', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -157 -> 128 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '128')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '128', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -86 -> 50 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-86')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-86', '50')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '50', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 13_sp -> 134 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_sp', '134')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '134', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 13_sp -> 128 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_sp', '128')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '128', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 13_sp -> 5_fp -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_sp', '5_fp')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '5_fp', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 13_sp -> 77 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_sp', '77')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '77', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 13_sp -> 160 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_sp', '160')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '160', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 13_sp -> 50 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_sp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_sp', '50')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '50', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 157 -> -158 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '157')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '157', '-158')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-158', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 157 -> 88 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '157')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '157', '88')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '88', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 157 -> 93 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '157')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '157', '93')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '93', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -230 -> 160 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '160')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '160', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 5_fp -> sg_fp -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', 'sg_fp')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'sg_fp', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 5_fp -> -158 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', '-158')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-158', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 88 -> sg_fp -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '88')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '88', 'sg_fp')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'sg_fp', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 88 -> 128 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '88')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '88', '128')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '128', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 88 -> -46 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '88')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '88', '-46')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-46', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 88 -> 50 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '88')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '88', '50')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '50', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 77 -> 88 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '77')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '77', '88')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '88', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 130 -> 93 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '130')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '130', '93')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '93', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -16 -> 160 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-16')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-16', '160')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '160', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -61 -> 160 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-61')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-61', '160')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '160', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -276 -> 77 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-276')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-276', '77')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '77', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -276 -> 160 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-276')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-276', '160')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '160', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -97 -> 128 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-97')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-97', '128')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '128', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 81 -> -158 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '81')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '81', '-158')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-158', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 4_sg -> sg_fp -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_sg')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_sg', 'sg_fp')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'sg_fp', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 40 -> 134 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '40')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '40', '134')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '134', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 40 -> sg_fp -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '40')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '40', 'sg_fp')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'sg_fp', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 40 -> 128 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '40')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '40', '128')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '128', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 40 -> -46 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '40')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '40', '-46')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-46', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 40 -> 5_fp -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '40')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '40', '5_fp')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '5_fp', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 40 -> 77 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '40')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '40', '77')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '77', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 40 -> 130 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '40')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '40', '130')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '130', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 40 -> 160 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '40')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '40', '160')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '160', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 40 -> 8_fp -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '40')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '40', '8_fp')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '8_fp', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 40 -> 50 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '40')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '40', '50')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '50', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -283 -> 128 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-283')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-283', '128')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '128', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -283 -> 160 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-283')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-283', '160')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '160', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -283 -> 50 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-283')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-283', '50')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '50', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 155 -> 50 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '155')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '155', '50')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '50', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 160 -> 151 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '160')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '160', '151')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '151', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 160 -> sg_fp -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '160')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '160', 'sg_fp')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'sg_fp', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 160 -> 128 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '160')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '160', '128')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '128', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 160 -> -46 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '160')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '160', '-46')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-46', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 160 -> -158 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '160')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '160', '-158')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-158', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 160 -> 5_fp -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '160')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '160', '5_fp')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '5_fp', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 160 -> 77 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '160')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '160', '77')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '77', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 160 -> 50 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '160')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '160', '50')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '50', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 93 -> 134 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '93')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '93', '134')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '134', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 93 -> 128 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '93')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '93', '128')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '128', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 93 -> 160 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '93')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '93', '160')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '160', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 93 -> 93 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '93')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '93', '93')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '93', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 50 -> -158 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '50')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '50', '-158')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-158', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 50 -> 93 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '50')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '50', '93')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '93', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 50 -> 8_fp -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '50')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '50', '8_fp')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '8_fp', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 50 -> 50 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '50')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '50', '50')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '50', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 139 -> 160 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '139')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '139', '160')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '160', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -271 -> 134 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-271')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-271', '134')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '134', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -271 -> 128 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-271')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-271', '128')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '128', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -271 -> 5_fp -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-271')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-271', '5_fp')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '5_fp', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -271 -> 77 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-271')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-271', '77')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '77', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -271 -> 160 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-271')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-271', '160')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '160', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -271 -> 50 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-271')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-271', '50')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '50', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 170 -> sg_fp -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '170')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '170', 'sg_fp')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'sg_fp', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> 170 -> -158 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '170')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '170', '-158')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-158', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -166 -> sg_fp -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-166')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-166', 'sg_fp')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'sg_fp', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -43 -> 128 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-43')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-43', '128')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '128', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -43 -> -158 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-43')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-43', '-158')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-158', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -69 -> 134 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-69')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-69', '134')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '134', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -69 -> 5_fp -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-69')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-69', '5_fp')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '5_fp', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -69 -> 77 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-69')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-69', '77')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '77', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # L -> -69 -> 160 -> T
    pleft = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-69')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-69', '160')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '160', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 151 -> 77 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '151')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '151', '77')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '77', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> gen. sg. -> 160 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'gen. sg.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. sg.', '160')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> sg_fp -> 77 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'sg_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg_fp', '77')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '77', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> sg_fp -> 160 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'sg_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg_fp', '160')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 128 -> 134 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '128')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '128', '134')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '134', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 128 -> sg_fp -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '128')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '128', 'sg_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 128 -> 128 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '128')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '128', '128')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '128', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 128 -> -46 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '128')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '128', '-46')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-46', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 128 -> 77 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '128')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '128', '77')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '77', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 128 -> 160 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '128')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '128', '160')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 128 -> 50 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '128')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '128', '50')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '50', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> -46 -> sg_fp -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-46')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-46', 'sg_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> -46 -> 128 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-46')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-46', '128')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '128', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> -46 -> 160 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-46')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-46', '160')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> -158 -> 77 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-158')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-158', '77')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '77', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> -158 -> 160 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-158')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-158', '160')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 88 -> 77 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '88')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '88', '77')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '77', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 88 -> 160 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '88')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '88', '160')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 77 -> 134 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '77')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '77', '134')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '134', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 77 -> 128 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '77')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '77', '128')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '128', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 77 -> 5_fp -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '77')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '77', '5_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 77 -> 77 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '77')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '77', '77')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '77', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 77 -> 160 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '77')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '77', '160')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 160 -> 134 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '160')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '160', '134')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '134', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 160 -> sg_fp -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '160')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '160', 'sg_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 160 -> 128 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '160')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '160', '128')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '128', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 160 -> 5_fp -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '160')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '160', '5_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 160 -> 77 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '160')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '160', '77')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '77', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 160 -> 160 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '160')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '160', '160')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 160 -> 8_fp -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '160')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '160', '8_fp')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 160 -> 50 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '160')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '160', '50')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '50', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 93 -> 134 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '93')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '93', '134')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '134', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 93 -> 160 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '93')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '93', '160')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 50 -> 77 -> L
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '50')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '50', '77')
    pright = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '77', node2.lemma)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 151 -> 93 -> T
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '151')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '151', '93')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '93', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 13_fp -> -158 -> T
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '13_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_fp', '-158')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-158', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 28_tp -> 50 -> T
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '28_tp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28_tp', '50')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '50', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> sg_fp -> 128 -> T
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'sg_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg_fp', '128')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '128', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> sg_fp -> -158 -> T
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'sg_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg_fp', '-158')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-158', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 128 -> -158 -> T
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '128')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '128', '-158')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-158', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> loc. pl. -> sg_fp -> T
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'loc. pl.')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. pl.', 'sg_fp')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'sg_fp', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> -158 -> 93 -> T
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-158')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-158', '93')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '93', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 157 -> -158 -> T
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '157')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '157', '-158')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-158', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 5_fp -> 93 -> T
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '5_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', '93')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '93', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 5_fp -> 50 -> T
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '5_fp')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', '50')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '50', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 77 -> sg_fp -> T
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '77')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '77', 'sg_fp')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'sg_fp', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 77 -> -46 -> T
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '77')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '77', '-46')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-46', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 77 -> -158 -> T
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '77')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '77', '-158')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-158', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> -97 -> -158 -> T
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-97')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-97', '-158')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-158', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 160 -> 128 -> T
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '160')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '160', '128')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '128', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 160 -> 160 -> T
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '160')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '160', '160')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '160', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 160 -> 93 -> T
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '160')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '160', '93')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '93', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 93 -> 93 -> T
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '93')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '93', '93')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '93', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    # T -> 50 -> sg_fp -> T
    pleft = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '50')
    pmid = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '50', 'sg_fp')
    pright = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'sg_fp', node2.tup)
    feats[fIndex] = pleft * pmid * pright; fIndex += 1
    
    feats[feats < 1e-25] = 1e-25
    return -np.log10(feats)

## END OF SUPER LARGE GET_FEATURS FUNCTION
    