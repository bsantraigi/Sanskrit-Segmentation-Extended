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
prob_lem2lem_linear = None
prob_lem2cng_linear = None
prob_lem2tup_linear = None

prob_cng2lem_linear = None
prob_cng2cng_linear = None
prob_cng2tup_linear = None

prob_tup2lem_linear = None
prob_tup2cng_linear = None
prob_tup2tup_linear = None

_cg_count = None
_edge_vector_dim = None
_full_cnglist = None

def word_definite_extInit(matDB):
    global prob_lem2lem_linear, prob_lem2cng_linear, prob_lem2tup_linear, \
     prob_cng2lem_linear, prob_cng2cng_linear, prob_cng2tup_linear, \
     prob_tup2lem_linear, prob_tup2cng_linear, prob_tup2tup_linear, \
     mat_lemCount_1D, mat_cngCount_1D, mat_tupCount_1D, \
     _cg_count, _full_cnglist, _edge_vector_dim

    prob_lem2lem_linear = matDB.prob_lem2lem_linear
    prob_lem2cng_linear = matDB.prob_lem2cng_linear
    prob_lem2tup_linear = matDB.prob_lem2tup_linear

    prob_cng2lem_linear = matDB.prob_cng2lem_linear
    prob_cng2cng_linear = matDB.prob_cng2cng_linear
    prob_cng2tup_linear = matDB.prob_cng2tup_linear
    
    prob_tup2lem_linear = matDB.prob_tup2lem_linear
    prob_tup2cng_linear = matDB.prob_tup2cng_linear
    prob_tup2tup_linear = matDB.prob_tup2tup_linear

    # TODO: Change to actual value
    # _cg_count = len(mat_cngCount_1D)
    # _edge_vector_dim = 9*_cg_count**2 + 9 * _cg_count + 8
    # _edge_vector_dim = 4*_cg_count**2 + 9*_cg_count + 9
    _edge_vector_dim = 1000
    # _full_cnglist = list(mat_cngCount_1D)
    # print(_cg_count)

"""
################################################################################################
######################  CREATING FEATURE MATRICES ##############################################
################################################################################################
"""
def tryProb_catchZero(mat2D, key1):
    # try:
    #     v = float(mat2D[key1][key2])/mat1D[key1];
    # except KeyError:
    #     v = 0

    if key1 in mat2D:
        return mat2D[key1]

    return 0
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
    feats[0] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' +  node2.lemma);
            
    # L->-90->L
    feats[1] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-90') * tryProb_catchZero(prob_cng2lem_linear, '-90^' +  node2.lemma);
            
    # L->-57->L
    feats[2] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-57') * tryProb_catchZero(prob_cng2lem_linear, '-57^' +  node2.lemma);
            
    # L->8_fp->L
    feats[3] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' +  node2.lemma);
            
    # L->169->L
    feats[4] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^169') * tryProb_catchZero(prob_cng2lem_linear, '169^' +  node2.lemma);
            
    # L->54->L
    feats[5] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' +  node2.lemma);
            
    # L->153->L
    feats[6] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' +  node2.lemma);
            
    # L->-24->L
    feats[7] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-24') * tryProb_catchZero(prob_cng2lem_linear, '-24^' +  node2.lemma);
            
    # L->137->L
    feats[8] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' +  node2.lemma);
            
    # L->-112->L
    feats[9] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' +  node2.lemma);
            
    # L->101->L
    feats[10] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' +  node2.lemma);
            
    # L->91->L
    feats[11] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^91') * tryProb_catchZero(prob_cng2lem_linear, '91^' +  node2.lemma);
            
    # L->3_du->L
    feats[12] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' +  node2.lemma);
            
    # L->91->T
    feats[13] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^91') * tryProb_catchZero(prob_cng2tup_linear, '91^' +  node2.tup);
            
    # C->3_du->T
    feats[14] = tryProb_catchZero(prob_cng2cng_linear, node1.cng + '^3_du') * tryProb_catchZero(prob_cng2tup_linear, '3_du^' +  node2.tup);
            
    # T->8_fp->L
    feats[15] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' +  node2.lemma);
            
    # T->54->L
    feats[16] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' +  node2.lemma);
            
    # T->153->L
    feats[17] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' +  node2.lemma);
            
    # T->137->L
    feats[18] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' +  node2.lemma);
            
    # T->-112->L
    feats[19] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' +  node2.lemma);
            
    # T->101->L
    feats[20] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' +  node2.lemma);
            
    # T->91->L
    feats[21] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^91') * tryProb_catchZero(prob_cng2lem_linear, '91^' +  node2.lemma);
            
    # L -> 2_sp -> dat. du. -> L
    feats[22] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^2_sp') * tryProb_catchZero(prob_cng2cng_linear, '2_sp^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # L -> 2_sp -> 8_fp -> L
    feats[23] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^2_sp') * tryProb_catchZero(prob_cng2cng_linear, '2_sp^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # L -> 2_sp -> 54 -> L
    feats[24] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^2_sp') * tryProb_catchZero(prob_cng2cng_linear, '2_sp^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> 2_sp -> pl -> L
    feats[25] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^2_sp') * tryProb_catchZero(prob_cng2cng_linear, '2_sp^pl') * tryProb_catchZero(prob_cng2lem_linear, 'pl^' + node2.lemma);
    
    # L -> 2_sp -> 28 -> L
    feats[26] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^2_sp') * tryProb_catchZero(prob_cng2cng_linear, '2_sp^28') * tryProb_catchZero(prob_cng2lem_linear, '28^' + node2.lemma);
    
    # L -> 2_sp -> 137 -> L
    feats[27] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^2_sp') * tryProb_catchZero(prob_cng2cng_linear, '2_sp^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> 2_sp -> 101 -> L
    feats[28] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^2_sp') * tryProb_catchZero(prob_cng2cng_linear, '2_sp^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> 2_sp -> 3_du -> L
    feats[29] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^2_sp') * tryProb_catchZero(prob_cng2cng_linear, '2_sp^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> 5_du -> dat. du. -> L
    feats[30] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^5_du') * tryProb_catchZero(prob_cng2cng_linear, '5_du^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # L -> 5_du -> -57 -> L
    feats[31] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^5_du') * tryProb_catchZero(prob_cng2cng_linear, '5_du^-57') * tryProb_catchZero(prob_cng2lem_linear, '-57^' + node2.lemma);
    
    # L -> 5_du -> 8_fp -> L
    feats[32] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^5_du') * tryProb_catchZero(prob_cng2cng_linear, '5_du^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # L -> 5_du -> 54 -> L
    feats[33] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^5_du') * tryProb_catchZero(prob_cng2cng_linear, '5_du^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> 5_du -> pl -> L
    feats[34] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^5_du') * tryProb_catchZero(prob_cng2cng_linear, '5_du^pl') * tryProb_catchZero(prob_cng2lem_linear, 'pl^' + node2.lemma);
    
    # L -> 5_du -> 153 -> L
    feats[35] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^5_du') * tryProb_catchZero(prob_cng2cng_linear, '5_du^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # L -> 5_du -> 28 -> L
    feats[36] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^5_du') * tryProb_catchZero(prob_cng2cng_linear, '5_du^28') * tryProb_catchZero(prob_cng2lem_linear, '28^' + node2.lemma);
    
    # L -> 5_du -> -24 -> L
    feats[37] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^5_du') * tryProb_catchZero(prob_cng2cng_linear, '5_du^-24') * tryProb_catchZero(prob_cng2lem_linear, '-24^' + node2.lemma);
    
    # L -> 5_du -> 137 -> L
    feats[38] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^5_du') * tryProb_catchZero(prob_cng2cng_linear, '5_du^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> 5_du -> -112 -> L
    feats[39] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^5_du') * tryProb_catchZero(prob_cng2cng_linear, '5_du^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # L -> 5_du -> 101 -> L
    feats[40] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^5_du') * tryProb_catchZero(prob_cng2cng_linear, '5_du^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> 5_du -> 91 -> L
    feats[41] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^5_du') * tryProb_catchZero(prob_cng2cng_linear, '5_du^91') * tryProb_catchZero(prob_cng2lem_linear, '91^' + node2.lemma);
    
    # L -> 5_du -> -35 -> L
    feats[42] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^5_du') * tryProb_catchZero(prob_cng2cng_linear, '5_du^-35') * tryProb_catchZero(prob_cng2lem_linear, '-35^' + node2.lemma);
    
    # L -> 5_du -> 3_du -> L
    feats[43] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^5_du') * tryProb_catchZero(prob_cng2cng_linear, '5_du^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> gen. sg. -> dat. du. -> L
    feats[44] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^gen. sg.') * tryProb_catchZero(prob_cng2cng_linear, 'gen. sg.^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # L -> gen. sg. -> -90 -> L
    feats[45] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^gen. sg.') * tryProb_catchZero(prob_cng2cng_linear, 'gen. sg.^-90') * tryProb_catchZero(prob_cng2lem_linear, '-90^' + node2.lemma);
    
    # L -> gen. sg. -> -57 -> L
    feats[46] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^gen. sg.') * tryProb_catchZero(prob_cng2cng_linear, 'gen. sg.^-57') * tryProb_catchZero(prob_cng2lem_linear, '-57^' + node2.lemma);
    
    # L -> gen. sg. -> 8_fp -> L
    feats[47] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^gen. sg.') * tryProb_catchZero(prob_cng2cng_linear, 'gen. sg.^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # L -> gen. sg. -> 54 -> L
    feats[48] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^gen. sg.') * tryProb_catchZero(prob_cng2cng_linear, 'gen. sg.^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> gen. sg. -> 153 -> L
    feats[49] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^gen. sg.') * tryProb_catchZero(prob_cng2cng_linear, 'gen. sg.^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # L -> gen. sg. -> -24 -> L
    feats[50] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^gen. sg.') * tryProb_catchZero(prob_cng2cng_linear, 'gen. sg.^-24') * tryProb_catchZero(prob_cng2lem_linear, '-24^' + node2.lemma);
    
    # L -> gen. sg. -> 137 -> L
    feats[51] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^gen. sg.') * tryProb_catchZero(prob_cng2cng_linear, 'gen. sg.^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> gen. sg. -> -112 -> L
    feats[52] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^gen. sg.') * tryProb_catchZero(prob_cng2cng_linear, 'gen. sg.^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # L -> gen. sg. -> 101 -> L
    feats[53] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^gen. sg.') * tryProb_catchZero(prob_cng2cng_linear, 'gen. sg.^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> gen. sg. -> 91 -> L
    feats[54] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^gen. sg.') * tryProb_catchZero(prob_cng2cng_linear, 'gen. sg.^91') * tryProb_catchZero(prob_cng2lem_linear, '91^' + node2.lemma);
    
    # L -> gen. sg. -> 3_du -> L
    feats[55] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^gen. sg.') * tryProb_catchZero(prob_cng2cng_linear, 'gen. sg.^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> 15_sp -> 101 -> L
    feats[56] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^15_sp') * tryProb_catchZero(prob_cng2cng_linear, '15_sp^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> 5_sp -> dat. du. -> L
    feats[57] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^5_sp') * tryProb_catchZero(prob_cng2cng_linear, '5_sp^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # L -> 5_sp -> -151 -> L
    feats[58] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^5_sp') * tryProb_catchZero(prob_cng2cng_linear, '5_sp^-151') * tryProb_catchZero(prob_cng2lem_linear, '-151^' + node2.lemma);
    
    # L -> 5_sp -> -57 -> L
    feats[59] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^5_sp') * tryProb_catchZero(prob_cng2cng_linear, '5_sp^-57') * tryProb_catchZero(prob_cng2lem_linear, '-57^' + node2.lemma);
    
    # L -> 5_sp -> 8_fp -> L
    feats[60] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^5_sp') * tryProb_catchZero(prob_cng2cng_linear, '5_sp^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # L -> 5_sp -> 169 -> L
    feats[61] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^5_sp') * tryProb_catchZero(prob_cng2cng_linear, '5_sp^169') * tryProb_catchZero(prob_cng2lem_linear, '169^' + node2.lemma);
    
    # L -> 5_sp -> pl -> L
    feats[62] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^5_sp') * tryProb_catchZero(prob_cng2cng_linear, '5_sp^pl') * tryProb_catchZero(prob_cng2lem_linear, 'pl^' + node2.lemma);
    
    # L -> 5_sp -> 153 -> L
    feats[63] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^5_sp') * tryProb_catchZero(prob_cng2cng_linear, '5_sp^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # L -> 5_sp -> -24 -> L
    feats[64] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^5_sp') * tryProb_catchZero(prob_cng2cng_linear, '5_sp^-24') * tryProb_catchZero(prob_cng2lem_linear, '-24^' + node2.lemma);
    
    # L -> 5_sp -> 137 -> L
    feats[65] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^5_sp') * tryProb_catchZero(prob_cng2cng_linear, '5_sp^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> 5_sp -> -112 -> L
    feats[66] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^5_sp') * tryProb_catchZero(prob_cng2cng_linear, '5_sp^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # L -> 5_sp -> 101 -> L
    feats[67] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^5_sp') * tryProb_catchZero(prob_cng2cng_linear, '5_sp^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> 5_sp -> 91 -> L
    feats[68] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^5_sp') * tryProb_catchZero(prob_cng2cng_linear, '5_sp^91') * tryProb_catchZero(prob_cng2lem_linear, '91^' + node2.lemma);
    
    # L -> 5_sp -> 3_du -> L
    feats[69] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^5_sp') * tryProb_catchZero(prob_cng2cng_linear, '5_sp^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> dat. du. -> 5_sp -> L
    feats[70] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^dat. du.') * tryProb_catchZero(prob_cng2cng_linear, 'dat. du.^5_sp') * tryProb_catchZero(prob_cng2lem_linear, '5_sp^' + node2.lemma);
    
    # L -> dat. du. -> dat. du. -> L
    feats[71] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^dat. du.') * tryProb_catchZero(prob_cng2cng_linear, 'dat. du.^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # L -> dat. du. -> -151 -> L
    feats[72] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^dat. du.') * tryProb_catchZero(prob_cng2cng_linear, 'dat. du.^-151') * tryProb_catchZero(prob_cng2lem_linear, '-151^' + node2.lemma);
    
    # L -> dat. du. -> -90 -> L
    feats[73] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^dat. du.') * tryProb_catchZero(prob_cng2cng_linear, 'dat. du.^-90') * tryProb_catchZero(prob_cng2lem_linear, '-90^' + node2.lemma);
    
    # L -> dat. du. -> 8_fp -> L
    feats[74] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^dat. du.') * tryProb_catchZero(prob_cng2cng_linear, 'dat. du.^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # L -> dat. du. -> 54 -> L
    feats[75] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^dat. du.') * tryProb_catchZero(prob_cng2cng_linear, 'dat. du.^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> dat. du. -> pl -> L
    feats[76] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^dat. du.') * tryProb_catchZero(prob_cng2cng_linear, 'dat. du.^pl') * tryProb_catchZero(prob_cng2lem_linear, 'pl^' + node2.lemma);
    
    # L -> dat. du. -> 153 -> L
    feats[77] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^dat. du.') * tryProb_catchZero(prob_cng2cng_linear, 'dat. du.^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # L -> dat. du. -> -24 -> L
    feats[78] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^dat. du.') * tryProb_catchZero(prob_cng2cng_linear, 'dat. du.^-24') * tryProb_catchZero(prob_cng2lem_linear, '-24^' + node2.lemma);
    
    # L -> dat. du. -> 137 -> L
    feats[79] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^dat. du.') * tryProb_catchZero(prob_cng2cng_linear, 'dat. du.^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> dat. du. -> -112 -> L
    feats[80] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^dat. du.') * tryProb_catchZero(prob_cng2cng_linear, 'dat. du.^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # L -> dat. du. -> 101 -> L
    feats[81] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^dat. du.') * tryProb_catchZero(prob_cng2cng_linear, 'dat. du.^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> dat. du. -> -35 -> L
    feats[82] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^dat. du.') * tryProb_catchZero(prob_cng2cng_linear, 'dat. du.^-35') * tryProb_catchZero(prob_cng2lem_linear, '-35^' + node2.lemma);
    
    # L -> dat. du. -> 3_du -> L
    feats[83] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^dat. du.') * tryProb_catchZero(prob_cng2cng_linear, 'dat. du.^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> -37 -> 54 -> L
    feats[84] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-37') * tryProb_catchZero(prob_cng2cng_linear, '-37^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> -37 -> 153 -> L
    feats[85] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-37') * tryProb_catchZero(prob_cng2cng_linear, '-37^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # L -> -37 -> 137 -> L
    feats[86] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-37') * tryProb_catchZero(prob_cng2cng_linear, '-37^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> -37 -> -112 -> L
    feats[87] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-37') * tryProb_catchZero(prob_cng2cng_linear, '-37^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # L -> -37 -> 101 -> L
    feats[88] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-37') * tryProb_catchZero(prob_cng2cng_linear, '-37^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> 15_du -> dat. du. -> L
    feats[89] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^15_du') * tryProb_catchZero(prob_cng2cng_linear, '15_du^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # L -> 15_du -> -151 -> L
    feats[90] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^15_du') * tryProb_catchZero(prob_cng2cng_linear, '15_du^-151') * tryProb_catchZero(prob_cng2lem_linear, '-151^' + node2.lemma);
    
    # L -> 15_du -> -90 -> L
    feats[91] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^15_du') * tryProb_catchZero(prob_cng2cng_linear, '15_du^-90') * tryProb_catchZero(prob_cng2lem_linear, '-90^' + node2.lemma);
    
    # L -> 15_du -> -57 -> L
    feats[92] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^15_du') * tryProb_catchZero(prob_cng2cng_linear, '15_du^-57') * tryProb_catchZero(prob_cng2lem_linear, '-57^' + node2.lemma);
    
    # L -> 15_du -> 8_fp -> L
    feats[93] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^15_du') * tryProb_catchZero(prob_cng2cng_linear, '15_du^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # L -> 15_du -> 169 -> L
    feats[94] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^15_du') * tryProb_catchZero(prob_cng2cng_linear, '15_du^169') * tryProb_catchZero(prob_cng2lem_linear, '169^' + node2.lemma);
    
    # L -> 15_du -> 54 -> L
    feats[95] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^15_du') * tryProb_catchZero(prob_cng2cng_linear, '15_du^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> 15_du -> pl -> L
    feats[96] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^15_du') * tryProb_catchZero(prob_cng2cng_linear, '15_du^pl') * tryProb_catchZero(prob_cng2lem_linear, 'pl^' + node2.lemma);
    
    # L -> 15_du -> 153 -> L
    feats[97] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^15_du') * tryProb_catchZero(prob_cng2cng_linear, '15_du^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # L -> 15_du -> 28 -> L
    feats[98] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^15_du') * tryProb_catchZero(prob_cng2cng_linear, '15_du^28') * tryProb_catchZero(prob_cng2lem_linear, '28^' + node2.lemma);
    
    # L -> 15_du -> -24 -> L
    feats[99] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^15_du') * tryProb_catchZero(prob_cng2cng_linear, '15_du^-24') * tryProb_catchZero(prob_cng2lem_linear, '-24^' + node2.lemma);
    
    # L -> 15_du -> 137 -> L
    feats[100] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^15_du') * tryProb_catchZero(prob_cng2cng_linear, '15_du^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> 15_du -> -112 -> L
    feats[101] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^15_du') * tryProb_catchZero(prob_cng2cng_linear, '15_du^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # L -> 15_du -> 101 -> L
    feats[102] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^15_du') * tryProb_catchZero(prob_cng2cng_linear, '15_du^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> 15_du -> 91 -> L
    feats[103] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^15_du') * tryProb_catchZero(prob_cng2cng_linear, '15_du^91') * tryProb_catchZero(prob_cng2lem_linear, '91^' + node2.lemma);
    
    # L -> 15_du -> -35 -> L
    feats[104] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^15_du') * tryProb_catchZero(prob_cng2cng_linear, '15_du^-35') * tryProb_catchZero(prob_cng2lem_linear, '-35^' + node2.lemma);
    
    # L -> 15_du -> 3_du -> L
    feats[105] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^15_du') * tryProb_catchZero(prob_cng2cng_linear, '15_du^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> -151 -> dat. du. -> L
    feats[106] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-151') * tryProb_catchZero(prob_cng2cng_linear, '-151^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # L -> -151 -> -90 -> L
    feats[107] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-151') * tryProb_catchZero(prob_cng2cng_linear, '-151^-90') * tryProb_catchZero(prob_cng2lem_linear, '-90^' + node2.lemma);
    
    # L -> -151 -> 8_fp -> L
    feats[108] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-151') * tryProb_catchZero(prob_cng2cng_linear, '-151^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # L -> -151 -> 169 -> L
    feats[109] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-151') * tryProb_catchZero(prob_cng2cng_linear, '-151^169') * tryProb_catchZero(prob_cng2lem_linear, '169^' + node2.lemma);
    
    # L -> -151 -> 54 -> L
    feats[110] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-151') * tryProb_catchZero(prob_cng2cng_linear, '-151^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> -151 -> pl -> L
    feats[111] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-151') * tryProb_catchZero(prob_cng2cng_linear, '-151^pl') * tryProb_catchZero(prob_cng2lem_linear, 'pl^' + node2.lemma);
    
    # L -> -151 -> 153 -> L
    feats[112] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-151') * tryProb_catchZero(prob_cng2cng_linear, '-151^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # L -> -151 -> 28 -> L
    feats[113] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-151') * tryProb_catchZero(prob_cng2cng_linear, '-151^28') * tryProb_catchZero(prob_cng2lem_linear, '28^' + node2.lemma);
    
    # L -> -151 -> -24 -> L
    feats[114] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-151') * tryProb_catchZero(prob_cng2cng_linear, '-151^-24') * tryProb_catchZero(prob_cng2lem_linear, '-24^' + node2.lemma);
    
    # L -> -151 -> 137 -> L
    feats[115] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-151') * tryProb_catchZero(prob_cng2cng_linear, '-151^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> -151 -> -112 -> L
    feats[116] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-151') * tryProb_catchZero(prob_cng2cng_linear, '-151^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # L -> -151 -> 101 -> L
    feats[117] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-151') * tryProb_catchZero(prob_cng2cng_linear, '-151^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> -151 -> 91 -> L
    feats[118] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-151') * tryProb_catchZero(prob_cng2cng_linear, '-151^91') * tryProb_catchZero(prob_cng2lem_linear, '91^' + node2.lemma);
    
    # L -> -151 -> -35 -> L
    feats[119] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-151') * tryProb_catchZero(prob_cng2cng_linear, '-151^-35') * tryProb_catchZero(prob_cng2lem_linear, '-35^' + node2.lemma);
    
    # L -> -151 -> 3_du -> L
    feats[120] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-151') * tryProb_catchZero(prob_cng2cng_linear, '-151^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> 39 -> dat. du. -> L
    feats[121] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^39') * tryProb_catchZero(prob_cng2cng_linear, '39^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # L -> 39 -> 8_fp -> L
    feats[122] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^39') * tryProb_catchZero(prob_cng2cng_linear, '39^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # L -> 39 -> 54 -> L
    feats[123] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^39') * tryProb_catchZero(prob_cng2cng_linear, '39^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> 39 -> 153 -> L
    feats[124] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^39') * tryProb_catchZero(prob_cng2cng_linear, '39^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # L -> 39 -> 137 -> L
    feats[125] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^39') * tryProb_catchZero(prob_cng2cng_linear, '39^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> 39 -> -112 -> L
    feats[126] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^39') * tryProb_catchZero(prob_cng2cng_linear, '39^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # L -> 39 -> 101 -> L
    feats[127] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^39') * tryProb_catchZero(prob_cng2cng_linear, '39^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> 39 -> 91 -> L
    feats[128] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^39') * tryProb_catchZero(prob_cng2cng_linear, '39^91') * tryProb_catchZero(prob_cng2lem_linear, '91^' + node2.lemma);
    
    # L -> 39 -> 3_du -> L
    feats[129] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^39') * tryProb_catchZero(prob_cng2cng_linear, '39^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> 12_fp -> 101 -> L
    feats[130] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^12_fp') * tryProb_catchZero(prob_cng2cng_linear, '12_fp^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> -306 -> 5_sp -> L
    feats[131] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-306') * tryProb_catchZero(prob_cng2cng_linear, '-306^5_sp') * tryProb_catchZero(prob_cng2lem_linear, '5_sp^' + node2.lemma);
    
    # L -> -306 -> dat. du. -> L
    feats[132] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-306') * tryProb_catchZero(prob_cng2cng_linear, '-306^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # L -> -306 -> -151 -> L
    feats[133] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-306') * tryProb_catchZero(prob_cng2cng_linear, '-306^-151') * tryProb_catchZero(prob_cng2lem_linear, '-151^' + node2.lemma);
    
    # L -> -306 -> -90 -> L
    feats[134] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-306') * tryProb_catchZero(prob_cng2cng_linear, '-306^-90') * tryProb_catchZero(prob_cng2lem_linear, '-90^' + node2.lemma);
    
    # L -> -306 -> -57 -> L
    feats[135] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-306') * tryProb_catchZero(prob_cng2cng_linear, '-306^-57') * tryProb_catchZero(prob_cng2lem_linear, '-57^' + node2.lemma);
    
    # L -> -306 -> 8_fp -> L
    feats[136] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-306') * tryProb_catchZero(prob_cng2cng_linear, '-306^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # L -> -306 -> 169 -> L
    feats[137] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-306') * tryProb_catchZero(prob_cng2cng_linear, '-306^169') * tryProb_catchZero(prob_cng2lem_linear, '169^' + node2.lemma);
    
    # L -> -306 -> 54 -> L
    feats[138] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-306') * tryProb_catchZero(prob_cng2cng_linear, '-306^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> -306 -> pl -> L
    feats[139] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-306') * tryProb_catchZero(prob_cng2cng_linear, '-306^pl') * tryProb_catchZero(prob_cng2lem_linear, 'pl^' + node2.lemma);
    
    # L -> -306 -> 153 -> L
    feats[140] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-306') * tryProb_catchZero(prob_cng2cng_linear, '-306^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # L -> -306 -> 28 -> L
    feats[141] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-306') * tryProb_catchZero(prob_cng2cng_linear, '-306^28') * tryProb_catchZero(prob_cng2lem_linear, '28^' + node2.lemma);
    
    # L -> -306 -> -24 -> L
    feats[142] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-306') * tryProb_catchZero(prob_cng2cng_linear, '-306^-24') * tryProb_catchZero(prob_cng2lem_linear, '-24^' + node2.lemma);
    
    # L -> -306 -> 137 -> L
    feats[143] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-306') * tryProb_catchZero(prob_cng2cng_linear, '-306^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> -306 -> -112 -> L
    feats[144] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-306') * tryProb_catchZero(prob_cng2cng_linear, '-306^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # L -> -306 -> 101 -> L
    feats[145] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-306') * tryProb_catchZero(prob_cng2cng_linear, '-306^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> -306 -> 91 -> L
    feats[146] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-306') * tryProb_catchZero(prob_cng2cng_linear, '-306^91') * tryProb_catchZero(prob_cng2lem_linear, '91^' + node2.lemma);
    
    # L -> -306 -> -35 -> L
    feats[147] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-306') * tryProb_catchZero(prob_cng2cng_linear, '-306^-35') * tryProb_catchZero(prob_cng2lem_linear, '-35^' + node2.lemma);
    
    # L -> -306 -> 3_du -> L
    feats[148] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-306') * tryProb_catchZero(prob_cng2cng_linear, '-306^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> -64 -> dat. du. -> L
    feats[149] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-64') * tryProb_catchZero(prob_cng2cng_linear, '-64^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # L -> -64 -> 54 -> L
    feats[150] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-64') * tryProb_catchZero(prob_cng2cng_linear, '-64^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> -64 -> 137 -> L
    feats[151] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-64') * tryProb_catchZero(prob_cng2cng_linear, '-64^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> -64 -> 101 -> L
    feats[152] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-64') * tryProb_catchZero(prob_cng2cng_linear, '-64^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> -64 -> 3_du -> L
    feats[153] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-64') * tryProb_catchZero(prob_cng2cng_linear, '-64^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> -90 -> -151 -> L
    feats[154] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-90') * tryProb_catchZero(prob_cng2cng_linear, '-90^-151') * tryProb_catchZero(prob_cng2lem_linear, '-151^' + node2.lemma);
    
    # L -> -90 -> -90 -> L
    feats[155] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-90') * tryProb_catchZero(prob_cng2cng_linear, '-90^-90') * tryProb_catchZero(prob_cng2lem_linear, '-90^' + node2.lemma);
    
    # L -> -90 -> 8_fp -> L
    feats[156] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-90') * tryProb_catchZero(prob_cng2cng_linear, '-90^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # L -> -90 -> 169 -> L
    feats[157] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-90') * tryProb_catchZero(prob_cng2cng_linear, '-90^169') * tryProb_catchZero(prob_cng2lem_linear, '169^' + node2.lemma);
    
    # L -> -90 -> 54 -> L
    feats[158] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-90') * tryProb_catchZero(prob_cng2cng_linear, '-90^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> -90 -> pl -> L
    feats[159] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-90') * tryProb_catchZero(prob_cng2cng_linear, '-90^pl') * tryProb_catchZero(prob_cng2lem_linear, 'pl^' + node2.lemma);
    
    # L -> -90 -> 153 -> L
    feats[160] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-90') * tryProb_catchZero(prob_cng2cng_linear, '-90^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # L -> -90 -> 28 -> L
    feats[161] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-90') * tryProb_catchZero(prob_cng2cng_linear, '-90^28') * tryProb_catchZero(prob_cng2lem_linear, '28^' + node2.lemma);
    
    # L -> -90 -> -24 -> L
    feats[162] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-90') * tryProb_catchZero(prob_cng2cng_linear, '-90^-24') * tryProb_catchZero(prob_cng2lem_linear, '-24^' + node2.lemma);
    
    # L -> -90 -> 137 -> L
    feats[163] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-90') * tryProb_catchZero(prob_cng2cng_linear, '-90^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> -90 -> -112 -> L
    feats[164] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-90') * tryProb_catchZero(prob_cng2cng_linear, '-90^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # L -> -90 -> 101 -> L
    feats[165] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-90') * tryProb_catchZero(prob_cng2cng_linear, '-90^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> -90 -> 91 -> L
    feats[166] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-90') * tryProb_catchZero(prob_cng2cng_linear, '-90^91') * tryProb_catchZero(prob_cng2lem_linear, '91^' + node2.lemma);
    
    # L -> -90 -> -35 -> L
    feats[167] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-90') * tryProb_catchZero(prob_cng2cng_linear, '-90^-35') * tryProb_catchZero(prob_cng2lem_linear, '-35^' + node2.lemma);
    
    # L -> -90 -> 3_du -> L
    feats[168] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-90') * tryProb_catchZero(prob_cng2cng_linear, '-90^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> -210 -> 8_fp -> L
    feats[169] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-210') * tryProb_catchZero(prob_cng2cng_linear, '-210^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # L -> -210 -> 54 -> L
    feats[170] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-210') * tryProb_catchZero(prob_cng2cng_linear, '-210^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> -210 -> 153 -> L
    feats[171] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-210') * tryProb_catchZero(prob_cng2cng_linear, '-210^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # L -> -210 -> -112 -> L
    feats[172] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-210') * tryProb_catchZero(prob_cng2cng_linear, '-210^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # L -> -210 -> 101 -> L
    feats[173] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-210') * tryProb_catchZero(prob_cng2cng_linear, '-210^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> -210 -> 91 -> L
    feats[174] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-210') * tryProb_catchZero(prob_cng2cng_linear, '-210^91') * tryProb_catchZero(prob_cng2lem_linear, '91^' + node2.lemma);
    
    # L -> -210 -> 3_du -> L
    feats[175] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-210') * tryProb_catchZero(prob_cng2cng_linear, '-210^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> -159 -> dat. du. -> L
    feats[176] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-159') * tryProb_catchZero(prob_cng2cng_linear, '-159^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # L -> -159 -> 8_fp -> L
    feats[177] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-159') * tryProb_catchZero(prob_cng2cng_linear, '-159^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # L -> -159 -> 54 -> L
    feats[178] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-159') * tryProb_catchZero(prob_cng2cng_linear, '-159^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> -159 -> pl -> L
    feats[179] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-159') * tryProb_catchZero(prob_cng2cng_linear, '-159^pl') * tryProb_catchZero(prob_cng2lem_linear, 'pl^' + node2.lemma);
    
    # L -> -159 -> 153 -> L
    feats[180] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-159') * tryProb_catchZero(prob_cng2cng_linear, '-159^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # L -> -159 -> -24 -> L
    feats[181] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-159') * tryProb_catchZero(prob_cng2cng_linear, '-159^-24') * tryProb_catchZero(prob_cng2lem_linear, '-24^' + node2.lemma);
    
    # L -> -159 -> 137 -> L
    feats[182] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-159') * tryProb_catchZero(prob_cng2cng_linear, '-159^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> -159 -> -112 -> L
    feats[183] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-159') * tryProb_catchZero(prob_cng2cng_linear, '-159^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # L -> -159 -> 101 -> L
    feats[184] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-159') * tryProb_catchZero(prob_cng2cng_linear, '-159^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> -159 -> 91 -> L
    feats[185] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-159') * tryProb_catchZero(prob_cng2cng_linear, '-159^91') * tryProb_catchZero(prob_cng2lem_linear, '91^' + node2.lemma);
    
    # L -> -159 -> 3_du -> L
    feats[186] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-159') * tryProb_catchZero(prob_cng2cng_linear, '-159^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> -57 -> dat. du. -> L
    feats[187] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-57') * tryProb_catchZero(prob_cng2cng_linear, '-57^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # L -> -57 -> 54 -> L
    feats[188] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-57') * tryProb_catchZero(prob_cng2cng_linear, '-57^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> -57 -> pl -> L
    feats[189] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-57') * tryProb_catchZero(prob_cng2cng_linear, '-57^pl') * tryProb_catchZero(prob_cng2lem_linear, 'pl^' + node2.lemma);
    
    # L -> -57 -> 153 -> L
    feats[190] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-57') * tryProb_catchZero(prob_cng2cng_linear, '-57^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # L -> -57 -> -24 -> L
    feats[191] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-57') * tryProb_catchZero(prob_cng2cng_linear, '-57^-24') * tryProb_catchZero(prob_cng2lem_linear, '-24^' + node2.lemma);
    
    # L -> -57 -> 137 -> L
    feats[192] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-57') * tryProb_catchZero(prob_cng2cng_linear, '-57^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> -57 -> -112 -> L
    feats[193] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-57') * tryProb_catchZero(prob_cng2cng_linear, '-57^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # L -> -57 -> 101 -> L
    feats[194] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-57') * tryProb_catchZero(prob_cng2cng_linear, '-57^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> -57 -> 91 -> L
    feats[195] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-57') * tryProb_catchZero(prob_cng2cng_linear, '-57^91') * tryProb_catchZero(prob_cng2lem_linear, '91^' + node2.lemma);
    
    # L -> -57 -> 3_du -> L
    feats[196] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-57') * tryProb_catchZero(prob_cng2cng_linear, '-57^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> 10_sp -> 5_sp -> L
    feats[197] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^10_sp') * tryProb_catchZero(prob_cng2cng_linear, '10_sp^5_sp') * tryProb_catchZero(prob_cng2lem_linear, '5_sp^' + node2.lemma);
    
    # L -> 10_sp -> dat. du. -> L
    feats[198] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^10_sp') * tryProb_catchZero(prob_cng2cng_linear, '10_sp^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # L -> 10_sp -> 8_fp -> L
    feats[199] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^10_sp') * tryProb_catchZero(prob_cng2cng_linear, '10_sp^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # L -> 10_sp -> 54 -> L
    feats[200] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^10_sp') * tryProb_catchZero(prob_cng2cng_linear, '10_sp^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> 10_sp -> pl -> L
    feats[201] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^10_sp') * tryProb_catchZero(prob_cng2cng_linear, '10_sp^pl') * tryProb_catchZero(prob_cng2lem_linear, 'pl^' + node2.lemma);
    
    # L -> 10_sp -> 153 -> L
    feats[202] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^10_sp') * tryProb_catchZero(prob_cng2cng_linear, '10_sp^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # L -> 10_sp -> 28 -> L
    feats[203] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^10_sp') * tryProb_catchZero(prob_cng2cng_linear, '10_sp^28') * tryProb_catchZero(prob_cng2lem_linear, '28^' + node2.lemma);
    
    # L -> 10_sp -> 137 -> L
    feats[204] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^10_sp') * tryProb_catchZero(prob_cng2cng_linear, '10_sp^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> 10_sp -> -112 -> L
    feats[205] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^10_sp') * tryProb_catchZero(prob_cng2cng_linear, '10_sp^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # L -> 10_sp -> 101 -> L
    feats[206] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^10_sp') * tryProb_catchZero(prob_cng2cng_linear, '10_sp^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> 10_sp -> -35 -> L
    feats[207] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^10_sp') * tryProb_catchZero(prob_cng2cng_linear, '10_sp^-35') * tryProb_catchZero(prob_cng2lem_linear, '-35^' + node2.lemma);
    
    # L -> 10_sp -> 3_du -> L
    feats[208] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^10_sp') * tryProb_catchZero(prob_cng2cng_linear, '10_sp^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> 6_sg -> dat. du. -> L
    feats[209] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^6_sg') * tryProb_catchZero(prob_cng2cng_linear, '6_sg^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # L -> 6_sg -> 8_fp -> L
    feats[210] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^6_sg') * tryProb_catchZero(prob_cng2cng_linear, '6_sg^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # L -> 6_sg -> 54 -> L
    feats[211] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^6_sg') * tryProb_catchZero(prob_cng2cng_linear, '6_sg^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> 6_sg -> pl -> L
    feats[212] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^6_sg') * tryProb_catchZero(prob_cng2cng_linear, '6_sg^pl') * tryProb_catchZero(prob_cng2lem_linear, 'pl^' + node2.lemma);
    
    # L -> 6_sg -> 28 -> L
    feats[213] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^6_sg') * tryProb_catchZero(prob_cng2cng_linear, '6_sg^28') * tryProb_catchZero(prob_cng2lem_linear, '28^' + node2.lemma);
    
    # L -> 6_sg -> 137 -> L
    feats[214] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^6_sg') * tryProb_catchZero(prob_cng2cng_linear, '6_sg^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> 6_sg -> 101 -> L
    feats[215] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^6_sg') * tryProb_catchZero(prob_cng2cng_linear, '6_sg^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> 6_sg -> 3_du -> L
    feats[216] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^6_sg') * tryProb_catchZero(prob_cng2cng_linear, '6_sg^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> 8_fp -> dat. du. -> L
    feats[217] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^8_fp') * tryProb_catchZero(prob_cng2cng_linear, '8_fp^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # L -> 8_fp -> -151 -> L
    feats[218] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^8_fp') * tryProb_catchZero(prob_cng2cng_linear, '8_fp^-151') * tryProb_catchZero(prob_cng2lem_linear, '-151^' + node2.lemma);
    
    # L -> 8_fp -> -57 -> L
    feats[219] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^8_fp') * tryProb_catchZero(prob_cng2cng_linear, '8_fp^-57') * tryProb_catchZero(prob_cng2lem_linear, '-57^' + node2.lemma);
    
    # L -> 8_fp -> 8_fp -> L
    feats[220] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^8_fp') * tryProb_catchZero(prob_cng2cng_linear, '8_fp^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # L -> 8_fp -> 169 -> L
    feats[221] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^8_fp') * tryProb_catchZero(prob_cng2cng_linear, '8_fp^169') * tryProb_catchZero(prob_cng2lem_linear, '169^' + node2.lemma);
    
    # L -> 8_fp -> 54 -> L
    feats[222] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^8_fp') * tryProb_catchZero(prob_cng2cng_linear, '8_fp^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> 8_fp -> pl -> L
    feats[223] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^8_fp') * tryProb_catchZero(prob_cng2cng_linear, '8_fp^pl') * tryProb_catchZero(prob_cng2lem_linear, 'pl^' + node2.lemma);
    
    # L -> 8_fp -> 153 -> L
    feats[224] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^8_fp') * tryProb_catchZero(prob_cng2cng_linear, '8_fp^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # L -> 8_fp -> 28 -> L
    feats[225] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^8_fp') * tryProb_catchZero(prob_cng2cng_linear, '8_fp^28') * tryProb_catchZero(prob_cng2lem_linear, '28^' + node2.lemma);
    
    # L -> 8_fp -> -24 -> L
    feats[226] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^8_fp') * tryProb_catchZero(prob_cng2cng_linear, '8_fp^-24') * tryProb_catchZero(prob_cng2lem_linear, '-24^' + node2.lemma);
    
    # L -> 8_fp -> 137 -> L
    feats[227] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^8_fp') * tryProb_catchZero(prob_cng2cng_linear, '8_fp^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> 8_fp -> -112 -> L
    feats[228] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^8_fp') * tryProb_catchZero(prob_cng2cng_linear, '8_fp^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # L -> 8_fp -> 101 -> L
    feats[229] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^8_fp') * tryProb_catchZero(prob_cng2cng_linear, '8_fp^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> 8_fp -> 91 -> L
    feats[230] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^8_fp') * tryProb_catchZero(prob_cng2cng_linear, '8_fp^91') * tryProb_catchZero(prob_cng2lem_linear, '91^' + node2.lemma);
    
    # L -> 8_fp -> -35 -> L
    feats[231] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^8_fp') * tryProb_catchZero(prob_cng2cng_linear, '8_fp^-35') * tryProb_catchZero(prob_cng2lem_linear, '-35^' + node2.lemma);
    
    # L -> 8_fp -> 3_du -> L
    feats[232] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^8_fp') * tryProb_catchZero(prob_cng2cng_linear, '8_fp^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> 4_pl -> dat. du. -> L
    feats[233] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^4_pl') * tryProb_catchZero(prob_cng2cng_linear, '4_pl^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # L -> 4_pl -> 8_fp -> L
    feats[234] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^4_pl') * tryProb_catchZero(prob_cng2cng_linear, '4_pl^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # L -> 4_pl -> 54 -> L
    feats[235] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^4_pl') * tryProb_catchZero(prob_cng2cng_linear, '4_pl^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> 4_pl -> pl -> L
    feats[236] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^4_pl') * tryProb_catchZero(prob_cng2cng_linear, '4_pl^pl') * tryProb_catchZero(prob_cng2lem_linear, 'pl^' + node2.lemma);
    
    # L -> 4_pl -> 153 -> L
    feats[237] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^4_pl') * tryProb_catchZero(prob_cng2cng_linear, '4_pl^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # L -> 4_pl -> 28 -> L
    feats[238] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^4_pl') * tryProb_catchZero(prob_cng2cng_linear, '4_pl^28') * tryProb_catchZero(prob_cng2lem_linear, '28^' + node2.lemma);
    
    # L -> 4_pl -> -24 -> L
    feats[239] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^4_pl') * tryProb_catchZero(prob_cng2cng_linear, '4_pl^-24') * tryProb_catchZero(prob_cng2lem_linear, '-24^' + node2.lemma);
    
    # L -> 4_pl -> 137 -> L
    feats[240] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^4_pl') * tryProb_catchZero(prob_cng2cng_linear, '4_pl^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> 4_pl -> -112 -> L
    feats[241] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^4_pl') * tryProb_catchZero(prob_cng2cng_linear, '4_pl^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # L -> 4_pl -> 101 -> L
    feats[242] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^4_pl') * tryProb_catchZero(prob_cng2cng_linear, '4_pl^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> 4_pl -> 91 -> L
    feats[243] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^4_pl') * tryProb_catchZero(prob_cng2cng_linear, '4_pl^91') * tryProb_catchZero(prob_cng2lem_linear, '91^' + node2.lemma);
    
    # L -> 4_pl -> -35 -> L
    feats[244] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^4_pl') * tryProb_catchZero(prob_cng2cng_linear, '4_pl^-35') * tryProb_catchZero(prob_cng2lem_linear, '-35^' + node2.lemma);
    
    # L -> 4_pl -> 3_du -> L
    feats[245] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^4_pl') * tryProb_catchZero(prob_cng2cng_linear, '4_pl^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> 169 -> dat. du. -> L
    feats[246] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^169') * tryProb_catchZero(prob_cng2cng_linear, '169^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # L -> 169 -> -151 -> L
    feats[247] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^169') * tryProb_catchZero(prob_cng2cng_linear, '169^-151') * tryProb_catchZero(prob_cng2lem_linear, '-151^' + node2.lemma);
    
    # L -> 169 -> -90 -> L
    feats[248] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^169') * tryProb_catchZero(prob_cng2cng_linear, '169^-90') * tryProb_catchZero(prob_cng2lem_linear, '-90^' + node2.lemma);
    
    # L -> 169 -> -57 -> L
    feats[249] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^169') * tryProb_catchZero(prob_cng2cng_linear, '169^-57') * tryProb_catchZero(prob_cng2lem_linear, '-57^' + node2.lemma);
    
    # L -> 169 -> 8_fp -> L
    feats[250] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^169') * tryProb_catchZero(prob_cng2cng_linear, '169^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # L -> 169 -> 169 -> L
    feats[251] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^169') * tryProb_catchZero(prob_cng2cng_linear, '169^169') * tryProb_catchZero(prob_cng2lem_linear, '169^' + node2.lemma);
    
    # L -> 169 -> 54 -> L
    feats[252] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^169') * tryProb_catchZero(prob_cng2cng_linear, '169^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> 169 -> 153 -> L
    feats[253] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^169') * tryProb_catchZero(prob_cng2cng_linear, '169^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # L -> 169 -> 28 -> L
    feats[254] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^169') * tryProb_catchZero(prob_cng2cng_linear, '169^28') * tryProb_catchZero(prob_cng2lem_linear, '28^' + node2.lemma);
    
    # L -> 169 -> -24 -> L
    feats[255] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^169') * tryProb_catchZero(prob_cng2cng_linear, '169^-24') * tryProb_catchZero(prob_cng2lem_linear, '-24^' + node2.lemma);
    
    # L -> 169 -> 137 -> L
    feats[256] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^169') * tryProb_catchZero(prob_cng2cng_linear, '169^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> 169 -> -112 -> L
    feats[257] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^169') * tryProb_catchZero(prob_cng2cng_linear, '169^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # L -> 169 -> 101 -> L
    feats[258] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^169') * tryProb_catchZero(prob_cng2cng_linear, '169^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> 169 -> 3_du -> L
    feats[259] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^169') * tryProb_catchZero(prob_cng2cng_linear, '169^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> 54 -> 5_sp -> L
    feats[260] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^54') * tryProb_catchZero(prob_cng2cng_linear, '54^5_sp') * tryProb_catchZero(prob_cng2lem_linear, '5_sp^' + node2.lemma);
    
    # L -> 54 -> dat. du. -> L
    feats[261] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^54') * tryProb_catchZero(prob_cng2cng_linear, '54^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # L -> 54 -> -151 -> L
    feats[262] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^54') * tryProb_catchZero(prob_cng2cng_linear, '54^-151') * tryProb_catchZero(prob_cng2lem_linear, '-151^' + node2.lemma);
    
    # L -> 54 -> -90 -> L
    feats[263] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^54') * tryProb_catchZero(prob_cng2cng_linear, '54^-90') * tryProb_catchZero(prob_cng2lem_linear, '-90^' + node2.lemma);
    
    # L -> 54 -> -57 -> L
    feats[264] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^54') * tryProb_catchZero(prob_cng2cng_linear, '54^-57') * tryProb_catchZero(prob_cng2lem_linear, '-57^' + node2.lemma);
    
    # L -> 54 -> 8_fp -> L
    feats[265] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^54') * tryProb_catchZero(prob_cng2cng_linear, '54^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # L -> 54 -> 169 -> L
    feats[266] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^54') * tryProb_catchZero(prob_cng2cng_linear, '54^169') * tryProb_catchZero(prob_cng2lem_linear, '169^' + node2.lemma);
    
    # L -> 54 -> 54 -> L
    feats[267] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^54') * tryProb_catchZero(prob_cng2cng_linear, '54^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> 54 -> 114 -> L
    feats[268] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^54') * tryProb_catchZero(prob_cng2cng_linear, '54^114') * tryProb_catchZero(prob_cng2lem_linear, '114^' + node2.lemma);
    
    # L -> 54 -> pl -> L
    feats[269] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^54') * tryProb_catchZero(prob_cng2cng_linear, '54^pl') * tryProb_catchZero(prob_cng2lem_linear, 'pl^' + node2.lemma);
    
    # L -> 54 -> 153 -> L
    feats[270] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^54') * tryProb_catchZero(prob_cng2cng_linear, '54^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # L -> 54 -> -24 -> L
    feats[271] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^54') * tryProb_catchZero(prob_cng2cng_linear, '54^-24') * tryProb_catchZero(prob_cng2lem_linear, '-24^' + node2.lemma);
    
    # L -> 54 -> 137 -> L
    feats[272] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^54') * tryProb_catchZero(prob_cng2cng_linear, '54^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> 54 -> -112 -> L
    feats[273] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^54') * tryProb_catchZero(prob_cng2cng_linear, '54^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # L -> 54 -> 101 -> L
    feats[274] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^54') * tryProb_catchZero(prob_cng2cng_linear, '54^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> 54 -> 91 -> L
    feats[275] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^54') * tryProb_catchZero(prob_cng2cng_linear, '54^91') * tryProb_catchZero(prob_cng2lem_linear, '91^' + node2.lemma);
    
    # L -> 54 -> -35 -> L
    feats[276] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^54') * tryProb_catchZero(prob_cng2cng_linear, '54^-35') * tryProb_catchZero(prob_cng2lem_linear, '-35^' + node2.lemma);
    
    # L -> 54 -> 3_du -> L
    feats[277] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^54') * tryProb_catchZero(prob_cng2cng_linear, '54^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> 114 -> dat. du. -> L
    feats[278] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^114') * tryProb_catchZero(prob_cng2cng_linear, '114^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # L -> 114 -> -151 -> L
    feats[279] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^114') * tryProb_catchZero(prob_cng2cng_linear, '114^-151') * tryProb_catchZero(prob_cng2lem_linear, '-151^' + node2.lemma);
    
    # L -> 114 -> -90 -> L
    feats[280] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^114') * tryProb_catchZero(prob_cng2cng_linear, '114^-90') * tryProb_catchZero(prob_cng2lem_linear, '-90^' + node2.lemma);
    
    # L -> 114 -> 8_fp -> L
    feats[281] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^114') * tryProb_catchZero(prob_cng2cng_linear, '114^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # L -> 114 -> 169 -> L
    feats[282] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^114') * tryProb_catchZero(prob_cng2cng_linear, '114^169') * tryProb_catchZero(prob_cng2lem_linear, '169^' + node2.lemma);
    
    # L -> 114 -> 54 -> L
    feats[283] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^114') * tryProb_catchZero(prob_cng2cng_linear, '114^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> 114 -> pl -> L
    feats[284] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^114') * tryProb_catchZero(prob_cng2cng_linear, '114^pl') * tryProb_catchZero(prob_cng2lem_linear, 'pl^' + node2.lemma);
    
    # L -> 114 -> 153 -> L
    feats[285] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^114') * tryProb_catchZero(prob_cng2cng_linear, '114^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # L -> 114 -> -24 -> L
    feats[286] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^114') * tryProb_catchZero(prob_cng2cng_linear, '114^-24') * tryProb_catchZero(prob_cng2lem_linear, '-24^' + node2.lemma);
    
    # L -> 114 -> 137 -> L
    feats[287] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^114') * tryProb_catchZero(prob_cng2cng_linear, '114^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> 114 -> -112 -> L
    feats[288] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^114') * tryProb_catchZero(prob_cng2cng_linear, '114^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # L -> 114 -> 101 -> L
    feats[289] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^114') * tryProb_catchZero(prob_cng2cng_linear, '114^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> 114 -> -35 -> L
    feats[290] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^114') * tryProb_catchZero(prob_cng2cng_linear, '114^-35') * tryProb_catchZero(prob_cng2lem_linear, '-35^' + node2.lemma);
    
    # L -> 114 -> 3_du -> L
    feats[291] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^114') * tryProb_catchZero(prob_cng2cng_linear, '114^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> loc. sg. -> 54 -> L
    feats[292] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^loc. sg.') * tryProb_catchZero(prob_cng2cng_linear, 'loc. sg.^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> -139 -> 5_sp -> L
    feats[293] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-139') * tryProb_catchZero(prob_cng2cng_linear, '-139^5_sp') * tryProb_catchZero(prob_cng2lem_linear, '5_sp^' + node2.lemma);
    
    # L -> -139 -> dat. du. -> L
    feats[294] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-139') * tryProb_catchZero(prob_cng2cng_linear, '-139^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # L -> -139 -> -151 -> L
    feats[295] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-139') * tryProb_catchZero(prob_cng2cng_linear, '-139^-151') * tryProb_catchZero(prob_cng2lem_linear, '-151^' + node2.lemma);
    
    # L -> -139 -> -90 -> L
    feats[296] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-139') * tryProb_catchZero(prob_cng2cng_linear, '-139^-90') * tryProb_catchZero(prob_cng2lem_linear, '-90^' + node2.lemma);
    
    # L -> -139 -> -57 -> L
    feats[297] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-139') * tryProb_catchZero(prob_cng2cng_linear, '-139^-57') * tryProb_catchZero(prob_cng2lem_linear, '-57^' + node2.lemma);
    
    # L -> -139 -> 8_fp -> L
    feats[298] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-139') * tryProb_catchZero(prob_cng2cng_linear, '-139^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # L -> -139 -> 169 -> L
    feats[299] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-139') * tryProb_catchZero(prob_cng2cng_linear, '-139^169') * tryProb_catchZero(prob_cng2lem_linear, '169^' + node2.lemma);
    
    # L -> -139 -> 54 -> L
    feats[300] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-139') * tryProb_catchZero(prob_cng2cng_linear, '-139^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> -139 -> pl -> L
    feats[301] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-139') * tryProb_catchZero(prob_cng2cng_linear, '-139^pl') * tryProb_catchZero(prob_cng2lem_linear, 'pl^' + node2.lemma);
    
    # L -> -139 -> 153 -> L
    feats[302] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-139') * tryProb_catchZero(prob_cng2cng_linear, '-139^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # L -> -139 -> 28 -> L
    feats[303] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-139') * tryProb_catchZero(prob_cng2cng_linear, '-139^28') * tryProb_catchZero(prob_cng2lem_linear, '28^' + node2.lemma);
    
    # L -> -139 -> -24 -> L
    feats[304] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-139') * tryProb_catchZero(prob_cng2cng_linear, '-139^-24') * tryProb_catchZero(prob_cng2lem_linear, '-24^' + node2.lemma);
    
    # L -> -139 -> 137 -> L
    feats[305] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-139') * tryProb_catchZero(prob_cng2cng_linear, '-139^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> -139 -> -112 -> L
    feats[306] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-139') * tryProb_catchZero(prob_cng2cng_linear, '-139^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # L -> -139 -> 101 -> L
    feats[307] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-139') * tryProb_catchZero(prob_cng2cng_linear, '-139^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> -139 -> 91 -> L
    feats[308] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-139') * tryProb_catchZero(prob_cng2cng_linear, '-139^91') * tryProb_catchZero(prob_cng2lem_linear, '91^' + node2.lemma);
    
    # L -> -139 -> -35 -> L
    feats[309] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-139') * tryProb_catchZero(prob_cng2cng_linear, '-139^-35') * tryProb_catchZero(prob_cng2lem_linear, '-35^' + node2.lemma);
    
    # L -> -139 -> 3_du -> L
    feats[310] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-139') * tryProb_catchZero(prob_cng2cng_linear, '-139^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> 8_du -> dat. du. -> L
    feats[311] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^8_du') * tryProb_catchZero(prob_cng2cng_linear, '8_du^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # L -> 8_du -> 8_fp -> L
    feats[312] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^8_du') * tryProb_catchZero(prob_cng2cng_linear, '8_du^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # L -> 8_du -> 54 -> L
    feats[313] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^8_du') * tryProb_catchZero(prob_cng2cng_linear, '8_du^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> 8_du -> pl -> L
    feats[314] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^8_du') * tryProb_catchZero(prob_cng2cng_linear, '8_du^pl') * tryProb_catchZero(prob_cng2lem_linear, 'pl^' + node2.lemma);
    
    # L -> 8_du -> 28 -> L
    feats[315] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^8_du') * tryProb_catchZero(prob_cng2cng_linear, '8_du^28') * tryProb_catchZero(prob_cng2lem_linear, '28^' + node2.lemma);
    
    # L -> 8_du -> 137 -> L
    feats[316] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^8_du') * tryProb_catchZero(prob_cng2cng_linear, '8_du^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> 8_du -> 101 -> L
    feats[317] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^8_du') * tryProb_catchZero(prob_cng2cng_linear, '8_du^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> 8_du -> -35 -> L
    feats[318] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^8_du') * tryProb_catchZero(prob_cng2cng_linear, '8_du^-35') * tryProb_catchZero(prob_cng2lem_linear, '-35^' + node2.lemma);
    
    # L -> 8_du -> 3_du -> L
    feats[319] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^8_du') * tryProb_catchZero(prob_cng2cng_linear, '8_du^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> -166 -> dat. du. -> L
    feats[320] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-166') * tryProb_catchZero(prob_cng2cng_linear, '-166^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # L -> -166 -> 8_fp -> L
    feats[321] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-166') * tryProb_catchZero(prob_cng2cng_linear, '-166^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # L -> -166 -> 54 -> L
    feats[322] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-166') * tryProb_catchZero(prob_cng2cng_linear, '-166^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> -166 -> pl -> L
    feats[323] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-166') * tryProb_catchZero(prob_cng2cng_linear, '-166^pl') * tryProb_catchZero(prob_cng2lem_linear, 'pl^' + node2.lemma);
    
    # L -> -166 -> 153 -> L
    feats[324] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-166') * tryProb_catchZero(prob_cng2cng_linear, '-166^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # L -> -166 -> 137 -> L
    feats[325] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-166') * tryProb_catchZero(prob_cng2cng_linear, '-166^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> -166 -> -112 -> L
    feats[326] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-166') * tryProb_catchZero(prob_cng2cng_linear, '-166^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # L -> -166 -> 101 -> L
    feats[327] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-166') * tryProb_catchZero(prob_cng2cng_linear, '-166^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> -166 -> -35 -> L
    feats[328] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-166') * tryProb_catchZero(prob_cng2cng_linear, '-166^-35') * tryProb_catchZero(prob_cng2lem_linear, '-35^' + node2.lemma);
    
    # L -> -166 -> 3_du -> L
    feats[329] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-166') * tryProb_catchZero(prob_cng2cng_linear, '-166^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> 139 -> 5_sp -> L
    feats[330] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^139') * tryProb_catchZero(prob_cng2cng_linear, '139^5_sp') * tryProb_catchZero(prob_cng2lem_linear, '5_sp^' + node2.lemma);
    
    # L -> 139 -> dat. du. -> L
    feats[331] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^139') * tryProb_catchZero(prob_cng2cng_linear, '139^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # L -> 139 -> -90 -> L
    feats[332] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^139') * tryProb_catchZero(prob_cng2cng_linear, '139^-90') * tryProb_catchZero(prob_cng2lem_linear, '-90^' + node2.lemma);
    
    # L -> 139 -> -57 -> L
    feats[333] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^139') * tryProb_catchZero(prob_cng2cng_linear, '139^-57') * tryProb_catchZero(prob_cng2lem_linear, '-57^' + node2.lemma);
    
    # L -> 139 -> 8_fp -> L
    feats[334] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^139') * tryProb_catchZero(prob_cng2cng_linear, '139^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # L -> 139 -> 169 -> L
    feats[335] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^139') * tryProb_catchZero(prob_cng2cng_linear, '139^169') * tryProb_catchZero(prob_cng2lem_linear, '169^' + node2.lemma);
    
    # L -> 139 -> 54 -> L
    feats[336] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^139') * tryProb_catchZero(prob_cng2cng_linear, '139^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> 139 -> pl -> L
    feats[337] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^139') * tryProb_catchZero(prob_cng2cng_linear, '139^pl') * tryProb_catchZero(prob_cng2lem_linear, 'pl^' + node2.lemma);
    
    # L -> 139 -> 153 -> L
    feats[338] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^139') * tryProb_catchZero(prob_cng2cng_linear, '139^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # L -> 139 -> 135 -> L
    feats[339] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^139') * tryProb_catchZero(prob_cng2cng_linear, '139^135') * tryProb_catchZero(prob_cng2lem_linear, '135^' + node2.lemma);
    
    # L -> 139 -> 28 -> L
    feats[340] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^139') * tryProb_catchZero(prob_cng2cng_linear, '139^28') * tryProb_catchZero(prob_cng2lem_linear, '28^' + node2.lemma);
    
    # L -> 139 -> -24 -> L
    feats[341] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^139') * tryProb_catchZero(prob_cng2cng_linear, '139^-24') * tryProb_catchZero(prob_cng2lem_linear, '-24^' + node2.lemma);
    
    # L -> 139 -> 137 -> L
    feats[342] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^139') * tryProb_catchZero(prob_cng2cng_linear, '139^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> 139 -> -112 -> L
    feats[343] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^139') * tryProb_catchZero(prob_cng2cng_linear, '139^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # L -> 139 -> 101 -> L
    feats[344] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^139') * tryProb_catchZero(prob_cng2cng_linear, '139^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> 139 -> 91 -> L
    feats[345] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^139') * tryProb_catchZero(prob_cng2cng_linear, '139^91') * tryProb_catchZero(prob_cng2lem_linear, '91^' + node2.lemma);
    
    # L -> 139 -> -35 -> L
    feats[346] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^139') * tryProb_catchZero(prob_cng2cng_linear, '139^-35') * tryProb_catchZero(prob_cng2lem_linear, '-35^' + node2.lemma);
    
    # L -> 139 -> 3_du -> L
    feats[347] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^139') * tryProb_catchZero(prob_cng2cng_linear, '139^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> 92 -> 101 -> L
    feats[348] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^92') * tryProb_catchZero(prob_cng2cng_linear, '92^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> pl -> dat. du. -> L
    feats[349] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^pl') * tryProb_catchZero(prob_cng2cng_linear, 'pl^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # L -> pl -> 54 -> L
    feats[350] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^pl') * tryProb_catchZero(prob_cng2cng_linear, 'pl^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> pl -> 101 -> L
    feats[351] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^pl') * tryProb_catchZero(prob_cng2cng_linear, 'pl^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> pl -> 3_du -> L
    feats[352] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^pl') * tryProb_catchZero(prob_cng2cng_linear, 'pl^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> -91 -> dat. du. -> L
    feats[353] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-91') * tryProb_catchZero(prob_cng2cng_linear, '-91^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # L -> -91 -> -151 -> L
    feats[354] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-91') * tryProb_catchZero(prob_cng2cng_linear, '-91^-151') * tryProb_catchZero(prob_cng2lem_linear, '-151^' + node2.lemma);
    
    # L -> -91 -> -90 -> L
    feats[355] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-91') * tryProb_catchZero(prob_cng2cng_linear, '-91^-90') * tryProb_catchZero(prob_cng2lem_linear, '-90^' + node2.lemma);
    
    # L -> -91 -> 8_fp -> L
    feats[356] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-91') * tryProb_catchZero(prob_cng2cng_linear, '-91^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # L -> -91 -> 54 -> L
    feats[357] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-91') * tryProb_catchZero(prob_cng2cng_linear, '-91^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> -91 -> pl -> L
    feats[358] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-91') * tryProb_catchZero(prob_cng2cng_linear, '-91^pl') * tryProb_catchZero(prob_cng2lem_linear, 'pl^' + node2.lemma);
    
    # L -> -91 -> 153 -> L
    feats[359] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-91') * tryProb_catchZero(prob_cng2cng_linear, '-91^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # L -> -91 -> 28 -> L
    feats[360] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-91') * tryProb_catchZero(prob_cng2cng_linear, '-91^28') * tryProb_catchZero(prob_cng2lem_linear, '28^' + node2.lemma);
    
    # L -> -91 -> -24 -> L
    feats[361] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-91') * tryProb_catchZero(prob_cng2cng_linear, '-91^-24') * tryProb_catchZero(prob_cng2lem_linear, '-24^' + node2.lemma);
    
    # L -> -91 -> 137 -> L
    feats[362] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-91') * tryProb_catchZero(prob_cng2cng_linear, '-91^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> -91 -> 101 -> L
    feats[363] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-91') * tryProb_catchZero(prob_cng2cng_linear, '-91^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> -91 -> 91 -> L
    feats[364] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-91') * tryProb_catchZero(prob_cng2cng_linear, '-91^91') * tryProb_catchZero(prob_cng2lem_linear, '91^' + node2.lemma);
    
    # L -> -91 -> -35 -> L
    feats[365] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-91') * tryProb_catchZero(prob_cng2cng_linear, '-91^-35') * tryProb_catchZero(prob_cng2lem_linear, '-35^' + node2.lemma);
    
    # L -> -91 -> 3_du -> L
    feats[366] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-91') * tryProb_catchZero(prob_cng2cng_linear, '-91^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> 28_sg -> 5_sp -> L
    feats[367] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^28_sg') * tryProb_catchZero(prob_cng2cng_linear, '28_sg^5_sp') * tryProb_catchZero(prob_cng2lem_linear, '5_sp^' + node2.lemma);
    
    # L -> 28_sg -> dat. du. -> L
    feats[368] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^28_sg') * tryProb_catchZero(prob_cng2cng_linear, '28_sg^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # L -> 28_sg -> 8_fp -> L
    feats[369] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^28_sg') * tryProb_catchZero(prob_cng2cng_linear, '28_sg^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # L -> 28_sg -> 54 -> L
    feats[370] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^28_sg') * tryProb_catchZero(prob_cng2cng_linear, '28_sg^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> 28_sg -> pl -> L
    feats[371] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^28_sg') * tryProb_catchZero(prob_cng2cng_linear, '28_sg^pl') * tryProb_catchZero(prob_cng2lem_linear, 'pl^' + node2.lemma);
    
    # L -> 28_sg -> 153 -> L
    feats[372] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^28_sg') * tryProb_catchZero(prob_cng2cng_linear, '28_sg^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # L -> 28_sg -> 28 -> L
    feats[373] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^28_sg') * tryProb_catchZero(prob_cng2cng_linear, '28_sg^28') * tryProb_catchZero(prob_cng2lem_linear, '28^' + node2.lemma);
    
    # L -> 28_sg -> -24 -> L
    feats[374] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^28_sg') * tryProb_catchZero(prob_cng2cng_linear, '28_sg^-24') * tryProb_catchZero(prob_cng2lem_linear, '-24^' + node2.lemma);
    
    # L -> 28_sg -> 137 -> L
    feats[375] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^28_sg') * tryProb_catchZero(prob_cng2cng_linear, '28_sg^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> 28_sg -> -112 -> L
    feats[376] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^28_sg') * tryProb_catchZero(prob_cng2cng_linear, '28_sg^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # L -> 28_sg -> 101 -> L
    feats[377] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^28_sg') * tryProb_catchZero(prob_cng2cng_linear, '28_sg^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> 28_sg -> 91 -> L
    feats[378] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^28_sg') * tryProb_catchZero(prob_cng2cng_linear, '28_sg^91') * tryProb_catchZero(prob_cng2lem_linear, '91^' + node2.lemma);
    
    # L -> 28_sg -> -35 -> L
    feats[379] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^28_sg') * tryProb_catchZero(prob_cng2cng_linear, '28_sg^-35') * tryProb_catchZero(prob_cng2lem_linear, '-35^' + node2.lemma);
    
    # L -> 28_sg -> 3_du -> L
    feats[380] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^28_sg') * tryProb_catchZero(prob_cng2cng_linear, '28_sg^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> 148 -> 5_sp -> L
    feats[381] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^148') * tryProb_catchZero(prob_cng2cng_linear, '148^5_sp') * tryProb_catchZero(prob_cng2lem_linear, '5_sp^' + node2.lemma);
    
    # L -> 148 -> dat. du. -> L
    feats[382] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^148') * tryProb_catchZero(prob_cng2cng_linear, '148^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # L -> 148 -> 8_fp -> L
    feats[383] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^148') * tryProb_catchZero(prob_cng2cng_linear, '148^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # L -> 148 -> 54 -> L
    feats[384] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^148') * tryProb_catchZero(prob_cng2cng_linear, '148^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> 148 -> pl -> L
    feats[385] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^148') * tryProb_catchZero(prob_cng2cng_linear, '148^pl') * tryProb_catchZero(prob_cng2lem_linear, 'pl^' + node2.lemma);
    
    # L -> 148 -> 153 -> L
    feats[386] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^148') * tryProb_catchZero(prob_cng2cng_linear, '148^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # L -> 148 -> 28 -> L
    feats[387] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^148') * tryProb_catchZero(prob_cng2cng_linear, '148^28') * tryProb_catchZero(prob_cng2lem_linear, '28^' + node2.lemma);
    
    # L -> 148 -> -24 -> L
    feats[388] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^148') * tryProb_catchZero(prob_cng2cng_linear, '148^-24') * tryProb_catchZero(prob_cng2lem_linear, '-24^' + node2.lemma);
    
    # L -> 148 -> 137 -> L
    feats[389] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^148') * tryProb_catchZero(prob_cng2cng_linear, '148^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> 148 -> -112 -> L
    feats[390] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^148') * tryProb_catchZero(prob_cng2cng_linear, '148^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # L -> 148 -> 101 -> L
    feats[391] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^148') * tryProb_catchZero(prob_cng2cng_linear, '148^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> 148 -> 91 -> L
    feats[392] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^148') * tryProb_catchZero(prob_cng2cng_linear, '148^91') * tryProb_catchZero(prob_cng2lem_linear, '91^' + node2.lemma);
    
    # L -> 148 -> -35 -> L
    feats[393] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^148') * tryProb_catchZero(prob_cng2cng_linear, '148^-35') * tryProb_catchZero(prob_cng2lem_linear, '-35^' + node2.lemma);
    
    # L -> 148 -> 3_du -> L
    feats[394] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^148') * tryProb_catchZero(prob_cng2cng_linear, '148^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> 5_sg -> -90 -> L
    feats[395] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^5_sg') * tryProb_catchZero(prob_cng2cng_linear, '5_sg^-90') * tryProb_catchZero(prob_cng2lem_linear, '-90^' + node2.lemma);
    
    # L -> 5_sg -> -57 -> L
    feats[396] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^5_sg') * tryProb_catchZero(prob_cng2cng_linear, '5_sg^-57') * tryProb_catchZero(prob_cng2lem_linear, '-57^' + node2.lemma);
    
    # L -> 5_sg -> 54 -> L
    feats[397] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^5_sg') * tryProb_catchZero(prob_cng2cng_linear, '5_sg^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> 5_sg -> -112 -> L
    feats[398] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^5_sg') * tryProb_catchZero(prob_cng2cng_linear, '5_sg^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # L -> 5_sg -> 101 -> L
    feats[399] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^5_sg') * tryProb_catchZero(prob_cng2cng_linear, '5_sg^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> 5_sg -> 91 -> L
    feats[400] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^5_sg') * tryProb_catchZero(prob_cng2cng_linear, '5_sg^91') * tryProb_catchZero(prob_cng2lem_linear, '91^' + node2.lemma);
    
    # L -> 153 -> 5_sp -> L
    feats[401] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^153') * tryProb_catchZero(prob_cng2cng_linear, '153^5_sp') * tryProb_catchZero(prob_cng2lem_linear, '5_sp^' + node2.lemma);
    
    # L -> 153 -> dat. du. -> L
    feats[402] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^153') * tryProb_catchZero(prob_cng2cng_linear, '153^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # L -> 153 -> -151 -> L
    feats[403] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^153') * tryProb_catchZero(prob_cng2cng_linear, '153^-151') * tryProb_catchZero(prob_cng2lem_linear, '-151^' + node2.lemma);
    
    # L -> 153 -> -90 -> L
    feats[404] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^153') * tryProb_catchZero(prob_cng2cng_linear, '153^-90') * tryProb_catchZero(prob_cng2lem_linear, '-90^' + node2.lemma);
    
    # L -> 153 -> -57 -> L
    feats[405] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^153') * tryProb_catchZero(prob_cng2cng_linear, '153^-57') * tryProb_catchZero(prob_cng2lem_linear, '-57^' + node2.lemma);
    
    # L -> 153 -> 8_fp -> L
    feats[406] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^153') * tryProb_catchZero(prob_cng2cng_linear, '153^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # L -> 153 -> 169 -> L
    feats[407] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^153') * tryProb_catchZero(prob_cng2cng_linear, '153^169') * tryProb_catchZero(prob_cng2lem_linear, '169^' + node2.lemma);
    
    # L -> 153 -> 54 -> L
    feats[408] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^153') * tryProb_catchZero(prob_cng2cng_linear, '153^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> 153 -> pl -> L
    feats[409] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^153') * tryProb_catchZero(prob_cng2cng_linear, '153^pl') * tryProb_catchZero(prob_cng2lem_linear, 'pl^' + node2.lemma);
    
    # L -> 153 -> 153 -> L
    feats[410] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^153') * tryProb_catchZero(prob_cng2cng_linear, '153^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # L -> 153 -> 135 -> L
    feats[411] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^153') * tryProb_catchZero(prob_cng2cng_linear, '153^135') * tryProb_catchZero(prob_cng2lem_linear, '135^' + node2.lemma);
    
    # L -> 153 -> 28 -> L
    feats[412] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^153') * tryProb_catchZero(prob_cng2cng_linear, '153^28') * tryProb_catchZero(prob_cng2lem_linear, '28^' + node2.lemma);
    
    # L -> 153 -> -24 -> L
    feats[413] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^153') * tryProb_catchZero(prob_cng2cng_linear, '153^-24') * tryProb_catchZero(prob_cng2lem_linear, '-24^' + node2.lemma);
    
    # L -> 153 -> 137 -> L
    feats[414] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^153') * tryProb_catchZero(prob_cng2cng_linear, '153^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> 153 -> -112 -> L
    feats[415] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^153') * tryProb_catchZero(prob_cng2cng_linear, '153^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # L -> 153 -> 89 -> L
    feats[416] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^153') * tryProb_catchZero(prob_cng2cng_linear, '153^89') * tryProb_catchZero(prob_cng2lem_linear, '89^' + node2.lemma);
    
    # L -> 153 -> 101 -> L
    feats[417] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^153') * tryProb_catchZero(prob_cng2cng_linear, '153^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> 153 -> 91 -> L
    feats[418] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^153') * tryProb_catchZero(prob_cng2cng_linear, '153^91') * tryProb_catchZero(prob_cng2lem_linear, '91^' + node2.lemma);
    
    # L -> 153 -> -35 -> L
    feats[419] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^153') * tryProb_catchZero(prob_cng2cng_linear, '153^-35') * tryProb_catchZero(prob_cng2lem_linear, '-35^' + node2.lemma);
    
    # L -> 153 -> 3_du -> L
    feats[420] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^153') * tryProb_catchZero(prob_cng2cng_linear, '153^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> -30 -> 101 -> L
    feats[421] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-30') * tryProb_catchZero(prob_cng2cng_linear, '-30^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> 11_sp -> dat. du. -> L
    feats[422] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^11_sp') * tryProb_catchZero(prob_cng2cng_linear, '11_sp^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # L -> 11_sp -> 8_fp -> L
    feats[423] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^11_sp') * tryProb_catchZero(prob_cng2cng_linear, '11_sp^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # L -> 11_sp -> 54 -> L
    feats[424] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^11_sp') * tryProb_catchZero(prob_cng2cng_linear, '11_sp^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> 11_sp -> pl -> L
    feats[425] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^11_sp') * tryProb_catchZero(prob_cng2cng_linear, '11_sp^pl') * tryProb_catchZero(prob_cng2lem_linear, 'pl^' + node2.lemma);
    
    # L -> 11_sp -> 153 -> L
    feats[426] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^11_sp') * tryProb_catchZero(prob_cng2cng_linear, '11_sp^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # L -> 11_sp -> 28 -> L
    feats[427] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^11_sp') * tryProb_catchZero(prob_cng2cng_linear, '11_sp^28') * tryProb_catchZero(prob_cng2lem_linear, '28^' + node2.lemma);
    
    # L -> 11_sp -> 137 -> L
    feats[428] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^11_sp') * tryProb_catchZero(prob_cng2cng_linear, '11_sp^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> 11_sp -> -112 -> L
    feats[429] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^11_sp') * tryProb_catchZero(prob_cng2cng_linear, '11_sp^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # L -> 11_sp -> 101 -> L
    feats[430] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^11_sp') * tryProb_catchZero(prob_cng2cng_linear, '11_sp^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> 11_sp -> -35 -> L
    feats[431] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^11_sp') * tryProb_catchZero(prob_cng2cng_linear, '11_sp^-35') * tryProb_catchZero(prob_cng2lem_linear, '-35^' + node2.lemma);
    
    # L -> 11_sp -> 3_du -> L
    feats[432] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^11_sp') * tryProb_catchZero(prob_cng2cng_linear, '11_sp^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> 30_sp -> 8_fp -> L
    feats[433] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^30_sp') * tryProb_catchZero(prob_cng2cng_linear, '30_sp^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # L -> 30_sp -> 54 -> L
    feats[434] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^30_sp') * tryProb_catchZero(prob_cng2cng_linear, '30_sp^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> 30_sp -> -24 -> L
    feats[435] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^30_sp') * tryProb_catchZero(prob_cng2cng_linear, '30_sp^-24') * tryProb_catchZero(prob_cng2lem_linear, '-24^' + node2.lemma);
    
    # L -> 30_sp -> -112 -> L
    feats[436] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^30_sp') * tryProb_catchZero(prob_cng2cng_linear, '30_sp^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # L -> 30_sp -> 101 -> L
    feats[437] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^30_sp') * tryProb_catchZero(prob_cng2cng_linear, '30_sp^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> 30_sp -> 3_du -> L
    feats[438] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^30_sp') * tryProb_catchZero(prob_cng2cng_linear, '30_sp^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> 135 -> dat. du. -> L
    feats[439] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^135') * tryProb_catchZero(prob_cng2cng_linear, '135^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # L -> 135 -> -151 -> L
    feats[440] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^135') * tryProb_catchZero(prob_cng2cng_linear, '135^-151') * tryProb_catchZero(prob_cng2lem_linear, '-151^' + node2.lemma);
    
    # L -> 135 -> -90 -> L
    feats[441] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^135') * tryProb_catchZero(prob_cng2cng_linear, '135^-90') * tryProb_catchZero(prob_cng2lem_linear, '-90^' + node2.lemma);
    
    # L -> 135 -> 54 -> L
    feats[442] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^135') * tryProb_catchZero(prob_cng2cng_linear, '135^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> 135 -> -24 -> L
    feats[443] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^135') * tryProb_catchZero(prob_cng2cng_linear, '135^-24') * tryProb_catchZero(prob_cng2lem_linear, '-24^' + node2.lemma);
    
    # L -> 135 -> 137 -> L
    feats[444] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^135') * tryProb_catchZero(prob_cng2cng_linear, '135^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> 135 -> -112 -> L
    feats[445] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^135') * tryProb_catchZero(prob_cng2cng_linear, '135^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # L -> 135 -> 101 -> L
    feats[446] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^135') * tryProb_catchZero(prob_cng2cng_linear, '135^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> -86 -> dat. du. -> L
    feats[447] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-86') * tryProb_catchZero(prob_cng2cng_linear, '-86^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # L -> -86 -> -90 -> L
    feats[448] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-86') * tryProb_catchZero(prob_cng2cng_linear, '-86^-90') * tryProb_catchZero(prob_cng2lem_linear, '-90^' + node2.lemma);
    
    # L -> -86 -> 8_fp -> L
    feats[449] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-86') * tryProb_catchZero(prob_cng2cng_linear, '-86^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # L -> -86 -> 54 -> L
    feats[450] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-86') * tryProb_catchZero(prob_cng2cng_linear, '-86^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> -86 -> 153 -> L
    feats[451] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-86') * tryProb_catchZero(prob_cng2cng_linear, '-86^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # L -> -86 -> -24 -> L
    feats[452] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-86') * tryProb_catchZero(prob_cng2cng_linear, '-86^-24') * tryProb_catchZero(prob_cng2lem_linear, '-24^' + node2.lemma);
    
    # L -> -86 -> 137 -> L
    feats[453] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-86') * tryProb_catchZero(prob_cng2cng_linear, '-86^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> -86 -> -112 -> L
    feats[454] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-86') * tryProb_catchZero(prob_cng2cng_linear, '-86^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # L -> -86 -> 101 -> L
    feats[455] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-86') * tryProb_catchZero(prob_cng2cng_linear, '-86^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> -86 -> 91 -> L
    feats[456] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-86') * tryProb_catchZero(prob_cng2cng_linear, '-86^91') * tryProb_catchZero(prob_cng2lem_linear, '91^' + node2.lemma);
    
    # L -> -86 -> -35 -> L
    feats[457] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-86') * tryProb_catchZero(prob_cng2cng_linear, '-86^-35') * tryProb_catchZero(prob_cng2lem_linear, '-35^' + node2.lemma);
    
    # L -> -86 -> 3_du -> L
    feats[458] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-86') * tryProb_catchZero(prob_cng2cng_linear, '-86^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> 7_fp -> 5_sp -> L
    feats[459] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^7_fp') * tryProb_catchZero(prob_cng2cng_linear, '7_fp^5_sp') * tryProb_catchZero(prob_cng2lem_linear, '5_sp^' + node2.lemma);
    
    # L -> 7_fp -> dat. du. -> L
    feats[460] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^7_fp') * tryProb_catchZero(prob_cng2cng_linear, '7_fp^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # L -> 7_fp -> -57 -> L
    feats[461] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^7_fp') * tryProb_catchZero(prob_cng2cng_linear, '7_fp^-57') * tryProb_catchZero(prob_cng2lem_linear, '-57^' + node2.lemma);
    
    # L -> 7_fp -> 8_fp -> L
    feats[462] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^7_fp') * tryProb_catchZero(prob_cng2cng_linear, '7_fp^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # L -> 7_fp -> 54 -> L
    feats[463] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^7_fp') * tryProb_catchZero(prob_cng2cng_linear, '7_fp^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> 7_fp -> pl -> L
    feats[464] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^7_fp') * tryProb_catchZero(prob_cng2cng_linear, '7_fp^pl') * tryProb_catchZero(prob_cng2lem_linear, 'pl^' + node2.lemma);
    
    # L -> 7_fp -> 153 -> L
    feats[465] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^7_fp') * tryProb_catchZero(prob_cng2cng_linear, '7_fp^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # L -> 7_fp -> -24 -> L
    feats[466] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^7_fp') * tryProb_catchZero(prob_cng2cng_linear, '7_fp^-24') * tryProb_catchZero(prob_cng2lem_linear, '-24^' + node2.lemma);
    
    # L -> 7_fp -> 137 -> L
    feats[467] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^7_fp') * tryProb_catchZero(prob_cng2cng_linear, '7_fp^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> 7_fp -> -112 -> L
    feats[468] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^7_fp') * tryProb_catchZero(prob_cng2cng_linear, '7_fp^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # L -> 7_fp -> 101 -> L
    feats[469] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^7_fp') * tryProb_catchZero(prob_cng2cng_linear, '7_fp^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> 7_fp -> 91 -> L
    feats[470] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^7_fp') * tryProb_catchZero(prob_cng2cng_linear, '7_fp^91') * tryProb_catchZero(prob_cng2lem_linear, '91^' + node2.lemma);
    
    # L -> 7_fp -> 3_du -> L
    feats[471] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^7_fp') * tryProb_catchZero(prob_cng2cng_linear, '7_fp^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> 111 -> dat. du. -> L
    feats[472] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^111') * tryProb_catchZero(prob_cng2cng_linear, '111^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # L -> 111 -> 54 -> L
    feats[473] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^111') * tryProb_catchZero(prob_cng2cng_linear, '111^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> 111 -> -24 -> L
    feats[474] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^111') * tryProb_catchZero(prob_cng2cng_linear, '111^-24') * tryProb_catchZero(prob_cng2lem_linear, '-24^' + node2.lemma);
    
    # L -> 111 -> -112 -> L
    feats[475] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^111') * tryProb_catchZero(prob_cng2cng_linear, '111^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # L -> 111 -> 101 -> L
    feats[476] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^111') * tryProb_catchZero(prob_cng2cng_linear, '111^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> -92 -> 137 -> L
    feats[477] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-92') * tryProb_catchZero(prob_cng2cng_linear, '-92^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> -92 -> 101 -> L
    feats[478] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-92') * tryProb_catchZero(prob_cng2cng_linear, '-92^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> 177 -> dat. du. -> L
    feats[479] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^177') * tryProb_catchZero(prob_cng2cng_linear, '177^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # L -> 177 -> -90 -> L
    feats[480] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^177') * tryProb_catchZero(prob_cng2cng_linear, '177^-90') * tryProb_catchZero(prob_cng2lem_linear, '-90^' + node2.lemma);
    
    # L -> 177 -> 8_fp -> L
    feats[481] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^177') * tryProb_catchZero(prob_cng2cng_linear, '177^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # L -> 177 -> 54 -> L
    feats[482] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^177') * tryProb_catchZero(prob_cng2cng_linear, '177^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> 177 -> pl -> L
    feats[483] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^177') * tryProb_catchZero(prob_cng2cng_linear, '177^pl') * tryProb_catchZero(prob_cng2lem_linear, 'pl^' + node2.lemma);
    
    # L -> 177 -> 153 -> L
    feats[484] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^177') * tryProb_catchZero(prob_cng2cng_linear, '177^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # L -> 177 -> 28 -> L
    feats[485] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^177') * tryProb_catchZero(prob_cng2cng_linear, '177^28') * tryProb_catchZero(prob_cng2lem_linear, '28^' + node2.lemma);
    
    # L -> 177 -> -24 -> L
    feats[486] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^177') * tryProb_catchZero(prob_cng2cng_linear, '177^-24') * tryProb_catchZero(prob_cng2lem_linear, '-24^' + node2.lemma);
    
    # L -> 177 -> 137 -> L
    feats[487] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^177') * tryProb_catchZero(prob_cng2cng_linear, '177^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> 177 -> -112 -> L
    feats[488] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^177') * tryProb_catchZero(prob_cng2cng_linear, '177^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # L -> 177 -> 101 -> L
    feats[489] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^177') * tryProb_catchZero(prob_cng2cng_linear, '177^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> 177 -> 91 -> L
    feats[490] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^177') * tryProb_catchZero(prob_cng2cng_linear, '177^91') * tryProb_catchZero(prob_cng2lem_linear, '91^' + node2.lemma);
    
    # L -> 177 -> -35 -> L
    feats[491] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^177') * tryProb_catchZero(prob_cng2cng_linear, '177^-35') * tryProb_catchZero(prob_cng2lem_linear, '-35^' + node2.lemma);
    
    # L -> 177 -> 3_du -> L
    feats[492] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^177') * tryProb_catchZero(prob_cng2cng_linear, '177^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> neutr -> dat. du. -> L
    feats[493] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^neutr') * tryProb_catchZero(prob_cng2cng_linear, 'neutr^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # L -> neutr -> -151 -> L
    feats[494] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^neutr') * tryProb_catchZero(prob_cng2cng_linear, 'neutr^-151') * tryProb_catchZero(prob_cng2lem_linear, '-151^' + node2.lemma);
    
    # L -> neutr -> 54 -> L
    feats[495] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^neutr') * tryProb_catchZero(prob_cng2cng_linear, 'neutr^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> neutr -> pl -> L
    feats[496] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^neutr') * tryProb_catchZero(prob_cng2cng_linear, 'neutr^pl') * tryProb_catchZero(prob_cng2lem_linear, 'pl^' + node2.lemma);
    
    # L -> neutr -> 153 -> L
    feats[497] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^neutr') * tryProb_catchZero(prob_cng2cng_linear, 'neutr^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # L -> neutr -> -24 -> L
    feats[498] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^neutr') * tryProb_catchZero(prob_cng2cng_linear, 'neutr^-24') * tryProb_catchZero(prob_cng2lem_linear, '-24^' + node2.lemma);
    
    # L -> neutr -> 137 -> L
    feats[499] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^neutr') * tryProb_catchZero(prob_cng2cng_linear, 'neutr^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> neutr -> -112 -> L
    feats[500] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^neutr') * tryProb_catchZero(prob_cng2cng_linear, 'neutr^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # L -> neutr -> 101 -> L
    feats[501] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^neutr') * tryProb_catchZero(prob_cng2cng_linear, 'neutr^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> neutr -> 3_du -> L
    feats[502] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^neutr') * tryProb_catchZero(prob_cng2cng_linear, 'neutr^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> 59 -> dat. du. -> L
    feats[503] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^59') * tryProb_catchZero(prob_cng2cng_linear, '59^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # L -> 59 -> 137 -> L
    feats[504] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^59') * tryProb_catchZero(prob_cng2cng_linear, '59^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> 59 -> 101 -> L
    feats[505] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^59') * tryProb_catchZero(prob_cng2cng_linear, '59^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> 59 -> 3_du -> L
    feats[506] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^59') * tryProb_catchZero(prob_cng2cng_linear, '59^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> 29 -> 5_sp -> L
    feats[507] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^29') * tryProb_catchZero(prob_cng2cng_linear, '29^5_sp') * tryProb_catchZero(prob_cng2lem_linear, '5_sp^' + node2.lemma);
    
    # L -> 29 -> dat. du. -> L
    feats[508] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^29') * tryProb_catchZero(prob_cng2cng_linear, '29^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # L -> 29 -> 8_fp -> L
    feats[509] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^29') * tryProb_catchZero(prob_cng2cng_linear, '29^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # L -> 29 -> 54 -> L
    feats[510] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^29') * tryProb_catchZero(prob_cng2cng_linear, '29^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> 29 -> pl -> L
    feats[511] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^29') * tryProb_catchZero(prob_cng2cng_linear, '29^pl') * tryProb_catchZero(prob_cng2lem_linear, 'pl^' + node2.lemma);
    
    # L -> 29 -> 153 -> L
    feats[512] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^29') * tryProb_catchZero(prob_cng2cng_linear, '29^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # L -> 29 -> 28 -> L
    feats[513] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^29') * tryProb_catchZero(prob_cng2cng_linear, '29^28') * tryProb_catchZero(prob_cng2lem_linear, '28^' + node2.lemma);
    
    # L -> 29 -> 137 -> L
    feats[514] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^29') * tryProb_catchZero(prob_cng2cng_linear, '29^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> 29 -> -112 -> L
    feats[515] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^29') * tryProb_catchZero(prob_cng2cng_linear, '29^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # L -> 29 -> 101 -> L
    feats[516] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^29') * tryProb_catchZero(prob_cng2cng_linear, '29^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> 29 -> -35 -> L
    feats[517] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^29') * tryProb_catchZero(prob_cng2cng_linear, '29^-35') * tryProb_catchZero(prob_cng2lem_linear, '-35^' + node2.lemma);
    
    # L -> 29 -> 3_du -> L
    feats[518] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^29') * tryProb_catchZero(prob_cng2cng_linear, '29^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> 156 -> 137 -> L
    feats[519] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^156') * tryProb_catchZero(prob_cng2cng_linear, '156^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> 156 -> 101 -> L
    feats[520] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^156') * tryProb_catchZero(prob_cng2cng_linear, '156^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> -153 -> 5_sp -> L
    feats[521] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-153') * tryProb_catchZero(prob_cng2cng_linear, '-153^5_sp') * tryProb_catchZero(prob_cng2lem_linear, '5_sp^' + node2.lemma);
    
    # L -> -153 -> dat. du. -> L
    feats[522] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-153') * tryProb_catchZero(prob_cng2cng_linear, '-153^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # L -> -153 -> -151 -> L
    feats[523] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-153') * tryProb_catchZero(prob_cng2cng_linear, '-153^-151') * tryProb_catchZero(prob_cng2lem_linear, '-151^' + node2.lemma);
    
    # L -> -153 -> -90 -> L
    feats[524] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-153') * tryProb_catchZero(prob_cng2cng_linear, '-153^-90') * tryProb_catchZero(prob_cng2lem_linear, '-90^' + node2.lemma);
    
    # L -> -153 -> -57 -> L
    feats[525] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-153') * tryProb_catchZero(prob_cng2cng_linear, '-153^-57') * tryProb_catchZero(prob_cng2lem_linear, '-57^' + node2.lemma);
    
    # L -> -153 -> 8_fp -> L
    feats[526] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-153') * tryProb_catchZero(prob_cng2cng_linear, '-153^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # L -> -153 -> 54 -> L
    feats[527] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-153') * tryProb_catchZero(prob_cng2cng_linear, '-153^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> -153 -> pl -> L
    feats[528] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-153') * tryProb_catchZero(prob_cng2cng_linear, '-153^pl') * tryProb_catchZero(prob_cng2lem_linear, 'pl^' + node2.lemma);
    
    # L -> -153 -> 153 -> L
    feats[529] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-153') * tryProb_catchZero(prob_cng2cng_linear, '-153^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # L -> -153 -> -24 -> L
    feats[530] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-153') * tryProb_catchZero(prob_cng2cng_linear, '-153^-24') * tryProb_catchZero(prob_cng2lem_linear, '-24^' + node2.lemma);
    
    # L -> -153 -> 137 -> L
    feats[531] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-153') * tryProb_catchZero(prob_cng2cng_linear, '-153^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> -153 -> -112 -> L
    feats[532] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-153') * tryProb_catchZero(prob_cng2cng_linear, '-153^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # L -> -153 -> 101 -> L
    feats[533] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-153') * tryProb_catchZero(prob_cng2cng_linear, '-153^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> -153 -> 91 -> L
    feats[534] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-153') * tryProb_catchZero(prob_cng2cng_linear, '-153^91') * tryProb_catchZero(prob_cng2lem_linear, '91^' + node2.lemma);
    
    # L -> -153 -> 3_du -> L
    feats[535] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-153') * tryProb_catchZero(prob_cng2cng_linear, '-153^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> 8_tp -> 5_sp -> L
    feats[536] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^8_tp') * tryProb_catchZero(prob_cng2cng_linear, '8_tp^5_sp') * tryProb_catchZero(prob_cng2lem_linear, '5_sp^' + node2.lemma);
    
    # L -> 8_tp -> dat. du. -> L
    feats[537] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^8_tp') * tryProb_catchZero(prob_cng2cng_linear, '8_tp^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # L -> 8_tp -> 8_fp -> L
    feats[538] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^8_tp') * tryProb_catchZero(prob_cng2cng_linear, '8_tp^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # L -> 8_tp -> 54 -> L
    feats[539] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^8_tp') * tryProb_catchZero(prob_cng2cng_linear, '8_tp^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> 8_tp -> pl -> L
    feats[540] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^8_tp') * tryProb_catchZero(prob_cng2cng_linear, '8_tp^pl') * tryProb_catchZero(prob_cng2lem_linear, 'pl^' + node2.lemma);
    
    # L -> 8_tp -> 153 -> L
    feats[541] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^8_tp') * tryProb_catchZero(prob_cng2cng_linear, '8_tp^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # L -> 8_tp -> 28 -> L
    feats[542] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^8_tp') * tryProb_catchZero(prob_cng2cng_linear, '8_tp^28') * tryProb_catchZero(prob_cng2lem_linear, '28^' + node2.lemma);
    
    # L -> 8_tp -> 137 -> L
    feats[543] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^8_tp') * tryProb_catchZero(prob_cng2cng_linear, '8_tp^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> 8_tp -> -112 -> L
    feats[544] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^8_tp') * tryProb_catchZero(prob_cng2cng_linear, '8_tp^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # L -> 8_tp -> 101 -> L
    feats[545] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^8_tp') * tryProb_catchZero(prob_cng2cng_linear, '8_tp^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> 8_tp -> -35 -> L
    feats[546] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^8_tp') * tryProb_catchZero(prob_cng2cng_linear, '8_tp^-35') * tryProb_catchZero(prob_cng2lem_linear, '-35^' + node2.lemma);
    
    # L -> 8_tp -> 3_du -> L
    feats[547] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^8_tp') * tryProb_catchZero(prob_cng2cng_linear, '8_tp^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> 4_tp -> 101 -> L
    feats[548] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^4_tp') * tryProb_catchZero(prob_cng2cng_linear, '4_tp^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> 28 -> dat. du. -> L
    feats[549] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^28') * tryProb_catchZero(prob_cng2cng_linear, '28^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # L -> 28 -> 54 -> L
    feats[550] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^28') * tryProb_catchZero(prob_cng2cng_linear, '28^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> 28 -> 137 -> L
    feats[551] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^28') * tryProb_catchZero(prob_cng2cng_linear, '28^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> 28 -> -112 -> L
    feats[552] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^28') * tryProb_catchZero(prob_cng2cng_linear, '28^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # L -> 28 -> 101 -> L
    feats[553] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^28') * tryProb_catchZero(prob_cng2cng_linear, '28^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> 28 -> 3_du -> L
    feats[554] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^28') * tryProb_catchZero(prob_cng2cng_linear, '28^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> -24 -> 5_sp -> L
    feats[555] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-24') * tryProb_catchZero(prob_cng2cng_linear, '-24^5_sp') * tryProb_catchZero(prob_cng2lem_linear, '5_sp^' + node2.lemma);
    
    # L -> -24 -> dat. du. -> L
    feats[556] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-24') * tryProb_catchZero(prob_cng2cng_linear, '-24^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # L -> -24 -> -151 -> L
    feats[557] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-24') * tryProb_catchZero(prob_cng2cng_linear, '-24^-151') * tryProb_catchZero(prob_cng2lem_linear, '-151^' + node2.lemma);
    
    # L -> -24 -> -90 -> L
    feats[558] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-24') * tryProb_catchZero(prob_cng2cng_linear, '-24^-90') * tryProb_catchZero(prob_cng2lem_linear, '-90^' + node2.lemma);
    
    # L -> -24 -> -57 -> L
    feats[559] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-24') * tryProb_catchZero(prob_cng2cng_linear, '-24^-57') * tryProb_catchZero(prob_cng2lem_linear, '-57^' + node2.lemma);
    
    # L -> -24 -> 8_fp -> L
    feats[560] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-24') * tryProb_catchZero(prob_cng2cng_linear, '-24^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # L -> -24 -> 54 -> L
    feats[561] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-24') * tryProb_catchZero(prob_cng2cng_linear, '-24^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> -24 -> pl -> L
    feats[562] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-24') * tryProb_catchZero(prob_cng2cng_linear, '-24^pl') * tryProb_catchZero(prob_cng2lem_linear, 'pl^' + node2.lemma);
    
    # L -> -24 -> 153 -> L
    feats[563] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-24') * tryProb_catchZero(prob_cng2cng_linear, '-24^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # L -> -24 -> 28 -> L
    feats[564] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-24') * tryProb_catchZero(prob_cng2cng_linear, '-24^28') * tryProb_catchZero(prob_cng2lem_linear, '28^' + node2.lemma);
    
    # L -> -24 -> -24 -> L
    feats[565] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-24') * tryProb_catchZero(prob_cng2cng_linear, '-24^-24') * tryProb_catchZero(prob_cng2lem_linear, '-24^' + node2.lemma);
    
    # L -> -24 -> 137 -> L
    feats[566] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-24') * tryProb_catchZero(prob_cng2cng_linear, '-24^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> -24 -> -112 -> L
    feats[567] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-24') * tryProb_catchZero(prob_cng2cng_linear, '-24^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # L -> -24 -> 101 -> L
    feats[568] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-24') * tryProb_catchZero(prob_cng2cng_linear, '-24^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> -24 -> 91 -> L
    feats[569] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-24') * tryProb_catchZero(prob_cng2cng_linear, '-24^91') * tryProb_catchZero(prob_cng2lem_linear, '91^' + node2.lemma);
    
    # L -> -24 -> -35 -> L
    feats[570] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-24') * tryProb_catchZero(prob_cng2cng_linear, '-24^-35') * tryProb_catchZero(prob_cng2lem_linear, '-35^' + node2.lemma);
    
    # L -> -24 -> 3_du -> L
    feats[571] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-24') * tryProb_catchZero(prob_cng2cng_linear, '-24^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> 137 -> dat. du. -> L
    feats[572] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^137') * tryProb_catchZero(prob_cng2cng_linear, '137^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # L -> 137 -> -151 -> L
    feats[573] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^137') * tryProb_catchZero(prob_cng2cng_linear, '137^-151') * tryProb_catchZero(prob_cng2lem_linear, '-151^' + node2.lemma);
    
    # L -> 137 -> -90 -> L
    feats[574] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^137') * tryProb_catchZero(prob_cng2cng_linear, '137^-90') * tryProb_catchZero(prob_cng2lem_linear, '-90^' + node2.lemma);
    
    # L -> 137 -> 8_fp -> L
    feats[575] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^137') * tryProb_catchZero(prob_cng2cng_linear, '137^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # L -> 137 -> 54 -> L
    feats[576] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^137') * tryProb_catchZero(prob_cng2cng_linear, '137^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> 137 -> 153 -> L
    feats[577] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^137') * tryProb_catchZero(prob_cng2cng_linear, '137^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # L -> 137 -> -24 -> L
    feats[578] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^137') * tryProb_catchZero(prob_cng2cng_linear, '137^-24') * tryProb_catchZero(prob_cng2lem_linear, '-24^' + node2.lemma);
    
    # L -> 137 -> 137 -> L
    feats[579] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^137') * tryProb_catchZero(prob_cng2cng_linear, '137^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> 137 -> -112 -> L
    feats[580] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^137') * tryProb_catchZero(prob_cng2cng_linear, '137^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # L -> 137 -> 101 -> L
    feats[581] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^137') * tryProb_catchZero(prob_cng2cng_linear, '137^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> 137 -> 91 -> L
    feats[582] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^137') * tryProb_catchZero(prob_cng2cng_linear, '137^91') * tryProb_catchZero(prob_cng2lem_linear, '91^' + node2.lemma);
    
    # L -> 137 -> 3_du -> L
    feats[583] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^137') * tryProb_catchZero(prob_cng2cng_linear, '137^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> -112 -> 5_sp -> L
    feats[584] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-112') * tryProb_catchZero(prob_cng2cng_linear, '-112^5_sp') * tryProb_catchZero(prob_cng2lem_linear, '5_sp^' + node2.lemma);
    
    # L -> -112 -> dat. du. -> L
    feats[585] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-112') * tryProb_catchZero(prob_cng2cng_linear, '-112^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # L -> -112 -> -151 -> L
    feats[586] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-112') * tryProb_catchZero(prob_cng2cng_linear, '-112^-151') * tryProb_catchZero(prob_cng2lem_linear, '-151^' + node2.lemma);
    
    # L -> -112 -> -57 -> L
    feats[587] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-112') * tryProb_catchZero(prob_cng2cng_linear, '-112^-57') * tryProb_catchZero(prob_cng2lem_linear, '-57^' + node2.lemma);
    
    # L -> -112 -> 8_fp -> L
    feats[588] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-112') * tryProb_catchZero(prob_cng2cng_linear, '-112^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # L -> -112 -> 169 -> L
    feats[589] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-112') * tryProb_catchZero(prob_cng2cng_linear, '-112^169') * tryProb_catchZero(prob_cng2lem_linear, '169^' + node2.lemma);
    
    # L -> -112 -> 54 -> L
    feats[590] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-112') * tryProb_catchZero(prob_cng2cng_linear, '-112^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> -112 -> pl -> L
    feats[591] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-112') * tryProb_catchZero(prob_cng2cng_linear, '-112^pl') * tryProb_catchZero(prob_cng2lem_linear, 'pl^' + node2.lemma);
    
    # L -> -112 -> 153 -> L
    feats[592] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-112') * tryProb_catchZero(prob_cng2cng_linear, '-112^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # L -> -112 -> 28 -> L
    feats[593] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-112') * tryProb_catchZero(prob_cng2cng_linear, '-112^28') * tryProb_catchZero(prob_cng2lem_linear, '28^' + node2.lemma);
    
    # L -> -112 -> -24 -> L
    feats[594] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-112') * tryProb_catchZero(prob_cng2cng_linear, '-112^-24') * tryProb_catchZero(prob_cng2lem_linear, '-24^' + node2.lemma);
    
    # L -> -112 -> 137 -> L
    feats[595] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-112') * tryProb_catchZero(prob_cng2cng_linear, '-112^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> -112 -> -112 -> L
    feats[596] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-112') * tryProb_catchZero(prob_cng2cng_linear, '-112^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # L -> -112 -> 101 -> L
    feats[597] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-112') * tryProb_catchZero(prob_cng2cng_linear, '-112^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> -112 -> 91 -> L
    feats[598] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-112') * tryProb_catchZero(prob_cng2cng_linear, '-112^91') * tryProb_catchZero(prob_cng2lem_linear, '91^' + node2.lemma);
    
    # L -> -112 -> -35 -> L
    feats[599] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-112') * tryProb_catchZero(prob_cng2cng_linear, '-112^-35') * tryProb_catchZero(prob_cng2lem_linear, '-35^' + node2.lemma);
    
    # L -> -112 -> 3_du -> L
    feats[600] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-112') * tryProb_catchZero(prob_cng2cng_linear, '-112^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> 89 -> -24 -> L
    feats[601] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^89') * tryProb_catchZero(prob_cng2cng_linear, '89^-24') * tryProb_catchZero(prob_cng2lem_linear, '-24^' + node2.lemma);
    
    # L -> 89 -> 3_du -> L
    feats[602] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^89') * tryProb_catchZero(prob_cng2cng_linear, '89^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> -161 -> 8_fp -> L
    feats[603] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-161') * tryProb_catchZero(prob_cng2cng_linear, '-161^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # L -> -161 -> 54 -> L
    feats[604] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-161') * tryProb_catchZero(prob_cng2cng_linear, '-161^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> -161 -> pl -> L
    feats[605] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-161') * tryProb_catchZero(prob_cng2cng_linear, '-161^pl') * tryProb_catchZero(prob_cng2lem_linear, 'pl^' + node2.lemma);
    
    # L -> -161 -> 153 -> L
    feats[606] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-161') * tryProb_catchZero(prob_cng2cng_linear, '-161^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # L -> -161 -> 137 -> L
    feats[607] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-161') * tryProb_catchZero(prob_cng2cng_linear, '-161^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> -161 -> -112 -> L
    feats[608] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-161') * tryProb_catchZero(prob_cng2cng_linear, '-161^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # L -> -161 -> 101 -> L
    feats[609] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-161') * tryProb_catchZero(prob_cng2cng_linear, '-161^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> -161 -> -35 -> L
    feats[610] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-161') * tryProb_catchZero(prob_cng2cng_linear, '-161^-35') * tryProb_catchZero(prob_cng2lem_linear, '-35^' + node2.lemma);
    
    # L -> -161 -> 3_du -> L
    feats[611] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-161') * tryProb_catchZero(prob_cng2cng_linear, '-161^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> -96 -> dat. du. -> L
    feats[612] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-96') * tryProb_catchZero(prob_cng2cng_linear, '-96^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # L -> -96 -> 8_fp -> L
    feats[613] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-96') * tryProb_catchZero(prob_cng2cng_linear, '-96^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # L -> -96 -> 54 -> L
    feats[614] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-96') * tryProb_catchZero(prob_cng2cng_linear, '-96^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> -96 -> pl -> L
    feats[615] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-96') * tryProb_catchZero(prob_cng2cng_linear, '-96^pl') * tryProb_catchZero(prob_cng2lem_linear, 'pl^' + node2.lemma);
    
    # L -> -96 -> 153 -> L
    feats[616] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-96') * tryProb_catchZero(prob_cng2cng_linear, '-96^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # L -> -96 -> -24 -> L
    feats[617] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-96') * tryProb_catchZero(prob_cng2cng_linear, '-96^-24') * tryProb_catchZero(prob_cng2lem_linear, '-24^' + node2.lemma);
    
    # L -> -96 -> 137 -> L
    feats[618] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-96') * tryProb_catchZero(prob_cng2cng_linear, '-96^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> -96 -> -112 -> L
    feats[619] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-96') * tryProb_catchZero(prob_cng2cng_linear, '-96^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # L -> -96 -> 101 -> L
    feats[620] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-96') * tryProb_catchZero(prob_cng2cng_linear, '-96^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> -96 -> 91 -> L
    feats[621] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-96') * tryProb_catchZero(prob_cng2cng_linear, '-96^91') * tryProb_catchZero(prob_cng2lem_linear, '91^' + node2.lemma);
    
    # L -> -96 -> 3_du -> L
    feats[622] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-96') * tryProb_catchZero(prob_cng2cng_linear, '-96^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> 174 -> 8_fp -> L
    feats[623] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^174') * tryProb_catchZero(prob_cng2cng_linear, '174^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # L -> 174 -> 54 -> L
    feats[624] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^174') * tryProb_catchZero(prob_cng2cng_linear, '174^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> 174 -> 101 -> L
    feats[625] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^174') * tryProb_catchZero(prob_cng2cng_linear, '174^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> 162 -> dat. du. -> L
    feats[626] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^162') * tryProb_catchZero(prob_cng2cng_linear, '162^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # L -> 162 -> 8_fp -> L
    feats[627] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^162') * tryProb_catchZero(prob_cng2cng_linear, '162^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # L -> 162 -> 54 -> L
    feats[628] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^162') * tryProb_catchZero(prob_cng2cng_linear, '162^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> 162 -> pl -> L
    feats[629] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^162') * tryProb_catchZero(prob_cng2cng_linear, '162^pl') * tryProb_catchZero(prob_cng2lem_linear, 'pl^' + node2.lemma);
    
    # L -> 162 -> 28 -> L
    feats[630] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^162') * tryProb_catchZero(prob_cng2cng_linear, '162^28') * tryProb_catchZero(prob_cng2lem_linear, '28^' + node2.lemma);
    
    # L -> 162 -> 137 -> L
    feats[631] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^162') * tryProb_catchZero(prob_cng2cng_linear, '162^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> 162 -> 101 -> L
    feats[632] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^162') * tryProb_catchZero(prob_cng2cng_linear, '162^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> 162 -> 3_du -> L
    feats[633] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^162') * tryProb_catchZero(prob_cng2cng_linear, '162^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> gen. pl. -> 137 -> L
    feats[634] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^gen. pl.') * tryProb_catchZero(prob_cng2cng_linear, 'gen. pl.^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> gen. pl. -> 101 -> L
    feats[635] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^gen. pl.') * tryProb_catchZero(prob_cng2cng_linear, 'gen. pl.^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> pl_tp -> 5_sp -> L
    feats[636] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^pl_tp') * tryProb_catchZero(prob_cng2cng_linear, 'pl_tp^5_sp') * tryProb_catchZero(prob_cng2lem_linear, '5_sp^' + node2.lemma);
    
    # L -> pl_tp -> dat. du. -> L
    feats[637] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^pl_tp') * tryProb_catchZero(prob_cng2cng_linear, 'pl_tp^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # L -> pl_tp -> -57 -> L
    feats[638] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^pl_tp') * tryProb_catchZero(prob_cng2cng_linear, 'pl_tp^-57') * tryProb_catchZero(prob_cng2lem_linear, '-57^' + node2.lemma);
    
    # L -> pl_tp -> 8_fp -> L
    feats[639] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^pl_tp') * tryProb_catchZero(prob_cng2cng_linear, 'pl_tp^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # L -> pl_tp -> 54 -> L
    feats[640] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^pl_tp') * tryProb_catchZero(prob_cng2cng_linear, 'pl_tp^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> pl_tp -> pl -> L
    feats[641] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^pl_tp') * tryProb_catchZero(prob_cng2cng_linear, 'pl_tp^pl') * tryProb_catchZero(prob_cng2lem_linear, 'pl^' + node2.lemma);
    
    # L -> pl_tp -> 153 -> L
    feats[642] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^pl_tp') * tryProb_catchZero(prob_cng2cng_linear, 'pl_tp^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # L -> pl_tp -> -24 -> L
    feats[643] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^pl_tp') * tryProb_catchZero(prob_cng2cng_linear, 'pl_tp^-24') * tryProb_catchZero(prob_cng2lem_linear, '-24^' + node2.lemma);
    
    # L -> pl_tp -> 137 -> L
    feats[644] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^pl_tp') * tryProb_catchZero(prob_cng2cng_linear, 'pl_tp^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> pl_tp -> -112 -> L
    feats[645] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^pl_tp') * tryProb_catchZero(prob_cng2cng_linear, 'pl_tp^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # L -> pl_tp -> 101 -> L
    feats[646] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^pl_tp') * tryProb_catchZero(prob_cng2cng_linear, 'pl_tp^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> pl_tp -> 91 -> L
    feats[647] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^pl_tp') * tryProb_catchZero(prob_cng2cng_linear, 'pl_tp^91') * tryProb_catchZero(prob_cng2lem_linear, '91^' + node2.lemma);
    
    # L -> pl_tp -> -35 -> L
    feats[648] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^pl_tp') * tryProb_catchZero(prob_cng2cng_linear, 'pl_tp^-35') * tryProb_catchZero(prob_cng2lem_linear, '-35^' + node2.lemma);
    
    # L -> pl_tp -> 3_du -> L
    feats[649] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^pl_tp') * tryProb_catchZero(prob_cng2cng_linear, 'pl_tp^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> 101 -> 5_sp -> L
    feats[650] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^101') * tryProb_catchZero(prob_cng2cng_linear, '101^5_sp') * tryProb_catchZero(prob_cng2lem_linear, '5_sp^' + node2.lemma);
    
    # L -> 101 -> dat. du. -> L
    feats[651] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^101') * tryProb_catchZero(prob_cng2cng_linear, '101^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # L -> 101 -> -151 -> L
    feats[652] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^101') * tryProb_catchZero(prob_cng2cng_linear, '101^-151') * tryProb_catchZero(prob_cng2lem_linear, '-151^' + node2.lemma);
    
    # L -> 101 -> -90 -> L
    feats[653] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^101') * tryProb_catchZero(prob_cng2cng_linear, '101^-90') * tryProb_catchZero(prob_cng2lem_linear, '-90^' + node2.lemma);
    
    # L -> 101 -> -57 -> L
    feats[654] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^101') * tryProb_catchZero(prob_cng2cng_linear, '101^-57') * tryProb_catchZero(prob_cng2lem_linear, '-57^' + node2.lemma);
    
    # L -> 101 -> 8_fp -> L
    feats[655] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^101') * tryProb_catchZero(prob_cng2cng_linear, '101^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # L -> 101 -> 169 -> L
    feats[656] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^101') * tryProb_catchZero(prob_cng2cng_linear, '101^169') * tryProb_catchZero(prob_cng2lem_linear, '169^' + node2.lemma);
    
    # L -> 101 -> 54 -> L
    feats[657] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^101') * tryProb_catchZero(prob_cng2cng_linear, '101^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> 101 -> pl -> L
    feats[658] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^101') * tryProb_catchZero(prob_cng2cng_linear, '101^pl') * tryProb_catchZero(prob_cng2lem_linear, 'pl^' + node2.lemma);
    
    # L -> 101 -> 153 -> L
    feats[659] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^101') * tryProb_catchZero(prob_cng2cng_linear, '101^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # L -> 101 -> 135 -> L
    feats[660] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^101') * tryProb_catchZero(prob_cng2cng_linear, '101^135') * tryProb_catchZero(prob_cng2lem_linear, '135^' + node2.lemma);
    
    # L -> 101 -> 28 -> L
    feats[661] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^101') * tryProb_catchZero(prob_cng2cng_linear, '101^28') * tryProb_catchZero(prob_cng2lem_linear, '28^' + node2.lemma);
    
    # L -> 101 -> -24 -> L
    feats[662] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^101') * tryProb_catchZero(prob_cng2cng_linear, '101^-24') * tryProb_catchZero(prob_cng2lem_linear, '-24^' + node2.lemma);
    
    # L -> 101 -> 137 -> L
    feats[663] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^101') * tryProb_catchZero(prob_cng2cng_linear, '101^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> 101 -> -112 -> L
    feats[664] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^101') * tryProb_catchZero(prob_cng2cng_linear, '101^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # L -> 101 -> 101 -> L
    feats[665] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^101') * tryProb_catchZero(prob_cng2cng_linear, '101^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> 101 -> 91 -> L
    feats[666] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^101') * tryProb_catchZero(prob_cng2cng_linear, '101^91') * tryProb_catchZero(prob_cng2lem_linear, '91^' + node2.lemma);
    
    # L -> 101 -> -35 -> L
    feats[667] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^101') * tryProb_catchZero(prob_cng2cng_linear, '101^-35') * tryProb_catchZero(prob_cng2lem_linear, '-35^' + node2.lemma);
    
    # L -> 101 -> 3_du -> L
    feats[668] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^101') * tryProb_catchZero(prob_cng2cng_linear, '101^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> 91 -> 5_sp -> L
    feats[669] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^91') * tryProb_catchZero(prob_cng2cng_linear, '91^5_sp') * tryProb_catchZero(prob_cng2lem_linear, '5_sp^' + node2.lemma);
    
    # L -> 91 -> dat. du. -> L
    feats[670] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^91') * tryProb_catchZero(prob_cng2cng_linear, '91^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # L -> 91 -> -151 -> L
    feats[671] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^91') * tryProb_catchZero(prob_cng2cng_linear, '91^-151') * tryProb_catchZero(prob_cng2lem_linear, '-151^' + node2.lemma);
    
    # L -> 91 -> -90 -> L
    feats[672] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^91') * tryProb_catchZero(prob_cng2cng_linear, '91^-90') * tryProb_catchZero(prob_cng2lem_linear, '-90^' + node2.lemma);
    
    # L -> 91 -> -57 -> L
    feats[673] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^91') * tryProb_catchZero(prob_cng2cng_linear, '91^-57') * tryProb_catchZero(prob_cng2lem_linear, '-57^' + node2.lemma);
    
    # L -> 91 -> 8_fp -> L
    feats[674] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^91') * tryProb_catchZero(prob_cng2cng_linear, '91^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # L -> 91 -> 169 -> L
    feats[675] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^91') * tryProb_catchZero(prob_cng2cng_linear, '91^169') * tryProb_catchZero(prob_cng2lem_linear, '169^' + node2.lemma);
    
    # L -> 91 -> 54 -> L
    feats[676] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^91') * tryProb_catchZero(prob_cng2cng_linear, '91^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> 91 -> pl -> L
    feats[677] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^91') * tryProb_catchZero(prob_cng2cng_linear, '91^pl') * tryProb_catchZero(prob_cng2lem_linear, 'pl^' + node2.lemma);
    
    # L -> 91 -> 153 -> L
    feats[678] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^91') * tryProb_catchZero(prob_cng2cng_linear, '91^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # L -> 91 -> 135 -> L
    feats[679] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^91') * tryProb_catchZero(prob_cng2cng_linear, '91^135') * tryProb_catchZero(prob_cng2lem_linear, '135^' + node2.lemma);
    
    # L -> 91 -> 28 -> L
    feats[680] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^91') * tryProb_catchZero(prob_cng2cng_linear, '91^28') * tryProb_catchZero(prob_cng2lem_linear, '28^' + node2.lemma);
    
    # L -> 91 -> -24 -> L
    feats[681] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^91') * tryProb_catchZero(prob_cng2cng_linear, '91^-24') * tryProb_catchZero(prob_cng2lem_linear, '-24^' + node2.lemma);
    
    # L -> 91 -> 137 -> L
    feats[682] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^91') * tryProb_catchZero(prob_cng2cng_linear, '91^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> 91 -> -112 -> L
    feats[683] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^91') * tryProb_catchZero(prob_cng2cng_linear, '91^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # L -> 91 -> 101 -> L
    feats[684] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^91') * tryProb_catchZero(prob_cng2cng_linear, '91^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> 91 -> 91 -> L
    feats[685] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^91') * tryProb_catchZero(prob_cng2cng_linear, '91^91') * tryProb_catchZero(prob_cng2lem_linear, '91^' + node2.lemma);
    
    # L -> 91 -> -35 -> L
    feats[686] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^91') * tryProb_catchZero(prob_cng2cng_linear, '91^-35') * tryProb_catchZero(prob_cng2lem_linear, '-35^' + node2.lemma);
    
    # L -> 91 -> 3_du -> L
    feats[687] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^91') * tryProb_catchZero(prob_cng2cng_linear, '91^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> -35 -> -90 -> L
    feats[688] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-35') * tryProb_catchZero(prob_cng2cng_linear, '-35^-90') * tryProb_catchZero(prob_cng2lem_linear, '-90^' + node2.lemma);
    
    # L -> -35 -> 54 -> L
    feats[689] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-35') * tryProb_catchZero(prob_cng2cng_linear, '-35^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> -35 -> 101 -> L
    feats[690] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-35') * tryProb_catchZero(prob_cng2cng_linear, '-35^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> 11_fp -> 5_sp -> L
    feats[691] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^11_fp') * tryProb_catchZero(prob_cng2cng_linear, '11_fp^5_sp') * tryProb_catchZero(prob_cng2lem_linear, '5_sp^' + node2.lemma);
    
    # L -> 11_fp -> dat. du. -> L
    feats[692] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^11_fp') * tryProb_catchZero(prob_cng2cng_linear, '11_fp^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # L -> 11_fp -> -57 -> L
    feats[693] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^11_fp') * tryProb_catchZero(prob_cng2cng_linear, '11_fp^-57') * tryProb_catchZero(prob_cng2lem_linear, '-57^' + node2.lemma);
    
    # L -> 11_fp -> 8_fp -> L
    feats[694] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^11_fp') * tryProb_catchZero(prob_cng2cng_linear, '11_fp^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # L -> 11_fp -> 54 -> L
    feats[695] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^11_fp') * tryProb_catchZero(prob_cng2cng_linear, '11_fp^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> 11_fp -> pl -> L
    feats[696] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^11_fp') * tryProb_catchZero(prob_cng2cng_linear, '11_fp^pl') * tryProb_catchZero(prob_cng2lem_linear, 'pl^' + node2.lemma);
    
    # L -> 11_fp -> 153 -> L
    feats[697] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^11_fp') * tryProb_catchZero(prob_cng2cng_linear, '11_fp^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # L -> 11_fp -> 28 -> L
    feats[698] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^11_fp') * tryProb_catchZero(prob_cng2cng_linear, '11_fp^28') * tryProb_catchZero(prob_cng2lem_linear, '28^' + node2.lemma);
    
    # L -> 11_fp -> -24 -> L
    feats[699] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^11_fp') * tryProb_catchZero(prob_cng2cng_linear, '11_fp^-24') * tryProb_catchZero(prob_cng2lem_linear, '-24^' + node2.lemma);
    
    # L -> 11_fp -> 137 -> L
    feats[700] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^11_fp') * tryProb_catchZero(prob_cng2cng_linear, '11_fp^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> 11_fp -> -112 -> L
    feats[701] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^11_fp') * tryProb_catchZero(prob_cng2cng_linear, '11_fp^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # L -> 11_fp -> 101 -> L
    feats[702] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^11_fp') * tryProb_catchZero(prob_cng2cng_linear, '11_fp^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> 11_fp -> 91 -> L
    feats[703] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^11_fp') * tryProb_catchZero(prob_cng2cng_linear, '11_fp^91') * tryProb_catchZero(prob_cng2lem_linear, '91^' + node2.lemma);
    
    # L -> 11_fp -> -35 -> L
    feats[704] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^11_fp') * tryProb_catchZero(prob_cng2cng_linear, '11_fp^-35') * tryProb_catchZero(prob_cng2lem_linear, '-35^' + node2.lemma);
    
    # L -> 11_fp -> 3_du -> L
    feats[705] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^11_fp') * tryProb_catchZero(prob_cng2cng_linear, '11_fp^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> 3_du -> 5_sp -> L
    feats[706] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^3_du') * tryProb_catchZero(prob_cng2cng_linear, '3_du^5_sp') * tryProb_catchZero(prob_cng2lem_linear, '5_sp^' + node2.lemma);
    
    # L -> 3_du -> dat. du. -> L
    feats[707] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^3_du') * tryProb_catchZero(prob_cng2cng_linear, '3_du^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # L -> 3_du -> -151 -> L
    feats[708] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^3_du') * tryProb_catchZero(prob_cng2cng_linear, '3_du^-151') * tryProb_catchZero(prob_cng2lem_linear, '-151^' + node2.lemma);
    
    # L -> 3_du -> -90 -> L
    feats[709] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^3_du') * tryProb_catchZero(prob_cng2cng_linear, '3_du^-90') * tryProb_catchZero(prob_cng2lem_linear, '-90^' + node2.lemma);
    
    # L -> 3_du -> -57 -> L
    feats[710] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^3_du') * tryProb_catchZero(prob_cng2cng_linear, '3_du^-57') * tryProb_catchZero(prob_cng2lem_linear, '-57^' + node2.lemma);
    
    # L -> 3_du -> 8_fp -> L
    feats[711] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^3_du') * tryProb_catchZero(prob_cng2cng_linear, '3_du^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # L -> 3_du -> 169 -> L
    feats[712] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^3_du') * tryProb_catchZero(prob_cng2cng_linear, '3_du^169') * tryProb_catchZero(prob_cng2lem_linear, '169^' + node2.lemma);
    
    # L -> 3_du -> 54 -> L
    feats[713] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^3_du') * tryProb_catchZero(prob_cng2cng_linear, '3_du^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> 3_du -> 114 -> L
    feats[714] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^3_du') * tryProb_catchZero(prob_cng2cng_linear, '3_du^114') * tryProb_catchZero(prob_cng2lem_linear, '114^' + node2.lemma);
    
    # L -> 3_du -> pl -> L
    feats[715] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^3_du') * tryProb_catchZero(prob_cng2cng_linear, '3_du^pl') * tryProb_catchZero(prob_cng2lem_linear, 'pl^' + node2.lemma);
    
    # L -> 3_du -> 153 -> L
    feats[716] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^3_du') * tryProb_catchZero(prob_cng2cng_linear, '3_du^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # L -> 3_du -> 28 -> L
    feats[717] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^3_du') * tryProb_catchZero(prob_cng2cng_linear, '3_du^28') * tryProb_catchZero(prob_cng2lem_linear, '28^' + node2.lemma);
    
    # L -> 3_du -> -24 -> L
    feats[718] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^3_du') * tryProb_catchZero(prob_cng2cng_linear, '3_du^-24') * tryProb_catchZero(prob_cng2lem_linear, '-24^' + node2.lemma);
    
    # L -> 3_du -> 137 -> L
    feats[719] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^3_du') * tryProb_catchZero(prob_cng2cng_linear, '3_du^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> 3_du -> -112 -> L
    feats[720] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^3_du') * tryProb_catchZero(prob_cng2cng_linear, '3_du^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # L -> 3_du -> 89 -> L
    feats[721] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^3_du') * tryProb_catchZero(prob_cng2cng_linear, '3_du^89') * tryProb_catchZero(prob_cng2lem_linear, '89^' + node2.lemma);
    
    # L -> 3_du -> 101 -> L
    feats[722] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^3_du') * tryProb_catchZero(prob_cng2cng_linear, '3_du^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> 3_du -> 91 -> L
    feats[723] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^3_du') * tryProb_catchZero(prob_cng2cng_linear, '3_du^91') * tryProb_catchZero(prob_cng2lem_linear, '91^' + node2.lemma);
    
    # L -> 3_du -> -35 -> L
    feats[724] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^3_du') * tryProb_catchZero(prob_cng2cng_linear, '3_du^-35') * tryProb_catchZero(prob_cng2lem_linear, '-35^' + node2.lemma);
    
    # L -> 3_du -> 3_du -> L
    feats[725] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^3_du') * tryProb_catchZero(prob_cng2cng_linear, '3_du^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> 8_sp -> 137 -> L
    feats[726] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^8_sp') * tryProb_catchZero(prob_cng2cng_linear, '8_sp^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> 8_sp -> 101 -> L
    feats[727] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^8_sp') * tryProb_catchZero(prob_cng2cng_linear, '8_sp^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> -296 -> dat. du. -> L
    feats[728] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-296') * tryProb_catchZero(prob_cng2cng_linear, '-296^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # L -> -296 -> 54 -> L
    feats[729] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-296') * tryProb_catchZero(prob_cng2cng_linear, '-296^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> -296 -> 28 -> L
    feats[730] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-296') * tryProb_catchZero(prob_cng2cng_linear, '-296^28') * tryProb_catchZero(prob_cng2lem_linear, '28^' + node2.lemma);
    
    # L -> -296 -> -24 -> L
    feats[731] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-296') * tryProb_catchZero(prob_cng2cng_linear, '-296^-24') * tryProb_catchZero(prob_cng2lem_linear, '-24^' + node2.lemma);
    
    # L -> -296 -> 137 -> L
    feats[732] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-296') * tryProb_catchZero(prob_cng2cng_linear, '-296^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> -296 -> -112 -> L
    feats[733] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-296') * tryProb_catchZero(prob_cng2cng_linear, '-296^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # L -> -296 -> 101 -> L
    feats[734] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-296') * tryProb_catchZero(prob_cng2cng_linear, '-296^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> -296 -> 91 -> L
    feats[735] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-296') * tryProb_catchZero(prob_cng2cng_linear, '-296^91') * tryProb_catchZero(prob_cng2lem_linear, '91^' + node2.lemma);
    
    # L -> -296 -> -35 -> L
    feats[736] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-296') * tryProb_catchZero(prob_cng2cng_linear, '-296^-35') * tryProb_catchZero(prob_cng2lem_linear, '-35^' + node2.lemma);
    
    # L -> -296 -> 3_du -> L
    feats[737] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-296') * tryProb_catchZero(prob_cng2cng_linear, '-296^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> 3_sp -> 5_sp -> L
    feats[738] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^3_sp') * tryProb_catchZero(prob_cng2cng_linear, '3_sp^5_sp') * tryProb_catchZero(prob_cng2lem_linear, '5_sp^' + node2.lemma);
    
    # L -> 3_sp -> dat. du. -> L
    feats[739] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^3_sp') * tryProb_catchZero(prob_cng2cng_linear, '3_sp^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # L -> 3_sp -> 8_fp -> L
    feats[740] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^3_sp') * tryProb_catchZero(prob_cng2cng_linear, '3_sp^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # L -> 3_sp -> 54 -> L
    feats[741] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^3_sp') * tryProb_catchZero(prob_cng2cng_linear, '3_sp^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> 3_sp -> pl -> L
    feats[742] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^3_sp') * tryProb_catchZero(prob_cng2cng_linear, '3_sp^pl') * tryProb_catchZero(prob_cng2lem_linear, 'pl^' + node2.lemma);
    
    # L -> 3_sp -> 153 -> L
    feats[743] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^3_sp') * tryProb_catchZero(prob_cng2cng_linear, '3_sp^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # L -> 3_sp -> 28 -> L
    feats[744] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^3_sp') * tryProb_catchZero(prob_cng2cng_linear, '3_sp^28') * tryProb_catchZero(prob_cng2lem_linear, '28^' + node2.lemma);
    
    # L -> 3_sp -> -24 -> L
    feats[745] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^3_sp') * tryProb_catchZero(prob_cng2cng_linear, '3_sp^-24') * tryProb_catchZero(prob_cng2lem_linear, '-24^' + node2.lemma);
    
    # L -> 3_sp -> 137 -> L
    feats[746] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^3_sp') * tryProb_catchZero(prob_cng2cng_linear, '3_sp^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> 3_sp -> -112 -> L
    feats[747] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^3_sp') * tryProb_catchZero(prob_cng2cng_linear, '3_sp^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # L -> 3_sp -> 101 -> L
    feats[748] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^3_sp') * tryProb_catchZero(prob_cng2cng_linear, '3_sp^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> 3_sp -> 91 -> L
    feats[749] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^3_sp') * tryProb_catchZero(prob_cng2cng_linear, '3_sp^91') * tryProb_catchZero(prob_cng2lem_linear, '91^' + node2.lemma);
    
    # L -> 3_sp -> -35 -> L
    feats[750] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^3_sp') * tryProb_catchZero(prob_cng2cng_linear, '3_sp^-35') * tryProb_catchZero(prob_cng2lem_linear, '-35^' + node2.lemma);
    
    # L -> 3_sp -> 3_du -> L
    feats[751] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^3_sp') * tryProb_catchZero(prob_cng2cng_linear, '3_sp^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> 120 -> 101 -> L
    feats[752] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^120') * tryProb_catchZero(prob_cng2cng_linear, '120^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> 94 -> 5_sp -> L
    feats[753] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^94') * tryProb_catchZero(prob_cng2cng_linear, '94^5_sp') * tryProb_catchZero(prob_cng2lem_linear, '5_sp^' + node2.lemma);
    
    # L -> 94 -> dat. du. -> L
    feats[754] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^94') * tryProb_catchZero(prob_cng2cng_linear, '94^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # L -> 94 -> 8_fp -> L
    feats[755] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^94') * tryProb_catchZero(prob_cng2cng_linear, '94^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # L -> 94 -> 54 -> L
    feats[756] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^94') * tryProb_catchZero(prob_cng2cng_linear, '94^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> 94 -> 153 -> L
    feats[757] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^94') * tryProb_catchZero(prob_cng2cng_linear, '94^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # L -> 94 -> -24 -> L
    feats[758] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^94') * tryProb_catchZero(prob_cng2cng_linear, '94^-24') * tryProb_catchZero(prob_cng2lem_linear, '-24^' + node2.lemma);
    
    # L -> 94 -> 137 -> L
    feats[759] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^94') * tryProb_catchZero(prob_cng2cng_linear, '94^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> 94 -> -112 -> L
    feats[760] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^94') * tryProb_catchZero(prob_cng2cng_linear, '94^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # L -> 94 -> 101 -> L
    feats[761] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^94') * tryProb_catchZero(prob_cng2cng_linear, '94^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> 94 -> 91 -> L
    feats[762] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^94') * tryProb_catchZero(prob_cng2cng_linear, '94^91') * tryProb_catchZero(prob_cng2lem_linear, '91^' + node2.lemma);
    
    # L -> 94 -> 3_du -> L
    feats[763] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^94') * tryProb_catchZero(prob_cng2cng_linear, '94^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> -200 -> 5_sp -> L
    feats[764] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-200') * tryProb_catchZero(prob_cng2cng_linear, '-200^5_sp') * tryProb_catchZero(prob_cng2lem_linear, '5_sp^' + node2.lemma);
    
    # L -> -200 -> dat. du. -> L
    feats[765] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-200') * tryProb_catchZero(prob_cng2cng_linear, '-200^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # L -> -200 -> 8_fp -> L
    feats[766] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-200') * tryProb_catchZero(prob_cng2cng_linear, '-200^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # L -> -200 -> 54 -> L
    feats[767] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-200') * tryProb_catchZero(prob_cng2cng_linear, '-200^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> -200 -> pl -> L
    feats[768] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-200') * tryProb_catchZero(prob_cng2cng_linear, '-200^pl') * tryProb_catchZero(prob_cng2lem_linear, 'pl^' + node2.lemma);
    
    # L -> -200 -> 153 -> L
    feats[769] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-200') * tryProb_catchZero(prob_cng2cng_linear, '-200^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # L -> -200 -> 28 -> L
    feats[770] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-200') * tryProb_catchZero(prob_cng2cng_linear, '-200^28') * tryProb_catchZero(prob_cng2lem_linear, '28^' + node2.lemma);
    
    # L -> -200 -> 137 -> L
    feats[771] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-200') * tryProb_catchZero(prob_cng2cng_linear, '-200^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> -200 -> -112 -> L
    feats[772] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-200') * tryProb_catchZero(prob_cng2cng_linear, '-200^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # L -> -200 -> 101 -> L
    feats[773] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-200') * tryProb_catchZero(prob_cng2cng_linear, '-200^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> -200 -> -35 -> L
    feats[774] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-200') * tryProb_catchZero(prob_cng2cng_linear, '-200^-35') * tryProb_catchZero(prob_cng2lem_linear, '-35^' + node2.lemma);
    
    # L -> -200 -> 3_du -> L
    feats[775] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-200') * tryProb_catchZero(prob_cng2cng_linear, '-200^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> -81 -> 5_sp -> L
    feats[776] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-81') * tryProb_catchZero(prob_cng2cng_linear, '-81^5_sp') * tryProb_catchZero(prob_cng2lem_linear, '5_sp^' + node2.lemma);
    
    # L -> -81 -> dat. du. -> L
    feats[777] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-81') * tryProb_catchZero(prob_cng2cng_linear, '-81^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # L -> -81 -> 8_fp -> L
    feats[778] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-81') * tryProb_catchZero(prob_cng2cng_linear, '-81^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # L -> -81 -> 54 -> L
    feats[779] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-81') * tryProb_catchZero(prob_cng2cng_linear, '-81^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # L -> -81 -> pl -> L
    feats[780] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-81') * tryProb_catchZero(prob_cng2cng_linear, '-81^pl') * tryProb_catchZero(prob_cng2lem_linear, 'pl^' + node2.lemma);
    
    # L -> -81 -> 153 -> L
    feats[781] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-81') * tryProb_catchZero(prob_cng2cng_linear, '-81^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # L -> -81 -> 28 -> L
    feats[782] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-81') * tryProb_catchZero(prob_cng2cng_linear, '-81^28') * tryProb_catchZero(prob_cng2lem_linear, '28^' + node2.lemma);
    
    # L -> -81 -> 137 -> L
    feats[783] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-81') * tryProb_catchZero(prob_cng2cng_linear, '-81^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # L -> -81 -> -112 -> L
    feats[784] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-81') * tryProb_catchZero(prob_cng2cng_linear, '-81^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # L -> -81 -> 101 -> L
    feats[785] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-81') * tryProb_catchZero(prob_cng2cng_linear, '-81^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # L -> -81 -> -35 -> L
    feats[786] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-81') * tryProb_catchZero(prob_cng2cng_linear, '-81^-35') * tryProb_catchZero(prob_cng2lem_linear, '-35^' + node2.lemma);
    
    # L -> -81 -> 3_du -> L
    feats[787] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-81') * tryProb_catchZero(prob_cng2cng_linear, '-81^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # L -> tp -> dat. du. -> T
    feats[788] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^tp') * tryProb_catchZero(prob_cng2cng_linear, 'tp^dat. du.') * tryProb_catchZero(prob_cng2tup_linear, 'dat. du.^' + node2.tup);
    
    # L -> tp -> 101 -> T
    feats[789] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^tp') * tryProb_catchZero(prob_cng2cng_linear, 'tp^101') * tryProb_catchZero(prob_cng2tup_linear, '101^' + node2.tup);
    
    # L -> 5_du -> 153 -> T
    feats[790] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^5_du') * tryProb_catchZero(prob_cng2cng_linear, '5_du^153') * tryProb_catchZero(prob_cng2tup_linear, '153^' + node2.tup);
    
    # L -> 157 -> 54 -> T
    feats[791] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^157') * tryProb_catchZero(prob_cng2cng_linear, '157^54') * tryProb_catchZero(prob_cng2tup_linear, '54^' + node2.tup);
    
    # L -> 157 -> 137 -> T
    feats[792] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^157') * tryProb_catchZero(prob_cng2cng_linear, '157^137') * tryProb_catchZero(prob_cng2tup_linear, '137^' + node2.tup);
    
    # L -> 157 -> 101 -> T
    feats[793] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^157') * tryProb_catchZero(prob_cng2cng_linear, '157^101') * tryProb_catchZero(prob_cng2tup_linear, '101^' + node2.tup);
    
    # L -> 157 -> 3_du -> T
    feats[794] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^157') * tryProb_catchZero(prob_cng2cng_linear, '157^3_du') * tryProb_catchZero(prob_cng2tup_linear, '3_du^' + node2.tup);
    
    # L -> 5_sp -> -90 -> T
    feats[795] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^5_sp') * tryProb_catchZero(prob_cng2cng_linear, '5_sp^-90') * tryProb_catchZero(prob_cng2tup_linear, '-90^' + node2.tup);
    
    # L -> 5_sp -> 91 -> T
    feats[796] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^5_sp') * tryProb_catchZero(prob_cng2cng_linear, '5_sp^91') * tryProb_catchZero(prob_cng2tup_linear, '91^' + node2.tup);
    
    # L -> -151 -> -24 -> T
    feats[797] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-151') * tryProb_catchZero(prob_cng2cng_linear, '-151^-24') * tryProb_catchZero(prob_cng2tup_linear, '-24^' + node2.tup);
    
    # L -> acc. sg. -> 101 -> T
    feats[798] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^acc. sg.') * tryProb_catchZero(prob_cng2cng_linear, 'acc. sg.^101') * tryProb_catchZero(prob_cng2tup_linear, '101^' + node2.tup);
    
    # L -> -25 -> 137 -> T
    feats[799] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-25') * tryProb_catchZero(prob_cng2cng_linear, '-25^137') * tryProb_catchZero(prob_cng2tup_linear, '137^' + node2.tup);
    
    # L -> -25 -> 101 -> T
    feats[800] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-25') * tryProb_catchZero(prob_cng2cng_linear, '-25^101') * tryProb_catchZero(prob_cng2tup_linear, '101^' + node2.tup);
    
    # L -> -57 -> 91 -> T
    feats[801] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-57') * tryProb_catchZero(prob_cng2cng_linear, '-57^91') * tryProb_catchZero(prob_cng2tup_linear, '91^' + node2.tup);
    
    # L -> 8_fp -> 91 -> T
    feats[802] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^8_fp') * tryProb_catchZero(prob_cng2cng_linear, '8_fp^91') * tryProb_catchZero(prob_cng2tup_linear, '91^' + node2.tup);
    
    # L -> 54 -> 8_fp -> T
    feats[803] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^54') * tryProb_catchZero(prob_cng2cng_linear, '54^8_fp') * tryProb_catchZero(prob_cng2tup_linear, '8_fp^' + node2.tup);
    
    # L -> 55 -> 101 -> T
    feats[804] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^55') * tryProb_catchZero(prob_cng2cng_linear, '55^101') * tryProb_catchZero(prob_cng2tup_linear, '101^' + node2.tup);
    
    # L -> -156 -> 101 -> T
    feats[805] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-156') * tryProb_catchZero(prob_cng2cng_linear, '-156^101') * tryProb_catchZero(prob_cng2tup_linear, '101^' + node2.tup);
    
    # L -> -29 -> 54 -> T
    feats[806] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-29') * tryProb_catchZero(prob_cng2cng_linear, '-29^54') * tryProb_catchZero(prob_cng2tup_linear, '54^' + node2.tup);
    
    # L -> -29 -> 137 -> T
    feats[807] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-29') * tryProb_catchZero(prob_cng2cng_linear, '-29^137') * tryProb_catchZero(prob_cng2tup_linear, '137^' + node2.tup);
    
    # L -> -29 -> 101 -> T
    feats[808] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-29') * tryProb_catchZero(prob_cng2cng_linear, '-29^101') * tryProb_catchZero(prob_cng2tup_linear, '101^' + node2.tup);
    
    # L -> 153 -> -24 -> T
    feats[809] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^153') * tryProb_catchZero(prob_cng2cng_linear, '153^-24') * tryProb_catchZero(prob_cng2tup_linear, '-24^' + node2.tup);
    
    # L -> -24 -> pl -> T
    feats[810] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^-24') * tryProb_catchZero(prob_cng2cng_linear, '-24^pl') * tryProb_catchZero(prob_cng2tup_linear, 'pl^' + node2.tup);
    
    # L -> acc -> dat. du. -> T
    feats[811] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^acc') * tryProb_catchZero(prob_cng2cng_linear, 'acc^dat. du.') * tryProb_catchZero(prob_cng2tup_linear, 'dat. du.^' + node2.tup);
    
    # L -> acc -> 101 -> T
    feats[812] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^acc') * tryProb_catchZero(prob_cng2cng_linear, 'acc^101') * tryProb_catchZero(prob_cng2tup_linear, '101^' + node2.tup);
    
    # L -> 13_pl -> 54 -> T
    feats[813] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^13_pl') * tryProb_catchZero(prob_cng2cng_linear, '13_pl^54') * tryProb_catchZero(prob_cng2tup_linear, '54^' + node2.tup);
    
    # L -> 101 -> -151 -> T
    feats[814] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^101') * tryProb_catchZero(prob_cng2cng_linear, '101^-151') * tryProb_catchZero(prob_cng2tup_linear, '-151^' + node2.tup);
    
    # L -> 101 -> 8_fp -> T
    feats[815] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^101') * tryProb_catchZero(prob_cng2cng_linear, '101^8_fp') * tryProb_catchZero(prob_cng2tup_linear, '8_fp^' + node2.tup);
    
    # L -> 91 -> 91 -> T
    feats[816] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^91') * tryProb_catchZero(prob_cng2cng_linear, '91^91') * tryProb_catchZero(prob_cng2tup_linear, '91^' + node2.tup);
    
    # L -> 3_du -> 91 -> T
    feats[817] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^3_du') * tryProb_catchZero(prob_cng2cng_linear, '3_du^91') * tryProb_catchZero(prob_cng2tup_linear, '91^' + node2.tup);
    
    # L -> 35 -> 8_fp -> T
    feats[818] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^35') * tryProb_catchZero(prob_cng2cng_linear, '35^8_fp') * tryProb_catchZero(prob_cng2tup_linear, '8_fp^' + node2.tup);
    
    # L -> 35 -> 54 -> T
    feats[819] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^35') * tryProb_catchZero(prob_cng2cng_linear, '35^54') * tryProb_catchZero(prob_cng2tup_linear, '54^' + node2.tup);
    
    # L -> 35 -> pl -> T
    feats[820] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^35') * tryProb_catchZero(prob_cng2cng_linear, '35^pl') * tryProb_catchZero(prob_cng2tup_linear, 'pl^' + node2.tup);
    
    # L -> 35 -> 101 -> T
    feats[821] = tryProb_catchZero(prob_lem2cng_linear, node1.lemma + '^35') * tryProb_catchZero(prob_cng2cng_linear, '35^101') * tryProb_catchZero(prob_cng2tup_linear, '101^' + node2.tup);
    
    # T -> gen. sg. -> 137 -> L
    feats[822] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^gen. sg.') * tryProb_catchZero(prob_cng2cng_linear, 'gen. sg.^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # T -> gen. sg. -> 101 -> L
    feats[823] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^gen. sg.') * tryProb_catchZero(prob_cng2cng_linear, 'gen. sg.^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # T -> 5_sp -> 54 -> L
    feats[824] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^5_sp') * tryProb_catchZero(prob_cng2cng_linear, '5_sp^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # T -> 5_sp -> -24 -> L
    feats[825] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^5_sp') * tryProb_catchZero(prob_cng2cng_linear, '5_sp^-24') * tryProb_catchZero(prob_cng2lem_linear, '-24^' + node2.lemma);
    
    # T -> 5_sp -> 137 -> L
    feats[826] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^5_sp') * tryProb_catchZero(prob_cng2cng_linear, '5_sp^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # T -> 5_sp -> -112 -> L
    feats[827] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^5_sp') * tryProb_catchZero(prob_cng2cng_linear, '5_sp^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # T -> 5_sp -> 101 -> L
    feats[828] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^5_sp') * tryProb_catchZero(prob_cng2cng_linear, '5_sp^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # T -> 5_sp -> 91 -> L
    feats[829] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^5_sp') * tryProb_catchZero(prob_cng2cng_linear, '5_sp^91') * tryProb_catchZero(prob_cng2lem_linear, '91^' + node2.lemma);
    
    # T -> 5_sp -> 3_du -> L
    feats[830] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^5_sp') * tryProb_catchZero(prob_cng2cng_linear, '5_sp^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # T -> dat. du. -> 54 -> L
    feats[831] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^dat. du.') * tryProb_catchZero(prob_cng2cng_linear, 'dat. du.^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # T -> 15_du -> dat. du. -> L
    feats[832] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^15_du') * tryProb_catchZero(prob_cng2cng_linear, '15_du^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # T -> 15_du -> 8_fp -> L
    feats[833] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^15_du') * tryProb_catchZero(prob_cng2cng_linear, '15_du^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # T -> 15_du -> 54 -> L
    feats[834] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^15_du') * tryProb_catchZero(prob_cng2cng_linear, '15_du^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # T -> 15_du -> pl -> L
    feats[835] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^15_du') * tryProb_catchZero(prob_cng2cng_linear, '15_du^pl') * tryProb_catchZero(prob_cng2lem_linear, 'pl^' + node2.lemma);
    
    # T -> 15_du -> 153 -> L
    feats[836] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^15_du') * tryProb_catchZero(prob_cng2cng_linear, '15_du^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # T -> 15_du -> 28 -> L
    feats[837] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^15_du') * tryProb_catchZero(prob_cng2cng_linear, '15_du^28') * tryProb_catchZero(prob_cng2lem_linear, '28^' + node2.lemma);
    
    # T -> 15_du -> 137 -> L
    feats[838] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^15_du') * tryProb_catchZero(prob_cng2cng_linear, '15_du^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # T -> 15_du -> -112 -> L
    feats[839] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^15_du') * tryProb_catchZero(prob_cng2cng_linear, '15_du^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # T -> 15_du -> 101 -> L
    feats[840] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^15_du') * tryProb_catchZero(prob_cng2cng_linear, '15_du^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # T -> 15_du -> 3_du -> L
    feats[841] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^15_du') * tryProb_catchZero(prob_cng2cng_linear, '15_du^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # T -> -151 -> dat. du. -> L
    feats[842] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^-151') * tryProb_catchZero(prob_cng2cng_linear, '-151^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # T -> -151 -> 8_fp -> L
    feats[843] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^-151') * tryProb_catchZero(prob_cng2cng_linear, '-151^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # T -> -151 -> 54 -> L
    feats[844] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^-151') * tryProb_catchZero(prob_cng2cng_linear, '-151^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # T -> -151 -> pl -> L
    feats[845] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^-151') * tryProb_catchZero(prob_cng2cng_linear, '-151^pl') * tryProb_catchZero(prob_cng2lem_linear, 'pl^' + node2.lemma);
    
    # T -> -151 -> 153 -> L
    feats[846] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^-151') * tryProb_catchZero(prob_cng2cng_linear, '-151^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # T -> -151 -> 28 -> L
    feats[847] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^-151') * tryProb_catchZero(prob_cng2cng_linear, '-151^28') * tryProb_catchZero(prob_cng2lem_linear, '28^' + node2.lemma);
    
    # T -> -151 -> 137 -> L
    feats[848] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^-151') * tryProb_catchZero(prob_cng2cng_linear, '-151^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # T -> -151 -> -112 -> L
    feats[849] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^-151') * tryProb_catchZero(prob_cng2cng_linear, '-151^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # T -> -151 -> 101 -> L
    feats[850] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^-151') * tryProb_catchZero(prob_cng2cng_linear, '-151^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # T -> -151 -> -35 -> L
    feats[851] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^-151') * tryProb_catchZero(prob_cng2cng_linear, '-151^-35') * tryProb_catchZero(prob_cng2lem_linear, '-35^' + node2.lemma);
    
    # T -> -151 -> 3_du -> L
    feats[852] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^-151') * tryProb_catchZero(prob_cng2cng_linear, '-151^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # T -> -306 -> 137 -> L
    feats[853] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^-306') * tryProb_catchZero(prob_cng2cng_linear, '-306^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # T -> -306 -> 101 -> L
    feats[854] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^-306') * tryProb_catchZero(prob_cng2cng_linear, '-306^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # T -> -90 -> dat. du. -> L
    feats[855] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^-90') * tryProb_catchZero(prob_cng2cng_linear, '-90^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # T -> -90 -> 8_fp -> L
    feats[856] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^-90') * tryProb_catchZero(prob_cng2cng_linear, '-90^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # T -> -90 -> 54 -> L
    feats[857] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^-90') * tryProb_catchZero(prob_cng2cng_linear, '-90^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # T -> -90 -> 153 -> L
    feats[858] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^-90') * tryProb_catchZero(prob_cng2cng_linear, '-90^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # T -> -90 -> 137 -> L
    feats[859] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^-90') * tryProb_catchZero(prob_cng2cng_linear, '-90^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # T -> -90 -> -112 -> L
    feats[860] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^-90') * tryProb_catchZero(prob_cng2cng_linear, '-90^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # T -> -90 -> 101 -> L
    feats[861] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^-90') * tryProb_catchZero(prob_cng2cng_linear, '-90^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # T -> -90 -> 3_du -> L
    feats[862] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^-90') * tryProb_catchZero(prob_cng2cng_linear, '-90^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # T -> -57 -> 8_fp -> L
    feats[863] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^-57') * tryProb_catchZero(prob_cng2cng_linear, '-57^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # T -> -57 -> 54 -> L
    feats[864] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^-57') * tryProb_catchZero(prob_cng2cng_linear, '-57^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # T -> -57 -> pl -> L
    feats[865] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^-57') * tryProb_catchZero(prob_cng2cng_linear, '-57^pl') * tryProb_catchZero(prob_cng2lem_linear, 'pl^' + node2.lemma);
    
    # T -> -57 -> 153 -> L
    feats[866] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^-57') * tryProb_catchZero(prob_cng2cng_linear, '-57^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # T -> -57 -> 137 -> L
    feats[867] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^-57') * tryProb_catchZero(prob_cng2cng_linear, '-57^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # T -> -57 -> -112 -> L
    feats[868] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^-57') * tryProb_catchZero(prob_cng2cng_linear, '-57^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # T -> -57 -> 101 -> L
    feats[869] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^-57') * tryProb_catchZero(prob_cng2cng_linear, '-57^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # T -> 8_fp -> dat. du. -> L
    feats[870] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^8_fp') * tryProb_catchZero(prob_cng2cng_linear, '8_fp^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # T -> 8_fp -> 8_fp -> L
    feats[871] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^8_fp') * tryProb_catchZero(prob_cng2cng_linear, '8_fp^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # T -> 8_fp -> 54 -> L
    feats[872] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^8_fp') * tryProb_catchZero(prob_cng2cng_linear, '8_fp^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # T -> 8_fp -> pl -> L
    feats[873] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^8_fp') * tryProb_catchZero(prob_cng2cng_linear, '8_fp^pl') * tryProb_catchZero(prob_cng2lem_linear, 'pl^' + node2.lemma);
    
    # T -> 8_fp -> 153 -> L
    feats[874] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^8_fp') * tryProb_catchZero(prob_cng2cng_linear, '8_fp^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # T -> 8_fp -> 28 -> L
    feats[875] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^8_fp') * tryProb_catchZero(prob_cng2cng_linear, '8_fp^28') * tryProb_catchZero(prob_cng2lem_linear, '28^' + node2.lemma);
    
    # T -> 8_fp -> -24 -> L
    feats[876] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^8_fp') * tryProb_catchZero(prob_cng2cng_linear, '8_fp^-24') * tryProb_catchZero(prob_cng2lem_linear, '-24^' + node2.lemma);
    
    # T -> 8_fp -> 137 -> L
    feats[877] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^8_fp') * tryProb_catchZero(prob_cng2cng_linear, '8_fp^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # T -> 8_fp -> -112 -> L
    feats[878] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^8_fp') * tryProb_catchZero(prob_cng2cng_linear, '8_fp^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # T -> 8_fp -> 101 -> L
    feats[879] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^8_fp') * tryProb_catchZero(prob_cng2cng_linear, '8_fp^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # T -> 8_fp -> 91 -> L
    feats[880] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^8_fp') * tryProb_catchZero(prob_cng2cng_linear, '8_fp^91') * tryProb_catchZero(prob_cng2lem_linear, '91^' + node2.lemma);
    
    # T -> 8_fp -> -35 -> L
    feats[881] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^8_fp') * tryProb_catchZero(prob_cng2cng_linear, '8_fp^-35') * tryProb_catchZero(prob_cng2lem_linear, '-35^' + node2.lemma);
    
    # T -> 8_fp -> 3_du -> L
    feats[882] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^8_fp') * tryProb_catchZero(prob_cng2cng_linear, '8_fp^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # T -> 169 -> dat. du. -> L
    feats[883] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^169') * tryProb_catchZero(prob_cng2cng_linear, '169^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # T -> 169 -> 8_fp -> L
    feats[884] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^169') * tryProb_catchZero(prob_cng2cng_linear, '169^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # T -> 54 -> dat. du. -> L
    feats[885] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^54') * tryProb_catchZero(prob_cng2cng_linear, '54^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # T -> 54 -> -151 -> L
    feats[886] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^54') * tryProb_catchZero(prob_cng2cng_linear, '54^-151') * tryProb_catchZero(prob_cng2lem_linear, '-151^' + node2.lemma);
    
    # T -> 54 -> -90 -> L
    feats[887] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^54') * tryProb_catchZero(prob_cng2cng_linear, '54^-90') * tryProb_catchZero(prob_cng2lem_linear, '-90^' + node2.lemma);
    
    # T -> 54 -> 8_fp -> L
    feats[888] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^54') * tryProb_catchZero(prob_cng2cng_linear, '54^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # T -> 54 -> 54 -> L
    feats[889] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^54') * tryProb_catchZero(prob_cng2cng_linear, '54^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # T -> 54 -> pl -> L
    feats[890] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^54') * tryProb_catchZero(prob_cng2cng_linear, '54^pl') * tryProb_catchZero(prob_cng2lem_linear, 'pl^' + node2.lemma);
    
    # T -> 54 -> 153 -> L
    feats[891] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^54') * tryProb_catchZero(prob_cng2cng_linear, '54^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # T -> 54 -> 135 -> L
    feats[892] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^54') * tryProb_catchZero(prob_cng2cng_linear, '54^135') * tryProb_catchZero(prob_cng2lem_linear, '135^' + node2.lemma);
    
    # T -> 54 -> -24 -> L
    feats[893] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^54') * tryProb_catchZero(prob_cng2cng_linear, '54^-24') * tryProb_catchZero(prob_cng2lem_linear, '-24^' + node2.lemma);
    
    # T -> 54 -> 137 -> L
    feats[894] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^54') * tryProb_catchZero(prob_cng2cng_linear, '54^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # T -> 54 -> -112 -> L
    feats[895] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^54') * tryProb_catchZero(prob_cng2cng_linear, '54^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # T -> 54 -> 101 -> L
    feats[896] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^54') * tryProb_catchZero(prob_cng2cng_linear, '54^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # T -> 54 -> 91 -> L
    feats[897] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^54') * tryProb_catchZero(prob_cng2cng_linear, '54^91') * tryProb_catchZero(prob_cng2lem_linear, '91^' + node2.lemma);
    
    # T -> 54 -> -35 -> L
    feats[898] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^54') * tryProb_catchZero(prob_cng2cng_linear, '54^-35') * tryProb_catchZero(prob_cng2lem_linear, '-35^' + node2.lemma);
    
    # T -> 54 -> 3_du -> L
    feats[899] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^54') * tryProb_catchZero(prob_cng2cng_linear, '54^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # T -> 114 -> 54 -> L
    feats[900] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^114') * tryProb_catchZero(prob_cng2cng_linear, '114^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # T -> 114 -> 101 -> L
    feats[901] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^114') * tryProb_catchZero(prob_cng2cng_linear, '114^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # T -> -139 -> 137 -> L
    feats[902] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^-139') * tryProb_catchZero(prob_cng2cng_linear, '-139^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # T -> -139 -> 101 -> L
    feats[903] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^-139') * tryProb_catchZero(prob_cng2cng_linear, '-139^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # T -> 139 -> dat. du. -> L
    feats[904] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^139') * tryProb_catchZero(prob_cng2cng_linear, '139^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # T -> 139 -> 8_fp -> L
    feats[905] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^139') * tryProb_catchZero(prob_cng2cng_linear, '139^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # T -> 139 -> 54 -> L
    feats[906] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^139') * tryProb_catchZero(prob_cng2cng_linear, '139^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # T -> 139 -> pl -> L
    feats[907] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^139') * tryProb_catchZero(prob_cng2cng_linear, '139^pl') * tryProb_catchZero(prob_cng2lem_linear, 'pl^' + node2.lemma);
    
    # T -> 139 -> 153 -> L
    feats[908] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^139') * tryProb_catchZero(prob_cng2cng_linear, '139^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # T -> 139 -> 28 -> L
    feats[909] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^139') * tryProb_catchZero(prob_cng2cng_linear, '139^28') * tryProb_catchZero(prob_cng2lem_linear, '28^' + node2.lemma);
    
    # T -> 139 -> 137 -> L
    feats[910] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^139') * tryProb_catchZero(prob_cng2cng_linear, '139^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # T -> 139 -> -112 -> L
    feats[911] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^139') * tryProb_catchZero(prob_cng2cng_linear, '139^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # T -> 139 -> 101 -> L
    feats[912] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^139') * tryProb_catchZero(prob_cng2cng_linear, '139^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # T -> 139 -> -35 -> L
    feats[913] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^139') * tryProb_catchZero(prob_cng2cng_linear, '139^-35') * tryProb_catchZero(prob_cng2lem_linear, '-35^' + node2.lemma);
    
    # T -> 139 -> 3_du -> L
    feats[914] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^139') * tryProb_catchZero(prob_cng2cng_linear, '139^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # T -> pl -> -90 -> L
    feats[915] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^pl') * tryProb_catchZero(prob_cng2cng_linear, 'pl^-90') * tryProb_catchZero(prob_cng2lem_linear, '-90^' + node2.lemma);
    
    # T -> pl -> 8_fp -> L
    feats[916] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^pl') * tryProb_catchZero(prob_cng2cng_linear, 'pl^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # T -> pl -> 54 -> L
    feats[917] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^pl') * tryProb_catchZero(prob_cng2cng_linear, 'pl^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # T -> pl -> -24 -> L
    feats[918] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^pl') * tryProb_catchZero(prob_cng2cng_linear, 'pl^-24') * tryProb_catchZero(prob_cng2lem_linear, '-24^' + node2.lemma);
    
    # T -> pl -> -112 -> L
    feats[919] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^pl') * tryProb_catchZero(prob_cng2cng_linear, 'pl^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # T -> pl -> 101 -> L
    feats[920] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^pl') * tryProb_catchZero(prob_cng2cng_linear, 'pl^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # T -> 153 -> dat. du. -> L
    feats[921] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^153') * tryProb_catchZero(prob_cng2cng_linear, '153^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # T -> 153 -> 8_fp -> L
    feats[922] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^153') * tryProb_catchZero(prob_cng2cng_linear, '153^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # T -> 153 -> 54 -> L
    feats[923] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^153') * tryProb_catchZero(prob_cng2cng_linear, '153^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # T -> 153 -> pl -> L
    feats[924] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^153') * tryProb_catchZero(prob_cng2cng_linear, '153^pl') * tryProb_catchZero(prob_cng2lem_linear, 'pl^' + node2.lemma);
    
    # T -> 153 -> 153 -> L
    feats[925] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^153') * tryProb_catchZero(prob_cng2cng_linear, '153^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # T -> 153 -> 28 -> L
    feats[926] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^153') * tryProb_catchZero(prob_cng2cng_linear, '153^28') * tryProb_catchZero(prob_cng2lem_linear, '28^' + node2.lemma);
    
    # T -> 153 -> -24 -> L
    feats[927] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^153') * tryProb_catchZero(prob_cng2cng_linear, '153^-24') * tryProb_catchZero(prob_cng2lem_linear, '-24^' + node2.lemma);
    
    # T -> 153 -> 137 -> L
    feats[928] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^153') * tryProb_catchZero(prob_cng2cng_linear, '153^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # T -> 153 -> -112 -> L
    feats[929] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^153') * tryProb_catchZero(prob_cng2cng_linear, '153^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # T -> 153 -> 101 -> L
    feats[930] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^153') * tryProb_catchZero(prob_cng2cng_linear, '153^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # T -> 153 -> 91 -> L
    feats[931] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^153') * tryProb_catchZero(prob_cng2cng_linear, '153^91') * tryProb_catchZero(prob_cng2lem_linear, '91^' + node2.lemma);
    
    # T -> 153 -> -35 -> L
    feats[932] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^153') * tryProb_catchZero(prob_cng2cng_linear, '153^-35') * tryProb_catchZero(prob_cng2lem_linear, '-35^' + node2.lemma);
    
    # T -> 153 -> 3_du -> L
    feats[933] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^153') * tryProb_catchZero(prob_cng2cng_linear, '153^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # T -> 135 -> 54 -> L
    feats[934] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^135') * tryProb_catchZero(prob_cng2cng_linear, '135^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # T -> 135 -> pl -> L
    feats[935] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^135') * tryProb_catchZero(prob_cng2cng_linear, '135^pl') * tryProb_catchZero(prob_cng2lem_linear, 'pl^' + node2.lemma);
    
    # T -> 135 -> 137 -> L
    feats[936] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^135') * tryProb_catchZero(prob_cng2cng_linear, '135^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # T -> 135 -> 3_du -> L
    feats[937] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^135') * tryProb_catchZero(prob_cng2cng_linear, '135^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # T -> 7_fp -> 101 -> L
    feats[938] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^7_fp') * tryProb_catchZero(prob_cng2cng_linear, '7_fp^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # T -> neutr -> 137 -> L
    feats[939] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^neutr') * tryProb_catchZero(prob_cng2cng_linear, 'neutr^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # T -> neutr -> 101 -> L
    feats[940] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^neutr') * tryProb_catchZero(prob_cng2cng_linear, 'neutr^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # T -> -153 -> dat. du. -> L
    feats[941] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^-153') * tryProb_catchZero(prob_cng2cng_linear, '-153^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # T -> -153 -> 8_fp -> L
    feats[942] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^-153') * tryProb_catchZero(prob_cng2cng_linear, '-153^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # T -> -153 -> 54 -> L
    feats[943] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^-153') * tryProb_catchZero(prob_cng2cng_linear, '-153^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # T -> -153 -> pl -> L
    feats[944] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^-153') * tryProb_catchZero(prob_cng2cng_linear, '-153^pl') * tryProb_catchZero(prob_cng2lem_linear, 'pl^' + node2.lemma);
    
    # T -> -153 -> 28 -> L
    feats[945] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^-153') * tryProb_catchZero(prob_cng2cng_linear, '-153^28') * tryProb_catchZero(prob_cng2lem_linear, '28^' + node2.lemma);
    
    # T -> -153 -> 137 -> L
    feats[946] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^-153') * tryProb_catchZero(prob_cng2cng_linear, '-153^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # T -> -153 -> 101 -> L
    feats[947] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^-153') * tryProb_catchZero(prob_cng2cng_linear, '-153^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # T -> -153 -> -35 -> L
    feats[948] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^-153') * tryProb_catchZero(prob_cng2cng_linear, '-153^-35') * tryProb_catchZero(prob_cng2lem_linear, '-35^' + node2.lemma);
    
    # T -> -153 -> 3_du -> L
    feats[949] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^-153') * tryProb_catchZero(prob_cng2cng_linear, '-153^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # T -> -24 -> 8_fp -> L
    feats[950] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^-24') * tryProb_catchZero(prob_cng2cng_linear, '-24^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # T -> -24 -> 54 -> L
    feats[951] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^-24') * tryProb_catchZero(prob_cng2cng_linear, '-24^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # T -> -24 -> pl -> L
    feats[952] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^-24') * tryProb_catchZero(prob_cng2cng_linear, '-24^pl') * tryProb_catchZero(prob_cng2lem_linear, 'pl^' + node2.lemma);
    
    # T -> -24 -> 137 -> L
    feats[953] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^-24') * tryProb_catchZero(prob_cng2cng_linear, '-24^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # T -> -24 -> -112 -> L
    feats[954] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^-24') * tryProb_catchZero(prob_cng2cng_linear, '-24^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # T -> -24 -> 101 -> L
    feats[955] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^-24') * tryProb_catchZero(prob_cng2cng_linear, '-24^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # T -> 137 -> dat. du. -> L
    feats[956] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^137') * tryProb_catchZero(prob_cng2cng_linear, '137^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # T -> 137 -> -151 -> L
    feats[957] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^137') * tryProb_catchZero(prob_cng2cng_linear, '137^-151') * tryProb_catchZero(prob_cng2lem_linear, '-151^' + node2.lemma);
    
    # T -> 137 -> 54 -> L
    feats[958] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^137') * tryProb_catchZero(prob_cng2cng_linear, '137^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # T -> 137 -> 153 -> L
    feats[959] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^137') * tryProb_catchZero(prob_cng2cng_linear, '137^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # T -> 137 -> 28 -> L
    feats[960] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^137') * tryProb_catchZero(prob_cng2cng_linear, '137^28') * tryProb_catchZero(prob_cng2lem_linear, '28^' + node2.lemma);
    
    # T -> 137 -> -24 -> L
    feats[961] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^137') * tryProb_catchZero(prob_cng2cng_linear, '137^-24') * tryProb_catchZero(prob_cng2lem_linear, '-24^' + node2.lemma);
    
    # T -> 137 -> 137 -> L
    feats[962] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^137') * tryProb_catchZero(prob_cng2cng_linear, '137^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # T -> 137 -> -112 -> L
    feats[963] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^137') * tryProb_catchZero(prob_cng2cng_linear, '137^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # T -> 137 -> 101 -> L
    feats[964] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^137') * tryProb_catchZero(prob_cng2cng_linear, '137^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # T -> 137 -> 3_du -> L
    feats[965] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^137') * tryProb_catchZero(prob_cng2cng_linear, '137^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # T -> -112 -> 153 -> L
    feats[966] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^-112') * tryProb_catchZero(prob_cng2cng_linear, '-112^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # T -> -112 -> -112 -> L
    feats[967] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^-112') * tryProb_catchZero(prob_cng2cng_linear, '-112^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # T -> 101 -> dat. du. -> L
    feats[968] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^101') * tryProb_catchZero(prob_cng2cng_linear, '101^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # T -> 101 -> -151 -> L
    feats[969] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^101') * tryProb_catchZero(prob_cng2cng_linear, '101^-151') * tryProb_catchZero(prob_cng2lem_linear, '-151^' + node2.lemma);
    
    # T -> 101 -> -90 -> L
    feats[970] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^101') * tryProb_catchZero(prob_cng2cng_linear, '101^-90') * tryProb_catchZero(prob_cng2lem_linear, '-90^' + node2.lemma);
    
    # T -> 101 -> -57 -> L
    feats[971] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^101') * tryProb_catchZero(prob_cng2cng_linear, '101^-57') * tryProb_catchZero(prob_cng2lem_linear, '-57^' + node2.lemma);
    
    # T -> 101 -> 8_fp -> L
    feats[972] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^101') * tryProb_catchZero(prob_cng2cng_linear, '101^8_fp') * tryProb_catchZero(prob_cng2lem_linear, '8_fp^' + node2.lemma);
    
    # T -> 101 -> 169 -> L
    feats[973] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^101') * tryProb_catchZero(prob_cng2cng_linear, '101^169') * tryProb_catchZero(prob_cng2lem_linear, '169^' + node2.lemma);
    
    # T -> 101 -> 54 -> L
    feats[974] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^101') * tryProb_catchZero(prob_cng2cng_linear, '101^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # T -> 101 -> 114 -> L
    feats[975] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^101') * tryProb_catchZero(prob_cng2cng_linear, '101^114') * tryProb_catchZero(prob_cng2lem_linear, '114^' + node2.lemma);
    
    # T -> 101 -> pl -> L
    feats[976] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^101') * tryProb_catchZero(prob_cng2cng_linear, '101^pl') * tryProb_catchZero(prob_cng2lem_linear, 'pl^' + node2.lemma);
    
    # T -> 101 -> 153 -> L
    feats[977] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^101') * tryProb_catchZero(prob_cng2cng_linear, '101^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # T -> 101 -> 28 -> L
    feats[978] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^101') * tryProb_catchZero(prob_cng2cng_linear, '101^28') * tryProb_catchZero(prob_cng2lem_linear, '28^' + node2.lemma);
    
    # T -> 101 -> -24 -> L
    feats[979] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^101') * tryProb_catchZero(prob_cng2cng_linear, '101^-24') * tryProb_catchZero(prob_cng2lem_linear, '-24^' + node2.lemma);
    
    # T -> 101 -> 137 -> L
    feats[980] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^101') * tryProb_catchZero(prob_cng2cng_linear, '101^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # T -> 101 -> -112 -> L
    feats[981] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^101') * tryProb_catchZero(prob_cng2cng_linear, '101^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # T -> 101 -> 101 -> L
    feats[982] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^101') * tryProb_catchZero(prob_cng2cng_linear, '101^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # T -> 101 -> 91 -> L
    feats[983] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^101') * tryProb_catchZero(prob_cng2cng_linear, '101^91') * tryProb_catchZero(prob_cng2lem_linear, '91^' + node2.lemma);
    
    # T -> 101 -> -35 -> L
    feats[984] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^101') * tryProb_catchZero(prob_cng2cng_linear, '101^-35') * tryProb_catchZero(prob_cng2lem_linear, '-35^' + node2.lemma);
    
    # T -> 101 -> 3_du -> L
    feats[985] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^101') * tryProb_catchZero(prob_cng2cng_linear, '101^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # T -> 91 -> dat. du. -> L
    feats[986] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^91') * tryProb_catchZero(prob_cng2cng_linear, '91^dat. du.') * tryProb_catchZero(prob_cng2lem_linear, 'dat. du.^' + node2.lemma);
    
    # T -> 91 -> 54 -> L
    feats[987] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^91') * tryProb_catchZero(prob_cng2cng_linear, '91^54') * tryProb_catchZero(prob_cng2lem_linear, '54^' + node2.lemma);
    
    # T -> 91 -> 153 -> L
    feats[988] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^91') * tryProb_catchZero(prob_cng2cng_linear, '91^153') * tryProb_catchZero(prob_cng2lem_linear, '153^' + node2.lemma);
    
    # T -> 91 -> -24 -> L
    feats[989] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^91') * tryProb_catchZero(prob_cng2cng_linear, '91^-24') * tryProb_catchZero(prob_cng2lem_linear, '-24^' + node2.lemma);
    
    # T -> 91 -> 137 -> L
    feats[990] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^91') * tryProb_catchZero(prob_cng2cng_linear, '91^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # T -> 91 -> -112 -> L
    feats[991] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^91') * tryProb_catchZero(prob_cng2cng_linear, '91^-112') * tryProb_catchZero(prob_cng2lem_linear, '-112^' + node2.lemma);
    
    # T -> 91 -> 101 -> L
    feats[992] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^91') * tryProb_catchZero(prob_cng2cng_linear, '91^101') * tryProb_catchZero(prob_cng2lem_linear, '101^' + node2.lemma);
    
    # T -> 91 -> 91 -> L
    feats[993] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^91') * tryProb_catchZero(prob_cng2cng_linear, '91^91') * tryProb_catchZero(prob_cng2lem_linear, '91^' + node2.lemma);
    
    # T -> 91 -> -35 -> L
    feats[994] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^91') * tryProb_catchZero(prob_cng2cng_linear, '91^-35') * tryProb_catchZero(prob_cng2lem_linear, '-35^' + node2.lemma);
    
    # T -> 91 -> 3_du -> L
    feats[995] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^91') * tryProb_catchZero(prob_cng2cng_linear, '91^3_du') * tryProb_catchZero(prob_cng2lem_linear, '3_du^' + node2.lemma);
    
    # T -> 3_du -> 137 -> L
    feats[996] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^3_du') * tryProb_catchZero(prob_cng2cng_linear, '3_du^137') * tryProb_catchZero(prob_cng2lem_linear, '137^' + node2.lemma);
    
    # T -> 153 -> 91 -> T
    feats[997] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^153') * tryProb_catchZero(prob_cng2cng_linear, '153^91') * tryProb_catchZero(prob_cng2tup_linear, '91^' + node2.tup);
    
    # T -> 28 -> 91 -> T
    feats[998] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^28') * tryProb_catchZero(prob_cng2cng_linear, '28^91') * tryProb_catchZero(prob_cng2tup_linear, '91^' + node2.tup);
    
    # T -> -296 -> 8_fp -> T
    feats[999] = tryProb_catchZero(prob_tup2cng_linear, node1.tup + '^-296') * tryProb_catchZero(prob_cng2cng_linear, '-296^8_fp') * tryProb_catchZero(prob_cng2tup_linear, '8_fp^' + node2.tup);
    
    feats[feats < 1e-25] = 1e-25
    return -np.log10(feats)
    