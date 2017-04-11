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
    _edge_vector_dim = 1500
    _full_cnglist = list(mat_cngCount_1D)
    print(_cg_count)

"""
################################################################################################
######################  CREATING FEATURE MATRICES ##############################################
################################################################################################
"""
def tryProb_catchZero(mat2D, mat1D, key1, key2):
    if key1 in mat1D and key1 in mat2D:
        if key2 in mat2D[key1]:
            return float(mat2D[key1][key2])/mat1D[key1]
    return 0.0

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
    
    # L->L
    feats[0] = tryProb_catchZero(mat_lem2lem_countonly, mat_lemCount_1D, node1.lemma, node2.lemma)
            
    # L->T
    feats[1] = tryProb_catchZero(mat_lem2tup_countonly, mat_lemCount_1D, node1.lemma, node2.tup)
            
    # C->C
    feats[2] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, node2.cng)
            
    # L -> 3 -> -210 -> L
    feats[3] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3', '-210') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-210', node2.lemma)
    
    # L -> 80 -> fem -> L
    feats[4] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '80') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '80', 'fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'fem', node2.lemma)
    
    # L -> 3 -> -307 -> L
    feats[5] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3', '-307') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-307', node2.lemma)
    
    # L -> -123 -> -161 -> L
    feats[6] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-123') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-123', '-161') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-161', node2.lemma)
    
    # L -> 48 -> -269 -> L
    feats[7] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '48', '-269') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-269', node2.lemma)
    
    # L -> 48 -> 28_sg -> L
    feats[8] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '48', '28_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28_sg', node2.lemma)
    
    # L -> 9_pl -> 6_tp -> L
    feats[9] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_pl', '6_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_tp', node2.lemma)
    
    # L -> acc -> 6_tp -> L
    feats[10] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc', '6_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_tp', node2.lemma)
    
    # L -> acc -> 7_sp -> L
    feats[11] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc', '7_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_sp', node2.lemma)
    
    # L -> 9_pl -> -15 -> L
    feats[12] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_pl', '-15') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-15', node2.lemma)
    
    # L -> acc -> -15 -> L
    feats[13] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc', '-15') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-15', node2.lemma)
    
    # L -> 9_pl -> -76 -> L
    feats[14] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_pl', '-76') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-76', node2.lemma)
    
    # L->-98->L
    feats[15] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-98') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-98', node2.lemma)
            
    # L -> 9_pl -> nom. sg. -> L
    feats[16] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_pl', 'nom. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. sg.', node2.lemma)
    
    # L->13_fp->L
    feats[17] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_fp', node2.lemma)
            
    # L->6_sg->L
    feats[18] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_sg', node2.lemma)
            
    # L -> -123 -> 10_sp -> L
    feats[19] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-123') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-123', '10_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_sp', node2.lemma)
    
    # L -> acc -> loc. sg. -> L
    feats[20] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc', 'loc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc. sg.', node2.lemma)
    
    # L -> 80 -> -92 -> L
    feats[21] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '80') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '80', '-92') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-92', node2.lemma)
    
    # L -> 9_pl -> 13_fp -> L
    feats[22] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_pl', '13_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_fp', node2.lemma)
    
    # L -> acc -> 42 -> L
    feats[23] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc', '42') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '42', node2.lemma)
    
    # L -> acc -> 6_sg -> L
    feats[24] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc', '6_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_sg', node2.lemma)
    
    # L -> 9_pl -> pl -> L
    feats[25] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_pl', 'pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    
    # L -> acc -> 11_pl -> L
    feats[26] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc', '11_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_pl', node2.lemma)
    
    # L -> acc -> 153 -> L
    feats[27] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc', '153') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    
    # L -> acc -> du_tp -> L
    feats[28] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc', 'du_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_tp', node2.lemma)
    
    # L -> 9_pl -> -36 -> L
    feats[29] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_pl', '-36') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-36', node2.lemma)
    
    # L -> 80 -> 29_sg -> L
    feats[30] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '80') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '80', '29_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '29_sg', node2.lemma)
    
    # L -> 80 -> -61 -> L
    feats[31] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '80') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '80', '-61') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-61', node2.lemma)
    
    # L -> 9_pl -> 11_du -> L
    feats[32] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_pl', '11_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_du', node2.lemma)
    
    # L -> 9_pl -> 3_fp -> L
    feats[33] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_pl', '3_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_fp', node2.lemma)
    
    # L -> 9_pl -> nom. du. -> L
    feats[34] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_pl', 'nom. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. du.', node2.lemma)
    
    # L -> -123 -> -119 -> L
    feats[35] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-123') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-123', '-119') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-119', node2.lemma)
    
    # L -> acc -> 3_fp -> L
    feats[36] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc', '3_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_fp', node2.lemma)
    
    # L -> acc -> nom. du. -> L
    feats[37] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc', 'nom. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. du.', node2.lemma)
    
    # L -> -141 -> 11_fp -> L
    feats[38] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-141') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-141', '11_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_fp', node2.lemma)
    
    # L -> 148 -> 5_pl -> L
    feats[39] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '148') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '148', '5_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_pl', node2.lemma)
    
    # L -> -33 -> -153 -> L
    feats[40] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-153') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-153', node2.lemma)
    
    # L -> -33 -> -81 -> L
    feats[41] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-81') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-81', node2.lemma)
    
    # L -> 9_pl -> 10_sp -> L
    feats[42] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_pl', '10_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_sp', node2.lemma)
    
    # L -> -33 -> 8_fp -> L
    feats[43] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '8_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    
    # L -> -141 -> adj -> L
    feats[44] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-141') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-141', 'adj') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'adj', node2.lemma)
    
    # L -> -33 -> 176 -> L
    feats[45] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '176') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '176', node2.lemma)
    
    # L -> -33 -> -62 -> L
    feats[46] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-62') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-62', node2.lemma)
    
    # L -> -33 -> -98 -> L
    feats[47] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-98') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-98', node2.lemma)
    
    # L -> 148 -> -50 -> L
    feats[48] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '148') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '148', '-50') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-50', node2.lemma)
    
    # L -> 148 -> -157 -> L
    feats[49] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '148') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '148', '-157') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-157', node2.lemma)
    
    # L -> -33 -> instr. sg. -> L
    feats[50] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', 'instr. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. sg.', node2.lemma)
    
    # L -> -33 -> 8_sp -> L
    feats[51] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '8_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_sp', node2.lemma)
    
    # L -> -33 -> 8_sg -> L
    feats[52] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '8_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_sg', node2.lemma)
    
    # L -> -33 -> 9_pl -> L
    feats[53] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '9_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_pl', node2.lemma)
    
    # L -> -33 -> -79 -> L
    feats[54] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-79') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-79', node2.lemma)
    
    # L -> -33 -> -50 -> L
    feats[55] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-50') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-50', node2.lemma)
    
    # L -> -33 -> -157 -> L
    feats[56] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-157') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-157', node2.lemma)
    
    # L -> -33 -> -297 -> L
    feats[57] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-297') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-297', node2.lemma)
    
    # L -> -33 -> 141 -> L
    feats[58] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '141') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '141', node2.lemma)
    
    # L -> -33 -> 79 -> L
    feats[59] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '79') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '79', node2.lemma)
    
    # L -> -33 -> 114 -> L
    feats[60] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '114') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '114', node2.lemma)
    
    # L -> -33 -> 159 -> L
    feats[61] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '159') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '159', node2.lemma)
    
    # L -> -33 -> -56 -> L
    feats[62] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-56') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-56', node2.lemma)
    
    # L -> -33 -> -18 -> L
    feats[63] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-18') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-18', node2.lemma)
    
    # L -> -33 -> -10 -> L
    feats[64] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-10') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-10', node2.lemma)
    
    # L -> 9_pl -> nom. fem -> L
    feats[65] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_pl', 'nom. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. fem', node2.lemma)
    
    # L -> -33 -> -41 -> L
    feats[66] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-41') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-41', node2.lemma)
    
    # L -> -33 -> -147 -> L
    feats[67] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-147') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-147', node2.lemma)
    
    # L -> -33 -> gen. sg. -> L
    feats[68] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', 'gen. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen. sg.', node2.lemma)
    
    # L -> -33 -> -166 -> L
    feats[69] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-166') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-166', node2.lemma)
    
    # L -> -33 -> 15_tp -> L
    feats[70] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '15_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_tp', node2.lemma)
    
    # L -> -33 -> 142 -> L
    feats[71] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '142') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '142', node2.lemma)
    
    # L -> -33 -> 11_du -> L
    feats[72] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '11_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_du', node2.lemma)
    
    # L -> -141 -> -303 -> L
    feats[73] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-141') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-141', '-303') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-303', node2.lemma)
    
    # L -> -33 -> gen -> L
    feats[74] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', 'gen') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen', node2.lemma)
    
    # L -> -33 -> -102 -> L
    feats[75] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-102') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-102', node2.lemma)
    
    # L -> -33 -> 171 -> L
    feats[76] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '171') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '171', node2.lemma)
    
    # L -> -33 -> instr. fem -> L
    feats[77] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', 'instr. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. fem', node2.lemma)
    
    # L -> -33 -> -27 -> L
    feats[78] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-27') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-27', node2.lemma)
    
    # L -> 148 -> 182 -> L
    feats[79] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '148') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '148', '182') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '182', node2.lemma)
    
    # L -> -33 -> 9_tp -> L
    feats[80] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '9_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_tp', node2.lemma)
    
    # L -> -33 -> -241 -> L
    feats[81] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-241') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-241', node2.lemma)
    
    # L -> -33 -> 6_du -> L
    feats[82] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '6_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_du', node2.lemma)
    
    # L -> 11_fp -> -104 -> L
    feats[83] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_fp', '-104') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-104', node2.lemma)
    
    # L -> -33 -> 173 -> L
    feats[84] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '173') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '173', node2.lemma)
    
    # L -> -33 -> nom. masc. -> L
    feats[85] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', 'nom. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. masc.', node2.lemma)
    
    # L -> -33 -> 2_sp -> L
    feats[86] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '2_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_sp', node2.lemma)
    
    # L -> -33 -> sg_fp -> L
    feats[87] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', 'sg_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_fp', node2.lemma)
    
    # L -> -33 -> 182 -> L
    feats[88] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '182') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '182', node2.lemma)
    
    # L -> -33 -> -59 -> L
    feats[89] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-59') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-59', node2.lemma)
    
    # L -> 56 -> 30_tp -> L
    feats[90] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '56') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '56', '30_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_tp', node2.lemma)
    
    # L -> -33 -> -299 -> L
    feats[91] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-299') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-299', node2.lemma)
    
    # L -> -33 -> 180 -> L
    feats[92] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '180') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '180', node2.lemma)
    
    # L -> -33 -> 14_sp -> L
    feats[93] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '14_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_sp', node2.lemma)
    
    # L -> -141 -> -83 -> L
    feats[94] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-141') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-141', '-83') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-83', node2.lemma)
    
    # L -> 56 -> -43 -> L
    feats[95] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '56') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '56', '-43') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-43', node2.lemma)
    
    # L -> 56 -> 81 -> L
    feats[96] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '56') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '56', '81') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '81', node2.lemma)
    
    # L -> 148 -> 36 -> L
    feats[97] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '148') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '148', '36') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '36', node2.lemma)
    
    # L -> 148 -> 14_sg -> L
    feats[98] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '148') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '148', '14_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_sg', node2.lemma)
    
    # L -> -141 -> 7_sg -> L
    feats[99] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-141') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-141', '7_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_sg', node2.lemma)
    
    # L -> -141 -> -44 -> L
    feats[100] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-141') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-141', '-44') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-44', node2.lemma)
    
    # L -> 11_fp -> 98 -> L
    feats[101] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_fp', '98') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '98', node2.lemma)
    
    # L -> -141 -> -132 -> L
    feats[102] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-141') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-141', '-132') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-132', node2.lemma)
    
    # L -> 11_fp -> 8_sp -> L
    feats[103] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_fp', '8_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_sp', node2.lemma)
    
    # L -> 56 -> 5_tp -> L
    feats[104] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '56') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '56', '5_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_tp', node2.lemma)
    
    # L -> -141 -> pl_fp -> L
    feats[105] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-141') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-141', 'pl_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_fp', node2.lemma)
    
    # L -> 56 -> -150 -> L
    feats[106] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '56') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '56', '-150') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-150', node2.lemma)
    
    # L -> -141 -> masc -> L
    feats[107] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-141') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-141', 'masc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'masc', node2.lemma)
    
    # L -> -141 -> 175 -> L
    feats[108] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-141') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-141', '175') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '175', node2.lemma)
    
    # L -> 56 -> 16_fp -> L
    feats[109] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '56') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '56', '16_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_fp', node2.lemma)
    
    # L -> -141 -> du_sp -> L
    feats[110] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-141') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-141', 'du_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_sp', node2.lemma)
    
    # L -> -141 -> -271 -> L
    feats[111] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-141') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-141', '-271') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-271', node2.lemma)
    
    # L -> -141 -> 2_pl -> L
    feats[112] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-141') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-141', '2_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_pl', node2.lemma)
    
    # L -> -141 -> 139 -> L
    feats[113] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-141') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-141', '139') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '139', node2.lemma)
    
    # L -> -141 -> 132 -> L
    feats[114] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-141') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-141', '132') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '132', node2.lemma)
    
    # L -> -79 -> -123 -> L
    feats[115] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-79', '-123') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-123', node2.lemma)
    
    # L -> -141 -> 11_sg -> L
    feats[116] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-141') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-141', '11_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_sg', node2.lemma)
    
    # L -> 56 -> pl -> L
    feats[117] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '56') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '56', 'pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    
    # L -> 56 -> 160 -> L
    feats[118] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '56') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '56', '160') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    
    # L -> -141 -> 4_sp -> L
    feats[119] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-141') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-141', '4_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_sp', node2.lemma)
    
    # L -> -141 -> 4_sg -> L
    feats[120] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-141') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-141', '4_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_sg', node2.lemma)
    
    # L -> -141 -> 97 -> L
    feats[121] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-141') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-141', '97') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '97', node2.lemma)
    
    # L -> -141 -> 10_tp -> L
    feats[122] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-141') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-141', '10_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_tp', node2.lemma)
    
    # L -> 56 -> -10 -> L
    feats[123] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '56') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '56', '-10') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-10', node2.lemma)
    
    # L -> 56 -> du_tp -> L
    feats[124] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '56') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '56', 'du_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_tp', node2.lemma)
    
    # L -> -141 -> 168 -> L
    feats[125] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-141') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-141', '168') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '168', node2.lemma)
    
    # L -> -141 -> -119 -> L
    feats[126] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-141') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-141', '-119') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-119', node2.lemma)
    
    # L -> 56 -> 128 -> L
    feats[127] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '56') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '56', '128') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '128', node2.lemma)
    
    # L -> 56 -> gen. sg. -> L
    feats[128] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '56') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '56', 'gen. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen. sg.', node2.lemma)
    
    # L -> -90 -> tp -> L
    feats[129] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-90', 'tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'tp', node2.lemma)
    
    # L -> -90 -> 179 -> L
    feats[130] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-90', '179') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '179', node2.lemma)
    
    # L -> 9_sp -> 10_fp -> L
    feats[131] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '10_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_fp', node2.lemma)
    
    # L -> -90 -> 61 -> L
    feats[132] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-90', '61') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '61', node2.lemma)
    
    # L->-113->C
    feats[133] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-113') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-113', node2.cng)
            
    # L -> 11_fp -> voc -> L
    feats[134] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_fp', 'voc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc', node2.lemma)
    
    # L -> -90 -> 27_du -> L
    feats[135] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-90', '27_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_du', node2.lemma)
    
    # L -> -90 -> 6_tp -> L
    feats[136] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-90', '6_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_tp', node2.lemma)
    
    # L -> -90 -> acc. neutr. -> L
    feats[137] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-90', 'acc. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. neutr.', node2.lemma)
    
    # L -> -90 -> 134 -> L
    feats[138] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-90', '134') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '134', node2.lemma)
    
    # L -> -90 -> voc. fem -> L
    feats[139] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-90', 'voc. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. fem', node2.lemma)
    
    # L->5_du->C
    feats[140] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_du', node2.cng)
            
    # L -> -90 -> -104 -> L
    feats[141] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-90', '-104') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-104', node2.lemma)
    
    # L->instr. neutr.->C
    feats[142] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. neutr.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. neutr.', node2.cng)
            
    # L -> -90 -> -84 -> L
    feats[143] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-90', '-84') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-84', node2.lemma)
    
    # L -> -90 -> -139 -> L
    feats[144] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-90', '-139') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-139', node2.lemma)
    
    # L -> -90 -> -159 -> L
    feats[145] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-90', '-159') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-159', node2.lemma)
    
    # L -> -90 -> 13_tp -> L
    feats[146] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-90', '13_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_tp', node2.lemma)
    
    # L -> -90 -> 115 -> L
    feats[147] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-90', '115') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '115', node2.lemma)
    
    # L -> -90 -> -57 -> L
    feats[148] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-90', '-57') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-57', node2.lemma)
    
    # L -> -90 -> -43 -> L
    feats[149] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-90', '-43') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-43', node2.lemma)
    
    # L -> -90 -> -71 -> L
    feats[150] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-90', '-71') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-71', node2.lemma)
    
    # L -> -90 -> -262 -> L
    feats[151] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-90', '-262') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-262', node2.lemma)
    
    # L -> -90 -> 59 -> L
    feats[152] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-90', '59') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '59', node2.lemma)
    
    # L -> 9_sp -> fp -> L
    feats[153] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', 'fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'fp', node2.lemma)
    
    # L -> 9_sp -> 176 -> L
    feats[154] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '176') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '176', node2.lemma)
    
    # L -> 11_fp -> 13_pl -> L
    feats[155] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_fp', '13_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_pl', node2.lemma)
    
    # L -> voc. sg. -> sg_sp -> L
    feats[156] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. sg.', 'sg_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_sp', node2.lemma)
    
    # L -> 56 -> 73 -> L
    feats[157] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '56') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '56', '73') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '73', node2.lemma)
    
    # L -> -79 -> 15_fp -> L
    feats[158] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-79', '15_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_fp', node2.lemma)
    
    # L -> voc. sg. -> 27_pl -> L
    feats[159] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. sg.', '27_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_pl', node2.lemma)
    
    # L -> voc. sg. -> 153 -> L
    feats[160] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. sg.', '153') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    
    # L->91->C
    feats[161] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '91') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '91', node2.cng)
            
    # L -> 56 -> -131 -> L
    feats[162] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '56') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '56', '-131') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-131', node2.lemma)
    
    # L -> 9_sp -> -102 -> L
    feats[163] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '-102') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-102', node2.lemma)
    
    # L -> 11_fp -> -29 -> L
    feats[164] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_fp', '-29') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-29', node2.lemma)
    
    # L->-190->T
    feats[165] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-190') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-190', node2.tup)
            
    # L -> -79 -> -29 -> L
    feats[166] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-79', '-29') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-29', node2.lemma)
    
    # L->157->T
    feats[167] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '157') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '157', node2.tup)
            
    # L->-37->T
    feats[168] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-37') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-37', node2.tup)
            
    # L -> -18 -> -109 -> L
    feats[169] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '-109') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-109', node2.lemma)
    
    # L->95->T
    feats[170] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '95') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '95', node2.tup)
            
    # L -> -18 -> 134 -> L
    feats[171] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '134') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '134', node2.lemma)
    
    # L->3_tp->T
    feats[172] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_tp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '3_tp', node2.tup)
            
    # L -> voc. sg. -> sp -> L
    feats[173] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. sg.', 'sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sp', node2.lemma)
    
    # L->81->T
    feats[174] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '81') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '81', node2.tup)
            
    # L -> 9_sp -> 28_tp -> L
    feats[175] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '28_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28_tp', node2.lemma)
    
    # L->-246->T
    feats[176] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-246') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-246', node2.tup)
            
    # L -> 157 -> 10_pl -> L
    feats[177] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '157', '10_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_pl', node2.lemma)
    
    # L -> 157 -> loc -> L
    feats[178] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '157', 'loc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc', node2.lemma)
    
    # L -> -18 -> -152 -> L
    feats[179] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '-152') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-152', node2.lemma)
    
    # L -> -18 -> 3_tp -> L
    feats[180] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '3_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_tp', node2.lemma)
    
    # L -> -90 -> -296 -> L
    feats[181] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-90', '-296') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-296', node2.lemma)
    
    # L->15_du->T
    feats[182] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_du') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '15_du', node2.tup)
            
    # L->-28->T
    feats[183] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-28') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-28', node2.tup)
            
    # L->8_sp->T
    feats[184] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_sp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '8_sp', node2.tup)
            
    # L->-169->T
    feats[185] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-169') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-169', node2.tup)
            
    # L -> -18 -> 30_sg -> L
    feats[186] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '30_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_sg', node2.lemma)
    
    # L -> 10_sg -> -35 -> L
    feats[187] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sg', '-35') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    
    # L -> -18 -> 15_du -> L
    feats[188] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '15_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_du', node2.lemma)
    
    # L->170->T
    feats[189] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '170') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '170', node2.tup)
            
    # L->5_fp->T
    feats[190] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '5_fp', node2.tup)
            
    # L -> 10_sg -> 42 -> L
    feats[191] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sg', '42') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '42', node2.lemma)
    
    # L -> 9_sp -> -240 -> L
    feats[192] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '-240') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-240', node2.lemma)
    
    # L -> 9_sp -> 11_sp -> L
    feats[193] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '11_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_sp', node2.lemma)
    
    # L -> voc. sg. -> 93 -> L
    feats[194] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. sg.', '93') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '93', node2.lemma)
    
    # L -> -18 -> 13_fp -> L
    feats[195] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '13_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_fp', node2.lemma)
    
    # L -> -18 -> -157 -> L
    feats[196] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '-157') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-157', node2.lemma)
    
    # L -> -18 -> abl. sg. -> L
    feats[197] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', 'abl. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. sg.', node2.lemma)
    
    # L -> -90 -> -64 -> L
    feats[198] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-90', '-64') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-64', node2.lemma)
    
    # L -> -90 -> 132 -> L
    feats[199] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-90', '132') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '132', node2.lemma)
    
    # L -> 9_sp -> -112 -> L
    feats[200] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '-112') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    
    # L -> 10_sg -> -147 -> L
    feats[201] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sg', '-147') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-147', node2.lemma)
    
    # L -> 10_fp -> 13_tp -> L
    feats[202] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_fp', '13_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_tp', node2.lemma)
    
    # L -> 10_fp -> -53 -> L
    feats[203] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_fp', '-53') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-53', node2.lemma)
    
    # L -> -18 -> du_fp -> L
    feats[204] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', 'du_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_fp', node2.lemma)
    
    # L -> 157 -> 1 -> L
    feats[205] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '157', '1') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '1', node2.lemma)
    
    # L -> 157 -> -20 -> L
    feats[206] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '157', '-20') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-20', node2.lemma)
    
    # L->99->T
    feats[207] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '99') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '99', node2.tup)
            
    # L -> 10_sg -> -49 -> L
    feats[208] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sg', '-49') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-49', node2.lemma)
    
    # L -> 114 -> 60 -> L
    feats[209] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '114') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '114', '60') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '60', node2.lemma)
    
    # L->-132->T
    feats[210] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-132') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-132', node2.tup)
            
    # L->voc. du.->T
    feats[211] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. du.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'voc. du.', node2.tup)
            
    # L -> -18 -> -309 -> L
    feats[212] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '-309') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-309', node2.lemma)
    
    # L -> -18 -> 8_du -> L
    feats[213] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '8_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_du', node2.lemma)
    
    # L->13_sp->T
    feats[214] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_sp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '13_sp', node2.tup)
            
    # L->11_sp->T
    feats[215] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_sp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '11_sp', node2.tup)
            
    # L -> 10_fp -> 37 -> L
    feats[216] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_fp', '37') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '37', node2.lemma)
    
    # L -> 157 -> -64 -> L
    feats[217] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '157', '-64') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-64', node2.lemma)
    
    # L -> 114 -> -19 -> L
    feats[218] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '114') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '114', '-19') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-19', node2.lemma)
    
    # L -> 10_sg -> -271 -> L
    feats[219] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sg', '-271') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-271', node2.lemma)
    
    # C->-158->L
    feats[220] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-158') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-158', node2.lemma)
            
    # C->35->L
    feats[221] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '35') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '35', node2.lemma)
            
    # C->92->L
    feats[222] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '92') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '92', node2.lemma)
            
    # C->5_sp->L
    feats[223] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '5_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sp', node2.lemma)
            
    # C->fem->L
    feats[224] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'fem', node2.lemma)
            
    # L -> 114 -> -137 -> L
    feats[225] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '114') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '114', '-137') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-137', node2.lemma)
    
    # L -> -156 -> -158 -> L
    feats[226] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-156', '-158') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-158', node2.lemma)
    
    # L -> 114 -> -303 -> L
    feats[227] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '114') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '114', '-303') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-303', node2.lemma)
    
    # L -> 10_fp -> -269 -> L
    feats[228] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_fp', '-269') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-269', node2.lemma)
    
    # L -> 114 -> loc. pl. -> L
    feats[229] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '114') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '114', 'loc. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc. pl.', node2.lemma)
    
    # L -> dat. pl. -> -46 -> L
    feats[230] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', '-46') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-46', node2.lemma)
    
    # L -> dat. pl. -> -53 -> L
    feats[231] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', '-53') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-53', node2.lemma)
    
    # L -> dat. pl. -> 5_pl -> L
    feats[232] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', '5_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_pl', node2.lemma)
    
    # L -> 114 -> 28_sg -> L
    feats[233] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '114') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '114', '28_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28_sg', node2.lemma)
    
    # L -> dat. pl. -> 115 -> L
    feats[234] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', '115') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '115', node2.lemma)
    
    # L -> dat. pl. -> 136 -> L
    feats[235] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', '136') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '136', node2.lemma)
    
    # L -> dat. pl. -> instr -> L
    feats[236] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', 'instr') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr', node2.lemma)
    
    # L -> 29 -> -149 -> L
    feats[237] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '29') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '29', '-149') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-149', node2.lemma)
    
    # L -> dat. pl. -> voc. masc. -> L
    feats[238] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', 'voc. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. masc.', node2.lemma)
    
    # L -> 114 -> 180 -> L
    feats[239] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '114') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '114', '180') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '180', node2.lemma)
    
    # L -> dat. pl. -> 14_fp -> L
    feats[240] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', '14_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_fp', node2.lemma)
    
    # L -> dat. pl. -> -153 -> L
    feats[241] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', '-153') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-153', node2.lemma)
    
    # L -> dat. pl. -> -43 -> L
    feats[242] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', '-43') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-43', node2.lemma)
    
    # L -> 13_sg -> -139 -> L
    feats[243] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_sg', '-139') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-139', node2.lemma)
    
    # L -> dat. pl. -> sg_sp -> L
    feats[244] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', 'sg_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_sp', node2.lemma)
    
    # L -> 114 -> -82 -> L
    feats[245] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '114') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '114', '-82') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-82', node2.lemma)
    
    # L -> dat. pl. -> 2 -> L
    feats[246] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', '2') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2', node2.lemma)
    
    # L -> dat. pl. -> 81 -> L
    feats[247] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', '81') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '81', node2.lemma)
    
    # L -> dat. pl. -> 59 -> L
    feats[248] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', '59') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '59', node2.lemma)
    
    # L -> dat. pl. -> -31 -> L
    feats[249] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', '-31') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-31', node2.lemma)
    
    # L -> dat. pl. -> instr. adj. -> L
    feats[250] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', 'instr. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. adj.', node2.lemma)
    
    # L -> dat. pl. -> dat -> L
    feats[251] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', 'dat') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat', node2.lemma)
    
    # L -> dat. pl. -> 3_sg -> L
    feats[252] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', '3_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_sg', node2.lemma)
    
    # L -> dat. pl. -> loc. sg. -> L
    feats[253] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', 'loc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc. sg.', node2.lemma)
    
    # L -> dat. pl. -> nom. sg. -> L
    feats[254] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', 'nom. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. sg.', node2.lemma)
    
    # L -> dat. pl. -> -22 -> L
    feats[255] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', '-22') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-22', node2.lemma)
    
    # L -> 114 -> -25 -> L
    feats[256] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '114') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '114', '-25') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-25', node2.lemma)
    
    # L -> dat. pl. -> 30_sg -> L
    feats[257] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', '30_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_sg', node2.lemma)
    
    # L -> 114 -> -21 -> L
    feats[258] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '114') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '114', '-21') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-21', node2.lemma)
    
    # L -> 114 -> voc. neutr. -> L
    feats[259] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '114') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '114', 'voc. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. neutr.', node2.lemma)
    
    # L -> dat. pl. -> fp -> L
    feats[260] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', 'fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'fp', node2.lemma)
    
    # L -> dat. pl. -> 5_tp -> L
    feats[261] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', '5_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_tp', node2.lemma)
    
    # L -> dat. pl. -> -35 -> L
    feats[262] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', '-35') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    
    # L -> dat. pl. -> 117 -> L
    feats[263] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', '117') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '117', node2.lemma)
    
    # L -> -156 -> 15_pl -> L
    feats[264] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-156', '15_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_pl', node2.lemma)
    
    # L -> -156 -> 59 -> L
    feats[265] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-156', '59') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '59', node2.lemma)
    
    # L -> dat. pl. -> -144 -> L
    feats[266] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', '-144') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-144', node2.lemma)
    
    # L -> dat. pl. -> -133 -> L
    feats[267] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', '-133') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-133', node2.lemma)
    
    # L -> dat. pl. -> -150 -> L
    feats[268] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', '-150') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-150', node2.lemma)
    
    # L -> dat. pl. -> adj -> L
    feats[269] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', 'adj') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'adj', node2.lemma)
    
    # L -> dat. pl. -> instr. sg. -> L
    feats[270] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', 'instr. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. sg.', node2.lemma)
    
    # L -> dat. pl. -> acc. adj. -> L
    feats[271] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', 'acc. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. adj.', node2.lemma)
    
    # L -> dat. pl. -> -28 -> L
    feats[272] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', '-28') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-28', node2.lemma)
    
    # L -> dat. pl. -> 27_sg -> L
    feats[273] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', '27_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_sg', node2.lemma)
    
    # L -> dat. pl. -> abl. pl. -> L
    feats[274] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', 'abl. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. pl.', node2.lemma)
    
    # L -> dat. pl. -> 31 -> L
    feats[275] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', '31') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '31', node2.lemma)
    
    # L -> dat. pl. -> 29_tp -> L
    feats[276] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', '29_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '29_tp', node2.lemma)
    
    # L -> 114 -> 121 -> L
    feats[277] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '114') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '114', '121') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '121', node2.lemma)
    
    # L -> dat. pl. -> 13_fp -> L
    feats[278] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', '13_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_fp', node2.lemma)
    
    # C->6_sg->L
    feats[279] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '6_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_sg', node2.lemma)
            
    # L -> dat. pl. -> -92 -> L
    feats[280] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', '-92') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-92', node2.lemma)
    
    # L -> dat. pl. -> 98 -> L
    feats[281] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', '98') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '98', node2.lemma)
    
    # L -> dat. pl. -> -169 -> L
    feats[282] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', '-169') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-169', node2.lemma)
    
    # L -> 114 -> -73 -> L
    feats[283] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '114') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '114', '-73') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-73', node2.lemma)
    
    # L -> dat. pl. -> -249 -> L
    feats[284] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', '-249') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-249', node2.lemma)
    
    # L -> dat. pl. -> -247 -> L
    feats[285] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', '-247') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-247', node2.lemma)
    
    # L -> dat. pl. -> 8_sg -> L
    feats[286] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', '8_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_sg', node2.lemma)
    
    # L -> dat. pl. -> 94 -> L
    feats[287] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', '94') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '94', node2.lemma)
    
    # L -> 114 -> -309 -> L
    feats[288] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '114') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '114', '-309') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-309', node2.lemma)
    
    # L -> dat. pl. -> -79 -> L
    feats[289] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', '-79') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-79', node2.lemma)
    
    # L -> dat. pl. -> -24 -> L
    feats[290] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', '-24') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    
    # L -> 114 -> -279 -> L
    feats[291] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '114') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '114', '-279') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-279', node2.lemma)
    
    # L -> dat. pl. -> -157 -> L
    feats[292] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', '-157') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-157', node2.lemma)
    
    # L -> dat. pl. -> abl. sg. -> L
    feats[293] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', 'abl. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. sg.', node2.lemma)
    
    # L -> 114 -> -44 -> L
    feats[294] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '114') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '114', '-44') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-44', node2.lemma)
    
    # L -> 29 -> -86 -> L
    feats[295] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '29') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '29', '-86') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-86', node2.lemma)
    
    # L -> 29 -> masc -> L
    feats[296] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '29') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '29', 'masc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'masc', node2.lemma)
    
    # L -> dat. pl. -> 5_fp -> L
    feats[297] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', '5_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_fp', node2.lemma)
    
    # L -> dat. pl. -> 160 -> L
    feats[298] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', '160') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    
    # L -> dat. pl. -> -19 -> L
    feats[299] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', '-19') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-19', node2.lemma)
    
    # L -> dat. pl. -> pl -> L
    feats[300] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', 'pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    
    # L -> dat. pl. -> neutr -> L
    feats[301] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', 'neutr') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'neutr', node2.lemma)
    
    # L -> dat. pl. -> 141 -> L
    feats[302] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', '141') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '141', node2.lemma)
    
    # L -> dat. pl. -> 170 -> L
    feats[303] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', '170') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '170', node2.lemma)
    
    # L -> dat. pl. -> -123 -> L
    feats[304] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', '-123') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-123', node2.lemma)
    
    # L -> 114 -> -11 -> L
    feats[305] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '114') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '114', '-11') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-11', node2.lemma)
    
    # L -> dat. pl. -> voc. sg. -> L
    feats[306] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', 'voc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. sg.', node2.lemma)
    
    # L -> dat. pl. -> 114 -> L
    feats[307] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', '114') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '114', node2.lemma)
    
    # L -> 13_sg -> -79 -> L
    feats[308] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_sg', '-79') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-79', node2.lemma)
    
    # L -> dat. pl. -> 172 -> L
    feats[309] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', '172') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '172', node2.lemma)
    
    # L -> dat. pl. -> 27_pl -> L
    feats[310] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', '27_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_pl', node2.lemma)
    
    # L -> dat. pl. -> 153 -> L
    feats[311] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', '153') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    
    # L -> dat. pl. -> 58 -> L
    feats[312] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', '58') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '58', node2.lemma)
    
    # C->118->L
    feats[313] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '118') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '118', node2.lemma)
            
    # C->instr. fem->L
    feats[314] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'instr. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. fem', node2.lemma)
            
    # L -> -190 -> -190 -> L
    feats[315] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-190') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-190', '-190') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-190', node2.lemma)
    
    # L -> -190 -> 38 -> L
    feats[316] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-190') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-190', '38') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '38', node2.lemma)
    
    # L -> 161 -> fp -> L
    feats[317] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '161') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '161', 'fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'fp', node2.lemma)
    
    # L -> 11_pl -> -31 -> L
    feats[318] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_pl', '-31') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-31', node2.lemma)
    
    # L -> dat. pl. -> voc. du. -> L
    feats[319] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', 'voc. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. du.', node2.lemma)
    
    # L -> 11_pl -> 27_sg -> L
    feats[320] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_pl', '27_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_sg', node2.lemma)
    
    # C->2_pl->L
    feats[321] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '2_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_pl', node2.lemma)
            
    # C->-131->L
    feats[322] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-131') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-131', node2.lemma)
            
    # L -> -156 -> 119 -> L
    feats[323] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-156', '119') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '119', node2.lemma)
    
    # C->4_fp->L
    feats[324] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '4_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_fp', node2.lemma)
            
    # C->-34->L
    feats[325] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-34') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-34', node2.lemma)
            
    # C->30_du->L
    feats[326] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '30_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_du', node2.lemma)
            
    # C->nom. fem->L
    feats[327] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'nom. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. fem', node2.lemma)
            
    # C->72->L
    feats[328] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '72') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '72', node2.lemma)
            
    # L -> 11_pl -> -123 -> L
    feats[329] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_pl', '-123') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-123', node2.lemma)
    
    # C->32->L
    feats[330] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '32') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '32', node2.lemma)
            
    # L -> 11_pl -> 11_pl -> L
    feats[331] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_pl', '11_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_pl', node2.lemma)
    
    # C->10_tp->L
    feats[332] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '10_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_tp', node2.lemma)
            
    # C->168->L
    feats[333] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '168') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '168', node2.lemma)
            
    # C->-119->L
    feats[334] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-119') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-119', node2.lemma)
            
    # C->-242->L
    feats[335] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-242') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-242', node2.lemma)
            
    # C->-89->L
    feats[336] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-89') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-89', node2.lemma)
            
    # L -> -190 -> abl. du. -> L
    feats[337] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-190') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-190', 'abl. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. du.', node2.lemma)
    
    # L -> -190 -> -91 -> L
    feats[338] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-190') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-190', '-91') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-91', node2.lemma)
    
    # C->158->L
    feats[339] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '158') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '158', node2.lemma)
            
    # C->-307->L
    feats[340] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-307') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-307', node2.lemma)
            
    # C->130->L
    feats[341] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '130') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '130', node2.lemma)
            
    # L -> 11_pl -> 135 -> L
    feats[342] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_pl', '135') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '135', node2.lemma)
    
    # C->-29->L
    feats[343] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-29') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-29', node2.lemma)
            
    # C->12_tp->L
    feats[344] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '12_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_tp', node2.lemma)
            
    # L -> 11_pl -> -291 -> L
    feats[345] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_pl', '-291') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-291', node2.lemma)
    
    # C->-158->C
    feats[346] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-158') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-158', node2.cng)
            
    # L -> 11_pl -> -36 -> L
    feats[347] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_pl', '-36') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-36', node2.lemma)
    
    # L -> 11_pl -> 7_tp -> L
    feats[348] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_pl', '7_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_tp', node2.lemma)
    
    # C->179->C
    feats[349] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '179') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '179', node2.cng)
            
    # C->-276->C
    feats[350] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-276') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-276', node2.cng)
            
    # C->-126->C
    feats[351] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-126') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-126', node2.cng)
            
    # C->9_sp->C
    feats[352] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', node2.cng)
            
    # L -> 13_sg -> -242 -> L
    feats[353] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_sg', '-242') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-242', node2.lemma)
    
    # C->5_sp->C
    feats[354] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '5_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', node2.cng)
            
    # C->fem->C
    feats[355] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'fem', node2.cng)
            
    # C->9_du->C
    feats[356] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', node2.cng)
            
    # C->7_du->C
    feats[357] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '7_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_du', node2.cng)
            
    # C->6_tp->C
    feats[358] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '6_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_tp', node2.cng)
            
    # C->134->C
    feats[359] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '134') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '134', node2.cng)
            
    # C->voc. fem->C
    feats[360] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'voc. fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. fem', node2.cng)
            
    # C->-273->C
    feats[361] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', node2.cng)
            
    # L -> 11_pl -> -61 -> L
    feats[362] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_pl', '-61') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-61', node2.lemma)
    
    # L -> 11_pl -> 110 -> L
    feats[363] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_pl', '110') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '110', node2.lemma)
    
    # C->30_tp->C
    feats[364] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '30_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_tp', node2.cng)
            
    # L -> -156 -> 155 -> L
    feats[365] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-156', '155') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '155', node2.lemma)
    
    # C->2_fp->C
    feats[366] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', node2.cng)
            
    # L -> 11_pl -> -303 -> L
    feats[367] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_pl', '-303') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-303', node2.lemma)
    
    # C->-37->C
    feats[368] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-37', node2.cng)
            
    # L -> 11_pl -> 111 -> L
    feats[369] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_pl', '111') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '111', node2.lemma)
    
    # C->-142->C
    feats[370] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-142') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-142', node2.cng)
            
    # C->-30->C
    feats[371] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', node2.cng)
            
    # C->102->C
    feats[372] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '102') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '102', node2.cng)
            
    # C->-302->C
    feats[373] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', node2.cng)
            
    # C->29->C
    feats[374] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '29') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '29', node2.cng)
            
    # C->161->C
    feats[375] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '161') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '161', node2.cng)
            
    # C->-139->C
    feats[376] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', node2.cng)
            
    # C->13_tp->C
    feats[377] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '13_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_tp', node2.cng)
            
    # C->-53->C
    feats[378] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-53') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-53', node2.cng)
            
    # C->5_pl->C
    feats[379] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '5_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_pl', node2.cng)
            
    # C->instr->C
    feats[380] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', node2.cng)
            
    # L -> 11_pl -> 28_sg -> L
    feats[381] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_pl', '28_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28_sg', node2.lemma)
    
    # C->-57->C
    feats[382] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-57') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-57', node2.cng)
            
    # C->77->C
    feats[383] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '77') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '77', node2.cng)
            
    # C->-153->C
    feats[384] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-153') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-153', node2.cng)
            
    # C->60->C
    feats[385] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '60') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '60', node2.cng)
            
    # C->81->C
    feats[386] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '81', node2.cng)
            
    # C->15_pl->C
    feats[387] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '15_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_pl', node2.cng)
            
    # C->59->C
    feats[388] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '59') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '59', node2.cng)
            
    # L -> 11_pl -> instr. masc. -> L
    feats[389] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_pl', 'instr. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. masc.', node2.lemma)
    
    # C->instr. adj.->C
    feats[390] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'instr. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. adj.', node2.cng)
            
    # C->loc->C
    feats[391] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'loc') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc', node2.cng)
            
    # C->-103->C
    feats[392] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-103') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-103', node2.cng)
            
    # C->nom. sg.->C
    feats[393] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'nom. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. sg.', node2.cng)
            
    # C->-81->C
    feats[394] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', node2.cng)
            
    # C->-22->C
    feats[395] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-22') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-22', node2.cng)
            
    # C->162->C
    feats[396] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '162') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '162', node2.cng)
            
    # L -> -190 -> -83 -> L
    feats[397] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-190') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-190', '-83') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-83', node2.lemma)
    
    # C->177->C
    feats[398] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '177') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '177', node2.cng)
            
    # C->fp->C
    feats[399] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'fp', node2.cng)
            
    # L -> -24 -> -262 -> L
    feats[400] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-24') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-24', '-262') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-262', node2.lemma)
    
    # L -> 96 -> 6_sg -> L
    feats[401] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '96') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '96', '6_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_sg', node2.lemma)
    
    # L -> 11_pl -> 36 -> L
    feats[402] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_pl', '36') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '36', node2.lemma)
    
    # L -> 96 -> 141 -> L
    feats[403] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '96') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '96', '141') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '141', node2.lemma)
    
    # L -> 96 -> 34 -> L
    feats[404] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '96') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '96', '34') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '34', node2.lemma)
    
    # L -> -190 -> 30_du -> L
    feats[405] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-190') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-190', '30_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_du', node2.lemma)
    
    # L -> -190 -> acc. pl. -> L
    feats[406] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-190') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-190', 'acc. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. pl.', node2.lemma)
    
    # L -> 96 -> -156 -> L
    feats[407] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '96') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '96', '-156') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-156', node2.lemma)
    
    # L -> 11_pl -> 2_sg -> L
    feats[408] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_pl', '2_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_sg', node2.lemma)
    
    # L -> 11_pl -> -99 -> L
    feats[409] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_pl', '-99') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-99', node2.lemma)
    
    # L -> 11_pl -> 132 -> L
    feats[410] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_pl', '132') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '132', node2.lemma)
    
    # L -> 96 -> -115 -> L
    feats[411] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '96') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '96', '-115') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-115', node2.lemma)
    
    # L -> 11_pl -> acc. pl. -> L
    feats[412] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_pl', 'acc. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. pl.', node2.lemma)
    
    # L -> -24 -> -101 -> L
    feats[413] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-24') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-24', '-101') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-101', node2.lemma)
    
    # L -> 161 -> 16_tp -> L
    feats[414] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '161') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '161', '16_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_tp', node2.lemma)
    
    # L -> 11_pl -> 51 -> L
    feats[415] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_pl', '51') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '51', node2.lemma)
    
    # L -> 11_pl -> 10_tp -> L
    feats[416] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_pl', '10_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_tp', node2.lemma)
    
    # L -> 11_pl -> 91 -> L
    feats[417] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_pl', '91') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    
    # L -> 96 -> 15_sp -> L
    feats[418] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '96') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '96', '15_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_sp', node2.lemma)
    
    # L -> 96 -> voc -> L
    feats[419] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '96') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '96', 'voc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc', node2.lemma)
    
    # L -> 172 -> dat. du. -> L
    feats[420] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '172') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '172', 'dat. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    
    # L -> 172 -> 10_fp -> L
    feats[421] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '172') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '172', '10_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_fp', node2.lemma)
    
    # L -> 96 -> 5_du -> L
    feats[422] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '96') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '96', '5_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_du', node2.lemma)
    
    # L -> 96 -> 116 -> L
    feats[423] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '96') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '96', '116') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '116', node2.lemma)
    
    # L -> 96 -> -93 -> L
    feats[424] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '96') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '96', '-93') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-93', node2.lemma)
    
    # L -> 172 -> 157 -> L
    feats[425] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '172') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '172', '157') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '157', node2.lemma)
    
    # C->150->C
    feats[426] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '150') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '150', node2.cng)
            
    # C->-261->C
    feats[427] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-261') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-261', node2.cng)
            
    # L -> 96 -> 1 -> L
    feats[428] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '96') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '96', '1') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '1', node2.lemma)
    
    # L -> 96 -> -20 -> L
    feats[429] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '96') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '96', '-20') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-20', node2.lemma)
    
    # L -> 38 -> 77 -> L
    feats[430] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '38') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '38', '77') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '77', node2.lemma)
    
    # L -> -24 -> -66 -> L
    feats[431] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-24') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-24', '-66') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-66', node2.lemma)
    
    # L -> -117 -> 3_du -> L
    feats[432] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-117') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-117', '3_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    
    # L -> -117 -> sg_sp -> L
    feats[433] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-117') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-117', 'sg_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_sp', node2.lemma)
    
    # L -> 38 -> 114 -> L
    feats[434] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '38') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '38', '114') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '114', node2.lemma)
    
    # L -> 38 -> 11_pl -> L
    feats[435] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '38') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '38', '11_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_pl', node2.lemma)
    
    # C->4_fp->C
    feats[436] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '4_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_fp', node2.cng)
            
    # L -> 96 -> -64 -> L
    feats[437] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '96') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '96', '-64') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-64', node2.lemma)
    
    # L -> -117 -> -41 -> L
    feats[438] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-117') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-117', '-41') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-41', node2.lemma)
    
    # L -> 172 -> instr. fem -> L
    feats[439] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '172') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '172', 'instr. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. fem', node2.lemma)
    
    # L -> 172 -> -68 -> L
    feats[440] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '172') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '172', '-68') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-68', node2.lemma)
    
    # C->-104->T
    feats[441] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-104') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-104', node2.tup)
            
    # C->48->T
    feats[442] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '48') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '48', node2.tup)
            
    # L -> 172 -> loc. pl. -> L
    feats[443] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '172') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '172', 'loc. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc. pl.', node2.lemma)
    
    # L -> 30_tp -> -46 -> L
    feats[444] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_tp', '-46') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-46', node2.lemma)
    
    # L -> 172 -> -87 -> L
    feats[445] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '172') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '172', '-87') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-87', node2.lemma)
    
    # L -> -117 -> -82 -> L
    feats[446] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-117') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-117', '-82') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-82', node2.lemma)
    
    # L -> -117 -> -163 -> L
    feats[447] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-117') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-117', '-163') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-163', node2.lemma)
    
    # L -> 172 -> -269 -> L
    feats[448] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '172') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '172', '-269') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-269', node2.lemma)
    
    # L -> 172 -> 82 -> L
    feats[449] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '172') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '172', '82') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '82', node2.lemma)
    
    # L -> 172 -> -25 -> L
    feats[450] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '172') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '172', '-25') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-25', node2.lemma)
    
    # L -> 172 -> 90 -> L
    feats[451] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '172') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '172', '90') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '90', node2.lemma)
    
    # L -> 172 -> 11_tp -> L
    feats[452] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '172') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '172', '11_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_tp', node2.lemma)
    
    # L -> 172 -> 13_pl -> L
    feats[453] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '172') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '172', '13_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_pl', node2.lemma)
    
    # L -> 172 -> -73 -> L
    feats[454] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '172') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '172', '-73') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-73', node2.lemma)
    
    # L -> 38 -> 30_pl -> L
    feats[455] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '38') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '38', '30_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_pl', node2.lemma)
    
    # L -> 172 -> 174 -> L
    feats[456] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '172') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '172', '174') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '174', node2.lemma)
    
    # C->79->T
    feats[457] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '79') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '79', node2.tup)
            
    # L -> 30_tp -> -19 -> L
    feats[458] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_tp', '-19') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-19', node2.lemma)
    
    # L -> 29_tp -> 88 -> L
    feats[459] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '29_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '29_tp', '88') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '88', node2.lemma)
    
    # L -> 172 -> -132 -> L
    feats[460] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '172') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '172', '-132') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-132', node2.lemma)
    
    # L -> 172 -> 16_sg -> L
    feats[461] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '172') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '172', '16_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_sg', node2.lemma)
    
    # L -> 172 -> 2_tp -> L
    feats[462] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '172') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '172', '2_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_tp', node2.lemma)
    
    # L -> 172 -> masc -> L
    feats[463] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '172') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '172', 'masc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'masc', node2.lemma)
    
    # L -> -117 -> -271 -> L
    feats[464] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-117') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-117', '-271') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-271', node2.lemma)
    
    # L -> 172 -> nom -> L
    feats[465] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '172') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '172', 'nom') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom', node2.lemma)
    
    # L -> 172 -> -240 -> L
    feats[466] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '172') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '172', '-240') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-240', node2.lemma)
    
    # L -> 172 -> -78 -> L
    feats[467] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '172') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '172', '-78') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-78', node2.lemma)
    
    # L -> 29_tp -> 4_sg -> L
    feats[468] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '29_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '29_tp', '4_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_sg', node2.lemma)
    
    # C->16_pl->T
    feats[469] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '16_pl') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '16_pl', node2.tup)
            
    # C->nom. masc.->T
    feats[470] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'nom. masc.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'nom. masc.', node2.tup)
            
    # L -> 30_tp -> -230 -> L
    feats[471] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_tp', '-230') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-230', node2.lemma)
    
    # L -> 27_du -> 80 -> L
    feats[472] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', '80') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '80', node2.lemma)
    
    # L -> -15 -> 2_fp -> L
    feats[473] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-15') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-15', '2_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_fp', node2.lemma)
    
    # L -> 27_pl -> 102 -> L
    feats[474] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_pl', '102') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '102', node2.lemma)
    
    # L -> 27_du -> adj -> L
    feats[475] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', 'adj') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'adj', node2.lemma)
    
    # L -> 27_du -> -263 -> L
    feats[476] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', '-263') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-263', node2.lemma)
    
    # C->88->T
    feats[477] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '88') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '88', node2.tup)
            
    # L -> 30_tp -> 30_pl -> L
    feats[478] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_tp', '30_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_pl', node2.lemma)
    
    # L -> 30_tp -> 119 -> L
    feats[479] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_tp', '119') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '119', node2.lemma)
    
    # L -> 27_pl -> -28 -> L
    feats[480] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_pl', '-28') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-28', node2.lemma)
    
    # L -> 30_tp -> 9_sg -> L
    feats[481] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_tp', '9_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_sg', node2.lemma)
    
    # L -> -15 -> voc. pl. -> L
    feats[482] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-15') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-15', 'voc. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. pl.', node2.lemma)
    
    # L -> 30_tp -> -271 -> L
    feats[483] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_tp', '-271') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-271', node2.lemma)
    
    # L -> 30_tp -> 4_fp -> L
    feats[484] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_tp', '4_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_fp', node2.lemma)
    
    # L -> 30_tp -> nom. fem -> L
    feats[485] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_tp', 'nom. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. fem', node2.lemma)
    
    # L -> 30_tp -> sg -> L
    feats[486] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_tp', 'sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg', node2.lemma)
    
    # L -> 30_tp -> 4_sp -> L
    feats[487] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_tp', '4_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_sp', node2.lemma)
    
    # L -> 30_tp -> 32 -> L
    feats[488] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_tp', '32') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '32', node2.lemma)
    
    # L -> 30_tp -> 137 -> L
    feats[489] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_tp', '137') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    
    # L -> 30_tp -> -210 -> L
    feats[490] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_tp', '-210') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-210', node2.lemma)
    
    # L -> 30_tp -> 33 -> L
    feats[491] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_tp', '33') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '33', node2.lemma)
    
    # L -> 42 -> 135 -> L
    feats[492] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '42') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '42', '135') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '135', node2.lemma)
    
    # L -> 30_tp -> -119 -> L
    feats[493] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_tp', '-119') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-119', node2.lemma)
    
    # T->-158->L
    feats[494] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-158') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-158', node2.lemma)
            
    # T->35->L
    feats[495] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '35') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '35', node2.lemma)
            
    # L -> 30_tp -> -29 -> L
    feats[496] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_tp', '-29') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-29', node2.lemma)
    
    # L -> 14_tp -> 35 -> L
    feats[497] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_tp', '35') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '35', node2.lemma)
    
    # L -> 27_du -> -303 -> L
    feats[498] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', '-303') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-303', node2.lemma)
    
    # L -> 14_tp -> 4_pl -> L
    feats[499] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_tp', '4_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_pl', node2.lemma)
    
    # L -> 14_tp -> 92 -> L
    feats[500] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_tp', '92') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '92', node2.lemma)
    
    # L -> -15 -> 69 -> L
    feats[501] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-15') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-15', '69') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '69', node2.lemma)
    
    # L -> -15 -> 70 -> L
    feats[502] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-15') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-15', '70') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '70', node2.lemma)
    
    # L -> 14_tp -> -273 -> L
    feats[503] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_tp', '-273') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-273', node2.lemma)
    
    # L -> 14_tp -> 2_fp -> L
    feats[504] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_tp', '2_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_fp', node2.lemma)
    
    # L -> 14_tp -> 12_pl -> L
    feats[505] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_tp', '12_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_pl', node2.lemma)
    
    # L -> 14_tp -> -30 -> L
    feats[506] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_tp', '-30') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-30', node2.lemma)
    
    # L -> 14_tp -> -143 -> L
    feats[507] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_tp', '-143') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-143', node2.lemma)
    
    # L -> 14_tp -> -104 -> L
    feats[508] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_tp', '-104') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-104', node2.lemma)
    
    # L -> 14_tp -> -141 -> L
    feats[509] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_tp', '-141') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-141', node2.lemma)
    
    # L -> 42 -> -154 -> L
    feats[510] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '42') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '42', '-154') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-154', node2.lemma)
    
    # L -> 14_tp -> -293 -> L
    feats[511] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_tp', '-293') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-293', node2.lemma)
    
    # L -> 14_tp -> 3_du -> L
    feats[512] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_tp', '3_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    
    # L -> 14_tp -> sg_sp -> L
    feats[513] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_tp', 'sg_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_sp', node2.lemma)
    
    # L -> 14_tp -> -262 -> L
    feats[514] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_tp', '-262') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-262', node2.lemma)
    
    # L -> 14_tp -> 15_pl -> L
    feats[515] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_tp', '15_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_pl', node2.lemma)
    
    # L -> 14_tp -> -31 -> L
    feats[516] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_tp', '-31') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-31', node2.lemma)
    
    # L -> -157 -> -76 -> L
    feats[517] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '-76') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-76', node2.lemma)
    
    # L -> 27_du -> 140 -> L
    feats[518] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', '140') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '140', node2.lemma)
    
    # L -> 14_tp -> -81 -> L
    feats[519] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_tp', '-81') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-81', node2.lemma)
    
    # L -> 14_tp -> 8_fp -> L
    feats[520] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_tp', '8_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    
    # L -> 27_pl -> 121 -> L
    feats[521] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_pl', '121') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '121', node2.lemma)
    
    # L -> -15 -> 50 -> L
    feats[522] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-15') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-15', '50') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '50', node2.lemma)
    
    # L -> -157 -> 79 -> L
    feats[523] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '79') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '79', node2.lemma)
    
    # T->-166->L
    feats[524] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-166') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-166', node2.lemma)
            
    # L -> 27_du -> 9_fp -> L
    feats[525] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', '9_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_fp', node2.lemma)
    
    # L -> -15 -> -307 -> L
    feats[526] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-15') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-15', '-307') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-307', node2.lemma)
    
    # L -> -15 -> 130 -> L
    feats[527] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-15') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-15', '130') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '130', node2.lemma)
    
    # L -> -157 -> 16_pl -> L
    feats[528] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '16_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_pl', node2.lemma)
    
    # L -> -109 -> -71 -> L
    feats[529] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-109') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-109', '-71') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-71', node2.lemma)
    
    # L -> -109 -> -260 -> L
    feats[530] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-109') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-109', '-260') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-260', node2.lemma)
    
    # L -> -157 -> 73 -> L
    feats[531] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '73') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '73', node2.lemma)
    
    # L -> 153 -> 117 -> L
    feats[532] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '153') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '153', '117') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '117', node2.lemma)
    
    # L -> 14_tp -> -86 -> L
    feats[533] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_tp', '-86') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-86', node2.lemma)
    
    # T->-99->L
    feats[534] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-99') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-99', node2.lemma)
            
    # T->-271->L
    feats[535] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-271') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-271', node2.lemma)
            
    # L -> 13_fp -> -157 -> L
    feats[536] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_fp', '-157') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-157', node2.lemma)
    
    # T->178->L
    feats[537] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '178') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '178', node2.lemma)
            
    # T->acc. sg.->L
    feats[538] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'acc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. sg.', node2.lemma)
            
    # T->4_sg->L
    feats[539] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '4_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_sg', node2.lemma)
            
    # T->91->L
    feats[540] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '91') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
            
    # L -> 13_fp -> 135 -> L
    feats[541] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_fp', '135') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '135', node2.lemma)
    
    # T->158->L
    feats[542] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '158') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '158', node2.lemma)
            
    # L -> -157 -> -306 -> L
    feats[543] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '-306') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-306', node2.lemma)
    
    # T->-276->C
    feats[544] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-276') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-276', node2.cng)
            
    # L -> 153 -> -12 -> L
    feats[545] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '153') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '153', '-12') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-12', node2.lemma)
    
    # T->148->C
    feats[546] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '148') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '148', node2.cng)
            
    # T->fem->C
    feats[547] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'fem', node2.cng)
            
    # L -> 2_fp -> 10_fp -> L
    feats[548] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', '10_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_fp', node2.lemma)
    
    # L -> 2_fp -> -190 -> L
    feats[549] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', '-190') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-190', node2.lemma)
    
    # L -> 13_fp -> instr. fem -> L
    feats[550] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_fp', 'instr. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. fem', node2.lemma)
    
    # T->7_sp->C
    feats[551] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '7_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_sp', node2.cng)
            
    # T->157->C
    feats[552] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '157', node2.cng)
            
    # T->2_fp->C
    feats[553] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', node2.cng)
            
    # T->-104->C
    feats[554] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', node2.cng)
            
    # T->29->C
    feats[555] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '29') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '29', node2.cng)
            
    # T->-84->C
    feats[556] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-84') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-84', node2.cng)
            
    # T->-159->C
    feats[557] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-159', node2.cng)
            
    # L -> -109 -> instr. masc. -> L
    feats[558] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-109') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-109', 'instr. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. masc.', node2.lemma)
    
    # L -> 13_fp -> -54 -> L
    feats[559] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_fp', '-54') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-54', node2.lemma)
    
    # L -> 13_fp -> -129 -> L
    feats[560] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_fp', '-129') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-129', node2.lemma)
    
    # T->-31->C
    feats[561] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-31') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-31', node2.cng)
            
    # L -> abl. sg. -> -159 -> L
    feats[562] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'abl. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'abl. sg.', '-159') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-159', node2.lemma)
    
    # L -> abl. sg. -> 13_tp -> L
    feats[563] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'abl. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'abl. sg.', '13_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_tp', node2.lemma)
    
    # T->-103->C
    feats[564] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-103') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-103', node2.cng)
            
    # T->-246->C
    feats[565] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-246') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-246', node2.cng)
            
    # T->loc. sg.->C
    feats[566] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'loc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', node2.cng)
            
    # T->-81->C
    feats[567] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', node2.cng)
            
    # T->acc. masc.->C
    feats[568] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'acc. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. masc.', node2.cng)
            
    # L -> 95 -> -73 -> L
    feats[569] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '95') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '95', '-73') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-73', node2.lemma)
    
    # T->-28->C
    feats[570] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-28') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-28', node2.cng)
            
    # T->27_sg->C
    feats[571] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', node2.cng)
            
    # T->abl. pl.->C
    feats[572] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'abl. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'abl. pl.', node2.cng)
            
    # L -> 153 -> 11_tp -> L
    feats[573] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '153') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '153', '11_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_tp', node2.lemma)
    
    # T->-92->C
    feats[574] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', node2.cng)
            
    # T->-247->C
    feats[575] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-247') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-247', node2.cng)
            
    # T->9_pl->C
    feats[576] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '9_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_pl', node2.cng)
            
    # T->-79->C
    feats[577] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-79', node2.cng)
            
    # L -> abl. sg. -> 79 -> L
    feats[578] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'abl. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'abl. sg.', '79') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '79', node2.lemma)
    
    # L -> 13_fp -> -32 -> L
    feats[579] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_fp', '-32') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-32', node2.lemma)
    
    # L -> 153 -> -34 -> L
    feats[580] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '153') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '153', '-34') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-34', node2.lemma)
    
    # L -> 2_fp -> -102 -> L
    feats[581] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', '-102') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-102', node2.lemma)
    
    # L -> 2_fp -> 171 -> L
    feats[582] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', '171') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '171', node2.lemma)
    
    # L -> 2_fp -> -241 -> L
    feats[583] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', '-241') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-241', node2.lemma)
    
    # L -> abl. sg. -> 9_tp -> L
    feats[584] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'abl. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'abl. sg.', '9_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_tp', node2.lemma)
    
    # L -> 2_fp -> 28_sg -> L
    feats[585] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', '28_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28_sg', node2.lemma)
    
    # L -> 6_sg -> 134 -> L
    feats[586] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', '134') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '134', node2.lemma)
    
    # L -> 6_sg -> -17 -> L
    feats[587] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', '-17') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-17', node2.lemma)
    
    # L -> -84 -> -302 -> L
    feats[588] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-84') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-84', '-302') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-302', node2.lemma)
    
    # L -> 2_fp -> sp -> L
    feats[589] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', 'sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sp', node2.lemma)
    
    # L -> 2_fp -> 40 -> L
    feats[590] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', '40') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '40', node2.lemma)
    
    # L -> 2_fp -> 76 -> L
    feats[591] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', '76') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '76', node2.lemma)
    
    # L -> 2_fp -> 6_fp -> L
    feats[592] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', '6_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_fp', node2.lemma)
    
    # L -> 2_fp -> 30_sp -> L
    feats[593] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', '30_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_sp', node2.lemma)
    
    # L -> 2_fp -> voc. neutr. -> L
    feats[594] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', 'voc. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. neutr.', node2.lemma)
    
    # L -> 2_fp -> 129 -> L
    feats[595] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', '129') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '129', node2.lemma)
    
    # L -> 2_fp -> 149 -> L
    feats[596] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', '149') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '149', node2.lemma)
    
    # L -> 58 -> -57 -> L
    feats[597] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '58', '-57') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-57', node2.lemma)
    
    # L -> 2_fp -> 120 -> L
    feats[598] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', '120') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '120', node2.lemma)
    
    # L -> 2_fp -> -16 -> L
    feats[599] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', '-16') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-16', node2.lemma)
    
    # L -> 5_sp -> 27_sg -> L
    feats[600] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', '27_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_sg', node2.lemma)
    
    # L -> 2_fp -> -132 -> L
    feats[601] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', '-132') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-132', node2.lemma)
    
    # L -> 2_fp -> -55 -> L
    feats[602] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', '-55') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-55', node2.lemma)
    
    # L -> 2_fp -> acc. du. -> L
    feats[603] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', 'acc. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. du.', node2.lemma)
    
    # L -> 2_fp -> gen. pl. -> L
    feats[604] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', 'gen. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen. pl.', node2.lemma)
    
    # L -> 6_sg -> 6_sg -> L
    feats[605] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', '6_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_sg', node2.lemma)
    
    # T->13_sp->C
    feats[606] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '13_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_sp', node2.cng)
            
    # L -> -84 -> -19 -> L
    feats[607] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-84') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-84', '-19') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-19', node2.lemma)
    
    # L -> 2_fp -> 11_sp -> L
    feats[608] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', '11_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_sp', node2.lemma)
    
    # L -> 2_fp -> -67 -> L
    feats[609] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', '-67') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-67', node2.lemma)
    
    # L -> 2_fp -> du_sp -> L
    feats[610] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', 'du_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_sp', node2.lemma)
    
    # L -> 2_fp -> -271 -> L
    feats[611] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', '-271') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-271', node2.lemma)
    
    # L -> 2_fp -> 132 -> L
    feats[612] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', '132') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '132', node2.lemma)
    
    # L -> 2_fp -> 11_sg -> L
    feats[613] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', '11_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_sg', node2.lemma)
    
    # L -> 58 -> -33 -> L
    feats[614] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '58', '-33') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-33', node2.lemma)
    
    # T->168->C
    feats[615] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '168') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '168', node2.cng)
            
    # T->8_pl->C
    feats[616] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '8_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_pl', node2.cng)
            
    # L -> 5_sp -> gen -> L
    feats[617] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', 'gen') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen', node2.lemma)
    
    # L -> 6_sg -> nom. du. -> L
    feats[618] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', 'nom. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. du.', node2.lemma)
    
    # L -> 5_sp -> 131 -> L
    feats[619] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', '131') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '131', node2.lemma)
    
    # L -> abl. sg. -> -114 -> L
    feats[620] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'abl. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'abl. sg.', '-114') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-114', node2.lemma)
    
    # L -> 5_sp -> -241 -> L
    feats[621] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', '-241') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-241', node2.lemma)
    
    # L -> -84 -> -230 -> L
    feats[622] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-84') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-84', '-230') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-230', node2.lemma)
    
    # L -> 5_sp -> -87 -> L
    feats[623] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', '-87') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-87', node2.lemma)
    
    # L -> -84 -> -241 -> L
    feats[624] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-84') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-84', '-241') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-241', node2.lemma)
    
    # L -> 5_sp -> -299 -> L
    feats[625] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', '-299') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-299', node2.lemma)
    
    # L -> 5_sp -> -82 -> L
    feats[626] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', '-82') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-82', node2.lemma)
    
    # L -> 5_sp -> 74 -> L
    feats[627] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', '74') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '74', node2.lemma)
    
    # L -> 5_sp -> 1 -> L
    feats[628] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', '1') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '1', node2.lemma)
    
    # L -> 5_sp -> sg_tp -> L
    feats[629] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', 'sg_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_tp', node2.lemma)
    
    # L -> 12_pl -> 136 -> L
    feats[630] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_pl', '136') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '136', node2.lemma)
    
    # L -> 5_sp -> 28_tp -> L
    feats[631] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', '28_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28_tp', node2.lemma)
    
    # L -> 5_sp -> -296 -> L
    feats[632] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', '-296') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-296', node2.lemma)
    
    # T->loc. sg.->T
    feats[633] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'loc. sg.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'loc. sg.', node2.tup)
            
    # L -> voc. pl. -> -57 -> L
    feats[634] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. pl.', '-57') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-57', node2.lemma)
    
    # L -> 5_sp -> 10_sp -> L
    feats[635] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', '10_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_sp', node2.lemma)
    
    # L -> 5_sp -> 5_sg -> L
    feats[636] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', '5_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sg', node2.lemma)
    
    # L -> 5_sp -> -73 -> L
    feats[637] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', '-73') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-73', node2.lemma)
    
    # L -> 5_sp -> 50 -> L
    feats[638] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', '50') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '50', node2.lemma)
    
    # L -> 5_sp -> 93 -> L
    feats[639] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', '93') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '93', node2.lemma)
    
    # L -> voc. pl. -> -35 -> L
    feats[640] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. pl.', '-35') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    
    # L -> 5_sp -> voc. du. -> L
    feats[641] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', 'voc. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. du.', node2.lemma)
    
    # L -> 5_sp -> 88 -> L
    feats[642] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', '88') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '88', node2.lemma)
    
    # L -> voc. pl. -> 117 -> L
    feats[643] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. pl.', '117') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '117', node2.lemma)
    
    # L -> 5_sp -> abl -> L
    feats[644] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', 'abl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl', node2.lemma)
    
    # L -> 12_pl -> 160 -> L
    feats[645] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_pl', '160') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    
    # L -> 12_pl -> -19 -> L
    feats[646] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_pl', '-19') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-19', node2.lemma)
    
    # L -> 5_sp -> du_sp -> L
    feats[647] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', 'du_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_sp', node2.lemma)
    
    # L -> 5_sp -> 2_sg -> L
    feats[648] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', '2_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_sg', node2.lemma)
    
    # L -> voc. pl. -> 9_pl -> L
    feats[649] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. pl.', '9_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_pl', node2.lemma)
    
    # L -> 5_sp -> 2_pl -> L
    feats[650] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', '2_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_pl', node2.lemma)
    
    # T->du_tp->T
    feats[651] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'du_tp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'du_tp', node2.tup)
            
    # L -> 5_sp -> 4_fp -> L
    feats[652] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', '4_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_fp', node2.lemma)
    
    # L -> 5_sp -> -34 -> L
    feats[653] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', '-34') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-34', node2.lemma)
    
    # L -> voc. pl. -> -123 -> L
    feats[654] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. pl.', '-123') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-123', node2.lemma)
    
    # L -> voc. pl. -> -33 -> L
    feats[655] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. pl.', '-33') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-33', node2.lemma)
    
    # L -> -84 -> 51 -> L
    feats[656] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-84') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-84', '51') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '51', node2.lemma)
    
    # L -> 5_sp -> 8_pl -> L
    feats[657] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', '8_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_pl', node2.lemma)
    
    # L -> 5_sp -> -119 -> L
    feats[658] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', '-119') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-119', node2.lemma)
    
    # L -> 5_sp -> 158 -> L
    feats[659] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', '158') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '158', node2.lemma)
    
    # L -> fem -> -158 -> L
    feats[660] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'fem', '-158') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-158', node2.lemma)
    
    # L -> fem -> 10_fp -> L
    feats[661] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'fem', '10_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_fp', node2.lemma)
    
    # L -> fem -> 5_sp -> L
    feats[662] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'fem', '5_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sp', node2.lemma)
    
    # L -> fem -> -17 -> L
    feats[663] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'fem', '-17') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-17', node2.lemma)
    
    # L -> fem -> 11_fp -> L
    feats[664] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'fem', '11_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_fp', node2.lemma)
    
    # L -> fem -> 12_pl -> L
    feats[665] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'fem', '12_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_pl', node2.lemma)
    
    # L -> fem -> -46 -> L
    feats[666] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'fem', '-46') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-46', node2.lemma)
    
    # L -> 12_pl -> 4_tp -> L
    feats[667] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_pl', '4_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_tp', node2.lemma)
    
    # L -> fem -> 161 -> L
    feats[668] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'fem', '161') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '161', node2.lemma)
    
    # L -> fem -> -84 -> L
    feats[669] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'fem', '-84') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-84', node2.lemma)
    
    # L -> fem -> 5_pl -> L
    feats[670] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'fem', '5_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_pl', node2.lemma)
    
    # L -> fem -> -152 -> L
    feats[671] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'fem', '-152') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-152', node2.lemma)
    
    # L -> fem -> -293 -> L
    feats[672] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'fem', '-293') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-293', node2.lemma)
    
    # T->10_du->T
    feats[673] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '10_du') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '10_du', node2.tup)
            
    # L -> fem -> -43 -> L
    feats[674] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'fem', '-43') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-43', node2.lemma)
    
    # L -> fem -> 3_du -> L
    feats[675] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'fem', '3_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    
    # L -> -139 -> 60 -> L
    feats[676] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '60') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '60', node2.lemma)
    
    # L -> fem -> dat -> L
    feats[677] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'fem', 'dat') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat', node2.lemma)
    
    # L -> fem -> 3_sg -> L
    feats[678] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'fem', '3_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_sg', node2.lemma)
    
    # L -> fem -> -246 -> L
    feats[679] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'fem', '-246') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-246', node2.lemma)
    
    # L -> fem -> -301 -> L
    feats[680] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'fem', '-301') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-301', node2.lemma)
    
    # L -> fem -> -263 -> L
    feats[681] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'fem', '-263') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-263', node2.lemma)
    
    # L -> fem -> instr. sg. -> L
    feats[682] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'fem', 'instr. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. sg.', node2.lemma)
    
    # L -> fem -> 42 -> L
    feats[683] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'fem', '42') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '42', node2.lemma)
    
    # L -> 12_pl -> -55 -> L
    feats[684] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_pl', '-55') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-55', node2.lemma)
    
    # L -> 12_pl -> acc. du. -> L
    feats[685] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_pl', 'acc. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. du.', node2.lemma)
    
    # L -> fem -> 8_sg -> L
    feats[686] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'fem', '8_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_sg', node2.lemma)
    
    # T->-67->T
    feats[687] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-67') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-67', node2.tup)
            
    # L -> -139 -> 159 -> L
    feats[688] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '159') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '159', node2.lemma)
    
    # L -> -72 -> 15_tp -> L
    feats[689] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-72') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-72', '15_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_tp', node2.lemma)
    
    # L -> voc. pl. -> -47 -> L
    feats[690] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. pl.', '-47') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-47', node2.lemma)
    
    # L -> fem -> 16_pl -> L
    feats[691] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'fem', '16_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_pl', node2.lemma)
    
    # L -> -92 -> 69 -> L
    feats[692] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '69') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '69', node2.lemma)
    
    # L -> -139 -> 76 -> L
    feats[693] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '76') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '76', node2.lemma)
    
    # L -> -297 -> -71 -> L
    feats[694] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-297') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-297', '-71') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-71', node2.lemma)
    
    # L -> 68 -> -51 -> L
    feats[695] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '68') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '68', '-51') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-51', node2.lemma)
    
    # L -> 68 -> 5_tp -> L
    feats[696] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '68') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '68', '5_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_tp', node2.lemma)
    
    # L -> fem -> 3_sp -> L
    feats[697] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'fem', '3_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_sp', node2.lemma)
    
    # L -> -92 -> -96 -> L
    feats[698] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '-96') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-96', node2.lemma)
    
    # L -> -139 -> 3_pl -> L
    feats[699] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '3_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_pl', node2.lemma)
    
    # L -> -139 -> 2_tp -> L
    feats[700] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '2_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_tp', node2.lemma)
    
    # L -> 68 -> -10 -> L
    feats[701] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '68') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '68', '-10') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-10', node2.lemma)
    
    # L -> 68 -> 89 -> L
    feats[702] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '68') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '68', '89') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '89', node2.lemma)
    
    # L -> 68 -> 7_pl -> L
    feats[703] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '68') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '68', '7_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_pl', node2.lemma)
    
    # L -> 68 -> 29_sg -> L
    feats[704] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '68') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '68', '29_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '29_sg', node2.lemma)
    
    # L -> 9_du -> dat. du. -> L
    feats[705] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', 'dat. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    
    # L -> -158 -> 131 -> L
    feats[706] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-158') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-158', '131') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '131', node2.lemma)
    
    # L -> -92 -> 158 -> L
    feats[707] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '158') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '158', node2.lemma)
    
    # L -> -72 -> -29 -> L
    feats[708] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-72') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-72', '-29') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-29', node2.lemma)
    
    # L -> 68 -> 182 -> L
    feats[709] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '68') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '68', '182') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '182', node2.lemma)
    
    # L -> 68 -> -292 -> L
    feats[710] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '68') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '68', '-292') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-292', node2.lemma)
    
    # L -> 9_du -> 77 -> L
    feats[711] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '77') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '77', node2.lemma)
    
    # L -> 68 -> 138 -> L
    feats[712] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '68') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '68', '138') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '138', node2.lemma)
    
    # L -> 68 -> 75 -> L
    feats[713] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '68') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '68', '75') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '75', node2.lemma)
    
    # L -> 68 -> -49 -> L
    feats[714] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '68') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '68', '-49') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-49', node2.lemma)
    
    # L -> 9_du -> -103 -> L
    feats[715] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '-103') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-103', node2.lemma)
    
    # L -> 9_du -> 3_sg -> L
    feats[716] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '3_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_sg', node2.lemma)
    
    # L -> -158 -> 5_sg -> L
    feats[717] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-158') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-158', '5_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sg', node2.lemma)
    
    # L -> -158 -> 13_pl -> L
    feats[718] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-158') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-158', '13_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_pl', node2.lemma)
    
    # L -> 181 -> 15_pl -> L
    feats[719] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '181') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '181', '15_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_pl', node2.lemma)
    
    # L -> 98 -> 15_sg -> L
    feats[720] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '98', '15_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_sg', node2.lemma)
    
    # L -> 68 -> 137 -> L
    feats[721] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '68') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '68', '137') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    
    # L -> 9_du -> 7_fp -> L
    feats[722] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '7_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_fp', node2.lemma)
    
    # L -> 98 -> 7_tp -> L
    feats[723] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '98', '7_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_tp', node2.lemma)
    
    # L -> -297 -> 158 -> L
    feats[724] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-297') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-297', '158') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '158', node2.lemma)
    
    # L -> 35 -> 80 -> L
    feats[725] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', '80') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '80', node2.lemma)
    
    # L -> 35 -> 157 -> L
    feats[726] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', '157') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '157', node2.lemma)
    
    # L -> 35 -> -37 -> L
    feats[727] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', '-37') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-37', node2.lemma)
    
    # L -> 35 -> -30 -> L
    feats[728] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', '-30') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-30', node2.lemma)
    
    # L -> 35 -> 48 -> L
    feats[729] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', '48') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '48', node2.lemma)
    
    # L -> 35 -> -159 -> L
    feats[730] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', '-159') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-159', node2.lemma)
    
    # L -> 35 -> -76 -> L
    feats[731] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', '-76') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-76', node2.lemma)
    
    # L -> -37 -> 3_tp -> L
    feats[732] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-37', '3_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_tp', node2.lemma)
    
    # L -> 9_du -> sg_tp -> L
    feats[733] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', 'sg_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_tp', node2.lemma)
    
    # L -> -37 -> -308 -> L
    feats[734] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-37', '-308') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-308', node2.lemma)
    
    # L -> 9_du -> 30_sp -> L
    feats[735] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '30_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_sp', node2.lemma)
    
    # L -> 35 -> -43 -> L
    feats[736] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', '-43') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-43', node2.lemma)
    
    # L -> 35 -> 81 -> L
    feats[737] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', '81') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '81', node2.lemma)
    
    # L -> 35 -> 59 -> L
    feats[738] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', '59') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '59', node2.lemma)
    
    # L -> 35 -> instr. adj. -> L
    feats[739] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', 'instr. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. adj.', node2.lemma)
    
    # L -> 35 -> 10_pl -> L
    feats[740] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', '10_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_pl', node2.lemma)
    
    # L -> 35 -> nom. sg. -> L
    feats[741] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', 'nom. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. sg.', node2.lemma)
    
    # L -> -159 -> -49 -> L
    feats[742] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-159', '-49') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-49', node2.lemma)
    
    # L -> 35 -> -22 -> L
    feats[743] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', '-22') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-22', node2.lemma)
    
    # L -> 35 -> 30_sg -> L
    feats[744] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', '30_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_sg', node2.lemma)
    
    # L -> 35 -> -301 -> L
    feats[745] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', '-301') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-301', node2.lemma)
    
    # L -> 35 -> 27_fp -> L
    feats[746] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', '27_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_fp', node2.lemma)
    
    # L -> 35 -> adj -> L
    feats[747] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', 'adj') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'adj', node2.lemma)
    
    # L -> 35 -> instr. sg. -> L
    feats[748] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', 'instr. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. sg.', node2.lemma)
    
    # L -> 35 -> 6_sg -> L
    feats[749] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', '6_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_sg', node2.lemma)
    
    # L -> 35 -> -249 -> L
    feats[750] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', '-249') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-249', node2.lemma)
    
    # L -> 35 -> 8_sg -> L
    feats[751] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', '8_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_sg', node2.lemma)
    
    # L -> 35 -> abl. sg. -> L
    feats[752] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', 'abl. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. sg.', node2.lemma)
    
    # L -> 35 -> neutr -> L
    feats[753] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', 'neutr') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'neutr', node2.lemma)
    
    # L -> 35 -> 141 -> L
    feats[754] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', '141') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '141', node2.lemma)
    
    # L -> 35 -> -123 -> L
    feats[755] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', '-123') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-123', node2.lemma)
    
    # L -> 35 -> 58 -> L
    feats[756] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', '58') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '58', node2.lemma)
    
    # L -> 35 -> 159 -> L
    feats[757] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', '159') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '159', node2.lemma)
    
    # L -> 35 -> 71 -> L
    feats[758] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', '71') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '71', node2.lemma)
    
    # L -> 9_du -> 156 -> L
    feats[759] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '156') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '156', node2.lemma)
    
    # L -> -159 -> 139 -> L
    feats[760] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-159', '139') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '139', node2.lemma)
    
    # L -> 35 -> 89 -> L
    feats[761] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', '89') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '89', node2.lemma)
    
    # L -> -159 -> -64 -> L
    feats[762] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-159', '-64') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-64', node2.lemma)
    
    # L -> 35 -> -291 -> L
    feats[763] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', '-291') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-291', node2.lemma)
    
    # L -> 35 -> 135 -> L
    feats[764] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', '135') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '135', node2.lemma)
    
    # L -> 35 -> -101 -> L
    feats[765] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', '-101') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-101', node2.lemma)
    
    # L -> 35 -> -121 -> L
    feats[766] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', '-121') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-121', node2.lemma)
    
    # L -> 35 -> -166 -> L
    feats[767] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', '-166') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-166', node2.lemma)
    
    # L -> 35 -> -115 -> L
    feats[768] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', '-115') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-115', node2.lemma)
    
    # L -> 35 -> 4_du -> L
    feats[769] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', '4_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_du', node2.lemma)
    
    # L -> 35 -> -91 -> L
    feats[770] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', '-91') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-91', node2.lemma)
    
    # L -> -37 -> nom. du. -> L
    feats[771] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-37', 'nom. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. du.', node2.lemma)
    
    # L -> -37 -> -61 -> L
    feats[772] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-37', '-61') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-61', node2.lemma)
    
    # L -> -37 -> -137 -> L
    feats[773] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-37', '-137') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-137', node2.lemma)
    
    # L -> -37 -> -161 -> L
    feats[774] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-37', '-161') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-161', node2.lemma)
    
    # L -> 98 -> 169 -> L
    feats[775] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '98', '169') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '169', node2.lemma)
    
    # L -> -37 -> nom. masc. -> L
    feats[776] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-37', 'nom. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. masc.', node2.lemma)
    
    # L -> -37 -> 12_sg -> L
    feats[777] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-37', '12_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_sg', node2.lemma)
    
    # L -> -37 -> 28 -> L
    feats[778] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-37', '28') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
    
    # L -> 7_du -> -46 -> L
    feats[779] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_du', '-46') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-46', node2.lemma)
    
    # L -> 7_du -> -30 -> L
    feats[780] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_du', '-30') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-30', node2.lemma)
    
    # L -> -37 -> -299 -> L
    feats[781] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-37', '-299') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-299', node2.lemma)
    
    # L -> -37 -> -154 -> L
    feats[782] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-37', '-154') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-154', node2.lemma)
    
    # L -> -37 -> 78 -> L
    feats[783] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-37', '78') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '78', node2.lemma)
    
    # L -> -37 -> 152 -> L
    feats[784] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-37', '152') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '152', node2.lemma)
    
    # L -> -37 -> sp -> L
    feats[785] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-37', 'sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sp', node2.lemma)
    
    # L -> -37 -> -93 -> L
    feats[786] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-37', '-93') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-93', node2.lemma)
    
    # L -> 7_du -> instr -> L
    feats[787] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_du', 'instr') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr', node2.lemma)
    
    # L -> -37 -> 1 -> L
    feats[788] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-37', '1') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '1', node2.lemma)
    
    # L -> 13_tp -> 3_tp -> L
    feats[789] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_tp', '3_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_tp', node2.lemma)
    
    # L -> -37 -> sg_tp -> L
    feats[790] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-37', 'sg_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_tp', node2.lemma)
    
    # L -> 7_du -> -43 -> L
    feats[791] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_du', '-43') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-43', node2.lemma)
    
    # L -> -37 -> -25 -> L
    feats[792] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-37', '-25') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-25', node2.lemma)
    
    # L -> -37 -> 90 -> L
    feats[793] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-37', '90') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '90', node2.lemma)
    
    # L -> 7_du -> 60 -> L
    feats[794] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_du', '60') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '60', node2.lemma)
    
    # L -> 35 -> -111 -> L
    feats[795] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', '-111') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-111', node2.lemma)
    
    # L -> 7_du -> 3_sg -> L
    feats[796] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_du', '3_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_sg', node2.lemma)
    
    # L -> -37 -> -296 -> L
    feats[797] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-37', '-296') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-296', node2.lemma)
    
    # L -> 37 -> instr -> L
    feats[798] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '37', 'instr') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr', node2.lemma)
    
    # L -> -37 -> 100 -> L
    feats[799] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-37', '100') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '100', node2.lemma)
    
    # L -> 8_sp -> -262 -> L
    feats[800] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_sp', '-262') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-262', node2.lemma)
    
    # L -> 7_du -> 176 -> L
    feats[801] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_du', '176') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '176', node2.lemma)
    
    # L -> -37 -> 5_sg -> L
    feats[802] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-37', '5_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sg', node2.lemma)
    
    # L -> 7_du -> -35 -> L
    feats[803] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_du', '-35') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    
    # L -> -37 -> gen. du. -> L
    feats[804] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-37', 'gen. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen. du.', node2.lemma)
    
    # L -> -37 -> -42 -> L
    feats[805] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-37', '-42') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-42', node2.lemma)
    
    # L -> -37 -> -73 -> L
    feats[806] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-37', '-73') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-73', node2.lemma)
    
    # L -> 7_du -> acc. masc. -> L
    feats[807] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_du', 'acc. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. masc.', node2.lemma)
    
    # L -> 7_du -> 117 -> L
    feats[808] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_du', '117') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '117', node2.lemma)
    
    # L -> 7_du -> -150 -> L
    feats[809] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_du', '-150') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-150', node2.lemma)
    
    # L -> 7_du -> adj -> L
    feats[810] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_du', 'adj') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'adj', node2.lemma)
    
    # L -> 7_du -> instr. sg. -> L
    feats[811] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_du', 'instr. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. sg.', node2.lemma)
    
    # L -> 7_du -> 13_fp -> L
    feats[812] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_du', '13_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_fp', node2.lemma)
    
    # L -> -37 -> 50 -> L
    feats[813] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-37', '50') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '50', node2.lemma)
    
    # L -> 7_du -> -169 -> L
    feats[814] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_du', '-169') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-169', node2.lemma)
    
    # L -> -37 -> gen. pl. -> L
    feats[815] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-37', 'gen. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen. pl.', node2.lemma)
    
    # L -> -37 -> voc. du. -> L
    feats[816] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-37', 'voc. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. du.', node2.lemma)
    
    # L -> -37 -> pl_fp -> L
    feats[817] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-37', 'pl_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_fp', node2.lemma)
    
    # L -> 7_du -> 13_sg -> L
    feats[818] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_du', '13_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_sg', node2.lemma)
    
    # L -> -37 -> 15_fp -> L
    feats[819] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-37', '15_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_fp', node2.lemma)
    
    # L -> 7_du -> abl. sg. -> L
    feats[820] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_du', 'abl. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. sg.', node2.lemma)
    
    # L -> -37 -> -14 -> L
    feats[821] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-37', '-14') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-14', node2.lemma)
    
    # L -> 7_du -> 160 -> L
    feats[822] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_du', '160') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    
    # L -> 7_du -> pl -> L
    feats[823] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_du', 'pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    
    # L -> -37 -> -26 -> L
    feats[824] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-37', '-26') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-26', node2.lemma)
    
    # L -> -37 -> nom -> L
    feats[825] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-37', 'nom') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom', node2.lemma)
    
    # L -> 13_tp -> 160 -> L
    feats[826] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_tp', '160') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    
    # L -> 13_tp -> -19 -> L
    feats[827] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_tp', '-19') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-19', node2.lemma)
    
    # L -> 7_du -> 172 -> L
    feats[828] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_du', '172') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '172', node2.lemma)
    
    # L -> -37 -> -32 -> L
    feats[829] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-37', '-32') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-32', node2.lemma)
    
    # L -> 7_du -> 159 -> L
    feats[830] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_du', '159') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '159', node2.lemma)
    
    # L -> 7_du -> -69 -> L
    feats[831] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_du', '-69') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-69', node2.lemma)
    
    # L -> 37 -> -297 -> L
    feats[832] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '37', '-297') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-297', node2.lemma)
    
    # L -> 7_du -> 56 -> L
    feats[833] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_du', '56') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '56', node2.lemma)
    
    # L -> 35 -> -13 -> L
    feats[834] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', '-13') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-13', node2.lemma)
    
    # L -> 7_du -> du_tp -> L
    feats[835] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_du', 'du_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_tp', node2.lemma)
    
    # L -> 7_du -> 135 -> L
    feats[836] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_du', '135') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '135', node2.lemma)
    
    # L -> 7_du -> -101 -> L
    feats[837] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_du', '-101') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-101', node2.lemma)
    
    # L -> 7_du -> nom. adj. -> L
    feats[838] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_du', 'nom. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. adj.', node2.lemma)
    
    # L -> 7_du -> -36 -> L
    feats[839] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_du', '-36') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-36', node2.lemma)
    
    # L -> 7_du -> -12 -> L
    feats[840] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_du', '-12') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-12', node2.lemma)
    
    # L -> 7_du -> 3_fp -> L
    feats[841] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_du', '3_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_fp', node2.lemma)
    
    # L -> 170 -> instr. pl. -> L
    feats[842] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '170') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '170', 'instr. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. pl.', node2.lemma)
    
    # L -> 7_du -> 29_sg -> L
    feats[843] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_du', '29_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '29_sg', node2.lemma)
    
    # L -> 7_du -> -61 -> L
    feats[844] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_du', '-61') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-61', node2.lemma)
    
    # L -> 7_du -> 110 -> L
    feats[845] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_du', '110') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '110', node2.lemma)
    
    # L -> 7_du -> -68 -> L
    feats[846] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_du', '-68') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-68', node2.lemma)
    
    # L -> 7_du -> -137 -> L
    feats[847] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_du', '-137') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-137', node2.lemma)
    
    # L -> 7_du -> 6_sp -> L
    feats[848] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_du', '6_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_sp', node2.lemma)
    
    # L -> 7_du -> 9_tp -> L
    feats[849] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_du', '9_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_tp', node2.lemma)
    
    # L -> 7_du -> 131 -> L
    feats[850] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_du', '131') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '131', node2.lemma)
    
    # L -> 7_du -> 15_sp -> L
    feats[851] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_du', '15_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_sp', node2.lemma)
    
    # L -> 7_du -> 69 -> L
    feats[852] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_du', '69') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '69', node2.lemma)
    
    # L -> 7_du -> -113 -> L
    feats[853] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_du', '-113') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-113', node2.lemma)
    
    # L -> 7_du -> -230 -> L
    feats[854] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_du', '-230') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-230', node2.lemma)
    
    # L -> 7_du -> -241 -> L
    feats[855] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_du', '-241') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-241', node2.lemma)
    
    # L -> 13_tp -> du_fp -> L
    feats[856] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_tp', 'du_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_fp', node2.lemma)
    
    # L -> 102 -> -142 -> L
    feats[857] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '102') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '102', '-142') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-142', node2.lemma)
    
    # L -> 13_tp -> 78 -> L
    feats[858] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_tp', '78') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '78', node2.lemma)
    
    # L -> 37 -> 28 -> L
    feats[859] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '37', '28') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
    
    # L -> dat. du. -> -57 -> L
    feats[860] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. du.', '-57') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-57', node2.lemma)
    
    # L -> 5_fp -> -84 -> L
    feats[861] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', '-84') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-84', node2.lemma)
    
    # L -> 7_du -> -96 -> L
    feats[862] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_du', '-96') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-96', node2.lemma)
    
    # L -> 102 -> 29_tp -> L
    feats[863] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '102') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '102', '29_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '29_tp', node2.lemma)
    
    # L -> 37 -> -309 -> L
    feats[864] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '37', '-309') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-309', node2.lemma)
    
    # L -> 37 -> 8_du -> L
    feats[865] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '37', '8_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_du', node2.lemma)
    
    # L -> dat. du. -> neutr -> L
    feats[866] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. du.', 'neutr') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'neutr', node2.lemma)
    
    # L -> dat. du. -> 141 -> L
    feats[867] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. du.', '141') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '141', node2.lemma)
    
    # L -> 8_sp -> -240 -> L
    feats[868] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_sp', '-240') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-240', node2.lemma)
    
    # L -> 8_sp -> -64 -> L
    feats[869] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_sp', '-64') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-64', node2.lemma)
    
    # L -> 7_du -> 158 -> L
    feats[870] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_du', '158') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '158', node2.lemma)
    
    # L -> 7_du -> -307 -> L
    feats[871] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_du', '-307') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-307', node2.lemma)
    
    # L -> -169 -> tp -> L
    feats[872] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-169') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-169', 'tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'tp', node2.lemma)
    
    # L -> dat. du. -> 14_sp -> L
    feats[873] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. du.', '14_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_sp', node2.lemma)
    
    # L -> -169 -> 6_tp -> L
    feats[874] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-169') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-169', '6_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_tp', node2.lemma)
    
    # L -> 5_fp -> 2_sp -> L
    feats[875] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', '2_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_sp', node2.lemma)
    
    # L -> 159 -> 11_fp -> L
    feats[876] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '11_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_fp', node2.lemma)
    
    # L -> 102 -> 5_sg -> L
    feats[877] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '102') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '102', '5_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sg', node2.lemma)
    
    # L -> 102 -> 13_pl -> L
    feats[878] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '102') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '102', '13_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_pl', node2.lemma)
    
    # L -> 6_tp -> -62 -> L
    feats[879] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_tp', '-62') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-62', node2.lemma)
    
    # L -> dat. du. -> acc. du. -> L
    feats[880] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. du.', 'acc. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. du.', node2.lemma)
    
    # L -> 159 -> 27_fp -> L
    feats[881] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '27_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_fp', node2.lemma)
    
    # L -> dat. du. -> nom -> L
    feats[882] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. du.', 'nom') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom', node2.lemma)
    
    # L -> -53 -> -297 -> L
    feats[883] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-53') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-53', '-297') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-297', node2.lemma)
    
    # L -> dat. du. -> -243 -> L
    feats[884] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. du.', '-243') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-243', node2.lemma)
    
    # L -> 5_fp -> abl -> L
    feats[885] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', 'abl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl', node2.lemma)
    
    # L -> dat. du. -> -99 -> L
    feats[886] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. du.', '-99') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-99', node2.lemma)
    
    # L -> dat. du. -> 6_pl -> L
    feats[887] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. du.', '6_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_pl', node2.lemma)
    
    # L -> -53 -> voc. sg. -> L
    feats[888] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-53') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-53', 'voc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. sg.', node2.lemma)
    
    # L -> -53 -> -33 -> L
    feats[889] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-53') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-53', '-33') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-33', node2.lemma)
    
    # L -> dat. du. -> 4_fp -> L
    feats[890] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. du.', '4_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_fp', node2.lemma)
    
    # L -> dat. du. -> 30_du -> L
    feats[891] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. du.', '30_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_du', node2.lemma)
    
    # L -> 6_tp -> 135 -> L
    feats[892] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_tp', '135') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '135', node2.lemma)
    
    # L -> 6_tp -> -41 -> L
    feats[893] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_tp', '-41') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-41', node2.lemma)
    
    # L -> dat. du. -> -112 -> L
    feats[894] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. du.', '-112') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    
    # L -> dat. du. -> -210 -> L
    feats[895] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. du.', '-210') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-210', node2.lemma)
    
    # L -> dat. du. -> 97 -> L
    feats[896] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. du.', '97') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '97', node2.lemma)
    
    # L -> dat. du. -> 33 -> L
    feats[897] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. du.', '33') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '33', node2.lemma)
    
    # L -> dat. du. -> 9_fp -> L
    feats[898] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. du.', '9_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_fp', node2.lemma)
    
    # L -> dat. du. -> 8_pl -> L
    feats[899] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. du.', '8_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_pl', node2.lemma)
    
    # L -> dat. du. -> -119 -> L
    feats[900] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. du.', '-119') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-119', node2.lemma)
    
    # L -> dat. du. -> -307 -> L
    feats[901] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. du.', '-307') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-307', node2.lemma)
    
    # L -> 102 -> -114 -> L
    feats[902] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '102') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '102', '-114') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-114', node2.lemma)
    
    # L -> dat. du. -> -29 -> L
    feats[903] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. du.', '-29') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-29', node2.lemma)
    
    # L -> dat. du. -> 155 -> L
    feats[904] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. du.', '155') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '155', node2.lemma)
    
    # L -> tp -> -126 -> L
    feats[905] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '-126') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-126', node2.lemma)
    
    # L -> tp -> -190 -> L
    feats[906] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '-190') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-190', node2.lemma)
    
    # L -> tp -> fem -> L
    feats[907] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', 'fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'fem', node2.lemma)
    
    # L -> tp -> 7_du -> L
    feats[908] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '7_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_du', node2.lemma)
    
    # L -> -53 -> 69 -> L
    feats[909] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-53') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-53', '69') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '69', node2.lemma)
    
    # L -> 5_fp -> 12_tp -> L
    feats[910] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', '12_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_tp', node2.lemma)
    
    # L -> 159 -> -61 -> L
    feats[911] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '-61') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-61', node2.lemma)
    
    # L -> tp -> 12_pl -> L
    feats[912] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '12_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_pl', node2.lemma)
    
    # L -> tp -> -37 -> L
    feats[913] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '-37') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-37', node2.lemma)
    
    # L -> tp -> 48 -> L
    feats[914] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '48') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '48', node2.lemma)
    
    # L -> tp -> -141 -> L
    feats[915] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '-141') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-141', node2.lemma)
    
    # L -> tp -> -15 -> L
    feats[916] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '-15') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-15', node2.lemma)
    
    # L -> tp -> -84 -> L
    feats[917] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '-84') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-84', node2.lemma)
    
    # L -> tp -> -139 -> L
    feats[918] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '-139') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-139', node2.lemma)
    
    # L -> tp -> 13_tp -> L
    feats[919] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '13_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_tp', node2.lemma)
    
    # L -> -53 -> 82 -> L
    feats[920] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-53') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-53', '82') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '82', node2.lemma)
    
    # L -> tp -> -71 -> L
    feats[921] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '-71') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-71', node2.lemma)
    
    # L -> tp -> -31 -> L
    feats[922] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '-31') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-31', node2.lemma)
    
    # L -> tp -> instr. adj. -> L
    feats[923] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', 'instr. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. adj.', node2.lemma)
    
    # L -> tp -> 10_pl -> L
    feats[924] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '10_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_pl', node2.lemma)
    
    # L -> -142 -> 60 -> L
    feats[925] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-142') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-142', '60') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '60', node2.lemma)
    
    # L -> tp -> dat -> L
    feats[926] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', 'dat') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat', node2.lemma)
    
    # L -> tp -> 3_sg -> L
    feats[927] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '3_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_sg', node2.lemma)
    
    # L -> tp -> loc. sg. -> L
    feats[928] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', 'loc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc. sg.', node2.lemma)
    
    # L -> tp -> 30_sg -> L
    feats[929] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '30_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_sg', node2.lemma)
    
    # L -> tp -> fp -> L
    feats[930] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', 'fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'fp', node2.lemma)
    
    # L -> tp -> 176 -> L
    feats[931] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '176') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '176', node2.lemma)
    
    # L -> 159 -> 28_tp -> L
    feats[932] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '28_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28_tp', node2.lemma)
    
    # L -> 159 -> -151 -> L
    feats[933] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '-151') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-151', node2.lemma)
    
    # L -> 159 -> 121 -> L
    feats[934] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '121') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '121', node2.lemma)
    
    # L -> 159 -> instr. du. -> L
    feats[935] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', 'instr. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. du.', node2.lemma)
    
    # L -> 159 -> 99 -> L
    feats[936] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '99') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '99', node2.lemma)
    
    # L -> -53 -> 88 -> L
    feats[937] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-53') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-53', '88') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '88', node2.lemma)
    
    # L -> -53 -> pl_fp -> L
    feats[938] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-53') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-53', 'pl_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_fp', node2.lemma)
    
    # L -> 159 -> 7_sg -> L
    feats[939] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '7_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_sg', node2.lemma)
    
    # L -> 159 -> 14_sg -> L
    feats[940] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '14_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_sg', node2.lemma)
    
    # L -> 159 -> -132 -> L
    feats[941] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '-132') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-132', node2.lemma)
    
    # L -> 6_tp -> 2_sg -> L
    feats[942] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_tp', '2_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_sg', node2.lemma)
    
    # L -> 159 -> 109 -> L
    feats[943] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '109') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '109', node2.lemma)
    
    # L -> 159 -> -86 -> L
    feats[944] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '-86') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-86', node2.lemma)
    
    # L -> 159 -> -14 -> L
    feats[945] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '-14') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-14', node2.lemma)
    
    # L -> 159 -> 9_sg -> L
    feats[946] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '9_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_sg', node2.lemma)
    
    # L -> 159 -> 112 -> L
    feats[947] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '112') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '112', node2.lemma)
    
    # L -> 159 -> pl_sp -> L
    feats[948] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', 'pl_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_sp', node2.lemma)
    
    # L -> 159 -> -240 -> L
    feats[949] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '-240') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-240', node2.lemma)
    
    # L -> 159 -> 13_sp -> L
    feats[950] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '13_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_sp', node2.lemma)
    
    # L -> 159 -> 175 -> L
    feats[951] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '175') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '175', node2.lemma)
    
    # L -> 159 -> 2_sg -> L
    feats[952] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '2_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_sg', node2.lemma)
    
    # L -> -142 -> -38 -> L
    feats[953] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-142') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-142', '-38') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-38', node2.lemma)
    
    # L -> 159 -> 6_pl -> L
    feats[954] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '6_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_pl', node2.lemma)
    
    # L -> 159 -> 30_du -> L
    feats[955] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '30_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_du', node2.lemma)
    
    # L -> 159 -> 72 -> L
    feats[956] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '72') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '72', node2.lemma)
    
    # L -> 159 -> sg -> L
    feats[957] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', 'sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg', node2.lemma)
    
    # L -> 159 -> 8_tp -> L
    feats[958] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '8_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_tp', node2.lemma)
    
    # L -> 159 -> acc. sg. -> L
    feats[959] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', 'acc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. sg.', node2.lemma)
    
    # L -> 159 -> 97 -> L
    feats[960] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '97') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '97', node2.lemma)
    
    # L -> 159 -> 91 -> L
    feats[961] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '91') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    
    # L -> 159 -> 168 -> L
    feats[962] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '168') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '168', node2.lemma)
    
    # L -> 159 -> 8_pl -> L
    feats[963] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '8_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_pl', node2.lemma)
    
    # L -> 159 -> 158 -> L
    feats[964] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '158') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '158', node2.lemma)
    
    # L -> 159 -> -307 -> L
    feats[965] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '-307') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-307', node2.lemma)
    
    # L -> -69 -> 35 -> L
    feats[966] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-69') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-69', '35') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '35', node2.lemma)
    
    # L -> -69 -> dat. du. -> L
    feats[967] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-69') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-69', 'dat. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    
    # L -> -69 -> 179 -> L
    feats[968] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-69') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-69', '179') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '179', node2.lemma)
    
    # L -> -69 -> -109 -> L
    feats[969] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-69') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-69', '-109') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-109', node2.lemma)
    
    # L -> 7_sp -> -117 -> L
    feats[970] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_sp', '-117') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-117', node2.lemma)
    
    # L -> -142 -> instr. masc. -> L
    feats[971] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-142') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-142', 'instr. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. masc.', node2.lemma)
    
    # L -> -142 -> -93 -> L
    feats[972] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-142') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-142', '-93') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-93', node2.lemma)
    
    # L -> -69 -> dat. pl. -> L
    feats[973] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-69') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-69', 'dat. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. pl.', node2.lemma)
    
    # L -> -69 -> 96 -> L
    feats[974] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-69') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-69', '96') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '96', node2.lemma)
    
    # L -> 7_sp -> -308 -> L
    feats[975] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_sp', '-308') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-308', node2.lemma)
    
    # L -> -69 -> 30_tp -> L
    feats[976] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-69') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-69', '30_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_tp', node2.lemma)
    
    # L -> -69 -> 68 -> L
    feats[977] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-69') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-69', '68') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '68', node2.lemma)
    
    # L -> -69 -> 102 -> L
    feats[978] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-69') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-69', '102') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '102', node2.lemma)
    
    # L -> -69 -> -143 -> L
    feats[979] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-69') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-69', '-143') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-143', node2.lemma)
    
    # L -> -69 -> -104 -> L
    feats[980] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-69') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-69', '-104') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-104', node2.lemma)
    
    # L -> -69 -> -90 -> L
    feats[981] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-69') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-69', '-90') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-90', node2.lemma)
    
    # L -> -69 -> 29 -> L
    feats[982] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-69') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-69', '29') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '29', node2.lemma)
    
    # L -> -69 -> -53 -> L
    feats[983] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-69') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-69', '-53') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-53', node2.lemma)
    
    # L -> 160 -> voc. neutr. -> L
    feats[984] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '160') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '160', 'voc. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. neutr.', node2.lemma)
    
    # L -> tp -> -97 -> L
    feats[985] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '-97') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-97', node2.lemma)
    
    # L -> 160 -> -111 -> L
    feats[986] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '160') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '160', '-111') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-111', node2.lemma)
    
    # L -> 5_pl -> -144 -> L
    feats[987] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_pl', '-144') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-144', node2.lemma)
    
    # L -> 160 -> 121 -> L
    feats[988] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '160') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '160', '121') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '121', node2.lemma)
    
    # L -> 160 -> 30_fp -> L
    feats[989] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '160') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '160', '30_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_fp', node2.lemma)
    
    # L -> tp -> 109 -> L
    feats[990] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '109') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '109', node2.lemma)
    
    # L -> tp -> 119 -> L
    feats[991] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '119') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '119', node2.lemma)
    
    # L -> tp -> -243 -> L
    feats[992] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '-243') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-243', node2.lemma)
    
    # L -> 160 -> 109 -> L
    feats[993] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '160') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '160', '109') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '109', node2.lemma)
    
    # L -> 7_sp -> 181 -> L
    feats[994] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_sp', '181') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '181', node2.lemma)
    
    # L -> -142 -> 2_sg -> L
    feats[995] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-142') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-142', '2_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_sg', node2.lemma)
    
    # L -> tp -> -94 -> L
    feats[996] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '-94') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-94', node2.lemma)
    
    # L -> tp -> acc. pl. -> L
    feats[997] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', 'acc. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. pl.', node2.lemma)
    
    # L -> 160 -> 13_sp -> L
    feats[998] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '160') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '160', '13_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_sp', node2.lemma)
    
    # L -> 160 -> 175 -> L
    feats[999] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '160') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '160', '175') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '175', node2.lemma)
    
    # L -> tp -> sg -> L
    feats[1000] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', 'sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg', node2.lemma)
    
    # L -> tp -> 4_sp -> L
    feats[1001] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '4_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_sp', node2.lemma)
    
    # L -> 5_pl -> -27 -> L
    feats[1002] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_pl', '-27') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-27', node2.lemma)
    
    # L -> 7_sp -> -245 -> L
    feats[1003] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_sp', '-245') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-245', node2.lemma)
    
    # L -> 7_sp -> du_fp -> L
    feats[1004] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_sp', 'du_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_fp', node2.lemma)
    
    # L -> -69 -> -137 -> L
    feats[1005] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-69') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-69', '-137') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-137', node2.lemma)
    
    # L -> -46 -> dat. pl. -> L
    feats[1006] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-46') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-46', 'dat. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. pl.', node2.lemma)
    
    # L -> 5_pl -> 100 -> L
    feats[1007] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_pl', '100') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '100', node2.lemma)
    
    # L -> 16_fp -> 138 -> L
    feats[1008] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_fp', '138') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '138', node2.lemma)
    
    # L -> 7_sp -> 119 -> L
    feats[1009] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_sp', '119') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '119', node2.lemma)
    
    # L -> -46 -> -24 -> L
    feats[1010] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-46') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-46', '-24') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    
    # L -> 179 -> voc. sg. -> L
    feats[1011] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '179') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '179', 'voc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. sg.', node2.lemma)
    
    # L -> 16_fp -> 175 -> L
    feats[1012] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_fp', '175') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '175', node2.lemma)
    
    # L -> 5_pl -> 10_tp -> L
    feats[1013] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_pl', '10_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_tp', node2.lemma)
    
    # L -> -69 -> 4_sp -> L
    feats[1014] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-69') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-69', '4_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_sp', node2.lemma)
    
    # L -> 16_fp -> 97 -> L
    feats[1015] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_fp', '97') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '97', node2.lemma)
    
    # L -> acc. neutr. -> 179 -> L
    feats[1016] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. neutr.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. neutr.', '179') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '179', node2.lemma)
    
    # L -> acc. neutr. -> -276 -> L
    feats[1017] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. neutr.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. neutr.', '-276') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-276', node2.lemma)
    
    # L -> -46 -> -268 -> L
    feats[1018] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-46') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-46', '-268') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-268', node2.lemma)
    
    # L -> -249 -> 4_pl -> L
    feats[1019] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-249') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-249', '4_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_pl', node2.lemma)
    
    # L -> -249 -> 134 -> L
    feats[1020] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-249') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-249', '134') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '134', node2.lemma)
    
    # L -> 179 -> sp -> L
    feats[1021] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '179') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '179', 'sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sp', node2.lemma)
    
    # L -> -152 -> voc. masc. -> L
    feats[1022] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-152') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-152', 'voc. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. masc.', node2.lemma)
    
    # L -> -152 -> -76 -> L
    feats[1023] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-152') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-152', '-76') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-76', node2.lemma)
    
    # L -> acc. neutr. -> -150 -> L
    feats[1024] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. neutr.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. neutr.', '-150') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-150', node2.lemma)
    
    # L -> -19 -> 30 -> L
    feats[1025] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-19') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-19', '30') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30', node2.lemma)
    
    # L -> 179 -> 16_sg -> L
    feats[1026] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '179') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '179', '16_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_sg', node2.lemma)
    
    # L -> 179 -> 30_pl -> L
    feats[1027] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '179') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '179', '30_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_pl', node2.lemma)
    
    # L -> -200 -> 13_sg -> L
    feats[1028] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '13_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_sg', node2.lemma)
    
    # L -> -152 -> -18 -> L
    feats[1029] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-152') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-152', '-18') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-18', node2.lemma)
    
    # L -> acc. neutr. -> 15_tp -> L
    feats[1030] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. neutr.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. neutr.', '15_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_tp', node2.lemma)
    
    # L -> acc. neutr. -> -115 -> L
    feats[1031] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. neutr.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. neutr.', '-115') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-115', node2.lemma)
    
    # L -> -46 -> 91 -> L
    feats[1032] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-46') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-46', '91') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    
    # L -> -19 -> -112 -> L
    feats[1033] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-19') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-19', '-112') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    
    # L -> -152 -> 5_du -> L
    feats[1034] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-152') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-152', '5_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_du', node2.lemma)
    
    # L -> -152 -> 116 -> L
    feats[1035] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-152') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-152', '116') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '116', node2.lemma)
    
    # L -> -200 -> -245 -> L
    feats[1036] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '-245') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-245', node2.lemma)
    
    # L -> -30 -> -153 -> L
    feats[1037] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '-153') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-153', node2.lemma)
    
    # L -> -30 -> -43 -> L
    feats[1038] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '-43') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-43', node2.lemma)
    
    # L -> -276 -> -103 -> L
    feats[1039] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-276') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-276', '-103') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-103', node2.lemma)
    
    # L -> -276 -> nom. sg. -> L
    feats[1040] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-276') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-276', 'nom. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. sg.', node2.lemma)
    
    # L -> -276 -> -22 -> L
    feats[1041] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-276') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-276', '-22') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-22', node2.lemma)
    
    # L -> -276 -> 5_tp -> L
    feats[1042] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-276') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-276', '5_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_tp', node2.lemma)
    
    # L -> -276 -> -263 -> L
    feats[1043] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-276') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-276', '-263') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-263', node2.lemma)
    
    # L -> -276 -> 42 -> L
    feats[1044] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-276') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-276', '42') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '42', node2.lemma)
    
    # L -> -276 -> 13_fp -> L
    feats[1045] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-276') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-276', '13_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_fp', node2.lemma)
    
    # L -> -276 -> 8_sp -> L
    feats[1046] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-276') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-276', '8_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_sp', node2.lemma)
    
    # L -> -276 -> 16_fp -> L
    feats[1047] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-276') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-276', '16_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_fp', node2.lemma)
    
    # L -> -276 -> -247 -> L
    feats[1048] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-276') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-276', '-247') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-247', node2.lemma)
    
    # L -> acc. neutr. -> 30_pl -> L
    feats[1049] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. neutr.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. neutr.', '30_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_pl', node2.lemma)
    
    # L -> acc. neutr. -> 109 -> L
    feats[1050] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. neutr.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. neutr.', '109') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '109', node2.lemma)
    
    # L -> -276 -> -24 -> L
    feats[1051] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-276') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-276', '-24') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    
    # L -> -200 -> 8_du -> L
    feats[1052] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '8_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_du', node2.lemma)
    
    # L -> -276 -> -297 -> L
    feats[1053] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-276') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-276', '-297') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-297', node2.lemma)
    
    # L -> -200 -> -97 -> L
    feats[1054] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '-97') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-97', node2.lemma)
    
    # L -> -276 -> -33 -> L
    feats[1055] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-276') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-276', '-33') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-33', node2.lemma)
    
    # L -> -276 -> voc. sg. -> L
    feats[1056] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-276') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-276', 'voc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. sg.', node2.lemma)
    
    # L -> -276 -> 172 -> L
    feats[1057] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-276') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-276', '172') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '172', node2.lemma)
    
    # L -> -276 -> 181 -> L
    feats[1058] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-276') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-276', '181') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '181', node2.lemma)
    
    # L -> -276 -> 159 -> L
    feats[1059] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-276') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-276', '159') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '159', node2.lemma)
    
    # L -> -276 -> 89 -> L
    feats[1060] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-276') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-276', '89') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '89', node2.lemma)
    
    # L -> -276 -> -291 -> L
    feats[1061] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-276') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-276', '-291') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-291', node2.lemma)
    
    # L -> -276 -> -41 -> L
    feats[1062] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-276') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-276', '-41') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-41', node2.lemma)
    
    # L -> -276 -> -101 -> L
    feats[1063] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-276') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-276', '-101') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-101', node2.lemma)
    
    # L -> -276 -> -147 -> L
    feats[1064] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-276') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-276', '-147') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-147', node2.lemma)
    
    # L -> -276 -> 54 -> L
    feats[1065] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-276') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-276', '54') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    
    # L -> -276 -> 15_tp -> L
    feats[1066] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-276') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-276', '15_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_tp', node2.lemma)
    
    # L -> -276 -> 4_du -> L
    feats[1067] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-276') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-276', '4_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_du', node2.lemma)
    
    # L -> -276 -> 142 -> L
    feats[1068] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-276') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-276', '142') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '142', node2.lemma)
    
    # L -> -276 -> -91 -> L
    feats[1069] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-276') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-276', '-91') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-91', node2.lemma)
    
    # L -> -276 -> -27 -> L
    feats[1070] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-276') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-276', '-27') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-27', node2.lemma)
    
    # L -> -152 -> -29 -> L
    feats[1071] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-152') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-152', '-29') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-29', node2.lemma)
    
    # L -> -308 -> -158 -> L
    feats[1072] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '-158') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-158', node2.lemma)
    
    # L -> -276 -> -268 -> L
    feats[1073] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-276') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-276', '-268') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-268', node2.lemma)
    
    # L -> -276 -> -161 -> L
    feats[1074] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-276') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-276', '-161') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-161', node2.lemma)
    
    # L -> -276 -> loc. pl. -> L
    feats[1075] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-276') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-276', 'loc. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc. pl.', node2.lemma)
    
    # L -> -276 -> 69 -> L
    feats[1076] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-276') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-276', '69') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '69', node2.lemma)
    
    # L -> -276 -> voc -> L
    feats[1077] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-276') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-276', 'voc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc', node2.lemma)
    
    # L -> -30 -> -268 -> L
    feats[1078] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '-268') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-268', node2.lemma)
    
    # L -> -276 -> -230 -> L
    feats[1079] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-276') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-276', '-230') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-230', node2.lemma)
    
    # L -> -30 -> loc. pl. -> L
    feats[1080] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', 'loc. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc. pl.', node2.lemma)
    
    # L -> -276 -> -63 -> L
    feats[1081] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-276') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-276', '-63') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-63', node2.lemma)
    
    # L -> -30 -> -241 -> L
    feats[1082] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '-241') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-241', node2.lemma)
    
    # L -> -30 -> 6_du -> L
    feats[1083] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '6_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_du', node2.lemma)
    
    # L -> -276 -> 12_sg -> L
    feats[1084] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-276') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-276', '12_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_sg', node2.lemma)
    
    # L -> -276 -> sg_fp -> L
    feats[1085] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-276') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-276', 'sg_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_fp', node2.lemma)
    
    # L -> -30 -> nom. masc. -> L
    feats[1086] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', 'nom. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. masc.', node2.lemma)
    
    # L -> 134 -> 80 -> L
    feats[1087] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '134') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '134', '80') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '80', node2.lemma)
    
    # L -> -30 -> 2_sp -> L
    feats[1088] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '2_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_sp', node2.lemma)
    
    # L -> -30 -> -87 -> L
    feats[1089] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '-87') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-87', node2.lemma)
    
    # L -> 34 -> 4_pl -> L
    feats[1090] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '34', '4_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_pl', node2.lemma)
    
    # L -> -30 -> -299 -> L
    feats[1091] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '-299') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-299', node2.lemma)
    
    # L -> -308 -> -46 -> L
    feats[1092] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '-46') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-46', node2.lemma)
    
    # L -> -30 -> 5_du -> L
    feats[1093] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '5_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_du', node2.lemma)
    
    # L -> -30 -> -54 -> L
    feats[1094] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '-54') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-54', node2.lemma)
    
    # L -> -308 -> 29 -> L
    feats[1095] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '29') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '29', node2.lemma)
    
    # L -> pl -> du_fp -> L
    feats[1096] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl', 'du_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_fp', node2.lemma)
    
    # L -> -30 -> -82 -> L
    feats[1097] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '-82') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-82', node2.lemma)
    
    # L -> -30 -> -163 -> L
    feats[1098] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '-163') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-163', node2.lemma)
    
    # L -> -30 -> 76 -> L
    feats[1099] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '76') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '76', node2.lemma)
    
    # L -> -30 -> 74 -> L
    feats[1100] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '74') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '74', node2.lemma)
    
    # L -> -30 -> -149 -> L
    feats[1101] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '-149') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-149', node2.lemma)
    
    # L -> -30 -> 1 -> L
    feats[1102] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '1') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '1', node2.lemma)
    
    # L -> -30 -> -20 -> L
    feats[1103] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '-20') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-20', node2.lemma)
    
    # L -> -30 -> -292 -> L
    feats[1104] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '-292') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-292', node2.lemma)
    
    # L -> -30 -> sg_tp -> L
    feats[1105] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', 'sg_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_tp', node2.lemma)
    
    # L -> -30 -> 150 -> L
    feats[1106] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '150') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '150', node2.lemma)
    
    # L -> -30 -> -21 -> L
    feats[1107] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '-21') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-21', node2.lemma)
    
    # L -> -30 -> -52 -> L
    feats[1108] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '-52') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-52', node2.lemma)
    
    # L -> -30 -> 28_tp -> L
    feats[1109] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '28_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28_tp', node2.lemma)
    
    # L -> -30 -> 11_tp -> L
    feats[1110] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '11_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_tp', node2.lemma)
    
    # L -> -30 -> 100 -> L
    feats[1111] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '100') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '100', node2.lemma)
    
    # L -> -30 -> -111 -> L
    feats[1112] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '-111') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-111', node2.lemma)
    
    # L -> -30 -> 140 -> L
    feats[1113] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '140') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '140', node2.lemma)
    
    # L -> -30 -> 154 -> L
    feats[1114] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '154') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '154', node2.lemma)
    
    # L -> -30 -> 13_pl -> L
    feats[1115] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '13_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_pl', node2.lemma)
    
    # L -> -30 -> 30 -> L
    feats[1116] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '30') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30', node2.lemma)
    
    # L -> -30 -> 174 -> L
    feats[1117] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '174') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '174', node2.lemma)
    
    # L -> -30 -> -309 -> L
    feats[1118] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '-309') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-309', node2.lemma)
    
    # L -> -30 -> 8_du -> L
    feats[1119] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '8_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_du', node2.lemma)
    
    # L -> -30 -> -97 -> L
    feats[1120] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '-97') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-97', node2.lemma)
    
    # L -> 134 -> instr. sg. -> L
    feats[1121] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '134') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '134', 'instr. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. sg.', node2.lemma)
    
    # L -> -30 -> -279 -> L
    feats[1122] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '-279') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-279', node2.lemma)
    
    # L -> -30 -> -44 -> L
    feats[1123] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '-44') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-44', node2.lemma)
    
    # L -> 34 -> -22 -> L
    feats[1124] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '34', '-22') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-22', node2.lemma)
    
    # L -> 134 -> -263 -> L
    feats[1125] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '134') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '134', '-263') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-263', node2.lemma)
    
    # L -> -30 -> 50 -> L
    feats[1126] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '50') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '50', node2.lemma)
    
    # L -> -30 -> -132 -> L
    feats[1127] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '-132') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-132', node2.lemma)
    
    # L -> -30 -> gen. pl. -> L
    feats[1128] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', 'gen. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen. pl.', node2.lemma)
    
    # L -> -30 -> 88 -> L
    feats[1129] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '88') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '88', node2.lemma)
    
    # L -> -30 -> pl_fp -> L
    feats[1130] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', 'pl_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_fp', node2.lemma)
    
    # L -> -30 -> 16_sg -> L
    feats[1131] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '16_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_sg', node2.lemma)
    
    # L -> pl -> -44 -> L
    feats[1132] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl', '-44') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-44', node2.lemma)
    
    # L -> pl -> -77 -> L
    feats[1133] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl', '-77') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-77', node2.lemma)
    
    # L -> -30 -> masc -> L
    feats[1134] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', 'masc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'masc', node2.lemma)
    
    # L -> -30 -> 12_sp -> L
    feats[1135] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '12_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_sp', node2.lemma)
    
    # L -> -30 -> pl_sp -> L
    feats[1136] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', 'pl_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_sp', node2.lemma)
    
    # L -> -30 -> -240 -> L
    feats[1137] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '-240') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-240', node2.lemma)
    
    # L -> -30 -> 2_sg -> L
    feats[1138] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '2_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_sg', node2.lemma)
    
    # L -> -30 -> 151 -> L
    feats[1139] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '151') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '151', node2.lemma)
    
    # L -> -30 -> -99 -> L
    feats[1140] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '-99') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-99', node2.lemma)
    
    # L -> -30 -> -271 -> L
    feats[1141] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '-271') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-271', node2.lemma)
    
    # L -> -30 -> 6_pl -> L
    feats[1142] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '6_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_pl', node2.lemma)
    
    # L -> -30 -> 2_pl -> L
    feats[1143] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '2_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_pl', node2.lemma)
    
    # L -> -30 -> 139 -> L
    feats[1144] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '139') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '139', node2.lemma)
    
    # L -> -30 -> -94 -> L
    feats[1145] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '-94') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-94', node2.lemma)
    
    # L -> -30 -> 4_fp -> L
    feats[1146] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '4_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_fp', node2.lemma)
    
    # L -> -30 -> 30_du -> L
    feats[1147] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '30_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_du', node2.lemma)
    
    # L -> -30 -> acc. pl. -> L
    feats[1148] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', 'acc. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. pl.', node2.lemma)
    
    # L -> -30 -> nom. fem -> L
    feats[1149] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', 'nom. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. fem', node2.lemma)
    
    # L -> -30 -> instr. pl. -> L
    feats[1150] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', 'instr. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. pl.', node2.lemma)
    
    # L -> -30 -> sg -> L
    feats[1151] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', 'sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg', node2.lemma)
    
    # L -> -30 -> 32 -> L
    feats[1152] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '32') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '32', node2.lemma)
    
    # L -> -308 -> -121 -> L
    feats[1153] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '-121') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-121', node2.lemma)
    
    # L -> -308 -> 128 -> L
    feats[1154] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '128') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '128', node2.lemma)
    
    # L -> -30 -> loc. du. -> L
    feats[1155] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', 'loc. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc. du.', node2.lemma)
    
    # L -> -308 -> 15_tp -> L
    feats[1156] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '15_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_tp', node2.lemma)
    
    # L -> -308 -> -115 -> L
    feats[1157] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '-115') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-115', node2.lemma)
    
    # L -> -276 -> -307 -> L
    feats[1158] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-276') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-276', '-307') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-307', node2.lemma)
    
    # L -> -276 -> 130 -> L
    feats[1159] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-276') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-276', '130') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '130', node2.lemma)
    
    # L -> -30 -> -114 -> L
    feats[1160] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '-114') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-114', node2.lemma)
    
    # L -> -30 -> 16_tp -> L
    feats[1161] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '16_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_tp', node2.lemma)
    
    # L -> -308 -> gen -> L
    feats[1162] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', 'gen') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen', node2.lemma)
    
    # L -> -30 -> 169 -> L
    feats[1163] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '169') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '169', node2.lemma)
    
    # L -> -308 -> 7_pl -> L
    feats[1164] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '7_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_pl', node2.lemma)
    
    # L -> -308 -> 110 -> L
    feats[1165] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '110') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '110', node2.lemma)
    
    # L -> -302 -> 4_pl -> L
    feats[1166] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', '4_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_pl', node2.lemma)
    
    # L -> 34 -> abl. du. -> L
    feats[1167] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '34', 'abl. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. du.', node2.lemma)
    
    # L -> -308 -> 131 -> L
    feats[1168] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '131') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '131', node2.lemma)
    
    # L -> -308 -> 69 -> L
    feats[1169] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '69') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '69', node2.lemma)
    
    # L -> -308 -> -113 -> L
    feats[1170] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '-113') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-113', node2.lemma)
    
    # L -> -308 -> -230 -> L
    feats[1171] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '-230') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-230', node2.lemma)
    
    # L -> -302 -> acc. neutr. -> L
    feats[1172] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', 'acc. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. neutr.', node2.lemma)
    
    # L -> -302 -> -17 -> L
    feats[1173] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', '-17') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-17', node2.lemma)
    
    # L -> -308 -> 16_pl -> L
    feats[1174] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '16_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_pl', node2.lemma)
    
    # L -> -302 -> -273 -> L
    feats[1175] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', '-273') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-273', node2.lemma)
    
    # L -> -308 -> -63 -> L
    feats[1176] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '-63') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-63', node2.lemma)
    
    # L -> -308 -> -87 -> L
    feats[1177] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '-87') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-87', node2.lemma)
    
    # L -> -308 -> sg_fp -> L
    feats[1178] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', 'sg_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_fp', node2.lemma)
    
    # L -> -302 -> 30_tp -> L
    feats[1179] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', '30_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_tp', node2.lemma)
    
    # L -> -302 -> dat. pl. -> L
    feats[1180] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', 'dat. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. pl.', node2.lemma)
    
    # L -> -302 -> 2_fp -> L
    feats[1181] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', '2_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_fp', node2.lemma)
    
    # L -> -308 -> -269 -> L
    feats[1182] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '-269') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-269', node2.lemma)
    
    # L -> -302 -> -30 -> L
    feats[1183] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', '-30') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-30', node2.lemma)
    
    # L -> -308 -> 4_tp -> L
    feats[1184] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '4_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_tp', node2.lemma)
    
    # L -> -302 -> -141 -> L
    feats[1185] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', '-141') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-141', node2.lemma)
    
    # L -> -308 -> -129 -> L
    feats[1186] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '-129') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-129', node2.lemma)
    
    # L -> -308 -> -66 -> L
    feats[1187] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '-66') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-66', node2.lemma)
    
    # L -> -308 -> instr. masc. -> L
    feats[1188] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', 'instr. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. masc.', node2.lemma)
    
    # L -> -302 -> -159 -> L
    feats[1189] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', '-159') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-159', node2.lemma)
    
    # L -> -308 -> -82 -> L
    feats[1190] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '-82') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-82', node2.lemma)
    
    # L -> -302 -> -152 -> L
    feats[1191] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', '-152') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-152', node2.lemma)
    
    # L -> -308 -> -149 -> L
    feats[1192] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '-149') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-149', node2.lemma)
    
    # L -> -308 -> 6_fp -> L
    feats[1193] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '6_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_fp', node2.lemma)
    
    # L -> -302 -> instr -> L
    feats[1194] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', 'instr') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr', node2.lemma)
    
    # L -> -302 -> voc. masc. -> L
    feats[1195] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', 'voc. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. masc.', node2.lemma)
    
    # L -> -308 -> 150 -> L
    feats[1196] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '150') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '150', node2.lemma)
    
    # L -> -308 -> 10_du -> L
    feats[1197] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '10_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_du', node2.lemma)
    
    # L -> -302 -> -71 -> L
    feats[1198] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', '-71') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-71', node2.lemma)
    
    # L -> -308 -> 129 -> L
    feats[1199] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '129') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '129', node2.lemma)
    
    # L -> -308 -> 138 -> L
    feats[1200] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '138') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '138', node2.lemma)
    
    # L -> -302 -> 15_pl -> L
    feats[1201] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', '15_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_pl', node2.lemma)
    
    # L -> -308 -> -151 -> L
    feats[1202] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '-151') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-151', node2.lemma)
    
    # L -> -308 -> -23 -> L
    feats[1203] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '-23') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-23', node2.lemma)
    
    # L -> -302 -> instr. adj. -> L
    feats[1204] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', 'instr. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. adj.', node2.lemma)
    
    # L -> -308 -> -296 -> L
    feats[1205] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '-296') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-296', node2.lemma)
    
    # L -> -308 -> 11_tp -> L
    feats[1206] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '11_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_tp', node2.lemma)
    
    # L -> 134 -> 10_sp -> L
    feats[1207] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '134') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '134', '10_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_sp', node2.lemma)
    
    # L -> 134 -> -16 -> L
    feats[1208] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '134') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '134', '-16') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-16', node2.lemma)
    
    # L -> 4_pl -> -62 -> L
    feats[1209] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_pl', '-62') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-62', node2.lemma)
    
    # L -> -302 -> nom. sg. -> L
    feats[1210] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', 'nom. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. sg.', node2.lemma)
    
    # L -> -308 -> 120 -> L
    feats[1211] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '120') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '120', node2.lemma)
    
    # L -> 34 -> -21 -> L
    feats[1212] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '34', '-21') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-21', node2.lemma)
    
    # L -> -302 -> 15_du -> L
    feats[1213] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', '15_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_du', node2.lemma)
    
    # L -> -308 -> 10_sp -> L
    feats[1214] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '10_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_sp', node2.lemma)
    
    # L -> -302 -> fp -> L
    feats[1215] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', 'fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'fp', node2.lemma)
    
    # L -> -308 -> 5_sg -> L
    feats[1216] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '5_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sg', node2.lemma)
    
    # L -> -308 -> -16 -> L
    feats[1217] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '-16') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-16', node2.lemma)
    
    # L -> -302 -> -62 -> L
    feats[1218] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', '-62') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-62', node2.lemma)
    
    # L -> -302 -> 30_sg -> L
    feats[1219] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', '30_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_sg', node2.lemma)
    
    # L -> -308 -> -73 -> L
    feats[1220] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '-73') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-73', node2.lemma)
    
    # L -> -308 -> 99 -> L
    feats[1221] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '99') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '99', node2.lemma)
    
    # L -> -308 -> -97 -> L
    feats[1222] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '-97') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-97', node2.lemma)
    
    # L -> -302 -> acc. adj. -> L
    feats[1223] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', 'acc. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. adj.', node2.lemma)
    
    # L -> -302 -> 13_fp -> L
    feats[1224] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', '13_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_fp', node2.lemma)
    
    # L -> -308 -> -96 -> L
    feats[1225] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '-96') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-96', node2.lemma)
    
    # L -> -308 -> 50 -> L
    feats[1226] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '50') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '50', node2.lemma)
    
    # L -> -308 -> -132 -> L
    feats[1227] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '-132') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-132', node2.lemma)
    
    # L -> -302 -> 8_sp -> L
    feats[1228] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', '8_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_sp', node2.lemma)
    
    # L -> -308 -> -11 -> L
    feats[1229] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '-11') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-11', node2.lemma)
    
    # L -> -308 -> acc. du. -> L
    feats[1230] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', 'acc. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. du.', node2.lemma)
    
    # L -> -302 -> -247 -> L
    feats[1231] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', '-247') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-247', node2.lemma)
    
    # L -> -308 -> gen. pl. -> L
    feats[1232] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', 'gen. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen. pl.', node2.lemma)
    
    # L -> -308 -> 88 -> L
    feats[1233] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '88') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '88', node2.lemma)
    
    # L -> -302 -> 10_sg -> L
    feats[1234] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', '10_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_sg', node2.lemma)
    
    # L -> -308 -> 119 -> L
    feats[1235] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '119') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '119', node2.lemma)
    
    # L -> -302 -> -157 -> L
    feats[1236] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', '-157') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-157', node2.lemma)
    
    # L -> -302 -> -50 -> L
    feats[1237] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', '-50') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-50', node2.lemma)
    
    # L -> -308 -> -86 -> L
    feats[1238] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '-86') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-86', node2.lemma)
    
    # L -> -308 -> masc -> L
    feats[1239] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', 'masc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'masc', node2.lemma)
    
    # L -> -302 -> -297 -> L
    feats[1240] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', '-297') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-297', node2.lemma)
    
    # L -> -308 -> 9_sg -> L
    feats[1241] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '9_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_sg', node2.lemma)
    
    # L -> -308 -> 12_sp -> L
    feats[1242] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '12_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_sp', node2.lemma)
    
    # L -> -302 -> neutr -> L
    feats[1243] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', 'neutr') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'neutr', node2.lemma)
    
    # L -> -302 -> 71 -> L
    feats[1244] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', '71') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '71', node2.lemma)
    
    # L -> -302 -> 3 -> L
    feats[1245] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', '3') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3', node2.lemma)
    
    # L -> -302 -> du_tp -> L
    feats[1246] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', 'du_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_tp', node2.lemma)
    
    # L -> 4_pl -> nom. adj. -> L
    feats[1247] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_pl', 'nom. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. adj.', node2.lemma)
    
    # L -> -302 -> -36 -> L
    feats[1248] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', '-36') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-36', node2.lemma)
    
    # L -> -302 -> 4_du -> L
    feats[1249] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', '4_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_du', node2.lemma)
    
    # L -> 134 -> 8_pl -> L
    feats[1250] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '134') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '134', '8_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_pl', node2.lemma)
    
    # L -> -302 -> 3_fp -> L
    feats[1251] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', '3_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_fp', node2.lemma)
    
    # L -> -302 -> 108 -> L
    feats[1252] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', '108') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '108', node2.lemma)
    
    # L -> -302 -> -102 -> L
    feats[1253] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', '-102') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-102', node2.lemma)
    
    # L -> 4_pl -> 110 -> L
    feats[1254] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_pl', '110') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '110', node2.lemma)
    
    # L -> -302 -> 7_pl -> L
    feats[1255] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', '7_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_pl', node2.lemma)
    
    # L -> -302 -> 110 -> L
    feats[1256] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', '110') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '110', node2.lemma)
    
    # L -> -247 -> 9_fp -> L
    feats[1257] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-247') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-247', '9_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_fp', node2.lemma)
    
    # L -> -302 -> 131 -> L
    feats[1258] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', '131') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '131', node2.lemma)
    
    # L -> -302 -> -230 -> L
    feats[1259] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', '-230') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-230', node2.lemma)
    
    # L -> 4_pl -> -87 -> L
    feats[1260] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_pl', '-87') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-87', node2.lemma)
    
    # L -> -302 -> 16_pl -> L
    feats[1261] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', '16_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_pl', node2.lemma)
    
    # L -> 4_pl -> 28 -> L
    feats[1262] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_pl', '28') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
    
    # L -> -302 -> sg_fp -> L
    feats[1263] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', 'sg_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_fp', node2.lemma)
    
    # L -> 4_pl -> 28_sg -> L
    feats[1264] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_pl', '28_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28_sg', node2.lemma)
    
    # L -> -302 -> -245 -> L
    feats[1265] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', '-245') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-245', node2.lemma)
    
    # L -> 4_pl -> pl_tp -> L
    feats[1266] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_pl', 'pl_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_tp', node2.lemma)
    
    # L -> 4_pl -> 5_du -> L
    feats[1267] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_pl', '5_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_du', node2.lemma)
    
    # L -> -302 -> 180 -> L
    feats[1268] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', '180') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '180', node2.lemma)
    
    # L -> 4_pl -> instr. masc. -> L
    feats[1269] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_pl', 'instr. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. masc.', node2.lemma)
    
    # L -> -302 -> 116 -> L
    feats[1270] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', '116') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '116', node2.lemma)
    
    # L -> 4_pl -> 76 -> L
    feats[1271] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_pl', '76') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '76', node2.lemma)
    
    # L -> -302 -> -163 -> L
    feats[1272] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', '-163') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-163', node2.lemma)
    
    # L -> -302 -> 40 -> L
    feats[1273] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', '40') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '40', node2.lemma)
    
    # L -> 4_pl -> 150 -> L
    feats[1274] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_pl', '150') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '150', node2.lemma)
    
    # L -> -302 -> 6_fp -> L
    feats[1275] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', '6_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_fp', node2.lemma)
    
    # L -> -302 -> -20 -> L
    feats[1276] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', '-20') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-20', node2.lemma)
    
    # L -> -302 -> -292 -> L
    feats[1277] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', '-292') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-292', node2.lemma)
    
    # L -> 4_pl -> -25 -> L
    feats[1278] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_pl', '-25') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-25', node2.lemma)
    
    # L -> -302 -> sg_tp -> L
    feats[1279] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', 'sg_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_tp', node2.lemma)
    
    # L -> 4_pl -> -283 -> L
    feats[1280] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_pl', '-283') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-283', node2.lemma)
    
    # L -> -17 -> -260 -> L
    feats[1281] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-17') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-17', '-260') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-260', node2.lemma)
    
    # L -> -17 -> 3_du -> L
    feats[1282] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-17') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-17', '3_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    
    # L -> -56 -> 95 -> L
    feats[1283] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-56') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-56', '95') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '95', node2.lemma)
    
    # L -> -302 -> -49 -> L
    feats[1284] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', '-49') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-49', node2.lemma)
    
    # L -> 8_sg -> -293 -> L
    feats[1285] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_sg', '-293') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-293', node2.lemma)
    
    # L -> 8_sg -> voc. masc. -> L
    feats[1286] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_sg', 'voc. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. masc.', node2.lemma)
    
    # L -> -302 -> 140 -> L
    feats[1287] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', '140') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '140', node2.lemma)
    
    # L -> -17 -> 162 -> L
    feats[1288] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-17') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-17', '162') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '162', node2.lemma)
    
    # L -> -302 -> -83 -> L
    feats[1289] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', '-83') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-83', node2.lemma)
    
    # L -> -302 -> 120 -> L
    feats[1290] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', '120') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '120', node2.lemma)
    
    # L -> -17 -> 177 -> L
    feats[1291] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-17') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-17', '177') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '177', node2.lemma)
    
    # L -> 4_pl -> -42 -> L
    feats[1292] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_pl', '-42') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-42', node2.lemma)
    
    # L -> -302 -> 121 -> L
    feats[1293] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', '121') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '121', node2.lemma)
    
    # L -> -17 -> 5_tp -> L
    feats[1294] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-17') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-17', '5_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_tp', node2.lemma)
    
    # L -> -17 -> -35 -> L
    feats[1295] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-17') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-17', '-35') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    
    # L -> 4_pl -> 99 -> L
    feats[1296] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_pl', '99') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '99', node2.lemma)
    
    # L -> 4_pl -> 14_pl -> L
    feats[1297] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_pl', '14_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_pl', node2.lemma)
    
    # L -> 4_pl -> -309 -> L
    feats[1298] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_pl', '-309') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-309', node2.lemma)
    
    # L -> -17 -> acc. masc. -> L
    feats[1299] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-17') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-17', 'acc. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. masc.', node2.lemma)
    
    # L -> 4_pl -> -97 -> L
    feats[1300] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_pl', '-97') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-97', node2.lemma)
    
    # L -> -302 -> 174 -> L
    feats[1301] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', '174') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '174', node2.lemma)
    
    # L -> -302 -> 99 -> L
    feats[1302] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', '99') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '99', node2.lemma)
    
    # L -> -17 -> 29_tp -> L
    feats[1303] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-17') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-17', '29_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '29_tp', node2.lemma)
    
    # L -> -17 -> 94 -> L
    feats[1304] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-17') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-17', '94') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '94', node2.lemma)
    
    # L -> 4_pl -> 2_tp -> L
    feats[1305] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_pl', '2_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_tp', node2.lemma)
    
    # L -> 4_pl -> masc -> L
    feats[1306] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_pl', 'masc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'masc', node2.lemma)
    
    # L -> -17 -> 13_sg -> L
    feats[1307] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-17') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-17', '13_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_sg', node2.lemma)
    
    # L -> -17 -> -157 -> L
    feats[1308] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-17') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-17', '-157') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-157', node2.lemma)
    
    # L -> -17 -> abl. sg. -> L
    feats[1309] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-17') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-17', 'abl. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. sg.', node2.lemma)
    
    # L -> -17 -> 170 -> L
    feats[1310] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-17') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-17', '170') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '170', node2.lemma)
    
    # L -> -17 -> -19 -> L
    feats[1311] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-17') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-17', '-19') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-19', node2.lemma)
    
    # L -> -17 -> neutr -> L
    feats[1312] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-17') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-17', 'neutr') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'neutr', node2.lemma)
    
    # L -> -17 -> 141 -> L
    feats[1313] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-17') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-17', '141') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '141', node2.lemma)
    
    # L -> 4_pl -> -243 -> L
    feats[1314] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_pl', '-243') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-243', node2.lemma)
    
    # L -> 4_pl -> du_sp -> L
    feats[1315] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_pl', 'du_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_sp', node2.lemma)
    
    # L -> -17 -> 172 -> L
    feats[1316] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-17') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-17', '172') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '172', node2.lemma)
    
    # L -> -17 -> 181 -> L
    feats[1317] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-17') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-17', '181') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '181', node2.lemma)
    
    # L -> -17 -> 56 -> L
    feats[1318] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-17') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-17', '56') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '56', node2.lemma)
    
    # L -> -17 -> du_tp -> L
    feats[1319] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-17') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-17', 'du_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_tp', node2.lemma)
    
    # L -> -17 -> 135 -> L
    feats[1320] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-17') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-17', '135') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '135', node2.lemma)
    
    # L -> -17 -> -101 -> L
    feats[1321] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-17') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-17', '-101') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-101', node2.lemma)
    
    # L -> 8_sg -> 71 -> L
    feats[1322] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_sg', '71') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '71', node2.lemma)
    
    # L -> -17 -> gen. sg. -> L
    feats[1323] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-17') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-17', 'gen. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen. sg.', node2.lemma)
    
    # L -> -17 -> 11_du -> L
    feats[1324] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-17') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-17', '11_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_du', node2.lemma)
    
    # L -> -17 -> abl. du. -> L
    feats[1325] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-17') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-17', 'abl. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. du.', node2.lemma)
    
    # L -> -17 -> 3_fp -> L
    feats[1326] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-17') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-17', '3_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_fp', node2.lemma)
    
    # L -> -17 -> instr. fem -> L
    feats[1327] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-17') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-17', 'instr. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. fem', node2.lemma)
    
    # L -> -17 -> 7_pl -> L
    feats[1328] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-17') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-17', '7_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_pl', node2.lemma)
    
    # L -> -17 -> -61 -> L
    feats[1329] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-17') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-17', '-61') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-61', node2.lemma)
    
    # L -> -17 -> 39 -> L
    feats[1330] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-17') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-17', '39') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '39', node2.lemma)
    
    # L -> -302 -> -29 -> L
    feats[1331] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', '-29') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-29', node2.lemma)
    
    # L -> -17 -> -27 -> L
    feats[1332] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-17') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-17', '-27') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-27', node2.lemma)
    
    # L -> -17 -> -137 -> L
    feats[1333] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-17') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-17', '-137') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-137', node2.lemma)
    
    # L -> neutr -> 49 -> L
    feats[1334] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'neutr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'neutr', '49') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '49', node2.lemma)
    
    # L -> -17 -> 6_sp -> L
    feats[1335] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-17') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-17', '6_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_sp', node2.lemma)
    
    # L -> -17 -> -161 -> L
    feats[1336] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-17') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-17', '-161') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-161', node2.lemma)
    
    # L -> neutr -> 16_tp -> L
    feats[1337] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'neutr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'neutr', '16_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_tp', node2.lemma)
    
    # L -> -17 -> 70 -> L
    feats[1338] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-17') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-17', '70') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '70', node2.lemma)
    
    # L -> -56 -> 171 -> L
    feats[1339] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-56') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-56', '171') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '171', node2.lemma)
    
    # L -> -17 -> 6_du -> L
    feats[1340] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-17') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-17', '6_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_du', node2.lemma)
    
    # L -> -17 -> -245 -> L
    feats[1341] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-17') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-17', '-245') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-245', node2.lemma)
    
    # L -> -17 -> du_fp -> L
    feats[1342] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-17') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-17', 'du_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_fp', node2.lemma)
    
    # L -> -17 -> 28 -> L
    feats[1343] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-17') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-17', '28') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
    
    # L -> -17 -> 182 -> L
    feats[1344] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-17') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-17', '182') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '182', node2.lemma)
    
    # L -> 141 -> 92 -> L
    feats[1345] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '141') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '141', '92') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '92', node2.lemma)
    
    # L -> -17 -> -269 -> L
    feats[1346] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-17') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-17', '-269') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-269', node2.lemma)
    
    # L -> 141 -> -190 -> L
    feats[1347] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '141') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '141', '-190') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-190', node2.lemma)
    
    # L -> -17 -> -154 -> L
    feats[1348] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-17') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-17', '-154') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-154', node2.lemma)
    
    # L -> 141 -> 5_sp -> L
    feats[1349] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '141') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '141', '5_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sp', node2.lemma)
    
    # L -> -17 -> pl_tp -> L
    feats[1350] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-17') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-17', 'pl_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_tp', node2.lemma)
    
    # L -> -17 -> 5_du -> L
    feats[1351] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-17') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-17', '5_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_du', node2.lemma)
    
    # L -> 141 -> acc. neutr. -> L
    feats[1352] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '141') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '141', 'acc. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. neutr.', node2.lemma)
    
    # L -> -17 -> -54 -> L
    feats[1353] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-17') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-17', '-54') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-54', node2.lemma)
    
    # L -> -17 -> -129 -> L
    feats[1354] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-17') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-17', '-129') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-129', node2.lemma)
    
    # L -> 141 -> voc. fem -> L
    feats[1355] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '141') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '141', 'voc. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. fem', node2.lemma)
    
    # L -> -17 -> -82 -> L
    feats[1356] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-17') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-17', '-82') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-82', node2.lemma)
    
    # L -> 8_sg -> -59 -> L
    feats[1357] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_sg', '-59') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-59', node2.lemma)
    
    # L -> 141 -> -142 -> L
    feats[1358] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '141') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '141', '-142') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-142', node2.lemma)
    
    # L -> 141 -> -104 -> L
    feats[1359] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '141') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '141', '-104') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-104', node2.lemma)
    
    # L -> 141 -> 48 -> L
    feats[1360] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '141') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '141', '48') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '48', node2.lemma)
    
    # L -> -143 -> 2 -> L
    feats[1361] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-143') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-143', '2') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2', node2.lemma)
    
    # L -> -143 -> 60 -> L
    feats[1362] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-143') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-143', '60') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '60', node2.lemma)
    
    # L -> 141 -> 5_pl -> L
    feats[1363] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '141') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '141', '5_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_pl', node2.lemma)
    
    # L -> 141 -> -76 -> L
    feats[1364] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '141') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '141', '-76') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-76', node2.lemma)
    
    # L -> 141 -> -57 -> L
    feats[1365] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '141') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '141', '-57') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-57', node2.lemma)
    
    # L -> 141 -> 3_du -> L
    feats[1366] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '141') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '141', '3_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    
    # L -> 141 -> sg_sp -> L
    feats[1367] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '141') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '141', 'sg_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_sp', node2.lemma)
    
    # L -> 141 -> -262 -> L
    feats[1368] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '141') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '141', '-262') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-262', node2.lemma)
    
    # L -> 141 -> 15_du -> L
    feats[1369] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '141') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '141', '15_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_du', node2.lemma)
    
    # L -> 141 -> fp -> L
    feats[1370] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '141') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '141', 'fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'fp', node2.lemma)
    
    # L -> 141 -> 176 -> L
    feats[1371] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '141') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '141', '176') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '176', node2.lemma)
    
    # L -> 141 -> -301 -> L
    feats[1372] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '141') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '141', '-301') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-301', node2.lemma)
    
    # L -> 141 -> -150 -> L
    feats[1373] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '141') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '141', '-150') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-150', node2.lemma)
    
    # L -> 141 -> du -> L
    feats[1374] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '141') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '141', 'du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du', node2.lemma)
    
    # L -> 141 -> 31 -> L
    feats[1375] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '141') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '141', '31') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '31', node2.lemma)
    
    # L -> -17 -> -32 -> L
    feats[1376] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-17') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-17', '-32') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-32', node2.lemma)
    
    # L -> 8_sg -> masc -> L
    feats[1377] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_sg', 'masc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'masc', node2.lemma)
    
    # L -> 61 -> 173 -> L
    feats[1378] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '61') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '61', '173') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '173', node2.lemma)
    
    # L -> voc. fem -> 80 -> L
    feats[1379] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. fem', '80') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '80', node2.lemma)
    
    # L -> 94 -> -126 -> L
    feats[1380] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '94') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '94', '-126') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-126', node2.lemma)
    
    # L -> 94 -> 148 -> L
    feats[1381] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '94') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '94', '148') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '148', node2.lemma)
    
    # L -> 71 -> -126 -> L
    feats[1382] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '71') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '71', '-126') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-126', node2.lemma)
    
    # L -> voc. fem -> 13_tp -> L
    feats[1383] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. fem', '13_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_tp', node2.lemma)
    
    # L -> -143 -> -25 -> L
    feats[1384] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-143') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-143', '-25') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-25', node2.lemma)
    
    # L -> voc. fem -> 30_sg -> L
    feats[1385] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. fem', '30_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_sg', node2.lemma)
    
    # L -> 61 -> 30_fp -> L
    feats[1386] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '61') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '61', '30_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_fp', node2.lemma)
    
    # L -> 61 -> 174 -> L
    feats[1387] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '61') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '61', '174') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '174', node2.lemma)
    
    # L -> 71 -> nom. sg. -> L
    feats[1388] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '71') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '71', 'nom. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. sg.', node2.lemma)
    
    # L -> voc. fem -> 6_sg -> L
    feats[1389] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. fem', '6_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_sg', node2.lemma)
    
    # L -> voc. fem -> 13_sg -> L
    feats[1390] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. fem', '13_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_sg', node2.lemma)
    
    # L -> voc. fem -> -123 -> L
    feats[1391] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. fem', '-123') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-123', node2.lemma)
    
    # L -> voc. fem -> voc. sg. -> L
    feats[1392] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. fem', 'voc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. sg.', node2.lemma)
    
    # L -> voc. fem -> 114 -> L
    feats[1393] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. fem', '114') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '114', node2.lemma)
    
    # L -> -143 -> -32 -> L
    feats[1394] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-143') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-143', '-32') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-32', node2.lemma)
    
    # L -> voc. fem -> 3 -> L
    feats[1395] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. fem', '3') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3', node2.lemma)
    
    # L -> voc. fem -> -291 -> L
    feats[1396] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. fem', '-291') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-291', node2.lemma)
    
    # L -> voc. fem -> -121 -> L
    feats[1397] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. fem', '-121') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-121', node2.lemma)
    
    # L -> voc. fem -> 15_tp -> L
    feats[1398] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. fem', '15_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_tp', node2.lemma)
    
    # L -> voc. fem -> -115 -> L
    feats[1399] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. fem', '-115') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-115', node2.lemma)
    
    # L -> 141 -> acc. sg. -> L
    feats[1400] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '141') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '141', 'acc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. sg.', node2.lemma)
    
    # L -> 94 -> -245 -> L
    feats[1401] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '94') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '94', '-245') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-245', node2.lemma)
    
    # L -> 92 -> voc. masc. -> L
    feats[1402] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '92', 'voc. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. masc.', node2.lemma)
    
    # L -> 79 -> -152 -> L
    feats[1403] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '-152') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-152', node2.lemma)
    
    # L -> 71 -> -25 -> L
    feats[1404] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '71') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '71', '-25') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-25', node2.lemma)
    
    # L -> -104 -> -144 -> L
    feats[1405] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '-144') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-144', node2.lemma)
    
    # L -> 94 -> -42 -> L
    feats[1406] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '94') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '94', '-42') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-42', node2.lemma)
    
    # L -> 94 -> -73 -> L
    feats[1407] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '94') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '94', '-73') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-73', node2.lemma)
    
    # L -> 71 -> -32 -> L
    feats[1408] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '71') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '71', '-32') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-32', node2.lemma)
    
    # L -> 79 -> -200 -> L
    feats[1409] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '-200') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-200', node2.lemma)
    
    # L -> 10_fp -> acc. fem -> L
    feats[1410] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_fp', 'acc. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. fem', node2.lemma)
    
    # L -> -104 -> 39 -> L
    feats[1411] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '39') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '39', node2.lemma)
    
    # L -> 94 -> 130 -> L
    feats[1412] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '94') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '94', '130') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '130', node2.lemma)
    
    # L -> 92 -> 180 -> L
    feats[1413] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '92', '180') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '180', node2.lemma)
    
    # L -> 92 -> -154 -> L
    feats[1414] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '92', '-154') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-154', node2.lemma)
    
    # L -> 79 -> -59 -> L
    feats[1415] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '-59') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-59', node2.lemma)
    
    # L -> -273 -> 3_tp -> L
    feats[1416] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '3_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_tp', node2.lemma)
    
    # L -> 3 -> -143 -> L
    feats[1417] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3', '-143') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-143', node2.lemma)
    
    # L -> -104 -> 75 -> L
    feats[1418] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '75') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '75', node2.lemma)
    
    # L -> -104 -> -49 -> L
    feats[1419] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '-49') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-49', node2.lemma)
    
    # L -> -58 -> 59 -> L
    feats[1420] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-58', '59') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '59', node2.lemma)
    
    # L -> 92 -> 12_sp -> L
    feats[1421] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '92', '12_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_sp', node2.lemma)
    
    # L -> -273 -> 160 -> L
    feats[1422] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '160') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    
    # L -> 3 -> 42 -> L
    feats[1423] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3', '42') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '42', node2.lemma)
    
    # L -> -273 -> 37 -> L
    feats[1424] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '37') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '37', node2.lemma)
    
    # L -> -273 -> 159 -> L
    feats[1425] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '159') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '159', node2.lemma)
    
    # L -> -273 -> -56 -> L
    feats[1426] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '-56') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-56', node2.lemma)
    
    # L -> -273 -> 71 -> L
    feats[1427] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '71') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '71', node2.lemma)
    
    # L -> -273 -> 3 -> L
    feats[1428] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '3') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3', node2.lemma)
    
    # L -> -273 -> 56 -> L
    feats[1429] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '56') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '56', node2.lemma)
    
    # L -> -273 -> -10 -> L
    feats[1430] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '-10') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-10', node2.lemma)
    
    # L -> -273 -> -220 -> L
    feats[1431] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '-220') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-220', node2.lemma)
    
    # L -> -273 -> nom. neutr. -> L
    feats[1432] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', 'nom. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. neutr.', node2.lemma)
    
    # L -> -273 -> nom. adj. -> L
    feats[1433] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', 'nom. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. adj.', node2.lemma)
    
    # L -> -273 -> -36 -> L
    feats[1434] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '-36') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-36', node2.lemma)
    
    # L -> -273 -> -121 -> L
    feats[1435] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '-121') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-121', node2.lemma)
    
    # L -> -273 -> 7_tp -> L
    feats[1436] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '7_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_tp', node2.lemma)
    
    # L -> -273 -> 128 -> L
    feats[1437] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '128') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '128', node2.lemma)
    
    # L -> -273 -> -147 -> L
    feats[1438] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '-147') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-147', node2.lemma)
    
    # L -> -273 -> -166 -> L
    feats[1439] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '-166') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-166', node2.lemma)
    
    # L -> -273 -> -12 -> L
    feats[1440] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '-12') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-12', node2.lemma)
    
    # L -> -273 -> 7_fp -> L
    feats[1441] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '7_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_fp', node2.lemma)
    
    # L -> -273 -> -115 -> L
    feats[1442] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '-115') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-115', node2.lemma)
    
    # L -> -273 -> 4_du -> L
    feats[1443] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '4_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_du', node2.lemma)
    
    # L -> -273 -> -91 -> L
    feats[1444] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '-91') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-91', node2.lemma)
    
    # L -> -273 -> 3_fp -> L
    feats[1445] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '3_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_fp', node2.lemma)
    
    # L -> -58 -> -41 -> L
    feats[1446] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-58', '-41') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-41', node2.lemma)
    
    # L -> -273 -> 108 -> L
    feats[1447] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '108') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '108', node2.lemma)
    
    # L -> -273 -> 7_pl -> L
    feats[1448] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '7_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_pl', node2.lemma)
    
    # L -> -273 -> dat. sg. -> L
    feats[1449] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', 'dat. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. sg.', node2.lemma)
    
    # L -> -273 -> 110 -> L
    feats[1450] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '110') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '110', node2.lemma)
    
    # L -> -273 -> -27 -> L
    feats[1451] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '-27') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-27', node2.lemma)
    
    # L -> -273 -> -137 -> L
    feats[1452] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '-137') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-137', node2.lemma)
    
    # L -> -273 -> -122 -> L
    feats[1453] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '-122') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-122', node2.lemma)
    
    # L -> -273 -> -48 -> L
    feats[1454] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '-48') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-48', node2.lemma)
    
    # L -> -273 -> -303 -> L
    feats[1455] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '-303') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-303', node2.lemma)
    
    # L -> -273 -> 122 -> L
    feats[1456] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '122') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '122', node2.lemma)
    
    # L -> -273 -> 111 -> L
    feats[1457] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '111') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '111', node2.lemma)
    
    # L -> -273 -> 9_tp -> L
    feats[1458] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '9_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_tp', node2.lemma)
    
    # L -> -273 -> 131 -> L
    feats[1459] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '131') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '131', node2.lemma)
    
    # L -> -273 -> -268 -> L
    feats[1460] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '-268') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-268', node2.lemma)
    
    # L -> -273 -> -161 -> L
    feats[1461] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '-161') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-161', node2.lemma)
    
    # L -> -126 -> 5_sp -> L
    feats[1462] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-126') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-126', '5_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sp', node2.lemma)
    
    # L -> -273 -> loc. pl. -> L
    feats[1463] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', 'loc. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc. pl.', node2.lemma)
    
    # L -> -273 -> 69 -> L
    feats[1464] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '69') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '69', node2.lemma)
    
    # L -> -126 -> fem -> L
    feats[1465] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-126') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-126', 'fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'fem', node2.lemma)
    
    # L -> -273 -> voc -> L
    feats[1466] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', 'voc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc', node2.lemma)
    
    # L -> -273 -> -113 -> L
    feats[1467] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '-113') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-113', node2.lemma)
    
    # L -> -273 -> 173 -> L
    feats[1468] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '173') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '173', node2.lemma)
    
    # L -> -273 -> 16_pl -> L
    feats[1469] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '16_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_pl', node2.lemma)
    
    # L -> -273 -> nom. masc. -> L
    feats[1470] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', 'nom. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. masc.', node2.lemma)
    
    # L -> -273 -> -63 -> L
    feats[1471] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '-63') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-63', node2.lemma)
    
    # L -> -273 -> sg_fp -> L
    feats[1472] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', 'sg_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_fp', node2.lemma)
    
    # L -> -273 -> du_fp -> L
    feats[1473] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', 'du_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_fp', node2.lemma)
    
    # L -> -273 -> 182 -> L
    feats[1474] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '182') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '182', node2.lemma)
    
    # L -> -273 -> -59 -> L
    feats[1475] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '-59') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-59', node2.lemma)
    
    # L -> -273 -> -299 -> L
    feats[1476] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '-299') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-299', node2.lemma)
    
    # L -> -273 -> 180 -> L
    feats[1477] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '180') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '180', node2.lemma)
    
    # L -> -273 -> -154 -> L
    feats[1478] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '-154') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-154', node2.lemma)
    
    # L -> -273 -> 4_tp -> L
    feats[1479] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '4_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_tp', node2.lemma)
    
    # L -> -273 -> 78 -> L
    feats[1480] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '78') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '78', node2.lemma)
    
    # L -> -273 -> 152 -> L
    feats[1481] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '152') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '152', node2.lemma)
    
    # L -> -273 -> 14_sp -> L
    feats[1482] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '14_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_sp', node2.lemma)
    
    # L -> -273 -> 5_du -> L
    feats[1483] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '5_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_du', node2.lemma)
    
    # L -> -273 -> 116 -> L
    feats[1484] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '116') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '116', node2.lemma)
    
    # L -> -273 -> 82 -> L
    feats[1485] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '82') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '82', node2.lemma)
    
    # L -> -273 -> sp -> L
    feats[1486] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', 'sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sp', node2.lemma)
    
    # L -> -273 -> instr. masc. -> L
    feats[1487] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', 'instr. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. masc.', node2.lemma)
    
    # L -> -273 -> -93 -> L
    feats[1488] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '-93') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-93', node2.lemma)
    
    # L -> 48 -> 13_tp -> L
    feats[1489] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '48', '13_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_tp', node2.lemma)
    
    # L -> 48 -> -53 -> L
    feats[1490] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '48', '-53') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-53', node2.lemma)
    
    # L -> -58 -> pl_tp -> L
    feats[1491] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-58', 'pl_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_tp', node2.lemma)
    
    # L -> -58 -> 5_du -> L
    feats[1492] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-58', '5_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_du', node2.lemma)
    
    # L -> 3 -> 8_du -> L
    feats[1493] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3', '8_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_du', node2.lemma)
    
    # L -> 3 -> 14_sg -> L
    feats[1494] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3', '14_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_sg', node2.lemma)
    
    # L -> -273 -> -240 -> L
    feats[1495] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '-240') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-240', node2.lemma)
    
    # L -> 3 -> -132 -> L
    feats[1496] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3', '-132') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-132', node2.lemma)
    
    # L -> 48 -> 37 -> L
    feats[1497] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '48', '37') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '37', node2.lemma)
    
    # L -> -58 -> -240 -> L
    feats[1498] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-58', '-240') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-240', node2.lemma)
    
    # L -> 3 -> 11_sp -> L
    feats[1499] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3', '11_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_sp', node2.lemma)
    
    feats[feats < 1e-25] = 1e-25
    return -np.log10(feats)
    