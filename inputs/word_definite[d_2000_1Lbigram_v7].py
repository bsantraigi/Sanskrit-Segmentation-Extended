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

def word_definite_extInit(matDB, _edge_vector_dim_from_universe):
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
    _edge_vector_dim = _edge_vector_dim_from_universe
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
            
    # L -> 13_pl -> -150 -> L
    feats[1] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '-150') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-150', node2.lemma)
    
    # L->T
    feats[2] = tryProb_catchZero(mat_lem2tup_countonly, mat_lemCount_1D, node1.lemma, node2.tup)
            
    # C->C
    feats[3] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, node2.cng)
            
    # L -> 13_pl -> -293 -> L
    feats[4] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '-293') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-293', node2.lemma)
    
    # L -> 13_pl -> instr. pl. -> L
    feats[5] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', 'instr. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. pl.', node2.lemma)
    
    # L -> 13_pl -> 176 -> L
    feats[6] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '176') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '176', node2.lemma)
    
    # L -> 13_pl -> 160 -> L
    feats[7] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '160') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    
    # L -> 75 -> -26 -> L
    feats[8] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '75') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '75', '-26') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-26', node2.lemma)
    
    # L -> 75 -> fem -> L
    feats[9] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '75') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '75', 'fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'fem', node2.lemma)
    
    # L -> 13_pl -> 2_tp -> L
    feats[10] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '2_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_tp', node2.lemma)
    
    # L -> 77 -> sg_tp -> L
    feats[11] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '77') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '77', 'sg_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_tp', node2.lemma)
    
    # L -> 12_fp -> 3_du -> L
    feats[12] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_fp', '3_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    
    # L -> -132 -> 3_sg -> L
    feats[13] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-132') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-132', '3_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_sg', node2.lemma)
    
    # L -> -132 -> -47 -> L
    feats[14] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-132') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-132', '-47') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-47', node2.lemma)
    
    # L -> 35 -> 5_fp -> L
    feats[15] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', '5_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_fp', node2.lemma)
    
    # L -> 77 -> 5_fp -> L
    feats[16] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '77') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '77', '5_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_fp', node2.lemma)
    
    # L -> 77 -> 8_sp -> L
    feats[17] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '77') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '77', '8_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_sp', node2.lemma)
    
    # L->-69->L
    feats[18] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-69') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-69', node2.lemma)
            
    # L->72->L
    feats[19] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '72') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '72', node2.lemma)
            
    # L -> 35 -> dat. pl. -> L
    feats[20] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', 'dat. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. pl.', node2.lemma)
    
    # L -> 35 -> sp -> L
    feats[21] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', 'sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sp', node2.lemma)
    
    # L -> 77 -> dat. pl. -> L
    feats[22] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '77') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '77', 'dat. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. pl.', node2.lemma)
    
    # L -> 35 -> -69 -> L
    feats[23] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', '-69') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-69', node2.lemma)
    
    # L->-241->L
    feats[24] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-241') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-241', node2.lemma)
            
    # L -> 77 -> 14_pl -> L
    feats[25] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '77') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '77', '14_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_pl', node2.lemma)
    
    # L->-17->L
    feats[26] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-17') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-17', node2.lemma)
            
    # L -> 77 -> -90 -> L
    feats[27] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '77') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '77', '-90') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-90', node2.lemma)
    
    # L -> 35 -> -301 -> L
    feats[28] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', '-301') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-301', node2.lemma)
    
    # L->-99->L
    feats[29] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-99') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-99', node2.lemma)
            
    # L->117->L
    feats[30] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '117') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '117', node2.lemma)
            
    # L -> 12_fp -> 10_sg -> L
    feats[31] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_fp', '10_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_sg', node2.lemma)
    
    # L -> 77 -> -84 -> L
    feats[32] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '77') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '77', '-84') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-84', node2.lemma)
    
    # L -> 12_fp -> 3_fp -> L
    feats[33] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_fp', '3_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_fp', node2.lemma)
    
    # L -> 35 -> 13_tp -> L
    feats[34] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', '13_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_tp', node2.lemma)
    
    # L->-31->L
    feats[35] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-31') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-31', node2.lemma)
            
    # L -> 75 -> acc. adj. -> L
    feats[36] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '75') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '75', 'acc. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. adj.', node2.lemma)
    
    # L -> 35 -> -17 -> L
    feats[37] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', '-17') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-17', node2.lemma)
    
    # L -> 77 -> 13_tp -> L
    feats[38] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '77') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '77', '13_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_tp', node2.lemma)
    
    # L -> 35 -> -99 -> L
    feats[39] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', '-99') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-99', node2.lemma)
    
    # L -> 77 -> abl. sg. -> L
    feats[40] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '77') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '77', 'abl. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. sg.', node2.lemma)
    
    # L -> 77 -> 117 -> L
    feats[41] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '77') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '77', '117') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '117', node2.lemma)
    
    # L -> 77 -> -268 -> L
    feats[42] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '77') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '77', '-268') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-268', node2.lemma)
    
    # L -> 77 -> -302 -> L
    feats[43] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '77') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '77', '-302') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-302', node2.lemma)
    
    # L -> 35 -> 137 -> L
    feats[44] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', '137') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    
    # L -> 35 -> -308 -> L
    feats[45] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', '-308') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-308', node2.lemma)
    
    # L -> 77 -> gen. pl. -> L
    feats[46] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '77') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '77', 'gen. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen. pl.', node2.lemma)
    
    # L -> 77 -> -48 -> L
    feats[47] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '77') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '77', '-48') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-48', node2.lemma)
    
    # L -> 77 -> -308 -> L
    feats[48] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '77') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '77', '-308') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-308', node2.lemma)
    
    # L -> 77 -> -240 -> L
    feats[49] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '77') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '77', '-240') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-240', node2.lemma)
    
    # L -> 35 -> 179 -> L
    feats[50] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', '179') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '179', node2.lemma)
    
    # L -> 75 -> 97 -> L
    feats[51] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '75') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '75', '97') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '97', node2.lemma)
    
    # L -> 75 -> 10_fp -> L
    feats[52] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '75') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '75', '10_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_fp', node2.lemma)
    
    # L -> 35 -> 9_fp -> L
    feats[53] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', '9_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_fp', node2.lemma)
    
    # L -> 35 -> 58 -> L
    feats[54] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', '58') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '58', node2.lemma)
    
    # L -> 35 -> 10_sp -> L
    feats[55] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', '10_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_sp', node2.lemma)
    
    # L -> 12_fp -> -123 -> L
    feats[56] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_fp', '-123') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-123', node2.lemma)
    
    # L -> 77 -> 58 -> L
    feats[57] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '77') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '77', '58') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '58', node2.lemma)
    
    # L -> 77 -> 10_sp -> L
    feats[58] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '77') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '77', '10_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_sp', node2.lemma)
    
    # L -> nom. du. -> nom. adj. -> L
    feats[59] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. du.', 'nom. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. adj.', node2.lemma)
    
    # L -> -30 -> -52 -> L
    feats[60] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '-52') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-52', node2.lemma)
    
    # L -> -306 -> -169 -> L
    feats[61] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '-169') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-169', node2.lemma)
    
    # L -> -306 -> -87 -> L
    feats[62] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '-87') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-87', node2.lemma)
    
    # L -> 35 -> 10_sg -> L
    feats[63] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', '10_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_sg', node2.lemma)
    
    # L -> -306 -> -53 -> L
    feats[64] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '-53') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-53', node2.lemma)
    
    # L -> -306 -> -101 -> L
    feats[65] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '-101') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-101', node2.lemma)
    
    # L -> nom. du. -> -39 -> L
    feats[66] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. du.', '-39') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-39', node2.lemma)
    
    # L -> -306 -> 2_sp -> L
    feats[67] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '2_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_sp', node2.lemma)
    
    # L -> -306 -> 16_du -> L
    feats[68] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '16_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_du', node2.lemma)
    
    # L -> -306 -> -17 -> L
    feats[69] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '-17') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-17', node2.lemma)
    
    # L -> -306 -> -68 -> L
    feats[70] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '-68') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-68', node2.lemma)
    
    # L -> -30 -> -302 -> L
    feats[71] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '-302') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-302', node2.lemma)
    
    # L -> -30 -> -157 -> L
    feats[72] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '-157') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-157', node2.lemma)
    
    # L -> -306 -> 15_tp -> L
    feats[73] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '15_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_tp', node2.lemma)
    
    # L -> -306 -> 10_pl -> L
    feats[74] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '10_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_pl', node2.lemma)
    
    # L -> -306 -> -151 -> L
    feats[75] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '-151') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-151', node2.lemma)
    
    # L -> -306 -> -262 -> L
    feats[76] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '-262') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-262', node2.lemma)
    
    # L -> -306 -> 35 -> L
    feats[77] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '35') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '35', node2.lemma)
    
    # L -> -306 -> 3_tp -> L
    feats[78] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '3_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_tp', node2.lemma)
    
    # L -> -306 -> -302 -> L
    feats[79] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '-302') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-302', node2.lemma)
    
    # L -> -306 -> -157 -> L
    feats[80] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '-157') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-157', node2.lemma)
    
    # L -> -306 -> nom. neutr. -> L
    feats[81] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', 'nom. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. neutr.', node2.lemma)
    
    # L -> -306 -> 10_tp -> L
    feats[82] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '10_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_tp', node2.lemma)
    
    # L -> -306 -> 95 -> L
    feats[83] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '95') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '95', node2.lemma)
    
    # L -> -306 -> -23 -> L
    feats[84] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '-23') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-23', node2.lemma)
    
    # L -> -306 -> 70 -> L
    feats[85] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '70') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '70', node2.lemma)
    
    # L -> -306 -> -33 -> L
    feats[86] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '-33') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-33', node2.lemma)
    
    # L -> -306 -> 28 -> L
    feats[87] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '28') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
    
    # L -> -306 -> -144 -> L
    feats[88] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '-144') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-144', node2.lemma)
    
    # L -> -306 -> 77 -> L
    feats[89] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '77') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '77', node2.lemma)
    
    # L -> -306 -> -73 -> L
    feats[90] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '-73') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-73', node2.lemma)
    
    # L -> -306 -> -96 -> L
    feats[91] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '-96') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-96', node2.lemma)
    
    # L -> 35 -> -25 -> L
    feats[92] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', '-25') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-25', node2.lemma)
    
    # L -> -306 -> -147 -> L
    feats[93] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '-147') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-147', node2.lemma)
    
    # L -> -306 -> 140 -> L
    feats[94] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '140') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '140', node2.lemma)
    
    # L -> -306 -> 15_du -> L
    feats[95] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '15_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_du', node2.lemma)
    
    # L -> -306 -> -32 -> L
    feats[96] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '-32') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-32', node2.lemma)
    
    # L -> -306 -> pl -> L
    feats[97] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', 'pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    
    # L -> -306 -> 1 -> L
    feats[98] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '1') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '1', node2.lemma)
    
    # L -> -306 -> -19 -> L
    feats[99] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '-19') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-19', node2.lemma)
    
    # L -> -306 -> 68 -> L
    feats[100] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '68') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '68', node2.lemma)
    
    # L -> -306 -> 9_fp -> L
    feats[101] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '9_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_fp', node2.lemma)
    
    # L -> nom. du. -> 94 -> L
    feats[102] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. du.', '94') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '94', node2.lemma)
    
    # L -> -306 -> 58 -> L
    feats[103] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '58') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '58', node2.lemma)
    
    # L -> -306 -> -161 -> L
    feats[104] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '-161') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-161', node2.lemma)
    
    # L -> -306 -> 50 -> L
    feats[105] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '50') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '50', node2.lemma)
    
    # L -> -306 -> 182 -> L
    feats[106] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '182') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '182', node2.lemma)
    
    # L -> -306 -> 118 -> L
    feats[107] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '118') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '118', node2.lemma)
    
    # L -> -306 -> -109 -> L
    feats[108] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '-109') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-109', node2.lemma)
    
    # L -> -30 -> 150 -> L
    feats[109] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '150') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '150', node2.lemma)
    
    # L -> -306 -> 30_pl -> L
    feats[110] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '30_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_pl', node2.lemma)
    
    # L -> -306 -> -247 -> L
    feats[111] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '-247') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-247', node2.lemma)
    
    # L -> -306 -> 9_du -> L
    feats[112] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '9_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_du', node2.lemma)
    
    # L -> -306 -> acc. pl. -> L
    feats[113] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', 'acc. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. pl.', node2.lemma)
    
    # L -> -306 -> 6_pl -> L
    feats[114] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '6_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_pl', node2.lemma)
    
    # L -> -306 -> 114 -> L
    feats[115] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '114') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '114', node2.lemma)
    
    # L -> nom. adj. -> 12_sg -> L
    feats[116] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. adj.', '12_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_sg', node2.lemma)
    
    # L -> -306 -> 9_tp -> L
    feats[117] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '9_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_tp', node2.lemma)
    
    # L -> -306 -> -41 -> L
    feats[118] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '-41') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-41', node2.lemma)
    
    # L -> -306 -> nom. masc. -> L
    feats[119] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', 'nom. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. masc.', node2.lemma)
    
    # L -> -306 -> 153 -> L
    feats[120] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '153') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    
    # L -> -306 -> 150 -> L
    feats[121] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '150') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '150', node2.lemma)
    
    # L -> -306 -> 98 -> L
    feats[122] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '98') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '98', node2.lemma)
    
    # L -> -24 -> 169 -> L
    feats[123] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-24') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-24', '169') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '169', node2.lemma)
    
    # L -> -306 -> 30_tp -> L
    feats[124] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '30_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_tp', node2.lemma)
    
    # L -> -306 -> sg_fp -> L
    feats[125] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', 'sg_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_fp', node2.lemma)
    
    # L -> -306 -> 14_tp -> L
    feats[126] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '14_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_tp', node2.lemma)
    
    # L -> -306 -> -102 -> L
    feats[127] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '-102') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-102', node2.lemma)
    
    # L -> -306 -> -159 -> L
    feats[128] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '-159') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-159', node2.lemma)
    
    # L -> -306 -> 5_pl -> L
    feats[129] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', '5_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_pl', node2.lemma)
    
    # L -> 3_tp -> -169 -> L
    feats[130] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_tp', '-169') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-169', node2.lemma)
    
    # L -> nom. du. -> -266 -> L
    feats[131] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. du.', '-266') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-266', node2.lemma)
    
    # L -> nom. du. -> -94 -> L
    feats[132] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. du.', '-94') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-94', node2.lemma)
    
    # L -> -24 -> 130 -> L
    feats[133] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-24') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-24', '130') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '130', node2.lemma)
    
    # L -> -24 -> 27_pl -> L
    feats[134] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-24') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-24', '27_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_pl', node2.lemma)
    
    # L -> -24 -> 11_tp -> L
    feats[135] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-24') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-24', '11_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_tp', node2.lemma)
    
    # L -> -30 -> -126 -> L
    feats[136] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '-126') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-126', node2.lemma)
    
    # L -> -30 -> 82 -> L
    feats[137] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', '82') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '82', node2.lemma)
    
    # L -> nom. du. -> 119 -> L
    feats[138] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. du.', '119') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '119', node2.lemma)
    
    # L -> nom. du. -> dat. sg. -> L
    feats[139] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. du.', 'dat. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. sg.', node2.lemma)
    
    # L -> -24 -> 13_sp -> L
    feats[140] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-24') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-24', '13_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_sp', node2.lemma)
    
    # L -> nom. adj. -> du -> L
    feats[141] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. adj.', 'du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du', node2.lemma)
    
    # L -> nom. du. -> 139 -> L
    feats[142] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. du.', '139') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '139', node2.lemma)
    
    # L -> nom. adj. -> 10_pl -> L
    feats[143] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. adj.', '10_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_pl', node2.lemma)
    
    # L -> nom. du. -> 8_tp -> L
    feats[144] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. du.', '8_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_tp', node2.lemma)
    
    # L -> -24 -> -45 -> L
    feats[145] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-24') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-24', '-45') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-45', node2.lemma)
    
    # L -> nom. du. -> -37 -> L
    feats[146] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. du.', '-37') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-37', node2.lemma)
    
    # L -> -24 -> -68 -> L
    feats[147] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-24') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-24', '-68') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-68', node2.lemma)
    
    # L -> nom. du. -> -46 -> L
    feats[148] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. du.', '-46') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-46', node2.lemma)
    
    # L -> -24 -> 10_pl -> L
    feats[149] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-24') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-24', '10_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_pl', node2.lemma)
    
    # L -> nom. du. -> 51 -> L
    feats[150] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. du.', '51') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '51', node2.lemma)
    
    # L -> -24 -> -15 -> L
    feats[151] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-24') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-24', '-15') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-15', node2.lemma)
    
    # L -> -24 -> -137 -> L
    feats[152] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-24') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-24', '-137') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-137', node2.lemma)
    
    # L -> nom. du. -> -156 -> L
    feats[153] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. du.', '-156') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-156', node2.lemma)
    
    # L -> nom. du. -> 3 -> L
    feats[154] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. du.', '3') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3', node2.lemma)
    
    # L -> nom. du. -> -104 -> L
    feats[155] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. du.', '-104') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-104', node2.lemma)
    
    # L -> nom. du. -> 102 -> L
    feats[156] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. du.', '102') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '102', node2.lemma)
    
    # L -> nom. du. -> 11_fp -> L
    feats[157] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. du.', '11_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_fp', node2.lemma)
    
    # L -> nom. du. -> 55 -> L
    feats[158] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. du.', '55') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '55', node2.lemma)
    
    # L -> 3_tp -> 12_fp -> L
    feats[159] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_tp', '12_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_fp', node2.lemma)
    
    # L -> nom. du. -> 12_sp -> L
    feats[160] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. du.', '12_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_sp', node2.lemma)
    
    # L -> -24 -> 137 -> L
    feats[161] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-24') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-24', '137') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    
    # L -> -24 -> 14_sg -> L
    feats[162] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-24') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-24', '14_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_sg', node2.lemma)
    
    # L -> -24 -> -152 -> L
    feats[163] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-24') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-24', '-152') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-152', node2.lemma)
    
    # L -> nom. du. -> instr. du. -> L
    feats[164] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. du.', 'instr. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. du.', node2.lemma)
    
    # L -> nom. du. -> 7_fp -> L
    feats[165] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. du.', '7_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_fp', node2.lemma)
    
    # L -> nom. du. -> -293 -> L
    feats[166] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. du.', '-293') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-293', node2.lemma)
    
    # L -> -24 -> 28 -> L
    feats[167] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-24') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-24', '28') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
    
    # L -> nom. du. -> 161 -> L
    feats[168] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. du.', '161') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '161', node2.lemma)
    
    # L -> nom. du. -> 160 -> L
    feats[169] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. du.', '160') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    
    # L -> nom. du. -> 73 -> L
    feats[170] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. du.', '73') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '73', node2.lemma)
    
    # L -> -24 -> -96 -> L
    feats[171] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-24') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-24', '-96') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-96', node2.lemma)
    
    # L -> -24 -> -240 -> L
    feats[172] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-24') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-24', '-240') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-240', node2.lemma)
    
    # L -> nom. du. -> fp -> L
    feats[173] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. du.', 'fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'fp', node2.lemma)
    
    # L -> -24 -> acc. masc. -> L
    feats[174] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-24') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-24', 'acc. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. masc.', node2.lemma)
    
    # L -> nom. du. -> -123 -> L
    feats[175] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. du.', '-123') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-123', node2.lemma)
    
    # L -> nom. du. -> -210 -> L
    feats[176] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. du.', '-210') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-210', node2.lemma)
    
    # L -> nom. du. -> voc. fem -> L
    feats[177] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. du.', 'voc. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. fem', node2.lemma)
    
    # L -> -24 -> 121 -> L
    feats[178] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-24') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-24', '121') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '121', node2.lemma)
    
    # L -> -24 -> pl -> L
    feats[179] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-24') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-24', 'pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    
    # L -> -24 -> 1 -> L
    feats[180] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-24') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-24', '1') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '1', node2.lemma)
    
    # L -> 31 -> 79 -> L
    feats[181] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '31') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '31', '79') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '79', node2.lemma)
    
    # L -> 31 -> sg_tp -> L
    feats[182] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '31') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '31', 'sg_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_tp', node2.lemma)
    
    # L -> -93 -> voc. sg. -> L
    feats[183] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-93') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-93', 'voc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. sg.', node2.lemma)
    
    # L -> 31 -> 120 -> L
    feats[184] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '31') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '31', '120') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '120', node2.lemma)
    
    # L->-190->C
    feats[185] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-190') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-190', node2.cng)
            
    # L -> 31 -> -93 -> L
    feats[186] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '31') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '31', '-93') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-93', node2.lemma)
    
    # L -> nom. adj. -> 13_fp -> L
    feats[187] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. adj.', '13_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_fp', node2.lemma)
    
    # L -> 31 -> -51 -> L
    feats[188] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '31') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '31', '-51') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-51', node2.lemma)
    
    # L -> 31 -> 27_fp -> L
    feats[189] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '31') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '31', '27_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_fp', node2.lemma)
    
    # L->8_pl->C
    feats[190] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_pl', node2.cng)
            
    # L -> 31 -> 5_fp -> L
    feats[191] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '31') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '31', '5_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_fp', node2.lemma)
    
    # L -> 31 -> -36 -> L
    feats[192] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '31') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '31', '-36') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-36', node2.lemma)
    
    # L -> 31 -> 155 -> L
    feats[193] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '31') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '31', '155') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '155', node2.lemma)
    
    # L -> 31 -> 156 -> L
    feats[194] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '31') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '31', '156') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '156', node2.lemma)
    
    # L->-129->C
    feats[195] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-129') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-129', node2.cng)
            
    # L -> 31 -> -13 -> L
    feats[196] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '31') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '31', '-13') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-13', node2.lemma)
    
    # L -> 31 -> 12_sg -> L
    feats[197] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '31') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '31', '12_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_sg', node2.lemma)
    
    # L->loc. pl.->C
    feats[198] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. pl.', node2.cng)
            
    # L -> 31 -> -276 -> L
    feats[199] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '31') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '31', '-276') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-276', node2.lemma)
    
    # L -> 31 -> du_fp -> L
    feats[200] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '31') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '31', 'du_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_fp', node2.lemma)
    
    # L -> 31 -> acc. neutr. -> L
    feats[201] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '31') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '31', 'acc. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. neutr.', node2.lemma)
    
    # L -> 31 -> 5_du -> L
    feats[202] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '31') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '31', '5_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_du', node2.lemma)
    
    # L -> 31 -> -52 -> L
    feats[203] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '31') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '31', '-52') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-52', node2.lemma)
    
    # L -> 31 -> 116 -> L
    feats[204] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '31') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '31', '116') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '116', node2.lemma)
    
    # L -> 31 -> -69 -> L
    feats[205] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '31') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '31', '-69') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-69', node2.lemma)
    
    # L -> 31 -> 88 -> L
    feats[206] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '31') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '31', '88') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '88', node2.lemma)
    
    # L -> 31 -> 130 -> L
    feats[207] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '31') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '31', '130') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '130', node2.lemma)
    
    # L -> 31 -> 12_pl -> L
    feats[208] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '31') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '31', '12_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_pl', node2.lemma)
    
    # L -> 31 -> -78 -> L
    feats[209] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '31') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '31', '-78') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-78', node2.lemma)
    
    # L -> 31 -> 90 -> L
    feats[210] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '31') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '31', '90') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '90', node2.lemma)
    
    # L -> -93 -> -72 -> L
    feats[211] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-93') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-93', '-72') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-72', node2.lemma)
    
    # L -> -93 -> 2_sp -> L
    feats[212] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-93') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-93', '2_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_sp', node2.lemma)
    
    # L -> nom. adj. -> 8_du -> L
    feats[213] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. adj.', '8_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_du', node2.lemma)
    
    # L -> -291 -> 132 -> L
    feats[214] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-291') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-291', '132') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '132', node2.lemma)
    
    # L -> -24 -> 5_sp -> L
    feats[215] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-24') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-24', '5_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sp', node2.lemma)
    
    # L -> -291 -> -58 -> L
    feats[216] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-291') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-291', '-58') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-58', node2.lemma)
    
    # L -> 3_tp -> -77 -> L
    feats[217] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_tp', '-77') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-77', node2.lemma)
    
    # L -> -291 -> -230 -> L
    feats[218] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-291') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-291', '-230') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-230', node2.lemma)
    
    # L -> -291 -> -308 -> L
    feats[219] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-291') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-291', '-308') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-308', node2.lemma)
    
    # L->-142->C
    feats[220] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-142') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-142', node2.cng)
            
    # L->-67->C
    feats[221] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', node2.cng)
            
    # L -> -24 -> -86 -> L
    feats[222] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-24') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-24', '-86') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-86', node2.lemma)
    
    # L -> -93 -> 50 -> L
    feats[223] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-93') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-93', '50') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '50', node2.lemma)
    
    # L->27_du->T
    feats[224] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '27_du', node2.tup)
            
    # L -> nom. adj. -> acc -> L
    feats[225] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. adj.', 'acc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc', node2.lemma)
    
    # L->voc. neutr.->T
    feats[226] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. neutr.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'voc. neutr.', node2.tup)
            
    # L -> 3_tp -> acc -> L
    feats[227] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_tp', 'acc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc', node2.lemma)
    
    # L -> 3_tp -> 177 -> L
    feats[228] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_tp', '177') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '177', node2.lemma)
    
    # L->30_du->T
    feats[229] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_du') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '30_du', node2.tup)
            
    # L->-81->T
    feats[230] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-81', node2.tup)
            
    # L->149->T
    feats[231] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '149', node2.tup)
            
    # L -> -73 -> -93 -> L
    feats[232] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-73') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-73', '-93') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-93', node2.lemma)
    
    # L -> -73 -> 78 -> L
    feats[233] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-73') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-73', '78') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '78', node2.lemma)
    
    # L->-263->T
    feats[234] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-263') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-263', node2.tup)
            
    # L -> -73 -> 155 -> L
    feats[235] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-73') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-73', '155') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '155', node2.lemma)
    
    # L->173->T
    feats[236] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '173') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '173', node2.tup)
            
    # L -> -291 -> gen. du. -> L
    feats[237] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-291') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-291', 'gen. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen. du.', node2.lemma)
    
    # L->27_pl->T
    feats[238] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_pl') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '27_pl', node2.tup)
            
    # L -> 30_du -> instr. adj. -> L
    feats[239] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_du', 'instr. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. adj.', node2.lemma)
    
    # L -> -93 -> voc. pl. -> L
    feats[240] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-93') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-93', 'voc. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. pl.', node2.lemma)
    
    # L -> -93 -> -269 -> L
    feats[241] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-93') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-93', '-269') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-269', node2.lemma)
    
    # L -> -73 -> -111 -> L
    feats[242] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-73') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-73', '-111') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-111', node2.lemma)
    
    # L -> -73 -> du_fp -> L
    feats[243] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-73') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-73', 'du_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_fp', node2.lemma)
    
    # L->2_fp->T
    feats[244] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '2_fp', node2.tup)
            
    # L -> 30_du -> 74 -> L
    feats[245] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_du', '74') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '74', node2.lemma)
    
    # L -> 30_du -> 42 -> L
    feats[246] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_du', '42') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '42', node2.lemma)
    
    # L -> -73 -> 131 -> L
    feats[247] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-73') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-73', '131') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '131', node2.lemma)
    
    # L -> -73 -> 173 -> L
    feats[248] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-73') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-73', '173') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '173', node2.lemma)
    
    # L -> 31 -> -83 -> L
    feats[249] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '31') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '31', '-83') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-83', node2.lemma)
    
    # L->-57->T
    feats[250] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-57') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-57', node2.tup)
            
    # L->2_sp->T
    feats[251] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_sp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '2_sp', node2.tup)
            
    # L -> -93 -> -35 -> L
    feats[252] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-93') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-93', '-35') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    
    # L->9_pl->T
    feats[253] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_pl') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '9_pl', node2.tup)
            
    # L -> -93 -> dat. sg. -> L
    feats[254] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-93') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-93', 'dat. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. sg.', node2.lemma)
    
    # L->10_pl->T
    feats[255] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_pl') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '10_pl', node2.tup)
            
    # L->10_du->T
    feats[256] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_du') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '10_du', node2.tup)
            
    # L -> -73 -> -241 -> L
    feats[257] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-73') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-73', '-241') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-241', node2.lemma)
    
    # L -> -73 -> 13_sp -> L
    feats[258] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-73') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-73', '13_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_sp', node2.lemma)
    
    # L -> 11_pl -> 34 -> L
    feats[259] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_pl', '34') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '34', node2.lemma)
    
    # L -> -73 -> -57 -> L
    feats[260] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-73') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-73', '-57') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-57', node2.lemma)
    
    # L->35->T
    feats[261] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '35', node2.tup)
            
    # L->10_tp->T
    feats[262] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_tp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '10_tp', node2.tup)
            
    # L->-152->T
    feats[263] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-152') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-152', node2.tup)
            
    # L -> 11_pl -> abl. sg. -> L
    feats[264] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_pl', 'abl. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. sg.', node2.lemma)
    
    # L -> -93 -> 162 -> L
    feats[265] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-93') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-93', '162') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '162', node2.lemma)
    
    # L -> -93 -> -79 -> L
    feats[266] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-93') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-93', '-79') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-79', node2.lemma)
    
    # L -> -291 -> -307 -> L
    feats[267] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-291') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-291', '-307') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-307', node2.lemma)
    
    # L -> -73 -> -99 -> L
    feats[268] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-73') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-73', '-99') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-99', node2.lemma)
    
    # L -> -291 -> 122 -> L
    feats[269] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-291') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-291', '122') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '122', node2.lemma)
    
    # L -> -73 -> -157 -> L
    feats[270] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-73') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-73', '-157') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-157', node2.lemma)
    
    # L -> -73 -> 16_tp -> L
    feats[271] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-73') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-73', '16_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_tp', node2.lemma)
    
    # L -> 31 -> 152 -> L
    feats[272] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '31') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '31', '152') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '152', node2.lemma)
    
    # L -> 31 -> 55 -> L
    feats[273] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '31') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '31', '55') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '55', node2.lemma)
    
    # L -> -93 -> 2_pl -> L
    feats[274] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-93') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-93', '2_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_pl', node2.lemma)
    
    # L -> -93 -> 49 -> L
    feats[275] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-93') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-93', '49') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '49', node2.lemma)
    
    # L -> 30_du -> pl -> L
    feats[276] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_du', 'pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    
    # L -> -93 -> 48 -> L
    feats[277] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-93') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-93', '48') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '48', node2.lemma)
    
    # L -> 11_pl -> -32 -> L
    feats[278] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_pl', '-32') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-32', node2.lemma)
    
    # L -> voc. sg. -> -30 -> L
    feats[279] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. sg.', '-30') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-30', node2.lemma)
    
    # L -> voc. sg. -> -13 -> L
    feats[280] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. sg.', '-13') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-13', node2.lemma)
    
    # L -> voc. sg. -> 5_du -> L
    feats[281] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. sg.', '5_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_du', node2.lemma)
    
    # L -> voc. sg. -> sp -> L
    feats[282] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. sg.', 'sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sp', node2.lemma)
    
    # L -> -73 -> -154 -> L
    feats[283] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-73') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-73', '-154') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-154', node2.lemma)
    
    # L -> 30_du -> dat. du. -> L
    feats[284] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_du', 'dat. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    
    # L -> 30_du -> 36 -> L
    feats[285] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_du', '36') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '36', node2.lemma)
    
    # L->-34->T
    feats[286] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-34', node2.tup)
            
    # L->16_sg->T
    feats[287] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_sg') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '16_sg', node2.tup)
            
    # L -> 11_pl -> 32 -> L
    feats[288] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_pl', '32') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '32', node2.lemma)
    
    # L -> -97 -> 16_du -> L
    feats[289] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-97') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-97', '16_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_du', node2.lemma)
    
    # L -> 70 -> 72 -> L
    feats[290] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '70') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '70', '72') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '72', node2.lemma)
    
    # L->dat. sg.->T
    feats[291] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. sg.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'dat. sg.', node2.tup)
            
    # L->82->T
    feats[292] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '82') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '82', node2.tup)
            
    # L->139->T
    feats[293] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '139') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '139', node2.tup)
            
    # L -> -73 -> 10_sg -> L
    feats[294] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-73') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-73', '10_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_sg', node2.lemma)
    
    # L->15_pl->T
    feats[295] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_pl') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '15_pl', node2.tup)
            
    # L -> -73 -> 4_sg -> L
    feats[296] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-73') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-73', '4_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_sg', node2.lemma)
    
    # L -> -73 -> pl_tp -> L
    feats[297] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-73') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-73', 'pl_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_tp', node2.lemma)
    
    # L->-141->T
    feats[298] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-141') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-141', node2.tup)
            
    # L->-79->T
    feats[299] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-79') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-79', node2.tup)
            
    # L -> voc. sg. -> gen -> L
    feats[300] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. sg.', 'gen') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen', node2.lemma)
    
    # L->55->T
    feats[301] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '55') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '55', node2.tup)
            
    # L->-150->T
    feats[302] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-150') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-150', node2.tup)
            
    # L -> 30_du -> 152 -> L
    feats[303] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_du', '152') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '152', node2.lemma)
    
    # L -> 70 -> 38 -> L
    feats[304] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '70') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '70', '38') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '38', node2.lemma)
    
    # L->gen. sg.->T
    feats[305] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. sg.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'gen. sg.', node2.tup)
            
    # L->-121->T
    feats[306] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-121') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-121', node2.tup)
            
    # L->49->T
    feats[307] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '49') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '49', node2.tup)
            
    # L -> 11_pl -> -104 -> L
    feats[308] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_pl', '-104') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-104', node2.lemma)
    
    # L -> 70 -> 6_fp -> L
    feats[309] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '70') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '70', '6_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_fp', node2.lemma)
    
    # C->99->L
    feats[310] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '99') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '99', node2.lemma)
            
    # C->27_du->L
    feats[311] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '27_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_du', node2.lemma)
            
    # C->7_sg->L
    feats[312] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '7_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_sg', node2.lemma)
            
    # C->7_sp->L
    feats[313] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '7_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_sp', node2.lemma)
            
    # C->-10->L
    feats[314] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-10') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-10', node2.lemma)
            
    # C->-26->L
    feats[315] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-26') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-26', node2.lemma)
            
    # L -> 70 -> 10_fp -> L
    feats[316] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '70') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '70', '10_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_fp', node2.lemma)
    
    # L -> 70 -> 89 -> L
    feats[317] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '70') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '70', '89') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '89', node2.lemma)
    
    # L -> 110 -> 99 -> L
    feats[318] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '110') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '110', '99') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '99', node2.lemma)
    
    # L -> 70 -> 94 -> L
    feats[319] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '70') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '70', '94') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '94', node2.lemma)
    
    # L -> voc. sg. -> 3_sg -> L
    feats[320] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. sg.', '3_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_sg', node2.lemma)
    
    # L -> 70 -> pl_sp -> L
    feats[321] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '70') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '70', 'pl_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_sp', node2.lemma)
    
    # L -> 2_du -> -98 -> L
    feats[322] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '-98') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-98', node2.lemma)
    
    # L -> 2_du -> 7_pl -> L
    feats[323] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '7_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_pl', node2.lemma)
    
    # L -> 70 -> 3_pl -> L
    feats[324] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '70') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '70', '3_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_pl', node2.lemma)
    
    # L -> 2_du -> sp -> L
    feats[325] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', 'sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sp', node2.lemma)
    
    # L -> 2_du -> -52 -> L
    feats[326] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '-52') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-52', node2.lemma)
    
    # L -> 70 -> -47 -> L
    feats[327] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '70') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '70', '-47') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-47', node2.lemma)
    
    # L -> 2_du -> 116 -> L
    feats[328] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '116') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '116', node2.lemma)
    
    # L -> 2_du -> 30_sp -> L
    feats[329] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '30_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_sp', node2.lemma)
    
    # L -> 2_du -> -82 -> L
    feats[330] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '-82') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-82', node2.lemma)
    
    # L -> 2_du -> 16_pl -> L
    feats[331] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '16_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_pl', node2.lemma)
    
    # L -> 2_du -> 14_pl -> L
    feats[332] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '14_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_pl', node2.lemma)
    
    # L -> -261 -> tp -> L
    feats[333] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-261') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-261', 'tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'tp', node2.lemma)
    
    # L -> 70 -> sg_fp -> L
    feats[334] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '70') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '70', 'sg_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_fp', node2.lemma)
    
    # L -> 2_du -> 157 -> L
    feats[335] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '157') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '157', node2.lemma)
    
    # L -> 2_du -> -169 -> L
    feats[336] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '-169') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-169', node2.lemma)
    
    # L -> 2_du -> 130 -> L
    feats[337] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '130') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '130', node2.lemma)
    
    # L -> -268 -> du_fp -> L
    feats[338] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', 'du_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_fp', node2.lemma)
    
    # L -> 2_du -> 132 -> L
    feats[339] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '132') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '132', node2.lemma)
    
    # L -> 2_du -> -90 -> L
    feats[340] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '-90') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-90', node2.lemma)
    
    # L -> 70 -> 2 -> L
    feats[341] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '70') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '70', '2') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2', node2.lemma)
    
    # L -> 2_du -> 27_pl -> L
    feats[342] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '27_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_pl', node2.lemma)
    
    # L -> 2_du -> -18 -> L
    feats[343] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '-18') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-18', node2.lemma)
    
    # L -> 2_du -> 90 -> L
    feats[344] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '90') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '90', node2.lemma)
    
    # L -> 2_du -> -87 -> L
    feats[345] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '-87') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-87', node2.lemma)
    
    # L -> 2_du -> 11_tp -> L
    feats[346] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '11_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_tp', node2.lemma)
    
    # L -> 2_du -> 74 -> L
    feats[347] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '74') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '74', node2.lemma)
    
    # L -> 70 -> 40 -> L
    feats[348] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '70') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '70', '40') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '40', node2.lemma)
    
    # L -> 2_du -> -44 -> L
    feats[349] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '-44') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-44', node2.lemma)
    
    # L -> 2_du -> voc. du. -> L
    feats[350] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', 'voc. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. du.', node2.lemma)
    
    # L -> 2_du -> du_tp -> L
    feats[351] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', 'du_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_tp', node2.lemma)
    
    # L -> 70 -> 36 -> L
    feats[352] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '70') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '70', '36') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '36', node2.lemma)
    
    # L -> 2_du -> -84 -> L
    feats[353] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '-84') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-84', node2.lemma)
    
    # L -> 2_du -> -301 -> L
    feats[354] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '-301') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-301', node2.lemma)
    
    # L -> 2_du -> -58 -> L
    feats[355] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '-58') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-58', node2.lemma)
    
    # L -> 70 -> -22 -> L
    feats[356] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '70') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '70', '-22') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-22', node2.lemma)
    
    # L -> 2_du -> -241 -> L
    feats[357] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '-241') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-241', node2.lemma)
    
    # L -> 70 -> 141 -> L
    feats[358] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '70') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '70', '141') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '141', node2.lemma)
    
    # L -> 70 -> 154 -> L
    feats[359] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '70') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '70', '154') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '154', node2.lemma)
    
    # L -> 2_du -> -72 -> L
    feats[360] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '-72') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-72', node2.lemma)
    
    # L -> 2_du -> -45 -> L
    feats[361] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '-45') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-45', node2.lemma)
    
    # L -> 2_du -> 34 -> L
    feats[362] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '34') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '34', node2.lemma)
    
    # L -> 2_du -> 60 -> L
    feats[363] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '60') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '60', node2.lemma)
    
    # L -> 2_du -> -17 -> L
    feats[364] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '-17') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-17', node2.lemma)
    
    # L -> 2_du -> -292 -> L
    feats[365] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '-292') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-292', node2.lemma)
    
    # L -> 110 -> -18 -> L
    feats[366] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '110') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '110', '-18') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-18', node2.lemma)
    
    # L -> 110 -> 90 -> L
    feats[367] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '110') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '110', '90') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '90', node2.lemma)
    
    # L -> 2_du -> -143 -> L
    feats[368] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '-143') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-143', node2.lemma)
    
    # L -> 2_du -> 59 -> L
    feats[369] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '59') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '59', node2.lemma)
    
    # L -> 2_du -> -68 -> L
    feats[370] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '-68') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-68', node2.lemma)
    
    # L -> 2_du -> -39 -> L
    feats[371] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '-39') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-39', node2.lemma)
    
    # L -> 2_du -> -122 -> L
    feats[372] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '-122') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-122', node2.lemma)
    
    # L -> 2_du -> 15_tp -> L
    feats[373] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '15_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_tp', node2.lemma)
    
    # L -> 2_du -> -91 -> L
    feats[374] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '-91') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-91', node2.lemma)
    
    # L -> 2_du -> 9_pl -> L
    feats[375] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '9_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_pl', node2.lemma)
    
    # L -> 2_du -> -279 -> L
    feats[376] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '-279') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-279', node2.lemma)
    
    # L -> 2_du -> 101 -> L
    feats[377] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '101') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    
    # L -> 2_du -> 109 -> L
    feats[378] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '109') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '109', node2.lemma)
    
    # L -> 2_du -> du_sp -> L
    feats[379] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', 'du_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_sp', node2.lemma)
    
    # L -> 70 -> voc. masc. -> L
    feats[380] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '70') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '70', 'voc. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. masc.', node2.lemma)
    
    # L -> 2_du -> -99 -> L
    feats[381] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '-99') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-99', node2.lemma)
    
    # C->117->L
    feats[382] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '117') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '117', node2.lemma)
            
    # L -> 2_du -> acc. adj. -> L
    feats[383] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', 'acc. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. adj.', node2.lemma)
    
    # L -> 2_du -> du -> L
    feats[384] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', 'du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du', node2.lemma)
    
    # L -> -261 -> 8_tp -> L
    feats[385] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-261') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-261', '8_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_tp', node2.lemma)
    
    # L -> 2_du -> 10_du -> L
    feats[386] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '10_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_du', node2.lemma)
    
    # L -> 70 -> 93 -> L
    feats[387] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '70') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '70', '93') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '93', node2.lemma)
    
    # L -> 2_du -> -31 -> L
    feats[388] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '-31') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-31', node2.lemma)
    
    # L -> 2_du -> -137 -> L
    feats[389] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '-137') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-137', node2.lemma)
    
    # L -> 2_du -> -151 -> L
    feats[390] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '-151') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-151', node2.lemma)
    
    # L -> 2_du -> -29 -> L
    feats[391] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '-29') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-29', node2.lemma)
    
    # L -> 70 -> 4_sg -> L
    feats[392] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '70') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '70', '4_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_sg', node2.lemma)
    
    # L -> 2_du -> 3_tp -> L
    feats[393] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '3_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_tp', node2.lemma)
    
    # L -> 2_du -> -245 -> L
    feats[394] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '-245') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-245', node2.lemma)
    
    # L -> 70 -> -273 -> L
    feats[395] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '70') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '70', '-273') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-273', node2.lemma)
    
    # L -> 2_du -> -157 -> L
    feats[396] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '-157') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-157', node2.lemma)
    
    # L -> 2_du -> 16_tp -> L
    feats[397] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '16_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_tp', node2.lemma)
    
    # L -> 70 -> dat. sg. -> L
    feats[398] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '70') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '70', 'dat. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. sg.', node2.lemma)
    
    # L -> -261 -> abl. du. -> L
    feats[399] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-261') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-261', 'abl. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. du.', node2.lemma)
    
    # L -> -261 -> -46 -> L
    feats[400] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-261') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-261', '-46') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-46', node2.lemma)
    
    # L -> 2_du -> -152 -> L
    feats[401] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '-152') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-152', node2.lemma)
    
    # L -> 2_du -> 14_sg -> L
    feats[402] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '14_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_sg', node2.lemma)
    
    # L -> 2_du -> 38 -> L
    feats[403] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '38') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '38', node2.lemma)
    
    # L -> 2_du -> 137 -> L
    feats[404] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '137') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    
    # L -> 2_du -> 37 -> L
    feats[405] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '37') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '37', node2.lemma)
    
    # L -> 2_du -> 95 -> L
    feats[406] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '95') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '95', node2.lemma)
    
    # L -> 2_du -> 10_tp -> L
    feats[407] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '10_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_tp', node2.lemma)
    
    # L -> 2_du -> 12_fp -> L
    feats[408] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '12_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_fp', node2.lemma)
    
    # L -> 2_du -> -306 -> L
    feats[409] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '-306') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-306', node2.lemma)
    
    # L -> 2_du -> -291 -> L
    feats[410] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '-291') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-291', node2.lemma)
    
    # L -> 2_du -> 70 -> L
    feats[411] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '70') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '70', node2.lemma)
    
    # L -> -268 -> 3_tp -> L
    feats[412] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', '3_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_tp', node2.lemma)
    
    # L -> 2_du -> -48 -> L
    feats[413] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '-48') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-48', node2.lemma)
    
    # L -> 2_du -> -230 -> L
    feats[414] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '-230') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-230', node2.lemma)
    
    # L -> 2_du -> -308 -> L
    feats[415] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '-308') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-308', node2.lemma)
    
    # L -> 2_du -> voc -> L
    feats[416] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', 'voc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc', node2.lemma)
    
    # L -> -268 -> 11_pl -> L
    feats[417] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', '11_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_pl', node2.lemma)
    
    # C->-21->L
    feats[418] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-21') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-21', node2.lemma)
            
    # C->118->L
    feats[419] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '118') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '118', node2.lemma)
            
    # L -> voc. neutr. -> voc. neutr. -> L
    feats[420] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. neutr.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. neutr.', 'voc. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. neutr.', node2.lemma)
    
    # L -> voc. neutr. -> 27_fp -> L
    feats[421] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. neutr.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. neutr.', '27_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_fp', node2.lemma)
    
    # L -> gen. pl. -> 79 -> L
    feats[422] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. pl.', '79') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '79', node2.lemma)
    
    # L -> 2_du -> -309 -> L
    feats[423] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '-309') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-309', node2.lemma)
    
    # L -> 70 -> -12 -> L
    feats[424] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '70') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '70', '-12') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-12', node2.lemma)
    
    # L -> 30_fp -> -72 -> L
    feats[425] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_fp', '-72') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-72', node2.lemma)
    
    # L -> gen. pl. -> -87 -> L
    feats[426] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. pl.', '-87') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-87', node2.lemma)
    
    # L -> voc. neutr. -> -279 -> L
    feats[427] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. neutr.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. neutr.', '-279') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-279', node2.lemma)
    
    # L -> 110 -> instr. sg. -> L
    feats[428] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '110') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '110', 'instr. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. sg.', node2.lemma)
    
    # L -> 110 -> -66 -> L
    feats[429] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '110') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '110', '-66') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-66', node2.lemma)
    
    # L -> gen. pl. -> 60 -> L
    feats[430] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. pl.', '60') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '60', node2.lemma)
    
    # L -> 2_du -> 15_pl -> L
    feats[431] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '15_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_pl', node2.lemma)
    
    # L -> 2_du -> -62 -> L
    feats[432] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', '-62') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-62', node2.lemma)
    
    # L -> gen. pl. -> -243 -> L
    feats[433] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. pl.', '-243') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-243', node2.lemma)
    
    # L -> gen. pl. -> -279 -> L
    feats[434] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. pl.', '-279') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-279', node2.lemma)
    
    # C->102->L
    feats[435] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '102') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '102', node2.lemma)
            
    # C->-86->L
    feats[436] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-86') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-86', node2.lemma)
            
    # L -> 110 -> -14 -> L
    feats[437] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '110') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '110', '-14') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-14', node2.lemma)
    
    # L -> gen. pl. -> -302 -> L
    feats[438] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. pl.', '-302') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-302', node2.lemma)
    
    # C->-242->L
    feats[439] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-242') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-242', node2.lemma)
            
    # C->136->L
    feats[440] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '136') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '136', node2.lemma)
            
    # C->-303->L
    feats[441] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-303') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-303', node2.lemma)
            
    # C->-25->L
    feats[442] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-25') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-25', node2.lemma)
            
    # C->7_fp->L
    feats[443] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '7_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_fp', node2.lemma)
            
    # C->-293->L
    feats[444] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-293') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-293', node2.lemma)
            
    # L -> gen. pl. -> 12_fp -> L
    feats[445] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. pl.', '12_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_fp', node2.lemma)
    
    # C->-121->L
    feats[446] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-121') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-121', node2.lemma)
            
    # L -> gen. pl. -> gen. pl. -> L
    feats[447] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. pl.', 'gen. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen. pl.', node2.lemma)
    
    # C->54->L
    feats[448] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '54') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
            
    # C->160->L
    feats[449] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '160') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
            
    # L -> gen. pl. -> gen -> L
    feats[450] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. pl.', 'gen') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen', node2.lemma)
    
    # C->73->L
    feats[451] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '73') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '73', node2.lemma)
            
    # C->fp->L
    feats[452] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'fp', node2.lemma)
            
    # C->-16->L
    feats[453] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-16') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-16', node2.lemma)
            
    # C->-123->L
    feats[454] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-123') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-123', node2.lemma)
            
    # C->48->L
    feats[455] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '48') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '48', node2.lemma)
            
    # C->148->L
    feats[456] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '148') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '148', node2.lemma)
            
    # L -> voc. neutr. -> 39 -> L
    feats[457] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. neutr.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. neutr.', '39') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '39', node2.lemma)
    
    # L -> voc. neutr. -> -114 -> L
    feats[458] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. neutr.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. neutr.', '-114') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-114', node2.lemma)
    
    # C->112->L
    feats[459] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '112') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '112', node2.lemma)
            
    # C->-210->L
    feats[460] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-210') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-210', node2.lemma)
            
    # C->2_tp->L
    feats[461] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '2_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_tp', node2.lemma)
            
    # C->6_tp->L
    feats[462] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '6_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_tp', node2.lemma)
            
    # L -> gen. pl. -> loc. sg. -> L
    feats[463] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. pl.', 'loc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc. sg.', node2.lemma)
    
    # C->acc->L
    feats[464] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'acc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc', node2.lemma)
            
    # C->177->L
    feats[465] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '177') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '177', node2.lemma)
            
    # L -> gen. pl. -> -147 -> L
    feats[466] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. pl.', '-147') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-147', node2.lemma)
    
    # C->99->C
    feats[467] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '99') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '99', node2.cng)
            
    # L -> gen. pl. -> 179 -> L
    feats[468] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. pl.', '179') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '179', node2.lemma)
    
    # L -> gen. pl. -> 15_du -> L
    feats[469] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. pl.', '15_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_du', node2.lemma)
    
    # C->sg_tp->C
    feats[470] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'sg_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg_tp', node2.cng)
            
    # C->-200->C
    feats[471] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', node2.cng)
            
    # C->-97->C
    feats[472] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-97') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-97', node2.cng)
            
    # C->-59->C
    feats[473] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-59') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-59', node2.cng)
            
    # C->-93->C
    feats[474] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-93') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-93', node2.cng)
            
    # L -> -268 -> 48 -> L
    feats[475] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', '48') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '48', node2.lemma)
    
    # C->-10->C
    feats[476] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-10') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-10', node2.cng)
            
    # C->-26->C
    feats[477] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-26') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-26', node2.cng)
            
    # C->fem->C
    feats[478] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'fem', node2.cng)
            
    # C->-139->C
    feats[479] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', node2.cng)
            
    # C->5_fp->C
    feats[480] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '5_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', node2.cng)
            
    # L -> 110 -> -142 -> L
    feats[481] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '110') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '110', '-142') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-142', node2.lemma)
    
    # C->155->C
    feats[482] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '155') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '155', node2.cng)
            
    # C->156->C
    feats[483] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '156', node2.cng)
            
    # C->9_sp->C
    feats[484] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', node2.cng)
            
    # L -> gen. pl. -> 10_fp -> L
    feats[485] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. pl.', '10_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_fp', node2.lemma)
    
    # L -> gen. pl. -> 80 -> L
    feats[486] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. pl.', '80') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '80', node2.lemma)
    
    # C->2_du->C
    feats[487] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', node2.cng)
            
    # C->169->C
    feats[488] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '169') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '169', node2.cng)
            
    # L -> 110 -> -49 -> L
    feats[489] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '110') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '110', '-49') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-49', node2.lemma)
    
    # C->sg->C
    feats[490] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg', node2.cng)
            
    # L -> gen. pl. -> 30_pl -> L
    feats[491] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. pl.', '30_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_pl', node2.lemma)
    
    # C->-81->C
    feats[492] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', node2.cng)
            
    # C->-13->C
    feats[493] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-13') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-13', node2.cng)
            
    # C->30->C
    feats[494] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30', node2.cng)
            
    # L -> gen. pl. -> 94 -> L
    feats[495] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. pl.', '94') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '94', node2.lemma)
    
    # C->7_pl->C
    feats[496] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', node2.cng)
            
    # C->149->C
    feats[497] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', node2.cng)
            
    # L -> gen. pl. -> nom -> L
    feats[498] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. pl.', 'nom') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom', node2.lemma)
    
    # C->-132->C
    feats[499] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-132') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-132', node2.cng)
            
    # C->-261->C
    feats[500] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-261') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-261', node2.cng)
            
    # C->30_fp->C
    feats[501] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '30_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_fp', node2.cng)
            
    # C->du_fp->C
    feats[502] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'du_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_fp', node2.cng)
            
    # C->5_du->C
    feats[503] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '5_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_du', node2.cng)
            
    # C->sp->C
    feats[504] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', node2.cng)
            
    # C->-52->C
    feats[505] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-52') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-52', node2.cng)
            
    # C->-82->C
    feats[506] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-82') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-82', node2.cng)
            
    # L -> gen. pl. -> -47 -> L
    feats[507] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. pl.', '-47') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-47', node2.lemma)
    
    # C->88->C
    feats[508] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '88') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '88', node2.cng)
            
    # C->29_tp->C
    feats[509] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '29_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '29_tp', node2.cng)
            
    # C->-169->C
    feats[510] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-169') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-169', node2.cng)
            
    # C->-90->C
    feats[511] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-90', node2.cng)
            
    # C->72->C
    feats[512] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '72') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '72', node2.cng)
            
    # C->27_pl->C
    feats[513] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '27_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_pl', node2.cng)
            
    # C->-18->C
    feats[514] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', node2.cng)
            
    # C->90->C
    feats[515] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '90', node2.cng)
            
    # L -> gen. pl. -> pl_fp -> L
    feats[516] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. pl.', 'pl_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_fp', node2.lemma)
    
    # L -> gen. pl. -> 4_du -> L
    feats[517] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. pl.', '4_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_du', node2.lemma)
    
    # C->11_tp->C
    feats[518] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '11_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_tp', node2.cng)
            
    # C->42->C
    feats[519] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '42') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '42', node2.cng)
            
    # C->voc. du.->C
    feats[520] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'voc. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. du.', node2.cng)
            
    # C->-301->C
    feats[521] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-301') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-301', node2.cng)
            
    # C->-101->C
    feats[522] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-101') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-101', node2.cng)
            
    # C->-58->C
    feats[523] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-58', node2.cng)
            
    # C->13_sp->C
    feats[524] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '13_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_sp', node2.cng)
            
    # C->-241->C
    feats[525] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-241') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-241', node2.cng)
            
    # C->170->C
    feats[526] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '170') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '170', node2.cng)
            
    # L -> voc. neutr. -> -266 -> L
    feats[527] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. neutr.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. neutr.', '-266') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-266', node2.lemma)
    
    # C->-72->C
    feats[528] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-72') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-72', node2.cng)
            
    # C->nom. sg.->C
    feats[529] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'nom. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. sg.', node2.cng)
            
    # L -> -245 -> -78 -> L
    feats[530] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-245') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-245', '-78') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-78', node2.lemma)
    
    # L -> -245 -> -90 -> L
    feats[531] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-245') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-245', '-90') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-90', node2.lemma)
    
    # L -> gen. pl. -> -117 -> L
    feats[532] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. pl.', '-117') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-117', node2.lemma)
    
    # L -> 71 -> 117 -> L
    feats[533] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '71') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '71', '117') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '117', node2.lemma)
    
    # L -> 71 -> -31 -> L
    feats[534] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '71') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '71', '-31') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-31', node2.lemma)
    
    # L -> 71 -> -245 -> L
    feats[535] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '71') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '71', '-245') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-245', node2.lemma)
    
    # L -> gen. pl. -> -126 -> L
    feats[536] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. pl.', '-126') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-126', node2.lemma)
    
    # L -> 71 -> 95 -> L
    feats[537] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '71') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '71', '95') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '95', node2.lemma)
    
    # L -> gen. pl. -> -119 -> L
    feats[538] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. pl.', '-119') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-119', node2.lemma)
    
    # L -> 71 -> gen -> L
    feats[539] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '71') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '71', 'gen') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen', node2.lemma)
    
    # L -> 71 -> -56 -> L
    feats[540] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '71') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '71', '-56') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-56', node2.lemma)
    
    # L -> voc. neutr. -> -303 -> L
    feats[541] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. neutr.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. neutr.', '-303') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-303', node2.lemma)
    
    # L -> voc. neutr. -> instr. du. -> L
    feats[542] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. neutr.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. neutr.', 'instr. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. du.', node2.lemma)
    
    # L -> 71 -> 110 -> L
    feats[543] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '71') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '71', '110') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '110', node2.lemma)
    
    # L -> gen. pl. -> 11_sp -> L
    feats[544] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. pl.', '11_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_sp', node2.lemma)
    
    # L -> gen. pl. -> 3 -> L
    feats[545] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. pl.', '3') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3', node2.lemma)
    
    # C->9_fp->C
    feats[546] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '9_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_fp', node2.cng)
            
    # L -> gen. pl. -> 55 -> L
    feats[547] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. pl.', '55') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '55', node2.lemma)
    
    # L -> 71 -> 128 -> L
    feats[548] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '71') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '71', '128') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '128', node2.lemma)
    
    # L -> gen. pl. -> instr. du. -> L
    feats[549] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. pl.', 'instr. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. du.', node2.lemma)
    
    # L -> -245 -> 81 -> L
    feats[550] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-245') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-245', '81') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '81', node2.lemma)
    
    # L -> 30_fp -> instr. neutr. -> L
    feats[551] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_fp', 'instr. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. neutr.', node2.lemma)
    
    # L -> gen. pl. -> 49 -> L
    feats[552] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. pl.', '49') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '49', node2.lemma)
    
    # L -> 71 -> 80 -> L
    feats[553] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '71') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '71', '80') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '80', node2.lemma)
    
    # L -> gen. pl. -> 73 -> L
    feats[554] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. pl.', '73') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '73', node2.lemma)
    
    # L -> gen. pl. -> dat -> L
    feats[555] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. pl.', 'dat') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat', node2.lemma)
    
    # L -> gen. pl. -> -67 -> L
    feats[556] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. pl.', '-67') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-67', node2.lemma)
    
    # L -> 71 -> -64 -> L
    feats[557] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '71') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '71', '-64') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-64', node2.lemma)
    
    # L -> 71 -> 4_sp -> L
    feats[558] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '71') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '71', '4_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_sp', node2.lemma)
    
    # L -> 71 -> 13_fp -> L
    feats[559] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '71') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '71', '13_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_fp', node2.lemma)
    
    # L -> 71 -> -190 -> L
    feats[560] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '71') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '71', '-190') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-190', node2.lemma)
    
    # L -> -48 -> 7_sg -> L
    feats[561] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', '7_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_sg', node2.lemma)
    
    # L -> -48 -> voc. sg. -> L
    feats[562] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', 'voc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. sg.', node2.lemma)
    
    # L -> 71 -> -129 -> L
    feats[563] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '71') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '71', '-129') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-129', node2.lemma)
    
    # L -> 71 -> 13_sg -> L
    feats[564] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '71') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '71', '13_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_sg', node2.lemma)
    
    # L -> 71 -> 4_du -> L
    feats[565] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '71') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '71', '4_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_du', node2.lemma)
    
    # L -> -48 -> 30_du -> L
    feats[566] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', '30_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_du', node2.lemma)
    
    # L -> -48 -> 2_du -> L
    feats[567] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', '2_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_du', node2.lemma)
    
    # C->111->C
    feats[568] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '111') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '111', node2.cng)
            
    # C->-260->C
    feats[569] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-260') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-260', node2.cng)
            
    # L -> 71 -> dat. du. -> L
    feats[570] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '71') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '71', 'dat. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    
    # L -> 71 -> 36 -> L
    feats[571] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '71') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '71', '36') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '36', node2.lemma)
    
    # L -> -48 -> 30 -> L
    feats[572] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', '30') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30', node2.lemma)
    
    # L -> 27_fp -> 29_tp -> L
    feats[573] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_fp', '29_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '29_tp', node2.lemma)
    
    # L -> -245 -> 69 -> L
    feats[574] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-245') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-245', '69') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '69', node2.lemma)
    
    # L -> -111 -> instr. adj. -> L
    feats[575] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-111') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-111', 'instr. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. adj.', node2.lemma)
    
    # L -> -111 -> 132 -> L
    feats[576] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-111') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-111', '132') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '132', node2.lemma)
    
    # L -> -48 -> 9_pl -> L
    feats[577] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', '9_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_pl', node2.lemma)
    
    # L -> 27_fp -> 70 -> L
    feats[578] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_fp', '70') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '70', node2.lemma)
    
    # L -> 27_fp -> gen. pl. -> L
    feats[579] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_fp', 'gen. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen. pl.', node2.lemma)
    
    # C->-242->C
    feats[580] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-242') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-242', node2.cng)
            
    # L -> 71 -> 152 -> L
    feats[581] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '71') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '71', '152') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '152', node2.lemma)
    
    # L -> -111 -> 140 -> L
    feats[582] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-111') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-111', '140') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '140', node2.lemma)
    
    # L -> -48 -> 9_fp -> L
    feats[583] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', '9_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_fp', node2.lemma)
    
    # L -> -48 -> 39 -> L
    feats[584] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', '39') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '39', node2.lemma)
    
    # L -> -48 -> 118 -> L
    feats[585] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', '118') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '118', node2.lemma)
    
    # L -> -48 -> acc. sg. -> L
    feats[586] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', 'acc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. sg.', node2.lemma)
    
    # L -> -48 -> 9_du -> L
    feats[587] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', '9_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_du', node2.lemma)
    
    # C->12_sg->T
    feats[588] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '12_sg') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '12_sg', node2.tup)
            
    # C->-132->T
    feats[589] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-132') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-132', node2.tup)
            
    # L -> -48 -> pl_sp -> L
    feats[590] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', 'pl_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_sp', node2.lemma)
    
    # L -> 169 -> -98 -> L
    feats[591] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '169') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '169', '-98') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-98', node2.lemma)
    
    # L -> 169 -> 7_pl -> L
    feats[592] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '169') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '169', '7_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_pl', node2.lemma)
    
    # L -> 27_fp -> -163 -> L
    feats[593] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_fp', '-163') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-163', node2.lemma)
    
    # L -> -48 -> 9_tp -> L
    feats[594] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', '9_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_tp', node2.lemma)
    
    # L -> -48 -> 3_pl -> L
    feats[595] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', '3_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_pl', node2.lemma)
    
    # L -> -111 -> 2 -> L
    feats[596] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-111') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-111', '2') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2', node2.lemma)
    
    # L -> -111 -> -271 -> L
    feats[597] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-111') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-111', '-271') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-271', node2.lemma)
    
    # L -> -48 -> 3_sg -> L
    feats[598] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', '3_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_sg', node2.lemma)
    
    # L -> -48 -> 168 -> L
    feats[599] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', '168') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '168', node2.lemma)
    
    # L -> -48 -> 5_pl -> L
    feats[600] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', '5_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_pl', node2.lemma)
    
    # L -> -48 -> pl_fp -> L
    feats[601] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', 'pl_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_fp', node2.lemma)
    
    # L -> -48 -> -22 -> L
    feats[602] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', '-22') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-22', node2.lemma)
    
    # L -> -48 -> adj -> L
    feats[603] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', 'adj') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'adj', node2.lemma)
    
    # L -> -48 -> nom. pl. -> L
    feats[604] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', 'nom. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. pl.', node2.lemma)
    
    # L -> -48 -> 8_du -> L
    feats[605] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', '8_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_du', node2.lemma)
    
    # L -> -48 -> acc. fem -> L
    feats[606] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', 'acc. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. fem', node2.lemma)
    
    # L -> -48 -> 93 -> L
    feats[607] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', '93') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '93', node2.lemma)
    
    # L -> 27_fp -> 15_fp -> L
    feats[608] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_fp', '15_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_fp', node2.lemma)
    
    # L -> -48 -> -35 -> L
    feats[609] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', '-35') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    
    # C->-23->T
    feats[610] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-23') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-23', node2.tup)
            
    # L -> 169 -> 38 -> L
    feats[611] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '169') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '169', '38') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '38', node2.lemma)
    
    # L -> du_sp -> -62 -> L
    feats[612] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_sp', '-62') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-62', node2.lemma)
    
    # L -> -48 -> 139 -> L
    feats[613] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', '139') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '139', node2.lemma)
    
    # L -> -48 -> 6_sp -> L
    feats[614] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', '6_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_sp', node2.lemma)
    
    # L -> -48 -> 15_fp -> L
    feats[615] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', '15_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_fp', node2.lemma)
    
    # L -> -48 -> 15_sp -> L
    feats[616] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', '15_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_sp', node2.lemma)
    
    # L -> -111 -> -104 -> L
    feats[617] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-111') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-111', '-104') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-104', node2.lemma)
    
    # L -> -48 -> -46 -> L
    feats[618] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', '-46') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-46', node2.lemma)
    
    # L -> -302 -> -152 -> L
    feats[619] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', '-152') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-152', node2.lemma)
    
    # L -> -48 -> 92 -> L
    feats[620] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', '92') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '92', node2.lemma)
    
    # L -> -48 -> 162 -> L
    feats[621] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', '162') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '162', node2.lemma)
    
    # L -> -48 -> 172 -> L
    feats[622] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', '172') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '172', node2.lemma)
    
    # L -> du_sp -> 161 -> L
    feats[623] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_sp', '161') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '161', node2.lemma)
    
    # L -> du_sp -> 160 -> L
    feats[624] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_sp', '160') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    
    # C->abl->T
    feats[625] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'abl') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'abl', node2.tup)
            
    # C->-41->T
    feats[626] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-41') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-41', node2.tup)
            
    # L -> 169 -> acc. pl. -> L
    feats[627] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '169') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '169', 'acc. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. pl.', node2.lemma)
    
    # L -> -51 -> 75 -> L
    feats[628] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-51') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-51', '75') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '75', node2.lemma)
    
    # L -> dat. pl. -> sg -> L
    feats[629] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', 'sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg', node2.lemma)
    
    # L -> -230 -> -13 -> L
    feats[630] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '-13') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-13', node2.lemma)
    
    # L -> -51 -> -39 -> L
    feats[631] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-51') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-51', '-39') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-39', node2.lemma)
    
    # L -> -51 -> -122 -> L
    feats[632] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-51') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-51', '-122') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-122', node2.lemma)
    
    # C->-62->T
    feats[633] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-62') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-62', node2.tup)
            
    # L -> 169 -> 15_fp -> L
    feats[634] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '169') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '169', '15_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_fp', node2.lemma)
    
    # L -> 169 -> -14 -> L
    feats[635] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '169') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '169', '-14') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-14', node2.lemma)
    
    # L -> -230 -> 9_pl -> L
    feats[636] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '9_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_pl', node2.lemma)
    
    # C->151->T
    feats[637] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '151') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '151', node2.tup)
            
    # L -> 169 -> abl. pl. -> L
    feats[638] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '169') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '169', 'abl. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. pl.', node2.lemma)
    
    # L -> dat. pl. -> -92 -> L
    feats[639] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', '-92') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-92', node2.lemma)
    
    # L -> 169 -> 172 -> L
    feats[640] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '169') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '169', '172') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '172', node2.lemma)
    
    # L -> 169 -> -141 -> L
    feats[641] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '169') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '169', '-141') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-141', node2.lemma)
    
    # L -> 169 -> -104 -> L
    feats[642] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '169') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '169', '-104') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-104', node2.lemma)
    
    # L -> 169 -> -242 -> L
    feats[643] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '169') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '169', '-242') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-242', node2.lemma)
    
    # L -> 169 -> -25 -> L
    feats[644] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '169') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '169', '-25') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-25', node2.lemma)
    
    # L -> 169 -> -38 -> L
    feats[645] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '169') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '169', '-38') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-38', node2.lemma)
    
    # L -> 169 -> -293 -> L
    feats[646] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '169') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '169', '-293') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-293', node2.lemma)
    
    # L -> 169 -> -121 -> L
    feats[647] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '169') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '169', '-121') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-121', node2.lemma)
    
    # L -> 169 -> 54 -> L
    feats[648] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '169') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '169', '54') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    
    # L -> 169 -> instr. pl. -> L
    feats[649] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '169') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '169', 'instr. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. pl.', node2.lemma)
    
    # L -> 169 -> -142 -> L
    feats[650] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '169') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '169', '-142') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-142', node2.lemma)
    
    # L -> abl. sg. -> loc. sg. -> L
    feats[651] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'abl. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'abl. sg.', 'loc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc. sg.', node2.lemma)
    
    # L -> 169 -> -50 -> L
    feats[652] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '169') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '169', '-50') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-50', node2.lemma)
    
    # L -> 169 -> -123 -> L
    feats[653] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '169') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '169', '-123') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-123', node2.lemma)
    
    # T->99->L
    feats[654] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '99') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '99', node2.lemma)
            
    # T->27_du->L
    feats[655] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '27_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_du', node2.lemma)
            
    # L -> 169 -> acc -> L
    feats[656] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '169') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '169', 'acc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc', node2.lemma)
    
    # L -> 29_sg -> 27_du -> L
    feats[657] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '29_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '29_sg', '27_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_du', node2.lemma)
    
    # L -> -51 -> 94 -> L
    feats[658] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-51') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-51', '94') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '94', node2.lemma)
    
    # L -> 29_sg -> -97 -> L
    feats[659] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '29_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '29_sg', '-97') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-97', node2.lemma)
    
    # L -> 29_sg -> 7_sp -> L
    feats[660] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '29_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '29_sg', '7_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_sp', node2.lemma)
    
    # L -> dat. pl. -> 12_tp -> L
    feats[661] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', '12_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_tp', node2.lemma)
    
    # L -> dat. pl. -> 4_sp -> L
    feats[662] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', '4_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_sp', node2.lemma)
    
    # L -> 29_sg -> instr -> L
    feats[663] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '29_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '29_sg', 'instr') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr', node2.lemma)
    
    # L -> 29_sg -> 9_sp -> L
    feats[664] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '29_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '29_sg', '9_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_sp', node2.lemma)
    
    # L -> 29_sg -> 75 -> L
    feats[665] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '29_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '29_sg', '75') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '75', node2.lemma)
    
    # L -> 29_sg -> sg -> L
    feats[666] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '29_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '29_sg', 'sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg', node2.lemma)
    
    # L -> 29_sg -> 9_sg -> L
    feats[667] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '29_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '29_sg', '9_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_sg', node2.lemma)
    
    # L -> -230 -> 171 -> L
    feats[668] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '171') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '171', node2.lemma)
    
    # L -> 29_sg -> 7_pl -> L
    feats[669] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '29_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '29_sg', '7_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_pl', node2.lemma)
    
    # L -> 29_sg -> 14_fp -> L
    feats[670] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '29_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '29_sg', '14_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_fp', node2.lemma)
    
    # L -> 29_sg -> 12_sg -> L
    feats[671] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '29_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '29_sg', '12_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_sg', node2.lemma)
    
    # L -> 29_sg -> nom. du. -> L
    feats[672] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '29_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '29_sg', 'nom. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. du.', node2.lemma)
    
    # L -> abl. sg. -> 168 -> L
    feats[673] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'abl. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'abl. sg.', '168') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '168', node2.lemma)
    
    # L -> 29_sg -> 131 -> L
    feats[674] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '29_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '29_sg', '131') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '131', node2.lemma)
    
    # L -> 29_sg -> 16_pl -> L
    feats[675] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '29_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '29_sg', '16_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_pl', node2.lemma)
    
    # L -> 29_sg -> instr. adj. -> L
    feats[676] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '29_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '29_sg', 'instr. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. adj.', node2.lemma)
    
    # L -> 29_sg -> 132 -> L
    feats[677] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '29_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '29_sg', '132') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '132', node2.lemma)
    
    # L -> 29_sg -> -78 -> L
    feats[678] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '29_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '29_sg', '-78') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-78', node2.lemma)
    
    # L -> 29_sg -> -18 -> L
    feats[679] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '29_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '29_sg', '-18') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-18', node2.lemma)
    
    # L -> 29_sg -> -87 -> L
    feats[680] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '29_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '29_sg', '-87') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-87', node2.lemma)
    
    # L -> -157 -> -69 -> L
    feats[681] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '-69') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-69', node2.lemma)
    
    # L -> -51 -> instr. sg. -> L
    feats[682] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-51') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-51', 'instr. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. sg.', node2.lemma)
    
    # L -> 29_sg -> voc. du. -> L
    feats[683] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '29_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '29_sg', 'voc. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. du.', node2.lemma)
    
    # L -> 29_sg -> du_tp -> L
    feats[684] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '29_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '29_sg', 'du_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_tp', node2.lemma)
    
    # L -> 29_sg -> -101 -> L
    feats[685] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '29_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '29_sg', '-101') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-101', node2.lemma)
    
    # L -> 29_sg -> -53 -> L
    feats[686] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '29_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '29_sg', '-53') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-53', node2.lemma)
    
    # L -> 29_sg -> -241 -> L
    feats[687] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '29_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '29_sg', '-241') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-241', node2.lemma)
    
    # L -> 29_sg -> 170 -> L
    feats[688] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '29_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '29_sg', '170') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '170', node2.lemma)
    
    # L -> -230 -> voc. masc. -> L
    feats[689] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', 'voc. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. masc.', node2.lemma)
    
    # L -> dat. pl. -> 8_tp -> L
    feats[690] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', '8_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_tp', node2.lemma)
    
    # L -> -157 -> -23 -> L
    feats[691] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '-23') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-23', node2.lemma)
    
    # L -> -157 -> 12_fp -> L
    feats[692] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '12_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_fp', node2.lemma)
    
    # T->1->L
    feats[693] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '1') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '1', node2.lemma)
            
    # L -> -51 -> -50 -> L
    feats[694] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-51') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-51', '-50') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-50', node2.lemma)
    
    # L -> dat. pl. -> 2_tp -> L
    feats[695] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', '2_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_tp', node2.lemma)
    
    # L -> dat. pl. -> 6_tp -> L
    feats[696] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', '6_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_tp', node2.lemma)
    
    # L -> 29_sg -> 3_sg -> L
    feats[697] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '29_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '29_sg', '3_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_sg', node2.lemma)
    
    # L -> -157 -> abl -> L
    feats[698] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', 'abl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl', node2.lemma)
    
    # L -> -263 -> 173 -> L
    feats[699] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-263') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-263', '173') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '173', node2.lemma)
    
    # L -> 78 -> 12_pl -> L
    feats[700] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '78') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '78', '12_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_pl', node2.lemma)
    
    # L -> 78 -> 4_tp -> L
    feats[701] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '78') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '78', '4_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_tp', node2.lemma)
    
    # T->32->L
    feats[702] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '32') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '32', node2.lemma)
            
    # L -> -263 -> -78 -> L
    feats[703] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-263') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-263', '-78') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-78', node2.lemma)
    
    # L -> -157 -> 5_sp -> L
    feats[704] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '5_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sp', node2.lemma)
    
    # L -> -157 -> masc -> L
    feats[705] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', 'masc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'masc', node2.lemma)
    
    # L -> -308 -> -292 -> L
    feats[706] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '-292') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-292', node2.lemma)
    
    # L -> -263 -> 3_tp -> L
    feats[707] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-263') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-263', '3_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_tp', node2.lemma)
    
    # L -> 29_sg -> abl. du. -> L
    feats[708] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '29_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '29_sg', 'abl. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. du.', node2.lemma)
    
    # T->3->L
    feats[709] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '3') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3', node2.lemma)
            
    # T->-104->L
    feats[710] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-104') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-104', node2.lemma)
            
    # L -> -99 -> -157 -> L
    feats[711] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-99') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-99', '-157') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-157', node2.lemma)
    
    # L -> -99 -> 37 -> L
    feats[712] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-99') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-99', '37') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '37', node2.lemma)
    
    # T->-293->L
    feats[713] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-293') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-293', node2.lemma)
            
    # L -> 78 -> -147 -> L
    feats[714] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '78') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '78', '-147') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-147', node2.lemma)
    
    # T->gen. sg.->L
    feats[715] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'gen. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen. sg.', node2.lemma)
            
    # L -> -157 -> -156 -> L
    feats[716] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '-156') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-156', node2.lemma)
    
    # T->176->L
    feats[717] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '176') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '176', node2.lemma)
            
    # T->161->L
    feats[718] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '161') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '161', node2.lemma)
            
    # T->-142->L
    feats[719] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-142') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-142', node2.lemma)
            
    # T->-67->L
    feats[720] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-67') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-67', node2.lemma)
            
    # L -> -99 -> loc. sg. -> L
    feats[721] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-99') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-99', 'loc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc. sg.', node2.lemma)
    
    # L -> -99 -> 81 -> L
    feats[722] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-99') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-99', '81') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '81', node2.lemma)
    
    # T->-210->L
    feats[723] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-210') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-210', node2.lemma)
            
    # T->voc. fem->L
    feats[724] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'voc. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. fem', node2.lemma)
            
    # L -> -157 -> 5_sg -> L
    feats[725] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '5_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sg', node2.lemma)
    
    # T->-200->C
    feats[726] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', node2.cng)
            
    # L -> -308 -> 28_tp -> L
    feats[727] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '28_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28_tp', node2.lemma)
    
    # T->-30->C
    feats[728] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', node2.cng)
            
    # T->-26->C
    feats[729] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-26') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-26', node2.cng)
            
    # L -> sg -> voc. sg. -> L
    feats[730] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg', 'voc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. sg.', node2.lemma)
    
    # L -> sg -> voc. neutr. -> L
    feats[731] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg', 'voc. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. neutr.', node2.lemma)
    
    # L -> -99 -> 118 -> L
    feats[732] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-99') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-99', '118') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '118', node2.lemma)
    
    # T->8_sp->C
    feats[733] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '8_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_sp', node2.cng)
            
    # T->30_du->C
    feats[734] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '30_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_du', node2.cng)
            
    # T->sg->C
    feats[735] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg', node2.cng)
            
    # T->12_sg->C
    feats[736] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', node2.cng)
            
    # T->-132->C
    feats[737] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-132') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-132', node2.cng)
            
    # T->-261->C
    feats[738] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-261') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-261', node2.cng)
            
    # T->dat. pl.->C
    feats[739] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', node2.cng)
            
    # T->-276->C
    feats[740] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-276') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-276', node2.cng)
            
    # T->acc. neutr.->C
    feats[741] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'acc. neutr.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. neutr.', node2.cng)
            
    # L -> 78 -> pl_fp -> L
    feats[742] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '78') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '78', 'pl_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_fp', node2.lemma)
    
    # T->-82->C
    feats[743] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-82') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-82', node2.cng)
            
    # T->29_tp->C
    feats[744] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '29_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '29_tp', node2.cng)
            
    # T->130->C
    feats[745] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '130') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '130', node2.cng)
            
    # L -> -99 -> -163 -> L
    feats[746] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-99') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-99', '-163') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-163', node2.lemma)
    
    # L -> -99 -> 4_pl -> L
    feats[747] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-99') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-99', '4_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_pl', node2.lemma)
    
    # T->instr. adj.->C
    feats[748] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'instr. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. adj.', node2.cng)
            
    # T->72->C
    feats[749] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '72') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '72', node2.cng)
            
    # T->90->C
    feats[750] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '90', node2.cng)
            
    # T->-87->C
    feats[751] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-87') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-87', node2.cng)
            
    # L -> 16_tp -> acc. neutr. -> L
    feats[752] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', 'acc. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. neutr.', node2.lemma)
    
    # L -> 16_tp -> 5_du -> L
    feats[753] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', '5_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_du', node2.lemma)
    
    # T->voc. du.->C
    feats[754] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'voc. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. du.', node2.cng)
            
    # T->2_fp->C
    feats[755] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', node2.cng)
            
    # T->-84->C
    feats[756] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-84') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-84', node2.cng)
            
    # T->-101->C
    feats[757] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-101') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-101', node2.cng)
            
    # T->2_sp->C
    feats[758] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '2_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_sp', node2.cng)
            
    # T->-63->C
    feats[759] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-63') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-63', node2.cng)
            
    # L -> -263 -> 2_sg -> L
    feats[760] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-263') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-263', '2_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_sg', node2.lemma)
    
    # L -> -263 -> 93 -> L
    feats[761] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-263') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-263', '93') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '93', node2.lemma)
    
    # T->9_pl->C
    feats[762] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '9_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_pl', node2.cng)
            
    # T->-279->C
    feats[763] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-279') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-279', node2.cng)
            
    # T->101->C
    feats[764] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '101') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '101', node2.cng)
            
    # L -> -308 -> nom. pl. -> L
    feats[765] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', 'nom. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. pl.', node2.lemma)
    
    # T->du_sp->C
    feats[766] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'du_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_sp', node2.cng)
            
    # T->acc. adj.->C
    feats[767] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'acc. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. adj.', node2.cng)
            
    # L -> sg -> 101 -> L
    feats[768] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg', '101') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    
    # T->-31->C
    feats[769] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-31') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-31', node2.cng)
            
    # T->-137->C
    feats[770] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-137') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-137', node2.cng)
            
    # T->35->C
    feats[771] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', node2.cng)
            
    # T->3_tp->C
    feats[772] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '3_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_tp', node2.cng)
            
    # T->11_pl->C
    feats[773] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '11_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_pl', node2.cng)
            
    # L -> 110 -> -150 -> L
    feats[774] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '110') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '110', '-150') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-150', node2.lemma)
    
    # L -> 16_tp -> -23 -> L
    feats[775] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', '-23') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-23', node2.lemma)
    
    # L -> -99 -> 129 -> L
    feats[776] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-99') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-99', '129') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '129', node2.lemma)
    
    # L -> -308 -> 136 -> L
    feats[777] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '136') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '136', node2.lemma)
    
    # L -> -308 -> 178 -> L
    feats[778] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '178') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '178', node2.lemma)
    
    # L -> sg -> 50 -> L
    feats[779] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg', '50') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '50', node2.lemma)
    
    # L -> sg -> 182 -> L
    feats[780] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg', '182') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '182', node2.lemma)
    
    # L -> sg -> 171 -> L
    feats[781] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg', '171') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '171', node2.lemma)
    
    # L -> sg -> 3_du -> L
    feats[782] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg', '3_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    
    # L -> sg -> acc. pl. -> L
    feats[783] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg', 'acc. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. pl.', node2.lemma)
    
    # L -> sg -> 6_pl -> L
    feats[784] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg', '6_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_pl', node2.lemma)
    
    # L -> 16_tp -> -247 -> L
    feats[785] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', '-247') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-247', node2.lemma)
    
    # L -> 16_tp -> 9_du -> L
    feats[786] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', '9_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_du', node2.lemma)
    
    # L -> sg -> -47 -> L
    feats[787] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg', '-47') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-47', node2.lemma)
    
    # L -> 117 -> 155 -> L
    feats[788] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '117') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '117', '155') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '155', node2.lemma)
    
    # L -> 117 -> instr -> L
    feats[789] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '117') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '117', 'instr') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr', node2.lemma)
    
    # L -> -276 -> 7_pl -> L
    feats[790] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-276') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-276', '7_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_pl', node2.lemma)
    
    # L -> -276 -> 149 -> L
    feats[791] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-276') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-276', '149') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '149', node2.lemma)
    
    # L -> sg -> gen. du. -> L
    feats[792] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg', 'gen. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen. du.', node2.lemma)
    
    # L -> sg -> 4_du -> L
    feats[793] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg', '4_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_du', node2.lemma)
    
    # L -> sg -> 40 -> L
    feats[794] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg', '40') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '40', node2.lemma)
    
    # L -> sg -> 135 -> L
    feats[795] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg', '135') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '135', node2.lemma)
    
    # L -> sg -> neutr -> L
    feats[796] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg', 'neutr') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'neutr', node2.lemma)
    
    # L -> sg -> 36 -> L
    feats[797] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg', '36') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '36', node2.lemma)
    
    # L -> sg -> -249 -> L
    feats[798] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg', '-249') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-249', node2.lemma)
    
    # L -> sg -> 154 -> L
    feats[799] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg', '154') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '154', node2.lemma)
    
    # L -> sg -> -113 -> L
    feats[800] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg', '-113') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-113', node2.lemma)
    
    # L -> sg -> -34 -> L
    feats[801] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg', '-34') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-34', node2.lemma)
    
    # L -> voc -> 88 -> L
    feats[802] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc', '88') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '88', node2.lemma)
    
    # L -> sg -> 96 -> L
    feats[803] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg', '96') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '96', node2.lemma)
    
    # L -> sg -> -94 -> L
    feats[804] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg', '-94') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-94', node2.lemma)
    
    # L -> sg -> -35 -> L
    feats[805] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg', '-35') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    
    # L -> sg -> 11_du -> L
    feats[806] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg', '11_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_du', node2.lemma)
    
    # L -> -10 -> -279 -> L
    feats[807] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-10') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-10', '-279') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-279', node2.lemma)
    
    # L -> sg -> loc -> L
    feats[808] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg', 'loc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc', node2.lemma)
    
    # L -> sg -> -126 -> L
    feats[809] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg', '-126') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-126', node2.lemma)
    
    # L -> sg -> 139 -> L
    feats[810] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg', '139') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '139', node2.lemma)
    
    # L -> sg -> 91 -> L
    feats[811] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg', '91') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    
    # L -> sg -> 7_du -> L
    feats[812] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg', '7_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_du', node2.lemma)
    
    # L -> sg -> 122 -> L
    feats[813] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg', '122') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '122', node2.lemma)
    
    # L -> 117 -> 117 -> L
    feats[814] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '117') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '117', '117') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '117', node2.lemma)
    
    # L -> sg -> abl. pl. -> L
    feats[815] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg', 'abl. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. pl.', node2.lemma)
    
    # T->-141->C
    feats[816] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-141') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-141', node2.cng)
            
    # L -> -276 -> 38 -> L
    feats[817] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-276') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-276', '38') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '38', node2.lemma)
    
    # L -> sg -> -79 -> L
    feats[818] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg', '-79') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-79', node2.lemma)
    
    # L -> sg -> -149 -> L
    feats[819] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg', '-149') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-149', node2.lemma)
    
    # L -> sg -> 129 -> L
    feats[820] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg', '129') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '129', node2.lemma)
    
    # L -> sg -> -156 -> L
    feats[821] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg', '-156') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-156', node2.lemma)
    
    # L -> 117 -> 16_tp -> L
    feats[822] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '117') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '117', '16_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_tp', node2.lemma)
    
    # L -> sg -> -297 -> L
    feats[823] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg', '-297') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-297', node2.lemma)
    
    # L -> sg -> -104 -> L
    feats[824] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg', '-104') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-104', node2.lemma)
    
    # L -> sg -> 55 -> L
    feats[825] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg', '55') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '55', node2.lemma)
    
    # L -> sg -> 12_sp -> L
    feats[826] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg', '12_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_sp', node2.lemma)
    
    # L -> voc -> -306 -> L
    feats[827] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc', '-306') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-306', node2.lemma)
    
    # T->fp->C
    feats[828] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'fp', node2.cng)
            
    # T->-16->C
    feats[829] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-16') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-16', node2.cng)
            
    # L -> -10 -> -161 -> L
    feats[830] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-10') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-10', '-161') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-161', node2.lemma)
    
    # L -> -10 -> 180 -> L
    feats[831] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-10') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-10', '180') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '180', node2.lemma)
    
    # L -> 117 -> 10_sp -> L
    feats[832] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '117') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '117', '10_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_sp', node2.lemma)
    
    # L -> -10 -> 9_du -> L
    feats[833] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-10') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-10', '9_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_du', node2.lemma)
    
    # L -> -10 -> 3_du -> L
    feats[834] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-10') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-10', '3_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    
    # L -> 16_tp -> -115 -> L
    feats[835] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', '-115') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-115', node2.lemma)
    
    # L -> -10 -> 6_pl -> L
    feats[836] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-10') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-10', '6_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_pl', node2.lemma)
    
    # L -> -276 -> acc. pl. -> L
    feats[837] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-276') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-276', 'acc. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. pl.', node2.lemma)
    
    # L -> -10 -> 3_pl -> L
    feats[838] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-10') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-10', '3_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_pl', node2.lemma)
    
    # L -> -276 -> 6_pl -> L
    feats[839] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-276') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-276', '6_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_pl', node2.lemma)
    
    # L -> -10 -> 30_tp -> L
    feats[840] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-10') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-10', '30_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_tp', node2.lemma)
    
    # L -> -10 -> -163 -> L
    feats[841] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-10') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-10', '-163') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-163', node2.lemma)
    
    # L -> -10 -> 2 -> L
    feats[842] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-10') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-10', '2') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2', node2.lemma)
    
    # L -> -92 -> 169 -> L
    feats[843] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '169') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '169', node2.lemma)
    
    # L -> -10 -> 27_tp -> L
    feats[844] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-10') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-10', '27_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_tp', node2.lemma)
    
    # L -> -10 -> dat. du. -> L
    feats[845] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-10') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-10', 'dat. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    
    # L -> -10 -> 16_fp -> L
    feats[846] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-10') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-10', '16_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_fp', node2.lemma)
    
    # L -> -10 -> 142 -> L
    feats[847] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-10') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-10', '142') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '142', node2.lemma)
    
    # L -> -10 -> 111 -> L
    feats[848] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-10') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-10', '111') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '111', node2.lemma)
    
    # L -> 9_sg -> 30_sp -> L
    feats[849] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sg', '30_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_sp', node2.lemma)
    
    # L -> -10 -> 141 -> L
    feats[850] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-10') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-10', '141') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '141', node2.lemma)
    
    # L -> -10 -> -269 -> L
    feats[851] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-10') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-10', '-269') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-269', node2.lemma)
    
    # L -> -10 -> -83 -> L
    feats[852] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-10') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-10', '-83') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-83', node2.lemma)
    
    # L -> -10 -> nom. pl. -> L
    feats[853] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-10') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-10', 'nom. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. pl.', node2.lemma)
    
    # T->-84->T
    feats[854] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-84') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-84', node2.tup)
            
    # L -> -92 -> 88 -> L
    feats[855] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '88') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '88', node2.lemma)
    
    # L -> -10 -> 10_sg -> L
    feats[856] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-10') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-10', '10_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_sg', node2.lemma)
    
    # L -> -10 -> -42 -> L
    feats[857] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-10') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-10', '-42') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-42', node2.lemma)
    
    # L -> -10 -> 93 -> L
    feats[858] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-10') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-10', '93') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '93', node2.lemma)
    
    # L -> -10 -> pl_tp -> L
    feats[859] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-10') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-10', 'pl_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_tp', node2.lemma)
    
    # L -> -92 -> -241 -> L
    feats[860] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '-241') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-241', node2.lemma)
    
    # L -> -10 -> 8_tp -> L
    feats[861] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-10') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-10', '8_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_tp', node2.lemma)
    
    # L -> -10 -> -307 -> L
    feats[862] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-10') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-10', '-307') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-307', node2.lemma)
    
    # L -> -92 -> 34 -> L
    feats[863] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '34') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '34', node2.lemma)
    
    # L -> -10 -> 15_pl -> L
    feats[864] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-10') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-10', '15_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_pl', node2.lemma)
    
    # L -> -10 -> -62 -> L
    feats[865] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-10') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-10', '-62') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-62', node2.lemma)
    
    # L -> -92 -> -292 -> L
    feats[866] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '-292') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-292', node2.lemma)
    
    # L -> -10 -> instr. masc. -> L
    feats[867] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-10') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-10', 'instr. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. masc.', node2.lemma)
    
    # L -> -10 -> 15_sp -> L
    feats[868] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-10') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-10', '15_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_sp', node2.lemma)
    
    # L -> -92 -> acc. adj. -> L
    feats[869] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', 'acc. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. adj.', node2.lemma)
    
    # L -> 9_sg -> 14_sg -> L
    feats[870] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sg', '14_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_sg', node2.lemma)
    
    # L -> 9_sg -> 38 -> L
    feats[871] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sg', '38') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '38', node2.lemma)
    
    # L -> -10 -> -156 -> L
    feats[872] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-10') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-10', '-156') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-156', node2.lemma)
    
    # L -> -10 -> 11_sp -> L
    feats[873] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-10') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-10', '11_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_sp', node2.lemma)
    
    # L -> -92 -> 35 -> L
    feats[874] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '35') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '35', node2.lemma)
    
    # L -> -10 -> -28 -> L
    feats[875] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-10') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-10', '-28') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-28', node2.lemma)
    
    # L -> -10 -> 102 -> L
    feats[876] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-10') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-10', '102') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '102', node2.lemma)
    
    # T->-240->T
    feats[877] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-240') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-240', node2.tup)
            
    # L -> -10 -> -242 -> L
    feats[878] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-10') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-10', '-242') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-242', node2.lemma)
    
    # L -> -10 -> 136 -> L
    feats[879] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-10') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-10', '136') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '136', node2.lemma)
    
    # L -> -92 -> 12_fp -> L
    feats[880] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '12_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_fp', node2.lemma)
    
    # L -> -92 -> -306 -> L
    feats[881] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '-306') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-306', node2.lemma)
    
    # L -> -276 -> 27_sg -> L
    feats[882] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-276') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-276', '27_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_sg', node2.lemma)
    
    # L -> -276 -> -121 -> L
    feats[883] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-276') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-276', '-121') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-121', node2.lemma)
    
    # L -> -276 -> 49 -> L
    feats[884] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-276') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-276', '49') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '49', node2.lemma)
    
    # L -> -10 -> -296 -> L
    feats[885] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-10') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-10', '-296') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-296', node2.lemma)
    
    # L -> -10 -> -16 -> L
    feats[886] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-10') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-10', '-16') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-16', node2.lemma)
    
    # L -> -10 -> -123 -> L
    feats[887] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-10') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-10', '-123') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-123', node2.lemma)
    
    # L -> -10 -> -210 -> L
    feats[888] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-10') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-10', '-210') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-210', node2.lemma)
    
    # L -> -26 -> 99 -> L
    feats[889] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-26') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-26', '99') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '99', node2.lemma)
    
    # L -> -26 -> voc. sg. -> L
    feats[890] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-26') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-26', 'voc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. sg.', node2.lemma)
    
    # L -> -26 -> -10 -> L
    feats[891] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-26') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-26', '-10') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-10', node2.lemma)
    
    # L -> -26 -> 5_fp -> L
    feats[892] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-26') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-26', '5_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_fp', node2.lemma)
    
    # L -> -26 -> instr -> L
    feats[893] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-26') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-26', 'instr') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr', node2.lemma)
    
    # L -> -26 -> 9_sp -> L
    feats[894] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-26') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-26', '9_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_sp', node2.lemma)
    
    # L -> -26 -> nom. adj. -> L
    feats[895] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-26') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-26', 'nom. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. adj.', node2.lemma)
    
    # L -> -26 -> 30_du -> L
    feats[896] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-26') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-26', '30_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_du', node2.lemma)
    
    # L -> -26 -> 9_sg -> L
    feats[897] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-26') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-26', '9_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_sg', node2.lemma)
    
    # L -> -26 -> -98 -> L
    feats[898] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-26') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-26', '-98') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-98', node2.lemma)
    
    # L -> 9_sg -> -43 -> L
    feats[899] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sg', '-43') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-43', node2.lemma)
    
    # L -> -26 -> 30_fp -> L
    feats[900] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-26') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-26', '30_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_fp', node2.lemma)
    
    # L -> -26 -> -276 -> L
    feats[901] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-26') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-26', '-276') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-276', node2.lemma)
    
    # L -> -26 -> -52 -> L
    feats[902] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-26') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-26', '-52') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-52', node2.lemma)
    
    # L -> -26 -> 131 -> L
    feats[903] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-26') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-26', '131') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '131', node2.lemma)
    
    # L -> -220 -> 71 -> L
    feats[904] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-220') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-220', '71') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '71', node2.lemma)
    
    # L -> -26 -> 16_pl -> L
    feats[905] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-26') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-26', '16_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_pl', node2.lemma)
    
    # T->-54->T
    feats[906] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-54') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-54', node2.tup)
            
    # L -> -26 -> 130 -> L
    feats[907] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-26') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-26', '130') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '130', node2.lemma)
    
    # L -> -26 -> 4_tp -> L
    feats[908] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-26') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-26', '4_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_tp', node2.lemma)
    
    # L -> -26 -> instr. adj. -> L
    feats[909] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-26') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-26', 'instr. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. adj.', node2.lemma)
    
    # L -> du_fp -> 72 -> L
    feats[910] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_fp', '72') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '72', node2.lemma)
    
    # L -> -26 -> -44 -> L
    feats[911] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-26') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-26', '-44') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-44', node2.lemma)
    
    # L -> -26 -> du_tp -> L
    feats[912] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-26') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-26', 'du_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_tp', node2.lemma)
    
    # L -> -26 -> 2_fp -> L
    feats[913] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-26') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-26', '2_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_fp', node2.lemma)
    
    # L -> -26 -> -84 -> L
    feats[914] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-26') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-26', '-84') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-84', node2.lemma)
    
    # L -> -26 -> 60 -> L
    feats[915] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-26') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-26', '60') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '60', node2.lemma)
    
    # L -> -26 -> -122 -> L
    feats[916] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-26') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-26', '-122') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-122', node2.lemma)
    
    # L -> -26 -> 15_tp -> L
    feats[917] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-26') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-26', '15_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_tp', node2.lemma)
    
    # L -> -26 -> -91 -> L
    feats[918] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-26') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-26', '-91') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-91', node2.lemma)
    
    # L -> -26 -> abl. sg. -> L
    feats[919] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-26') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-26', 'abl. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. sg.', node2.lemma)
    
    # L -> -26 -> 10_pl -> L
    feats[920] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-26') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-26', '10_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_pl', node2.lemma)
    
    # L -> -26 -> 10_du -> L
    feats[921] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-26') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-26', '10_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_du', node2.lemma)
    
    # L -> 9_sg -> 91 -> L
    feats[922] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sg', '91') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    
    # L -> 9_sg -> 7_du -> L
    feats[923] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sg', '7_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_du', node2.lemma)
    
    # L -> -26 -> -137 -> L
    feats[924] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-26') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-26', '-137') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-137', node2.lemma)
    
    # L -> -26 -> -151 -> L
    feats[925] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-26') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-26', '-151') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-151', node2.lemma)
    
    # T->-149->T
    feats[926] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-149') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-149', node2.tup)
            
    # L -> -220 -> acc. adj. -> L
    feats[927] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-220') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-220', 'acc. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. adj.', node2.lemma)
    
    # L -> du_fp -> -33 -> L
    feats[928] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_fp', '-33') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-33', node2.lemma)
    
    # L -> acc. adj. -> gen -> L
    feats[929] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. adj.', 'gen') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen', node2.lemma)
    
    # L -> -220 -> -19 -> L
    feats[930] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-220') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-220', '-19') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-19', node2.lemma)
    
    # L -> -92 -> 112 -> L
    feats[931] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '112') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '112', node2.lemma)
    
    # L -> loc. du. -> fem -> L
    feats[932] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. du.', 'fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'fem', node2.lemma)
    
    # L -> -26 -> abl -> L
    feats[933] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-26') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-26', 'abl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl', node2.lemma)
    
    # L -> acc. adj. -> 12_tp -> L
    feats[934] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. adj.', '12_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_tp', node2.lemma)
    
    # L -> du_fp -> 135 -> L
    feats[935] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_fp', '135') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '135', node2.lemma)
    
    # L -> nom. neutr. -> 12_pl -> L
    feats[936] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. neutr.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. neutr.', '12_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_pl', node2.lemma)
    
    # L -> loc. du. -> nom. sg. -> L
    feats[937] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. du.', 'nom. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. sg.', node2.lemma)
    
    # L -> loc. du. -> -45 -> L
    feats[938] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. du.', '-45') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-45', node2.lemma)
    
    # L -> -26 -> 11_du -> L
    feats[939] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-26') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-26', '11_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_du', node2.lemma)
    
    # L -> -26 -> -27 -> L
    feats[940] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-26') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-26', '-27') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-27', node2.lemma)
    
    # L -> loc. du. -> -99 -> L
    feats[941] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. du.', '-99') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-99', node2.lemma)
    
    # L -> -26 -> abl. du. -> L
    feats[942] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-26') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-26', 'abl. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. du.', node2.lemma)
    
    # L -> acc. adj. -> 14_sp -> L
    feats[943] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. adj.', '14_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_sp', node2.lemma)
    
    # L -> loc. du. -> -157 -> L
    feats[944] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. du.', '-157') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-157', node2.lemma)
    
    # L -> du_fp -> 138 -> L
    feats[945] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_fp', '138') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '138', node2.lemma)
    
    # L -> du_fp -> 15_sp -> L
    feats[946] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_fp', '15_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_sp', node2.lemma)
    
    # L -> -26 -> instr. du. -> L
    feats[947] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-26') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-26', 'instr. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. du.', node2.lemma)
    
    # L -> loc. du. -> -96 -> L
    feats[948] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. du.', '-96') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-96', node2.lemma)
    
    # L -> loc. du. -> 41 -> L
    feats[949] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. du.', '41') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '41', node2.lemma)
    
    # L -> loc. du. -> pl -> L
    feats[950] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. du.', 'pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    
    # L -> loc. du. -> 180 -> L
    feats[951] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. du.', '180') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '180', node2.lemma)
    
    # L -> loc. du. -> 5_tp -> L
    feats[952] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. du.', '5_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_tp', node2.lemma)
    
    # L -> loc. du. -> 97 -> L
    feats[953] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. du.', '97') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '97', node2.lemma)
    
    # L -> fem -> 7_sg -> L
    feats[954] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'fem', '7_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_sg', node2.lemma)
    
    # L -> 99 -> 9_du -> L
    feats[955] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '99') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '99', '9_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_du', node2.lemma)
    
    # L -> acc. adj. -> -210 -> L
    feats[956] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. adj.', '-210') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-210', node2.lemma)
    
    # L -> loc. du. -> pl_sp -> L
    feats[957] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. du.', 'pl_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_sp', node2.lemma)
    
    # L -> fem -> 5_fp -> L
    feats[958] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'fem', '5_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_fp', node2.lemma)
    
    # L -> -220 -> acc -> L
    feats[959] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-220') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-220', 'acc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc', node2.lemma)
    
    # L -> loc. du. -> -41 -> L
    feats[960] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. du.', '-41') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-41', node2.lemma)
    
    # L -> loc. du. -> 150 -> L
    feats[961] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. du.', '150') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '150', node2.lemma)
    
    # L -> loc. du. -> 3_sg -> L
    feats[962] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. du.', '3_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_sg', node2.lemma)
    
    # L -> fem -> 30 -> L
    feats[963] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'fem', '30') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30', node2.lemma)
    
    # L -> loc. du. -> 16_fp -> L
    feats[964] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. du.', '16_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_fp', node2.lemma)
    
    # L -> loc. du. -> -249 -> L
    feats[965] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. du.', '-249') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-249', node2.lemma)
    
    # L -> fem -> 29_tp -> L
    feats[966] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'fem', '29_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '29_tp', node2.lemma)
    
    # L -> fem -> -169 -> L
    feats[967] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'fem', '-169') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-169', node2.lemma)
    
    # L -> fem -> 130 -> L
    feats[968] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'fem', '130') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '130', node2.lemma)
    
    # L -> loc. du. -> -299 -> L
    feats[969] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. du.', '-299') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-299', node2.lemma)
    
    # L -> loc. du. -> 76 -> L
    feats[970] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. du.', '76') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '76', node2.lemma)
    
    # L -> loc. du. -> 32 -> L
    feats[971] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. du.', '32') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '32', node2.lemma)
    
    # L -> fem -> voc. du. -> L
    feats[972] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'fem', 'voc. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. du.', node2.lemma)
    
    # L -> fem -> du_tp -> L
    feats[973] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'fem', 'du_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_tp', node2.lemma)
    
    # L -> 99 -> -42 -> L
    feats[974] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '99') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '99', '-42') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-42', node2.lemma)
    
    # L -> 99 -> 8_du -> L
    feats[975] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '99') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '99', '8_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_du', node2.lemma)
    
    # L -> 8_sg -> -18 -> L
    feats[976] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_sg', '-18') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-18', node2.lemma)
    
    # L -> du -> 13_tp -> L
    feats[977] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du', '13_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_tp', node2.lemma)
    
    # L -> nom. neutr. -> -14 -> L
    feats[978] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. neutr.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. neutr.', '-14') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-14', node2.lemma)
    
    # L -> loc. du. -> 54 -> L
    feats[979] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. du.', '54') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    
    # L -> fem -> -112 -> L
    feats[980] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'fem', '-112') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    
    # L -> du -> 15_du -> L
    feats[981] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du', '15_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_du', node2.lemma)
    
    # L -> nom. neutr. -> -210 -> L
    feats[982] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. neutr.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. neutr.', '-210') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-210', node2.lemma)
    
    # L -> 27_du -> 75 -> L
    feats[983] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', '75') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '75', node2.lemma)
    
    # L -> 27_du -> 30_du -> L
    feats[984] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', '30_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_du', node2.lemma)
    
    # L -> 27_du -> -81 -> L
    feats[985] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', '-81') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-81', node2.lemma)
    
    # L -> 27_du -> 7_pl -> L
    feats[986] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', '7_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_pl', node2.lemma)
    
    # L -> 27_du -> -132 -> L
    feats[987] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', '-132') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-132', node2.lemma)
    
    # L -> 27_du -> acc. neutr. -> L
    feats[988] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', 'acc. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. neutr.', node2.lemma)
    
    # L -> 27_du -> -69 -> L
    feats[989] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', '-69') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-69', node2.lemma)
    
    # L -> -81 -> 173 -> L
    feats[990] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', '173') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '173', node2.lemma)
    
    # L -> fem -> 142 -> L
    feats[991] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'fem', '142') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '142', node2.lemma)
    
    # L -> -81 -> 159 -> L
    feats[992] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', '159') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '159', node2.lemma)
    
    # L -> fem -> -249 -> L
    feats[993] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'fem', '-249') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-249', node2.lemma)
    
    # L -> 27_du -> 130 -> L
    feats[994] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', '130') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '130', node2.lemma)
    
    # L -> 27_du -> 4_tp -> L
    feats[995] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', '4_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_tp', node2.lemma)
    
    # L -> 27_du -> instr. adj. -> L
    feats[996] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', 'instr. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. adj.', node2.lemma)
    
    # L -> 27_du -> 132 -> L
    feats[997] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', '132') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '132', node2.lemma)
    
    # L -> 27_du -> 27_pl -> L
    feats[998] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', '27_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_pl', node2.lemma)
    
    # L -> 27_du -> 90 -> L
    feats[999] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', '90') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '90', node2.lemma)
    
    # L -> 27_du -> 11_tp -> L
    feats[1000] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', '11_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_tp', node2.lemma)
    
    # L -> 27_du -> 74 -> L
    feats[1001] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', '74') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '74', node2.lemma)
    
    # L -> 27_du -> -301 -> L
    feats[1002] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', '-301') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-301', node2.lemma)
    
    # L -> acc. neutr. -> 32 -> L
    feats[1003] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. neutr.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. neutr.', '32') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '32', node2.lemma)
    
    # L -> 27_du -> -58 -> L
    feats[1004] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', '-58') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-58', node2.lemma)
    
    # L -> 27_du -> -241 -> L
    feats[1005] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', '-241') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-241', node2.lemma)
    
    # L -> 27_du -> 60 -> L
    feats[1006] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', '60') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '60', node2.lemma)
    
    # L -> 27_du -> -17 -> L
    feats[1007] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', '-17') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-17', node2.lemma)
    
    # L -> 27_du -> -292 -> L
    feats[1008] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', '-292') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-292', node2.lemma)
    
    # L -> 27_du -> -63 -> L
    feats[1009] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', '-63') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-63', node2.lemma)
    
    # L -> 27_du -> 7_tp -> L
    feats[1010] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', '7_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_tp', node2.lemma)
    
    # L -> 27_du -> -68 -> L
    feats[1011] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', '-68') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-68', node2.lemma)
    
    # L -> 27_du -> -39 -> L
    feats[1012] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', '-39') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-39', node2.lemma)
    
    # L -> 27_du -> 15_tp -> L
    feats[1013] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', '15_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_tp', node2.lemma)
    
    # L -> 27_du -> 117 -> L
    feats[1014] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', '117') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '117', node2.lemma)
    
    # L -> 27_du -> -31 -> L
    feats[1015] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', '-31') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-31', node2.lemma)
    
    # L -> 27_du -> -137 -> L
    feats[1016] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', '-137') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-137', node2.lemma)
    
    # L -> 27_du -> -151 -> L
    feats[1017] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', '-151') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-151', node2.lemma)
    
    # L -> 27_du -> -262 -> L
    feats[1018] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', '-262') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-262', node2.lemma)
    
    # L -> 27_du -> -245 -> L
    feats[1019] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', '-245') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-245', node2.lemma)
    
    # L -> 27_du -> 16_tp -> L
    feats[1020] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', '16_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_tp', node2.lemma)
    
    # L -> 27_du -> -92 -> L
    feats[1021] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', '-92') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-92', node2.lemma)
    
    # L -> 27_du -> 38 -> L
    feats[1022] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', '38') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '38', node2.lemma)
    
    # L -> 27_du -> 37 -> L
    feats[1023] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', '37') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '37', node2.lemma)
    
    # L -> 27_du -> 95 -> L
    feats[1024] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', '95') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '95', node2.lemma)
    
    # L -> 27_du -> 12_fp -> L
    feats[1025] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', '12_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_fp', node2.lemma)
    
    # L -> 27_du -> voc -> L
    feats[1026] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', 'voc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc', node2.lemma)
    
    # L -> 27_du -> -33 -> L
    feats[1027] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', '-33') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-33', node2.lemma)
    
    # L -> 27_du -> 33 -> L
    feats[1028] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', '33') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '33', node2.lemma)
    
    # L -> fem -> 27_sg -> L
    feats[1029] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'fem', '27_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_sg', node2.lemma)
    
    # L -> acc. neutr. -> 11_fp -> L
    feats[1030] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. neutr.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. neutr.', '11_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_fp', node2.lemma)
    
    # L -> 27_du -> 41 -> L
    feats[1031] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', '41') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '41', node2.lemma)
    
    # L -> acc. neutr. -> 152 -> L
    feats[1032] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. neutr.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. neutr.', '152') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '152', node2.lemma)
    
    # L -> 27_du -> -147 -> L
    feats[1033] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', '-147') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-147', node2.lemma)
    
    # L -> 27_du -> loc. sg. -> L
    feats[1034] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', 'loc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc. sg.', node2.lemma)
    
    # L -> 27_du -> 81 -> L
    feats[1035] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', '81') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '81', node2.lemma)
    
    # L -> 27_du -> -158 -> L
    feats[1036] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', '-158') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-158', node2.lemma)
    
    # L -> 27_du -> 15_sg -> L
    feats[1037] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', '15_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_sg', node2.lemma)
    
    # L -> 27_du -> pl -> L
    feats[1038] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', 'pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    
    # L -> 27_du -> 1 -> L
    feats[1039] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', '1') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '1', node2.lemma)
    
    # L -> 27_du -> 128 -> L
    feats[1040] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', '128') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '128', node2.lemma)
    
    # L -> 27_du -> -283 -> L
    feats[1041] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', '-283') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-283', node2.lemma)
    
    # L -> 27_du -> -114 -> L
    feats[1042] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', '-114') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-114', node2.lemma)
    
    # L -> -81 -> 10_sp -> L
    feats[1043] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', '10_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_sp', node2.lemma)
    
    # L -> -81 -> 97 -> L
    feats[1044] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', '97') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '97', node2.lemma)
    
    # L -> -81 -> 10_fp -> L
    feats[1045] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', '10_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_fp', node2.lemma)
    
    # L -> -81 -> 89 -> L
    feats[1046] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', '89') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '89', node2.lemma)
    
    # L -> -81 -> 3_du -> L
    feats[1047] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', '3_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    
    # L -> du -> voc. fem -> L
    feats[1048] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du', 'voc. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. fem', node2.lemma)
    
    # L -> -81 -> -41 -> L
    feats[1049] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', '-41') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-41', node2.lemma)
    
    # L -> -81 -> 175 -> L
    feats[1050] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', '175') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '175', node2.lemma)
    
    # L -> -81 -> -154 -> L
    feats[1051] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', '-154') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-154', node2.lemma)
    
    # L -> -81 -> 134 -> L
    feats[1052] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', '134') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '134', node2.lemma)
    
    # L -> -81 -> 150 -> L
    feats[1053] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', '150') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '150', node2.lemma)
    
    # L -> -81 -> 30_tp -> L
    feats[1054] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', '30_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_tp', node2.lemma)
    
    # L -> -139 -> 7_pl -> L
    feats[1055] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '7_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_pl', node2.lemma)
    
    # L -> -81 -> 168 -> L
    feats[1056] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', '168') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '168', node2.lemma)
    
    # L -> -139 -> -98 -> L
    feats[1057] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-98') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-98', node2.lemma)
    
    # L -> -81 -> 14_tp -> L
    feats[1058] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', '14_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_tp', node2.lemma)
    
    # L -> -81 -> -102 -> L
    feats[1059] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', '-102') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-102', node2.lemma)
    
    # L -> -81 -> gen. du. -> L
    feats[1060] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', 'gen. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen. du.', node2.lemma)
    
    # L -> -81 -> 4_du -> L
    feats[1061] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', '4_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_du', node2.lemma)
    
    # L -> -139 -> 116 -> L
    feats[1062] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '116') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '116', node2.lemma)
    
    # L -> -139 -> -82 -> L
    feats[1063] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-82') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-82', node2.lemma)
    
    # L -> -81 -> dat. du. -> L
    feats[1064] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', 'dat. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    
    # L -> -81 -> 36 -> L
    feats[1065] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', '36') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '36', node2.lemma)
    
    # L -> 5_du -> 173 -> L
    feats[1066] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_du', '173') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '173', node2.lemma)
    
    # L -> -81 -> 142 -> L
    feats[1067] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', '142') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '142', node2.lemma)
    
    # L -> -139 -> -169 -> L
    feats[1068] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-169') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-169', node2.lemma)
    
    # L -> -139 -> 130 -> L
    feats[1069] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '130') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '130', node2.lemma)
    
    # L -> -81 -> -22 -> L
    feats[1070] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', '-22') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-22', node2.lemma)
    
    # L -> -81 -> adj -> L
    feats[1071] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', 'adj') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'adj', node2.lemma)
    
    # L -> -139 -> 72 -> L
    feats[1072] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '72') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '72', node2.lemma)
    
    # L -> 27_du -> -55 -> L
    feats[1073] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', '-55') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-55', node2.lemma)
    
    # L -> -139 -> du_tp -> L
    feats[1074] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', 'du_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_tp', node2.lemma)
    
    # L -> -81 -> masc -> L
    feats[1075] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', 'masc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'masc', node2.lemma)
    
    # L -> -81 -> -83 -> L
    feats[1076] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', '-83') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-83', node2.lemma)
    
    # L -> gen -> -82 -> L
    feats[1077] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen', '-82') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-82', node2.lemma)
    
    # L -> -81 -> nom. fem -> L
    feats[1078] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', 'nom. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. fem', node2.lemma)
    
    # L -> -139 -> -241 -> L
    feats[1079] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-241') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-241', node2.lemma)
    
    # L -> 10_pl -> -78 -> L
    feats[1080] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_pl', '-78') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-78', node2.lemma)
    
    # L -> -81 -> -266 -> L
    feats[1081] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', '-266') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-266', node2.lemma)
    
    # L -> 10_pl -> -90 -> L
    feats[1082] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_pl', '-90') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-90', node2.lemma)
    
    # L -> -139 -> 2_sp -> L
    feats[1083] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '2_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_sp', node2.lemma)
    
    # L -> -81 -> -42 -> L
    feats[1084] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', '-42') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-42', node2.lemma)
    
    # L -> -81 -> 3_fp -> L
    feats[1085] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', '3_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_fp', node2.lemma)
    
    # L -> -139 -> 34 -> L
    feats[1086] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '34') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '34', node2.lemma)
    
    # L -> -81 -> 56 -> L
    feats[1087] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', '56') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '56', node2.lemma)
    
    # L -> -139 -> -17 -> L
    feats[1088] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-17') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-17', node2.lemma)
    
    # L -> -81 -> 2_sg -> L
    feats[1089] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', '2_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_sg', node2.lemma)
    
    # L -> -81 -> 93 -> L
    feats[1090] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', '93') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '93', node2.lemma)
    
    # L -> -139 -> -63 -> L
    feats[1091] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-63') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-63', node2.lemma)
    
    # L -> -81 -> 28_sg -> L
    feats[1092] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', '28_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28_sg', node2.lemma)
    
    # L -> -139 -> -292 -> L
    feats[1093] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-292') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-292', node2.lemma)
    
    # L -> -139 -> -68 -> L
    feats[1094] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-68') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-68', node2.lemma)
    
    # L -> -139 -> -39 -> L
    feats[1095] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-39') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-39', node2.lemma)
    
    # L -> -139 -> 15_tp -> L
    feats[1096] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '15_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_tp', node2.lemma)
    
    # L -> -139 -> -99 -> L
    feats[1097] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-99') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-99', node2.lemma)
    
    # L -> -81 -> 8_tp -> L
    feats[1098] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', '8_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_tp', node2.lemma)
    
    # L -> -139 -> 10_du -> L
    feats[1099] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '10_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_du', node2.lemma)
    
    # L -> -81 -> 122 -> L
    feats[1100] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', '122') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '122', node2.lemma)
    
    # L -> -81 -> 15_pl -> L
    feats[1101] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', '15_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_pl', node2.lemma)
    
    # L -> -81 -> -62 -> L
    feats[1102] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', '-62') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-62', node2.lemma)
    
    # L -> -81 -> -37 -> L
    feats[1103] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', '-37') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-37', node2.lemma)
    
    # L -> -81 -> 6_sp -> L
    feats[1104] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', '6_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_sp', node2.lemma)
    
    # L -> -139 -> -268 -> L
    feats[1105] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-268') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-268', node2.lemma)
    
    # L -> -81 -> -119 -> L
    feats[1106] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', '-119') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-119', node2.lemma)
    
    # L -> -139 -> -245 -> L
    feats[1107] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-245') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-245', node2.lemma)
    
    # L -> -81 -> -77 -> L
    feats[1108] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', '-77') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-77', node2.lemma)
    
    # L -> -139 -> 16_tp -> L
    feats[1109] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '16_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_tp', node2.lemma)
    
    # L -> -139 -> -157 -> L
    feats[1110] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-157') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-157', node2.lemma)
    
    # L -> -81 -> 100 -> L
    feats[1111] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', '100') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '100', node2.lemma)
    
    # L -> -81 -> 115 -> L
    feats[1112] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', '115') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '115', node2.lemma)
    
    # L -> -139 -> 14_sg -> L
    feats[1113] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '14_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_sg', node2.lemma)
    
    # L -> -139 -> 137 -> L
    feats[1114] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '137') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    
    # L -> -81 -> 11_sg -> L
    feats[1115] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', '11_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_sg', node2.lemma)
    
    # L -> -81 -> 92 -> L
    feats[1116] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', '92') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '92', node2.lemma)
    
    # L -> 5_du -> 14_sg -> L
    feats[1117] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_du', '14_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_sg', node2.lemma)
    
    # L -> 5_du -> 38 -> L
    feats[1118] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_du', '38') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '38', node2.lemma)
    
    # L -> -139 -> -291 -> L
    feats[1119] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-291') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-291', node2.lemma)
    
    # L -> -81 -> 162 -> L
    feats[1120] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', '162') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '162', node2.lemma)
    
    # L -> -139 -> -48 -> L
    feats[1121] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-48') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-48', node2.lemma)
    
    # L -> -81 -> -141 -> L
    feats[1122] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', '-141') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-141', node2.lemma)
    
    # L -> -81 -> 129 -> L
    feats[1123] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', '129') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '129', node2.lemma)
    
    # L -> -139 -> gen -> L
    feats[1124] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', 'gen') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen', node2.lemma)
    
    # L -> -139 -> -33 -> L
    feats[1125] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-33') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-33', node2.lemma)
    
    # L -> -139 -> 28 -> L
    feats[1126] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '28') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
    
    # L -> -139 -> -56 -> L
    feats[1127] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-56') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-56', node2.lemma)
    
    # L -> -139 -> -144 -> L
    feats[1128] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-144') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-144', node2.lemma)
    
    # L -> gen -> nom. neutr. -> L
    feats[1129] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen', 'nom. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. neutr.', node2.lemma)
    
    # L -> -139 -> -24 -> L
    feats[1130] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-24') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    
    # L -> 27_du -> -303 -> L
    feats[1131] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', '-303') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-303', node2.lemma)
    
    # L -> 27_du -> 178 -> L
    feats[1132] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', '178') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '178', node2.lemma)
    
    # L -> -139 -> -96 -> L
    feats[1133] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-96') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-96', node2.lemma)
    
    # L -> -139 -> -240 -> L
    feats[1134] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-240') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-240', node2.lemma)
    
    # L -> -139 -> 41 -> L
    feats[1135] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '41') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '41', node2.lemma)
    
    # L -> -139 -> loc. sg. -> L
    feats[1136] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', 'loc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc. sg.', node2.lemma)
    
    # L -> -139 -> 81 -> L
    feats[1137] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '81') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '81', node2.lemma)
    
    # L -> -139 -> 6_fp -> L
    feats[1138] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '6_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_fp', node2.lemma)
    
    # L -> -139 -> 179 -> L
    feats[1139] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '179') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '179', node2.lemma)
    
    # L -> -139 -> 28_tp -> L
    feats[1140] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '28_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28_tp', node2.lemma)
    
    # L -> -139 -> 58 -> L
    feats[1141] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '58') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '58', node2.lemma)
    
    # L -> 10_tp -> -103 -> L
    feats[1142] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_tp', '-103') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-103', node2.lemma)
    
    # L -> gen -> 41 -> L
    feats[1143] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen', '41') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '41', node2.lemma)
    
    # L -> -139 -> -131 -> L
    feats[1144] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-131') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-131', node2.lemma)
    
    # L -> -139 -> 97 -> L
    feats[1145] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '97') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '97', node2.lemma)
    
    # L -> -139 -> 10_fp -> L
    feats[1146] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '10_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_fp', node2.lemma)
    
    # L -> -139 -> 80 -> L
    feats[1147] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '80') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '80', node2.lemma)
    
    # L -> -139 -> acc. sg. -> L
    feats[1148] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', 'acc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. sg.', node2.lemma)
    
    # L -> -139 -> 89 -> L
    feats[1149] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '89') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '89', node2.lemma)
    
    # L -> -139 -> 29 -> L
    feats[1150] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '29') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '29', node2.lemma)
    
    # L -> -139 -> -247 -> L
    feats[1151] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-247') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-247', node2.lemma)
    
    # L -> -139 -> 9_du -> L
    feats[1152] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '9_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_du', node2.lemma)
    
    # L -> -139 -> -64 -> L
    feats[1153] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-64') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-64', node2.lemma)
    
    # L -> -139 -> 12_tp -> L
    feats[1154] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '12_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_tp', node2.lemma)
    
    # L -> -139 -> 4_sp -> L
    feats[1155] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '4_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_sp', node2.lemma)
    
    # L -> -139 -> -190 -> L
    feats[1156] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-190') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-190', node2.lemma)
    
    # L -> -139 -> acc. pl. -> L
    feats[1157] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', 'acc. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. pl.', node2.lemma)
    
    # L -> -139 -> 6_pl -> L
    feats[1158] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '6_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_pl', node2.lemma)
    
    # L -> 5_du -> -154 -> L
    feats[1159] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_du', '-154') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-154', node2.lemma)
    
    # L -> -13 -> 30 -> L
    feats[1160] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-13') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-13', '30') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30', node2.lemma)
    
    # L -> 5_du -> -43 -> L
    feats[1161] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_du', '-43') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-43', node2.lemma)
    
    # L -> 5_du -> 14_tp -> L
    feats[1162] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_du', '14_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_tp', node2.lemma)
    
    # L -> gen -> 134 -> L
    feats[1163] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen', '134') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '134', node2.lemma)
    
    # L -> 7_sg -> 88 -> L
    feats[1164] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_sg', '88') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '88', node2.lemma)
    
    # L -> -152 -> -276 -> L
    feats[1165] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-152') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-152', '-276') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-276', node2.lemma)
    
    # L -> 10_pl -> 8_fp -> L
    feats[1166] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_pl', '8_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    
    # L -> -139 -> 14_sp -> L
    feats[1167] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '14_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_sp', node2.lemma)
    
    # L -> -13 -> du_sp -> L
    feats[1168] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-13') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-13', 'du_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_sp', node2.lemma)
    
    # L -> 7_sg -> -151 -> L
    feats[1169] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_sg', '-151') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-151', node2.lemma)
    
    # L -> gen -> 4_sg -> L
    feats[1170] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen', '4_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_sg', node2.lemma)
    
    # L -> gen -> pl_tp -> L
    feats[1171] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen', 'pl_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_tp', node2.lemma)
    
    # L -> 7_sg -> 37 -> L
    feats[1172] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_sg', '37') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '37', node2.lemma)
    
    # L -> 7_sg -> 95 -> L
    feats[1173] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_sg', '95') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '95', node2.lemma)
    
    # L -> 10_pl -> 162 -> L
    feats[1174] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_pl', '162') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '162', node2.lemma)
    
    # L -> 10_pl -> 152 -> L
    feats[1175] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_pl', '152') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '152', node2.lemma)
    
    # L -> -139 -> -210 -> L
    feats[1176] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-210') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-210', node2.lemma)
    
    # L -> -139 -> 2_tp -> L
    feats[1177] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '2_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_tp', node2.lemma)
    
    # L -> -13 -> -64 -> L
    feats[1178] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-13') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-13', '-64') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-64', node2.lemma)
    
    # L -> 10_du -> 79 -> L
    feats[1179] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_du', '79') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '79', node2.lemma)
    
    # L -> 7_sg -> -159 -> L
    feats[1180] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_sg', '-159') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-159', node2.lemma)
    
    # L -> 10_du -> 5_fp -> L
    feats[1181] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_du', '5_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_fp', node2.lemma)
    
    # L -> 10_du -> 8_sp -> L
    feats[1182] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_du', '8_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_sp', node2.lemma)
    
    # L -> -152 -> nom. masc. -> L
    feats[1183] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-152') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-152', 'nom. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. masc.', node2.lemma)
    
    # L -> -33 -> nom. adj. -> L
    feats[1184] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', 'nom. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. adj.', node2.lemma)
    
    # L -> sp -> 130 -> L
    feats[1185] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', '130') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '130', node2.lemma)
    
    # L -> -13 -> -42 -> L
    feats[1186] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-13') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-13', '-42') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-42', node2.lemma)
    
    # L -> -13 -> 8_du -> L
    feats[1187] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-13') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-13', '8_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_du', node2.lemma)
    
    # L -> 5_fp -> 16_du -> L
    feats[1188] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', '16_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_du', node2.lemma)
    
    # L -> 7_sg -> 91 -> L
    feats[1189] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_sg', '91') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    
    # L -> 7_sg -> 7_du -> L
    feats[1190] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_sg', '7_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_du', node2.lemma)
    
    # L -> -33 -> 7_tp -> L
    feats[1191] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '7_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_tp', node2.lemma)
    
    # L -> 7_sg -> 92 -> L
    feats[1192] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_sg', '92') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '92', node2.lemma)
    
    # L -> 7_sg -> 151 -> L
    feats[1193] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_sg', '151') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '151', node2.lemma)
    
    # L -> sp -> nom. neutr. -> L
    feats[1194] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', 'nom. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. neutr.', node2.lemma)
    
    # L -> 7_sg -> -61 -> L
    feats[1195] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_sg', '-61') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-61', node2.lemma)
    
    # L -> -152 -> instr. masc. -> L
    feats[1196] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-152') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-152', 'instr. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. masc.', node2.lemma)
    
    # L -> 7_sg -> 3 -> L
    feats[1197] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_sg', '3') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3', node2.lemma)
    
    # L -> 7_sg -> -28 -> L
    feats[1198] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_sg', '-28') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-28', node2.lemma)
    
    # L -> sp -> -291 -> L
    feats[1199] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', '-291') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-291', node2.lemma)
    
    # L -> sp -> -306 -> L
    feats[1200] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', '-306') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-306', node2.lemma)
    
    # L -> 7_sg -> -242 -> L
    feats[1201] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_sg', '-242') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-242', node2.lemma)
    
    # L -> 7_sg -> -303 -> L
    feats[1202] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_sg', '-303') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-303', node2.lemma)
    
    # L -> 5_fp -> loc. sg. -> L
    feats[1203] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', 'loc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc. sg.', node2.lemma)
    
    # L -> 5_fp -> 140 -> L
    feats[1204] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', '140') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '140', node2.lemma)
    
    # L -> 7_sg -> 2_pl -> L
    feats[1205] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_sg', '2_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_pl', node2.lemma)
    
    # L -> 7_sg -> instr. pl. -> L
    feats[1206] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_sg', 'instr. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. pl.', node2.lemma)
    
    # L -> 7_sg -> 160 -> L
    feats[1207] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_sg', '160') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    
    # L -> 7_sg -> -296 -> L
    feats[1208] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_sg', '-296') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-296', node2.lemma)
    
    # L -> 7_sg -> -142 -> L
    feats[1209] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_sg', '-142') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-142', node2.lemma)
    
    # L -> 7_sg -> -50 -> L
    feats[1210] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_sg', '-50') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-50', node2.lemma)
    
    # L -> 7_sg -> -16 -> L
    feats[1211] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_sg', '-16') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-16', node2.lemma)
    
    # L -> 7_sg -> -123 -> L
    feats[1212] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_sg', '-123') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-123', node2.lemma)
    
    # L -> sp -> 1 -> L
    feats[1213] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', '1') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '1', node2.lemma)
    
    # L -> 7_sg -> 2_tp -> L
    feats[1214] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_sg', '2_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_tp', node2.lemma)
    
    # L -> -13 -> -115 -> L
    feats[1215] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-13') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-13', '-115') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-115', node2.lemma)
    
    # L -> 7_sg -> acc -> L
    feats[1216] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_sg', 'acc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc', node2.lemma)
    
    # L -> 7_sg -> -49 -> L
    feats[1217] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_sg', '-49') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-49', node2.lemma)
    
    # L -> 79 -> -59 -> L
    feats[1218] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '-59') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-59', node2.lemma)
    
    # L -> 79 -> -30 -> L
    feats[1219] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '-30') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-30', node2.lemma)
    
    # L -> 79 -> voc. neutr. -> L
    feats[1220] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', 'voc. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. neutr.', node2.lemma)
    
    # L -> 79 -> 78 -> L
    feats[1221] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '78') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '78', node2.lemma)
    
    # L -> 79 -> -26 -> L
    feats[1222] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '-26') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-26', node2.lemma)
    
    # L -> 79 -> -139 -> L
    feats[1223] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '-139') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-139', node2.lemma)
    
    # L -> sp -> 171 -> L
    feats[1224] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', '171') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '171', node2.lemma)
    
    # L -> 79 -> instr -> L
    feats[1225] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', 'instr') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr', node2.lemma)
    
    # L -> sp -> 12_tp -> L
    feats[1226] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', '12_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_tp', node2.lemma)
    
    # L -> -152 -> 177 -> L
    feats[1227] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-152') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-152', '177') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '177', node2.lemma)
    
    # L -> 79 -> 75 -> L
    feats[1228] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '75') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '75', node2.lemma)
    
    # L -> -33 -> 10_fp -> L
    feats[1229] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '10_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_fp', node2.lemma)
    
    # L -> -152 -> -49 -> L
    feats[1230] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-152') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-152', '-49') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-49', node2.lemma)
    
    # L -> 79 -> 9_sg -> L
    feats[1231] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '9_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_sg', node2.lemma)
    
    # L -> 79 -> -81 -> L
    feats[1232] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '-81') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-81', node2.lemma)
    
    # L -> 79 -> -13 -> L
    feats[1233] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '-13') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-13', node2.lemma)
    
    # L -> 79 -> -132 -> L
    feats[1234] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '-132') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-132', node2.lemma)
    
    # L -> 79 -> nom. du. -> L
    feats[1235] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', 'nom. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. du.', node2.lemma)
    
    # L -> sp -> sg_fp -> L
    feats[1236] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', 'sg_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_fp', node2.lemma)
    
    # L -> 79 -> 30_fp -> L
    feats[1237] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '30_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_fp', node2.lemma)
    
    # L -> 79 -> dat. pl. -> L
    feats[1238] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', 'dat. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. pl.', node2.lemma)
    
    # L -> 79 -> -276 -> L
    feats[1239] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '-276') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-276', node2.lemma)
    
    # L -> 79 -> du_fp -> L
    feats[1240] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', 'du_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_fp', node2.lemma)
    
    # L -> 79 -> 5_du -> L
    feats[1241] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '5_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_du', node2.lemma)
    
    # L -> sp -> 5_pl -> L
    feats[1242] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', '5_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_pl', node2.lemma)
    
    # L -> sp -> -163 -> L
    feats[1243] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', '-163') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-163', node2.lemma)
    
    # L -> 79 -> 88 -> L
    feats[1244] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '88') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '88', node2.lemma)
    
    # L -> 79 -> 12_pl -> L
    feats[1245] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '12_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_pl', node2.lemma)
    
    # L -> 79 -> 90 -> L
    feats[1246] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '90') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '90', node2.lemma)
    
    # L -> 79 -> -87 -> L
    feats[1247] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '-87') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-87', node2.lemma)
    
    # L -> 79 -> 11_tp -> L
    feats[1248] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '11_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_tp', node2.lemma)
    
    # L -> 79 -> 74 -> L
    feats[1249] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '74') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '74', node2.lemma)
    
    # L -> 30 -> 72 -> L
    feats[1250] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30', '72') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '72', node2.lemma)
    
    # L -> 79 -> -44 -> L
    feats[1251] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '-44') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-44', node2.lemma)
    
    # L -> 79 -> du_tp -> L
    feats[1252] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', 'du_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_tp', node2.lemma)
    
    # L -> 79 -> -84 -> L
    feats[1253] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '-84') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-84', node2.lemma)
    
    # L -> sp -> -34 -> L
    feats[1254] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', '-34') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-34', node2.lemma)
    
    # L -> 79 -> 13_sp -> L
    feats[1255] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '13_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_sp', node2.lemma)
    
    # L -> 79 -> -241 -> L
    feats[1256] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '-241') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-241', node2.lemma)
    
    # L -> 79 -> -72 -> L
    feats[1257] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '-72') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-72', node2.lemma)
    
    # L -> 79 -> 2_sp -> L
    feats[1258] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '2_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_sp', node2.lemma)
    
    # L -> 79 -> -45 -> L
    feats[1259] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '-45') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-45', node2.lemma)
    
    # L -> -33 -> -269 -> L
    feats[1260] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-269') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-269', node2.lemma)
    
    # L -> -33 -> voc. pl. -> L
    feats[1261] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', 'voc. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. pl.', node2.lemma)
    
    # L -> sp -> 119 -> L
    feats[1262] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', '119') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '119', node2.lemma)
    
    # L -> sp -> -273 -> L
    feats[1263] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', '-273') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-273', node2.lemma)
    
    # L -> -33 -> voc. masc. -> L
    feats[1264] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', 'voc. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. masc.', node2.lemma)
    
    # L -> -33 -> 158 -> L
    feats[1265] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '158') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '158', node2.lemma)
    
    # L -> -33 -> 16_sg -> L
    feats[1266] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '16_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_sg', node2.lemma)
    
    # L -> sp -> -62 -> L
    feats[1267] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', '-62') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-62', node2.lemma)
    
    # L -> sp -> -37 -> L
    feats[1268] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', '-37') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-37', node2.lemma)
    
    # L -> -33 -> 4_sg -> L
    feats[1269] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '4_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_sg', node2.lemma)
    
    # L -> -33 -> 119 -> L
    feats[1270] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '119') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '119', node2.lemma)
    
    # L -> -33 -> 82 -> L
    feats[1271] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '82') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '82', node2.lemma)
    
    # L -> -33 -> 139 -> L
    feats[1272] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '139') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '139', node2.lemma)
    
    # L -> 5_fp -> 11_sp -> L
    feats[1273] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', '11_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_sp', node2.lemma)
    
    # L -> -33 -> -37 -> L
    feats[1274] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-37') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-37', node2.lemma)
    
    # L -> -33 -> -119 -> L
    feats[1275] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-119') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-119', node2.lemma)
    
    # L -> -33 -> abl. du. -> L
    feats[1276] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', 'abl. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. du.', node2.lemma)
    
    # L -> -33 -> -46 -> L
    feats[1277] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-46') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-46', node2.lemma)
    
    # L -> -33 -> 115 -> L
    feats[1278] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '115') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '115', node2.lemma)
    
    # L -> -33 -> abl. pl. -> L
    feats[1279] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', 'abl. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. pl.', node2.lemma)
    
    # L -> -33 -> 92 -> L
    feats[1280] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '92') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '92', node2.lemma)
    
    # L -> -33 -> -11 -> L
    feats[1281] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-11') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-11', node2.lemma)
    
    # L -> -33 -> 151 -> L
    feats[1282] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '151') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '151', node2.lemma)
    
    # L -> -33 -> 162 -> L
    feats[1283] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '162') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '162', node2.lemma)
    
    # L -> -33 -> -141 -> L
    feats[1284] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-141') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-141', node2.lemma)
    
    # L -> -33 -> 51 -> L
    feats[1285] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '51') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '51', node2.lemma)
    
    # L -> -33 -> 11_sp -> L
    feats[1286] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '11_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_sp', node2.lemma)
    
    # L -> 30 -> -158 -> L
    feats[1287] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30', '-158') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-158', node2.lemma)
    
    # L -> -33 -> -89 -> L
    feats[1288] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-89') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-89', node2.lemma)
    
    # L -> -33 -> -28 -> L
    feats[1289] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-28') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-28', node2.lemma)
    
    # L -> -33 -> -150 -> L
    feats[1290] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-150') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-150', node2.lemma)
    
    # L -> -33 -> -242 -> L
    feats[1291] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-242') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-242', node2.lemma)
    
    # L -> -33 -> -303 -> L
    feats[1292] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-303') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-303', node2.lemma)
    
    # L -> -33 -> 7_fp -> L
    feats[1293] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '7_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_fp', node2.lemma)
    
    # L -> -33 -> -38 -> L
    feats[1294] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-38') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-38', node2.lemma)
    
    # L -> -33 -> 174 -> L
    feats[1295] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '174') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '174', node2.lemma)
    
    # L -> -33 -> 5_sg -> L
    feats[1296] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '5_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sg', node2.lemma)
    
    # L -> -33 -> 176 -> L
    feats[1297] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '176') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '176', node2.lemma)
    
    # L -> -33 -> 160 -> L
    feats[1298] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '160') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    
    # L -> -33 -> -67 -> L
    feats[1299] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-67') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-67', node2.lemma)
    
    # L -> -33 -> fp -> L
    feats[1300] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', 'fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'fp', node2.lemma)
    
    # L -> -33 -> -16 -> L
    feats[1301] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-16') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-16', node2.lemma)
    
    # L -> -33 -> -115 -> L
    feats[1302] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-115') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-115', node2.lemma)
    
    # L -> -33 -> -210 -> L
    feats[1303] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-210') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-210', node2.lemma)
    
    # L -> -33 -> 2_tp -> L
    feats[1304] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '2_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_tp', node2.lemma)
    
    # L -> -33 -> acc -> L
    feats[1305] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', 'acc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc', node2.lemma)
    
    # L -> -33 -> -49 -> L
    feats[1306] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-49') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-49', node2.lemma)
    
    # L -> 28 -> 27_du -> L
    feats[1307] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28', '27_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_du', node2.lemma)
    
    # L -> 28 -> 7_sg -> L
    feats[1308] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28', '7_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_sg', node2.lemma)
    
    # L -> 28 -> sg_tp -> L
    feats[1309] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28', 'sg_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_tp', node2.lemma)
    
    # L -> 28 -> voc. sg. -> L
    feats[1310] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28', 'voc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. sg.', node2.lemma)
    
    # L -> 28 -> 78 -> L
    feats[1311] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28', '78') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '78', node2.lemma)
    
    # L -> 28 -> -139 -> L
    feats[1312] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28', '-139') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-139', node2.lemma)
    
    # L -> 8_sp -> -111 -> L
    feats[1313] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_sp', '-111') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-111', node2.lemma)
    
    # L -> 30 -> pl_fp -> L
    feats[1314] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30', 'pl_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_fp', node2.lemma)
    
    # L -> 30 -> 4_du -> L
    feats[1315] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30', '4_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_du', node2.lemma)
    
    # L -> 28 -> 2_du -> L
    feats[1316] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28', '2_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_du', node2.lemma)
    
    # L -> 28 -> 71 -> L
    feats[1317] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28', '71') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '71', node2.lemma)
    
    # L -> 8_sp -> 159 -> L
    feats[1318] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_sp', '159') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '159', node2.lemma)
    
    # L -> 28 -> 169 -> L
    feats[1319] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28', '169') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '169', node2.lemma)
    
    # L -> 28 -> sg -> L
    feats[1320] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28', 'sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg', node2.lemma)
    
    # L -> 28 -> loc. du. -> L
    feats[1321] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28', 'loc. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc. du.', node2.lemma)
    
    # L -> 14_sg -> 14_tp -> L
    feats[1322] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_sg', '14_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_tp', node2.lemma)
    
    # L -> 28 -> -13 -> L
    feats[1323] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28', '-13') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-13', node2.lemma)
    
    # L -> 28 -> 14_fp -> L
    feats[1324] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28', '14_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_fp', node2.lemma)
    
    # L -> 28 -> 12_sg -> L
    feats[1325] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28', '12_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_sg', node2.lemma)
    
    # L -> 28 -> 31 -> L
    feats[1326] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28', '31') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '31', node2.lemma)
    
    # L -> 28 -> -261 -> L
    feats[1327] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28', '-261') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-261', node2.lemma)
    
    # L -> 28 -> dat. pl. -> L
    feats[1328] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28', 'dat. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. pl.', node2.lemma)
    
    # L -> 79 -> 6_du -> L
    feats[1329] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '6_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_du', node2.lemma)
    
    # L -> 28 -> -263 -> L
    feats[1330] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28', '-263') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-263', node2.lemma)
    
    # L -> 28 -> -276 -> L
    feats[1331] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28', '-276') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-276', node2.lemma)
    
    # L -> 79 -> nom. pl. -> L
    feats[1332] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', 'nom. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. pl.', node2.lemma)
    
    # L -> 28 -> sp -> L
    feats[1333] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28', 'sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sp', node2.lemma)
    
    # L -> 14_sg -> 111 -> L
    feats[1334] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_sg', '111') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '111', node2.lemma)
    
    # L -> 79 -> 10_sg -> L
    feats[1335] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '10_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_sg', node2.lemma)
    
    # L -> 14_sg -> 154 -> L
    feats[1336] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_sg', '154') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '154', node2.lemma)
    
    # L -> 79 -> 158 -> L
    feats[1337] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '158') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '158', node2.lemma)
    
    # L -> 79 -> acc. fem -> L
    feats[1338] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', 'acc. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. fem', node2.lemma)
    
    # L -> 79 -> 11_du -> L
    feats[1339] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '11_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_du', node2.lemma)
    
    # L -> 79 -> -153 -> L
    feats[1340] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '-153') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-153', node2.lemma)
    
    # L -> 14_sg -> -55 -> L
    feats[1341] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_sg', '-55') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-55', node2.lemma)
    
    # L -> 79 -> -246 -> L
    feats[1342] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '-246') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-246', node2.lemma)
    
    # L -> -52 -> 7_tp -> L
    feats[1343] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-52') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-52', '7_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_tp', node2.lemma)
    
    # L -> -52 -> -143 -> L
    feats[1344] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-52') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-52', '-143') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-143', node2.lemma)
    
    # L -> 14_sg -> voc. masc. -> L
    feats[1345] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_sg', 'voc. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. masc.', node2.lemma)
    
    # L -> 14_sg -> 28_sg -> L
    feats[1346] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_sg', '28_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28_sg', node2.lemma)
    
    # L -> 14_sg -> -35 -> L
    feats[1347] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_sg', '-35') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    
    # L -> 79 -> -119 -> L
    feats[1348] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '-119') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-119', node2.lemma)
    
    # L -> 79 -> -14 -> L
    feats[1349] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '-14') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-14', node2.lemma)
    
    # L -> 14_sg -> 11_du -> L
    feats[1350] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_sg', '11_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_du', node2.lemma)
    
    # L -> 79 -> 129 -> L
    feats[1351] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '129') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '129', node2.lemma)
    
    # L -> 79 -> -61 -> L
    feats[1352] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '-61') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-61', node2.lemma)
    
    # L -> 79 -> -156 -> L
    feats[1353] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '-156') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-156', node2.lemma)
    
    # L -> 14_sg -> -119 -> L
    feats[1354] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_sg', '-119') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-119', node2.lemma)
    
    # L -> 8_sp -> 8_sg -> L
    feats[1355] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_sp', '8_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_sg', node2.lemma)
    
    # L -> 30 -> 11_sp -> L
    feats[1356] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30', '11_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_sp', node2.lemma)
    
    # L -> 79 -> -150 -> L
    feats[1357] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '-150') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-150', node2.lemma)
    
    # L -> 79 -> 136 -> L
    feats[1358] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '136') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '136', node2.lemma)
    
    # L -> 14_sg -> -11 -> L
    feats[1359] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_sg', '-11') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-11', node2.lemma)
    
    # L -> 79 -> instr. du. -> L
    feats[1360] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', 'instr. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. du.', node2.lemma)
    
    # L -> 14_sg -> -141 -> L
    feats[1361] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_sg', '-141') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-141', node2.lemma)
    
    # L -> 14_sg -> 51 -> L
    feats[1362] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_sg', '51') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '51', node2.lemma)
    
    # L -> 79 -> -38 -> L
    feats[1363] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '-38') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-38', node2.lemma)
    
    # L -> 79 -> -293 -> L
    feats[1364] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '-293') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-293', node2.lemma)
    
    # L -> -52 -> -109 -> L
    feats[1365] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-52') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-52', '-109') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-109', node2.lemma)
    
    # L -> -52 -> 6_pl -> L
    feats[1366] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-52') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-52', '6_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_pl', node2.lemma)
    
    # L -> 8_sp -> -166 -> L
    feats[1367] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_sp', '-166') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-166', node2.lemma)
    
    # L -> 8_sp -> -154 -> L
    feats[1368] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_sp', '-154') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-154', node2.lemma)
    
    # L -> 28 -> 89 -> L
    feats[1369] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28', '89') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '89', node2.lemma)
    
    # L -> -98 -> 2_du -> L
    feats[1370] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-98', '2_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_du', node2.lemma)
    
    # L -> -98 -> 71 -> L
    feats[1371] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-98', '71') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '71', node2.lemma)
    
    # L -> sg_tp -> -169 -> L
    feats[1372] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg_tp', '-169') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-169', node2.lemma)
    
    # L -> -52 -> nom. fem -> L
    feats[1373] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-52') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-52', 'nom. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. fem', node2.lemma)
    
    # L -> -52 -> -55 -> L
    feats[1374] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-52') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-52', '-55') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-55', node2.lemma)
    
    # L -> -15 -> -299 -> L
    feats[1375] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-15') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-15', '-299') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-299', node2.lemma)
    
    # L -> -15 -> sg_sp -> L
    feats[1376] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-15') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-15', 'sg_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_sp', node2.lemma)
    
    # L -> 8_sp -> -14 -> L
    feats[1377] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_sp', '-14') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-14', node2.lemma)
    
    # L -> -98 -> -245 -> L
    feats[1378] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-98', '-245') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-245', node2.lemma)
    
    # L -> sg_tp -> -291 -> L
    feats[1379] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg_tp', '-291') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-291', node2.lemma)
    
    # L -> sg_tp -> 70 -> L
    feats[1380] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg_tp', '70') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '70', node2.lemma)
    
    # L -> -15 -> 51 -> L
    feats[1381] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-15') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-15', '51') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '51', node2.lemma)
    
    # L -> -15 -> 27_sg -> L
    feats[1382] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-15') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-15', '27_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_sg', node2.lemma)
    
    # L -> -52 -> 73 -> L
    feats[1383] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-52') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-52', '73') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '73', node2.lemma)
    
    # L -> 28 -> -293 -> L
    feats[1384] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28', '-293') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-293', node2.lemma)
    
    # L -> -15 -> 176 -> L
    feats[1385] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-15') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-15', '176') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '176', node2.lemma)
    
    # L -> -15 -> 160 -> L
    feats[1386] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-15') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-15', '160') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    
    # L -> -36 -> sg_tp -> L
    feats[1387] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-36') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-36', 'sg_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_tp', node2.lemma)
    
    # L -> -36 -> -200 -> L
    feats[1388] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-36') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-36', '-200') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-200', node2.lemma)
    
    # L -> -98 -> 171 -> L
    feats[1389] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-98', '171') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '171', node2.lemma)
    
    # L -> -31 -> -97 -> L
    feats[1390] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-31') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-31', '-97') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-97', node2.lemma)
    
    # L -> -31 -> 27_fp -> L
    feats[1391] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-31') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-31', '27_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_fp', node2.lemma)
    
    # L -> -31 -> 155 -> L
    feats[1392] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-31') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-31', '155') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '155', node2.lemma)
    
    # L -> sg_tp -> gen. du. -> L
    feats[1393] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg_tp', 'gen. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen. du.', node2.lemma)
    
    # L -> -31 -> instr -> L
    feats[1394] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-31') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-31', 'instr') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr', node2.lemma)
    
    # L -> 131 -> 14_pl -> L
    feats[1395] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '131') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '131', '14_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_pl', node2.lemma)
    
    # L -> 131 -> -69 -> L
    feats[1396] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '131') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '131', '-69') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-69', node2.lemma)
    
    # L -> -36 -> -68 -> L
    feats[1397] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-36') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-36', '-68') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-68', node2.lemma)
    
    # L -> 38 -> sg_sp -> L
    feats[1398] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '38') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '38', 'sg_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_sp', node2.lemma)
    
    # L -> sg_tp -> 6_sp -> L
    feats[1399] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg_tp', '6_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_sp', node2.lemma)
    
    # L -> sg_tp -> 15_fp -> L
    feats[1400] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg_tp', '15_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_fp', node2.lemma)
    
    # L -> 38 -> -35 -> L
    feats[1401] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '38') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '38', '-35') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    
    # L -> 108 -> -268 -> L
    feats[1402] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '108') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '108', '-268') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-268', node2.lemma)
    
    # L -> 131 -> -73 -> L
    feats[1403] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '131') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '131', '-73') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-73', node2.lemma)
    
    # L -> -36 -> -19 -> L
    feats[1404] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-36') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-36', '-19') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-19', node2.lemma)
    
    # L -> -36 -> 128 -> L
    feats[1405] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-36') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-36', '128') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '128', node2.lemma)
    
    # L -> -98 -> -67 -> L
    feats[1406] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-98', '-67') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-67', node2.lemma)
    
    # L -> 38 -> 2_pl -> L
    feats[1407] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '38') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '38', '2_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_pl', node2.lemma)
    
    # L -> 38 -> -121 -> L
    feats[1408] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '38') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '38', '-121') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-121', node2.lemma)
    
    # L -> 131 -> 4_sp -> L
    feats[1409] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '131') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '131', '4_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_sp', node2.lemma)
    
    # L -> -36 -> -159 -> L
    feats[1410] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-36') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-36', '-159') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-159', node2.lemma)
    
    # L -> 131 -> -129 -> L
    feats[1411] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '131') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '131', '-129') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-129', node2.lemma)
    
    # L -> 131 -> 13_sg -> L
    feats[1412] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '131') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '131', '13_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_sg', node2.lemma)
    
    # L -> 108 -> -166 -> L
    feats[1413] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '108') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '108', '-166') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-166', node2.lemma)
    
    # L -> -200 -> -90 -> L
    feats[1414] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '-90') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-90', node2.lemma)
    
    # L -> 7_pl -> 130 -> L
    feats[1415] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '130') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '130', node2.lemma)
    
    # L -> -200 -> 27_pl -> L
    feats[1416] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '27_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_pl', node2.lemma)
    
    # L -> 7_pl -> -169 -> L
    feats[1417] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '-169') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-169', node2.lemma)
    
    # L -> -200 -> -44 -> L
    feats[1418] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '-44') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-44', node2.lemma)
    
    # L -> -200 -> voc. du. -> L
    feats[1419] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', 'voc. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. du.', node2.lemma)
    
    # L -> -200 -> -301 -> L
    feats[1420] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '-301') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-301', node2.lemma)
    
    # L -> -200 -> -58 -> L
    feats[1421] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '-58') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-58', node2.lemma)
    
    # L -> -36 -> voc. masc. -> L
    feats[1422] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-36') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-36', 'voc. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. masc.', node2.lemma)
    
    # L -> -200 -> -45 -> L
    feats[1423] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '-45') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-45', node2.lemma)
    
    # L -> -36 -> 158 -> L
    feats[1424] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-36') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-36', '158') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '158', node2.lemma)
    
    # L -> -200 -> 59 -> L
    feats[1425] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '59') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '59', node2.lemma)
    
    # L -> -200 -> -122 -> L
    feats[1426] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '-122') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-122', node2.lemma)
    
    # L -> -200 -> 101 -> L
    feats[1427] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '101') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    
    # L -> -200 -> abl. sg. -> L
    feats[1428] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', 'abl. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. sg.', node2.lemma)
    
    # L -> -200 -> -99 -> L
    feats[1429] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '-99') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-99', node2.lemma)
    
    # L -> -200 -> 10_pl -> L
    feats[1430] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '10_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_pl', node2.lemma)
    
    # L -> -200 -> -15 -> L
    feats[1431] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '-15') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-15', node2.lemma)
    
    # L -> -200 -> -137 -> L
    feats[1432] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '-137') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-137', node2.lemma)
    
    # L -> -200 -> -151 -> L
    feats[1433] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '-151') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-151', node2.lemma)
    
    # L -> -36 -> 15_fp -> L
    feats[1434] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-36') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-36', '15_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_fp', node2.lemma)
    
    # L -> -36 -> -119 -> L
    feats[1435] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-36') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-36', '-119') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-119', node2.lemma)
    
    # L -> -200 -> -245 -> L
    feats[1436] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '-245') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-245', node2.lemma)
    
    # L -> -200 -> -157 -> L
    feats[1437] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '-157') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-157', node2.lemma)
    
    # L -> 108 -> pl_tp -> L
    feats[1438] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '108') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '108', 'pl_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_tp', node2.lemma)
    
    # L -> -200 -> nom. neutr. -> L
    feats[1439] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', 'nom. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. neutr.', node2.lemma)
    
    # L -> 108 -> -153 -> L
    feats[1440] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '108') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '108', '-153') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-153', node2.lemma)
    
    # L -> -200 -> -306 -> L
    feats[1441] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '-306') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-306', node2.lemma)
    
    # L -> -200 -> -291 -> L
    feats[1442] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '-291') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-291', node2.lemma)
    
    # L -> -200 -> 70 -> L
    feats[1443] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '70') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '70', node2.lemma)
    
    # L -> -200 -> -48 -> L
    feats[1444] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '-48') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-48', node2.lemma)
    
    # L -> -200 -> 8_sg -> L
    feats[1445] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '8_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_sg', node2.lemma)
    
    # L -> -200 -> -33 -> L
    feats[1446] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '-33') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-33', node2.lemma)
    
    # L -> -200 -> 33 -> L
    feats[1447] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '33') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '33', node2.lemma)
    
    # L -> -200 -> -24 -> L
    feats[1448] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '-24') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    
    # L -> -200 -> 110 -> L
    feats[1449] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '110') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '110', node2.lemma)
    
    # L -> -200 -> 41 -> L
    feats[1450] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '41') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '41', node2.lemma)
    
    # L -> -200 -> -147 -> L
    feats[1451] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '-147') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-147', node2.lemma)
    
    # L -> -200 -> 140 -> L
    feats[1452] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '140') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '140', node2.lemma)
    
    # L -> -200 -> 81 -> L
    feats[1453] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '81') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '81', node2.lemma)
    
    # L -> -200 -> -32 -> L
    feats[1454] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '-32') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-32', node2.lemma)
    
    # L -> -200 -> 4_fp -> L
    feats[1455] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '4_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_fp', node2.lemma)
    
    # L -> -200 -> -19 -> L
    feats[1456] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '-19') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-19', node2.lemma)
    
    # L -> -200 -> -283 -> L
    feats[1457] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '-283') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-283', node2.lemma)
    
    # L -> -200 -> 68 -> L
    feats[1458] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '68') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '68', node2.lemma)
    
    # L -> -200 -> -114 -> L
    feats[1459] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '-114') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-114', node2.lemma)
    
    # L -> 131 -> 160 -> L
    feats[1460] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '131') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '131', '160') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    
    # L -> 131 -> 3_sp -> L
    feats[1461] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '131') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '131', '3_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_sp', node2.lemma)
    
    # L -> 131 -> -16 -> L
    feats[1462] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '131') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '131', '-16') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-16', node2.lemma)
    
    # L -> -200 -> -109 -> L
    feats[1463] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '-109') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-109', node2.lemma)
    
    # L -> 131 -> acc -> L
    feats[1464] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '131') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '131', 'acc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc', node2.lemma)
    
    # L -> 7_pl -> -109 -> L
    feats[1465] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '-109') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-109', node2.lemma)
    
    # L -> 159 -> 99 -> L
    feats[1466] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '99') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '99', node2.lemma)
    
    # L -> -200 -> 171 -> L
    feats[1467] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '171') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '171', node2.lemma)
    
    # L -> -200 -> 3_du -> L
    feats[1468] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '3_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    
    # L -> -200 -> pl_sp -> L
    feats[1469] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', 'pl_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_sp', node2.lemma)
    
    # L -> -200 -> 9_du -> L
    feats[1470] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '9_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_du', node2.lemma)
    
    # L -> -200 -> 12_tp -> L
    feats[1471] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '12_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_tp', node2.lemma)
    
    # L -> -200 -> 13_fp -> L
    feats[1472] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '13_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_fp', node2.lemma)
    
    # L -> 7_pl -> 171 -> L
    feats[1473] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '171') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '171', node2.lemma)
    
    # L -> -200 -> acc. pl. -> L
    feats[1474] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', 'acc. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. pl.', node2.lemma)
    
    # L -> 7_pl -> pl_sp -> L
    feats[1475] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', 'pl_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_sp', node2.lemma)
    
    # L -> 159 -> -93 -> L
    feats[1476] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '-93') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-93', node2.lemma)
    
    # L -> 7_pl -> 13_fp -> L
    feats[1477] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '13_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_fp', node2.lemma)
    
    # L -> -200 -> 8_pl -> L
    feats[1478] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '8_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_pl', node2.lemma)
    
    # L -> 7_pl -> 6_pl -> L
    feats[1479] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '6_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_pl', node2.lemma)
    
    # L -> 7_pl -> 114 -> L
    feats[1480] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '114') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '114', node2.lemma)
    
    # L -> -200 -> 175 -> L
    feats[1481] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '175') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '175', node2.lemma)
    
    # L -> -200 -> 153 -> L
    feats[1482] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '153') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    
    # L -> 7_pl -> -41 -> L
    feats[1483] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '-41') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-41', node2.lemma)
    
    # L -> 159 -> 5_fp -> L
    feats[1484] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '5_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_fp', node2.lemma)
    
    # L -> 7_pl -> nom. masc. -> L
    feats[1485] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', 'nom. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. masc.', node2.lemma)
    
    # L -> 7_pl -> 3_pl -> L
    feats[1486] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '3_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_pl', node2.lemma)
    
    # L -> 155 -> 75 -> L
    feats[1487] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '155') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '155', '75') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '75', node2.lemma)
    
    # L -> -56 -> -97 -> L
    feats[1488] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-56') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-56', '-97') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-97', node2.lemma)
    
    # L -> 159 -> nom. adj. -> L
    feats[1489] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', 'nom. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. adj.', node2.lemma)
    
    # L -> 159 -> 30_du -> L
    feats[1490] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '30_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_du', node2.lemma)
    
    # L -> 159 -> 2_du -> L
    feats[1491] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '2_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_du', node2.lemma)
    
    # L -> 7_pl -> 30_tp -> L
    feats[1492] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '30_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_tp', node2.lemma)
    
    # L -> -137 -> 5_fp -> L
    feats[1493] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-137') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-137', '5_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_fp', node2.lemma)
    
    # L -> 7_pl -> -43 -> L
    feats[1494] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '-43') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-43', node2.lemma)
    
    # L -> 159 -> -98 -> L
    feats[1495] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '-98') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-98', node2.lemma)
    
    # L -> 7_pl -> -129 -> L
    feats[1496] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '-129') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-129', node2.lemma)
    
    # L -> 159 -> 12_sg -> L
    feats[1497] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '12_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_sg', node2.lemma)
    
    # L -> 7_pl -> -163 -> L
    feats[1498] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '-163') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-163', node2.lemma)
    
    # L -> 159 -> -261 -> L
    feats[1499] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '-261') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-261', node2.lemma)
    
    # L -> 7_pl -> 4_du -> L
    feats[1500] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '4_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_du', node2.lemma)
    
    # L -> 7_pl -> 2 -> L
    feats[1501] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '2') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2', node2.lemma)
    
    # L -> 7_pl -> -271 -> L
    feats[1502] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '-271') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-271', node2.lemma)
    
    # L -> 137 -> -154 -> L
    feats[1503] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '137') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '137', '-154') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-154', node2.lemma)
    
    # L -> 7_pl -> 135 -> L
    feats[1504] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '135') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '135', node2.lemma)
    
    # L -> 7_pl -> 27_tp -> L
    feats[1505] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '27_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_tp', node2.lemma)
    
    # L -> 7_pl -> tp -> L
    feats[1506] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', 'tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'tp', node2.lemma)
    
    # L -> 7_pl -> dat. du. -> L
    feats[1507] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', 'dat. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    
    # L -> 7_pl -> 36 -> L
    feats[1508] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '36') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '36', node2.lemma)
    
    # L -> 7_pl -> 16_fp -> L
    feats[1509] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '16_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_fp', node2.lemma)
    
    # L -> 7_pl -> 142 -> L
    feats[1510] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '142') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '142', node2.lemma)
    
    # L -> 7_pl -> 111 -> L
    feats[1511] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '111') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '111', node2.lemma)
    
    # L -> 7_pl -> -54 -> L
    feats[1512] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '-54') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-54', node2.lemma)
    
    # L -> 7_pl -> 141 -> L
    feats[1513] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '141') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '141', node2.lemma)
    
    # L -> 7_pl -> 8_fp -> L
    feats[1514] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '8_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    
    # L -> 7_pl -> -113 -> L
    feats[1515] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '-113') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-113', node2.lemma)
    
    # L -> 7_pl -> -299 -> L
    feats[1516] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '-299') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-299', node2.lemma)
    
    # L -> 7_pl -> -269 -> L
    feats[1517] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '-269') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-269', node2.lemma)
    
    # L -> 7_pl -> nom. pl. -> L
    feats[1518] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', 'nom. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. pl.', node2.lemma)
    
    # L -> 7_pl -> nom. fem -> L
    feats[1519] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', 'nom. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. fem', node2.lemma)
    
    # L -> 7_pl -> -55 -> L
    feats[1520] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '-55') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-55', node2.lemma)
    
    # L -> 7_pl -> instr. sg. -> L
    feats[1521] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', 'instr. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. sg.', node2.lemma)
    
    # L -> 137 -> 154 -> L
    feats[1522] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '137') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '137', '154') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '154', node2.lemma)
    
    # L -> 7_pl -> -94 -> L
    feats[1523] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '-94') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-94', node2.lemma)
    
    # L -> 7_pl -> 3_fp -> L
    feats[1524] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '3_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_fp', node2.lemma)
    
    # L -> 7_pl -> 8_du -> L
    feats[1525] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '8_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_du', node2.lemma)
    
    # L -> 7_pl -> acc. fem -> L
    feats[1526] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', 'acc. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. fem', node2.lemma)
    
    # L -> 7_pl -> sg_sp -> L
    feats[1527] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', 'sg_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_sp', node2.lemma)
    
    # L -> 7_pl -> -35 -> L
    feats[1528] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '-35') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    
    # L -> 7_pl -> 4_sg -> L
    feats[1529] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '4_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_sg', node2.lemma)
    
    # L -> 7_pl -> pl_tp -> L
    feats[1530] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', 'pl_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_tp', node2.lemma)
    
    # L -> 7_pl -> -153 -> L
    feats[1531] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '-153') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-153', node2.lemma)
    
    # L -> 155 -> 15_tp -> L
    feats[1532] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '155') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '155', '15_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_tp', node2.lemma)
    
    # L -> 7_pl -> -273 -> L
    feats[1533] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '-273') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-273', node2.lemma)
    
    # L -> 7_pl -> dat. sg. -> L
    feats[1534] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', 'dat. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. sg.', node2.lemma)
    
    # L -> -56 -> -58 -> L
    feats[1535] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-56') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-56', '-58') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-58', node2.lemma)
    
    # L -> -56 -> 13_sp -> L
    feats[1536] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-56') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-56', '13_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_sp', node2.lemma)
    
    # L -> 155 -> -122 -> L
    feats[1537] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '155') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '155', '-122') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-122', node2.lemma)
    
    # L -> 7_pl -> 82 -> L
    feats[1538] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '82') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '82', node2.lemma)
    
    # L -> 137 -> voc. masc. -> L
    feats[1539] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '137') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '137', 'voc. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. masc.', node2.lemma)
    
    # L -> 7_pl -> 8_tp -> L
    feats[1540] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '8_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_tp', node2.lemma)
    
    # L -> 7_pl -> 139 -> L
    feats[1541] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '139') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '139', node2.lemma)
    
    # L -> 7_pl -> 7_du -> L
    feats[1542] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '7_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_du', node2.lemma)
    
    # L -> 7_pl -> 122 -> L
    feats[1543] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '122') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '122', node2.lemma)
    
    # L -> 7_pl -> -62 -> L
    feats[1544] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '-62') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-62', node2.lemma)
    
    # L -> 7_pl -> -37 -> L
    feats[1545] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '-37') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-37', node2.lemma)
    
    # L -> 7_pl -> 6_sp -> L
    feats[1546] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '6_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_sp', node2.lemma)
    
    # L -> 7_pl -> -14 -> L
    feats[1547] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '-14') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-14', node2.lemma)
    
    # L -> 137 -> dat. sg. -> L
    feats[1548] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '137') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '137', 'dat. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. sg.', node2.lemma)
    
    # L -> 137 -> -246 -> L
    feats[1549] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '137') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '137', '-246') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-246', node2.lemma)
    
    # L -> 7_pl -> -46 -> L
    feats[1550] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '-46') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-46', node2.lemma)
    
    # L -> 7_pl -> -20 -> L
    feats[1551] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '-20') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-20', node2.lemma)
    
    # L -> 7_pl -> 151 -> L
    feats[1552] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '151') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '151', node2.lemma)
    
    # L -> 7_pl -> 162 -> L
    feats[1553] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '162') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '162', node2.lemma)
    
    # L -> 7_pl -> 11_sp -> L
    feats[1554] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '11_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_sp', node2.lemma)
    
    # L -> 7_pl -> -297 -> L
    feats[1555] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '-297') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-297', node2.lemma)
    
    # L -> 7_pl -> 3 -> L
    feats[1556] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '3') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3', node2.lemma)
    
    # L -> 7_pl -> -104 -> L
    feats[1557] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '-104') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-104', node2.lemma)
    
    # L -> 7_pl -> -28 -> L
    feats[1558] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '-28') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-28', node2.lemma)
    
    # L -> 7_pl -> 102 -> L
    feats[1559] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '102') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '102', node2.lemma)
    
    # L -> 7_pl -> 11_fp -> L
    feats[1560] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '11_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_fp', node2.lemma)
    
    # L -> 7_pl -> -150 -> L
    feats[1561] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '-150') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-150', node2.lemma)
    
    # L -> 7_pl -> -242 -> L
    feats[1562] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '-242') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-242', node2.lemma)
    
    # L -> 7_pl -> -303 -> L
    feats[1563] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '-303') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-303', node2.lemma)
    
    # L -> 7_pl -> instr. du. -> L
    feats[1564] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', 'instr. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. du.', node2.lemma)
    
    # L -> 7_pl -> -25 -> L
    feats[1565] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '-25') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-25', node2.lemma)
    
    # L -> 7_pl -> -103 -> L
    feats[1566] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '-103') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-103', node2.lemma)
    
    # L -> 7_pl -> -38 -> L
    feats[1567] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '-38') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-38', node2.lemma)
    
    # L -> 7_pl -> -121 -> L
    feats[1568] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '-121') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-121', node2.lemma)
    
    # L -> 159 -> 15_sg -> L
    feats[1569] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '15_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_sg', node2.lemma)
    
    # L -> 159 -> 121 -> L
    feats[1570] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '121') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '121', node2.lemma)
    
    # L -> 7_pl -> dat -> L
    feats[1571] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', 'dat') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat', node2.lemma)
    
    # L -> 159 -> -19 -> L
    feats[1572] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '-19') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-19', node2.lemma)
    
    # L -> 159 -> 128 -> L
    feats[1573] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '128') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '128', node2.lemma)
    
    # L -> -200 -> 2_tp -> L
    feats[1574] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '2_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_tp', node2.lemma)
    
    # L -> -200 -> 6_tp -> L
    feats[1575] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', '6_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_tp', node2.lemma)
    
    # L -> 7_pl -> -115 -> L
    feats[1576] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', '-115') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-115', node2.lemma)
    
    # L -> 7_pl -> instr. neutr. -> L
    feats[1577] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', 'instr. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. neutr.', node2.lemma)
    
    # L -> 159 -> -161 -> L
    feats[1578] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '-161') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-161', node2.lemma)
    
    # L -> 7_pl -> voc. fem -> L
    feats[1579] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', 'voc. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. fem', node2.lemma)
    
    # L -> 159 -> 5_tp -> L
    feats[1580] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '5_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_tp', node2.lemma)
    
    # L -> 159 -> 80 -> L
    feats[1581] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '80') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '80', node2.lemma)
    
    # L -> 155 -> -76 -> L
    feats[1582] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '155') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '155', '-76') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-76', node2.lemma)
    
    # L -> 149 -> -97 -> L
    feats[1583] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '-97') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-97', node2.lemma)
    
    # L -> 149 -> -59 -> L
    feats[1584] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '-59') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-59', node2.lemma)
    
    # L -> -56 -> 39 -> L
    feats[1585] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-56') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-56', '39') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '39', node2.lemma)
    
    # L -> 159 -> 9_du -> L
    feats[1586] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '9_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_du', node2.lemma)
    
    # L -> 149 -> -26 -> L
    feats[1587] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '-26') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-26', node2.lemma)
    
    # L -> 159 -> 12_tp -> L
    feats[1588] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '12_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_tp', node2.lemma)
    
    # L -> 159 -> 4_sp -> L
    feats[1589] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '4_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_sp', node2.lemma)
    
    # L -> 159 -> 13_fp -> L
    feats[1590] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '13_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_fp', node2.lemma)
    
    # L -> 159 -> -190 -> L
    feats[1591] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '-190') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-190', node2.lemma)
    
    # L -> 159 -> acc. pl. -> L
    feats[1592] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', 'acc. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. pl.', node2.lemma)
    
    # L -> 149 -> -36 -> L
    feats[1593] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '-36') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-36', node2.lemma)
    
    # L -> 149 -> instr -> L
    feats[1594] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', 'instr') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr', node2.lemma)
    
    # L -> 159 -> abl -> L
    feats[1595] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', 'abl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl', node2.lemma)
    
    # L -> 149 -> 9_sp -> L
    feats[1596] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '9_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_sp', node2.lemma)
    
    # L -> 159 -> 8_pl -> L
    feats[1597] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '8_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_pl', node2.lemma)
    
    # L -> 159 -> 3_pl -> L
    feats[1598] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '3_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_pl', node2.lemma)
    
    # L -> 159 -> 153 -> L
    feats[1599] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '153') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    
    # L -> 149 -> 169 -> L
    feats[1600] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '169') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '169', node2.lemma)
    
    # L -> 149 -> 2_du -> L
    feats[1601] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '2_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_du', node2.lemma)
    
    # L -> 149 -> sg -> L
    feats[1602] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', 'sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg', node2.lemma)
    
    # L -> 149 -> 9_sg -> L
    feats[1603] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '9_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_sg', node2.lemma)
    
    # L -> 159 -> 3_sg -> L
    feats[1604] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '3_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_sg', node2.lemma)
    
    # L -> 149 -> 7_pl -> L
    feats[1605] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '7_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_pl', node2.lemma)
    
    # L -> 159 -> 168 -> L
    feats[1606] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '168') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '168', node2.lemma)
    
    # L -> 159 -> -43 -> L
    feats[1607] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '-43') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-43', node2.lemma)
    
    # L -> 159 -> sg_fp -> L
    feats[1608] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', 'sg_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_fp', node2.lemma)
    
    # L -> 149 -> nom. du. -> L
    feats[1609] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', 'nom. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. du.', node2.lemma)
    
    # L -> 159 -> 4_pl -> L
    feats[1610] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '4_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_pl', node2.lemma)
    
    # L -> 159 -> 69 -> L
    feats[1611] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '69') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '69', node2.lemma)
    
    # L -> 159 -> pl_fp -> L
    feats[1612] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', 'pl_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_fp', node2.lemma)
    
    # L -> 149 -> acc. neutr. -> L
    feats[1613] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', 'acc. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. neutr.', node2.lemma)
    
    # L -> 159 -> 2 -> L
    feats[1614] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '2') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2', node2.lemma)
    
    # L -> 149 -> 5_du -> L
    feats[1615] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '5_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_du', node2.lemma)
    
    # L -> 149 -> 131 -> L
    feats[1616] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '131') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '131', node2.lemma)
    
    # L -> 159 -> tp -> L
    feats[1617] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', 'tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'tp', node2.lemma)
    
    # L -> 159 -> neutr -> L
    feats[1618] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', 'neutr') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'neutr', node2.lemma)
    
    # L -> 149 -> -82 -> L
    feats[1619] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '-82') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-82', node2.lemma)
    
    # L -> 149 -> 14_pl -> L
    feats[1620] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '14_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_pl', node2.lemma)
    
    # L -> 159 -> 111 -> L
    feats[1621] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '111') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '111', node2.lemma)
    
    # L -> 159 -> -54 -> L
    feats[1622] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '-54') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-54', node2.lemma)
    
    # L -> 149 -> 12_pl -> L
    feats[1623] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '12_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_pl', node2.lemma)
    
    # L -> 159 -> -117 -> L
    feats[1624] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '-117') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-117', node2.lemma)
    
    # L -> 159 -> -113 -> L
    feats[1625] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '-113') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-113', node2.lemma)
    
    # L -> 159 -> -299 -> L
    feats[1626] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '-299') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-299', node2.lemma)
    
    # L -> 149 -> -18 -> L
    feats[1627] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '-18') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-18', node2.lemma)
    
    # L -> 159 -> voc. pl. -> L
    feats[1628] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', 'voc. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. pl.', node2.lemma)
    
    # L -> 159 -> 30_sg -> L
    feats[1629] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '30_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_sg', node2.lemma)
    
    # L -> 149 -> 11_tp -> L
    feats[1630] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '11_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_tp', node2.lemma)
    
    # L -> 149 -> -87 -> L
    feats[1631] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '-87') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-87', node2.lemma)
    
    # L -> 149 -> voc. du. -> L
    feats[1632] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', 'voc. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. du.', node2.lemma)
    
    # L -> 159 -> -34 -> L
    feats[1633] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '-34') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-34', node2.lemma)
    
    # L -> 159 -> -83 -> L
    feats[1634] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '-83') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-83', node2.lemma)
    
    # L -> 159 -> nom. pl. -> L
    feats[1635] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', 'nom. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. pl.', node2.lemma)
    
    # L -> 155 -> 10_sg -> L
    feats[1636] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '155') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '155', '10_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_sg', node2.lemma)
    
    # L -> 155 -> -94 -> L
    feats[1637] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '155') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '155', '-94') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-94', node2.lemma)
    
    # L -> 149 -> -84 -> L
    feats[1638] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '-84') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-84', node2.lemma)
    
    # L -> 149 -> -301 -> L
    feats[1639] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '-301') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-301', node2.lemma)
    
    # L -> 159 -> 96 -> L
    feats[1640] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '96') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '96', node2.lemma)
    
    # L -> 149 -> -241 -> L
    feats[1641] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '-241') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-241', node2.lemma)
    
    # L -> 149 -> -57 -> L
    feats[1642] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '-57') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-57', node2.lemma)
    
    # L -> 159 -> 10_sg -> L
    feats[1643] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '10_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_sg', node2.lemma)
    
    # L -> 149 -> -72 -> L
    feats[1644] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '-72') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-72', node2.lemma)
    
    # L -> 159 -> -42 -> L
    feats[1645] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '-42') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-42', node2.lemma)
    
    # L -> 159 -> -94 -> L
    feats[1646] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '-94') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-94', node2.lemma)
    
    # L -> 149 -> 16_du -> L
    feats[1647] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '16_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_du', node2.lemma)
    
    # L -> -56 -> 141 -> L
    feats[1648] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-56') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-56', '141') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '141', node2.lemma)
    
    # L -> 159 -> 93 -> L
    feats[1649] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '93') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '93', node2.lemma)
    
    # L -> 159 -> 16_sg -> L
    feats[1650] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '16_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_sg', node2.lemma)
    
    # L -> 149 -> -68 -> L
    feats[1651] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '-68') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-68', node2.lemma)
    
    # L -> 159 -> pl_tp -> L
    feats[1652] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', 'pl_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_tp', node2.lemma)
    
    # L -> 159 -> -153 -> L
    feats[1653] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '-153') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-153', node2.lemma)
    
    # L -> 149 -> -91 -> L
    feats[1654] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '-91') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-91', node2.lemma)
    
    # L -> 149 -> -99 -> L
    feats[1655] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '-99') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-99', node2.lemma)
    
    # L -> 159 -> 14_sp -> L
    feats[1656] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '14_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_sp', node2.lemma)
    
    # L -> 159 -> 8_tp -> L
    feats[1657] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '8_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_tp', node2.lemma)
    
    # L -> 159 -> 139 -> L
    feats[1658] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '139') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '139', node2.lemma)
    
    # L -> 149 -> 10_pl -> L
    feats[1659] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '10_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_pl', node2.lemma)
    
    # L -> 159 -> -12 -> L
    feats[1660] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '-12') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-12', node2.lemma)
    
    # L -> 159 -> 7_du -> L
    feats[1661] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '7_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_du', node2.lemma)
    
    # L -> 149 -> -137 -> L
    feats[1662] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '-137') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-137', node2.lemma)
    
    # L -> 159 -> 122 -> L
    feats[1663] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '122') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '122', node2.lemma)
    
    # L -> 149 -> acc. adj. -> L
    feats[1664] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', 'acc. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. adj.', node2.lemma)
    
    # L -> 159 -> -62 -> L
    feats[1665] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '-62') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-62', node2.lemma)
    
    # L -> 149 -> -15 -> L
    feats[1666] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '-15') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-15', node2.lemma)
    
    # L -> 149 -> 11_pl -> L
    feats[1667] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '11_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_pl', node2.lemma)
    
    # L -> 159 -> -119 -> L
    feats[1668] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '-119') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-119', node2.lemma)
    
    # L -> 159 -> -14 -> L
    feats[1669] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '-14') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-14', node2.lemma)
    
    # L -> 149 -> -157 -> L
    feats[1670] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '-157') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-157', node2.lemma)
    
    # L -> 149 -> -302 -> L
    feats[1671] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '-302') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-302', node2.lemma)
    
    # L -> 159 -> abl. du. -> L
    feats[1672] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', 'abl. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. du.', node2.lemma)
    
    # L -> 159 -> -46 -> L
    feats[1673] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '-46') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-46', node2.lemma)
    
    # L -> 149 -> nom. neutr. -> L
    feats[1674] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', 'nom. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. neutr.', node2.lemma)
    
    # L -> 159 -> abl. pl. -> L
    feats[1675] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', 'abl. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. pl.', node2.lemma)
    
    # L -> 159 -> -20 -> L
    feats[1676] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '-20') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-20', node2.lemma)
    
    # L -> 159 -> 11_sg -> L
    feats[1677] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '11_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_sg', node2.lemma)
    
    # L -> 149 -> 37 -> L
    feats[1678] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '37') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '37', node2.lemma)
    
    # L -> 159 -> 162 -> L
    feats[1679] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '162') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '162', node2.lemma)
    
    # L -> 149 -> 8_sg -> L
    feats[1680] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '8_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_sg', node2.lemma)
    
    # L -> 149 -> 33 -> L
    feats[1681] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '33') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '33', node2.lemma)
    
    # L -> 149 -> 13_pl -> L
    feats[1682] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '13_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_pl', node2.lemma)
    
    # L -> 149 -> -96 -> L
    feats[1683] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '-96') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-96', node2.lemma)
    
    # L -> 149 -> -240 -> L
    feats[1684] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '-240') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-240', node2.lemma)
    
    # L -> -97 -> 6_fp -> L
    feats[1685] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-97') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-97', '6_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_fp', node2.lemma)
    
    # L -> 149 -> 179 -> L
    feats[1686] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '179') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '179', node2.lemma)
    
    # L -> 149 -> -283 -> L
    feats[1687] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '-283') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-283', node2.lemma)
    
    # L -> -97 -> -161 -> L
    feats[1688] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-97') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-97', '-161') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-161', node2.lemma)
    
    # L -> 155 -> -16 -> L
    feats[1689] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '155') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '155', '-16') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-16', node2.lemma)
    
    # L -> 149 -> 58 -> L
    feats[1690] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '58') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '58', node2.lemma)
    
    # L -> 149 -> 180 -> L
    feats[1691] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '180') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '180', node2.lemma)
    
    # L -> 149 -> 50 -> L
    feats[1692] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '50') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '50', node2.lemma)
    
    # L -> -97 -> 80 -> L
    feats[1693] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-97') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-97', '80') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '80', node2.lemma)
    
    # L -> 149 -> 5_tp -> L
    feats[1694] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '5_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_tp', node2.lemma)
    
    # L -> 149 -> 80 -> L
    feats[1695] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '80') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '80', node2.lemma)
    
    # L -> -97 -> 30_pl -> L
    feats[1696] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-97') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-97', '30_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_pl', node2.lemma)
    
    # L -> -137 -> -50 -> L
    feats[1697] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-137') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-137', '-50') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-50', node2.lemma)
    
    # L -> -97 -> 12_tp -> L
    feats[1698] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-97') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-97', '12_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_tp', node2.lemma)
    
    # L -> 149 -> 9_du -> L
    feats[1699] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '9_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_du', node2.lemma)
    
    # L -> 149 -> acc. pl. -> L
    feats[1700] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', 'acc. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. pl.', node2.lemma)
    
    # L -> -97 -> 3_pl -> L
    feats[1701] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-97') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-97', '3_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_pl', node2.lemma)
    
    # L -> 149 -> abl -> L
    feats[1702] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', 'abl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl', node2.lemma)
    
    # L -> -97 -> 134 -> L
    feats[1703] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-97') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-97', '134') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '134', node2.lemma)
    
    # L -> 149 -> 153 -> L
    feats[1704] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '153') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    
    # L -> -97 -> -47 -> L
    feats[1705] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-97') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-97', '-47') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-47', node2.lemma)
    
    # L -> 149 -> -166 -> L
    feats[1706] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '-166') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-166', node2.lemma)
    
    # L -> -97 -> -43 -> L
    feats[1707] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-97') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-97', '-43') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-43', node2.lemma)
    
    # L -> 149 -> -47 -> L
    feats[1708] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '-47') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-47', node2.lemma)
    
    # L -> -97 -> -309 -> L
    feats[1709] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-97') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-97', '-309') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-309', node2.lemma)
    
    # L -> -97 -> -129 -> L
    feats[1710] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-97') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-97', '-129') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-129', node2.lemma)
    
    # L -> 149 -> sg_fp -> L
    feats[1711] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', 'sg_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_fp', node2.lemma)
    
    # L -> -97 -> pl_fp -> L
    feats[1712] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-97') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-97', 'pl_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_fp', node2.lemma)
    
    # L -> 149 -> 13_sg -> L
    feats[1713] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '13_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_sg', node2.lemma)
    
    # L -> -97 -> 135 -> L
    feats[1714] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-97') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-97', '135') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '135', node2.lemma)
    
    # L -> 149 -> -271 -> L
    feats[1715] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '-271') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-271', node2.lemma)
    
    # L -> 149 -> 40 -> L
    feats[1716] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '40') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '40', node2.lemma)
    
    # L -> -97 -> 111 -> L
    feats[1717] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-97') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-97', '111') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '111', node2.lemma)
    
    # L -> 149 -> neutr -> L
    feats[1718] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', 'neutr') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'neutr', node2.lemma)
    
    # L -> 149 -> 36 -> L
    feats[1719] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '36') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '36', node2.lemma)
    
    # L -> 149 -> 16_fp -> L
    feats[1720] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '16_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_fp', node2.lemma)
    
    # L -> -97 -> -22 -> L
    feats[1721] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-97') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-97', '-22') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-22', node2.lemma)
    
    # L -> 149 -> 142 -> L
    feats[1722] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '142') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '142', node2.lemma)
    
    # L -> -97 -> -117 -> L
    feats[1723] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-97') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-97', '-117') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-117', node2.lemma)
    
    # L -> -97 -> -299 -> L
    feats[1724] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-97') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-97', '-299') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-299', node2.lemma)
    
    # L -> instr -> instr. adj. -> L
    feats[1725] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', 'instr. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. adj.', node2.lemma)
    
    # L -> instr -> 4_tp -> L
    feats[1726] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', '4_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_tp', node2.lemma)
    
    # L -> -144 -> -263 -> L
    feats[1727] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-144') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-144', '-263') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-263', node2.lemma)
    
    # L -> -97 -> masc -> L
    feats[1728] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-97') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-97', 'masc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'masc', node2.lemma)
    
    # L -> 149 -> 32 -> L
    feats[1729] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '32') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '32', node2.lemma)
    
    # L -> -151 -> 16_pl -> L
    feats[1730] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-151') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-151', '16_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_pl', node2.lemma)
    
    # L -> -151 -> 14_pl -> L
    feats[1731] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-151') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-151', '14_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_pl', node2.lemma)
    
    # L -> -97 -> 10_sg -> L
    feats[1732] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-97') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-97', '10_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_sg', node2.lemma)
    
    # L -> -97 -> 3_fp -> L
    feats[1733] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-97') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-97', '3_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_fp', node2.lemma)
    
    # L -> 149 -> instr. sg. -> L
    feats[1734] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', 'instr. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. sg.', node2.lemma)
    
    # L -> instr -> 13_sp -> L
    feats[1735] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', '13_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_sp', node2.lemma)
    
    # L -> 149 -> -266 -> L
    feats[1736] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '-266') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-266', node2.lemma)
    
    # L -> 149 -> 96 -> L
    feats[1737] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '96') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '96', node2.lemma)
    
    # L -> instr -> 170 -> L
    feats[1738] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', '170') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '170', node2.lemma)
    
    # L -> 149 -> 10_sg -> L
    feats[1739] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '10_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_sg', node2.lemma)
    
    # L -> -97 -> 2_sg -> L
    feats[1740] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-97') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-97', '2_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_sg', node2.lemma)
    
    # L -> 149 -> voc. masc. -> L
    feats[1741] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', 'voc. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. masc.', node2.lemma)
    
    # L -> instr -> -45 -> L
    feats[1742] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', '-45') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-45', node2.lemma)
    
    # L -> instr -> 34 -> L
    feats[1743] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', '34') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '34', node2.lemma)
    
    # L -> -97 -> 16_sg -> L
    feats[1744] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-97') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-97', '16_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_sg', node2.lemma)
    
    # L -> -97 -> -27 -> L
    feats[1745] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-97') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-97', '-27') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-27', node2.lemma)
    
    # L -> -97 -> 4_sg -> L
    feats[1746] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-97') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-97', '4_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_sg', node2.lemma)
    
    # L -> instr -> -63 -> L
    feats[1747] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', '-63') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-63', node2.lemma)
    
    # L -> -97 -> -153 -> L
    feats[1748] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-97') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-97', '-153') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-153', node2.lemma)
    
    # L -> 149 -> -35 -> L
    feats[1749] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '-35') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    
    # L -> 149 -> 16_sg -> L
    feats[1750] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '16_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_sg', node2.lemma)
    
    # L -> instr -> -68 -> L
    feats[1751] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', '-68') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-68', node2.lemma)
    
    # L -> 149 -> 93 -> L
    feats[1752] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', '93') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '93', node2.lemma)
    
    # L -> 173 -> -68 -> L
    feats[1753] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '173') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '173', '-68') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-68', node2.lemma)
    
    # L -> instr -> du_sp -> L
    feats[1754] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', 'du_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_sp', node2.lemma)
    
    # L -> -97 -> -12 -> L
    feats[1755] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-97') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-97', '-12') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-12', node2.lemma)
    
    # L -> instr -> -15 -> L
    feats[1756] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', '-15') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-15', node2.lemma)
    
    # L -> instr -> -29 -> L
    feats[1757] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', '-29') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-29', node2.lemma)
    
    # L -> -97 -> 15_sp -> L
    feats[1758] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-97') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-97', '15_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_sp', node2.lemma)
    
    # L -> -97 -> abl. du. -> L
    feats[1759] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-97') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-97', 'abl. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. du.', node2.lemma)
    
    # L -> -97 -> -46 -> L
    feats[1760] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-97') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-97', '-46') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-46', node2.lemma)
    
    # L -> instr -> -268 -> L
    feats[1761] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', '-268') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-268', node2.lemma)
    
    # L -> instr -> -157 -> L
    feats[1762] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', '-157') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-157', node2.lemma)
    
    # L -> instr -> 16_tp -> L
    feats[1763] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', '16_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_tp', node2.lemma)
    
    # L -> -97 -> -20 -> L
    feats[1764] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-97') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-97', '-20') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-20', node2.lemma)
    
    # L -> -97 -> 92 -> L
    feats[1765] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-97') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-97', '92') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '92', node2.lemma)
    
    # L -> instr -> 10_tp -> L
    feats[1766] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', '10_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_tp', node2.lemma)
    
    # L -> instr -> -152 -> L
    feats[1767] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', '-152') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-152', node2.lemma)
    
    # L -> instr -> 38 -> L
    feats[1768] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', '38') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '38', node2.lemma)
    
    # L -> instr -> 37 -> L
    feats[1769] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', '37') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '37', node2.lemma)
    
    # L -> instr -> 95 -> L
    feats[1770] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', '95') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '95', node2.lemma)
    
    # L -> -97 -> -61 -> L
    feats[1771] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-97') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-97', '-61') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-61', node2.lemma)
    
    # L -> -97 -> -156 -> L
    feats[1772] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-97') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-97', '-156') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-156', node2.lemma)
    
    # L -> instr -> -48 -> L
    feats[1773] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', '-48') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-48', node2.lemma)
    
    # L -> instr -> 8_sg -> L
    feats[1774] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', '8_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_sg', node2.lemma)
    
    # L -> instr -> -24 -> L
    feats[1775] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', '-24') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    
    # L -> instr -> -240 -> L
    feats[1776] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', '-240') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-240', node2.lemma)
    
    # L -> instr -> loc. sg. -> L
    feats[1777] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', 'loc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc. sg.', node2.lemma)
    
    # L -> instr -> 81 -> L
    feats[1778] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', '81') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '81', node2.lemma)
    
    # L -> -151 -> 33 -> L
    feats[1779] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-151') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-151', '33') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '33', node2.lemma)
    
    # L -> instr -> pl -> L
    feats[1780] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', 'pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    
    # L -> instr -> -283 -> L
    feats[1781] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', '-283') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-283', node2.lemma)
    
    # L -> instr -> 68 -> L
    feats[1782] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', '68') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '68', node2.lemma)
    
    # L -> instr -> 9_fp -> L
    feats[1783] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', '9_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_fp', node2.lemma)
    
    # L -> instr -> 39 -> L
    feats[1784] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', '39') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '39', node2.lemma)
    
    # L -> instr -> 58 -> L
    feats[1785] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', '58') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '58', node2.lemma)
    
    # L -> instr -> 50 -> L
    feats[1786] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', '50') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '50', node2.lemma)
    
    # L -> instr -> 118 -> L
    feats[1787] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', '118') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '118', node2.lemma)
    
    # L -> instr -> 5_tp -> L
    feats[1788] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', '5_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_tp', node2.lemma)
    
    # L -> 37 -> -121 -> L
    feats[1789] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '37', '-121') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-121', node2.lemma)
    
    # L -> instr -> 10_fp -> L
    feats[1790] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', '10_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_fp', node2.lemma)
    
    # L -> instr -> 61 -> L
    feats[1791] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', '61') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '61', node2.lemma)
    
    # L -> 149 -> acc -> L
    feats[1792] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', 'acc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc', node2.lemma)
    
    # L -> instr -> -109 -> L
    feats[1793] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', '-109') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-109', node2.lemma)
    
    # L -> instr -> 89 -> L
    feats[1794] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', '89') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '89', node2.lemma)
    
    # L -> 37 -> 3_sp -> L
    feats[1795] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '37', '3_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_sp', node2.lemma)
    
    # L -> instr -> 29 -> L
    feats[1796] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', '29') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '29', node2.lemma)
    
    # L -> 37 -> -296 -> L
    feats[1797] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '37', '-296') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-296', node2.lemma)
    
    # L -> instr -> 3_du -> L
    feats[1798] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', '3_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    
    # L -> 120 -> -26 -> L
    feats[1799] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '120') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '120', '-26') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-26', node2.lemma)
    
    # L -> 37 -> instr. neutr. -> L
    feats[1800] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '37', 'instr. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. neutr.', node2.lemma)
    
    # L -> instr -> 4_sp -> L
    feats[1801] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', '4_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_sp', node2.lemma)
    
    # L -> -144 -> 182 -> L
    feats[1802] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-144') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-144', '182') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '182', node2.lemma)
    
    # L -> instr -> 114 -> L
    feats[1803] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', '114') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '114', node2.lemma)
    
    # L -> instr -> -166 -> L
    feats[1804] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', '-166') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-166', node2.lemma)
    
    # L -> instr -> -154 -> L
    feats[1805] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', '-154') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-154', node2.lemma)
    
    # L -> instr -> 134 -> L
    feats[1806] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', '134') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '134', node2.lemma)
    
    # L -> instr -> 150 -> L
    feats[1807] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', '150') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '150', node2.lemma)
    
    # L -> 95 -> 7_sp -> L
    feats[1808] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '95') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '95', '7_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_sp', node2.lemma)
    
    # L -> instr -> 3_sg -> L
    feats[1809] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', '3_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_sg', node2.lemma)
    
    # L -> instr -> 30_tp -> L
    feats[1810] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', '30_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_tp', node2.lemma)
    
    # L -> 95 -> voc. neutr. -> L
    feats[1811] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '95') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '95', 'voc. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. neutr.', node2.lemma)
    
    # L -> instr -> 168 -> L
    feats[1812] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', '168') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '168', node2.lemma)
    
    # L -> 95 -> -10 -> L
    feats[1813] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '95') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '95', '-10') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-10', node2.lemma)
    
    # L -> instr -> -309 -> L
    feats[1814] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', '-309') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-309', node2.lemma)
    
    # L -> instr -> -129 -> L
    feats[1815] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', '-129') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-129', node2.lemma)
    
    # L -> 95 -> -36 -> L
    feats[1816] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '95') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '95', '-36') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-36', node2.lemma)
    
    # L -> instr -> -163 -> L
    feats[1817] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', '-163') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-163', node2.lemma)
    
    # L -> instr -> 4_pl -> L
    feats[1818] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', '4_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_pl', node2.lemma)
    
    # L -> 95 -> 156 -> L
    feats[1819] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '95') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '95', '156') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '156', node2.lemma)
    
    # L -> 95 -> 75 -> L
    feats[1820] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '95') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '95', '75') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '75', node2.lemma)
    
    # L -> instr -> 4_du -> L
    feats[1821] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', '4_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_du', node2.lemma)
    
    # L -> instr -> 2 -> L
    feats[1822] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', '2') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2', node2.lemma)
    
    # L -> -151 -> 98 -> L
    feats[1823] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-151') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-151', '98') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '98', node2.lemma)
    
    # L -> 95 -> 30 -> L
    feats[1824] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '95') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '95', '30') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30', node2.lemma)
    
    # L -> 95 -> 12_sg -> L
    feats[1825] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '95') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '95', '12_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_sg', node2.lemma)
    
    # L -> 95 -> -132 -> L
    feats[1826] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '95') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '95', '-132') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-132', node2.lemma)
    
    # L -> 95 -> dat. pl. -> L
    feats[1827] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '95') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '95', 'dat. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. pl.', node2.lemma)
    
    # L -> 14_fp -> -90 -> L
    feats[1828] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_fp', '-90') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-90', node2.lemma)
    
    # L -> 14_fp -> 72 -> L
    feats[1829] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_fp', '72') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '72', node2.lemma)
    
    # L -> 95 -> -52 -> L
    feats[1830] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '95') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '95', '-52') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-52', node2.lemma)
    
    # L -> 95 -> 173 -> L
    feats[1831] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '95') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '95', '173') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '173', node2.lemma)
    
    # L -> 95 -> -69 -> L
    feats[1832] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '95') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '95', '-69') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-69', node2.lemma)
    
    # L -> 95 -> 88 -> L
    feats[1833] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '95') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '95', '88') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '88', node2.lemma)
    
    # L -> 95 -> 12_pl -> L
    feats[1834] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '95') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '95', '12_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_pl', node2.lemma)
    
    # L -> 95 -> instr. adj. -> L
    feats[1835] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '95') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '95', 'instr. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. adj.', node2.lemma)
    
    # L -> 95 -> 132 -> L
    feats[1836] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '95') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '95', '132') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '132', node2.lemma)
    
    # L -> 95 -> -78 -> L
    feats[1837] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '95') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '95', '-78') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-78', node2.lemma)
    
    # L -> 95 -> -87 -> L
    feats[1838] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '95') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '95', '-87') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-87', node2.lemma)
    
    # L -> 95 -> -57 -> L
    feats[1839] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '95') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '95', '-57') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-57', node2.lemma)
    
    # L -> 95 -> -72 -> L
    feats[1840] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '95') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '95', '-72') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-72', node2.lemma)
    
    # L -> 95 -> 2_sp -> L
    feats[1841] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '95') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '95', '2_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_sp', node2.lemma)
    
    # L -> 95 -> -45 -> L
    feats[1842] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '95') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '95', '-45') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-45', node2.lemma)
    
    # L -> 95 -> 60 -> L
    feats[1843] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '95') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '95', '60') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '60', node2.lemma)
    
    # L -> 95 -> -68 -> L
    feats[1844] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '95') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '95', '-68') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-68', node2.lemma)
    
    # L -> 95 -> -243 -> L
    feats[1845] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '95') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '95', '-243') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-243', node2.lemma)
    
    # L -> 95 -> 109 -> L
    feats[1846] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '95') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '95', '109') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '109', node2.lemma)
    
    # L -> instr -> 129 -> L
    feats[1847] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', '129') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '129', node2.lemma)
    
    # L -> -151 -> -46 -> L
    feats[1848] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-151') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-151', '-46') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-46', node2.lemma)
    
    # L -> -151 -> -297 -> L
    feats[1849] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-151') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-151', '-297') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-297', node2.lemma)
    
    # L -> 120 -> 9_tp -> L
    feats[1850] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '120') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '120', '9_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_tp', node2.lemma)
    
    # L -> 156 -> 75 -> L
    feats[1851] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '156', '75') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '75', node2.lemma)
    
    # L -> 156 -> nom. adj. -> L
    feats[1852] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '156', 'nom. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. adj.', node2.lemma)
    
    # L -> -29 -> -59 -> L
    feats[1853] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-29') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-29', '-59') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-59', node2.lemma)
    
    # L -> -29 -> -30 -> L
    feats[1854] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-29') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-29', '-30') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-30', node2.lemma)
    
    # L -> 33 -> -59 -> L
    feats[1855] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '33', '-59') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-59', node2.lemma)
    
    # L -> 156 -> 31 -> L
    feats[1856] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '156', '31') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '31', node2.lemma)
    
    # L -> 156 -> 5_du -> L
    feats[1857] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '156', '5_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_du', node2.lemma)
    
    # L -> 14_fp -> -22 -> L
    feats[1858] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_fp', '-22') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-22', node2.lemma)
    
    # L -> 156 -> 27_pl -> L
    feats[1859] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '156', '27_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_pl', node2.lemma)
    
    # L -> 156 -> -101 -> L
    feats[1860] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '156', '-101') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-101', node2.lemma)
    
    # L -> 156 -> -241 -> L
    feats[1861] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '156', '-241') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-241', node2.lemma)
    
    # L -> 156 -> nom. sg. -> L
    feats[1862] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '156', 'nom. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. sg.', node2.lemma)
    
    # L -> 120 -> 28_sg -> L
    feats[1863] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '120') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '120', '28_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28_sg', node2.lemma)
    
    # L -> 120 -> -35 -> L
    feats[1864] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '120') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '120', '-35') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    
    # L -> 156 -> -292 -> L
    feats[1865] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '156', '-292') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-292', node2.lemma)
    
    # L -> 95 -> -34 -> L
    feats[1866] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '95') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '95', '-34') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-34', node2.lemma)
    
    # L -> 156 -> 15_tp -> L
    feats[1867] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '156', '15_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_tp', node2.lemma)
    
    # L -> 33 -> -301 -> L
    feats[1868] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '33', '-301') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-301', node2.lemma)
    
    # L -> 156 -> 117 -> L
    feats[1869] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '156', '117') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '117', node2.lemma)
    
    # L -> 156 -> 11_pl -> L
    feats[1870] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '156', '11_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_pl', node2.lemma)
    
    # L -> 156 -> -268 -> L
    feats[1871] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '156', '-268') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-268', node2.lemma)
    
    # L -> 156 -> 14_sg -> L
    feats[1872] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '156', '14_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_sg', node2.lemma)
    
    # L -> 156 -> 12_fp -> L
    feats[1873] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '156', '12_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_fp', node2.lemma)
    
    # L -> 156 -> -291 -> L
    feats[1874] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '156', '-291') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-291', node2.lemma)
    
    # L -> 156 -> 70 -> L
    feats[1875] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '156', '70') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '70', node2.lemma)
    
    # L -> 156 -> gen. pl. -> L
    feats[1876] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '156', 'gen. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen. pl.', node2.lemma)
    
    # L -> 14_fp -> 129 -> L
    feats[1877] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_fp', '129') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '129', node2.lemma)
    
    # L -> 156 -> 13_pl -> L
    feats[1878] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '156', '13_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_pl', node2.lemma)
    
    # L -> 156 -> -147 -> L
    feats[1879] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '156', '-147') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-147', node2.lemma)
    
    # L -> 156 -> loc. sg. -> L
    feats[1880] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '156', 'loc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc. sg.', node2.lemma)
    
    # L -> 156 -> 15_sg -> L
    feats[1881] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '156', '15_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_sg', node2.lemma)
    
    # L -> 156 -> -19 -> L
    feats[1882] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '156', '-19') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-19', node2.lemma)
    
    # L -> 156 -> 128 -> L
    feats[1883] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '156', '128') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '128', node2.lemma)
    
    # L -> -29 -> -240 -> L
    feats[1884] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-29') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-29', '-240') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-240', node2.lemma)
    
    # L -> 95 -> 176 -> L
    feats[1885] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '95') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '95', '176') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '176', node2.lemma)
    
    # L -> -29 -> 97 -> L
    feats[1886] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-29') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-29', '97') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '97', node2.lemma)
    
    # L -> -29 -> -166 -> L
    feats[1887] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-29') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-29', '-166') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-166', node2.lemma)
    
    # L -> 7_sp -> 14_pl -> L
    feats[1888] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_sp', '14_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_pl', node2.lemma)
    
    # L -> -23 -> 131 -> L
    feats[1889] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-23') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-23', '131') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '131', node2.lemma)
    
    # L -> -23 -> 159 -> L
    feats[1890] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-23') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-23', '159') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '159', node2.lemma)
    
    # L -> 33 -> -22 -> L
    feats[1891] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '33', '-22') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-22', node2.lemma)
    
    # L -> 12_sg -> 7_tp -> L
    feats[1892] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', '7_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_tp', node2.lemma)
    
    # L -> 12_sg -> -143 -> L
    feats[1893] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', '-143') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-143', node2.lemma)
    
    # L -> -29 -> 2_sg -> L
    feats[1894] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-29') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-29', '2_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_sg', node2.lemma)
    
    # L -> -29 -> 93 -> L
    feats[1895] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-29') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-29', '93') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '93', node2.lemma)
    
    # L -> 7_sp -> -144 -> L
    feats[1896] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_sp', '-144') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-144', node2.lemma)
    
    # L -> 156 -> 152 -> L
    feats[1897] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '156', '152') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '152', node2.lemma)
    
    # L -> 33 -> 129 -> L
    feats[1898] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '33', '129') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '129', node2.lemma)
    
    # L -> -23 -> 108 -> L
    feats[1899] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-23') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-23', '108') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '108', node2.lemma)
    
    # L -> voc. sg. -> 100 -> L
    feats[1900] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. sg.', '100') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '100', node2.lemma)
    
    # L -> 12_sg -> 61 -> L
    feats[1901] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', '61') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '61', node2.lemma)
    
    # L -> -29 -> 6_tp -> L
    feats[1902] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-29') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-29', '6_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_tp', node2.lemma)
    
    # L -> 7_sp -> sg_fp -> L
    feats[1903] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_sp', 'sg_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_fp', node2.lemma)
    
    # L -> 7_sp -> 168 -> L
    feats[1904] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_sp', '168') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '168', node2.lemma)
    
    # L -> -23 -> 98 -> L
    feats[1905] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-23') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-23', '98') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '98', node2.lemma)
    
    # L -> -23 -> 3_sg -> L
    feats[1906] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-23') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-23', '3_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_sg', node2.lemma)
    
    # L -> 9_sp -> 173 -> L
    feats[1907] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '173') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '173', node2.lemma)
    
    # L -> 13_pl -> 14_fp -> L
    feats[1908] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '14_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_fp', node2.lemma)
    
    # L -> 12_sg -> 76 -> L
    feats[1909] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', '76') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '76', node2.lemma)
    
    # L -> 12_sg -> 32 -> L
    feats[1910] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', '32') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '32', node2.lemma)
    
    # L -> -262 -> 90 -> L
    feats[1911] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-262') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-262', '90') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '90', node2.lemma)
    
    # L -> 7_sp -> -20 -> L
    feats[1912] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_sp', '-20') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-20', node2.lemma)
    
    # L -> 9_sp -> 14_sg -> L
    feats[1913] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '14_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_sg', node2.lemma)
    
    # L -> 13_pl -> abl. sg. -> L
    feats[1914] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', 'abl. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. sg.', node2.lemma)
    
    # L -> 9_sp -> gen -> L
    feats[1915] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', 'gen') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen', node2.lemma)
    
    # L -> 9_sp -> -33 -> L
    feats[1916] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '-33') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-33', node2.lemma)
    
    # L -> 9_sp -> -144 -> L
    feats[1917] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '-144') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-144', node2.lemma)
    
    # L -> 9_sp -> 33 -> L
    feats[1918] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '33') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '33', node2.lemma)
    
    # L -> 9_sp -> 13_pl -> L
    feats[1919] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '13_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_pl', node2.lemma)
    
    # L -> 9_sp -> 77 -> L
    feats[1920] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '77') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '77', node2.lemma)
    
    # L -> 9_sp -> -24 -> L
    feats[1921] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '-24') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    
    # L -> 9_sp -> -96 -> L
    feats[1922] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '-96') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-96', node2.lemma)
    
    # L -> 9_sp -> acc. masc. -> L
    feats[1923] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', 'acc. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. masc.', node2.lemma)
    
    # L -> 9_sp -> -133 -> L
    feats[1924] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '-133') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-133', node2.lemma)
    
    # L -> 9_sp -> 140 -> L
    feats[1925] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '140') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '140', node2.lemma)
    
    # L -> 9_sp -> 6_fp -> L
    feats[1926] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '6_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_fp', node2.lemma)
    
    # L -> 9_sp -> 179 -> L
    feats[1927] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '179') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '179', node2.lemma)
    
    # L -> 9_sp -> 15_sg -> L
    feats[1928] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '15_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_sg', node2.lemma)
    
    # L -> 9_sp -> 15_du -> L
    feats[1929] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '15_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_du', node2.lemma)
    
    # L -> 9_sp -> 121 -> L
    feats[1930] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '121') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '121', node2.lemma)
    
    # L -> 9_sp -> -32 -> L
    feats[1931] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '-32') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-32', node2.lemma)
    
    # L -> 9_sp -> 1 -> L
    feats[1932] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '1') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '1', node2.lemma)
    
    # L -> 9_sp -> 28_tp -> L
    feats[1933] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '28_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28_tp', node2.lemma)
    
    # L -> 9_sp -> -112 -> L
    feats[1934] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '-112') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    
    # L -> 9_sp -> 128 -> L
    feats[1935] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '128') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '128', node2.lemma)
    
    # L -> 9_sp -> -283 -> L
    feats[1936] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '-283') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-283', node2.lemma)
    
    # L -> 9_sp -> -114 -> L
    feats[1937] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '-114') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-114', node2.lemma)
    
    # L -> 9_sp -> 58 -> L
    feats[1938] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '58') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '58', node2.lemma)
    
    # L -> -262 -> 140 -> L
    feats[1939] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-262') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-262', '140') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '140', node2.lemma)
    
    # L -> 9_sp -> 180 -> L
    feats[1940] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '180') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '180', node2.lemma)
    
    # L -> 9_sp -> 5_tp -> L
    feats[1941] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '5_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_tp', node2.lemma)
    
    # L -> 9_sp -> -131 -> L
    feats[1942] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '-131') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-131', node2.lemma)
    
    # L -> 9_sp -> 97 -> L
    feats[1943] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '97') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '97', node2.lemma)
    
    # L -> 9_sp -> 10_fp -> L
    feats[1944] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '10_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_fp', node2.lemma)
    
    # L -> 9_sp -> 80 -> L
    feats[1945] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '80') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '80', node2.lemma)
    
    # L -> 9_sp -> -109 -> L
    feats[1946] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '-109') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-109', node2.lemma)
    
    # L -> 9_sp -> 89 -> L
    feats[1947] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '89') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '89', node2.lemma)
    
    # L -> 9_sp -> instr. fem -> L
    feats[1948] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', 'instr. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. fem', node2.lemma)
    
    # L -> 9_sp -> 30_pl -> L
    feats[1949] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '30_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_pl', node2.lemma)
    
    # L -> 9_sp -> 94 -> L
    feats[1950] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '94') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '94', node2.lemma)
    
    # L -> 9_sp -> -76 -> L
    feats[1951] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '-76') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-76', node2.lemma)
    
    # L -> 9_sp -> nom -> L
    feats[1952] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', 'nom') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom', node2.lemma)
    
    # L -> 9_sp -> -247 -> L
    feats[1953] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '-247') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-247', node2.lemma)
    
    # L -> 9_sp -> 9_du -> L
    feats[1954] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '9_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_du', node2.lemma)
    
    # L -> 9_sp -> 171 -> L
    feats[1955] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '171') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '171', node2.lemma)
    
    # L -> 9_sp -> 3_du -> L
    feats[1956] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '3_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    
    # L -> 9_sp -> pl_sp -> L
    feats[1957] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', 'pl_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_sp', node2.lemma)
    
    # L -> -59 -> -26 -> L
    feats[1958] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-59') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-59', '-26') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-26', node2.lemma)
    
    # L -> 9_sp -> 12_tp -> L
    feats[1959] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '12_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_tp', node2.lemma)
    
    # L -> -59 -> -10 -> L
    feats[1960] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-59') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-59', '-10') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-10', node2.lemma)
    
    # L -> 9_sp -> 13_fp -> L
    feats[1961] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '13_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_fp', node2.lemma)
    
    # L -> 9_sp -> -190 -> L
    feats[1962] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '-190') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-190', node2.lemma)
    
    # L -> 9_sp -> 6_pl -> L
    feats[1963] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '6_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_pl', node2.lemma)
    
    # L -> 9_sp -> 9_tp -> L
    feats[1964] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '9_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_tp', node2.lemma)
    
    # L -> 9_sp -> abl -> L
    feats[1965] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', 'abl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl', node2.lemma)
    
    # L -> 9_sp -> -41 -> L
    feats[1966] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '-41') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-41', node2.lemma)
    
    # L -> 9_sp -> 8_pl -> L
    feats[1967] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '8_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_pl', node2.lemma)
    
    # L -> 9_sp -> 153 -> L
    feats[1968] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '153') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    
    # L -> 9_sp -> -154 -> L
    feats[1969] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '-154') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-154', node2.lemma)
    
    # L -> 9_sp -> 150 -> L
    feats[1970] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '150') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '150', node2.lemma)
    
    # L -> 9_sp -> 98 -> L
    feats[1971] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '98') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '98', node2.lemma)
    
    # L -> 9_sp -> 30_tp -> L
    feats[1972] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '30_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_tp', node2.lemma)
    
    # L -> 9_sp -> sg_fp -> L
    feats[1973] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', 'sg_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_fp', node2.lemma)
    
    # L -> 9_sp -> 168 -> L
    feats[1974] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '168') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '168', node2.lemma)
    
    # L -> 9_sp -> -43 -> L
    feats[1975] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '-43') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-43', node2.lemma)
    
    # L -> 9_sp -> 14_tp -> L
    feats[1976] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '14_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_tp', node2.lemma)
    
    # L -> 9_sp -> -102 -> L
    feats[1977] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '-102') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-102', node2.lemma)
    
    # L -> 9_sp -> -159 -> L
    feats[1978] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '-159') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-159', node2.lemma)
    
    # L -> 9_sp -> -129 -> L
    feats[1979] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '-129') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-129', node2.lemma)
    
    # L -> 9_sp -> 13_sg -> L
    feats[1980] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '13_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_sg', node2.lemma)
    
    # L -> 9_sp -> 5_pl -> L
    feats[1981] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '5_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_pl', node2.lemma)
    
    # L -> 9_sp -> gen. du. -> L
    feats[1982] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', 'gen. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen. du.', node2.lemma)
    
    # L -> 9_sp -> -163 -> L
    feats[1983] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '-163') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-163', node2.lemma)
    
    # L -> 9_sp -> 4_pl -> L
    feats[1984] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '4_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_pl', node2.lemma)
    
    # L -> 9_sp -> 69 -> L
    feats[1985] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '69') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '69', node2.lemma)
    
    # L -> 9_sp -> pl_fp -> L
    feats[1986] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', 'pl_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_fp', node2.lemma)
    
    # L -> 9_sp -> 4_du -> L
    feats[1987] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '4_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_du', node2.lemma)
    
    # L -> -132 -> 5_du -> L
    feats[1988] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-132') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-132', '5_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_du', node2.lemma)
    
    # L -> -132 -> sp -> L
    feats[1989] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-132') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-132', 'sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sp', node2.lemma)
    
    # L -> -262 -> -309 -> L
    feats[1990] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-262') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-262', '-309') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-309', node2.lemma)
    
    # L -> -262 -> -129 -> L
    feats[1991] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-262') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-262', '-129') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-129', node2.lemma)
    
    # L -> 13_pl -> -83 -> L
    feats[1992] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '-83') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-83', node2.lemma)
    
    # L -> 13_pl -> pl_tp -> L
    feats[1993] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', 'pl_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_tp', node2.lemma)
    
    # L -> 13_pl -> 82 -> L
    feats[1994] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '82') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '82', node2.lemma)
    
    # L -> 9_sp -> 162 -> L
    feats[1995] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', '162') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '162', node2.lemma)
    
    # L -> 13_pl -> 139 -> L
    feats[1996] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '139') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '139', node2.lemma)
    
    # L -> -132 -> gen -> L
    feats[1997] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-132') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-132', 'gen') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen', node2.lemma)
    
    # L -> -262 -> 162 -> L
    feats[1998] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-262') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-262', '162') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '162', node2.lemma)
    
    # L -> 13_pl -> -79 -> L
    feats[1999] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '-79') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-79', node2.lemma)
    
    feats[feats < 1e-25] = 1e-25
    return -np.log10(feats)
    