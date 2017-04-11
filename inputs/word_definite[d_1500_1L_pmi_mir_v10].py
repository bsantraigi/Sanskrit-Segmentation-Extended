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
    
    # L->instr. sg.->L
    feats[0] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. sg.', node2.lemma)
            
    # L->dat. du.->L
    feats[1] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
            
    # L->-266->L
    feats[2] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-266') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-266', node2.lemma)
            
    # L->40->L
    feats[3] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '40') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '40', node2.lemma)
            
    # L->12_tp->L
    feats[4] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_tp', node2.lemma)
            
    # L->31->C
    feats[5] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '31') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '31', node2.cng)
            
    # L->94->C
    feats[6] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '94') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '94', node2.cng)
            
    # L->-14->C
    feats[7] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-14') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-14', node2.cng)
            
    # L->voc. fem->C
    feats[8] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. fem', node2.cng)
            
    # L->8_du->C
    feats[9] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_du', node2.cng)
            
    # L->5_sg->C
    feats[10] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sg', node2.cng)
            
    # L->15_sg->C
    feats[11] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_sg', node2.cng)
            
    # L->-241->C
    feats[12] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-241') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-241', node2.cng)
            
    # L->6_sg->C
    feats[13] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', node2.cng)
            
    # L->-308->C
    feats[14] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', node2.cng)
            
    # L->-77->C
    feats[15] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-77') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-77', node2.cng)
            
    # L->30_sp->C
    feats[16] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_sp', node2.cng)
            
    # L->151->C
    feats[17] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '151') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '151', node2.cng)
            
    # L->8_sp->C
    feats[18] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_sp', node2.cng)
            
    # L->132->C
    feats[19] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '132') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '132', node2.cng)
            
    # L->nom. du.->C
    feats[20] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. du.', node2.cng)
            
    # L->38->C
    feats[21] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '38') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '38', node2.cng)
            
    # L->42->C
    feats[22] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '42') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '42', node2.cng)
            
    # L->9_fp->C
    feats[23] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_fp', node2.cng)
            
    # L->3_pl->C
    feats[24] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_pl', node2.cng)
            
    # L->10_fp->C
    feats[25] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_fp', node2.cng)
            
    # L->99->C
    feats[26] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '99') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '99', node2.cng)
            
    # L->142->C
    feats[27] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '142') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '142', node2.cng)
            
    # L->gen. du.->C
    feats[28] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. du.', node2.cng)
            
    # L->-246->C
    feats[29] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-246') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-246', node2.cng)
            
    # L->9_tp->C
    feats[30] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_tp', node2.cng)
            
    # L->voc. masc.->C
    feats[31] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. masc.', node2.cng)
            
    # L->32->C
    feats[32] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '32') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '32', node2.cng)
            
    # L->-220->C
    feats[33] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-220') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-220', node2.cng)
            
    # L -> 28_sg -> 134 -> L
    feats[34] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28_sg', '134') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '134', node2.lemma)
    
    # L -> 28_sg -> 34 -> L
    feats[35] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28_sg', '34') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '34', node2.lemma)
    
    # L -> 28_sg -> instr. neutr. -> L
    feats[36] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28_sg', 'instr. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. neutr.', node2.lemma)
    
    # L -> 28_sg -> 100 -> L
    feats[37] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28_sg', '100') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '100', node2.lemma)
    
    # L -> 28_sg -> -41 -> L
    feats[38] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28_sg', '-41') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-41', node2.lemma)
    
    # L -> 28_sg -> 12_tp -> L
    feats[39] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28_sg', '12_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_tp', node2.lemma)
    
    # L -> 28_sg -> dat. pl. -> L
    feats[40] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28_sg', 'dat. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. pl.', node2.lemma)
    
    # L -> 28_sg -> -273 -> L
    feats[41] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28_sg', '-273') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-273', node2.lemma)
    
    # L -> 28_sg -> 4_sg -> L
    feats[42] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28_sg', '4_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_sg', node2.lemma)
    
    # L -> 28_sg -> 30_tp -> L
    feats[43] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28_sg', '30_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_tp', node2.lemma)
    
    # L -> -12 -> 156 -> L
    feats[44] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-12') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-12', '156') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '156', node2.lemma)
    
    # L -> -12 -> 5_sg -> L
    feats[45] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-12') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-12', '5_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sg', node2.lemma)
    
    # L -> -12 -> -269 -> L
    feats[46] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-12') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-12', '-269') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-269', node2.lemma)
    
    # L -> -12 -> 13_sg -> L
    feats[47] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-12') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-12', '13_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_sg', node2.lemma)
    
    # L -> -12 -> -43 -> L
    feats[48] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-12') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-12', '-43') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-43', node2.lemma)
    
    # L -> -12 -> -308 -> L
    feats[49] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-12') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-12', '-308') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-308', node2.lemma)
    
    # L -> -12 -> 11_pl -> L
    feats[50] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-12') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-12', '11_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_pl', node2.lemma)
    
    # L -> -12 -> -92 -> L
    feats[51] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-12') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-12', '-92') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-92', node2.lemma)
    
    # L -> -12 -> -17 -> L
    feats[52] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-12') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-12', '-17') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-17', node2.lemma)
    
    # L -> -12 -> 8_tp -> L
    feats[53] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-12') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-12', '8_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_tp', node2.lemma)
    
    # L -> -12 -> -83 -> L
    feats[54] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-12') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-12', '-83') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-83', node2.lemma)
    
    # L -> -12 -> sg_sp -> L
    feats[55] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-12') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-12', 'sg_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_sp', node2.lemma)
    
    # L -> -12 -> 42 -> L
    feats[56] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-12') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-12', '42') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '42', node2.lemma)
    
    # L -> -12 -> 96 -> L
    feats[57] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-12') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-12', '96') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '96', node2.lemma)
    
    # L -> -12 -> 134 -> L
    feats[58] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-12') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-12', '134') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '134', node2.lemma)
    
    # L -> -12 -> voc. neutr. -> L
    feats[59] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-12') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-12', 'voc. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. neutr.', node2.lemma)
    
    # L -> -12 -> 139 -> L
    feats[60] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-12') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-12', '139') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '139', node2.lemma)
    
    # L -> -12 -> 136 -> L
    feats[61] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-12') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-12', '136') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '136', node2.lemma)
    
    # L -> -12 -> 174 -> L
    feats[62] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-12') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-12', '174') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '174', node2.lemma)
    
    # L -> -12 -> -133 -> L
    feats[63] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-12') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-12', '-133') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-133', node2.lemma)
    
    # L -> -12 -> 120 -> L
    feats[64] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-12') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-12', '120') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '120', node2.lemma)
    
    # L -> -241 -> 80 -> L
    feats[65] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-241') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-241', '80') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '80', node2.lemma)
    
    # L -> -241 -> 8_pl -> L
    feats[66] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-241') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-241', '8_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_pl', node2.lemma)
    
    # L -> -241 -> 15_sp -> L
    feats[67] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-241') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-241', '15_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_sp', node2.lemma)
    
    # L -> -241 -> 16_du -> L
    feats[68] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-241') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-241', '16_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_du', node2.lemma)
    
    # L -> -241 -> 112 -> L
    feats[69] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-241') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-241', '112') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '112', node2.lemma)
    
    # L -> -241 -> 10_sg -> L
    feats[70] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-241') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-241', '10_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_sg', node2.lemma)
    
    # L -> -241 -> -49 -> L
    feats[71] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-241') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-241', '-49') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-49', node2.lemma)
    
    # L -> -241 -> 8_sp -> L
    feats[72] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-241') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-241', '8_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_sp', node2.lemma)
    
    # L -> -241 -> -56 -> L
    feats[73] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-241') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-241', '-56') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-56', node2.lemma)
    
    # L -> -241 -> 38 -> L
    feats[74] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-241') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-241', '38') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '38', node2.lemma)
    
    # L -> -241 -> 2_du -> L
    feats[75] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-241') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-241', '2_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_du', node2.lemma)
    
    # L -> -241 -> -154 -> L
    feats[76] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-241') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-241', '-154') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-154', node2.lemma)
    
    # L -> -241 -> 96 -> L
    feats[77] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-241') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-241', '96') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '96', node2.lemma)
    
    # L -> -241 -> 90 -> L
    feats[78] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-241') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-241', '90') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '90', node2.lemma)
    
    # L -> -241 -> 8_sg -> L
    feats[79] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-241') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-241', '8_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_sg', node2.lemma)
    
    # L -> -241 -> 13_sp -> L
    feats[80] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-241') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-241', '13_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_sp', node2.lemma)
    
    # L -> -241 -> 136 -> L
    feats[81] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-241') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-241', '136') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '136', node2.lemma)
    
    # L -> -241 -> -133 -> L
    feats[82] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-241') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-241', '-133') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-133', node2.lemma)
    
    # L -> -241 -> 75 -> L
    feats[83] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-241') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-241', '75') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '75', node2.lemma)
    
    # L -> -241 -> 181 -> L
    feats[84] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-241') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-241', '181') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '181', node2.lemma)
    
    # L -> -241 -> -89 -> L
    feats[85] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-241') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-241', '-89') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-89', node2.lemma)
    
    # L -> -241 -> -97 -> L
    feats[86] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-241') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-241', '-97') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-97', node2.lemma)
    
    # L -> -241 -> 14_tp -> L
    feats[87] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-241') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-241', '14_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_tp', node2.lemma)
    
    # L -> -241 -> 92 -> L
    feats[88] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-241') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-241', '92') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '92', node2.lemma)
    
    # L -> -241 -> 6_fp -> L
    feats[89] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-241') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-241', '6_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_fp', node2.lemma)
    
    # L -> -293 -> 11_fp -> L
    feats[90] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '11_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_fp', node2.lemma)
    
    # L -> -293 -> -241 -> L
    feats[91] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '-241') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-241', node2.lemma)
    
    # L -> -293 -> 55 -> L
    feats[92] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '55') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '55', node2.lemma)
    
    # L -> -293 -> -131 -> L
    feats[93] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '-131') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-131', node2.lemma)
    
    # L -> -293 -> 40 -> L
    feats[94] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '40') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '40', node2.lemma)
    
    # L -> -293 -> 6_sp -> L
    feats[95] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '6_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_sp', node2.lemma)
    
    # L -> -293 -> 89 -> L
    feats[96] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '89') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '89', node2.lemma)
    
    # L -> -210 -> -139 -> L
    feats[97] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '-139') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-139', node2.lemma)
    
    # L -> -210 -> 2 -> L
    feats[98] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '2') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2', node2.lemma)
    
    # L -> -210 -> voc. sg. -> L
    feats[99] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', 'voc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. sg.', node2.lemma)
    
    # L -> -210 -> -54 -> L
    feats[100] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '-54') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-54', node2.lemma)
    
    # L -> -210 -> -38 -> L
    feats[101] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '-38') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-38', node2.lemma)
    
    # L -> -210 -> 5_pl -> L
    feats[102] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '5_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_pl', node2.lemma)
    
    # L -> -210 -> 29_tp -> L
    feats[103] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '29_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '29_tp', node2.lemma)
    
    # L -> -210 -> 4_fp -> L
    feats[104] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '4_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_fp', node2.lemma)
    
    # L -> -210 -> 116 -> L
    feats[105] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '116') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '116', node2.lemma)
    
    # L -> -210 -> instr -> L
    feats[106] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', 'instr') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr', node2.lemma)
    
    # L -> -210 -> -86 -> L
    feats[107] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '-86') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-86', node2.lemma)
    
    # L -> -210 -> 61 -> L
    feats[108] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '61') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '61', node2.lemma)
    
    # L -> -210 -> voc. du. -> L
    feats[109] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', 'voc. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. du.', node2.lemma)
    
    # L -> -210 -> 6_fp -> L
    feats[110] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '6_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_fp', node2.lemma)
    
    # L -> -35 -> 11_du -> L
    feats[111] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-35', '11_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_du', node2.lemma)
    
    # L -> -35 -> -276 -> L
    feats[112] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-35', '-276') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-276', node2.lemma)
    
    # L -> -35 -> dat -> L
    feats[113] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-35', 'dat') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat', node2.lemma)
    
    # L -> -35 -> 10_sp -> L
    feats[114] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-35', '10_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_sp', node2.lemma)
    
    # L -> -35 -> 72 -> L
    feats[115] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-35', '72') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '72', node2.lemma)
    
    # L -> -35 -> 28 -> L
    feats[116] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-35', '28') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
    
    # L -> -35 -> 50 -> L
    feats[117] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-35', '50') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '50', node2.lemma)
    
    # L -> -35 -> 2_sg -> L
    feats[118] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-35', '2_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_sg', node2.lemma)
    
    # L -> -35 -> 27_fp -> L
    feats[119] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-35', '27_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_fp', node2.lemma)
    
    # L -> -35 -> 8_du -> L
    feats[120] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-35', '8_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_du', node2.lemma)
    
    # L -> -35 -> fp -> L
    feats[121] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-35', 'fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'fp', node2.lemma)
    
    # L -> -35 -> 12_sp -> L
    feats[122] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-35', '12_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_sp', node2.lemma)
    
    # L -> -35 -> -10 -> L
    feats[123] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-35', '-10') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-10', node2.lemma)
    
    # L -> -35 -> 170 -> L
    feats[124] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-35', '170') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '170', node2.lemma)
    
    # L -> -35 -> 15_sg -> L
    feats[125] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-35', '15_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_sg', node2.lemma)
    
    # L -> -35 -> acc. sg. -> L
    feats[126] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-35', 'acc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. sg.', node2.lemma)
    
    # L -> -35 -> 28_sg -> L
    feats[127] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-35', '28_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28_sg', node2.lemma)
    
    # L -> -35 -> 13_sg -> L
    feats[128] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-35', '13_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_sg', node2.lemma)
    
    # L -> -35 -> 157 -> L
    feats[129] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-35', '157') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '157', node2.lemma)
    
    # L -> -35 -> 15_sp -> L
    feats[130] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-35', '15_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_sp', node2.lemma)
    
    # L -> -35 -> 128 -> L
    feats[131] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-35', '128') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '128', node2.lemma)
    
    # L -> -35 -> 16_du -> L
    feats[132] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-35', '16_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_du', node2.lemma)
    
    # L -> -35 -> 41 -> L
    feats[133] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-35', '41') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '41', node2.lemma)
    
    # L -> -35 -> pl_sp -> L
    feats[134] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-35', 'pl_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_sp', node2.lemma)
    
    # L -> -35 -> 112 -> L
    feats[135] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-35', '112') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '112', node2.lemma)
    
    # L -> -35 -> 16_sg -> L
    feats[136] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-35', '16_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_sg', node2.lemma)
    
    # L -> -35 -> du -> L
    feats[137] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-35', 'du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du', node2.lemma)
    
    # L -> -25 -> 174 -> L
    feats[138] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-25') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-25', '174') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '174', node2.lemma)
    
    # L -> -150 -> 9_sp -> L
    feats[139] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-150') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-150', '9_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_sp', node2.lemma)
    
    # L -> -150 -> 6_tp -> L
    feats[140] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-150') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-150', '6_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_tp', node2.lemma)
    
    # L -> -150 -> -31 -> L
    feats[141] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-150') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-150', '-31') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-31', node2.lemma)
    
    # L -> -150 -> -220 -> L
    feats[142] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-150') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-150', '-220') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-220', node2.lemma)
    
    # L -> -150 -> 2 -> L
    feats[143] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-150') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-150', '2') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2', node2.lemma)
    
    # L -> -150 -> 114 -> L
    feats[144] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-150') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-150', '114') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '114', node2.lemma)
    
    # L -> -150 -> 13_tp -> L
    feats[145] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-150') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-150', '13_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_tp', node2.lemma)
    
    # L -> -150 -> nom -> L
    feats[146] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-150') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-150', 'nom') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom', node2.lemma)
    
    # L -> -150 -> 7_tp -> L
    feats[147] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-150') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-150', '7_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_tp', node2.lemma)
    
    # L -> -150 -> -112 -> L
    feats[148] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-150') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-150', '-112') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    
    # L -> -150 -> -132 -> L
    feats[149] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-150') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-150', '-132') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-132', node2.lemma)
    
    # L -> -150 -> -115 -> L
    feats[150] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-150') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-150', '-115') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-115', node2.lemma)
    
    # L -> -150 -> -261 -> L
    feats[151] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-150') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-150', '-261') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-261', node2.lemma)
    
    # L -> -150 -> 68 -> L
    feats[152] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-150') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-150', '68') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '68', node2.lemma)
    
    # L -> -150 -> du -> L
    feats[153] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-150') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-150', 'du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du', node2.lemma)
    
    # L -> -150 -> 116 -> L
    feats[154] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-150') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-150', '116') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '116', node2.lemma)
    
    # L -> -150 -> du_fp -> L
    feats[155] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-150') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-150', 'du_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_fp', node2.lemma)
    
    # L -> -150 -> tp -> L
    feats[156] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-150') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-150', 'tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'tp', node2.lemma)
    
    # L -> -150 -> 61 -> L
    feats[157] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-150') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-150', '61') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '61', node2.lemma)
    
    # L -> -150 -> 160 -> L
    feats[158] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-150') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-150', '160') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    
    # L -> 6_sg -> nom. neutr. -> L
    feats[159] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', 'nom. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. neutr.', node2.lemma)
    
    # L -> 6_sg -> loc. du. -> L
    feats[160] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', 'loc. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc. du.', node2.lemma)
    
    # L -> 6_sg -> -79 -> L
    feats[161] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', '-79') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-79', node2.lemma)
    
    # L -> 6_sg -> 14_fp -> L
    feats[162] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', '14_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_fp', node2.lemma)
    
    # L -> 6_sg -> -46 -> L
    feats[163] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', '-46') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-46', node2.lemma)
    
    # L -> 6_sg -> 11_fp -> L
    feats[164] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', '11_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_fp', node2.lemma)
    
    # L -> 6_sg -> -242 -> L
    feats[165] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', '-242') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-242', node2.lemma)
    
    # L -> 6_sg -> 60 -> L
    feats[166] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', '60') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '60', node2.lemma)
    
    # L -> 6_sg -> 50 -> L
    feats[167] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', '50') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '50', node2.lemma)
    
    # L -> 6_sg -> -14 -> L
    feats[168] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', '-14') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-14', node2.lemma)
    
    # L -> 6_sg -> -109 -> L
    feats[169] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', '-109') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-109', node2.lemma)
    
    # L -> 6_sg -> 2_sg -> L
    feats[170] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', '2_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_sg', node2.lemma)
    
    # L -> 6_sg -> voc. fem -> L
    feats[171] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', 'voc. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. fem', node2.lemma)
    
    # L -> 6_sg -> 35 -> L
    feats[172] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', '35') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '35', node2.lemma)
    
    # L -> 6_sg -> -271 -> L
    feats[173] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', '-271') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-271', node2.lemma)
    
    # L -> 6_sg -> 12_sp -> L
    feats[174] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', '12_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_sp', node2.lemma)
    
    # L -> 6_sg -> 5_sg -> L
    feats[175] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', '5_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sg', node2.lemma)
    
    # L -> 6_sg -> 54 -> L
    feats[176] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', '54') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    
    # L -> 6_sg -> -10 -> L
    feats[177] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', '-10') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-10', node2.lemma)
    
    # L -> 6_sg -> 15_sg -> L
    feats[178] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', '15_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_sg', node2.lemma)
    
    # L -> 6_sg -> -241 -> L
    feats[179] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', '-241') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-241', node2.lemma)
    
    # L -> 6_sg -> -35 -> L
    feats[180] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', '-35') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    
    # L -> 6_sg -> -25 -> L
    feats[181] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', '-25') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-25', node2.lemma)
    
    # L -> 6_sg -> 110 -> L
    feats[182] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', '110') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '110', node2.lemma)
    
    # L -> 6_sg -> -262 -> L
    feats[183] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', '-262') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-262', node2.lemma)
    
    # L -> 6_sg -> 88 -> L
    feats[184] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', '88') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '88', node2.lemma)
    
    # L -> 6_sg -> -308 -> L
    feats[185] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', '-308') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-308', node2.lemma)
    
    # L -> 6_sg -> 15_sp -> L
    feats[186] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', '15_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_sp', node2.lemma)
    
    # L -> 6_sg -> 55 -> L
    feats[187] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', '55') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '55', node2.lemma)
    
    # L -> 6_sg -> -156 -> L
    feats[188] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', '-156') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-156', node2.lemma)
    
    # L -> 6_sg -> -123 -> L
    feats[189] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', '-123') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-123', node2.lemma)
    
    # L -> 6_sg -> 41 -> L
    feats[190] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', '41') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '41', node2.lemma)
    
    # L -> 6_sg -> pl_sp -> L
    feats[191] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', 'pl_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_sp', node2.lemma)
    
    # L -> 6_sg -> -17 -> L
    feats[192] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', '-17') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-17', node2.lemma)
    
    # L -> 6_sg -> -18 -> L
    feats[193] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', '-18') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-18', node2.lemma)
    
    # L -> 6_sg -> acc. masc. -> L
    feats[194] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', 'acc. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. masc.', node2.lemma)
    
    # L -> 6_sg -> 174 -> L
    feats[195] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', '174') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '174', node2.lemma)
    
    # L -> 6_sg -> -89 -> L
    feats[196] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', '-89') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-89', node2.lemma)
    
    # L -> 6_sg -> -113 -> L
    feats[197] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', '-113') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-113', node2.lemma)
    
    # L -> 6_sg -> 5_fp -> L
    feats[198] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', '5_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_fp', node2.lemma)
    
    # L -> 6_sg -> -266 -> L
    feats[199] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', '-266') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-266', node2.lemma)
    
    # L -> 6_sg -> 91 -> L
    feats[200] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', '91') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    
    # L -> 6_sg -> 30_pl -> L
    feats[201] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', '30_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_pl', node2.lemma)
    
    # L -> 6_sg -> nom. fem -> L
    feats[202] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', 'nom. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. fem', node2.lemma)
    
    # L -> 6_sg -> -21 -> L
    feats[203] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', '-21') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-21', node2.lemma)
    
    # L -> 6_sg -> 9_sg -> L
    feats[204] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', '9_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_sg', node2.lemma)
    
    # L -> loc. sg. -> -309 -> L
    feats[205] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', '-309') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-309', node2.lemma)
    
    # L -> loc. sg. -> -296 -> L
    feats[206] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', '-296') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-296', node2.lemma)
    
    # L -> loc. sg. -> -102 -> L
    feats[207] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', '-102') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-102', node2.lemma)
    
    # L -> loc. sg. -> instr. du. -> L
    feats[208] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', 'instr. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. du.', node2.lemma)
    
    # L -> loc. sg. -> -121 -> L
    feats[209] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', '-121') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-121', node2.lemma)
    
    # L -> loc. sg. -> -33 -> L
    feats[210] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', '-33') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-33', node2.lemma)
    
    # L -> loc. sg. -> instr. masc. -> L
    feats[211] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', 'instr. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. masc.', node2.lemma)
    
    # L -> loc. sg. -> 152 -> L
    feats[212] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', '152') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '152', node2.lemma)
    
    # L -> loc. sg. -> 71 -> L
    feats[213] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', '71') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '71', node2.lemma)
    
    # L -> loc. sg. -> -143 -> L
    feats[214] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', '-143') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-143', node2.lemma)
    
    # L -> loc. sg. -> -152 -> L
    feats[215] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', '-152') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-152', node2.lemma)
    
    # L -> loc. sg. -> -306 -> L
    feats[216] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', '-306') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-306', node2.lemma)
    
    # L -> loc. sg. -> 40 -> L
    feats[217] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', '40') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '40', node2.lemma)
    
    # L -> loc. sg. -> -41 -> L
    feats[218] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', '-41') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-41', node2.lemma)
    
    # L -> loc. sg. -> 12_tp -> L
    feats[219] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', '12_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_tp', node2.lemma)
    
    # L -> loc. sg. -> dat. pl. -> L
    feats[220] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', 'dat. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. pl.', node2.lemma)
    
    # L -> loc. sg. -> -50 -> L
    feats[221] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', '-50') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-50', node2.lemma)
    
    # L -> loc. sg. -> -11 -> L
    feats[222] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', '-11') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-11', node2.lemma)
    
    # L -> loc. sg. -> 82 -> L
    feats[223] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', '82') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '82', node2.lemma)
    
    # L -> loc. sg. -> 10_tp -> L
    feats[224] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', '10_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_tp', node2.lemma)
    
    # L -> loc. sg. -> -78 -> L
    feats[225] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', '-78') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-78', node2.lemma)
    
    # L -> loc. sg. -> 9_sg -> L
    feats[226] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', '9_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_sg', node2.lemma)
    
    # L -> loc. sg. -> 4_sg -> L
    feats[227] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', '4_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_sg', node2.lemma)
    
    # L -> loc. sg. -> -240 -> L
    feats[228] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', '-240') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-240', node2.lemma)
    
    # L -> loc. sg. -> -51 -> L
    feats[229] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', '-51') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-51', node2.lemma)
    
    # L -> loc. sg. -> 12_pl -> L
    feats[230] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', '12_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_pl', node2.lemma)
    
    # L -> loc. sg. -> 93 -> L
    feats[231] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', '93') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '93', node2.lemma)
    
    # L -> loc. sg. -> 138 -> L
    feats[232] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', '138') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '138', node2.lemma)
    
    # L -> loc. sg. -> 101 -> L
    feats[233] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', '101') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    
    # L -> loc. sg. -> -72 -> L
    feats[234] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', '-72') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-72', node2.lemma)
    
    # L -> loc. sg. -> 29_sg -> L
    feats[235] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', '29_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '29_sg', node2.lemma)
    
    # L -> loc. sg. -> -15 -> L
    feats[236] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', '-15') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-15', node2.lemma)
    
    # L -> loc. sg. -> 51 -> L
    feats[237] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', '51') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '51', node2.lemma)
    
    # L -> loc. sg. -> -91 -> L
    feats[238] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', '-91') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-91', node2.lemma)
    
    # L -> loc. sg. -> 59 -> L
    feats[239] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', '59') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '59', node2.lemma)
    
    # L -> loc. sg. -> -103 -> L
    feats[240] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', '-103') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-103', node2.lemma)
    
    # L -> loc. sg. -> 49 -> L
    feats[241] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', '49') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '49', node2.lemma)
    
    # L -> loc. sg. -> 5_du -> L
    feats[242] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', '5_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_du', node2.lemma)
    
    # L -> loc. sg. -> -55 -> L
    feats[243] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', '-55') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-55', node2.lemma)
    
    # L -> loc. sg. -> nom. sg. -> L
    feats[244] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', 'nom. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. sg.', node2.lemma)
    
    # L -> loc. sg. -> -147 -> L
    feats[245] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', '-147') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-147', node2.lemma)
    
    # L -> loc. sg. -> -87 -> L
    feats[246] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', '-87') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-87', node2.lemma)
    
    # L -> loc. sg. -> 56 -> L
    feats[247] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', '56') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '56', node2.lemma)
    
    # L -> loc. sg. -> -117 -> L
    feats[248] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', '-117') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-117', node2.lemma)
    
    # L -> loc. sg. -> 118 -> L
    feats[249] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', '118') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '118', node2.lemma)
    
    # L -> loc. sg. -> -200 -> L
    feats[250] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', '-200') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-200', node2.lemma)
    
    # L -> loc. sg. -> -303 -> L
    feats[251] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', '-303') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-303', node2.lemma)
    
    # L -> loc. sg. -> 11_sp -> L
    feats[252] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', '11_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_sp', node2.lemma)
    
    # L -> loc. sg. -> -69 -> L
    feats[253] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', '-69') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-69', node2.lemma)
    
    # L -> loc. sg. -> -42 -> L
    feats[254] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', '-42') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-42', node2.lemma)
    
    # L -> loc. sg. -> -297 -> L
    feats[255] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', '-297') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-297', node2.lemma)
    
    # L -> loc. sg. -> 9_tp -> L
    feats[256] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', '9_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_tp', node2.lemma)
    
    # L -> 13_sg -> 7_fp -> L
    feats[257] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_sg', '7_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_fp', node2.lemma)
    
    # L -> -149 -> 30_sg -> L
    feats[258] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-149', '30_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_sg', node2.lemma)
    
    # L -> -149 -> 15_sp -> L
    feats[259] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-149', '15_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_sp', node2.lemma)
    
    # L -> -149 -> 128 -> L
    feats[260] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-149', '128') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '128', node2.lemma)
    
    # L -> -149 -> -200 -> L
    feats[261] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-149', '-200') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-200', node2.lemma)
    
    # L -> -149 -> -57 -> L
    feats[262] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-149', '-57') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-57', node2.lemma)
    
    # L -> -149 -> -144 -> L
    feats[263] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-149', '-144') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-144', node2.lemma)
    
    # L -> -149 -> -69 -> L
    feats[264] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-149', '-69') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-69', node2.lemma)
    
    # L -> -149 -> 9_tp -> L
    feats[265] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-149', '9_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_tp', node2.lemma)
    
    # L -> -149 -> 89 -> L
    feats[266] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-149', '89') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '89', node2.lemma)
    
    # L -> -149 -> -90 -> L
    feats[267] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-149', '-90') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-90', node2.lemma)
    
    # L -> -149 -> masc -> L
    feats[268] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-149', 'masc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'masc', node2.lemma)
    
    # L -> -149 -> 27_tp -> L
    feats[269] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-149', '27_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_tp', node2.lemma)
    
    # L -> -149 -> -119 -> L
    feats[270] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-149', '-119') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-119', node2.lemma)
    
    # L -> -149 -> 29 -> L
    feats[271] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-149', '29') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '29', node2.lemma)
    
    # L -> -149 -> 14_pl -> L
    feats[272] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-149', '14_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_pl', node2.lemma)
    
    # L -> -149 -> -31 -> L
    feats[273] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-149', '-31') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-31', node2.lemma)
    
    # L -> -149 -> 173 -> L
    feats[274] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-149', '173') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '173', node2.lemma)
    
    # L -> -149 -> -58 -> L
    feats[275] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-149', '-58') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-58', node2.lemma)
    
    # L -> -149 -> 114 -> L
    feats[276] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-149', '114') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '114', node2.lemma)
    
    # L -> -149 -> acc. neutr. -> L
    feats[277] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-149', 'acc. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. neutr.', node2.lemma)
    
    # L -> -149 -> -159 -> L
    feats[278] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-149', '-159') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-159', node2.lemma)
    
    # L -> -149 -> -73 -> L
    feats[279] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-149', '-73') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-73', node2.lemma)
    
    # L -> -149 -> 81 -> L
    feats[280] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-149', '81') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '81', node2.lemma)
    
    # L -> -149 -> 5_pl -> L
    feats[281] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-149', '5_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_pl', node2.lemma)
    
    # L -> -149 -> -142 -> L
    feats[282] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-149', '-142') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-142', node2.lemma)
    
    # L -> -149 -> 178 -> L
    feats[283] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-149', '178') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '178', node2.lemma)
    
    # L -> -149 -> -112 -> L
    feats[284] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-149', '-112') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    
    # L -> -149 -> -115 -> L
    feats[285] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-149', '-115') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-115', node2.lemma)
    
    # L -> -149 -> -292 -> L
    feats[286] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-149', '-292') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-292', node2.lemma)
    
    # L -> -149 -> 68 -> L
    feats[287] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-149', '68') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '68', node2.lemma)
    
    # L -> -149 -> du_fp -> L
    feats[288] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-149', 'du_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_fp', node2.lemma)
    
    # L -> -149 -> -169 -> L
    feats[289] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-149', '-169') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-169', node2.lemma)
    
    # L -> -149 -> 7_fp -> L
    feats[290] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-149', '7_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_fp', node2.lemma)
    
    # L -> -149 -> voc. du. -> L
    feats[291] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-149', 'voc. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. du.', node2.lemma)
    
    # L -> -149 -> 160 -> L
    feats[292] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-149', '160') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    
    # L -> -149 -> 6_fp -> L
    feats[293] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-149', '6_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_fp', node2.lemma)
    
    # L -> -20 -> nom. neutr. -> L
    feats[294] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-20') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-20', 'nom. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. neutr.', node2.lemma)
    
    # L -> -20 -> -283 -> L
    feats[295] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-20') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-20', '-283') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-283', node2.lemma)
    
    # L -> -20 -> 155 -> L
    feats[296] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-20') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-20', '155') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '155', node2.lemma)
    
    # L -> -20 -> 6_pl -> L
    feats[297] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-20') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-20', '6_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_pl', node2.lemma)
    
    # L -> -20 -> 14_fp -> L
    feats[298] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-20') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-20', '14_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_fp', node2.lemma)
    
    # L -> -20 -> -151 -> L
    feats[299] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-20') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-20', '-151') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-151', node2.lemma)
    
    # L -> -20 -> -46 -> L
    feats[300] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-20') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-20', '-46') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-46', node2.lemma)
    
    # L -> -20 -> instr. pl. -> L
    feats[301] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-20') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-20', 'instr. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. pl.', node2.lemma)
    
    # L -> -20 -> 94 -> L
    feats[302] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-20') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-20', '94') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '94', node2.lemma)
    
    # L -> -20 -> -291 -> L
    feats[303] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-20') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-20', '-291') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-291', node2.lemma)
    
    # L -> -20 -> -14 -> L
    feats[304] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-20') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-20', '-14') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-14', node2.lemma)
    
    # L -> -20 -> -247 -> L
    feats[305] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-20') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-20', '-247') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-247', node2.lemma)
    
    # L -> -20 -> 35 -> L
    feats[306] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-20') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-20', '35') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '35', node2.lemma)
    
    # L -> -20 -> -47 -> L
    feats[307] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-20') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-20', '-47') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-47', node2.lemma)
    
    # L -> -20 -> 153 -> L
    feats[308] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-20') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-20', '153') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    
    # L -> -20 -> 11_sg -> L
    feats[309] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-20') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-20', '11_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_sg', node2.lemma)
    
    # L -> -20 -> 54 -> L
    feats[310] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-20') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-20', '54') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    
    # L -> -20 -> -10 -> L
    feats[311] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-20') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-20', '-10') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-10', node2.lemma)
    
    # L -> -20 -> -71 -> L
    feats[312] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-20') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-20', '-71') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-71', node2.lemma)
    
    # L -> -20 -> 170 -> L
    feats[313] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-20') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-20', '170') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '170', node2.lemma)
    
    # L -> -20 -> 30_sg -> L
    feats[314] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-20') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-20', '30_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_sg', node2.lemma)
    
    # L -> -20 -> 15_sg -> L
    feats[315] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-20') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-20', '15_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_sg', node2.lemma)
    
    # L -> -20 -> acc. sg. -> L
    feats[316] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-20') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-20', 'acc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. sg.', node2.lemma)
    
    # L -> -20 -> -112 -> L
    feats[317] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-20') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-20', '-112') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    
    # L -> -20 -> -132 -> L
    feats[318] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-20') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-20', '-132') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-132', node2.lemma)
    
    # L -> -43 -> 70 -> L
    feats[319] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-43') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-43', '70') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '70', node2.lemma)
    
    # L -> -43 -> -26 -> L
    feats[320] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-43') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-43', '-26') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-26', node2.lemma)
    
    # L -> -43 -> 10_sp -> L
    feats[321] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-43') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-43', '10_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_sp', node2.lemma)
    
    # L -> -43 -> -247 -> L
    feats[322] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-43') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-43', '-247') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-247', node2.lemma)
    
    # L -> -43 -> 12_sp -> L
    feats[323] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-43') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-43', '12_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_sp', node2.lemma)
    
    # L -> -43 -> 88 -> L
    feats[324] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-43') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-43', '88') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '88', node2.lemma)
    
    # L -> -43 -> 80 -> L
    feats[325] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-43') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-43', '80') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '80', node2.lemma)
    
    # L -> -43 -> 8_pl -> L
    feats[326] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-43') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-43', '8_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_pl', node2.lemma)
    
    # L -> -43 -> -77 -> L
    feats[327] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-43') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-43', '-77') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-77', node2.lemma)
    
    # L -> -43 -> 16_sg -> L
    feats[328] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-43') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-43', '16_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_sg', node2.lemma)
    
    # L -> -43 -> -98 -> L
    feats[329] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-43') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-43', '-98') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-98', node2.lemma)
    
    # L -> -43 -> 8_tp -> L
    feats[330] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-43') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-43', '8_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_tp', node2.lemma)
    
    # L -> -43 -> -104 -> L
    feats[331] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-43') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-43', '-104') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-104', node2.lemma)
    
    # L -> -43 -> 90 -> L
    feats[332] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-43') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-43', '90') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '90', node2.lemma)
    
    # L -> -43 -> 30_du -> L
    feats[333] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-43') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-43', '30_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_du', node2.lemma)
    
    # L -> -43 -> -41 -> L
    feats[334] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-43') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-43', '-41') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-41', node2.lemma)
    
    # L -> -43 -> instr. adj. -> L
    feats[335] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-43') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-43', 'instr. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. adj.', node2.lemma)
    
    # L -> -43 -> -240 -> L
    feats[336] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-43') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-43', '-240') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-240', node2.lemma)
    
    # L -> 110 -> 4_fp -> L
    feats[337] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '110') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '110', '4_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_fp', node2.lemma)
    
    # L -> -262 -> sg_tp -> L
    feats[338] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-262') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-262', 'sg_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_tp', node2.lemma)
    
    # L -> -262 -> nom. pl. -> L
    feats[339] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-262') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-262', 'nom. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. pl.', node2.lemma)
    
    # L -> -262 -> 11_sg -> L
    feats[340] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-262') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-262', '11_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_sg', node2.lemma)
    
    # L -> -262 -> 12_sp -> L
    feats[341] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-262') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-262', '12_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_sp', node2.lemma)
    
    # L -> -262 -> -20 -> L
    feats[342] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-262') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-262', '-20') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-20', node2.lemma)
    
    # L -> -262 -> -43 -> L
    feats[343] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-262') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-262', '-43') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-43', node2.lemma)
    
    # L -> -262 -> -77 -> L
    feats[344] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-262') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-262', '-77') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-77', node2.lemma)
    
    # L -> 88 -> dat. du. -> L
    feats[345] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '88') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '88', 'dat. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    
    # L -> 88 -> 16_tp -> L
    feats[346] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '88') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '88', '16_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_tp', node2.lemma)
    
    # L -> 88 -> -91 -> L
    feats[347] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '88') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '88', '-91') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-91', node2.lemma)
    
    # L -> 88 -> -38 -> L
    feats[348] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '88') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '88', '-38') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-38', node2.lemma)
    
    # L -> -308 -> 72 -> L
    feats[349] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '72') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '72', node2.lemma)
    
    # L -> -308 -> 12_sp -> L
    feats[350] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '12_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_sp', node2.lemma)
    
    # L -> -308 -> 10_tp -> L
    feats[351] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '10_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_tp', node2.lemma)
    
    # L -> -308 -> -21 -> L
    feats[352] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '-21') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-21', node2.lemma)
    
    # L -> -308 -> -78 -> L
    feats[353] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '-78') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-78', node2.lemma)
    
    # L -> -308 -> 74 -> L
    feats[354] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '74') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '74', node2.lemma)
    
    # L -> -308 -> -240 -> L
    feats[355] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '-240') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-240', node2.lemma)
    
    # L -> -308 -> 177 -> L
    feats[356] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '177') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '177', node2.lemma)
    
    # L -> -308 -> -15 -> L
    feats[357] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '-15') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-15', node2.lemma)
    
    # L -> -308 -> fem -> L
    feats[358] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', 'fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'fem', node2.lemma)
    
    # L -> -308 -> -53 -> L
    feats[359] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '-53') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-53', node2.lemma)
    
    # L -> -308 -> 14_sp -> L
    feats[360] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '14_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_sp', node2.lemma)
    
    # L -> -308 -> 49 -> L
    feats[361] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '49') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '49', node2.lemma)
    
    # L -> -308 -> voc. pl. -> L
    feats[362] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', 'voc. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. pl.', node2.lemma)
    
    # L -> -308 -> -114 -> L
    feats[363] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '-114') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-114', node2.lemma)
    
    # L -> -308 -> -246 -> L
    feats[364] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '-246') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-246', node2.lemma)
    
    # L -> -308 -> nom. sg. -> L
    feats[365] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', 'nom. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. sg.', node2.lemma)
    
    # L -> -308 -> 56 -> L
    feats[366] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '56') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '56', node2.lemma)
    
    # L -> -308 -> 115 -> L
    feats[367] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '115') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '115', node2.lemma)
    
    # L -> -308 -> acc. du. -> L
    feats[368] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', 'acc. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. du.', node2.lemma)
    
    # L -> -308 -> -230 -> L
    feats[369] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '-230') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-230', node2.lemma)
    
    # L -> -308 -> -42 -> L
    feats[370] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '-42') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-42', node2.lemma)
    
    # L -> -308 -> 76 -> L
    feats[371] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '76') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '76', node2.lemma)
    
    # L -> -308 -> gen. sg. -> L
    feats[372] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', 'gen. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen. sg.', node2.lemma)
    
    # L -> -308 -> 4_tp -> L
    feats[373] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '4_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_tp', node2.lemma)
    
    # L -> -308 -> 9_sp -> L
    feats[374] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '9_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_sp', node2.lemma)
    
    # L -> -308 -> 27_tp -> L
    feats[375] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '27_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_tp', node2.lemma)
    
    # L -> -308 -> -166 -> L
    feats[376] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '-166') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-166', node2.lemma)
    
    # L -> -308 -> 32 -> L
    feats[377] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '32') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '32', node2.lemma)
    
    # L -> -308 -> 114 -> L
    feats[378] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '114') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '114', node2.lemma)
    
    # L -> -308 -> -64 -> L
    feats[379] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '-64') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-64', node2.lemma)
    
    # L -> -308 -> -142 -> L
    feats[380] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '-142') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-142', node2.lemma)
    
    # L -> -308 -> 179 -> L
    feats[381] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '179') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '179', node2.lemma)
    
    # L -> -308 -> 137 -> L
    feats[382] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '137') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    
    # L -> -308 -> -111 -> L
    feats[383] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '-111') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-111', node2.lemma)
    
    # L -> -308 -> 117 -> L
    feats[384] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '117') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '117', node2.lemma)
    
    # L -> -308 -> 116 -> L
    feats[385] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', '116') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '116', node2.lemma)
    
    # L -> 80 -> 170 -> L
    feats[386] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '80') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '80', '170') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '170', node2.lemma)
    
    # L -> 80 -> 73 -> L
    feats[387] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '80') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '80', '73') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '73', node2.lemma)
    
    # L -> 80 -> 1 -> L
    feats[388] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '80') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '80', '1') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '1', node2.lemma)
    
    # L -> 8_pl -> 9_sp -> L
    feats[389] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_pl', '9_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_sp', node2.lemma)
    
    # L -> 8_pl -> -38 -> L
    feats[390] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_pl', '-38') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-38', node2.lemma)
    
    # L -> -77 -> -283 -> L
    feats[391] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-77') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-77', '-283') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-283', node2.lemma)
    
    # L -> -77 -> 94 -> L
    feats[392] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-77') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-77', '94') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '94', node2.lemma)
    
    # L -> -77 -> 50 -> L
    feats[393] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-77') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-77', '50') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '50', node2.lemma)
    
    # L -> -77 -> 2_sg -> L
    feats[394] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-77') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-77', '2_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_sg', node2.lemma)
    
    # L -> -77 -> 29_sg -> L
    feats[395] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-77') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-77', '29_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '29_sg', node2.lemma)
    
    # L -> -77 -> 98 -> L
    feats[396] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-77') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-77', '98') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '98', node2.lemma)
    
    # L -> -77 -> 28_tp -> L
    feats[397] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-77') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-77', '28_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28_tp', node2.lemma)
    
    # L -> -77 -> 2 -> L
    feats[398] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-77') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-77', '2') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2', node2.lemma)
    
    # L -> 11_pl -> -96 -> L
    feats[399] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_pl', '-96') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-96', node2.lemma)
    
    # L -> 11_pl -> -121 -> L
    feats[400] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_pl', '-121') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-121', node2.lemma)
    
    # L -> 11_pl -> -307 -> L
    feats[401] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_pl', '-307') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-307', node2.lemma)
    
    # L -> 11_pl -> 71 -> L
    feats[402] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_pl', '71') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '71', node2.lemma)
    
    # L -> 11_pl -> dat. pl. -> L
    feats[403] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_pl', 'dat. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. pl.', node2.lemma)
    
    # L -> 11_pl -> -51 -> L
    feats[404] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_pl', '-51') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-51', node2.lemma)
    
    # L -> 11_pl -> 30_tp -> L
    feats[405] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_pl', '30_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_tp', node2.lemma)
    
    # L -> 11_pl -> 169 -> L
    feats[406] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_pl', '169') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '169', node2.lemma)
    
    # L -> 11_pl -> 97 -> L
    feats[407] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_pl', '97') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '97', node2.lemma)
    
    # L -> 11_pl -> -114 -> L
    feats[408] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_pl', '-114') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-114', node2.lemma)
    
    # L -> 11_pl -> 148 -> L
    feats[409] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_pl', '148') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '148', node2.lemma)
    
    # L -> 11_pl -> 27_sg -> L
    feats[410] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_pl', '27_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_sg', node2.lemma)
    
    # L -> 157 -> 28 -> L
    feats[411] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '157', '28') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
    
    # L -> 157 -> 54 -> L
    feats[412] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '157', '54') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    
    # L -> 157 -> 117 -> L
    feats[413] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '157', '117') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '117', node2.lemma)
    
    # L -> 157 -> 68 -> L
    feats[414] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '157', '68') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '68', node2.lemma)
    
    # L -> 157 -> -169 -> L
    feats[415] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '157', '-169') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-169', node2.lemma)
    
    # L -> 157 -> 3_du -> L
    feats[416] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '157', '3_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    
    # L -> 157 -> 6_fp -> L
    feats[417] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '157', '6_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_fp', node2.lemma)
    
    # L -> 15_sp -> -276 -> L
    feats[418] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_sp', '-276') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-276', node2.lemma)
    
    # L -> 15_sp -> -151 -> L
    feats[419] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_sp', '-151') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-151', node2.lemma)
    
    # L -> 15_sp -> 11_fp -> L
    feats[420] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_sp', '11_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_fp', node2.lemma)
    
    # L -> 15_sp -> -242 -> L
    feats[421] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_sp', '-242') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-242', node2.lemma)
    
    # L -> 15_sp -> 94 -> L
    feats[422] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_sp', '94') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '94', node2.lemma)
    
    # L -> 15_sp -> -299 -> L
    feats[423] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_sp', '-299') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-299', node2.lemma)
    
    # L -> 15_sp -> -14 -> L
    feats[424] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_sp', '-14') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-14', node2.lemma)
    
    # L -> 15_sp -> -94 -> L
    feats[425] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_sp', '-94') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-94', node2.lemma)
    
    # L -> 15_sp -> -247 -> L
    feats[426] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_sp', '-247') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-247', node2.lemma)
    
    # L -> 15_sp -> -271 -> L
    feats[427] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_sp', '-271') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-271', node2.lemma)
    
    # L -> 15_sp -> 54 -> L
    feats[428] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_sp', '54') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    
    # L -> 15_sp -> 30_sg -> L
    feats[429] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_sp', '30_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_sg', node2.lemma)
    
    # L -> 15_sp -> 28_sg -> L
    feats[430] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_sp', '28_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28_sg', node2.lemma)
    
    # L -> 15_sp -> -150 -> L
    feats[431] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_sp', '-150') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-150', node2.lemma)
    
    # L -> 15_sp -> 88 -> L
    feats[432] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_sp', '88') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '88', node2.lemma)
    
    # L -> 15_sp -> -308 -> L
    feats[433] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_sp', '-308') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-308', node2.lemma)
    
    # L -> 15_sp -> -77 -> L
    feats[434] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_sp', '-77') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-77', node2.lemma)
    
    # L -> 15_sp -> -156 -> L
    feats[435] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_sp', '-156') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-156', node2.lemma)
    
    # L -> 15_sp -> -268 -> L
    feats[436] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_sp', '-268') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-268', node2.lemma)
    
    # L -> 15_sp -> 37 -> L
    feats[437] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_sp', '37') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '37', node2.lemma)
    
    # L -> 15_sp -> -123 -> L
    feats[438] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_sp', '-123') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-123', node2.lemma)
    
    # L -> 15_sp -> 5_sp -> L
    feats[439] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_sp', '5_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sp', node2.lemma)
    
    # L -> 15_sp -> 132 -> L
    feats[440] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_sp', '132') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '132', node2.lemma)
    
    # L -> 15_sp -> nom. du. -> L
    feats[441] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_sp', 'nom. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. du.', node2.lemma)
    
    # L -> 15_sp -> 58 -> L
    feats[442] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_sp', '58') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '58', node2.lemma)
    
    # L -> 15_sp -> abl. du. -> L
    feats[443] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_sp', 'abl. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. du.', node2.lemma)
    
    # L -> 15_sp -> voc -> L
    feats[444] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_sp', 'voc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc', node2.lemma)
    
    # L -> 15_sp -> dat. pl. -> L
    feats[445] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_sp', 'dat. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. pl.', node2.lemma)
    
    # L -> 15_sp -> 4_sp -> L
    feats[446] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_sp', '4_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_sp', node2.lemma)
    
    # L -> 15_sp -> 4_sg -> L
    feats[447] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_sp', '4_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_sg', node2.lemma)
    
    # L -> 15_sp -> 59 -> L
    feats[448] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_sp', '59') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '59', node2.lemma)
    
    # L -> 15_sp -> -69 -> L
    feats[449] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_sp', '-69') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-69', node2.lemma)
    
    # L -> 15_sp -> masc -> L
    feats[450] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_sp', 'masc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'masc', node2.lemma)
    
    # L -> 15_sp -> -132 -> L
    feats[451] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_sp', '-132') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-132', node2.lemma)
    
    # L -> 128 -> sg_tp -> L
    feats[452] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '128') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '128', 'sg_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_tp', node2.lemma)
    
    # L -> 128 -> 155 -> L
    feats[453] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '128') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '128', '155') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '155', node2.lemma)
    
    # L -> 128 -> 5_sg -> L
    feats[454] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '128') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '128', '5_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sg', node2.lemma)
    
    # L -> 128 -> 37 -> L
    feats[455] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '128') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '128', '37') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '37', node2.lemma)
    
    # L -> 128 -> -56 -> L
    feats[456] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '128') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '128', '-56') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-56', node2.lemma)
    
    # L -> 128 -> 58 -> L
    feats[457] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '128') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '128', '58') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '58', node2.lemma)
    
    # L -> 128 -> 122 -> L
    feats[458] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '128') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '128', '122') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '122', node2.lemma)
    
    # L -> 128 -> 8_sg -> L
    feats[459] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '128') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '128', '8_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_sg', node2.lemma)
    
    # L -> 128 -> 34 -> L
    feats[460] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '128') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '128', '34') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '34', node2.lemma)
    
    # L -> 128 -> 161 -> L
    feats[461] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '128') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '128', '161') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '161', node2.lemma)
    
    # L -> 128 -> 14_tp -> L
    feats[462] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '128') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '128', '14_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_tp', node2.lemma)
    
    # L -> -92 -> 176 -> L
    feats[463] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '176') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '176', node2.lemma)
    
    # L -> -92 -> 7_sp -> L
    feats[464] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '7_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_sp', node2.lemma)
    
    # L -> -92 -> 96 -> L
    feats[465] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '96') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '96', node2.lemma)
    
    # L -> -92 -> 140 -> L
    feats[466] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '140') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '140', node2.lemma)
    
    # L -> -92 -> -84 -> L
    feats[467] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '-84') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-84', node2.lemma)
    
    # L -> -92 -> 122 -> L
    feats[468] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '122') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '122', node2.lemma)
    
    # L -> -92 -> 174 -> L
    feats[469] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '174') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '174', node2.lemma)
    
    # L -> -92 -> -23 -> L
    feats[470] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '-23') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-23', node2.lemma)
    
    # L -> -92 -> 120 -> L
    feats[471] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '120') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '120', node2.lemma)
    
    # L -> -92 -> dat. sg. -> L
    feats[472] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', 'dat. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. sg.', node2.lemma)
    
    # L -> -92 -> -131 -> L
    feats[473] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '-131') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-131', node2.lemma)
    
    # L -> -92 -> 34 -> L
    feats[474] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '34') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '34', node2.lemma)
    
    # L -> -92 -> 102 -> L
    feats[475] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '102') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '102', node2.lemma)
    
    # L -> -92 -> du_tp -> L
    feats[476] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', 'du_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_tp', node2.lemma)
    
    # L -> -92 -> pl -> L
    feats[477] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', 'pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    
    # L -> -92 -> 119 -> L
    feats[478] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '119') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '119', node2.lemma)
    
    # L -> -92 -> -296 -> L
    feats[479] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '-296') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-296', node2.lemma)
    
    # L -> -92 -> -67 -> L
    feats[480] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '-67') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-67', node2.lemma)
    
    # L -> -92 -> -102 -> L
    feats[481] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '-102') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-102', node2.lemma)
    
    # L -> -92 -> -28 -> L
    feats[482] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '-28') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-28', node2.lemma)
    
    # L -> -92 -> -33 -> L
    feats[483] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '-33') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-33', node2.lemma)
    
    # L -> -92 -> 182 -> L
    feats[484] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '182') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '182', node2.lemma)
    
    # L -> -92 -> 159 -> L
    feats[485] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '159') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '159', node2.lemma)
    
    # L -> -92 -> -36 -> L
    feats[486] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '-36') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-36', node2.lemma)
    
    # L -> -92 -> 5_fp -> L
    feats[487] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '5_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_fp', node2.lemma)
    
    # L -> -92 -> -39 -> L
    feats[488] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '-39') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-39', node2.lemma)
    
    # L -> -92 -> -249 -> L
    feats[489] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '-249') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-249', node2.lemma)
    
    # L -> -92 -> -306 -> L
    feats[490] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '-306') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-306', node2.lemma)
    
    # L -> -92 -> 40 -> L
    feats[491] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '40') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '40', node2.lemma)
    
    # L -> -92 -> 15_fp -> L
    feats[492] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '15_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_fp', node2.lemma)
    
    # L -> -92 -> -153 -> L
    feats[493] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '-153') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-153', node2.lemma)
    
    # L -> -92 -> 30_pl -> L
    feats[494] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '30_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_pl', node2.lemma)
    
    # L -> -92 -> 6_du -> L
    feats[495] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '6_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_du', node2.lemma)
    
    # L -> -92 -> 4_sp -> L
    feats[496] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '4_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_sp', node2.lemma)
    
    # L -> -92 -> 6_sp -> L
    feats[497] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '6_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_sp', node2.lemma)
    
    # L -> -92 -> instr. adj. -> L
    feats[498] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', 'instr. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. adj.', node2.lemma)
    
    # L -> -92 -> -82 -> L
    feats[499] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '-82') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-82', node2.lemma)
    
    # L -> -92 -> -51 -> L
    feats[500] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '-51') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-51', node2.lemma)
    
    # L -> -92 -> 30_tp -> L
    feats[501] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '30_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_tp', node2.lemma)
    
    # L -> -92 -> 177 -> L
    feats[502] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '177') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '177', node2.lemma)
    
    # L -> -92 -> -302 -> L
    feats[503] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '-302') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-302', node2.lemma)
    
    # L -> -92 -> voc. pl. -> L
    feats[504] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', 'voc. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. pl.', node2.lemma)
    
    # L -> -92 -> -279 -> L
    feats[505] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '-279') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-279', node2.lemma)
    
    # L -> -92 -> 9_tp -> L
    feats[506] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '9_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_tp', node2.lemma)
    
    # L -> -92 -> 172 -> L
    feats[507] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '172') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '172', node2.lemma)
    
    # L -> -92 -> -119 -> L
    feats[508] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '-119') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-119', node2.lemma)
    
    # L -> -92 -> -58 -> L
    feats[509] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '-58') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-58', node2.lemma)
    
    # L -> -92 -> acc. neutr. -> L
    feats[510] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', 'acc. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. neutr.', node2.lemma)
    
    # L -> -92 -> 68 -> L
    feats[511] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '68') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '68', node2.lemma)
    
    # L -> -92 -> du -> L
    feats[512] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', 'du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du', node2.lemma)
    
    # L -> -92 -> 116 -> L
    feats[513] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-92', '116') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '116', node2.lemma)
    
    # L -> 30_sp -> -79 -> L
    feats[514] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_sp', '-79') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-79', node2.lemma)
    
    # L -> 30_sp -> -37 -> L
    feats[515] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_sp', '-37') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-37', node2.lemma)
    
    # L -> 30_sp -> 11_fp -> L
    feats[516] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_sp', '11_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_fp', node2.lemma)
    
    # L -> 30_sp -> 50 -> L
    feats[517] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_sp', '50') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '50', node2.lemma)
    
    # L -> 30_sp -> -299 -> L
    feats[518] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_sp', '-299') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-299', node2.lemma)
    
    # L -> 55 -> -34 -> L
    feats[519] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '55') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '55', '-34') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-34', node2.lemma)
    
    # L -> 55 -> -142 -> L
    feats[520] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '55') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '55', '-142') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-142', node2.lemma)
    
    # L -> -156 -> -48 -> L
    feats[521] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-156', '-48') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-48', node2.lemma)
    
    # L -> -156 -> 173 -> L
    feats[522] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-156', '173') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '173', node2.lemma)
    
    # L -> -156 -> 95 -> L
    feats[523] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-156', '95') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '95', node2.lemma)
    
    # L -> -156 -> -159 -> L
    feats[524] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-156', '-159') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-159', node2.lemma)
    
    # L -> -156 -> -73 -> L
    feats[525] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-156', '-73') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-73', node2.lemma)
    
    # L -> -156 -> 8_fp -> L
    feats[526] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-156', '8_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    
    # L -> -156 -> 81 -> L
    feats[527] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-156', '81') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '81', node2.lemma)
    
    # L -> -156 -> -24 -> L
    feats[528] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-156', '-24') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    
    # L -> -156 -> -132 -> L
    feats[529] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-156', '-132') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-132', node2.lemma)
    
    # L -> -156 -> -292 -> L
    feats[530] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-156', '-292') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-292', node2.lemma)
    
    # L -> -156 -> du -> L
    feats[531] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-156', 'du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du', node2.lemma)
    
    # L -> -156 -> 116 -> L
    feats[532] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-156', '116') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '116', node2.lemma)
    
    # L -> -156 -> -22 -> L
    feats[533] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-156', '-22') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-22', node2.lemma)
    
    # L -> -156 -> du_fp -> L
    feats[534] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-156', 'du_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_fp', node2.lemma)
    
    # L -> -156 -> instr -> L
    feats[535] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-156', 'instr') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr', node2.lemma)
    
    # L -> -156 -> -86 -> L
    feats[536] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-156', '-86') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-86', node2.lemma)
    
    # L -> -156 -> 61 -> L
    feats[537] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-156', '61') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '61', node2.lemma)
    
    # L -> -156 -> voc. du. -> L
    feats[538] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-156', 'voc. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. du.', node2.lemma)
    
    # L -> -156 -> 6_fp -> L
    feats[539] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-156', '6_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_fp', node2.lemma)
    
    # L -> 16_du -> 11_du -> L
    feats[540] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_du', '11_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_du', node2.lemma)
    
    # L -> 16_du -> -276 -> L
    feats[541] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_du', '-276') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-276', node2.lemma)
    
    # L -> 16_du -> dat -> L
    feats[542] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_du', 'dat') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat', node2.lemma)
    
    # L -> 16_du -> 11_fp -> L
    feats[543] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_du', '11_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_fp', node2.lemma)
    
    # L -> 16_du -> nom. adj. -> L
    feats[544] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_du', 'nom. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. adj.', node2.lemma)
    
    # L -> 16_du -> -242 -> L
    feats[545] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_du', '-242') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-242', node2.lemma)
    
    # L -> 16_du -> 60 -> L
    feats[546] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_du', '60') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '60', node2.lemma)
    
    # L -> 16_du -> -299 -> L
    feats[547] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_du', '-299') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-299', node2.lemma)
    
    # L -> 16_du -> -109 -> L
    feats[548] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_du', '-109') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-109', node2.lemma)
    
    # L -> 16_du -> voc. fem -> L
    feats[549] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_du', 'voc. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. fem', node2.lemma)
    
    # L -> 16_du -> -247 -> L
    feats[550] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_du', '-247') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-247', node2.lemma)
    
    # L -> 16_du -> 153 -> L
    feats[551] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_du', '153') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    
    # L -> 16_du -> -271 -> L
    feats[552] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_du', '-271') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-271', node2.lemma)
    
    # L -> 16_du -> fp -> L
    feats[553] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_du', 'fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'fp', node2.lemma)
    
    # L -> 16_du -> nom. pl. -> L
    feats[554] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_du', 'nom. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. pl.', node2.lemma)
    
    # L -> 16_du -> 156 -> L
    feats[555] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_du', '156') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '156', node2.lemma)
    
    # L -> 16_du -> acc -> L
    feats[556] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_du', 'acc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc', node2.lemma)
    
    # L -> 16_du -> 30_sg -> L
    feats[557] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_du', '30_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_sg', node2.lemma)
    
    # L -> 16_du -> acc. sg. -> L
    feats[558] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_du', 'acc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. sg.', node2.lemma)
    
    # L -> 16_du -> 28_sg -> L
    feats[559] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_du', '28_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28_sg', node2.lemma)
    
    # L -> 16_du -> -12 -> L
    feats[560] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_du', '-12') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-12', node2.lemma)
    
    # L -> 16_du -> -20 -> L
    feats[561] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_du', '-20') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-20', node2.lemma)
    
    # L -> 16_du -> -262 -> L
    feats[562] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_du', '-262') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-262', node2.lemma)
    
    # L -> 16_du -> 88 -> L
    feats[563] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_du', '88') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '88', node2.lemma)
    
    # L -> 16_du -> 128 -> L
    feats[564] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_du', '128') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '128', node2.lemma)
    
    # L -> 16_du -> -92 -> L
    feats[565] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_du', '-92') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-92', node2.lemma)
    
    # L -> 16_du -> 30_sp -> L
    feats[566] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_du', '30_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_sp', node2.lemma)
    
    # L -> 16_du -> 129 -> L
    feats[567] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_du', '129') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '129', node2.lemma)
    
    # L -> 16_du -> -123 -> L
    feats[568] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_du', '-123') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-123', node2.lemma)
    
    # L -> 16_du -> -139 -> L
    feats[569] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_du', '-139') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-139', node2.lemma)
    
    # L -> 16_du -> 112 -> L
    feats[570] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_du', '112') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '112', node2.lemma)
    
    # L -> 16_du -> 8_tp -> L
    feats[571] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_du', '8_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_tp', node2.lemma)
    
    # L -> 16_du -> -49 -> L
    feats[572] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_du', '-49') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-49', node2.lemma)
    
    # L -> 16_du -> 9_du -> L
    feats[573] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_du', '9_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_du', node2.lemma)
    
    # L -> 16_du -> 2_tp -> L
    feats[574] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_du', '2_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_tp', node2.lemma)
    
    # L -> 16_du -> -45 -> L
    feats[575] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_du', '-45') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-45', node2.lemma)
    
    # L -> 16_du -> -122 -> L
    feats[576] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_du', '-122') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-122', node2.lemma)
    
    # L -> 16_du -> -61 -> L
    feats[577] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_du', '-61') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-61', node2.lemma)
    
    # L -> 16_du -> -152 -> L
    feats[578] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_du', '-152') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-152', node2.lemma)
    
    # L -> 16_du -> 30_fp -> L
    feats[579] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_du', '30_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_fp', node2.lemma)
    
    # L -> 16_du -> 10_du -> L
    feats[580] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_du', '10_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_du', node2.lemma)
    
    # L -> 16_du -> -273 -> L
    feats[581] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_du', '-273') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-273', node2.lemma)
    
    # L -> 16_du -> acc. du. -> L
    feats[582] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_du', 'acc. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. du.', node2.lemma)
    
    # L -> 16_du -> 27_sg -> L
    feats[583] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_du', '27_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_sg', node2.lemma)
    
    # L -> 16_du -> -111 -> L
    feats[584] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_du', '-111') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-111', node2.lemma)
    
    # L -> -268 -> instr. pl. -> L
    feats[585] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', 'instr. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. pl.', node2.lemma)
    
    # L -> -268 -> 30_sg -> L
    feats[586] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', '30_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_sg', node2.lemma)
    
    # L -> -268 -> -308 -> L
    feats[587] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', '-308') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-308', node2.lemma)
    
    # L -> -268 -> 38 -> L
    feats[588] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', '38') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '38', node2.lemma)
    
    # L -> -268 -> 140 -> L
    feats[589] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', '140') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '140', node2.lemma)
    
    # L -> -268 -> acc. adj. -> L
    feats[590] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', 'acc. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. adj.', node2.lemma)
    
    # L -> -268 -> 91 -> L
    feats[591] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', '91') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    
    # L -> -268 -> -69 -> L
    feats[592] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', '-69') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-69', node2.lemma)
    
    # L -> -268 -> 2_sp -> L
    feats[593] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', '2_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_sp', node2.lemma)
    
    # L -> -268 -> 33 -> L
    feats[594] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', '33') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '33', node2.lemma)
    
    # L -> -268 -> loc. pl. -> L
    feats[595] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', 'loc. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc. pl.', node2.lemma)
    
    # L -> -268 -> -90 -> L
    feats[596] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', '-90') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-90', node2.lemma)
    
    # L -> -268 -> -245 -> L
    feats[597] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', '-245') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-245', node2.lemma)
    
    # L -> -268 -> -126 -> L
    feats[598] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', '-126') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-126', node2.lemma)
    
    # L -> -268 -> 9_sp -> L
    feats[599] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', '9_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_sp', node2.lemma)
    
    # L -> -268 -> -166 -> L
    feats[600] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', '-166') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-166', node2.lemma)
    
    # L -> -268 -> 14_pl -> L
    feats[601] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', '14_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_pl', node2.lemma)
    
    # L -> -268 -> -31 -> L
    feats[602] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', '-31') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-31', node2.lemma)
    
    # L -> -268 -> 95 -> L
    feats[603] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', '95') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '95', node2.lemma)
    
    # L -> -268 -> acc. neutr. -> L
    feats[604] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', 'acc. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. neutr.', node2.lemma)
    
    # L -> -268 -> -73 -> L
    feats[605] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', '-73') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-73', node2.lemma)
    
    # L -> -268 -> 27_sg -> L
    feats[606] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', '27_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_sg', node2.lemma)
    
    # L -> -268 -> -129 -> L
    feats[607] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', '-129') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-129', node2.lemma)
    
    # L -> -268 -> -142 -> L
    feats[608] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', '-142') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-142', node2.lemma)
    
    # L -> -268 -> 150 -> L
    feats[609] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', '150') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '150', node2.lemma)
    
    # L -> -268 -> -132 -> L
    feats[610] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', '-132') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-132', node2.lemma)
    
    # L -> -268 -> -261 -> L
    feats[611] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', '-261') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-261', node2.lemma)
    
    # L -> -268 -> -292 -> L
    feats[612] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', '-292') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-292', node2.lemma)
    
    # L -> -268 -> 4_fp -> L
    feats[613] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', '4_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_fp', node2.lemma)
    
    # L -> -268 -> 117 -> L
    feats[614] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', '117') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '117', node2.lemma)
    
    # L -> -268 -> -86 -> L
    feats[615] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', '-86') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-86', node2.lemma)
    
    # L -> -268 -> 61 -> L
    feats[616] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', '61') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '61', node2.lemma)
    
    # L -> -268 -> voc. du. -> L
    feats[617] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', 'voc. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. du.', node2.lemma)
    
    # L -> 37 -> 11_du -> L
    feats[618] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '37', '11_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_du', node2.lemma)
    
    # L -> 37 -> -299 -> L
    feats[619] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '37', '-299') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-299', node2.lemma)
    
    # L -> 37 -> -291 -> L
    feats[620] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '37', '-291') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-291', node2.lemma)
    
    # L -> 37 -> -71 -> L
    feats[621] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '37', '-71') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-71', node2.lemma)
    
    # L -> 37 -> 28_sg -> L
    feats[622] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '37', '28_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28_sg', node2.lemma)
    
    # L -> 37 -> -210 -> L
    feats[623] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '37', '-210') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-210', node2.lemma)
    
    # L -> 37 -> 6_sg -> L
    feats[624] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '37', '6_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_sg', node2.lemma)
    
    # L -> 37 -> -20 -> L
    feats[625] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '37', '-20') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-20', node2.lemma)
    
    # L -> 37 -> 8_pl -> L
    feats[626] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '37', '8_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_pl', node2.lemma)
    
    # L -> 37 -> 128 -> L
    feats[627] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '37', '128') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '128', node2.lemma)
    
    # L -> 37 -> -123 -> L
    feats[628] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '37', '-123') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-123', node2.lemma)
    
    # L -> 37 -> -260 -> L
    feats[629] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '37', '-260') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-260', node2.lemma)
    
    # L -> 37 -> 176 -> L
    feats[630] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '37', '176') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '176', node2.lemma)
    
    # L -> 37 -> 141 -> L
    feats[631] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '37', '141') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '141', node2.lemma)
    
    # L -> 37 -> 36 -> L
    feats[632] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '37', '36') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '36', node2.lemma)
    
    # L -> 37 -> 9_fp -> L
    feats[633] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '37', '9_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_fp', node2.lemma)
    
    # L -> 37 -> 13_sp -> L
    feats[634] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '37', '13_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_sp', node2.lemma)
    
    # L -> 37 -> 168 -> L
    feats[635] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '37', '168') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '168', node2.lemma)
    
    # L -> 37 -> 13_fp -> L
    feats[636] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '37', '13_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_fp', node2.lemma)
    
    # L -> 37 -> 181 -> L
    feats[637] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '37', '181') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '181', node2.lemma)
    
    # L -> 37 -> 4_pl -> L
    feats[638] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '37', '4_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_pl', node2.lemma)
    
    # L -> 37 -> -11 -> L
    feats[639] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '37', '-11') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-11', node2.lemma)
    
    # L -> 37 -> 9_sg -> L
    feats[640] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '37', '9_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_sg', node2.lemma)
    
    # L -> 37 -> -114 -> L
    feats[641] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '37', '-114') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-114', node2.lemma)
    
    # L -> 37 -> -16 -> L
    feats[642] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '37', '-16') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-16', node2.lemma)
    
    # L -> 37 -> 98 -> L
    feats[643] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '37', '98') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '98', node2.lemma)
    
    # L -> 37 -> acc. du. -> L
    feats[644] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '37', 'acc. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. du.', node2.lemma)
    
    # L -> 37 -> -200 -> L
    feats[645] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '37', '-200') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-200', node2.lemma)
    
    # L -> 37 -> -76 -> L
    feats[646] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '37', '-76') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-76', node2.lemma)
    
    # L -> 41 -> -10 -> L
    feats[647] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '41') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '41', '-10') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-10', node2.lemma)
    
    # L -> 41 -> -20 -> L
    feats[648] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '41') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '41', '-20') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-20', node2.lemma)
    
    # L -> 41 -> 132 -> L
    feats[649] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '41') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '41', '132') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '132', node2.lemma)
    
    # L -> 41 -> 140 -> L
    feats[650] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '41') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '41', '140') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '140', node2.lemma)
    
    # L -> pl_sp -> 12_fp -> L
    feats[651] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_sp', '12_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_fp', node2.lemma)
    
    # L -> pl_sp -> 3 -> L
    feats[652] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_sp', '3') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3', node2.lemma)
    
    # L -> pl_sp -> 13_sp -> L
    feats[653] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_sp', '13_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_sp', node2.lemma)
    
    # L -> pl_sp -> 174 -> L
    feats[654] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_sp', '174') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '174', node2.lemma)
    
    # L -> pl_sp -> -190 -> L
    feats[655] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_sp', '-190') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-190', node2.lemma)
    
    # L -> pl_sp -> 120 -> L
    feats[656] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_sp', '120') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '120', node2.lemma)
    
    # L -> pl_sp -> -122 -> L
    feats[657] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_sp', '-122') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-122', node2.lemma)
    
    # L -> pl_sp -> -131 -> L
    feats[658] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_sp', '-131') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-131', node2.lemma)
    
    # L -> pl_sp -> 102 -> L
    feats[659] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_sp', '102') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '102', node2.lemma)
    
    # L -> pl_sp -> -89 -> L
    feats[660] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_sp', '-89') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-89', node2.lemma)
    
    # L -> pl_sp -> -29 -> L
    feats[661] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_sp', '-29') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-29', node2.lemma)
    
    # L -> pl_sp -> 119 -> L
    feats[662] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_sp', '119') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '119', node2.lemma)
    
    # L -> pl_sp -> -309 -> L
    feats[663] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_sp', '-309') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-309', node2.lemma)
    
    # L -> pl_sp -> 100 -> L
    feats[664] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_sp', '100') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '100', node2.lemma)
    
    # L -> pl_sp -> abl. pl. -> L
    feats[665] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_sp', 'abl. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. pl.', node2.lemma)
    
    # L -> pl_sp -> -33 -> L
    feats[666] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_sp', '-33') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-33', node2.lemma)
    
    # L -> pl_sp -> 182 -> L
    feats[667] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_sp', '182') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '182', node2.lemma)
    
    # L -> pl_sp -> sp -> L
    feats[668] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_sp', 'sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sp', node2.lemma)
    
    # L -> pl_sp -> 154 -> L
    feats[669] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_sp', '154') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '154', node2.lemma)
    
    # L -> pl_sp -> 14_sg -> L
    feats[670] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_sp', '14_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_sg', node2.lemma)
    
    # L -> pl_sp -> 10_pl -> L
    feats[671] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_sp', '10_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_pl', node2.lemma)
    
    # L -> pl_sp -> instr. fem -> L
    feats[672] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_sp', 'instr. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. fem', node2.lemma)
    
    # L -> pl_sp -> 152 -> L
    feats[673] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_sp', '152') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '152', node2.lemma)
    
    # L -> pl_sp -> 71 -> L
    feats[674] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_sp', '71') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '71', node2.lemma)
    
    # L -> -18 -> -279 -> L
    feats[675] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '-279') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-279', node2.lemma)
    
    # L -> -18 -> 56 -> L
    feats[676] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '56') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '56', node2.lemma)
    
    # L -> -18 -> -34 -> L
    feats[677] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '-34') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-34', node2.lemma)
    
    # L -> -18 -> -263 -> L
    feats[678] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '-263') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-263', node2.lemma)
    
    # L -> -18 -> 11_sp -> L
    feats[679] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '11_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_sp', node2.lemma)
    
    # L -> -18 -> -68 -> L
    feats[680] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '-68') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-68', node2.lemma)
    
    # L -> -18 -> -69 -> L
    feats[681] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '-69') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-69', node2.lemma)
    
    # L -> -18 -> 9_tp -> L
    feats[682] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '9_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_tp', node2.lemma)
    
    # L -> -18 -> gen. sg. -> L
    feats[683] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', 'gen. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen. sg.', node2.lemma)
    
    # L -> -18 -> -126 -> L
    feats[684] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '-126') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-126', node2.lemma)
    
    # L -> -18 -> 9_sp -> L
    feats[685] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '9_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_sp', node2.lemma)
    
    # L -> -18 -> 32 -> L
    feats[686] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '32') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '32', node2.lemma)
    
    # L -> -18 -> -119 -> L
    feats[687] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '-119') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-119', node2.lemma)
    
    # L -> -18 -> 29 -> L
    feats[688] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '29') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '29', node2.lemma)
    
    # L -> -18 -> 14_pl -> L
    feats[689] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '14_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_pl', node2.lemma)
    
    # L -> -18 -> -31 -> L
    feats[690] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '-31') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-31', node2.lemma)
    
    # L -> -18 -> -220 -> L
    feats[691] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '-220') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-220', node2.lemma)
    
    # L -> -18 -> 2 -> L
    feats[692] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '2') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2', node2.lemma)
    
    # L -> -18 -> voc. sg. -> L
    feats[693] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', 'voc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. sg.', node2.lemma)
    
    # L -> -18 -> -58 -> L
    feats[694] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '-58') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-58', node2.lemma)
    
    # L -> -18 -> 114 -> L
    feats[695] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '114') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '114', node2.lemma)
    
    # L -> -18 -> -158 -> L
    feats[696] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '-158') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-158', node2.lemma)
    
    # L -> -18 -> -159 -> L
    feats[697] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '-159') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-159', node2.lemma)
    
    # L -> -18 -> 27_sg -> L
    feats[698] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '27_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_sg', node2.lemma)
    
    # L -> -18 -> 81 -> L
    feats[699] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '81') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '81', node2.lemma)
    
    # L -> -18 -> -129 -> L
    feats[700] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '-129') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-129', node2.lemma)
    
    # L -> -18 -> 5_pl -> L
    feats[701] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '5_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_pl', node2.lemma)
    
    # L -> -18 -> -64 -> L
    feats[702] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '-64') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-64', node2.lemma)
    
    # L -> -18 -> 178 -> L
    feats[703] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '178') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '178', node2.lemma)
    
    # L -> -18 -> -112 -> L
    feats[704] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '-112') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    
    # L -> -18 -> -115 -> L
    feats[705] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '-115') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-115', node2.lemma)
    
    # L -> -18 -> -292 -> L
    feats[706] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '-292') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-292', node2.lemma)
    
    # L -> -18 -> 117 -> L
    feats[707] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '117') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '117', node2.lemma)
    
    # L -> -18 -> 68 -> L
    feats[708] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '68') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '68', node2.lemma)
    
    # L -> -18 -> tp -> L
    feats[709] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', 'tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'tp', node2.lemma)
    
    # L -> -18 -> instr -> L
    feats[710] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', 'instr') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr', node2.lemma)
    
    # L -> -18 -> 3_du -> L
    feats[711] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '3_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    
    # L -> -18 -> -86 -> L
    feats[712] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '-86') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-86', node2.lemma)
    
    # L -> -18 -> 160 -> L
    feats[713] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '160') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    
    # L -> -18 -> 6_fp -> L
    feats[714] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '6_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_fp', node2.lemma)
    
    # L -> -139 -> 11_du -> L
    feats[715] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '11_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_du', node2.lemma)
    
    # L -> -139 -> -283 -> L
    feats[716] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-283') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-283', node2.lemma)
    
    # L -> -139 -> -301 -> L
    feats[717] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-301') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-301', node2.lemma)
    
    # L -> -139 -> -276 -> L
    feats[718] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-276') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-276', node2.lemma)
    
    # L -> -139 -> 6_pl -> L
    feats[719] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '6_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_pl', node2.lemma)
    
    # L -> -139 -> loc. du. -> L
    feats[720] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', 'loc. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc. du.', node2.lemma)
    
    # L -> -139 -> -37 -> L
    feats[721] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-37') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-37', node2.lemma)
    
    # L -> -139 -> -26 -> L
    feats[722] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-26') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-26', node2.lemma)
    
    # L -> -139 -> -46 -> L
    feats[723] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-46') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-46', node2.lemma)
    
    # L -> -139 -> instr. pl. -> L
    feats[724] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', 'instr. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. pl.', node2.lemma)
    
    # L -> -139 -> 10_sp -> L
    feats[725] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '10_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_sp', node2.lemma)
    
    # L -> -139 -> 28 -> L
    feats[726] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '28') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
    
    # L -> -139 -> -299 -> L
    feats[727] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-299') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-299', node2.lemma)
    
    # L -> -139 -> -109 -> L
    feats[728] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-109') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-109', node2.lemma)
    
    # L -> -139 -> -81 -> L
    feats[729] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-81') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-81', node2.lemma)
    
    # L -> -139 -> -59 -> L
    feats[730] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-59') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-59', node2.lemma)
    
    # L -> -139 -> -271 -> L
    feats[731] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-271') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-271', node2.lemma)
    
    # L -> -139 -> 8_du -> L
    feats[732] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '8_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_du', node2.lemma)
    
    # L -> -139 -> fp -> L
    feats[733] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', 'fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'fp', node2.lemma)
    
    # L -> -139 -> 11_sg -> L
    feats[734] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '11_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_sg', node2.lemma)
    
    # L -> -139 -> 156 -> L
    feats[735] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '156') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '156', node2.lemma)
    
    # L -> -139 -> 54 -> L
    feats[736] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '54') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    
    # L -> -139 -> -71 -> L
    feats[737] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-71') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-71', node2.lemma)
    
    # L -> -139 -> 170 -> L
    feats[738] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '170') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '170', node2.lemma)
    
    # L -> -139 -> acc. sg. -> L
    feats[739] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', 'acc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. sg.', node2.lemma)
    
    # L -> -139 -> 28_sg -> L
    feats[740] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '28_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28_sg', node2.lemma)
    
    # L -> -139 -> -12 -> L
    feats[741] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-12') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-12', node2.lemma)
    
    # L -> -139 -> -241 -> L
    feats[742] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-241') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-241', node2.lemma)
    
    # L -> -139 -> -293 -> L
    feats[743] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-293') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-293', node2.lemma)
    
    # L -> -139 -> -210 -> L
    feats[744] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-210') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-210', node2.lemma)
    
    # L -> -139 -> 6_sp -> L
    feats[745] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '6_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_sp', node2.lemma)
    
    # L -> -139 -> 138 -> L
    feats[746] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '138') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '138', node2.lemma)
    
    # L -> -139 -> 131 -> L
    feats[747] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '131') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '131', node2.lemma)
    
    # L -> -139 -> 89 -> L
    feats[748] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '89') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '89', node2.lemma)
    
    # L -> 112 -> 9_fp -> L
    feats[749] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '112') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '112', '9_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_fp', node2.lemma)
    
    # L -> 112 -> 3_pl -> L
    feats[750] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '112') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '112', '3_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_pl', node2.lemma)
    
    # L -> 112 -> 139 -> L
    feats[751] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '112') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '112', '139') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '139', node2.lemma)
    
    # L -> 112 -> 90 -> L
    feats[752] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '112') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '112', '90') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '90', node2.lemma)
    
    # L -> 112 -> 99 -> L
    feats[753] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '112') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '112', '99') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '99', node2.lemma)
    
    # L -> 112 -> 12_fp -> L
    feats[754] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '112') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '112', '12_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_fp', node2.lemma)
    
    # L -> 112 -> -133 -> L
    feats[755] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '112') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '112', '-133') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-133', node2.lemma)
    
    # L -> 112 -> -23 -> L
    feats[756] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '112') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '112', '-23') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-23', node2.lemma)
    
    # L -> 112 -> dat. sg. -> L
    feats[757] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '112') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '112', 'dat. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. sg.', node2.lemma)
    
    # L -> 112 -> -131 -> L
    feats[758] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '112') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '112', '-131') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-131', node2.lemma)
    
    # L -> 112 -> 34 -> L
    feats[759] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '112') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '112', '34') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '34', node2.lemma)
    
    # L -> 112 -> 102 -> L
    feats[760] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '112') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '112', '102') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '102', node2.lemma)
    
    # L -> 112 -> du_tp -> L
    feats[761] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '112') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '112', 'du_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_tp', node2.lemma)
    
    # L -> 112 -> pl -> L
    feats[762] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '112') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '112', 'pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    
    # L -> 112 -> 119 -> L
    feats[763] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '112') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '112', '119') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '119', node2.lemma)
    
    # L -> 112 -> -309 -> L
    feats[764] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '112') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '112', '-309') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-309', node2.lemma)
    
    # L -> 112 -> abl. pl. -> L
    feats[765] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '112') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '112', 'abl. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. pl.', node2.lemma)
    
    # L -> 112 -> -307 -> L
    feats[766] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '112') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '112', '-307') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-307', node2.lemma)
    
    # L -> 112 -> 159 -> L
    feats[767] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '112') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '112', '159') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '159', node2.lemma)
    
    # L -> 112 -> 69 -> L
    feats[768] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '112') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '112', '69') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '69', node2.lemma)
    
    # L -> 112 -> instr. fem -> L
    feats[769] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '112') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '112', 'instr. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. fem', node2.lemma)
    
    # L -> 112 -> 71 -> L
    feats[770] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '112') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '112', '71') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '71', node2.lemma)
    
    # L -> 112 -> -152 -> L
    feats[771] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '112') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '112', '-152') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-152', node2.lemma)
    
    # L -> 112 -> -306 -> L
    feats[772] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '112') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '112', '-306') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-306', node2.lemma)
    
    # L -> 112 -> 2_pl -> L
    feats[773] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '112') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '112', '2_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_pl', node2.lemma)
    
    # L -> 112 -> 10_tp -> L
    feats[774] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '112') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '112', '10_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_tp', node2.lemma)
    
    # L -> 112 -> -21 -> L
    feats[775] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '112') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '112', '-21') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-21', node2.lemma)
    
    # L -> 112 -> -30 -> L
    feats[776] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '112') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '112', '-30') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-30', node2.lemma)
    
    # L -> 112 -> instr. adj. -> L
    feats[777] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '112') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '112', 'instr. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. adj.', node2.lemma)
    
    # L -> 112 -> 5_tp -> L
    feats[778] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '112') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '112', '5_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_tp', node2.lemma)
    
    # L -> 151 -> 55 -> L
    feats[779] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '151') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '151', '55') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '55', node2.lemma)
    
    # L -> 10_sg -> 155 -> L
    feats[780] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sg', '155') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '155', node2.lemma)
    
    # L -> 10_sg -> -276 -> L
    feats[781] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sg', '-276') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-276', node2.lemma)
    
    # L -> 10_sg -> -37 -> L
    feats[782] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sg', '-37') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-37', node2.lemma)
    
    # L -> 10_sg -> 14_fp -> L
    feats[783] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sg', '14_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_fp', node2.lemma)
    
    # L -> 10_sg -> 11_fp -> L
    feats[784] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sg', '11_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_fp', node2.lemma)
    
    # L -> 10_sg -> 94 -> L
    feats[785] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sg', '94') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '94', node2.lemma)
    
    # L -> 10_sg -> 35 -> L
    feats[786] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sg', '35') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '35', node2.lemma)
    
    # L -> 10_sg -> 8_du -> L
    feats[787] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sg', '8_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_du', node2.lemma)
    
    # L -> 10_sg -> 12_sp -> L
    feats[788] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sg', '12_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_sp', node2.lemma)
    
    # L -> 10_sg -> -10 -> L
    feats[789] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sg', '-10') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-10', node2.lemma)
    
    # L -> 10_sg -> -269 -> L
    feats[790] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sg', '-269') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-269', node2.lemma)
    
    # L -> 10_sg -> -71 -> L
    feats[791] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sg', '-71') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-71', node2.lemma)
    
    # L -> 10_sg -> 170 -> L
    feats[792] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sg', '170') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '170', node2.lemma)
    
    # L -> 10_sg -> -241 -> L
    feats[793] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sg', '-241') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-241', node2.lemma)
    
    # L -> 10_sg -> -150 -> L
    feats[794] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sg', '-150') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-150', node2.lemma)
    
    # L -> 10_sg -> 88 -> L
    feats[795] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sg', '88') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '88', node2.lemma)
    
    # L -> 10_sg -> -77 -> L
    feats[796] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sg', '-77') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-77', node2.lemma)
    
    # L -> 10_sg -> 15_sp -> L
    feats[797] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sg', '15_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_sp', node2.lemma)
    
    # L -> 10_sg -> 30_sp -> L
    feats[798] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sg', '30_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_sp', node2.lemma)
    
    # L -> 10_sg -> 41 -> L
    feats[799] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sg', '41') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '41', node2.lemma)
    
    # L -> 10_sg -> -98 -> L
    feats[800] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sg', '-98') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-98', node2.lemma)
    
    # L -> 10_sg -> 8_tp -> L
    feats[801] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sg', '8_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_tp', node2.lemma)
    
    # L -> 10_sg -> 2_tp -> L
    feats[802] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sg', '2_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_tp', node2.lemma)
    
    # L -> 10_sg -> abl -> L
    feats[803] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sg', 'abl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl', node2.lemma)
    
    # L -> 10_sg -> 132 -> L
    feats[804] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sg', '132') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '132', node2.lemma)
    
    # L -> 10_sg -> 15_tp -> L
    feats[805] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sg', '15_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_tp', node2.lemma)
    
    # L -> 10_sg -> 38 -> L
    feats[806] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sg', '38') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '38', node2.lemma)
    
    # L -> 10_sg -> sg_sp -> L
    feats[807] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sg', 'sg_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_sp', node2.lemma)
    
    # L -> 10_sg -> -154 -> L
    feats[808] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sg', '-154') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-154', node2.lemma)
    
    # L -> 10_sg -> 139 -> L
    feats[809] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sg', '139') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '139', node2.lemma)
    
    # L -> 10_sg -> 90 -> L
    feats[810] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sg', '90') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '90', node2.lemma)
    
    # L -> 10_sg -> 1 -> L
    feats[811] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sg', '1') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '1', node2.lemma)
    
    # L -> 10_sg -> 3_sp -> L
    feats[812] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sg', '3_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_sp', node2.lemma)
    
    # L -> 10_sg -> 152 -> L
    feats[813] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sg', '152') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '152', node2.lemma)
    
    # L -> 10_sg -> 30_pl -> L
    feats[814] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sg', '30_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_pl', node2.lemma)
    
    # L -> 10_sg -> 3_tp -> L
    feats[815] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sg', '3_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_tp', node2.lemma)
    
    # L -> 10_sg -> 148 -> L
    feats[816] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sg', '148') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '148', node2.lemma)
    
    # L -> 10_sg -> nom. sg. -> L
    feats[817] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sg', 'nom. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. sg.', node2.lemma)
    
    # L -> 16_sg -> -14 -> L
    feats[818] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_sg', '-14') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-14', node2.lemma)
    
    # L -> 16_sg -> 35 -> L
    feats[819] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_sg', '35') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '35', node2.lemma)
    
    # L -> 16_sg -> -47 -> L
    feats[820] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_sg', '-47') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-47', node2.lemma)
    
    # L -> 16_sg -> 8_du -> L
    feats[821] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_sg', '8_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_du', node2.lemma)
    
    # L -> 16_sg -> 54 -> L
    feats[822] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_sg', '54') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    
    # L -> 16_sg -> 28_sg -> L
    feats[823] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_sg', '28_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28_sg', node2.lemma)
    
    # L -> 16_sg -> -210 -> L
    feats[824] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_sg', '-210') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-210', node2.lemma)
    
    # L -> 16_sg -> 6_sg -> L
    feats[825] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_sg', '6_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_sg', node2.lemma)
    
    # L -> 16_sg -> -20 -> L
    feats[826] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_sg', '-20') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-20', node2.lemma)
    
    # L -> 16_sg -> 110 -> L
    feats[827] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_sg', '110') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '110', node2.lemma)
    
    # L -> 16_sg -> 88 -> L
    feats[828] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_sg', '88') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '88', node2.lemma)
    
    # L -> 16_sg -> -308 -> L
    feats[829] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_sg', '-308') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-308', node2.lemma)
    
    # L -> 16_sg -> 15_sp -> L
    feats[830] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_sg', '15_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_sp', node2.lemma)
    
    # L -> 16_sg -> 55 -> L
    feats[831] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_sg', '55') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '55', node2.lemma)
    
    # L -> 16_sg -> -17 -> L
    feats[832] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_sg', '-17') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-17', node2.lemma)
    
    # L -> 16_sg -> 151 -> L
    feats[833] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_sg', '151') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '151', node2.lemma)
    
    # L -> 16_sg -> 9_du -> L
    feats[834] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_sg', '9_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_du', node2.lemma)
    
    # L -> 16_sg -> 15_tp -> L
    feats[835] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_sg', '15_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_tp', node2.lemma)
    
    # L -> 16_sg -> sg_sp -> L
    feats[836] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_sg', 'sg_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_sp', node2.lemma)
    
    # L -> 16_sg -> -260 -> L
    feats[837] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_sg', '-260') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-260', node2.lemma)
    
    # L -> 16_sg -> 176 -> L
    feats[838] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_sg', '176') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '176', node2.lemma)
    
    # L -> 16_sg -> 141 -> L
    feats[839] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_sg', '141') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '141', node2.lemma)
    
    # L -> 16_sg -> 36 -> L
    feats[840] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_sg', '36') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '36', node2.lemma)
    
    # L -> 16_sg -> 9_fp -> L
    feats[841] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_sg', '9_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_fp', node2.lemma)
    
    # L -> 16_sg -> abl. du. -> L
    feats[842] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_sg', 'abl. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. du.', node2.lemma)
    
    # L -> 16_sg -> 12_fp -> L
    feats[843] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_sg', '12_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_fp', node2.lemma)
    
    # L -> 16_sg -> 3 -> L
    feats[844] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_sg', '3') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3', node2.lemma)
    
    # L -> 16_sg -> 13_fp -> L
    feats[845] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_sg', '13_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_fp', node2.lemma)
    
    # L -> 16_sg -> dat. sg. -> L
    feats[846] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_sg', 'dat. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. sg.', node2.lemma)
    
    # L -> 16_sg -> 2_pl -> L
    feats[847] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_sg', '2_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_pl', node2.lemma)
    
    # L -> 8_tp -> 6_sg -> L
    feats[848] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_tp', '6_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_sg', node2.lemma)
    
    # L -> 8_tp -> 176 -> L
    feats[849] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_tp', '176') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '176', node2.lemma)
    
    # L -> 8_tp -> 3 -> L
    feats[850] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_tp', '3') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3', node2.lemma)
    
    # L -> 8_tp -> dat. pl. -> L
    feats[851] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_tp', 'dat. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. pl.', node2.lemma)
    
    # L -> 8_tp -> 5_du -> L
    feats[852] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_tp', '5_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_du', node2.lemma)
    
    # L -> 8_tp -> 115 -> L
    feats[853] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_tp', '115') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '115', node2.lemma)
    
    # L -> -49 -> -150 -> L
    feats[854] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-49') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-49', '-150') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-150', node2.lemma)
    
    # L -> -49 -> 89 -> L
    feats[855] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-49') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-49', '89') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '89', node2.lemma)
    
    # L -> -49 -> 172 -> L
    feats[856] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-49') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-49', '172') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '172', node2.lemma)
    
    # L -> -49 -> 175 -> L
    feats[857] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-49') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-49', '175') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '175', node2.lemma)
    
    # L -> -49 -> 7_pl -> L
    feats[858] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-49') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-49', '7_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_pl', node2.lemma)
    
    # L -> -49 -> masc -> L
    feats[859] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-49') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-49', 'masc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'masc', node2.lemma)
    
    # L -> -49 -> -119 -> L
    feats[860] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-49') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-49', '-119') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-119', node2.lemma)
    
    # L -> -49 -> 114 -> L
    feats[861] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-49') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-49', '114') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '114', node2.lemma)
    
    # L -> -49 -> -73 -> L
    feats[862] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-49') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-49', '-73') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-73', node2.lemma)
    
    # L -> -49 -> -38 -> L
    feats[863] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-49') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-49', '-38') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-38', node2.lemma)
    
    # L -> -49 -> 13_tp -> L
    feats[864] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-49') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-49', '13_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_tp', node2.lemma)
    
    # L -> -49 -> 5_pl -> L
    feats[865] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-49') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-49', '5_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_pl', node2.lemma)
    
    # L -> -49 -> -13 -> L
    feats[866] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-49') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-49', '-13') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-13', node2.lemma)
    
    # L -> -49 -> -64 -> L
    feats[867] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-49') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-49', '-64') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-64', node2.lemma)
    
    # L -> -49 -> 179 -> L
    feats[868] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-49') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-49', '179') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '179', node2.lemma)
    
    # L -> -49 -> -137 -> L
    feats[869] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-49') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-49', '-137') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-137', node2.lemma)
    
    # L -> -49 -> 178 -> L
    feats[870] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-49') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-49', '178') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '178', node2.lemma)
    
    # L -> -49 -> -111 -> L
    feats[871] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-49') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-49', '-111') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-111', node2.lemma)
    
    # L -> -49 -> -22 -> L
    feats[872] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-49') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-49', '-22') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-22', node2.lemma)
    
    # L -> -49 -> -169 -> L
    feats[873] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-49') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-49', '-169') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-169', node2.lemma)
    
    # L -> -49 -> -86 -> L
    feats[874] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-49') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-49', '-86') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-86', node2.lemma)
    
    # L -> -49 -> 160 -> L
    feats[875] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-49') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-49', '160') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    
    # L -> 9_du -> 155 -> L
    feats[876] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '155') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '155', node2.lemma)
    
    # L -> 9_du -> loc. du. -> L
    feats[877] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', 'loc. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc. du.', node2.lemma)
    
    # L -> 9_du -> 14_fp -> L
    feats[878] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '14_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_fp', node2.lemma)
    
    # L -> 9_du -> -26 -> L
    feats[879] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '-26') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-26', node2.lemma)
    
    # L -> 9_du -> instr. pl. -> L
    feats[880] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', 'instr. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. pl.', node2.lemma)
    
    # L -> 9_du -> nom. adj. -> L
    feats[881] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', 'nom. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. adj.', node2.lemma)
    
    # L -> 9_du -> 60 -> L
    feats[882] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '60') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '60', node2.lemma)
    
    # L -> 9_du -> 50 -> L
    feats[883] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '50') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '50', node2.lemma)
    
    # L -> 9_du -> -299 -> L
    feats[884] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '-299') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-299', node2.lemma)
    
    # L -> 9_du -> -291 -> L
    feats[885] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '-291') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-291', node2.lemma)
    
    # L -> 9_du -> -47 -> L
    feats[886] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '-47') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-47', node2.lemma)
    
    # L -> 9_du -> 54 -> L
    feats[887] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '54') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    
    # L -> 9_du -> -10 -> L
    feats[888] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '-10') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-10', node2.lemma)
    
    # L -> 9_du -> acc -> L
    feats[889] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', 'acc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc', node2.lemma)
    
    # L -> 9_du -> -269 -> L
    feats[890] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '-269') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-269', node2.lemma)
    
    # L -> 9_du -> acc. sg. -> L
    feats[891] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', 'acc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. sg.', node2.lemma)
    
    # L -> 9_du -> -293 -> L
    feats[892] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '-293') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-293', node2.lemma)
    
    # L -> 9_du -> -210 -> L
    feats[893] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '-210') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-210', node2.lemma)
    
    # L -> 9_du -> -149 -> L
    feats[894] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '-149') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-149', node2.lemma)
    
    # L -> 9_du -> -92 -> L
    feats[895] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '-92') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-92', node2.lemma)
    
    # L -> 9_du -> -123 -> L
    feats[896] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '-123') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-123', node2.lemma)
    
    # L -> 9_du -> 41 -> L
    feats[897] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '41') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '41', node2.lemma)
    
    # L -> 9_du -> -139 -> L
    feats[898] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '-139') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-139', node2.lemma)
    
    # L -> 9_du -> -49 -> L
    feats[899] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '-49') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-49', node2.lemma)
    
    # L -> 9_du -> 5_sp -> L
    feats[900] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '5_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sp', node2.lemma)
    
    # L -> 9_du -> 38 -> L
    feats[901] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '38') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '38', node2.lemma)
    
    # L -> 9_du -> 58 -> L
    feats[902] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '58') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '58', node2.lemma)
    
    # L -> 9_du -> 3_fp -> L
    feats[903] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '3_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_fp', node2.lemma)
    
    # L -> 9_du -> 12_sg -> L
    feats[904] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '12_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_sg', node2.lemma)
    
    # L -> 9_du -> -154 -> L
    feats[905] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '-154') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-154', node2.lemma)
    
    # L -> 9_du -> -45 -> L
    feats[906] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '-45') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-45', node2.lemma)
    
    # L -> 9_du -> 3_pl -> L
    feats[907] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '3_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_pl', node2.lemma)
    
    # L -> 9_du -> 134 -> L
    feats[908] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '134') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '134', node2.lemma)
    
    # L -> 9_du -> 99 -> L
    feats[909] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '99') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '99', node2.lemma)
    
    # L -> 9_du -> 122 -> L
    feats[910] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '122') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '122', node2.lemma)
    
    # L -> 9_du -> 27_pl -> L
    feats[911] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '27_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_pl', node2.lemma)
    
    # L -> 9_du -> 120 -> L
    feats[912] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '120') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '120', node2.lemma)
    
    # L -> 9_du -> -44 -> L
    feats[913] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '-44') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-44', node2.lemma)
    
    # L -> 9_du -> 161 -> L
    feats[914] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '161') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '161', node2.lemma)
    
    # L -> 9_du -> -89 -> L
    feats[915] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '-89') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-89', node2.lemma)
    
    # L -> 9_du -> 14_tp -> L
    feats[916] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '14_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_tp', node2.lemma)
    
    # L -> 9_du -> pl -> L
    feats[917] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', 'pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    
    # L -> 9_du -> voc -> L
    feats[918] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', 'voc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc', node2.lemma)
    
    # L -> 9_du -> 149 -> L
    feats[919] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '149') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '149', node2.lemma)
    
    # L -> 9_du -> instr. du. -> L
    feats[920] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', 'instr. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. du.', node2.lemma)
    
    # L -> 9_du -> abl. pl. -> L
    feats[921] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', 'abl. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. pl.', node2.lemma)
    
    # L -> 9_du -> -28 -> L
    feats[922] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '-28') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-28', node2.lemma)
    
    # L -> 9_du -> -307 -> L
    feats[923] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '-307') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-307', node2.lemma)
    
    # L -> 9_du -> 159 -> L
    feats[924] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '159') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '159', node2.lemma)
    
    # L -> 9_du -> 5_fp -> L
    feats[925] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '5_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_fp', node2.lemma)
    
    # L -> 9_du -> -11 -> L
    feats[926] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '-11') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-11', node2.lemma)
    
    # L -> 9_du -> -101 -> L
    feats[927] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '-101') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-101', node2.lemma)
    
    # L -> 9_du -> 82 -> L
    feats[928] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '82') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '82', node2.lemma)
    
    # L -> 9_du -> 16_tp -> L
    feats[929] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '16_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_tp', node2.lemma)
    
    # L -> 9_du -> nom. masc. -> L
    feats[930] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', 'nom. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. masc.', node2.lemma)
    
    # L -> 9_du -> pl_fp -> L
    feats[931] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', 'pl_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_fp', node2.lemma)
    
    # L -> 9_du -> -15 -> L
    feats[932] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '-15') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-15', node2.lemma)
    
    # L -> 9_du -> fem -> L
    feats[933] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', 'fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'fem', node2.lemma)
    
    # L -> 9_du -> 14_sp -> L
    feats[934] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '14_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_sp', node2.lemma)
    
    # L -> 9_du -> 158 -> L
    feats[935] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '158') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '158', node2.lemma)
    
    # L -> 9_du -> -200 -> L
    feats[936] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '-200') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-200', node2.lemma)
    
    # L -> 9_du -> -303 -> L
    feats[937] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '-303') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-303', node2.lemma)
    
    # L -> 9_du -> -68 -> L
    feats[938] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '-68') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-68', node2.lemma)
    
    # L -> 9_du -> -76 -> L
    feats[939] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '-76') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-76', node2.lemma)
    
    # L -> 9_du -> -144 -> L
    feats[940] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '-144') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-144', node2.lemma)
    
    # L -> 9_du -> 180 -> L
    feats[941] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '180') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '180', node2.lemma)
    
    # L -> 9_du -> masc -> L
    feats[942] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', 'masc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'masc', node2.lemma)
    
    # L -> 9_du -> 32 -> L
    feats[943] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', '32') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '32', node2.lemma)
    
    # L -> 121 -> 14_fp -> L
    feats[944] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '121') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '121', '14_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_fp', node2.lemma)
    
    # L -> 121 -> loc. sg. -> L
    feats[945] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '121') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '121', 'loc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc. sg.', node2.lemma)
    
    # L -> 121 -> 55 -> L
    feats[946] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '121') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '121', '55') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '55', node2.lemma)
    
    # L -> -104 -> 15_sp -> L
    feats[947] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '15_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_sp', node2.lemma)
    
    # L -> -104 -> 13_pl -> L
    feats[948] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '13_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_pl', node2.lemma)
    
    # L -> -104 -> 100 -> L
    feats[949] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '100') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '100', node2.lemma)
    
    # L -> -104 -> instr. du. -> L
    feats[950] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', 'instr. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. du.', node2.lemma)
    
    # L -> -104 -> abl. pl. -> L
    feats[951] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', 'abl. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. pl.', node2.lemma)
    
    # L -> -104 -> 78 -> L
    feats[952] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '78') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '78', node2.lemma)
    
    # L -> -104 -> 182 -> L
    feats[953] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '182') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '182', node2.lemma)
    
    # L -> -104 -> 14_sg -> L
    feats[954] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '14_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_sg', node2.lemma)
    
    # L -> -104 -> 2_fp -> L
    feats[955] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '2_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_fp', node2.lemma)
    
    # L -> -104 -> 10_pl -> L
    feats[956] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '10_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_pl', node2.lemma)
    
    # L -> -104 -> -306 -> L
    feats[957] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '-306') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-306', node2.lemma)
    
    # L -> -104 -> -153 -> L
    feats[958] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '-153') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-153', node2.lemma)
    
    # L -> -104 -> dat. pl. -> L
    feats[959] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', 'dat. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. pl.', node2.lemma)
    
    # L -> -104 -> 7_du -> L
    feats[960] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '7_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_du', node2.lemma)
    
    # L -> -104 -> -50 -> L
    feats[961] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '-50') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-50', node2.lemma)
    
    # L -> -104 -> -11 -> L
    feats[962] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '-11') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-11', node2.lemma)
    
    # L -> -104 -> -101 -> L
    feats[963] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '-101') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-101', node2.lemma)
    
    # L -> -104 -> 77 -> L
    feats[964] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '77') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '77', node2.lemma)
    
    # L -> -104 -> nom. fem -> L
    feats[965] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', 'nom. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. fem', node2.lemma)
    
    # L -> -104 -> 82 -> L
    feats[966] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '82') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '82', node2.lemma)
    
    # L -> -104 -> -48 -> L
    feats[967] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '-48') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-48', node2.lemma)
    
    # L -> -104 -> 6_du -> L
    feats[968] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '6_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_du', node2.lemma)
    
    # L -> -104 -> du_sp -> L
    feats[969] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', 'du_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_sp', node2.lemma)
    
    # L -> -104 -> 10_tp -> L
    feats[970] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '10_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_tp', node2.lemma)
    
    # L -> -104 -> 6_sp -> L
    feats[971] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '6_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_sp', node2.lemma)
    
    # L -> -104 -> instr. adj. -> L
    feats[972] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', 'instr. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. adj.', node2.lemma)
    
    # L -> -104 -> -78 -> L
    feats[973] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '-78') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-78', node2.lemma)
    
    # L -> -104 -> 74 -> L
    feats[974] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '74') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '74', node2.lemma)
    
    # L -> -104 -> -273 -> L
    feats[975] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '-273') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-273', node2.lemma)
    
    # L -> -104 -> 130 -> L
    feats[976] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '130') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '130', node2.lemma)
    
    # L -> -104 -> abl. sg. -> L
    feats[977] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', 'abl. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. sg.', node2.lemma)
    
    # L -> -104 -> 12_pl -> L
    feats[978] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '12_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_pl', node2.lemma)
    
    # L -> -104 -> nom. masc. -> L
    feats[979] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', 'nom. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. masc.', node2.lemma)
    
    # L -> -104 -> 138 -> L
    feats[980] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '138') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '138', node2.lemma)
    
    # L -> -104 -> 111 -> L
    feats[981] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '111') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '111', node2.lemma)
    
    # L -> -104 -> 101 -> L
    feats[982] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '101') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    
    # L -> -104 -> -15 -> L
    feats[983] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '-15') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-15', node2.lemma)
    
    # L -> -104 -> -52 -> L
    feats[984] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '-52') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-52', node2.lemma)
    
    # L -> -104 -> 51 -> L
    feats[985] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '51') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '51', node2.lemma)
    
    # L -> -104 -> fem -> L
    feats[986] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', 'fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'fem', node2.lemma)
    
    # L -> -104 -> 15_pl -> L
    feats[987] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '15_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_pl', node2.lemma)
    
    # L -> -104 -> -53 -> L
    feats[988] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '-53') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-53', node2.lemma)
    
    # L -> -104 -> 169 -> L
    feats[989] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '169') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '169', node2.lemma)
    
    # L -> -104 -> 59 -> L
    feats[990] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '59') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '59', node2.lemma)
    
    # L -> -104 -> 14_sp -> L
    feats[991] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '14_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_sp', node2.lemma)
    
    # L -> -104 -> 5_du -> L
    feats[992] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '5_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_du', node2.lemma)
    
    # L -> -104 -> -16 -> L
    feats[993] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '-16') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-16', node2.lemma)
    
    # L -> -104 -> 98 -> L
    feats[994] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '98') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '98', node2.lemma)
    
    # L -> -104 -> 158 -> L
    feats[995] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '158') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '158', node2.lemma)
    
    # L -> -104 -> 11_tp -> L
    feats[996] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '11_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_tp', node2.lemma)
    
    # L -> -104 -> adj -> L
    feats[997] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', 'adj') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'adj', node2.lemma)
    
    # L -> -104 -> 148 -> L
    feats[998] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '148') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '148', node2.lemma)
    
    # L -> -104 -> nom. sg. -> L
    feats[999] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', 'nom. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. sg.', node2.lemma)
    
    # L -> -104 -> 118 -> L
    feats[1000] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '118') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '118', node2.lemma)
    
    # L -> -104 -> -263 -> L
    feats[1001] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '-263') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-263', node2.lemma)
    
    # L -> -104 -> acc. du. -> L
    feats[1002] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', 'acc. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. du.', node2.lemma)
    
    # L -> -104 -> -76 -> L
    feats[1003] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '-76') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-76', node2.lemma)
    
    # L -> -104 -> -144 -> L
    feats[1004] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '-144') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-144', node2.lemma)
    
    # L -> -104 -> -69 -> L
    feats[1005] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '-69') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-69', node2.lemma)
    
    # L -> -104 -> -42 -> L
    feats[1006] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '-42') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-42', node2.lemma)
    
    # L -> -104 -> 9_tp -> L
    feats[1007] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '9_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_tp', node2.lemma)
    
    # L -> -104 -> 33 -> L
    feats[1008] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '33') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '33', node2.lemma)
    
    # L -> -104 -> 172 -> L
    feats[1009] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '172') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '172', node2.lemma)
    
    # L -> -104 -> 7_sg -> L
    feats[1010] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '7_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_sg', node2.lemma)
    
    # L -> -104 -> -90 -> L
    feats[1011] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '-90') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-90', node2.lemma)
    
    # L -> -104 -> voc. masc. -> L
    feats[1012] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', 'voc. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. masc.', node2.lemma)
    
    # L -> -104 -> -245 -> L
    feats[1013] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '-245') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-245', node2.lemma)
    
    # L -> -104 -> 4_tp -> L
    feats[1014] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '4_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_tp', node2.lemma)
    
    # L -> -104 -> -126 -> L
    feats[1015] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '-126') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-126', node2.lemma)
    
    # L -> -104 -> 180 -> L
    feats[1016] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', '180') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '180', node2.lemma)
    
    # L -> -104 -> masc -> L
    feats[1017] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', 'masc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'masc', node2.lemma)
    
    # L -> -56 -> 132 -> L
    feats[1018] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-56') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-56', '132') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '132', node2.lemma)
    
    # L -> -56 -> 13_pl -> L
    feats[1019] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-56') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-56', '13_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_pl', node2.lemma)
    
    # L -> -56 -> 138 -> L
    feats[1020] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-56') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-56', '138') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '138', node2.lemma)
    
    # L -> -56 -> du_fp -> L
    feats[1021] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-56') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-56', 'du_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_fp', node2.lemma)
    
    # L -> -83 -> -109 -> L
    feats[1022] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-83') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-83', '-109') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-109', node2.lemma)
    
    # L -> -83 -> -17 -> L
    feats[1023] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-83') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-83', '-17') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-17', node2.lemma)
    
    # L -> -83 -> -18 -> L
    feats[1024] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-83') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-83', '-18') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-18', node2.lemma)
    
    # L -> -83 -> -139 -> L
    feats[1025] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-83') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-83', '-139') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-139', node2.lemma)
    
    # L -> -83 -> 58 -> L
    feats[1026] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-83') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-83', '58') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '58', node2.lemma)
    
    # L -> -83 -> 2_du -> L
    feats[1027] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-83') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-83', '2_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_du', node2.lemma)
    
    # L -> -83 -> -154 -> L
    feats[1028] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-83') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-83', '-154') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-154', node2.lemma)
    
    # L -> -83 -> 181 -> L
    feats[1029] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-83') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-83', '181') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '181', node2.lemma)
    
    # L -> -83 -> -102 -> L
    feats[1030] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-83') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-83', '-102') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-102', node2.lemma)
    
    # L -> -83 -> -121 -> L
    feats[1031] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-83') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-83', '-121') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-121', node2.lemma)
    
    # L -> -83 -> instr. fem -> L
    feats[1032] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-83') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-83', 'instr. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. fem', node2.lemma)
    
    # L -> -83 -> -153 -> L
    feats[1033] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-83') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-83', '-153') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-153', node2.lemma)
    
    # L -> -83 -> acc. pl. -> L
    feats[1034] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-83') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-83', 'acc. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. pl.', node2.lemma)
    
    # L -> -83 -> -240 -> L
    feats[1035] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-83') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-83', '-240') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-240', node2.lemma)
    
    # L -> -83 -> -62 -> L
    feats[1036] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-83') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-83', '-62') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-62', node2.lemma)
    
    # L -> -83 -> -103 -> L
    feats[1037] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-83') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-83', '-103') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-103', node2.lemma)
    
    # L -> -83 -> -54 -> L
    feats[1038] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-83') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-83', '-54') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-54', node2.lemma)
    
    # L -> 5_sp -> dat. du. -> L
    feats[1039] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', 'dat. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    
    # L -> 5_sp -> 77 -> L
    feats[1040] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', '77') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '77', node2.lemma)
    
    # L -> 5_sp -> -273 -> L
    feats[1041] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', '-273') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-273', node2.lemma)
    
    # L -> 5_sp -> voc. pl. -> L
    feats[1042] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', 'voc. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. pl.', node2.lemma)
    
    # L -> 5_sp -> -114 -> L
    feats[1043] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', '-114') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-114', node2.lemma)
    
    # L -> 5_sp -> 16_fp -> L
    feats[1044] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', '16_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_fp', node2.lemma)
    
    # L -> 5_sp -> 115 -> L
    feats[1045] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', '115') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '115', node2.lemma)
    
    # L -> 5_sp -> 16_pl -> L
    feats[1046] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', '16_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_pl', node2.lemma)
    
    # L -> 132 -> 73 -> L
    feats[1047] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '132') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '132', '73') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '73', node2.lemma)
    
    # L -> nom. du. -> 110 -> L
    feats[1048] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. du.', '110') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '110', node2.lemma)
    
    # L -> nom. du. -> 151 -> L
    feats[1049] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. du.', '151') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '151', node2.lemma)
    
    # L -> 15_tp -> 12_fp -> L
    feats[1050] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_tp', '12_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_fp', node2.lemma)
    
    # L -> 15_tp -> 3_sg -> L
    feats[1051] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_tp', '3_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_sg', node2.lemma)
    
    # L -> 15_tp -> -87 -> L
    feats[1052] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_tp', '-87') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-87', node2.lemma)
    
    # L -> 15_tp -> 16_fp -> L
    feats[1053] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_tp', '16_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_fp', node2.lemma)
    
    # L -> 15_tp -> 172 -> L
    feats[1054] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_tp', '172') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '172', node2.lemma)
    
    # L -> 15_tp -> 4_tp -> L
    feats[1055] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_tp', '4_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_tp', node2.lemma)
    
    # L -> 15_tp -> 81 -> L
    feats[1056] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_tp', '81') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '81', node2.lemma)
    
    # L -> 15_tp -> 5_pl -> L
    feats[1057] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_tp', '5_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_pl', node2.lemma)
    
    # L -> 15_tp -> -13 -> L
    feats[1058] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_tp', '-13') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-13', node2.lemma)
    
    # L -> 15_tp -> -64 -> L
    feats[1059] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_tp', '-64') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-64', node2.lemma)
    
    # L -> 15_tp -> instr -> L
    feats[1060] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_tp', 'instr') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr', node2.lemma)
    
    # L -> 15_tp -> -169 -> L
    feats[1061] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_tp', '-169') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-169', node2.lemma)
    
    # L -> 15_tp -> 6_fp -> L
    feats[1062] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_tp', '6_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_fp', node2.lemma)
    
    # L -> 38 -> 60 -> L
    feats[1063] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '38') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '38', '60') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '60', node2.lemma)
    
    # L -> 38 -> 5_sg -> L
    feats[1064] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '38') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '38', '5_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sg', node2.lemma)
    
    # L -> 38 -> -12 -> L
    feats[1065] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '38') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '38', '-12') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-12', node2.lemma)
    
    # L -> 38 -> 55 -> L
    feats[1066] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '38') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '38', '55') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '55', node2.lemma)
    
    # L -> 38 -> 112 -> L
    feats[1067] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '38') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '38', '112') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '112', node2.lemma)
    
    # L -> 38 -> 151 -> L
    feats[1068] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '38') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '38', '151') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '151', node2.lemma)
    
    # L -> 38 -> 8_tp -> L
    feats[1069] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '38') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '38', '8_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_tp', node2.lemma)
    
    # L -> 38 -> -30 -> L
    feats[1070] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '38') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '38', '-30') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-30', node2.lemma)
    
    # L -> 38 -> 30_fp -> L
    feats[1071] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '38') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '38', '30_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_fp', node2.lemma)
    
    # L -> 38 -> 5_tp -> L
    feats[1072] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '38') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '38', '5_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_tp', node2.lemma)
    
    # L -> 38 -> 111 -> L
    feats[1073] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '38') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '38', '111') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '111', node2.lemma)
    
    # L -> 38 -> 29_sg -> L
    feats[1074] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '38') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '38', '29_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '29_sg', node2.lemma)
    
    # L -> 38 -> 51 -> L
    feats[1075] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '38') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '38', '51') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '51', node2.lemma)
    
    # L -> 38 -> -91 -> L
    feats[1076] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '38') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '38', '-91') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-91', node2.lemma)
    
    # L -> 38 -> 97 -> L
    feats[1077] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '38') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '38', '97') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '97', node2.lemma)
    
    # L -> 38 -> -114 -> L
    feats[1078] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '38') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '38', '-114') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-114', node2.lemma)
    
    # L -> 38 -> 158 -> L
    feats[1079] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '38') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '38', '158') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '158', node2.lemma)
    
    # L -> 38 -> -27 -> L
    feats[1080] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '38') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '38', '-27') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-27', node2.lemma)
    
    # L -> 38 -> 115 -> L
    feats[1081] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '38') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '38', '115') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '115', node2.lemma)
    
    # L -> 38 -> -245 -> L
    feats[1082] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '38') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '38', '-245') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-245', node2.lemma)
    
    # L -> 38 -> -126 -> L
    feats[1083] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '38') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '38', '-126') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-126', node2.lemma)
    
    # L -> 38 -> 180 -> L
    feats[1084] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '38') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '38', '180') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '180', node2.lemma)
    
    # L -> 38 -> masc -> L
    feats[1085] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '38') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '38', 'masc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'masc', node2.lemma)
    
    # L -> 38 -> 13_tp -> L
    feats[1086] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '38') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '38', '13_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_tp', node2.lemma)
    
    # L -> 38 -> 27_sg -> L
    feats[1087] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '38') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '38', '27_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_sg', node2.lemma)
    
    # L -> 38 -> 179 -> L
    feats[1088] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '38') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '38', '179') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '179', node2.lemma)
    
    # L -> sg_sp -> -276 -> L
    feats[1089] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg_sp', '-276') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-276', node2.lemma)
    
    # L -> sg_sp -> 70 -> L
    feats[1090] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg_sp', '70') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '70', node2.lemma)
    
    # L -> sg_sp -> -37 -> L
    feats[1091] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg_sp', '-37') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-37', node2.lemma)
    
    # L -> sg_sp -> -94 -> L
    feats[1092] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg_sp', '-94') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-94', node2.lemma)
    
    # L -> sg_sp -> 27_fp -> L
    feats[1093] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg_sp', '27_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_fp', node2.lemma)
    
    # L -> sg_sp -> 153 -> L
    feats[1094] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg_sp', '153') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    
    # L -> sg_sp -> fp -> L
    feats[1095] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg_sp', 'fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'fp', node2.lemma)
    
    # L -> sg_sp -> 11_sg -> L
    feats[1096] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg_sp', '11_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_sg', node2.lemma)
    
    # L -> sg_sp -> 156 -> L
    feats[1097] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg_sp', '156') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '156', node2.lemma)
    
    # L -> sg_sp -> -71 -> L
    feats[1098] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg_sp', '-71') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-71', node2.lemma)
    
    # L -> sg_sp -> 28_sg -> L
    feats[1099] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg_sp', '28_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28_sg', node2.lemma)
    
    # L -> sg_sp -> loc. sg. -> L
    feats[1100] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg_sp', 'loc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc. sg.', node2.lemma)
    
    # L -> sg_sp -> 8_pl -> L
    feats[1101] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg_sp', '8_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_pl', node2.lemma)
    
    # L -> sg_sp -> 8_sp -> L
    feats[1102] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg_sp', '8_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_sp', node2.lemma)
    
    # L -> sg_sp -> abl -> L
    feats[1103] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg_sp', 'abl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl', node2.lemma)
    
    # L -> sg_sp -> 171 -> L
    feats[1104] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg_sp', '171') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '171', node2.lemma)
    
    # L -> sg_sp -> -56 -> L
    feats[1105] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg_sp', '-56') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-56', node2.lemma)
    
    # L -> sg_sp -> 3_pl -> L
    feats[1106] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg_sp', '3_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_pl', node2.lemma)
    
    # L -> sg_sp -> pl_tp -> L
    feats[1107] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg_sp', 'pl_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_tp', node2.lemma)
    
    # L -> sg_sp -> 90 -> L
    feats[1108] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg_sp', '90') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '90', node2.lemma)
    
    # L -> 73 -> -291 -> L
    feats[1109] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '73') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '73', '-291') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-291', node2.lemma)
    
    # L -> 73 -> sg_sp -> L
    feats[1110] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '73') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '73', 'sg_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_sp', node2.lemma)
    
    # L -> 73 -> acc. fem -> L
    feats[1111] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '73') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '73', 'acc. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. fem', node2.lemma)
    
    # L -> 73 -> 2 -> L
    feats[1112] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '73') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '73', '2') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2', node2.lemma)
    
    # L -> acc. masc. -> -67 -> L
    feats[1113] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. masc.', '-67') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-67', node2.lemma)
    
    # L -> acc. masc. -> 79 -> L
    feats[1114] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. masc.', '79') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '79', node2.lemma)
    
    # L -> acc. masc. -> -163 -> L
    feats[1115] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. masc.', '-163') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-163', node2.lemma)
    
    # L -> acc. masc. -> masc -> L
    feats[1116] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. masc.', 'masc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'masc', node2.lemma)
    
    # L -> acc. masc. -> 14_pl -> L
    feats[1117] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. masc.', '14_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_pl', node2.lemma)
    
    # L -> acc. masc. -> -31 -> L
    feats[1118] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. masc.', '-31') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-31', node2.lemma)
    
    # L -> acc. masc. -> 114 -> L
    feats[1119] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. masc.', '114') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '114', node2.lemma)
    
    # L -> acc. masc. -> acc. neutr. -> L
    feats[1120] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. masc.', 'acc. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. neutr.', node2.lemma)
    
    # L -> acc. masc. -> -73 -> L
    feats[1121] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. masc.', '-73') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-73', node2.lemma)
    
    # L -> acc. masc. -> 8_fp -> L
    feats[1122] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. masc.', '8_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    
    # L -> acc. masc. -> -142 -> L
    feats[1123] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. masc.', '-142') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-142', node2.lemma)
    
    # L -> acc. masc. -> 178 -> L
    feats[1124] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. masc.', '178') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '178', node2.lemma)
    
    # L -> acc. masc. -> -115 -> L
    feats[1125] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. masc.', '-115') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-115', node2.lemma)
    
    # L -> acc. masc. -> 4_fp -> L
    feats[1126] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. masc.', '4_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_fp', node2.lemma)
    
    # L -> acc. masc. -> 117 -> L
    feats[1127] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. masc.', '117') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '117', node2.lemma)
    
    # L -> acc. masc. -> 68 -> L
    feats[1128] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. masc.', '68') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '68', node2.lemma)
    
    # L -> acc. masc. -> du -> L
    feats[1129] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. masc.', 'du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du', node2.lemma)
    
    # L -> acc. masc. -> -22 -> L
    feats[1130] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. masc.', '-22') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-22', node2.lemma)
    
    # L -> acc. masc. -> instr -> L
    feats[1131] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. masc.', 'instr') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr', node2.lemma)
    
    # L -> acc. masc. -> -169 -> L
    feats[1132] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. masc.', '-169') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-169', node2.lemma)
    
    # L -> acc. masc. -> 3_du -> L
    feats[1133] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. masc.', '3_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    
    # L -> acc. masc. -> 7_fp -> L
    feats[1134] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. masc.', '7_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_fp', node2.lemma)
    
    # L -> acc. masc. -> 61 -> L
    feats[1135] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. masc.', '61') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '61', node2.lemma)
    
    # L -> acc. masc. -> voc. du. -> L
    feats[1136] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. masc.', 'voc. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. du.', node2.lemma)
    
    # L -> 58 -> -301 -> L
    feats[1137] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '58', '-301') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-301', node2.lemma)
    
    # L -> 58 -> sg_tp -> L
    feats[1138] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '58', 'sg_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_tp', node2.lemma)
    
    # L -> 58 -> 155 -> L
    feats[1139] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '58', '155') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '155', node2.lemma)
    
    # L -> 58 -> gen. pl. -> L
    feats[1140] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '58', 'gen. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen. pl.', node2.lemma)
    
    # L -> 58 -> -46 -> L
    feats[1141] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '58', '-46') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-46', node2.lemma)
    
    # L -> 58 -> 10_sp -> L
    feats[1142] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '58', '10_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_sp', node2.lemma)
    
    # L -> 58 -> 31 -> L
    feats[1143] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '58', '31') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '31', node2.lemma)
    
    # L -> 58 -> -242 -> L
    feats[1144] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '58', '-242') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-242', node2.lemma)
    
    # L -> 58 -> 60 -> L
    feats[1145] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '58', '60') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '60', node2.lemma)
    
    # L -> 58 -> -291 -> L
    feats[1146] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '58', '-291') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-291', node2.lemma)
    
    # L -> 58 -> -14 -> L
    feats[1147] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '58', '-14') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-14', node2.lemma)
    
    # L -> 58 -> -109 -> L
    feats[1148] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '58', '-109') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-109', node2.lemma)
    
    # L -> 58 -> 2_sg -> L
    feats[1149] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '58', '2_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_sg', node2.lemma)
    
    # L -> 58 -> -247 -> L
    feats[1150] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '58', '-247') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-247', node2.lemma)
    
    # L -> 58 -> 35 -> L
    feats[1151] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '58', '35') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '35', node2.lemma)
    
    # L -> 58 -> -47 -> L
    feats[1152] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '58', '-47') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-47', node2.lemma)
    
    # L -> 58 -> -59 -> L
    feats[1153] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '58', '-59') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-59', node2.lemma)
    
    # L -> 58 -> -271 -> L
    feats[1154] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '58', '-271') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-271', node2.lemma)
    
    # L -> 58 -> nom. pl. -> L
    feats[1155] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '58', 'nom. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. pl.', node2.lemma)
    
    # L -> 58 -> 5_sg -> L
    feats[1156] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '58', '5_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sg', node2.lemma)
    
    # L -> 58 -> acc -> L
    feats[1157] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '58', 'acc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc', node2.lemma)
    
    # L -> 58 -> -141 -> L
    feats[1158] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '58', '-141') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-141', node2.lemma)
    
    # L -> 58 -> 170 -> L
    feats[1159] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '58', '170') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '170', node2.lemma)
    
    # L -> 58 -> 30_sg -> L
    feats[1160] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '58', '30_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_sg', node2.lemma)
    
    # L -> 58 -> 15_sg -> L
    feats[1161] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '58', '15_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_sg', node2.lemma)
    
    # L -> 58 -> loc. sg. -> L
    feats[1162] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '58', 'loc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc. sg.', node2.lemma)
    
    # L -> 58 -> -20 -> L
    feats[1163] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '58', '-20') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-20', node2.lemma)
    
    # L -> 58 -> -43 -> L
    feats[1164] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '58', '-43') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-43', node2.lemma)
    
    # L -> 58 -> 11_pl -> L
    feats[1165] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '58', '11_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_pl', node2.lemma)
    
    # L -> 58 -> 157 -> L
    feats[1166] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '58', '157') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '157', node2.lemma)
    
    # L -> 58 -> 15_sp -> L
    feats[1167] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '58', '15_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_sp', node2.lemma)
    
    # L -> 58 -> 128 -> L
    feats[1168] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '58', '128') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '128', node2.lemma)
    
    # L -> 58 -> 16_du -> L
    feats[1169] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '58', '16_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_du', node2.lemma)
    
    # L -> 58 -> pl_sp -> L
    feats[1170] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '58', 'pl_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_sp', node2.lemma)
    
    # L -> 58 -> -17 -> L
    feats[1171] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '58', '-17') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-17', node2.lemma)
    
    # L -> 58 -> 151 -> L
    feats[1172] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '58', '151') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '151', node2.lemma)
    
    # L -> 58 -> 10_sg -> L
    feats[1173] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '58', '10_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_sg', node2.lemma)
    
    # L -> 58 -> 16_sg -> L
    feats[1174] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '58', '16_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_sg', node2.lemma)
    
    # L -> 58 -> -98 -> L
    feats[1175] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '58', '-98') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-98', node2.lemma)
    
    # L -> 58 -> 8_tp -> L
    feats[1176] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '58', '8_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_tp', node2.lemma)
    
    # L -> -131 -> -271 -> L
    feats[1177] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-131') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-131', '-271') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-271', node2.lemma)
    
    # L -> -131 -> 8_tp -> L
    feats[1178] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-131') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-131', '8_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_tp', node2.lemma)
    
    # L -> -131 -> -260 -> L
    feats[1179] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-131') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-131', '-260') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-260', node2.lemma)
    
    # L -> -131 -> dat. du. -> L
    feats[1180] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-131') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-131', 'dat. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    
    # L -> -131 -> -240 -> L
    feats[1181] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-131') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-131', '-240') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-240', node2.lemma)
    
    # L -> -131 -> fem -> L
    feats[1182] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-131') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-131', 'fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'fem', node2.lemma)
    
    # L -> 39 -> nom. sg. -> L
    feats[1183] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', 'nom. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. sg.', node2.lemma)
    
    # L -> 39 -> 89 -> L
    feats[1184] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '89') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '89', node2.lemma)
    
    # L -> 39 -> 6_tp -> L
    feats[1185] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '6_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_tp', node2.lemma)
    
    # L -> 39 -> 7_tp -> L
    feats[1186] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '7_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_tp', node2.lemma)
    
    # L -> 39 -> 7_fp -> L
    feats[1187] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '7_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_fp', node2.lemma)
    
    # L -> 39 -> -86 -> L
    feats[1188] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '-86') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-86', node2.lemma)
    
    # L -> 34 -> nom. neutr. -> L
    feats[1189] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '34', 'nom. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. neutr.', node2.lemma)
    
    # L -> 161 -> 5_pl -> L
    feats[1190] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '161') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '161', '5_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_pl', node2.lemma)
    
    # L -> 102 -> 70 -> L
    feats[1191] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '102') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '102', '70') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '70', node2.lemma)
    
    # L -> 102 -> 13_sg -> L
    feats[1192] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '102') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '102', '13_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_sg', node2.lemma)
    
    # L -> 102 -> -149 -> L
    feats[1193] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '102') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '102', '-149') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-149', node2.lemma)
    
    # L -> 102 -> -17 -> L
    feats[1194] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '102') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '102', '-17') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-17', node2.lemma)
    
    # L -> 102 -> 8_tp -> L
    feats[1195] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '102') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '102', '8_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_tp', node2.lemma)
    
    # L -> du_tp -> 16_pl -> L
    feats[1196] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', '16_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_pl', node2.lemma)
    
    # L -> du_tp -> -58 -> L
    feats[1197] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', '-58') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-58', node2.lemma)
    
    # L -> instr. neutr. -> -151 -> L
    feats[1198] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. neutr.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. neutr.', '-151') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-151', node2.lemma)
    
    # L -> instr. neutr. -> 153 -> L
    feats[1199] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. neutr.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. neutr.', '153') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    
    # L -> instr. neutr. -> -59 -> L
    feats[1200] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. neutr.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. neutr.', '-59') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-59', node2.lemma)
    
    # L -> instr. neutr. -> 156 -> L
    feats[1201] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. neutr.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. neutr.', '156') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '156', node2.lemma)
    
    # L -> -89 -> -21 -> L
    feats[1202] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-89') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-89', '-21') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-21', node2.lemma)
    
    # L -> -89 -> 12_pl -> L
    feats[1203] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-89') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-89', '12_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_pl', node2.lemma)
    
    # L -> -89 -> 148 -> L
    feats[1204] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-89') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-89', '148') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '148', node2.lemma)
    
    # L -> -89 -> nom. sg. -> L
    feats[1205] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-89') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-89', 'nom. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. sg.', node2.lemma)
    
    # L -> -89 -> -147 -> L
    feats[1206] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-89') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-89', '-147') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-147', node2.lemma)
    
    # L -> -89 -> -144 -> L
    feats[1207] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-89') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-89', '-144') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-144', node2.lemma)
    
    # L -> -89 -> -69 -> L
    feats[1208] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-89') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-89', '-69') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-69', node2.lemma)
    
    # L -> -89 -> 172 -> L
    feats[1209] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-89') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-89', '172') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '172', node2.lemma)
    
    # L -> -97 -> -76 -> L
    feats[1210] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-97') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-97', '-76') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-76', node2.lemma)
    
    # L -> 14_tp -> 29 -> L
    feats[1211] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_tp', '29') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '29', node2.lemma)
    
    # L -> 14_tp -> 14_pl -> L
    feats[1212] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_tp', '14_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_pl', node2.lemma)
    
    # L -> 14_tp -> 2 -> L
    feats[1213] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_tp', '2') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2', node2.lemma)
    
    # L -> 14_tp -> -64 -> L
    feats[1214] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_tp', '-64') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-64', node2.lemma)
    
    # L -> 14_tp -> -24 -> L
    feats[1215] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_tp', '-24') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    
    # L -> 14_tp -> 137 -> L
    feats[1216] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_tp', '137') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    
    # L -> 14_tp -> -132 -> L
    feats[1217] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_tp', '-132') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-132', node2.lemma)
    
    # L -> 14_tp -> -292 -> L
    feats[1218] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_tp', '-292') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-292', node2.lemma)
    
    # L -> 14_tp -> du_fp -> L
    feats[1219] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_tp', 'du_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_fp', node2.lemma)
    
    # L -> 14_tp -> 7_fp -> L
    feats[1220] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_tp', '7_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_fp', node2.lemma)
    
    # L -> pl -> -283 -> L
    feats[1221] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl', '-283') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-283', node2.lemma)
    
    # L -> pl -> -79 -> L
    feats[1222] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl', '-79') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-79', node2.lemma)
    
    # L -> pl -> instr. pl. -> L
    feats[1223] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl', 'instr. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. pl.', node2.lemma)
    
    # L -> pl -> -247 -> L
    feats[1224] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl', '-247') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-247', node2.lemma)
    
    # L -> pl -> -47 -> L
    feats[1225] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl', '-47') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-47', node2.lemma)
    
    # L -> pl -> 153 -> L
    feats[1226] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl', '153') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    
    # L -> pl -> -59 -> L
    feats[1227] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl', '-59') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-59', node2.lemma)
    
    # L -> pl -> -12 -> L
    feats[1228] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl', '-12') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-12', node2.lemma)
    
    # L -> pl -> -241 -> L
    feats[1229] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl', '-241') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-241', node2.lemma)
    
    # L -> pl -> -293 -> L
    feats[1230] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl', '-293') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-293', node2.lemma)
    
    # L -> pl -> 6_sg -> L
    feats[1231] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl', '6_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_sg', node2.lemma)
    
    # L -> pl -> -49 -> L
    feats[1232] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl', '-49') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-49', node2.lemma)
    
    # L -> pl -> -84 -> L
    feats[1233] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl', '-84') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-84', node2.lemma)
    
    # L -> pl -> -102 -> L
    feats[1234] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl', '-102') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-102', node2.lemma)
    
    # L -> pl -> 152 -> L
    feats[1235] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl', '152') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '152', node2.lemma)
    
    # L -> voc -> fem -> L
    feats[1236] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc', 'fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'fem', node2.lemma)
    
    # L -> voc -> -230 -> L
    feats[1237] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc', '-230') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-230', node2.lemma)
    
    # L -> voc -> 175 -> L
    feats[1238] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc', '175') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '175', node2.lemma)
    
    # L -> voc -> 29_tp -> L
    feats[1239] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc', '29_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '29_tp', node2.lemma)
    
    # L -> voc -> nom -> L
    feats[1240] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc', 'nom') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom', node2.lemma)
    
    # L -> voc -> 150 -> L
    feats[1241] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc', '150') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '150', node2.lemma)
    
    # L -> voc -> 160 -> L
    feats[1242] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc', '160') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    
    # L -> voc -> 6_fp -> L
    feats[1243] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc', '6_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_fp', node2.lemma)
    
    # L -> -29 -> -276 -> L
    feats[1244] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-29') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-29', '-276') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-276', node2.lemma)
    
    # L -> -29 -> -26 -> L
    feats[1245] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-29') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-29', '-26') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-26', node2.lemma)
    
    # L -> -66 -> 10_sg -> L
    feats[1246] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-66') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-66', '10_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_sg', node2.lemma)
    
    # L -> -66 -> 12_sg -> L
    feats[1247] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-66') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-66', '12_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_sg', node2.lemma)
    
    # L -> -66 -> 7_sp -> L
    feats[1248] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-66') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-66', '7_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_sp', node2.lemma)
    
    # L -> -66 -> 140 -> L
    feats[1249] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-66') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-66', '140') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '140', node2.lemma)
    
    # L -> -66 -> voc. neutr. -> L
    feats[1250] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-66') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-66', 'voc. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. neutr.', node2.lemma)
    
    # L -> -66 -> 90 -> L
    feats[1251] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-66') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-66', '90') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '90', node2.lemma)
    
    # L -> -66 -> 3_sp -> L
    feats[1252] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-66') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-66', '3_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_sp', node2.lemma)
    
    # L -> -66 -> 174 -> L
    feats[1253] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-66') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-66', '174') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '174', node2.lemma)
    
    # L -> -66 -> -19 -> L
    feats[1254] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-66') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-66', '-19') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-19', node2.lemma)
    
    # L -> -66 -> acc. fem -> L
    feats[1255] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-66') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-66', 'acc. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. fem', node2.lemma)
    
    # L -> -66 -> -102 -> L
    feats[1256] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-66') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-66', '-102') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-102', node2.lemma)
    
    # L -> -66 -> instr. du. -> L
    feats[1257] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-66') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-66', 'instr. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. du.', node2.lemma)
    
    # L -> -66 -> 152 -> L
    feats[1258] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-66') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-66', '152') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '152', node2.lemma)
    
    # L -> -66 -> -249 -> L
    feats[1259] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-66') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-66', '-249') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-249', node2.lemma)
    
    # L -> -66 -> 40 -> L
    feats[1260] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-66') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-66', '40') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '40', node2.lemma)
    
    # L -> 119 -> 157 -> L
    feats[1261] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '119') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '119', '157') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '157', node2.lemma)
    
    # L -> 119 -> 15_sp -> L
    feats[1262] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '119') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '119', '15_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_sp', node2.lemma)
    
    # L -> 119 -> -18 -> L
    feats[1263] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '119') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '119', '-18') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-18', node2.lemma)
    
    # L -> 119 -> acc. masc. -> L
    feats[1264] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '119') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '119', 'acc. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. masc.', node2.lemma)
    
    # L -> 119 -> 2_du -> L
    feats[1265] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '119') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '119', '2_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_du', node2.lemma)
    
    # L -> 119 -> 4_du -> L
    feats[1266] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '119') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '119', '4_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_du', node2.lemma)
    
    # L -> 119 -> 1 -> L
    feats[1267] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '119') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '119', '1') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '1', node2.lemma)
    
    # L -> 119 -> abl. du. -> L
    feats[1268] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '119') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '119', 'abl. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. du.', node2.lemma)
    
    # L -> 119 -> -97 -> L
    feats[1269] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '119') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '119', '-97') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-97', node2.lemma)
    
    # L -> 119 -> 14_tp -> L
    feats[1270] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '119') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '119', '14_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_tp', node2.lemma)
    
    # L -> 119 -> pl -> L
    feats[1271] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '119') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '119', 'pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    
    # L -> 119 -> 108 -> L
    feats[1272] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '119') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '119', '108') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '108', node2.lemma)
    
    # L -> 119 -> dat. du. -> L
    feats[1273] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '119') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '119', 'dat. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    
    # L -> 119 -> 14_sg -> L
    feats[1274] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '119') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '119', '14_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_sg', node2.lemma)
    
    # L -> 119 -> 91 -> L
    feats[1275] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '119') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '119', '91') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    
    # L -> -63 -> -78 -> L
    feats[1276] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-63') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-63', '-78') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-78', node2.lemma)
    
    # L -> -63 -> 15_pl -> L
    feats[1277] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-63') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-63', '15_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_pl', node2.lemma)
    
    # L -> -63 -> -55 -> L
    feats[1278] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-63') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-63', '-55') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-55', node2.lemma)
    
    # L -> -63 -> gen. sg. -> L
    feats[1279] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-63') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-63', 'gen. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen. sg.', node2.lemma)
    
    # L -> -63 -> 14_pl -> L
    feats[1280] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-63') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-63', '14_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_pl', node2.lemma)
    
    # L -> -63 -> 68 -> L
    feats[1281] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-63') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-63', '68') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '68', node2.lemma)
    
    # L -> -61 -> 94 -> L
    feats[1282] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-61') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-61', '94') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '94', node2.lemma)
    
    # L -> -61 -> 12_sp -> L
    feats[1283] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-61') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-61', '12_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_sp', node2.lemma)
    
    # L -> -61 -> 16_du -> L
    feats[1284] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-61') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-61', '16_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_du', node2.lemma)
    
    # L -> -61 -> -268 -> L
    feats[1285] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-61') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-61', '-268') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-268', node2.lemma)
    
    # L -> -61 -> 37 -> L
    feats[1286] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-61') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-61', '37') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '37', node2.lemma)
    
    # L -> -61 -> nom. du. -> L
    feats[1287] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-61') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-61', 'nom. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. du.', node2.lemma)
    
    # L -> -61 -> 15_tp -> L
    feats[1288] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-61') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-61', '15_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_tp', node2.lemma)
    
    # L -> -61 -> gen -> L
    feats[1289] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-61') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-61', 'gen') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen', node2.lemma)
    
    # L -> -61 -> 7_sp -> L
    feats[1290] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-61') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-61', '7_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_sp', node2.lemma)
    
    # L -> 30 -> -28 -> L
    feats[1291] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30', '-28') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-28', node2.lemma)
    
    # L -> 30 -> acc. pl. -> L
    feats[1292] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30', 'acc. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. pl.', node2.lemma)
    
    # L -> 30 -> fem -> L
    feats[1293] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30', 'fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'fem', node2.lemma)
    
    # L -> 13_pl -> -276 -> L
    feats[1294] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '-276') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-276', node2.lemma)
    
    # L -> 13_pl -> -25 -> L
    feats[1295] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '-25') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-25', node2.lemma)
    
    # L -> 13_pl -> -156 -> L
    feats[1296] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '-156') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-156', node2.lemma)
    
    # L -> -309 -> 169 -> L
    feats[1297] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-309') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-309', '169') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '169', node2.lemma)
    
    # L -> -309 -> -76 -> L
    feats[1298] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-309') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-309', '-76') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-76', node2.lemma)
    
    # L -> -309 -> 33 -> L
    feats[1299] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-309') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-309', '33') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '33', node2.lemma)
    
    # L -> -309 -> -126 -> L
    feats[1300] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-309') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-309', '-126') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-126', node2.lemma)
    
    # L -> -309 -> -119 -> L
    feats[1301] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-309') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-309', '-119') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-119', node2.lemma)
    
    # L -> -309 -> 7_tp -> L
    feats[1302] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-309') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-309', '7_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_tp', node2.lemma)
    
    # L -> -309 -> -112 -> L
    feats[1303] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-309') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-309', '-112') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    
    # L -> -309 -> -132 -> L
    feats[1304] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-309') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-309', '-132') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-132', node2.lemma)
    
    # L -> -296 -> -301 -> L
    feats[1305] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-296') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-296', '-301') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-301', node2.lemma)
    
    # L -> -296 -> sg_tp -> L
    feats[1306] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-296') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-296', 'sg_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_tp', node2.lemma)
    
    # L -> -296 -> dat -> L
    feats[1307] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-296') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-296', 'dat') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat', node2.lemma)
    
    # L -> -296 -> 11_fp -> L
    feats[1308] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-296') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-296', '11_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_fp', node2.lemma)
    
    # L -> -296 -> 31 -> L
    feats[1309] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-296') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-296', '31') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '31', node2.lemma)
    
    # L -> sg_fp -> 28_sg -> L
    feats[1310] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg_fp', '28_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28_sg', node2.lemma)
    
    # L -> -67 -> -262 -> L
    feats[1311] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', '-262') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-262', node2.lemma)
    
    # L -> -67 -> 88 -> L
    feats[1312] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', '88') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '88', node2.lemma)
    
    # L -> -67 -> 15_sp -> L
    feats[1313] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', '15_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_sp', node2.lemma)
    
    # L -> -67 -> 128 -> L
    feats[1314] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', '128') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '128', node2.lemma)
    
    # L -> -67 -> 30_sp -> L
    feats[1315] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', '30_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_sp', node2.lemma)
    
    # L -> -67 -> -156 -> L
    feats[1316] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', '-156') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-156', node2.lemma)
    
    # L -> -67 -> -17 -> L
    feats[1317] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', '-17') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-17', node2.lemma)
    
    # L -> -67 -> -98 -> L
    feats[1318] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', '-98') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-98', node2.lemma)
    
    # L -> -67 -> 9_du -> L
    feats[1319] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', '9_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_du', node2.lemma)
    
    # L -> -67 -> 8_sp -> L
    feats[1320] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', '8_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_sp', node2.lemma)
    
    # L -> -67 -> -104 -> L
    feats[1321] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', '-104') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-104', node2.lemma)
    
    # L -> -67 -> abl -> L
    feats[1322] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', 'abl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl', node2.lemma)
    
    # L -> -67 -> 171 -> L
    feats[1323] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', '171') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '171', node2.lemma)
    
    # L -> -67 -> -83 -> L
    feats[1324] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', '-83') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-83', node2.lemma)
    
    # L -> -67 -> 15_tp -> L
    feats[1325] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', '15_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_tp', node2.lemma)
    
    # L -> -67 -> 38 -> L
    feats[1326] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', '38') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '38', node2.lemma)
    
    # L -> -67 -> sg_sp -> L
    feats[1327] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', 'sg_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_sp', node2.lemma)
    
    # L -> -67 -> acc. masc. -> L
    feats[1328] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', 'acc. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. masc.', node2.lemma)
    
    # L -> -67 -> 58 -> L
    feats[1329] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', '58') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '58', node2.lemma)
    
    # L -> -67 -> 12_sg -> L
    feats[1330] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', '12_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_sg', node2.lemma)
    
    # L -> -67 -> -154 -> L
    feats[1331] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', '-154') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-154', node2.lemma)
    
    # L -> -67 -> 7_sp -> L
    feats[1332] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', '7_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_sp', node2.lemma)
    
    # L -> -67 -> 140 -> L
    feats[1333] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', '140') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '140', node2.lemma)
    
    # L -> -67 -> voc. neutr. -> L
    feats[1334] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', 'voc. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. neutr.', node2.lemma)
    
    # L -> -67 -> 90 -> L
    feats[1335] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', '90') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '90', node2.lemma)
    
    # L -> -67 -> -161 -> L
    feats[1336] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', '-161') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-161', node2.lemma)
    
    # L -> -67 -> 142 -> L
    feats[1337] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', '142') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '142', node2.lemma)
    
    # L -> -67 -> 122 -> L
    feats[1338] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', '122') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '122', node2.lemma)
    
    # L -> -67 -> 8_sg -> L
    feats[1339] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', '8_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_sg', node2.lemma)
    
    # L -> -67 -> 12_fp -> L
    feats[1340] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', '12_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_fp', node2.lemma)
    
    # L -> -67 -> gen. du. -> L
    feats[1341] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', 'gen. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen. du.', node2.lemma)
    
    # L -> -67 -> 168 -> L
    feats[1342] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', '168') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '168', node2.lemma)
    
    # L -> -67 -> sg -> L
    feats[1343] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', 'sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg', node2.lemma)
    
    # L -> -67 -> 174 -> L
    feats[1344] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', '174') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '174', node2.lemma)
    
    # L -> -67 -> -133 -> L
    feats[1345] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', '-133') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-133', node2.lemma)
    
    # L -> -67 -> -190 -> L
    feats[1346] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', '-190') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-190', node2.lemma)
    
    # L -> -67 -> -44 -> L
    feats[1347] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', '-44') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-44', node2.lemma)
    
    # L -> -67 -> instr. sg. -> L
    feats[1348] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', 'instr. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. sg.', node2.lemma)
    
    # L -> -67 -> 39 -> L
    feats[1349] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', '39') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '39', node2.lemma)
    
    # L -> -67 -> 34 -> L
    feats[1350] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', '34') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '34', node2.lemma)
    
    # L -> -67 -> 161 -> L
    feats[1351] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', '161') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '161', node2.lemma)
    
    # L -> -67 -> voc -> L
    feats[1352] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', 'voc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc', node2.lemma)
    
    # L -> -67 -> -66 -> L
    feats[1353] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', '-66') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-66', node2.lemma)
    
    # L -> -67 -> 92 -> L
    feats[1354] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', '92') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '92', node2.lemma)
    
    # L -> -67 -> 13_pl -> L
    feats[1355] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', '13_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_pl', node2.lemma)
    
    # L -> -67 -> -309 -> L
    feats[1356] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', '-309') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-309', node2.lemma)
    
    # L -> -67 -> -296 -> L
    feats[1357] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', '-296') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-296', node2.lemma)
    
    # L -> -67 -> sg_fp -> L
    feats[1358] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', 'sg_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_fp', node2.lemma)
    
    # L -> -67 -> acc. fem -> L
    feats[1359] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', 'acc. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. fem', node2.lemma)
    
    # L -> -67 -> -28 -> L
    feats[1360] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', '-28') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-28', node2.lemma)
    
    # L -> -67 -> -121 -> L
    feats[1361] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', '-121') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-121', node2.lemma)
    
    # L -> -67 -> -307 -> L
    feats[1362] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', '-307') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-307', node2.lemma)
    
    # L -> -67 -> 182 -> L
    feats[1363] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', '182') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '182', node2.lemma)
    
    # L -> -67 -> sp -> L
    feats[1364] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', 'sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sp', node2.lemma)
    
    # L -> -67 -> -112 -> L
    feats[1365] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', '-112') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    
    # L -> 4_pl -> -47 -> L
    feats[1366] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_pl', '-47') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-47', node2.lemma)
    
    # L -> 4_pl -> -35 -> L
    feats[1367] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_pl', '-35') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    
    # L -> 4_pl -> -84 -> L
    feats[1368] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_pl', '-84') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-84', node2.lemma)
    
    # L -> 4_pl -> 4_pl -> L
    feats[1369] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_pl', '4_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_pl', node2.lemma)
    
    # L -> 4_pl -> -266 -> L
    feats[1370] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_pl', '-266') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-266', node2.lemma)
    
    # L -> -102 -> 11_du -> L
    feats[1371] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-102') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-102', '11_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_du', node2.lemma)
    
    # L -> -102 -> nom. neutr. -> L
    feats[1372] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-102') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-102', 'nom. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. neutr.', node2.lemma)
    
    # L -> -102 -> -276 -> L
    feats[1373] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-102') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-102', '-276') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-276', node2.lemma)
    
    # L -> -102 -> -37 -> L
    feats[1374] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-102') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-102', '-37') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-37', node2.lemma)
    
    # L -> -102 -> 72 -> L
    feats[1375] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-102') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-102', '72') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '72', node2.lemma)
    
    # L -> -102 -> 60 -> L
    feats[1376] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-102') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-102', '60') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '60', node2.lemma)
    
    # L -> -102 -> -291 -> L
    feats[1377] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-102') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-102', '-291') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-291', node2.lemma)
    
    # L -> -102 -> -81 -> L
    feats[1378] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-102') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-102', '-81') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-81', node2.lemma)
    
    # L -> -102 -> -247 -> L
    feats[1379] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-102') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-102', '-247') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-247', node2.lemma)
    
    # L -> -102 -> fp -> L
    feats[1380] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-102') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-102', 'fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'fp', node2.lemma)
    
    # L -> -102 -> 5_sg -> L
    feats[1381] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-102') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-102', '5_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sg', node2.lemma)
    
    # L -> -102 -> 30_sg -> L
    feats[1382] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-102') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-102', '30_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_sg', node2.lemma)
    
    # L -> -102 -> -150 -> L
    feats[1383] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-102') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-102', '-150') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-150', node2.lemma)
    
    # L -> -102 -> 15_sp -> L
    feats[1384] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-102') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-102', '15_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_sp', node2.lemma)
    
    # L -> -102 -> pl_sp -> L
    feats[1385] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-102') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-102', 'pl_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_sp', node2.lemma)
    
    # L -> -102 -> -18 -> L
    feats[1386] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-102') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-102', '-18') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-18', node2.lemma)
    
    # L -> -102 -> -139 -> L
    feats[1387] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-102') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-102', '-139') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-139', node2.lemma)
    
    # L -> -102 -> 112 -> L
    feats[1388] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-102') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-102', '112') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '112', node2.lemma)
    
    # L -> -102 -> 73 -> L
    feats[1389] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-102') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-102', '73') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '73', node2.lemma)
    
    # L -> -102 -> acc. masc. -> L
    feats[1390] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-102') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-102', 'acc. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. masc.', node2.lemma)
    
    # L -> -102 -> 58 -> L
    feats[1391] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-102') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-102', '58') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '58', node2.lemma)
    
    # L -> -102 -> -260 -> L
    feats[1392] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-102') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-102', '-260') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-260', node2.lemma)
    
    # L -> -102 -> 136 -> L
    feats[1393] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-102') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-102', '136') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '136', node2.lemma)
    
    # L -> -102 -> 4_pl -> L
    feats[1394] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-102') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-102', '4_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_pl', node2.lemma)
    
    # L -> -102 -> 74 -> L
    feats[1395] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-102') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-102', '74') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '74', node2.lemma)
    
    # L -> -102 -> -15 -> L
    feats[1396] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-102') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-102', '-15') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-15', node2.lemma)
    
    # L -> -102 -> -52 -> L
    feats[1397] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-102') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-102', '-52') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-52', node2.lemma)
    
    # L -> -102 -> 15_pl -> L
    feats[1398] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-102') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-102', '15_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_pl', node2.lemma)
    
    # L -> abl. pl. -> 7_du -> L
    feats[1399] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'abl. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'abl. pl.', '7_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_du', node2.lemma)
    
    # L -> abl. pl. -> voc. pl. -> L
    feats[1400] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'abl. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'abl. pl.', 'voc. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. pl.', node2.lemma)
    
    # L -> abl. pl. -> -34 -> L
    feats[1401] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'abl. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'abl. pl.', '-34') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-34', node2.lemma)
    
    # L -> abl. pl. -> 7_fp -> L
    feats[1402] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'abl. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'abl. pl.', '7_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_fp', node2.lemma)
    
    # L -> -28 -> -26 -> L
    feats[1403] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-28') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-28', '-26') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-26', node2.lemma)
    
    # L -> -28 -> 30_sg -> L
    feats[1404] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-28') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-28', '30_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_sg', node2.lemma)
    
    # L -> -28 -> 15_sp -> L
    feats[1405] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-28') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-28', '15_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_sp', node2.lemma)
    
    # L -> -28 -> -268 -> L
    feats[1406] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-28') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-28', '-268') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-268', node2.lemma)
    
    # L -> -28 -> 120 -> L
    feats[1407] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-28') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-28', '120') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '120', node2.lemma)
    
    # L -> -28 -> 15_du -> L
    feats[1408] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-28') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-28', '15_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_du', node2.lemma)
    
    # L -> -33 -> 91 -> L
    feats[1409] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '91') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    
    # L -> -33 -> -39 -> L
    feats[1410] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-39') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-39', node2.lemma)
    
    # L -> -33 -> -143 -> L
    feats[1411] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-143') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-143', node2.lemma)
    
    # L -> -33 -> -152 -> L
    feats[1412] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-152') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-152', node2.lemma)
    
    # L -> -33 -> 2_pl -> L
    feats[1413] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '2_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_pl', node2.lemma)
    
    # L -> -33 -> -153 -> L
    feats[1414] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-153') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-153', node2.lemma)
    
    # L -> -33 -> -11 -> L
    feats[1415] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-11') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-11', node2.lemma)
    
    # L -> -33 -> du_sp -> L
    feats[1416] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', 'du_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_sp', node2.lemma)
    
    # L -> -33 -> -30 -> L
    feats[1417] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-30') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-30', node2.lemma)
    
    # L -> -33 -> 30_fp -> L
    feats[1418] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '30_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_fp', node2.lemma)
    
    # L -> -33 -> 10_du -> L
    feats[1419] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '10_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_du', node2.lemma)
    
    # L -> -33 -> 4_sp -> L
    feats[1420] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '4_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_sp', node2.lemma)
    
    # L -> -33 -> instr. adj. -> L
    feats[1421] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', 'instr. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. adj.', node2.lemma)
    
    # L -> -33 -> -273 -> L
    feats[1422] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-273') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-273', node2.lemma)
    
    # L -> -33 -> 4_sg -> L
    feats[1423] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '4_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_sg', node2.lemma)
    
    # L -> -33 -> -82 -> L
    feats[1424] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-82') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-82', node2.lemma)
    
    # L -> -33 -> -240 -> L
    feats[1425] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-240') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-240', node2.lemma)
    
    # L -> -33 -> 16_tp -> L
    feats[1426] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '16_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_tp', node2.lemma)
    
    # L -> -33 -> 12_pl -> L
    feats[1427] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '12_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_pl', node2.lemma)
    
    # L -> -33 -> nom. masc. -> L
    feats[1428] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', 'nom. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. masc.', node2.lemma)
    
    # L -> -33 -> 138 -> L
    feats[1429] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '138') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '138', node2.lemma)
    
    # L -> -33 -> -15 -> L
    feats[1430] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-15') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-15', node2.lemma)
    
    # L -> -33 -> -243 -> L
    feats[1431] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-243') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-243', node2.lemma)
    
    # L -> -33 -> -91 -> L
    feats[1432] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-91') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-91', node2.lemma)
    
    # L -> -33 -> 169 -> L
    feats[1433] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '169') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '169', node2.lemma)
    
    # L -> -33 -> 14_sp -> L
    feats[1434] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '14_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_sp', node2.lemma)
    
    # L -> -33 -> -103 -> L
    feats[1435] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-103') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-103', node2.lemma)
    
    # L -> -33 -> 5_du -> L
    feats[1436] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '5_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_du', node2.lemma)
    
    # L -> -33 -> -114 -> L
    feats[1437] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-114') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-114', node2.lemma)
    
    # L -> -33 -> loc -> L
    feats[1438] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', 'loc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc', node2.lemma)
    
    # L -> -33 -> 158 -> L
    feats[1439] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '158') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '158', node2.lemma)
    
    # L -> -33 -> 11_tp -> L
    feats[1440] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '11_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_tp', node2.lemma)
    
    # L -> -33 -> -55 -> L
    feats[1441] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-55') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-55', node2.lemma)
    
    # L -> -33 -> nom. sg. -> L
    feats[1442] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', 'nom. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. sg.', node2.lemma)
    
    # L -> -33 -> -279 -> L
    feats[1443] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-279') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-279', node2.lemma)
    
    # L -> -33 -> 56 -> L
    feats[1444] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '56') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '56', node2.lemma)
    
    # L -> -33 -> 118 -> L
    feats[1445] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '118') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '118', node2.lemma)
    
    # L -> -33 -> 16_fp -> L
    feats[1446] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '16_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_fp', node2.lemma)
    
    # L -> -33 -> 115 -> L
    feats[1447] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '115') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '115', node2.lemma)
    
    # L -> -33 -> -68 -> L
    feats[1448] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-68') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-68', node2.lemma)
    
    # L -> -33 -> -69 -> L
    feats[1449] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-69') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-69', node2.lemma)
    
    # L -> -33 -> -42 -> L
    feats[1450] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-42') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-42', node2.lemma)
    
    # L -> -33 -> loc. pl. -> L
    feats[1451] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', 'loc. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc. pl.', node2.lemma)
    
    # L -> -33 -> 28_tp -> L
    feats[1452] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '28_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28_tp', node2.lemma)
    
    # L -> -33 -> 76 -> L
    feats[1453] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '76') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '76', node2.lemma)
    
    # L -> -33 -> gen. sg. -> L
    feats[1454] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', 'gen. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen. sg.', node2.lemma)
    
    # L -> -33 -> voc. masc. -> L
    feats[1455] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', 'voc. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. masc.', node2.lemma)
    
    # L -> -33 -> 180 -> L
    feats[1456] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '180') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '180', node2.lemma)
    
    # L -> -33 -> masc -> L
    feats[1457] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', 'masc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'masc', node2.lemma)
    
    # L -> -33 -> -166 -> L
    feats[1458] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-166') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-166', node2.lemma)
    
    # L -> -33 -> -119 -> L
    feats[1459] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '-119') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-119', node2.lemma)
    
    # L -> -33 -> 29 -> L
    feats[1460] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', '29') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '29', node2.lemma)
    
    # L -> 78 -> 75 -> L
    feats[1461] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '78') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '78', '75') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '75', node2.lemma)
    
    # L -> 78 -> -28 -> L
    feats[1462] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '78') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '78', '-28') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-28', node2.lemma)
    
    # L -> 78 -> 10_pl -> L
    feats[1463] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '78') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '78', '10_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_pl', node2.lemma)
    
    # L -> 78 -> 15_fp -> L
    feats[1464] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '78') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '78', '15_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_fp', node2.lemma)
    
    # L -> 78 -> -82 -> L
    feats[1465] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '78') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '78', '-82') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-82', node2.lemma)
    
    # L -> 78 -> acc. pl. -> L
    feats[1466] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '78') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '78', 'acc. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. pl.', node2.lemma)
    
    # L -> 78 -> -240 -> L
    feats[1467] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '78') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '78', '-240') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-240', node2.lemma)
    
    # L -> 78 -> fem -> L
    feats[1468] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '78') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '78', 'fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'fem', node2.lemma)
    
    # L -> 78 -> -62 -> L
    feats[1469] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '78') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '78', '-62') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-62', node2.lemma)
    
    # L -> 78 -> 14_sp -> L
    feats[1470] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '78') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '78', '14_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_sp', node2.lemma)
    
    # L -> 78 -> -158 -> L
    feats[1471] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '78') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '78', '-158') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-158', node2.lemma)
    
    # L -> 182 -> 10_sp -> L
    feats[1472] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '182') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '182', '10_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_sp', node2.lemma)
    
    # L -> 182 -> 8_du -> L
    feats[1473] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '182') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '182', '8_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_du', node2.lemma)
    
    # L -> 154 -> -142 -> L
    feats[1474] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '154') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '154', '-142') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-142', node2.lemma)
    
    # L -> 14_sg -> 54 -> L
    feats[1475] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_sg', '54') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    
    # L -> 14_sg -> 41 -> L
    feats[1476] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_sg', '41') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '41', node2.lemma)
    
    # L -> 14_sg -> -45 -> L
    feats[1477] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_sg', '-45') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-45', node2.lemma)
    
    # L -> 14_sg -> -29 -> L
    feats[1478] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_sg', '-29') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-29', node2.lemma)
    
    # L -> 14_sg -> -121 -> L
    feats[1479] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_sg', '-121') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-121', node2.lemma)
    
    # L -> 14_sg -> sp -> L
    feats[1480] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_sg', 'sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sp', node2.lemma)
    
    # L -> instr. masc. -> 15_sp -> L
    feats[1481] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. masc.', '15_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_sp', node2.lemma)
    
    # L -> instr. masc. -> 128 -> L
    feats[1482] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. masc.', '128') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '128', node2.lemma)
    
    # L -> instr. masc. -> -49 -> L
    feats[1483] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. masc.', '-49') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-49', node2.lemma)
    
    # L -> instr. masc. -> nom. du. -> L
    feats[1484] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. masc.', 'nom. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. du.', node2.lemma)
    
    # L -> instr. masc. -> 73 -> L
    feats[1485] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. masc.', '73') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '73', node2.lemma)
    
    # L -> instr. masc. -> 96 -> L
    feats[1486] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. masc.', '96') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '96', node2.lemma)
    
    # L -> instr. masc. -> 99 -> L
    feats[1487] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. masc.', '99') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '99', node2.lemma)
    
    # L -> instr. masc. -> abl. du. -> L
    feats[1488] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. masc.', 'abl. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. du.', node2.lemma)
    
    # L -> instr. masc. -> 174 -> L
    feats[1489] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. masc.', '174') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '174', node2.lemma)
    
    # L -> instr. masc. -> -19 -> L
    feats[1490] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. masc.', '-19') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-19', node2.lemma)
    
    # L -> instr. masc. -> -131 -> L
    feats[1491] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. masc.', '-131') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-131', node2.lemma)
    
    # L -> instr. masc. -> 39 -> L
    feats[1492] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. masc.', '39') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '39', node2.lemma)
    
    # L -> instr. masc. -> 34 -> L
    feats[1493] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. masc.', '34') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '34', node2.lemma)
    
    # L -> instr. masc. -> -309 -> L
    feats[1494] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. masc.', '-309') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-309', node2.lemma)
    
    # L -> instr. masc. -> -296 -> L
    feats[1495] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. masc.', '-296') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-296', node2.lemma)
    
    # L -> instr. masc. -> acc. fem -> L
    feats[1496] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. masc.', 'acc. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. fem', node2.lemma)
    
    # L -> 2_fp -> 115 -> L
    feats[1497] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', '115') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '115', node2.lemma)
    
    # L -> 2_fp -> 27_sg -> L
    feats[1498] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', '27_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_sg', node2.lemma)
    
    # L -> 2_fp -> -292 -> L
    feats[1499] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', '-292') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-292', node2.lemma)
    
    feats[feats < 1e-25] = 1e-25
    return -np.log10(feats)
    