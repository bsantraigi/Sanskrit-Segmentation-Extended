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
                    s = _neuralnet.Get_Energy(featVMat[i][j])
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
                    s = _neuralnet.Get_Energy(featVMat[i][j])
                    WScalarMat[i, j] = s   
        return WScalarMat



########################## ARTIFICIALLY GENERATED CODE FOR FEATURE GENERATION ###########################
############################# FOLLOWING FEATURES WERE SELECTED BASED ON A ###############################
############################# MUTUAL INFORMATION BASED SELECTION CRITERIA ###############################


def Get_Features(node1, node2):
    feats = np.zeros((_edge_vector_dim, 1))
    fIndex = 0
    
    # L->sg_sp->T
    feats[0] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg_sp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'sg_sp', node2.tup)
            
    # L->-266->T
    feats[1] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-266') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-266', node2.tup)
            
    # L->5_pl->T
    feats[2] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_pl') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '5_pl', node2.tup)
            
    # L->115->T
    feats[3] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '115') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '115', node2.tup)
            
    # L->-137->T
    feats[4] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-137') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-137', node2.tup)
            
    # L->voc. neutr.->T
    feats[5] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. neutr.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'voc. neutr.', node2.tup)
            
    # L->-77->T
    feats[6] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-77') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-77', node2.tup)
            
    # L->acc->T
    feats[7] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'acc', node2.tup)
            
    # L->149->T
    feats[8] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '149', node2.tup)
            
    # L->9_pl->T
    feats[9] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_pl') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '9_pl', node2.tup)
            
    # L->instr. du.->T
    feats[10] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. du.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'instr. du.', node2.tup)
            
    # L->-10->T
    feats[11] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-10') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-10', node2.tup)
            
    # L->160->T
    feats[12] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '160') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '160', node2.tup)
            
    # L->dat->T
    feats[13] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'dat', node2.tup)
            
    # L->48->T
    feats[14] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '48') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '48', node2.tup)
            
    # L->-161->T
    feats[15] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-161') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-161', node2.tup)
            
    # C->-268->L
    feats[16] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-268') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-268', node2.lemma)
            
    # C->36->L
    feats[17] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '36') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '36', node2.lemma)
            
    # T->-39->L
    feats[18] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-39') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-39', node2.lemma)
            
    # T->-71->L
    feats[19] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-71') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-71', node2.lemma)
            
    # T->nom. adj.->L
    feats[20] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'nom. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. adj.', node2.lemma)
            
    # T->94->L
    feats[21] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '94') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '94', node2.lemma)
            
    # T->-32->L
    feats[22] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-32') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-32', node2.lemma)
            
    # T->-121->L
    feats[23] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-121') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-121', node2.lemma)
            
    # T->-144->L
    feats[24] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-144') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-144', node2.lemma)
            
    # T->28->L
    feats[25] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '28') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
            
    # T->10_sp->C
    feats[26] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '10_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sp', node2.cng)
            
    # T->nom. masc.->C
    feats[27] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'nom. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. masc.', node2.cng)
            
    # T->-262->C
    feats[28] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-262') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-262', node2.cng)
            
    # T->34->C
    feats[29] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '34', node2.cng)
            
    # T->-83->C
    feats[30] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-83') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-83', node2.cng)
            
    # T->-78->C
    feats[31] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-78') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-78', node2.cng)
            
    # T->-126->C
    feats[32] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-126') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-126', node2.cng)
            
    # T->3_du->C
    feats[33] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '3_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_du', node2.cng)
            
    # T->81->C
    feats[34] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '81', node2.cng)
            
    # T->-48->C
    feats[35] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', node2.cng)
            
    # T->-297->C
    feats[36] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-297') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-297', node2.cng)
            
    # T->5_du->C
    feats[37] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '5_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_du', node2.cng)
            
    # T->90->C
    feats[38] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '90', node2.cng)
            
    # T->-25->C
    feats[39] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-25') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-25', node2.cng)
            
    # T->voc. sg.->C
    feats[40] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'voc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. sg.', node2.cng)
            
    # T->175->C
    feats[41] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '175') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '175', node2.cng)
            
    # T->119->C
    feats[42] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '119') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '119', node2.cng)
            
    # T->16_pl->C
    feats[43] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '16_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_pl', node2.cng)
            
    # T->4_fp->C
    feats[44] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '4_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_fp', node2.cng)
            
    # L -> 4_du -> 12_pl -> L
    feats[45] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_du', '12_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_pl', node2.lemma)
    
    # L -> 4_du -> 30_du -> L
    feats[46] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_du', '30_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_du', node2.lemma)
    
    # L -> 4_du -> -262 -> L
    feats[47] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_du', '-262') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-262', node2.lemma)
    
    # L -> 4_du -> -109 -> L
    feats[48] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_du', '-109') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-109', node2.lemma)
    
    # L -> 4_du -> tp -> L
    feats[49] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_du', 'tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'tp', node2.lemma)
    
    # L -> 4_du -> 34 -> L
    feats[50] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_du', '34') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '34', node2.lemma)
    
    # L -> 4_du -> 16_tp -> L
    feats[51] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_du', '16_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_tp', node2.lemma)
    
    # L -> 4_du -> 7_pl -> L
    feats[52] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_du', '7_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_pl', node2.lemma)
    
    # L -> 4_du -> voc. du. -> L
    feats[53] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_du', 'voc. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. du.', node2.lemma)
    
    # L -> 4_du -> 59 -> L
    feats[54] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_du', '59') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '59', node2.lemma)
    
    # L -> 4_du -> abl. sg. -> L
    feats[55] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_du', 'abl. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. sg.', node2.lemma)
    
    # L -> 4_du -> -126 -> L
    feats[56] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_du', '-126') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-126', node2.lemma)
    
    # L -> 4_du -> -29 -> L
    feats[57] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_du', '-29') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-29', node2.lemma)
    
    # L -> 4_du -> -122 -> L
    feats[58] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_du', '-122') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-122', node2.lemma)
    
    # L -> 4_du -> 14_sg -> L
    feats[59] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_du', '14_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_sg', node2.lemma)
    
    # L -> 4_du -> 9_tp -> L
    feats[60] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_du', '9_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_tp', node2.lemma)
    
    # L -> 4_du -> 37 -> L
    feats[61] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_du', '37') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '37', node2.lemma)
    
    # L -> 4_du -> -283 -> L
    feats[62] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_du', '-283') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-283', node2.lemma)
    
    # L -> 4_du -> 14_tp -> L
    feats[63] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_du', '14_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_tp', node2.lemma)
    
    # L -> 4_du -> 1 -> L
    feats[64] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_du', '1') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '1', node2.lemma)
    
    # L -> 4_du -> 16_pl -> L
    feats[65] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_du', '16_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_pl', node2.lemma)
    
    # L -> 4_du -> 118 -> L
    feats[66] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_du', '118') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '118', node2.lemma)
    
    # L -> 4_du -> -68 -> L
    feats[67] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_du', '-68') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-68', node2.lemma)
    
    # L -> 4_du -> 91 -> L
    feats[68] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_du', '91') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    
    # L -> 4_du -> 93 -> L
    feats[69] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_du', '93') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '93', node2.lemma)
    
    # L -> 4_du -> 38 -> L
    feats[70] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_du', '38') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '38', node2.lemma)
    
    # L -> 4_du -> 4_sp -> L
    feats[71] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_du', '4_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_sp', node2.lemma)
    
    # L -> 4_du -> 13_tp -> L
    feats[72] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_du', '13_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_tp', node2.lemma)
    
    # L -> 4_du -> 88 -> L
    feats[73] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_du', '88') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '88', node2.lemma)
    
    # L -> 4_du -> -17 -> L
    feats[74] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_du', '-17') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-17', node2.lemma)
    
    # L -> 4_du -> acc. fem -> L
    feats[75] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_du', 'acc. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. fem', node2.lemma)
    
    # L -> 4_du -> 12_sg -> L
    feats[76] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_du', '12_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_sg', node2.lemma)
    
    # L -> 4_du -> instr -> L
    feats[77] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_du', 'instr') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr', node2.lemma)
    
    # L -> 4_du -> 27_pl -> L
    feats[78] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_du', '27_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_pl', node2.lemma)
    
    # L -> 4_du -> 58 -> L
    feats[79] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_du', '58') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '58', node2.lemma)
    
    # L -> 4_du -> 8_sp -> L
    feats[80] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_du', '8_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_sp', node2.lemma)
    
    # L -> 30 -> 7_fp -> L
    feats[81] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30', '7_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_fp', node2.lemma)
    
    # L -> 30 -> -269 -> L
    feats[82] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30', '-269') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-269', node2.lemma)
    
    # L -> 30 -> dat. sg. -> L
    feats[83] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30', 'dat. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. sg.', node2.lemma)
    
    # L -> 30 -> 134 -> L
    feats[84] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30', '134') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '134', node2.lemma)
    
    # L -> 30 -> 9_sg -> L
    feats[85] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30', '9_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_sg', node2.lemma)
    
    # L -> 30 -> -119 -> L
    feats[86] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30', '-119') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-119', node2.lemma)
    
    # L -> 30 -> 32 -> L
    feats[87] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30', '32') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '32', node2.lemma)
    
    # L -> 30 -> 129 -> L
    feats[88] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30', '129') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '129', node2.lemma)
    
    # L -> 30 -> -21 -> L
    feats[89] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30', '-21') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-21', node2.lemma)
    
    # L -> 30 -> -37 -> L
    feats[90] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30', '-37') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-37', node2.lemma)
    
    # L -> 30 -> 94 -> L
    feats[91] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30', '94') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '94', node2.lemma)
    
    # L -> 30 -> -149 -> L
    feats[92] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30', '-149') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-149', node2.lemma)
    
    # L -> 30 -> -33 -> L
    feats[93] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30', '-33') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-33', node2.lemma)
    
    # L -> 30 -> -32 -> L
    feats[94] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30', '-32') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-32', node2.lemma)
    
    # L -> 30 -> -73 -> L
    feats[95] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30', '-73') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-73', node2.lemma)
    
    # L -> 30 -> -151 -> L
    feats[96] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30', '-151') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-151', node2.lemma)
    
    # L -> 30 -> 75 -> L
    feats[97] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30', '75') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '75', node2.lemma)
    
    # L -> 30 -> 11_pl -> L
    feats[98] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30', '11_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_pl', node2.lemma)
    
    # L -> 30 -> -190 -> L
    feats[99] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30', '-190') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-190', node2.lemma)
    
    # L -> 30 -> 51 -> L
    feats[100] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30', '51') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '51', node2.lemma)
    
    # L -> 30 -> -144 -> L
    feats[101] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30', '-144') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-144', node2.lemma)
    
    # L -> -268 -> 30 -> L
    feats[102] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', '30') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30', node2.lemma)
    
    # L -> -268 -> 14_pl -> L
    feats[103] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', '14_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_pl', node2.lemma)
    
    # L -> -268 -> -293 -> L
    feats[104] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', '-293') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-293', node2.lemma)
    
    # L -> -268 -> masc -> L
    feats[105] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', 'masc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'masc', node2.lemma)
    
    # L -> -268 -> 10_sp -> L
    feats[106] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', '10_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_sp', node2.lemma)
    
    # L -> -268 -> 12_pl -> L
    feats[107] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', '12_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_pl', node2.lemma)
    
    # L -> -268 -> 162 -> L
    feats[108] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', '162') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '162', node2.lemma)
    
    # L -> -268 -> 6_sg -> L
    feats[109] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', '6_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_sg', node2.lemma)
    
    # L -> -268 -> -28 -> L
    feats[110] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', '-28') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-28', node2.lemma)
    
    # L -> -268 -> tp -> L
    feats[111] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', 'tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'tp', node2.lemma)
    
    # L -> -268 -> -98 -> L
    feats[112] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', '-98') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-98', node2.lemma)
    
    # L -> -268 -> gen. du. -> L
    feats[113] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', 'gen. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen. du.', node2.lemma)
    
    # L -> -268 -> 16_tp -> L
    feats[114] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', '16_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_tp', node2.lemma)
    
    # L -> -268 -> -78 -> L
    feats[115] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', '-78') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-78', node2.lemma)
    
    # L -> -268 -> abl. pl. -> L
    feats[116] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', 'abl. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. pl.', node2.lemma)
    
    # L -> -268 -> du_tp -> L
    feats[117] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', 'du_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_tp', node2.lemma)
    
    # L -> -268 -> 117 -> L
    feats[118] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', '117') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '117', node2.lemma)
    
    # L -> -268 -> 138 -> L
    feats[119] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', '138') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '138', node2.lemma)
    
    # L -> -268 -> 90 -> L
    feats[120] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', '90') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '90', node2.lemma)
    
    # L -> -268 -> 136 -> L
    feats[121] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', '136') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '136', node2.lemma)
    
    # L -> -268 -> 9_sp -> L
    feats[122] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', '9_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_sp', node2.lemma)
    
    # L -> -268 -> voc. sg. -> L
    feats[123] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-268') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-268', 'voc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. sg.', node2.lemma)
    
    # L -> 14_pl -> du_tp -> L
    feats[124] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_pl', 'du_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_tp', node2.lemma)
    
    # L -> 14_pl -> -139 -> L
    feats[125] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_pl', '-139') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-139', node2.lemma)
    
    # L -> 14_pl -> -48 -> L
    feats[126] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_pl', '-48') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-48', node2.lemma)
    
    # L -> 14_pl -> -44 -> L
    feats[127] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_pl', '-44') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-44', node2.lemma)
    
    # L -> 14_pl -> -122 -> L
    feats[128] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_pl', '-122') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-122', node2.lemma)
    
    # L -> 14_pl -> -13 -> L
    feats[129] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_pl', '-13') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-13', node2.lemma)
    
    # L -> 14_pl -> 117 -> L
    feats[130] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_pl', '117') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '117', node2.lemma)
    
    # L -> 14_pl -> -91 -> L
    feats[131] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_pl', '-91') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-91', node2.lemma)
    
    # L -> 14_pl -> -84 -> L
    feats[132] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_pl', '-84') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-84', node2.lemma)
    
    # L -> 14_pl -> -242 -> L
    feats[133] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_pl', '-242') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-242', node2.lemma)
    
    # L -> 14_pl -> 14_tp -> L
    feats[134] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_pl', '14_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_tp', node2.lemma)
    
    # L -> 14_pl -> 74 -> L
    feats[135] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_pl', '74') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '74', node2.lemma)
    
    # L -> 14_pl -> 175 -> L
    feats[136] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_pl', '175') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '175', node2.lemma)
    
    # L -> 14_pl -> -307 -> L
    feats[137] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_pl', '-307') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-307', node2.lemma)
    
    # L -> 14_pl -> 96 -> L
    feats[138] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_pl', '96') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '96', node2.lemma)
    
    # L -> 14_pl -> 2_sg -> L
    feats[139] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_pl', '2_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_sg', node2.lemma)
    
    # L -> 14_pl -> -247 -> L
    feats[140] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_pl', '-247') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-247', node2.lemma)
    
    # L -> 14_pl -> acc. du. -> L
    feats[141] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_pl', 'acc. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. du.', node2.lemma)
    
    # L -> 14_pl -> 28_tp -> L
    feats[142] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_pl', '28_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28_tp', node2.lemma)
    
    # L -> 14_pl -> -111 -> L
    feats[143] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_pl', '-111') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-111', node2.lemma)
    
    # L -> 14_pl -> nom. fem -> L
    feats[144] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_pl', 'nom. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. fem', node2.lemma)
    
    # L -> 14_pl -> 1 -> L
    feats[145] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_pl', '1') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '1', node2.lemma)
    
    # L -> 14_pl -> fp -> L
    feats[146] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_pl', 'fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'fp', node2.lemma)
    
    # L -> 14_pl -> -220 -> L
    feats[147] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_pl', '-220') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-220', node2.lemma)
    
    # L -> 14_pl -> 6_tp -> L
    feats[148] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_pl', '6_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_tp', node2.lemma)
    
    # L -> 14_pl -> 93 -> L
    feats[149] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_pl', '93') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '93', node2.lemma)
    
    # L -> 14_pl -> 69 -> L
    feats[150] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_pl', '69') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '69', node2.lemma)
    
    # L -> 14_pl -> -92 -> L
    feats[151] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_pl', '-92') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-92', node2.lemma)
    
    # L -> 14_pl -> 38 -> L
    feats[152] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_pl', '38') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '38', node2.lemma)
    
    # L -> 14_pl -> instr. masc. -> L
    feats[153] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_pl', 'instr. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. masc.', node2.lemma)
    
    # L -> 14_pl -> 88 -> L
    feats[154] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_pl', '88') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '88', node2.lemma)
    
    # L -> 14_pl -> 30_fp -> L
    feats[155] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_pl', '30_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_fp', node2.lemma)
    
    # L -> 14_pl -> adj -> L
    feats[156] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_pl', 'adj') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'adj', node2.lemma)
    
    # L -> 14_pl -> -72 -> L
    feats[157] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_pl', '-72') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-72', node2.lemma)
    
    # L -> 14_pl -> du_fp -> L
    feats[158] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_pl', 'du_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_fp', node2.lemma)
    
    # L -> 14_pl -> 156 -> L
    feats[159] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_pl', '156') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '156', node2.lemma)
    
    # L -> 14_pl -> 30_pl -> L
    feats[160] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_pl', '30_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_pl', node2.lemma)
    
    # L -> 14_pl -> 15_sg -> L
    feats[161] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_pl', '15_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_sg', node2.lemma)
    
    # L -> 14_pl -> 50 -> L
    feats[162] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_pl', '50') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '50', node2.lemma)
    
    # L -> 14_pl -> nom. pl. -> L
    feats[163] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_pl', 'nom. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. pl.', node2.lemma)
    
    # L -> 14_pl -> -15 -> L
    feats[164] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_pl', '-15') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-15', node2.lemma)
    
    # L -> 14_pl -> 12_tp -> L
    feats[165] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_pl', '12_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_tp', node2.lemma)
    
    # L -> 14_pl -> 8_pl -> L
    feats[166] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_pl', '8_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_pl', node2.lemma)
    
    # L -> 14_pl -> -249 -> L
    feats[167] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_pl', '-249') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-249', node2.lemma)
    
    # L -> 14_pl -> -245 -> L
    feats[168] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_pl', '-245') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-245', node2.lemma)
    
    # L -> 14_pl -> -158 -> L
    feats[169] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_pl', '-158') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-158', node2.lemma)
    
    # L -> 14_pl -> 4_pl -> L
    feats[170] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_pl', '4_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_pl', node2.lemma)
    
    # L -> 14_pl -> -24 -> L
    feats[171] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_pl', '-24') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    
    # L -> 14_pl -> loc. du. -> L
    feats[172] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_pl', 'loc. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc. du.', node2.lemma)
    
    # L -> 14_pl -> -104 -> L
    feats[173] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_pl', '-104') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-104', node2.lemma)
    
    # L -> -293 -> -240 -> L
    feats[174] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '-240') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-240', node2.lemma)
    
    # L -> -293 -> 121 -> L
    feats[175] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '121') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '121', node2.lemma)
    
    # L -> -293 -> 11_sp -> L
    feats[176] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '11_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_sp', node2.lemma)
    
    # L -> -293 -> 110 -> L
    feats[177] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '110') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '110', node2.lemma)
    
    # L -> -293 -> 77 -> L
    feats[178] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '77') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '77', node2.lemma)
    
    # L -> -293 -> 8_pl -> L
    feats[179] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '8_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_pl', node2.lemma)
    
    # L -> -293 -> -249 -> L
    feats[180] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '-249') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-249', node2.lemma)
    
    # L -> -293 -> nom. du. -> L
    feats[181] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', 'nom. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. du.', node2.lemma)
    
    # L -> -293 -> 101 -> L
    feats[182] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '101') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    
    # L -> -293 -> 10_sg -> L
    feats[183] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '10_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_sg', node2.lemma)
    
    # L -> -293 -> 173 -> L
    feats[184] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '173') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '173', node2.lemma)
    
    # L -> -293 -> 5_tp -> L
    feats[185] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '5_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_tp', node2.lemma)
    
    # L -> -293 -> abl -> L
    feats[186] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', 'abl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl', node2.lemma)
    
    # L -> -293 -> -101 -> L
    feats[187] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '-101') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-101', node2.lemma)
    
    # L -> -293 -> 95 -> L
    feats[188] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '95') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '95', node2.lemma)
    
    # L -> -293 -> 111 -> L
    feats[189] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '111') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '111', node2.lemma)
    
    # L -> -293 -> -271 -> L
    feats[190] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '-271') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-271', node2.lemma)
    
    # L -> -293 -> -46 -> L
    feats[191] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '-46') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-46', node2.lemma)
    
    # L -> -293 -> -35 -> L
    feats[192] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '-35') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    
    # L -> -293 -> abl. du. -> L
    feats[193] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', 'abl. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. du.', node2.lemma)
    
    # L -> -293 -> 9_du -> L
    feats[194] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '9_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_du', node2.lemma)
    
    # L -> -293 -> -119 -> L
    feats[195] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '-119') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-119', node2.lemma)
    
    # L -> masc -> 13_pl -> L
    feats[196] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'masc') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'masc', '13_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_pl', node2.lemma)
    
    # L -> masc -> -141 -> L
    feats[197] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'masc') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'masc', '-141') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-141', node2.lemma)
    
    # L -> masc -> 16_fp -> L
    feats[198] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'masc') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'masc', '16_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_fp', node2.lemma)
    
    # L -> masc -> tp -> L
    feats[199] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'masc') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'masc', 'tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'tp', node2.lemma)
    
    # L -> masc -> -308 -> L
    feats[200] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'masc') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'masc', '-308') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-308', node2.lemma)
    
    # L -> masc -> -147 -> L
    feats[201] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'masc') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'masc', '-147') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-147', node2.lemma)
    
    # L -> masc -> 13_sp -> L
    feats[202] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'masc') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'masc', '13_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_sp', node2.lemma)
    
    # L -> masc -> 177 -> L
    feats[203] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'masc') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'masc', '177') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '177', node2.lemma)
    
    # L -> masc -> instr. du. -> L
    feats[204] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'masc') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'masc', 'instr. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. du.', node2.lemma)
    
    # L -> 13_pl -> 10_tp -> L
    feats[205] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '10_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_tp', node2.lemma)
    
    # L -> 13_pl -> 33 -> L
    feats[206] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '33') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '33', node2.lemma)
    
    # L -> 13_pl -> nom. du. -> L
    feats[207] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', 'nom. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. du.', node2.lemma)
    
    # L -> 13_pl -> -291 -> L
    feats[208] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '-291') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-291', node2.lemma)
    
    # L -> 13_pl -> -166 -> L
    feats[209] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '-166') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-166', node2.lemma)
    
    # L -> 13_pl -> instr. sg. -> L
    feats[210] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', 'instr. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. sg.', node2.lemma)
    
    # L -> 13_pl -> 101 -> L
    feats[211] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '101') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    
    # L -> 13_pl -> -261 -> L
    feats[212] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '-261') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-261', node2.lemma)
    
    # L -> 13_pl -> 109 -> L
    feats[213] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '109') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '109', node2.lemma)
    
    # L -> 13_pl -> sg_fp -> L
    feats[214] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', 'sg_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_fp', node2.lemma)
    
    # L -> 13_pl -> 29_tp -> L
    feats[215] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '29_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '29_tp', node2.lemma)
    
    # L -> 13_pl -> 10_sg -> L
    feats[216] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '10_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_sg', node2.lemma)
    
    # L -> 13_pl -> sg_sp -> L
    feats[217] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', 'sg_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_sp', node2.lemma)
    
    # L -> 13_pl -> 5_pl -> L
    feats[218] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '5_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_pl', node2.lemma)
    
    # L -> 13_pl -> -45 -> L
    feats[219] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '-45') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-45', node2.lemma)
    
    # L -> 13_pl -> -27 -> L
    feats[220] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '-27') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-27', node2.lemma)
    
    # L -> 13_pl -> 12_sp -> L
    feats[221] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '12_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_sp', node2.lemma)
    
    # L -> 13_pl -> 40 -> L
    feats[222] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '40') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '40', node2.lemma)
    
    # L -> 13_pl -> -156 -> L
    feats[223] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '-156') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-156', node2.lemma)
    
    # L -> 13_pl -> 89 -> L
    feats[224] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '89') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '89', node2.lemma)
    
    # L -> 13_pl -> 78 -> L
    feats[225] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '78') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '78', node2.lemma)
    
    # L -> 13_pl -> -246 -> L
    feats[226] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '-246') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-246', node2.lemma)
    
    # L -> 13_pl -> 115 -> L
    feats[227] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '115') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '115', node2.lemma)
    
    # L -> 13_pl -> 173 -> L
    feats[228] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '173') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '173', node2.lemma)
    
    # L -> 13_pl -> pl_fp -> L
    feats[229] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', 'pl_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_fp', node2.lemma)
    
    # L -> 13_pl -> 5_tp -> L
    feats[230] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '5_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_tp', node2.lemma)
    
    # L -> 13_pl -> abl -> L
    feats[231] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', 'abl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl', node2.lemma)
    
    # L -> 13_pl -> 7_sp -> L
    feats[232] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '7_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_sp', node2.lemma)
    
    # L -> 13_pl -> acc. neutr. -> L
    feats[233] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', 'acc. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. neutr.', node2.lemma)
    
    # L -> 13_pl -> 11_du -> L
    feats[234] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '11_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_du', node2.lemma)
    
    # L -> 13_pl -> 137 -> L
    feats[235] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '137') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    
    # L -> 13_pl -> 70 -> L
    feats[236] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '70') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '70', node2.lemma)
    
    # L -> 13_pl -> 169 -> L
    feats[237] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '169') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '169', node2.lemma)
    
    # L -> 13_pl -> 95 -> L
    feats[238] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '95') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '95', node2.lemma)
    
    # L -> 13_pl -> acc. adj. -> L
    feats[239] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', 'acc. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. adj.', node2.lemma)
    
    # L -> 13_pl -> -11 -> L
    feats[240] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '-11') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-11', node2.lemma)
    
    # L -> 13_pl -> 139 -> L
    feats[241] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '139') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '139', node2.lemma)
    
    # L -> 13_pl -> acc -> L
    feats[242] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', 'acc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc', node2.lemma)
    
    # L -> 13_pl -> 102 -> L
    feats[243] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '102') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '102', node2.lemma)
    
    # L -> 13_pl -> -22 -> L
    feats[244] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '-22') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-22', node2.lemma)
    
    # L -> 13_pl -> 15_fp -> L
    feats[245] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '15_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_fp', node2.lemma)
    
    # L -> 13_pl -> 3_pl -> L
    feats[246] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '3_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_pl', node2.lemma)
    
    # L -> 13_pl -> 16_du -> L
    feats[247] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '16_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_du', node2.lemma)
    
    # L -> 13_pl -> 15_pl -> L
    feats[248] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '15_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_pl', node2.lemma)
    
    # L -> 13_pl -> 7_tp -> L
    feats[249] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '7_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_tp', node2.lemma)
    
    # L -> 13_pl -> -46 -> L
    feats[250] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '-46') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-46', node2.lemma)
    
    # L -> 13_pl -> -58 -> L
    feats[251] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '-58') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-58', node2.lemma)
    
    # L -> 13_pl -> -269 -> L
    feats[252] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '-269') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-269', node2.lemma)
    
    # L -> 13_pl -> dat. sg. -> L
    feats[253] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', 'dat. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. sg.', node2.lemma)
    
    # L -> 13_pl -> instr. du. -> L
    feats[254] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', 'instr. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. du.', node2.lemma)
    
    # L -> 13_pl -> -82 -> L
    feats[255] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '-82') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-82', node2.lemma)
    
    # L -> 13_pl -> 9_du -> L
    feats[256] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '9_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_du', node2.lemma)
    
    # L -> 13_pl -> 61 -> L
    feats[257] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '61') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '61', node2.lemma)
    
    # L -> 13_pl -> -87 -> L
    feats[258] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '-87') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-87', node2.lemma)
    
    # L -> 13_pl -> 140 -> L
    feats[259] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '140') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '140', node2.lemma)
    
    # L -> 13_pl -> 6_pl -> L
    feats[260] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '6_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_pl', node2.lemma)
    
    # L -> 13_pl -> -133 -> L
    feats[261] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '-133') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-133', node2.lemma)
    
    # L -> 13_pl -> -169 -> L
    feats[262] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '-169') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-169', node2.lemma)
    
    # L -> 13_pl -> 3 -> L
    feats[263] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '3') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3', node2.lemma)
    
    # L -> 13_pl -> 2_fp -> L
    feats[264] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '2_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_fp', node2.lemma)
    
    # L -> 13_pl -> 56 -> L
    feats[265] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '56') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '56', node2.lemma)
    
    # L -> 13_pl -> 15_tp -> L
    feats[266] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '15_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_tp', node2.lemma)
    
    # L -> 13_pl -> -64 -> L
    feats[267] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '-64') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-64', node2.lemma)
    
    # L -> 13_pl -> -10 -> L
    feats[268] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '-10') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-10', node2.lemma)
    
    # L -> 13_pl -> -123 -> L
    feats[269] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '-123') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-123', node2.lemma)
    
    # L -> 13_pl -> -21 -> L
    feats[270] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '-21') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-21', node2.lemma)
    
    # L -> 13_pl -> nom. neutr. -> L
    feats[271] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', 'nom. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. neutr.', node2.lemma)
    
    # L -> 13_pl -> 8_sg -> L
    feats[272] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '8_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_sg', node2.lemma)
    
    # L -> 13_pl -> 49 -> L
    feats[273] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '49') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '49', node2.lemma)
    
    # L -> 13_pl -> -39 -> L
    feats[274] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', '-39') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-39', node2.lemma)
    
    # L -> 132 -> 157 -> L
    feats[275] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '132') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '132', '157') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '157', node2.lemma)
    
    # L -> 132 -> 2_tp -> L
    feats[276] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '132') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '132', '2_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_tp', node2.lemma)
    
    # L -> 132 -> 11_du -> L
    feats[277] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '132') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '132', '11_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_du', node2.lemma)
    
    # L -> 132 -> -137 -> L
    feats[278] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '132') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '132', '-137') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-137', node2.lemma)
    
    # L -> 132 -> 176 -> L
    feats[279] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '132') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '132', '176') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '176', node2.lemma)
    
    # L -> 132 -> voc. neutr. -> L
    feats[280] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '132') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '132', 'voc. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. neutr.', node2.lemma)
    
    # L -> 132 -> 35 -> L
    feats[281] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '132') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '132', '35') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '35', node2.lemma)
    
    # L -> 132 -> acc. adj. -> L
    feats[282] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '132') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '132', 'acc. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. adj.', node2.lemma)
    
    # L -> 132 -> -22 -> L
    feats[283] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '132') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '132', '-22') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-22', node2.lemma)
    
    # L -> 132 -> 15_pl -> L
    feats[284] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '132') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '132', '15_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_pl', node2.lemma)
    
    # L -> 132 -> 6_du -> L
    feats[285] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '132') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '132', '6_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_du', node2.lemma)
    
    # L -> 132 -> -46 -> L
    feats[286] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '132') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '132', '-46') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-46', node2.lemma)
    
    # L -> 132 -> fem -> L
    feats[287] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '132') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '132', 'fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'fem', node2.lemma)
    
    # L -> 132 -> 149 -> L
    feats[288] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '132') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '132', '149') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '149', node2.lemma)
    
    # L -> 132 -> 9_sg -> L
    feats[289] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '132') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '132', '9_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_sg', node2.lemma)
    
    # L -> 132 -> 140 -> L
    feats[290] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '132') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '132', '140') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '140', node2.lemma)
    
    # L -> 132 -> -133 -> L
    feats[291] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '132') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '132', '-133') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-133', node2.lemma)
    
    # L -> 132 -> -169 -> L
    feats[292] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '132') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '132', '-169') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-169', node2.lemma)
    
    # L -> 132 -> 129 -> L
    feats[293] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '132') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '132', '129') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '129', node2.lemma)
    
    # L -> 132 -> -10 -> L
    feats[294] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '132') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '132', '-10') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-10', node2.lemma)
    
    # L -> 132 -> 100 -> L
    feats[295] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '132') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '132', '100') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '100', node2.lemma)
    
    # L -> 132 -> nom. neutr. -> L
    feats[296] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '132') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '132', 'nom. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. neutr.', node2.lemma)
    
    # L -> 132 -> 82 -> L
    feats[297] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '132') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '132', '82') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '82', node2.lemma)
    
    # L -> 132 -> -39 -> L
    feats[298] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '132') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '132', '-39') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-39', node2.lemma)
    
    # L -> 132 -> instr. pl. -> L
    feats[299] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '132') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '132', 'instr. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. pl.', node2.lemma)
    
    # L -> 132 -> -121 -> L
    feats[300] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '132') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '132', '-121') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-121', node2.lemma)
    
    # L -> 132 -> 11_pl -> L
    feats[301] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '132') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '132', '11_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_pl', node2.lemma)
    
    # L -> 132 -> 48 -> L
    feats[302] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '132') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '132', '48') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '48', node2.lemma)
    
    # L -> 132 -> 168 -> L
    feats[303] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '132') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '132', '168') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '168', node2.lemma)
    
    # L -> 132 -> -161 -> L
    feats[304] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '132') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '132', '-161') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-161', node2.lemma)
    
    # L -> 132 -> 28 -> L
    feats[305] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '132') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '132', '28') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
    
    # L -> 36 -> voc. masc. -> L
    feats[306] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '36') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '36', 'voc. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. masc.', node2.lemma)
    
    # L -> 36 -> 4_du -> L
    feats[307] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '36') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '36', '4_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_du', node2.lemma)
    
    # L -> 36 -> -293 -> L
    feats[308] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '36') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '36', '-293') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-293', node2.lemma)
    
    # L -> 36 -> 3_sp -> L
    feats[309] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '36') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '36', '3_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_sp', node2.lemma)
    
    # L -> 36 -> 5_sg -> L
    feats[310] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '36') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '36', '5_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sg', node2.lemma)
    
    # L -> 36 -> 6_sg -> L
    feats[311] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '36') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '36', '6_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_sg', node2.lemma)
    
    # L -> 36 -> 34 -> L
    feats[312] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '36') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '36', '34') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '34', node2.lemma)
    
    # L -> 36 -> -18 -> L
    feats[313] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '36') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '36', '-18') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-18', node2.lemma)
    
    # L -> 36 -> 16_tp -> L
    feats[314] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '36') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '36', '16_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_tp', node2.lemma)
    
    # L -> 36 -> 11_sg -> L
    feats[315] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '36') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '36', '11_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_sg', node2.lemma)
    
    # L -> 36 -> 59 -> L
    feats[316] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '36') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '36', '59') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '59', node2.lemma)
    
    # L -> 36 -> -302 -> L
    feats[317] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '36') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '36', '-302') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-302', node2.lemma)
    
    # L -> nom. masc. -> -245 -> L
    feats[318] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. masc.', '-245') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-245', node2.lemma)
    
    # L -> nom. masc. -> 97 -> L
    feats[319] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. masc.', '97') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '97', node2.lemma)
    
    # L -> nom. masc. -> -166 -> L
    feats[320] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. masc.', '-166') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-166', node2.lemma)
    
    # L -> nom. masc. -> 161 -> L
    feats[321] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. masc.', '161') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '161', node2.lemma)
    
    # L -> nom. masc. -> dat. du. -> L
    feats[322] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. masc.', 'dat. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    
    # L -> nom. masc. -> 169 -> L
    feats[323] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. masc.', '169') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '169', node2.lemma)
    
    # L -> nom. masc. -> 16_du -> L
    feats[324] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. masc.', '16_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_du', node2.lemma)
    
    # L -> nom. masc. -> -87 -> L
    feats[325] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. masc.', '-87') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-87', node2.lemma)
    
    # L -> nom. masc. -> -169 -> L
    feats[326] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. masc.', '-169') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-169', node2.lemma)
    
    # L -> nom. masc. -> -37 -> L
    feats[327] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. masc.', '-37') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-37', node2.lemma)
    
    # L -> nom. masc. -> 92 -> L
    feats[328] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. masc.', '92') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '92', node2.lemma)
    
    # L -> nom. masc. -> -94 -> L
    feats[329] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. masc.', '-94') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-94', node2.lemma)
    
    # L -> nom. masc. -> 168 -> L
    feats[330] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. masc.', '168') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '168', node2.lemma)
    
    # L -> nom. masc. -> 16_sg -> L
    feats[331] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. masc.', '16_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_sg', node2.lemma)
    
    # L -> loc -> -293 -> L
    feats[332] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc', '-293') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-293', node2.lemma)
    
    # L -> loc -> instr. adj. -> L
    feats[333] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc', 'instr. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. adj.', node2.lemma)
    
    # L -> loc -> 7_pl -> L
    feats[334] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc', '7_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_pl', node2.lemma)
    
    # L -> loc -> 27_sg -> L
    feats[335] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc', '27_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_sg', node2.lemma)
    
    # L -> loc -> -29 -> L
    feats[336] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc', '-29') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-29', node2.lemma)
    
    # L -> loc -> 135 -> L
    feats[337] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc', '135') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '135', node2.lemma)
    
    # L -> loc -> 117 -> L
    feats[338] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc', '117') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '117', node2.lemma)
    
    # L -> loc -> -283 -> L
    feats[339] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc', '-283') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-283', node2.lemma)
    
    # L -> loc -> voc. sg. -> L
    feats[340] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc', 'voc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. sg.', node2.lemma)
    
    # L -> loc -> 30_sp -> L
    feats[341] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc', '30_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_sp', node2.lemma)
    
    # L -> loc -> 16_pl -> L
    feats[342] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc', '16_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_pl', node2.lemma)
    
    # L -> -296 -> 12_fp -> L
    feats[343] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-296') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-296', '12_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_fp', node2.lemma)
    
    # L -> -296 -> 8_du -> L
    feats[344] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-296') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-296', '8_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_du', node2.lemma)
    
    # L -> -296 -> 27_sg -> L
    feats[345] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-296') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-296', '27_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_sg', node2.lemma)
    
    # L -> -296 -> du_tp -> L
    feats[346] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-296') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-296', 'du_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_tp', node2.lemma)
    
    # L -> -296 -> -273 -> L
    feats[347] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-296') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-296', '-273') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-273', node2.lemma)
    
    # L -> -296 -> -91 -> L
    feats[348] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-296') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-296', '-91') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-91', node2.lemma)
    
    # L -> -296 -> 90 -> L
    feats[349] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-296') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-296', '90') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '90', node2.lemma)
    
    # L -> -296 -> sp -> L
    feats[350] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-296') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-296', 'sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sp', node2.lemma)
    
    # L -> -296 -> 9_sp -> L
    feats[351] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-296') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-296', '9_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_sp', node2.lemma)
    
    # L -> -296 -> -307 -> L
    feats[352] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-296') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-296', '-307') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-307', node2.lemma)
    
    # L -> -296 -> 3_sg -> L
    feats[353] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-296') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-296', '3_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_sg', node2.lemma)
    
    # L -> -296 -> 118 -> L
    feats[354] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-296') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-296', '118') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '118', node2.lemma)
    
    # L -> -296 -> 13_tp -> L
    feats[355] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-296') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-296', '13_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_tp', node2.lemma)
    
    # L -> -296 -> -19 -> L
    feats[356] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-296') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-296', '-19') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-19', node2.lemma)
    
    # L -> -296 -> 30_fp -> L
    feats[357] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-296') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-296', '30_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_fp', node2.lemma)
    
    # L -> -296 -> gen. sg. -> L
    feats[358] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-296') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-296', 'gen. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen. sg.', node2.lemma)
    
    # L -> -296 -> -31 -> L
    feats[359] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-296') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-296', '-31') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-31', node2.lemma)
    
    # L -> -296 -> -72 -> L
    feats[360] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-296') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-296', '-72') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-72', node2.lemma)
    
    # L -> -296 -> du_fp -> L
    feats[361] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-296') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-296', 'du_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_fp', node2.lemma)
    
    # L -> -296 -> 156 -> L
    feats[362] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-296') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-296', '156') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '156', node2.lemma)
    
    # L -> -296 -> -96 -> L
    feats[363] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-296') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-296', '-96') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-96', node2.lemma)
    
    # L -> -296 -> -62 -> L
    feats[364] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-296') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-296', '-62') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-62', node2.lemma)
    
    # L -> -296 -> 121 -> L
    feats[365] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-296') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-296', '121') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '121', node2.lemma)
    
    # L -> -296 -> 58 -> L
    feats[366] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-296') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-296', '58') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '58', node2.lemma)
    
    # L -> -296 -> 41 -> L
    feats[367] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-296') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-296', '41') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '41', node2.lemma)
    
    # L -> -296 -> 30_tp -> L
    feats[368] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-296') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-296', '30_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_tp', node2.lemma)
    
    # L -> -296 -> 110 -> L
    feats[369] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-296') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-296', '110') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '110', node2.lemma)
    
    # L -> -296 -> -14 -> L
    feats[370] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-296') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-296', '-14') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-14', node2.lemma)
    
    # L -> -296 -> 13_sg -> L
    feats[371] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-296') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-296', '13_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_sg', node2.lemma)
    
    # L -> -210 -> 138 -> L
    feats[372] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '138') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '138', node2.lemma)
    
    # L -> -210 -> 5_du -> L
    feats[373] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '5_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_du', node2.lemma)
    
    # L -> -210 -> -84 -> L
    feats[374] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '-84') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-84', node2.lemma)
    
    # L -> -210 -> 90 -> L
    feats[375] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '90') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '90', node2.lemma)
    
    # L -> -210 -> -52 -> L
    feats[376] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '-52') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-52', node2.lemma)
    
    # L -> -210 -> 9_sp -> L
    feats[377] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '9_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_sp', node2.lemma)
    
    # L -> -210 -> 14_tp -> L
    feats[378] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '14_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_tp', node2.lemma)
    
    # L -> -210 -> 74 -> L
    feats[379] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '74') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '74', node2.lemma)
    
    # L -> -210 -> 175 -> L
    feats[380] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '175') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '175', node2.lemma)
    
    # L -> -210 -> -132 -> L
    feats[381] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '-132') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-132', node2.lemma)
    
    # L -> -210 -> 2_sg -> L
    feats[382] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '2_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_sg', node2.lemma)
    
    # L -> -210 -> 119 -> L
    feats[383] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '119') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '119', node2.lemma)
    
    # L -> -210 -> nom. fem -> L
    feats[384] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', 'nom. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. fem', node2.lemma)
    
    # L -> -210 -> -220 -> L
    feats[385] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '-220') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-220', node2.lemma)
    
    # L -> -210 -> 39 -> L
    feats[386] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '39') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '39', node2.lemma)
    
    # L -> -210 -> -68 -> L
    feats[387] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '-68') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-68', node2.lemma)
    
    # L -> -210 -> 182 -> L
    feats[388] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '182') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '182', node2.lemma)
    
    # L -> -210 -> pl -> L
    feats[389] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', 'pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    
    # L -> -210 -> -92 -> L
    feats[390] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '-92') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-92', node2.lemma)
    
    # L -> -210 -> 3_tp -> L
    feats[391] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '3_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_tp', node2.lemma)
    
    # L -> -210 -> 13_tp -> L
    feats[392] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '13_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_tp', node2.lemma)
    
    # L -> -210 -> -23 -> L
    feats[393] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '-23') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-23', node2.lemma)
    
    # L -> -210 -> sg_tp -> L
    feats[394] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', 'sg_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_tp', node2.lemma)
    
    # L -> -210 -> 88 -> L
    feats[395] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '88') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '88', node2.lemma)
    
    # L -> -210 -> sg -> L
    feats[396] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', 'sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg', node2.lemma)
    
    # L -> -210 -> gen. sg. -> L
    feats[397] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', 'gen. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen. sg.', node2.lemma)
    
    # L -> -210 -> du_fp -> L
    feats[398] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', 'du_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_fp', node2.lemma)
    
    # L -> -210 -> -96 -> L
    feats[399] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '-96') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-96', node2.lemma)
    
    # L -> -210 -> 15_sg -> L
    feats[400] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '15_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_sg', node2.lemma)
    
    # L -> -210 -> 10_du -> L
    feats[401] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '10_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_du', node2.lemma)
    
    # L -> -210 -> 30_tp -> L
    feats[402] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '30_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_tp', node2.lemma)
    
    # L -> -210 -> nom. pl. -> L
    feats[403] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', 'nom. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. pl.', node2.lemma)
    
    # L -> -210 -> -15 -> L
    feats[404] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '-15') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-15', node2.lemma)
    
    # L -> -210 -> -55 -> L
    feats[405] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '-55') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-55', node2.lemma)
    
    # L -> -210 -> 2 -> L
    feats[406] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '2') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2', node2.lemma)
    
    # L -> -210 -> 120 -> L
    feats[407] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '120') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '120', node2.lemma)
    
    # L -> -210 -> 181 -> L
    feats[408] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '181') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '181', node2.lemma)
    
    # L -> -210 -> -113 -> L
    feats[409] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '-113') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-113', node2.lemma)
    
    # L -> -210 -> -158 -> L
    feats[410] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '-158') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-158', node2.lemma)
    
    # L -> -210 -> gen. pl. -> L
    feats[411] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', 'gen. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen. pl.', node2.lemma)
    
    # L -> -210 -> -263 -> L
    feats[412] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '-263') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-263', node2.lemma)
    
    # L -> -210 -> nom -> L
    feats[413] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', 'nom') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom', node2.lemma)
    
    # L -> -210 -> 151 -> L
    feats[414] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '151') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '151', node2.lemma)
    
    # L -> -210 -> 71 -> L
    feats[415] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '71') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '71', node2.lemma)
    
    # L -> -210 -> -153 -> L
    feats[416] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '-153') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-153', node2.lemma)
    
    # L -> -210 -> -53 -> L
    feats[417] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '-53') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-53', node2.lemma)
    
    # L -> -210 -> 14_sp -> L
    feats[418] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '14_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_sp', node2.lemma)
    
    # L -> -210 -> -261 -> L
    feats[419] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', '-261') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-261', node2.lemma)
    
    # L -> -210 -> gen -> L
    feats[420] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', 'gen') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen', node2.lemma)
    
    # L -> 30_du -> pl_sp -> L
    feats[421] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_du', 'pl_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_sp', node2.lemma)
    
    # L -> 30_du -> 8_fp -> L
    feats[422] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_du', '8_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    
    # L -> 30_du -> abl. sg. -> L
    feats[423] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_du', 'abl. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. sg.', node2.lemma)
    
    # L -> 30_du -> 3_du -> L
    feats[424] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_du', '3_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    
    # L -> 30_du -> -29 -> L
    feats[425] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_du', '-29') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-29', node2.lemma)
    
    # L -> 30_du -> 2_sg -> L
    feats[426] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_du', '2_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_sg', node2.lemma)
    
    # L -> 30_du -> -111 -> L
    feats[427] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_du', '-111') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-111', node2.lemma)
    
    # L -> 30_du -> -220 -> L
    feats[428] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_du', '-220') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-220', node2.lemma)
    
    # L -> 30_du -> 91 -> L
    feats[429] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_du', '91') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    
    # L -> 30_du -> -41 -> L
    feats[430] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_du', '-41') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-41', node2.lemma)
    
    # L -> 30_du -> 12_sg -> L
    feats[431] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_du', '12_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_sg', node2.lemma)
    
    # L -> 30_du -> 112 -> L
    feats[432] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_du', '112') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '112', node2.lemma)
    
    # L -> -262 -> -103 -> L
    feats[433] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-262') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-262', '-103') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-103', node2.lemma)
    
    # L -> -262 -> 3_fp -> L
    feats[434] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-262') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-262', '3_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_fp', node2.lemma)
    
    # L -> -262 -> 160 -> L
    feats[435] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-262') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-262', '160') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    
    # L -> -262 -> nom. adj. -> L
    feats[436] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-262') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-262', 'nom. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. adj.', node2.lemma)
    
    # L -> -262 -> 168 -> L
    feats[437] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-262') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-262', '168') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '168', node2.lemma)
    
    # L -> -262 -> -161 -> L
    feats[438] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-262') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-262', '-161') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-161', node2.lemma)
    
    # L -> -262 -> 180 -> L
    feats[439] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-262') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-262', '180') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '180', node2.lemma)
    
    # L -> 162 -> -268 -> L
    feats[440] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '162') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '162', '-268') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-268', node2.lemma)
    
    # L -> 162 -> 14_pl -> L
    feats[441] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '162') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '162', '14_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_pl', node2.lemma)
    
    # L -> 162 -> -293 -> L
    feats[442] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '162') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '162', '-293') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-293', node2.lemma)
    
    # L -> 162 -> 10_sp -> L
    feats[443] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '162') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '162', '10_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_sp', node2.lemma)
    
    # L -> 162 -> 13_pl -> L
    feats[444] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '162') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '162', '13_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_pl', node2.lemma)
    
    # L -> 162 -> -152 -> L
    feats[445] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '162') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '162', '-152') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-152', node2.lemma)
    
    # L -> 162 -> 132 -> L
    feats[446] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '162') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '162', '132') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '132', node2.lemma)
    
    # L -> 162 -> loc -> L
    feats[447] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '162') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '162', 'loc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc', node2.lemma)
    
    # L -> 6_sg -> 81 -> L
    feats[448] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', '81') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '81', node2.lemma)
    
    # L -> 6_sg -> -44 -> L
    feats[449] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', '-44') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-44', node2.lemma)
    
    # L -> 6_sg -> -242 -> L
    feats[450] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', '-242') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-242', node2.lemma)
    
    # L -> 6_sg -> 2_pl -> L
    feats[451] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', '2_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_pl', node2.lemma)
    
    # L -> 55 -> 6_du -> L
    feats[452] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '55') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '55', '6_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_du', node2.lemma)
    
    # L -> 55 -> abl. du. -> L
    feats[453] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '55') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '55', 'abl. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. du.', node2.lemma)
    
    # L -> 55 -> -143 -> L
    feats[454] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '55') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '55', '-143') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-143', node2.lemma)
    
    # L -> 55 -> 149 -> L
    feats[455] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '55') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '55', '149') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '149', node2.lemma)
    
    # L -> 55 -> 177 -> L
    feats[456] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '55') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '55', '177') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '177', node2.lemma)
    
    # L -> 55 -> 9_pl -> L
    feats[457] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '55') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '55', '9_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_pl', node2.lemma)
    
    # L -> 55 -> instr. du. -> L
    feats[458] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '55') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '55', 'instr. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. du.', node2.lemma)
    
    # L -> 55 -> 9_sg -> L
    feats[459] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '55') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '55', '9_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_sg', node2.lemma)
    
    # L -> 55 -> -119 -> L
    feats[460] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '55') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '55', '-119') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-119', node2.lemma)
    
    # L -> 55 -> -169 -> L
    feats[461] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '55') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '55', '-169') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-169', node2.lemma)
    
    # L -> 55 -> -123 -> L
    feats[462] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '55') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '55', '-123') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-123', node2.lemma)
    
    # L -> 55 -> nom. neutr. -> L
    feats[463] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '55') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '55', 'nom. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. neutr.', node2.lemma)
    
    # L -> 55 -> 7_du -> L
    feats[464] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '55') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '55', '7_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_du', node2.lemma)
    
    # L -> 55 -> -39 -> L
    feats[465] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '55') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '55', '-39') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-39', node2.lemma)
    
    # L -> 55 -> -33 -> L
    feats[466] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '55') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '55', '-33') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-33', node2.lemma)
    
    # L -> 55 -> 168 -> L
    feats[467] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '55') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '55', '168') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '168', node2.lemma)
    
    # L -> 55 -> 180 -> L
    feats[468] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '55') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '55', '180') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '180', node2.lemma)
    
    # L -> 55 -> 16_sg -> L
    feats[469] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '55') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '55', '16_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_sg', node2.lemma)
    
    # L -> 16_fp -> 5_sg -> L
    feats[470] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_fp', '5_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sg', node2.lemma)
    
    # L -> 16_fp -> 12_pl -> L
    feats[471] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_fp', '12_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_pl', node2.lemma)
    
    # L -> 16_fp -> 6_sg -> L
    feats[472] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_fp', '6_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_sg', node2.lemma)
    
    # L -> 16_fp -> -28 -> L
    feats[473] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_fp', '-28') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-28', node2.lemma)
    
    # L -> 16_fp -> tp -> L
    feats[474] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_fp', 'tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'tp', node2.lemma)
    
    # L -> 16_fp -> -98 -> L
    feats[475] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_fp', '-98') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-98', node2.lemma)
    
    # L -> 16_fp -> -18 -> L
    feats[476] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_fp', '-18') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-18', node2.lemma)
    
    # L -> 16_fp -> 5_fp -> L
    feats[477] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_fp', '5_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_fp', node2.lemma)
    
    # L -> 16_fp -> -126 -> L
    feats[478] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_fp', '-126') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-126', node2.lemma)
    
    # L -> 16_fp -> 135 -> L
    feats[479] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_fp', '135') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '135', node2.lemma)
    
    # L -> 16_fp -> -273 -> L
    feats[480] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_fp', '-273') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-273', node2.lemma)
    
    # L -> 16_fp -> -13 -> L
    feats[481] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_fp', '-13') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-13', node2.lemma)
    
    # L -> 16_fp -> 117 -> L
    feats[482] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_fp', '117') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '117', node2.lemma)
    
    # L -> 16_fp -> -297 -> L
    feats[483] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_fp', '-297') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-297', node2.lemma)
    
    # L -> 16_fp -> 9_tp -> L
    feats[484] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_fp', '9_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_tp', node2.lemma)
    
    # L -> 16_fp -> -84 -> L
    feats[485] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_fp', '-84') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-84', node2.lemma)
    
    # L -> 16_fp -> 76 -> L
    feats[486] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_fp', '76') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '76', node2.lemma)
    
    # L -> tp -> -279 -> L
    feats[487] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '-279') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-279', node2.lemma)
    
    # L -> tp -> 8_pl -> L
    feats[488] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '8_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_pl', node2.lemma)
    
    # L -> tp -> -61 -> L
    feats[489] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '-61') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-61', node2.lemma)
    
    # L -> tp -> -113 -> L
    feats[490] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '-113') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-113', node2.lemma)
    
    # L -> tp -> gen. pl. -> L
    feats[491] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', 'gen. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen. pl.', node2.lemma)
    
    # L -> tp -> du -> L
    feats[492] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', 'du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du', node2.lemma)
    
    # L -> tp -> nom -> L
    feats[493] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', 'nom') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom', node2.lemma)
    
    # L -> tp -> nom. du. -> L
    feats[494] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', 'nom. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. du.', node2.lemma)
    
    # L -> tp -> -291 -> L
    feats[495] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '-291') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-291', node2.lemma)
    
    # L -> tp -> -166 -> L
    feats[496] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '-166') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-166', node2.lemma)
    
    # L -> tp -> 148 -> L
    feats[497] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '148') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '148', node2.lemma)
    
    # L -> tp -> 101 -> L
    feats[498] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '101') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    
    # L -> tp -> 109 -> L
    feats[499] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '109') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '109', node2.lemma)
    
    # L -> tp -> 170 -> L
    feats[500] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '170') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '170', node2.lemma)
    
    # L -> tp -> -42 -> L
    feats[501] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '-42') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-42', node2.lemma)
    
    # L -> tp -> 10_sg -> L
    feats[502] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '10_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_sg', node2.lemma)
    
    # L -> tp -> 12_sp -> L
    feats[503] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '12_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_sp', node2.lemma)
    
    # L -> tp -> -26 -> L
    feats[504] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '-26') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-26', node2.lemma)
    
    # L -> tp -> 80 -> L
    feats[505] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '80') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '80', node2.lemma)
    
    # L -> tp -> voc -> L
    feats[506] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', 'voc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc', node2.lemma)
    
    # L -> tp -> 78 -> L
    feats[507] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '78') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '78', node2.lemma)
    
    # L -> tp -> 131 -> L
    feats[508] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '131') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '131', node2.lemma)
    
    # L -> tp -> 2_tp -> L
    feats[509] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '2_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_tp', node2.lemma)
    
    # L -> tp -> 122 -> L
    feats[510] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '122') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '122', node2.lemma)
    
    # L -> tp -> 176 -> L
    feats[511] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '176') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '176', node2.lemma)
    
    # L -> tp -> -101 -> L
    feats[512] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '-101') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-101', node2.lemma)
    
    # L -> tp -> 70 -> L
    feats[513] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '70') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '70', node2.lemma)
    
    # L -> tp -> voc. neutr. -> L
    feats[514] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', 'voc. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. neutr.', node2.lemma)
    
    # L -> tp -> -115 -> L
    feats[515] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '-115') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-115', node2.lemma)
    
    # L -> tp -> 95 -> L
    feats[516] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '95') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '95', node2.lemma)
    
    # L -> tp -> 35 -> L
    feats[517] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '35') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '35', node2.lemma)
    
    # L -> tp -> acc -> L
    feats[518] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', 'acc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc', node2.lemma)
    
    # L -> tp -> -22 -> L
    feats[519] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '-22') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-22', node2.lemma)
    
    # L -> tp -> -117 -> L
    feats[520] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '-117') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-117', node2.lemma)
    
    # L -> tp -> 15_pl -> L
    feats[521] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '15_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_pl', node2.lemma)
    
    # L -> tp -> -131 -> L
    feats[522] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '-131') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-131', node2.lemma)
    
    # L -> tp -> -35 -> L
    feats[523] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '-35') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    
    # L -> tp -> -58 -> L
    feats[524] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '-58') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-58', node2.lemma)
    
    # L -> tp -> 79 -> L
    feats[525] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '79') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '79', node2.lemma)
    
    # L -> tp -> 7_fp -> L
    feats[526] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '7_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_fp', node2.lemma)
    
    # L -> tp -> 9_pl -> L
    feats[527] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '9_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_pl', node2.lemma)
    
    # L -> tp -> instr. du. -> L
    feats[528] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', 'instr. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. du.', node2.lemma)
    
    # L -> tp -> -82 -> L
    feats[529] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '-82') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-82', node2.lemma)
    
    # L -> tp -> 140 -> L
    feats[530] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '140') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '140', node2.lemma)
    
    # L -> tp -> 2_fp -> L
    feats[531] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '2_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_fp', node2.lemma)
    
    # L -> tp -> -64 -> L
    feats[532] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '-64') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-64', node2.lemma)
    
    # L -> tp -> -10 -> L
    feats[533] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '-10') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-10', node2.lemma)
    
    # L -> tp -> -123 -> L
    feats[534] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '-123') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-123', node2.lemma)
    
    # L -> tp -> 8_tp -> L
    feats[535] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '8_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_tp', node2.lemma)
    
    # L -> tp -> 7_du -> L
    feats[536] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '7_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_du', node2.lemma)
    
    # L -> tp -> 141 -> L
    feats[537] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '141') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '141', node2.lemma)
    
    # L -> -98 -> 119 -> L
    feats[538] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-98', '119') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '119', node2.lemma)
    
    # L -> -98 -> 154 -> L
    feats[539] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-98', '154') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '154', node2.lemma)
    
    # L -> -98 -> 128 -> L
    feats[540] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-98', '128') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '128', node2.lemma)
    
    # L -> -98 -> 39 -> L
    feats[541] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-98', '39') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '39', node2.lemma)
    
    # L -> -98 -> 118 -> L
    feats[542] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-98', '118') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '118', node2.lemma)
    
    # L -> -98 -> -68 -> L
    feats[543] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-98', '-68') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-68', node2.lemma)
    
    # L -> -98 -> 182 -> L
    feats[544] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-98', '182') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '182', node2.lemma)
    
    # L -> -98 -> 93 -> L
    feats[545] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-98', '93') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '93', node2.lemma)
    
    # L -> -98 -> 3_tp -> L
    feats[546] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-98', '3_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_tp', node2.lemma)
    
    # L -> -98 -> 11_tp -> L
    feats[547] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-98', '11_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_tp', node2.lemma)
    
    # L -> -98 -> -154 -> L
    feats[548] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-98', '-154') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-154', node2.lemma)
    
    # L -> -98 -> -41 -> L
    feats[549] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-98', '-41') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-41', node2.lemma)
    
    # L -> -98 -> sg_tp -> L
    feats[550] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-98', 'sg_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_tp', node2.lemma)
    
    # L -> -98 -> -17 -> L
    feats[551] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-98', '-17') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-17', node2.lemma)
    
    # L -> -98 -> -12 -> L
    feats[552] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-98', '-12') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-12', node2.lemma)
    
    # L -> -98 -> -31 -> L
    feats[553] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-98', '-31') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-31', node2.lemma)
    
    # L -> -98 -> 12_sg -> L
    feats[554] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-98', '12_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_sg', node2.lemma)
    
    # L -> -98 -> -62 -> L
    feats[555] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-98', '-62') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-62', node2.lemma)
    
    # L -> -98 -> 10_du -> L
    feats[556] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-98', '10_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_du', node2.lemma)
    
    # L -> -98 -> 27_pl -> L
    feats[557] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-98', '27_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_pl', node2.lemma)
    
    # L -> -98 -> -79 -> L
    feats[558] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-98', '-79') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-79', node2.lemma)
    
    # L -> -98 -> 50 -> L
    feats[559] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-98', '50') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '50', node2.lemma)
    
    # L -> -98 -> nom. pl. -> L
    feats[560] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-98', 'nom. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. pl.', node2.lemma)
    
    # L -> -98 -> -15 -> L
    feats[561] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-98', '-15') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-15', node2.lemma)
    
    # L -> -98 -> -14 -> L
    feats[562] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-98', '-14') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-14', node2.lemma)
    
    # L -> -98 -> 12_tp -> L
    feats[563] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-98', '12_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_tp', node2.lemma)
    
    # L -> -98 -> 10_pl -> L
    feats[564] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-98', '10_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_pl', node2.lemma)
    
    # L -> -98 -> instr. fem -> L
    feats[565] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-98', 'instr. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. fem', node2.lemma)
    
    # L -> -98 -> 120 -> L
    feats[566] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-98', '120') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '120', node2.lemma)
    
    # L -> -98 -> -113 -> L
    feats[567] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-98', '-113') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-113', node2.lemma)
    
    # L -> -98 -> 4_pl -> L
    feats[568] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-98', '4_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_pl', node2.lemma)
    
    # L -> -98 -> -104 -> L
    feats[569] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-98', '-104') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-104', node2.lemma)
    
    # L -> -98 -> 29 -> L
    feats[570] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-98', '29') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '29', node2.lemma)
    
    # L -> -98 -> -153 -> L
    feats[571] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-98', '-153') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-153', node2.lemma)
    
    # L -> -98 -> -166 -> L
    feats[572] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-98', '-166') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-166', node2.lemma)
    
    # L -> -98 -> 148 -> L
    feats[573] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-98', '148') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '148', node2.lemma)
    
    # L -> -98 -> sg_fp -> L
    feats[574] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-98', 'sg_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_fp', node2.lemma)
    
    # L -> -98 -> 142 -> L
    feats[575] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-98', '142') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '142', node2.lemma)
    
    # L -> -98 -> 170 -> L
    feats[576] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-98', '170') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '170', node2.lemma)
    
    # L -> -98 -> sg_sp -> L
    feats[577] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-98', 'sg_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_sp', node2.lemma)
    
    # L -> -98 -> 5_pl -> L
    feats[578] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-98', '5_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_pl', node2.lemma)
    
    # L -> -98 -> -49 -> L
    feats[579] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-98', '-49') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-49', node2.lemma)
    
    # L -> -18 -> 15_sg -> L
    feats[580] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '15_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_sg', node2.lemma)
    
    # L -> -18 -> 41 -> L
    feats[581] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '41') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '41', node2.lemma)
    
    # L -> -18 -> 50 -> L
    feats[582] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '50') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '50', node2.lemma)
    
    # L -> -18 -> instr. fem -> L
    feats[583] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', 'instr. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. fem', node2.lemma)
    
    # L -> -18 -> 120 -> L
    feats[584] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '120') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '120', node2.lemma)
    
    # L -> -18 -> -113 -> L
    feats[585] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '-113') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-113', node2.lemma)
    
    # L -> -18 -> -292 -> L
    feats[586] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '-292') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-292', node2.lemma)
    
    # L -> -18 -> 116 -> L
    feats[587] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '116') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '116', node2.lemma)
    
    # L -> -18 -> -53 -> L
    feats[588] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '-53') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-53', node2.lemma)
    
    # L -> -18 -> 33 -> L
    feats[589] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '33') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '33', node2.lemma)
    
    # L -> -18 -> 101 -> L
    feats[590] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '101') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    
    # L -> -18 -> voc -> L
    feats[591] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', 'voc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc', node2.lemma)
    
    # L -> -18 -> 78 -> L
    feats[592] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '78') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '78', node2.lemma)
    
    # L -> -18 -> -246 -> L
    feats[593] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '-246') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-246', node2.lemma)
    
    # L -> -18 -> acc. neutr. -> L
    feats[594] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', 'acc. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. neutr.', node2.lemma)
    
    # L -> -18 -> -137 -> L
    feats[595] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '-137') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-137', node2.lemma)
    
    # L -> -18 -> 122 -> L
    feats[596] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '122') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '122', node2.lemma)
    
    # L -> -18 -> acc. adj. -> L
    feats[597] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', 'acc. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. adj.', node2.lemma)
    
    # L -> -18 -> -77 -> L
    feats[598] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '-77') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-77', node2.lemma)
    
    # L -> -18 -> -117 -> L
    feats[599] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '-117') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-117', node2.lemma)
    
    # L -> -18 -> 15_fp -> L
    feats[600] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '15_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_fp', node2.lemma)
    
    # L -> -18 -> 15_pl -> L
    feats[601] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '15_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_pl', node2.lemma)
    
    # L -> -18 -> 6_du -> L
    feats[602] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '6_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_du', node2.lemma)
    
    # L -> -18 -> -131 -> L
    feats[603] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '-131') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-131', node2.lemma)
    
    # L -> -18 -> abl. du. -> L
    feats[604] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', 'abl. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. du.', node2.lemma)
    
    # L -> -18 -> 9_pl -> L
    feats[605] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '9_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_pl', node2.lemma)
    
    # L -> -18 -> 9_du -> L
    feats[606] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '9_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_du', node2.lemma)
    
    # L -> -18 -> 9_sg -> L
    feats[607] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '9_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_sg', node2.lemma)
    
    # L -> -18 -> -87 -> L
    feats[608] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '-87') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-87', node2.lemma)
    
    # L -> -18 -> 140 -> L
    feats[609] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '140') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '140', node2.lemma)
    
    # L -> -18 -> 6_pl -> L
    feats[610] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '6_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_pl', node2.lemma)
    
    # L -> -18 -> -169 -> L
    feats[611] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', '-169') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-169', node2.lemma)
    
    # L -> 16_tp -> 37 -> L
    feats[612] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', '37') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '37', node2.lemma)
    
    # L -> 16_tp -> -242 -> L
    feats[613] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', '-242') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-242', node2.lemma)
    
    # L -> 16_tp -> 54 -> L
    feats[614] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', '54') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    
    # L -> 16_tp -> -111 -> L
    feats[615] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', '-111') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-111', node2.lemma)
    
    # L -> 16_tp -> 3_sg -> L
    feats[616] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', '3_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_sg', node2.lemma)
    
    # L -> 16_tp -> -68 -> L
    feats[617] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', '-68') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-68', node2.lemma)
    
    # L -> 16_tp -> 69 -> L
    feats[618] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', '69') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '69', node2.lemma)
    
    # L -> 16_tp -> 4_tp -> L
    feats[619] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', '4_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_tp', node2.lemma)
    
    # L -> 16_tp -> -23 -> L
    feats[620] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', '-23') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-23', node2.lemma)
    
    # L -> 16_tp -> -17 -> L
    feats[621] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', '-17') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-17', node2.lemma)
    
    # L -> 16_tp -> acc. fem -> L
    feats[622] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', 'acc. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. fem', node2.lemma)
    
    # L -> 16_tp -> 12_sg -> L
    feats[623] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', '12_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_sg', node2.lemma)
    
    # L -> 16_tp -> -93 -> L
    feats[624] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', '-93') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-93', node2.lemma)
    
    # L -> 16_tp -> instr -> L
    feats[625] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', 'instr') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr', node2.lemma)
    
    # L -> 16_tp -> 112 -> L
    feats[626] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', '112') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '112', node2.lemma)
    
    # L -> 16_tp -> 14_fp -> L
    feats[627] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', '14_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_fp', node2.lemma)
    
    # L -> 16_tp -> 50 -> L
    feats[628] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', '50') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '50', node2.lemma)
    
    # L -> 16_tp -> -90 -> L
    feats[629] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', '-90') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-90', node2.lemma)
    
    # L -> 16_tp -> 77 -> L
    feats[630] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', '77') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '77', node2.lemma)
    
    # L -> 16_tp -> -15 -> L
    feats[631] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', '-15') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-15', node2.lemma)
    
    # L -> 16_tp -> 10_pl -> L
    feats[632] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', '10_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_pl', node2.lemma)
    
    # L -> 16_tp -> instr. fem -> L
    feats[633] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', 'instr. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. fem', node2.lemma)
    
    # L -> 16_tp -> -249 -> L
    feats[634] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', '-249') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-249', node2.lemma)
    
    # L -> 16_tp -> -113 -> L
    feats[635] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', '-113') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-113', node2.lemma)
    
    # L -> 16_tp -> -200 -> L
    feats[636] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', '-200') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-200', node2.lemma)
    
    # L -> 16_tp -> -245 -> L
    feats[637] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', '-245') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-245', node2.lemma)
    
    # L -> 16_tp -> 97 -> L
    feats[638] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', '97') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '97', node2.lemma)
    
    # L -> 16_tp -> 71 -> L
    feats[639] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', '71') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '71', node2.lemma)
    
    # L -> 16_tp -> nom. du. -> L
    feats[640] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', 'nom. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. du.', node2.lemma)
    
    # L -> 16_tp -> 14_sp -> L
    feats[641] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', '14_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_sp', node2.lemma)
    
    # L -> 16_tp -> gen -> L
    feats[642] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', 'gen') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen', node2.lemma)
    
    # L -> 16_tp -> sg_fp -> L
    feats[643] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', 'sg_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_fp', node2.lemma)
    
    # L -> 16_tp -> 115 -> L
    feats[644] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', '115') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '115', node2.lemma)
    
    # L -> 16_tp -> 15_du -> L
    feats[645] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', '15_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_du', node2.lemma)
    
    # L -> 16_tp -> 137 -> L
    feats[646] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', '137') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    
    # L -> 16_tp -> 176 -> L
    feats[647] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', '176') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '176', node2.lemma)
    
    # L -> 16_tp -> 70 -> L
    feats[648] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', '70') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '70', node2.lemma)
    
    # L -> 16_tp -> 95 -> L
    feats[649] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', '95') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '95', node2.lemma)
    
    # L -> 16_tp -> -11 -> L
    feats[650] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', '-11') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-11', node2.lemma)
    
    # L -> 16_tp -> -22 -> L
    feats[651] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', '-22') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-22', node2.lemma)
    
    # L -> 16_tp -> 15_pl -> L
    feats[652] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', '15_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_pl', node2.lemma)
    
    # L -> 16_tp -> -271 -> L
    feats[653] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', '-271') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-271', node2.lemma)
    
    # L -> 16_tp -> -46 -> L
    feats[654] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', '-46') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-46', node2.lemma)
    
    # L -> 16_tp -> 7_fp -> L
    feats[655] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', '7_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_fp', node2.lemma)
    
    # L -> 16_tp -> -82 -> L
    feats[656] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', '-82') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-82', node2.lemma)
    
    # L -> 16_tp -> -157 -> L
    feats[657] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', '-157') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-157', node2.lemma)
    
    # L -> 16_tp -> -87 -> L
    feats[658] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', '-87') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-87', node2.lemma)
    
    # L -> 16_tp -> -10 -> L
    feats[659] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', '-10') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-10', node2.lemma)
    
    # L -> 16_tp -> 8_sg -> L
    feats[660] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', '8_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_sg', node2.lemma)
    
    # L -> 16_tp -> 98 -> L
    feats[661] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', '98') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '98', node2.lemma)
    
    # L -> 16_tp -> 94 -> L
    feats[662] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', '94') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '94', node2.lemma)
    
    # L -> 16_tp -> dat -> L
    feats[663] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', 'dat') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat', node2.lemma)
    
    # L -> 16_tp -> 60 -> L
    feats[664] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', '60') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '60', node2.lemma)
    
    # L -> 16_tp -> nom. sg. -> L
    feats[665] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', 'nom. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. sg.', node2.lemma)
    
    # L -> 16_tp -> 51 -> L
    feats[666] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', '51') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '51', node2.lemma)
    
    # L -> 16_tp -> 178 -> L
    feats[667] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', '178') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '178', node2.lemma)
    
    # L -> 5_fp -> 14_pl -> L
    feats[668] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', '14_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_pl', node2.lemma)
    
    # L -> 5_fp -> masc -> L
    feats[669] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', 'masc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'masc', node2.lemma)
    
    # L -> 5_fp -> 36 -> L
    feats[670] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', '36') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '36', node2.lemma)
    
    # L -> 5_fp -> -66 -> L
    feats[671] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', '-66') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-66', node2.lemma)
    
    # L -> 5_fp -> -141 -> L
    feats[672] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', '-141') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-141', node2.lemma)
    
    # L -> 5_fp -> instr. adj. -> L
    feats[673] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', 'instr. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. adj.', node2.lemma)
    
    # L -> 5_fp -> 11_sg -> L
    feats[674] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', '11_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_sg', node2.lemma)
    
    # L -> 5_fp -> -78 -> L
    feats[675] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', '-78') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-78', node2.lemma)
    
    # L -> 5_fp -> -308 -> L
    feats[676] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', '-308') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-308', node2.lemma)
    
    # L -> 5_fp -> -147 -> L
    feats[677] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', '-147') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-147', node2.lemma)
    
    # L -> 5_fp -> -122 -> L
    feats[678] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', '-122') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-122', node2.lemma)
    
    # L -> 5_fp -> 117 -> L
    feats[679] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', '117') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '117', node2.lemma)
    
    # L -> 5_fp -> -84 -> L
    feats[680] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', '-84') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-84', node2.lemma)
    
    # L -> 5_fp -> -52 -> L
    feats[681] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', '-52') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-52', node2.lemma)
    
    # L -> 5_fp -> 29_sg -> L
    feats[682] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', '29_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '29_sg', node2.lemma)
    
    # L -> 5_fp -> 179 -> L
    feats[683] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', '179') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '179', node2.lemma)
    
    # L -> 5_fp -> 119 -> L
    feats[684] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', '119') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '119', node2.lemma)
    
    # L -> 5_fp -> 3_sg -> L
    feats[685] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', '3_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_sg', node2.lemma)
    
    # L -> 5_fp -> 1 -> L
    feats[686] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', '1') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '1', node2.lemma)
    
    # L -> 5_fp -> 39 -> L
    feats[687] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', '39') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '39', node2.lemma)
    
    # L -> 5_fp -> 3_tp -> L
    feats[688] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', '3_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_tp', node2.lemma)
    
    # L -> 5_fp -> 174 -> L
    feats[689] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', '174') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '174', node2.lemma)
    
    # L -> 5_fp -> 4_tp -> L
    feats[690] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', '4_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_tp', node2.lemma)
    
    # L -> 5_fp -> 88 -> L
    feats[691] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', '88') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '88', node2.lemma)
    
    # L -> 5_fp -> 156 -> L
    feats[692] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', '156') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '156', node2.lemma)
    
    # L -> 5_fp -> 112 -> L
    feats[693] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', '112') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '112', node2.lemma)
    
    # L -> 5_fp -> acc. masc. -> L
    feats[694] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', 'acc. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. masc.', node2.lemma)
    
    # L -> 5_fp -> 10_du -> L
    feats[695] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', '10_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_du', node2.lemma)
    
    # L -> 5_fp -> 41 -> L
    feats[696] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', '41') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '41', node2.lemma)
    
    # L -> 5_fp -> voc. pl. -> L
    feats[697] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', 'voc. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. pl.', node2.lemma)
    
    # L -> 5_fp -> 13_sg -> L
    feats[698] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', '13_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_sg', node2.lemma)
    
    # L -> 5_fp -> 120 -> L
    feats[699] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', '120') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '120', node2.lemma)
    
    # L -> 5_fp -> -163 -> L
    feats[700] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', '-163') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-163', node2.lemma)
    
    # L -> 5_fp -> du -> L
    feats[701] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', 'du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du', node2.lemma)
    
    # L -> 5_fp -> -153 -> L
    feats[702] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', '-153') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-153', node2.lemma)
    
    # L -> 5_fp -> 33 -> L
    feats[703] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', '33') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '33', node2.lemma)
    
    # L -> 5_fp -> 68 -> L
    feats[704] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', '68') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '68', node2.lemma)
    
    # L -> 5_fp -> 158 -> L
    feats[705] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', '158') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '158', node2.lemma)
    
    # L -> 5_fp -> 14_sp -> L
    feats[706] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', '14_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_sp', node2.lemma)
    
    # L -> 5_fp -> -230 -> L
    feats[707] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', '-230') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-230', node2.lemma)
    
    # L -> 5_fp -> -129 -> L
    feats[708] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', '-129') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-129', node2.lemma)
    
    # L -> 5_fp -> 176 -> L
    feats[709] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', '176') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '176', node2.lemma)
    
    # L -> 5_fp -> 99 -> L
    feats[710] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', '99') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '99', node2.lemma)
    
    # L -> 5_fp -> -11 -> L
    feats[711] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', '-11') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-11', node2.lemma)
    
    # L -> -83 -> 100 -> L
    feats[712] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-83') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-83', '100') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '100', node2.lemma)
    
    # L -> -83 -> 7_du -> L
    feats[713] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-83') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-83', '7_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_du', node2.lemma)
    
    # L -> -83 -> -260 -> L
    feats[714] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-83') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-83', '-260') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-260', node2.lemma)
    
    # L -> -83 -> 60 -> L
    feats[715] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-83') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-83', '60') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '60', node2.lemma)
    
    # L -> 8_fp -> 14_pl -> L
    feats[716] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_fp', '14_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_pl', node2.lemma)
    
    # L -> 8_fp -> -243 -> L
    feats[717] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_fp', '-243') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-243', node2.lemma)
    
    # L -> 8_fp -> 16_fp -> L
    feats[718] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_fp', '16_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_fp', node2.lemma)
    
    # L -> 8_fp -> -16 -> L
    feats[719] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_fp', '-16') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-16', node2.lemma)
    
    # L -> 8_fp -> 7_pl -> L
    feats[720] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_fp', '7_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_pl', node2.lemma)
    
    # L -> 8_fp -> -302 -> L
    feats[721] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_fp', '-302') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-302', node2.lemma)
    
    # L -> 8_fp -> 81 -> L
    feats[722] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_fp', '81') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '81', node2.lemma)
    
    # L -> 8_fp -> -44 -> L
    feats[723] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_fp', '-44') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-44', node2.lemma)
    
    # L -> 8_fp -> -122 -> L
    feats[724] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_fp', '-122') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-122', node2.lemma)
    
    # L -> 8_fp -> 117 -> L
    feats[725] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_fp', '117') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '117', node2.lemma)
    
    # L -> 8_fp -> -84 -> L
    feats[726] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_fp', '-84') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-84', node2.lemma)
    
    # L -> 8_fp -> 153 -> L
    feats[727] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_fp', '153') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    
    # L -> 8_fp -> 14_tp -> L
    feats[728] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_fp', '14_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_tp', node2.lemma)
    
    # L -> -16 -> 70 -> L
    feats[729] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-16') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-16', '70') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '70', node2.lemma)
    
    # L -> -16 -> 95 -> L
    feats[730] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-16') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-16', '95') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '95', node2.lemma)
    
    # L -> -16 -> -77 -> L
    feats[731] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-16') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-16', '-77') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-77', node2.lemma)
    
    # L -> -16 -> acc -> L
    feats[732] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-16') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-16', 'acc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc', node2.lemma)
    
    # L -> -16 -> 7_tp -> L
    feats[733] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-16') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-16', '7_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_tp', node2.lemma)
    
    # L -> -16 -> -143 -> L
    feats[734] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-16') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-16', '-143') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-143', node2.lemma)
    
    # L -> -16 -> fem -> L
    feats[735] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-16') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-16', 'fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'fem', node2.lemma)
    
    # L -> -16 -> 149 -> L
    feats[736] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-16') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-16', '149') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '149', node2.lemma)
    
    # L -> -16 -> -82 -> L
    feats[737] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-16') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-16', '-82') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-82', node2.lemma)
    
    # L -> -16 -> 140 -> L
    feats[738] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-16') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-16', '140') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '140', node2.lemma)
    
    # L -> -16 -> -133 -> L
    feats[739] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-16') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-16', '-133') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-133', node2.lemma)
    
    # L -> -16 -> 56 -> L
    feats[740] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-16') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-16', '56') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '56', node2.lemma)
    
    # L -> -16 -> 8_tp -> L
    feats[741] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-16') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-16', '8_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_tp', node2.lemma)
    
    # L -> -16 -> 7_du -> L
    feats[742] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-16') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-16', '7_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_du', node2.lemma)
    
    # L -> -16 -> 82 -> L
    feats[743] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-16') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-16', '82') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '82', node2.lemma)
    
    # L -> -16 -> 160 -> L
    feats[744] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-16') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-16', '160') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    
    # L -> -16 -> -32 -> L
    feats[745] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-16') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-16', '-32') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-32', node2.lemma)
    
    # L -> -16 -> -73 -> L
    feats[746] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-16') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-16', '-73') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-73', node2.lemma)
    
    # L -> -16 -> nom. sg. -> L
    feats[747] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-16') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-16', 'nom. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. sg.', node2.lemma)
    
    # L -> -16 -> acc. pl. -> L
    feats[748] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-16') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-16', 'acc. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. pl.', node2.lemma)
    
    # L -> -16 -> 168 -> L
    feats[749] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-16') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-16', '168') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '168', node2.lemma)
    
    # L -> -16 -> 51 -> L
    feats[750] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-16') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-16', '51') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '51', node2.lemma)
    
    # L -> 28_sg -> voc. masc. -> L
    feats[751] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28_sg', 'voc. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. masc.', node2.lemma)
    
    # L -> 28_sg -> 2_sp -> L
    feats[752] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28_sg', '2_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_sp', node2.lemma)
    
    # L -> 28_sg -> -152 -> L
    feats[753] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28_sg', '-152') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-152', node2.lemma)
    
    # L -> 28_sg -> nom. masc. -> L
    feats[754] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28_sg', 'nom. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. masc.', node2.lemma)
    
    # L -> 28_sg -> du_sp -> L
    feats[755] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28_sg', 'du_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_sp', node2.lemma)
    
    # L -> 28_sg -> -262 -> L
    feats[756] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28_sg', '-262') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-262', node2.lemma)
    
    # L -> 28_sg -> 16_fp -> L
    feats[757] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28_sg', '16_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_fp', node2.lemma)
    
    # L -> 28_sg -> 16_tp -> L
    feats[758] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28_sg', '16_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_tp', node2.lemma)
    
    # L -> 28_sg -> -67 -> L
    feats[759] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28_sg', '-67') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-67', node2.lemma)
    
    # L -> 28_sg -> 59 -> L
    feats[760] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28_sg', '59') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '59', node2.lemma)
    
    # L -> 28_sg -> -276 -> L
    feats[761] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28_sg', '-276') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-276', node2.lemma)
    
    # L -> 28_sg -> -44 -> L
    feats[762] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28_sg', '-44') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-44', node2.lemma)
    
    # L -> -78 -> masc -> L
    feats[763] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-78') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-78', 'masc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'masc', node2.lemma)
    
    # L -> -78 -> 30_du -> L
    feats[764] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-78') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-78', '30_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_du', node2.lemma)
    
    # L -> -78 -> -243 -> L
    feats[765] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-78') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-78', '-243') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-243', node2.lemma)
    
    # L -> -78 -> 6_sg -> L
    feats[766] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-78') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-78', '6_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_sg', node2.lemma)
    
    # L -> -78 -> 59 -> L
    feats[767] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-78') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-78', '59') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '59', node2.lemma)
    
    # L -> -78 -> abl. sg. -> L
    feats[768] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-78') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-78', 'abl. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. sg.', node2.lemma)
    
    # L -> -78 -> 3_du -> L
    feats[769] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-78') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-78', '3_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    
    # L -> -78 -> -139 -> L
    feats[770] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-78') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-78', '-139') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-139', node2.lemma)
    
    # L -> -78 -> -44 -> L
    feats[771] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-78') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-78', '-44') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-44', node2.lemma)
    
    # L -> -78 -> -13 -> L
    feats[772] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-78') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-78', '-13') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-13', node2.lemma)
    
    # L -> -78 -> 155 -> L
    feats[773] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-78') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-78', '155') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '155', node2.lemma)
    
    # L -> -78 -> 37 -> L
    feats[774] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-78') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-78', '37') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '37', node2.lemma)
    
    # L -> -78 -> -89 -> L
    feats[775] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-78') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-78', '-89') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-89', node2.lemma)
    
    # L -> -78 -> 74 -> L
    feats[776] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-78') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-78', '74') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '74', node2.lemma)
    
    # L -> -78 -> 1 -> L
    feats[777] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-78') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-78', '1') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '1', node2.lemma)
    
    # L -> -78 -> -220 -> L
    feats[778] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-78') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-78', '-220') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-220', node2.lemma)
    
    # L -> -78 -> 118 -> L
    feats[779] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-78') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-78', '118') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '118', node2.lemma)
    
    # L -> voc. du. -> -97 -> L
    feats[780] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. du.', '-97') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-97', node2.lemma)
    
    # L -> voc. du. -> acc. fem -> L
    feats[781] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. du.', 'acc. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. fem', node2.lemma)
    
    # L -> voc. du. -> neutr -> L
    feats[782] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. du.', 'neutr') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'neutr', node2.lemma)
    
    # L -> voc. du. -> 121 -> L
    feats[783] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. du.', '121') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '121', node2.lemma)
    
    # L -> voc. du. -> -86 -> L
    feats[784] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. du.', '-86') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-86', node2.lemma)
    
    # L -> 27_sg -> loc -> L
    feats[785] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', 'loc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc', node2.lemma)
    
    # L -> 27_sg -> -210 -> L
    feats[786] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', '-210') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-210', node2.lemma)
    
    # L -> 27_sg -> 9_tp -> L
    feats[787] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', '9_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_tp', node2.lemma)
    
    # L -> 27_sg -> 2_pl -> L
    feats[788] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', '2_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_pl', node2.lemma)
    
    # L -> 27_sg -> -34 -> L
    feats[789] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', '-34') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-34', node2.lemma)
    
    # L -> 27_sg -> 154 -> L
    feats[790] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', '154') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '154', node2.lemma)
    
    # L -> 27_sg -> 1 -> L
    feats[791] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', '1') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '1', node2.lemma)
    
    # L -> 27_sg -> 182 -> L
    feats[792] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', '182') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '182', node2.lemma)
    
    # L -> 27_sg -> 6_tp -> L
    feats[793] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', '6_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_tp', node2.lemma)
    
    # L -> 27_sg -> 35 -> L
    feats[794] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', '35') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '35', node2.lemma)
    
    # L -> abl. sg. -> neutr -> L
    feats[795] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'abl. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'abl. sg.', 'neutr') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'neutr', node2.lemma)
    
    # L -> abl. sg. -> 15_sp -> L
    feats[796] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'abl. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'abl. sg.', '15_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_sp', node2.lemma)
    
    # L -> abl. sg. -> du_fp -> L
    feats[797] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'abl. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'abl. sg.', 'du_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_fp', node2.lemma)
    
    # L -> abl. sg. -> 42 -> L
    feats[798] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'abl. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'abl. sg.', '42') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '42', node2.lemma)
    
    # L -> abl. sg. -> -96 -> L
    feats[799] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'abl. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'abl. sg.', '-96') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-96', node2.lemma)
    
    # L -> abl. sg. -> -114 -> L
    feats[800] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'abl. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'abl. sg.', '-114') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-114', node2.lemma)
    
    # L -> abl. sg. -> 27_pl -> L
    feats[801] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'abl. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'abl. sg.', '27_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_pl', node2.lemma)
    
    # L -> abl. sg. -> -79 -> L
    feats[802] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'abl. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'abl. sg.', '-79') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-79', node2.lemma)
    
    # L -> abl. sg. -> -240 -> L
    feats[803] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'abl. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'abl. sg.', '-240') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-240', node2.lemma)
    
    # L -> abl. sg. -> 8_sp -> L
    feats[804] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'abl. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'abl. sg.', '8_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_sp', node2.lemma)
    
    # L -> abl. sg. -> 11_fp -> L
    feats[805] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'abl. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'abl. sg.', '11_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_fp', node2.lemma)
    
    # L -> abl. sg. -> -90 -> L
    feats[806] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'abl. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'abl. sg.', '-90') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-90', node2.lemma)
    
    # L -> abl. sg. -> 13_sg -> L
    feats[807] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'abl. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'abl. sg.', '13_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_sg', node2.lemma)
    
    # L -> abl. sg. -> instr. fem -> L
    feats[808] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'abl. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'abl. sg.', 'instr. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. fem', node2.lemma)
    
    # L -> abl. sg. -> 2 -> L
    feats[809] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'abl. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'abl. sg.', '2') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2', node2.lemma)
    
    # L -> abl. sg. -> -292 -> L
    feats[810] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'abl. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'abl. sg.', '-292') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-292', node2.lemma)
    
    # L -> abl. sg. -> 158 -> L
    feats[811] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'abl. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'abl. sg.', '158') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '158', node2.lemma)
    
    # L -> abl. sg. -> -27 -> L
    feats[812] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'abl. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'abl. sg.', '-27') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-27', node2.lemma)
    
    # L -> abl. sg. -> 6_fp -> L
    feats[813] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'abl. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'abl. sg.', '6_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_fp', node2.lemma)
    
    # L -> abl. sg. -> -49 -> L
    feats[814] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'abl. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'abl. sg.', '-49') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-49', node2.lemma)
    
    # L -> abl. sg. -> 173 -> L
    feats[815] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'abl. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'abl. sg.', '173') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '173', node2.lemma)
    
    # L -> abl. sg. -> -137 -> L
    feats[816] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'abl. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'abl. sg.', '-137') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-137', node2.lemma)
    
    # L -> abl. sg. -> 70 -> L
    feats[817] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'abl. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'abl. sg.', '70') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '70', node2.lemma)
    
    # L -> abl. sg. -> 111 -> L
    feats[818] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'abl. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'abl. sg.', '111') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '111', node2.lemma)
    
    # L -> abl. sg. -> -271 -> L
    feats[819] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'abl. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'abl. sg.', '-271') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-271', node2.lemma)
    
    # L -> abl. sg. -> 7_tp -> L
    feats[820] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'abl. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'abl. sg.', '7_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_tp', node2.lemma)
    
    # L -> abl. sg. -> -58 -> L
    feats[821] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'abl. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'abl. sg.', '-58') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-58', node2.lemma)
    
    # L -> abl. sg. -> 7_fp -> L
    feats[822] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'abl. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'abl. sg.', '7_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_fp', node2.lemma)
    
    # L -> abl. sg. -> -82 -> L
    feats[823] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'abl. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'abl. sg.', '-82') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-82', node2.lemma)
    
    # L -> -43 -> 34 -> L
    feats[824] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-43') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-43', '34') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '34', node2.lemma)
    
    # L -> -43 -> -43 -> L
    feats[825] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-43') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-43', '-43') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-43', node2.lemma)
    
    # L -> -43 -> 5_du -> L
    feats[826] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-43') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-43', '5_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_du', node2.lemma)
    
    # L -> -43 -> -84 -> L
    feats[827] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-43') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-43', '-84') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-84', node2.lemma)
    
    # L -> -43 -> 153 -> L
    feats[828] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-43') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-43', '153') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    
    # L -> -43 -> 174 -> L
    feats[829] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-43') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-43', '174') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '174', node2.lemma)
    
    # L -> -43 -> 4_sp -> L
    feats[830] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-43') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-43', '4_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_sp', node2.lemma)
    
    # L -> -43 -> -142 -> L
    feats[831] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-43') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-43', '-142') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-142', node2.lemma)
    
    # L -> -43 -> -36 -> L
    feats[832] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-43') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-43', '-36') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-36', node2.lemma)
    
    # L -> -43 -> 88 -> L
    feats[833] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-43') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-43', '88') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '88', node2.lemma)
    
    # L -> -126 -> 3_fp -> L
    feats[834] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-126') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-126', '3_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_fp', node2.lemma)
    
    # L -> -126 -> 30_sg -> L
    feats[835] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-126') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-126', '30_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_sg', node2.lemma)
    
    # L -> 31 -> 4_du -> L
    feats[836] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '31') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '31', '4_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_du', node2.lemma)
    
    # L -> 31 -> -50 -> L
    feats[837] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '31') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '31', '-50') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-50', node2.lemma)
    
    # L -> 31 -> 5_sg -> L
    feats[838] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '31') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '31', '5_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sg', node2.lemma)
    
    # L -> 31 -> loc -> L
    feats[839] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '31') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '31', 'loc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc', node2.lemma)
    
    # L -> 31 -> 28_sg -> L
    feats[840] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '31') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '31', '28_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28_sg', node2.lemma)
    
    # L -> 31 -> -308 -> L
    feats[841] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '31') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '31', '-308') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-308', node2.lemma)
    
    # L -> 31 -> 135 -> L
    feats[842] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '31') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '31', '135') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '135', node2.lemma)
    
    # L -> 31 -> 13_fp -> L
    feats[843] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '31') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '31', '13_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_fp', node2.lemma)
    
    # L -> 31 -> 9_tp -> L
    feats[844] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '31') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '31', '9_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_tp', node2.lemma)
    
    # L -> 31 -> 138 -> L
    feats[845] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '31') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '31', '138') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '138', node2.lemma)
    
    # L -> 81 -> -299 -> L
    feats[846] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '81', '-299') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-299', node2.lemma)
    
    # L -> 81 -> -157 -> L
    feats[847] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '81', '-157') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-157', node2.lemma)
    
    # L -> 81 -> -119 -> L
    feats[848] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '81', '-119') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-119', node2.lemma)
    
    # L -> 81 -> 56 -> L
    feats[849] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '81', '56') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '56', node2.lemma)
    
    # L -> 81 -> -123 -> L
    feats[850] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '81', '-123') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-123', node2.lemma)
    
    # L -> 81 -> -21 -> L
    feats[851] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '81', '-21') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-21', node2.lemma)
    
    # L -> 81 -> -71 -> L
    feats[852] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '81', '-71') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-71', node2.lemma)
    
    # L -> 81 -> 160 -> L
    feats[853] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '81', '160') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    
    # L -> 81 -> nom. adj. -> L
    feats[854] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '81', 'nom. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. adj.', node2.lemma)
    
    # L -> 81 -> -149 -> L
    feats[855] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '81', '-149') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-149', node2.lemma)
    
    # L -> 81 -> -32 -> L
    feats[856] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '81', '-32') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-32', node2.lemma)
    
    # L -> 81 -> 152 -> L
    feats[857] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '81', '152') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '152', node2.lemma)
    
    # L -> 81 -> 60 -> L
    feats[858] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '81', '60') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '60', node2.lemma)
    
    # L -> 81 -> nom. sg. -> L
    feats[859] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '81', 'nom. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. sg.', node2.lemma)
    
    # L -> 81 -> 75 -> L
    feats[860] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '81', '75') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '75', node2.lemma)
    
    # L -> 81 -> -94 -> L
    feats[861] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '81', '-94') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-94', node2.lemma)
    
    # L -> 81 -> instr. pl. -> L
    feats[862] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '81', 'instr. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. pl.', node2.lemma)
    
    # L -> 81 -> 11_pl -> L
    feats[863] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '81', '11_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_pl', node2.lemma)
    
    # L -> 81 -> -190 -> L
    feats[864] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '81', '-190') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-190', node2.lemma)
    
    # L -> -29 -> -268 -> L
    feats[865] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-29') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-29', '-268') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-268', node2.lemma)
    
    # L -> -29 -> -293 -> L
    feats[866] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-29') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-29', '-293') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-293', node2.lemma)
    
    # L -> -29 -> 10_sp -> L
    feats[867] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-29') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-29', '10_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_sp', node2.lemma)
    
    # L -> -29 -> 13_pl -> L
    feats[868] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-29') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-29', '13_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_pl', node2.lemma)
    
    # L -> -29 -> 132 -> L
    feats[869] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-29') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-29', '132') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '132', node2.lemma)
    
    # L -> -29 -> 3_sp -> L
    feats[870] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-29') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-29', '3_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_sp', node2.lemma)
    
    # L -> -29 -> -141 -> L
    feats[871] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-29') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-29', '-141') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-141', node2.lemma)
    
    # L -> -29 -> -69 -> L
    feats[872] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-29') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-29', '-69') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-69', node2.lemma)
    
    # L -> -29 -> 30_du -> L
    feats[873] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-29') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-29', '30_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_du', node2.lemma)
    
    # L -> -29 -> 6_sg -> L
    feats[874] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-29') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-29', '6_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_sg', node2.lemma)
    
    # L -> -29 -> -28 -> L
    feats[875] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-29') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-29', '-28') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-28', node2.lemma)
    
    # L -> -29 -> 55 -> L
    feats[876] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-29') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-29', '55') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '55', node2.lemma)
    
    # L -> -29 -> tp -> L
    feats[877] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-29') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-29', 'tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'tp', node2.lemma)
    
    # L -> -29 -> -98 -> L
    feats[878] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-29') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-29', '-98') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-98', node2.lemma)
    
    # L -> -29 -> 11_sg -> L
    feats[879] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-29') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-29', '11_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_sg', node2.lemma)
    
    # L -> -29 -> 27_tp -> L
    feats[880] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-29') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-29', '27_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_tp', node2.lemma)
    
    # L -> -29 -> 59 -> L
    feats[881] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-29') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-29', '59') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '59', node2.lemma)
    
    # L -> -29 -> -276 -> L
    feats[882] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-29') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-29', '-276') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-276', node2.lemma)
    
    # L -> -29 -> abl. pl. -> L
    feats[883] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-29') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-29', 'abl. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. pl.', node2.lemma)
    
    # L -> -29 -> -29 -> L
    feats[884] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-29') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-29', '-29') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-29', node2.lemma)
    
    # L -> -29 -> du_tp -> L
    feats[885] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-29') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-29', 'du_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_tp', node2.lemma)
    
    # L -> -29 -> -273 -> L
    feats[886] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-29') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-29', '-273') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-273', node2.lemma)
    
    # L -> -29 -> loc. sg. -> L
    feats[887] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-29') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-29', 'loc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc. sg.', node2.lemma)
    
    # L -> -29 -> loc. pl. -> L
    feats[888] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-29') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-29', 'loc. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc. pl.', node2.lemma)
    
    # L -> -29 -> -122 -> L
    feats[889] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-29') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-29', '-122') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-122', node2.lemma)
    
    # L -> -29 -> 13_fp -> L
    feats[890] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-29') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-29', '13_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_fp', node2.lemma)
    
    # L -> -29 -> 14_sg -> L
    feats[891] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-29') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-29', '14_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_sg', node2.lemma)
    
    # L -> -29 -> 138 -> L
    feats[892] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-29') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-29', '138') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '138', node2.lemma)
    
    # L -> -29 -> 5_du -> L
    feats[893] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-29') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-29', '5_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_du', node2.lemma)
    
    # L -> -29 -> -242 -> L
    feats[894] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-29') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-29', '-242') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-242', node2.lemma)
    
    # L -> -29 -> -25 -> L
    feats[895] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-29') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-29', '-25') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-25', node2.lemma)
    
    # L -> -29 -> 153 -> L
    feats[896] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-29') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-29', '153') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    
    # L -> -29 -> 2_pl -> L
    feats[897] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-29') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-29', '2_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_pl', node2.lemma)
    
    # L -> du_tp -> -96 -> L
    feats[898] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', '-96') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-96', node2.lemma)
    
    # L -> du_tp -> -62 -> L
    feats[899] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', '-62') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-62', node2.lemma)
    
    # L -> du_tp -> 27_pl -> L
    feats[900] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', '27_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_pl', node2.lemma)
    
    # L -> du_tp -> -79 -> L
    feats[901] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', '-79') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-79', node2.lemma)
    
    # L -> du_tp -> 41 -> L
    feats[902] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', '41') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '41', node2.lemma)
    
    # L -> du_tp -> 150 -> L
    feats[903] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', '150') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '150', node2.lemma)
    
    # L -> du_tp -> 30_tp -> L
    feats[904] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', '30_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_tp', node2.lemma)
    
    # L -> du_tp -> 14_fp -> L
    feats[905] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', '14_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_fp', node2.lemma)
    
    # L -> du_tp -> -112 -> L
    feats[906] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', '-112') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    
    # L -> du_tp -> nom. pl. -> L
    feats[907] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', 'nom. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. pl.', node2.lemma)
    
    # L -> du_tp -> -15 -> L
    feats[908] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', '-15') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-15', node2.lemma)
    
    # L -> du_tp -> -279 -> L
    feats[909] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', '-279') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-279', node2.lemma)
    
    # L -> du_tp -> 12_tp -> L
    feats[910] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', '12_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_tp', node2.lemma)
    
    # L -> du_tp -> 181 -> L
    feats[911] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', '181') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '181', node2.lemma)
    
    # L -> du_tp -> -249 -> L
    feats[912] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', '-249') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-249', node2.lemma)
    
    # L -> du_tp -> -200 -> L
    feats[913] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', '-200') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-200', node2.lemma)
    
    # L -> du_tp -> -245 -> L
    feats[914] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', '-245') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-245', node2.lemma)
    
    # L -> du_tp -> gen. pl. -> L
    feats[915] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', 'gen. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen. pl.', node2.lemma)
    
    # L -> du_tp -> -57 -> L
    feats[916] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', '-57') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-57', node2.lemma)
    
    # L -> du_tp -> -306 -> L
    feats[917] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', '-306') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-306', node2.lemma)
    
    # L -> du_tp -> -38 -> L
    feats[918] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', '-38') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-38', node2.lemma)
    
    # L -> du_tp -> -24 -> L
    feats[919] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', '-24') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    
    # L -> du_tp -> 116 -> L
    feats[920] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', '116') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '116', node2.lemma)
    
    # L -> du_tp -> -263 -> L
    feats[921] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', '-263') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-263', node2.lemma)
    
    # L -> du_tp -> du -> L
    feats[922] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', 'du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du', node2.lemma)
    
    # L -> du_tp -> loc. du. -> L
    feats[923] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', 'loc. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc. du.', node2.lemma)
    
    # L -> du_tp -> -153 -> L
    feats[924] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', '-153') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-153', node2.lemma)
    
    # L -> du_tp -> 33 -> L
    feats[925] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', '33') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '33', node2.lemma)
    
    # L -> du_tp -> -166 -> L
    feats[926] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', '-166') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-166', node2.lemma)
    
    # L -> du_tp -> 101 -> L
    feats[927] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', '101') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    
    # L -> du_tp -> 14_sp -> L
    feats[928] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', '14_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_sp', node2.lemma)
    
    # L -> du_tp -> -42 -> L
    feats[929] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', '-42') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-42', node2.lemma)
    
    # L -> du_tp -> sg_sp -> L
    feats[930] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', 'sg_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_sp', node2.lemma)
    
    # L -> du_tp -> 40 -> L
    feats[931] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', '40') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '40', node2.lemma)
    
    # L -> du_tp -> 80 -> L
    feats[932] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', '80') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '80', node2.lemma)
    
    # L -> du_tp -> -246 -> L
    feats[933] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', '-246') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-246', node2.lemma)
    
    # L -> du_tp -> 115 -> L
    feats[934] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', '115') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '115', node2.lemma)
    
    # L -> du_tp -> 131 -> L
    feats[935] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', '131') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '131', node2.lemma)
    
    # L -> du_tp -> 2_tp -> L
    feats[936] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', '2_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_tp', node2.lemma)
    
    # L -> du_tp -> 11_du -> L
    feats[937] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', '11_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_du', node2.lemma)
    
    # L -> du_tp -> 137 -> L
    feats[938] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', '137') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    
    # L -> du_tp -> -129 -> L
    feats[939] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', '-129') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-129', node2.lemma)
    
    # L -> du_tp -> -101 -> L
    feats[940] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', '-101') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-101', node2.lemma)
    
    # L -> du_tp -> 70 -> L
    feats[941] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', '70') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '70', node2.lemma)
    
    # L -> du_tp -> voc. neutr. -> L
    feats[942] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', 'voc. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. neutr.', node2.lemma)
    
    # L -> du_tp -> -115 -> L
    feats[943] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', '-115') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-115', node2.lemma)
    
    # L -> du_tp -> 35 -> L
    feats[944] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', '35') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '35', node2.lemma)
    
    # L -> du_tp -> acc. adj. -> L
    feats[945] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', 'acc. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. adj.', node2.lemma)
    
    # L -> du_tp -> acc -> L
    feats[946] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', 'acc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc', node2.lemma)
    
    # L -> du_tp -> -22 -> L
    feats[947] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', '-22') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-22', node2.lemma)
    
    # L -> du_tp -> 6_du -> L
    feats[948] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', '6_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_du', node2.lemma)
    
    # L -> du_tp -> abl. du. -> L
    feats[949] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', 'abl. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. du.', node2.lemma)
    
    # L -> du_tp -> 177 -> L
    feats[950] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', '177') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '177', node2.lemma)
    
    # L -> du_tp -> dat. sg. -> L
    feats[951] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', 'dat. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. sg.', node2.lemma)
    
    # L -> 135 -> 102 -> L
    feats[952] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '135') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', '102') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '102', node2.lemma)
    
    # L -> 135 -> 15_fp -> L
    feats[953] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '135') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', '15_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_fp', node2.lemma)
    
    # L -> 135 -> -46 -> L
    feats[954] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '135') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', '-46') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-46', node2.lemma)
    
    # L -> 135 -> abl. du. -> L
    feats[955] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '135') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', 'abl. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. du.', node2.lemma)
    
    # L -> 135 -> fem -> L
    feats[956] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '135') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', 'fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'fem', node2.lemma)
    
    # L -> 135 -> 79 -> L
    feats[957] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '135') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', '79') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '79', node2.lemma)
    
    # L -> 135 -> -269 -> L
    feats[958] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '135') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', '-269') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-269', node2.lemma)
    
    # L -> 135 -> -82 -> L
    feats[959] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '135') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', '-82') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-82', node2.lemma)
    
    # L -> 135 -> 61 -> L
    feats[960] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '135') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', '61') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '61', node2.lemma)
    
    # L -> 135 -> -157 -> L
    feats[961] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '135') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', '-157') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-157', node2.lemma)
    
    # L -> 135 -> -119 -> L
    feats[962] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '135') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', '-119') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-119', node2.lemma)
    
    # L -> 135 -> 32 -> L
    feats[963] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '135') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', '32') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '32', node2.lemma)
    
    # L -> 135 -> -133 -> L
    feats[964] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '135') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', '-133') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-133', node2.lemma)
    
    # L -> 135 -> 3 -> L
    feats[965] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '135') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', '3') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3', node2.lemma)
    
    # L -> 135 -> 56 -> L
    feats[966] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '135') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', '56') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '56', node2.lemma)
    
    # L -> 135 -> -64 -> L
    feats[967] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '135') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', '-64') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-64', node2.lemma)
    
    # L -> 135 -> -103 -> L
    feats[968] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '135') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', '-103') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-103', node2.lemma)
    
    # L -> 135 -> -123 -> L
    feats[969] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '135') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', '-123') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-123', node2.lemma)
    
    # L -> 135 -> 100 -> L
    feats[970] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '135') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', '100') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '100', node2.lemma)
    
    # L -> 135 -> nom. neutr. -> L
    feats[971] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '135') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', 'nom. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. neutr.', node2.lemma)
    
    # L -> 135 -> 49 -> L
    feats[972] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '135') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', '49') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '49', node2.lemma)
    
    # L -> 135 -> 160 -> L
    feats[973] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '135') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', '160') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    
    # L -> 135 -> -33 -> L
    feats[974] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '135') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', '-33') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-33', node2.lemma)
    
    # L -> 135 -> -260 -> L
    feats[975] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '135') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', '-260') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-260', node2.lemma)
    
    # L -> 135 -> -30 -> L
    feats[976] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '135') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', '-30') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-30', node2.lemma)
    
    # L -> 135 -> -121 -> L
    feats[977] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '135') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', '-121') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-121', node2.lemma)
    
    # L -> 135 -> 168 -> L
    feats[978] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '135') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', '168') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '168', node2.lemma)
    
    # L -> 135 -> 51 -> L
    feats[979] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '135') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', '51') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '51', node2.lemma)
    
    # L -> 135 -> 16_sg -> L
    feats[980] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '135') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', '16_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_sg', node2.lemma)
    
    # L -> 135 -> -144 -> L
    feats[981] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '135') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', '-144') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-144', node2.lemma)
    
    # L -> -48 -> voc. masc. -> L
    feats[982] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', 'voc. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. masc.', node2.lemma)
    
    # L -> -48 -> 2_sp -> L
    feats[983] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', '2_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_sp', node2.lemma)
    
    # L -> -48 -> 10_sp -> L
    feats[984] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', '10_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_sp', node2.lemma)
    
    # L -> -48 -> -152 -> L
    feats[985] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', '-152') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-152', node2.lemma)
    
    # L -> -48 -> -66 -> L
    feats[986] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', '-66') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-66', node2.lemma)
    
    # L -> -48 -> 5_sg -> L
    feats[987] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', '5_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sg', node2.lemma)
    
    # L -> -48 -> loc -> L
    feats[988] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', 'loc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc', node2.lemma)
    
    # L -> -48 -> -243 -> L
    feats[989] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', '-243') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-243', node2.lemma)
    
    # L -> -48 -> -109 -> L
    feats[990] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', '-109') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-109', node2.lemma)
    
    # L -> -48 -> 6_sg -> L
    feats[991] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', '6_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_sg', node2.lemma)
    
    # L -> -48 -> 55 -> L
    feats[992] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', '55') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '55', node2.lemma)
    
    # L -> -48 -> tp -> L
    feats[993] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', 'tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'tp', node2.lemma)
    
    # L -> -48 -> -98 -> L
    feats[994] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', '-98') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-98', node2.lemma)
    
    # L -> -48 -> gen. du. -> L
    feats[995] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', 'gen. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen. du.', node2.lemma)
    
    # L -> -48 -> 5_fp -> L
    feats[996] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', '5_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_fp', node2.lemma)
    
    # L -> -48 -> 8_fp -> L
    feats[997] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', '8_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    
    # L -> -48 -> 27_du -> L
    feats[998] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', '27_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_du', node2.lemma)
    
    # L -> -48 -> 171 -> L
    feats[999] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', '171') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '171', node2.lemma)
    
    # L -> -48 -> 27_tp -> L
    feats[1000] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', '27_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_tp', node2.lemma)
    
    # L -> -48 -> 27_sg -> L
    feats[1001] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', '27_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_sg', node2.lemma)
    
    # L -> -48 -> abl. sg. -> L
    feats[1002] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', 'abl. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. sg.', node2.lemma)
    
    # L -> -48 -> -276 -> L
    feats[1003] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', '-276') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-276', node2.lemma)
    
    # L -> loc. pl. -> loc. du. -> L
    feats[1004] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. pl.', 'loc. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc. du.', node2.lemma)
    
    # L -> loc. pl. -> nom. du. -> L
    feats[1005] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. pl.', 'nom. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. du.', node2.lemma)
    
    # L -> loc. pl. -> -166 -> L
    feats[1006] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. pl.', '-166') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-166', node2.lemma)
    
    # L -> loc. pl. -> 6_fp -> L
    feats[1007] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. pl.', '6_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_fp', node2.lemma)
    
    # L -> loc. pl. -> 89 -> L
    feats[1008] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. pl.', '89') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '89', node2.lemma)
    
    # L -> loc. pl. -> 80 -> L
    feats[1009] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. pl.', '80') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '80', node2.lemma)
    
    # L -> loc. pl. -> -49 -> L
    feats[1010] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. pl.', '-49') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-49', node2.lemma)
    
    # L -> loc. pl. -> -63 -> L
    feats[1011] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. pl.', '-63') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-63', node2.lemma)
    
    # L -> loc. pl. -> -230 -> L
    feats[1012] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. pl.', '-230') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-230', node2.lemma)
    
    # L -> loc. pl. -> 157 -> L
    feats[1013] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. pl.', '157') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '157', node2.lemma)
    
    # L -> loc. pl. -> 5_tp -> L
    feats[1014] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. pl.', '5_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_tp', node2.lemma)
    
    # L -> loc. pl. -> 122 -> L
    feats[1015] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. pl.', '122') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '122', node2.lemma)
    
    # L -> loc. pl. -> -129 -> L
    feats[1016] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. pl.', '-129') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-129', node2.lemma)
    
    # L -> loc. pl. -> -241 -> L
    feats[1017] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. pl.', '-241') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-241', node2.lemma)
    
    # L -> loc. pl. -> -11 -> L
    feats[1018] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. pl.', '-11') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-11', node2.lemma)
    
    # L -> loc. pl. -> -271 -> L
    feats[1019] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. pl.', '-271') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-271', node2.lemma)
    
    # L -> loc. pl. -> -269 -> L
    feats[1020] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. pl.', '-269') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-269', node2.lemma)
    
    # L -> loc. pl. -> -309 -> L
    feats[1021] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. pl.', '-309') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-309', node2.lemma)
    
    # L -> loc. pl. -> -87 -> L
    feats[1022] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. pl.', '-87') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-87', node2.lemma)
    
    # L -> loc. pl. -> -133 -> L
    feats[1023] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. pl.', '-133') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-133', node2.lemma)
    
    # L -> loc. pl. -> 3 -> L
    feats[1024] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. pl.', '3') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3', node2.lemma)
    
    # L -> loc. pl. -> -103 -> L
    feats[1025] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. pl.', '-103') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-103', node2.lemma)
    
    # L -> loc. pl. -> nom. neutr. -> L
    feats[1026] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. pl.', 'nom. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. neutr.', node2.lemma)
    
    # L -> loc. pl. -> 49 -> L
    feats[1027] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. pl.', '49') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '49', node2.lemma)
    
    # L -> loc. pl. -> -39 -> L
    feats[1028] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. pl.', '-39') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-39', node2.lemma)
    
    # L -> loc. pl. -> 160 -> L
    feats[1029] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. pl.', '160') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    
    # L -> loc. pl. -> nom. adj. -> L
    feats[1030] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. pl.', 'nom. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. adj.', node2.lemma)
    
    # L -> loc. pl. -> -73 -> L
    feats[1031] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. pl.', '-73') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-73', node2.lemma)
    
    # L -> loc. pl. -> 92 -> L
    feats[1032] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. pl.', '92') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '92', node2.lemma)
    
    # L -> loc. pl. -> -86 -> L
    feats[1033] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. pl.', '-86') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-86', node2.lemma)
    
    # L -> loc. pl. -> 11_pl -> L
    feats[1034] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. pl.', '11_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_pl', node2.lemma)
    
    # L -> loc. pl. -> -161 -> L
    feats[1035] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. pl.', '-161') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-161', node2.lemma)
    
    # L -> 117 -> 111 -> L
    feats[1036] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '117') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '117', '111') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '111', node2.lemma)
    
    # L -> 155 -> voc. masc. -> L
    feats[1037] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '155') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '155', 'voc. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. masc.', node2.lemma)
    
    # L -> 155 -> nom. sg. -> L
    feats[1038] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '155') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '155', 'nom. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. sg.', node2.lemma)
    
    # L -> 155 -> -94 -> L
    feats[1039] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '155') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '155', '-94') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-94', node2.lemma)
    
    # L -> -297 -> -262 -> L
    feats[1040] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-297') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-297', '-262') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-262', node2.lemma)
    
    # L -> -297 -> 171 -> L
    feats[1041] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-297') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-297', '171') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '171', node2.lemma)
    
    # L -> -297 -> 9_tp -> L
    feats[1042] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-297') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-297', '9_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_tp', node2.lemma)
    
    # L -> -297 -> -84 -> L
    feats[1043] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-297') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-297', '-84') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-84', node2.lemma)
    
    # L -> -297 -> 175 -> L
    feats[1044] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-297') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-297', '175') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '175', node2.lemma)
    
    # L -> -297 -> acc. du. -> L
    feats[1045] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-297') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-297', 'acc. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. du.', node2.lemma)
    
    # L -> 9_tp -> 4_fp -> L
    feats[1046] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_tp', '4_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_fp', node2.lemma)
    
    # L -> 9_tp -> 30_fp -> L
    feats[1047] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_tp', '30_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_fp', node2.lemma)
    
    # L -> 9_tp -> -61 -> L
    feats[1048] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_tp', '-61') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-61', node2.lemma)
    
    # L -> 9_tp -> 71 -> L
    feats[1049] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_tp', '71') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '71', node2.lemma)
    
    # L -> 9_tp -> instr. sg. -> L
    feats[1050] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_tp', 'instr. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. sg.', node2.lemma)
    
    # L -> -91 -> -133 -> L
    feats[1051] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-91') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-91', '-133') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-133', node2.lemma)
    
    # L -> -91 -> -169 -> L
    feats[1052] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-91') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-91', '-169') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-169', node2.lemma)
    
    # L -> -91 -> 129 -> L
    feats[1053] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-91') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-91', '129') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '129', node2.lemma)
    
    # L -> -91 -> -64 -> L
    feats[1054] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-91') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-91', '-64') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-64', node2.lemma)
    
    # L -> -91 -> -103 -> L
    feats[1055] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-91') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-91', '-103') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-103', node2.lemma)
    
    # L -> -91 -> 100 -> L
    feats[1056] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-91') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-91', '100') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '100', node2.lemma)
    
    # L -> -91 -> 49 -> L
    feats[1057] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-91') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-91', '49') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '49', node2.lemma)
    
    # L -> -91 -> -39 -> L
    feats[1058] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-91') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-91', '-39') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-39', node2.lemma)
    
    # L -> -91 -> 141 -> L
    feats[1059] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-91') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-91', '141') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '141', node2.lemma)
    
    # L -> -91 -> nom. adj. -> L
    feats[1060] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-91') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-91', 'nom. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. adj.', node2.lemma)
    
    # L -> -91 -> -20 -> L
    feats[1061] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-91') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-91', '-20') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-20', node2.lemma)
    
    # L -> -91 -> dat -> L
    feats[1062] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-91') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-91', 'dat') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat', node2.lemma)
    
    # L -> -91 -> -73 -> L
    feats[1063] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-91') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-91', '-73') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-73', node2.lemma)
    
    # L -> -91 -> nom. sg. -> L
    feats[1064] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-91') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-91', 'nom. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. sg.', node2.lemma)
    
    # L -> -91 -> 75 -> L
    feats[1065] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-91') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-91', '75') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '75', node2.lemma)
    
    # L -> -91 -> 30_sg -> L
    feats[1066] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-91') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-91', '30_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_sg', node2.lemma)
    
    # L -> -91 -> instr. pl. -> L
    feats[1067] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-91') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-91', 'instr. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. pl.', node2.lemma)
    
    # L -> -91 -> acc. pl. -> L
    feats[1068] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-91') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-91', 'acc. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. pl.', node2.lemma)
    
    # L -> -91 -> -161 -> L
    feats[1069] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-91') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-91', '-161') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-161', node2.lemma)
    
    # L -> -91 -> 180 -> L
    feats[1070] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-91') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-91', '180') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '180', node2.lemma)
    
    # L -> -91 -> 28 -> L
    feats[1071] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-91') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-91', '28') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
    
    # L -> -84 -> -50 -> L
    feats[1072] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-84') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-84', '-50') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-50', node2.lemma)
    
    # L -> -84 -> 5_sg -> L
    feats[1073] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-84') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-84', '5_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sg', node2.lemma)
    
    # L -> -84 -> loc -> L
    feats[1074] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-84') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-84', 'loc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc', node2.lemma)
    
    # L -> -84 -> -141 -> L
    feats[1075] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-84') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-84', '-141') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-141', node2.lemma)
    
    # L -> -84 -> -210 -> L
    feats[1076] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-84') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-84', '-210') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-210', node2.lemma)
    
    # L -> -84 -> -98 -> L
    feats[1077] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-84') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-84', '-98') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-98', node2.lemma)
    
    # L -> -84 -> 5_fp -> L
    feats[1078] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-84') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-84', '5_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_fp', node2.lemma)
    
    # L -> -84 -> -16 -> L
    feats[1079] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-84') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-84', '-16') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-16', node2.lemma)
    
    # L -> -84 -> -78 -> L
    feats[1080] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-84') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-84', '-78') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-78', node2.lemma)
    
    # L -> -84 -> 12_fp -> L
    feats[1081] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-84') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-84', '12_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_fp', node2.lemma)
    
    # L -> -84 -> 59 -> L
    feats[1082] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-84') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-84', '59') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '59', node2.lemma)
    
    # L -> -84 -> -302 -> L
    feats[1083] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-84') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-84', '-302') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-302', node2.lemma)
    
    # L -> -84 -> -29 -> L
    feats[1084] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-84') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-84', '-29') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-29', node2.lemma)
    
    # L -> -84 -> -273 -> L
    feats[1085] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-84') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-84', '-273') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-273', node2.lemma)
    
    # L -> -84 -> -13 -> L
    feats[1086] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-84') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-84', '-13') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-13', node2.lemma)
    
    # L -> -84 -> 5_du -> L
    feats[1087] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-84') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-84', '5_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_du', node2.lemma)
    
    # L -> -84 -> -84 -> L
    feats[1088] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-84') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-84', '-84') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-84', node2.lemma)
    
    # L -> -84 -> -89 -> L
    feats[1089] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-84') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-84', '-89') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-89', node2.lemma)
    
    # L -> -84 -> 136 -> L
    feats[1090] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-84') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-84', '136') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '136', node2.lemma)
    
    # L -> -84 -> 54 -> L
    feats[1091] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-84') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-84', '54') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    
    # L -> -242 -> -50 -> L
    feats[1092] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-242') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-242', '-50') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-50', node2.lemma)
    
    # L -> -242 -> -66 -> L
    feats[1093] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-242') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-242', '-66') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-66', node2.lemma)
    
    # L -> -242 -> nom. masc. -> L
    feats[1094] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-242') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-242', 'nom. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. masc.', node2.lemma)
    
    # L -> -242 -> -141 -> L
    feats[1095] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-242') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-242', '-141') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-141', node2.lemma)
    
    # L -> -242 -> -69 -> L
    feats[1096] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-242') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-242', '-69') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-69', node2.lemma)
    
    # L -> -242 -> 30_du -> L
    feats[1097] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-242') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-242', '30_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_du', node2.lemma)
    
    # L -> -242 -> 162 -> L
    feats[1098] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-242') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-242', '162') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '162', node2.lemma)
    
    # L -> -242 -> -243 -> L
    feats[1099] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-242') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-242', '-243') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-243', node2.lemma)
    
    # L -> -242 -> 6_sg -> L
    feats[1100] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-242') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-242', '6_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_sg', node2.lemma)
    
    # L -> -242 -> 55 -> L
    feats[1101] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-242') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-242', '55') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '55', node2.lemma)
    
    # L -> -242 -> 16_fp -> L
    feats[1102] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-242') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-242', '16_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_fp', node2.lemma)
    
    # L -> -242 -> tp -> L
    feats[1103] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-242') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-242', 'tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'tp', node2.lemma)
    
    # L -> -242 -> pl_sp -> L
    feats[1104] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-242') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-242', 'pl_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_sp', node2.lemma)
    
    # L -> -242 -> -83 -> L
    feats[1105] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-242') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-242', '-83') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-83', node2.lemma)
    
    # L -> -242 -> 27_du -> L
    feats[1106] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-242') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-242', '27_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_du', node2.lemma)
    
    # L -> -242 -> 171 -> L
    feats[1107] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-242') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-242', '171') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '171', node2.lemma)
    
    # L -> -242 -> 12_fp -> L
    feats[1108] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-242') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-242', '12_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_fp', node2.lemma)
    
    # L -> -242 -> 3_du -> L
    feats[1109] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-242') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-242', '3_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    
    # L -> -242 -> -139 -> L
    feats[1110] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-242') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-242', '-139') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-139', node2.lemma)
    
    # L -> -242 -> -122 -> L
    feats[1111] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-242') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-242', '-122') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-122', node2.lemma)
    
    # L -> -242 -> 155 -> L
    feats[1112] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-242') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-242', '155') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '155', node2.lemma)
    
    # L -> -242 -> 37 -> L
    feats[1113] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-242') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-242', '37') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '37', node2.lemma)
    
    # L -> -242 -> 154 -> L
    feats[1114] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-242') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-242', '154') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '154', node2.lemma)
    
    # L -> -242 -> 108 -> L
    feats[1115] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-242') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-242', '108') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '108', node2.lemma)
    
    # L -> -242 -> pl -> L
    feats[1116] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-242') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-242', 'pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    
    # L -> -242 -> 4_tp -> L
    feats[1117] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-242') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-242', '4_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_tp', node2.lemma)
    
    # L -> -242 -> -41 -> L
    feats[1118] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-242') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-242', '-41') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-41', node2.lemma)
    
    # L -> -242 -> -97 -> L
    feats[1119] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-242') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-242', '-97') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-97', node2.lemma)
    
    # L -> -52 -> -77 -> L
    feats[1120] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-52') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-52', '-77') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-77', node2.lemma)
    
    # L -> -52 -> -271 -> L
    feats[1121] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-52') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-52', '-271') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-271', node2.lemma)
    
    # L -> -52 -> -131 -> L
    feats[1122] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-52') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-52', '-131') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-131', node2.lemma)
    
    # L -> -52 -> -143 -> L
    feats[1123] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-52') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-52', '-143') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-143', node2.lemma)
    
    # L -> -52 -> 9_pl -> L
    feats[1124] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-52') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-52', '9_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_pl', node2.lemma)
    
    # L -> -52 -> -82 -> L
    feats[1125] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-52') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-52', '-82') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-82', node2.lemma)
    
    # L -> -52 -> 9_sg -> L
    feats[1126] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-52') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-52', '9_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_sg', node2.lemma)
    
    # L -> -52 -> 140 -> L
    feats[1127] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-52') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-52', '140') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '140', node2.lemma)
    
    # L -> -52 -> -119 -> L
    feats[1128] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-52') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-52', '-119') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-119', node2.lemma)
    
    # L -> -52 -> 32 -> L
    feats[1129] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-52') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-52', '32') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '32', node2.lemma)
    
    # L -> -52 -> 15_tp -> L
    feats[1130] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-52') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-52', '15_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_tp', node2.lemma)
    
    # L -> -52 -> -64 -> L
    feats[1131] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-52') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-52', '-64') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-64', node2.lemma)
    
    # L -> -52 -> -10 -> L
    feats[1132] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-52') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-52', '-10') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-10', node2.lemma)
    
    # L -> -52 -> -123 -> L
    feats[1133] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-52') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-52', '-123') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-123', node2.lemma)
    
    # L -> -52 -> -37 -> L
    feats[1134] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-52') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-52', '-37') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-37', node2.lemma)
    
    # L -> -52 -> 100 -> L
    feats[1135] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-52') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-52', '100') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '100', node2.lemma)
    
    # L -> -52 -> 8_sg -> L
    feats[1136] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-52') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-52', '8_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_sg', node2.lemma)
    
    # L -> -52 -> -39 -> L
    feats[1137] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-52') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-52', '-39') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-39', node2.lemma)
    
    # L -> -52 -> 98 -> L
    feats[1138] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-52') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-52', '98') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '98', node2.lemma)
    
    # L -> -52 -> -149 -> L
    feats[1139] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-52') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-52', '-149') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-149', node2.lemma)
    
    # L -> -52 -> dat -> L
    feats[1140] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-52') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-52', 'dat') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat', node2.lemma)
    
    # L -> -52 -> -30 -> L
    feats[1141] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-52') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-52', '-30') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-30', node2.lemma)
    
    # L -> -52 -> -151 -> L
    feats[1142] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-52') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-52', '-151') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-151', node2.lemma)
    
    # L -> -52 -> instr. pl. -> L
    feats[1143] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-52') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-52', 'instr. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. pl.', node2.lemma)
    
    # L -> -52 -> -121 -> L
    feats[1144] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-52') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-52', '-121') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-121', node2.lemma)
    
    # L -> 153 -> 13_pl -> L
    feats[1145] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '153') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '153', '13_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_pl', node2.lemma)
    
    # L -> 153 -> -152 -> L
    feats[1146] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '153') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '153', '-152') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-152', node2.lemma)
    
    # L -> 153 -> 132 -> L
    feats[1147] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '153') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '153', '132') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '132', node2.lemma)
    
    # L -> 153 -> -50 -> L
    feats[1148] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '153') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '153', '-50') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-50', node2.lemma)
    
    # L -> 153 -> -210 -> L
    feats[1149] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '153') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '153', '-210') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-210', node2.lemma)
    
    # L -> 153 -> 6_sg -> L
    feats[1150] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '153') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '153', '6_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_sg', node2.lemma)
    
    # L -> 153 -> -28 -> L
    feats[1151] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '153') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '153', '-28') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-28', node2.lemma)
    
    # L -> 153 -> -98 -> L
    feats[1152] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '153') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '153', '-98') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-98', node2.lemma)
    
    # L -> 153 -> 16_tp -> L
    feats[1153] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '153') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '153', '16_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_tp', node2.lemma)
    
    # L -> 153 -> 7_pl -> L
    feats[1154] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '153') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '153', '7_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_pl', node2.lemma)
    
    # L -> 153 -> -78 -> L
    feats[1155] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '153') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '153', '-78') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-78', node2.lemma)
    
    # L -> 153 -> 12_fp -> L
    feats[1156] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '153') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '153', '12_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_fp', node2.lemma)
    
    # L -> 153 -> -43 -> L
    feats[1157] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '153') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '153', '-43') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-43', node2.lemma)
    
    # L -> 153 -> -302 -> L
    feats[1158] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '153') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '153', '-302') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-302', node2.lemma)
    
    # L -> 153 -> du_tp -> L
    feats[1159] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '153') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '153', 'du_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_tp', node2.lemma)
    
    # L -> 76 -> 128 -> L
    feats[1160] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '76') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '76', '128') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '128', node2.lemma)
    
    # L -> 76 -> 159 -> L
    feats[1161] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '76') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '76', '159') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '159', node2.lemma)
    
    # L -> 76 -> 182 -> L
    feats[1162] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '76') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '76', '182') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '182', node2.lemma)
    
    # L -> 76 -> 6_tp -> L
    feats[1163] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '76') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '76', '6_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_tp', node2.lemma)
    
    # L -> 76 -> 69 -> L
    feats[1164] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '76') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '76', '69') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '69', node2.lemma)
    
    # L -> 76 -> -142 -> L
    feats[1165] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '76') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '76', '-142') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-142', node2.lemma)
    
    # L -> 76 -> -17 -> L
    feats[1166] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '76') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '76', '-17') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-17', node2.lemma)
    
    # L -> 76 -> 15_sp -> L
    feats[1167] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '76') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '76', '15_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_sp', node2.lemma)
    
    # L -> 76 -> -79 -> L
    feats[1168] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '76') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '76', '-79') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-79', node2.lemma)
    
    # L -> 76 -> -90 -> L
    feats[1169] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '76') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '76', '-90') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-90', node2.lemma)
    
    # L -> 76 -> -292 -> L
    feats[1170] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '76') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '76', '-292') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-292', node2.lemma)
    
    # L -> 76 -> 151 -> L
    feats[1171] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '76') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '76', '151') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '151', node2.lemma)
    
    # L -> 76 -> nom. du. -> L
    feats[1172] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '76') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '76', 'nom. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. du.', node2.lemma)
    
    # L -> 76 -> -291 -> L
    feats[1173] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '76') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '76', '-291') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-291', node2.lemma)
    
    # L -> 76 -> -42 -> L
    feats[1174] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '76') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '76', '-42') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-42', node2.lemma)
    
    # L -> 2_pl -> 51 -> L
    feats[1175] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_pl', '51') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '51', node2.lemma)
    
    # L -> 2_pl -> 16_sg -> L
    feats[1176] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_pl', '16_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_sg', node2.lemma)
    
    # L -> 2_pl -> -144 -> L
    feats[1177] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_pl', '-144') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-144', node2.lemma)
    
    # L -> 2_pl -> 28 -> L
    feats[1178] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_pl', '28') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
    
    # L -> sp -> voc. masc. -> L
    feats[1179] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', 'voc. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. masc.', node2.lemma)
    
    # L -> sp -> 4_du -> L
    feats[1180] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', '4_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_du', node2.lemma)
    
    # L -> sp -> -268 -> L
    feats[1181] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', '-268') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-268', node2.lemma)
    
    # L -> sp -> 2_sp -> L
    feats[1182] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', '2_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_sp', node2.lemma)
    
    # L -> sp -> -293 -> L
    feats[1183] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', '-293') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-293', node2.lemma)
    
    # L -> sp -> 13_pl -> L
    feats[1184] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', '13_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_pl', node2.lemma)
    
    # L -> sp -> 132 -> L
    feats[1185] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', '132') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '132', node2.lemma)
    
    # L -> sp -> -50 -> L
    feats[1186] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', '-50') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-50', node2.lemma)
    
    # L -> sp -> 3_sp -> L
    feats[1187] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', '3_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_sp', node2.lemma)
    
    # L -> sp -> -262 -> L
    feats[1188] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', '-262') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-262', node2.lemma)
    
    # L -> sp -> 162 -> L
    feats[1189] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', '162') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '162', node2.lemma)
    
    # L -> sp -> 16_fp -> L
    feats[1190] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', '16_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_fp', node2.lemma)
    
    # L -> sp -> -18 -> L
    feats[1191] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', '-18') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-18', node2.lemma)
    
    # L -> sp -> pl_sp -> L
    feats[1192] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', 'pl_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_sp', node2.lemma)
    
    # L -> sp -> 7_pl -> L
    feats[1193] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', '7_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_pl', node2.lemma)
    
    # L -> sp -> 171 -> L
    feats[1194] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', '171') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '171', node2.lemma)
    
    # L -> sp -> -67 -> L
    feats[1195] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', '-67') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-67', node2.lemma)
    
    # L -> sp -> 12_fp -> L
    feats[1196] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', '12_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_fp', node2.lemma)
    
    # L -> sp -> -43 -> L
    feats[1197] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', '-43') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-43', node2.lemma)
    
    # L -> sp -> -302 -> L
    feats[1198] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', '-302') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-302', node2.lemma)
    
    # L -> sp -> -126 -> L
    feats[1199] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', '-126') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-126', node2.lemma)
    
    # L -> sp -> 9_tp -> L
    feats[1200] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', '9_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_tp', node2.lemma)
    
    # L -> sp -> 37 -> L
    feats[1201] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', '37') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '37', node2.lemma)
    
    # L -> sp -> 138 -> L
    feats[1202] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', '138') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '138', node2.lemma)
    
    # L -> sp -> 136 -> L
    feats[1203] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', '136') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '136', node2.lemma)
    
    # L -> sp -> 76 -> L
    feats[1204] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', '76') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '76', node2.lemma)
    
    # L -> sp -> sp -> L
    feats[1205] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', 'sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sp', node2.lemma)
    
    # L -> sp -> 29_sg -> L
    feats[1206] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', '29_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '29_sg', node2.lemma)
    
    # L -> sp -> -132 -> L
    feats[1207] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', '-132') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-132', node2.lemma)
    
    # L -> sp -> -307 -> L
    feats[1208] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', '-307') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-307', node2.lemma)
    
    # L -> sp -> pl_tp -> L
    feats[1209] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', 'pl_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_tp', node2.lemma)
    
    # L -> sp -> instr. neutr. -> L
    feats[1210] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', 'instr. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. neutr.', node2.lemma)
    
    # L -> sp -> 2_sg -> L
    feats[1211] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', '2_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_sg', node2.lemma)
    
    # L -> sp -> 72 -> L
    feats[1212] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', '72') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '72', node2.lemma)
    
    # L -> sp -> 130 -> L
    feats[1213] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', '130') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '130', node2.lemma)
    
    # L -> sp -> 119 -> L
    feats[1214] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', '119') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '119', node2.lemma)
    
    # L -> sp -> 27_fp -> L
    feats[1215] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', '27_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_fp', node2.lemma)
    
    # L -> sp -> 128 -> L
    feats[1216] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', '128') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '128', node2.lemma)
    
    # L -> sp -> -301 -> L
    feats[1217] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', '-301') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-301', node2.lemma)
    
    # L -> sp -> 159 -> L
    feats[1218] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', '159') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '159', node2.lemma)
    
    # L -> sp -> 91 -> L
    feats[1219] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', '91') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    
    # L -> sp -> -92 -> L
    feats[1220] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', '-92') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-92', node2.lemma)
    
    # L -> sp -> -303 -> L
    feats[1221] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', '-303') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-303', node2.lemma)
    
    # L -> sp -> 11_tp -> L
    feats[1222] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', '11_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_tp', node2.lemma)
    
    # L -> sp -> -154 -> L
    feats[1223] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', '-154') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-154', node2.lemma)
    
    # L -> sp -> dat. pl. -> L
    feats[1224] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', 'dat. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. pl.', node2.lemma)
    
    # L -> sp -> -41 -> L
    feats[1225] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', '-41') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-41', node2.lemma)
    
    # L -> 14_tp -> -143 -> L
    feats[1226] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_tp', '-143') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-143', node2.lemma)
    
    # L -> 14_tp -> 149 -> L
    feats[1227] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_tp', '149') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '149', node2.lemma)
    
    # L -> 14_tp -> -299 -> L
    feats[1228] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_tp', '-299') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-299', node2.lemma)
    
    # L -> 14_tp -> -309 -> L
    feats[1229] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_tp', '-309') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-309', node2.lemma)
    
    # L -> 14_tp -> 2_fp -> L
    feats[1230] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_tp', '2_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_fp', node2.lemma)
    
    # L -> 14_tp -> -39 -> L
    feats[1231] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_tp', '-39') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-39', node2.lemma)
    
    # L -> 14_tp -> 98 -> L
    feats[1232] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_tp', '98') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '98', node2.lemma)
    
    # L -> 14_tp -> dat -> L
    feats[1233] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_tp', 'dat') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat', node2.lemma)
    
    # L -> 14_tp -> -150 -> L
    feats[1234] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_tp', '-150') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-150', node2.lemma)
    
    # L -> 114 -> 2_sp -> L
    feats[1235] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '114') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '114', '2_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_sp', node2.lemma)
    
    # L -> 114 -> -141 -> L
    feats[1236] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '114') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '114', '-141') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-141', node2.lemma)
    
    # L -> 114 -> 162 -> L
    feats[1237] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '114') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '114', '162') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '162', node2.lemma)
    
    # L -> 114 -> 171 -> L
    feats[1238] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '114') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '114', '171') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '171', node2.lemma)
    
    # L -> 114 -> 27_sg -> L
    feats[1239] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '114') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '114', '27_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_sg', node2.lemma)
    
    # L -> 175 -> 27_pl -> L
    feats[1240] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '175') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '175', '27_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_pl', node2.lemma)
    
    # L -> 175 -> 121 -> L
    feats[1241] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '175') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '175', '121') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '121', node2.lemma)
    
    # L -> 175 -> 8_sp -> L
    feats[1242] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '175') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '175', '8_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_sp', node2.lemma)
    
    # L -> 175 -> 50 -> L
    feats[1243] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '175') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '175', '50') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '50', node2.lemma)
    
    # L -> 175 -> -15 -> L
    feats[1244] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '175') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '175', '-15') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-15', node2.lemma)
    
    # L -> 175 -> 120 -> L
    feats[1245] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '175') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '175', '120') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '120', node2.lemma)
    
    # L -> 175 -> 8_pl -> L
    feats[1246] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '175') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '175', '8_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_pl', node2.lemma)
    
    # L -> 175 -> -249 -> L
    feats[1247] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '175') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '175', '-249') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-249', node2.lemma)
    
    # L -> 175 -> -292 -> L
    feats[1248] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '175') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '175', '-292') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-292', node2.lemma)
    
    # L -> 175 -> -57 -> L
    feats[1249] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '175') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '175', '-57') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-57', node2.lemma)
    
    # L -> 175 -> -306 -> L
    feats[1250] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '175') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '175', '-306') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-306', node2.lemma)
    
    # L -> 175 -> -38 -> L
    feats[1251] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '175') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '175', '-38') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-38', node2.lemma)
    
    # L -> 175 -> 73 -> L
    feats[1252] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '175') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '175', '73') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '73', node2.lemma)
    
    # L -> 175 -> du -> L
    feats[1253] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '175') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '175', 'du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du', node2.lemma)
    
    # L -> 175 -> loc. du. -> L
    feats[1254] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '175') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '175', 'loc. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc. du.', node2.lemma)
    
    # L -> 175 -> 97 -> L
    feats[1255] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '175') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '175', '97') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '97', node2.lemma)
    
    # L -> 175 -> -104 -> L
    feats[1256] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '175') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '175', '-104') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-104', node2.lemma)
    
    # L -> 175 -> 151 -> L
    feats[1257] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '175') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '175', '151') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '151', node2.lemma)
    
    # L -> 175 -> 71 -> L
    feats[1258] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '175') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '175', '71') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '71', node2.lemma)
    
    # L -> 175 -> 10_tp -> L
    feats[1259] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '175') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '175', '10_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_tp', node2.lemma)
    
    # L -> 175 -> nom. du. -> L
    feats[1260] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '175') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '175', 'nom. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. du.', node2.lemma)
    
    # L -> 175 -> instr. sg. -> L
    feats[1261] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '175') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '175', 'instr. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. sg.', node2.lemma)
    
    # L -> 175 -> -42 -> L
    feats[1262] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '175') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '175', '-42') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-42', node2.lemma)
    
    # L -> 175 -> -266 -> L
    feats[1263] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '175') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '175', '-266') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-266', node2.lemma)
    
    # L -> 175 -> 5_pl -> L
    feats[1264] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '175') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '175', '5_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_pl', node2.lemma)
    
    # L -> 175 -> 6_fp -> L
    feats[1265] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '175') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '175', '6_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_fp', node2.lemma)
    
    # L -> 175 -> 12_sp -> L
    feats[1266] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '175') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '175', '12_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_sp', node2.lemma)
    
    # L -> 175 -> 40 -> L
    feats[1267] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '175') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '175', '40') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '40', node2.lemma)
    
    # L -> 175 -> 80 -> L
    feats[1268] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '175') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '175', '80') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '80', node2.lemma)
    
    # L -> 175 -> 7_sp -> L
    feats[1269] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '175') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '175', '7_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_sp', node2.lemma)
    
    # L -> 175 -> 11_du -> L
    feats[1270] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '175') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '175', '11_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_du', node2.lemma)
    
    # L -> 175 -> 15_du -> L
    feats[1271] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '175') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '175', '15_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_du', node2.lemma)
    
    # L -> 175 -> 122 -> L
    feats[1272] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '175') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '175', '122') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '122', node2.lemma)
    
    # L -> 175 -> 169 -> L
    feats[1273] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '175') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '175', '169') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '169', node2.lemma)
    
    # L -> 175 -> acc -> L
    feats[1274] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '175') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '175', 'acc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc', node2.lemma)
    
    # L -> 175 -> 102 -> L
    feats[1275] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '175') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '175', '102') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '102', node2.lemma)
    
    # L -> 175 -> 16_du -> L
    feats[1276] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '175') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '175', '16_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_du', node2.lemma)
    
    # L -> 175 -> -131 -> L
    feats[1277] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '175') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '175', '-131') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-131', node2.lemma)
    
    # L -> 175 -> -46 -> L
    feats[1278] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '175') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '175', '-46') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-46', node2.lemma)
    
    # L -> 175 -> 79 -> L
    feats[1279] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '175') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '175', '79') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '79', node2.lemma)
    
    # L -> 175 -> -269 -> L
    feats[1280] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '175') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '175', '-269') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-269', node2.lemma)
    
    # L -> 175 -> -82 -> L
    feats[1281] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '175') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '175', '-82') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-82', node2.lemma)
    
    # L -> 175 -> -299 -> L
    feats[1282] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '175') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '175', '-299') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-299', node2.lemma)
    
    # L -> 175 -> -157 -> L
    feats[1283] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '175') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '175', '-157') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-157', node2.lemma)
    
    # L -> 175 -> -309 -> L
    feats[1284] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '175') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '175', '-309') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-309', node2.lemma)
    
    # L -> 175 -> -169 -> L
    feats[1285] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '175') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '175', '-169') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-169', node2.lemma)
    
    # L -> 72 -> 116 -> L
    feats[1286] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '72') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '72', '116') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '116', node2.lemma)
    
    # L -> 72 -> -263 -> L
    feats[1287] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '72') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '72', '-263') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-263', node2.lemma)
    
    # L -> 72 -> 97 -> L
    feats[1288] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '72') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '72', '97') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '97', node2.lemma)
    
    # L -> 72 -> nom. du. -> L
    feats[1289] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '72') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '72', 'nom. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. du.', node2.lemma)
    
    # L -> 72 -> 170 -> L
    feats[1290] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '72') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '72', '170') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '170', node2.lemma)
    
    # L -> 72 -> -266 -> L
    feats[1291] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '72') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '72', '-266') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-266', node2.lemma)
    
    # L -> 72 -> 161 -> L
    feats[1292] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '72') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '72', '161') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '161', node2.lemma)
    
    # L -> 72 -> -246 -> L
    feats[1293] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '72') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '72', '-246') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-246', node2.lemma)
    
    # L -> 72 -> 131 -> L
    feats[1294] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '72') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '72', '131') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '131', node2.lemma)
    
    # L -> 72 -> 157 -> L
    feats[1295] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '72') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '72', '157') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '157', node2.lemma)
    
    # L -> 72 -> 7_sp -> L
    feats[1296] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '72') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '72', '7_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_sp', node2.lemma)
    
    # L -> 72 -> 15_du -> L
    feats[1297] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '72') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '72', '15_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_du', node2.lemma)
    
    # L -> 72 -> 176 -> L
    feats[1298] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '72') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '72', '176') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '176', node2.lemma)
    
    # L -> 72 -> -157 -> L
    feats[1299] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '72') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '72', '-157') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-157', node2.lemma)
    
    # L -> 72 -> -10 -> L
    feats[1300] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '72') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '72', '-10') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-10', node2.lemma)
    
    # L -> 72 -> -37 -> L
    feats[1301] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '72') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '72', '-37') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-37', node2.lemma)
    
    # L -> 72 -> 98 -> L
    feats[1302] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '72') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '72', '98') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '98', node2.lemma)
    
    # L -> 72 -> nom. adj. -> L
    feats[1303] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '72') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '72', 'nom. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. adj.', node2.lemma)
    
    # L -> 28_tp -> -119 -> L
    feats[1304] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28_tp', '-119') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-119', node2.lemma)
    
    # L -> 28_tp -> 2_fp -> L
    feats[1305] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28_tp', '2_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_fp', node2.lemma)
    
    # L -> 28_tp -> 180 -> L
    feats[1306] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28_tp', '180') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '180', node2.lemma)
    
    # L -> 28_tp -> -144 -> L
    feats[1307] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28_tp', '-144') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-144', node2.lemma)
    
    # L -> 179 -> 4_du -> L
    feats[1308] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '179') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '179', '4_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_du', node2.lemma)
    
    # L -> 179 -> 5_sg -> L
    feats[1309] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '179') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '179', '5_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sg', node2.lemma)
    
    # L -> 179 -> -296 -> L
    feats[1310] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '179') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '179', '-296') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-296', node2.lemma)
    
    # L -> 179 -> -18 -> L
    feats[1311] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '179') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '179', '-18') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-18', node2.lemma)
    
    # L -> 179 -> voc. du. -> L
    feats[1312] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '179') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '179', 'voc. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. du.', node2.lemma)
    
    # L -> 179 -> -147 -> L
    feats[1313] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '179') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '179', '-147') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-147', node2.lemma)
    
    # L -> 179 -> 9_tp -> L
    feats[1314] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '179') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '179', '9_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_tp', node2.lemma)
    
    # L -> 179 -> 136 -> L
    feats[1315] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '179') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '179', '136') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '136', node2.lemma)
    
    # L -> 179 -> -52 -> L
    feats[1316] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '179') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '179', '-52') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-52', node2.lemma)
    
    # L -> nom. fem -> instr. neutr. -> L
    feats[1317] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. fem', 'instr. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. neutr.', node2.lemma)
    
    # L -> nom. fem -> 42 -> L
    feats[1318] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. fem', '42') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '42', node2.lemma)
    
    # L -> nom. fem -> -55 -> L
    feats[1319] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. fem', '-55') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-55', node2.lemma)
    
    # L -> 128 -> -101 -> L
    feats[1320] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '128') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '128', '-101') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-101', node2.lemma)
    
    # L -> 39 -> 16_fp -> L
    feats[1321] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '16_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_fp', node2.lemma)
    
    # L -> 39 -> tp -> L
    feats[1322] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', 'tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'tp', node2.lemma)
    
    # L -> 39 -> 34 -> L
    feats[1323] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '34') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '34', node2.lemma)
    
    # L -> 39 -> pl_sp -> L
    feats[1324] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', 'pl_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_sp', node2.lemma)
    
    # L -> 39 -> 8_fp -> L
    feats[1325] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '8_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    
    # L -> 39 -> 28_sg -> L
    feats[1326] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '28_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28_sg', node2.lemma)
    
    # L -> 39 -> -78 -> L
    feats[1327] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '-78') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-78', node2.lemma)
    
    # L -> 39 -> 12_fp -> L
    feats[1328] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '12_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_fp', node2.lemma)
    
    # L -> 39 -> 27_sg -> L
    feats[1329] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '27_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_sg', node2.lemma)
    
    # L -> 39 -> -276 -> L
    feats[1330] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '-276') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-276', node2.lemma)
    
    # L -> 39 -> -43 -> L
    feats[1331] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '-43') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-43', node2.lemma)
    
    # L -> 39 -> 31 -> L
    feats[1332] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '31') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '31', node2.lemma)
    
    # L -> 39 -> 3_du -> L
    feats[1333] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '3_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    
    # L -> 39 -> -147 -> L
    feats[1334] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '-147') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-147', node2.lemma)
    
    # L -> 39 -> -122 -> L
    feats[1335] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '-122') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-122', node2.lemma)
    
    # L -> 39 -> -13 -> L
    feats[1336] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '-13') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-13', node2.lemma)
    
    # L -> 39 -> 14_sg -> L
    feats[1337] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '14_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_sg', node2.lemma)
    
    # L -> 39 -> 138 -> L
    feats[1338] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '138') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '138', node2.lemma)
    
    # L -> 39 -> -91 -> L
    feats[1339] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '-91') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-91', node2.lemma)
    
    # L -> 39 -> -242 -> L
    feats[1340] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '-242') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-242', node2.lemma)
    
    # L -> 39 -> 76 -> L
    feats[1341] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '76') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '76', node2.lemma)
    
    # L -> 39 -> 2_pl -> L
    feats[1342] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '2_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_pl', node2.lemma)
    
    # L -> 39 -> -283 -> L
    feats[1343] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '-283') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-283', node2.lemma)
    
    # L -> 39 -> 114 -> L
    feats[1344] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '114') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '114', node2.lemma)
    
    # L -> 39 -> pl_tp -> L
    feats[1345] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', 'pl_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_tp', node2.lemma)
    
    # L -> 39 -> instr. neutr. -> L
    feats[1346] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', 'instr. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. neutr.', node2.lemma)
    
    # L -> 39 -> 119 -> L
    feats[1347] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '119') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '119', node2.lemma)
    
    # L -> 39 -> -76 -> L
    feats[1348] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '-76') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-76', node2.lemma)
    
    # L -> 39 -> nom. fem -> L
    feats[1349] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', 'nom. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. fem', node2.lemma)
    
    # L -> 39 -> 128 -> L
    feats[1350] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '128') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '128', node2.lemma)
    
    # L -> 39 -> -220 -> L
    feats[1351] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '-220') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-220', node2.lemma)
    
    # L -> 39 -> 39 -> L
    feats[1352] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '39') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '39', node2.lemma)
    
    # L -> 39 -> -301 -> L
    feats[1353] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '-301') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-301', node2.lemma)
    
    # L -> 39 -> 6_tp -> L
    feats[1354] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '6_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_tp', node2.lemma)
    
    # L -> 39 -> 93 -> L
    feats[1355] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '93') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '93', node2.lemma)
    
    # L -> 39 -> 38 -> L
    feats[1356] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '38') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '38', node2.lemma)
    
    # L -> 39 -> 174 -> L
    feats[1357] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '174') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '174', node2.lemma)
    
    # L -> 39 -> -154 -> L
    feats[1358] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '-154') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-154', node2.lemma)
    
    # L -> 39 -> acc. sg. -> L
    feats[1359] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', 'acc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. sg.', node2.lemma)
    
    # L -> 39 -> 13_tp -> L
    feats[1360] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '13_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_tp', node2.lemma)
    
    # L -> 39 -> -41 -> L
    feats[1361] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '-41') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-41', node2.lemma)
    
    # L -> 39 -> -23 -> L
    feats[1362] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '-23') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-23', node2.lemma)
    
    # L -> 39 -> -19 -> L
    feats[1363] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '-19') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-19', node2.lemma)
    
    # L -> 39 -> instr. masc. -> L
    feats[1364] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', 'instr. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. masc.', node2.lemma)
    
    # L -> 39 -> adj -> L
    feats[1365] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', 'adj') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'adj', node2.lemma)
    
    # L -> 39 -> acc. fem -> L
    feats[1366] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', 'acc. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. fem', node2.lemma)
    
    # L -> 39 -> -31 -> L
    feats[1367] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '-31') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-31', node2.lemma)
    
    # L -> 39 -> -72 -> L
    feats[1368] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '-72') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-72', node2.lemma)
    
    # L -> 39 -> 10_fp -> L
    feats[1369] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '10_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_fp', node2.lemma)
    
    # L -> 39 -> 30_pl -> L
    feats[1370] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '30_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_pl', node2.lemma)
    
    # L -> 39 -> -81 -> L
    feats[1371] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '-81') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-81', node2.lemma)
    
    # L -> 39 -> -96 -> L
    feats[1372] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '-96') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-96', node2.lemma)
    
    # L -> 39 -> 15_sg -> L
    feats[1373] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '15_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_sg', node2.lemma)
    
    # L -> 39 -> acc. masc. -> L
    feats[1374] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', 'acc. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. masc.', node2.lemma)
    
    # L -> 39 -> 10_du -> L
    feats[1375] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '10_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_du', node2.lemma)
    
    # L -> 39 -> -79 -> L
    feats[1376] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '-79') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-79', node2.lemma)
    
    # L -> 39 -> -240 -> L
    feats[1377] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '-240') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-240', node2.lemma)
    
    # L -> 39 -> -54 -> L
    feats[1378] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '-54') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-54', node2.lemma)
    
    # L -> 39 -> 11_sp -> L
    feats[1379] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '11_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_sp', node2.lemma)
    
    # L -> 39 -> -112 -> L
    feats[1380] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', '-112') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    
    # L -> 159 -> 97 -> L
    feats[1381] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '159', '97') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '97', node2.lemma)
    
    # L -> 6_tp -> 82 -> L
    feats[1382] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_tp', '82') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '82', node2.lemma)
    
    # L -> 6_tp -> 98 -> L
    feats[1383] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_tp', '98') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '98', node2.lemma)
    
    # L -> 6_tp -> acc. pl. -> L
    feats[1384] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_tp', 'acc. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. pl.', node2.lemma)
    
    # L -> 6_tp -> 168 -> L
    feats[1385] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_tp', '168') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '168', node2.lemma)
    
    # L -> pl -> 14_pl -> L
    feats[1386] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl', '14_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_pl', node2.lemma)
    
    # L -> -303 -> -297 -> L
    feats[1387] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', '-297') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-297', node2.lemma)
    
    # L -> -303 -> 14_sg -> L
    feats[1388] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', '14_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_sg', node2.lemma)
    
    # L -> -303 -> 9_tp -> L
    feats[1389] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', '9_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_tp', node2.lemma)
    
    # L -> -303 -> -84 -> L
    feats[1390] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', '-84') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-84', node2.lemma)
    
    # L -> -303 -> -242 -> L
    feats[1391] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', '-242') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-242', node2.lemma)
    
    # L -> -303 -> 136 -> L
    feats[1392] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', '136') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '136', node2.lemma)
    
    # L -> -303 -> 54 -> L
    feats[1393] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', '54') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    
    # L -> -303 -> 76 -> L
    feats[1394] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', '76') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '76', node2.lemma)
    
    # L -> -303 -> sp -> L
    feats[1395] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', 'sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sp', node2.lemma)
    
    # L -> -303 -> 9_sp -> L
    feats[1396] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', '9_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_sp', node2.lemma)
    
    # L -> -303 -> -283 -> L
    feats[1397] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', '-283') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-283', node2.lemma)
    
    # L -> -303 -> -132 -> L
    feats[1398] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', '-132') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-132', node2.lemma)
    
    # L -> -303 -> -307 -> L
    feats[1399] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', '-307') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-307', node2.lemma)
    
    # L -> -303 -> 2_sg -> L
    feats[1400] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', '2_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_sg', node2.lemma)
    
    # L -> -303 -> -247 -> L
    feats[1401] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', '-247') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-247', node2.lemma)
    
    # L -> -303 -> acc. du. -> L
    feats[1402] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', 'acc. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. du.', node2.lemma)
    
    # L -> -303 -> 130 -> L
    feats[1403] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', '130') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '130', node2.lemma)
    
    # L -> -303 -> nom. fem -> L
    feats[1404] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', 'nom. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. fem', node2.lemma)
    
    # L -> -303 -> 1 -> L
    feats[1405] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', '1') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '1', node2.lemma)
    
    # L -> -303 -> 27_fp -> L
    feats[1406] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', '27_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_fp', node2.lemma)
    
    # L -> -303 -> fp -> L
    feats[1407] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', 'fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'fp', node2.lemma)
    
    # L -> -303 -> 39 -> L
    feats[1408] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', '39') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '39', node2.lemma)
    
    # L -> -303 -> -99 -> L
    feats[1409] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', '-99') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-99', node2.lemma)
    
    # L -> -303 -> 4_fp -> L
    feats[1410] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', '4_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_fp', node2.lemma)
    
    # L -> -303 -> 13_sp -> L
    feats[1411] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', '13_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_sp', node2.lemma)
    
    # L -> -303 -> -301 -> L
    feats[1412] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', '-301') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-301', node2.lemma)
    
    # L -> -303 -> pl -> L
    feats[1413] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', 'pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    
    # L -> -303 -> 91 -> L
    feats[1414] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', '91') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    
    # L -> -303 -> 93 -> L
    feats[1415] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', '93') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '93', node2.lemma)
    
    # L -> -303 -> 38 -> L
    feats[1416] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', '38') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '38', node2.lemma)
    
    # L -> -303 -> 3_tp -> L
    feats[1417] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', '3_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_tp', node2.lemma)
    
    # L -> -303 -> dat. pl. -> L
    feats[1418] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', 'dat. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. pl.', node2.lemma)
    
    # L -> -303 -> 4_tp -> L
    feats[1419] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', '4_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_tp', node2.lemma)
    
    # L -> -303 -> -23 -> L
    feats[1420] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', '-23') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-23', node2.lemma)
    
    # L -> -303 -> -97 -> L
    feats[1421] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', '-97') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-97', node2.lemma)
    
    # L -> -303 -> sg_tp -> L
    feats[1422] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', 'sg_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_tp', node2.lemma)
    
    # L -> -303 -> 30_fp -> L
    feats[1423] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', '30_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_fp', node2.lemma)
    
    # L -> -303 -> sg -> L
    feats[1424] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', 'sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg', node2.lemma)
    
    # L -> -303 -> 5_sp -> L
    feats[1425] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', '5_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sp', node2.lemma)
    
    # L -> -303 -> adj -> L
    feats[1426] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', 'adj') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'adj', node2.lemma)
    
    # L -> -303 -> neutr -> L
    feats[1427] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', 'neutr') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'neutr', node2.lemma)
    
    # L -> -303 -> 12_sg -> L
    feats[1428] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', '12_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_sg', node2.lemma)
    
    # L -> -303 -> 10_fp -> L
    feats[1429] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', '10_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_fp', node2.lemma)
    
    # L -> -303 -> 42 -> L
    feats[1430] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', '42') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '42', node2.lemma)
    
    # L -> -303 -> -93 -> L
    feats[1431] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', '-93') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-93', node2.lemma)
    
    # L -> -303 -> -96 -> L
    feats[1432] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', '-96') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-96', node2.lemma)
    
    # L -> -303 -> 15_sg -> L
    feats[1433] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', '15_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_sg', node2.lemma)
    
    # L -> -303 -> 10_du -> L
    feats[1434] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', '10_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_du', node2.lemma)
    
    # L -> -303 -> 27_pl -> L
    feats[1435] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', '27_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_pl', node2.lemma)
    
    # L -> -303 -> -240 -> L
    feats[1436] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', '-240') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-240', node2.lemma)
    
    # L -> -303 -> -54 -> L
    feats[1437] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', '-54') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-54', node2.lemma)
    
    # L -> -303 -> 150 -> L
    feats[1438] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', '150') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '150', node2.lemma)
    
    # L -> -303 -> voc. pl. -> L
    feats[1439] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', 'voc. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. pl.', node2.lemma)
    
    # L -> -303 -> 11_fp -> L
    feats[1440] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', '11_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_fp', node2.lemma)
    
    # L -> -303 -> -61 -> L
    feats[1441] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', '-61') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-61', node2.lemma)
    
    # L -> -303 -> -306 -> L
    feats[1442] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', '-306') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-306', node2.lemma)
    
    # L -> -303 -> 116 -> L
    feats[1443] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', '116') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '116', node2.lemma)
    
    # L -> -303 -> 4_sg -> L
    feats[1444] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', '4_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_sg', node2.lemma)
    
    # L -> 3_tp -> 15_tp -> L
    feats[1445] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_tp', '15_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_tp', node2.lemma)
    
    # L -> 3_tp -> -30 -> L
    feats[1446] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_tp', '-30') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-30', node2.lemma)
    
    # L -> 174 -> -114 -> L
    feats[1447] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '174') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '174', '-114') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-114', node2.lemma)
    
    # L -> 174 -> -113 -> L
    feats[1448] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '174') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '174', '-113') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-113', node2.lemma)
    
    # L -> 174 -> 116 -> L
    feats[1449] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '174') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '174', '116') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '116', node2.lemma)
    
    # L -> 174 -> 71 -> L
    feats[1450] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '174') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '174', '71') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '71', node2.lemma)
    
    # L -> 174 -> -159 -> L
    feats[1451] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '174') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '174', '-159') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-159', node2.lemma)
    
    # L -> acc. sg. -> 109 -> L
    feats[1452] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. sg.', '109') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '109', node2.lemma)
    
    # L -> acc. sg. -> 11_du -> L
    feats[1453] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. sg.', '11_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_du', node2.lemma)
    
    # L -> dat. pl. -> -276 -> L
    feats[1454] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', '-276') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-276', node2.lemma)
    
    # L -> dat. pl. -> 5_sp -> L
    feats[1455] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', '5_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sp', node2.lemma)
    
    # L -> 4_tp -> -25 -> L
    feats[1456] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_tp', '-25') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-25', node2.lemma)
    
    # L -> 4_tp -> 153 -> L
    feats[1457] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_tp', '153') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    
    # L -> 4_tp -> -34 -> L
    feats[1458] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_tp', '-34') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-34', node2.lemma)
    
    # L -> 4_tp -> 29_sg -> L
    feats[1459] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_tp', '29_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '29_sg', node2.lemma)
    
    # L -> 4_tp -> fp -> L
    feats[1460] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_tp', 'fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'fp', node2.lemma)
    
    # L -> 4_tp -> pl -> L
    feats[1461] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_tp', 'pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
    
    # L -> 4_tp -> 11_tp -> L
    feats[1462] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_tp', '11_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_tp', node2.lemma)
    
    # L -> 4_tp -> -90 -> L
    feats[1463] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_tp', '-90') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-90', node2.lemma)
    
    # L -> -41 -> -157 -> L
    feats[1464] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-41') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-41', '-157') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-157', node2.lemma)
    
    # L -> -41 -> -64 -> L
    feats[1465] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-41') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-41', '-64') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-64', node2.lemma)
    
    # L -> -41 -> 92 -> L
    feats[1466] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-41') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-41', '92') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '92', node2.lemma)
    
    # L -> -23 -> 3_sp -> L
    feats[1467] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-23') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-23', '3_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_sp', node2.lemma)
    
    # L -> -23 -> -18 -> L
    feats[1468] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-23') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-23', '-18') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-18', node2.lemma)
    
    # L -> -23 -> 5_fp -> L
    feats[1469] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-23') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-23', '5_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_fp', node2.lemma)
    
    # L -> -23 -> -83 -> L
    feats[1470] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-23') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-23', '-83') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-83', node2.lemma)
    
    # L -> 42 -> 8_fp -> L
    feats[1471] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '42') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '42', '8_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    
    # L -> 42 -> 81 -> L
    feats[1472] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '42') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '42', '81') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '81', node2.lemma)
    
    # L -> 58 -> 41 -> L
    feats[1473] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '58', '41') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '41', node2.lemma)
    
    # L -> -54 -> -86 -> L
    feats[1474] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-54') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-54', '-86') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-86', node2.lemma)
    
    # L -> -54 -> -144 -> L
    feats[1475] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-54') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-54', '-144') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-144', node2.lemma)
    
    # L -> 8_sp -> 14_pl -> L
    feats[1476] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_sp', '14_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_pl', node2.lemma)
    
    # L -> 8_sp -> 3_sp -> L
    feats[1477] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_sp', '3_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_sp', node2.lemma)
    
    # L -> 8_sp -> loc -> L
    feats[1478] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_sp', 'loc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc', node2.lemma)
    
    # L -> 8_sp -> -141 -> L
    feats[1479] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_sp', '-141') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-141', node2.lemma)
    
    # L -> 8_sp -> 34 -> L
    feats[1480] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_sp', '34') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '34', node2.lemma)
    
    # L -> 8_sp -> -83 -> L
    feats[1481] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_sp', '-83') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-83', node2.lemma)
    
    # L -> 8_sp -> 8_fp -> L
    feats[1482] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_sp', '8_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    
    # L -> 8_sp -> 28_sg -> L
    feats[1483] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_sp', '28_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28_sg', node2.lemma)
    
    # L -> 8_sp -> -276 -> L
    feats[1484] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_sp', '-276') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-276', node2.lemma)
    
    # L -> 8_sp -> abl. pl. -> L
    feats[1485] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_sp', 'abl. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. pl.', node2.lemma)
    
    # L -> 8_sp -> 31 -> L
    feats[1486] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_sp', '31') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '31', node2.lemma)
    
    # L -> 8_sp -> 3_du -> L
    feats[1487] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_sp', '3_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    
    # L -> 8_sp -> -48 -> L
    feats[1488] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_sp', '-48') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-48', node2.lemma)
    
    # L -> 8_sp -> -122 -> L
    feats[1489] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_sp', '-122') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-122', node2.lemma)
    
    # L -> 8_sp -> -89 -> L
    feats[1490] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_sp', '-89') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-89', node2.lemma)
    
    # L -> 8_sp -> -242 -> L
    feats[1491] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_sp', '-242') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-242', node2.lemma)
    
    # L -> 8_sp -> 9_sp -> L
    feats[1492] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_sp', '9_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_sp', node2.lemma)
    
    # L -> 8_sp -> -247 -> L
    feats[1493] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_sp', '-247') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-247', node2.lemma)
    
    # L -> 8_sp -> nom. fem -> L
    feats[1494] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_sp', 'nom. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. fem', node2.lemma)
    
    # L -> 8_sp -> -220 -> L
    feats[1495] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_sp', '-220') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-220', node2.lemma)
    
    # L -> 8_sp -> 2_du -> L
    feats[1496] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_sp', '2_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_du', node2.lemma)
    
    # L -> 8_sp -> 13_sp -> L
    feats[1497] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_sp', '13_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_sp', node2.lemma)
    
    # L -> 8_sp -> -301 -> L
    feats[1498] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_sp', '-301') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-301', node2.lemma)
    
    # L -> 50 -> 134 -> L
    feats[1499] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '50') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '50', '134') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '134', node2.lemma)
    
    feats[feats < 1e-25] = 1e-25
    return -np.log10(feats)
    