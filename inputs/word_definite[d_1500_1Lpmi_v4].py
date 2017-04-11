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
    
    # L->-154->L
    feats[0] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-154') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-154', node2.lemma)
            
    # L->acc. masc.->L
    feats[1] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. masc.', node2.lemma)
            
    # L->loc. pl.->L
    feats[2] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc. pl.', node2.lemma)
            
    # L->11_fp->L
    feats[3] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_fp', node2.lemma)
            
    # L->-79->L
    feats[4] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-79') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-79', node2.lemma)
            
    # L->-94->C
    feats[5] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-94') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-94', node2.cng)
            
    # L->-122->C
    feats[6] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-122') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-122', node2.cng)
            
    # L->110->C
    feats[7] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '110') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '110', node2.cng)
            
    # L->16_pl->C
    feats[8] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_pl', node2.cng)
            
    # L->nom. fem->C
    feats[9] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. fem', node2.cng)
            
    # L->loc. sg.->C
    feats[10] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', node2.cng)
            
    # L->-17->C
    feats[11] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-17') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-17', node2.cng)
            
    # L->14_sp->C
    feats[12] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_sp', node2.cng)
            
    # L->6_pl->C
    feats[13] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_pl', node2.cng)
            
    # L->10_fp->C
    feats[14] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_fp', node2.cng)
            
    # L->-142->C
    feats[15] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-142') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-142', node2.cng)
            
    # L->-22->C
    feats[16] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-22') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-22', node2.cng)
            
    # L->177->C
    feats[17] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '177') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '177', node2.cng)
            
    # L->-126->C
    feats[18] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-126') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-126', node2.cng)
            
    # L->75->C
    feats[19] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '75') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '75', node2.cng)
            
    # L->sg_fp->C
    feats[20] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg_fp', node2.cng)
            
    # L->90->C
    feats[21] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '90', node2.cng)
            
    # L->115->C
    feats[22] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '115') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '115', node2.cng)
            
    # L->149->C
    feats[23] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', node2.cng)
            
    # L->-14->C
    feats[24] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-14') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-14', node2.cng)
            
    # L->140->C
    feats[25] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '140') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '140', node2.cng)
            
    # L->13_tp->C
    feats[26] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_tp', node2.cng)
            
    # L->92->C
    feats[27] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '92', node2.cng)
            
    # L->3_du->C
    feats[28] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_du', node2.cng)
            
    # L->6_fp->C
    feats[29] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_fp', node2.cng)
            
    # L->-18->C
    feats[30] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', node2.cng)
            
    # L->instr. adj.->C
    feats[31] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. adj.', node2.cng)
            
    # L->voc. pl.->C
    feats[32] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. pl.', node2.cng)
            
    # L->4_tp->C
    feats[33] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_tp', node2.cng)
            
    # L -> -161 -> voc. sg. -> L
    feats[34] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-161') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-161', 'voc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. sg.', node2.lemma)
    
    # L -> -161 -> -113 -> L
    feats[35] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-161') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-161', '-113') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-113', node2.lemma)
    
    # L -> -161 -> dat. du. -> L
    feats[36] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-161') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-161', 'dat. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
    
    # L -> -161 -> -24 -> L
    feats[37] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-161') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-161', '-24') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    
    # L -> -161 -> -247 -> L
    feats[38] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-161') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-161', '-247') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-247', node2.lemma)
    
    # L -> -161 -> -79 -> L
    feats[39] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-161') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-161', '-79') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-79', node2.lemma)
    
    # L -> -161 -> 152 -> L
    feats[40] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-161') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-161', '152') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '152', node2.lemma)
    
    # L -> -161 -> 16_fp -> L
    feats[41] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-161') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-161', '16_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_fp', node2.lemma)
    
    # L -> -161 -> -133 -> L
    feats[42] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-161') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-161', '-133') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-133', node2.lemma)
    
    # L -> -161 -> 58 -> L
    feats[43] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-161') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-161', '58') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '58', node2.lemma)
    
    # L -> pl_fp -> -71 -> L
    feats[44] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_fp', '-71') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-71', node2.lemma)
    
    # L -> pl_fp -> loc. sg. -> L
    feats[45] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_fp', 'loc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc. sg.', node2.lemma)
    
    # L -> pl_fp -> 170 -> L
    feats[46] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_fp', '170') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '170', node2.lemma)
    
    # L -> pl_fp -> -59 -> L
    feats[47] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_fp', '-59') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-59', node2.lemma)
    
    # L -> pl_fp -> 138 -> L
    feats[48] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_fp', '138') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '138', node2.lemma)
    
    # L -> pl_fp -> 10_fp -> L
    feats[49] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_fp', '10_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_fp', node2.lemma)
    
    # L -> pl_fp -> -291 -> L
    feats[50] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_fp', '-291') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-291', node2.lemma)
    
    # L -> pl_fp -> 27_sg -> L
    feats[51] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_fp', '27_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_sg', node2.lemma)
    
    # L -> pl_fp -> -93 -> L
    feats[52] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_fp', '-93') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-93', node2.lemma)
    
    # L -> pl_fp -> 12_tp -> L
    feats[53] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_fp', '12_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_tp', node2.lemma)
    
    # L -> pl_fp -> -58 -> L
    feats[54] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_fp', '-58') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-58', node2.lemma)
    
    # L -> pl_fp -> 2_pl -> L
    feats[55] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_fp', '2_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_pl', node2.lemma)
    
    # L -> pl_fp -> 115 -> L
    feats[56] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_fp', '115') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '115', node2.lemma)
    
    # L -> pl_fp -> 30_sg -> L
    feats[57] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_fp', '30_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_sg', node2.lemma)
    
    # L -> pl_fp -> voc. sg. -> L
    feats[58] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_fp', 'voc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. sg.', node2.lemma)
    
    # L -> pl_fp -> 136 -> L
    feats[59] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_fp', '136') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '136', node2.lemma)
    
    # L -> pl_fp -> 171 -> L
    feats[60] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_fp', '171') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '171', node2.lemma)
    
    # L -> pl_fp -> 37 -> L
    feats[61] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_fp', '37') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '37', node2.lemma)
    
    # L -> pl_fp -> 114 -> L
    feats[62] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_fp', '114') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '114', node2.lemma)
    
    # L -> pl_fp -> 7_du -> L
    feats[63] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_fp', '7_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_du', node2.lemma)
    
    # L -> pl_fp -> -27 -> L
    feats[64] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_fp', '-27') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-27', node2.lemma)
    
    # L -> 14_sp -> 132 -> L
    feats[65] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_sp', '132') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '132', node2.lemma)
    
    # L -> 14_sp -> -297 -> L
    feats[66] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_sp', '-297') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-297', node2.lemma)
    
    # L -> 14_sp -> -271 -> L
    feats[67] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_sp', '-271') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-271', node2.lemma)
    
    # L -> 14_sp -> -245 -> L
    feats[68] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_sp', '-245') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-245', node2.lemma)
    
    # L -> 14_sp -> 98 -> L
    feats[69] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_sp', '98') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '98', node2.lemma)
    
    # L -> 14_sp -> 2_fp -> L
    feats[70] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_sp', '2_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_fp', node2.lemma)
    
    # L -> 14_sp -> 7_tp -> L
    feats[71] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_sp', '7_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_tp', node2.lemma)
    
    # L -> 14_sp -> -126 -> L
    feats[72] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_sp', '-126') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-126', node2.lemma)
    
    # L -> 14_sp -> -132 -> L
    feats[73] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_sp', '-132') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-132', node2.lemma)
    
    # L -> 14_sp -> 90 -> L
    feats[74] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_sp', '90') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '90', node2.lemma)
    
    # L -> 14_sp -> 120 -> L
    feats[75] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_sp', '120') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '120', node2.lemma)
    
    # L -> 14_sp -> 156 -> L
    feats[76] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_sp', '156') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '156', node2.lemma)
    
    # L -> 14_sp -> 30_sg -> L
    feats[77] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_sp', '30_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_sg', node2.lemma)
    
    # L -> 14_sp -> fp -> L
    feats[78] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_sp', 'fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'fp', node2.lemma)
    
    # L -> 14_sp -> -279 -> L
    feats[79] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_sp', '-279') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-279', node2.lemma)
    
    # L -> 14_sp -> -263 -> L
    feats[80] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_sp', '-263') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-263', node2.lemma)
    
    # L -> 14_sp -> 37 -> L
    feats[81] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_sp', '37') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '37', node2.lemma)
    
    # L -> 14_sp -> 7_du -> L
    feats[82] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_sp', '7_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_du', node2.lemma)
    
    # L -> 14_sp -> -36 -> L
    feats[83] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_sp', '-36') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-36', node2.lemma)
    
    # L -> 14_sp -> -156 -> L
    feats[84] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_sp', '-156') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-156', node2.lemma)
    
    # L -> 14_sp -> -143 -> L
    feats[85] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_sp', '-143') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-143', node2.lemma)
    
    # L -> 14_sp -> 142 -> L
    feats[86] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_sp', '142') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '142', node2.lemma)
    
    # L -> 14_sp -> 9_fp -> L
    feats[87] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_sp', '9_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_fp', node2.lemma)
    
    # L -> 14_sp -> 3 -> L
    feats[88] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_sp', '3') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3', node2.lemma)
    
    # L -> 14_sp -> -303 -> L
    feats[89] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_sp', '-303') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-303', node2.lemma)
    
    # L -> 77 -> 2_sg -> L
    feats[90] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '77') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '77', '2_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_sg', node2.lemma)
    
    # L -> 77 -> 14_sp -> L
    feats[91] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '77') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '77', '14_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_sp', node2.lemma)
    
    # L -> 77 -> 10_tp -> L
    feats[92] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '77') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '77', '10_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_tp', node2.lemma)
    
    # L -> 77 -> 30 -> L
    feats[93] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '77') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '77', '30') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30', node2.lemma)
    
    # L -> 77 -> 11_fp -> L
    feats[94] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '77') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '77', '11_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_fp', node2.lemma)
    
    # L -> 77 -> 11_du -> L
    feats[95] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '77') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '77', '11_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_du', node2.lemma)
    
    # L -> 77 -> 7_sp -> L
    feats[96] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '77') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '77', '7_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_sp', node2.lemma)
    
    # L -> 79 -> 151 -> L
    feats[97] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '151') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '151', node2.lemma)
    
    # L -> 79 -> 6_sg -> L
    feats[98] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '6_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_sg', node2.lemma)
    
    # L -> 79 -> 29 -> L
    feats[99] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '29') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '29', node2.lemma)
    
    # L -> 79 -> 33 -> L
    feats[100] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '33') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '33', node2.lemma)
    
    # L -> 79 -> -68 -> L
    feats[101] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '-68') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-68', node2.lemma)
    
    # L -> 79 -> 3_pl -> L
    feats[102] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '3_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_pl', node2.lemma)
    
    # L -> 79 -> 162 -> L
    feats[103] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '162') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '162', node2.lemma)
    
    # L -> 79 -> 12_sp -> L
    feats[104] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '12_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_sp', node2.lemma)
    
    # L -> 79 -> 15_fp -> L
    feats[105] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '15_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_fp', node2.lemma)
    
    # L -> 79 -> 176 -> L
    feats[106] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '176') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '176', node2.lemma)
    
    # L -> 79 -> nom. masc. -> L
    feats[107] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', 'nom. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. masc.', node2.lemma)
    
    # L -> 79 -> 8_pl -> L
    feats[108] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '8_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_pl', node2.lemma)
    
    # L -> 79 -> -28 -> L
    feats[109] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '-28') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-28', node2.lemma)
    
    # L -> 79 -> -303 -> L
    feats[110] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', '-303') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-303', node2.lemma)
    
    # L -> 34 -> abl. du. -> L
    feats[111] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '34', 'abl. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. du.', node2.lemma)
    
    # L -> 34 -> -147 -> L
    feats[112] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '34', '-147') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-147', node2.lemma)
    
    # L -> 34 -> -159 -> L
    feats[113] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '34', '-159') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-159', node2.lemma)
    
    # L -> 34 -> -296 -> L
    feats[114] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '34', '-296') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-296', node2.lemma)
    
    # L -> 34 -> -86 -> L
    feats[115] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '34', '-86') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-86', node2.lemma)
    
    # L -> 34 -> 27_fp -> L
    feats[116] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '34', '27_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_fp', node2.lemma)
    
    # L -> 34 -> 13_sg -> L
    feats[117] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '34', '13_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_sg', node2.lemma)
    
    # L -> 34 -> du_fp -> L
    feats[118] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '34', 'du_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_fp', node2.lemma)
    
    # L -> 34 -> -309 -> L
    feats[119] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '34', '-309') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-309', node2.lemma)
    
    # L -> 34 -> nom. fem -> L
    feats[120] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '34', 'nom. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. fem', node2.lemma)
    
    # L -> 34 -> 128 -> L
    feats[121] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '34', '128') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '128', node2.lemma)
    
    # L -> 34 -> -48 -> L
    feats[122] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '34', '-48') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-48', node2.lemma)
    
    # L -> 34 -> -261 -> L
    feats[123] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '34', '-261') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-261', node2.lemma)
    
    # L -> 34 -> 5_sp -> L
    feats[124] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '34', '5_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sp', node2.lemma)
    
    # L -> 34 -> -17 -> L
    feats[125] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '34', '-17') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-17', node2.lemma)
    
    # L -> 34 -> 150 -> L
    feats[126] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '34', '150') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '150', node2.lemma)
    
    # L -> 34 -> -161 -> L
    feats[127] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '34', '-161') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-161', node2.lemma)
    
    # L -> 34 -> -59 -> L
    feats[128] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '34', '-59') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-59', node2.lemma)
    
    # L -> 34 -> -149 -> L
    feats[129] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '34', '-149') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-149', node2.lemma)
    
    # L -> 34 -> -271 -> L
    feats[130] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '34', '-271') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-271', node2.lemma)
    
    # L -> 34 -> 10_sg -> L
    feats[131] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '34', '10_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_sg', node2.lemma)
    
    # L -> 34 -> -245 -> L
    feats[132] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '34', '-245') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-245', node2.lemma)
    
    # L -> 34 -> dat -> L
    feats[133] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '34', 'dat') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat', node2.lemma)
    
    # L -> 34 -> -25 -> L
    feats[134] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '34', '-25') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-25', node2.lemma)
    
    # L -> 34 -> 98 -> L
    feats[135] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '34', '98') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '98', node2.lemma)
    
    # L -> 34 -> 5_pl -> L
    feats[136] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '34', '5_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_pl', node2.lemma)
    
    # L -> 34 -> 2_sp -> L
    feats[137] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '34', '2_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_sp', node2.lemma)
    
    # L -> -151 -> 114 -> L
    feats[138] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-151') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-151', '114') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '114', node2.lemma)
    
    # L -> 179 -> 180 -> L
    feats[139] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '179') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '179', '180') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '180', node2.lemma)
    
    # L -> 179 -> 2_tp -> L
    feats[140] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '179') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '179', '2_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_tp', node2.lemma)
    
    # L -> 179 -> 91 -> L
    feats[141] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '179') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '179', '91') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    
    # L -> 179 -> 4_tp -> L
    feats[142] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '179') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '179', '4_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_tp', node2.lemma)
    
    # L -> 179 -> 6_sg -> L
    feats[143] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '179') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '179', '6_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_sg', node2.lemma)
    
    # L -> 179 -> -37 -> L
    feats[144] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '179') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '179', '-37') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-37', node2.lemma)
    
    # L -> 179 -> -131 -> L
    feats[145] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '179') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '179', '-131') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-131', node2.lemma)
    
    # L -> 179 -> 172 -> L
    feats[146] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '179') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '179', '172') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '172', node2.lemma)
    
    # L -> 179 -> 122 -> L
    feats[147] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '179') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '179', '122') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '122', node2.lemma)
    
    # L -> 179 -> 3_sp -> L
    feats[148] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '179') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '179', '3_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_sp', node2.lemma)
    
    # L -> 179 -> 74 -> L
    feats[149] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '179') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '179', '74') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '74', node2.lemma)
    
    # L -> 179 -> 55 -> L
    feats[150] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '179') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '179', '55') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '55', node2.lemma)
    
    # L -> 179 -> -78 -> L
    feats[151] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '179') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '179', '-78') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-78', node2.lemma)
    
    # L -> 179 -> -141 -> L
    feats[152] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '179') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '179', '-141') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-141', node2.lemma)
    
    # L -> 179 -> 2_sp -> L
    feats[153] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '179') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '179', '2_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_sp', node2.lemma)
    
    # L -> 179 -> 15_fp -> L
    feats[154] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '179') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '179', '15_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_fp', node2.lemma)
    
    # L -> 179 -> 13_fp -> L
    feats[155] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '179') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '179', '13_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_fp', node2.lemma)
    
    # L -> 179 -> -64 -> L
    feats[156] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '179') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '179', '-64') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-64', node2.lemma)
    
    # L -> 179 -> 8_pl -> L
    feats[157] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '179') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '179', '8_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_pl', node2.lemma)
    
    # L -> 179 -> 1 -> L
    feats[158] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '179') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '179', '1') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '1', node2.lemma)
    
    # L -> 6_pl -> -150 -> L
    feats[159] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_pl', '-150') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-150', node2.lemma)
    
    # L -> 6_pl -> 4_pl -> L
    feats[160] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_pl', '4_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_pl', node2.lemma)
    
    # L -> 6_pl -> 3_tp -> L
    feats[161] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_pl', '3_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_tp', node2.lemma)
    
    # L -> 6_pl -> -246 -> L
    feats[162] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_pl', '-246') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-246', node2.lemma)
    
    # L -> 6_pl -> voc. masc. -> L
    feats[163] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_pl', 'voc. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. masc.', node2.lemma)
    
    # L -> 6_pl -> 2_sg -> L
    feats[164] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_pl', '2_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_sg', node2.lemma)
    
    # L -> 6_pl -> 11_sp -> L
    feats[165] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_pl', '11_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_sp', node2.lemma)
    
    # L -> 6_pl -> du -> L
    feats[166] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_pl', 'du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du', node2.lemma)
    
    # L -> 6_pl -> 13_sg -> L
    feats[167] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_pl', '13_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_sg', node2.lemma)
    
    # L -> 6_pl -> 110 -> L
    feats[168] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_pl', '110') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '110', node2.lemma)
    
    # L -> 6_pl -> -10 -> L
    feats[169] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_pl', '-10') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-10', node2.lemma)
    
    # L -> 6_pl -> du_fp -> L
    feats[170] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_pl', 'du_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_fp', node2.lemma)
    
    # L -> 6_pl -> 16_pl -> L
    feats[171] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_pl', '16_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_pl', node2.lemma)
    
    # L -> 6_pl -> -220 -> L
    feats[172] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_pl', '-220') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-220', node2.lemma)
    
    # L -> 6_pl -> -169 -> L
    feats[173] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_pl', '-169') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-169', node2.lemma)
    
    # L -> 6_pl -> -48 -> L
    feats[174] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_pl', '-48') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-48', node2.lemma)
    
    # L -> 6_pl -> loc. sg. -> L
    feats[175] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_pl', 'loc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc. sg.', node2.lemma)
    
    # L -> 6_pl -> sg_tp -> L
    feats[176] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_pl', 'sg_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_tp', node2.lemma)
    
    # L -> 6_pl -> -261 -> L
    feats[177] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_pl', '-261') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-261', node2.lemma)
    
    # L -> 6_pl -> -17 -> L
    feats[178] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_pl', '-17') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-17', node2.lemma)
    
    # L -> 6_pl -> 14_sp -> L
    feats[179] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_pl', '14_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_sp', node2.lemma)
    
    # L -> 6_pl -> 34 -> L
    feats[180] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_pl', '34') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '34', node2.lemma)
    
    # L -> 6_pl -> -151 -> L
    feats[181] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_pl', '-151') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-151', node2.lemma)
    
    # L -> 6_pl -> -152 -> L
    feats[182] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_pl', '-152') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-152', node2.lemma)
    
    # L -> 6_pl -> acc. neutr. -> L
    feats[183] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_pl', 'acc. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. neutr.', node2.lemma)
    
    # L -> 6_pl -> tp -> L
    feats[184] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_pl', 'tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'tp', node2.lemma)
    
    # L -> 6_pl -> 10_fp -> L
    feats[185] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_pl', '10_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_fp', node2.lemma)
    
    # L -> 6_pl -> -271 -> L
    feats[186] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_pl', '-271') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-271', node2.lemma)
    
    # L -> 6_pl -> 10_tp -> L
    feats[187] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_pl', '10_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_tp', node2.lemma)
    
    # L -> 6_pl -> -77 -> L
    feats[188] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_pl', '-77') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-77', node2.lemma)
    
    # L -> 6_pl -> -166 -> L
    feats[189] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_pl', '-166') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-166', node2.lemma)
    
    # L -> 6_pl -> dat -> L
    feats[190] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_pl', 'dat') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat', node2.lemma)
    
    # L -> 6_pl -> -25 -> L
    feats[191] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_pl', '-25') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-25', node2.lemma)
    
    # L -> 6_pl -> -93 -> L
    feats[192] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_pl', '-93') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-93', node2.lemma)
    
    # L -> 6_pl -> -139 -> L
    feats[193] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_pl', '-139') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-139', node2.lemma)
    
    # L -> 6_pl -> instr. sg. -> L
    feats[194] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_pl', 'instr. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. sg.', node2.lemma)
    
    # L -> 6_pl -> 114 -> L
    feats[195] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_pl', '114') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '114', node2.lemma)
    
    # L -> 6_pl -> -143 -> L
    feats[196] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_pl', '-143') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-143', node2.lemma)
    
    # L -> 6_pl -> nom. du. -> L
    feats[197] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_pl', 'nom. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. du.', node2.lemma)
    
    # L -> 6_pl -> 11_tp -> L
    feats[198] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_pl', '11_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_tp', node2.lemma)
    
    # L -> 6_pl -> loc. pl. -> L
    feats[199] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_pl', 'loc. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc. pl.', node2.lemma)
    
    # L -> 6_pl -> -47 -> L
    feats[200] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_pl', '-47') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-47', node2.lemma)
    
    # L -> 6_pl -> instr. masc. -> L
    feats[201] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_pl', 'instr. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. masc.', node2.lemma)
    
    # L -> 6_pl -> 119 -> L
    feats[202] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_pl', '119') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '119', node2.lemma)
    
    # L -> 6_pl -> 51 -> L
    feats[203] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_pl', '51') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '51', node2.lemma)
    
    # L -> 6_pl -> -26 -> L
    feats[204] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_pl', '-26') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-26', node2.lemma)
    
    # L -> -157 -> 81 -> L
    feats[205] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '81') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '81', node2.lemma)
    
    # L -> -157 -> 95 -> L
    feats[206] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '95') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '95', node2.lemma)
    
    # L -> -157 -> -50 -> L
    feats[207] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '-50') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-50', node2.lemma)
    
    # L -> -157 -> -200 -> L
    feats[208] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '-200') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-200', node2.lemma)
    
    # L -> -157 -> -99 -> L
    feats[209] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '-99') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-99', node2.lemma)
    
    # L -> -157 -> -293 -> L
    feats[210] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '-293') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-293', node2.lemma)
    
    # L -> -157 -> 5_tp -> L
    feats[211] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '5_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_tp', node2.lemma)
    
    # L -> -157 -> 89 -> L
    feats[212] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '89') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '89', node2.lemma)
    
    # L -> -157 -> 7_fp -> L
    feats[213] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '7_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_fp', node2.lemma)
    
    # L -> -157 -> 5_fp -> L
    feats[214] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '5_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_fp', node2.lemma)
    
    # L -> -157 -> 139 -> L
    feats[215] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '139') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '139', node2.lemma)
    
    # L -> -157 -> -276 -> L
    feats[216] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '-276') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-276', node2.lemma)
    
    # L -> -157 -> 11_fp -> L
    feats[217] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '11_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_fp', node2.lemma)
    
    # L -> -157 -> -247 -> L
    feats[218] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '-247') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-247', node2.lemma)
    
    # L -> -157 -> -79 -> L
    feats[219] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '-79') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-79', node2.lemma)
    
    # L -> -157 -> 152 -> L
    feats[220] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '152') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '152', node2.lemma)
    
    # L -> -157 -> -102 -> L
    feats[221] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '-102') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-102', node2.lemma)
    
    # L -> -157 -> 97 -> L
    feats[222] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '97') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '97', node2.lemma)
    
    # L -> -157 -> -153 -> L
    feats[223] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '-153') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-153', node2.lemma)
    
    # L -> -157 -> -76 -> L
    feats[224] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '-76') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-76', node2.lemma)
    
    # L -> -157 -> acc -> L
    feats[225] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', 'acc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc', node2.lemma)
    
    # L -> -157 -> -26 -> L
    feats[226] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '-26') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-26', node2.lemma)
    
    # L -> -157 -> -133 -> L
    feats[227] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '-133') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-133', node2.lemma)
    
    # L -> -157 -> 78 -> L
    feats[228] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '78') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '78', node2.lemma)
    
    # L -> -157 -> -96 -> L
    feats[229] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '-96') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-96', node2.lemma)
    
    # L -> -157 -> 109 -> L
    feats[230] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '109') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '109', node2.lemma)
    
    # L -> -157 -> acc. pl. -> L
    feats[231] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', 'acc. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. pl.', node2.lemma)
    
    # L -> -157 -> gen. sg. -> L
    feats[232] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', 'gen. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen. sg.', node2.lemma)
    
    # L -> -157 -> instr. fem -> L
    feats[233] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', 'instr. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. fem', node2.lemma)
    
    # L -> -157 -> 102 -> L
    feats[234] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '102') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '102', node2.lemma)
    
    # L -> -157 -> -66 -> L
    feats[235] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '-66') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-66', node2.lemma)
    
    # L -> -157 -> 8_sp -> L
    feats[236] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '8_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_sp', node2.lemma)
    
    # L -> -157 -> 148 -> L
    feats[237] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '148') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '148', node2.lemma)
    
    # L -> -157 -> 60 -> L
    feats[238] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '60') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '60', node2.lemma)
    
    # L -> -157 -> 181 -> L
    feats[239] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '181') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '181', node2.lemma)
    
    # L -> -157 -> -210 -> L
    feats[240] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '-210') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-210', node2.lemma)
    
    # L -> -157 -> -55 -> L
    feats[241] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '-55') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-55', node2.lemma)
    
    # L -> -157 -> 68 -> L
    feats[242] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '68') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '68', node2.lemma)
    
    # L -> -157 -> 93 -> L
    feats[243] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '93') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '93', node2.lemma)
    
    # L -> -157 -> 6_du -> L
    feats[244] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '6_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_du', node2.lemma)
    
    # L -> -157 -> 71 -> L
    feats[245] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '71') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '71', node2.lemma)
    
    # L -> -157 -> 15_tp -> L
    feats[246] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '15_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_tp', node2.lemma)
    
    # L -> -157 -> -32 -> L
    feats[247] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '-32') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-32', node2.lemma)
    
    # L -> -157 -> instr. neutr. -> L
    feats[248] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', 'instr. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. neutr.', node2.lemma)
    
    # L -> -157 -> -54 -> L
    feats[249] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '-54') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-54', node2.lemma)
    
    # L -> -157 -> -81 -> L
    feats[250] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '-81') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-81', node2.lemma)
    
    # L -> -157 -> -190 -> L
    feats[251] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '-190') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-190', node2.lemma)
    
    # L -> -157 -> 157 -> L
    feats[252] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '157') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '157', node2.lemma)
    
    # L -> -157 -> -299 -> L
    feats[253] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '-299') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-299', node2.lemma)
    
    # L -> -157 -> 4_sg -> L
    feats[254] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '4_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_sg', node2.lemma)
    
    # L -> -157 -> 35 -> L
    feats[255] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '35') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '35', node2.lemma)
    
    # L -> -157 -> -18 -> L
    feats[256] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-157', '-18') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-18', node2.lemma)
    
    # L -> -59 -> 100 -> L
    feats[257] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-59') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-59', '100') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '100', node2.lemma)
    
    # L -> -273 -> 54 -> L
    feats[258] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '54') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    
    # L -> -273 -> -271 -> L
    feats[259] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '-271') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-271', node2.lemma)
    
    # L -> -273 -> 10_sg -> L
    feats[260] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '10_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_sg', node2.lemma)
    
    # L -> -273 -> -81 -> L
    feats[261] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '-81') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-81', node2.lemma)
    
    # L -> -273 -> -43 -> L
    feats[262] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '-43') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-43', node2.lemma)
    
    # L -> -273 -> 30_fp -> L
    feats[263] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '30_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_fp', node2.lemma)
    
    # L -> -273 -> -299 -> L
    feats[264] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '-299') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-299', node2.lemma)
    
    # L -> -273 -> -18 -> L
    feats[265] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '-18') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-18', node2.lemma)
    
    # L -> -273 -> 7_sp -> L
    feats[266] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '7_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_sp', node2.lemma)
    
    # L -> -273 -> 155 -> L
    feats[267] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '155') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '155', node2.lemma)
    
    # L -> -273 -> 108 -> L
    feats[268] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '108') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '108', node2.lemma)
    
    # L -> -273 -> 141 -> L
    feats[269] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '141') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '141', node2.lemma)
    
    # L -> -273 -> -35 -> L
    feats[270] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '-35') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    
    # L -> -273 -> -117 -> L
    feats[271] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '-117') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-117', node2.lemma)
    
    # L -> -273 -> 9_pl -> L
    feats[272] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '9_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_pl', node2.lemma)
    
    # L -> -273 -> 91 -> L
    feats[273] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '91') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    
    # L -> -273 -> 8_du -> L
    feats[274] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '8_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_du', node2.lemma)
    
    # L -> -273 -> -308 -> L
    feats[275] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '-308') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-308', node2.lemma)
    
    # L -> -273 -> -37 -> L
    feats[276] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '-37') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-37', node2.lemma)
    
    # L -> -273 -> 131 -> L
    feats[277] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '131') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '131', node2.lemma)
    
    # L -> -273 -> -69 -> L
    feats[278] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '-69') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-69', node2.lemma)
    
    # L -> -273 -> 101 -> L
    feats[279] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '101') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    
    # L -> -273 -> 160 -> L
    feats[280] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '160') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    
    # L -> -273 -> 3_pl -> L
    feats[281] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '3_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_pl', node2.lemma)
    
    # L -> -273 -> -90 -> L
    feats[282] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '-90') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-90', node2.lemma)
    
    # L -> -273 -> 10_pl -> L
    feats[283] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '10_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_pl', node2.lemma)
    
    # L -> -273 -> 3_sp -> L
    feats[284] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '3_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_sp', node2.lemma)
    
    # L -> -273 -> 55 -> L
    feats[285] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '55') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '55', node2.lemma)
    
    # L -> -273 -> 27_du -> L
    feats[286] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '27_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_du', node2.lemma)
    
    # L -> -273 -> -141 -> L
    feats[287] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '-141') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-141', node2.lemma)
    
    # L -> -273 -> 13_fp -> L
    feats[288] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '13_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_fp', node2.lemma)
    
    # L -> -273 -> -129 -> L
    feats[289] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '-129') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-129', node2.lemma)
    
    # L -> -273 -> 100 -> L
    feats[290] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '100') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '100', node2.lemma)
    
    # L -> -273 -> -28 -> L
    feats[291] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '-28') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-28', node2.lemma)
    
    # L -> -273 -> 1 -> L
    feats[292] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '1') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '1', node2.lemma)
    
    # L -> -273 -> -303 -> L
    feats[293] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', '-303') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-303', node2.lemma)
    
    # L -> voc. fem -> -150 -> L
    feats[294] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. fem', '-150') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-150', node2.lemma)
    
    # L -> voc. fem -> 28_sg -> L
    feats[295] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. fem', '28_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28_sg', node2.lemma)
    
    # L -> voc. fem -> abl. pl. -> L
    feats[296] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. fem', 'abl. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. pl.', node2.lemma)
    
    # L -> voc. fem -> -15 -> L
    feats[297] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. fem', '-15') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-15', node2.lemma)
    
    # L -> voc. fem -> -246 -> L
    feats[298] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. fem', '-246') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-246', node2.lemma)
    
    # L -> voc. fem -> voc. du. -> L
    feats[299] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. fem', 'voc. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. du.', node2.lemma)
    
    # L -> voc. fem -> voc. masc. -> L
    feats[300] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. fem', 'voc. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. masc.', node2.lemma)
    
    # L -> voc. fem -> -163 -> L
    feats[301] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. fem', '-163') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-163', node2.lemma)
    
    # L -> voc. fem -> -122 -> L
    feats[302] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. fem', '-122') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-122', node2.lemma)
    
    # L -> voc. fem -> 159 -> L
    feats[303] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. fem', '159') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '159', node2.lemma)
    
    # L -> voc. fem -> 110 -> L
    feats[304] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. fem', '110') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '110', node2.lemma)
    
    # L -> voc. fem -> 12_fp -> L
    feats[305] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. fem', '12_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_fp', node2.lemma)
    
    # L -> voc. fem -> -220 -> L
    feats[306] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. fem', '-220') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-220', node2.lemma)
    
    # L -> voc. fem -> 137 -> L
    feats[307] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. fem', '137') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    
    # L -> voc. fem -> 76 -> L
    feats[308] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. fem', '76') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '76', node2.lemma)
    
    # L -> voc. fem -> -11 -> L
    feats[309] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. fem', '-11') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-11', node2.lemma)
    
    # L -> voc. fem -> sg_tp -> L
    feats[310] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. fem', 'sg_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_tp', node2.lemma)
    
    # L -> voc. fem -> -261 -> L
    feats[311] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. fem', '-261') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-261', node2.lemma)
    
    # L -> voc. fem -> 9_sp -> L
    feats[312] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. fem', '9_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_sp', node2.lemma)
    
    # L -> voc. fem -> 5_sp -> L
    feats[313] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. fem', '5_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sp', node2.lemma)
    
    # L -> voc. fem -> 54 -> L
    feats[314] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. fem', '54') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    
    # L -> voc. fem -> -17 -> L
    feats[315] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. fem', '-17') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-17', node2.lemma)
    
    # L -> voc. fem -> 150 -> L
    feats[316] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. fem', '150') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '150', node2.lemma)
    
    # L -> voc. fem -> 3_sp -> L
    feats[317] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. fem', '3_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_sp', node2.lemma)
    
    # L -> voc. fem -> 74 -> L
    feats[318] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. fem', '74') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '74', node2.lemma)
    
    # L -> 138 -> -114 -> L
    feats[319] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '138') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '138', '-114') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-114', node2.lemma)
    
    # L -> 138 -> du_tp -> L
    feats[320] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '138') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '138', 'du_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_tp', node2.lemma)
    
    # L -> 138 -> -296 -> L
    feats[321] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '138') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '138', '-296') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-296', node2.lemma)
    
    # L -> 138 -> 12_fp -> L
    feats[322] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '138') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '138', '12_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_fp', node2.lemma)
    
    # L -> 138 -> -48 -> L
    feats[323] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '138') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '138', '-48') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-48', node2.lemma)
    
    # L -> 138 -> tp -> L
    feats[324] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '138') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '138', 'tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'tp', node2.lemma)
    
    # L -> 138 -> 132 -> L
    feats[325] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '138') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '138', '132') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '132', node2.lemma)
    
    # L -> 138 -> -297 -> L
    feats[326] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '138') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '138', '-297') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-297', node2.lemma)
    
    # L -> 138 -> -142 -> L
    feats[327] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '138') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '138', '-142') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-142', node2.lemma)
    
    # L -> 138 -> 5_pl -> L
    feats[328] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '138') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '138', '5_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_pl', node2.lemma)
    
    # L -> 138 -> adj -> L
    feats[329] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '138') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '138', 'adj') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'adj', node2.lemma)
    
    # L -> 138 -> 12_tp -> L
    feats[330] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '138') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '138', '12_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_tp', node2.lemma)
    
    # L -> 138 -> -34 -> L
    feats[331] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '138') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '138', '-34') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-34', node2.lemma)
    
    # L -> 138 -> fp -> L
    feats[332] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '138') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '138', 'fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'fp', node2.lemma)
    
    # L -> 138 -> -13 -> L
    feats[333] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '138') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '138', '-13') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-13', node2.lemma)
    
    # L -> 138 -> -247 -> L
    feats[334] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '138') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '138', '-247') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-247', node2.lemma)
    
    # L -> 138 -> 8_fp -> L
    feats[335] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '138') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '138', '8_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    
    # L -> 138 -> 78 -> L
    feats[336] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '138') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '138', '78') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '78', node2.lemma)
    
    # L -> -152 -> 12_sp -> L
    feats[337] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-152') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-152', '12_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_sp', node2.lemma)
    
    # L -> acc. neutr. -> 16_sg -> L
    feats[338] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. neutr.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. neutr.', '16_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_sg', node2.lemma)
    
    # L -> acc. neutr. -> -123 -> L
    feats[339] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. neutr.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. neutr.', '-123') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-123', node2.lemma)
    
    # L -> acc. neutr. -> -11 -> L
    feats[340] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. neutr.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. neutr.', '-11') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-11', node2.lemma)
    
    # L -> acc. neutr. -> -48 -> L
    feats[341] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. neutr.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. neutr.', '-48') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-48', node2.lemma)
    
    # L -> acc. neutr. -> voc. fem -> L
    feats[342] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. neutr.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. neutr.', 'voc. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. fem', node2.lemma)
    
    # L -> acc. neutr. -> 138 -> L
    feats[343] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. neutr.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. neutr.', '138') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '138', node2.lemma)
    
    # L -> acc. neutr. -> -142 -> L
    feats[344] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. neutr.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. neutr.', '-142') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-142', node2.lemma)
    
    # L -> tp -> acc. masc. -> L
    feats[345] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', 'acc. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. masc.', node2.lemma)
    
    # L -> tp -> 169 -> L
    feats[346] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '169') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '169', node2.lemma)
    
    # L -> tp -> 60 -> L
    feats[347] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '60') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '60', node2.lemma)
    
    # L -> tp -> -68 -> L
    feats[348] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', '-68') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-68', node2.lemma)
    
    # L -> 10_fp -> -86 -> L
    feats[349] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_fp', '-86') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-86', node2.lemma)
    
    # L -> 10_fp -> -48 -> L
    feats[350] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_fp', '-48') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-48', node2.lemma)
    
    # L -> 10_fp -> -76 -> L
    feats[351] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_fp', '-76') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-76', node2.lemma)
    
    # L -> 10_fp -> 51 -> L
    feats[352] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_fp', '51') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '51', node2.lemma)
    
    # L -> 10_fp -> acc -> L
    feats[353] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_fp', 'acc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc', node2.lemma)
    
    # L -> 10_fp -> 80 -> L
    feats[354] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_fp', '80') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '80', node2.lemma)
    
    # L -> 10_fp -> 78 -> L
    feats[355] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_fp', '78') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '78', node2.lemma)
    
    # L -> 10_fp -> 153 -> L
    feats[356] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_fp', '153') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    
    # L -> 10_fp -> 8_sp -> L
    feats[357] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_fp', '8_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_sp', node2.lemma)
    
    # L -> 10_fp -> 112 -> L
    feats[358] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_fp', '112') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '112', node2.lemma)
    
    # L -> 10_fp -> -119 -> L
    feats[359] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_fp', '-119') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-119', node2.lemma)
    
    # L -> 10_fp -> voc -> L
    feats[360] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_fp', 'voc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc', node2.lemma)
    
    # L -> 10_fp -> -55 -> L
    feats[361] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_fp', '-55') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-55', node2.lemma)
    
    # L -> 10_fp -> 14_sg -> L
    feats[362] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_fp', '14_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_sg', node2.lemma)
    
    # L -> 10_fp -> acc. du. -> L
    feats[363] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_fp', 'acc. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. du.', node2.lemma)
    
    # L -> 10_fp -> 6_fp -> L
    feats[364] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_fp', '6_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_fp', node2.lemma)
    
    # L -> 10_fp -> 6_du -> L
    feats[365] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_fp', '6_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_du', node2.lemma)
    
    # L -> 10_fp -> -32 -> L
    feats[366] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_fp', '-32') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-32', node2.lemma)
    
    # L -> 10_fp -> masc -> L
    feats[367] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_fp', 'masc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'masc', node2.lemma)
    
    # L -> 10_fp -> 134 -> L
    feats[368] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_fp', '134') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '134', node2.lemma)
    
    # L -> 10_fp -> 16_du -> L
    feats[369] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_fp', '16_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_du', node2.lemma)
    
    # L -> 10_fp -> 4_sg -> L
    feats[370] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_fp', '4_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_sg', node2.lemma)
    
    # L -> 10_fp -> 31 -> L
    feats[371] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_fp', '31') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '31', node2.lemma)
    
    # L -> 10_fp -> 8_sg -> L
    feats[372] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_fp', '8_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_sg', node2.lemma)
    
    # L -> 10_fp -> 168 -> L
    feats[373] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_fp', '168') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '168', node2.lemma)
    
    # L -> 10_fp -> 180 -> L
    feats[374] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_fp', '180') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '180', node2.lemma)
    
    # L -> 10_fp -> 141 -> L
    feats[375] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_fp', '141') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '141', node2.lemma)
    
    # L -> 10_fp -> -56 -> L
    feats[376] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_fp', '-56') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-56', node2.lemma)
    
    # L -> 10_fp -> voc. pl. -> L
    feats[377] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_fp', 'voc. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. pl.', node2.lemma)
    
    # L -> 10_fp -> -37 -> L
    feats[378] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_fp', '-37') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-37', node2.lemma)
    
    # L -> 10_fp -> 30_du -> L
    feats[379] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_fp', '30_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_du', node2.lemma)
    
    # L -> 10_fp -> -90 -> L
    feats[380] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_fp', '-90') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-90', node2.lemma)
    
    # L -> 10_fp -> sg -> L
    feats[381] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_fp', 'sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg', node2.lemma)
    
    # L -> 10_fp -> -16 -> L
    feats[382] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_fp', '-16') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-16', node2.lemma)
    
    # L -> 10_fp -> -84 -> L
    feats[383] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_fp', '-84') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-84', node2.lemma)
    
    # L -> 10_fp -> -72 -> L
    feats[384] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_fp', '-72') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-72', node2.lemma)
    
    # L -> 10_fp -> 15_fp -> L
    feats[385] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_fp', '15_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_fp', node2.lemma)
    
    # L -> 132 -> 5_sp -> L
    feats[386] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '132') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '132', '5_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sp', node2.lemma)
    
    # L -> 132 -> -121 -> L
    feats[387] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '132') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '132', '-121') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-121', node2.lemma)
    
    # L -> 132 -> -23 -> L
    feats[388] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '132') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '132', '-23') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-23', node2.lemma)
    
    # L -> -297 -> 180 -> L
    feats[389] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-297') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-297', '180') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '180', node2.lemma)
    
    # L -> -297 -> -68 -> L
    feats[390] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-297') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-297', '-68') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-68', node2.lemma)
    
    # L -> -142 -> 28_sg -> L
    feats[391] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-142') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-142', '28_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28_sg', node2.lemma)
    
    # L -> -142 -> -122 -> L
    feats[392] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-142') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-142', '-122') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-122', node2.lemma)
    
    # L -> -142 -> 13_sg -> L
    feats[393] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-142') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-142', '13_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_sg', node2.lemma)
    
    # L -> -142 -> du_fp -> L
    feats[394] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-142') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-142', 'du_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_fp', node2.lemma)
    
    # L -> -142 -> -66 -> L
    feats[395] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-142') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-142', '-66') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-66', node2.lemma)
    
    # L -> -142 -> -83 -> L
    feats[396] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-142') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-142', '-83') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-83', node2.lemma)
    
    # L -> -142 -> 27_pl -> L
    feats[397] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-142') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-142', '27_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_pl', node2.lemma)
    
    # L -> -142 -> 6_sg -> L
    feats[398] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-142') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-142', '6_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_sg', node2.lemma)
    
    # L -> -291 -> -53 -> L
    feats[399] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-291') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-291', '-53') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-53', node2.lemma)
    
    # L -> -291 -> -99 -> L
    feats[400] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-291') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-291', '-99') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-99', node2.lemma)
    
    # L -> -291 -> 6_sp -> L
    feats[401] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-291') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-291', '6_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_sp', node2.lemma)
    
    # L -> -291 -> 7_fp -> L
    feats[402] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-291') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-291', '7_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_fp', node2.lemma)
    
    # L -> -291 -> 152 -> L
    feats[403] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-291') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-291', '152') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '152', node2.lemma)
    
    # L -> -291 -> -96 -> L
    feats[404] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-291') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-291', '-96') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-96', node2.lemma)
    
    # L -> -291 -> 58 -> L
    feats[405] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-291') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-291', '58') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '58', node2.lemma)
    
    # L -> -291 -> 5_sg -> L
    feats[406] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-291') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-291', '5_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sg', node2.lemma)
    
    # L -> -291 -> 14_pl -> L
    feats[407] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-291') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-291', '14_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_pl', node2.lemma)
    
    # L -> -291 -> acc. du. -> L
    feats[408] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-291') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-291', 'acc. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. du.', node2.lemma)
    
    # L -> -291 -> 16_tp -> L
    feats[409] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-291') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-291', '16_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_tp', node2.lemma)
    
    # L -> -291 -> 158 -> L
    feats[410] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-291') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-291', '158') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '158', node2.lemma)
    
    # L -> -149 -> 27_fp -> L
    feats[411] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-149', '27_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_fp', node2.lemma)
    
    # L -> -149 -> sg_tp -> L
    feats[412] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-149', 'sg_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_tp', node2.lemma)
    
    # L -> -149 -> -72 -> L
    feats[413] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-149', '-72') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-72', node2.lemma)
    
    # L -> -149 -> -141 -> L
    feats[414] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-149', '-141') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-141', node2.lemma)
    
    # L -> -149 -> -129 -> L
    feats[415] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-149', '-129') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-129', node2.lemma)
    
    # L -> -149 -> 30_sp -> L
    feats[416] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-149', '30_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_sp', node2.lemma)
    
    # L -> -149 -> -303 -> L
    feats[417] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-149', '-303') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-303', node2.lemma)
    
    # L -> -271 -> -147 -> L
    feats[418] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-271') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-271', '-147') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-147', node2.lemma)
    
    # L -> -271 -> voc. du. -> L
    feats[419] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-271') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-271', 'voc. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. du.', node2.lemma)
    
    # L -> -271 -> 2_sg -> L
    feats[420] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-271') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-271', '2_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_sg', node2.lemma)
    
    # L -> -271 -> 11_sp -> L
    feats[421] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-271') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-271', '11_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_sp', node2.lemma)
    
    # L -> -271 -> -122 -> L
    feats[422] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-271') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-271', '-122') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-122', node2.lemma)
    
    # L -> -271 -> -97 -> L
    feats[423] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-271') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-271', '-97') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-97', node2.lemma)
    
    # L -> -271 -> 110 -> L
    feats[424] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-271') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-271', '110') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '110', node2.lemma)
    
    # L -> -271 -> 72 -> L
    feats[425] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-271') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-271', '72') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '72', node2.lemma)
    
    # L -> -271 -> 12_fp -> L
    feats[426] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-271') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-271', '12_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_fp', node2.lemma)
    
    # L -> -271 -> -169 -> L
    feats[427] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-271') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-271', '-169') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-169', node2.lemma)
    
    # L -> -271 -> sg_tp -> L
    feats[428] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-271') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-271', 'sg_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_tp', node2.lemma)
    
    # L -> -271 -> 54 -> L
    feats[429] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-271') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-271', '54') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    
    # L -> -271 -> -161 -> L
    feats[430] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-271') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-271', '-161') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-161', node2.lemma)
    
    # L -> -271 -> 179 -> L
    feats[431] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-271') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-271', '179') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '179', node2.lemma)
    
    # L -> -271 -> tp -> L
    feats[432] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-271') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-271', 'tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'tp', node2.lemma)
    
    # L -> -271 -> 10_fp -> L
    feats[433] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-271') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-271', '10_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_fp', node2.lemma)
    
    # L -> -271 -> -142 -> L
    feats[434] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-271') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-271', '-142') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-142', node2.lemma)
    
    # L -> -271 -> -77 -> L
    feats[435] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-271') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-271', '-77') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-77', node2.lemma)
    
    # L -> -271 -> instr. du. -> L
    feats[436] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-271') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-271', 'instr. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. du.', node2.lemma)
    
    # L -> -271 -> instr. pl. -> L
    feats[437] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-271') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-271', 'instr. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. pl.', node2.lemma)
    
    # L -> -271 -> -166 -> L
    feats[438] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-271') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-271', '-166') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-166', node2.lemma)
    
    # L -> -271 -> -52 -> L
    feats[439] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-271') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-271', '-52') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-52', node2.lemma)
    
    # L -> -271 -> 75 -> L
    feats[440] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-271') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-271', '75') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '75', node2.lemma)
    
    # L -> -271 -> sg_fp -> L
    feats[441] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-271') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-271', 'sg_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_fp', node2.lemma)
    
    # L -> -271 -> acc. adj. -> L
    feats[442] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-271') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-271', 'acc. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. adj.', node2.lemma)
    
    # L -> -271 -> 49 -> L
    feats[443] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-271') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-271', '49') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '49', node2.lemma)
    
    # L -> -271 -> 3_fp -> L
    feats[444] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-271') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-271', '3_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_fp', node2.lemma)
    
    # L -> -271 -> 152 -> L
    feats[445] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-271') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-271', '152') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '152', node2.lemma)
    
    # L -> -271 -> -283 -> L
    feats[446] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-271') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-271', '-283') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-283', node2.lemma)
    
    # L -> -271 -> -133 -> L
    feats[447] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-271') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-271', '-133') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-133', node2.lemma)
    
    # L -> -271 -> 181 -> L
    feats[448] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-271') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-271', '181') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '181', node2.lemma)
    
    # L -> -271 -> -299 -> L
    feats[449] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-271') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-271', '-299') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-299', node2.lemma)
    
    # L -> -271 -> 108 -> L
    feats[450] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-271') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-271', '108') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '108', node2.lemma)
    
    # L -> -271 -> 74 -> L
    feats[451] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-271') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-271', '74') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '74', node2.lemma)
    
    # L -> 10_sg -> 16_sg -> L
    feats[452] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sg', '16_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_sg', node2.lemma)
    
    # L -> 10_sg -> abl. pl. -> L
    feats[453] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sg', 'abl. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. pl.', node2.lemma)
    
    # L -> 10_sg -> loc. sg. -> L
    feats[454] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sg', 'loc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc. sg.', node2.lemma)
    
    # L -> 10_sg -> instr. pl. -> L
    feats[455] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sg', 'instr. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. pl.', node2.lemma)
    
    # L -> 10_sg -> -132 -> L
    feats[456] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sg', '-132') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-132', node2.lemma)
    
    # L -> 10_sg -> acc. adj. -> L
    feats[457] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sg', 'acc. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. adj.', node2.lemma)
    
    # L -> 10_sg -> -73 -> L
    feats[458] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sg', '-73') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-73', node2.lemma)
    
    # L -> 10_sg -> -279 -> L
    feats[459] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sg', '-279') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-279', node2.lemma)
    
    # L -> 10_sg -> -113 -> L
    feats[460] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sg', '-113') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-113', node2.lemma)
    
    # L -> 10_sg -> 4_sp -> L
    feats[461] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sg', '4_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_sp', node2.lemma)
    
    # L -> 10_sg -> 9_fp -> L
    feats[462] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sg', '9_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_fp', node2.lemma)
    
    # L -> 27_sg -> -112 -> L
    feats[463] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', '-112') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    
    # L -> 27_sg -> acc. fem -> L
    feats[464] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', 'acc. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. fem', node2.lemma)
    
    # L -> 27_sg -> 30_sg -> L
    feats[465] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', '30_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_sg', node2.lemma)
    
    # L -> 27_sg -> -302 -> L
    feats[466] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', '-302') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-302', node2.lemma)
    
    # L -> 27_sg -> -104 -> L
    feats[467] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', '-104') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-104', node2.lemma)
    
    # L -> 27_sg -> -73 -> L
    feats[468] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', '-73') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-73', node2.lemma)
    
    # L -> 27_sg -> 114 -> L
    feats[469] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', '114') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '114', node2.lemma)
    
    # L -> 27_sg -> 70 -> L
    feats[470] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', '70') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '70', node2.lemma)
    
    # L -> 27_sg -> -27 -> L
    feats[471] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', '-27') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-27', node2.lemma)
    
    # L -> 27_sg -> -92 -> L
    feats[472] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', '-92') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-92', node2.lemma)
    
    # L -> 27_sg -> 30 -> L
    feats[473] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', '30') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30', node2.lemma)
    
    # L -> 27_sg -> -113 -> L
    feats[474] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', '-113') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-113', node2.lemma)
    
    # L -> 27_sg -> 117 -> L
    feats[475] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', '117') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '117', node2.lemma)
    
    # L -> 27_sg -> -158 -> L
    feats[476] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', '-158') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-158', node2.lemma)
    
    # L -> 27_sg -> 135 -> L
    feats[477] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', '135') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '135', node2.lemma)
    
    # L -> 27_sg -> -307 -> L
    feats[478] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', '-307') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-307', node2.lemma)
    
    # L -> 27_sg -> 95 -> L
    feats[479] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', '95') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '95', node2.lemma)
    
    # L -> 27_sg -> 12_sg -> L
    feats[480] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', '12_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_sg', node2.lemma)
    
    # L -> 27_sg -> -50 -> L
    feats[481] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', '-50') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-50', node2.lemma)
    
    # L -> 27_sg -> -101 -> L
    feats[482] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', '-101') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-101', node2.lemma)
    
    # L -> 27_sg -> -293 -> L
    feats[483] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', '-293') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-293', node2.lemma)
    
    # L -> 27_sg -> -63 -> L
    feats[484] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', '-63') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-63', node2.lemma)
    
    # L -> 27_sg -> 42 -> L
    feats[485] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', '42') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '42', node2.lemma)
    
    # L -> 27_sg -> gen. pl. -> L
    feats[486] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', 'gen. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen. pl.', node2.lemma)
    
    # L -> 27_sg -> 11_tp -> L
    feats[487] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', '11_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_tp', node2.lemma)
    
    # L -> 27_sg -> neutr -> L
    feats[488] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', 'neutr') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'neutr', node2.lemma)
    
    # L -> 27_sg -> -243 -> L
    feats[489] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', '-243') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-243', node2.lemma)
    
    # L -> 27_sg -> -276 -> L
    feats[490] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', '-276') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-276', node2.lemma)
    
    # L -> 27_sg -> 11_fp -> L
    feats[491] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', '11_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_fp', node2.lemma)
    
    # L -> 27_sg -> -301 -> L
    feats[492] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', '-301') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-301', node2.lemma)
    
    # L -> 27_sg -> 73 -> L
    feats[493] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', '73') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '73', node2.lemma)
    
    # L -> 27_sg -> instr. masc. -> L
    feats[494] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', 'instr. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. masc.', node2.lemma)
    
    # L -> 27_sg -> gen. du. -> L
    feats[495] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', 'gen. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen. du.', node2.lemma)
    
    # L -> 27_sg -> -283 -> L
    feats[496] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', '-283') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-283', node2.lemma)
    
    # L -> 27_sg -> 11_du -> L
    feats[497] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', '11_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_du', node2.lemma)
    
    # L -> 27_sg -> 8_fp -> L
    feats[498] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', '8_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    
    # L -> 27_sg -> 129 -> L
    feats[499] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', '129') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '129', node2.lemma)
    
    # L -> 27_sg -> -96 -> L
    feats[500] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', '-96') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-96', node2.lemma)
    
    # L -> 27_sg -> 58 -> L
    feats[501] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', '58') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '58', node2.lemma)
    
    # L -> 27_sg -> 153 -> L
    feats[502] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', '153') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
    
    # L -> 27_sg -> 178 -> L
    feats[503] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', '178') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '178', node2.lemma)
    
    # L -> 27_sg -> 14_sg -> L
    feats[504] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', '14_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_sg', node2.lemma)
    
    # L -> 27_sg -> -29 -> L
    feats[505] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', '-29') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-29', node2.lemma)
    
    # L -> 27_sg -> -18 -> L
    feats[506] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', '-18') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-18', node2.lemma)
    
    # L -> 27_sg -> pl_sp -> L
    feats[507] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', 'pl_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_sp', node2.lemma)
    
    # L -> 27_sg -> -35 -> L
    feats[508] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', '-35') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    
    # L -> 27_sg -> -308 -> L
    feats[509] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', '-308') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-308', node2.lemma)
    
    # L -> 27_sg -> 131 -> L
    feats[510] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', '131') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '131', node2.lemma)
    
    # L -> 27_sg -> -141 -> L
    feats[511] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', '-141') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-141', node2.lemma)
    
    # L -> 27_sg -> 2_sp -> L
    feats[512] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', '2_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_sp', node2.lemma)
    
    # L -> 27_sg -> 15_fp -> L
    feats[513] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_sg', '15_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_fp', node2.lemma)
    
    # L -> -22 -> 3_tp -> L
    feats[514] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-22') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-22', '3_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_tp', node2.lemma)
    
    # L -> -22 -> 39 -> L
    feats[515] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-22') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-22', '39') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '39', node2.lemma)
    
    # L -> -22 -> 2_sg -> L
    feats[516] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-22') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-22', '2_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_sg', node2.lemma)
    
    # L -> -22 -> 13_sg -> L
    feats[517] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-22') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-22', '13_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_sg', node2.lemma)
    
    # L -> -22 -> -97 -> L
    feats[518] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-22') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-22', '-97') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-97', node2.lemma)
    
    # L -> 10_tp -> -12 -> L
    feats[519] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_tp', '-12') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-12', node2.lemma)
    
    # L -> 10_tp -> -90 -> L
    feats[520] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_tp', '-90') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-90', node2.lemma)
    
    # L -> -77 -> voc. neutr. -> L
    feats[521] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-77') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-77', 'voc. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. neutr.', node2.lemma)
    
    # L -> -77 -> 8_du -> L
    feats[522] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-77') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-77', '8_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_du', node2.lemma)
    
    # L -> -77 -> -45 -> L
    feats[523] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-77') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-77', '-45') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-45', node2.lemma)
    
    # L -> -77 -> -69 -> L
    feats[524] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-77') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-77', '-69') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-69', node2.lemma)
    
    # L -> -77 -> 101 -> L
    feats[525] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-77') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-77', '101') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    
    # L -> -77 -> -240 -> L
    feats[526] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-77') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-77', '-240') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-240', node2.lemma)
    
    # L -> -77 -> 160 -> L
    feats[527] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-77') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-77', '160') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    
    # L -> -77 -> 94 -> L
    feats[528] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-77') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-77', '94') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '94', node2.lemma)
    
    # L -> -77 -> 74 -> L
    feats[529] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-77') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-77', '74') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '74', node2.lemma)
    
    # L -> -77 -> 27_du -> L
    feats[530] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-77') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-77', '27_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_du', node2.lemma)
    
    # L -> -77 -> 2_sp -> L
    feats[531] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-77') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-77', '2_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_sp', node2.lemma)
    
    # L -> -77 -> 15_fp -> L
    feats[532] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-77') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-77', '15_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_fp', node2.lemma)
    
    # L -> -77 -> nom -> L
    feats[533] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-77') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-77', 'nom') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom', node2.lemma)
    
    # L -> -77 -> 13_fp -> L
    feats[534] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-77') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-77', '13_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_fp', node2.lemma)
    
    # L -> -77 -> 176 -> L
    feats[535] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-77') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-77', '176') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '176', node2.lemma)
    
    # L -> -77 -> nom. masc. -> L
    feats[536] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-77') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-77', 'nom. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. masc.', node2.lemma)
    
    # L -> -77 -> 8_pl -> L
    feats[537] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-77') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-77', '8_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_pl', node2.lemma)
    
    # L -> -77 -> -28 -> L
    feats[538] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-77') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-77', '-28') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-28', node2.lemma)
    
    # L -> -77 -> -303 -> L
    feats[539] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-77') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-77', '-303') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-303', node2.lemma)
    
    # L -> -245 -> abl. du. -> L
    feats[540] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-245') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-245', 'abl. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. du.', node2.lemma)
    
    # L -> -245 -> -147 -> L
    feats[541] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-245') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-245', '-147') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-147', node2.lemma)
    
    # L -> -245 -> -159 -> L
    feats[542] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-245') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-245', '-159') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-159', node2.lemma)
    
    # L -> -245 -> 2_sg -> L
    feats[543] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-245') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-245', '2_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_sg', node2.lemma)
    
    # L -> -245 -> 88 -> L
    feats[544] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-245') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-245', '88') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '88', node2.lemma)
    
    # L -> -245 -> 11_sp -> L
    feats[545] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-245') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-245', '11_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_sp', node2.lemma)
    
    # L -> -245 -> du -> L
    feats[546] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-245') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-245', 'du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du', node2.lemma)
    
    # L -> -245 -> -97 -> L
    feats[547] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-245') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-245', '-97') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-97', node2.lemma)
    
    # L -> -245 -> -10 -> L
    feats[548] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-245') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-245', '-10') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-10', node2.lemma)
    
    # L -> -245 -> 16_pl -> L
    feats[549] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-245') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-245', '16_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_pl', node2.lemma)
    
    # L -> -245 -> 12_fp -> L
    feats[550] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-245') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-245', '12_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_fp', node2.lemma)
    
    # L -> -245 -> 76 -> L
    feats[551] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-245') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-245', '76') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '76', node2.lemma)
    
    # L -> -245 -> -169 -> L
    feats[552] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-245') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-245', '-169') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-169', node2.lemma)
    
    # L -> -245 -> 128 -> L
    feats[553] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-245') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-245', '128') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '128', node2.lemma)
    
    # L -> -245 -> -123 -> L
    feats[554] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-245') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-245', '-123') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-123', node2.lemma)
    
    # L -> -245 -> -71 -> L
    feats[555] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-245') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-245', '-71') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-71', node2.lemma)
    
    # L -> -245 -> 56 -> L
    feats[556] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-245') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-245', '56') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '56', node2.lemma)
    
    # L -> -245 -> 54 -> L
    feats[557] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-245') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-245', '54') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    
    # L -> -245 -> 150 -> L
    feats[558] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-245') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-245', '150') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '150', node2.lemma)
    
    # L -> -245 -> -161 -> L
    feats[559] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-245') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-245', '-161') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-161', node2.lemma)
    
    # L -> -245 -> pl_fp -> L
    feats[560] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-245') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-245', 'pl_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_fp', node2.lemma)
    
    # L -> -245 -> voc. fem -> L
    feats[561] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-245') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-245', 'voc. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. fem', node2.lemma)
    
    # L -> -245 -> acc. neutr. -> L
    feats[562] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-245') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-245', 'acc. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. neutr.', node2.lemma)
    
    # L -> -245 -> tp -> L
    feats[563] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-245') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-245', 'tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'tp', node2.lemma)
    
    # L -> -245 -> 10_sg -> L
    feats[564] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-245') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-245', '10_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_sg', node2.lemma)
    
    # L -> -245 -> 27_sg -> L
    feats[565] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-245') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-245', '27_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_sg', node2.lemma)
    
    # L -> -245 -> -22 -> L
    feats[566] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-245') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-245', '-22') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-22', node2.lemma)
    
    # L -> -245 -> 10_sp -> L
    feats[567] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-245') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-245', '10_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_sp', node2.lemma)
    
    # L -> -245 -> -166 -> L
    feats[568] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-245') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-245', '-166') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-166', node2.lemma)
    
    # L -> -245 -> 151 -> L
    feats[569] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-245') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-245', '151') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '151', node2.lemma)
    
    # L -> -245 -> 98 -> L
    feats[570] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-245') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-245', '98') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '98', node2.lemma)
    
    # L -> -245 -> 12_tp -> L
    feats[571] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-245') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-245', '12_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_tp', node2.lemma)
    
    # L -> -245 -> 7_tp -> L
    feats[572] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-245') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-245', '7_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_tp', node2.lemma)
    
    # L -> -245 -> -230 -> L
    feats[573] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-245') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-245', '-230') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-230', node2.lemma)
    
    # L -> -245 -> 9_tp -> L
    feats[574] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-245') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-245', '9_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_tp', node2.lemma)
    
    # L -> -245 -> 118 -> L
    feats[575] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-245') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-245', '118') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '118', node2.lemma)
    
    # L -> -245 -> -111 -> L
    feats[576] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-245') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-245', '-111') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-111', node2.lemma)
    
    # L -> -245 -> 8_tp -> L
    feats[577] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-245') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-245', '8_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_tp', node2.lemma)
    
    # L -> -245 -> 139 -> L
    feats[578] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-245') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-245', '139') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '139', node2.lemma)
    
    # L -> -245 -> 36 -> L
    feats[579] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-245') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-245', '36') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '36', node2.lemma)
    
    # L -> -245 -> -269 -> L
    feats[580] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-245') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-245', '-269') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-269', node2.lemma)
    
    # L -> -245 -> 16_fp -> L
    feats[581] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-245') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-245', '16_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_fp', node2.lemma)
    
    # L -> -245 -> 134 -> L
    feats[582] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-245') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-245', '134') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '134', node2.lemma)
    
    # L -> -245 -> 158 -> L
    feats[583] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-245') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-245', '158') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '158', node2.lemma)
    
    # L -> -245 -> -84 -> L
    feats[584] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-245') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-245', '-84') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-84', node2.lemma)
    
    # L -> instr. du. -> -163 -> L
    feats[585] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. du.', '-163') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-163', node2.lemma)
    
    # L -> instr. du. -> 54 -> L
    feats[586] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. du.', '54') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    
    # L -> instr. du. -> 10_fp -> L
    feats[587] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. du.', '10_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_fp', node2.lemma)
    
    # L -> instr. du. -> 90 -> L
    feats[588] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. du.', '90') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '90', node2.lemma)
    
    # L -> instr. du. -> -302 -> L
    feats[589] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. du.', '-302') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-302', node2.lemma)
    
    # L -> instr. du. -> -98 -> L
    feats[590] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. du.', '-98') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-98', node2.lemma)
    
    # L -> instr. du. -> -47 -> L
    feats[591] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. du.', '-47') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-47', node2.lemma)
    
    # L -> instr. du. -> -299 -> L
    feats[592] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. du.', '-299') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-299', node2.lemma)
    
    # L -> instr. du. -> 2 -> L
    feats[593] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. du.', '2') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2', node2.lemma)
    
    # L -> instr. du. -> 61 -> L
    feats[594] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. du.', '61') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '61', node2.lemma)
    
    # L -> instr. du. -> 5_du -> L
    feats[595] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. du.', '5_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_du', node2.lemma)
    
    # L -> instr. du. -> 155 -> L
    feats[596] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. du.', '155') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '155', node2.lemma)
    
    # L -> instr. du. -> 82 -> L
    feats[597] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. du.', '82') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '82', node2.lemma)
    
    # L -> instr. du. -> -57 -> L
    feats[598] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. du.', '-57') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-57', node2.lemma)
    
    # L -> instr. du. -> 180 -> L
    feats[599] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. du.', '180') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '180', node2.lemma)
    
    # L -> instr. du. -> -56 -> L
    feats[600] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. du.', '-56') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-56', node2.lemma)
    
    # L -> instr. du. -> 9_pl -> L
    feats[601] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. du.', '9_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_pl', node2.lemma)
    
    # L -> instr. du. -> 91 -> L
    feats[602] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. du.', '91') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    
    # L -> instr. du. -> -45 -> L
    feats[603] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. du.', '-45') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-45', node2.lemma)
    
    # L -> instr. du. -> 131 -> L
    feats[604] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. du.', '131') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '131', node2.lemma)
    
    # L -> instr. du. -> 101 -> L
    feats[605] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. du.', '101') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    
    # L -> instr. du. -> 158 -> L
    feats[606] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. du.', '158') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '158', node2.lemma)
    
    # L -> instr. du. -> -91 -> L
    feats[607] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. du.', '-91') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-91', node2.lemma)
    
    # L -> instr. du. -> -90 -> L
    feats[608] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. du.', '-90') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-90', node2.lemma)
    
    # L -> instr. du. -> -30 -> L
    feats[609] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. du.', '-30') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-30', node2.lemma)
    
    # L -> instr. du. -> 74 -> L
    feats[610] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. du.', '74') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '74', node2.lemma)
    
    # L -> instr. du. -> -78 -> L
    feats[611] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. du.', '-78') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-78', node2.lemma)
    
    # L -> instr. du. -> 27_du -> L
    feats[612] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. du.', '27_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_du', node2.lemma)
    
    # L -> instr. du. -> 12_sp -> L
    feats[613] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. du.', '12_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_sp', node2.lemma)
    
    # L -> instr. du. -> -72 -> L
    feats[614] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. du.', '-72') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-72', node2.lemma)
    
    # L -> instr. du. -> nom. masc. -> L
    feats[615] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. du.', 'nom. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. masc.', node2.lemma)
    
    # L -> instr. du. -> 8_pl -> L
    feats[616] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. du.', '8_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_pl', node2.lemma)
    
    # L -> instr. du. -> -28 -> L
    feats[617] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. du.', '-28') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-28', node2.lemma)
    
    # L -> instr. pl. -> abl. du. -> L
    feats[618] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. pl.', 'abl. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. du.', node2.lemma)
    
    # L -> instr. pl. -> -97 -> L
    feats[619] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. pl.', '-97') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-97', node2.lemma)
    
    # L -> instr. pl. -> 159 -> L
    feats[620] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. pl.', '159') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '159', node2.lemma)
    
    # L -> instr. pl. -> 9_sp -> L
    feats[621] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. pl.', '9_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_sp', node2.lemma)
    
    # L -> instr. pl. -> -161 -> L
    feats[622] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. pl.', '-161') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-161', node2.lemma)
    
    # L -> instr. pl. -> 79 -> L
    feats[623] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. pl.', '79') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '79', node2.lemma)
    
    # L -> instr. pl. -> 6_pl -> L
    feats[624] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. pl.', '6_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_pl', node2.lemma)
    
    # L -> instr. pl. -> voc. fem -> L
    feats[625] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. pl.', 'voc. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. fem', node2.lemma)
    
    # L -> instr. pl. -> -297 -> L
    feats[626] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. pl.', '-297') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-297', node2.lemma)
    
    # L -> instr. pl. -> 10_sg -> L
    feats[627] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. pl.', '10_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_sg', node2.lemma)
    
    # L -> instr. pl. -> -166 -> L
    feats[628] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. pl.', '-166') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-166', node2.lemma)
    
    # L -> instr. pl. -> -82 -> L
    feats[629] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. pl.', '-82') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-82', node2.lemma)
    
    # L -> instr. pl. -> -112 -> L
    feats[630] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. pl.', '-112') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    
    # L -> instr. pl. -> -266 -> L
    feats[631] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. pl.', '-266') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-266', node2.lemma)
    
    # L -> instr. pl. -> loc -> L
    feats[632] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. pl.', 'loc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc', node2.lemma)
    
    # L -> instr. pl. -> 149 -> L
    feats[633] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. pl.', '149') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '149', node2.lemma)
    
    # L -> instr. pl. -> -263 -> L
    feats[634] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. pl.', '-263') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-263', node2.lemma)
    
    # L -> instr. pl. -> 14_tp -> L
    feats[635] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. pl.', '14_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_tp', node2.lemma)
    
    # L -> instr. pl. -> -292 -> L
    feats[636] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. pl.', '-292') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-292', node2.lemma)
    
    # L -> instr. pl. -> -156 -> L
    feats[637] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. pl.', '-156') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-156', node2.lemma)
    
    # L -> instr. pl. -> -42 -> L
    feats[638] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. pl.', '-42') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-42', node2.lemma)
    
    # L -> instr. pl. -> 97 -> L
    feats[639] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. pl.', '97') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '97', node2.lemma)
    
    # L -> instr. pl. -> -26 -> L
    feats[640] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. pl.', '-26') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-26', node2.lemma)
    
    # L -> instr. pl. -> acc. du. -> L
    feats[641] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. pl.', 'acc. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. du.', node2.lemma)
    
    # L -> instr. pl. -> 15_pl -> L
    feats[642] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. pl.', '15_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_pl', node2.lemma)
    
    # L -> instr. pl. -> -83 -> L
    feats[643] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. pl.', '-83') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-83', node2.lemma)
    
    # L -> instr. pl. -> 134 -> L
    feats[644] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. pl.', '134') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '134', node2.lemma)
    
    # L -> instr. pl. -> -81 -> L
    feats[645] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. pl.', '-81') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-81', node2.lemma)
    
    # L -> instr. pl. -> nom. pl. -> L
    feats[646] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. pl.', 'nom. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. pl.', node2.lemma)
    
    # L -> dat -> -261 -> L
    feats[647] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat', '-261') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-261', node2.lemma)
    
    # L -> dat -> voc. fem -> L
    feats[648] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat', 'voc. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. fem', node2.lemma)
    
    # L -> dat -> 75 -> L
    feats[649] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat', '75') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '75', node2.lemma)
    
    # L -> dat -> -302 -> L
    feats[650] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat', '-302') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-302', node2.lemma)
    
    # L -> -25 -> -44 -> L
    feats[651] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-25') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-25', '-44') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-44', node2.lemma)
    
    # L -> -25 -> -115 -> L
    feats[652] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-25') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-25', '-115') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-115', node2.lemma)
    
    # L -> -25 -> -263 -> L
    feats[653] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-25') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-25', '-263') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-263', node2.lemma)
    
    # L -> -25 -> 114 -> L
    feats[654] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-25') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-25', '114') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '114', node2.lemma)
    
    # L -> -25 -> -49 -> L
    feats[655] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-25') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-25', '-49') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-49', node2.lemma)
    
    # L -> -25 -> -27 -> L
    feats[656] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-25') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-25', '-27') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-27', node2.lemma)
    
    # L -> -25 -> -111 -> L
    feats[657] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-25') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-25', '-111') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-111', node2.lemma)
    
    # L -> -25 -> 30 -> L
    feats[658] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-25') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-25', '30') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30', node2.lemma)
    
    # L -> -25 -> 117 -> L
    feats[659] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-25') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-25', '117') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '117', node2.lemma)
    
    # L -> -25 -> -143 -> L
    feats[660] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-25') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-25', '-143') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-143', node2.lemma)
    
    # L -> -25 -> gen -> L
    feats[661] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-25') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-25', 'gen') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen', node2.lemma)
    
    # L -> -25 -> -307 -> L
    feats[662] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-25') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-25', '-307') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-307', node2.lemma)
    
    # L -> -25 -> 81 -> L
    feats[663] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-25') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-25', '81') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '81', node2.lemma)
    
    # L -> -25 -> -24 -> L
    feats[664] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-25') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-25', '-24') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    
    # L -> -25 -> -51 -> L
    feats[665] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-25') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-25', '-51') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-51', node2.lemma)
    
    # L -> -25 -> -293 -> L
    feats[666] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-25') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-25', '-293') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-293', node2.lemma)
    
    # L -> -25 -> -63 -> L
    feats[667] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-25') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-25', '-63') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-63', node2.lemma)
    
    # L -> -25 -> acc. sg. -> L
    feats[668] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-25') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-25', 'acc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. sg.', node2.lemma)
    
    # L -> -25 -> -89 -> L
    feats[669] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-25') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-25', '-89') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-89', node2.lemma)
    
    # L -> -25 -> sg_sp -> L
    feats[670] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-25') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-25', 'sg_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_sp', node2.lemma)
    
    # L -> -25 -> 15_sp -> L
    feats[671] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-25') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-25', '15_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_sp', node2.lemma)
    
    # L -> -25 -> dat. pl. -> L
    feats[672] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-25') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-25', 'dat. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. pl.', node2.lemma)
    
    # L -> -25 -> 89 -> L
    feats[673] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-25') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-25', '89') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '89', node2.lemma)
    
    # L -> -25 -> 7_fp -> L
    feats[674] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-25') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-25', '7_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_fp', node2.lemma)
    
    # L -> -139 -> -29 -> L
    feats[675] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-29') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-29', node2.lemma)
    
    # L -> -139 -> -32 -> L
    feats[676] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-32') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-32', node2.lemma)
    
    # L -> -139 -> -12 -> L
    feats[677] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-12') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-12', node2.lemma)
    
    # L -> -139 -> -41 -> L
    feats[678] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-41') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-41', node2.lemma)
    
    # L -> -139 -> 157 -> L
    feats[679] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '157') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '157', node2.lemma)
    
    # L -> -139 -> 29_tp -> L
    feats[680] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '29_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '29_tp', node2.lemma)
    
    # L -> -139 -> -299 -> L
    feats[681] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-299') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-299', node2.lemma)
    
    # L -> -139 -> -18 -> L
    feats[682] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-18') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-18', node2.lemma)
    
    # L -> -139 -> 8_sg -> L
    feats[683] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '8_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_sg', node2.lemma)
    
    # L -> -139 -> -57 -> L
    feats[684] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-57') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-57', node2.lemma)
    
    # L -> -139 -> 180 -> L
    feats[685] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '180') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '180', node2.lemma)
    
    # L -> -139 -> voc. pl. -> L
    feats[686] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', 'voc. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. pl.', node2.lemma)
    
    # L -> -139 -> -35 -> L
    feats[687] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-35') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    
    # L -> -139 -> -117 -> L
    feats[688] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-117') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-117', node2.lemma)
    
    # L -> -139 -> 9_pl -> L
    feats[689] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '9_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_pl', node2.lemma)
    
    # L -> -139 -> 91 -> L
    feats[690] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '91') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    
    # L -> -139 -> 4_tp -> L
    feats[691] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '4_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_tp', node2.lemma)
    
    # L -> -139 -> 6_sg -> L
    feats[692] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '6_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_sg', node2.lemma)
    
    # L -> -139 -> 29 -> L
    feats[693] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '29') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '29', node2.lemma)
    
    # L -> -139 -> -308 -> L
    feats[694] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-308') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-308', node2.lemma)
    
    # L -> -139 -> -37 -> L
    feats[695] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-37') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-37', node2.lemma)
    
    # L -> -139 -> 130 -> L
    feats[696] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '130') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '130', node2.lemma)
    
    # L -> -139 -> -69 -> L
    feats[697] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-69') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-69', node2.lemma)
    
    # L -> -139 -> 158 -> L
    feats[698] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '158') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '158', node2.lemma)
    
    # L -> -139 -> 160 -> L
    feats[699] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '160') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    
    # L -> -139 -> -91 -> L
    feats[700] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-91') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-91', node2.lemma)
    
    # L -> -139 -> 3_pl -> L
    feats[701] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '3_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_pl', node2.lemma)
    
    # L -> -139 -> 30_du -> L
    feats[702] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '30_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_du', node2.lemma)
    
    # L -> -139 -> 10_pl -> L
    feats[703] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '10_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_pl', node2.lemma)
    
    # L -> -139 -> 3_sp -> L
    feats[704] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '3_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_sp', node2.lemma)
    
    # L -> -139 -> 55 -> L
    feats[705] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '55') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '55', node2.lemma)
    
    # L -> -139 -> 27_du -> L
    feats[706] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '27_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_du', node2.lemma)
    
    # L -> -139 -> -72 -> L
    feats[707] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-72') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-72', node2.lemma)
    
    # L -> -139 -> -141 -> L
    feats[708] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-141') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-141', node2.lemma)
    
    # L -> -139 -> -64 -> L
    feats[709] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-64') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-64', node2.lemma)
    
    # L -> -139 -> 176 -> L
    feats[710] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '176') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '176', node2.lemma)
    
    # L -> -139 -> 30_sp -> L
    feats[711] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '30_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_sp', node2.lemma)
    
    # L -> -139 -> nom. masc. -> L
    feats[712] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', 'nom. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. masc.', node2.lemma)
    
    # L -> -139 -> 1 -> L
    feats[713] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '1') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '1', node2.lemma)
    
    # L -> -139 -> -303 -> L
    feats[714] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', '-303') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-303', node2.lemma)
    
    # L -> 151 -> abl. du. -> L
    feats[715] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '151') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '151', 'abl. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. du.', node2.lemma)
    
    # L -> 151 -> 28_sg -> L
    feats[716] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '151') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '151', '28_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28_sg', node2.lemma)
    
    # L -> 151 -> 11_sg -> L
    feats[717] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '151') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '151', '11_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_sg', node2.lemma)
    
    # L -> 151 -> -147 -> L
    feats[718] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '151') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '151', '-147') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-147', node2.lemma)
    
    # L -> 151 -> -15 -> L
    feats[719] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '151') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '151', '-15') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-15', node2.lemma)
    
    # L -> 151 -> 4_pl -> L
    feats[720] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '151') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '151', '4_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_pl', node2.lemma)
    
    # L -> 151 -> 39 -> L
    feats[721] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '151') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '151', '39') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '39', node2.lemma)
    
    # L -> 151 -> du_tp -> L
    feats[722] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '151') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '151', 'du_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_tp', node2.lemma)
    
    # L -> 151 -> voc. masc. -> L
    feats[723] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '151') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '151', 'voc. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. masc.', node2.lemma)
    
    # L -> 151 -> -163 -> L
    feats[724] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '151') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '151', '-163') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-163', node2.lemma)
    
    # L -> 151 -> -296 -> L
    feats[725] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '151') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '151', '-296') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-296', node2.lemma)
    
    # L -> 151 -> 27_fp -> L
    feats[726] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '151') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '151', '27_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_fp', node2.lemma)
    
    # L -> 151 -> -97 -> L
    feats[727] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '151') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '151', '-97') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-97', node2.lemma)
    
    # L -> 151 -> -10 -> L
    feats[728] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '151') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '151', '-10') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-10', node2.lemma)
    
    # L -> 151 -> -241 -> L
    feats[729] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '151') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '151', '-241') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-241', node2.lemma)
    
    # L -> 151 -> -33 -> L
    feats[730] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '151') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '151', '-33') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-33', node2.lemma)
    
    # L -> 151 -> -169 -> L
    feats[731] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '151') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '151', '-169') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-169', node2.lemma)
    
    # L -> 151 -> nom. fem -> L
    feats[732] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '151') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '151', 'nom. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. fem', node2.lemma)
    
    # L -> 151 -> 128 -> L
    feats[733] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '151') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '151', '128') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '128', node2.lemma)
    
    # L -> 151 -> -11 -> L
    feats[734] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '151') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '151', '-11') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-11', node2.lemma)
    
    # L -> 151 -> -71 -> L
    feats[735] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '151') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '151', '-71') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-71', node2.lemma)
    
    # L -> 151 -> sg_tp -> L
    feats[736] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '151') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '151', 'sg_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_tp', node2.lemma)
    
    # L -> 151 -> 9_sp -> L
    feats[737] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '151') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '151', '9_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_sp', node2.lemma)
    
    # L -> 151 -> 5_sp -> L
    feats[738] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '151') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '151', '5_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sp', node2.lemma)
    
    # L -> 151 -> 150 -> L
    feats[739] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '151') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '151', '150') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '150', node2.lemma)
    
    # L -> 151 -> -161 -> L
    feats[740] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '151') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '151', '-161') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-161', node2.lemma)
    
    # L -> 151 -> pl_fp -> L
    feats[741] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '151') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '151', 'pl_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_fp', node2.lemma)
    
    # L -> 151 -> 14_sp -> L
    feats[742] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '151') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '151', '14_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_sp', node2.lemma)
    
    # L -> 151 -> 77 -> L
    feats[743] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '151') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '151', '77') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '77', node2.lemma)
    
    # L -> 151 -> 79 -> L
    feats[744] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '151') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '151', '79') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '79', node2.lemma)
    
    # L -> 151 -> 11_du -> L
    feats[745] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '151') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '151', '11_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_du', node2.lemma)
    
    # L -> 151 -> gen. sg. -> L
    feats[746] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '151') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '151', 'gen. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen. sg.', node2.lemma)
    
    # L -> 151 -> 10_du -> L
    feats[747] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '151') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '151', '10_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_du', node2.lemma)
    
    # L -> 151 -> 7_sp -> L
    feats[748] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '151') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '151', '7_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_sp', node2.lemma)
    
    # L -> 98 -> 149 -> L
    feats[749] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '98', '149') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '149', node2.lemma)
    
    # L -> 98 -> -14 -> L
    feats[750] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '98', '-14') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-14', node2.lemma)
    
    # L -> 98 -> 171 -> L
    feats[751] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '98', '171') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '171', node2.lemma)
    
    # L -> 98 -> fp -> L
    feats[752] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '98', 'fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'fp', node2.lemma)
    
    # L -> 98 -> 13_tp -> L
    feats[753] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '98', '13_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_tp', node2.lemma)
    
    # L -> 98 -> -44 -> L
    feats[754] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '98', '-44') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-44', node2.lemma)
    
    # L -> 98 -> 7_du -> L
    feats[755] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '98', '7_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_du', node2.lemma)
    
    # L -> 98 -> 70 -> L
    feats[756] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '98', '70') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '70', node2.lemma)
    
    # L -> 98 -> -92 -> L
    feats[757] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '98', '-92') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-92', node2.lemma)
    
    # L -> 98 -> 30 -> L
    feats[758] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '98', '30') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30', node2.lemma)
    
    # L -> 98 -> -113 -> L
    feats[759] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '98', '-113') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-113', node2.lemma)
    
    # L -> 98 -> 117 -> L
    feats[760] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '98', '117') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '117', node2.lemma)
    
    # L -> 98 -> -158 -> L
    feats[761] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '98', '-158') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-158', node2.lemma)
    
    # L -> 98 -> 135 -> L
    feats[762] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '98', '135') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '135', node2.lemma)
    
    # L -> 98 -> -307 -> L
    feats[763] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '98', '-307') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-307', node2.lemma)
    
    # L -> 98 -> 81 -> L
    feats[764] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '98', '81') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '81', node2.lemma)
    
    # L -> 98 -> -51 -> L
    feats[765] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '98', '-51') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-51', node2.lemma)
    
    # L -> 98 -> 6_sp -> L
    feats[766] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '98', '6_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_sp', node2.lemma)
    
    # L -> 98 -> 42 -> L
    feats[767] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '98', '42') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '42', node2.lemma)
    
    # L -> 98 -> 32 -> L
    feats[768] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '98', '32') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '32', node2.lemma)
    
    # L -> 98 -> dat. pl. -> L
    feats[769] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '98', 'dat. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. pl.', node2.lemma)
    
    # L -> 98 -> 7_fp -> L
    feats[770] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '98', '7_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_fp', node2.lemma)
    
    # L -> 98 -> 139 -> L
    feats[771] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '98', '139') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '139', node2.lemma)
    
    # L -> 98 -> -276 -> L
    feats[772] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '98', '-276') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-276', node2.lemma)
    
    # L -> 98 -> 6_tp -> L
    feats[773] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '98', '6_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_tp', node2.lemma)
    
    # L -> 98 -> -76 -> L
    feats[774] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '98', '-76') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-76', node2.lemma)
    
    # L -> 98 -> 51 -> L
    feats[775] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '98', '51') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '51', node2.lemma)
    
    # L -> 98 -> 99 -> L
    feats[776] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '98', '99') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '99', node2.lemma)
    
    # L -> 98 -> 8_fp -> L
    feats[777] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '98', '8_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    
    # L -> 98 -> -61 -> L
    feats[778] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '98', '-61') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-61', node2.lemma)
    
    # L -> 177 -> 10_tp -> L
    feats[779] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '177') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '177', '10_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_tp', node2.lemma)
    
    # L -> 2_fp -> abl. pl. -> L
    feats[780] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', 'abl. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. pl.', node2.lemma)
    
    # L -> 2_fp -> -147 -> L
    feats[781] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', '-147') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-147', node2.lemma)
    
    # L -> 2_fp -> 39 -> L
    feats[782] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', '39') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '39', node2.lemma)
    
    # L -> 2_fp -> -246 -> L
    feats[783] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', '-246') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-246', node2.lemma)
    
    # L -> 2_fp -> 2_sg -> L
    feats[784] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', '2_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_sg', node2.lemma)
    
    # L -> 2_fp -> -122 -> L
    feats[785] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', '-122') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-122', node2.lemma)
    
    # L -> 2_fp -> -220 -> L
    feats[786] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', '-220') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-220', node2.lemma)
    
    # L -> 2_fp -> nom. fem -> L
    feats[787] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', 'nom. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. fem', node2.lemma)
    
    # L -> 2_fp -> -48 -> L
    feats[788] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', '-48') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-48', node2.lemma)
    
    # L -> 2_fp -> -261 -> L
    feats[789] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', '-261') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-261', node2.lemma)
    
    # L -> 2_fp -> 170 -> L
    feats[790] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', '170') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '170', node2.lemma)
    
    # L -> 2_fp -> 9_sp -> L
    feats[791] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', '9_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_sp', node2.lemma)
    
    # L -> 2_fp -> 5_sp -> L
    feats[792] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', '5_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sp', node2.lemma)
    
    # L -> 2_fp -> 14_sp -> L
    feats[793] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', '14_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_sp', node2.lemma)
    
    # L -> 2_fp -> 179 -> L
    feats[794] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', '179') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '179', node2.lemma)
    
    # L -> 2_fp -> tp -> L
    feats[795] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', 'tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'tp', node2.lemma)
    
    # L -> 2_fp -> -142 -> L
    feats[796] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', '-142') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-142', node2.lemma)
    
    # L -> 2_fp -> -271 -> L
    feats[797] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', '-271') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-271', node2.lemma)
    
    # L -> 2_fp -> -22 -> L
    feats[798] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', '-22') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-22', node2.lemma)
    
    # L -> 2_fp -> dat -> L
    feats[799] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', 'dat') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat', node2.lemma)
    
    # L -> 2_fp -> adj -> L
    feats[800] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', 'adj') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'adj', node2.lemma)
    
    # L -> 2_fp -> 12_tp -> L
    feats[801] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', '12_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_tp', node2.lemma)
    
    # L -> 2_fp -> 9_tp -> L
    feats[802] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', '9_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_tp', node2.lemma)
    
    # L -> 2_fp -> 29_sg -> L
    feats[803] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', '29_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '29_sg', node2.lemma)
    
    # L -> 2_fp -> 75 -> L
    feats[804] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', '75') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '75', node2.lemma)
    
    # L -> 2_fp -> nom. sg. -> L
    feats[805] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', 'nom. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. sg.', node2.lemma)
    
    # L -> 2_fp -> 90 -> L
    feats[806] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', '90') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '90', node2.lemma)
    
    # L -> 2_fp -> 2_pl -> L
    feats[807] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', '2_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_pl', node2.lemma)
    
    # L -> 2_fp -> 156 -> L
    feats[808] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', '156') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '156', node2.lemma)
    
    # L -> 2_fp -> 171 -> L
    feats[809] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', '171') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '171', node2.lemma)
    
    # L -> 2_fp -> fp -> L
    feats[810] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', 'fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'fp', node2.lemma)
    
    # L -> 2_fp -> -23 -> L
    feats[811] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', '-23') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-23', node2.lemma)
    
    # L -> 2_fp -> -21 -> L
    feats[812] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', '-21') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-21', node2.lemma)
    
    # L -> 2_fp -> 89 -> L
    feats[813] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', '89') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '89', node2.lemma)
    
    # L -> 2_fp -> instr. masc. -> L
    feats[814] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', 'instr. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. masc.', node2.lemma)
    
    # L -> 2_fp -> -306 -> L
    feats[815] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', '-306') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-306', node2.lemma)
    
    # L -> 2_fp -> 16_tp -> L
    feats[816] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', '16_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_tp', node2.lemma)
    
    # L -> 2_fp -> 6_du -> L
    feats[817] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', '6_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_du', node2.lemma)
    
    # L -> 5_pl -> 110 -> L
    feats[818] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_pl', '110') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '110', node2.lemma)
    
    # L -> 5_pl -> -220 -> L
    feats[819] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_pl', '-220') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-220', node2.lemma)
    
    # L -> 5_pl -> 137 -> L
    feats[820] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_pl', '137') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    
    # L -> 5_pl -> nom. fem -> L
    feats[821] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_pl', 'nom. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. fem', node2.lemma)
    
    # L -> 5_pl -> sg_tp -> L
    feats[822] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_pl', 'sg_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_tp', node2.lemma)
    
    # L -> 5_pl -> -161 -> L
    feats[823] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_pl', '-161') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-161', node2.lemma)
    
    # L -> 5_pl -> 79 -> L
    feats[824] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_pl', '79') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '79', node2.lemma)
    
    # L -> 5_pl -> 6_pl -> L
    feats[825] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_pl', '6_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_pl', node2.lemma)
    
    # L -> 5_pl -> voc. fem -> L
    feats[826] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_pl', 'voc. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. fem', node2.lemma)
    
    # L -> 5_pl -> -152 -> L
    feats[827] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_pl', '-152') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-152', node2.lemma)
    
    # L -> 5_pl -> tp -> L
    feats[828] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_pl', 'tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'tp', node2.lemma)
    
    # L -> 5_pl -> 10_fp -> L
    feats[829] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_pl', '10_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_fp', node2.lemma)
    
    # L -> 5_pl -> -271 -> L
    feats[830] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_pl', '-271') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-271', node2.lemma)
    
    # L -> 5_pl -> 10_tp -> L
    feats[831] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_pl', '10_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_tp', node2.lemma)
    
    # L -> 5_pl -> -93 -> L
    feats[832] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_pl', '-93') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-93', node2.lemma)
    
    # L -> 5_pl -> 177 -> L
    feats[833] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_pl', '177') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '177', node2.lemma)
    
    # L -> 5_pl -> -230 -> L
    feats[834] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_pl', '-230') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-230', node2.lemma)
    
    # L -> 5_pl -> nom. sg. -> L
    feats[835] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_pl', 'nom. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. sg.', node2.lemma)
    
    # L -> 5_pl -> 2_pl -> L
    feats[836] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_pl', '2_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_pl', node2.lemma)
    
    # L -> 5_pl -> -82 -> L
    feats[837] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_pl', '-82') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-82', node2.lemma)
    
    # L -> 5_pl -> -112 -> L
    feats[838] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_pl', '-112') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    
    # L -> 5_pl -> -266 -> L
    feats[839] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_pl', '-266') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-266', node2.lemma)
    
    # L -> 5_pl -> loc -> L
    feats[840] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_pl', 'loc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc', node2.lemma)
    
    # L -> 5_pl -> 149 -> L
    feats[841] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_pl', '149') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '149', node2.lemma)
    
    # L -> 5_pl -> 49 -> L
    feats[842] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_pl', '49') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '49', node2.lemma)
    
    # L -> 5_pl -> -44 -> L
    feats[843] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_pl', '-44') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-44', node2.lemma)
    
    # L -> 5_pl -> -115 -> L
    feats[844] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_pl', '-115') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-115', node2.lemma)
    
    # L -> 5_pl -> -292 -> L
    feats[845] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_pl', '-292') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-292', node2.lemma)
    
    # L -> 5_pl -> -92 -> L
    feats[846] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_pl', '-92') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-92', node2.lemma)
    
    # L -> 5_pl -> 6_tp -> L
    feats[847] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_pl', '6_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_tp', node2.lemma)
    
    # L -> 12_tp -> 6_pl -> L
    feats[848] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_tp', '6_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_pl', node2.lemma)
    
    # L -> 12_tp -> -112 -> L
    feats[849] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_tp', '-112') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
    
    # L -> 12_tp -> -115 -> L
    feats[850] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_tp', '-115') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-115', node2.lemma)
    
    # L -> 12_tp -> 152 -> L
    feats[851] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_tp', '152') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '152', node2.lemma)
    
    # L -> 12_tp -> 68 -> L
    feats[852] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_tp', '68') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '68', node2.lemma)
    
    # L -> 12_tp -> masc -> L
    feats[853] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_tp', 'masc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'masc', node2.lemma)
    
    # L -> 7_tp -> 179 -> L
    feats[854] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_tp', '179') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '179', node2.lemma)
    
    # L -> 7_tp -> 7_sp -> L
    feats[855] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_tp', '7_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_sp', node2.lemma)
    
    # L -> 7_tp -> pl_sp -> L
    feats[856] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_tp', 'pl_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_sp', node2.lemma)
    
    # L -> 7_tp -> 154 -> L
    feats[857] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_tp', '154') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '154', node2.lemma)
    
    # L -> 7_tp -> abl -> L
    feats[858] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_tp', 'abl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl', node2.lemma)
    
    # L -> 7_tp -> 108 -> L
    feats[859] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_tp', '108') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '108', node2.lemma)
    
    # L -> 7_tp -> -35 -> L
    feats[860] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_tp', '-35') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    
    # L -> 7_tp -> -37 -> L
    feats[861] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_tp', '-37') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-37', node2.lemma)
    
    # L -> 7_tp -> 101 -> L
    feats[862] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_tp', '101') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    
    # L -> 7_tp -> -68 -> L
    feats[863] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_tp', '-68') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-68', node2.lemma)
    
    # L -> 7_tp -> -131 -> L
    feats[864] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_tp', '-131') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-131', node2.lemma)
    
    # L -> 7_tp -> 3_pl -> L
    feats[865] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_tp', '3_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_pl', node2.lemma)
    
    # L -> 7_tp -> 48 -> L
    feats[866] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_tp', '48') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '48', node2.lemma)
    
    # L -> 7_tp -> 30_du -> L
    feats[867] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_tp', '30_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_du', node2.lemma)
    
    # L -> 7_tp -> sg -> L
    feats[868] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_tp', 'sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg', node2.lemma)
    
    # L -> 7_tp -> 161 -> L
    feats[869] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_tp', '161') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '161', node2.lemma)
    
    # L -> 7_tp -> 10_pl -> L
    feats[870] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_tp', '10_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_pl', node2.lemma)
    
    # L -> 7_tp -> -84 -> L
    feats[871] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_tp', '-84') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-84', node2.lemma)
    
    # L -> 7_tp -> nom -> L
    feats[872] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_tp', 'nom') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom', node2.lemma)
    
    # L -> 7_tp -> -129 -> L
    feats[873] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_tp', '-129') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-129', node2.lemma)
    
    # L -> 7_tp -> nom. masc. -> L
    feats[874] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_tp', 'nom. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. masc.', node2.lemma)
    
    # L -> 7_tp -> 1 -> L
    feats[875] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_tp', '1') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '1', node2.lemma)
    
    # L -> -230 -> abl. pl. -> L
    feats[876] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', 'abl. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. pl.', node2.lemma)
    
    # L -> -230 -> 4_pl -> L
    feats[877] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '4_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_pl', node2.lemma)
    
    # L -> -230 -> -246 -> L
    feats[878] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '-246') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-246', node2.lemma)
    
    # L -> -230 -> du_tp -> L
    feats[879] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', 'du_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_tp', node2.lemma)
    
    # L -> -230 -> -163 -> L
    feats[880] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '-163') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-163', node2.lemma)
    
    # L -> -230 -> 88 -> L
    feats[881] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '88') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '88', node2.lemma)
    
    # L -> -230 -> du -> L
    feats[882] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', 'du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du', node2.lemma)
    
    # L -> -230 -> 13_sg -> L
    feats[883] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '13_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_sg', node2.lemma)
    
    # L -> -230 -> -97 -> L
    feats[884] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '-97') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-97', node2.lemma)
    
    # L -> -230 -> 159 -> L
    feats[885] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '159') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '159', node2.lemma)
    
    # L -> -230 -> 137 -> L
    feats[886] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '137') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    
    # L -> -230 -> sg_tp -> L
    feats[887] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', 'sg_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_tp', node2.lemma)
    
    # L -> -230 -> -261 -> L
    feats[888] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '-261') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-261', node2.lemma)
    
    # L -> -230 -> 56 -> L
    feats[889] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '56') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '56', node2.lemma)
    
    # L -> -230 -> 170 -> L
    feats[890] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '170') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '170', node2.lemma)
    
    # L -> -230 -> 150 -> L
    feats[891] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '150') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '150', node2.lemma)
    
    # L -> -230 -> 77 -> L
    feats[892] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '77') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '77', node2.lemma)
    
    # L -> -230 -> 79 -> L
    feats[893] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '79') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '79', node2.lemma)
    
    # L -> -230 -> -273 -> L
    feats[894] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '-273') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-273', node2.lemma)
    
    # L -> -230 -> 27_sg -> L
    feats[895] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '27_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_sg', node2.lemma)
    
    # L -> -230 -> -166 -> L
    feats[896] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '-166') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-166', node2.lemma)
    
    # L -> -230 -> dat -> L
    feats[897] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', 'dat') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat', node2.lemma)
    
    # L -> -230 -> 151 -> L
    feats[898] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '151') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '151', node2.lemma)
    
    # L -> -230 -> 7_tp -> L
    feats[899] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '7_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_tp', node2.lemma)
    
    # L -> -230 -> -52 -> L
    feats[900] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '-52') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-52', node2.lemma)
    
    # L -> -230 -> 90 -> L
    feats[901] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '90') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '90', node2.lemma)
    
    # L -> -230 -> acc. adj. -> L
    feats[902] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', 'acc. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. adj.', node2.lemma)
    
    # L -> -230 -> 4_du -> L
    feats[903] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '4_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_du', node2.lemma)
    
    # L -> -230 -> 11_pl -> L
    feats[904] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '11_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_pl', node2.lemma)
    
    # L -> -230 -> 156 -> L
    feats[905] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '156') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '156', node2.lemma)
    
    # L -> -230 -> 118 -> L
    feats[906] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '118') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '118', node2.lemma)
    
    # L -> -230 -> -14 -> L
    feats[907] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '-14') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-14', node2.lemma)
    
    # L -> -230 -> voc. sg. -> L
    feats[908] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', 'voc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. sg.', node2.lemma)
    
    # L -> -230 -> 13_tp -> L
    feats[909] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '13_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_tp', node2.lemma)
    
    # L -> -230 -> -73 -> L
    feats[910] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '-73') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-73', node2.lemma)
    
    # L -> -230 -> 12_pl -> L
    feats[911] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '12_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_pl', node2.lemma)
    
    # L -> -230 -> -27 -> L
    feats[912] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '-27') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-27', node2.lemma)
    
    # L -> -230 -> -144 -> L
    feats[913] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '-144') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-144', node2.lemma)
    
    # L -> -230 -> 4_sp -> L
    feats[914] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '4_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_sp', node2.lemma)
    
    # L -> -230 -> -143 -> L
    feats[915] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '-143') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-143', node2.lemma)
    
    # L -> -230 -> 9_fp -> L
    feats[916] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '9_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_fp', node2.lemma)
    
    # L -> -230 -> 135 -> L
    feats[917] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '135') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '135', node2.lemma)
    
    # L -> -230 -> 3_fp -> L
    feats[918] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '3_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_fp', node2.lemma)
    
    # L -> -230 -> -260 -> L
    feats[919] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '-260') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-260', node2.lemma)
    
    # L -> -230 -> -200 -> L
    feats[920] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '-200') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-200', node2.lemma)
    
    # L -> -230 -> -51 -> L
    feats[921] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '-51') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-51', node2.lemma)
    
    # L -> -230 -> -101 -> L
    feats[922] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '-101') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-101', node2.lemma)
    
    # L -> -230 -> 6_sp -> L
    feats[923] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '6_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_sp', node2.lemma)
    
    # L -> -230 -> 42 -> L
    feats[924] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '42') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '42', node2.lemma)
    
    # L -> -230 -> 11_tp -> L
    feats[925] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '11_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_tp', node2.lemma)
    
    # L -> -230 -> 97 -> L
    feats[926] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '97') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '97', node2.lemma)
    
    # L -> -230 -> -268 -> L
    feats[927] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '-268') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-268', node2.lemma)
    
    # L -> -230 -> -153 -> L
    feats[928] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '-153') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-153', node2.lemma)
    
    # L -> -230 -> 169 -> L
    feats[929] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '169') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '169', node2.lemma)
    
    # L -> -230 -> 3_sg -> L
    feats[930] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '3_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_sg', node2.lemma)
    
    # L -> -230 -> -262 -> L
    feats[931] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '-262') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-262', node2.lemma)
    
    # L -> -230 -> 8_sp -> L
    feats[932] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '8_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_sp', node2.lemma)
    
    # L -> -230 -> 112 -> L
    feats[933] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '112') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '112', node2.lemma)
    
    # L -> -230 -> voc -> L
    feats[934] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', 'voc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc', node2.lemma)
    
    # L -> -230 -> -249 -> L
    feats[935] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '-249') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-249', node2.lemma)
    
    # L -> -230 -> -81 -> L
    feats[936] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '-81') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-81', node2.lemma)
    
    # L -> -230 -> -190 -> L
    feats[937] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '-190') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-190', node2.lemma)
    
    # L -> -230 -> 29_tp -> L
    feats[938] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '29_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '29_tp', node2.lemma)
    
    # L -> -230 -> nom. pl. -> L
    feats[939] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', 'nom. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. pl.', node2.lemma)
    
    # L -> -230 -> 30_fp -> L
    feats[940] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '30_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_fp', node2.lemma)
    
    # L -> -230 -> -46 -> L
    feats[941] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '-46') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-46', node2.lemma)
    
    # L -> -230 -> 108 -> L
    feats[942] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', '108') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '108', node2.lemma)
    
    # L -> -230 -> voc. pl. -> L
    feats[943] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', 'voc. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. pl.', node2.lemma)
    
    # L -> loc. du. -> -246 -> L
    feats[944] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. du.', '-246') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-246', node2.lemma)
    
    # L -> loc. du. -> -157 -> L
    feats[945] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. du.', '-157') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-157', node2.lemma)
    
    # L -> loc. du. -> 10_tp -> L
    feats[946] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. du.', '10_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_tp', node2.lemma)
    
    # L -> -34 -> -271 -> L
    feats[947] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '-271') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-271', node2.lemma)
    
    # L -> -34 -> 116 -> L
    feats[948] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '116') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '116', node2.lemma)
    
    # L -> -34 -> -24 -> L
    feats[949] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '-24') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
    
    # L -> -34 -> -200 -> L
    feats[950] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '-200') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-200', node2.lemma)
    
    # L -> -34 -> -51 -> L
    feats[951] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '-51') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-51', node2.lemma)
    
    # L -> -34 -> 30_pl -> L
    feats[952] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '30_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_pl', node2.lemma)
    
    # L -> -34 -> -63 -> L
    feats[953] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '-63') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-63', node2.lemma)
    
    # L -> -34 -> sg_sp -> L
    feats[954] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', 'sg_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_sp', node2.lemma)
    
    # L -> -34 -> 27_tp -> L
    feats[955] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '27_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_tp', node2.lemma)
    
    # L -> -34 -> 15_sp -> L
    feats[956] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '15_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_sp', node2.lemma)
    
    # L -> -34 -> -276 -> L
    feats[957] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '-276') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-276', node2.lemma)
    
    # L -> -34 -> 73 -> L
    feats[958] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '73') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '73', node2.lemma)
    
    # L -> -34 -> 152 -> L
    feats[959] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '152') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '152', node2.lemma)
    
    # L -> -34 -> 7_pl -> L
    feats[960] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '7_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_pl', node2.lemma)
    
    # L -> -34 -> -102 -> L
    feats[961] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '-102') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-102', node2.lemma)
    
    # L -> -34 -> 97 -> L
    feats[962] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '97') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '97', node2.lemma)
    
    # L -> -34 -> -268 -> L
    feats[963] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '-268') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-268', node2.lemma)
    
    # L -> -34 -> -39 -> L
    feats[964] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '-39') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-39', node2.lemma)
    
    # L -> -34 -> 119 -> L
    feats[965] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '119') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '119', node2.lemma)
    
    # L -> -34 -> -153 -> L
    feats[966] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '-153') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-153', node2.lemma)
    
    # L -> -34 -> voc. neutr. -> L
    feats[967] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', 'voc. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. neutr.', node2.lemma)
    
    # L -> -34 -> gen. du. -> L
    feats[968] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', 'gen. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen. du.', node2.lemma)
    
    # L -> -34 -> -62 -> L
    feats[969] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '-62') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-62', node2.lemma)
    
    # L -> -34 -> -76 -> L
    feats[970] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '-76') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-76', node2.lemma)
    
    # L -> -34 -> 11_du -> L
    feats[971] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '11_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_du', node2.lemma)
    
    # L -> -34 -> 8_fp -> L
    feats[972] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '8_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    
    # L -> -34 -> acc -> L
    feats[973] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', 'acc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc', node2.lemma)
    
    # L -> -34 -> 80 -> L
    feats[974] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '80') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '80', node2.lemma)
    
    # L -> -34 -> 16_fp -> L
    feats[975] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '16_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_fp', node2.lemma)
    
    # L -> -34 -> 59 -> L
    feats[976] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '59') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '59', node2.lemma)
    
    # L -> -34 -> 15_sg -> L
    feats[977] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '15_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_sg', node2.lemma)
    
    # L -> -34 -> 109 -> L
    feats[978] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '109') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '109', node2.lemma)
    
    # L -> -34 -> 3_sg -> L
    feats[979] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '3_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_sg', node2.lemma)
    
    # L -> -34 -> gen. sg. -> L
    feats[980] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', 'gen. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen. sg.', node2.lemma)
    
    # L -> -34 -> sp -> L
    feats[981] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', 'sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sp', node2.lemma)
    
    # L -> -34 -> instr. fem -> L
    feats[982] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', 'instr. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. fem', node2.lemma)
    
    # L -> -34 -> 8_sp -> L
    feats[983] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '8_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_sp', node2.lemma)
    
    # L -> -34 -> 38 -> L
    feats[984] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '38') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '38', node2.lemma)
    
    # L -> -34 -> 148 -> L
    feats[985] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '148') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '148', node2.lemma)
    
    # L -> -34 -> 112 -> L
    feats[986] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '112') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '112', node2.lemma)
    
    # L -> -34 -> -137 -> L
    feats[987] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '-137') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-137', node2.lemma)
    
    # L -> -34 -> -119 -> L
    feats[988] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '-119') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-119', node2.lemma)
    
    # L -> -34 -> 5_sg -> L
    feats[989] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '5_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sg', node2.lemma)
    
    # L -> -34 -> 181 -> L
    feats[990] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '181') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '181', node2.lemma)
    
    # L -> -34 -> voc -> L
    feats[991] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', 'voc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc', node2.lemma)
    
    # L -> -34 -> 68 -> L
    feats[992] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '68') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '68', node2.lemma)
    
    # L -> -34 -> 15_pl -> L
    feats[993] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '15_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_pl', node2.lemma)
    
    # L -> -34 -> -83 -> L
    feats[994] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '-83') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-83', node2.lemma)
    
    # L -> -34 -> -249 -> L
    feats[995] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '-249') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-249', node2.lemma)
    
    # L -> -34 -> 13_pl -> L
    feats[996] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '13_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_pl', node2.lemma)
    
    # L -> -34 -> nom. adj. -> L
    feats[997] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', 'nom. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. adj.', node2.lemma)
    
    # L -> -34 -> 16_tp -> L
    feats[998] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '16_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_tp', node2.lemma)
    
    # L -> -34 -> 6_du -> L
    feats[999] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '6_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_du', node2.lemma)
    
    # L -> -34 -> -54 -> L
    feats[1000] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '-54') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-54', node2.lemma)
    
    # L -> -34 -> -41 -> L
    feats[1001] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '-41') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-41', node2.lemma)
    
    # L -> -34 -> 134 -> L
    feats[1002] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '134') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '134', node2.lemma)
    
    # L -> -34 -> nom. pl. -> L
    feats[1003] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', 'nom. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. pl.', node2.lemma)
    
    # L -> -34 -> 30_fp -> L
    feats[1004] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '30_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_fp', node2.lemma)
    
    # L -> -34 -> -299 -> L
    feats[1005] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '-299') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-299', node2.lemma)
    
    # L -> -34 -> 4_sg -> L
    feats[1006] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '4_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_sg', node2.lemma)
    
    # L -> -34 -> -18 -> L
    feats[1007] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '-18') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-18', node2.lemma)
    
    # L -> -34 -> 61 -> L
    feats[1008] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '61') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '61', node2.lemma)
    
    # L -> -34 -> pl_sp -> L
    feats[1009] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', 'pl_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_sp', node2.lemma)
    
    # L -> -34 -> 69 -> L
    feats[1010] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '69') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '69', node2.lemma)
    
    # L -> -34 -> 155 -> L
    feats[1011] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '155') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '155', node2.lemma)
    
    # L -> -34 -> instr. adj. -> L
    feats[1012] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', 'instr. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. adj.', node2.lemma)
    
    # L -> -34 -> 82 -> L
    feats[1013] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '82') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '82', node2.lemma)
    
    # L -> -34 -> 168 -> L
    feats[1014] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '168') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '168', node2.lemma)
    
    # L -> -34 -> -57 -> L
    feats[1015] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '-57') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-57', node2.lemma)
    
    # L -> -34 -> -46 -> L
    feats[1016] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '-46') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-46', node2.lemma)
    
    # L -> -34 -> 108 -> L
    feats[1017] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', '108') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '108', node2.lemma)
    
    # L -> -132 -> 75 -> L
    feats[1018] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-132') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-132', '75') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '75', node2.lemma)
    
    # L -> -132 -> 116 -> L
    feats[1019] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-132') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-132', '116') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '116', node2.lemma)
    
    # L -> -132 -> gen. sg. -> L
    feats[1020] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-132') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-132', 'gen. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen. sg.', node2.lemma)
    
    # L -> -132 -> 13_fp -> L
    feats[1021] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-132') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-132', '13_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_fp', node2.lemma)
    
    # L -> -58 -> -10 -> L
    feats[1022] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-58', '-10') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-10', node2.lemma)
    
    # L -> -58 -> -93 -> L
    feats[1023] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-58', '-93') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-93', node2.lemma)
    
    # L -> -58 -> -139 -> L
    feats[1024] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-58', '-139') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-139', node2.lemma)
    
    # L -> -58 -> 151 -> L
    feats[1025] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-58', '151') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '151', node2.lemma)
    
    # L -> -58 -> acc. adj. -> L
    feats[1026] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-58', 'acc. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. adj.', node2.lemma)
    
    # L -> -58 -> 120 -> L
    feats[1027] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-58', '120') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '120', node2.lemma)
    
    # L -> -58 -> 156 -> L
    feats[1028] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-58', '156') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '156', node2.lemma)
    
    # L -> -58 -> -156 -> L
    feats[1029] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-58', '-156') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-156', node2.lemma)
    
    # L -> -58 -> -50 -> L
    feats[1030] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-58', '-50') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-50', node2.lemma)
    
    # L -> -58 -> -99 -> L
    feats[1031] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-58', '-99') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-99', node2.lemma)
    
    # L -> -58 -> dat. pl. -> L
    feats[1032] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-58', 'dat. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. pl.', node2.lemma)
    
    # L -> -58 -> 73 -> L
    feats[1033] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-58', '73') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '73', node2.lemma)
    
    # L -> -58 -> 173 -> L
    feats[1034] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-58', '173') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '173', node2.lemma)
    
    # L -> -58 -> 78 -> L
    feats[1035] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-58', '78') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '78', node2.lemma)
    
    # L -> -58 -> 28 -> L
    feats[1036] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-58', '28') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
    
    # L -> -58 -> -210 -> L
    feats[1037] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-58', '-210') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-210', node2.lemma)
    
    # L -> -58 -> 33 -> L
    feats[1038] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-58', '33') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '33', node2.lemma)
    
    # L -> -52 -> acc. masc. -> L
    feats[1039] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-52') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-52', 'acc. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. masc.', node2.lemma)
    
    # L -> -52 -> -39 -> L
    feats[1040] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-52') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-52', '-39') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-39', node2.lemma)
    
    # L -> -52 -> 16_fp -> L
    feats[1041] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-52') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-52', '16_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_fp', node2.lemma)
    
    # L -> -52 -> 14_sg -> L
    feats[1042] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-52') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-52', '14_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_sg', node2.lemma)
    
    # L -> -52 -> acc. du. -> L
    feats[1043] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-52') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-52', 'acc. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. du.', node2.lemma)
    
    # L -> -52 -> 175 -> L
    feats[1044] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-52') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-52', '175') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '175', node2.lemma)
    
    # L -> -52 -> masc -> L
    feats[1045] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-52') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-52', 'masc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'masc', node2.lemma)
    
    # L -> -52 -> -38 -> L
    feats[1046] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-52') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-52', '-38') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-38', node2.lemma)
    
    # L -> 75 -> -121 -> L
    feats[1047] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '75') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '75', '-121') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-121', node2.lemma)
    
    # L -> sg_fp -> -152 -> L
    feats[1048] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg_fp', '-152') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-152', node2.lemma)
    
    # L -> sg_fp -> 177 -> L
    feats[1049] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg_fp', '177') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '177', node2.lemma)
    
    # L -> nom. sg. -> -44 -> L
    feats[1050] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. sg.', '-44') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-44', node2.lemma)
    
    # L -> nom. sg. -> pl_tp -> L
    feats[1051] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. sg.', 'pl_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_tp', node2.lemma)
    
    # L -> nom. sg. -> 15_tp -> L
    feats[1052] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. sg.', '15_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_tp', node2.lemma)
    
    # L -> nom. sg. -> 175 -> L
    feats[1053] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. sg.', '175') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '175', node2.lemma)
    
    # L -> nom. sg. -> pl_sp -> L
    feats[1054] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. sg.', 'pl_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_sp', node2.lemma)
    
    # L -> nom. sg. -> 168 -> L
    feats[1055] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. sg.', '168') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '168', node2.lemma)
    
    # L -> nom. sg. -> 160 -> L
    feats[1056] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. sg.', '160') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
    
    # L -> nom. sg. -> 3_pl -> L
    feats[1057] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. sg.', '3_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_pl', node2.lemma)
    
    # L -> nom. sg. -> 48 -> L
    feats[1058] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. sg.', '48') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '48', node2.lemma)
    
    # L -> nom. sg. -> 30_du -> L
    feats[1059] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. sg.', '30_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_du', node2.lemma)
    
    # L -> nom. sg. -> 176 -> L
    feats[1060] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. sg.', '176') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '176', node2.lemma)
    
    # L -> nom. sg. -> -129 -> L
    feats[1061] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. sg.', '-129') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-129', node2.lemma)
    
    # L -> nom. sg. -> -303 -> L
    feats[1062] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. sg.', '-303') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-303', node2.lemma)
    
    # L -> 90 -> du -> L
    feats[1063] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '90', 'du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du', node2.lemma)
    
    # L -> 90 -> loc. sg. -> L
    feats[1064] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '90', 'loc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc. sg.', node2.lemma)
    
    # L -> 90 -> pl_fp -> L
    feats[1065] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '90', 'pl_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_fp', node2.lemma)
    
    # L -> 90 -> 10_tp -> L
    feats[1066] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '90', '10_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_tp', node2.lemma)
    
    # L -> 90 -> 98 -> L
    feats[1067] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '90', '98') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '98', node2.lemma)
    
    # L -> 90 -> 177 -> L
    feats[1068] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '90', '177') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '177', node2.lemma)
    
    # L -> 90 -> 12_tp -> L
    feats[1069] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '90', '12_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_tp', node2.lemma)
    
    # L -> 90 -> 99 -> L
    feats[1070] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '90', '99') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '99', node2.lemma)
    
    # L -> 90 -> 36 -> L
    feats[1071] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '90', '36') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '36', node2.lemma)
    
    # L -> 90 -> -61 -> L
    feats[1072] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '90', '-61') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-61', node2.lemma)
    
    # L -> 90 -> sp -> L
    feats[1073] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '90', 'sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sp', node2.lemma)
    
    # L -> 90 -> -66 -> L
    feats[1074] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '90', '-66') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-66', node2.lemma)
    
    # L -> 90 -> 148 -> L
    feats[1075] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '90', '148') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '148', node2.lemma)
    
    # L -> 90 -> 60 -> L
    feats[1076] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '90', '60') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '60', node2.lemma)
    
    # L -> 90 -> 14_pl -> L
    feats[1077] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '90', '14_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_pl', node2.lemma)
    
    # L -> 90 -> acc. du. -> L
    feats[1078] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '90', 'acc. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. du.', node2.lemma)
    
    # L -> 90 -> -249 -> L
    feats[1079] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '90', '-249') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-249', node2.lemma)
    
    # L -> 90 -> 15_du -> L
    feats[1080] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '90', '15_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_du', node2.lemma)
    
    # L -> 90 -> masc -> L
    feats[1081] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '90', 'masc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'masc', node2.lemma)
    
    # L -> 90 -> 82 -> L
    feats[1082] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '90', '82') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '82', node2.lemma)
    
    # L -> 90 -> -57 -> L
    feats[1083] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '90', '-57') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-57', node2.lemma)
    
    # L -> 90 -> -46 -> L
    feats[1084] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '90', '-46') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-46', node2.lemma)
    
    # L -> 90 -> 108 -> L
    feats[1085] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '90', '108') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '108', node2.lemma)
    
    # L -> 90 -> -131 -> L
    feats[1086] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '90', '-131') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-131', node2.lemma)
    
    # L -> 90 -> 158 -> L
    feats[1087] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '90', '158') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '158', node2.lemma)
    
    # L -> 90 -> sg -> L
    feats[1088] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '90', 'sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg', node2.lemma)
    
    # L -> 2_pl -> -147 -> L
    feats[1089] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_pl', '-147') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-147', node2.lemma)
    
    # L -> 2_pl -> -114 -> L
    feats[1090] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_pl', '-114') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-114', node2.lemma)
    
    # L -> 2_pl -> 39 -> L
    feats[1091] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_pl', '39') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '39', node2.lemma)
    
    # L -> 2_pl -> 72 -> L
    feats[1092] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_pl', '72') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '72', node2.lemma)
    
    # L -> 2_pl -> -309 -> L
    feats[1093] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_pl', '-309') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-309', node2.lemma)
    
    # L -> 2_pl -> 76 -> L
    feats[1094] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_pl', '76') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '76', node2.lemma)
    
    # L -> 2_pl -> 128 -> L
    feats[1095] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_pl', '128') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '128', node2.lemma)
    
    # L -> 2_pl -> -11 -> L
    feats[1096] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_pl', '-11') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-11', node2.lemma)
    
    # L -> 2_pl -> -71 -> L
    feats[1097] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_pl', '-71') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-71', node2.lemma)
    
    # L -> 2_pl -> 9_sp -> L
    feats[1098] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_pl', '9_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_sp', node2.lemma)
    
    # L -> 2_pl -> -161 -> L
    feats[1099] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_pl', '-161') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-161', node2.lemma)
    
    # L -> 2_pl -> -157 -> L
    feats[1100] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_pl', '-157') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-157', node2.lemma)
    
    # L -> 2_pl -> -297 -> L
    feats[1101] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_pl', '-297') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-297', node2.lemma)
    
    # L -> 2_pl -> -126 -> L
    feats[1102] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_pl', '-126') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-126', node2.lemma)
    
    # L -> 2_pl -> 29_sg -> L
    feats[1103] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_pl', '29_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '29_sg', node2.lemma)
    
    # L -> 2_pl -> -103 -> L
    feats[1104] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_pl', '-103') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-103', node2.lemma)
    
    # L -> 2_pl -> -132 -> L
    feats[1105] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_pl', '-132') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-132', node2.lemma)
    
    # L -> 2_pl -> -14 -> L
    feats[1106] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_pl', '-14') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-14', node2.lemma)
    
    # L -> 2_pl -> 121 -> L
    feats[1107] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_pl', '121') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '121', node2.lemma)
    
    # L -> 2_pl -> fp -> L
    feats[1108] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_pl', 'fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'fp', node2.lemma)
    
    # L -> -121 -> 159 -> L
    feats[1109] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-121') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-121', '159') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '159', node2.lemma)
    
    # L -> -121 -> 2_pl -> L
    feats[1110] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-121') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-121', '2_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_pl', node2.lemma)
    
    # L -> -121 -> 7_sg -> L
    feats[1111] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-121') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-121', '7_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_sg', node2.lemma)
    
    # L -> -121 -> 6_sg -> L
    feats[1112] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-121') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-121', '6_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_sg', node2.lemma)
    
    # L -> instr. sg. -> 12_sg -> L
    feats[1113] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. sg.', '12_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_sg', node2.lemma)
    
    # L -> instr. sg. -> 9_sg -> L
    feats[1114] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. sg.', '9_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_sg', node2.lemma)
    
    # L -> instr. sg. -> 4_fp -> L
    feats[1115] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. sg.', '4_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_fp', node2.lemma)
    
    # L -> instr. sg. -> 108 -> L
    feats[1116] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. sg.', '108') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '108', node2.lemma)
    
    # L -> instr. sg. -> 9_pl -> L
    feats[1117] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. sg.', '9_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_pl', node2.lemma)
    
    # L -> instr. sg. -> 91 -> L
    feats[1118] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. sg.', '91') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
    
    # L -> instr. sg. -> -37 -> L
    feats[1119] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. sg.', '-37') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-37', node2.lemma)
    
    # L -> instr. sg. -> 131 -> L
    feats[1120] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. sg.', '131') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '131', node2.lemma)
    
    # L -> instr. sg. -> 101 -> L
    feats[1121] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. sg.', '101') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
    
    # L -> instr. sg. -> -240 -> L
    feats[1122] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. sg.', '-240') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-240', node2.lemma)
    
    # L -> instr. sg. -> -90 -> L
    feats[1123] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. sg.', '-90') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-90', node2.lemma)
    
    # L -> instr. sg. -> 10_pl -> L
    feats[1124] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. sg.', '10_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_pl', node2.lemma)
    
    # L -> instr. sg. -> 55 -> L
    feats[1125] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. sg.', '55') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '55', node2.lemma)
    
    # L -> instr. sg. -> 12_sp -> L
    feats[1126] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. sg.', '12_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_sp', node2.lemma)
    
    # L -> instr. sg. -> -72 -> L
    feats[1127] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. sg.', '-72') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-72', node2.lemma)
    
    # L -> instr. sg. -> -141 -> L
    feats[1128] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. sg.', '-141') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-141', node2.lemma)
    
    # L -> instr. sg. -> 2_sp -> L
    feats[1129] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. sg.', '2_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_sp', node2.lemma)
    
    # L -> instr. sg. -> nom -> L
    feats[1130] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. sg.', 'nom') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom', node2.lemma)
    
    # L -> instr. sg. -> 176 -> L
    feats[1131] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. sg.', '176') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '176', node2.lemma)
    
    # L -> instr. sg. -> -129 -> L
    feats[1132] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. sg.', '-129') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-129', node2.lemma)
    
    # L -> instr. sg. -> 30_sp -> L
    feats[1133] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. sg.', '30_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_sp', node2.lemma)
    
    # L -> instr. sg. -> 100 -> L
    feats[1134] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. sg.', '100') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '100', node2.lemma)
    
    # L -> instr. sg. -> 8_pl -> L
    feats[1135] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. sg.', '8_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_pl', node2.lemma)
    
    # L -> instr. sg. -> -28 -> L
    feats[1136] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. sg.', '-28') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-28', node2.lemma)
    
    # L -> acc. adj. -> 11_sg -> L
    feats[1137] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. adj.', '11_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_sg', node2.lemma)
    
    # L -> acc. adj. -> 16_sg -> L
    feats[1138] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. adj.', '16_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_sg', node2.lemma)
    
    # L -> acc. adj. -> abl. pl. -> L
    feats[1139] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. adj.', 'abl. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. pl.', node2.lemma)
    
    # L -> acc. adj. -> fem -> L
    feats[1140] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. adj.', 'fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'fem', node2.lemma)
    
    # L -> acc. adj. -> voc. masc. -> L
    feats[1141] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. adj.', 'voc. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. masc.', node2.lemma)
    
    # L -> acc. adj. -> -296 -> L
    feats[1142] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. adj.', '-296') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-296', node2.lemma)
    
    # L -> acc. adj. -> -94 -> L
    feats[1143] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. adj.', '-94') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-94', node2.lemma)
    
    # L -> acc. adj. -> 11_sp -> L
    feats[1144] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. adj.', '11_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_sp', node2.lemma)
    
    # L -> acc. adj. -> du -> L
    feats[1145] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. adj.', 'du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du', node2.lemma)
    
    # L -> acc. adj. -> 159 -> L
    feats[1146] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. adj.', '159') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '159', node2.lemma)
    
    # L -> acc. adj. -> 110 -> L
    feats[1147] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. adj.', '110') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '110', node2.lemma)
    
    # L -> acc. adj. -> -10 -> L
    feats[1148] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. adj.', '-10') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-10', node2.lemma)
    
    # L -> acc. adj. -> du_fp -> L
    feats[1149] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. adj.', 'du_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_fp', node2.lemma)
    
    # L -> acc. adj. -> 12_fp -> L
    feats[1150] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. adj.', '12_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_fp', node2.lemma)
    
    # L -> acc. adj. -> -220 -> L
    feats[1151] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. adj.', '-220') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-220', node2.lemma)
    
    # L -> acc. adj. -> 137 -> L
    feats[1152] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. adj.', '137') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    
    # L -> acc. adj. -> -33 -> L
    feats[1153] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. adj.', '-33') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-33', node2.lemma)
    
    # L -> acc. adj. -> -169 -> L
    feats[1154] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. adj.', '-169') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-169', node2.lemma)
    
    # L -> acc. adj. -> -123 -> L
    feats[1155] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. adj.', '-123') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-123', node2.lemma)
    
    # L -> acc. adj. -> loc. sg. -> L
    feats[1156] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. adj.', 'loc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc. sg.', node2.lemma)
    
    # L -> acc. adj. -> 56 -> L
    feats[1157] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. adj.', '56') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '56', node2.lemma)
    
    # L -> acc. adj. -> 41 -> L
    feats[1158] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. adj.', '41') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '41', node2.lemma)
    
    # L -> acc. adj. -> 5_sp -> L
    feats[1159] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. adj.', '5_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sp', node2.lemma)
    
    # L -> acc. adj. -> 54 -> L
    feats[1160] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. adj.', '54') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    
    # L -> acc. adj. -> -17 -> L
    feats[1161] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. adj.', '-17') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-17', node2.lemma)
    
    # L -> acc. adj. -> -157 -> L
    feats[1162] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. adj.', '-157') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-157', node2.lemma)
    
    # L -> acc. adj. -> voc. fem -> L
    feats[1163] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. adj.', 'voc. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. fem', node2.lemma)
    
    # L -> acc. adj. -> 138 -> L
    feats[1164] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. adj.', '138') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '138', node2.lemma)
    
    # L -> acc. adj. -> -291 -> L
    feats[1165] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. adj.', '-291') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-291', node2.lemma)
    
    # L -> acc. adj. -> -149 -> L
    feats[1166] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. adj.', '-149') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-149', node2.lemma)
    
    # L -> acc. adj. -> -271 -> L
    feats[1167] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. adj.', '-271') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-271', node2.lemma)
    
    # L -> acc. adj. -> 10_sg -> L
    feats[1168] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. adj.', '10_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_sg', node2.lemma)
    
    # L -> acc. adj. -> -245 -> L
    feats[1169] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. adj.', '-245') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-245', node2.lemma)
    
    # L -> acc. adj. -> -25 -> L
    feats[1170] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. adj.', '-25') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-25', node2.lemma)
    
    # L -> acc. adj. -> -93 -> L
    feats[1171] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. adj.', '-93') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-93', node2.lemma)
    
    # L -> acc. adj. -> 177 -> L
    feats[1172] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. adj.', '177') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '177', node2.lemma)
    
    # L -> acc. adj. -> 2_fp -> L
    feats[1173] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. adj.', '2_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_fp', node2.lemma)
    
    # L -> acc. adj. -> 5_pl -> L
    feats[1174] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. adj.', '5_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_pl', node2.lemma)
    
    # L -> acc. adj. -> adj -> L
    feats[1175] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. adj.', 'adj') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'adj', node2.lemma)
    
    # L -> acc. adj. -> 12_tp -> L
    feats[1176] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. adj.', '12_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_tp', node2.lemma)
    
    # L -> 30 -> -169 -> L
    feats[1177] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30', '-169') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-169', node2.lemma)
    
    # L -> 30 -> 12_tp -> L
    feats[1178] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30', '12_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_tp', node2.lemma)
    
    # L -> 30 -> -82 -> L
    feats[1179] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30', '-82') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-82', node2.lemma)
    
    # L -> 30 -> acc. masc. -> L
    feats[1180] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30', 'acc. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. masc.', node2.lemma)
    
    # L -> 30 -> 78 -> L
    feats[1181] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30', '78') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '78', node2.lemma)
    
    # L -> 30 -> 112 -> L
    feats[1182] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30', '112') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '112', node2.lemma)
    
    # L -> -242 -> 6_du -> L
    feats[1183] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-242') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-242', '6_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_du', node2.lemma)
    
    # L -> -242 -> 7_sp -> L
    feats[1184] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-242') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-242', '7_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_sp', node2.lemma)
    
    # L -> -242 -> 2_tp -> L
    feats[1185] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-242') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-242', '2_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_tp', node2.lemma)
    
    # L -> -242 -> 122 -> L
    feats[1186] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-242') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-242', '122') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '122', node2.lemma)
    
    # L -> -242 -> 100 -> L
    feats[1187] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-242') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-242', '100') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '100', node2.lemma)
    
    # L -> -242 -> nom. masc. -> L
    feats[1188] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-242') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-242', 'nom. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. masc.', node2.lemma)
    
    # L -> -113 -> -150 -> L
    feats[1189] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-113') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-113', '-150') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-150', node2.lemma)
    
    # L -> 4_sp -> 3_pl -> L
    feats[1190] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_sp', '3_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_pl', node2.lemma)
    
    # L -> 117 -> -114 -> L
    feats[1191] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '117') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '117', '-114') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-114', node2.lemma)
    
    # L -> 117 -> -59 -> L
    feats[1192] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '117') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '117', '-59') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-59', node2.lemma)
    
    # L -> 117 -> -273 -> L
    feats[1193] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '117') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '117', '-273') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-273', node2.lemma)
    
    # L -> 117 -> -93 -> L
    feats[1194] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '117') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '117', '-93') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-93', node2.lemma)
    
    # L -> 117 -> 12_tp -> L
    feats[1195] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '117') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '117', '12_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_tp', node2.lemma)
    
    # L -> -158 -> -38 -> L
    feats[1196] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-158') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-158', '-38') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-38', node2.lemma)
    
    # L -> -158 -> -308 -> L
    feats[1197] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-158') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-158', '-308') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-308', node2.lemma)
    
    # L -> dat. du. -> voc. du. -> L
    feats[1198] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. du.', 'voc. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. du.', node2.lemma)
    
    # L -> dat. du. -> 76 -> L
    feats[1199] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. du.', '76') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '76', node2.lemma)
    
    # L -> dat. du. -> -33 -> L
    feats[1200] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. du.', '-33') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-33', node2.lemma)
    
    # L -> dat. du. -> -71 -> L
    feats[1201] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. du.', '-71') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-71', node2.lemma)
    
    # L -> -143 -> 51 -> L
    feats[1202] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-143') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-143', '51') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '51', node2.lemma)
    
    # L -> -143 -> 109 -> L
    feats[1203] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-143') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-143', '109') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '109', node2.lemma)
    
    # L -> -143 -> 16_tp -> L
    feats[1204] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-143') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-143', '16_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_tp', node2.lemma)
    
    # L -> -143 -> 6_du -> L
    feats[1205] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-143') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-143', '6_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_du', node2.lemma)
    
    # L -> -143 -> 71 -> L
    feats[1206] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-143') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-143', '71') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '71', node2.lemma)
    
    # L -> -143 -> 30_fp -> L
    feats[1207] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-143') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-143', '30_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_fp', node2.lemma)
    
    # L -> -143 -> -299 -> L
    feats[1208] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-143') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-143', '-299') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-299', node2.lemma)
    
    # L -> -143 -> pl_sp -> L
    feats[1209] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-143') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-143', 'pl_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_sp', node2.lemma)
    
    # L -> 142 -> nom. pl. -> L
    feats[1210] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '142') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '142', 'nom. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. pl.', node2.lemma)
    
    # L -> 9_fp -> -117 -> L
    feats[1211] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_fp', '-117') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-117', node2.lemma)
    
    # L -> 9_fp -> 9_pl -> L
    feats[1212] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_fp', '9_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_pl', node2.lemma)
    
    # L -> 9_fp -> 6_sg -> L
    feats[1213] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_fp', '6_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_sg', node2.lemma)
    
    # L -> 9_fp -> 30_du -> L
    feats[1214] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_fp', '30_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_du', node2.lemma)
    
    # L -> 9_fp -> 94 -> L
    feats[1215] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_fp', '94') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '94', node2.lemma)
    
    # L -> 9_fp -> -16 -> L
    feats[1216] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_fp', '-16') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-16', node2.lemma)
    
    # L -> 9_fp -> 74 -> L
    feats[1217] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_fp', '74') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '74', node2.lemma)
    
    # L -> 9_fp -> 27_du -> L
    feats[1218] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_fp', '27_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_du', node2.lemma)
    
    # L -> 9_fp -> 13_fp -> L
    feats[1219] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_fp', '13_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_fp', node2.lemma)
    
    # L -> 9_fp -> 100 -> L
    feats[1220] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_fp', '100') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '100', node2.lemma)
    
    # L -> 135 -> 28_sg -> L
    feats[1221] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '135') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', '28_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28_sg', node2.lemma)
    
    # L -> 135 -> 3_tp -> L
    feats[1222] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '135') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', '3_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_tp', node2.lemma)
    
    # L -> 135 -> -163 -> L
    feats[1223] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '135') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', '-163') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-163', node2.lemma)
    
    # L -> 135 -> 12_fp -> L
    feats[1224] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '135') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', '12_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_fp', node2.lemma)
    
    # L -> 135 -> 137 -> L
    feats[1225] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '135') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', '137') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    
    # L -> 135 -> 76 -> L
    feats[1226] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '135') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', '76') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '76', node2.lemma)
    
    # L -> 135 -> -33 -> L
    feats[1227] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '135') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', '-33') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-33', node2.lemma)
    
    # L -> 135 -> pl_fp -> L
    feats[1228] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '135') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', 'pl_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_fp', node2.lemma)
    
    # L -> 135 -> 14_sp -> L
    feats[1229] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '135') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', '14_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_sp', node2.lemma)
    
    # L -> 135 -> 77 -> L
    feats[1230] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '135') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', '77') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '77', node2.lemma)
    
    # L -> 135 -> 6_pl -> L
    feats[1231] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '135') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', '6_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_pl', node2.lemma)
    
    # L -> 135 -> 7_tp -> L
    feats[1232] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '135') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', '7_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_tp', node2.lemma)
    
    # L -> 135 -> -104 -> L
    feats[1233] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '135') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', '-104') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-104', node2.lemma)
    
    # L -> 135 -> -50 -> L
    feats[1234] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '135') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', '-50') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-50', node2.lemma)
    
    # L -> 135 -> 89 -> L
    feats[1235] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '135') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', '89') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '89', node2.lemma)
    
    # L -> 3_fp -> 112 -> L
    feats[1236] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_fp', '112') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '112', node2.lemma)
    
    # L -> 3_fp -> 16_du -> L
    feats[1237] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_fp', '16_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_du', node2.lemma)
    
    # L -> 3_fp -> 154 -> L
    feats[1238] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_fp', '154') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '154', node2.lemma)
    
    # L -> 3_fp -> 162 -> L
    feats[1239] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_fp', '162') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '162', node2.lemma)
    
    # L -> 3_fp -> 172 -> L
    feats[1240] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_fp', '172') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '172', node2.lemma)
    
    # L -> 3_fp -> -30 -> L
    feats[1241] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_fp', '-30') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-30', node2.lemma)
    
    # L -> 3_fp -> 1 -> L
    feats[1242] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_fp', '1') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '1', node2.lemma)
    
    # L -> 3_fp -> -303 -> L
    feats[1243] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_fp', '-303') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-303', node2.lemma)
    
    # L -> gen -> -147 -> L
    feats[1244] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen', '-147') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-147', node2.lemma)
    
    # L -> gen -> du_tp -> L
    feats[1245] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen', 'du_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_tp', node2.lemma)
    
    # L -> 14_fp -> 2_fp -> L
    feats[1246] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_fp', '2_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_fp', node2.lemma)
    
    # L -> 14_fp -> 11_pl -> L
    feats[1247] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_fp', '11_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_pl', node2.lemma)
    
    # L -> 14_fp -> acc. fem -> L
    feats[1248] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_fp', 'acc. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. fem', node2.lemma)
    
    # L -> 14_fp -> -302 -> L
    feats[1249] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_fp', '-302') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-302', node2.lemma)
    
    # L -> 14_fp -> 136 -> L
    feats[1250] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_fp', '136') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '136', node2.lemma)
    
    # L -> 14_fp -> fp -> L
    feats[1251] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_fp', 'fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'fp', node2.lemma)
    
    # L -> 14_fp -> -21 -> L
    feats[1252] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_fp', '-21') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-21', node2.lemma)
    
    # L -> 14_fp -> 114 -> L
    feats[1253] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_fp', '114') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '114', node2.lemma)
    
    # L -> 14_fp -> 96 -> L
    feats[1254] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_fp', '96') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '96', node2.lemma)
    
    # L -> 14_fp -> 7_sg -> L
    feats[1255] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_fp', '7_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_sg', node2.lemma)
    
    # L -> 14_fp -> -50 -> L
    feats[1256] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_fp', '-50') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-50', node2.lemma)
    
    # L -> 14_fp -> -200 -> L
    feats[1257] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_fp', '-200') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-200', node2.lemma)
    
    # L -> 14_fp -> 89 -> L
    feats[1258] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_fp', '89') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '89', node2.lemma)
    
    # L -> 14_fp -> -243 -> L
    feats[1259] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_fp', '-243') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-243', node2.lemma)
    
    # L -> 14_fp -> 11_fp -> L
    feats[1260] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_fp', '11_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_fp', node2.lemma)
    
    # L -> -307 -> -149 -> L
    feats[1261] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-307') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-307', '-149') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-149', node2.lemma)
    
    # L -> -307 -> -271 -> L
    feats[1262] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-307') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-307', '-271') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-271', node2.lemma)
    
    # L -> -307 -> -139 -> L
    feats[1263] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-307') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-307', '-139') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-139', node2.lemma)
    
    # L -> -307 -> instr. sg. -> L
    feats[1264] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-307') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-307', 'instr. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. sg.', node2.lemma)
    
    # L -> -307 -> 120 -> L
    feats[1265] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-307') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-307', '120') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '120', node2.lemma)
    
    # L -> -307 -> 28_tp -> L
    feats[1266] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-307') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-307', '28_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28_tp', node2.lemma)
    
    # L -> -307 -> -23 -> L
    feats[1267] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-307') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-307', '-23') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-23', node2.lemma)
    
    # L -> -307 -> 49 -> L
    feats[1268] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-307') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-307', '49') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '49', node2.lemma)
    
    # L -> -307 -> 142 -> L
    feats[1269] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-307') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-307', '142') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '142', node2.lemma)
    
    # L -> -307 -> 9_fp -> L
    feats[1270] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-307') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-307', '9_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_fp', node2.lemma)
    
    # L -> -307 -> 135 -> L
    feats[1271] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-307') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-307', '135') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '135', node2.lemma)
    
    # L -> -307 -> abl. sg. -> L
    feats[1272] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-307') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-307', 'abl. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. sg.', node2.lemma)
    
    # L -> -307 -> acc. masc. -> L
    feats[1273] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-307') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-307', 'acc. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. masc.', node2.lemma)
    
    # L -> -307 -> sg_sp -> L
    feats[1274] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-307') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-307', 'sg_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_sp', node2.lemma)
    
    # L -> -307 -> -47 -> L
    feats[1275] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-307') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-307', '-47') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-47', node2.lemma)
    
    # L -> -20 -> acc -> L
    feats[1276] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-20') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-20', 'acc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc', node2.lemma)
    
    # L -> -20 -> -137 -> L
    feats[1277] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-20') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-20', '-137') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-137', node2.lemma)
    
    # L -> -20 -> 93 -> L
    feats[1278] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-20') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-20', '93') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '93', node2.lemma)
    
    # L -> -20 -> 8_sg -> L
    feats[1279] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-20') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-20', '8_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_sg', node2.lemma)
    
    # L -> -20 -> 9_pl -> L
    feats[1280] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-20') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-20', '9_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_pl', node2.lemma)
    
    # L -> -20 -> -141 -> L
    feats[1281] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-20') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-20', '-141') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-141', node2.lemma)
    
    # L -> 8_tp -> -122 -> L
    feats[1282] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_tp', '-122') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-122', node2.lemma)
    
    # L -> 8_tp -> -48 -> L
    feats[1283] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_tp', '-48') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-48', node2.lemma)
    
    # L -> 8_tp -> -245 -> L
    feats[1284] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_tp', '-245') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-245', node2.lemma)
    
    # L -> 8_tp -> instr. du. -> L
    feats[1285] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_tp', 'instr. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. du.', node2.lemma)
    
    # L -> 8_tp -> instr. pl. -> L
    feats[1286] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_tp', 'instr. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. pl.', node2.lemma)
    
    # L -> 8_tp -> sg_fp -> L
    feats[1287] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_tp', 'sg_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_fp', node2.lemma)
    
    # L -> 8_tp -> nom. sg. -> L
    feats[1288] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_tp', 'nom. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. sg.', node2.lemma)
    
    # L -> 8_tp -> 182 -> L
    feats[1289] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_tp', '182') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '182', node2.lemma)
    
    # L -> 8_tp -> acc. fem -> L
    feats[1290] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_tp', 'acc. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. fem', node2.lemma)
    
    # L -> du_sp -> -101 -> L
    feats[1291] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_sp', '-101') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-101', node2.lemma)
    
    # L -> du_sp -> 173 -> L
    feats[1292] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_sp', '173') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '173', node2.lemma)
    
    # L -> du_sp -> 112 -> L
    feats[1293] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_sp', '112') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '112', node2.lemma)
    
    # L -> 116 -> -147 -> L
    feats[1294] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '116') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '116', '-147') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-147', node2.lemma)
    
    # L -> 116 -> -151 -> L
    feats[1295] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '116') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '116', '-151') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-151', node2.lemma)
    
    # L -> 116 -> -77 -> L
    feats[1296] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '116') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '116', '-77') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-77', node2.lemma)
    
    # L -> 81 -> 5_sg -> L
    feats[1297] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '81', '5_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sg', node2.lemma)
    
    # L -> 81 -> nom. pl. -> L
    feats[1298] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '81', 'nom. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. pl.', node2.lemma)
    
    # L -> 81 -> 61 -> L
    feats[1299] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '81', '61') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '61', node2.lemma)
    
    # L -> 81 -> -57 -> L
    feats[1300] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '81', '-57') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-57', node2.lemma)
    
    # L -> 81 -> -35 -> L
    feats[1301] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '81', '-35') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    
    # L -> 81 -> 122 -> L
    feats[1302] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '81', '122') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '122', node2.lemma)
    
    # L -> 81 -> 3_sp -> L
    feats[1303] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '81', '3_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_sp', node2.lemma)
    
    # L -> 81 -> 74 -> L
    feats[1304] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '81', '74') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '74', node2.lemma)
    
    # L -> 95 -> 11_sg -> L
    feats[1305] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '95') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '95', '11_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_sg', node2.lemma)
    
    # L -> 95 -> 16_sg -> L
    feats[1306] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '95') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '95', '16_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_sg', node2.lemma)
    
    # L -> 95 -> -159 -> L
    feats[1307] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '95') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '95', '-159') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-159', node2.lemma)
    
    # L -> 95 -> 2_sg -> L
    feats[1308] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '95') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '95', '2_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_sg', node2.lemma)
    
    # L -> 95 -> -94 -> L
    feats[1309] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '95') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '95', '-94') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-94', node2.lemma)
    
    # L -> -109 -> -161 -> L
    feats[1310] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-109') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-109', '-161') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-161', node2.lemma)
    
    # L -> 12_sg -> acc. neutr. -> L
    feats[1311] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', 'acc. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. neutr.', node2.lemma)
    
    # L -> 12_sg -> tp -> L
    feats[1312] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', 'tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'tp', node2.lemma)
    
    # L -> 12_sg -> -271 -> L
    feats[1313] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', '-271') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-271', node2.lemma)
    
    # L -> 12_sg -> 10_sg -> L
    feats[1314] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', '10_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_sg', node2.lemma)
    
    # L -> 12_sg -> -22 -> L
    feats[1315] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', '-22') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-22', node2.lemma)
    
    # L -> 12_sg -> -77 -> L
    feats[1316] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', '-77') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-77', node2.lemma)
    
    # L -> 12_sg -> -93 -> L
    feats[1317] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', '-93') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-93', node2.lemma)
    
    # L -> 12_sg -> adj -> L
    feats[1318] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', 'adj') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'adj', node2.lemma)
    
    # L -> 12_sg -> -230 -> L
    feats[1319] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', '-230') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-230', node2.lemma)
    
    # L -> 12_sg -> -126 -> L
    feats[1320] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', '-126') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-126', node2.lemma)
    
    # L -> 12_sg -> -34 -> L
    feats[1321] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', '-34') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-34', node2.lemma)
    
    # L -> 12_sg -> 29_sg -> L
    feats[1322] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', '29_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '29_sg', node2.lemma)
    
    # L -> 12_sg -> -103 -> L
    feats[1323] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', '-103') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-103', node2.lemma)
    
    # L -> 12_sg -> -58 -> L
    feats[1324] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', '-58') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-58', node2.lemma)
    
    # L -> 12_sg -> nom. sg. -> L
    feats[1325] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', 'nom. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. sg.', node2.lemma)
    
    # L -> 12_sg -> 90 -> L
    feats[1326] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', '90') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '90', node2.lemma)
    
    # L -> 12_sg -> 2_pl -> L
    feats[1327] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', '2_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_pl', node2.lemma)
    
    # L -> 12_sg -> instr. sg. -> L
    feats[1328] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', 'instr. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. sg.', node2.lemma)
    
    # L -> 12_sg -> acc. adj. -> L
    feats[1329] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', 'acc. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. adj.', node2.lemma)
    
    # L -> 12_sg -> 11_pl -> L
    feats[1330] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', '11_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_pl', node2.lemma)
    
    # L -> 12_sg -> 156 -> L
    feats[1331] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', '156') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '156', node2.lemma)
    
    # L -> 12_sg -> acc. fem -> L
    feats[1332] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', 'acc. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. fem', node2.lemma)
    
    # L -> 12_sg -> -302 -> L
    feats[1333] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', '-302') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-302', node2.lemma)
    
    # L -> 12_sg -> 136 -> L
    feats[1334] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', '136') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '136', node2.lemma)
    
    # L -> 12_sg -> fp -> L
    feats[1335] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', 'fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'fp', node2.lemma)
    
    # L -> 12_sg -> 2_du -> L
    feats[1336] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', '2_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_du', node2.lemma)
    
    # L -> 12_sg -> 92 -> L
    feats[1337] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', '92') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '92', node2.lemma)
    
    # L -> 12_sg -> -73 -> L
    feats[1338] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', '-73') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-73', node2.lemma)
    
    # L -> 12_sg -> -279 -> L
    feats[1339] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', '-279') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-279', node2.lemma)
    
    # L -> 12_sg -> -44 -> L
    feats[1340] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', '-44') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-44', node2.lemma)
    
    # L -> 12_sg -> 3_du -> L
    feats[1341] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', '3_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
    
    # L -> 12_sg -> 14_tp -> L
    feats[1342] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', '14_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_tp', node2.lemma)
    
    # L -> 12_sg -> 9_du -> L
    feats[1343] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', '9_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_du', node2.lemma)
    
    # L -> 12_sg -> 114 -> L
    feats[1344] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', '114') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '114', node2.lemma)
    
    # L -> 12_sg -> 7_du -> L
    feats[1345] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', '7_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_du', node2.lemma)
    
    # L -> 12_sg -> -49 -> L
    feats[1346] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', '-49') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-49', node2.lemma)
    
    # L -> 12_sg -> -144 -> L
    feats[1347] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', '-144') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-144', node2.lemma)
    
    # L -> 12_sg -> -154 -> L
    feats[1348] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', '-154') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-154', node2.lemma)
    
    # L -> 12_sg -> -242 -> L
    feats[1349] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', '-242') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-242', node2.lemma)
    
    # L -> 12_sg -> -113 -> L
    feats[1350] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', '-113') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-113', node2.lemma)
    
    # L -> 12_sg -> 4_sp -> L
    feats[1351] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', '4_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_sp', node2.lemma)
    
    # L -> 12_sg -> 3_fp -> L
    feats[1352] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', '3_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_fp', node2.lemma)
    
    # L -> 12_sg -> 14_fp -> L
    feats[1353] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', '14_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_fp', node2.lemma)
    
    # L -> 12_sg -> 3 -> L
    feats[1354] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', '3') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3', node2.lemma)
    
    # L -> 12_sg -> 116 -> L
    feats[1355] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', '116') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '116', node2.lemma)
    
    # L -> 12_sg -> 81 -> L
    feats[1356] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', '81') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '81', node2.lemma)
    
    # L -> 12_sg -> 95 -> L
    feats[1357] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', '95') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '95', node2.lemma)
    
    # L -> 12_sg -> -109 -> L
    feats[1358] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', '-109') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-109', node2.lemma)
    
    # L -> 12_sg -> 7_sg -> L
    feats[1359] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', '7_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_sg', node2.lemma)
    
    # L -> 12_sg -> -101 -> L
    feats[1360] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', '-101') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-101', node2.lemma)
    
    # L -> 12_sg -> -99 -> L
    feats[1361] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', '-99') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-99', node2.lemma)
    
    # L -> 12_sg -> 6_sp -> L
    feats[1362] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', '6_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_sp', node2.lemma)
    
    # L -> 12_sg -> -63 -> L
    feats[1363] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', '-63') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-63', node2.lemma)
    
    # L -> 12_sg -> acc. sg. -> L
    feats[1364] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', 'acc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. sg.', node2.lemma)
    
    # L -> 12_sg -> 3_sp -> L
    feats[1365] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', '3_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_sp', node2.lemma)
    
    # L -> -42 -> 137 -> L
    feats[1366] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-42') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-42', '137') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
    
    # L -> -42 -> 34 -> L
    feats[1367] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-42') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-42', '34') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '34', node2.lemma)
    
    # L -> -42 -> -104 -> L
    feats[1368] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-42') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-42', '-104') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-104', node2.lemma)
    
    # L -> -42 -> -42 -> L
    feats[1369] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-42') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-42', '-42') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-42', node2.lemma)
    
    # L -> -42 -> loc. pl. -> L
    feats[1370] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-42') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-42', 'loc. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc. pl.', node2.lemma)
    
    # L -> -50 -> abl. du. -> L
    feats[1371] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-50') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-50', 'abl. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. du.', node2.lemma)
    
    # L -> -50 -> -150 -> L
    feats[1372] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-50') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-50', '-150') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-150', node2.lemma)
    
    # L -> -50 -> -147 -> L
    feats[1373] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-50') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-50', '-147') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-147', node2.lemma)
    
    # L -> -50 -> 39 -> L
    feats[1374] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-50') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-50', '39') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '39', node2.lemma)
    
    # L -> -50 -> -86 -> L
    feats[1375] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-50') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-50', '-86') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-86', node2.lemma)
    
    # L -> -50 -> du -> L
    feats[1376] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-50') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-50', 'du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du', node2.lemma)
    
    # L -> -50 -> 159 -> L
    feats[1377] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-50') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-50', '159') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '159', node2.lemma)
    
    # L -> -50 -> -241 -> L
    feats[1378] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-50') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-50', '-241') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-241', node2.lemma)
    
    # L -> -50 -> 12_fp -> L
    feats[1379] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-50') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-50', '12_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_fp', node2.lemma)
    
    # L -> -50 -> 128 -> L
    feats[1380] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-50') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-50', '128') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '128', node2.lemma)
    
    # L -> -50 -> loc. sg. -> L
    feats[1381] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-50') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-50', 'loc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc. sg.', node2.lemma)
    
    # L -> -50 -> 54 -> L
    feats[1382] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-50') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-50', '54') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    
    # L -> -50 -> 179 -> L
    feats[1383] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-50') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-50', '179') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '179', node2.lemma)
    
    # L -> -50 -> -271 -> L
    feats[1384] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-50') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-50', '-271') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-271', node2.lemma)
    
    # L -> -50 -> -25 -> L
    feats[1385] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-50') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-50', '-25') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-25', node2.lemma)
    
    # L -> -50 -> -139 -> L
    feats[1386] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-50') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-50', '-139') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-139', node2.lemma)
    
    # L -> -50 -> 151 -> L
    feats[1387] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-50') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-50', '151') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '151', node2.lemma)
    
    # L -> -50 -> 98 -> L
    feats[1388] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-50') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-50', '98') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '98', node2.lemma)
    
    # L -> -50 -> -121 -> L
    feats[1389] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-50') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-50', '-121') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-121', node2.lemma)
    
    # L -> -50 -> instr. sg. -> L
    feats[1390] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-50') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-50', 'instr. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. sg.', node2.lemma)
    
    # L -> -50 -> acc. adj. -> L
    feats[1391] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-50') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-50', 'acc. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. adj.', node2.lemma)
    
    # L -> -50 -> -82 -> L
    feats[1392] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-50') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-50', '-82') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-82', node2.lemma)
    
    # L -> -50 -> 37 -> L
    feats[1393] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-50') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-50', '37') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '37', node2.lemma)
    
    # L -> -50 -> -42 -> L
    feats[1394] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-50') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-50', '-42') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-42', node2.lemma)
    
    # L -> -50 -> 80 -> L
    feats[1395] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-50') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-50', '80') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '80', node2.lemma)
    
    # L -> -50 -> 8_sp -> L
    feats[1396] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-50') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-50', '8_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_sp', node2.lemma)
    
    # L -> -50 -> 38 -> L
    feats[1397] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-50') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-50', '38') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '38', node2.lemma)
    
    # L -> -50 -> -137 -> L
    feats[1398] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-50') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-50', '-137') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-137', node2.lemma)
    
    # L -> -51 -> 7_pl -> L
    feats[1399] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-51') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-51', '7_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_pl', node2.lemma)
    
    # L -> -51 -> 14_sg -> L
    feats[1400] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-51') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-51', '14_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_sg', node2.lemma)
    
    # L -> -51 -> -12 -> L
    feats[1401] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-51') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-51', '-12') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-12', node2.lemma)
    
    # L -> -51 -> 100 -> L
    feats[1402] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-51') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-51', '100') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '100', node2.lemma)
    
    # L -> -101 -> du_tp -> L
    feats[1403] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-101') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-101', 'du_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_tp', node2.lemma)
    
    # L -> -101 -> 54 -> L
    feats[1404] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-101') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-101', '54') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
    
    # L -> -101 -> -271 -> L
    feats[1405] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-101') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-101', '-271') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-271', node2.lemma)
    
    # L -> -101 -> instr. du. -> L
    feats[1406] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-101') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-101', 'instr. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. du.', node2.lemma)
    
    # L -> -101 -> -27 -> L
    feats[1407] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-101') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-101', '-27') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-27', node2.lemma)
    
    # L -> -101 -> -67 -> L
    feats[1408] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-101') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-101', '-67') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-67', node2.lemma)
    
    # L -> -293 -> -47 -> L
    feats[1409] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '-47') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-47', node2.lemma)
    
    # L -> -293 -> neutr -> L
    feats[1410] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', 'neutr') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'neutr', node2.lemma)
    
    # L -> -293 -> 5_fp -> L
    feats[1411] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '5_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_fp', node2.lemma)
    
    # L -> -293 -> 139 -> L
    feats[1412] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '139') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '139', node2.lemma)
    
    # L -> -293 -> 6_tp -> L
    feats[1413] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '6_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_tp', node2.lemma)
    
    # L -> -293 -> 73 -> L
    feats[1414] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '73') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '73', node2.lemma)
    
    # L -> -293 -> 97 -> L
    feats[1415] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '97') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '97', node2.lemma)
    
    # L -> -293 -> -62 -> L
    feats[1416] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '-62') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-62', node2.lemma)
    
    # L -> -293 -> 99 -> L
    feats[1417] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '99') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '99', node2.lemma)
    
    # L -> -293 -> 36 -> L
    feats[1418] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '36') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '36', node2.lemma)
    
    # L -> -293 -> -269 -> L
    feats[1419] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '-269') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-269', node2.lemma)
    
    # L -> -293 -> -283 -> L
    feats[1420] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '-283') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-283', node2.lemma)
    
    # L -> -293 -> 8_fp -> L
    feats[1421] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '8_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
    
    # L -> -293 -> 16_fp -> L
    feats[1422] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '16_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_fp', node2.lemma)
    
    # L -> -293 -> -133 -> L
    feats[1423] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '-133') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-133', node2.lemma)
    
    # L -> -293 -> 129 -> L
    feats[1424] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '129') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '129', node2.lemma)
    
    # L -> -293 -> 78 -> L
    feats[1425] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '78') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '78', node2.lemma)
    
    # L -> -293 -> 169 -> L
    feats[1426] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '169') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '169', node2.lemma)
    
    # L -> -293 -> 109 -> L
    feats[1427] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '109') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '109', node2.lemma)
    
    # L -> -293 -> 3_sg -> L
    feats[1428] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '3_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_sg', node2.lemma)
    
    # L -> -293 -> gen. sg. -> L
    feats[1429] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', 'gen. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen. sg.', node2.lemma)
    
    # L -> -293 -> 8_sp -> L
    feats[1430] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '8_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_sp', node2.lemma)
    
    # L -> -293 -> 174 -> L
    feats[1431] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '174') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '174', node2.lemma)
    
    # L -> -293 -> 60 -> L
    feats[1432] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '60') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '60', node2.lemma)
    
    # L -> -293 -> 5_sg -> L
    feats[1433] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '5_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sg', node2.lemma)
    
    # L -> -293 -> voc -> L
    feats[1434] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', 'voc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc', node2.lemma)
    
    # L -> -293 -> -210 -> L
    feats[1435] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '-210') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-210', node2.lemma)
    
    # L -> -293 -> 68 -> L
    feats[1436] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '68') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '68', node2.lemma)
    
    # L -> -293 -> acc. du. -> L
    feats[1437] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', 'acc. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. du.', node2.lemma)
    
    # L -> -293 -> 30_tp -> L
    feats[1438] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '30_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_tp', node2.lemma)
    
    # L -> -293 -> -249 -> L
    feats[1439] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '-249') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-249', node2.lemma)
    
    # L -> -293 -> 13_pl -> L
    feats[1440] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '13_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_pl', node2.lemma)
    
    # L -> -293 -> 93 -> L
    feats[1441] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '93') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '93', node2.lemma)
    
    # L -> -293 -> 6_du -> L
    feats[1442] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '6_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_du', node2.lemma)
    
    # L -> -293 -> -29 -> L
    feats[1443] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '-29') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-29', node2.lemma)
    
    # L -> -293 -> -32 -> L
    feats[1444] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '-32') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-32', node2.lemma)
    
    # L -> -293 -> -54 -> L
    feats[1445] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '-54') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-54', node2.lemma)
    
    # L -> -293 -> 175 -> L
    feats[1446] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '175') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '175', node2.lemma)
    
    # L -> -293 -> masc -> L
    feats[1447] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', 'masc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'masc', node2.lemma)
    
    # L -> -293 -> 29_tp -> L
    feats[1448] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '29_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '29_tp', node2.lemma)
    
    # L -> -293 -> -299 -> L
    feats[1449] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '-299') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-299', node2.lemma)
    
    # L -> -293 -> 4_sg -> L
    feats[1450] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '4_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_sg', node2.lemma)
    
    # L -> -293 -> 5_du -> L
    feats[1451] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '5_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_du', node2.lemma)
    
    # L -> -293 -> 27_pl -> L
    feats[1452] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '27_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_pl', node2.lemma)
    
    # L -> -293 -> 31 -> L
    feats[1453] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '31') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '31', node2.lemma)
    
    # L -> -293 -> 8_sg -> L
    feats[1454] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '8_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_sg', node2.lemma)
    
    # L -> -293 -> instr. adj. -> L
    feats[1455] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', 'instr. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. adj.', node2.lemma)
    
    # L -> -293 -> -46 -> L
    feats[1456] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '-46') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-46', node2.lemma)
    
    # L -> -293 -> 108 -> L
    feats[1457] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '108') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '108', node2.lemma)
    
    # L -> -293 -> -56 -> L
    feats[1458] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '-56') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-56', node2.lemma)
    
    # L -> -293 -> -35 -> L
    feats[1459] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '-35') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
    
    # L -> -293 -> -117 -> L
    feats[1460] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', '-117') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-117', node2.lemma)
    
    # L -> 30_pl -> -36 -> L
    feats[1461] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_pl', '-36') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-36', node2.lemma)
    
    # L -> 30_pl -> -101 -> L
    feats[1462] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_pl', '-101') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-101', node2.lemma)
    
    # L -> 30_pl -> 15_sp -> L
    feats[1463] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_pl', '15_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_sp', node2.lemma)
    
    # L -> 30_pl -> -301 -> L
    feats[1464] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_pl', '-301') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-301', node2.lemma)
    
    # L -> 30_pl -> 129 -> L
    feats[1465] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_pl', '129') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '129', node2.lemma)
    
    # L -> 30_pl -> 173 -> L
    feats[1466] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_pl', '173') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '173', node2.lemma)
    
    # L -> 30_pl -> 78 -> L
    feats[1467] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_pl', '78') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '78', node2.lemma)
    
    # L -> 30_pl -> 112 -> L
    feats[1468] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_pl', '112') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '112', node2.lemma)
    
    # L -> 30_pl -> 28 -> L
    feats[1469] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_pl', '28') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
    
    # L -> 30_pl -> voc -> L
    feats[1470] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_pl', 'voc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc', node2.lemma)
    
    # L -> 30_pl -> 130 -> L
    feats[1471] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_pl', '130') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '130', node2.lemma)
    
    # L -> -63 -> -296 -> L
    feats[1472] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-63') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-63', '-296') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-296', node2.lemma)
    
    # L -> -63 -> nom. fem -> L
    feats[1473] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-63') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-63', 'nom. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. fem', node2.lemma)
    
    # L -> -89 -> -90 -> L
    feats[1474] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-89') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-89', '-90') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-90', node2.lemma)
    
    # L -> sg_sp -> sg_tp -> L
    feats[1475] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg_sp', 'sg_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_tp', node2.lemma)
    
    # L -> sg_sp -> dat -> L
    feats[1476] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg_sp', 'dat') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat', node2.lemma)
    
    # L -> sg_sp -> 118 -> L
    feats[1477] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg_sp', '118') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '118', node2.lemma)
    
    # L -> sg_sp -> gen -> L
    feats[1478] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg_sp', 'gen') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen', node2.lemma)
    
    # L -> sg_sp -> -99 -> L
    feats[1479] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg_sp', '-99') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-99', node2.lemma)
    
    # L -> sg_sp -> acc. sg. -> L
    feats[1480] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg_sp', 'acc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. sg.', node2.lemma)
    
    # L -> 5_tp -> -271 -> L
    feats[1481] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_tp', '-271') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-271', node2.lemma)
    
    # L -> 5_tp -> 10_sg -> L
    feats[1482] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_tp', '10_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_sg', node2.lemma)
    
    # L -> 5_tp -> 7_tp -> L
    feats[1483] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_tp', '7_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_tp', node2.lemma)
    
    # L -> 5_tp -> sg_fp -> L
    feats[1484] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_tp', 'sg_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_fp', node2.lemma)
    
    # L -> 5_tp -> -121 -> L
    feats[1485] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_tp', '-121') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-121', node2.lemma)
    
    # L -> 5_tp -> 30_sg -> L
    feats[1486] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_tp', '30_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_sg', node2.lemma)
    
    # L -> 5_tp -> 13_tp -> L
    feats[1487] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_tp', '13_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_tp', node2.lemma)
    
    # L -> 5_tp -> 49 -> L
    feats[1488] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_tp', '49') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '49', node2.lemma)
    
    # L -> 5_tp -> 114 -> L
    feats[1489] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_tp', '114') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '114', node2.lemma)
    
    # L -> 5_tp -> 96 -> L
    feats[1490] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_tp', '96') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '96', node2.lemma)
    
    # L -> 5_tp -> 30 -> L
    feats[1491] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_tp', '30') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30', node2.lemma)
    
    # L -> 5_tp -> -242 -> L
    feats[1492] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_tp', '-242') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-242', node2.lemma)
    
    # L -> 5_tp -> -113 -> L
    feats[1493] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_tp', '-113') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-113', node2.lemma)
    
    # L -> 5_tp -> 81 -> L
    feats[1494] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_tp', '81') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '81', node2.lemma)
    
    # L -> 5_tp -> 95 -> L
    feats[1495] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_tp', '95') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '95', node2.lemma)
    
    # L -> 5_tp -> 7_sg -> L
    feats[1496] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_tp', '7_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_sg', node2.lemma)
    
    # L -> 27_tp -> masc -> L
    feats[1497] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_tp', 'masc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'masc', node2.lemma)
    
    # L -> 27_tp -> 158 -> L
    feats[1498] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_tp', '158') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '158', node2.lemma)
    
    # L -> 27_tp -> 27_du -> L
    feats[1499] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_tp', '27_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_du', node2.lemma)
    
    feats[feats < 1e-25] = 1e-25
    return -np.log10(feats)
    