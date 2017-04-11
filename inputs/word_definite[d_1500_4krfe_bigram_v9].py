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
    
    # L->121->L
    feats[0] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '121') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '121', node2.lemma)
            
    # L->-15->L
    feats[1] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-15') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-15', node2.lemma)
            
    # L->-121->L
    feats[2] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-121') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-121', node2.lemma)
            
    # L->117->L
    feats[3] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '117') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '117', node2.lemma)
            
    # L->10_sg->L
    feats[4] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_sg', node2.lemma)
            
    # L->-240->L
    feats[5] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-240') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-240', node2.lemma)
            
    # L->55->L
    feats[6] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '55') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '55', node2.lemma)
            
    # L->-245->L
    feats[7] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-245') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-245', node2.lemma)
            
    # L->171->L
    feats[8] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '171') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '171', node2.lemma)
            
    # L->154->L
    feats[9] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '154') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '154', node2.lemma)
            
    # L->14_sp->L
    feats[10] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_sp', node2.lemma)
            
    # L->38->L
    feats[11] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '38') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '38', node2.lemma)
            
    # L->-151->L
    feats[12] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-151') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-151', node2.lemma)
            
    # L->du_tp->L
    feats[13] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_tp', node2.lemma)
            
    # L->-299->L
    feats[14] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-299') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-299', node2.lemma)
            
    # L->61->L
    feats[15] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '61') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '61', node2.lemma)
            
    # L->-79->L
    feats[16] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-79') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-79', node2.lemma)
            
    # L->9_fp->L
    feats[17] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_fp', node2.lemma)
            
    # L->12_tp->L
    feats[18] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_tp', node2.lemma)
            
    # L->29_tp->L
    feats[19] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '29_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '29_tp', node2.lemma)
            
    # L->pl_sp->L
    feats[20] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_sp', node2.lemma)
            
    # L->135->L
    feats[21] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '135') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '135', node2.lemma)
            
    # L->-122->L
    feats[22] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-122') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-122', node2.lemma)
            
    # L->acc. adj.->L
    feats[23] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. adj.', node2.lemma)
            
    # L->102->L
    feats[24] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '102') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '102', node2.lemma)
            
    # L->-152->L
    feats[25] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-152') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-152', node2.lemma)
            
    # L->3_pl->L
    feats[26] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_pl', node2.lemma)
            
    # L->151->L
    feats[27] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '151') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '151', node2.lemma)
            
    # L->15_du->L
    feats[28] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_du', node2.lemma)
            
    # L->82->L
    feats[29] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '82') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '82', node2.lemma)
            
    # L->instr. fem->L
    feats[30] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. fem', node2.lemma)
            
    # L->-98->L
    feats[31] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-98') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-98', node2.lemma)
            
    # L->-242->L
    feats[32] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-242') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-242', node2.lemma)
            
    # L->27_du->L
    feats[33] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_du', node2.lemma)
            
    # L->92->L
    feats[34] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '92') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '92', node2.lemma)
            
    # L->gen. pl.->L
    feats[35] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen. pl.', node2.lemma)
            
    # L->74->L
    feats[36] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '74') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '74', node2.lemma)
            
    # L->-86->L
    feats[37] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-86') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-86', node2.lemma)
            
    # L->8_du->L
    feats[38] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_du', node2.lemma)
            
    # L->97->L
    feats[39] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '97') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '97', node2.lemma)
            
    # L->15_sg->L
    feats[40] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_sg', node2.lemma)
            
    # L->7_du->L
    feats[41] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_du', node2.lemma)
            
    # L->11_sp->L
    feats[42] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_sp', node2.lemma)
            
    # L->89->L
    feats[43] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '89') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '89', node2.lemma)
            
    # L->-261->L
    feats[44] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-261') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-261', node2.lemma)
            
    # L->16_fp->L
    feats[45] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_fp', node2.lemma)
            
    # L->178->L
    feats[46] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '178') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '178', node2.lemma)
            
    # L->-200->L
    feats[47] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-200', node2.lemma)
            
    # L->4_sp->L
    feats[48] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_sp', node2.lemma)
            
    # L->adj->L
    feats[49] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'adj') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'adj', node2.lemma)
            
    # L->10_sp->L
    feats[50] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_sp', node2.lemma)
            
    # L->-22->L
    feats[51] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-22') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-22', node2.lemma)
            
    # L->-46->L
    feats[52] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-46') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-46', node2.lemma)
            
    # L->30_fp->L
    feats[53] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_fp', node2.lemma)
            
    # L->instr. pl.->L
    feats[54] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. pl.', node2.lemma)
            
    # L->nom. adj.->L
    feats[55] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. adj.', node2.lemma)
            
    # L->-157->L
    feats[56] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-157', node2.lemma)
            
    # L->-27->L
    feats[57] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-27') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-27', node2.lemma)
            
    # L->-35->L
    feats[58] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-35') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
            
    # L->-243->L
    feats[59] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-243') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-243', node2.lemma)
            
    # L->-18->L
    feats[60] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-18', node2.lemma)
            
    # L->voc. sg.->L
    feats[61] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. sg.', node2.lemma)
            
    # L->49->L
    feats[62] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '49') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '49', node2.lemma)
            
    # L->-67->L
    feats[63] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-67', node2.lemma)
            
    # L->-89->L
    feats[64] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-89') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-89', node2.lemma)
            
    # L->-302->L
    feats[65] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-302', node2.lemma)
            
    # L->11_tp->L
    feats[66] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_tp', node2.lemma)
            
    # L->-139->L
    feats[67] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-139', node2.lemma)
            
    # L->28->L
    feats[68] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
            
    # L->-56->L
    feats[69] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-56') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-56', node2.lemma)
            
    # L->acc. masc.->L
    feats[70] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. masc.', node2.lemma)
            
    # L->-47->L
    feats[71] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-47') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-47', node2.lemma)
            
    # L->-17->L
    feats[72] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-17') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-17', node2.lemma)
            
    # L->5_du->L
    feats[73] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_du', node2.lemma)
            
    # L->98->L
    feats[74] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '98') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '98', node2.lemma)
            
    # L->81->L
    feats[75] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '81') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '81', node2.lemma)
            
    # L->-113->L
    feats[76] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-113') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-113', node2.lemma)
            
    # L->sg_tp->L
    feats[77] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_tp', node2.lemma)
            
    # L->-169->L
    feats[78] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-169') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-169', node2.lemma)
            
    # L->-283->L
    feats[79] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-283') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-283', node2.lemma)
            
    # L->-210->L
    feats[80] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-210', node2.lemma)
            
    # L->-25->L
    feats[81] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-25') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-25', node2.lemma)
            
    # L->162->L
    feats[82] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '162') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '162', node2.lemma)
            
    # L->-101->L
    feats[83] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-101') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-101', node2.lemma)
            
    # L->instr. adj.->L
    feats[84] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. adj.', node2.lemma)
            
    # L->du_fp->L
    feats[85] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_fp', node2.lemma)
            
    # L->-303->L
    feats[86] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-303', node2.lemma)
            
    # L->5_tp->L
    feats[87] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_tp', node2.lemma)
            
    # L->175->L
    feats[88] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '175') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '175', node2.lemma)
            
    # L->11_du->L
    feats[89] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_du', node2.lemma)
            
    # L->15_sp->L
    feats[90] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_sp', node2.lemma)
            
    # L->30_du->L
    feats[91] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_du', node2.lemma)
            
    # L->128->L
    feats[92] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '128') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '128', node2.lemma)
            
    # L->-273->L
    feats[93] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-273', node2.lemma)
            
    # L->7_pl->L
    feats[94] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_pl', node2.lemma)
            
    # L->acc. pl.->L
    feats[95] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. pl.', node2.lemma)
            
    # L->-103->L
    feats[96] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-103') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-103', node2.lemma)
            
    # L->71->L
    feats[97] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '71') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '71', node2.lemma)
            
    # L->136->L
    feats[98] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '136') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '136', node2.lemma)
            
    # L->-306->L
    feats[99] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-306', node2.lemma)
            
    # L->dat. du.->L
    feats[100] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
            
    # L->109->L
    feats[101] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '109') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '109', node2.lemma)
            
    # L->-291->L
    feats[102] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-291') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-291', node2.lemma)
            
    # L->-57->L
    feats[103] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-57') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-57', node2.lemma)
            
    # L->-161->L
    feats[104] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-161') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-161', node2.lemma)
            
    # L->13_fp->L
    feats[105] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_fp', node2.lemma)
            
    # L->-266->L
    feats[106] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-266') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-266', node2.lemma)
            
    # L->5_sp->L
    feats[107] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sp', node2.lemma)
            
    # L->181->L
    feats[108] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '181') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '181', node2.lemma)
            
    # L->-111->L
    feats[109] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-111') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-111', node2.lemma)
            
    # L->3_du->L
    feats[110] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
            
    # L->4_tp->L
    feats[111] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_tp', node2.lemma)
            
    # L->112->L
    feats[112] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '112') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '112', node2.lemma)
            
    # L->pl_fp->L
    feats[113] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_fp', node2.lemma)
            
    # L->28_sg->L
    feats[114] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28_sg', node2.lemma)
            
    # L->138->L
    feats[115] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '138') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '138', node2.lemma)
            
    # L->-78->L
    feats[116] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-78') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-78', node2.lemma)
            
    # L->-54->L
    feats[117] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-54') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-54', node2.lemma)
            
    # L->176->L
    feats[118] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '176') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '176', node2.lemma)
            
    # L->2_tp->L
    feats[119] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_tp', node2.lemma)
            
    # L->100->L
    feats[120] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '100') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '100', node2.lemma)
            
    # L->-99->L
    feats[121] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-99') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-99', node2.lemma)
            
    # L->14_pl->L
    feats[122] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_pl', node2.lemma)
            
    # L->37->L
    feats[123] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '37') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '37', node2.lemma)
            
    # L->-44->L
    feats[124] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-44') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-44', node2.lemma)
            
    # L->95->L
    feats[125] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '95') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '95', node2.lemma)
            
    # L->-83->L
    feats[126] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-83') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-83', node2.lemma)
            
    # L->76->L
    feats[127] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '76') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '76', node2.lemma)
            
    # L->36->L
    feats[128] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '36') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '36', node2.lemma)
            
    # L->-68->L
    feats[129] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-68') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-68', node2.lemma)
            
    # L->du->L
    feats[130] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du', node2.lemma)
            
    # L->voc->L
    feats[131] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc', node2.lemma)
            
    # L->-53->L
    feats[132] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-53') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-53', node2.lemma)
            
    # L->48->L
    feats[133] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '48') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '48', node2.lemma)
            
    # L->156->L
    feats[134] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '156') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '156', node2.lemma)
            
    # L->-72->L
    feats[135] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-72') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-72', node2.lemma)
            
    # L->-158->L
    feats[136] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-158') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-158', node2.lemma)
            
    # L->-293->L
    feats[137] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-293', node2.lemma)
            
    # L->34->L
    feats[138] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '34') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '34', node2.lemma)
            
    # L->instr. masc.->L
    feats[139] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. masc.', node2.lemma)
            
    # L->7_sp->L
    feats[140] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_sp', node2.lemma)
            
    # L->2_fp->L
    feats[141] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_fp', node2.lemma)
            
    # L->152->L
    feats[142] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '152') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '152', node2.lemma)
            
    # L->-43->L
    feats[143] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-43') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-43', node2.lemma)
            
    # L->-149->L
    feats[144] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-149') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-149', node2.lemma)
            
    # L->nom. du.->L
    feats[145] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. du.', node2.lemma)
            
    # L->30_sp->L
    feats[146] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_sp', node2.lemma)
            
    # L->-143->L
    feats[147] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-143') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-143', node2.lemma)
            
    # L->-62->L
    feats[148] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-62') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-62', node2.lemma)
            
    # L->180->L
    feats[149] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '180') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '180', node2.lemma)
            
    # L->16_sg->L
    feats[150] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_sg', node2.lemma)
            
    # L->108->L
    feats[151] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '108') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '108', node2.lemma)
            
    # L->16_du->L
    feats[152] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_du', node2.lemma)
            
    # L->5_sg->L
    feats[153] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sg', node2.lemma)
            
    # L->32->L
    feats[154] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '32') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '32', node2.lemma)
            
    # L->27_tp->L
    feats[155] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_tp', node2.lemma)
            
    # L->54->L
    feats[156] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '54') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
            
    # L->acc. du.->L
    feats[157] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. du.', node2.lemma)
            
    # L->-48->L
    feats[158] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-48', node2.lemma)
            
    # L->-126->L
    feats[159] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-126') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-126', node2.lemma)
            
    # L->10_du->L
    feats[160] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_du', node2.lemma)
            
    # L->-61->L
    feats[161] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-61') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-61', node2.lemma)
            
    # L->16_tp->L
    feats[162] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_tp', node2.lemma)
            
    # L->-90->L
    feats[163] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-90') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-90', node2.lemma)
            
    # L->abl. sg.->L
    feats[164] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'abl. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. sg.', node2.lemma)
            
    # L->acc. sg.->L
    feats[165] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. sg.', node2.lemma)
            
    # L->4_fp->L
    feats[166] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_fp', node2.lemma)
            
    # L->-71->L
    feats[167] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-71') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-71', node2.lemma)
            
    # L->-42->L
    feats[168] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-42') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-42', node2.lemma)
            
    # L->6_sp->L
    feats[169] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_sp', node2.lemma)
            
    # L->3_sp->L
    feats[170] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_sp', node2.lemma)
            
    # L->51->L
    feats[171] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '51') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '51', node2.lemma)
            
    # L->dat. pl.->L
    feats[172] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. pl.', node2.lemma)
            
    # L->du_sp->L
    feats[173] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_sp', node2.lemma)
            
    # L->179->L
    feats[174] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '179') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '179', node2.lemma)
            
    # L->pl_tp->L
    feats[175] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_tp', node2.lemma)
            
    # L->79->L
    feats[176] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '79', node2.lemma)
            
    # L->30_tp->L
    feats[177] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_tp', node2.lemma)
            
    # L->-133->L
    feats[178] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-133') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-133', node2.lemma)
            
    # L->5_pl->L
    feats[179] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_pl', node2.lemma)
            
    # L->-144->L
    feats[180] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-144') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-144', node2.lemma)
            
    # L->7_sg->L
    feats[181] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_sg', node2.lemma)
            
    # L->29->L
    feats[182] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '29') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '29', node2.lemma)
            
    # L->30_pl->L
    feats[183] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_pl', node2.lemma)
            
    # L->149->L
    feats[184] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '149', node2.lemma)
            
    # L->182->L
    feats[185] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '182') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '182', node2.lemma)
            
    # L->-292->L
    feats[186] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-292') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-292', node2.lemma)
            
    # L->6_pl->L
    feats[187] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_pl', node2.lemma)
            
    # L->-269->L
    feats[188] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-269') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-269', node2.lemma)
            
    # L->72->L
    feats[189] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '72') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '72', node2.lemma)
            
    # L->instr. du.->L
    feats[190] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. du.', node2.lemma)
            
    # L->sp->L
    feats[191] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sp', node2.lemma)
            
    # L->173->L
    feats[192] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '173') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '173', node2.lemma)
            
    # L->29_sg->L
    feats[193] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '29_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '29_sg', node2.lemma)
            
    # L->7_fp->L
    feats[194] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_fp', node2.lemma)
            
    # L->-114->L
    feats[195] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-114') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-114', node2.lemma)
            
    # L->80->L
    feats[196] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '80') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '80', node2.lemma)
            
    # L->-119->L
    feats[197] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-119') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-119', node2.lemma)
            
    # L->-279->L
    feats[198] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-279') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-279', node2.lemma)
            
    # L->-36->L
    feats[199] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-36') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-36', node2.lemma)
            
    # L->-34->L
    feats[200] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-34', node2.lemma)
            
    # L->9_du->L
    feats[201] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_du', node2.lemma)
            
    # L->155->L
    feats[202] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '155') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '155', node2.lemma)
            
    # L->168->L
    feats[203] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '168') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '168', node2.lemma)
            
    # L->-307->L
    feats[204] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-307') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-307', node2.lemma)
            
    # L->110->L
    feats[205] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '110') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '110', node2.lemma)
            
    # L->69->L
    feats[206] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '69') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '69', node2.lemma)
            
    # L->3_sg->L
    feats[207] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_sg', node2.lemma)
            
    # L->28_tp->L
    feats[208] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28_tp', node2.lemma)
            
    # L->4_pl->L
    feats[209] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_pl', node2.lemma)
            
    # L->14_fp->L
    feats[210] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_fp', node2.lemma)
            
    # L->-64->L
    feats[211] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-64') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-64', node2.lemma)
            
    # L->-97->L
    feats[212] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-97') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-97', node2.lemma)
            
    # L->77->L
    feats[213] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '77') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '77', node2.lemma)
            
    # L->148->L
    feats[214] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '148') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '148', node2.lemma)
            
    # L->-41->L
    feats[215] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-41') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-41', node2.lemma)
            
    # L->pl->L
    feats[216] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
            
    # L->122->L
    feats[217] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '122') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '122', node2.lemma)
            
    # L->-117->L
    feats[218] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-117') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-117', node2.lemma)
            
    # L->-141->L
    feats[219] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-141') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-141', node2.lemma)
            
    # L->-142->L
    feats[220] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-142') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-142', node2.lemma)
            
    # L->134->L
    feats[221] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '134') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '134', node2.lemma)
            
    # L->2_sp->L
    feats[222] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_sp', node2.lemma)
            
    # L->voc. du.->L
    feats[223] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. du.', node2.lemma)
            
    # L->12_pl->L
    feats[224] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_pl', node2.lemma)
            
    # L->-301->L
    feats[225] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-301') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-301', node2.lemma)
            
    # L->-262->L
    feats[226] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-262') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-262', node2.lemma)
            
    # L->-150->L
    feats[227] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-150') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-150', node2.lemma)
            
    # L->6_tp->L
    feats[228] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_tp', node2.lemma)
            
    # L->instr. neutr.->L
    feats[229] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. neutr.', node2.lemma)
            
    # L->2_pl->L
    feats[230] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_pl', node2.lemma)
            
    # L->15_tp->L
    feats[231] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_tp', node2.lemma)
            
    # L->-93->L
    feats[232] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-93') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-93', node2.lemma)
            
    # L->neutr->L
    feats[233] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'neutr') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'neutr', node2.lemma)
            
    # L->73->L
    feats[234] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '73') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '73', node2.lemma)
            
    # L->voc. neutr.->L
    feats[235] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. neutr.', node2.lemma)
            
    # L->-246->L
    feats[236] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-246') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-246', node2.lemma)
            
    # L->3_fp->L
    feats[237] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_fp', node2.lemma)
            
    # L->88->L
    feats[238] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '88') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '88', node2.lemma)
            
    # L->96->L
    feats[239] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '96') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '96', node2.lemma)
            
    # L->59->L
    feats[240] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '59') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '59', node2.lemma)
            
    # L->loc->L
    feats[241] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc', node2.lemma)
            
    # L->-153->L
    feats[242] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-153') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-153', node2.lemma)
            
    # L->-23->L
    feats[243] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-23') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-23', node2.lemma)
            
    # L->-14->L
    feats[244] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-14') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-14', node2.lemma)
            
    # L->-45->L
    feats[245] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-45') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-45', node2.lemma)
            
    # L->nom. neutr.->L
    feats[246] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. neutr.', node2.lemma)
            
    # L->-159->L
    feats[247] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-159') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-159', node2.lemma)
            
    # L->-309->L
    feats[248] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-309') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-309', node2.lemma)
            
    # L->158->L
    feats[249] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '158') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '158', node2.lemma)
            
    # L->90->L
    feats[250] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '90') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '90', node2.lemma)
            
    # L->39->L
    feats[251] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '39', node2.lemma)
            
    # L->137->L
    feats[252] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '137') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
            
    # L->-10->L
    feats[253] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-10') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-10', node2.lemma)
            
    # L->12_fp->L
    feats[254] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_fp', node2.lemma)
            
    # L->-26->L
    feats[255] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-26') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-26', node2.lemma)
            
    # L->132->L
    feats[256] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '132') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '132', node2.lemma)
            
    # L->gen->L
    feats[257] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen', node2.lemma)
            
    # L->instr->L
    feats[258] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr', node2.lemma)
            
    # L->-59->L
    feats[259] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-59') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-59', node2.lemma)
            
    # L->10_fp->L
    feats[260] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_fp', node2.lemma)
            
    # L->3_tp->L
    feats[261] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_tp', node2.lemma)
            
    # L->114->L
    feats[262] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '114') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '114', node2.lemma)
            
    # L->8_sg->L
    feats[263] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_sg', node2.lemma)
            
    # L->-76->L
    feats[264] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-76') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-76', node2.lemma)
            
    # L->-24->L
    feats[265] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-24') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
            
    # L->2->L
    feats[266] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2', node2.lemma)
            
    # L->6_sg->L
    feats[267] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_sg', node2.lemma)
            
    # L->masc->L
    feats[268] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'masc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'masc', node2.lemma)
            
    # L->-96->L
    feats[269] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-96') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-96', node2.lemma)
            
    # L->-82->L
    feats[270] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-82') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-82', node2.lemma)
            
    # L->9_tp->L
    feats[271] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_tp', node2.lemma)
            
    # L->-109->L
    feats[272] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-109') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-109', node2.lemma)
            
    # L->abl->L
    feats[273] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'abl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl', node2.lemma)
            
    # L->4_du->L
    feats[274] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_du', node2.lemma)
            
    # L->abl. du.->L
    feats[275] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'abl. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. du.', node2.lemma)
            
    # L->4_sg->L
    feats[276] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_sg', node2.lemma)
            
    # L->60->L
    feats[277] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '60') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '60', node2.lemma)
            
    # L->50->L
    feats[278] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '50') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '50', node2.lemma)
            
    # L->-91->L
    feats[279] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-91') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-91', node2.lemma)
            
    # L->-263->L
    feats[280] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-263') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-263', node2.lemma)
            
    # L->voc. pl.->L
    feats[281] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. pl.', node2.lemma)
            
    # L->-31->L
    feats[282] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-31') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-31', node2.lemma)
            
    # L->6_du->L
    feats[283] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_du', node2.lemma)
            
    # L->30->L
    feats[284] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30', node2.lemma)
            
    # L->41->L
    feats[285] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '41') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '41', node2.lemma)
            
    # L->91->L
    feats[286] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '91') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
            
    # L->94->L
    feats[287] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '94') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '94', node2.lemma)
            
    # L->-69->L
    feats[288] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-69') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-69', node2.lemma)
            
    # L->131->L
    feats[289] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '131') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '131', node2.lemma)
            
    # L->10_tp->L
    feats[290] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_tp', node2.lemma)
            
    # L->140->L
    feats[291] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '140') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '140', node2.lemma)
            
    # L->tp->L
    feats[292] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'tp', node2.lemma)
            
    # L->159->L
    feats[293] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '159', node2.lemma)
            
    # L->gen. du.->L
    feats[294] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen. du.', node2.lemma)
            
    # L->-87->L
    feats[295] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-87') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-87', node2.lemma)
            
    # L->150->L
    feats[296] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '150') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '150', node2.lemma)
            
    # L->8_tp->L
    feats[297] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_tp', node2.lemma)
            
    # L->3->L
    feats[298] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3', node2.lemma)
            
    # L->-308->L
    feats[299] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-308', node2.lemma)
            
    # L->33->L
    feats[300] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '33') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '33', node2.lemma)
            
    # L->acc. fem->L
    feats[301] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. fem', node2.lemma)
            
    # L->9_sg->L
    feats[302] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_sg', node2.lemma)
            
    # L->129->L
    feats[303] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '129') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '129', node2.lemma)
            
    # L->-94->L
    feats[304] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-94') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-94', node2.lemma)
            
    # L->nom. sg.->L
    feats[305] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. sg.', node2.lemma)
            
    # L->142->L
    feats[306] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '142') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '142', node2.lemma)
            
    # L->nom. fem->L
    feats[307] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. fem', node2.lemma)
            
    # L->-276->L
    feats[308] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-276') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-276', node2.lemma)
            
    # L->111->L
    feats[309] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '111') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '111', node2.lemma)
            
    # L->fp->L
    feats[310] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'fp', node2.lemma)
            
    # L->nom. masc.->L
    feats[311] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. masc.', node2.lemma)
            
    # L->-271->L
    feats[312] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-271') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-271', node2.lemma)
            
    # L->5_fp->L
    feats[313] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_fp', node2.lemma)
            
    # L->174->L
    feats[314] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '174') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '174', node2.lemma)
            
    # L->99->L
    feats[315] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '99') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '99', node2.lemma)
            
    # L->-30->L
    feats[316] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-30', node2.lemma)
            
    # L->voc. fem->L
    feats[317] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. fem', node2.lemma)
            
    # L->-73->L
    feats[318] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-73') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-73', node2.lemma)
            
    # L->nom->L
    feats[319] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom', node2.lemma)
            
    # L->78->L
    feats[320] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '78') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '78', node2.lemma)
            
    # L->-13->L
    feats[321] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-13') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-13', node2.lemma)
            
    # L->115->L
    feats[322] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '115') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '115', node2.lemma)
            
    # L->15_pl->L
    feats[323] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_pl', node2.lemma)
            
    # L->70->L
    feats[324] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '70') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '70', node2.lemma)
            
    # L->35->L
    feats[325] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '35', node2.lemma)
            
    # L->-247->L
    feats[326] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-247') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-247', node2.lemma)
            
    # L->157->L
    feats[327] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '157') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '157', node2.lemma)
            
    # L->1->L
    feats[328] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '1') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '1', node2.lemma)
            
    # L->-249->L
    feats[329] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-249') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-249', node2.lemma)
            
    # L->-38->L
    feats[330] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-38') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-38', node2.lemma)
            
    # L->161->L
    feats[331] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '161') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '161', node2.lemma)
            
    # L->27_sg->L
    feats[332] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_sg', node2.lemma)
            
    # L->-55->L
    feats[333] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-55') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-55', node2.lemma)
            
    # L->acc->L
    feats[334] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc', node2.lemma)
            
    # L->-51->L
    feats[335] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-51') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-51', node2.lemma)
            
    # L->-16->L
    feats[336] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-16') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-16', node2.lemma)
            
    # L->-52->L
    feats[337] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-52') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-52', node2.lemma)
            
    # L->loc. pl.->L
    feats[338] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc. pl.', node2.lemma)
            
    # L->42->L
    feats[339] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '42') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '42', node2.lemma)
            
    # L->-84->L
    feats[340] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-84') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-84', node2.lemma)
            
    # L->27_fp->L
    feats[341] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_fp', node2.lemma)
            
    # L->-20->L
    feats[342] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-20') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-20', node2.lemma)
            
    # L->141->L
    feats[343] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '141') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '141', node2.lemma)
            
    # L->139->L
    feats[344] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '139') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '139', node2.lemma)
            
    # L->-129->L
    feats[345] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-129') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-129', node2.lemma)
            
    # L->120->L
    feats[346] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '120') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '120', node2.lemma)
            
    # L->-39->L
    feats[347] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-39') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-39', node2.lemma)
            
    # L->-166->L
    feats[348] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-166') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-166', node2.lemma)
            
    # L->-28->L
    feats[349] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-28') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-28', node2.lemma)
            
    # L->170->L
    feats[350] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '170') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '170', node2.lemma)
            
    # L->-131->L
    feats[351] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-131') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-131', node2.lemma)
            
    # L->-19->L
    feats[352] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-19') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-19', node2.lemma)
            
    # L->2_du->L
    feats[353] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_du', node2.lemma)
            
    # L->-112->L
    feats[354] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-112') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
            
    # L->12_sp->L
    feats[355] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_sp', node2.lemma)
            
    # L->-104->L
    feats[356] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-104', node2.lemma)
            
    # L->-12->L
    feats[357] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-12') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-12', node2.lemma)
            
    # L->acc. neutr.->L
    feats[358] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. neutr.', node2.lemma)
            
    # L->14_tp->L
    feats[359] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_tp', node2.lemma)
            
    # L->-29->L
    feats[360] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-29') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-29', node2.lemma)
            
    # L->sg->L
    feats[361] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg', node2.lemma)
            
    # L->177->L
    feats[362] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '177') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '177', node2.lemma)
            
    # L->7_tp->L
    feats[363] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_tp', node2.lemma)
            
    # L->-81->L
    feats[364] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-81', node2.lemma)
            
    # L->-92->L
    feats[365] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-92', node2.lemma)
            
    # L->-230->L
    feats[366] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-230', node2.lemma)
            
    # L->-297->L
    feats[367] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-297') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-297', node2.lemma)
            
    # L->voc. masc.->L
    feats[368] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. masc.', node2.lemma)
            
    # L->56->L
    feats[369] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '56') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '56', node2.lemma)
            
    # L->gen. sg.->L
    feats[370] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen. sg.', node2.lemma)
            
    # L->172->L
    feats[371] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '172') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '172', node2.lemma)
            
    # L->11_pl->L
    feats[372] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_pl', node2.lemma)
            
    # L->-11->L
    feats[373] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-11') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-11', node2.lemma)
            
    # L->13_tp->L
    feats[374] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_tp', node2.lemma)
            
    # L->160->L
    feats[375] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '160') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
            
    # L->-115->L
    feats[376] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-115') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-115', node2.lemma)
            
    # L->-50->L
    feats[377] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-50') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-50', node2.lemma)
            
    # L->6_fp->L
    feats[378] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_fp', node2.lemma)
            
    # L->68->L
    feats[379] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '68') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '68', node2.lemma)
            
    # L->-15->C
    feats[380] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-15') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-15', node2.cng)
            
    # L->-121->C
    feats[381] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-121') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-121', node2.cng)
            
    # L->117->C
    feats[382] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '117') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '117', node2.cng)
            
    # L->-240->C
    feats[383] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-240') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-240', node2.cng)
            
    # L->55->C
    feats[384] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '55') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '55', node2.cng)
            
    # L->-245->C
    feats[385] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-245') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-245', node2.cng)
            
    # L->171->C
    feats[386] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '171') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '171', node2.cng)
            
    # L->154->C
    feats[387] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '154') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '154', node2.cng)
            
    # L->14_sp->C
    feats[388] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_sp', node2.cng)
            
    # L->du_tp->C
    feats[389] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', node2.cng)
            
    # L->40->C
    feats[390] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '40') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '40', node2.cng)
            
    # L->-79->C
    feats[391] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-79', node2.cng)
            
    # L->9_fp->C
    feats[392] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_fp', node2.cng)
            
    # L->12_tp->C
    feats[393] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_tp', node2.cng)
            
    # L->29_tp->C
    feats[394] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '29_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '29_tp', node2.cng)
            
    # L->pl_sp->C
    feats[395] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_sp', node2.cng)
            
    # L->-132->C
    feats[396] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-132') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-132', node2.cng)
            
    # L->-122->C
    feats[397] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-122') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-122', node2.cng)
            
    # L->102->C
    feats[398] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '102') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '102', node2.cng)
            
    # L->-152->C
    feats[399] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-152') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-152', node2.cng)
            
    # L->3_pl->C
    feats[400] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_pl', node2.cng)
            
    # L->151->C
    feats[401] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '151') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '151', node2.cng)
            
    # L->15_du->C
    feats[402] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_du', node2.cng)
            
    # L->82->C
    feats[403] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '82') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '82', node2.cng)
            
    # L->-242->C
    feats[404] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-242') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-242', node2.cng)
            
    # L->92->C
    feats[405] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '92', node2.cng)
            
    # L->74->C
    feats[406] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '74') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '74', node2.cng)
            
    # L->97->C
    feats[407] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '97') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '97', node2.cng)
            
    # L->15_sg->C
    feats[408] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_sg', node2.cng)
            
    # L->7_du->C
    feats[409] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_du', node2.cng)
            
    # L->11_sp->C
    feats[410] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_sp', node2.cng)
            
    # L->89->C
    feats[411] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '89') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '89', node2.cng)
            
    # L->16_fp->C
    feats[412] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_fp', node2.cng)
            
    # L->178->C
    feats[413] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '178') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '178', node2.cng)
            
    # L->sg_fp->C
    feats[414] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg_fp', node2.cng)
            
    # L->4_sp->C
    feats[415] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_sp', node2.cng)
            
    # L->adj->C
    feats[416] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'adj') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'adj', node2.cng)
            
    # L->10_sp->C
    feats[417] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sp', node2.cng)
            
    # L->-22->C
    feats[418] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-22') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-22', node2.cng)
            
    # L->nom. adj.->C
    feats[419] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. adj.', node2.cng)
            
    # L->-27->C
    feats[420] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-27') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-27', node2.cng)
            
    # L->-35->C
    feats[421] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-35', node2.cng)
            
    # L->-243->C
    feats[422] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-243') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-243', node2.cng)
            
    # L->voc. sg.->C
    feats[423] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. sg.', node2.cng)
            
    # L->49->C
    feats[424] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '49') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '49', node2.cng)
            
    # L->-89->C
    feats[425] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-89') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-89', node2.cng)
            
    # L->-302->C
    feats[426] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', node2.cng)
            
    # L->-139->C
    feats[427] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', node2.cng)
            
    # L->28->C
    feats[428] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28', node2.cng)
            
    # L->-56->C
    feats[429] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-56') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-56', node2.cng)
            
    # L->-47->C
    feats[430] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-47') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-47', node2.cng)
            
    # L->-210->C
    feats[431] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', node2.cng)
            
    # L->-25->C
    feats[432] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-25') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-25', node2.cng)
            
    # L->-101->C
    feats[433] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-101') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-101', node2.cng)
            
    # L->instr. adj.->C
    feats[434] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. adj.', node2.cng)
            
    # L->30_du->C
    feats[435] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_du', node2.cng)
            
    # L->128->C
    feats[436] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '128') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '128', node2.cng)
            
    # L->-273->C
    feats[437] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', node2.cng)
            
    # L->acc. pl.->C
    feats[438] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. pl.', node2.cng)
            
    # L->-306->C
    feats[439] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', node2.cng)
            
    # L->dat. du.->C
    feats[440] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. du.', node2.cng)
            
    # L->109->C
    feats[441] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '109') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '109', node2.cng)
            
    # L->13_fp->C
    feats[442] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_fp', node2.cng)
            
    # L->181->C
    feats[443] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '181') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '181', node2.cng)
            
    # L->-111->C
    feats[444] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-111') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-111', node2.cng)
            
    # L->3_du->C
    feats[445] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_du', node2.cng)
            
    # L->-37->C
    feats[446] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-37', node2.cng)
            
    # L->28_sg->C
    feats[447] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28_sg', node2.cng)
            
    # L->138->C
    feats[448] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '138') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '138', node2.cng)
            
    # L->-78->C
    feats[449] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-78') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-78', node2.cng)
            
    # L->-54->C
    feats[450] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-54') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-54', node2.cng)
            
    # L->-44->C
    feats[451] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-44') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-44', node2.cng)
            
    # L->-83->C
    feats[452] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-83') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-83', node2.cng)
            
    # L->36->C
    feats[453] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '36') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '36', node2.cng)
            
    # L->-68->C
    feats[454] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-68') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-68', node2.cng)
            
    # L->du->C
    feats[455] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du', node2.cng)
            
    # L->voc->C
    feats[456] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc', node2.cng)
            
    # L->-53->C
    feats[457] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-53') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-53', node2.cng)
            
    # L->48->C
    feats[458] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '48', node2.cng)
            
    # L->156->C
    feats[459] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '156', node2.cng)
            
    # L->-72->C
    feats[460] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-72') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-72', node2.cng)
            
    # L->-158->C
    feats[461] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-158') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-158', node2.cng)
            
    # L->34->C
    feats[462] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '34', node2.cng)
            
    # L->instr. masc.->C
    feats[463] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. masc.', node2.cng)
            
    # L->-149->C
    feats[464] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-149', node2.cng)
            
    # L->30_sp->C
    feats[465] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_sp', node2.cng)
            
    # L->108->C
    feats[466] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '108') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '108', node2.cng)
            
    # L->16_du->C
    feats[467] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_du', node2.cng)
            
    # L->5_sg->C
    feats[468] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sg', node2.cng)
            
    # L->27_tp->C
    feats[469] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_tp', node2.cng)
            
    # L->54->C
    feats[470] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '54') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '54', node2.cng)
            
    # L->-48->C
    feats[471] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-48', node2.cng)
            
    # L->10_du->C
    feats[472] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_du', node2.cng)
            
    # L->16_tp->C
    feats[473] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', node2.cng)
            
    # L->-90->C
    feats[474] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-90', node2.cng)
            
    # L->abl. sg.->C
    feats[475] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'abl. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'abl. sg.', node2.cng)
            
    # L->4_fp->C
    feats[476] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_fp', node2.cng)
            
    # L->-42->C
    feats[477] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-42') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-42', node2.cng)
            
    # L->6_sp->C
    feats[478] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sp', node2.cng)
            
    # L->dat. pl.->C
    feats[479] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. pl.', node2.cng)
            
    # L->du_sp->C
    feats[480] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_sp', node2.cng)
            
    # L->pl_tp->C
    feats[481] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_tp', node2.cng)
            
    # L->-133->C
    feats[482] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-133') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-133', node2.cng)
            
    # L->7_sg->C
    feats[483] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_sg', node2.cng)
            
    # L->30_pl->C
    feats[484] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_pl', node2.cng)
            
    # L->182->C
    feats[485] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '182') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '182', node2.cng)
            
    # L->-269->C
    feats[486] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-269') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-269', node2.cng)
            
    # L->7_fp->C
    feats[487] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_fp', node2.cng)
            
    # L->-114->C
    feats[488] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-114') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-114', node2.cng)
            
    # L->-279->C
    feats[489] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-279') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-279', node2.cng)
            
    # L->-36->C
    feats[490] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-36') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-36', node2.cng)
            
    # L->-34->C
    feats[491] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', node2.cng)
            
    # L->155->C
    feats[492] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '155') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '155', node2.cng)
            
    # L->-307->C
    feats[493] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-307') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-307', node2.cng)
            
    # L->110->C
    feats[494] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '110') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '110', node2.cng)
            
    # L->3_sg->C
    feats[495] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_sg', node2.cng)
            
    # L->28_tp->C
    feats[496] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28_tp', node2.cng)
            
    # L->4_pl->C
    feats[497] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_pl', node2.cng)
            
    # L->14_fp->C
    feats[498] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_fp', node2.cng)
            
    # L->148->C
    feats[499] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '148') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '148', node2.cng)
            
    # L->122->C
    feats[500] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '122') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '122', node2.cng)
            
    # L->-117->C
    feats[501] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-117') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-117', node2.cng)
            
    # L->-141->C
    feats[502] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-141') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-141', node2.cng)
            
    # L->-142->C
    feats[503] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-142') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-142', node2.cng)
            
    # L->2_sp->C
    feats[504] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_sp', node2.cng)
            
    # L->-301->C
    feats[505] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-301') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-301', node2.cng)
            
    # L->-262->C
    feats[506] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-262') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-262', node2.cng)
            
    # L->-150->C
    feats[507] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-150') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-150', node2.cng)
            
    # L->2_pl->C
    feats[508] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_pl', node2.cng)
            
    # L->neutr->C
    feats[509] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'neutr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'neutr', node2.cng)
            
    # L->73->C
    feats[510] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '73') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '73', node2.cng)
            
    # L->voc. neutr.->C
    feats[511] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. neutr.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. neutr.', node2.cng)
            
    # L->-246->C
    feats[512] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-246') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-246', node2.cng)
            
    # L->88->C
    feats[513] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '88') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '88', node2.cng)
            
    # L->59->C
    feats[514] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '59') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '59', node2.cng)
            
    # L->-45->C
    feats[515] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-45') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-45', node2.cng)
            
    # L->13_sp->C
    feats[516] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_sp', node2.cng)
            
    # L->2_sg->C
    feats[517] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_sg', node2.cng)
            
    # L->-159->C
    feats[518] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-159', node2.cng)
            
    # L->-309->C
    feats[519] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-309') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-309', node2.cng)
            
    # L->158->C
    feats[520] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '158') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '158', node2.cng)
            
    # L->39->C
    feats[521] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', node2.cng)
            
    # L->132->C
    feats[522] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '132') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '132', node2.cng)
            
    # L->gen->C
    feats[523] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen', node2.cng)
            
    # L->instr->C
    feats[524] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', node2.cng)
            
    # L->-59->C
    feats[525] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-59') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-59', node2.cng)
            
    # L->10_fp->C
    feats[526] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_fp', node2.cng)
            
    # L->8_sg->C
    feats[527] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_sg', node2.cng)
            
    # L->-76->C
    feats[528] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-76') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-76', node2.cng)
            
    # L->-24->C
    feats[529] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-24') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-24', node2.cng)
            
    # L->6_sg->C
    feats[530] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', node2.cng)
            
    # L->masc->C
    feats[531] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'masc') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'masc', node2.cng)
            
    # L->-82->C
    feats[532] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-82') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-82', node2.cng)
            
    # L->abl->C
    feats[533] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'abl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'abl', node2.cng)
            
    # L->abl. du.->C
    feats[534] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'abl. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'abl. du.', node2.cng)
            
    # L->4_sg->C
    feats[535] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_sg', node2.cng)
            
    # L->50->C
    feats[536] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '50') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '50', node2.cng)
            
    # L->-91->C
    feats[537] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-91') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-91', node2.cng)
            
    # L->-263->C
    feats[538] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-263') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-263', node2.cng)
            
    # L->voc. pl.->C
    feats[539] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. pl.', node2.cng)
            
    # L->-31->C
    feats[540] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-31') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-31', node2.cng)
            
    # L->abl. pl.->C
    feats[541] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'abl. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'abl. pl.', node2.cng)
            
    # L->6_du->C
    feats[542] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_du', node2.cng)
            
    # L->41->C
    feats[543] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '41') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '41', node2.cng)
            
    # L->94->C
    feats[544] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '94') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '94', node2.cng)
            
    # L->-69->C
    feats[545] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-69') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-69', node2.cng)
            
    # L->-87->C
    feats[546] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-87') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-87', node2.cng)
            
    # L->8_tp->C
    feats[547] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_tp', node2.cng)
            
    # L->-308->C
    feats[548] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', node2.cng)
            
    # L->33->C
    feats[549] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '33', node2.cng)
            
    # L->acc. fem->C
    feats[550] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. fem', node2.cng)
            
    # L->-94->C
    feats[551] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-94') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-94', node2.cng)
            
    # L->nom. sg.->C
    feats[552] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. sg.', node2.cng)
            
    # L->nom. fem->C
    feats[553] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. fem', node2.cng)
            
    # L->-276->C
    feats[554] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-276') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-276', node2.cng)
            
    # L->111->C
    feats[555] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '111') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '111', node2.cng)
            
    # L->fp->C
    feats[556] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'fp', node2.cng)
            
    # L->nom. masc.->C
    feats[557] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. masc.', node2.cng)
            
    # L->5_fp->C
    feats[558] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', node2.cng)
            
    # L->99->C
    feats[559] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '99') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '99', node2.cng)
            
    # L->-30->C
    feats[560] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', node2.cng)
            
    # L->voc. fem->C
    feats[561] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. fem', node2.cng)
            
    # L->-73->C
    feats[562] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-73') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-73', node2.cng)
            
    # L->115->C
    feats[563] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '115') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '115', node2.cng)
            
    # L->15_pl->C
    feats[564] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_pl', node2.cng)
            
    # L->70->C
    feats[565] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '70') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '70', node2.cng)
            
    # L->35->C
    feats[566] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', node2.cng)
            
    # L->1->C
    feats[567] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '1') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '1', node2.cng)
            
    # L->-38->C
    feats[568] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-38') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-38', node2.cng)
            
    # L->161->C
    feats[569] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '161') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '161', node2.cng)
            
    # L->-55->C
    feats[570] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-55') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-55', node2.cng)
            
    # L->acc->C
    feats[571] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc', node2.cng)
            
    # L->-51->C
    feats[572] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-51') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-51', node2.cng)
            
    # L->-52->C
    feats[573] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-52') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-52', node2.cng)
            
    # L->loc. pl.->C
    feats[574] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. pl.', node2.cng)
            
    # L->42->C
    feats[575] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '42') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '42', node2.cng)
            
    # L->27_fp->C
    feats[576] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_fp', node2.cng)
            
    # L->-20->C
    feats[577] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-20') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-20', node2.cng)
            
    # L->139->C
    feats[578] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '139', node2.cng)
            
    # L->-166->C
    feats[579] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-166') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-166', node2.cng)
            
    # L->-28->C
    feats[580] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-28') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-28', node2.cng)
            
    # L->170->C
    feats[581] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '170') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '170', node2.cng)
            
    # L->-19->C
    feats[582] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-19') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-19', node2.cng)
            
    # L->-104->C
    feats[583] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', node2.cng)
            
    # L->-12->C
    feats[584] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-12') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-12', node2.cng)
            
    # L->acc. neutr.->C
    feats[585] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. neutr.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. neutr.', node2.cng)
            
    # L->-29->C
    feats[586] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-29') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-29', node2.cng)
            
    # L->sg->C
    feats[587] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg', node2.cng)
            
    # L->7_tp->C
    feats[588] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_tp', node2.cng)
            
    # L->-81->C
    feats[589] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', node2.cng)
            
    # L->-230->C
    feats[590] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', node2.cng)
            
    # L->-297->C
    feats[591] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-297') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-297', node2.cng)
            
    # L->voc. masc.->C
    feats[592] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. masc.', node2.cng)
            
    # L->-11->C
    feats[593] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-11') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-11', node2.cng)
            
    # L->160->C
    feats[594] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '160') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '160', node2.cng)
            
    # L->-115->C
    feats[595] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-115') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-115', node2.cng)
            
    # L->6_fp->C
    feats[596] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_fp', node2.cng)
            
    # L->121->T
    feats[597] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '121') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '121', node2.tup)
            
    # L->-15->T
    feats[598] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-15') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-15', node2.tup)
            
    # L->-121->T
    feats[599] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-121') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-121', node2.tup)
            
    # L->117->T
    feats[600] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '117') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '117', node2.tup)
            
    # L->10_sg->T
    feats[601] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sg') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '10_sg', node2.tup)
            
    # L->-240->T
    feats[602] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-240') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-240', node2.tup)
            
    # L->55->T
    feats[603] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '55') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '55', node2.tup)
            
    # L->-245->T
    feats[604] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-245') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-245', node2.tup)
            
    # L->171->T
    feats[605] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '171') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '171', node2.tup)
            
    # L->154->T
    feats[606] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '154') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '154', node2.tup)
            
    # L->14_sp->T
    feats[607] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_sp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '14_sp', node2.tup)
            
    # L->38->T
    feats[608] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '38') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '38', node2.tup)
            
    # L->-151->T
    feats[609] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-151') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-151', node2.tup)
            
    # L->du_tp->T
    feats[610] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_tp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'du_tp', node2.tup)
            
    # L->-299->T
    feats[611] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-299') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-299', node2.tup)
            
    # L->40->T
    feats[612] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '40') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '40', node2.tup)
            
    # L->61->T
    feats[613] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '61') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '61', node2.tup)
            
    # L->-79->T
    feats[614] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-79') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-79', node2.tup)
            
    # L->9_fp->T
    feats[615] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_fp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '9_fp', node2.tup)
            
    # L->12_tp->T
    feats[616] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_tp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '12_tp', node2.tup)
            
    # L->29_tp->T
    feats[617] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '29_tp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '29_tp', node2.tup)
            
    # L->pl_sp->T
    feats[618] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_sp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'pl_sp', node2.tup)
            
    # L->135->T
    feats[619] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '135') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '135', node2.tup)
            
    # L->-122->T
    feats[620] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-122') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-122', node2.tup)
            
    # L->acc. adj.->T
    feats[621] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. adj.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'acc. adj.', node2.tup)
            
    # L->-123->T
    feats[622] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-123') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-123', node2.tup)
            
    # L->102->T
    feats[623] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '102') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '102', node2.tup)
            
    # L->-152->T
    feats[624] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-152') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-152', node2.tup)
            
    # L->3_pl->T
    feats[625] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_pl') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '3_pl', node2.tup)
            
    # L->instr. sg.->T
    feats[626] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. sg.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'instr. sg.', node2.tup)
            
    # L->151->T
    feats[627] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '151') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '151', node2.tup)
            
    # L->15_du->T
    feats[628] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_du') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '15_du', node2.tup)
            
    # L->82->T
    feats[629] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '82') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '82', node2.tup)
            
    # L->instr. fem->T
    feats[630] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. fem') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'instr. fem', node2.tup)
            
    # L->-98->T
    feats[631] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-98') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-98', node2.tup)
            
    # L->-242->T
    feats[632] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-242') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-242', node2.tup)
            
    # L->27_du->T
    feats[633] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '27_du', node2.tup)
            
    # L->-156->T
    feats[634] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-156') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-156', node2.tup)
            
    # L->92->T
    feats[635] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '92') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '92', node2.tup)
            
    # L->gen. pl.->T
    feats[636] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. pl.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'gen. pl.', node2.tup)
            
    # L->74->T
    feats[637] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '74') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '74', node2.tup)
            
    # L->-86->T
    feats[638] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-86') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-86', node2.tup)
            
    # L->8_du->T
    feats[639] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_du') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '8_du', node2.tup)
            
    # L->loc. sg.->T
    feats[640] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'loc. sg.', node2.tup)
            
    # L->97->T
    feats[641] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '97') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '97', node2.tup)
            
    # L->7_du->T
    feats[642] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_du') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '7_du', node2.tup)
            
    # L->11_sp->T
    feats[643] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_sp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '11_sp', node2.tup)
            
    # L->89->T
    feats[644] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '89') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '89', node2.tup)
            
    # L->-261->T
    feats[645] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-261') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-261', node2.tup)
            
    # L->16_fp->T
    feats[646] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_fp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '16_fp', node2.tup)
            
    # L->178->T
    feats[647] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '178') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '178', node2.tup)
            
    # L->12_sg->T
    feats[648] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '12_sg', node2.tup)
            
    # L->sg_fp->T
    feats[649] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg_fp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'sg_fp', node2.tup)
            
    # L->-200->T
    feats[650] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-200') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-200', node2.tup)
            
    # L->4_sp->T
    feats[651] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_sp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '4_sp', node2.tup)
            
    # L->adj->T
    feats[652] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'adj') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'adj', node2.tup)
            
    # L->10_sp->T
    feats[653] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_sp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '10_sp', node2.tup)
            
    # L->-22->T
    feats[654] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-22') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-22', node2.tup)
            
    # L->-46->T
    feats[655] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-46') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-46', node2.tup)
            
    # L->30_fp->T
    feats[656] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_fp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '30_fp', node2.tup)
            
    # L->instr. pl.->T
    feats[657] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. pl.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'instr. pl.', node2.tup)
            
    # L->nom. adj.->T
    feats[658] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. adj.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'nom. adj.', node2.tup)
            
    # L->-157->T
    feats[659] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-157') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-157', node2.tup)
            
    # L->-27->T
    feats[660] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-27') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-27', node2.tup)
            
    # L->-35->T
    feats[661] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-35') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-35', node2.tup)
            
    # L->-243->T
    feats[662] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-243') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-243', node2.tup)
            
    # L->-18->T
    feats[663] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-18') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-18', node2.tup)
            
    # L->voc. sg.->T
    feats[664] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. sg.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'voc. sg.', node2.tup)
            
    # L->49->T
    feats[665] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '49') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '49', node2.tup)
            
    # L->-67->T
    feats[666] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-67', node2.tup)
            
    # L->-89->T
    feats[667] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-89') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-89', node2.tup)
            
    # L->-302->T
    feats[668] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-302', node2.tup)
            
    # L->11_tp->T
    feats[669] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_tp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '11_tp', node2.tup)
            
    # L->-139->T
    feats[670] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-139') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-139', node2.tup)
            
    # L->28->T
    feats[671] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '28', node2.tup)
            
    # L->-56->T
    feats[672] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-56') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-56', node2.tup)
            
    # L->acc. masc.->T
    feats[673] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. masc.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'acc. masc.', node2.tup)
            
    # L->-47->T
    feats[674] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-47') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-47', node2.tup)
            
    # L->-17->T
    feats[675] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-17') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-17', node2.tup)
            
    # L->98->T
    feats[676] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '98') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '98', node2.tup)
            
    # L->81->T
    feats[677] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '81') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '81', node2.tup)
            
    # L->-113->T
    feats[678] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-113') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-113', node2.tup)
            
    # L->sg_tp->T
    feats[679] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg_tp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'sg_tp', node2.tup)
            
    # L->-169->T
    feats[680] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-169') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-169', node2.tup)
            
    # L->-283->T
    feats[681] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-283') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-283', node2.tup)
            
    # L->-210->T
    feats[682] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-210') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-210', node2.tup)
            
    # L->-25->T
    feats[683] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-25') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-25', node2.tup)
            
    # L->162->T
    feats[684] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '162') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '162', node2.tup)
            
    # L->-101->T
    feats[685] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-101') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-101', node2.tup)
            
    # L->instr. adj.->T
    feats[686] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. adj.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'instr. adj.', node2.tup)
            
    # L->du_fp->T
    feats[687] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_fp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'du_fp', node2.tup)
            
    # L->-303->T
    feats[688] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-303') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-303', node2.tup)
            
    # L->5_tp->T
    feats[689] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_tp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '5_tp', node2.tup)
            
    # L->175->T
    feats[690] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '175') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '175', node2.tup)
            
    # L->11_du->T
    feats[691] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_du') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '11_du', node2.tup)
            
    # L->15_sp->T
    feats[692] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_sp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '15_sp', node2.tup)
            
    # L->30_du->T
    feats[693] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_du') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '30_du', node2.tup)
            
    # L->128->T
    feats[694] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '128') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '128', node2.tup)
            
    # L->-273->T
    feats[695] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-273') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-273', node2.tup)
            
    # L->7_pl->T
    feats[696] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '7_pl', node2.tup)
            
    # L->acc. pl.->T
    feats[697] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. pl.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'acc. pl.', node2.tup)
            
    # L->-103->T
    feats[698] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-103') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-103', node2.tup)
            
    # L->71->T
    feats[699] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '71') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '71', node2.tup)
            
    # L->136->T
    feats[700] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '136') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '136', node2.tup)
            
    # L->sg_sp->T
    feats[701] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg_sp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'sg_sp', node2.tup)
            
    # L->-306->T
    feats[702] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-306') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-306', node2.tup)
            
    # L->dat. du.->T
    feats[703] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. du.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'dat. du.', node2.tup)
            
    # L->109->T
    feats[704] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '109') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '109', node2.tup)
            
    # L->-291->T
    feats[705] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-291') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-291', node2.tup)
            
    # L->-57->T
    feats[706] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-57') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-57', node2.tup)
            
    # L->-161->T
    feats[707] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-161') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-161', node2.tup)
            
    # L->13_fp->T
    feats[708] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_fp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '13_fp', node2.tup)
            
    # L->-266->T
    feats[709] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-266') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-266', node2.tup)
            
    # L->30_sg->T
    feats[710] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_sg') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '30_sg', node2.tup)
            
    # L->5_sp->T
    feats[711] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '5_sp', node2.tup)
            
    # L->181->T
    feats[712] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '181') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '181', node2.tup)
            
    # L->-111->T
    feats[713] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-111') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-111', node2.tup)
            
    # L->3_du->T
    feats[714] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_du') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '3_du', node2.tup)
            
    # L->119->T
    feats[715] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '119') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '119', node2.tup)
            
    # L->4_tp->T
    feats[716] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_tp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '4_tp', node2.tup)
            
    # L->-33->T
    feats[717] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-33') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-33', node2.tup)
            
    # L->-37->T
    feats[718] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-37') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-37', node2.tup)
            
    # L->112->T
    feats[719] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '112') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '112', node2.tup)
            
    # L->pl_fp->T
    feats[720] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_fp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'pl_fp', node2.tup)
            
    # L->28_sg->T
    feats[721] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_sg') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '28_sg', node2.tup)
            
    # L->fem->T
    feats[722] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'fem') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'fem', node2.tup)
            
    # L->-78->T
    feats[723] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-78') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-78', node2.tup)
            
    # L->-54->T
    feats[724] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-54') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-54', node2.tup)
            
    # L->176->T
    feats[725] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '176') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '176', node2.tup)
            
    # L->2_tp->T
    feats[726] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_tp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '2_tp', node2.tup)
            
    # L->100->T
    feats[727] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '100') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '100', node2.tup)
            
    # L->-99->T
    feats[728] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-99') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-99', node2.tup)
            
    # L->-163->T
    feats[729] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-163') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-163', node2.tup)
            
    # L->14_pl->T
    feats[730] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_pl') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '14_pl', node2.tup)
            
    # L->11_fp->T
    feats[731] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_fp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '11_fp', node2.tup)
            
    # L->37->T
    feats[732] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '37') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '37', node2.tup)
            
    # L->-44->T
    feats[733] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-44') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-44', node2.tup)
            
    # L->95->T
    feats[734] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '95') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '95', node2.tup)
            
    # L->-83->T
    feats[735] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-83') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-83', node2.tup)
            
    # L->76->T
    feats[736] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '76') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '76', node2.tup)
            
    # L->36->T
    feats[737] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '36') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '36', node2.tup)
            
    # L->-68->T
    feats[738] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-68') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-68', node2.tup)
            
    # L->du->T
    feats[739] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'du', node2.tup)
            
    # L->voc->T
    feats[740] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'voc', node2.tup)
            
    # L->-53->T
    feats[741] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-53') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-53', node2.tup)
            
    # L->48->T
    feats[742] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '48') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '48', node2.tup)
            
    # L->156->T
    feats[743] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '156') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '156', node2.tup)
            
    # L->-158->T
    feats[744] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-158') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-158', node2.tup)
            
    # L->-293->T
    feats[745] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-293', node2.tup)
            
    # L->34->T
    feats[746] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '34') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '34', node2.tup)
            
    # L->instr. masc.->T
    feats[747] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. masc.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'instr. masc.', node2.tup)
            
    # L->7_sp->T
    feats[748] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_sp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '7_sp', node2.tup)
            
    # L->2_fp->T
    feats[749] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '2_fp', node2.tup)
            
    # L->152->T
    feats[750] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '152') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '152', node2.tup)
            
    # L->-43->T
    feats[751] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-43') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-43', node2.tup)
            
    # L->-149->T
    feats[752] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-149') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-149', node2.tup)
            
    # L->nom. du.->T
    feats[753] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. du.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'nom. du.', node2.tup)
            
    # L->30_sp->T
    feats[754] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_sp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '30_sp', node2.tup)
            
    # L->-143->T
    feats[755] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-143') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-143', node2.tup)
            
    # L->-62->T
    feats[756] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-62') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-62', node2.tup)
            
    # L->180->T
    feats[757] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '180') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '180', node2.tup)
            
    # L->16_sg->T
    feats[758] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_sg') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '16_sg', node2.tup)
            
    # L->108->T
    feats[759] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '108') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '108', node2.tup)
            
    # L->16_du->T
    feats[760] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_du') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '16_du', node2.tup)
            
    # L->5_sg->T
    feats[761] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sg') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '5_sg', node2.tup)
            
    # L->32->T
    feats[762] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '32') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '32', node2.tup)
            
    # L->27_tp->T
    feats[763] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_tp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '27_tp', node2.tup)
            
    # L->54->T
    feats[764] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '54') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '54', node2.tup)
            
    # L->acc. du.->T
    feats[765] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. du.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'acc. du.', node2.tup)
            
    # L->-48->T
    feats[766] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-48') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-48', node2.tup)
            
    # L->-126->T
    feats[767] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-126') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-126', node2.tup)
            
    # L->10_du->T
    feats[768] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_du') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '10_du', node2.tup)
            
    # L->-61->T
    feats[769] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-61') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-61', node2.tup)
            
    # L->16_tp->T
    feats[770] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_tp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '16_tp', node2.tup)
            
    # L->-90->T
    feats[771] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-90') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-90', node2.tup)
            
    # L->abl. sg.->T
    feats[772] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'abl. sg.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'abl. sg.', node2.tup)
            
    # L->acc. sg.->T
    feats[773] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. sg.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'acc. sg.', node2.tup)
            
    # L->4_fp->T
    feats[774] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_fp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '4_fp', node2.tup)
            
    # L->-71->T
    feats[775] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-71') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-71', node2.tup)
            
    # L->14_sg->T
    feats[776] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_sg') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '14_sg', node2.tup)
            
    # L->-42->T
    feats[777] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-42') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-42', node2.tup)
            
    # L->6_sp->T
    feats[778] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '6_sp', node2.tup)
            
    # L->3_sp->T
    feats[779] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_sp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '3_sp', node2.tup)
            
    # L->51->T
    feats[780] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '51') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '51', node2.tup)
            
    # L->dat. pl.->T
    feats[781] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. pl.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'dat. pl.', node2.tup)
            
    # L->du_sp->T
    feats[782] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'du_sp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'du_sp', node2.tup)
            
    # L->179->T
    feats[783] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '179') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '179', node2.tup)
            
    # L->pl_tp->T
    feats[784] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_tp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'pl_tp', node2.tup)
            
    # L->79->T
    feats[785] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '79') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '79', node2.tup)
            
    # L->30_tp->T
    feats[786] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_tp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '30_tp', node2.tup)
            
    # L->-133->T
    feats[787] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-133') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-133', node2.tup)
            
    # L->5_pl->T
    feats[788] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_pl') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '5_pl', node2.tup)
            
    # L->-144->T
    feats[789] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-144') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-144', node2.tup)
            
    # L->7_sg->T
    feats[790] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_sg') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '7_sg', node2.tup)
            
    # L->29->T
    feats[791] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '29') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '29', node2.tup)
            
    # L->30_pl->T
    feats[792] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_pl') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '30_pl', node2.tup)
            
    # L->149->T
    feats[793] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '149') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '149', node2.tup)
            
    # L->182->T
    feats[794] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '182') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '182', node2.tup)
            
    # L->-292->T
    feats[795] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-292') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-292', node2.tup)
            
    # L->6_pl->T
    feats[796] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_pl') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '6_pl', node2.tup)
            
    # L->-269->T
    feats[797] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-269') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-269', node2.tup)
            
    # L->72->T
    feats[798] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '72') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '72', node2.tup)
            
    # L->instr. du.->T
    feats[799] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. du.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'instr. du.', node2.tup)
            
    # L->sp->T
    feats[800] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'sp', node2.tup)
            
    # L->173->T
    feats[801] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '173') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '173', node2.tup)
            
    # L->29_sg->T
    feats[802] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '29_sg') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '29_sg', node2.tup)
            
    # L->7_fp->T
    feats[803] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_fp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '7_fp', node2.tup)
            
    # L->-114->T
    feats[804] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-114') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-114', node2.tup)
            
    # L->80->T
    feats[805] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '80') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '80', node2.tup)
            
    # L->-119->T
    feats[806] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-119') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-119', node2.tup)
            
    # L->-279->T
    feats[807] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-279') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-279', node2.tup)
            
    # L->-36->T
    feats[808] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-36') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-36', node2.tup)
            
    # L->-34->T
    feats[809] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-34', node2.tup)
            
    # L->9_du->T
    feats[810] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_du') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '9_du', node2.tup)
            
    # L->155->T
    feats[811] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '155') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '155', node2.tup)
            
    # L->168->T
    feats[812] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '168') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '168', node2.tup)
            
    # L->-307->T
    feats[813] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-307') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-307', node2.tup)
            
    # L->110->T
    feats[814] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '110') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '110', node2.tup)
            
    # L->69->T
    feats[815] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '69') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '69', node2.tup)
            
    # L->3_sg->T
    feats[816] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_sg') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '3_sg', node2.tup)
            
    # L->28_tp->T
    feats[817] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_tp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '28_tp', node2.tup)
            
    # L->4_pl->T
    feats[818] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_pl') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '4_pl', node2.tup)
            
    # L->14_fp->T
    feats[819] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_fp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '14_fp', node2.tup)
            
    # L->-64->T
    feats[820] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-64') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-64', node2.tup)
            
    # L->-97->T
    feats[821] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-97') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-97', node2.tup)
            
    # L->77->T
    feats[822] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '77') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '77', node2.tup)
            
    # L->148->T
    feats[823] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '148') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '148', node2.tup)
            
    # L->-41->T
    feats[824] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-41') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-41', node2.tup)
            
    # L->pl->T
    feats[825] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'pl', node2.tup)
            
    # L->122->T
    feats[826] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '122') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '122', node2.tup)
            
    # L->-117->T
    feats[827] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-117') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-117', node2.tup)
            
    # L->-141->T
    feats[828] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-141') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-141', node2.tup)
            
    # L->-142->T
    feats[829] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-142') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-142', node2.tup)
            
    # L->134->T
    feats[830] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '134') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '134', node2.tup)
            
    # L->2_sp->T
    feats[831] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_sp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '2_sp', node2.tup)
            
    # L->voc. du.->T
    feats[832] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. du.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'voc. du.', node2.tup)
            
    # L->12_pl->T
    feats[833] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_pl') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '12_pl', node2.tup)
            
    # L->-301->T
    feats[834] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-301') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-301', node2.tup)
            
    # L->-262->T
    feats[835] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-262') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-262', node2.tup)
            
    # L->-150->T
    feats[836] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-150') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-150', node2.tup)
            
    # L->6_tp->T
    feats[837] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_tp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '6_tp', node2.tup)
            
    # L->instr. neutr.->T
    feats[838] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. neutr.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'instr. neutr.', node2.tup)
            
    # L->2_pl->T
    feats[839] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_pl') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '2_pl', node2.tup)
            
    # L->15_tp->T
    feats[840] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_tp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '15_tp', node2.tup)
            
    # L->-93->T
    feats[841] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-93') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-93', node2.tup)
            
    # L->neutr->T
    feats[842] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'neutr') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'neutr', node2.tup)
            
    # L->73->T
    feats[843] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '73') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '73', node2.tup)
            
    # L->voc. neutr.->T
    feats[844] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. neutr.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'voc. neutr.', node2.tup)
            
    # L->-246->T
    feats[845] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-246') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-246', node2.tup)
            
    # L->3_fp->T
    feats[846] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_fp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '3_fp', node2.tup)
            
    # L->88->T
    feats[847] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '88') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '88', node2.tup)
            
    # L->-49->T
    feats[848] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-49') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-49', node2.tup)
            
    # L->96->T
    feats[849] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '96') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '96', node2.tup)
            
    # L->59->T
    feats[850] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '59') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '59', node2.tup)
            
    # L->loc->T
    feats[851] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'loc', node2.tup)
            
    # L->-153->T
    feats[852] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-153') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-153', node2.tup)
            
    # L->-23->T
    feats[853] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-23') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-23', node2.tup)
            
    # L->-14->T
    feats[854] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-14') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-14', node2.tup)
            
    # L->-45->T
    feats[855] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-45') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-45', node2.tup)
            
    # L->nom. neutr.->T
    feats[856] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. neutr.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'nom. neutr.', node2.tup)
            
    # L->2_sg->T
    feats[857] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_sg') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '2_sg', node2.tup)
            
    # L->-159->T
    feats[858] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-159') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-159', node2.tup)
            
    # L->-309->T
    feats[859] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-309') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-309', node2.tup)
            
    # L->158->T
    feats[860] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '158') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '158', node2.tup)
            
    # L->90->T
    feats[861] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '90') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '90', node2.tup)
            
    # L->39->T
    feats[862] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '39', node2.tup)
            
    # L->137->T
    feats[863] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '137') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '137', node2.tup)
            
    # L->-10->T
    feats[864] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-10') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-10', node2.tup)
            
    # L->12_fp->T
    feats[865] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_fp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '12_fp', node2.tup)
            
    # L->-26->T
    feats[866] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-26') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-26', node2.tup)
            
    # L->132->T
    feats[867] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '132') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '132', node2.tup)
            
    # L->gen->T
    feats[868] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'gen', node2.tup)
            
    # L->instr->T
    feats[869] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'instr', node2.tup)
            
    # L->-59->T
    feats[870] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-59') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-59', node2.tup)
            
    # L->10_fp->T
    feats[871] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_fp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '10_fp', node2.tup)
            
    # L->3_tp->T
    feats[872] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_tp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '3_tp', node2.tup)
            
    # L->114->T
    feats[873] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '114') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '114', node2.tup)
            
    # L->8_sg->T
    feats[874] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_sg') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '8_sg', node2.tup)
            
    # L->-76->T
    feats[875] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-76') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-76', node2.tup)
            
    # L->-24->T
    feats[876] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-24') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-24', node2.tup)
            
    # L->2->T
    feats[877] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '2', node2.tup)
            
    # L->6_sg->T
    feats[878] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sg') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '6_sg', node2.tup)
            
    # L->masc->T
    feats[879] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'masc') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'masc', node2.tup)
            
    # L->-96->T
    feats[880] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-96') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-96', node2.tup)
            
    # L->-82->T
    feats[881] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-82') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-82', node2.tup)
            
    # L->9_tp->T
    feats[882] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_tp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '9_tp', node2.tup)
            
    # L->-109->T
    feats[883] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-109') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-109', node2.tup)
            
    # L->abl->T
    feats[884] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'abl') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'abl', node2.tup)
            
    # L->4_du->T
    feats[885] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_du') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '4_du', node2.tup)
            
    # L->abl. du.->T
    feats[886] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'abl. du.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'abl. du.', node2.tup)
            
    # L->4_sg->T
    feats[887] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_sg') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '4_sg', node2.tup)
            
    # L->60->T
    feats[888] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '60') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '60', node2.tup)
            
    # L->50->T
    feats[889] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '50') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '50', node2.tup)
            
    # L->-91->T
    feats[890] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-91') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-91', node2.tup)
            
    # L->-263->T
    feats[891] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-263') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-263', node2.tup)
            
    # L->voc. pl.->T
    feats[892] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. pl.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'voc. pl.', node2.tup)
            
    # L->-31->T
    feats[893] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-31') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-31', node2.tup)
            
    # L->6_du->T
    feats[894] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_du') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '6_du', node2.tup)
            
    # L->30->T
    feats[895] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '30', node2.tup)
            
    # L->41->T
    feats[896] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '41') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '41', node2.tup)
            
    # L->91->T
    feats[897] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '91') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '91', node2.tup)
            
    # L->94->T
    feats[898] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '94') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '94', node2.tup)
            
    # L->-69->T
    feats[899] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-69') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-69', node2.tup)
            
    # L->131->T
    feats[900] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '131') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '131', node2.tup)
            
    # L->10_tp->T
    feats[901] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_tp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '10_tp', node2.tup)
            
    # L->140->T
    feats[902] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '140') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '140', node2.tup)
            
    # L->tp->T
    feats[903] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'tp', node2.tup)
            
    # L->159->T
    feats[904] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '159') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '159', node2.tup)
            
    # L->gen. du.->T
    feats[905] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. du.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'gen. du.', node2.tup)
            
    # L->-87->T
    feats[906] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-87') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-87', node2.tup)
            
    # L->150->T
    feats[907] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '150') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '150', node2.tup)
            
    # L->8_tp->T
    feats[908] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_tp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '8_tp', node2.tup)
            
    # L->3->T
    feats[909] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '3', node2.tup)
            
    # L->-308->T
    feats[910] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-308') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-308', node2.tup)
            
    # L->33->T
    feats[911] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '33') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '33', node2.tup)
            
    # L->acc. fem->T
    feats[912] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. fem') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'acc. fem', node2.tup)
            
    # L->9_sg->T
    feats[913] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sg') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '9_sg', node2.tup)
            
    # L->129->T
    feats[914] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '129') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '129', node2.tup)
            
    # L->-94->T
    feats[915] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-94') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-94', node2.tup)
            
    # L->nom. sg.->T
    feats[916] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. sg.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'nom. sg.', node2.tup)
            
    # L->142->T
    feats[917] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '142') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '142', node2.tup)
            
    # L->nom. fem->T
    feats[918] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. fem') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'nom. fem', node2.tup)
            
    # L->-276->T
    feats[919] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-276') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-276', node2.tup)
            
    # L->111->T
    feats[920] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '111') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '111', node2.tup)
            
    # L->fp->T
    feats[921] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'fp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'fp', node2.tup)
            
    # L->nom. masc.->T
    feats[922] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. masc.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'nom. masc.', node2.tup)
            
    # L->-271->T
    feats[923] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-271') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-271', node2.tup)
            
    # L->5_fp->T
    feats[924] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_fp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '5_fp', node2.tup)
            
    # L->174->T
    feats[925] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '174') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '174', node2.tup)
            
    # L->99->T
    feats[926] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '99') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '99', node2.tup)
            
    # L->-30->T
    feats[927] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-30', node2.tup)
            
    # L->voc. fem->T
    feats[928] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. fem') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'voc. fem', node2.tup)
            
    # L->-73->T
    feats[929] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-73') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-73', node2.tup)
            
    # L->nom->T
    feats[930] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'nom', node2.tup)
            
    # L->78->T
    feats[931] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '78') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '78', node2.tup)
            
    # L->-13->T
    feats[932] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-13') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-13', node2.tup)
            
    # L->115->T
    feats[933] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '115') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '115', node2.tup)
            
    # L->15_pl->T
    feats[934] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_pl') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '15_pl', node2.tup)
            
    # L->70->T
    feats[935] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '70') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '70', node2.tup)
            
    # L->35->T
    feats[936] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '35') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '35', node2.tup)
            
    # L->-247->T
    feats[937] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-247') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-247', node2.tup)
            
    # L->157->T
    feats[938] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '157') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '157', node2.tup)
            
    # L->1->T
    feats[939] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '1') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '1', node2.tup)
            
    # L->-249->T
    feats[940] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-249') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-249', node2.tup)
            
    # L->-38->T
    feats[941] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-38') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-38', node2.tup)
            
    # L->161->T
    feats[942] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '161') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '161', node2.tup)
            
    # L->27_sg->T
    feats[943] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_sg') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '27_sg', node2.tup)
            
    # L->-55->T
    feats[944] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-55') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-55', node2.tup)
            
    # L->acc->T
    feats[945] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'acc', node2.tup)
            
    # L->-51->T
    feats[946] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-51') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-51', node2.tup)
            
    # L->-16->T
    feats[947] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-16') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-16', node2.tup)
            
    # L->-52->T
    feats[948] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-52') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-52', node2.tup)
            
    # L->loc. pl.->T
    feats[949] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. pl.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'loc. pl.', node2.tup)
            
    # L->42->T
    feats[950] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '42') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '42', node2.tup)
            
    # L->-84->T
    feats[951] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-84') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-84', node2.tup)
            
    # L->27_fp->T
    feats[952] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_fp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '27_fp', node2.tup)
            
    # L->-20->T
    feats[953] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-20') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-20', node2.tup)
            
    # L->141->T
    feats[954] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '141') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '141', node2.tup)
            
    # L->139->T
    feats[955] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '139') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '139', node2.tup)
            
    # L->-129->T
    feats[956] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-129') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-129', node2.tup)
            
    # L->120->T
    feats[957] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '120') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '120', node2.tup)
            
    # L->-39->T
    feats[958] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-39') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-39', node2.tup)
            
    # L->-166->T
    feats[959] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-166') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-166', node2.tup)
            
    # L->-28->T
    feats[960] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-28') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-28', node2.tup)
            
    # L->170->T
    feats[961] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '170') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '170', node2.tup)
            
    # L->-131->T
    feats[962] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-131') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-131', node2.tup)
            
    # L->-19->T
    feats[963] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-19') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-19', node2.tup)
            
    # L->2_du->T
    feats[964] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '2_du', node2.tup)
            
    # L->-112->T
    feats[965] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-112') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-112', node2.tup)
            
    # L->-104->T
    feats[966] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-104', node2.tup)
            
    # L->-12->T
    feats[967] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-12') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-12', node2.tup)
            
    # L->acc. neutr.->T
    feats[968] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. neutr.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'acc. neutr.', node2.tup)
            
    # L->14_tp->T
    feats[969] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_tp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '14_tp', node2.tup)
            
    # L->-29->T
    feats[970] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-29') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-29', node2.tup)
            
    # L->sg->T
    feats[971] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'sg', node2.tup)
            
    # L->177->T
    feats[972] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '177') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '177', node2.tup)
            
    # L->7_tp->T
    feats[973] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_tp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '7_tp', node2.tup)
            
    # L->-81->T
    feats[974] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-81') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-81', node2.tup)
            
    # L->-92->T
    feats[975] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-92') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-92', node2.tup)
            
    # L->-230->T
    feats[976] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-230', node2.tup)
            
    # L->-297->T
    feats[977] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-297') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-297', node2.tup)
            
    # L->voc. masc.->T
    feats[978] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. masc.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'voc. masc.', node2.tup)
            
    # L->56->T
    feats[979] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '56') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '56', node2.tup)
            
    # L->gen. sg.->T
    feats[980] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. sg.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'gen. sg.', node2.tup)
            
    # L->172->T
    feats[981] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '172') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '172', node2.tup)
            
    # L->11_pl->T
    feats[982] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '11_pl') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '11_pl', node2.tup)
            
    # L->-11->T
    feats[983] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-11') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-11', node2.tup)
            
    # L->13_tp->T
    feats[984] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_tp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '13_tp', node2.tup)
            
    # L->160->T
    feats[985] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '160') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '160', node2.tup)
            
    # L->-115->T
    feats[986] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-115') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-115', node2.tup)
            
    # L->-50->T
    feats[987] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-50') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-50', node2.tup)
            
    # L->6_fp->T
    feats[988] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_fp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '6_fp', node2.tup)
            
    # L->68->T
    feats[989] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '68') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '68', node2.tup)
            
    # C->121->L
    feats[990] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '121') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '121', node2.lemma)
            
    # C->-15->L
    feats[991] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-15') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-15', node2.lemma)
            
    # C->-121->L
    feats[992] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-121') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-121', node2.lemma)
            
    # C->117->L
    feats[993] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '117') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '117', node2.lemma)
            
    # C->10_sg->L
    feats[994] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '10_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_sg', node2.lemma)
            
    # C->-240->L
    feats[995] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-240') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-240', node2.lemma)
            
    # C->27_pl->L
    feats[996] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '27_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_pl', node2.lemma)
            
    # C->55->L
    feats[997] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '55') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '55', node2.lemma)
            
    # C->-245->L
    feats[998] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-245') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-245', node2.lemma)
            
    # C->171->L
    feats[999] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '171') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '171', node2.lemma)
            
    # C->154->L
    feats[1000] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '154') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '154', node2.lemma)
            
    # C->14_sp->L
    feats[1001] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '14_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_sp', node2.lemma)
            
    # C->38->L
    feats[1002] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '38') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '38', node2.lemma)
            
    # C->-151->L
    feats[1003] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-151') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-151', node2.lemma)
            
    # C->du_tp->L
    feats[1004] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'du_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_tp', node2.lemma)
            
    # C->-299->L
    feats[1005] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-299') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-299', node2.lemma)
            
    # C->40->L
    feats[1006] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '40') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '40', node2.lemma)
            
    # C->61->L
    feats[1007] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '61') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '61', node2.lemma)
            
    # C->-79->L
    feats[1008] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-79') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-79', node2.lemma)
            
    # C->9_fp->L
    feats[1009] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '9_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_fp', node2.lemma)
            
    # C->12_tp->L
    feats[1010] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '12_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_tp', node2.lemma)
            
    # C->101->L
    feats[1011] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '101') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '101', node2.lemma)
            
    # C->29_tp->L
    feats[1012] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '29_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '29_tp', node2.lemma)
            
    # C->pl_sp->L
    feats[1013] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'pl_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_sp', node2.lemma)
            
    # C->135->L
    feats[1014] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '135') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '135', node2.lemma)
            
    # C->-132->L
    feats[1015] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-132') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-132', node2.lemma)
            
    # C->-122->L
    feats[1016] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-122') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-122', node2.lemma)
            
    # C->acc. adj.->L
    feats[1017] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'acc. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. adj.', node2.lemma)
            
    # C->-123->L
    feats[1018] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-123') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-123', node2.lemma)
            
    # C->102->L
    feats[1019] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '102') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '102', node2.lemma)
            
    # C->-152->L
    feats[1020] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-152') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-152', node2.lemma)
            
    # C->3_pl->L
    feats[1021] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '3_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_pl', node2.lemma)
            
    # C->instr. sg.->L
    feats[1022] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'instr. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. sg.', node2.lemma)
            
    # C->-66->L
    feats[1023] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-66') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-66', node2.lemma)
            
    # C->151->L
    feats[1024] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '151') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '151', node2.lemma)
            
    # C->15_du->L
    feats[1025] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '15_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_du', node2.lemma)
            
    # C->82->L
    feats[1026] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '82') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '82', node2.lemma)
            
    # C->instr. fem->L
    feats[1027] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'instr. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. fem', node2.lemma)
            
    # C->-98->L
    feats[1028] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-98') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-98', node2.lemma)
            
    # C->-242->L
    feats[1029] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-242') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-242', node2.lemma)
            
    # C->27_du->L
    feats[1030] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '27_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_du', node2.lemma)
            
    # C->-156->L
    feats[1031] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-156') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-156', node2.lemma)
            
    # C->92->L
    feats[1032] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '92') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '92', node2.lemma)
            
    # C->gen. pl.->L
    feats[1033] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'gen. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen. pl.', node2.lemma)
            
    # C->74->L
    feats[1034] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '74') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '74', node2.lemma)
            
    # C->-86->L
    feats[1035] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-86') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-86', node2.lemma)
            
    # C->8_du->L
    feats[1036] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '8_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_du', node2.lemma)
            
    # C->97->L
    feats[1037] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '97') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '97', node2.lemma)
            
    # C->15_sg->L
    feats[1038] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '15_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_sg', node2.lemma)
            
    # C->7_du->L
    feats[1039] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '7_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_du', node2.lemma)
            
    # C->11_sp->L
    feats[1040] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '11_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_sp', node2.lemma)
            
    # C->89->L
    feats[1041] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '89') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '89', node2.lemma)
            
    # C->-261->L
    feats[1042] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-261') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-261', node2.lemma)
            
    # C->16_fp->L
    feats[1043] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '16_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_fp', node2.lemma)
            
    # C->178->L
    feats[1044] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '178') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '178', node2.lemma)
            
    # C->12_sg->L
    feats[1045] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '12_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_sg', node2.lemma)
            
    # C->sg_fp->L
    feats[1046] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'sg_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_fp', node2.lemma)
            
    # C->-200->L
    feats[1047] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-200') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-200', node2.lemma)
            
    # C->4_sp->L
    feats[1048] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '4_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_sp', node2.lemma)
            
    # C->adj->L
    feats[1049] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'adj') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'adj', node2.lemma)
            
    # C->10_sp->L
    feats[1050] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '10_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_sp', node2.lemma)
            
    # C->-22->L
    feats[1051] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-22') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-22', node2.lemma)
            
    # C->-46->L
    feats[1052] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-46') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-46', node2.lemma)
            
    # C->30_fp->L
    feats[1053] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '30_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_fp', node2.lemma)
            
    # C->instr. pl.->L
    feats[1054] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'instr. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. pl.', node2.lemma)
            
    # C->nom. adj.->L
    feats[1055] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'nom. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. adj.', node2.lemma)
            
    # C->-27->L
    feats[1056] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-27') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-27', node2.lemma)
            
    # C->-35->L
    feats[1057] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-35') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
            
    # C->-243->L
    feats[1058] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-243') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-243', node2.lemma)
            
    # C->-18->L
    feats[1059] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-18') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-18', node2.lemma)
            
    # C->voc. sg.->L
    feats[1060] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'voc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. sg.', node2.lemma)
            
    # C->49->L
    feats[1061] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '49') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '49', node2.lemma)
            
    # C->-67->L
    feats[1062] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-67') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-67', node2.lemma)
            
    # C->-89->L
    feats[1063] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-89') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-89', node2.lemma)
            
    # C->-302->L
    feats[1064] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-302') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-302', node2.lemma)
            
    # C->-21->L
    feats[1065] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-21') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-21', node2.lemma)
            
    # C->11_tp->L
    feats[1066] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '11_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_tp', node2.lemma)
            
    # C->-139->L
    feats[1067] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-139') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-139', node2.lemma)
            
    # C->28->L
    feats[1068] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '28') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28', node2.lemma)
            
    # C->-56->L
    feats[1069] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-56') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-56', node2.lemma)
            
    # C->acc. masc.->L
    feats[1070] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'acc. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. masc.', node2.lemma)
            
    # C->-47->L
    feats[1071] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-47') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-47', node2.lemma)
            
    # C->-17->L
    feats[1072] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-17') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-17', node2.lemma)
            
    # C->5_du->L
    feats[1073] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '5_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_du', node2.lemma)
            
    # C->98->L
    feats[1074] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '98') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '98', node2.lemma)
            
    # C->81->L
    feats[1075] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '81') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '81', node2.lemma)
            
    # C->-113->L
    feats[1076] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-113') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-113', node2.lemma)
            
    # C->sg_tp->L
    feats[1077] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'sg_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_tp', node2.lemma)
            
    # C->-169->L
    feats[1078] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-169') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-169', node2.lemma)
            
    # C->-283->L
    feats[1079] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-283') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-283', node2.lemma)
            
    # C->-210->L
    feats[1080] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-210') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-210', node2.lemma)
            
    # C->-25->L
    feats[1081] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-25') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-25', node2.lemma)
            
    # C->162->L
    feats[1082] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '162') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '162', node2.lemma)
            
    # C->-101->L
    feats[1083] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-101') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-101', node2.lemma)
            
    # C->instr. adj.->L
    feats[1084] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'instr. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. adj.', node2.lemma)
            
    # C->du_fp->L
    feats[1085] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'du_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_fp', node2.lemma)
            
    # C->-303->L
    feats[1086] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-303') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-303', node2.lemma)
            
    # C->5_tp->L
    feats[1087] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '5_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_tp', node2.lemma)
            
    # C->175->L
    feats[1088] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '175') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '175', node2.lemma)
            
    # C->11_du->L
    feats[1089] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '11_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_du', node2.lemma)
            
    # C->15_sp->L
    feats[1090] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '15_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_sp', node2.lemma)
            
    # C->30_du->L
    feats[1091] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '30_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_du', node2.lemma)
            
    # C->128->L
    feats[1092] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '128') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '128', node2.lemma)
            
    # C->-273->L
    feats[1093] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-273') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-273', node2.lemma)
            
    # C->7_pl->L
    feats[1094] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '7_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_pl', node2.lemma)
            
    # C->-137->L
    feats[1095] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-137') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-137', node2.lemma)
            
    # C->9_pl->L
    feats[1096] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '9_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_pl', node2.lemma)
            
    # C->acc. pl.->L
    feats[1097] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'acc. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. pl.', node2.lemma)
            
    # C->-103->L
    feats[1098] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-103') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-103', node2.lemma)
            
    # C->71->L
    feats[1099] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '71') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '71', node2.lemma)
            
    # C->136->L
    feats[1100] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '136') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '136', node2.lemma)
            
    # C->sg_sp->L
    feats[1101] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'sg_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_sp', node2.lemma)
            
    # C->-306->L
    feats[1102] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-306') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-306', node2.lemma)
            
    # C->dat. du.->L
    feats[1103] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'dat. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
            
    # C->109->L
    feats[1104] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '109') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '109', node2.lemma)
            
    # C->-291->L
    feats[1105] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-291') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-291', node2.lemma)
            
    # C->-296->L
    feats[1106] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-296') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-296', node2.lemma)
            
    # C->11_sg->L
    feats[1107] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '11_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_sg', node2.lemma)
            
    # C->-57->L
    feats[1108] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-57') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-57', node2.lemma)
            
    # C->-77->L
    feats[1109] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-77') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-77', node2.lemma)
            
    # C->-161->L
    feats[1110] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-161') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-161', node2.lemma)
            
    # C->13_fp->L
    feats[1111] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '13_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_fp', node2.lemma)
            
    # C->-266->L
    feats[1112] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-266') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-266', node2.lemma)
            
    # C->30_sg->L
    feats[1113] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '30_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_sg', node2.lemma)
            
    # C->5_sp->L
    feats[1114] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '5_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sp', node2.lemma)
            
    # C->181->L
    feats[1115] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '181') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '181', node2.lemma)
            
    # C->-111->L
    feats[1116] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-111') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-111', node2.lemma)
            
    # C->3_du->L
    feats[1117] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '3_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_du', node2.lemma)
            
    # C->119->L
    feats[1118] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '119') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '119', node2.lemma)
            
    # C->4_tp->L
    feats[1119] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '4_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_tp', node2.lemma)
            
    # C->-33->L
    feats[1120] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-33') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-33', node2.lemma)
            
    # C->-37->L
    feats[1121] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-37') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-37', node2.lemma)
            
    # C->112->L
    feats[1122] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '112') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '112', node2.lemma)
            
    # C->pl_fp->L
    feats[1123] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'pl_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_fp', node2.lemma)
            
    # C->28_sg->L
    feats[1124] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '28_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28_sg', node2.lemma)
            
    # C->138->L
    feats[1125] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '138') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '138', node2.lemma)
            
    # C->fem->L
    feats[1126] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'fem', node2.lemma)
            
    # C->-78->L
    feats[1127] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-78') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-78', node2.lemma)
            
    # C->-54->L
    feats[1128] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-54') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-54', node2.lemma)
            
    # C->176->L
    feats[1129] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '176') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '176', node2.lemma)
            
    # C->2_tp->L
    feats[1130] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '2_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_tp', node2.lemma)
            
    # C->100->L
    feats[1131] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '100') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '100', node2.lemma)
            
    # C->-99->L
    feats[1132] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-99') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-99', node2.lemma)
            
    # C->-163->L
    feats[1133] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-163') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-163', node2.lemma)
            
    # C->14_pl->L
    feats[1134] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '14_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_pl', node2.lemma)
            
    # C->11_fp->L
    feats[1135] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '11_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_fp', node2.lemma)
            
    # C->37->L
    feats[1136] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '37') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '37', node2.lemma)
            
    # C->-44->L
    feats[1137] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-44') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-44', node2.lemma)
            
    # C->95->L
    feats[1138] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '95') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '95', node2.lemma)
            
    # C->-83->L
    feats[1139] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-83') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-83', node2.lemma)
            
    # C->76->L
    feats[1140] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '76') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '76', node2.lemma)
            
    # C->36->L
    feats[1141] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '36') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '36', node2.lemma)
            
    # C->-68->L
    feats[1142] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-68') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-68', node2.lemma)
            
    # C->du->L
    feats[1143] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du', node2.lemma)
            
    # C->voc->L
    feats[1144] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'voc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc', node2.lemma)
            
    # C->-53->L
    feats[1145] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-53') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-53', node2.lemma)
            
    # C->48->L
    feats[1146] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '48') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '48', node2.lemma)
            
    # C->156->L
    feats[1147] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '156') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '156', node2.lemma)
            
    # C->-72->L
    feats[1148] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-72') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-72', node2.lemma)
            
    # C->-158->L
    feats[1149] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-158') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-158', node2.lemma)
            
    # C->-293->L
    feats[1150] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-293') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-293', node2.lemma)
            
    # C->34->L
    feats[1151] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '34') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '34', node2.lemma)
            
    # C->instr. masc.->L
    feats[1152] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'instr. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. masc.', node2.lemma)
            
    # C->7_sp->L
    feats[1153] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '7_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_sp', node2.lemma)
            
    # C->2_fp->L
    feats[1154] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '2_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_fp', node2.lemma)
            
    # C->152->L
    feats[1155] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '152') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '152', node2.lemma)
            
    # C->-43->L
    feats[1156] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-43') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-43', node2.lemma)
            
    # C->-149->L
    feats[1157] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-149') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-149', node2.lemma)
            
    # C->nom. du.->L
    feats[1158] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'nom. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. du.', node2.lemma)
            
    # C->30_sp->L
    feats[1159] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '30_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_sp', node2.lemma)
            
    # C->-143->L
    feats[1160] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-143') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-143', node2.lemma)
            
    # C->180->L
    feats[1161] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '180') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '180', node2.lemma)
            
    # C->16_sg->L
    feats[1162] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '16_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_sg', node2.lemma)
            
    # C->108->L
    feats[1163] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '108') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '108', node2.lemma)
            
    # C->16_du->L
    feats[1164] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '16_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_du', node2.lemma)
            
    # C->5_sg->L
    feats[1165] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '5_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sg', node2.lemma)
            
    # C->32->L
    feats[1166] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '32') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '32', node2.lemma)
            
    # C->27_tp->L
    feats[1167] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '27_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_tp', node2.lemma)
            
    # C->54->L
    feats[1168] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '54') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '54', node2.lemma)
            
    # C->acc. du.->L
    feats[1169] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'acc. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. du.', node2.lemma)
            
    # C->-48->L
    feats[1170] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-48') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-48', node2.lemma)
            
    # C->-126->L
    feats[1171] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-126') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-126', node2.lemma)
            
    # C->10_du->L
    feats[1172] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '10_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_du', node2.lemma)
            
    # C->-61->L
    feats[1173] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-61') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-61', node2.lemma)
            
    # C->-102->L
    feats[1174] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-102') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-102', node2.lemma)
            
    # C->16_tp->L
    feats[1175] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '16_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_tp', node2.lemma)
            
    # C->-90->L
    feats[1176] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-90') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-90', node2.lemma)
            
    # C->abl. sg.->L
    feats[1177] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'abl. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. sg.', node2.lemma)
            
    # C->acc. sg.->L
    feats[1178] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'acc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. sg.', node2.lemma)
            
    # C->4_fp->L
    feats[1179] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '4_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_fp', node2.lemma)
            
    # C->-71->L
    feats[1180] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-71') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-71', node2.lemma)
            
    # C->14_sg->L
    feats[1181] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '14_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_sg', node2.lemma)
            
    # C->-42->L
    feats[1182] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-42') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-42', node2.lemma)
            
    # C->6_sp->L
    feats[1183] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '6_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_sp', node2.lemma)
            
    # C->3_sp->L
    feats[1184] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '3_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_sp', node2.lemma)
            
    # C->51->L
    feats[1185] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '51') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '51', node2.lemma)
            
    # C->dat. pl.->L
    feats[1186] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'dat. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. pl.', node2.lemma)
            
    # C->du_sp->L
    feats[1187] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'du_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du_sp', node2.lemma)
            
    # C->179->L
    feats[1188] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '179') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '179', node2.lemma)
            
    # C->pl_tp->L
    feats[1189] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'pl_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_tp', node2.lemma)
            
    # C->79->L
    feats[1190] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '79') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '79', node2.lemma)
            
    # C->30_tp->L
    feats[1191] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '30_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_tp', node2.lemma)
            
    # C->-133->L
    feats[1192] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-133') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-133', node2.lemma)
            
    # C->5_pl->L
    feats[1193] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '5_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_pl', node2.lemma)
            
    # C->-144->L
    feats[1194] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-144') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-144', node2.lemma)
            
    # C->7_sg->L
    feats[1195] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '7_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_sg', node2.lemma)
            
    # C->29->L
    feats[1196] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '29') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '29', node2.lemma)
            
    # C->153->L
    feats[1197] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '153') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
            
    # C->30_pl->L
    feats[1198] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '30_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_pl', node2.lemma)
            
    # C->149->L
    feats[1199] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '149') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '149', node2.lemma)
            
    # C->182->L
    feats[1200] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '182') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '182', node2.lemma)
            
    # C->6_pl->L
    feats[1201] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '6_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_pl', node2.lemma)
            
    # C->-269->L
    feats[1202] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-269') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-269', node2.lemma)
            
    # C->72->L
    feats[1203] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '72') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '72', node2.lemma)
            
    # C->instr. du.->L
    feats[1204] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'instr. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. du.', node2.lemma)
            
    # C->sp->L
    feats[1205] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sp', node2.lemma)
            
    # C->173->L
    feats[1206] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '173') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '173', node2.lemma)
            
    # C->29_sg->L
    feats[1207] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '29_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '29_sg', node2.lemma)
            
    # C->7_fp->L
    feats[1208] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '7_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_fp', node2.lemma)
            
    # C->-114->L
    feats[1209] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-114') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-114', node2.lemma)
            
    # C->-119->L
    feats[1210] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-119') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-119', node2.lemma)
            
    # C->-32->L
    feats[1211] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-32') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-32', node2.lemma)
            
    # C->-279->L
    feats[1212] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-279') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-279', node2.lemma)
            
    # C->-36->L
    feats[1213] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-36') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-36', node2.lemma)
            
    # C->-34->L
    feats[1214] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-34') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-34', node2.lemma)
            
    # C->9_du->L
    feats[1215] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '9_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_du', node2.lemma)
            
    # C->155->L
    feats[1216] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '155') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '155', node2.lemma)
            
    # C->168->L
    feats[1217] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '168') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '168', node2.lemma)
            
    # C->-307->L
    feats[1218] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-307') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-307', node2.lemma)
            
    # C->110->L
    feats[1219] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '110') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '110', node2.lemma)
            
    # C->69->L
    feats[1220] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '69') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '69', node2.lemma)
            
    # C->3_sg->L
    feats[1221] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '3_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_sg', node2.lemma)
            
    # C->28_tp->L
    feats[1222] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '28_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28_tp', node2.lemma)
            
    # C->4_pl->L
    feats[1223] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '4_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_pl', node2.lemma)
            
    # C->14_fp->L
    feats[1224] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '14_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_fp', node2.lemma)
            
    # C->-64->L
    feats[1225] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-64') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-64', node2.lemma)
            
    # C->-97->L
    feats[1226] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-97') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-97', node2.lemma)
            
    # C->75->L
    feats[1227] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '75') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '75', node2.lemma)
            
    # C->77->L
    feats[1228] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '77') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '77', node2.lemma)
            
    # C->148->L
    feats[1229] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '148') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '148', node2.lemma)
            
    # C->-41->L
    feats[1230] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-41') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-41', node2.lemma)
            
    # C->pl->L
    feats[1231] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
            
    # C->122->L
    feats[1232] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '122') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '122', node2.lemma)
            
    # C->-117->L
    feats[1233] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-117') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-117', node2.lemma)
            
    # C->-141->L
    feats[1234] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-141') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-141', node2.lemma)
            
    # C->-142->L
    feats[1235] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-142') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-142', node2.lemma)
            
    # C->134->L
    feats[1236] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '134') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '134', node2.lemma)
            
    # C->2_sp->L
    feats[1237] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '2_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_sp', node2.lemma)
            
    # C->voc. du.->L
    feats[1238] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'voc. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. du.', node2.lemma)
            
    # C->12_pl->L
    feats[1239] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '12_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_pl', node2.lemma)
            
    # C->-301->L
    feats[1240] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-301') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-301', node2.lemma)
            
    # C->-262->L
    feats[1241] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-262') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-262', node2.lemma)
            
    # C->-150->L
    feats[1242] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-150') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-150', node2.lemma)
            
    # C->6_tp->L
    feats[1243] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '6_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_tp', node2.lemma)
            
    # C->instr. neutr.->L
    feats[1244] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'instr. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. neutr.', node2.lemma)
            
    # C->2_pl->L
    feats[1245] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '2_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_pl', node2.lemma)
            
    # C->15_tp->L
    feats[1246] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '15_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_tp', node2.lemma)
            
    # C->-93->L
    feats[1247] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-93') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-93', node2.lemma)
            
    # C->neutr->L
    feats[1248] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'neutr') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'neutr', node2.lemma)
            
    # C->73->L
    feats[1249] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '73') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '73', node2.lemma)
            
    # C->voc. neutr.->L
    feats[1250] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'voc. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. neutr.', node2.lemma)
            
    # C->-246->L
    feats[1251] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-246') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-246', node2.lemma)
            
    # C->3_fp->L
    feats[1252] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '3_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_fp', node2.lemma)
            
    # C->88->L
    feats[1253] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '88') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '88', node2.lemma)
            
    # C->-49->L
    feats[1254] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-49') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-49', node2.lemma)
            
    # C->59->L
    feats[1255] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '59') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '59', node2.lemma)
            
    # C->loc->L
    feats[1256] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'loc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc', node2.lemma)
            
    # C->-153->L
    feats[1257] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-153') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-153', node2.lemma)
            
    # C->-23->L
    feats[1258] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-23') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-23', node2.lemma)
            
    # C->-14->L
    feats[1259] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-14') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-14', node2.lemma)
            
    # C->-45->L
    feats[1260] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-45') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-45', node2.lemma)
            
    # C->13_sp->L
    feats[1261] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '13_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_sp', node2.lemma)
            
    # C->2_sg->L
    feats[1262] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '2_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_sg', node2.lemma)
            
    # C->-159->L
    feats[1263] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-159') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-159', node2.lemma)
            
    # C->-309->L
    feats[1264] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-309') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-309', node2.lemma)
            
    # C->158->L
    feats[1265] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '158') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '158', node2.lemma)
            
    # C->90->L
    feats[1266] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '90') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '90', node2.lemma)
            
    # C->39->L
    feats[1267] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '39') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '39', node2.lemma)
            
    # C->-10->L
    feats[1268] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-10') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-10', node2.lemma)
            
    # C->12_fp->L
    feats[1269] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '12_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_fp', node2.lemma)
            
    # C->-26->L
    feats[1270] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-26') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-26', node2.lemma)
            
    # C->132->L
    feats[1271] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '132') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '132', node2.lemma)
            
    # C->gen->L
    feats[1272] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'gen') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen', node2.lemma)
            
    # C->instr->L
    feats[1273] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'instr') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr', node2.lemma)
            
    # C->-59->L
    feats[1274] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-59') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-59', node2.lemma)
            
    # C->10_fp->L
    feats[1275] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '10_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_fp', node2.lemma)
            
    # C->3_tp->L
    feats[1276] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '3_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_tp', node2.lemma)
            
    # C->114->L
    feats[1277] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '114') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '114', node2.lemma)
            
    # C->8_sg->L
    feats[1278] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '8_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_sg', node2.lemma)
            
    # C->-76->L
    feats[1279] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-76') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-76', node2.lemma)
            
    # C->-24->L
    feats[1280] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-24') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
            
    # C->2->L
    feats[1281] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '2') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2', node2.lemma)
            
    # C->6_sg->L
    feats[1282] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '6_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_sg', node2.lemma)
            
    # C->masc->L
    feats[1283] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'masc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'masc', node2.lemma)
            
    # C->-96->L
    feats[1284] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-96') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-96', node2.lemma)
            
    # C->-82->L
    feats[1285] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-82') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-82', node2.lemma)
            
    # C->9_tp->L
    feats[1286] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '9_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_tp', node2.lemma)
            
    # C->-109->L
    feats[1287] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-109') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-109', node2.lemma)
            
    # C->abl->L
    feats[1288] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'abl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl', node2.lemma)
            
    # C->4_du->L
    feats[1289] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '4_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_du', node2.lemma)
            
    # C->abl. du.->L
    feats[1290] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'abl. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. du.', node2.lemma)
            
    # C->4_sg->L
    feats[1291] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '4_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_sg', node2.lemma)
            
    # C->60->L
    feats[1292] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '60') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '60', node2.lemma)
            
    # C->50->L
    feats[1293] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '50') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '50', node2.lemma)
            
    # C->-91->L
    feats[1294] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-91') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-91', node2.lemma)
            
    # C->-263->L
    feats[1295] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-263') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-263', node2.lemma)
            
    # C->voc. pl.->L
    feats[1296] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'voc. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. pl.', node2.lemma)
            
    # C->-31->L
    feats[1297] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-31') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-31', node2.lemma)
            
    # C->abl. pl.->L
    feats[1298] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'abl. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. pl.', node2.lemma)
            
    # C->6_du->L
    feats[1299] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '6_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_du', node2.lemma)
            
    # C->30->L
    feats[1300] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '30') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30', node2.lemma)
            
    # C->41->L
    feats[1301] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '41') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '41', node2.lemma)
            
    # C->91->L
    feats[1302] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '91') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
            
    # C->94->L
    feats[1303] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '94') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '94', node2.lemma)
            
    # C->-69->L
    feats[1304] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-69') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-69', node2.lemma)
            
    # C->131->L
    feats[1305] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '131') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '131', node2.lemma)
            
    # C->10_tp->L
    feats[1306] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '10_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_tp', node2.lemma)
            
    # C->140->L
    feats[1307] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '140') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '140', node2.lemma)
            
    # C->tp->L
    feats[1308] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'tp', node2.lemma)
            
    # C->159->L
    feats[1309] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '159') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '159', node2.lemma)
            
    # C->gen. du.->L
    feats[1310] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'gen. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen. du.', node2.lemma)
            
    # C->-87->L
    feats[1311] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-87') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-87', node2.lemma)
            
    # C->150->L
    feats[1312] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '150') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '150', node2.lemma)
            
    # C->8_tp->L
    feats[1313] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '8_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_tp', node2.lemma)
            
    # C->9_sp->L
    feats[1314] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '9_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_sp', node2.lemma)
            
    # C->3->L
    feats[1315] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '3') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3', node2.lemma)
            
    # C->-308->L
    feats[1316] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-308') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-308', node2.lemma)
            
    # C->33->L
    feats[1317] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '33') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '33', node2.lemma)
            
    # C->acc. fem->L
    feats[1318] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'acc. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. fem', node2.lemma)
            
    # C->9_sg->L
    feats[1319] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '9_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_sg', node2.lemma)
            
    # C->16_pl->L
    feats[1320] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '16_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_pl', node2.lemma)
            
    # C->129->L
    feats[1321] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '129') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '129', node2.lemma)
            
    # C->-94->L
    feats[1322] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-94') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-94', node2.lemma)
            
    # C->nom. sg.->L
    feats[1323] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'nom. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. sg.', node2.lemma)
            
    # C->142->L
    feats[1324] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '142') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '142', node2.lemma)
            
    # C->nom. fem->L
    feats[1325] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'nom. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. fem', node2.lemma)
            
    # C->-276->L
    feats[1326] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-276') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-276', node2.lemma)
            
    # C->111->L
    feats[1327] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '111') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '111', node2.lemma)
            
    # C->15_fp->L
    feats[1328] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '15_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_fp', node2.lemma)
            
    # C->fp->L
    feats[1329] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'fp', node2.lemma)
            
    # C->nom. masc.->L
    feats[1330] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'nom. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. masc.', node2.lemma)
            
    # C->-271->L
    feats[1331] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-271') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-271', node2.lemma)
            
    # C->130->L
    feats[1332] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '130') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '130', node2.lemma)
            
    # C->5_fp->L
    feats[1333] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '5_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_fp', node2.lemma)
            
    # C->174->L
    feats[1334] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '174') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '174', node2.lemma)
            
    # C->99->L
    feats[1335] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '99') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '99', node2.lemma)
            
    # C->-30->L
    feats[1336] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-30') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-30', node2.lemma)
            
    # C->voc. fem->L
    feats[1337] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'voc. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. fem', node2.lemma)
            
    # C->-73->L
    feats[1338] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-73') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-73', node2.lemma)
            
    # C->nom->L
    feats[1339] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'nom') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom', node2.lemma)
            
    # C->78->L
    feats[1340] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '78') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '78', node2.lemma)
            
    # C->-13->L
    feats[1341] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-13') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-13', node2.lemma)
            
    # C->115->L
    feats[1342] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '115') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '115', node2.lemma)
            
    # C->15_pl->L
    feats[1343] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '15_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_pl', node2.lemma)
            
    # C->70->L
    feats[1344] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '70') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '70', node2.lemma)
            
    # C->35->L
    feats[1345] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '35') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '35', node2.lemma)
            
    # C->-247->L
    feats[1346] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-247') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-247', node2.lemma)
            
    # C->157->L
    feats[1347] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '157') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '157', node2.lemma)
            
    # C->1->L
    feats[1348] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '1') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '1', node2.lemma)
            
    # C->-249->L
    feats[1349] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-249') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-249', node2.lemma)
            
    # C->-38->L
    feats[1350] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-38') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-38', node2.lemma)
            
    # C->161->L
    feats[1351] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '161') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '161', node2.lemma)
            
    # C->27_sg->L
    feats[1352] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '27_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_sg', node2.lemma)
            
    # C->-55->L
    feats[1353] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-55') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-55', node2.lemma)
            
    # C->acc->L
    feats[1354] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'acc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc', node2.lemma)
            
    # C->-51->L
    feats[1355] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-51') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-51', node2.lemma)
            
    # C->-268->L
    feats[1356] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-268') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-268', node2.lemma)
            
    # C->-16->L
    feats[1357] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-16') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-16', node2.lemma)
            
    # C->-52->L
    feats[1358] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-52') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-52', node2.lemma)
            
    # C->loc. pl.->L
    feats[1359] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'loc. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc. pl.', node2.lemma)
            
    # C->13_pl->L
    feats[1360] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '13_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_pl', node2.lemma)
            
    # C->42->L
    feats[1361] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '42') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '42', node2.lemma)
            
    # C->-84->L
    feats[1362] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-84') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-84', node2.lemma)
            
    # C->27_fp->L
    feats[1363] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '27_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_fp', node2.lemma)
            
    # C->-20->L
    feats[1364] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-20') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-20', node2.lemma)
            
    # C->141->L
    feats[1365] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '141') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '141', node2.lemma)
            
    # C->139->L
    feats[1366] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '139') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '139', node2.lemma)
            
    # C->-129->L
    feats[1367] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-129') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-129', node2.lemma)
            
    # C->120->L
    feats[1368] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '120') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '120', node2.lemma)
            
    # C->-39->L
    feats[1369] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-39') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-39', node2.lemma)
            
    # C->-166->L
    feats[1370] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-166') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-166', node2.lemma)
            
    # C->-28->L
    feats[1371] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-28') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-28', node2.lemma)
            
    # C->170->L
    feats[1372] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '170') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '170', node2.lemma)
            
    # C->-131->L
    feats[1373] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-131') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-131', node2.lemma)
            
    # C->-19->L
    feats[1374] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-19') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-19', node2.lemma)
            
    # C->2_du->L
    feats[1375] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '2_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_du', node2.lemma)
            
    # C->-112->L
    feats[1376] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-112') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-112', node2.lemma)
            
    # C->12_sp->L
    feats[1377] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '12_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_sp', node2.lemma)
            
    # C->-104->L
    feats[1378] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-104') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-104', node2.lemma)
            
    # C->-12->L
    feats[1379] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-12') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-12', node2.lemma)
            
    # C->acc. neutr.->L
    feats[1380] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'acc. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. neutr.', node2.lemma)
            
    # C->14_tp->L
    feats[1381] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '14_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_tp', node2.lemma)
            
    # C->-29->L
    feats[1382] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-29') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-29', node2.lemma)
            
    # C->sg->L
    feats[1383] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg', node2.lemma)
            
    # C->177->L
    feats[1384] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '177') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '177', node2.lemma)
            
    # C->7_tp->L
    feats[1385] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '7_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_tp', node2.lemma)
            
    # C->-81->L
    feats[1386] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-81') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-81', node2.lemma)
            
    # C->-230->L
    feats[1387] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-230') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-230', node2.lemma)
            
    # C->-297->L
    feats[1388] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-297') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-297', node2.lemma)
            
    # C->voc. masc.->L
    feats[1389] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'voc. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. masc.', node2.lemma)
            
    # C->172->L
    feats[1390] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '172') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '172', node2.lemma)
            
    # C->11_pl->L
    feats[1391] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '11_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_pl', node2.lemma)
            
    # C->-11->L
    feats[1392] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-11') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-11', node2.lemma)
            
    # C->13_tp->L
    feats[1393] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '13_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_tp', node2.lemma)
            
    # C->160->L
    feats[1394] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '160') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
            
    # C->-115->L
    feats[1395] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-115') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-115', node2.lemma)
            
    # C->-50->L
    feats[1396] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-50') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-50', node2.lemma)
            
    # C->6_fp->L
    feats[1397] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '6_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_fp', node2.lemma)
            
    # C->68->L
    feats[1398] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '68') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '68', node2.lemma)
            
    # C->121->C
    feats[1399] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '121') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '121', node2.cng)
            
    # C->-15->C
    feats[1400] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-15') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-15', node2.cng)
            
    # C->-121->C
    feats[1401] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-121') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-121', node2.cng)
            
    # C->117->C
    feats[1402] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '117') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '117', node2.cng)
            
    # C->10_sg->C
    feats[1403] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '10_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sg', node2.cng)
            
    # C->-240->C
    feats[1404] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-240') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-240', node2.cng)
            
    # C->55->C
    feats[1405] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '55') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '55', node2.cng)
            
    # C->-245->C
    feats[1406] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-245') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-245', node2.cng)
            
    # C->171->C
    feats[1407] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '171') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '171', node2.cng)
            
    # C->154->C
    feats[1408] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '154') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '154', node2.cng)
            
    # C->14_sp->C
    feats[1409] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '14_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_sp', node2.cng)
            
    # C->-151->C
    feats[1410] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-151') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-151', node2.cng)
            
    # C->du_tp->C
    feats[1411] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'du_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_tp', node2.cng)
            
    # C->40->C
    feats[1412] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '40') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '40', node2.cng)
            
    # C->-79->C
    feats[1413] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-79', node2.cng)
            
    # C->12_tp->C
    feats[1414] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '12_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_tp', node2.cng)
            
    # C->101->C
    feats[1415] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '101') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '101', node2.cng)
            
    # C->29_tp->C
    feats[1416] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '29_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '29_tp', node2.cng)
            
    # C->pl_sp->C
    feats[1417] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'pl_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_sp', node2.cng)
            
    # C->135->C
    feats[1418] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '135') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', node2.cng)
            
    # C->-132->C
    feats[1419] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-132') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-132', node2.cng)
            
    # C->-122->C
    feats[1420] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-122') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-122', node2.cng)
            
    # C->acc. adj.->C
    feats[1421] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'acc. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. adj.', node2.cng)
            
    # C->102->C
    feats[1422] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '102') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '102', node2.cng)
            
    # C->-152->C
    feats[1423] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-152') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-152', node2.cng)
            
    # C->3_pl->C
    feats[1424] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '3_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_pl', node2.cng)
            
    # C->-66->C
    feats[1425] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-66') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-66', node2.cng)
            
    # C->151->C
    feats[1426] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '151') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '151', node2.cng)
            
    # C->15_du->C
    feats[1427] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '15_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_du', node2.cng)
            
    # C->82->C
    feats[1428] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '82') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '82', node2.cng)
            
    # C->instr. fem->C
    feats[1429] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'instr. fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. fem', node2.cng)
            
    # C->-242->C
    feats[1430] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-242') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-242', node2.cng)
            
    # C->27_du->C
    feats[1431] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', node2.cng)
            
    # C->92->C
    feats[1432] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '92', node2.cng)
            
    # C->gen. pl.->C
    feats[1433] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'gen. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. pl.', node2.cng)
            
    # C->74->C
    feats[1434] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '74') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '74', node2.cng)
            
    # C->-86->C
    feats[1435] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-86') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-86', node2.cng)
            
    # C->8_du->C
    feats[1436] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '8_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_du', node2.cng)
            
    # C->97->C
    feats[1437] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '97') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '97', node2.cng)
            
    # C->15_sg->C
    feats[1438] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '15_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_sg', node2.cng)
            
    # C->7_du->C
    feats[1439] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '7_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_du', node2.cng)
            
    # C->11_sp->C
    feats[1440] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '11_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_sp', node2.cng)
            
    # C->89->C
    feats[1441] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '89') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '89', node2.cng)
            
    # C->-261->C
    feats[1442] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-261') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-261', node2.cng)
            
    # C->16_fp->C
    feats[1443] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '16_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_fp', node2.cng)
            
    # C->178->C
    feats[1444] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '178') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '178', node2.cng)
            
    # C->-200->C
    feats[1445] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-200') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-200', node2.cng)
            
    # C->4_sp->C
    feats[1446] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '4_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_sp', node2.cng)
            
    # C->adj->C
    feats[1447] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'adj') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'adj', node2.cng)
            
    # C->10_sp->C
    feats[1448] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '10_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sp', node2.cng)
            
    # C->-22->C
    feats[1449] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-22') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-22', node2.cng)
            
    # C->-46->C
    feats[1450] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-46') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-46', node2.cng)
            
    # C->30_fp->C
    feats[1451] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '30_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_fp', node2.cng)
            
    # C->instr. pl.->C
    feats[1452] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'instr. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. pl.', node2.cng)
            
    # C->nom. adj.->C
    feats[1453] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'nom. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. adj.', node2.cng)
            
    # C->-27->C
    feats[1454] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-27') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-27', node2.cng)
            
    # C->-35->C
    feats[1455] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-35', node2.cng)
            
    # C->-243->C
    feats[1456] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-243') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-243', node2.cng)
            
    # C->-18->C
    feats[1457] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', node2.cng)
            
    # C->voc. sg.->C
    feats[1458] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'voc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. sg.', node2.cng)
            
    # C->49->C
    feats[1459] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '49') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '49', node2.cng)
            
    # C->-89->C
    feats[1460] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-89') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-89', node2.cng)
            
    # C->-302->C
    feats[1461] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', node2.cng)
            
    # C->11_tp->C
    feats[1462] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '11_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_tp', node2.cng)
            
    # C->-139->C
    feats[1463] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', node2.cng)
            
    # C->28->C
    feats[1464] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '28') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28', node2.cng)
            
    # C->-56->C
    feats[1465] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-56') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-56', node2.cng)
            
    # C->acc. masc.->C
    feats[1466] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'acc. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. masc.', node2.cng)
            
    # C->-47->C
    feats[1467] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-47') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-47', node2.cng)
            
    # C->-17->C
    feats[1468] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-17') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-17', node2.cng)
            
    # C->5_du->C
    feats[1469] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '5_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_du', node2.cng)
            
    # C->98->C
    feats[1470] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '98', node2.cng)
            
    # C->81->C
    feats[1471] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '81', node2.cng)
            
    # C->sg_tp->C
    feats[1472] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'sg_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg_tp', node2.cng)
            
    # C->-169->C
    feats[1473] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-169') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-169', node2.cng)
            
    # C->-210->C
    feats[1474] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', node2.cng)
            
    # C->-25->C
    feats[1475] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-25') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-25', node2.cng)
            
    # C->162->C
    feats[1476] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '162') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '162', node2.cng)
            
    # C->-101->C
    feats[1477] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-101') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-101', node2.cng)
            
    # C->instr. adj.->C
    feats[1478] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'instr. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. adj.', node2.cng)
            
    # C->du_fp->C
    feats[1479] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'du_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_fp', node2.cng)
            
    # C->-303->C
    feats[1480] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', node2.cng)
            
    # C->5_tp->C
    feats[1481] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '5_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_tp', node2.cng)
            
    # C->11_du->C
    feats[1482] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '11_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_du', node2.cng)
            
    # C->15_sp->C
    feats[1483] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '15_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_sp', node2.cng)
            
    # C->30_du->C
    feats[1484] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '30_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_du', node2.cng)
            
    # C->128->C
    feats[1485] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '128') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '128', node2.cng)
            
    # C->-273->C
    feats[1486] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', node2.cng)
            
    # C->-137->C
    feats[1487] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-137') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-137', node2.cng)
            
    # C->9_pl->C
    feats[1488] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '9_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_pl', node2.cng)
            
    # C->acc. pl.->C
    feats[1489] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'acc. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. pl.', node2.cng)
            
    # C->-103->C
    feats[1490] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-103') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-103', node2.cng)
            
    # C->71->C
    feats[1491] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '71') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '71', node2.cng)
            
    # C->-306->C
    feats[1492] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', node2.cng)
            
    # C->dat. du.->C
    feats[1493] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'dat. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. du.', node2.cng)
            
    # C->109->C
    feats[1494] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '109') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '109', node2.cng)
            
    # C->-291->C
    feats[1495] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-291') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-291', node2.cng)
            
    # C->-296->C
    feats[1496] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-296') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-296', node2.cng)
            
    # C->11_sg->C
    feats[1497] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '11_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_sg', node2.cng)
            
    # C->-57->C
    feats[1498] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-57') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-57', node2.cng)
            
    # C->-77->C
    feats[1499] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-77') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-77', node2.cng)
            
    feats[feats < 1e-25] = 1e-25
    return -np.log10(feats)
    