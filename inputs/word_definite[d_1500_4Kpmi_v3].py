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
    
    # L->C
    feats[0] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, node2.cng)
            
    # C->L
    feats[1] = tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, node1.cng, node2.lemma)
            
    # C->C
    feats[2] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, node2.cng)
            
    # C->T
    feats[3] = tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, node1.cng, node2.tup)
            
    # T->C
    feats[4] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, node2.cng)
            
    # T->T
    feats[5] = tryProb_catchZero(mat_tup2tup_countonly, mat_tupCount_1D, node1.tup, node2.tup)
            
    # L->abl. du.->L
    feats[6] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'abl. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. du.', node2.lemma)
            
    # L->-150->L
    feats[7] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-150') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-150', node2.lemma)
            
    # L->3_tp->L
    feats[8] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_tp', node2.lemma)
            
    # L->39->L
    feats[9] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '39', node2.lemma)
            
    # L->voc. masc.->L
    feats[10] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. masc.', node2.lemma)
            
    # L->-94->L
    feats[11] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-94') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-94', node2.lemma)
            
    # L->-31->L
    feats[12] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-31') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-31', node2.lemma)
            
    # L->-122->L
    feats[13] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-122') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-122', node2.lemma)
            
    # L->13_sg->L
    feats[14] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_sg', node2.lemma)
            
    # L->72->L
    feats[15] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '72') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '72', node2.lemma)
            
    # L->12_fp->L
    feats[16] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_fp', node2.lemma)
            
    # L->137->L
    feats[17] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '137') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
            
    # L->-169->L
    feats[18] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-169') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-169', node2.lemma)
            
    # L->-11->L
    feats[19] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-11') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-11', node2.lemma)
            
    # L->loc. sg.->L
    feats[20] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc. sg.', node2.lemma)
            
    # L->170->L
    feats[21] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '170') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '170', node2.lemma)
            
    # L->9_sp->L
    feats[22] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_sp', node2.lemma)
            
    # L->5_sp->L
    feats[23] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sp', node2.lemma)
            
    # L->34->L
    feats[24] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '34') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '34', node2.lemma)
            
    # L->179->L
    feats[25] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '179') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '179', node2.lemma)
            
    # L->-152->L
    feats[26] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-152') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-152', node2.lemma)
            
    # L->tp->L
    feats[27] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'tp', node2.lemma)
            
    # L->-149->L
    feats[28] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-149') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-149', node2.lemma)
            
    # L->10_tp->L
    feats[29] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_tp', node2.lemma)
            
    # L->-77->L
    feats[30] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-77') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-77', node2.lemma)
            
    # L->-245->L
    feats[31] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-245') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-245', node2.lemma)
            
    # L->instr. du.->L
    feats[32] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. du.', node2.lemma)
            
    # L->151->L
    feats[33] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '151') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '151', node2.lemma)
            
    # L->2_fp->L
    feats[34] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_fp', node2.lemma)
            
    # L->5_pl->L
    feats[35] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_pl', node2.lemma)
            
    # L->-230->L
    feats[36] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-230', node2.lemma)
            
    # L->-132->L
    feats[37] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-132') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-132', node2.lemma)
            
    # L->acc. adj.->L
    feats[38] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. adj.', node2.lemma)
            
    # L->120->L
    feats[39] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '120') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '120', node2.lemma)
            
    # L->182->L
    feats[40] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '182') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '182', node2.lemma)
            
    # L->28_tp->L
    feats[41] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28_tp', node2.lemma)
            
    # L->115->L
    feats[42] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '115') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '115', node2.lemma)
            
    # L->156->L
    feats[43] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '156') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '156', node2.lemma)
            
    # L->118->L
    feats[44] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '118') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '118', node2.lemma)
            
    # L->acc. fem->L
    feats[45] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. fem', node2.lemma)
            
    # L->121->L
    feats[46] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '121') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '121', node2.lemma)
            
    # L->-302->L
    feats[47] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-302', node2.lemma)
            
    # L->voc. sg.->L
    feats[48] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. sg.', node2.lemma)
            
    # L->fp->L
    feats[49] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'fp', node2.lemma)
            
    # L->-23->L
    feats[50] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-23') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-23', node2.lemma)
            
    # L->-104->L
    feats[51] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-104', node2.lemma)
            
    # L->2_du->L
    feats[52] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_du', node2.lemma)
            
    # L->pl->L
    feats[53] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
            
    # L->-44->L
    feats[54] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-44') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-44', node2.lemma)
            
    # L->-292->L
    feats[55] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-292') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-292', node2.lemma)
            
    # L->114->L
    feats[56] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '114') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '114', node2.lemma)
            
    # L->-36->L
    feats[57] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-36') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-36', node2.lemma)
            
    # L->-156->L
    feats[58] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-156') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-156', node2.lemma)
            
    # L->70->L
    feats[59] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '70') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '70', node2.lemma)
            
    # L->-27->L
    feats[60] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-27') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-27', node2.lemma)
            
    # L->-144->L
    feats[61] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-144') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-144', node2.lemma)
            
    # L->116->L
    feats[62] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '116') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '116', node2.lemma)
            
    # L->-24->L
    feats[63] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-24') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
            
    # L->12_sg->L
    feats[64] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_sg', node2.lemma)
            
    # L->-51->L
    feats[65] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-51') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-51', node2.lemma)
            
    # L->abl. sg.->L
    feats[66] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'abl. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. sg.', node2.lemma)
            
    # L->-293->L
    feats[67] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-293', node2.lemma)
            
    # L->6_sp->L
    feats[68] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_sp', node2.lemma)
            
    # L->-63->L
    feats[69] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-63') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-63', node2.lemma)
            
    # L->acc. sg.->L
    feats[70] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. sg.', node2.lemma)
            
    # L->-89->L
    feats[71] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-89') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-89', node2.lemma)
            
    # L->sg_sp->L
    feats[72] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_sp', node2.lemma)
            
    # L->42->L
    feats[73] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '42') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '42', node2.lemma)
            
    # L->7_fp->L
    feats[74] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_fp', node2.lemma)
            
    # L->139->L
    feats[75] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '139') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '139', node2.lemma)
            
    # L->6_tp->L
    feats[76] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_tp', node2.lemma)
            
    # L->-301->L
    feats[77] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-301') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-301', node2.lemma)
            
    # L->73->L
    feats[78] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '73') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '73', node2.lemma)
            
    # L->-247->L
    feats[79] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-247') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-247', node2.lemma)
            
    # L->7_pl->L
    feats[80] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_pl', node2.lemma)
            
    # L->50->L
    feats[81] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '50') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '50', node2.lemma)
            
    # L->instr. masc.->L
    feats[82] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. masc.', node2.lemma)
            
    # L->119->L
    feats[83] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '119') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '119', node2.lemma)
            
    # L->-153->L
    feats[84] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-153') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-153', node2.lemma)
            
    # L->-67->L
    feats[85] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-67', node2.lemma)
            
    # L->-76->L
    feats[86] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-76') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-76', node2.lemma)
            
    # L->51->L
    feats[87] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '51') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '51', node2.lemma)
            
    # L->8_fp->L
    feats[88] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
            
    # L->80->L
    feats[89] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '80') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '80', node2.lemma)
            
    # L->16_fp->L
    feats[90] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_fp', node2.lemma)
            
    # L->-133->L
    feats[91] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-133') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-133', node2.lemma)
            
    # L->173->L
    feats[92] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '173') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '173', node2.lemma)
            
    # L->-61->L
    feats[93] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-61') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-61', node2.lemma)
            
    # L->58->L
    feats[94] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '58') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '58', node2.lemma)
            
    # L->3_sg->L
    feats[95] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_sg', node2.lemma)
            
    # L->gen. sg.->L
    feats[96] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen. sg.', node2.lemma)
            
    # L->instr. fem->L
    feats[97] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. fem', node2.lemma)
            
    # L->pl_tp->L
    feats[98] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_tp', node2.lemma)
            
    # L->8_sp->L
    feats[99] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_sp', node2.lemma)
            
    # L->60->L
    feats[100] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '60') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '60', node2.lemma)
            
    # L->14_pl->L
    feats[101] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_pl', node2.lemma)
            
    # L->14_sg->L
    feats[102] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_sg', node2.lemma)
            
    # L->15_pl->L
    feats[103] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_pl', node2.lemma)
            
    # L->-83->L
    feats[104] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-83') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-83', node2.lemma)
            
    # L->93->L
    feats[105] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '93') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '93', node2.lemma)
            
    # L->15_du->L
    feats[106] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_du', node2.lemma)
            
    # L->10_du->L
    feats[107] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_du', node2.lemma)
            
    # L->15_tp->L
    feats[108] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_tp', node2.lemma)
            
    # L->-32->L
    feats[109] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-32') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-32', node2.lemma)
            
    # L->134->L
    feats[110] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '134') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '134', node2.lemma)
            
    # L->-43->L
    feats[111] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-43') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-43', node2.lemma)
            
    # L->-299->L
    feats[112] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-299') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-299', node2.lemma)
            
    # L->27_pl->L
    feats[113] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_pl', node2.lemma)
            
    # L->8_sg->L
    feats[114] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_sg', node2.lemma)
            
    # L->69->L
    feats[115] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '69') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '69', node2.lemma)
            
    # L->155->L
    feats[116] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '155') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '155', node2.lemma)
            
    # L->154->L
    feats[117] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '154') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '154', node2.lemma)
            
    # L->-57->L
    feats[118] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-57') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-57', node2.lemma)
            
    # L->108->L
    feats[119] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '108') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '108', node2.lemma)
            
    # L->2_tp->L
    feats[120] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_tp', node2.lemma)
            
    # L->-35->L
    feats[121] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-35') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
            
    # L->-240->L
    feats[122] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-240') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-240', node2.lemma)
            
    # L->-90->L
    feats[123] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-90') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-90', node2.lemma)
            
    # L->161->L
    feats[124] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '161') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '161', node2.lemma)
            
    # L->10_pl->L
    feats[125] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_pl', node2.lemma)
            
    # L->27_du->L
    feats[126] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_du', node2.lemma)
            
    # L->15_fp->L
    feats[127] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_fp', node2.lemma)
            
    # L->13_fp->L
    feats[128] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_fp', node2.lemma)
            
    # L->-64->L
    feats[129] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-64') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-64', node2.lemma)
            
    # L->176->L
    feats[130] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '176') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '176', node2.lemma)
            
    # L->100->L
    feats[131] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '100') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '100', node2.lemma)
            
    # L->-28->L
    feats[132] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-28') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-28', node2.lemma)
            
    # L->abl. du.->C
    feats[133] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'abl. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'abl. du.', node2.cng)
            
    # L->-150->C
    feats[134] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-150') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-150', node2.cng)
            
    # L->3_tp->C
    feats[135] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_tp', node2.cng)
            
    # L->voc. masc.->C
    feats[136] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. masc.', node2.cng)
            
    # L->-94->C
    feats[137] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-94') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-94', node2.cng)
            
    # L->-31->C
    feats[138] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-31') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-31', node2.cng)
            
    # L->13_sg->C
    feats[139] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_sg', node2.cng)
            
    # L->137->C
    feats[140] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '137') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '137', node2.cng)
            
    # L->-169->C
    feats[141] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-169') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-169', node2.cng)
            
    # L->-11->C
    feats[142] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-11') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-11', node2.cng)
            
    # L->loc. sg.->C
    feats[143] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', node2.cng)
            
    # L->170->C
    feats[144] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '170') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '170', node2.cng)
            
    # L->9_sp->C
    feats[145] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', node2.cng)
            
    # L->34->C
    feats[146] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '34', node2.cng)
            
    # L->179->C
    feats[147] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '179') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '179', node2.cng)
            
    # L->-152->C
    feats[148] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-152') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-152', node2.cng)
            
    # L->acc. neutr.->C
    feats[149] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. neutr.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. neutr.', node2.cng)
            
    # L->tp->C
    feats[150] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', node2.cng)
            
    # L->-149->C
    feats[151] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-149', node2.cng)
            
    # L->10_tp->C
    feats[152] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_tp', node2.cng)
            
    # L->-245->C
    feats[153] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-245') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-245', node2.cng)
            
    # L->instr. du.->C
    feats[154] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. du.', node2.cng)
            
    # L->151->C
    feats[155] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '151') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '151', node2.cng)
            
    # L->2_fp->C
    feats[156] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', node2.cng)
            
    # L->5_pl->C
    feats[157] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_pl', node2.cng)
            
    # L->-132->C
    feats[158] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-132') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-132', node2.cng)
            
    # L->acc. adj.->C
    feats[159] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. adj.', node2.cng)
            
    # L->182->C
    feats[160] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '182') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '182', node2.cng)
            
    # L->28_tp->C
    feats[161] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28_tp', node2.cng)
            
    # L->156->C
    feats[162] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '156', node2.cng)
            
    # L->121->C
    feats[163] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '121') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '121', node2.cng)
            
    # L->-302->C
    feats[164] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', node2.cng)
            
    # L->-23->C
    feats[165] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-23') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-23', node2.cng)
            
    # L->-104->C
    feats[166] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', node2.cng)
            
    # L->2_du->C
    feats[167] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', node2.cng)
            
    # L->pl->C
    feats[168] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl', node2.cng)
            
    # L->-44->C
    feats[169] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-44') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-44', node2.cng)
            
    # L->12_pl->C
    feats[170] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_pl', node2.cng)
            
    # L->-292->C
    feats[171] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-292') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-292', node2.cng)
            
    # L->114->C
    feats[172] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '114') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '114', node2.cng)
            
    # L->-36->C
    feats[173] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-36') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-36', node2.cng)
            
    # L->-156->C
    feats[174] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-156', node2.cng)
            
    # L->-27->C
    feats[175] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-27') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-27', node2.cng)
            
    # L->-144->C
    feats[176] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-144') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-144', node2.cng)
            
    # L->116->C
    feats[177] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '116') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '116', node2.cng)
            
    # L->-24->C
    feats[178] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-24') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-24', node2.cng)
            
    # L->-51->C
    feats[179] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-51') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-51', node2.cng)
            
    # L->abl. sg.->C
    feats[180] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'abl. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'abl. sg.', node2.cng)
            
    # L->-293->C
    feats[181] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', node2.cng)
            
    # L->6_sp->C
    feats[182] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sp', node2.cng)
            
    # L->-63->C
    feats[183] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-63') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-63', node2.cng)
            
    # L->-89->C
    feats[184] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-89') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-89', node2.cng)
            
    # L->sg_sp->C
    feats[185] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg_sp', node2.cng)
            
    # L->nom. neutr.->C
    feats[186] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. neutr.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. neutr.', node2.cng)
            
    # L->42->C
    feats[187] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '42') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '42', node2.cng)
            
    # L->7_fp->C
    feats[188] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_fp', node2.cng)
            
    # L->6_tp->C
    feats[189] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_tp', node2.cng)
            
    # L->73->C
    feats[190] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '73') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '73', node2.cng)
            
    # L->7_pl->C
    feats[191] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', node2.cng)
            
    # L->50->C
    feats[192] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '50') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '50', node2.cng)
            
    # L->instr. masc.->C
    feats[193] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. masc.', node2.cng)
            
    # L->-67->C
    feats[194] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', node2.cng)
            
    # L->-76->C
    feats[195] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-76') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-76', node2.cng)
            
    # L->8_fp->C
    feats[196] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_fp', node2.cng)
            
    # L->16_fp->C
    feats[197] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_fp', node2.cng)
            
    # L->-133->C
    feats[198] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-133') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-133', node2.cng)
            
    # L->173->C
    feats[199] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '173') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '173', node2.cng)
            
    # L->-61->C
    feats[200] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-61') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-61', node2.cng)
            
    # L->58->C
    feats[201] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '58', node2.cng)
            
    # L->3_sg->C
    feats[202] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_sg', node2.cng)
            
    # L->gen. sg.->C
    feats[203] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. sg.', node2.cng)
            
    # L->instr. fem->C
    feats[204] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. fem', node2.cng)
            
    # L->8_sp->C
    feats[205] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_sp', node2.cng)
            
    # L->60->C
    feats[206] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '60') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '60', node2.cng)
            
    # L->14_pl->C
    feats[207] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_pl', node2.cng)
            
    # L->14_sg->C
    feats[208] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_sg', node2.cng)
            
    # L->15_pl->C
    feats[209] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_pl', node2.cng)
            
    # L->-83->C
    feats[210] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-83') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-83', node2.cng)
            
    # L->93->C
    feats[211] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '93') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '93', node2.cng)
            
    # L->6_du->C
    feats[212] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_du', node2.cng)
            
    # L->111->C
    feats[213] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '111') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '111', node2.cng)
            
    # L->10_du->C
    feats[214] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_du', node2.cng)
            
    # L->instr. neutr.->C
    feats[215] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. neutr.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. neutr.', node2.cng)
            
    # L->175->C
    feats[216] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '175') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '175', node2.cng)
            
    # L->134->C
    feats[217] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '134') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '134', node2.cng)
            
    # L->-43->C
    feats[218] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-43') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-43', node2.cng)
            
    # L->-299->C
    feats[219] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-299') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-299', node2.cng)
            
    # L->27_pl->C
    feats[220] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_pl', node2.cng)
            
    # L->8_sg->C
    feats[221] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_sg', node2.cng)
            
    # L->69->C
    feats[222] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '69') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '69', node2.cng)
            
    # L->155->C
    feats[223] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '155') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '155', node2.cng)
            
    # L->154->C
    feats[224] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '154') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '154', node2.cng)
            
    # L->-57->C
    feats[225] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-57') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-57', node2.cng)
            
    # L->108->C
    feats[226] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '108') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '108', node2.cng)
            
    # L->voc. pl.->C
    feats[227] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. pl.', node2.cng)
            
    # L->-35->C
    feats[228] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-35', node2.cng)
            
    # L->-240->C
    feats[229] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-240') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-240', node2.cng)
            
    # L->-90->C
    feats[230] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-90', node2.cng)
            
    # L->161->C
    feats[231] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '161') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '161', node2.cng)
            
    # L->10_pl->C
    feats[232] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_pl', node2.cng)
            
    # L->27_du->C
    feats[233] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', node2.cng)
            
    # L->15_fp->C
    feats[234] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_fp', node2.cng)
            
    # L->13_fp->C
    feats[235] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_fp', node2.cng)
            
    # L->-64->C
    feats[236] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-64') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-64', node2.cng)
            
    # L->176->C
    feats[237] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '176') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '176', node2.cng)
            
    # L->100->C
    feats[238] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '100') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '100', node2.cng)
            
    # L->-28->C
    feats[239] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-28') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-28', node2.cng)
            
    # L->abl. du.->T
    feats[240] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'abl. du.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'abl. du.', node2.tup)
            
    # L->-150->T
    feats[241] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-150') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-150', node2.tup)
            
    # L->-159->T
    feats[242] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-159') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-159', node2.tup)
            
    # L->3_tp->T
    feats[243] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_tp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '3_tp', node2.tup)
            
    # L->39->T
    feats[244] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '39') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '39', node2.tup)
            
    # L->voc. masc.->T
    feats[245] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. masc.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'voc. masc.', node2.tup)
            
    # L->-94->T
    feats[246] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-94') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-94', node2.tup)
            
    # L->-31->T
    feats[247] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-31') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-31', node2.tup)
            
    # L->-122->T
    feats[248] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-122') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-122', node2.tup)
            
    # L->13_sg->T
    feats[249] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_sg') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '13_sg', node2.tup)
            
    # L->72->T
    feats[250] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '72') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '72', node2.tup)
            
    # L->12_fp->T
    feats[251] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_fp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '12_fp', node2.tup)
            
    # L->137->T
    feats[252] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '137') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '137', node2.tup)
            
    # L->76->T
    feats[253] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '76') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '76', node2.tup)
            
    # L->-169->T
    feats[254] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-169') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-169', node2.tup)
            
    # L->-123->T
    feats[255] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-123') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-123', node2.tup)
            
    # L->-11->T
    feats[256] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-11') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-11', node2.tup)
            
    # L->loc. sg.->T
    feats[257] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'loc. sg.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'loc. sg.', node2.tup)
            
    # L->170->T
    feats[258] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '170') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '170', node2.tup)
            
    # L->9_sp->T
    feats[259] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '9_sp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '9_sp', node2.tup)
            
    # L->5_sp->T
    feats[260] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_sp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '5_sp', node2.tup)
            
    # L->14_sp->T
    feats[261] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_sp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '14_sp', node2.tup)
            
    # L->34->T
    feats[262] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '34') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '34', node2.tup)
            
    # L->179->T
    feats[263] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '179') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '179', node2.tup)
            
    # L->-152->T
    feats[264] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-152') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-152', node2.tup)
            
    # L->acc. neutr.->T
    feats[265] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. neutr.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'acc. neutr.', node2.tup)
            
    # L->tp->T
    feats[266] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'tp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'tp', node2.tup)
            
    # L->10_fp->T
    feats[267] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_fp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '10_fp', node2.tup)
            
    # L->-149->T
    feats[268] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-149') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-149', node2.tup)
            
    # L->10_tp->T
    feats[269] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_tp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '10_tp', node2.tup)
            
    # L->-77->T
    feats[270] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-77') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-77', node2.tup)
            
    # L->-245->T
    feats[271] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-245') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-245', node2.tup)
            
    # L->instr. du.->T
    feats[272] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. du.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'instr. du.', node2.tup)
            
    # L->-166->T
    feats[273] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-166') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-166', node2.tup)
            
    # L->151->T
    feats[274] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '151') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '151', node2.tup)
            
    # L->98->T
    feats[275] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '98') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '98', node2.tup)
            
    # L->177->T
    feats[276] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '177') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '177', node2.tup)
            
    # L->2_fp->T
    feats[277] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_fp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '2_fp', node2.tup)
            
    # L->5_pl->T
    feats[278] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '5_pl') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '5_pl', node2.tup)
            
    # L->-230->T
    feats[279] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-230') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-230', node2.tup)
            
    # L->-34->T
    feats[280] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-34') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-34', node2.tup)
            
    # L->-103->T
    feats[281] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-103') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-103', node2.tup)
            
    # L->-132->T
    feats[282] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-132') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-132', node2.tup)
            
    # L->acc. adj.->T
    feats[283] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. adj.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'acc. adj.', node2.tup)
            
    # L->120->T
    feats[284] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '120') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '120', node2.tup)
            
    # L->182->T
    feats[285] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '182') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '182', node2.tup)
            
    # L->28_tp->T
    feats[286] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '28_tp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '28_tp', node2.tup)
            
    # L->115->T
    feats[287] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '115') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '115', node2.tup)
            
    # L->156->T
    feats[288] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '156') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '156', node2.tup)
            
    # L->118->T
    feats[289] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '118') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '118', node2.tup)
            
    # L->acc. fem->T
    feats[290] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. fem') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'acc. fem', node2.tup)
            
    # L->-266->T
    feats[291] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-266') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-266', node2.tup)
            
    # L->121->T
    feats[292] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '121') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '121', node2.tup)
            
    # L->-302->T
    feats[293] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-302') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-302', node2.tup)
            
    # L->voc. sg.->T
    feats[294] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. sg.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'voc. sg.', node2.tup)
            
    # L->fp->T
    feats[295] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'fp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'fp', node2.tup)
            
    # L->-23->T
    feats[296] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-23') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-23', node2.tup)
            
    # L->-104->T
    feats[297] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-104') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-104', node2.tup)
            
    # L->2_du->T
    feats[298] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_du') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '2_du', node2.tup)
            
    # L->92->T
    feats[299] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '92') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '92', node2.tup)
            
    # L->pl->T
    feats[300] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'pl', node2.tup)
            
    # L->-44->T
    feats[301] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-44') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-44', node2.tup)
            
    # L->-292->T
    feats[302] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-292') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-292', node2.tup)
            
    # L->114->T
    feats[303] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '114') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '114', node2.tup)
            
    # L->-36->T
    feats[304] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-36') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-36', node2.tup)
            
    # L->-156->T
    feats[305] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-156') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-156', node2.tup)
            
    # L->70->T
    feats[306] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '70') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '70', node2.tup)
            
    # L->-27->T
    feats[307] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-27') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-27', node2.tup)
            
    # L->-144->T
    feats[308] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-144') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-144', node2.tup)
            
    # L->-111->T
    feats[309] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-111') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-111', node2.tup)
            
    # L->117->T
    feats[310] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '117') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '117', node2.tup)
            
    # L->dat. du.->T
    feats[311] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. du.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'dat. du.', node2.tup)
            
    # L->-307->T
    feats[312] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-307') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-307', node2.tup)
            
    # L->nom. du.->T
    feats[313] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. du.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'nom. du.', node2.tup)
            
    # L->116->T
    feats[314] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '116') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '116', node2.tup)
            
    # L->81->T
    feats[315] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '81') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '81', node2.tup)
            
    # L->95->T
    feats[316] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '95') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '95', node2.tup)
            
    # L->-24->T
    feats[317] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-24') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-24', node2.tup)
            
    # L->12_sg->T
    feats[318] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '12_sg') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '12_sg', node2.tup)
            
    # L->-51->T
    feats[319] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-51') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-51', node2.tup)
            
    # L->abl. sg.->T
    feats[320] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'abl. sg.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'abl. sg.', node2.tup)
            
    # L->-293->T
    feats[321] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-293') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-293', node2.tup)
            
    # L->6_sp->T
    feats[322] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_sp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '6_sp', node2.tup)
            
    # L->-63->T
    feats[323] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-63') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-63', node2.tup)
            
    # L->acc. sg.->T
    feats[324] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'acc. sg.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'acc. sg.', node2.tup)
            
    # L->-89->T
    feats[325] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-89') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-89', node2.tup)
            
    # L->sg_sp->T
    feats[326] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'sg_sp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'sg_sp', node2.tup)
            
    # L->nom. neutr.->T
    feats[327] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'nom. neutr.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'nom. neutr.', node2.tup)
            
    # L->42->T
    feats[328] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '42') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '42', node2.tup)
            
    # L->7_fp->T
    feats[329] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_fp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '7_fp', node2.tup)
            
    # L->139->T
    feats[330] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '139') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '139', node2.tup)
            
    # L->6_tp->T
    feats[331] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_tp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '6_tp', node2.tup)
            
    # L->-301->T
    feats[332] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-301') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-301', node2.tup)
            
    # L->73->T
    feats[333] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '73') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '73', node2.tup)
            
    # L->-247->T
    feats[334] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-247') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-247', node2.tup)
            
    # L->7_pl->T
    feats[335] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '7_pl') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '7_pl', node2.tup)
            
    # L->50->T
    feats[336] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '50') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '50', node2.tup)
            
    # L->119->T
    feats[337] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '119') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '119', node2.tup)
            
    # L->-153->T
    feats[338] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-153') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-153', node2.tup)
            
    # L->-67->T
    feats[339] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-67') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-67', node2.tup)
            
    # L->-76->T
    feats[340] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-76') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-76', node2.tup)
            
    # L->51->T
    feats[341] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '51') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '51', node2.tup)
            
    # L->8_fp->T
    feats[342] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_fp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '8_fp', node2.tup)
            
    # L->80->T
    feats[343] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '80') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '80', node2.tup)
            
    # L->16_fp->T
    feats[344] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '16_fp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '16_fp', node2.tup)
            
    # L->-133->T
    feats[345] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-133') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-133', node2.tup)
            
    # L->173->T
    feats[346] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '173') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '173', node2.tup)
            
    # L->-61->T
    feats[347] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-61') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-61', node2.tup)
            
    # L->58->T
    feats[348] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '58') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '58', node2.tup)
            
    # L->15_sg->T
    feats[349] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_sg') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '15_sg', node2.tup)
            
    # L->3_sg->T
    feats[350] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '3_sg') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '3_sg', node2.tup)
            
    # L->gen. sg.->T
    feats[351] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'gen. sg.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'gen. sg.', node2.tup)
            
    # L->-262->T
    feats[352] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-262') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-262', node2.tup)
            
    # L->instr. fem->T
    feats[353] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'instr. fem') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'instr. fem', node2.tup)
            
    # L->pl_tp->T
    feats[354] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'pl_tp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'pl_tp', node2.tup)
            
    # L->8_sp->T
    feats[355] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_sp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '8_sp', node2.tup)
            
    # L->60->T
    feats[356] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '60') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '60', node2.tup)
            
    # L->14_pl->T
    feats[357] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_pl') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '14_pl', node2.tup)
            
    # L->14_sg->T
    feats[358] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '14_sg') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '14_sg', node2.tup)
            
    # L->15_pl->T
    feats[359] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_pl') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '15_pl', node2.tup)
            
    # L->-83->T
    feats[360] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-83') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-83', node2.tup)
            
    # L->93->T
    feats[361] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '93') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '93', node2.tup)
            
    # L->6_du->T
    feats[362] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '6_du') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '6_du', node2.tup)
            
    # L->111->T
    feats[363] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '111') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '111', node2.tup)
            
    # L->15_du->T
    feats[364] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_du') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '15_du', node2.tup)
            
    # L->10_du->T
    feats[365] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_du') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '10_du', node2.tup)
            
    # L->dat. sg.->T
    feats[366] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'dat. sg.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'dat. sg.', node2.tup)
            
    # L->15_tp->T
    feats[367] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_tp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '15_tp', node2.tup)
            
    # L->-32->T
    feats[368] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-32') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-32', node2.tup)
            
    # L->-54->T
    feats[369] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-54') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-54', node2.tup)
            
    # L->134->T
    feats[370] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '134') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '134', node2.tup)
            
    # L->-43->T
    feats[371] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-43') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-43', node2.tup)
            
    # L->-190->T
    feats[372] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-190') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-190', node2.tup)
            
    # L->30_fp->T
    feats[373] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '30_fp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '30_fp', node2.tup)
            
    # L->-299->T
    feats[374] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-299') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-299', node2.tup)
            
    # L->2->T
    feats[375] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '2', node2.tup)
            
    # L->27_pl->T
    feats[376] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_pl') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '27_pl', node2.tup)
            
    # L->8_sg->T
    feats[377] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_sg') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '8_sg', node2.tup)
            
    # L->69->T
    feats[378] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '69') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '69', node2.tup)
            
    # L->155->T
    feats[379] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '155') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '155', node2.tup)
            
    # L->154->T
    feats[380] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '154') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '154', node2.tup)
            
    # L->-57->T
    feats[381] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-57') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-57', node2.tup)
            
    # L->108->T
    feats[382] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '108') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '108', node2.tup)
            
    # L->-87->T
    feats[383] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-87') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-87', node2.tup)
            
    # L->2_tp->T
    feats[384] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '2_tp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '2_tp', node2.tup)
            
    # L->voc. pl.->T
    feats[385] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, 'voc. pl.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'voc. pl.', node2.tup)
            
    # L->-35->T
    feats[386] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-35') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-35', node2.tup)
            
    # L->4_tp->T
    feats[387] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '4_tp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '4_tp', node2.tup)
            
    # L->-240->T
    feats[388] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-240') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-240', node2.tup)
            
    # L->-90->T
    feats[389] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-90') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-90', node2.tup)
            
    # L->161->T
    feats[390] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '161') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '161', node2.tup)
            
    # L->-30->T
    feats[391] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-30') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-30', node2.tup)
            
    # L->10_pl->T
    feats[392] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '10_pl') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '10_pl', node2.tup)
            
    # L->27_du->T
    feats[393] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '27_du') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '27_du', node2.tup)
            
    # L->-141->T
    feats[394] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-141') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-141', node2.tup)
            
    # L->15_fp->T
    feats[395] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '15_fp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '15_fp', node2.tup)
            
    # L->13_fp->T
    feats[396] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '13_fp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '13_fp', node2.tup)
            
    # L->-64->T
    feats[397] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-64') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-64', node2.tup)
            
    # L->176->T
    feats[398] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '176') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '176', node2.tup)
            
    # L->100->T
    feats[399] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '100') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '100', node2.tup)
            
    # L->8_pl->T
    feats[400] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '8_pl') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '8_pl', node2.tup)
            
    # L->-28->T
    feats[401] = tryProb_catchZero(mat_lem2cng_countonly, mat_lemCount_1D, node1.lemma, '-28') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-28', node2.tup)
            
    # C->abl. du.->L
    feats[402] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'abl. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. du.', node2.lemma)
            
    # C->-150->L
    feats[403] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-150') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-150', node2.lemma)
            
    # C->16_sg->L
    feats[404] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '16_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_sg', node2.lemma)
            
    # C->fem->L
    feats[405] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'fem', node2.lemma)
            
    # C->-15->L
    feats[406] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-15') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-15', node2.lemma)
            
    # C->-159->L
    feats[407] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-159') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-159', node2.lemma)
            
    # C->3_tp->L
    feats[408] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '3_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_tp', node2.lemma)
            
    # C->39->L
    feats[409] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '39') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '39', node2.lemma)
            
    # C->voc. masc.->L
    feats[410] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'voc. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. masc.', node2.lemma)
            
    # C->-296->L
    feats[411] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-296') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-296', node2.lemma)
            
    # C->-94->L
    feats[412] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-94') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-94', node2.lemma)
            
    # C->-31->L
    feats[413] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-31') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-31', node2.lemma)
            
    # C->du->L
    feats[414] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'du', node2.lemma)
            
    # C->-122->L
    feats[415] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-122') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-122', node2.lemma)
            
    # C->13_sg->L
    feats[416] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '13_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_sg', node2.lemma)
            
    # C->-97->L
    feats[417] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-97') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-97', node2.lemma)
            
    # C->72->L
    feats[418] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '72') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '72', node2.lemma)
            
    # C->-241->L
    feats[419] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-241') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-241', node2.lemma)
            
    # C->16_pl->L
    feats[420] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '16_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_pl', node2.lemma)
            
    # C->12_fp->L
    feats[421] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '12_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_fp', node2.lemma)
            
    # C->137->L
    feats[422] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '137') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
            
    # C->76->L
    feats[423] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '76') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '76', node2.lemma)
            
    # C->-33->L
    feats[424] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-33') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-33', node2.lemma)
            
    # C->-169->L
    feats[425] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-169') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-169', node2.lemma)
            
    # C->nom. fem->L
    feats[426] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'nom. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. fem', node2.lemma)
            
    # C->128->L
    feats[427] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '128') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '128', node2.lemma)
            
    # C->-123->L
    feats[428] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-123') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-123', node2.lemma)
            
    # C->-11->L
    feats[429] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-11') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-11', node2.lemma)
            
    # C->loc. sg.->L
    feats[430] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'loc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc. sg.', node2.lemma)
            
    # C->-261->L
    feats[431] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-261') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-261', node2.lemma)
            
    # C->56->L
    feats[432] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '56') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '56', node2.lemma)
            
    # C->170->L
    feats[433] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '170') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '170', node2.lemma)
            
    # C->41->L
    feats[434] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '41') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '41', node2.lemma)
            
    # C->9_sp->L
    feats[435] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '9_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_sp', node2.lemma)
            
    # C->5_sp->L
    feats[436] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '5_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sp', node2.lemma)
            
    # C->-17->L
    feats[437] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-17') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-17', node2.lemma)
            
    # C->pl_fp->L
    feats[438] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'pl_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_fp', node2.lemma)
            
    # C->14_sp->L
    feats[439] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '14_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_sp', node2.lemma)
            
    # C->34->L
    feats[440] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '34') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '34', node2.lemma)
            
    # C->179->L
    feats[441] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '179') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '179', node2.lemma)
            
    # C->voc. fem->L
    feats[442] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'voc. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. fem', node2.lemma)
            
    # C->-152->L
    feats[443] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-152') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-152', node2.lemma)
            
    # C->acc. neutr.->L
    feats[444] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'acc. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. neutr.', node2.lemma)
            
    # C->tp->L
    feats[445] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'tp', node2.lemma)
            
    # C->10_fp->L
    feats[446] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '10_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_fp', node2.lemma)
            
    # C->132->L
    feats[447] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '132') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '132', node2.lemma)
            
    # C->-297->L
    feats[448] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-297') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-297', node2.lemma)
            
    # C->-149->L
    feats[449] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-149') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-149', node2.lemma)
            
    # C->10_sg->L
    feats[450] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '10_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_sg', node2.lemma)
            
    # C->-22->L
    feats[451] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-22') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-22', node2.lemma)
            
    # C->10_sp->L
    feats[452] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '10_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_sp', node2.lemma)
            
    # C->10_tp->L
    feats[453] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '10_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_tp', node2.lemma)
            
    # C->-77->L
    feats[454] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-77') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-77', node2.lemma)
            
    # C->-245->L
    feats[455] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-245') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-245', node2.lemma)
            
    # C->instr. du.->L
    feats[456] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'instr. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. du.', node2.lemma)
            
    # C->-166->L
    feats[457] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-166') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-166', node2.lemma)
            
    # C->dat->L
    feats[458] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'dat') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat', node2.lemma)
            
    # C->-139->L
    feats[459] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-139') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-139', node2.lemma)
            
    # C->151->L
    feats[460] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '151') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '151', node2.lemma)
            
    # C->98->L
    feats[461] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '98') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '98', node2.lemma)
            
    # C->2_fp->L
    feats[462] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '2_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_fp', node2.lemma)
            
    # C->5_pl->L
    feats[463] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '5_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_pl', node2.lemma)
            
    # C->12_tp->L
    feats[464] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '12_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_tp', node2.lemma)
            
    # C->-230->L
    feats[465] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-230') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-230', node2.lemma)
            
    # C->-34->L
    feats[466] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-34') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-34', node2.lemma)
            
    # C->29_sg->L
    feats[467] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '29_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '29_sg', node2.lemma)
            
    # C->-103->L
    feats[468] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-103') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-103', node2.lemma)
            
    # C->-132->L
    feats[469] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-132') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-132', node2.lemma)
            
    # C->-58->L
    feats[470] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-58') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-58', node2.lemma)
            
    # C->sg_fp->L
    feats[471] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'sg_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_fp', node2.lemma)
            
    # C->-121->L
    feats[472] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-121') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-121', node2.lemma)
            
    # C->acc. adj.->L
    feats[473] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'acc. adj.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. adj.', node2.lemma)
            
    # C->120->L
    feats[474] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '120') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '120', node2.lemma)
            
    # C->182->L
    feats[475] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '182') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '182', node2.lemma)
            
    # C->28_tp->L
    feats[476] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '28_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28_tp', node2.lemma)
            
    # C->-82->L
    feats[477] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-82') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-82', node2.lemma)
            
    # C->115->L
    feats[478] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '115') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '115', node2.lemma)
            
    # C->156->L
    feats[479] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '156') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '156', node2.lemma)
            
    # C->118->L
    feats[480] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '118') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '118', node2.lemma)
            
    # C->acc. fem->L
    feats[481] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'acc. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. fem', node2.lemma)
            
    # C->-266->L
    feats[482] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-266') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-266', node2.lemma)
            
    # C->loc->L
    feats[483] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'loc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc', node2.lemma)
            
    # C->149->L
    feats[484] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '149') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '149', node2.lemma)
            
    # C->121->L
    feats[485] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '121') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '121', node2.lemma)
            
    # C->140->L
    feats[486] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '140') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '140', node2.lemma)
            
    # C->30_sg->L
    feats[487] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '30_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_sg', node2.lemma)
            
    # C->-302->L
    feats[488] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-302') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-302', node2.lemma)
            
    # C->voc. sg.->L
    feats[489] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'voc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. sg.', node2.lemma)
            
    # C->fp->L
    feats[490] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'fp', node2.lemma)
            
    # C->-23->L
    feats[491] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-23') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-23', node2.lemma)
            
    # C->-104->L
    feats[492] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-104') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-104', node2.lemma)
            
    # C->2_du->L
    feats[493] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '2_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_du', node2.lemma)
            
    # C->-98->L
    feats[494] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-98') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-98', node2.lemma)
            
    # C->92->L
    feats[495] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '92') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '92', node2.lemma)
            
    # C->-21->L
    feats[496] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-21') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-21', node2.lemma)
            
    # C->pl->L
    feats[497] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
            
    # C->-279->L
    feats[498] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-279') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-279', node2.lemma)
            
    # C->-44->L
    feats[499] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-44') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-44', node2.lemma)
            
    # C->-115->L
    feats[500] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-115') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-115', node2.lemma)
            
    # C->12_pl->L
    feats[501] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '12_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_pl', node2.lemma)
            
    # C->-292->L
    feats[502] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-292') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-292', node2.lemma)
            
    # C->9_du->L
    feats[503] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '9_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_du', node2.lemma)
            
    # C->114->L
    feats[504] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '114') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '114', node2.lemma)
            
    # C->-36->L
    feats[505] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-36') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-36', node2.lemma)
            
    # C->-156->L
    feats[506] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-156') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-156', node2.lemma)
            
    # C->70->L
    feats[507] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '70') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '70', node2.lemma)
            
    # C->-27->L
    feats[508] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-27') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-27', node2.lemma)
            
    # C->-144->L
    feats[509] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-144') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-144', node2.lemma)
            
    # C->-111->L
    feats[510] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-111') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-111', node2.lemma)
            
    # C->-154->L
    feats[511] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-154') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-154', node2.lemma)
            
    # C->4_sp->L
    feats[512] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '4_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_sp', node2.lemma)
            
    # C->117->L
    feats[513] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '117') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '117', node2.lemma)
            
    # C->dat. du.->L
    feats[514] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'dat. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. du.', node2.lemma)
            
    # C->-307->L
    feats[515] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-307') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-307', node2.lemma)
            
    # C->-20->L
    feats[516] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-20') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-20', node2.lemma)
            
    # C->nom. du.->L
    feats[517] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'nom. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. du.', node2.lemma)
            
    # C->116->L
    feats[518] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '116') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '116', node2.lemma)
            
    # C->81->L
    feats[519] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '81') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '81', node2.lemma)
            
    # C->95->L
    feats[520] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '95') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '95', node2.lemma)
            
    # C->-24->L
    feats[521] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-24') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
            
    # C->12_sg->L
    feats[522] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '12_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_sg', node2.lemma)
            
    # C->-42->L
    feats[523] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-42') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-42', node2.lemma)
            
    # C->-51->L
    feats[524] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-51') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-51', node2.lemma)
            
    # C->abl. sg.->L
    feats[525] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'abl. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. sg.', node2.lemma)
            
    # C->-293->L
    feats[526] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-293') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-293', node2.lemma)
            
    # C->6_sp->L
    feats[527] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '6_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_sp', node2.lemma)
            
    # C->-63->L
    feats[528] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-63') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-63', node2.lemma)
            
    # C->acc. sg.->L
    feats[529] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'acc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. sg.', node2.lemma)
            
    # C->-89->L
    feats[530] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-89') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-89', node2.lemma)
            
    # C->sg_sp->L
    feats[531] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'sg_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_sp', node2.lemma)
            
    # C->nom. neutr.->L
    feats[532] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'nom. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom. neutr.', node2.lemma)
            
    # C->42->L
    feats[533] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '42') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '42', node2.lemma)
            
    # C->32->L
    feats[534] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '32') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '32', node2.lemma)
            
    # C->neutr->L
    feats[535] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'neutr') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'neutr', node2.lemma)
            
    # C->89->L
    feats[536] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '89') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '89', node2.lemma)
            
    # C->9_sg->L
    feats[537] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '9_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_sg', node2.lemma)
            
    # C->7_fp->L
    feats[538] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '7_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_fp', node2.lemma)
            
    # C->139->L
    feats[539] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '139') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '139', node2.lemma)
            
    # C->6_tp->L
    feats[540] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '6_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_tp', node2.lemma)
            
    # C->-301->L
    feats[541] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-301') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-301', node2.lemma)
            
    # C->73->L
    feats[542] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '73') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '73', node2.lemma)
            
    # C->-247->L
    feats[543] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-247') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-247', node2.lemma)
            
    # C->152->L
    feats[544] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '152') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '152', node2.lemma)
            
    # C->7_pl->L
    feats[545] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '7_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_pl', node2.lemma)
            
    # C->97->L
    feats[546] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '97') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '97', node2.lemma)
            
    # C->50->L
    feats[547] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '50') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '50', node2.lemma)
            
    # C->instr. masc.->L
    feats[548] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'instr. masc.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. masc.', node2.lemma)
            
    # C->119->L
    feats[549] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '119') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '119', node2.lemma)
            
    # C->-153->L
    feats[550] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-153') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-153', node2.lemma)
            
    # C->voc. neutr.->L
    feats[551] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'voc. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. neutr.', node2.lemma)
            
    # C->-67->L
    feats[552] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-67') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-67', node2.lemma)
            
    # C->-62->L
    feats[553] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-62') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-62', node2.lemma)
            
    # C->-76->L
    feats[554] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-76') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-76', node2.lemma)
            
    # C->51->L
    feats[555] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '51') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '51', node2.lemma)
            
    # C->36->L
    feats[556] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '36') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '36', node2.lemma)
            
    # C->11_du->L
    feats[557] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '11_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '11_du', node2.lemma)
            
    # C->8_fp->L
    feats[558] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '8_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
            
    # C->80->L
    feats[559] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '80') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '80', node2.lemma)
            
    # C->16_fp->L
    feats[560] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '16_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_fp', node2.lemma)
            
    # C->-133->L
    feats[561] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-133') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-133', node2.lemma)
            
    # C->173->L
    feats[562] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '173') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '173', node2.lemma)
            
    # C->169->L
    feats[563] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '169') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '169', node2.lemma)
            
    # C->-61->L
    feats[564] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-61') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-61', node2.lemma)
            
    # C->58->L
    feats[565] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '58') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '58', node2.lemma)
            
    # C->15_sg->L
    feats[566] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '15_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_sg', node2.lemma)
            
    # C->153->L
    feats[567] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '153') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '153', node2.lemma)
            
    # C->3_sg->L
    feats[568] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '3_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_sg', node2.lemma)
            
    # C->gen. sg.->L
    feats[569] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'gen. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen. sg.', node2.lemma)
            
    # C->-262->L
    feats[570] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-262') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-262', node2.lemma)
            
    # C->instr. fem->L
    feats[571] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'instr. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. fem', node2.lemma)
            
    # C->pl_tp->L
    feats[572] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'pl_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_tp', node2.lemma)
            
    # C->instr->L
    feats[573] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'instr') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr', node2.lemma)
            
    # C->8_sp->L
    feats[574] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '8_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_sp', node2.lemma)
            
    # C->-306->L
    feats[575] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-306') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-306', node2.lemma)
            
    # C->60->L
    feats[576] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '60') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '60', node2.lemma)
            
    # C->-119->L
    feats[577] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-119') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-119', node2.lemma)
            
    # C->14_pl->L
    feats[578] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '14_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_pl', node2.lemma)
            
    # C->voc->L
    feats[579] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'voc') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc', node2.lemma)
            
    # C->68->L
    feats[580] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '68') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '68', node2.lemma)
            
    # C->14_sg->L
    feats[581] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '14_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_sg', node2.lemma)
            
    # C->acc. du.->L
    feats[582] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'acc. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. du.', node2.lemma)
            
    # C->15_pl->L
    feats[583] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '15_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_pl', node2.lemma)
            
    # C->-83->L
    feats[584] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-83') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-83', node2.lemma)
            
    # C->30_tp->L
    feats[585] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '30_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_tp', node2.lemma)
            
    # C->-249->L
    feats[586] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-249') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-249', node2.lemma)
            
    # C->6_fp->L
    feats[587] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '6_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_fp', node2.lemma)
            
    # C->93->L
    feats[588] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '93') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '93', node2.lemma)
            
    # C->6_du->L
    feats[589] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '6_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_du', node2.lemma)
            
    # C->71->L
    feats[590] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '71') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '71', node2.lemma)
            
    # C->111->L
    feats[591] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '111') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '111', node2.lemma)
            
    # C->15_du->L
    feats[592] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '15_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_du', node2.lemma)
            
    # C->-29->L
    feats[593] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-29') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-29', node2.lemma)
            
    # C->10_du->L
    feats[594] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '10_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_du', node2.lemma)
            
    # C->dat. sg.->L
    feats[595] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'dat. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'dat. sg.', node2.lemma)
            
    # C->15_tp->L
    feats[596] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '15_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_tp', node2.lemma)
            
    # C->-32->L
    feats[597] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-32') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-32', node2.lemma)
            
    # C->instr. neutr.->L
    feats[598] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'instr. neutr.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. neutr.', node2.lemma)
            
    # C->-54->L
    feats[599] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-54') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-54', node2.lemma)
            
    # C->175->L
    feats[600] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '175') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '175', node2.lemma)
            
    # C->134->L
    feats[601] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '134') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '134', node2.lemma)
            
    # C->-81->L
    feats[602] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-81') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-81', node2.lemma)
            
    # C->-43->L
    feats[603] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-43') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-43', node2.lemma)
            
    # C->-38->L
    feats[604] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-38') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-38', node2.lemma)
            
    # C->-190->L
    feats[605] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-190') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-190', node2.lemma)
            
    # C->30_fp->L
    feats[606] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '30_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_fp', node2.lemma)
            
    # C->-299->L
    feats[607] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-299') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-299', node2.lemma)
            
    # C->4_sg->L
    feats[608] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '4_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_sg', node2.lemma)
            
    # C->2->L
    feats[609] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '2') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2', node2.lemma)
            
    # C->pl_sp->L
    feats[610] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'pl_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl_sp', node2.lemma)
            
    # C->27_pl->L
    feats[611] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '27_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_pl', node2.lemma)
            
    # C->8_sg->L
    feats[612] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '8_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_sg', node2.lemma)
            
    # C->69->L
    feats[613] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '69') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '69', node2.lemma)
            
    # C->155->L
    feats[614] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '155') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '155', node2.lemma)
            
    # C->154->L
    feats[615] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '154') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '154', node2.lemma)
            
    # C->82->L
    feats[616] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '82') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '82', node2.lemma)
            
    # C->-57->L
    feats[617] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-57') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-57', node2.lemma)
            
    # C->108->L
    feats[618] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '108') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '108', node2.lemma)
            
    # C->-87->L
    feats[619] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-87') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-87', node2.lemma)
            
    # C->180->L
    feats[620] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '180') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '180', node2.lemma)
            
    # C->2_tp->L
    feats[621] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '2_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_tp', node2.lemma)
            
    # C->voc. pl.->L
    feats[622] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'voc. pl.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. pl.', node2.lemma)
            
    # C->-35->L
    feats[623] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-35') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
            
    # C->9_pl->L
    feats[624] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '9_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '9_pl', node2.lemma)
            
    # C->91->L
    feats[625] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '91') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '91', node2.lemma)
            
    # C->4_tp->L
    feats[626] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '4_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '4_tp', node2.lemma)
            
    # C->-45->L
    feats[627] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-45') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-45', node2.lemma)
            
    # C->-308->L
    feats[628] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-308') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-308', node2.lemma)
            
    # C->130->L
    feats[629] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '130') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '130', node2.lemma)
            
    # C->-240->L
    feats[630] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-240') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-240', node2.lemma)
            
    # C->160->L
    feats[631] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '160') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '160', node2.lemma)
            
    # C->-90->L
    feats[632] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-90') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-90', node2.lemma)
            
    # C->161->L
    feats[633] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '161') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '161', node2.lemma)
            
    # C->162->L
    feats[634] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '162') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '162', node2.lemma)
            
    # C->172->L
    feats[635] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '172') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '172', node2.lemma)
            
    # C->-30->L
    feats[636] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-30') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-30', node2.lemma)
            
    # C->-16->L
    feats[637] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-16') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-16', node2.lemma)
            
    # C->10_pl->L
    feats[638] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '10_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_pl', node2.lemma)
            
    # C->122->L
    feats[639] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '122') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '122', node2.lemma)
            
    # C->-78->L
    feats[640] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-78') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-78', node2.lemma)
            
    # C->27_du->L
    feats[641] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '27_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_du', node2.lemma)
            
    # C->-72->L
    feats[642] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-72') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-72', node2.lemma)
            
    # C->-141->L
    feats[643] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-141') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-141', node2.lemma)
            
    # C->15_fp->L
    feats[644] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '15_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_fp', node2.lemma)
            
    # C->nom->L
    feats[645] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'nom') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'nom', node2.lemma)
            
    # C->13_fp->L
    feats[646] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '13_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_fp', node2.lemma)
            
    # C->-64->L
    feats[647] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-64') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-64', node2.lemma)
            
    # C->176->L
    feats[648] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '176') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '176', node2.lemma)
            
    # C->30_sp->L
    feats[649] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '30_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '30_sp', node2.lemma)
            
    # C->100->L
    feats[650] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '100') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '100', node2.lemma)
            
    # C->8_pl->L
    feats[651] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '8_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_pl', node2.lemma)
            
    # C->-28->L
    feats[652] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-28') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-28', node2.lemma)
            
    # C->-303->L
    feats[653] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-303') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-303', node2.lemma)
            
    # C->abl. du.->C
    feats[654] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'abl. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'abl. du.', node2.cng)
            
    # C->-150->C
    feats[655] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-150') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-150', node2.cng)
            
    # C->28_sg->C
    feats[656] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '28_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28_sg', node2.cng)
            
    # C->16_sg->C
    feats[657] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '16_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_sg', node2.cng)
            
    # C->abl. pl.->C
    feats[658] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'abl. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'abl. pl.', node2.cng)
            
    # C->fem->C
    feats[659] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'fem', node2.cng)
            
    # C->-15->C
    feats[660] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-15') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-15', node2.cng)
            
    # C->4_pl->C
    feats[661] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '4_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_pl', node2.cng)
            
    # C->-159->C
    feats[662] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-159') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-159', node2.cng)
            
    # C->3_tp->C
    feats[663] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '3_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_tp', node2.cng)
            
    # C->39->C
    feats[664] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '39', node2.cng)
            
    # C->voc. du.->C
    feats[665] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'voc. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. du.', node2.cng)
            
    # C->voc. masc.->C
    feats[666] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'voc. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. masc.', node2.cng)
            
    # C->-296->C
    feats[667] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-296') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-296', node2.cng)
            
    # C->-94->C
    feats[668] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-94') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-94', node2.cng)
            
    # C->-31->C
    feats[669] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-31') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-31', node2.cng)
            
    # C->du->C
    feats[670] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du', node2.cng)
            
    # C->-122->C
    feats[671] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-122') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-122', node2.cng)
            
    # C->13_sg->C
    feats[672] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '13_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_sg', node2.cng)
            
    # C->72->C
    feats[673] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '72') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '72', node2.cng)
            
    # C->-241->C
    feats[674] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-241') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-241', node2.cng)
            
    # C->du_fp->C
    feats[675] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'du_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'du_fp', node2.cng)
            
    # C->16_pl->C
    feats[676] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '16_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_pl', node2.cng)
            
    # C->-309->C
    feats[677] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-309') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-309', node2.cng)
            
    # C->12_fp->C
    feats[678] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '12_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_fp', node2.cng)
            
    # C->-220->C
    feats[679] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-220') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-220', node2.cng)
            
    # C->137->C
    feats[680] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '137') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '137', node2.cng)
            
    # C->76->C
    feats[681] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '76') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '76', node2.cng)
            
    # C->-33->C
    feats[682] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-33') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-33', node2.cng)
            
    # C->-169->C
    feats[683] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-169') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-169', node2.cng)
            
    # C->nom. fem->C
    feats[684] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'nom. fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. fem', node2.cng)
            
    # C->128->C
    feats[685] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '128') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '128', node2.cng)
            
    # C->-123->C
    feats[686] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-123') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-123', node2.cng)
            
    # C->-11->C
    feats[687] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-11') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-11', node2.cng)
            
    # C->loc. sg.->C
    feats[688] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'loc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc. sg.', node2.cng)
            
    # C->-261->C
    feats[689] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-261') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-261', node2.cng)
            
    # C->56->C
    feats[690] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '56') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '56', node2.cng)
            
    # C->170->C
    feats[691] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '170') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '170', node2.cng)
            
    # C->41->C
    feats[692] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '41') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '41', node2.cng)
            
    # C->9_sp->C
    feats[693] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '9_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sp', node2.cng)
            
    # C->5_sp->C
    feats[694] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '5_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', node2.cng)
            
    # C->54->C
    feats[695] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '54') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '54', node2.cng)
            
    # C->-17->C
    feats[696] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-17') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-17', node2.cng)
            
    # C->150->C
    feats[697] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '150') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '150', node2.cng)
            
    # C->pl_fp->C
    feats[698] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'pl_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_fp', node2.cng)
            
    # C->14_sp->C
    feats[699] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '14_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_sp', node2.cng)
            
    # C->79->C
    feats[700] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '79') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '79', node2.cng)
            
    # C->34->C
    feats[701] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '34', node2.cng)
            
    # C->179->C
    feats[702] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '179') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '179', node2.cng)
            
    # C->6_pl->C
    feats[703] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '6_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_pl', node2.cng)
            
    # C->-273->C
    feats[704] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-273') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-273', node2.cng)
            
    # C->voc. fem->C
    feats[705] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'voc. fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. fem', node2.cng)
            
    # C->138->C
    feats[706] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '138') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '138', node2.cng)
            
    # C->-152->C
    feats[707] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-152') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-152', node2.cng)
            
    # C->acc. neutr.->C
    feats[708] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'acc. neutr.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. neutr.', node2.cng)
            
    # C->tp->C
    feats[709] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'tp', node2.cng)
            
    # C->10_fp->C
    feats[710] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '10_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_fp', node2.cng)
            
    # C->132->C
    feats[711] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '132') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '132', node2.cng)
            
    # C->-297->C
    feats[712] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-297') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-297', node2.cng)
            
    # C->-142->C
    feats[713] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-142') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-142', node2.cng)
            
    # C->-291->C
    feats[714] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-291') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-291', node2.cng)
            
    # C->-149->C
    feats[715] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-149', node2.cng)
            
    # C->10_sg->C
    feats[716] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '10_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sg', node2.cng)
            
    # C->-22->C
    feats[717] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-22') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-22', node2.cng)
            
    # C->10_sp->C
    feats[718] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '10_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_sp', node2.cng)
            
    # C->10_tp->C
    feats[719] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '10_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_tp', node2.cng)
            
    # C->-77->C
    feats[720] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-77') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-77', node2.cng)
            
    # C->-245->C
    feats[721] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-245') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-245', node2.cng)
            
    # C->instr. du.->C
    feats[722] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'instr. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. du.', node2.cng)
            
    # C->instr. pl.->C
    feats[723] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'instr. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. pl.', node2.cng)
            
    # C->-166->C
    feats[724] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-166') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-166', node2.cng)
            
    # C->dat->C
    feats[725] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'dat') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat', node2.cng)
            
    # C->-25->C
    feats[726] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-25') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-25', node2.cng)
            
    # C->-139->C
    feats[727] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-139', node2.cng)
            
    # C->151->C
    feats[728] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '151') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '151', node2.cng)
            
    # C->98->C
    feats[729] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '98', node2.cng)
            
    # C->2_fp->C
    feats[730] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', node2.cng)
            
    # C->5_pl->C
    feats[731] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '5_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_pl', node2.cng)
            
    # C->12_tp->C
    feats[732] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '12_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_tp', node2.cng)
            
    # C->-230->C
    feats[733] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', node2.cng)
            
    # C->-34->C
    feats[734] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-34', node2.cng)
            
    # C->29_sg->C
    feats[735] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '29_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '29_sg', node2.cng)
            
    # C->-103->C
    feats[736] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-103') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-103', node2.cng)
            
    # C->-132->C
    feats[737] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-132') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-132', node2.cng)
            
    # C->-58->C
    feats[738] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-58', node2.cng)
            
    # C->-52->C
    feats[739] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-52') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-52', node2.cng)
            
    # C->sg_fp->C
    feats[740] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'sg_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg_fp', node2.cng)
            
    # C->90->C
    feats[741] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '90', node2.cng)
            
    # C->-121->C
    feats[742] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-121') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-121', node2.cng)
            
    # C->acc. adj.->C
    feats[743] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'acc. adj.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. adj.', node2.cng)
            
    # C->120->C
    feats[744] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '120') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '120', node2.cng)
            
    # C->182->C
    feats[745] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '182') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '182', node2.cng)
            
    # C->28_tp->C
    feats[746] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '28_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28_tp', node2.cng)
            
    # C->-82->C
    feats[747] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-82') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-82', node2.cng)
            
    # C->115->C
    feats[748] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '115') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '115', node2.cng)
            
    # C->156->C
    feats[749] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '156', node2.cng)
            
    # C->118->C
    feats[750] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '118') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '118', node2.cng)
            
    # C->-112->C
    feats[751] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-112') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-112', node2.cng)
            
    # C->acc. fem->C
    feats[752] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'acc. fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. fem', node2.cng)
            
    # C->-266->C
    feats[753] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-266') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-266', node2.cng)
            
    # C->loc->C
    feats[754] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'loc') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'loc', node2.cng)
            
    # C->149->C
    feats[755] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '149', node2.cng)
            
    # C->121->C
    feats[756] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '121') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '121', node2.cng)
            
    # C->140->C
    feats[757] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '140') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '140', node2.cng)
            
    # C->30_sg->C
    feats[758] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '30_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_sg', node2.cng)
            
    # C->-302->C
    feats[759] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', node2.cng)
            
    # C->voc. sg.->C
    feats[760] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'voc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. sg.', node2.cng)
            
    # C->136->C
    feats[761] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '136') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '136', node2.cng)
            
    # C->171->C
    feats[762] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '171') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '171', node2.cng)
            
    # C->fp->C
    feats[763] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'fp', node2.cng)
            
    # C->-23->C
    feats[764] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-23') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-23', node2.cng)
            
    # C->-104->C
    feats[765] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', node2.cng)
            
    # C->2_du->C
    feats[766] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', node2.cng)
            
    # C->-98->C
    feats[767] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-98') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-98', node2.cng)
            
    # C->92->C
    feats[768] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '92') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '92', node2.cng)
            
    # C->-21->C
    feats[769] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-21') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-21', node2.cng)
            
    # C->pl->C
    feats[770] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl', node2.cng)
            
    # C->-73->C
    feats[771] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-73') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-73', node2.cng)
            
    # C->-279->C
    feats[772] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-279') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-279', node2.cng)
            
    # C->-44->C
    feats[773] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-44') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-44', node2.cng)
            
    # C->-115->C
    feats[774] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-115') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-115', node2.cng)
            
    # C->12_pl->C
    feats[775] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '12_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_pl', node2.cng)
            
    # C->-292->C
    feats[776] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-292') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-292', node2.cng)
            
    # C->9_du->C
    feats[777] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '9_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_du', node2.cng)
            
    # C->114->C
    feats[778] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '114') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '114', node2.cng)
            
    # C->-36->C
    feats[779] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-36') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-36', node2.cng)
            
    # C->-156->C
    feats[780] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-156', node2.cng)
            
    # C->70->C
    feats[781] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '70') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '70', node2.cng)
            
    # C->-27->C
    feats[782] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-27') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-27', node2.cng)
            
    # C->-53->C
    feats[783] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-53') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-53', node2.cng)
            
    # C->-144->C
    feats[784] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-144') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-144', node2.cng)
            
    # C->-111->C
    feats[785] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-111') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-111', node2.cng)
            
    # C->-154->C
    feats[786] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-154') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-154', node2.cng)
            
    # C->-242->C
    feats[787] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-242') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-242', node2.cng)
            
    # C->-113->C
    feats[788] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-113') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-113', node2.cng)
            
    # C->4_sp->C
    feats[789] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '4_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_sp', node2.cng)
            
    # C->117->C
    feats[790] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '117') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '117', node2.cng)
            
    # C->dat. du.->C
    feats[791] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'dat. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. du.', node2.cng)
            
    # C->-143->C
    feats[792] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-143') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-143', node2.cng)
            
    # C->135->C
    feats[793] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '135') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '135', node2.cng)
            
    # C->gen->C
    feats[794] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'gen') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen', node2.cng)
            
    # C->-307->C
    feats[795] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-307') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-307', node2.cng)
            
    # C->-20->C
    feats[796] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-20') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-20', node2.cng)
            
    # C->8_tp->C
    feats[797] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '8_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_tp', node2.cng)
            
    # C->nom. du.->C
    feats[798] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'nom. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. du.', node2.cng)
            
    # C->116->C
    feats[799] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '116') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '116', node2.cng)
            
    # C->81->C
    feats[800] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '81', node2.cng)
            
    # C->95->C
    feats[801] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '95') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '95', node2.cng)
            
    # C->-24->C
    feats[802] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-24') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-24', node2.cng)
            
    # C->12_sg->C
    feats[803] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '12_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sg', node2.cng)
            
    # C->-42->C
    feats[804] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-42') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-42', node2.cng)
            
    # C->7_sg->C
    feats[805] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '7_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_sg', node2.cng)
            
    # C->-50->C
    feats[806] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-50') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-50', node2.cng)
            
    # C->-51->C
    feats[807] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-51') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-51', node2.cng)
            
    # C->abl. sg.->C
    feats[808] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'abl. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'abl. sg.', node2.cng)
            
    # C->acc. masc.->C
    feats[809] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'acc. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. masc.', node2.cng)
            
    # C->-293->C
    feats[810] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', node2.cng)
            
    # C->6_sp->C
    feats[811] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '6_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sp', node2.cng)
            
    # C->-63->C
    feats[812] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-63') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-63', node2.cng)
            
    # C->acc. sg.->C
    feats[813] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'acc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. sg.', node2.cng)
            
    # C->-89->C
    feats[814] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-89') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-89', node2.cng)
            
    # C->sg_sp->C
    feats[815] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'sg_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg_sp', node2.cng)
            
    # C->5_tp->C
    feats[816] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '5_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_tp', node2.cng)
            
    # C->nom. neutr.->C
    feats[817] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'nom. neutr.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. neutr.', node2.cng)
            
    # C->27_tp->C
    feats[818] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '27_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_tp', node2.cng)
            
    # C->42->C
    feats[819] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '42') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '42', node2.cng)
            
    # C->32->C
    feats[820] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '32') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '32', node2.cng)
            
    # C->11_tp->C
    feats[821] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '11_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_tp', node2.cng)
            
    # C->-47->C
    feats[822] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-47') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-47', node2.cng)
            
    # C->15_sp->C
    feats[823] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '15_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_sp', node2.cng)
            
    # C->89->C
    feats[824] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '89') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '89', node2.cng)
            
    # C->9_sg->C
    feats[825] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '9_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_sg', node2.cng)
            
    # C->7_fp->C
    feats[826] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '7_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_fp', node2.cng)
            
    # C->5_fp->C
    feats[827] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '5_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_fp', node2.cng)
            
    # C->139->C
    feats[828] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '139') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '139', node2.cng)
            
    # C->6_tp->C
    feats[829] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '6_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_tp', node2.cng)
            
    # C->-301->C
    feats[830] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-301') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-301', node2.cng)
            
    # C->73->C
    feats[831] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '73') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '73', node2.cng)
            
    # C->-247->C
    feats[832] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-247') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-247', node2.cng)
            
    # C->152->C
    feats[833] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '152') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '152', node2.cng)
            
    # C->7_pl->C
    feats[834] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '7_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_pl', node2.cng)
            
    # C->97->C
    feats[835] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '97') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '97', node2.cng)
            
    # C->-39->C
    feats[836] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-39') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-39', node2.cng)
            
    # C->50->C
    feats[837] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '50') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '50', node2.cng)
            
    # C->instr. masc.->C
    feats[838] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'instr. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. masc.', node2.cng)
            
    # C->119->C
    feats[839] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '119') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '119', node2.cng)
            
    # C->-153->C
    feats[840] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-153') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-153', node2.cng)
            
    # C->voc. neutr.->C
    feats[841] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'voc. neutr.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. neutr.', node2.cng)
            
    # C->-67->C
    feats[842] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', node2.cng)
            
    # C->-62->C
    feats[843] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-62') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-62', node2.cng)
            
    # C->-76->C
    feats[844] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-76') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-76', node2.cng)
            
    # C->51->C
    feats[845] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '51') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '51', node2.cng)
            
    # C->36->C
    feats[846] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '36') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '36', node2.cng)
            
    # C->-269->C
    feats[847] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-269') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-269', node2.cng)
            
    # C->11_du->C
    feats[848] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '11_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '11_du', node2.cng)
            
    # C->8_fp->C
    feats[849] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '8_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_fp', node2.cng)
            
    # C->80->C
    feats[850] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '80') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '80', node2.cng)
            
    # C->16_fp->C
    feats[851] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '16_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_fp', node2.cng)
            
    # C->-133->C
    feats[852] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-133') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-133', node2.cng)
            
    # C->173->C
    feats[853] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '173') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '173', node2.cng)
            
    # C->169->C
    feats[854] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '169') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '169', node2.cng)
            
    # C->-61->C
    feats[855] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-61') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-61', node2.cng)
            
    # C->-96->C
    feats[856] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-96') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-96', node2.cng)
            
    # C->58->C
    feats[857] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '58', node2.cng)
            
    # C->15_sg->C
    feats[858] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '15_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_sg', node2.cng)
            
    # C->153->C
    feats[859] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '153') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '153', node2.cng)
            
    # C->3_sg->C
    feats[860] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '3_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_sg', node2.cng)
            
    # C->178->C
    feats[861] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '178') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '178', node2.cng)
            
    # C->gen. sg.->C
    feats[862] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'gen. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. sg.', node2.cng)
            
    # C->-262->C
    feats[863] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-262') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-262', node2.cng)
            
    # C->sp->C
    feats[864] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sp', node2.cng)
            
    # C->instr. fem->C
    feats[865] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'instr. fem') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. fem', node2.cng)
            
    # C->pl_tp->C
    feats[866] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'pl_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_tp', node2.cng)
            
    # C->102->C
    feats[867] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '102') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '102', node2.cng)
            
    # C->instr->C
    feats[868] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'instr') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr', node2.cng)
            
    # C->8_sp->C
    feats[869] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '8_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_sp', node2.cng)
            
    # C->38->C
    feats[870] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '38') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '38', node2.cng)
            
    # C->-306->C
    feats[871] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-306') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-306', node2.cng)
            
    # C->148->C
    feats[872] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '148') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '148', node2.cng)
            
    # C->174->C
    feats[873] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '174') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '174', node2.cng)
            
    # C->112->C
    feats[874] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '112') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '112', node2.cng)
            
    # C->60->C
    feats[875] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '60') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '60', node2.cng)
            
    # C->-119->C
    feats[876] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-119') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-119', node2.cng)
            
    # C->14_pl->C
    feats[877] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '14_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_pl', node2.cng)
            
    # C->181->C
    feats[878] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '181') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '181', node2.cng)
            
    # C->voc->C
    feats[879] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'voc') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc', node2.cng)
            
    # C->-210->C
    feats[880] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-210') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-210', node2.cng)
            
    # C->-55->C
    feats[881] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-55') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-55', node2.cng)
            
    # C->68->C
    feats[882] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '68') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '68', node2.cng)
            
    # C->14_sg->C
    feats[883] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '14_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_sg', node2.cng)
            
    # C->acc. du.->C
    feats[884] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'acc. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. du.', node2.cng)
            
    # C->15_pl->C
    feats[885] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '15_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_pl', node2.cng)
            
    # C->-83->C
    feats[886] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-83') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-83', node2.cng)
            
    # C->30_tp->C
    feats[887] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '30_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_tp', node2.cng)
            
    # C->-19->C
    feats[888] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-19') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-19', node2.cng)
            
    # C->-249->C
    feats[889] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-249') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-249', node2.cng)
            
    # C->6_fp->C
    feats[890] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '6_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_fp', node2.cng)
            
    # C->13_pl->C
    feats[891] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '13_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_pl', node2.cng)
            
    # C->93->C
    feats[892] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '93') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '93', node2.cng)
            
    # C->16_tp->C
    feats[893] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '16_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_tp', node2.cng)
            
    # C->6_du->C
    feats[894] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '6_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_du', node2.cng)
            
    # C->71->C
    feats[895] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '71') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '71', node2.cng)
            
    # C->111->C
    feats[896] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '111') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '111', node2.cng)
            
    # C->15_du->C
    feats[897] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '15_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_du', node2.cng)
            
    # C->10_du->C
    feats[898] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '10_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_du', node2.cng)
            
    # C->dat. sg.->C
    feats[899] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'dat. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. sg.', node2.cng)
            
    # C->15_tp->C
    feats[900] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '15_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_tp', node2.cng)
            
    # C->-32->C
    feats[901] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-32') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-32', node2.cng)
            
    # C->instr. neutr.->C
    feats[902] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'instr. neutr.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. neutr.', node2.cng)
            
    # C->-54->C
    feats[903] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-54') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-54', node2.cng)
            
    # C->175->C
    feats[904] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '175') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '175', node2.cng)
            
    # C->masc->C
    feats[905] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'masc') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'masc', node2.cng)
            
    # C->-41->C
    feats[906] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-41') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-41', node2.cng)
            
    # C->134->C
    feats[907] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '134') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '134', node2.cng)
            
    # C->-81->C
    feats[908] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-81', node2.cng)
            
    # C->-43->C
    feats[909] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-43') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-43', node2.cng)
            
    # C->-38->C
    feats[910] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-38') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-38', node2.cng)
            
    # C->-190->C
    feats[911] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-190') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-190', node2.cng)
            
    # C->nom. pl.->C
    feats[912] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'nom. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. pl.', node2.cng)
            
    # C->30_fp->C
    feats[913] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '30_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_fp', node2.cng)
            
    # C->-299->C
    feats[914] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-299') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-299', node2.cng)
            
    # C->4_sg->C
    feats[915] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '4_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_sg', node2.cng)
            
    # C->35->C
    feats[916] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '35', node2.cng)
            
    # C->-18->C
    feats[917] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-18') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-18', node2.cng)
            
    # C->2->C
    feats[918] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '2') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2', node2.cng)
            
    # C->61->C
    feats[919] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '61') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '61', node2.cng)
            
    # C->7_sp->C
    feats[920] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '7_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_sp', node2.cng)
            
    # C->pl_sp->C
    feats[921] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'pl_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl_sp', node2.cng)
            
    # C->27_pl->C
    feats[922] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '27_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_pl', node2.cng)
            
    # C->31->C
    feats[923] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '31') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '31', node2.cng)
            
    # C->8_sg->C
    feats[924] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '8_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_sg', node2.cng)
            
    # C->69->C
    feats[925] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '69') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '69', node2.cng)
            
    # C->155->C
    feats[926] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '155') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '155', node2.cng)
            
    # C->154->C
    feats[927] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '154') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '154', node2.cng)
            
    # C->82->C
    feats[928] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '82') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '82', node2.cng)
            
    # C->168->C
    feats[929] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '168') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '168', node2.cng)
            
    # C->-57->C
    feats[930] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-57') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-57', node2.cng)
            
    # C->108->C
    feats[931] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '108') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '108', node2.cng)
            
    # C->-87->C
    feats[932] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-87') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-87', node2.cng)
            
    # C->180->C
    feats[933] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '180') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '180', node2.cng)
            
    # C->2_tp->C
    feats[934] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '2_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_tp', node2.cng)
            
    # C->141->C
    feats[935] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '141') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '141', node2.cng)
            
    # C->-56->C
    feats[936] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-56') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-56', node2.cng)
            
    # C->voc. pl.->C
    feats[937] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'voc. pl.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'voc. pl.', node2.cng)
            
    # C->-35->C
    feats[938] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-35') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-35', node2.cng)
            
    # C->-117->C
    feats[939] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-117') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-117', node2.cng)
            
    # C->9_pl->C
    feats[940] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '9_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '9_pl', node2.cng)
            
    # C->91->C
    feats[941] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '91') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '91', node2.cng)
            
    # C->4_tp->C
    feats[942] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '4_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '4_tp', node2.cng)
            
    # C->-45->C
    feats[943] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-45') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-45', node2.cng)
            
    # C->6_sg->C
    feats[944] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '6_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sg', node2.cng)
            
    # C->29->C
    feats[945] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '29') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '29', node2.cng)
            
    # C->-308->C
    feats[946] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-308') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-308', node2.cng)
            
    # C->-37->C
    feats[947] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-37') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-37', node2.cng)
            
    # C->131->C
    feats[948] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '131') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '131', node2.cng)
            
    # C->130->C
    feats[949] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '130') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '130', node2.cng)
            
    # C->-68->C
    feats[950] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-68') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-68', node2.cng)
            
    # C->-240->C
    feats[951] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-240') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-240', node2.cng)
            
    # C->160->C
    feats[952] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '160') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '160', node2.cng)
            
    # C->-91->C
    feats[953] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-91') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-91', node2.cng)
            
    # C->48->C
    feats[954] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '48') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '48', node2.cng)
            
    # C->30_du->C
    feats[955] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '30_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '30_du', node2.cng)
            
    # C->-90->C
    feats[956] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-90', node2.cng)
            
    # C->161->C
    feats[957] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '161') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '161', node2.cng)
            
    # C->162->C
    feats[958] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '162') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '162', node2.cng)
            
    # C->172->C
    feats[959] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '172') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '172', node2.cng)
            
    # C->-30->C
    feats[960] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-30') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-30', node2.cng)
            
    # C->-16->C
    feats[961] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-16') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-16', node2.cng)
            
    # C->10_pl->C
    feats[962] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '10_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_pl', node2.cng)
            
    # C->55->C
    feats[963] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '55') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '55', node2.cng)
            
    # C->-78->C
    feats[964] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-78') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-78', node2.cng)
            
    # C->27_du->C
    feats[965] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', node2.cng)
            
    # C->12_sp->C
    feats[966] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '12_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_sp', node2.cng)
            
    # C->-141->C
    feats[967] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-141') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-141', node2.cng)
            
    # C->2_sp->C
    feats[968] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '2_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_sp', node2.cng)
            
    # C->15_fp->C
    feats[969] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '15_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_fp', node2.cng)
            
    # C->nom->C
    feats[970] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'nom') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom', node2.cng)
            
    # C->13_fp->C
    feats[971] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '13_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_fp', node2.cng)
            
    # C->-64->C
    feats[972] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-64') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-64', node2.cng)
            
    # C->176->C
    feats[973] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '176') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '176', node2.cng)
            
    # C->-129->C
    feats[974] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-129') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-129', node2.cng)
            
    # C->100->C
    feats[975] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '100') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '100', node2.cng)
            
    # C->nom. masc.->C
    feats[976] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'nom. masc.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'nom. masc.', node2.cng)
            
    # C->8_pl->C
    feats[977] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '8_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_pl', node2.cng)
            
    # C->-28->C
    feats[978] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-28') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-28', node2.cng)
            
    # C->-303->C
    feats[979] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-303') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-303', node2.cng)
            
    # C->abl. du.->T
    feats[980] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'abl. du.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'abl. du.', node2.tup)
            
    # C->-150->T
    feats[981] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-150') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-150', node2.tup)
            
    # C->16_sg->T
    feats[982] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '16_sg') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '16_sg', node2.tup)
            
    # C->-159->T
    feats[983] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-159') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-159', node2.tup)
            
    # C->3_tp->T
    feats[984] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '3_tp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '3_tp', node2.tup)
            
    # C->39->T
    feats[985] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '39') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '39', node2.tup)
            
    # C->voc. masc.->T
    feats[986] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'voc. masc.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'voc. masc.', node2.tup)
            
    # C->-296->T
    feats[987] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-296') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-296', node2.tup)
            
    # C->-94->T
    feats[988] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-94') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-94', node2.tup)
            
    # C->-31->T
    feats[989] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-31') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-31', node2.tup)
            
    # C->du->T
    feats[990] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'du') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'du', node2.tup)
            
    # C->-122->T
    feats[991] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-122') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-122', node2.tup)
            
    # C->13_sg->T
    feats[992] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '13_sg') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '13_sg', node2.tup)
            
    # C->-97->T
    feats[993] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-97') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-97', node2.tup)
            
    # C->72->T
    feats[994] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '72') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '72', node2.tup)
            
    # C->-241->T
    feats[995] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-241') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-241', node2.tup)
            
    # C->12_fp->T
    feats[996] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '12_fp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '12_fp', node2.tup)
            
    # C->137->T
    feats[997] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '137') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '137', node2.tup)
            
    # C->76->T
    feats[998] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '76') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '76', node2.tup)
            
    # C->-33->T
    feats[999] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-33') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-33', node2.tup)
            
    # C->-169->T
    feats[1000] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-169') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-169', node2.tup)
            
    # C->nom. fem->T
    feats[1001] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'nom. fem') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'nom. fem', node2.tup)
            
    # C->128->T
    feats[1002] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '128') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '128', node2.tup)
            
    # C->-123->T
    feats[1003] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-123') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-123', node2.tup)
            
    # C->-11->T
    feats[1004] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-11') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-11', node2.tup)
            
    # C->loc. sg.->T
    feats[1005] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'loc. sg.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'loc. sg.', node2.tup)
            
    # C->56->T
    feats[1006] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '56') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '56', node2.tup)
            
    # C->170->T
    feats[1007] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '170') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '170', node2.tup)
            
    # C->41->T
    feats[1008] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '41') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '41', node2.tup)
            
    # C->9_sp->T
    feats[1009] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '9_sp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '9_sp', node2.tup)
            
    # C->5_sp->T
    feats[1010] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '5_sp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '5_sp', node2.tup)
            
    # C->-17->T
    feats[1011] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-17') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-17', node2.tup)
            
    # C->pl_fp->T
    feats[1012] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'pl_fp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'pl_fp', node2.tup)
            
    # C->14_sp->T
    feats[1013] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '14_sp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '14_sp', node2.tup)
            
    # C->34->T
    feats[1014] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '34') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '34', node2.tup)
            
    # C->179->T
    feats[1015] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '179') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '179', node2.tup)
            
    # C->-152->T
    feats[1016] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-152') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-152', node2.tup)
            
    # C->acc. neutr.->T
    feats[1017] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'acc. neutr.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'acc. neutr.', node2.tup)
            
    # C->tp->T
    feats[1018] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'tp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'tp', node2.tup)
            
    # C->10_fp->T
    feats[1019] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '10_fp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '10_fp', node2.tup)
            
    # C->132->T
    feats[1020] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '132') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '132', node2.tup)
            
    # C->-297->T
    feats[1021] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-297') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-297', node2.tup)
            
    # C->-149->T
    feats[1022] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-149') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-149', node2.tup)
            
    # C->10_sg->T
    feats[1023] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '10_sg') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '10_sg', node2.tup)
            
    # C->10_sp->T
    feats[1024] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '10_sp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '10_sp', node2.tup)
            
    # C->10_tp->T
    feats[1025] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '10_tp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '10_tp', node2.tup)
            
    # C->-77->T
    feats[1026] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-77') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-77', node2.tup)
            
    # C->-245->T
    feats[1027] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-245') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-245', node2.tup)
            
    # C->instr. du.->T
    feats[1028] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'instr. du.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'instr. du.', node2.tup)
            
    # C->-166->T
    feats[1029] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-166') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-166', node2.tup)
            
    # C->dat->T
    feats[1030] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'dat') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'dat', node2.tup)
            
    # C->-139->T
    feats[1031] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-139') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-139', node2.tup)
            
    # C->151->T
    feats[1032] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '151') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '151', node2.tup)
            
    # C->98->T
    feats[1033] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '98') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '98', node2.tup)
            
    # C->177->T
    feats[1034] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '177') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '177', node2.tup)
            
    # C->2_fp->T
    feats[1035] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '2_fp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '2_fp', node2.tup)
            
    # C->5_pl->T
    feats[1036] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '5_pl') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '5_pl', node2.tup)
            
    # C->12_tp->T
    feats[1037] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '12_tp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '12_tp', node2.tup)
            
    # C->-230->T
    feats[1038] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-230') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-230', node2.tup)
            
    # C->-34->T
    feats[1039] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-34') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-34', node2.tup)
            
    # C->29_sg->T
    feats[1040] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '29_sg') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '29_sg', node2.tup)
            
    # C->-103->T
    feats[1041] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-103') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-103', node2.tup)
            
    # C->-132->T
    feats[1042] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-132') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-132', node2.tup)
            
    # C->sg_fp->T
    feats[1043] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'sg_fp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'sg_fp', node2.tup)
            
    # C->90->T
    feats[1044] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '90') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '90', node2.tup)
            
    # C->-121->T
    feats[1045] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-121') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-121', node2.tup)
            
    # C->acc. adj.->T
    feats[1046] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'acc. adj.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'acc. adj.', node2.tup)
            
    # C->120->T
    feats[1047] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '120') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '120', node2.tup)
            
    # C->182->T
    feats[1048] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '182') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '182', node2.tup)
            
    # C->28_tp->T
    feats[1049] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '28_tp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '28_tp', node2.tup)
            
    # C->-82->T
    feats[1050] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-82') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-82', node2.tup)
            
    # C->115->T
    feats[1051] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '115') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '115', node2.tup)
            
    # C->156->T
    feats[1052] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '156') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '156', node2.tup)
            
    # C->118->T
    feats[1053] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '118') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '118', node2.tup)
            
    # C->acc. fem->T
    feats[1054] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'acc. fem') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'acc. fem', node2.tup)
            
    # C->-266->T
    feats[1055] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-266') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-266', node2.tup)
            
    # C->loc->T
    feats[1056] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'loc') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'loc', node2.tup)
            
    # C->121->T
    feats[1057] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '121') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '121', node2.tup)
            
    # C->140->T
    feats[1058] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '140') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '140', node2.tup)
            
    # C->30_sg->T
    feats[1059] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '30_sg') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '30_sg', node2.tup)
            
    # C->-302->T
    feats[1060] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-302') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-302', node2.tup)
            
    # C->voc. sg.->T
    feats[1061] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'voc. sg.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'voc. sg.', node2.tup)
            
    # C->fp->T
    feats[1062] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'fp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'fp', node2.tup)
            
    # C->-23->T
    feats[1063] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-23') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-23', node2.tup)
            
    # C->-104->T
    feats[1064] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-104') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-104', node2.tup)
            
    # C->2_du->T
    feats[1065] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '2_du') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '2_du', node2.tup)
            
    # C->-98->T
    feats[1066] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-98') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-98', node2.tup)
            
    # C->92->T
    feats[1067] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '92') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '92', node2.tup)
            
    # C->-21->T
    feats[1068] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-21') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-21', node2.tup)
            
    # C->pl->T
    feats[1069] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'pl') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'pl', node2.tup)
            
    # C->-279->T
    feats[1070] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-279') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-279', node2.tup)
            
    # C->-44->T
    feats[1071] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-44') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-44', node2.tup)
            
    # C->-115->T
    feats[1072] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-115') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-115', node2.tup)
            
    # C->12_pl->T
    feats[1073] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '12_pl') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '12_pl', node2.tup)
            
    # C->-292->T
    feats[1074] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-292') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-292', node2.tup)
            
    # C->114->T
    feats[1075] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '114') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '114', node2.tup)
            
    # C->-36->T
    feats[1076] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-36') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-36', node2.tup)
            
    # C->-156->T
    feats[1077] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-156') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-156', node2.tup)
            
    # C->70->T
    feats[1078] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '70') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '70', node2.tup)
            
    # C->-27->T
    feats[1079] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-27') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-27', node2.tup)
            
    # C->-144->T
    feats[1080] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-144') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-144', node2.tup)
            
    # C->-92->T
    feats[1081] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-92') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-92', node2.tup)
            
    # C->-111->T
    feats[1082] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-111') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-111', node2.tup)
            
    # C->-154->T
    feats[1083] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-154') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-154', node2.tup)
            
    # C->4_sp->T
    feats[1084] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '4_sp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '4_sp', node2.tup)
            
    # C->117->T
    feats[1085] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '117') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '117', node2.tup)
            
    # C->dat. du.->T
    feats[1086] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'dat. du.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'dat. du.', node2.tup)
            
    # C->-307->T
    feats[1087] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-307') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-307', node2.tup)
            
    # C->nom. du.->T
    feats[1088] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'nom. du.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'nom. du.', node2.tup)
            
    # C->116->T
    feats[1089] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '116') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '116', node2.tup)
            
    # C->81->T
    feats[1090] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '81') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '81', node2.tup)
            
    # C->95->T
    feats[1091] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '95') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '95', node2.tup)
            
    # C->-24->T
    feats[1092] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-24') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-24', node2.tup)
            
    # C->12_sg->T
    feats[1093] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '12_sg') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '12_sg', node2.tup)
            
    # C->-42->T
    feats[1094] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-42') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-42', node2.tup)
            
    # C->-51->T
    feats[1095] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-51') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-51', node2.tup)
            
    # C->abl. sg.->T
    feats[1096] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'abl. sg.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'abl. sg.', node2.tup)
            
    # C->-293->T
    feats[1097] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-293') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-293', node2.tup)
            
    # C->6_sp->T
    feats[1098] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '6_sp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '6_sp', node2.tup)
            
    # C->-63->T
    feats[1099] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-63') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-63', node2.tup)
            
    # C->acc. sg.->T
    feats[1100] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'acc. sg.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'acc. sg.', node2.tup)
            
    # C->-89->T
    feats[1101] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-89') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-89', node2.tup)
            
    # C->sg_sp->T
    feats[1102] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'sg_sp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'sg_sp', node2.tup)
            
    # C->nom. neutr.->T
    feats[1103] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'nom. neutr.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'nom. neutr.', node2.tup)
            
    # C->42->T
    feats[1104] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '42') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '42', node2.tup)
            
    # C->32->T
    feats[1105] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '32') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '32', node2.tup)
            
    # C->neutr->T
    feats[1106] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'neutr') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'neutr', node2.tup)
            
    # C->15_sp->T
    feats[1107] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '15_sp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '15_sp', node2.tup)
            
    # C->89->T
    feats[1108] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '89') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '89', node2.tup)
            
    # C->9_sg->T
    feats[1109] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '9_sg') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '9_sg', node2.tup)
            
    # C->7_fp->T
    feats[1110] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '7_fp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '7_fp', node2.tup)
            
    # C->139->T
    feats[1111] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '139') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '139', node2.tup)
            
    # C->6_tp->T
    feats[1112] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '6_tp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '6_tp', node2.tup)
            
    # C->-301->T
    feats[1113] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-301') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-301', node2.tup)
            
    # C->73->T
    feats[1114] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '73') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '73', node2.tup)
            
    # C->-247->T
    feats[1115] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-247') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-247', node2.tup)
            
    # C->152->T
    feats[1116] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '152') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '152', node2.tup)
            
    # C->7_pl->T
    feats[1117] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '7_pl') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '7_pl', node2.tup)
            
    # C->97->T
    feats[1118] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '97') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '97', node2.tup)
            
    # C->50->T
    feats[1119] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '50') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '50', node2.tup)
            
    # C->instr. masc.->T
    feats[1120] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'instr. masc.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'instr. masc.', node2.tup)
            
    # C->119->T
    feats[1121] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '119') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '119', node2.tup)
            
    # C->-153->T
    feats[1122] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-153') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-153', node2.tup)
            
    # C->voc. neutr.->T
    feats[1123] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'voc. neutr.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'voc. neutr.', node2.tup)
            
    # C->-67->T
    feats[1124] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-67') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-67', node2.tup)
            
    # C->-62->T
    feats[1125] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-62') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-62', node2.tup)
            
    # C->-76->T
    feats[1126] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-76') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-76', node2.tup)
            
    # C->51->T
    feats[1127] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '51') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '51', node2.tup)
            
    # C->36->T
    feats[1128] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '36') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '36', node2.tup)
            
    # C->11_du->T
    feats[1129] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '11_du') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '11_du', node2.tup)
            
    # C->8_fp->T
    feats[1130] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '8_fp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '8_fp', node2.tup)
            
    # C->80->T
    feats[1131] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '80') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '80', node2.tup)
            
    # C->16_fp->T
    feats[1132] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '16_fp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '16_fp', node2.tup)
            
    # C->-133->T
    feats[1133] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-133') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-133', node2.tup)
            
    # C->173->T
    feats[1134] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '173') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '173', node2.tup)
            
    # C->169->T
    feats[1135] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '169') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '169', node2.tup)
            
    # C->-61->T
    feats[1136] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-61') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-61', node2.tup)
            
    # C->58->T
    feats[1137] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '58') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '58', node2.tup)
            
    # C->15_sg->T
    feats[1138] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '15_sg') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '15_sg', node2.tup)
            
    # C->153->T
    feats[1139] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '153') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '153', node2.tup)
            
    # C->3_sg->T
    feats[1140] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '3_sg') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '3_sg', node2.tup)
            
    # C->gen. sg.->T
    feats[1141] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'gen. sg.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'gen. sg.', node2.tup)
            
    # C->-262->T
    feats[1142] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-262') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-262', node2.tup)
            
    # C->sp->T
    feats[1143] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'sp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'sp', node2.tup)
            
    # C->instr. fem->T
    feats[1144] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'instr. fem') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'instr. fem', node2.tup)
            
    # C->pl_tp->T
    feats[1145] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'pl_tp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'pl_tp', node2.tup)
            
    # C->instr->T
    feats[1146] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'instr') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'instr', node2.tup)
            
    # C->8_sp->T
    feats[1147] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '8_sp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '8_sp', node2.tup)
            
    # C->60->T
    feats[1148] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '60') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '60', node2.tup)
            
    # C->-119->T
    feats[1149] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-119') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-119', node2.tup)
            
    # C->14_pl->T
    feats[1150] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '14_pl') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '14_pl', node2.tup)
            
    # C->voc->T
    feats[1151] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'voc') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'voc', node2.tup)
            
    # C->68->T
    feats[1152] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '68') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '68', node2.tup)
            
    # C->14_sg->T
    feats[1153] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '14_sg') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '14_sg', node2.tup)
            
    # C->15_pl->T
    feats[1154] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '15_pl') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '15_pl', node2.tup)
            
    # C->-83->T
    feats[1155] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-83') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-83', node2.tup)
            
    # C->-249->T
    feats[1156] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-249') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-249', node2.tup)
            
    # C->93->T
    feats[1157] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '93') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '93', node2.tup)
            
    # C->6_du->T
    feats[1158] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '6_du') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '6_du', node2.tup)
            
    # C->71->T
    feats[1159] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '71') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '71', node2.tup)
            
    # C->111->T
    feats[1160] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '111') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '111', node2.tup)
            
    # C->15_du->T
    feats[1161] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '15_du') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '15_du', node2.tup)
            
    # C->10_du->T
    feats[1162] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '10_du') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '10_du', node2.tup)
            
    # C->dat. sg.->T
    feats[1163] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'dat. sg.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'dat. sg.', node2.tup)
            
    # C->15_tp->T
    feats[1164] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '15_tp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '15_tp', node2.tup)
            
    # C->-32->T
    feats[1165] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-32') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-32', node2.tup)
            
    # C->instr. neutr.->T
    feats[1166] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'instr. neutr.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'instr. neutr.', node2.tup)
            
    # C->-54->T
    feats[1167] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-54') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-54', node2.tup)
            
    # C->175->T
    feats[1168] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '175') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '175', node2.tup)
            
    # C->134->T
    feats[1169] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '134') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '134', node2.tup)
            
    # C->-43->T
    feats[1170] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-43') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-43', node2.tup)
            
    # C->-38->T
    feats[1171] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-38') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-38', node2.tup)
            
    # C->-190->T
    feats[1172] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-190') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-190', node2.tup)
            
    # C->157->T
    feats[1173] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '157') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '157', node2.tup)
            
    # C->29_tp->T
    feats[1174] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '29_tp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '29_tp', node2.tup)
            
    # C->30_fp->T
    feats[1175] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '30_fp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '30_fp', node2.tup)
            
    # C->-299->T
    feats[1176] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-299') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-299', node2.tup)
            
    # C->2->T
    feats[1177] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '2') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '2', node2.tup)
            
    # C->pl_sp->T
    feats[1178] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'pl_sp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'pl_sp', node2.tup)
            
    # C->27_pl->T
    feats[1179] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '27_pl') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '27_pl', node2.tup)
            
    # C->8_sg->T
    feats[1180] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '8_sg') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '8_sg', node2.tup)
            
    # C->69->T
    feats[1181] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '69') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '69', node2.tup)
            
    # C->155->T
    feats[1182] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '155') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '155', node2.tup)
            
    # C->154->T
    feats[1183] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '154') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '154', node2.tup)
            
    # C->-57->T
    feats[1184] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-57') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-57', node2.tup)
            
    # C->108->T
    feats[1185] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '108') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '108', node2.tup)
            
    # C->-87->T
    feats[1186] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-87') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-87', node2.tup)
            
    # C->2_tp->T
    feats[1187] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '2_tp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '2_tp', node2.tup)
            
    # C->voc. pl.->T
    feats[1188] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, 'voc. pl.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'voc. pl.', node2.tup)
            
    # C->-35->T
    feats[1189] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-35') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-35', node2.tup)
            
    # C->91->T
    feats[1190] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '91') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '91', node2.tup)
            
    # C->4_tp->T
    feats[1191] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '4_tp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '4_tp', node2.tup)
            
    # C->-308->T
    feats[1192] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-308') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-308', node2.tup)
            
    # C->130->T
    feats[1193] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '130') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '130', node2.tup)
            
    # C->-68->T
    feats[1194] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-68') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-68', node2.tup)
            
    # C->-240->T
    feats[1195] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-240') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-240', node2.tup)
            
    # C->160->T
    feats[1196] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '160') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '160', node2.tup)
            
    # C->-90->T
    feats[1197] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-90') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-90', node2.tup)
            
    # C->161->T
    feats[1198] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '161') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '161', node2.tup)
            
    # C->162->T
    feats[1199] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '162') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '162', node2.tup)
            
    # C->172->T
    feats[1200] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '172') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '172', node2.tup)
            
    # C->-30->T
    feats[1201] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-30') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-30', node2.tup)
            
    # C->-16->T
    feats[1202] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-16') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-16', node2.tup)
            
    # C->10_pl->T
    feats[1203] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '10_pl') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '10_pl', node2.tup)
            
    # C->27_du->T
    feats[1204] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '27_du') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '27_du', node2.tup)
            
    # C->-72->T
    feats[1205] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-72') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-72', node2.tup)
            
    # C->-141->T
    feats[1206] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-141') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-141', node2.tup)
            
    # C->15_fp->T
    feats[1207] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '15_fp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '15_fp', node2.tup)
            
    # C->13_fp->T
    feats[1208] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '13_fp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '13_fp', node2.tup)
            
    # C->-64->T
    feats[1209] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-64') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-64', node2.tup)
            
    # C->176->T
    feats[1210] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '176') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '176', node2.tup)
            
    # C->-129->T
    feats[1211] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-129') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-129', node2.tup)
            
    # C->100->T
    feats[1212] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '100') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '100', node2.tup)
            
    # C->8_pl->T
    feats[1213] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '8_pl') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '8_pl', node2.tup)
            
    # C->-28->T
    feats[1214] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-28') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-28', node2.tup)
            
    # C->-303->T
    feats[1215] = tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, node1.cng, '-303') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-303', node2.tup)
            
    # T->abl. du.->L
    feats[1216] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'abl. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. du.', node2.lemma)
            
    # T->39->L
    feats[1217] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '39') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '39', node2.lemma)
            
    # T->-94->L
    feats[1218] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-94') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-94', node2.lemma)
            
    # T->13_sg->L
    feats[1219] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '13_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_sg', node2.lemma)
            
    # T->72->L
    feats[1220] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '72') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '72', node2.lemma)
            
    # T->12_fp->L
    feats[1221] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '12_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '12_fp', node2.lemma)
            
    # T->137->L
    feats[1222] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '137') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '137', node2.lemma)
            
    # T->-169->L
    feats[1223] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-169') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-169', node2.lemma)
            
    # T->-11->L
    feats[1224] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-11') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-11', node2.lemma)
            
    # T->loc. sg.->L
    feats[1225] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'loc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'loc. sg.', node2.lemma)
            
    # T->5_sp->L
    feats[1226] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '5_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_sp', node2.lemma)
            
    # T->34->L
    feats[1227] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '34') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '34', node2.lemma)
            
    # T->179->L
    feats[1228] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '179') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '179', node2.lemma)
            
    # T->-152->L
    feats[1229] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-152') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-152', node2.lemma)
            
    # T->-149->L
    feats[1230] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-149') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-149', node2.lemma)
            
    # T->instr. du.->L
    feats[1231] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'instr. du.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. du.', node2.lemma)
            
    # T->151->L
    feats[1232] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '151') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '151', node2.lemma)
            
    # T->2_fp->L
    feats[1233] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '2_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_fp', node2.lemma)
            
    # T->5_pl->L
    feats[1234] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '5_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '5_pl', node2.lemma)
            
    # T->-230->L
    feats[1235] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-230') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-230', node2.lemma)
            
    # T->-132->L
    feats[1236] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-132') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-132', node2.lemma)
            
    # T->120->L
    feats[1237] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '120') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '120', node2.lemma)
            
    # T->182->L
    feats[1238] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '182') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '182', node2.lemma)
            
    # T->28_tp->L
    feats[1239] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '28_tp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '28_tp', node2.lemma)
            
    # T->115->L
    feats[1240] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '115') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '115', node2.lemma)
            
    # T->156->L
    feats[1241] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '156') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '156', node2.lemma)
            
    # T->118->L
    feats[1242] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '118') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '118', node2.lemma)
            
    # T->121->L
    feats[1243] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '121') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '121', node2.lemma)
            
    # T->-302->L
    feats[1244] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-302') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-302', node2.lemma)
            
    # T->voc. sg.->L
    feats[1245] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'voc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'voc. sg.', node2.lemma)
            
    # T->fp->L
    feats[1246] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'fp', node2.lemma)
            
    # T->-104->L
    feats[1247] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-104') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-104', node2.lemma)
            
    # T->2_du->L
    feats[1248] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '2_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '2_du', node2.lemma)
            
    # T->pl->L
    feats[1249] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'pl', node2.lemma)
            
    # T->-292->L
    feats[1250] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-292') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-292', node2.lemma)
            
    # T->114->L
    feats[1251] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '114') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '114', node2.lemma)
            
    # T->-36->L
    feats[1252] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-36') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-36', node2.lemma)
            
    # T->-156->L
    feats[1253] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-156') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-156', node2.lemma)
            
    # T->70->L
    feats[1254] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '70') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '70', node2.lemma)
            
    # T->-27->L
    feats[1255] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-27') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-27', node2.lemma)
            
    # T->-144->L
    feats[1256] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-144') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-144', node2.lemma)
            
    # T->-24->L
    feats[1257] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-24') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-24', node2.lemma)
            
    # T->-51->L
    feats[1258] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-51') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-51', node2.lemma)
            
    # T->abl. sg.->L
    feats[1259] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'abl. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'abl. sg.', node2.lemma)
            
    # T->-293->L
    feats[1260] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-293') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-293', node2.lemma)
            
    # T->6_sp->L
    feats[1261] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '6_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '6_sp', node2.lemma)
            
    # T->-63->L
    feats[1262] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-63') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-63', node2.lemma)
            
    # T->acc. sg.->L
    feats[1263] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'acc. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'acc. sg.', node2.lemma)
            
    # T->-89->L
    feats[1264] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-89') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-89', node2.lemma)
            
    # T->sg_sp->L
    feats[1265] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'sg_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'sg_sp', node2.lemma)
            
    # T->42->L
    feats[1266] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '42') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '42', node2.lemma)
            
    # T->7_fp->L
    feats[1267] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '7_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '7_fp', node2.lemma)
            
    # T->-301->L
    feats[1268] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-301') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-301', node2.lemma)
            
    # T->-247->L
    feats[1269] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-247') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-247', node2.lemma)
            
    # T->119->L
    feats[1270] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '119') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '119', node2.lemma)
            
    # T->-67->L
    feats[1271] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-67') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-67', node2.lemma)
            
    # T->-76->L
    feats[1272] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-76') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-76', node2.lemma)
            
    # T->51->L
    feats[1273] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '51') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '51', node2.lemma)
            
    # T->8_fp->L
    feats[1274] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '8_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_fp', node2.lemma)
            
    # T->80->L
    feats[1275] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '80') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '80', node2.lemma)
            
    # T->16_fp->L
    feats[1276] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '16_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '16_fp', node2.lemma)
            
    # T->-61->L
    feats[1277] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-61') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-61', node2.lemma)
            
    # T->58->L
    feats[1278] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '58') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '58', node2.lemma)
            
    # T->3_sg->L
    feats[1279] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '3_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '3_sg', node2.lemma)
            
    # T->gen. sg.->L
    feats[1280] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'gen. sg.') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'gen. sg.', node2.lemma)
            
    # T->instr. fem->L
    feats[1281] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'instr. fem') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, 'instr. fem', node2.lemma)
            
    # T->8_sp->L
    feats[1282] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '8_sp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_sp', node2.lemma)
            
    # T->60->L
    feats[1283] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '60') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '60', node2.lemma)
            
    # T->14_pl->L
    feats[1284] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '14_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '14_pl', node2.lemma)
            
    # T->15_pl->L
    feats[1285] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '15_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '15_pl', node2.lemma)
            
    # T->-83->L
    feats[1286] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-83') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-83', node2.lemma)
            
    # T->93->L
    feats[1287] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '93') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '93', node2.lemma)
            
    # T->10_du->L
    feats[1288] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '10_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_du', node2.lemma)
            
    # T->-32->L
    feats[1289] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-32') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-32', node2.lemma)
            
    # T->134->L
    feats[1290] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '134') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '134', node2.lemma)
            
    # T->-43->L
    feats[1291] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-43') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-43', node2.lemma)
            
    # T->8_sg->L
    feats[1292] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '8_sg') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '8_sg', node2.lemma)
            
    # T->69->L
    feats[1293] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '69') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '69', node2.lemma)
            
    # T->155->L
    feats[1294] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '155') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '155', node2.lemma)
            
    # T->154->L
    feats[1295] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '154') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '154', node2.lemma)
            
    # T->108->L
    feats[1296] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '108') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '108', node2.lemma)
            
    # T->-35->L
    feats[1297] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-35') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-35', node2.lemma)
            
    # T->-240->L
    feats[1298] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-240') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-240', node2.lemma)
            
    # T->-90->L
    feats[1299] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-90') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-90', node2.lemma)
            
    # T->10_pl->L
    feats[1300] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '10_pl') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '10_pl', node2.lemma)
            
    # T->27_du->L
    feats[1301] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '27_du') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '27_du', node2.lemma)
            
    # T->13_fp->L
    feats[1302] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '13_fp') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '13_fp', node2.lemma)
            
    # T->-64->L
    feats[1303] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-64') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-64', node2.lemma)
            
    # T->176->L
    feats[1304] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '176') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '176', node2.lemma)
            
    # T->-28->L
    feats[1305] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-28') * tryProb_catchZero(mat_cng2lem_countonly, mat_cngCount_1D, '-28', node2.lemma)
            
    # T->abl. du.->C
    feats[1306] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'abl. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'abl. du.', node2.cng)
            
    # T->-94->C
    feats[1307] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-94') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-94', node2.cng)
            
    # T->13_sg->C
    feats[1308] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '13_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_sg', node2.cng)
            
    # T->72->C
    feats[1309] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '72') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '72', node2.cng)
            
    # T->12_fp->C
    feats[1310] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '12_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '12_fp', node2.cng)
            
    # T->137->C
    feats[1311] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '137') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '137', node2.cng)
            
    # T->-169->C
    feats[1312] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-169') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-169', node2.cng)
            
    # T->5_sp->C
    feats[1313] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '5_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_sp', node2.cng)
            
    # T->34->C
    feats[1314] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '34') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '34', node2.cng)
            
    # T->179->C
    feats[1315] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '179') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '179', node2.cng)
            
    # T->-152->C
    feats[1316] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-152') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-152', node2.cng)
            
    # T->-149->C
    feats[1317] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-149') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-149', node2.cng)
            
    # T->10_tp->C
    feats[1318] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '10_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_tp', node2.cng)
            
    # T->-245->C
    feats[1319] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-245') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-245', node2.cng)
            
    # T->instr. du.->C
    feats[1320] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'instr. du.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'instr. du.', node2.cng)
            
    # T->151->C
    feats[1321] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '151') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '151', node2.cng)
            
    # T->2_fp->C
    feats[1322] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '2_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_fp', node2.cng)
            
    # T->5_pl->C
    feats[1323] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '5_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '5_pl', node2.cng)
            
    # T->-230->C
    feats[1324] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-230') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-230', node2.cng)
            
    # T->-132->C
    feats[1325] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-132') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-132', node2.cng)
            
    # T->120->C
    feats[1326] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '120') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '120', node2.cng)
            
    # T->182->C
    feats[1327] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '182') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '182', node2.cng)
            
    # T->28_tp->C
    feats[1328] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '28_tp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '28_tp', node2.cng)
            
    # T->115->C
    feats[1329] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '115') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '115', node2.cng)
            
    # T->156->C
    feats[1330] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '156', node2.cng)
            
    # T->118->C
    feats[1331] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '118') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '118', node2.cng)
            
    # T->121->C
    feats[1332] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '121') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '121', node2.cng)
            
    # T->-302->C
    feats[1333] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-302') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-302', node2.cng)
            
    # T->fp->C
    feats[1334] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'fp', node2.cng)
            
    # T->-104->C
    feats[1335] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-104') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-104', node2.cng)
            
    # T->2_du->C
    feats[1336] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '2_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '2_du', node2.cng)
            
    # T->pl->C
    feats[1337] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'pl', node2.cng)
            
    # T->-292->C
    feats[1338] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-292') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-292', node2.cng)
            
    # T->114->C
    feats[1339] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '114') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '114', node2.cng)
            
    # T->-36->C
    feats[1340] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-36') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-36', node2.cng)
            
    # T->-156->C
    feats[1341] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-156') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-156', node2.cng)
            
    # T->70->C
    feats[1342] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '70') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '70', node2.cng)
            
    # T->-27->C
    feats[1343] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-27') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-27', node2.cng)
            
    # T->-144->C
    feats[1344] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-144') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-144', node2.cng)
            
    # T->-307->C
    feats[1345] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-307') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-307', node2.cng)
            
    # T->81->C
    feats[1346] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '81') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '81', node2.cng)
            
    # T->-24->C
    feats[1347] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-24') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-24', node2.cng)
            
    # T->abl. sg.->C
    feats[1348] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'abl. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'abl. sg.', node2.cng)
            
    # T->-293->C
    feats[1349] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-293') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-293', node2.cng)
            
    # T->6_sp->C
    feats[1350] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '6_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '6_sp', node2.cng)
            
    # T->-63->C
    feats[1351] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-63') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-63', node2.cng)
            
    # T->acc. sg.->C
    feats[1352] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'acc. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'acc. sg.', node2.cng)
            
    # T->sg_sp->C
    feats[1353] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'sg_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'sg_sp', node2.cng)
            
    # T->42->C
    feats[1354] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '42') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '42', node2.cng)
            
    # T->7_fp->C
    feats[1355] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '7_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '7_fp', node2.cng)
            
    # T->-301->C
    feats[1356] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-301') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-301', node2.cng)
            
    # T->119->C
    feats[1357] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '119') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '119', node2.cng)
            
    # T->-67->C
    feats[1358] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-67') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-67', node2.cng)
            
    # T->-76->C
    feats[1359] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-76') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-76', node2.cng)
            
    # T->51->C
    feats[1360] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '51') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '51', node2.cng)
            
    # T->8_fp->C
    feats[1361] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '8_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_fp', node2.cng)
            
    # T->80->C
    feats[1362] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '80') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '80', node2.cng)
            
    # T->16_fp->C
    feats[1363] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '16_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '16_fp', node2.cng)
            
    # T->-61->C
    feats[1364] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-61') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-61', node2.cng)
            
    # T->58->C
    feats[1365] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '58') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '58', node2.cng)
            
    # T->3_sg->C
    feats[1366] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '3_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '3_sg', node2.cng)
            
    # T->gen. sg.->C
    feats[1367] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'gen. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'gen. sg.', node2.cng)
            
    # T->-262->C
    feats[1368] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-262') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-262', node2.cng)
            
    # T->8_sp->C
    feats[1369] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '8_sp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_sp', node2.cng)
            
    # T->60->C
    feats[1370] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '60') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '60', node2.cng)
            
    # T->14_pl->C
    feats[1371] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '14_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '14_pl', node2.cng)
            
    # T->15_pl->C
    feats[1372] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '15_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '15_pl', node2.cng)
            
    # T->-83->C
    feats[1373] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-83') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-83', node2.cng)
            
    # T->93->C
    feats[1374] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '93') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '93', node2.cng)
            
    # T->10_du->C
    feats[1375] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '10_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_du', node2.cng)
            
    # T->dat. sg.->C
    feats[1376] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'dat. sg.') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, 'dat. sg.', node2.cng)
            
    # T->-32->C
    feats[1377] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-32') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-32', node2.cng)
            
    # T->134->C
    feats[1378] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '134') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '134', node2.cng)
            
    # T->-43->C
    feats[1379] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-43') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-43', node2.cng)
            
    # T->-299->C
    feats[1380] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-299') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-299', node2.cng)
            
    # T->8_sg->C
    feats[1381] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '8_sg') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '8_sg', node2.cng)
            
    # T->69->C
    feats[1382] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '69') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '69', node2.cng)
            
    # T->155->C
    feats[1383] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '155') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '155', node2.cng)
            
    # T->108->C
    feats[1384] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '108') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '108', node2.cng)
            
    # T->-240->C
    feats[1385] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-240') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-240', node2.cng)
            
    # T->-90->C
    feats[1386] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-90') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-90', node2.cng)
            
    # T->10_pl->C
    feats[1387] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '10_pl') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '10_pl', node2.cng)
            
    # T->27_du->C
    feats[1388] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '27_du') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '27_du', node2.cng)
            
    # T->13_fp->C
    feats[1389] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '13_fp') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '13_fp', node2.cng)
            
    # T->-64->C
    feats[1390] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-64') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-64', node2.cng)
            
    # T->176->C
    feats[1391] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '176') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '176', node2.cng)
            
    # T->-28->C
    feats[1392] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-28') * tryProb_catchZero(mat_cng2cng_countonly, mat_cngCount_1D, '-28', node2.cng)
            
    # T->abl. du.->T
    feats[1393] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'abl. du.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'abl. du.', node2.tup)
            
    # T->39->T
    feats[1394] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '39') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '39', node2.tup)
            
    # T->-94->T
    feats[1395] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-94') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-94', node2.tup)
            
    # T->13_sg->T
    feats[1396] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '13_sg') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '13_sg', node2.tup)
            
    # T->72->T
    feats[1397] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '72') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '72', node2.tup)
            
    # T->12_fp->T
    feats[1398] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '12_fp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '12_fp', node2.tup)
            
    # T->137->T
    feats[1399] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '137') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '137', node2.tup)
            
    # T->-169->T
    feats[1400] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-169') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-169', node2.tup)
            
    # T->-11->T
    feats[1401] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-11') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-11', node2.tup)
            
    # T->loc. sg.->T
    feats[1402] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'loc. sg.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'loc. sg.', node2.tup)
            
    # T->5_sp->T
    feats[1403] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '5_sp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '5_sp', node2.tup)
            
    # T->34->T
    feats[1404] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '34') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '34', node2.tup)
            
    # T->179->T
    feats[1405] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '179') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '179', node2.tup)
            
    # T->-152->T
    feats[1406] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-152') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-152', node2.tup)
            
    # T->-149->T
    feats[1407] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-149') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-149', node2.tup)
            
    # T->-245->T
    feats[1408] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-245') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-245', node2.tup)
            
    # T->instr. du.->T
    feats[1409] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'instr. du.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'instr. du.', node2.tup)
            
    # T->151->T
    feats[1410] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '151') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '151', node2.tup)
            
    # T->2_fp->T
    feats[1411] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '2_fp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '2_fp', node2.tup)
            
    # T->5_pl->T
    feats[1412] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '5_pl') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '5_pl', node2.tup)
            
    # T->-230->T
    feats[1413] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-230') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-230', node2.tup)
            
    # T->-103->T
    feats[1414] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-103') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-103', node2.tup)
            
    # T->-132->T
    feats[1415] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-132') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-132', node2.tup)
            
    # T->120->T
    feats[1416] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '120') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '120', node2.tup)
            
    # T->182->T
    feats[1417] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '182') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '182', node2.tup)
            
    # T->28_tp->T
    feats[1418] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '28_tp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '28_tp', node2.tup)
            
    # T->115->T
    feats[1419] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '115') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '115', node2.tup)
            
    # T->156->T
    feats[1420] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '156') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '156', node2.tup)
            
    # T->118->T
    feats[1421] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '118') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '118', node2.tup)
            
    # T->acc. fem->T
    feats[1422] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'acc. fem') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'acc. fem', node2.tup)
            
    # T->121->T
    feats[1423] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '121') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '121', node2.tup)
            
    # T->-302->T
    feats[1424] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-302') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-302', node2.tup)
            
    # T->voc. sg.->T
    feats[1425] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'voc. sg.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'voc. sg.', node2.tup)
            
    # T->fp->T
    feats[1426] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'fp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'fp', node2.tup)
            
    # T->-23->T
    feats[1427] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-23') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-23', node2.tup)
            
    # T->-104->T
    feats[1428] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-104') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-104', node2.tup)
            
    # T->2_du->T
    feats[1429] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '2_du') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '2_du', node2.tup)
            
    # T->92->T
    feats[1430] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '92') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '92', node2.tup)
            
    # T->pl->T
    feats[1431] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'pl') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'pl', node2.tup)
            
    # T->-292->T
    feats[1432] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-292') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-292', node2.tup)
            
    # T->114->T
    feats[1433] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '114') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '114', node2.tup)
            
    # T->-36->T
    feats[1434] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-36') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-36', node2.tup)
            
    # T->-156->T
    feats[1435] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-156') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-156', node2.tup)
            
    # T->70->T
    feats[1436] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '70') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '70', node2.tup)
            
    # T->-27->T
    feats[1437] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-27') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-27', node2.tup)
            
    # T->-144->T
    feats[1438] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-144') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-144', node2.tup)
            
    # T->117->T
    feats[1439] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '117') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '117', node2.tup)
            
    # T->-24->T
    feats[1440] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-24') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-24', node2.tup)
            
    # T->12_sg->T
    feats[1441] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '12_sg') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '12_sg', node2.tup)
            
    # T->-51->T
    feats[1442] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-51') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-51', node2.tup)
            
    # T->abl. sg.->T
    feats[1443] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'abl. sg.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'abl. sg.', node2.tup)
            
    # T->-293->T
    feats[1444] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-293') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-293', node2.tup)
            
    # T->6_sp->T
    feats[1445] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '6_sp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '6_sp', node2.tup)
            
    # T->-63->T
    feats[1446] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-63') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-63', node2.tup)
            
    # T->acc. sg.->T
    feats[1447] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'acc. sg.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'acc. sg.', node2.tup)
            
    # T->-89->T
    feats[1448] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-89') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-89', node2.tup)
            
    # T->sg_sp->T
    feats[1449] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'sg_sp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'sg_sp', node2.tup)
            
    # T->42->T
    feats[1450] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '42') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '42', node2.tup)
            
    # T->7_fp->T
    feats[1451] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '7_fp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '7_fp', node2.tup)
            
    # T->139->T
    feats[1452] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '139') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '139', node2.tup)
            
    # T->-301->T
    feats[1453] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-301') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-301', node2.tup)
            
    # T->-247->T
    feats[1454] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-247') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-247', node2.tup)
            
    # T->119->T
    feats[1455] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '119') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '119', node2.tup)
            
    # T->-153->T
    feats[1456] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-153') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-153', node2.tup)
            
    # T->-67->T
    feats[1457] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-67') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-67', node2.tup)
            
    # T->-76->T
    feats[1458] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-76') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-76', node2.tup)
            
    # T->51->T
    feats[1459] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '51') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '51', node2.tup)
            
    # T->8_fp->T
    feats[1460] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '8_fp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '8_fp', node2.tup)
            
    # T->80->T
    feats[1461] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '80') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '80', node2.tup)
            
    # T->16_fp->T
    feats[1462] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '16_fp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '16_fp', node2.tup)
            
    # T->-133->T
    feats[1463] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-133') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-133', node2.tup)
            
    # T->-61->T
    feats[1464] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-61') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-61', node2.tup)
            
    # T->58->T
    feats[1465] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '58') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '58', node2.tup)
            
    # T->3_sg->T
    feats[1466] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '3_sg') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '3_sg', node2.tup)
            
    # T->gen. sg.->T
    feats[1467] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'gen. sg.') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'gen. sg.', node2.tup)
            
    # T->-262->T
    feats[1468] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-262') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-262', node2.tup)
            
    # T->instr. fem->T
    feats[1469] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, 'instr. fem') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, 'instr. fem', node2.tup)
            
    # T->8_sp->T
    feats[1470] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '8_sp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '8_sp', node2.tup)
            
    # T->60->T
    feats[1471] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '60') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '60', node2.tup)
            
    # T->14_pl->T
    feats[1472] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '14_pl') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '14_pl', node2.tup)
            
    # T->15_pl->T
    feats[1473] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '15_pl') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '15_pl', node2.tup)
            
    # T->-83->T
    feats[1474] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-83') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-83', node2.tup)
            
    # T->93->T
    feats[1475] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '93') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '93', node2.tup)
            
    # T->15_du->T
    feats[1476] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '15_du') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '15_du', node2.tup)
            
    # T->10_du->T
    feats[1477] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '10_du') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '10_du', node2.tup)
            
    # T->-32->T
    feats[1478] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-32') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-32', node2.tup)
            
    # T->-54->T
    feats[1479] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-54') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-54', node2.tup)
            
    # T->134->T
    feats[1480] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '134') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '134', node2.tup)
            
    # T->-43->T
    feats[1481] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-43') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-43', node2.tup)
            
    # T->27_pl->T
    feats[1482] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '27_pl') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '27_pl', node2.tup)
            
    # T->8_sg->T
    feats[1483] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '8_sg') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '8_sg', node2.tup)
            
    # T->69->T
    feats[1484] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '69') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '69', node2.tup)
            
    # T->155->T
    feats[1485] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '155') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '155', node2.tup)
            
    # T->154->T
    feats[1486] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '154') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '154', node2.tup)
            
    # T->108->T
    feats[1487] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '108') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '108', node2.tup)
            
    # T->2_tp->T
    feats[1488] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '2_tp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '2_tp', node2.tup)
            
    # T->-35->T
    feats[1489] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-35') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-35', node2.tup)
            
    # T->4_tp->T
    feats[1490] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '4_tp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '4_tp', node2.tup)
            
    # T->-240->T
    feats[1491] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-240') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-240', node2.tup)
            
    # T->-90->T
    feats[1492] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-90') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-90', node2.tup)
            
    # T->10_pl->T
    feats[1493] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '10_pl') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '10_pl', node2.tup)
            
    # T->27_du->T
    feats[1494] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '27_du') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '27_du', node2.tup)
            
    # T->-141->T
    feats[1495] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-141') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-141', node2.tup)
            
    # T->13_fp->T
    feats[1496] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '13_fp') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '13_fp', node2.tup)
            
    # T->-64->T
    feats[1497] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-64') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-64', node2.tup)
            
    # T->176->T
    feats[1498] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '176') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '176', node2.tup)
            
    # T->-28->T
    feats[1499] = tryProb_catchZero(mat_tup2cng_countonly, mat_tupCount_1D, node1.tup, '-28') * tryProb_catchZero(mat_cng2tup_countonly, mat_cngCount_1D, '-28', node2.tup)
            
    feats[feats < 1e-25] = 1e-25
    return -np.log10(feats)
    