from sentences import *
from DCS import *
from romtoslp import *
from collections import defaultdict
import numpy as np

class word_definite:
    def __init__(self,derived, lemma, cng, pos, chunk_id):
        self.lemma = lemma
        self.derived = derived
        self.cng = cng
        # self.form = form
        self.pos = pos
        self.chunk_id = chunk_id        
        # Fields Required for heap
        self.dist = np.inf
        self.src = -1
        self.id = -1
        self.isConflicted = False
    def __str__(self):
        return 'WD_Node[C: %d, P: %d, %s @(%d) => %s]' %(self.chunk_id, self.pos, self.lemma, self.cng, self.derived)
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
################################################################################################
"""
mat_lem2cng_countonly = pickle.load(open('../NewData/ultimate_new_lastsem/mat_lem2cng_countonly.p', 'rb'))
mat_cng2lem_countonly = pickle.load(open('../NewData/ultimate_new_lastsem/mat_cng2lem_countonly.p', 'rb'))
mat_lem2cng_1D = pickle.load(open('../NewData/ultimate_new_lastsem/mat_lem2cng_1D.p', 'rb'))
mat_cng2lem_1D = pickle.load(open('../NewData/ultimate_new_lastsem/mat_cng2lem_1D.p', 'rb'))

_edge_vector_dim = len(mat_cng2lem_1D)
_full_cnglist = list(mat_cng2lem_1D)
 

"""
################################################################################################
######################  CREATING FEATURE MATRICES ##############################################
################################################################################################
"""
def Get_Feat_Vec_Matrix(nodelist_new, conflicts_Dict):
    nodesCount = len(nodelist_new)
    featVMat = [[None for _ in range(nodesCount)] for _ in range(nodesCount)]
    for i in range(nodesCount):
        for j in range(nodesCount):
            if j in conflicts_Dict[i] or i == j:                
                # featVMat[i][j][:] = 1e-35
                pass
            else:
                wd1 = nodelist_new[i]
                wd2 = nodelist_new[j]
                featVMat[i][j] = np.zeros((_edge_vector_dim, 1))
                # print('For ', wd1, wd2)
                for k in range(_edge_vector_dim):
                    cng_k = _full_cnglist[k]
                    # TODO: Some lemma's still missing
                    try:
                        pleft = float(mat_lem2cng_countonly[wd1.lemma][cng_k]) / mat_lem2cng_1D[wd1.lemma]
                    except KeyError:
                        pleft = 0
                    try:
                        pright = float(mat_cng2lem_countonly[cng_k][wd2.lemma]) / mat_cng2lem_1D[cng_k]
                    except KeyError:
                        pright = 0
                        
                    featVMat[i][j][k] = pleft * pright                
    return featVMat

# This function will be some neural network or a linear function or something of sorts
def Get_W_Scalar_Matrix_from_FeatVect_Matrix(featVMat, nodelist_new, conflicts_Dict, _neuralnet):
    nodesCount = len(nodelist_new)
    WScalarMat = np.zeros((nodesCount, nodesCount))   
    sigmoid_25 = 1/(1+np.exp(-25))
    for i in range(nodesCount):
        for j in range(nodesCount):
            if featVMat[i][j] is None:
                pass
            else:
                # Since s is output of a sigmoid gate, it will always be greater than zero
                (_, _, s) = _neuralnet.Forward_Prop(featVMat[i][j])
                WScalarMat[i, j] = np.minimum(s, sigmoid_25)
    toinf = (WScalarMat == 0)
    WScalarMat[WScalarMat > 0] = -np.log2(WScalarMat[WScalarMat > 0])
    WScalarMat[toinf] = np.inf
    return WScalarMat
