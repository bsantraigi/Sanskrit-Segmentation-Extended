import numpy as np

"""
################################################################################################
###################  METHODs: SIGMOID and DERIVATIVE OF SIGMOID ################################
################################################################################################
"""

def sigmoid(vec):
    evec = 1 + np.exp(-vec)
    return 1/evec
    
def d_sigmoid(output_of_gate):
    return output_of_gate*(1-output_of_gate)
    
"""
################################################################################################
###################  METHODs: ReLU AND DERIVATE OF ReLU ########################################
################################################################################################
"""

def relu(vec_x):
    relu_x = vec_x.copy()
    relu_x[vec_x < 0] = 0
    return relu_x
def d_relu(vec_x):
    d_relu_x = vec_x.copy()
    d_relu_x[vec_x > 0] = 1
    d_relu_x[vec_x <= 0] = 0
    return d_relu_x
    
"""
################################################################################################
##################  IMPLEMENTATION OF NEURAL NETWORK  ##########################################
################################################################################################
"""

class NN:
    def __init__(self, input_dimension, hidden_layer_size):
        # d: Input feature dimension i.e. the dimension of the edge feature vectors
        # n: Hidden layer size
        
        # TODO: Add Bias terms
        
        self.n = hidden_layer_size
        self.d = input_dimension
        # rand_init_range = np.sqrt(6.0/(self.n+self.d))
        rand_init_range = 10
        self.W = np.random.uniform(-rand_init_range, rand_init_range, (self.n, self.d))
        # rand_init_range = np.sqrt(6.0/(self.n))
        self.U = np.random.uniform(-rand_init_range, rand_init_range, (self.n, 1))

    def Forward_Prop(self, x):
        z2 = np.matmul(self.W, x)
        a2 = relu(z2)
        o = np.matmul(self.U.transpose(), a2)
        s = sigmoid(o)
        return (z2, a2, s)
    
    def Back_Prop(self, dLdS, nodeLen, featVMat):
        N = nodeLen
        dLdU = np.zeros(self.U.shape)
        dLdW = np.zeros(self.W.shape)
        eta = 100
        _actual_sample_count = 0
        for i in range(N):
            for j in range(N):
                if dLdS[i, j] != 0:
                    _actual_sample_count += 1
                    (z2, a2, s) = self.Forward_Prop(featVMat[i][j])                
                    dLdU += dLdS[i, j]*a2*d_sigmoid(s)
                    dRelu = d_relu(z2)
                    for k in range(self.n):
                        dLdW[k, :] = dLdW[k, :] + ((dLdS[i, j])*self.U[k]*dRelu[k]*d_sigmoid(s))*featVMat[i][j].transpose()
        
        self.W -= eta*dLdW
        self.U -= eta*dLdU
        
        
        
        
        
        
        
        
        
        
        
