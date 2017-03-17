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

def lrelu(vec_x):
    relu_x = vec_x.copy()
    relu_x[vec_x < 0] = relu_x[vec_x < 0]/100
    return relu_x

def d_relu(vec_x):
    d_relu_x = vec_x.copy()
    d_relu_x[vec_x > 0] = 1
    d_relu_x[vec_x <= 0] = 0
    return d_relu_x

def d_lrelu(vec_x):
    d_relu_x = vec_x.copy()
    d_relu_x[vec_x > 0] = 1
    d_relu_x[vec_x <= 0] = 0.01
    return d_relu_x
    
"""
################################################################################################
##################  IMPLEMENTATION OF NEURAL NETWORK  ##########################################
################################################################################################
"""

class NN:
    def __init__(self, input_dimension, hidden_layer_size, outer_relu = True):
        # d: Input feature dimension i.e. the dimension of the edge feature vectors
        # n: Hidden layer size
        
        # TODO: Add Bias terms
        
        self.n = hidden_layer_size
        self.d = input_dimension
        # rand_init_range = np.sqrt(6.0/(self.n+self.d))
        # print(rand_init_range)
        
        rand_init_range = 1e-2
        self.W = np.random.uniform(-rand_init_range, rand_init_range, (self.n, self.d))
        # rand_init_range = np.sqrt(6.0/(self.n))

        rand_init_range = 1e-1
        self.U = np.random.uniform(-rand_init_range, rand_init_range, (self.n, 1))

        # Apply relu or sigmoid at the output layer
        # If relu is applied it will be assumed that log is applied to the 
        #   feature before passing it to the network
        # Else in case of outer sigmoid 
        #   log is applied after the neural network
        self.outer_relu = outer_relu

    def Forward_Prop(self, x):
        z2 = np.matmul(self.W, x)
        a2 = lrelu(z2)
        o = np.matmul(self.U.transpose(), a2)
        if self.outer_relu:
            # s = relu(o)
            s = o
        else:
            s = sigmoid(o)
        return (z2, a2, s)
    
    # Back_Propagate gradient of Loss, L:  Assuming S is the direct output of the network
    def Back_Prop(self, dLdOut, nodeLen, featVMat):
        N = nodeLen
        dLdU = np.zeros(self.U.shape)
        dLdW = np.zeros(self.W.shape)
        if not self.outer_relu:
            etaW = 5e6
            etaU = 5e3
            batch_size = 0
            for i in range(N):
                for j in range(N):
                    if dLdOut[i, j] != 0 and (featVMat[i][j] is not None):
                        batch_size += 1
                        (z2, a2, s) = self.Forward_Prop(featVMat[i][j])
                        dLdU += dLdOut[i, j]*a2
                        dRelu = d_relu(z2)
                        for k in range(self.n):
                            if dRelu[k] != 0:
                                dLdW[k, :, None] += (dLdOut[i, j])*self.U[k]*dRelu[k]*(featVMat[i][j])
            
            self.W -= etaW*dLdW/(batch_size)
            self.U -= etaU*dLdU/(batch_size)
        else:
            etaW = 0.0001
            etaU = 0.00008
            batch_size = 0
            # print(dLdOut)
            for i in range(N):
                for j in range(N):
                    if dLdOut[i, j] != 0 and (featVMat[i][j] is not None):
                        batch_size += 1
                        (z2, a2, s) = self.Forward_Prop(featVMat[i][j])
                        # print(a2.transpose())
                        # print('o')
                        # print(np.matmul(self.U.transpose(), a2))
                        dLdU += dLdOut[i, j]*a2
                        dRelu = d_lrelu(z2)
                        dLdW += (dLdOut[i, j])*(self.U*dRelu)*featVMat[i][j].transpose()
                        # for k in range(self.n):
                        #     if dRelu[k] != 0:
                        #         dLdW[k, :, None] += (dLdOut[i, j])*self.U[k]*dRelu[k]*(featVMat[i][j])
            # print('dlDW:')
            # print(dLdW/(batch_size))
            # print('dlDU:')
            # print(dLdU/(batch_size))
            # print('Batch size: ', batch_size)
            if batch_size > 0:
                # delW = etaW*dLdW/(batch_size)
                # delU = etaU*dLdU/(batch_size)
                delW = etaW*dLdW/(batch_size)
                delU = etaU*dLdU/(batch_size)
                print('Max(delW): %10.6f\tMax(delU): %10.6f'%(np.max(np.abs(delW)), np.max(np.abs(delU))))
                self.W -= delW
                self.U -= delU

        
        
        
        
        
        
        
        
        
        
        
