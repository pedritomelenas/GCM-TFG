import numpy as np

def WeighProd(x, y, sigmas):
    m = x.transpose()*sigmas    #np.multiply(x, sigmas)
    n = m.transpose()*y     #np.multiply(m, y)
    return n.sum(axis=0)