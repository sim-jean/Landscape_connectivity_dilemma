# This file codes functions for new indicators

import numpy as np
import math

shape = np.array([1,2,1,1,0,0,1,2,2,1,0,2,1,0,2,1])
shape2 = np.array([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1])
shape3 = np.array([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0])
def shannon_index(shape):
    """
    Shannon index on land types with pseudo count to avoid 0 values
    :param shape: np.array()
    :return: value in [0,1]
    """
    store_ = []
    for i in range(3):
        if sum(shape == i) ==0 :
            store_.append(0)
        else:
            p_ = sum(shape == i)/len(shape)
            store_.append(p_*math.log2(p_))
    return -sum(store_)
def simpson_index(shape):
    """
    Simpson index on land types
    :param shape: np.array()
    :return: value in [0,1]
    """
    H = 1
    for i in range(3):
        H = H- ((sum(shape==i))/len(shape))**2
    return H

# It'd be good to have all the components parameters

