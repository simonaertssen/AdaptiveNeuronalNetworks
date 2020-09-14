import sys
import time

import numpy as np
from c_lib import c_functions
#c_functions = ctypes.CDLL('lib/c_functions.cpython-37m-darwin.so')
#c_functions.c_a_n.argtype = (ctypes.c_int)

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from scipy.stats import cauchy

def compute():
    # This main function has been adapted to use as many c functions where possible!
    # Other functions are gathered in the pxd declarator
    # Parameters
    start = time.time()
    pars = {}
    # for _ in range(1000000):
    #     c_functions.c_a_n(2)
    print(c_functions.factorial(4))
    # print(c_functions.c_a_n(2))
    # pars["a_n"] = c_functions.c_a_n(2)
    pars["N"] = 1

    print("time: ", time.time() - start)
    # Takes about 0.3462 seconds so not the fastest...

if __name__ == '__main__':
    compute()
