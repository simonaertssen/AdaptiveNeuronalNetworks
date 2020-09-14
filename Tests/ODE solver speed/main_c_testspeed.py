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
    c_functions.printing(23)
    start = time.time()
    pars = {}

    # Implementation does not work as there is a segmentation fault in the c file...
    # Takes about 0.3462 seconds so not the fastest...

if __name__ == '__main__':
    compute()
