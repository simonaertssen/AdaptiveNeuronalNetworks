import sys
import time

start = time.time()

import numpy as np
cimport numpy as np

from cython_a_n cimport cython_a_n

# For the pure c++ files we need a cython declaration:
cdef extern from "c_a_n.h":
    double c_a_n(int)

cpdef testtiming_cython():
  # Test speed improvement by calling cython functions
  start = time.time()
  cdef int n = 2
  for _ in range(1000000):
    cython_a_n(n)
  print(cython_a_n(n))
  print("time: ", time.time() - start)
  # Takes about 0.65s unoptimised
  # Takes about 0.038s using cpdef for a_n and cdef for factorial
  # Takes about 0.034s wrapping with @cython.boundscheck(False) and @cython.wraparound(False)
  # Takes about 0.045s using from libc.math cimport pow
  # Takes about 0.0063s using .pxd for all functions

cpdef testtiming():
  # Test speed improvement by calling c++ functions
  start = time.time()
  cdef int n = 2
  for _ in range(1000000):
    c_a_n(n)
  print(c_a_n(n))
  print("time: ", time.time() - start)
  # Takes about 0.0034s unoptimised

from cython_DOPRI import cython_DOPRI

cdef extern from "cython_c_functions.h":
    double c_a_n(int)

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from scipy.stats import cauchy

cpdef main():
  # This main function has been adapted to use as many c functions where possible!
  # Other functions are gathered in the pxd declarator
  # Parameters
  # pars = {}
  # cdef int n = 2
  # pars["a_n"] = c_a_n(n)
  # pars["N"] = 1
  #
  # cdef float t = 1.0
  # cdef double F = c_thetaneurons(t)
  # print(F)
  # tnow = 0
  # tend = 10
  # h = 0.005
  # IC = np.random.randn(pars["N"])*0.9 + 1
  #
  # pars["a_n"] = cython_a_n(2)
  # pars["eta0"] = 10.75
  # pars["delta"] = 0.5
  # pars["K"] = -9
  # seed = 0
  # pars["e"] = cauchy.rvs(random_state=seed, loc=pars["eta0"], scale=pars["delta"], size=pars["N"]);
  #
  # # t, x = cython_DOPRI(F, tnow, tend, IC, h, pars)
  # # tnew = np.vstack([t] * pars["N"])
  # # data = np.stack((tnew,x), axis=2)
  # #
  # # fig, ax = plt.subplots()
  # # ax.add_collection(LineCollection(data))
  # # ax.set_ylim([x.min(1).min(), x.max(1).max()])
  # #plt.show()
  # print(time.time() - start)

# Elapsed time is 21.666523933410645 seconds (changed nothing, just ran with .pyx instead of py)
