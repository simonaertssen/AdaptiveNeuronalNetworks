from cython cimport cdivision, boundscheck, wraparound, locals

import numpy as np
cimport numpy as np

np.import_array()

cimport scipy.linalg.cython_blas as blas
from libc.math cimport cos

cdef np.ndarray[np.float64_t, ndim=1] cython_pulse(np.ndarray[np.float64_t, ndim=1] theta):
  cdef int thetalength = theta.size
  # cdef tuple shape = (<object> theta).shape
  # cdef int dimension0 = shape[0]
  for i in range(thetalength):
    theta[i] = (1 - cos(theta[i]))**2
  return theta


cdef np.ndarray[np.float64_t, ndim=1] cython_thetaneurons(double t, np.ndarray[np.float64_t, ndim=1] x, np.ndarray[np.float64_t, ndim=1] e, double KdivN, double a):
  cdef int xlength = x.size, i
  cdef np.ndarray[np.float64_t, ndim=1] sumarray = cython_pulse(x)
  cdef double sum
  for i in range(xlength):
    sum += sumarray[i]

  cdef double mult = (a * KdivN * sum)
  for i in range(xlength):
    x[i] = (1 - cos(x[i])) + (1 + cos(x[i])) * (e[i] + mult)
  return x


# @cdivision(True)
cdef np.ndarray[np.float64_t, ndim=1] cython_DOPRIstep(double t, np.ndarray[np.float64_t, ndim=1] x, double h, np.ndarray[np.float64_t, ndim=1] e, double kn, double an):
  cdef np.ndarray[np.float64_t, ndim=1] K1 = cython_thetaneurons(t,x,e,kn,an).dot(h)
  cdef np.ndarray[np.float64_t, ndim=1] K2 = cython_thetaneurons(t+h/5,    x + K1.dot(0.2), e,kn,an).dot(h)
  cdef np.ndarray[np.float64_t, ndim=1] K3 = cython_thetaneurons(t+3*h/10, x + K1.dot(0.075) + K2.dot(0.225), e,kn,an).dot(h)
  cdef np.ndarray[np.float64_t, ndim=1] K4 = cython_thetaneurons(t+4*h/5,  x + K1.dot(44/45) - K2.dot(53/15) + K3.dot(32/9), e,kn,an).dot(h)
  cdef np.ndarray[np.float64_t, ndim=1] K5 = cython_thetaneurons(t+8*h/9,  x + K1.dot(19372/6561) - K2.dot(25360/2187) + K3.dot(64448/6561) - K4.dot(212/729), e,kn,an).dot(h)
  cdef np.ndarray[np.float64_t, ndim=1] K6 = cython_thetaneurons(t+h,      x + K1.dot(9017/3168)  - K2.dot(355/33)     + K3.dot(46732/5247) + K4.dot(49/176) - K5.dot(5103/18656), e,kn,an).dot(h)
  return x + K1.dot(35/384) + K3.dot(500/1113) + K4.dot(125/192) - K5.dot(2187/6784) + K6.dot(11/84)

# cdef np.ndarray[np.float64_t, ndim=1] cython_DOPRIstep(double t, np.ndarray[np.float64_t, ndim=1] x, double h, double e, double kn, double an):
#     K1 = h*cython_thetaneurons(t,x,e,kn,an)
#     K2 = h*cython_thetaneurons(t+h/5,x+K1/5,e,kn,an)
#     K3 = h*cython_thetaneurons(t+3*h/10,x+3*K1/40+9*K2/40,e,kn,an)
#     K4 = h*cython_thetaneurons(t+4*h/5,x+44*K1/45-56*K2/15+32*K3/9,e,kn,an)
#     K5 = h*cython_thetaneurons(t+8*h/9,x+19372*K1/6561-25360*K2/2187+64448*K3/6561-212*K4/729,e,kn,an)
#     K6 = h*cython_thetaneurons(t+h,x+9017*K1/3168-355*K2/33+46732*K3/5247+49*K4/176-5103*K5/18656,e,kn,an)
#     return x + 35*K1/384+500*K3/1113+125*K4/192-2187*K5/6784+11*K6/84


# cdef class FunctionHandle:
#     def __init__(self, handle, double e, double kn, double an):
#         super(FunctionHandle, self).__init__(handle, e, kn, an)
#         self.handle = handle
#         self.e  = e
#         self.kn = kn
#         self.an = an
#
#     # @boundscheck(False)
#     # @cdivision(True)
#     # @wraparound(False)
#     # @nonecheck(False)
#     cdef void eval(self, double t, double x):
#         self.handle(t, x, self.e, self.kn, self.an)


@cdivision(True)
cdef np.ndarray[np.float64_t, ndim=2] cython_DOPRI(double ta, double tb, np.ndarray[np.float64_t, ndim=1] x0, double h, dict pars):
    cdef int npts = int(round((tb - ta)/h + 1))
    cdef double htemp = (tb - ta)/(npts-1)
    if h != htemp:
      h = htemp
      print("Setting h to", h)
    cdef int xlength = x0.size
    cdef np.ndarray[np.float64_t, ndim=2] xout = np.zeros((xlength, npts))
    xout[:,0] = x0

    # cdef np.ndarray[np.float64_t, ndim=1] tout = linspace(ta,tb,npts);
    # Make new function handle to improve speed!
    #cdef double KN =p["K"]/p["N"]
    # cdef inline np.ndarray[np.float64_t, ndim=1] func(doub t, np.ndarray[np.float64_t, ndim=1] x, e=p["e"], KN=p["K"]/p["N"], a_n=p["a_n"]):
    #     originalfunc(t, x, e, KN, a_n)

    # cdef FunctionHandle handle = FunctionHandle(cython_thetaneurons, e=pars["e"], KN=pars["K"]/pars["N"], a_n=pars["a_n"])
    # cdef int i
    # for i in range(npts-1):
    #     xout[:,i+1] = cython_DOPRIstep(tout[i], xout[:,i], h, pars["e"], pars["K"]/pars["N"], pars["a_n"]);
    return xout
