from cython cimport cdivision, boundscheck, wraparound, locals

import numpy as np
cimport numpy as np
np.import_array()

cimport scipy.linalg.cython_blas as blas
from libc.math cimport cos

import time

cdef np.ndarray[np.float64_t, mode="fortran", ndim=1] DSCAL(np.ndarray[np.float64_t, mode="fortran", ndim=1] A, double scalar):
  cdef int Asize = A.size
  cdef double *alpha = &scalar
  cdef int increment = 1
  blas.dscal(&Asize, alpha, &A[0], &increment);
  return A

cdef np.ndarray[np.float64_t, mode="fortran", ndim=1] cython_pulse(np.ndarray[np.float64_t, mode="fortran", ndim=1] theta):
  cdef int thetalength = theta.size
  # cdef tuple shape = (<object> theta).shape
  # cdef int dimension0 = shape[0]
  for i in range(thetalength):
    theta[i] = (1 - cos(theta[i]))**2
  return theta


cdef np.ndarray[np.float64_t, mode="fortran", ndim=1] cython_thetaneurons(double t, np.ndarray[np.float64_t, mode="fortran", ndim=1] x, np.ndarray[np.float64_t, mode="fortran", ndim=1] e, double KdivN, double a):
  cdef int xlength = x.size, i
  cdef np.ndarray[np.float64_t, mode="fortran", ndim=1] sumarray = cython_pulse(x)
  cdef double sum
  for i in range(xlength):
    sum += sumarray[i]

  cdef double mult = (a * KdivN * sum)
  for i in range(xlength):
    x[i] = (1 - cos(x[i])) + (1 + cos(x[i])) * (e[i] + mult)
  return x


# @cdivision(True)
cdef np.ndarray[np.float64_t, mode="fortran", ndim=1] cython_DOPRIstep(double t, np.ndarray[np.float64_t, mode="fortran", ndim=1] x, double h, np.ndarray[np.float64_t, mode="fortran", ndim=1] e, double kn, double an):
  cdef np.ndarray[np.float64_t, mode="fortran", ndim=1] K1 = cython_thetaneurons(t,x,e,kn,an).dot(h)
  cdef np.ndarray[np.float64_t, mode="fortran", ndim=1] K2 = cython_thetaneurons(t+h/5,    x + K1.dot(0.2), e,kn,an).dot(h)
  cdef np.ndarray[np.float64_t, mode="fortran", ndim=1] K3 = cython_thetaneurons(t+3*h/10, x + K1.dot(0.075) + K2.dot(0.225), e,kn,an).dot(h)
  cdef np.ndarray[np.float64_t, mode="fortran", ndim=1] K4 = cython_thetaneurons(t+4*h/5,  x + K1.dot(44/45) - K2.dot(53/15) + K3.dot(32/9), e,kn,an).dot(h)
  cdef np.ndarray[np.float64_t, mode="fortran", ndim=1] K5 = cython_thetaneurons(t+8*h/9,  x + K1.dot(19372/6561) - K2.dot(25360/2187) + K3.dot(64448/6561) - K4.dot(212/729), e,kn,an).dot(h)
  cdef np.ndarray[np.float64_t, mode="fortran", ndim=1] K6 = cython_thetaneurons(t+h,      x + K1.dot(9017/3168)  - K2.dot(355/33)     + K3.dot(46732/5247) + K4.dot(49/176) - K5.dot(5103/18656), e,kn,an).dot(h)
  return x + K1.dot(35/384) + K3.dot(500/1113) + K4.dot(125/192) - K5.dot(2187/6784) + K6.dot(11/84)


@cdivision(True)
cdef np.ndarray[np.float64_t, mode="fortran", ndim=2] cython_DOPRI(double ta, double tb, np.ndarray[np.float64_t, mode="fortran", ndim=1] x0, double h, dict pars):
    cdef int npts = int(round((tb - ta)/h + 1))
    cdef double htemp = (tb - ta)/(npts-1)
    if h != htemp:
      h = htemp
      print("Setting h to", h)
    cdef int xlength = x0.size
    cdef np.ndarray[np.float64_t, mode="fortran", ndim=2] xout = np.zeros((xlength, npts), order='F')
    xout[:,0] = x0

    cdef np.ndarray[np.float64_t, mode="fortran", ndim=1] tout = np.linspace(ta,tb,npts);
    cdef int i
    cdef double timingresult = 0
    for i in range(npts-1):
        # start = time.time()
        xout[:,i+1] = cython_DOPRIstep(tout[i], xout[:,i], h, pars["e"], pars["K"]/pars["N"], pars["a_n"]);
        # timingresult = timingresult * (npts-1)/npts + (time.time() - start)/npts
    # print("Average time =", timingresult)
    return xout
