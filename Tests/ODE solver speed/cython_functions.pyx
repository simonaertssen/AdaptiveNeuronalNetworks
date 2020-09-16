from cython cimport cdivision, boundscheck, wraparound, locals

import numpy as np
cimport numpy as np
np.import_array()

from scipy.linalg.cython_blas cimport dscal, daxpy
from libc.math cimport cos

import time

cdef double[:] DSCAL(double scalar, double[:] A):
  cdef int Asize = A.size
  cdef double *alpha = &scalar
  cdef int increment = 1
  dscal(&Asize, alpha, &A[0], &increment)
  return A


cdef double[:] DAXPY(double scalar, double[:] x, double[:] y):
  cdef int xsize = x.size
  cdef double *alpha = &scalar
  cdef int increment = 1
  daxpy(&xsize, alpha, &x[0], &increment, &y[0], &increment)
  return y


cdef double[:] cython_pulse(double[:] theta):
  cdef int thetalength = theta.size
  for i in range(thetalength):
    theta[i] = (1 - cos(theta[i]))**2
  return theta


cdef double[:] cython_thetaneurons(double t, double[:] x, double[:] e, double KdivN, double a):
  cdef int xlength = x.size, i
  cdef double[:] sumarray = cython_pulse(x)
  cdef double sum
  for i in range(xlength):
    sum += sumarray[i]

  cdef double mult = (a * KdivN * sum)
  for i in range(xlength):
    x[i] = (1 - cos(x[i])) + (1 + cos(x[i])) * (e[i] + mult)
  return x


cdef double[:] cython_DOPRIstep(double t, double[:] x, double h, double[:] e, double kn, double an):
  cdef double[:] K1 = DSCAL(h,cython_thetaneurons(t,x,e,kn,an))
  cdef double[:] K2 = DSCAL(h,cython_thetaneurons(t+h/5,    DAXPY(0.2,K1,x),e,kn,an))
  cdef double[:] K3 = DSCAL(h,cython_thetaneurons(t+3*h/10, DAXPY(1, DAXPY(0.075,K1,x), DSCAL(0.225,K2)),e,kn,an))
  cdef double[:] K4 = DSCAL(h,cython_thetaneurons(t+4*h/5,  DAXPY(1, DAXPY(1, DAXPY(44/45,K1,x), DSCAL(-53/15, K2)), DSCAL(32/9,K3)),e,kn,an))
  cdef double[:] K5 = DSCAL(h,cython_thetaneurons(t+8*h/9,  DAXPY(1, DAXPY(1, DAXPY(1, DAXPY(19372/6561,K1,x), DSCAL(-25360/2187,K2)), DSCAL(64448/6561,K3)), DSCAL(-212/729,K4)),e,kn,an))
  cdef double[:] K6 = DSCAL(h,cython_thetaneurons(t+h,      DAXPY(1, DAXPY(1, DAXPY(1, DAXPY(1, DAXPY(9017/3168,K1,x), DSCAL(-355/33,K2)), DSCAL(46732/5247,K3)), DSCAL(49/176,K4)), DSCAL(-5103/18656,K5)), e,kn,an))
  return DAXPY(1, x, DAXPY(1, DAXPY(1, DAXPY(1, DAXPY(1,DSCAL(35/384,K1), DSCAL(500/1113,K3)), DSCAL(125/192,K4)), DSCAL(-2187/6784,K5)), DSCAL(11/84,K6)))


cdef double[:,:] cython_DOPRI(double ta, double tb, double[:] x0, double h, dict pars):
    cdef int npts = int(round((tb - ta)/h + 1))
    cdef double htemp = (tb - ta)/(npts-1)
    if h != htemp:
      h = htemp
      print("Setting h to", h)
    cdef int xlength = x0.size
    cdef double[:,:] xout = np.zeros((xlength, npts), order='F', dtype=np.double)
    cdef double[:,:] *xtest = np.zeros((xlength, npts), order='F', dtype=np.double)

    xout[:,0] = x0

    cdef double[:] tout = np.linspace(ta,tb,npts);
    cdef int i
    cdef double timingresult = 0
    for i in range(npts-1):
        xout[:,i+1] = cython_DOPRIstep(tout[i], xout[:,i], h, pars["e"], pars["K"]/pars["N"], pars["a_n"]);
    return xout
