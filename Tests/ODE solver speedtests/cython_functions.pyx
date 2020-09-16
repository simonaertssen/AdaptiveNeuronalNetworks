from cython cimport cdivision, boundscheck, wraparound, view
from cython.parallel cimport prange

cimport openmp
openmp.omp_set_dynamic(0)
openmp.omp_set_num_threads(4)

import ctypes

import numpy as np
cimport numpy as np
np.import_array()

from scipy.linalg.cython_blas cimport dscal, daxpy
from libc.math cimport cos
from libc.stdlib cimport malloc

import time

cdef double[:] DSCAL(double scalar, double[:] x):
  cdef int xsize = x.size
  cdef double *alpha = &scalar
  cdef int increment = 1
  dscal(&xsize, alpha, &x[0], &increment)
  return x


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
  for i in prange(xlength, nogil=True):
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

cdef double[:] cython_DOPRIstep_loop(double t, double[:] x, double h, double[:] e, double kn, double an):
  cdef double[:] K1, K2, K3, K4, K5, K6, Ktmp
  cdef int xlength = x.size, i

  K1 = cython_thetaneurons(t,x,e,kn,an)
  for i in range(xlength):
    K1[i] *= h

  Ktmp = K1
  for i in range(xlength): # x + K1/5
    Ktmp[i] = x[i] + 0.2*Ktmp[i]
  K2 = cython_thetaneurons(t,Ktmp,e,kn,an)
  for i in range(xlength):
    K2[i] *= h

  Ktmp = K2
  for i in range(xlength): # x + 3*K1/40 + 9*K2/40
    Ktmp[i] = x[i] + 3*K1[i]/40 + 9*K2[i]/40
  K3 = cython_thetaneurons(t,Ktmp,e,kn,an)
  for i in range(xlength):
    K3[i] *= h

  Ktmp = K3
  for i in range(xlength): # x + 44*K1/45 - 56*K2/15 + 32*K3/9
    Ktmp[i] = x[i] + 44*K1[i]/45 - 56*K2[i]/15 + 32*K3[i]/9
  K4 = cython_thetaneurons(t,Ktmp,e,kn,an)
  for i in range(xlength):
    K4[i] *= h

  Ktmp = K4
  for i in range(xlength): # x + 19372*K1/6561 - 25360*K2/2187 + 64448*K3/6561 - 212*K4/729
    Ktmp[i] = x[i] + 19372*K1[i]/6561 - 25360*K2[i]/2187 + 64448*K3[i]/6561 - 212*K4[i]/729
  K5 = cython_thetaneurons(t,Ktmp,e,kn,an)
  for i in range(xlength):
    K5[i] *= h

  Ktmp = K5
  for i in range(xlength): # x + 9017*K1/3168 - 355*K2/33 + 46732*K3/5247 + 49*K4/176 - 5103*K5/18656
    Ktmp[i] = x[i] + 9017*K1[i]/3168 - 355*K2[i]/33 + 46732*K3[i]/5247 + 49*K4[i]/176 - 5103*K5[i]/18656
  K6 = cython_thetaneurons(t,Ktmp,e,kn,an)
  for i in range(xlength):
    K6[i] *= h

  for i in range(xlength): # x + 35*K1/384 + 500*K3/1113 + 125*K4/192 - 2187*K5/6784 + 11*K6/84
    x[i] = x[i] + 35*K1[i]/384 + 500*K3[i]/1113 + 125*K4[i]/192 - 2187*K5[i]/6784 + 11*K6[i]/84

  return x

@boundscheck(False)
cdef double[:] cython_DOPRIstep_optimloop(double t, double[:] x, double h, double[:] e, double kn, double an):
  cdef double[:] K1=x, K2, K3, K4, K5, K6, Ktmp
  cdef double* xptr = &x[0]
  cdef int xlength = x.size, i

  Ktmp = cython_thetaneurons(t,x,e,kn,an)
  for i in prange(xlength, nogil=True): # x + K1/5
    Ktmp[i] *= h
    K1[i] = Ktmp[i]
    Ktmp[i] = xptr[i] + 0.2*Ktmp[i]

  K2 = Ktmp
  Ktmp = cython_thetaneurons(t,Ktmp,e,kn,an)
  for i in prange(xlength, nogil=True): # x + 3*K1/40 + 9*K2/40
    Ktmp[i] *= h
    K2[i] = Ktmp[i]
    Ktmp[i] = xptr[i] + 3*K1[i]/40 + 9*K2[i]/40

  K3 = Ktmp
  Ktmp = cython_thetaneurons(t,Ktmp,e,kn,an)
  for i in prange(xlength, nogil=True): # x + 44*K1/45 - 56*K2/15 + 32*K3/9
    Ktmp[i] *= h
    K2[i] = Ktmp[i]
    Ktmp[i] = xptr[i] + 44*K1[i]/45 - 56*K2[i]/15 + 32*K3[i]/9

  K4 = Ktmp
  Ktmp = cython_thetaneurons(t,Ktmp,e,kn,an)
  for i in prange(xlength, nogil=True): # x + 19372*K1/6561 - 25360*K2/2187 + 64448*K3/6561 - 212*K4/729
    Ktmp[i] *= h
    K2[i] = Ktmp[i]
    Ktmp[i] = xptr[i] + 19372*K1[i]/6561 - 25360*K2[i]/2187 + 64448*K3[i]/6561 - 212*K4[i]/729

  K5 = Ktmp
  Ktmp = cython_thetaneurons(t,Ktmp,e,kn,an)
  for i in prange(xlength, nogil=True): # x + 9017*K1/3168 - 355*K2/33 + 46732*K3/5247 + 49*K4/176 - 5103*K5/18656
    Ktmp[i] *= h
    K2[i] = Ktmp[i]
    Ktmp[i] = xptr[i] + 9017*K1[i]/3168 - 355*K2[i]/33 + 46732*K3[i]/5247 + 49*K4[i]/176 - 5103*K5[i]/18656

  K6 = Ktmp
  for i in prange(xlength, nogil=True): # x + 35*K1/384 + 500*K3/1113 + 125*K4/192 - 2187*K5/6784 + 11*K6/84
    x[i] = x[i] + 35*K1[i]/384 + 500*K3[i]/1113 + 125*K4[i]/192 - 2187*K5[i]/6784 + 11*K6[i]/84
  return x


cdef double[:,:] cython_DOPRI(double ta, double tb, double[:] x0, double h, dict pars):
    cdef int npts = int(round((tb - ta)/h + 1))
    cdef double htemp = (tb - ta)/(npts-1)
    if h != htemp:
      h = htemp
      print("Setting h to", h)
    cdef int xlength = x0.size
    cdef double[:] tout = np.linspace(ta,tb,npts)
    cdef double[:,:] xout = np.zeros((xlength, npts), order='F', dtype=np.double).data
    xout[:,0] = x0
    
    cdef int i
    for i in range(npts-1):
        xout[:,i+1] = cython_DOPRIstep_optimloop(tout[i], xout[:,i], h, pars["e"], pars["K"]/pars["N"], pars["a_n"]);
    return xout
