cimport numpy as np
import numpy as np

cdef double[:] DSCAL(double scalar, double[:] A)
cdef double[:] DAXPY(double scalar, double[:] x, double[:] y)

cdef double[:] cython_pulse(double[:] theta)
cdef double[:] cython_thetaneurons(double t, double[:] x, double[:] e, double KdivN, double a)
cdef double[:] cython_DOPRIstep(double t, double[:] x, double h, double[:] e, double kn, double an)
cdef double[:,:] cython_DOPRI(double tnow, double tend, double[:] IC, double h, dict pars)
