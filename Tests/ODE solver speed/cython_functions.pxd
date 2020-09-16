cdef class FunctionHandle:
    cdef double e
    cdef double kn
    cdef double an

cdef double[:] cython_pulse(double[:] theta)
cdef double[:] cython_thetaneurons(double t, double[:] x, double[:] e, double KdivN, double a)
cdef double[:] numpy_DOPRIstep(func, double t, double[:] x, double h)
cdef double[:,:] cython_DOPRI(F, double tnow, double tend, double[:] IC, double h, dict pars)
