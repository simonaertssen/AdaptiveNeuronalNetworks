
cimport numpy as np
import numpy as np

cdef np.ndarray[np.float64_t, mode="fortran", ndim=1] DSCAL(np.ndarray[np.float64_t, mode="fortran", ndim=1] A, double scalar)
cdef np.ndarray[np.float64_t, mode="fortran", ndim=1] cython_pulse(np.ndarray[np.float64_t, mode="fortran", ndim=1] theta)
cdef np.ndarray[np.float64_t, mode="fortran", ndim=1] cython_thetaneurons(double t, np.ndarray[np.float64_t, mode="fortran", ndim=1] x, np.ndarray[np.float64_t, mode="fortran", ndim=1] e, double KdivN, double a)
cdef np.ndarray[np.float64_t, mode="fortran", ndim=1] cython_DOPRIstep(double t, np.ndarray[np.float64_t, mode="fortran", ndim=1] x, double h, np.ndarray[np.float64_t, mode="fortran", ndim=1] e, double kn, double an)
cdef np.ndarray[np.float64_t, mode="fortran", ndim=2] cython_DOPRI(double tnow, double tend, np.ndarray[np.float64_t, mode="fortran", ndim=1] IC, double h, dict pars)
