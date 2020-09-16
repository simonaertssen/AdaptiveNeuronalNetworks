from cython cimport cdivision, boundscheck, wraparound

cimport numpy as np
from np import zeros, linspace
np.import_array()

from libc.math cimport cos

# @boundscheck(False)
# @wraparound(False)
# cdef np.ndarray[double, ndim=1] cython_pulse(np.ndarray[double, ndim=1] theta):
#   cdef int thetalength = theta.size
#   for i in range(thetalength):
#     theta[i] = (1 - cos(theta[i]))**2
#   return theta

cdef double[:] cython_pulse(double[:] theta):
  cdef int thetalength = theta.size
  # cdef tuple shape = (<object> theta).shape
  # cdef int dimension0 = shape[0]
  for i in range(thetalength):
    theta[i] = (1 - cos(theta[i]))**2
  return theta


cdef double[:] cython_thetaneurons(double t, double[:] x, double[:] e, double KdivN, double a):
  cdef int xlength = x.size
  cdef double mult = (a * KdivN * cython_pulse(x).sum())
  for i in range(xlength):
    x[i] = (1 - cos(x[i,j])) + (1 + cos(x[i])) * (e[i] + mult)
  return x


@cdivision(True)
cdef double[:] cython_DOPRIstep(func, double t, double[:] x, double h):
    cdef double[:] K1 = h*func(t,x)
    cdef double[:] K2 = h*func(t+h/5,x+K1/5)
    cdef double[:] K3 = h*func(t+3*h/10,x+3*K1/40+9*K2/40)
    cdef double[:] K4 = h*func(t+4*h/5,x+44*K1/45-56*K2/15+32*K3/9)
    cdef double[:] K5 = h*func(t+8*h/9,x+19372*K1/6561-25360*K2/2187+64448*K3/6561-212*K4/729)
    cdef double[:] K6 = h*func(t+h,x+9017*K1/3168-355*K2/33+46732*K3/5247+49*K4/176-5103*K5/18656)
    return x + 35*K1/384+500*K3/1113+125*K4/192-2187*K5/6784+11*K6/84


cdef class FunctionHandle:
    cdef double e
    cdef double kn
    cdef double an

    def __init__(self):
        super(FunctionHandle, self).__init__(handle, e, KN, a_n)
        self.handle = handle
        self.e  = e
        self.kn = kn
        self.an = an

    # @boundscheck(False)
    # @wraparound(False)
    cpdef inline void eval(self, t, x):
        self.handle(t, x, self.e, self.kn, self.an)


@cdivision(True)
cdef double[:,:] cython_DOPRI(originalfunc, double ta, double tb, np.ndarray[double] x0, double h, dict pars):
    cdef int npts = int(round((tb - ta)/h + 1))
    cdef double h = (tb - ta)/(npts-1)
    cdef tuple xshape = (<object> x0).shape
    cdef int dimension0 = xshape[0]
    cdef int dimension1 = xshape[1]
    cdef double[:,:] xout = zeros(dimension0, npts)
    xout[:,0] = x0

    cdef double[:] tout = linspace(ta,tb,npts);
    # Make new function handle to improve speed!
    #cdef double KN =p["K"]/p["N"]
    # cdef inline double[:] func(doub t, double[:] x, e=p["e"], KN=p["K"]/p["N"], a_n=p["a_n"]):
    #     originalfunc(t, x, e, KN, a_n)
    cdef FunctionHandle handle = FunctionHandle(e=p["e"], KN=p["K"]/p["N"], a_n=p["a_n"])
    for i in range(npts-1):
        xout[:,i+1] = cython_DOPRIstep(handle.eval,tout[i],xout[:,i],h);
    return tout, xout
