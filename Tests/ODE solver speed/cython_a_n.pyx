from cython import boundscheck
from cython import wraparound

#from libc.math cimport pow

cdef double factorial(int n):
    if n==1 or n==0:
        return 1
    else:
        return n * factorial(n-1)

cdef double factorial_it(int n):
  # SLOWER than recursive implementation for small numbers
    cdef double fac = 1
    for i in range(2, n):
        fac = fac * i
    return fac

@boundscheck(False)
@wraparound(False)
cdef double cython_a_n(int n):
    return (2**n)*(factorial(n)**2)/factorial(2*n)
