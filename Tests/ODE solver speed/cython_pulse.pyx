from numpy import cos, power

cdef cython_pulse(theta):
    return power(1 - cos(theta), 2)
