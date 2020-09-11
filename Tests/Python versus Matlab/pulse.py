from numpy import cos, power
from numba import njit, cfunc

@njit(fastmath=True)
def pulse(theta):
    return power(1 - cos(theta), 2)
