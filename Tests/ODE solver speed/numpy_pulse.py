from numpy import cos, power
from numba import njit

@njit(fastmath=True)
def numpy_pulse(theta):
    return power(1 - cos(theta), 2)
