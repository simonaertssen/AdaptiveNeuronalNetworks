from torch import cos, pow
from numba import njit

def pytorch_pulse(theta):
    return pow(1 - cos(theta), 2)
