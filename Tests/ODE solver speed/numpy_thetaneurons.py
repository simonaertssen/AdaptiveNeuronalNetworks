from numpy import cos, power
from numba import njit
from numpy_pulse import numpy_pulse

@njit(fastmath=True)
def numpy_thetaneurons(t, x, e, KdivN, a):
    # Model a network of oscillators
    I_sync = numpy_pulse(x).sum()
    return (1 - cos(x)) + (1 + cos(x)) * (e + a * KdivN * I_sync);
