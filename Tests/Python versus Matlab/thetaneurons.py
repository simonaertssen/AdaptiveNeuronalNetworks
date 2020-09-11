from numpy import cos, power
from numba import jit, njit, vectorize
from pulse import pulse

@njit
def thetaneurons(t, x, e, KdivN, a):
    # Model a network of oscillators
    I_sync = pulse(x).sum()
    return (1 - cos(x)) + (1 + cos(x)) * (e + a * KdivN * I_sync);
