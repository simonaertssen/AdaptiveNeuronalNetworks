from numpy import cos, power
from cython_pulse import cython_pulse

def cython_thetaneurons(t, x, e, KdivN, a):
    # Model a network of oscillators
    I_sync = cython_pulse(x).sum()
    return (1 - cos(x)) + (1 + cos(x)) * (e + a * KdivN * I_sync);
