from torch import cos
from pytorch_pulse import pytorch_pulse

def pytorch_thetaneurons(t, x, e, KdivN, a):
    # Model a network of oscillators
    I_sync = pytorch_pulse(x).sum()
    return (1 - cos(x)).add(1 + cos(x)).matmul(e + a * KdivN * I_sync);
