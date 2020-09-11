from torch import cos, pow

def pytorch_pulse(theta):
    return pow(1 - cos(theta), 2)
