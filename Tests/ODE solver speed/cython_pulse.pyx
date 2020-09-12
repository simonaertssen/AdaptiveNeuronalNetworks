from numpy import cos, power

def cython_pulse(theta):
    return power(1 - cos(theta), 2)
