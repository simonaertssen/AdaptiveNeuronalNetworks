from numba import jit, njit

def DOPRIstep(func,t,x,h):
    K1 = h*func(t,x)
    K2 = h*func(t+h/5,x+K1/5)
    K3 = h*func(t+3*h/10,x+3*K1/40+9*K2/40)
    K4 = h*func(t+4*h/5,x+44*K1/45-56*K2/15+32*K3/9)
    K5 = h*func(t+8*h/9,x+19372*K1/6561-25360*K2/2187+64448*K3/6561-212*K4/729)
    K6 = h*func(t+h,x+9017*K1/3168-355*K2/33+46732*K3/5247+49*K4/176-5103*K5/18656)
    return x + 35*K1/384+500*K3/1113+125*K4/192-2187*K5/6784+11*K6/84
