from torch import linspace, zeros
from pytorch_DOPRIstep import pytorch_DOPRIstep

def pytorch_DOPRI(originalfunc,ta,tb,x0,h,p):
    npts = round((tb - ta)/h + 1)
    h = (tb - ta)/(npts-1)
    dim = x0.shape
    xout = zeros((dim[0], npts))
    xout[:,0] = x0

    tout = linspace(ta,tb,npts);
    # Make new function handle to improve speed!
    def func(t, x, e=p["e"], KN=p["K"]/p["N"], a_n=p["a_n"]):
        return originalfunc(t, x, e, KN, a_n)

    for i in range(npts-1):
        xout[:,i+1] = pytorch_DOPRIstep(func,tout[i],xout[:,i],h);
    return tout, xout
