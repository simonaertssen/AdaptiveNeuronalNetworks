from torch import linspace, zeros, cuda, device
from pytorch_DOPRIstep import pytorch_DOPRIstep

def pytorch_DOPRI(originalfunc,ta,tb,x0,h,p):
    npts = int(round((tb - ta)/h + 1))
    h = (tb - ta)/(npts-1)
    dim = x0.shape

    cuda_available = cuda.is_available()
    if cuda_available:
        print("cuda session enabled")
        device = device("cuda")
    else:
        print("cpu session enabled")
        device = device("cpu")

    xout = zeros(dim[0], npts).to(device)
    xout[:,0] = x0

    tout = linspace(ta,tb,npts).to(device)
    # Make new function handle to improve speed!
    def func(t, x, e=p["e"], KN=p["K"]/p["N"], a_n=p["a_n"]):
        return originalfunc(t, x, e, KN, a_n)

    for i in range(npts-1):
        xout[:,i+1] = pytorch_DOPRIstep(func,tout[i],xout[:,i],h);
    return tout, xout
