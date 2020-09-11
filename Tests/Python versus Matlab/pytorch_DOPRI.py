from torch import linspace, zeros
from pytorch_DOPRIstep import pytorch_DOPRIstep

def pytorch_DOPRI(originalfunc,ta,tb,x0,h,p):
    npts = int(round((tb - ta)/h + 1))
    h = (tb - ta)/(npts-1)
    dim = x0.shape

    cuda = torch.cuda.is_available()
    if cuda:
        print("cuda session enabled")
        device = torch.device("cuda")
        torch.set_default_tensor_type('torch.cuda.FloatTensor')
    else:
        print("cpu session enabled")
        device = torch.device("cpu")
        torch.set_default_tensor_type('torch.FloatTensor')

    xout = zeros(dim[0], npts).to(device)
    xout[:,0] = x0

    tout = linspace(ta,tb,npts);
    # Make new function handle to improve speed!
    def func(t, x, e=p["e"], KN=p["K"]/p["N"], a_n=p["a_n"]):
        return originalfunc(t, x, e, KN, a_n)

    for i in range(npts-1):
        xout[:,i+1] = pytorch_DOPRIstep(func,tout[i],xout[:,i],h);
    return tout, xout
