import sys
import time

start = time.time()

import torch

cuda = torch.cuda.is_available()
if cuda:
    print("cuda session enabled")
    device = torch.device("cuda")
    torch.set_default_tensor_type('torch.cuda.FloatTensor')
else:
    print("cpu session enabled")
    device = torch.device("cpu")
    torch.set_default_tensor_type('torch.FloatTensor')

import numpy as np
from pytorch_a_n import pytorch_a_n
from pytorch_DOPRI import pytorch_DOPRI
from pytorch_thetaneurons import pytorch_thetaneurons

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from scipy.stats import cauchy

# Parameters
pars = {}
pars["a_n"] = pytorch_a_n(2)
pars["N"] = 10000
F = pytorch_thetaneurons
tnow = 0
tend = 10
h = 0.005
IC = torch.randn(pars["N"])*0.9 + 1

pars["a_n"] = pytorch_a_n(2)
pars["eta0"] = 10.75
pars["delta"] = 0.5
pars["K"] = -9
seed = 0
pars["e"] = torch.distributions.Cauchy(loc=pars["eta0"], scale=pars["delta"]).sample((pars["N"],))
#pars["e"] = torch.tensor(cauchy.rvs(random_state=seed, loc=pars["eta0"], scale=pars["delta"], size=pars["N"]));

t, x = pytorch_DOPRI(F, tnow, tend, IC, h, pars)
tnew = np.vstack([t.numpy()] * pars["N"])
data = np.stack((tnew,x.numpy()), axis=2)

fig, ax = plt.subplots()
ax.add_collection(LineCollection(data))
ax.set_ylim([x.min(), x.max()])
#plt.show()
print(time.time() - start)

# % Elapsed time is 8.17128300666809 seconds.
