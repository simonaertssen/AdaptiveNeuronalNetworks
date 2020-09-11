import sys
import time

start = time.time()

import numpy as np
from numpy_a_n import numpy_a_n
from numpy_DOPRI import numpy_DOPRI
from numpy_thetaneurons import numpy_thetaneurons

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from scipy.stats import cauchy

# Parameters
pars = {}
pars["a_n"] = numpy_a_n(2)
pars["N"] = 10000
F = numpy_thetaneurons
tnow = 0
tend = 10
h = 0.005
IC = np.random.randn(pars["N"])*0.9 + 1

pars["a_n"] = numpy_a_n(2)
pars["eta0"] = 10.75
pars["delta"] = 0.5
pars["K"] = -9
seed = 0
pars["e"] = cauchy.rvs(random_state=seed, loc=pars["eta0"], scale=pars["delta"], size=pars["N"]);

t, x = numpy_DOPRI(F, tnow, tend, IC, h, pars)
tnew = np.vstack([t] * pars["N"])
data = np.stack((tnew,x), axis=2)

fig, ax = plt.subplots()
ax.add_collection(LineCollection(data))
ax.set_ylim([x.min(1).min(), x.max(1).max()])
#plt.show()
print(time.time() - start)

# % Elapsed time is 11.177425146102905 seconds.
