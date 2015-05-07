import matplotlib.pyplot as plt
import numpy as np
import os
import shutil
import time
import timeit

import fcmodel
import mcmc


prex = 0.15
pex = 2.2
aex = 9.0
eex = 0.0
iex = 87.5
wex = 0
wwex = 0
t0ex = 100.0
iex = 87.5
wwex = 0
wex = 0

a1ex = 0.4786
a2ex = 0.375
a3ex = -0.5961
a4ex = 0.2407

metall = 0.01
teff = 6590
logg = 4.1
filter = 'J'

x = np.arange(t0ex - pex / 20, t0ex + pex / 20, 0.001)

y = fcmodel.model((a1ex, a2ex, a3ex, a4ex), prex, pex, aex, eex, iex, wex, t0ex, x, wwex)
y = y + np.random.normal(loc=0.0, scale=0.0005, size=len(y))

if os.path.isdir('mcmc_fit'):
	shutil.rmtree('mcmc_fit')

iter = 20000
burn = 5000

# the prior for each parameter is a uniform [ Value - var*Valur, Value, Value + var*Valur]
prvar = 0.1
avar  = 0.3
ivar  = 0.05


# use (input data x, input data y) as first argument
mcmc.model_fit((x,y),iter,burn,(a1ex, a2ex, a3ex, a4ex), prex, prvar, pex, aex, avar, eex, iex, ivar, wex, t0ex, wwex)

# alternatively 'file name' as first argument
# the file should contain input data x as first column and input data y as second column
# np.savetxt('test.txt',np.swapaxes([x,y],0,1))
# mcmc.model_fit('test.txt',iter,burn,(a1ex, a2ex, a3ex, a4ex), prex, prvar, pex, aex, avar, eex, iex, ivar, wex, t0ex, wwex)
# os.system("rm test.txt")



