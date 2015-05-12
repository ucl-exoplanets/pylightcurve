import matplotlib.pyplot as plt
import numpy as np
import time
import timeit

import fcmodel
import mcmc

import warnings
warnings.filterwarnings('ignore')


# given limb darkening coefficients
a1ex = 0.4786
a2ex = 0.375
a3ex = -0.5961
a4ex = 0.2407



# or search in claret's tables based on stellar parameters
# metallicity
metall = 0.01
# effective temperature
teff = 6590
# Log(g)
logg = 4.1
# filter
filter = 'J'
# ldcoeff function: returns a tuple ( a1, a2, a3, a4 )
a1ex, a2ex, a3ex, a4ex = fcmodel.ldcoeff(metall, teff, logg, filter)



# transit parameters
# period [ days ]
prex = 0.15
# planet radius / stellar radius [ no units]
pex = 2.2
# semi-major axis / stellar radius [ no units ]
aex = 9.0
# eccentricity [ no units ]
eex = 0.0
# inclination [ degrees ]
iex = 87.5
# argument of periastron [ degrees ]
wex = 0
# mid-transit point [ days, same as time array ]
t0ex = 100.0
# time array [ days, same as mid-transit point ]
ttex = np.arange(t0ex - pex / 20, t0ex + pex / 20, 0.0005)
# Big Omega, if not given default value = 0 
wwex = 0



# model function: returns an numpy array ( relative flux )
light_curve = fcmodel.model( ( a1ex, a2ex, a3ex, a4ex ), prex, pex, aex, eex, iex, wex, t0ex, ttex )



# test plot ant runtime
plt.plot((ttex-t0ex)/pex, light_curve, 'k-', lw=2)
plt.xlabel(r'$ phase $')
plt.ylabel(r'$ relative \, flux $')
plt.ylim((plt.ylim()[0], 1.002))

iter = 1000
benchtime = timeit.timeit("fcmodel.model((a1ex, a2ex, a3ex, a4ex), prex, pex, aex, eex, iex, wex, t0ex, ttex)", setup="from __main__ import *", number=iter)
print "\n{:.5}ms per lightcurve ( from {} iterations )\n".format(benchtime/iter*1000,iter)

plt.show()



# create a fake observation with noise to test fitting
noise = 0.0005
observation_time = np.arange(t0ex - pex / 20, t0ex + pex / 20, 0.0005)
observation_flux = fcmodel.model((a1ex, a2ex, a3ex, a4ex), prex, pex, aex, eex, iex, wex, t0ex, observation_time)
observation_flux = observation_flux + np.random.normal(loc=0.0, scale=noise, size=len(observation_flux))

# fitting parameters
# number of iterations
iter = 20000
# number of iterations to be excluded from the beginning 
burn = 5000
# planet radius / stellar radius range [ +/- % ] 
prvar = 0.1
# semi-major axis / stellar radius range [ +/- % ]
avar = 0.3
# inclination range [ +/- % ]
ivar = 0.05


# model_fit function: creates a directory 'mcmc_fit' with all the mcmc fitting results
# if there is such a directory the new one will be 'mcmc_fit_2' or 'mcmc_fit_3' etc ...
mcmc.model_fit( ( observation_time, observation_flux ), iter, burn, ( a1ex, a2ex, a3ex, a4ex ), prex, prvar, pex, aex, avar, eex, iex, ivar, wex, t0ex, wwex )

# alternatively, if the observation is stored in a file 'file_name'
# containing time as first column and relative flux as second column
# the model_fit function can be used as follows:
# mcmc.model_fit( 'test.txt', iter, burn, ( a1ex, a2ex, a3ex, a4ex ), prex, prvar, pex, aex, avar, eex, iex, ivar, wex, t0ex )
