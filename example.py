import matplotlib.pyplot as plt
import numpy as np

import fcmodel


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

x1 = np.arange(t0ex - pex / 20, t0ex + pex / 20, 0.0005)[::2]
x2 = np.arange(t0ex - pex / 20, t0ex + pex / 20, 0.0005)[1::2]
y1 = fcmodel.model((a1ex, a2ex, a3ex, a4ex), prex, pex, aex, eex, iex, wex, wwex, t0ex, x1)
y2 = fcmodel.model(fcmodel.ldcoeff(metall, teff, logg, filter), prex, pex, aex, eex, iex, wex, wwex, t0ex, x2)

plt.plot(x1, y1, 'ko')
plt.plot(x2, y2, 'ro')
plt.xlabel(r'$time\,(t)\,[days]$')
plt.ylabel(r'$relative\,flux\,(f(t))$')
plt.ylim((plt.ylim()[0], 1.002))
plt.xlim((x1[0], x2[-1]))
plt.show()

# wrong filter test
y2 = fcmodel.model(fcmodel.ldcoeff(metall, teff, logg, 'a'), prex, pex, aex, eex, iex, wex, wwex, t0ex, x2)
