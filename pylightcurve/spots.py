from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from ._1databases import *


def theta_to_phi(th, th0, ph0, a):

    def model(th_m, ph_m):

        return np.cos(a) - (np.cos(th_m) * np.sin(ph_m) * np.cos(th0) * np.sin(ph0) + np.sin(th_m) * np.sin(ph_m) * np.sin(th0) * np.sin(ph0) + np.cos(ph_m) * np.cos(ph0))

    ppov, pcov = curve_fit(model, th, 0)

    return ppov[0]

x=[]
y=[]
z=[]
for th in range(0, 360):
    ph = theta_to_phi(th, 0, np.pi/2, np.pi/10)
    x.append(np.cos(th) * np.sin(ph))
    y.append(np.sin(th) * np.sin(ph))
    z.append(np.cos(ph))

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x,y,z)
plt.show()
# plt.plot(y, z, 'o')
# plt.show()