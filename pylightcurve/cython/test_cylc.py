
import pylightcurve as plc
import time
import numpy as np

def testc(N, projected_distance):
    testcy = np.zeros(len(projected_distance))
    n = len(projected_distance)
    tpy, tcy = 0, 0
    for xx in range(N):
        t0=time.time()
        testpy = plc.transit_flux_drop('claret', [0.5, 0.4, 0.3, 0.2], 0.1, projected_distance)
        tpy += time.time() - t0
    for xx in range(N):
        t0 = time.time()
        plc.pyparallel.cc_transit_flux_drop_claret(testcy, 0.5, 0.4, 0.3, 0.2, 0.1, projected_distance, n)
        tcy += time.time() - t0
    print(100 * tcy/tpy, tcy, tpy)
    print(testpy-testcy)


def test(N, time_array):
    t11, t22 = 0, 0
    for xx in range(N):
        t0=time.time()
        position_vector = plc.exoplanet_orbit(5.0, 3.0, 0.0, 90.0, 0.0, 100, time_array)
        projected_distance = np.where(
        position_vector[0] < 0, 1.0 + 5.0 * 0.1,
        np.sqrt(position_vector[1] * position_vector[1] + position_vector[2] * position_vector[2]))
        t11 += time.time() - t0
        t0 = time.time()
        test = plc.transit_flux_drop('claret', [0.5, 0.4, 0.3, 0.2], 0.1, projected_distance)
        t22 += time.time() - t0
    print(100*t11/(t11+t22), 100*t22/(t11+t22))
    print(t22)
    return projected_distance

# print('')
# time_array = np.arange(100 + 10.0 * 5 - 0.11, 100 + 10.0 * 5 + 0.11, 10.0 / 60.0 / 60.0 / 24.0)
# print(len(time_array))
# print('')
# projected_distance = test(1000, time_array)
# print('')
# testc(1000, projected_distance)

import numpy as np
from pytransit import MandelAgol


t = np.linspace(0.8,1.2,500)
print(t)
k, tt0, p, a, i, e, w = 0.1, 1.01, 4, 8, 0.48*np.pi, 0.0, 0.0
u = [0.25,0.10]

tm = 0
tp = 0
for loop in range(100):
    t0 = time.time()
    m = MandelAgol()
    f1 = m.evaluate(t, k, u, tt0, p, a, i, e, w)
    tp += time.time() - t0
    t0 = time.time()
    f2 = plc.transit('quad', u, k, p, a, e, i*180/np.pi, w, tt0, t)
    tm += time.time() - t0
    f3 = plc.transit('quad', u, k, p, a, e, i * 180 / np.pi, w, tt0, t, precision=6)

f1=np.array([ff[0] for ff in f1])
print('')
print(tm, tp)
import matplotlib.pyplot as plt
plt.plot(t, f1, 'ro')
plt.plot(t, f2, 'bo')
plt.show()
plt.plot(t, f1-f2, 'ro')
plt.plot(t, f1-f3, 'bo')
plt.show()