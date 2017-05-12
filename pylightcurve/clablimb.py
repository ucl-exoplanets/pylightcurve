__all__ = ['clablimb']


import os
import glob
import numpy as np

data_base_location = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'clablimb_data_base')


def tri_linear(x, y, z, x0, x1, y0, y1, z0, z1, v000, v100, v010, v001, v101, v011, v110, v111):
    c0 = v000
    c1 = v100 - v000
    c2 = v010 - v000
    c3 = v001 - v000
    c4 = v110 - v010 - v100 + v000
    c5 = v011 - v001 - v010 + v000
    c6 = v101 - v001 - v100 + v000
    c7 = v111 - v011 - v101 - v110 + v100 + v001 + v010 - v000
    dx = (x - x0) / (x1 - x0)
    dy = (y - y0) / (y1 - y0)
    dz = (z - z0) / (z1 - z0)
    return c0 + c1 * dx + c2 * dy + c3 * dz + c4 * dx * dy + c5 * dy * dz + c6 * dz * dx + c7 * dx * dy * dz


def clablimb(method, stellar_logg, stellar_temperature, stellar_metallicity, photometric_filter, stellar_model='ATLAS'):

    if method != 'claret':
        raise ImportError('Limb darkening model not currently supported.')

    data_file = glob.glob(os.path.join(data_base_location, '*_' + stellar_model + '_' + photometric_filter + '.txt'))[0]
    data = np.loadtxt(data_file, usecols=[0, 1, 2, 3, 4, 5, 6, 7], unpack=True)

    x = np.unique(data[0])
    y = np.unique(data[1])
    z = np.unique(data[2])

    yin = float(stellar_temperature)
    zin = float(stellar_metallicity)
    xin = float(stellar_logg)

    test = np.argmin(np.abs(np.ones_like(x) * xin - x))
    xmin = x[max(0, test - 1)]
    xmax = x[min(len(x) - 1, test + 1)]
    test = np.argmin(np.abs(np.ones_like(y) * yin - y))
    ymin = y[max(0, test - 1)]
    ymax = y[min(len(y) - 1, test + 1)]
    test = np.argmin(np.abs(np.ones_like(z) * zin - z))
    zmin = z[max(0, test - 1)]
    zmax = z[min(len(z) - 1, test + 1)]

    final_coefficients = []

    for index in [4, 5, 6, 7]:

        v000 = data[index][np.where((data[0] == xmin) & (data[1] == ymin) & (data[2] == zmin))][0]
        v100 = data[index][np.where((data[0] == xmax) & (data[1] == ymin) & (data[2] == zmin))][0]
        v010 = data[index][np.where((data[0] == xmin) & (data[1] == ymax) & (data[2] == zmin))][0]
        v001 = data[index][np.where((data[0] == xmin) & (data[1] == ymin) & (data[2] == zmax))][0]
        v101 = data[index][np.where((data[0] == xmax) & (data[1] == ymin) & (data[2] == zmax))][0]
        v011 = data[index][np.where((data[0] == xmin) & (data[1] == ymax) & (data[2] == zmax))][0]
        v110 = data[index][np.where((data[0] == xmax) & (data[1] == ymax) & (data[2] == zmin))][0]
        v111 = data[index][np.where((data[0] == xmax) & (data[1] == ymax) & (data[2] == zmax))][0]

        res = tri_linear(xin, yin, zin, xmin, xmax, ymin, ymax, zmin, zmax,
                         v000, v100, v010, v001, v101, v011, v110, v111)

        final_coefficients.append(res)

    return np.array(final_coefficients)