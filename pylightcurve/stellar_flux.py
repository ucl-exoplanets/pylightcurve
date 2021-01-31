
from .tools_files import *


def clablimb(method, stellar_logg, stellar_temperature, stellar_metallicity, photometric_filter, stellar_model='ATLAS'):

    if method != 'claret':
        raise ImportError('Limb darkening model not currently supported.')

    data_file = glob.glob(os.path.join(databases.clablimb,
                                       '*_{0}_{1}.txt'.format(stellar_model, photometric_filter)))[0]
    data = np.loadtxt(data_file, usecols=[0, 1, 2, 3, 4, 5, 6, 7], unpack=True)

    xx = np.unique(data[0])
    yy = np.unique(data[1])
    zz = np.unique(data[2])

    yin = float(stellar_temperature)
    zin = float(stellar_metallicity)
    xin = float(stellar_logg)

    if xin in xx:
        xmin, xmax = xin, xin
    else:
        test = np.argmin(np.abs(np.ones_like(xx) * xin - xx))
        if xx[test] < xin:
            xmin = xx[test]
            xmax = xx[min(len(xx) - 1, test + 1)]
        else:
            xmin = xx[max(0, test - 1)]
            xmax = xx[test]
    if yin in yy:
        ymin, ymax = yin, yin
    else:
        test = np.argmin(np.abs(np.ones_like(yy) * yin - yy))
        if yy[test] < yin:
            ymin = yy[test]
            ymax = yy[min(len(yy) - 1, test + 1)]
        else:
            ymin = yy[max(0, test - 1)]
            ymax = yy[test]
    if zin in zz:
        zmin, zmax = zin, zin
    else:
        test = np.argmin(np.abs(np.ones_like(zz) * zin - zz))
        if zz[test] < zin:
            zmin = zz[test]
            zmax = zz[min(len(zz) - 1, test + 1)]
        else:
            zmin = zz[max(0, test - 1)]
            zmax = zz[test]

    def tri_linear(x, y, z, x0, x1, y0, y1, z0, z1, v000, v100, v010, v001, v101, v011, v110, v111):
        c0 = v000
        c1 = v100 - v000
        c2 = v010 - v000
        c3 = v001 - v000
        c4 = v110 - v010 - v100 + v000
        c5 = v011 - v001 - v010 + v000
        c6 = v101 - v001 - v100 + v000
        c7 = v111 - v011 - v101 - v110 + v100 + v001 + v010 - v000
        if x == x0 == x1:
            dx = 0
        else:
            dx = (x - x0) / (x1 - x0)
        if y == y0 == y1:
            dy = 0
        else:
            dy = (y - y0) / (y1 - y0)
        if z == z0 == z1:
            dz = 0
        else:
            dz = (z - z0) / (z1 - z0)
        return c0 + c1 * dx + c2 * dy + c3 * dz + c4 * dx * dy + c5 * dy * dz + c6 * dz * dx + c7 * dx * dy * dz

    final_coefficients = []

    for index in [4, 5, 6, 7]:

        vv000 = data[index][np.where((data[0] == xmin) & (data[1] == ymin) & (data[2] == zmin))][0]
        vv100 = data[index][np.where((data[0] == xmax) & (data[1] == ymin) & (data[2] == zmin))][0]
        vv010 = data[index][np.where((data[0] == xmin) & (data[1] == ymax) & (data[2] == zmin))][0]
        vv001 = data[index][np.where((data[0] == xmin) & (data[1] == ymin) & (data[2] == zmax))][0]
        vv101 = data[index][np.where((data[0] == xmax) & (data[1] == ymin) & (data[2] == zmax))][0]
        vv011 = data[index][np.where((data[0] == xmin) & (data[1] == ymax) & (data[2] == zmax))][0]
        vv110 = data[index][np.where((data[0] == xmax) & (data[1] == ymax) & (data[2] == zmin))][0]
        vv111 = data[index][np.where((data[0] == xmax) & (data[1] == ymax) & (data[2] == zmax))][0]

        res = tri_linear(xin, yin, zin, xmin, xmax, ymin, ymax, zmin, zmax,
                         vv000, vv100, vv010, vv001, vv101, vv011, vv110, vv111)

        final_coefficients.append(res)

    return np.array(final_coefficients)
