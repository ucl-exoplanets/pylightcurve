import math
import os
import glob

import numpy as np

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
pi = np.pi

def position(P, A, E, I, W, WW, T0, tt):
    #
    if W < pi / 2:
        AA = 1.0 * pi / 2 - W
    else:
        AA = 5.0 * pi / 2 - W
    BB = 2 * np.arctan(np.sqrt((1 - E) / (1 + E)) * np.tan(AA / 2))
    if BB < 0:
        BB = 2 * pi + BB
    T0 = T0 - (P / 2.0 / pi) * (BB - E * np.sin(BB))
    #
    M = (tt - T0 - np.int_((tt - T0) / P) * P) * 2.0 * pi / P
    u0 = M

    stop = False
    for i in xrange(10000):  # setting a limit of 1k iterations - arbitrary limit
        u1 = u0 - (u0 - E * np.sin(u0) - M) / (1 - E * np.cos(u0))
        stop = (np.abs(u1 - u0) < 10 ** (-7)).all()
        if stop:
            break
        else:
            u0 = u1

            if math.isnan(stop):
                raise ValueError("Nan produced in loop, check inputs")

    if not stop:
        raise PyLCOptimiseError("Failed to find a solution in 10000 loops")

    vv = 2 * np.arctan(np.sqrt((1 + E) / (1 - E)) * np.tan(u1 / 2))
    #
    rr = A * (1 - (E ** 2)) / (np.ones_like(vv) + E * np.cos(vv))
    AA = np.cos(vv + W)
    BB = np.sin(vv + W)
    X = rr * BB * np.sin(I)
    Y = rr * (-AA * np.cos(WW) + BB * np.sin(WW) * np.cos(I))
    Z = rr * (-AA * np.sin(WW) - BB * np.cos(WW) * np.cos(I))
    return [X, Y, Z]




def ldcoeff(Z, Teff, Logg, Filter="V"):
    """ Looks up the non quadtractic limb darkening coefficients in the Claret table

    :param Z: stellar metallicity (Fe/H)
    :param Teff: Effective temperature of the star (Kelvin)
    :param Logg: log(g) of the star
    :param Filter: which filter to retreive the coefficents for out of u, v, b, y, U, B, V, R, I, J, H, K

    :return: The 4 non quadratic limb darkening coefficients
    :rtype: (a1, a2 ,a3, a4)

    :raises PyLC_FilterError: If invalid filter is given
    """
    filterlist = (('u', 'v', 'b', 'y', 'U', 'B', 'V', 'R', 'I', 'J', 'H', 'K'),
                  ( 4,   5,   6,   7,   8,   9,   10,  11,  12,  13,  14,  15)
                  )
    if Filter not in filterlist[0]:
        raise PyLCFilterError("Invalid filter, got {} must be in {}".format(Filter, filterlist[0]))

    # This could probably all be cleaned up by importing to a pandas dataframe
    Filter = filterlist[1][filterlist[0].index(Filter)]
    tables, mett = np.loadtxt(glob.glob(__location__ + '/*claretinfo*')[0], usecols=(0, 4), unpack=True)
    Table = str(int(tables[np.argmin(abs(Z - mett))]))
    File = glob.glob(__location__ + '/*/TABLE' + Table)[0]
    logglist, tefflist = np.loadtxt(File, usecols=(1, 2), unpack=True, skiprows=5)
    Teff = tefflist[np.argmin(abs(Teff - tefflist))]
    Logg = logglist[np.argmin(abs(Logg - logglist))]
    ld_coeffs = []
    for i in open(File).readlines()[5:]:
        coef = float(i.split()[Filter])
        logg = float(i.split()[1])
        teff = float(i.split()[2])
        if (logg == Logg and teff == Teff):
            ld_coeffs.append(coef)
    return tuple(ld_coeffs)


class PyLCError(BaseException):
    pass

class PyLCOptimiseError(PyLCError):
    pass

class PyLCFilterError(PyLCError):
    pass

