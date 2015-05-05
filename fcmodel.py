import glob
import os
import math

import numpy as np
import scipy.integrate as spi

pi = np.pi
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))


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
    return ld_coeffs


class PyLCError(BaseException):
    pass


class PyLCFilterError(PyLCError):
    pass


class PyLCOptimiseError(PyLCError):
    pass


def INTR(a1, a2, a3, a4, r):
    a0 = 1.0 - a1 - a2 - a3 - a4
    RR = (1.0 - r ** 2) ** (1.0 / 4)
    AA4 = - (2.0 / 4) * a0 * (RR ** 4.0)
    AA5 = - (2.0 / 5) * a1 * (RR ** 5.0)
    AA6 = - (2.0 / 6) * a2 * (RR ** 6.0)
    AA7 = - (2.0 / 7) * a3 * (RR ** 7.0)
    AA8 = - (2.0 / 8) * a4 * (RR ** 8.0)
    return AA4 + AA5 + AA6 + AA7 + AA8


# def num(r, a1, a2, a3, a4, p, z):
#     return (1.0 - a1 * (1.0 - (1.0 - r ** 2) ** (1.0 / 4)) - a2 * (1.0 - (1.0 - r ** 2) ** (1.0 / 2)) - a3 * (
#     1.0 - (1.0 - r ** 2) ** (3.0 / 4)) - a4 * (1.0 - (1.0 - r ** 2)) ) * r * np.arccos(
#         ( -p ** 2 + z ** 2 + r ** 2 ) / ( 2.0 * z * r ))
# 
# 
# def intr0(a1, a2, a3, a4, p, z, r1, r2):
#     return spi.quad(num, r1, r2, args=(a1, a2, a3, a4, p, z))[0]
# 
# 
# intr = np.vectorize(intr0)

def num(r,a1,a2,a3,a4,p,z):
	arccos	=	(( -p**2 + z**2 +r**2 )/( 2.0*z*r ))
	arccos	=	np.where(arccos>1,1,arccos)
	return ( 1.0 - a1*(1.0 - (1.0 - r**2)**(1.0/4)) - a2*(1.0 - (1.0 - r**2)**(1.0/2)) - a3*(1.0 - (1.0 - r**2)**(3.0/4)) - a4*(1.0 - (1.0 - r**2)) )*r*np.arccos( arccos )


def intr(a1,a2,a3,a4,p,z,r1,r2):
	bins	=	30.0
	h		=	(r2-r1)/bins
	test1	=	np.meshgrid(r1,range(int(bins)))
	test1	=	test1[0] + h*test1[1]
	test2	=	test1 + h
	z		=	np.meshgrid(z,range(int(bins)))[0]
	result	=	((test2-test1)/6.0)*(num(test1,a1,a2,a3,a4,p,z)+4.0*num((test1+test2)/2.0,a1,a2,a3,a4,p,z)+num(test2,a1,a2,a3,a4,p,z))
	return np.sum(result,0)


def INTCENT(a1, a2, a3, a4, p, z, ww1, ww2):
    w1 = np.minimum(ww1, ww2)
    w2 = np.maximum(ww1, ww2)
    return ( INTR(a1, a2, a3, a4, p) - INTR(a1, a2, a3, a4, 0.0) ) * ( w2 - w1 )


def INTPLUS(a1, a2, a3, a4, p, z, ww1, rr1, ww2, rr2):
    if len(z) == 0:
        return []
    w1 = np.minimum(ww1, ww2)
    r1 = np.minimum(rr1, rr2)
    w2 = np.maximum(ww1, ww2)
    r2 = np.maximum(rr1, rr2)
    PARTA = INTR(a1, a2, a3, a4, 0.0) * ( w1 - w2 )
    PARTB = INTR(a1, a2, a3, a4, r1) * (  w2 )
    PARTC = INTR(a1, a2, a3, a4, r2) * ( -w1 )
    PARTD = intr(a1, a2, a3, a4, p, z, r1, r2)
    return PARTA + PARTB + PARTC + PARTD


def INTMINS(a1, a2, a3, a4, p, z, ww1, rr1, ww2, rr2):
    if len(z) == 0:
        return []
    w1 = np.minimum(ww1, ww2)
    r1 = np.minimum(rr1, rr2)
    w2 = np.maximum(ww1, ww2)
    r2 = np.maximum(rr1, rr2)
    PARTA = INTR(a1, a2, a3, a4, 0.0) * ( w1 - w2 )
    PARTB = INTR(a1, a2, a3, a4, r1) * ( -w1 )
    PARTC = INTR(a1, a2, a3, a4, r2) * (  w2 )
    PARTD = intr(a1, a2, a3, a4, p, z, r1, r2)
    return PARTA + PARTB + PARTC - PARTD


def position(P, A, E, I, W, T0, tt, WW=0):
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
    for i in xrange(1000):  # setting a limit of 1k iterations - arbitrary limit
        u1 = u0 - (u0 - E * np.sin(u0) - M) / (1 - E * np.cos(u0))
        stop = (np.abs(u1 - u0) < 10 ** (-7)).all()
        if stop:
            break
        else:
            u0 = u1

            if math.isnan(stop):
                raise ValueError("Nan produced in loop, check inputs")

    if not stop:
        raise PyLCOptimiseError("Failed to find a solution in 1000 loops")

    vv = 2 * np.arctan(np.sqrt((1 + E) / (1 - E)) * np.tan(u1 / 2))
    #
    rr = A * (1 - (E ** 2)) / (np.ones_like(vv) + E * np.cos(vv))
    AA = np.cos(vv + W)
    BB = np.sin(vv + W)
    X = rr * BB * np.sin(I)
    Y = rr * (-AA * np.cos(WW) + BB * np.sin(WW) * np.cos(I))
    Z = rr * (-AA * np.sin(WW) - BB * np.cos(WW) * np.cos(I))
    return [X, Y, Z]


def model(ldcoeffs, RpRs, P, a, e, i, W, T0, tt, WW=0):
    """ Generates the lightcurve model

    :param RpRs: [dimensionless]
    :param P: Period [days]
    :param a: Semi-major axis/ Rstar [dimensionless]
    :param e: eccentricity [no units]
    :param i: inclination [degrees]
    :param W: argument of periastron [degrees]
    :param T0: epoch [JD] (units must match tt)
    :param tt: time array [JD] (units must match T0)
    :param WW: Omega [degrees]
    :return: transit depth for each element tt
    """

    a1, a2, a3, a4 = ldcoeffs

    if math.isnan(W):
        W = 0.

    p =RpRs
    ## projected distance
    pos = position(P, a, e, i * pi / 180, W * pi / 180, T0, tt, WW * pi / 180)
    fx = pos[0]
    fy = pos[1]
    fz = pos[2]
    z = np.sqrt(fy ** 2 + fz ** 2)
    ## cases
    case1 = np.where((fx > 0) & (z == 0))
    case2 = np.where((fx > 0) & (z < p))
    case3 = np.where((fx > 0) & (z == p))
    case4 = np.where((fx > 0) & (z > p) & (z < 1 - p))
    case5 = np.where((fx > 0) & (z == 1 - p))
    case6 = np.where((fx > 0) & (z > 1 - p) & (z ** 2 - p ** 2 < 1))
    case7 = np.where((fx > 0) & (z ** 2 - p ** 2 == 1))
    case8 = np.where((fx > 0) & (z ** 2 - p ** 2 > 1) & (z < 1 + p))
    ## cross points
    th = np.arcsin(np.where(p / z > 1.0, 1.0, p / z))
    ro = np.sqrt(abs(z ** 2 - p ** 2))
    ph = np.arccos(
        np.where((1.0 - p ** 2 + z ** 2) / (2.0 * z) > 1.0, 1.0, (1.0 - p ** 2 + z ** 2) / (2.0 * z)))
    ## flux
    plusflux = np.zeros(len(z))
    plusflux[case1] = INTCENT(a1, a2, a3, a4, p, z[case1], 0.0, 2 * pi)
    plusflux[case2] = INTPLUS(a1, a2, a3, a4, p, z[case2], 0.0, p + z[case2], pi, p - z[case2])
    plusflux[case3] = INTPLUS(a1, a2, a3, a4, p, z[case3], 0.0, 2 * p, pi / 2, 0.0)
    plusflux[case4] = INTPLUS(a1, a2, a3, a4, p, z[case4], 0.0, p + z[case4], th[case4], ro[case4])
    plusflux[case5] = INTPLUS(a1, a2, a3, a4, p, z[case5], 0.0, 1.0, th[case5], ro[case5])
    plusflux[case6] = INTPLUS(a1, a2, a3, a4, p, z[case6], ph[case6], 1.0, th[case6], ro[case6])
    minsflux = np.zeros(len(z))
    minscase = np.concatenate((case4[0],case5[0],case6[0],case7[0]))
    minsflux[minscase] = INTMINS(a1, a2, a3, a4, p, z[minscase], 0.0, z[minscase] - p, th[minscase], ro[minscase])
    minsflux[case8] = INTMINS(a1, a2, a3, a4, p, z[case8], 0.0, z[case8] - p, ph[case8], 1.0)
    starflux = np.zeros(len(z))
    starcase = np.concatenate((case6[0],case7[0],case8[0]))
    starflux[starcase] = INTCENT(a1, a2, a3, a4, 1, z[starcase], 0.0, ph[starcase])
    F0 = INTCENT(a1, a2, a3, a4, 1, 0, 0.0, 2.0 * pi)
    return np.array(1 - 2.0 * ( plusflux + starflux - minsflux ) / F0)