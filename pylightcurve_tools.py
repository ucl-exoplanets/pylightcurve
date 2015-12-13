import glob
import os

import numpy as np

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
pi = np.pi


class PYLCError(BaseException):
    pass


class PYLCOptimiseError(PYLCError):
    pass
    

class PYLCFilterError(PYLCError):
    pass


def ldcoeff(metall, teff, logg, phot_filter):

    filterlist = (('u', 'v', 'b', 'y', 'U', 'B', 'V', 'R', 'I', 'J', 'H', 'K'),
                  (4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15))

    if phot_filter not in filterlist[0]:
        raise PYLCFilterError("Invalid filter, got {} must be in {}".format(phot_filter, filterlist[0]))

    # This could probably all be cleaned up by importing to a pandas dataframe
    phot_filter = filterlist[1][filterlist[0].index(phot_filter)]
    tables, mett = np.loadtxt(glob.glob(__location__ + '/claret/claretinfo.txt')[0], usecols=(0, 4), unpack=True)
    table = str(int(tables[np.argmin(abs(metall - mett))]))
    table_file = glob.glob(__location__ + '/claret/claret_tables/TABLE' + table)[0]
    logglist, tefflist = np.loadtxt(table_file, usecols=(1, 2), unpack=True, skiprows=5)
    teff0 = tefflist[np.argmin(abs(teff - tefflist))]
    logg0 = logglist[np.argmin(abs(logg - logglist))]
    ld_coeffs = []
    for i in open(table_file).readlines()[5:]:
        coef = float(i.split()[phot_filter])
        logg = float(i.split()[1])
        teff = float(i.split()[2])
        if logg == logg0 and teff == teff0:
            ld_coeffs.append(coef)
    return tuple(ld_coeffs)


def position(p, a, e, i, w, ww, t0, tt):
    if w < pi / 2:
        aa = 1.0 * pi / 2 - w
    else:
        aa = 5.0 * pi / 2 - w
    bb = 2 * np.arctan(np.sqrt((1 - e) / (1 + e)) * np.tan(aa / 2))
    if bb < 0:
        bb += 2 * pi
    t0 -= (p / 2.0 / pi) * (bb - e * np.sin(bb))
    m = (tt - t0 - np.int_((tt - t0) / p) * p) * 2.0 * pi / p
    u0 = m
    stop = False
    u1 = 0
    for ii in xrange(10000):  # setting a limit of 1k iterations - arbitrary limit
        u1 = u0 - (u0 - e * np.sin(u0) - m) / (1 - e * np.cos(u0))
        stop = (np.abs(u1 - u0) < 10 ** (-7)).all()
        if stop:
            break
        else:
            u0 = u1
    if not stop:
        raise PYLCOptimiseError("Failed to find a solution in 10000 loops")
    vv = 2 * np.arctan(np.sqrt((1 + e) / (1 - e)) * np.tan(u1 / 2))
    #
    rr = a * (1 - (e ** 2)) / (np.ones_like(vv) + e * np.cos(vv))
    aa = np.cos(vv + w)
    bb = np.sin(vv + w)
    x = rr * bb * np.sin(i)
    y = rr * (-aa * np.cos(ww) + bb * np.sin(ww) * np.cos(i))
    z = rr * (-aa * np.sin(ww) - bb * np.cos(ww) * np.cos(i))
    return [x, y, z]


gauss30 = [
    [0.1028526528935588, -0.0514718425553177],
    [0.1028526528935588, 0.0514718425553177],
    [0.1017623897484055, -0.1538699136085835],
    [0.1017623897484055, 0.1538699136085835],
    [0.0995934205867953, -0.2546369261678899],
    [0.0995934205867953, 0.2546369261678899],
    [0.0963687371746443, -0.3527047255308781],
    [0.0963687371746443, 0.3527047255308781],
    [0.0921225222377861, -0.4470337695380892],
    [0.0921225222377861, 0.4470337695380892],
    [0.0868997872010830, -0.5366241481420199],
    [0.0868997872010830, 0.5366241481420199],
    [0.0807558952294202, -0.6205261829892429],
    [0.0807558952294202, 0.6205261829892429],
    [0.0737559747377052, -0.6978504947933158],
    [0.0737559747377052, 0.6978504947933158],
    [0.0659742298821805, -0.7677774321048262],
    [0.0659742298821805, 0.7677774321048262],
    [0.0574931562176191, -0.8295657623827684],
    [0.0574931562176191, 0.8295657623827684],
    [0.0484026728305941, -0.8825605357920527],
    [0.0484026728305941, 0.8825605357920527],
    [0.0387991925696271, -0.9262000474292743],
    [0.0387991925696271, 0.9262000474292743],
    [0.0287847078833234, -0.9600218649683075],
    [0.0287847078833234, 0.9600218649683075],
    [0.0184664683110910, -0.9836681232797472],
    [0.0184664683110910, 0.9836681232797472],
    [0.0079681924961666, -0.9968934840746495],
    [0.0079681924961666, 0.9968934840746495]
]
gausstab = np.swapaxes(gauss30, 0, 1)


def integral_r(a1, a2, a3, a4, r):
    mu44 = 1.0 - r * r
    mu24 = np.sqrt(mu44)
    mu14 = np.sqrt(mu24)
    return - (2.0 * (1.0 - a1 - a2 - a3 - a4) / 4) * mu44 \
           - (2.0 * a1 / 5) * mu44 * mu14 \
           - (2.0 * a2 / 6) * mu44 * mu24 \
           - (2.0 * a3 / 7) * mu44 * mu24 * mu14 \
           - (2.0 * a4 / 8) * mu44 * mu44


def num(r, a1, a2, a3, a4, rprs, z):
    rsq = r * r
    mu44 = 1.0 - rsq
    mu24 = np.sqrt(mu44)
    mu14 = np.sqrt(mu24)
    return ((1.0 - a1 - a2 - a3 - a4) + a1 * mu14 + a2 * mu24 + a3 * mu24 * mu14 + a4 * mu44) \
        * r * np.arccos(np.minimum((-rprs ** 2 + z * z + rsq) / (2.0 * z * r), 1.0))


def integral_r_f(a1, a2, a3, a4, rprs, z, r1, r2):
    x1 = (r2 - r1) / 2
    x2 = (r2 + r1) / 2
    return x1 * np.sum(gausstab[0][:, None] * num(x1[None, :] * gausstab[1][:, None] + x2[None, :],
                                                  a1, a2, a3, a4, rprs, z), 0)


def integral_centred(a1, a2, a3, a4, rprs, ww1, ww2):
    return (integral_r(a1, a2, a3, a4, rprs) - integral_r(a1, a2, a3, a4, 0.0)) * np.abs(ww2 - ww1)


def integral_plus_core(a1, a2, a3, a4, rprs, z, ww1, ww2):
    if len(z) == 0:
        return z
    rr1 = z * np.cos(ww1) + np.sqrt(np.maximum(rprs ** 2 - (z * np.sin(ww1)) ** 2, 0))
    rr1 = np.clip(rr1, 0, 1)
    rr2 = z * np.cos(ww2) + np.sqrt(np.maximum(rprs ** 2 - (z * np.sin(ww2)) ** 2, 0))
    rr2 = np.clip(rr2, 0, 1)
    w1 = np.minimum(ww1, ww2)
    r1 = np.minimum(rr1, rr2)
    w2 = np.maximum(ww1, ww2)
    r2 = np.maximum(rr1, rr2)
    parta = integral_r(a1, a2, a3, a4, 0.0) * (w1 - w2)
    partb = integral_r(a1, a2, a3, a4, r1) * w2
    partc = integral_r(a1, a2, a3, a4, r2) * (-w1)
    partd = integral_r_f(a1, a2, a3, a4, rprs, z, r1, r2)
    return parta + partb + partc + partd


def integral_plus(a1, a2, a3, a4, rprs, z, ww1, ww2):
    split = np.where(ww1 * ww2 < 0)
    if len(split[0]) == 0:
        return integral_plus_core(a1, a2, a3, a4, rprs, z, np.abs(ww1), np.abs(ww2))
    else:
        w1 = np.minimum(ww1, ww2)
        w2 = np.maximum(ww1, ww2)
        intplus = integral_plus_core(a1, a2, a3, a4, rprs, z, np.where(w1 * w2 < 0, 0, np.abs(w1)), np.abs(w2))
        intplus[split] = intplus[split] + integral_plus_core(a1, a2, a3, a4, rprs, z[split], 0, np.abs(w1[split]))
        return intplus    


def integral_minus_core(a1, a2, a3, a4, rprs, z, ww1, ww2):
    if len(z) == 0:
        return z
    rr1 = z * np.cos(ww1) - np.sqrt(np.maximum(rprs ** 2 - (z * np.sin(ww1)) ** 2, 0))
    rr1 = np.clip(rr1, 0, 1)
    rr2 = z * np.cos(ww2) - np.sqrt(np.maximum(rprs ** 2 - (z * np.sin(ww2)) ** 2, 0))
    rr2 = np.clip(rr2, 0, 1)
    w1 = np.minimum(ww1, ww2)
    r1 = np.minimum(rr1, rr2)
    w2 = np.maximum(ww1, ww2)
    r2 = np.maximum(rr1, rr2)
    parta = integral_r(a1, a2, a3, a4, 0.0) * (w1 - w2)
    partb = integral_r(a1, a2, a3, a4, r1) * (-w1)
    partc = integral_r(a1, a2, a3, a4, r2) * w2
    partd = integral_r_f(a1, a2, a3, a4, rprs, z, r1, r2)
    return parta + partb + partc - partd


def integral_minus(a1, a2, a3, a4, rprs, z, ww1, ww2):
    split = np.where(ww1 * ww2 < 0)
    if len(split[0]) == 0:
        return integral_minus_core(a1, a2, a3, a4, rprs, z, np.abs(ww1), np.abs(ww2))
    else:
        w1 = np.minimum(ww1, ww2)
        w2 = np.maximum(ww1, ww2)
        intmins = integral_minus_core(a1, a2, a3, a4, rprs, z, np.where(w1 * w2 < 0, 0, np.abs(w1)), np.abs(w2))
        intmins[split] = intmins[split] + integral_minus_core(a1, a2, a3, a4, rprs, z[split], 0, np.abs(w1[split]))
        return intmins    


def single_model(ldcoeffs, rprs, xyz, tt):

    if len(tt) == 0:
        return np.array([])

    a1, a2, a3, a4 = ldcoeffs

    # projected distance
    fx, fy, fz = xyz
    z = np.sqrt(fy * fy + fz * fz)
    zsq = z * z

    # cases
    sum_z_rprs = z + rprs
    dif_z_rprs = rprs - z
    sqr_dif_z_rprs = zsq - rprs ** 2
    case0 = np.where((fx > 0) & (z == 0) & (rprs <= 1))
    case1 = np.where((fx > 0) & (z < rprs) & (sum_z_rprs <= 1))
    casea = np.where((fx > 0) & (z < rprs) & (sum_z_rprs > 1) & (dif_z_rprs < 1))
    caseb = np.where((fx > 0) & (z < rprs) & (sum_z_rprs > 1) & (dif_z_rprs > 1))
    case2 = np.where((fx > 0) & (z == rprs) & (sum_z_rprs <= 1))
    casec = np.where((fx > 0) & (z == rprs) & (sum_z_rprs > 1))
    case3 = np.where((fx > 0) & (z > rprs) & (sum_z_rprs < 1))
    case4 = np.where((fx > 0) & (z > rprs) & (sum_z_rprs == 1))
    case5 = np.where((fx > 0) & (z > rprs) & (sum_z_rprs > 1) & (sqr_dif_z_rprs < 1))
    case6 = np.where((fx > 0) & (z > rprs) & (sum_z_rprs > 1) & (sqr_dif_z_rprs == 1))
    case7 = np.where((fx > 0) & (z > rprs) & (sum_z_rprs > 1) & (sqr_dif_z_rprs > 1) & (-1 < dif_z_rprs))
    plus_case = np.concatenate((case1[0], case2[0], case3[0], case4[0], case5[0], casea[0], casec[0]))
    minus_case = np.concatenate((case3[0], case4[0], case5[0], case6[0], case7[0]))
    star_case = np.concatenate((case5[0], case6[0], case7[0], casea[0], casec[0]))

    # cross points
    th = np.arcsin(np.minimum(rprs / z, 1))
    ph = np.arccos(np.clip((1.0 - rprs ** 2 + zsq) / (2.0 * z), -1, 1))
    theta_1 = np.zeros(len(z))
    theta_1[case5] = ph[case5]
    theta_1[casea] = ph[casea]
    theta_1[casec] = ph[casec]
    theta_2 = np.full_like(th, th)
    theta_2[case1] = pi
    theta_2[case2] = pi / 2.0
    theta_2[casea] = pi
    theta_2[casec] = pi / 2.0
    theta_2[case7] = ph[case7]

    # flux_upper
    plusflux = np.zeros(len(z))
    plusflux[plus_case] = integral_plus_core(a1, a2, a3, a4, rprs, z[plus_case], theta_1[plus_case], theta_2[plus_case])
    if len(case0[0]) > 0:
        plusflux[case0] = integral_centred(a1, a2, a3, a4, rprs, 0.0, pi)
    if len(caseb[0]) > 0:
        plusflux[caseb] = integral_centred(a1, a2, a3, a4, 1, 0.0, pi)

    # flux_lower
    minsflux = np.zeros(len(z))
    minsflux[minus_case] = integral_minus_core(a1, a2, a3, a4, rprs, z[minus_case], 0.0, theta_2[minus_case])

    # flux_star
    starflux = np.zeros(len(z))
    starflux[star_case] = integral_centred(a1, a2, a3, a4, 1, 0.0, ph[star_case])

    # flux_total
    total_flux = integral_centred(a1, a2, a3, a4, 1, 0.0, 2.0 * pi)

    return 1 - (2.0 / total_flux) * (plusflux + starflux - minsflux)
