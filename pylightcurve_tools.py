__all__ = ['limb_darkening', 'position_vector', 'flux_drop']

import glob
import os

import numpy as np

from scipy.optimize import curve_fit

import matplotlib.pyplot as plt


__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
pi = np.pi


class PYLCError(BaseException):
    pass


class PYLCOptimiseError(PYLCError):
    pass
    

class PYLCFilterError(PYLCError):
    pass


def limb_darkening(metallicity, effective_temperature, logg, photometric_filter):

    filterlist = (('u', 'v', 'b', 'y', 'U', 'B', 'V', 'R', 'I', 'J', 'H', 'K'),
                  (4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15))

    if photometric_filter not in filterlist[0]:
        raise PYLCFilterError("Invalid filter, got {} must be in {}".format(photometric_filter, filterlist[0]))

    # This could probably all be cleaned up by importing to a pandas dataframe
    photometric_filter = filterlist[1][filterlist[0].index(photometric_filter)]
    tables, mett = np.loadtxt(glob.glob(__location__ + '/claret/claretinfo.txt')[0], usecols=(0, 4), unpack=True)
    table = str(int(tables[np.argmin(abs(metallicity - mett))]))
    table_file = glob.glob(__location__ + '/claret/claret_tables/TABLE' + table)[0]
    logglist, tefflist = np.loadtxt(table_file, usecols=(1, 2), unpack=True, skiprows=5)
    teff0 = tefflist[np.argmin(abs(effective_temperature - tefflist))]
    logg0 = logglist[np.argmin(abs(logg - logglist))]
    ld_coeffs = []
    for i in open(table_file).readlines()[5:]:
        coef = float(i.split()[photometric_filter])
        logg = float(i.split()[1])
        effective_temperature = float(i.split()[2])
        if logg == logg0 and effective_temperature == teff0:
            ld_coeffs.append(coef)
    return tuple(ld_coeffs)


def position_vector(period, sma_over_rs, eccentricity, inclination, periastron, mid_time, time_array, ww=0):
    if np.isnan(periastron):
        periastron = 0.
    inclination = inclination * pi / 180.0
    periastron = periastron * pi / 180.0
    ww = ww * pi / 180.0

    if eccentricity == 0 and ww == 0:
        vv = 2 * pi * (time_array - mid_time) / period
        bb = sma_over_rs * np.cos(vv)
        return [bb * np.sin(inclination), sma_over_rs * np.sin(vv), - bb * np.cos(inclination)]

    if periastron < pi / 2:
        aa = 1.0 * pi / 2 - periastron
    else:
        aa = 5.0 * pi / 2 - periastron
    bb = 2 * np.arctan(np.sqrt((1 - eccentricity) / (1 + eccentricity)) * np.tan(aa / 2))
    if bb < 0:
        bb += 2 * pi
    mid_time -= (period / 2.0 / pi) * (bb - eccentricity * np.sin(bb))
    m = (time_array - mid_time - np.int_((time_array - mid_time) / period) * period) * 2.0 * pi / period
    u0 = m
    stop = False
    u1 = 0
    for ii in xrange(10000):  # setting a limit of 1k iterations - arbitrary limit
        u1 = u0 - (u0 - eccentricity * np.sin(u0) - m) / (1 - eccentricity * np.cos(u0))
        stop = (np.abs(u1 - u0) < 10 ** (-7)).all()
        if stop:
            break
        else:
            u0 = u1
    if not stop:
        raise PYLCOptimiseError("Failed to find a solution in 10000 loops")
    vv = 2 * np.arctan(np.sqrt((1 + eccentricity) / (1 - eccentricity)) * np.tan(u1 / 2))
    #
    rr = sma_over_rs * (1 - (eccentricity ** 2)) / (np.ones_like(vv) + eccentricity * np.cos(vv))
    aa = np.cos(vv + periastron)
    bb = np.sin(vv + periastron)
    x = rr * bb * np.sin(inclination)
    y = rr * (-aa * np.cos(ww) + bb * np.sin(ww) * np.cos(inclination))
    z = rr * (-aa * np.sin(ww) - bb * np.cos(ww) * np.cos(inclination))
    return [x, y, z]


def plot_trajectory(p_vector):

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    xx, yy, zz = p_vector
    ax.plot(xx, yy, zz, color='k', label='trajectory')

    u = np.linspace(0, 2 * np.pi, 200)
    v = np.linspace(0, np.pi, 100)
    x = 1 * np.outer(np.cos(u), np.sin(v))
    y = 1 * np.outer(np.sin(u), np.sin(v))
    z = 1 * np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_surface(x, y, z, color='y', linewidth=0, antialiased=False)

    max_range = np.array([xx.max() - xx.min(), yy.max() - yy.min(), zz.max() - zz.min()]).max()
    xb = 0.5 * max_range * np.mgrid[-1:2:2, -1:2:2, -1:2:2][0].flatten() + 0.5 * (xx.max() + xx.min())
    yb = 0.5 * max_range * np.mgrid[-1:2:2, -1:2:2, -1:2:2][1].flatten() + 0.5 * (yy.max() + yy.min())
    zb = 0.5 * max_range * np.mgrid[-1:2:2, -1:2:2, -1:2:2][2].flatten() + 0.5 * (zz.max() + zz.min())
    for xb, yb, zb in zip(xb, yb, zb):
        ax.plot([xb], [yb], [zb], 'w')

    ax.plot([0, plt.xlim()[1]], [0, 0], [0, 0], color='r', label='line of sight')
    ax.legend()

    plt.show()

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


def flux_drop(limb_darkening_coefficients, rp_over_rs, z_over_rs):

    if len(z_over_rs) == 0:
        return np.array([])

    a1, a2, a3, a4 = limb_darkening_coefficients

    # cases
    zsq = z_over_rs * z_over_rs
    sum_z_rprs = z_over_rs + rp_over_rs
    dif_z_rprs = rp_over_rs - z_over_rs
    sqr_dif_z_rprs = zsq - rp_over_rs ** 2
    case0 = np.where((z_over_rs == 0) & (rp_over_rs <= 1))
    case1 = np.where((z_over_rs < rp_over_rs) & (sum_z_rprs <= 1))
    casea = np.where((z_over_rs < rp_over_rs) & (sum_z_rprs > 1) & (dif_z_rprs < 1))
    caseb = np.where((z_over_rs < rp_over_rs) & (sum_z_rprs > 1) & (dif_z_rprs > 1))
    case2 = np.where((z_over_rs == rp_over_rs) & (sum_z_rprs <= 1))
    casec = np.where((z_over_rs == rp_over_rs) & (sum_z_rprs > 1))
    case3 = np.where((z_over_rs > rp_over_rs) & (sum_z_rprs < 1))
    case4 = np.where((z_over_rs > rp_over_rs) & (sum_z_rprs == 1))
    case5 = np.where((z_over_rs > rp_over_rs) & (sum_z_rprs > 1) & (sqr_dif_z_rprs < 1))
    case6 = np.where((z_over_rs > rp_over_rs) & (sum_z_rprs > 1) & (sqr_dif_z_rprs == 1))
    case7 = np.where((z_over_rs > rp_over_rs) & (sum_z_rprs > 1) & (sqr_dif_z_rprs > 1) & (-1 < dif_z_rprs))
    plus_case = np.concatenate((case1[0], case2[0], case3[0], case4[0], case5[0], casea[0], casec[0]))
    minus_case = np.concatenate((case3[0], case4[0], case5[0], case6[0], case7[0]))
    star_case = np.concatenate((case5[0], case6[0], case7[0], casea[0], casec[0]))

    # cross points
    ph = np.arccos(np.clip((1.0 - rp_over_rs ** 2 + zsq) / (2.0 * z_over_rs), -1, 1))
    theta_1 = np.zeros(len(z_over_rs))
    ph_case = np.concatenate((case5[0], casea[0], casec[0]))
    theta_1[ph_case] = ph[ph_case]
    theta_2 = np.arcsin(np.minimum(rp_over_rs / z_over_rs, 1))
    theta_2[case1] = pi
    theta_2[case2] = pi / 2.0
    theta_2[casea] = pi
    theta_2[casec] = pi / 2.0
    theta_2[case7] = ph[case7]

    # flux_upper
    plusflux = np.zeros(len(z_over_rs))
    plusflux[plus_case] = integral_plus_core(a1, a2, a3, a4, rp_over_rs, z_over_rs[plus_case], theta_1[plus_case], theta_2[plus_case])
    if len(case0[0]) > 0:
        plusflux[case0] = integral_centred(a1, a2, a3, a4, rp_over_rs, 0.0, pi)
    if len(caseb[0]) > 0:
        plusflux[caseb] = integral_centred(a1, a2, a3, a4, 1, 0.0, pi)

    # flux_lower
    minsflux = np.zeros(len(z_over_rs))
    minsflux[minus_case] = integral_minus_core(a1, a2, a3, a4, rp_over_rs, z_over_rs[minus_case], 0.0, theta_2[minus_case])

    # flux_star
    starflux = np.zeros(len(z_over_rs))
    starflux[star_case] = integral_centred(a1, a2, a3, a4, 1, 0.0, ph[star_case])

    # flux_total
    total_flux = integral_centred(a1, a2, a3, a4, 1, 0.0, 2.0 * pi)

    return 1 - (2.0 / total_flux) * (plusflux + starflux - minsflux)


def distribution(data_i):

    def gauss(x, aa, x0, bb):
        return aa * np.exp(-(x - x0) ** 2 / (2 * bb ** 2))

    data = np.array(data_i)
    xstep = np.sqrt(np.median((data - np.median(data)) ** 2)) / 5.0
    xmin = min(data)
    xmax = max(data)
    x_size = int((xmax - xmin) / xstep) + 1

    distrx = xmin + np.arange(x_size) * xstep
    data = np.int_((data - xmin) / xstep)
    distr = np.bincount(data)
    distr = np.insert(distr, len(distr), np.zeros(x_size - len(distr)))

    pick = np.max(distr)
    mean = distrx[np.argmax(distr)]
    sigma = np.abs(distrx[np.argmin(np.abs(distr - pick / 2))] - mean)
    try:
        popt, pcov = curve_fit(gauss, distrx, distr, p0=[pick, mean, sigma])
    except RuntimeError:
        popt = [0, np.mean(data_i), np.std(data_i)]
    return popt


def plot_correlations(names, traces, results, errors):

    def od_distribution(data):

        data = np.array(data)
        xstep = np.sqrt(np.median((data - np.median(data)) ** 2)) / 5.0
        xmin = min(data)
        xmax = max(data)
        x_size = int((xmax - xmin) / xstep) + 1

        distrx = xmin + np.arange(x_size) * xstep
        data = np.int_((data - xmin) / xstep)
        distr = np.bincount(data)
        distr = np.insert(distr, len(distr), np.zeros(x_size - len(distr)))

        plt.step(distrx, distr, c='k')

    def td_distribution(datax, datay):

        datax = np.array(datax)
        median = np.median(datax)
        med = np.sqrt(np.median((datax - median) ** 2))
        xstep = med / 5.0
        xmin = min(datax)
        xmax = max(datax)
        x_size = int((xmax - xmin) / xstep) + 1
        datax = np.int_((datax - xmin) / xstep)
        datay = np.array(datay)
        median = np.median(datay)
        med = np.sqrt(np.median((datay - median) ** 2))
        ystep = med / 5.0
        ymin = min(datay)
        ymax = max(datay)
        y_size = int((ymax - ymin) / ystep) + 1
        datay = np.int_((datay - ymin) / ystep)

        yx_size = x_size * y_size
        yx = datay * x_size + datax

        yx = np.bincount(yx)
        yx = np.insert(yx, len(yx), np.zeros(yx_size - len(yx)))

        xx, yy = np.meshgrid(xmin + np.arange(x_size) * xstep, ymin + np.arange(y_size) * ystep)

        final = np.reshape(yx, (y_size, x_size))
        plt.imshow(np.where(final > 0, np.log(final), 0), extent=(np.min(xx), np.max(xx), np.min(yy), np.max(yy)),
                   cmap=plt.cm.Greys, origin='lower', aspect='auto')

    test = 0
    while test == 0:

        test = 1
        for var in range(len(names)):
            if len(traces[var]) == 1:
                del traces[var]
                del names[var]
                test = 0
                break

    all_var = len(traces)
    plt.figure(figsize=(2 * all_var, 2 * all_var))

    for var in range(len(names)):

        plt.subplot(all_var, all_var, all_var * var + var + 1)
        od_distribution(traces[var])

        plt.axvline(results[var], c='r')
        plt.axvline(results[var] - errors[var], c='r', ls='--', lw=0.5)
        plt.axvline(results[var] + errors[var], c='r', ls='--', lw=0.5)

        plt.xticks(plt.xticks()[0], np.ones_like(plt.yticks()[0]))
        plt.yticks(plt.yticks()[0], np.ones_like(plt.yticks()[0]))
        plt.tick_params(left='off', right='off', top='off', bottom='off', labelbottom='off', labelleft='off')
        plt.xlabel(r'${0}$'.format(names[var]) + '\n' +
                   r'${0:.{width}f}$'.format(results[var], width=abs(int(np.log10(errors[var]))) + 2) + '\n' +
                   r'$\pm{0:.{width}f}$'.format(errors[var], width=abs(int(np.log10(errors[var]))) + 2), fontsize=15)

        plt.xlim(results[var] - 6 * errors[var], results[var] + 6 * errors[var])

        for j in range(var + 1, all_var):

            plt.subplot(all_var, all_var, all_var * var + 1 + j)
            td_distribution(traces[j], traces[var])

            plt.axvline(results[j], c='r')
            plt.axvline(results[j] - errors[j], c='r', ls='--', lw=0.5)
            plt.axvline(results[j] + errors[j], c='r', ls='--', lw=0.5)

            plt.axhline(results[var], c='r')
            plt.axhline(results[var] - errors[var], c='r', ls='--', lw=0.5)
            plt.axhline(results[var] + errors[var], c='r', ls='--', lw=0.5)

            plt.yticks(plt.yticks()[0], np.arange(len(plt.yticks()[0])))
            plt.xticks(plt.xticks()[0], np.arange(len(plt.xticks()[0])))
            plt.tick_params(bottom='off', left='off', right='off', top='off', labelbottom='off',
                            labelleft='off', labelright='off', labeltop='off')

            plt.xlim(results[j] - 6 * errors[j], results[j] + 6 * errors[j])
            plt.ylim(results[var] - 6 * errors[var], results[var] + 6 * errors[var])

    plt.subplots_adjust(hspace=0, wspace=0)
    plt.savefig('traces_correlations.pdf', bbox_inches='tight', dpi=200)
    plt.close()


def plot_traces(names, traces, results, errors):

    test = 0
    while test == 0:

        test = 1
        for var in range(len(names)):
            if len(traces[var]) == 1:
                del traces[var]
                del names[var]
                del results[var]
                del errors[var]
                test = 0
                break

    all_var = len(traces)
    plt.figure(figsize=(7, 2 * all_var))

    for var in range(len(names)):

        plt.subplot(all_var, 1, var + 1)
        plt.plot(traces[var], 'k-', lw=0.5)
        plt.axhline(results[var], c='r')
        plt.axhline(results[var] - errors[var], ls='--', c='r', lw=0.5)
        plt.axhline(results[var] + errors[var], ls='--', c='r', lw=0.5)

        plt.yticks(plt.yticks()[0], np.ones_like(plt.yticks()[0]))
        plt.tick_params(left='off', right='off', labelleft='off')
        plt.ylabel(r'${0}$'.format(names[var]) + '\n' +
                   r'${0:.{width}f}$'.format(results[var], width=abs(int(np.log10(errors[var]))) + 2) + '\n' +
                   r'$\pm{0:.{width}f}$'.format(errors[var], width=abs(int(np.log10(errors[var]))) + 2), fontsize=15)

        if var != all_var - 1:
            plt.xticks(plt.xticks()[0], np.ones_like(plt.yticks()[0]))
            plt.tick_params(labelbottom='off')
        else:
            plt.xlabel(r'$\mathrm{iteration}$', fontsize=20)

        plt.ylim(results[var] - 8 * errors[var], results[var] + 8 * errors[var])

    plt.subplots_adjust(hspace=0, wspace=0)
    plt.savefig('traces_all.pdf', bbox_inches='tight', dpi=200)
    plt.close()


def plot_model(datax, datay, model):

    plt.subplot2grid((4, 1), (0, 0), rowspan=3)

    plt.plot(datax, datay, 'bo', mec='b')
    datax2 = np.arange(datax[0], datax[-1], (datax[1] - datax[0]) / 100)
    plt.plot(datax2, model(datax2), 'r-')

    di = plt.yticks()[0][1] - plt.yticks()[0][0]
    ticks = np.arange(plt.ylim()[0], plt.ylim()[1] + di / 2, di)
    plt.yticks(np.round(ticks, 4), np.round(ticks, 4))
    plt.ylim(plt.ylim()[0] - di / 2, plt.ylim()[1] + di / 2)
    plt.ylabel(r'$\mathrm{relative} \ \mathrm{flux}$', fontsize=20)
    plt.tick_params(labelbottom='off')

    plt.text(plt.xlim()[0] + 0.02 * (plt.xlim()[-1] - plt.xlim()[0]),
             plt.ylim()[0] + 0.02 * (plt.ylim()[-1] - plt.ylim()[0]),
             r'$\mathrm{rms}_\mathrm{res} = %.3e$' % np.std(datay - model(datax)))

    plt.subplot(4, 1, 4)
    plt.axhline(0, color='k')
    plt.plot(datax, datay - model(datax), 'ko', ms=3)

    di = plt.yticks()[0][1] - plt.yticks()[0][0]
    ticks = np.arange(plt.ylim()[0], plt.ylim()[1] + di / 2, di)
    plt.yticks(np.round(ticks, 4), np.round(ticks, 4))
    plt.ylim(plt.ylim()[0] - di / 2, plt.ylim()[1] + di / 2)

    plt.xlabel(r'$\mathrm{time} \ \mathrm{days}$', fontsize=20)
    plt.ylabel(r'$\mathrm{residuals}$', fontsize=20)

    plt.subplots_adjust(hspace=0.0)
    plt.savefig('model.pdf', bbox_inches='tight', dpi=200)


def save_results(names, results, errors):

    w = open('fitting_results.txt', 'w')

    for var in range(len(names)):
        w.write('{0}\t{1}\t{2}\n'.format(names[var], results[var], errors[var]))

    w.close()


def save_model(datax, datay, final_model):

    np.savetxt('model.txt', np.swapaxes([datax, datay, final_model(datax), datay - final_model(datax)], 0, 1))


def save_traces(names, traces):

    for var in range(len(names)):
        if len(traces[var]) > 1:
            np.savetxt('trace_' + names[var] + '.txt', traces[var])


