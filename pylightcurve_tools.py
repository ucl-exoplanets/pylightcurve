__all__ = ['limb_darkening', 'jd_to_hjd']

import glob
import os

import numpy as np

from scipy.optimize import curve_fit

import matplotlib.pyplot as plt

import ephem


__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
pi = np.pi


class PYLCError(BaseException):
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


def jd_to_hjd(jd, ra, dec):

    ra_target, dec_target = ra * np.pi / 180, dec * np.pi / 180

    try:
        hjd = []

        for julian_date in jd:
            # calculate the position of the sun on the sky for this date

            sun = ephem.Sun()
            sun.compute(ephem.date(julian_date - 2415020))
            ra_sun, dec_sun = float(sun.ra), float(sun.dec)

            # calculate the correction coefficients

            a = 149597870700.0 / ephem.c
            b = np.sin(dec_target) * np.sin(dec_sun)
            c = np.cos(dec_target) * np.cos(dec_sun) * np.cos(ra_target - ra_sun)

            # apply the correction and save the result as the heliocentric_julian_date keyword

            heliocentric_julian_date = julian_date - (a * (b + c)) / (24.0 * 60.0 * 60.0)

            hjd.append(heliocentric_julian_date)

        return np.array(hjd)

    except TypeError:

        julian_date = jd

        # calculate the position of the sun on the sky for this date

        sun = ephem.Sun()
        sun.compute(ephem.date(julian_date - 2415020))
        ra_sun, dec_sun = float(sun.ra), float(sun.dec)

        # calculate the correction coefficients

        a = 149597870700.0 / ephem.c
        b = np.sin(dec_target) * np.sin(dec_sun)
        c = np.cos(dec_target) * np.cos(dec_sun) * np.cos(ra_target - ra_sun)

        # apply the correction and save the result as the heliocentric_julian_date keyword

        heliocentric_julian_date = julian_date - (a * (b + c)) / (24.0 * 60.0 * 60.0)

        return heliocentric_julian_date


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
        digit = abs(int(np.log10(errors[var]))) + 1
        done = 1 / int(errors[var] * (10 ** digit))
        plt.xlabel(r'${0}$'.format(names[var]) + '\n' +
                   r'${0:.{width}f}$'.format(results[var], width=digit + done) + '\n' +
                   r'$\pm{0:.{width}f}$'.format(errors[var], width=digit + done), fontsize=15)

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
    plt.close('all')


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
        digit = abs(int(np.log10(errors[var]))) + 1
        done = 1 / int(errors[var] * (10 ** digit))
        plt.ylabel(r'${0}$'.format(names[var]) + '\n' +
                   r'${0:.{width}f}$'.format(results[var], width=digit + done) + '\n' +
                   r'$\pm{0:.{width}f}$'.format(errors[var], width=digit + done), fontsize=15)

        if var != all_var - 1:
            plt.xticks(plt.xticks()[0], np.ones_like(plt.xticks()[0]))
            plt.tick_params(labelbottom='off')
        else:
            plt.xlabel(r'$\mathrm{iteration}$', fontsize=20)

        plt.ylim(results[var] - 8 * errors[var], results[var] + 8 * errors[var])

    plt.subplots_adjust(hspace=0, wspace=0)
    plt.savefig('traces_all.pdf', bbox_inches='tight', dpi=200)
    plt.close('all')


def transit_duration(period, rp_rs, semimajoraxis, inclination, eccentricity, periastron):
    ww = periastron * np.pi / 180
    ii = inclination * np.pi / 180
    ee = eccentricity
    aa = semimajoraxis
    ro_pt = (1 - ee ** 2) / (1 + ee * np.sin(ww))
    b_pt = aa * ro_pt * np.cos(ii)
    if b_pt > 1:
        b_pt = 0.5
    s_ps = 1.0 + rp_rs
    df = np.arcsin(np.sqrt((s_ps ** 2 - b_pt ** 2) / ((aa ** 2) * (ro_pt ** 2) - b_pt ** 2)))
    abs_value = (period * (ro_pt ** 2)) / (np.pi * np.sqrt(1 - ee ** 2)) * df
    return abs_value


def plot_model(datax, datay, names, initial, results, errors, final_model, final_systematics_model, set_number):

    rp = results[names.index('rp')]
    rp_e = errors[names.index('rp')]
    mt = results[names.index('mt')] + \
        round((np.mean(datax) - results[names.index('mt')]) / results[names.index('P')]) * results[names.index('P')]
    mt_e = errors[names.index('mt')]
    period = results[names.index('P')]

    prediction = initial[names.index('mt')] + \
        round((np.mean(datax) - initial[names.index('mt')]) / initial[names.index('P')]) * initial[names.index('P')]

    duration = transit_duration(initial[names.index('P')], initial[names.index('rp')], initial[names.index('a')],
                                initial[names.index('i')], initial[names.index('e')], initial[names.index('w')])

    ingress = prediction - duration / 2
    egress = prediction + duration / 2

    # plot1

    plt.subplot2grid((4, 1), (0, 0), rowspan=3)

    plt.plot((datax - mt) / period, datay, 'ko', ms=3)
    datax2 = np.arange(datax[0], datax[-1], (datax[1] - datax[0]) / 100)
    plt.plot((datax2 - mt) / period, final_systematics_model(datax2) * final_model(datax2), 'r-')

    plt.yticks(plt.yticks()[0][1:])
    plt.ylabel(r'$\mathrm{relative} \ \mathrm{flux}$', fontsize=20)

    plt.ylim(min(datay), max(datay) + np.std(datay - final_systematics_model(datax) * final_model(datax)))
    plt.ylim(plt.ylim()[0] - 0.2 * (plt.ylim()[1] - plt.ylim()[0]), plt.ylim()[1])

    x_max = max(np.abs((datax - mt) / period)) + 0.05 * (max((datax - mt) / period) - min((datax - mt) / period))
    plt.xlim(-x_max, x_max)
    plt.tick_params(labelbottom='off')

    digit = abs(int(np.log10(rp_e))) + 1
    done = 1 / int(rp_e * (10 ** digit))
    rpstr = r'$R_\mathrm{p}/R_* = $' + r'${0:.{width}f}$'.format(rp, width=digit + done) + \
            r'$ \, \pm \, {0:.{width}f}$'.format(rp_e, width=digit + done)
    digit = abs(int(np.log10(mt_e))) + 1
    done = 1 / int(mt_e * (10 ** digit))
    mtstr = r'$T_\mathrm{HJD} = $' + r'${0:.{width}f}$'.format(mt, width=digit + done) + \
            r'$ \, \pm \, {0:.{width}f} \,$'.format(mt_e, width=digit + done)

    plt.text(plt.xlim()[0] + 0.5 * (plt.xlim()[-1] - plt.xlim()[0]),
             plt.ylim()[0] + 0.07 * (plt.ylim()[-1] - plt.ylim()[0]),
             rpstr + '\n' + mtstr, ha='center', va='center', fontsize=10)

    plt.axvline((ingress - mt) / period, 0.5, 1.0, ls='--', c='k')
    plt.text((ingress - mt) / period, 0.5 * (plt.ylim()[1] + plt.ylim()[0]),
             r'$\mathrm{predicted}$' + '\n' + r'$\mathrm{ingress start}$', ha='right', va='top', fontsize=10)
    plt.axvline((egress - mt) / period, 0.5, 1.0, ls='--', c='k')
    plt.text((egress - mt) / period, 0.5 * (plt.ylim()[1] + plt.ylim()[0]),
             r'$\mathrm{predicted}$' + '\n' + r'$\mathrm{egress end}$', ha='left', va='top', fontsize=10)

    plt.subplot(4, 1, 4)
    plt.cla()
    plt.plot((datax - mt) / period, np.zeros_like(datax), 'r-')
    plt.axhline(3 * np.std(datay / final_systematics_model(datax) - final_model(datax)), lw=0)
    plt.axhline(- 3 * np.std(datay / final_systematics_model(datax) - final_model(datax)), lw=0)
    plt.plot((datax - mt) / period, datay / final_systematics_model(datax) - final_model(datax), 'ko', ms=3)

    plt.xlabel(r'$\mathrm{phase}$', fontsize=20)
    plt.ylabel(r'$\mathrm{residuals}$', fontsize=20)

    plt.xlim(-x_max, x_max)

    plt.text(plt.xlim()[0] + 0.02 * (plt.xlim()[-1] - plt.xlim()[0]),
             plt.ylim()[0] + 0.07 * (plt.ylim()[-1] - plt.ylim()[0]),
             r'$\mathrm{rms}_\mathrm{res} = %.1e$'
             % np.std(datay - final_systematics_model(datax) * final_model(datax)), fontsize=10)

    plt.suptitle('dataset' + str(set_number), fontsize=20)
    plt.subplots_adjust(hspace=0.0)
    plt.savefig('model_dataset' + str(set_number) + '.pdf', bbox_inches='tight', dpi=200)
    plt.close('all')

    # plot2

    plt.subplot2grid((4, 1), (0, 0), rowspan=3)

    plt.plot((datax - mt) / period, datay / final_systematics_model(datax), 'ko', ms=3)
    datax2 = np.arange(datax[0], datax[-1], (datax[1] - datax[0]) / 100)
    plt.plot((datax2 - mt) / period, final_model(datax2), 'r-')

    plt.yticks(plt.yticks()[0][1:])
    plt.ylabel(r'$\mathrm{relative} \ \mathrm{flux}$', fontsize=20)

    plt.ylim(min(datay / final_systematics_model(datax)), max(datay / final_systematics_model(datax))
             + np.std(datay / final_systematics_model(datax) - final_model(datax)))
    plt.ylim(plt.ylim()[0] - 0.2 * (plt.ylim()[1] - plt.ylim()[0]), plt.ylim()[1])

    plt.xlim(-x_max, x_max)
    plt.tick_params(labelbottom='off')

    digit = abs(int(np.log10(rp_e))) + 1
    done = 1 / int(rp_e * (10 ** digit))
    rpstr = r'$R_\mathrm{p}/R_* = $' + r'${0:.{width}f}$'.format(rp, width=digit + done) + \
            r'$ \, \pm \, {0:.{width}f}$'.format(rp_e, width=digit + done)
    digit = abs(int(np.log10(mt_e))) + 1
    done = 1 / int(mt_e * (10 ** digit))
    mtstr = r'$T_\mathrm{HJD} = $' + r'${0:.{width}f}$'.format(mt, width=digit + done) + \
            r'$ \, \pm \, {0:.{width}f} \,$'.format(mt_e, width=digit + done)

    plt.text(plt.xlim()[0] + 0.5 * (plt.xlim()[1] - plt.xlim()[0]),
             plt.ylim()[0] + 0.07 * (plt.ylim()[1] - plt.ylim()[0]),
             rpstr + '\n' + mtstr, ha='center', va='center', fontsize=10)

    text_y = 2.0 - plt.ylim()[1]
    plt.axvline((ingress - mt) / period,
                (text_y - plt.ylim()[0]) / (plt.ylim()[1] - plt.ylim()[0]), 1.0, ls='--', c='k')
    plt.text((ingress - mt) / period, text_y,
             r'$\mathrm{predicted}$' + '\n' + r'$\mathrm{ingress}$', ha='right', va='top', fontsize=10)
    plt.axvline((egress - mt) / period,
                (text_y - plt.ylim()[0]) / (plt.ylim()[1] - plt.ylim()[0]), 1.0, ls='--', c='k')
    plt.text((egress - mt) / period, text_y,
             r'$\mathrm{predicted}$' + '\n' + r'$\mathrm{egress}$', ha='left', va='top', fontsize=10)

    plt.subplot(4, 1, 4)
    plt.cla()

    plt.plot((datax - mt) / period, np.zeros_like(datax), 'r-')
    plt.axhline(3 * np.std(datay / final_systematics_model(datax) - final_model(datax)), lw=0)
    plt.axhline(- 3 * np.std(datay / final_systematics_model(datax) - final_model(datax)), lw=0)
    plt.plot((datax - mt) / period, datay / final_systematics_model(datax) - final_model(datax), 'ko', ms=3)

    plt.xlabel(r'$\mathrm{phase}$', fontsize=20)
    plt.ylabel(r'$\mathrm{residuals}$', fontsize=20)

    plt.xlim(-x_max, x_max)

    plt.text(plt.xlim()[0] + 0.02 * (plt.xlim()[-1] - plt.xlim()[0]),
             plt.ylim()[0] + 0.07 * (plt.ylim()[-1] - plt.ylim()[0]),
             r'$\mathrm{rms}_\mathrm{res} = %.1e$'
             % np.std(datay / final_systematics_model(datax) - final_model(datax)), fontsize=10)

    plt.suptitle('dataset' + str(set_number), fontsize=20)
    plt.subplots_adjust(hspace=0.0)
    plt.savefig('model2_dataset' + str(set_number) + '.pdf', bbox_inches='tight', dpi=200)
    plt.close('all')


def save_results(names, initial, results, errors, limb_darkening_coefficients):

    w = open('fitting_results.txt', 'w')

    w.write('# variable\tinitial\tfinal\tuncertainty\n')

    for var in range(len(names)):
        w.write('{0}\t{1}\t{2}\t{3}\n'.format(names[var], initial[var], results[var], errors[var]))

    for i in range(4):
        w.write('# limb_darkening_coefficient_{0}\t{1}\n'.format(i + 1, limb_darkening_coefficients[i]))

    w.close()


def save_model(datax, datay, final_model, final_systematics_model, set_number):

    np.savetxt('model_dataset' + str(set_number) + '.txt',
               np.swapaxes([datax, datay, final_systematics_model(datax),
                            datay / final_systematics_model(datax), final_model(datax),
                            datay / final_systematics_model(datax) - final_model(datax)], 0, 1))


def save_traces(names, traces):

    for var in range(len(names)):
        if len(traces[var]) > 1:
            np.savetxt('trace_' + names[var] + '.txt', traces[var])
