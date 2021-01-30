
__all__ = ['planet_orbit', 'planet_star_projected_distance', 'planet_phase',
           'transit', 'transit_integrated', 'transit_depth', 'transit_duration',
           'eclipse', 'eclipse_integrated', 'eclipse_depth', 'eclipse_duration', 'eclipse_mid_time',
           'fp_over_fs', 'exotethys']


import os
import numpy as np

from pylightcurve.errors import *
from pylightcurve.__databases__ import plc_data
from pylightcurve.analysis.numerical_integration import gauss_numerical_integration
from pylightcurve.analysis.curve_fit import curve_fit
from pylightcurve.processes.files import open_dict


# orbit


def planet_orbit(period, sma_over_rs, eccentricity, inclination, periastron, mid_time, time_array, ww=0):

    inclination = inclination * np.pi / 180.0
    periastron = periastron * np.pi / 180.0
    ww = ww * np.pi / 180.0

    if eccentricity == 0 and ww == 0:
        vv = 2 * np.pi * (time_array - mid_time) / period
        bb = sma_over_rs * np.cos(vv)
        return [bb * np.sin(inclination), sma_over_rs * np.sin(vv), - bb * np.cos(inclination)]

    if periastron < np.pi / 2:
        aa = 1.0 * np.pi / 2 - periastron
    else:
        aa = 5.0 * np.pi / 2 - periastron
    bb = 2 * np.arctan(np.sqrt((1 - eccentricity) / (1 + eccentricity)) * np.tan(aa / 2))
    if bb < 0:
        bb += 2 * np.pi
    mid_time = float(mid_time) - (period / 2.0 / np.pi) * (bb - eccentricity * np.sin(bb))
    m = (time_array - mid_time - np.int_((time_array - mid_time) / period) * period) * 2.0 * np.pi / period
    u0 = m
    stop = False
    u1 = 0
    for ii in range(10000):  # setting a limit of 1k iterations - arbitrary limit
        u1 = u0 - (u0 - eccentricity * np.sin(u0) - m) / (1 - eccentricity * np.cos(u0))
        stop = (np.abs(u1 - u0) < 10 ** (-7)).all()
        if stop:
            break
        else:
            u0 = u1
    if not stop:
        raise RuntimeError('Failed to find a solution in 10000 loops')

    vv = 2 * np.arctan(np.sqrt((1 + eccentricity) / (1 - eccentricity)) * np.tan(u1 / 2))
    #
    rr = sma_over_rs * (1 - (eccentricity ** 2)) / (np.ones_like(vv) + eccentricity * np.cos(vv))
    aa = np.cos(vv + periastron)
    bb = np.sin(vv + periastron)
    x = rr * bb * np.sin(inclination)
    y = rr * (-aa * np.cos(ww) + bb * np.sin(ww) * np.cos(inclination))
    z = rr * (-aa * np.sin(ww) - bb * np.cos(ww) * np.cos(inclination))

    return [x, y, z]


def planet_star_projected_distance(period, sma_over_rs, eccentricity, inclination, periastron, mid_time, time_array):

    position_vector = planet_orbit(period, sma_over_rs, eccentricity, inclination, periastron, mid_time, time_array)

    return np.sqrt(position_vector[1] * position_vector[1] + position_vector[2] * position_vector[2])


def planet_phase(period, mid_time, time_array):
    return (time_array - mid_time)/period


# flux drop


def integral_r_claret(limb_darkening_coefficients, r):
    a1, a2, a3, a4 = limb_darkening_coefficients
    mu44 = 1.0 - r * r
    mu24 = np.sqrt(mu44)
    mu14 = np.sqrt(mu24)
    return - (2.0 * (1.0 - a1 - a2 - a3 - a4) / 4) * mu44 \
           - (2.0 * a1 / 5) * mu44 * mu14 \
           - (2.0 * a2 / 6) * mu44 * mu24 \
           - (2.0 * a3 / 7) * mu44 * mu24 * mu14 \
           - (2.0 * a4 / 8) * mu44 * mu44


def num_claret(r, limb_darkening_coefficients, rprs, z):
    a1, a2, a3, a4 = limb_darkening_coefficients
    rsq = r * r
    mu44 = 1.0 - rsq
    mu24 = np.sqrt(mu44)
    mu14 = np.sqrt(mu24)
    return ((1.0 - a1 - a2 - a3 - a4) + a1 * mu14 + a2 * mu24 + a3 * mu24 * mu14 + a4 * mu44) \
        * r * np.arccos(np.minimum((-rprs ** 2 + z * z + rsq) / (2.0 * z * r), 1.0))


def integral_r_f_claret(limb_darkening_coefficients, rprs, z, r1, r2, precision=3):
    return gauss_numerical_integration(num_claret, r1, r2, precision, limb_darkening_coefficients, rprs, z)


# integral definitions for zero method


def integral_r_zero(limb_darkening_coefficients, r):
    musq = 1 - r * r
    return (-1.0 / 6) * musq * 3.0


def num_zero(r, limb_darkening_coefficients, rprs, z):
    rsq = r * r
    return r * np.arccos(np.minimum((-rprs ** 2 + z * z + rsq) / (2.0 * z * r), 1.0))


def integral_r_f_zero(limb_darkening_coefficients, rprs, z, r1, r2, precision=3):
    return gauss_numerical_integration(num_zero, r1, r2, precision, limb_darkening_coefficients, rprs, z)


# integral definitions for linear method


def integral_r_linear(limb_darkening_coefficients, r):
    a1 = limb_darkening_coefficients[0]
    musq = 1 - r * r
    return (-1.0 / 6) * musq * (3.0 + a1 * (-3.0 + 2.0 * np.sqrt(musq)))


def num_linear(r, limb_darkening_coefficients, rprs, z):
    a1 = limb_darkening_coefficients[0]
    rsq = r * r
    return (1.0 - a1 * (1.0 - np.sqrt(1.0 - rsq))) \
        * r * np.arccos(np.minimum((-rprs ** 2 + z * z + rsq) / (2.0 * z * r), 1.0))


def integral_r_f_linear(limb_darkening_coefficients, rprs, z, r1, r2, precision=3):
    return gauss_numerical_integration(num_linear, r1, r2, precision, limb_darkening_coefficients, rprs, z)


# integral definitions for quadratic method


def integral_r_quad(limb_darkening_coefficients, r):
    a1, a2 = limb_darkening_coefficients[:2]
    musq = 1 - r * r
    mu = np.sqrt(musq)
    return (1.0 / 12) * (-4.0 * (a1 + 2.0 * a2) * mu * musq + 6.0 * (-1 + a1 + a2) * musq + 3.0 * a2 * musq * musq)


def num_quad(r, limb_darkening_coefficients, rprs, z):
    a1, a2 = limb_darkening_coefficients[:2]
    rsq = r * r
    cc = 1.0 - np.sqrt(1.0 - rsq)
    return (1.0 - a1 * cc - a2 * cc * cc) \
        * r * np.arccos(np.minimum((-rprs ** 2 + z * z + rsq) / (2.0 * z * r), 1.0))


def integral_r_f_quad(limb_darkening_coefficients, rprs, z, r1, r2, precision=3):
    return gauss_numerical_integration(num_quad, r1, r2, precision, limb_darkening_coefficients, rprs, z)


# integral definitions for square root method


def integral_r_sqrt(limb_darkening_coefficients, r):
    a1, a2 = limb_darkening_coefficients[:2]
    musq = 1 - r * r
    mu = np.sqrt(musq)
    return ((-2.0 / 5) * a2 * np.sqrt(mu) - (1.0 / 3) * a1 * mu + (1.0 / 2) * (-1 + a1 + a2)) * musq


def num_sqrt(r, limb_darkening_coefficients, rprs, z):
    a1, a2 = limb_darkening_coefficients[:2]
    rsq = r * r
    mu = np.sqrt(1.0 - rsq)
    return (1.0 - a1 * (1 - mu) - a2 * (1.0 - np.sqrt(mu))) \
        * r * np.arccos(np.minimum((-rprs ** 2 + z * z + rsq) / (2.0 * z * r), 1.0))


def integral_r_f_sqrt(limb_darkening_coefficients, rprs, z, r1, r2, precision=3):
    return gauss_numerical_integration(num_sqrt, r1, r2, precision, limb_darkening_coefficients, rprs, z)


# dictionaries containing the different methods,
# if you define a new method, include the functions in the dictionary as well

integral_r = {
    'claret': integral_r_claret,
    'linear': integral_r_linear,
    'quad': integral_r_quad,
    'sqrt': integral_r_sqrt,
    'zero': integral_r_zero
}

integral_r_f = {
    'claret': integral_r_f_claret,
    'linear': integral_r_f_linear,
    'quad': integral_r_f_quad,
    'sqrt': integral_r_f_sqrt,
    'zero': integral_r_f_zero,
}


def integral_centred(method, limb_darkening_coefficients, rprs, ww1, ww2):
    return (integral_r[method](limb_darkening_coefficients, rprs)
            - integral_r[method](limb_darkening_coefficients, 0.0)) * np.abs(ww2 - ww1)


def integral_plus_core(method, limb_darkening_coefficients, rprs, z, ww1, ww2, precision=3):
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
    parta = integral_r[method](limb_darkening_coefficients, 0.0) * (w1 - w2)
    partb = integral_r[method](limb_darkening_coefficients, r1) * w2
    partc = integral_r[method](limb_darkening_coefficients, r2) * (-w1)
    partd = integral_r_f[method](limb_darkening_coefficients, rprs, z, r1, r2, precision=precision)
    return parta + partb + partc + partd


def integral_minus_core(method, limb_darkening_coefficients, rprs, z, ww1, ww2, precision=3):
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
    parta = integral_r[method](limb_darkening_coefficients, 0.0) * (w1 - w2)
    partb = integral_r[method](limb_darkening_coefficients, r1) * (-w1)
    partc = integral_r[method](limb_darkening_coefficients, r2) * w2
    partd = integral_r_f[method](limb_darkening_coefficients, rprs, z, r1, r2, precision=precision)
    return parta + partb + partc - partd


def transit_flux_drop(limb_darkening_coefficients, rp_over_rs, z_over_rs, method='claret', precision=3):

    z_over_rs = np.where(z_over_rs < 0, 1.0 + 100.0 * rp_over_rs, z_over_rs)
    z_over_rs = np.maximum(z_over_rs, 10**(-10))

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
    theta_2[case1] = np.pi
    theta_2[case2] = np.pi / 2.0
    theta_2[casea] = np.pi
    theta_2[casec] = np.pi / 2.0
    theta_2[case7] = ph[case7]

    # flux_upper
    plusflux = np.zeros(len(z_over_rs))
    plusflux[plus_case] = integral_plus_core(method, limb_darkening_coefficients, rp_over_rs, z_over_rs[plus_case],
                                             theta_1[plus_case], theta_2[plus_case], precision=precision)
    if len(case0[0]) > 0:
        plusflux[case0] = integral_centred(method, limb_darkening_coefficients, rp_over_rs, 0.0, np.pi)
    if len(caseb[0]) > 0:
        plusflux[caseb] = integral_centred(method, limb_darkening_coefficients, 1, 0.0, np.pi)

    # flux_lower
    minsflux = np.zeros(len(z_over_rs))
    minsflux[minus_case] = integral_minus_core(method, limb_darkening_coefficients, rp_over_rs,
                                               z_over_rs[minus_case], 0.0, theta_2[minus_case], precision=precision)

    # flux_star
    starflux = np.zeros(len(z_over_rs))
    starflux[star_case] = integral_centred(method, limb_darkening_coefficients, 1, 0.0, ph[star_case])

    # flux_total
    total_flux = integral_centred(method, limb_darkening_coefficients, 1, 0.0, 2.0 * np.pi)

    return 1 - (2.0 / total_flux) * (plusflux + starflux - minsflux)


# transit

def transit(limb_darkening_coefficients, rp_over_rs, period, sma_over_rs, eccentricity, inclination, periastron,
            mid_time, time_array, method='claret', precision=3):

    position_vector = planet_orbit(period, sma_over_rs, eccentricity, inclination, periastron, mid_time, time_array)

    projected_distance = np.where(
        position_vector[0] < 0, 1.0 + 5.0 * rp_over_rs,
        np.sqrt(position_vector[1] * position_vector[1] + position_vector[2] * position_vector[2]))

    return transit_flux_drop(limb_darkening_coefficients, rp_over_rs, projected_distance,
                             method=method, precision=precision)


def transit_integrated(limb_darkening_coefficients, rp_over_rs, period, sma_over_rs, eccentricity, inclination,
                       periastron, mid_time, time_array, exp_time, max_sub_exp_time=10, method='claret', precision=3):

    time_factor = int(exp_time / max_sub_exp_time) + 1
    exp_time /= (60.0 * 60.0 * 24.0)

    time_array_hr = (time_array[:, None] + np.arange(-exp_time / 2 + exp_time / time_factor / 2, exp_time / 2,
                                                     exp_time / time_factor)).flatten()

    position_vector = planet_orbit(period, sma_over_rs, eccentricity,
                                   inclination, periastron, mid_time, time_array_hr)

    projected_distance = np.where(
        position_vector[0] < 0, 1.0 + 5.0 * rp_over_rs,
        np.sqrt(position_vector[1] * position_vector[1] + position_vector[2] * position_vector[2]))

    return np.mean(np.reshape(
        transit_flux_drop(limb_darkening_coefficients, rp_over_rs, projected_distance,
                          method=method, precision=precision),
        (len(time_array), time_factor)), 1)


def transit_duration(rp_over_rs, period, sma_over_rs, eccentricity, inclination, periastron):

    ww = periastron * np.pi / 180
    ii = inclination * np.pi / 180
    ee = eccentricity
    aa = sma_over_rs
    ro_pt = (1 - ee ** 2) / (1 + ee * np.sin(ww))
    b_pt = aa * ro_pt * np.cos(ii)
    if b_pt > 1:
        b_pt = 0.5
    s_ps = 1.0 + rp_over_rs
    df = np.arcsin(np.sqrt((s_ps ** 2 - b_pt ** 2) / ((aa ** 2) * (ro_pt ** 2) - b_pt ** 2)))
    aprox = (period * (ro_pt ** 2)) / (np.pi * np.sqrt(1 - ee ** 2)) * df * 60 * 60 * 24

    def function_to_fit(x, t):
        return planet_star_projected_distance(period, sma_over_rs, eccentricity, inclination, periastron,
                                              10000, np.array(10000 + t / 24 / 60 / 60))

    popt1, pcov1 = curve_fit(function_to_fit, [0], [1.0 + rp_over_rs], p0=[-aprox / 2])
    popt2, pcov2 = curve_fit(function_to_fit, [0], [1.0 + rp_over_rs], p0=[aprox / 2])

    return (popt2[0] - popt1[0]) / 24 / 60 / 60


def transit_depth(limb_darkening_coefficients, rp_over_rs, period, sma_over_rs, eccentricity, inclination,
                  periastron, method='claret', precision=6):

    return 1 - transit(limb_darkening_coefficients, rp_over_rs, period, sma_over_rs, eccentricity, inclination,
                       periastron, 10000, np.array([10000]), method=method, precision=precision)[0]


# eclipse

def eclipse(fp_over_fs, rp_over_rs, period, sma_over_rs, eccentricity, inclination, periastron, mid_time,
            time_array, precision=3):

    position_vector = planet_orbit(period, sma_over_rs / rp_over_rs, eccentricity, inclination, periastron + 180,
                                   mid_time, time_array)

    projected_distance = np.where(
        position_vector[0] < 0, 1.0 + 5.0 / rp_over_rs,
        np.sqrt(position_vector[1] * position_vector[1] + position_vector[2] * position_vector[2]))

    return (1.0 + fp_over_fs * transit_flux_drop([0, 0, 0, 0], 1 / rp_over_rs, projected_distance,
                                                 precision=precision, method='zero')) / (1.0 + fp_over_fs)


def eclipse_integrated(fp_over_fs, rp_over_rs, period, sma_over_rs, eccentricity, inclination, periastron,
                       mid_time, time_array, exp_time, max_sub_exp_time=10, precision=3):

    time_factor = int(exp_time / max_sub_exp_time) + 1
    exp_time /= (60.0 * 60.0 * 24.0)

    time_array_hr = (time_array[:, None] + np.arange(-exp_time / 2 + exp_time / time_factor / 2, exp_time / 2,
                                                     exp_time / time_factor)).flatten()

    position_vector = planet_orbit(period, sma_over_rs / rp_over_rs, eccentricity, inclination, periastron + 180,
                                   mid_time, time_array_hr)

    projected_distance = np.where(
        position_vector[0] < 0, 1.0 + 5.0 * rp_over_rs,
        np.sqrt(position_vector[1] * position_vector[1] + position_vector[2] * position_vector[2]))

    return np.mean(np.reshape(
        (1.0 + fp_over_fs * transit_flux_drop([0, 0, 0, 0], 1 / rp_over_rs, projected_distance, method='zero',
                                              precision=precision)) / (1.0 + fp_over_fs),
        (len(time_array), time_factor)), 1)


def eclipse_mid_time(period, sma_over_rs, eccentricity, inclination, periastron, mid_time):
    test_array = np.arange(0, period, 0.001)
    xx, yy, zz = planet_orbit(period, sma_over_rs, eccentricity, inclination, periastron, mid_time,
                              test_array + mid_time)

    test1 = np.where(xx < 0)
    yy = yy[test1]
    test_array = test_array[test1]

    aprox = test_array[np.argmin(np.abs(yy))]

    def function_to_fit(x, t):
        xx, yy, zz = planet_orbit(period, sma_over_rs, eccentricity, inclination, periastron, mid_time,
                            np.array(mid_time + t))
        return yy

    popt, pcov = curve_fit(function_to_fit, [0], [0], p0=[aprox])

    return mid_time + popt[0]


def eclipse_duration(rp_over_rs, period, sma_over_rs, eccentricity, inclination, periastron):

    ww = periastron * np.pi / 180
    ii = inclination * np.pi / 180
    ee = eccentricity
    aa = sma_over_rs
    ro_pt = (1 - ee ** 2) / (1 + ee * np.sin(ww))
    b_pt = aa * ro_pt * np.cos(ii)
    if b_pt > 1:
        b_pt = 0.5
    s_ps = 1.0 + rp_over_rs
    df = np.arcsin(np.sqrt((s_ps ** 2 - b_pt ** 2) / ((aa ** 2) * (ro_pt ** 2) - b_pt ** 2)))
    aprox = (period * (ro_pt ** 2)) / (np.pi * np.sqrt(1 - ee ** 2)) * df * 60 * 60 * 24

    emt = eclipse_mid_time(period, sma_over_rs, eccentricity, inclination, periastron, 10000)

    def function_to_fit(x, t):
        xx = planet_star_projected_distance(period, sma_over_rs, eccentricity, inclination, periastron,
                                              10000, np.array(emt + t / 24 / 60 / 60))
        return xx

    popt1, pcov1 = curve_fit(function_to_fit, [0], [1.0 + rp_over_rs], p0=[-aprox / 2])
    popt2, pcov2 = curve_fit(function_to_fit, [0], [1.0 + rp_over_rs], p0=[aprox / 2])

    return (popt2[0] - popt1[0]) / 24 / 60 / 60


def eclipse_depth(fp_over_fs, rp_over_rs, period, sma_over_rs, eccentricity, inclination, periastron, precision=6):

    return 1 - eclipse(fp_over_fs, rp_over_rs, period, sma_over_rs, eccentricity, inclination,
                       periastron, 10000, np.array([10000]), precision=precision)[0]


def _get_filter(photometric_filter):

    if photometric_filter not in plc_data.all_filters():
        raise PyLCInputError('{0} is not available. Available filters: {1}'.format(
            photometric_filter, ','.join(plc_data.all_filters())))


def fp_over_fs(rp_over_rs, sma_over_rs, albedo, emissivity, stellar_temperature, filter_name):

    _get_filter(filter_name)

    def _black_body(w, t):
        # w in mu
        w = w / (10 ** 6)
        h = 6.62607004 * (10 ** (-34))
        c = 3 * (10 ** 8)
        w5 = w ** 5
        k = 1.38064852 * (10 ** (-23))
        return (2 * h * c * c / w5) / (np.exp(h * c / w / k / t) - 1)

    planet_temperature = stellar_temperature * np.sqrt(0.5 / sma_over_rs) * (((1 - albedo) / emissivity) ** 0.25)

    wavelength_array, band = np.loadtxt(os.path.join(plc_data.photometry(), filter_name + '.pass'), unpack=True)
    binsedge = 0.5 * (wavelength_array[:-1] + wavelength_array[1:])
    binsedge1 = np.append(wavelength_array[0] - (binsedge[0] - wavelength_array[0]), binsedge)
    binsedge2 = np.append(binsedge, wavelength_array[-1] + (wavelength_array[-1] - binsedge[-1]))
    binswidth = binsedge2 - binsedge1

    weights = band * binswidth / 10000
    emission = ((rp_over_rs ** 2) *
                _black_body(wavelength_array / 10000, planet_temperature) /
                _black_body(wavelength_array / 10000, stellar_temperature))

    emission = np.sum(emission * weights) / np.sum(weights)

    reflection = albedo * (rp_over_rs ** 2) / (sma_over_rs ** 2)

    return reflection + emission


def exotethys(stellar_logg, stellar_temperature, stellar_metallicity, filter_name, method='claret',
              stellar_model='phoenix'):

    available_methods = ['claret', 'sqrt', 'quad', 'linear']

    if not isinstance(method, str) or method.lower() not in available_methods:
        raise PyLCInputError('{0} is not available. Available methods: {1}'.format(
            method, ','.join(available_methods)))

    method = method.lower()

    method_map = {
        'claret': 'claret4',
        'sqrt': 'square_root',
        'quad': 'quadratic',
        'linear': 'linear'
    }

    _get_filter(filter_name)

    available_stellar_models = ['atlas', 'phoenix']

    if not isinstance(stellar_model, str) or stellar_model.lower() not in available_stellar_models:
        raise PyLCInputError('{0} is not available. Available stellar models: {1}'.format(
            stellar_model, ','.join(available_stellar_models)))

    stellar_model = stellar_model.lower()

    if stellar_model == 'phoenix' and stellar_metallicity != 0:
        print('PHOENIX models are only computed for solar metallicity stars. Setting stellar_metallicity = 0.')
        stellar_metallicity = 0
    elif stellar_model == 'phoenix' and filter_name in ['irac1', 'irac2', 'irac3', 'irac4']:
        print('PHOENIX models are not available for this filter. Setting stellar_model = atlas')
        stellar_model = 'atlas'

    data = open_dict(os.path.join(plc_data.exotethys(), '{0}_{1}.pickle'.format(stellar_model, filter_name)))

    xin = float(stellar_logg)
    yin = float(stellar_temperature)
    zin = float(stellar_metallicity)

    xunique = np.unique(data['logg'])
    xunique = sorted(xunique, key=lambda x:(x - xin)**2)

    yunique = np.unique(data['teff'])
    yunique = sorted(yunique, key=lambda x:(x - yin)**2)

    zunique = np.unique(data['meta'])
    zunique = sorted(zunique, key=lambda x:(x - zin)**2)

    xmin = 0
    xmax = 0
    ymin = 0
    ymax = 0
    zmin = 0
    zmax = 0

    complete = False
    tsearch = 0

    while not complete and tsearch < len(yunique):

        ymin = min(yunique[:tsearch + 2])
        ymax = max(yunique[:tsearch + 2])

        tsearch += 1

        gsearch = 0

        while not complete and gsearch < len(xunique):

            xmin = min(xunique[:gsearch + 2])
            xmax = max(xunique[:gsearch + 2])

            gsearch += 1

            msearch = 0

            while not complete and msearch < len(zunique):

                zmin = min(zunique[:msearch + 2])
                zmax = max(zunique[:msearch + 2])

                msearch += 1

                # print(xmin, xmax, ymin, ymax, zmin, zmax)

                if all([
                    len(np.where((data['logg'] == ff[0]) & (data['teff'] == ff[1]) & (data['meta'] == ff[2]))[0]) > 0
                        for ff in [
                        [xmin, ymin, zmin],
                        [xmin, ymax, zmin],
                        [xmax, ymin, zmin],
                        [xmax, ymax, zmin],
                        [xmin, ymin, zmax],
                        [xmin, ymax, zmax],
                        [xmax, ymin, zmax],
                        [xmax, ymax, zmax],

                    ]]):

                    if xmin <= xin and xmax >= xin and ymin <= yin and ymax >= yin and zmin <= zin and zmax >= zin:

                        complete = True

    if not complete:
        raise PyLCProcessError('LDCs could not be estimated for the given set of stellar parameters: '
                               'Log(g): {0} to {1} (cgs), Teff: {2} to {3} K, MH: {4} to {5} (dex)'.format(
            xmin, xmax, ymin, ymax, zmin, zmax))
    # else:
    #     print('Trilinear interpolation between: Log(g): {0} to {1} (cgs), Teff: {2} to {3} K, MH: {4} to {5} (dex)'.format(
    #         xmin, xmax, ymin, ymax, zmin, zmax))

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

    for index in [1, 2, 3, 4]:

        vv000 = data[method_map[method]]['ldc' + str(index)][np.where((data['logg'] == xmin) & (data['teff'] == ymin) & (data['meta'] == zmin))][0]
        vv100 = data[method_map[method]]['ldc' + str(index)][np.where((data['logg'] == xmax) & (data['teff'] == ymin) & (data['meta'] == zmin))][0]
        vv010 = data[method_map[method]]['ldc' + str(index)][np.where((data['logg'] == xmin) & (data['teff'] == ymax) & (data['meta'] == zmin))][0]
        vv001 = data[method_map[method]]['ldc' + str(index)][np.where((data['logg'] == xmin) & (data['teff'] == ymin) & (data['meta'] == zmax))][0]
        vv101 = data[method_map[method]]['ldc' + str(index)][np.where((data['logg'] == xmax) & (data['teff'] == ymin) & (data['meta'] == zmax))][0]
        vv011 = data[method_map[method]]['ldc' + str(index)][np.where((data['logg'] == xmin) & (data['teff'] == ymax) & (data['meta'] == zmax))][0]
        vv110 = data[method_map[method]]['ldc' + str(index)][np.where((data['logg'] == xmax) & (data['teff'] == ymax) & (data['meta'] == zmin))][0]
        vv111 = data[method_map[method]]['ldc' + str(index)][np.where((data['logg'] == xmax) & (data['teff'] == ymax) & (data['meta'] == zmax))][0]

        res = tri_linear(xin, yin, zin, xmin, xmax, ymin, ymax, zmin, zmax,
                         vv000, vv100, vv010, vv001, vv101, vv011, vv110, vv111)

        final_coefficients.append(res)

    return np.array(final_coefficients)
