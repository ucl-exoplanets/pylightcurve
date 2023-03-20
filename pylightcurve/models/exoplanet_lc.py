
__all__ = ['planet_orbit', 'planet_star_projected_distance', 'planet_phase',
           'transit', 'transit_integrated', 'transit_depth', 'transit_duration',
           'eclipse', 'eclipse_integrated', 'eclipse_depth', 'eclipse_duration', 'eclipse_mid_time',
           'fp_over_fs','transit_t12', 'exotethys', 'convert_to_bjd_tdb', 'convert_to_jd_utc', 'convert_to_relflux']

import os
import glob
import warnings
import pickle
import numpy as np
import logging

from exotethys import sail
from scipy.optimize import curve_fit as scipy_curve_fit

from ..errors import *
from ..databases import plc_data
from ..analysis.numerical_integration import gauss_numerical_integration
from pylightcurve.processes.files import open_dict


# orbit

def curve_fit(*args, **kwargs):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore",
                                message='Covariance of the parameters could not be estimated')
        return scipy_curve_fit(*args, **kwargs)


def planet_orbit(period, sma_over_rs, eccentricity, inclination, periastron, mid_time, time_array, ww=0):

    inclination = inclination * np.pi / 180.0
    periastron = periastron * np.pi / 180.0
    ww = ww * np.pi / 180.0

    if eccentricity == 0 and ww == 0:

        e_t = (2 * np.pi / period) * (time_array - mid_time)
        cos_e_t = np.cos(e_t)
        sin_e_t = np.sin(e_t)

        x_t = (sma_over_rs * np.sin(inclination)) * cos_e_t
        y_t = sma_over_rs * sin_e_t
        z_t = (- sma_over_rs * np.cos(inclination)) * cos_e_t

    else:

        f_tmid = np.pi / 2 - periastron
        e_tmid = 2 * np.arctan(np.sqrt((1 - eccentricity) / (1 + eccentricity)) * np.tan(f_tmid / 2))
        if e_tmid < 0:
            e_tmid += 2 * np.pi
        tp = mid_time - (period / 2.0 / np.pi) * (e_tmid - eccentricity * np.sin(e_tmid))

        m = (time_array - tp - np.int_((time_array - tp) / period) * period) * 2.0 * np.pi / period
        e_t0 = m
        e_t = e_t0
        stop = False
        for ii in range(10000):  # setting a limit of 10k iterations - arbitrary limit
            e_t = e_t0 - (e_t0 - eccentricity * np.sin(e_t0) - m) / (1 - eccentricity * np.cos(e_t0))
            stop = (np.abs(e_t - e_t0) < 10 ** (-7)).all()
            if stop:
                break
            else:
                e_t0 = e_t
        if not stop:
            raise RuntimeError('Failed to find a solution in 10000 loops')

        f_t = 2 * np.arctan(np.sqrt((1 + eccentricity) / (1 - eccentricity)) * np.tan(e_t / 2))
        r_t = sma_over_rs * (1 - (eccentricity ** 2)) / (1 + eccentricity * np.cos(f_t))
        f_t_plus_periastron = f_t + periastron
        cos_f_t_plus_periastron = np.cos(f_t_plus_periastron)
        sin_f_t_plus_periastron = np.sin(f_t_plus_periastron)

        x_t = r_t * sin_f_t_plus_periastron * np.sin(inclination)

        if ww == 0:
            y_t = - r_t * cos_f_t_plus_periastron
            z_t = - r_t * sin_f_t_plus_periastron * np.cos(inclination)

        else:
            y_t = - r_t * (cos_f_t_plus_periastron * np.cos(ww) - sin_f_t_plus_periastron * np.sin(ww) * np.cos(inclination))
            z_t = - r_t * (cos_f_t_plus_periastron * np.sin(ww) + sin_f_t_plus_periastron * np.cos(ww) * np.cos(inclination))

    return [x_t, y_t, z_t]


def planet_star_projected_distance(period, sma_over_rs, eccentricity, inclination, periastron, mid_time, time_array):

    x_t, y_t, z_t = planet_orbit(period, sma_over_rs, eccentricity, inclination, periastron, mid_time, time_array)

    return np.sqrt(y_t * y_t + z_t * z_t)


def planet_phase(period, mid_time, time_array):

    phase = (time_array - mid_time)/period

    return phase - np.int_(phase)

# flux drop - new

class Integrals:

    def __init__(self):

        self.methods = ['linear', 'quad', 'claret']
        self.kappa = {
            'linear': self.kappa_linear,
            'quad': self.kappa_quad,
            'claret': self.kappa_claret,
            'eclipse': self.kappa_eclipse
        }

    def kappa_claret(self, limb_darkening_coefficients, r):
        a1, a2, a3, a4 = limb_darkening_coefficients

        mu44 = np.abs(1.0 - r * r)
        mu24 = np.sqrt(mu44)
        mu14 = np.sqrt(mu24)
        return - mu44 * (((1.0 - a1 - a2 - a3 - a4) / 2) + (2 * a1 / 5) * mu14 + (2 * a2 / 6) * mu24 + (2 * a3 / 7) * mu24 * mu14 + (2 * a4 / 8) * mu44)

    def kappa_quad(self, limb_darkening_coefficients, r):
        a1, a2 = limb_darkening_coefficients[:2]
        musq = np.abs(1.0 - r * r)
        mu = np.sqrt(musq)

        return (1.0 / 12) * (-4.0 * (a1 + 2.0 * a2) * mu * musq + 6.0 * (-1 + a1 + a2) * musq + 3.0 * a2 * musq * musq)

    def kappa_linear(self, limb_darkening_coefficients, r):
        a1 = limb_darkening_coefficients[0]
        musq = np.abs(1.0 - r * r)

        return (-1.0 / 6) * musq * (3.0 + a1 * (-3.0 + 2.0 * np.sqrt(musq)))

    def kappa_eclipse(self, limb_darkening_coefficients, r):
        musq = np.abs(1.0 - r * r)

        return (-1.0 / 6) * musq * 3.0

    def kappa_r(self, u, limb_darkening_coefficients, rprs, z, zsq, sign, method):

        a = z * np.cos(u)
        r = a + sign * np.sqrt(np.abs(rprs*rprs - zsq + a*a))

        return self.kappa[method](limb_darkening_coefficients, r)

    def integral(self, limb_darkening_coefficients, rprs, z, zsq, sign, method, u1, u2, precision):

        spoint = 0.99

        if precision == 'ref':

            u11 = u1
            u12 = u1 + spoint * (u2-u1)
            u21 = u1 + spoint * (u2-u1)
            u22 = u2

            result = np.zeros(len(z))
            du = (u12-u11)/1000
            for ii in np.arange(0, 1000, 1):
                result += gauss_numerical_integration(self.kappa_r, u11 + ii * du, u11 + (ii + 1) * du, 10, limb_darkening_coefficients, rprs, z, zsq, sign, method)
            du = (u22-u21)/1000
            for ii in np.arange(0, 1000, 1):
                result += gauss_numerical_integration(self.kappa_r, u21 + ii * du, u21 + (ii + 1) * du, 10, limb_darkening_coefficients, rprs, z, zsq, sign, method)

            return result

        else:

            return (gauss_numerical_integration(self.kappa_r, u1, u1 + spoint * (u2-u1), precision, limb_darkening_coefficients, rprs, z, zsq, sign, method)
                    + gauss_numerical_integration(self.kappa_r, u1 + spoint * (u2-u1), u2, precision, limb_darkening_coefficients, rprs, z, zsq, sign, method))

    def kappa_zero(self, limb_darkening_coefficients, method):
        return self.kappa[method](limb_darkening_coefficients, 0)


integral_class = Integrals()


def transit_flux_drop(limb_darkening_coefficients, rp_over_rs, z_over_rs, method, precision):

    len_z_over_rs = len(z_over_rs)

    zsq = z_over_rs * z_over_rs

    test1 = (z_over_rs > -1 + rp_over_rs) * (z_over_rs <= rp_over_rs)
    test2 = (z_over_rs > rp_over_rs) * (z_over_rs < 1 + rp_over_rs)
    test3 = (z_over_rs <= -1 + rp_over_rs)
    test4 = z_over_rs <= 1 - rp_over_rs
    test5 = zsq < rp_over_rs ** 2 + 1
    test6 = zsq >= rp_over_rs ** 2 - 1

    case1 = np.where(test1 * test4)
    case2 = np.where(test2 * test4)
    case3 = np.where(test2 * (~test4) * test5)
    case4 = np.where(test2 * (~test4) * (~test5))
    case5 = np.where(test1 * (~test4) * test6)
    case5b = np.where(test1 * (~test4) * (~test6))
    case6 = np.where(test3)

    # flux_total
    pi_kappa_zero = np.pi * integral_class.kappa_zero(limb_darkening_coefficients, method)

    semi_df = np.zeros(len_z_over_rs)
    u1 = np.zeros(len_z_over_rs)
    u2 = np.zeros(len_z_over_rs)

    u1[case3] = np.arccos((1.0 - rp_over_rs ** 2 + zsq[case3]) / (2.0 * z_over_rs[case3]))
    u1[case5] = 3*np.pi/2
    u1[case5b] = np.pi
    u2[case1] = np.pi/2
    u2[case2] = np.arcsin(rp_over_rs / z_over_rs[case2])
    u2[case3] = np.arcsin(rp_over_rs / z_over_rs[case3])
    u2[case4] = np.arccos((1.0 - rp_over_rs ** 2 + zsq[case4]) / (2.0 * z_over_rs[case4]))
    u2[case5] = 2*np.pi - np.arccos((1.0 - rp_over_rs ** 2 + zsq[case5]) / (2.0 * z_over_rs[case5]))
    u2[case5b] = 2*np.pi - np.arccos((1.0 - rp_over_rs ** 2 + zsq[case5b]) / (2.0 * z_over_rs[case5b]))

    case_plus = np.concatenate([case1[0], case2[0], case3[0], case5[0], case5b[0]])
    semi_df[case_plus] += integral_class.integral(limb_darkening_coefficients, rp_over_rs, z_over_rs[case_plus], zsq[case_plus], 1, method, u1[case_plus], u2[case_plus], precision)

    case_pi = np.concatenate([case1[0], case5[0]])
    semi_df[case_pi] -= integral_class.integral(limb_darkening_coefficients, rp_over_rs, z_over_rs[case_pi], zsq[case_pi], +1, method, np.pi, np.ones_like(case_pi) * np.pi/2, precision)

    case_minus = np.concatenate([case2[0], case3[0], case4[0]])
    semi_df[case_minus] -= integral_class.integral(limb_darkening_coefficients, rp_over_rs, z_over_rs[case_minus], zsq[case_minus], -1, method, 0, u2[case_minus], precision)

    case_star = np.concatenate([case1[0], case5[0], case5b[0], case6[0]])
    semi_df[case_star] -= pi_kappa_zero

    return 1 + (semi_df / pi_kappa_zero)


# transit

def transit(limb_darkening_coefficients, rp_over_rs, period, sma_over_rs, eccentricity, inclination, periastron,
            mid_time, time_array, method='claret', precision=2):

    position_vector = planet_orbit(period, sma_over_rs, eccentricity, inclination, periastron, mid_time, time_array)

    projected_distance = np.where(
        position_vector[0] < 0, 1.0 + 10.0 * rp_over_rs,
        np.sqrt(position_vector[1] * position_vector[1] + position_vector[2] * position_vector[2]))

    return transit_flux_drop(limb_darkening_coefficients, rp_over_rs, projected_distance, method, precision)


def transit_integrated(limb_darkening_coefficients, rp_over_rs, period, sma_over_rs, eccentricity, inclination,
                       periastron, mid_time, time_array, exp_time, max_sub_exp_time=10, method='claret', precision=2):

    time_factor = int(exp_time / max_sub_exp_time) + 1
    exp_time /= (60.0 * 60.0 * 24.0)

    time_array_hr = (time_array[:, None] + np.arange(-exp_time / 2 + exp_time / time_factor / 2, exp_time / 2,
                                                     exp_time / time_factor)).flatten()

    position_vector = planet_orbit(period, sma_over_rs, eccentricity,
                                   inclination, periastron, mid_time, time_array_hr)

    projected_distance = np.where(
        position_vector[0] < 0, 1.0 + 10.0 * rp_over_rs,
        np.sqrt(position_vector[1] * position_vector[1] + position_vector[2] * position_vector[2]))

    return np.mean(np.reshape(
        transit_flux_drop(limb_darkening_coefficients, rp_over_rs, projected_distance, method, precision),
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


def transit_t12(rp_over_rs, period, sma_over_rs, eccentricity, inclination, periastron):

    aprox = transit_duration(rp_over_rs, period, sma_over_rs, eccentricity, inclination, periastron) * 60 * 60 * 24

    def function_to_fit(x, t):
        return planet_star_projected_distance(period, sma_over_rs, eccentricity, inclination, periastron,
                                              10000, np.array(10000 + t / 24 / 60 / 60))

    popt1, pcov1 = curve_fit(function_to_fit, [0], [1.0 + rp_over_rs], p0=[-aprox / 2])
    popt2, pcov2 = curve_fit(function_to_fit, [0], [1.0 - rp_over_rs], p0=[-aprox / 2])

    return min((popt2[0] - popt1[0]) / 24 / 60 / 60,
               0.5*transit_duration(rp_over_rs, period, sma_over_rs, eccentricity, inclination, periastron)
               )


def transit_depth(limb_darkening_coefficients, rp_over_rs, period, sma_over_rs, eccentricity, inclination,
                  periastron, method='claret', precision=2):

    return 1 - transit(limb_darkening_coefficients, rp_over_rs, period, sma_over_rs, eccentricity, inclination,
                       periastron, 10000, np.array([10000]), method, precision)[0]


# eclipse

def eclipse(fp_over_fs, rp_over_rs, period, sma_over_rs, eccentricity, inclination, periastron, mid_time,
            time_array, precision=2):

    position_vector = planet_orbit(period, sma_over_rs / rp_over_rs, eccentricity, inclination, periastron + 180,
                                   mid_time, time_array)

    projected_distance = np.where(
        position_vector[0] < 0, 1.0 + 5.0 / rp_over_rs,
        np.sqrt(position_vector[1] * position_vector[1] + position_vector[2] * position_vector[2]))

    return (1.0 + fp_over_fs * transit_flux_drop([0, 0, 0, 0], 1 / rp_over_rs, projected_distance,
                                                 'eclipse', precision)) / (1.0 + fp_over_fs)


def eclipse_integrated(fp_over_fs, rp_over_rs, period, sma_over_rs, eccentricity, inclination, periastron,
                       mid_time, time_array, exp_time, max_sub_exp_time=10, precision=2):

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
        (1.0 + fp_over_fs * transit_flux_drop([0, 0, 0, 0], 1 / rp_over_rs, projected_distance, 'eclipse',
                                              precision)) / (1.0 + fp_over_fs),
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

    aprox = transit_duration(rp_over_rs, period, sma_over_rs, eccentricity, inclination, periastron) * 60 * 60 * 24

    emt = eclipse_mid_time(period, sma_over_rs, eccentricity, inclination, periastron, 10000)

    def function_to_fit(x, t):
        xx = planet_star_projected_distance(period, sma_over_rs, eccentricity, inclination, periastron,
                                            10000, np.array(emt + t / 24 / 60 / 60))
        return xx

    popt1, pcov1 = curve_fit(function_to_fit, [0], [1.0 + rp_over_rs], p0=[-aprox / 2])
    popt2, pcov2 = curve_fit(function_to_fit, [0], [1.0 + rp_over_rs], p0=[aprox / 2])

    return (popt2[0] - popt1[0]) / 24 / 60 / 60


def eclipse_depth(fp_over_fs, rp_over_rs, period, sma_over_rs, eccentricity, inclination, periastron, precision=2):

    return 1 - eclipse(fp_over_fs, rp_over_rs, period, sma_over_rs, eccentricity, inclination,
                       periastron, 10000, np.array([10000]), precision)[0]


def fp_over_fs(rp_over_rs, sma_over_rs, albedo, emissivity, stellar_temperature, filter_name, wlrange=None):

    def _black_body(w, t):
        # w in mu
        w = w / (10 ** 6)
        h = 6.62607004 * (10 ** (-34))
        c = 3 * (10 ** 8)
        w5 = w ** 5
        k = 1.38064852 * (10 ** (-23))
        return (2 * h * c * c / w5) / (np.exp(h * c / w / k / t) - 1)

    planet_temperature = stellar_temperature * np.sqrt(0.5 / sma_over_rs) * (((1 - albedo) / emissivity) ** 0.25)

    sub_passband = plc_data.get_filter(filter_name, wlrange)
    wavelength_array, band = np.swapaxes(sub_passband, 0, 1)

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


def exotethys(stellar_logg, stellar_temperature, stellar_metallicity, filter_name, wlrange=None, ldc_method='claret',
              stellar_model='Phoenix_2018'):
    """Calculates the limb-darkening coefficients for the host star, by calling the ExoTehys package.
     Not all the combinations of filters and stellar models are available. This will depend on the stellar temperature,
     the stellar model and the wavelength range of the filter.
     Please check the `ExoTETHyS <https://github.com/ucl-exoplanets/ExoTETHyS>`_,
     `Morello et al. 2020 <https://iopscience.iop.org/article/10.3847/1538-3881/ab63dc>`_
     documentation for all the available parameter space.

     :param filter_name: Observing filter
     :type filter_name: str
     :param wlrange: A 2-element list, defining the wavelength range (in Angstrom) within the filter, Useful for spectroscopic observations.
     :type wlrange: list
     :param ldc_method: Limb-darkening model to de used. Available methods: claret, quad, linear, sqrt
     :type ldc_method: str
     :param stellar_model: Stellar model to be used. Available methods: check ExoTethys documentation.
     :type stellar_model: str
     :return: Limb-darkening coefficients. Four coefficients are always returned (for methods that need less than four coefficents, the rest are 0).
     :rtype: numpy.ndarray
     """

    if 'phoenix' in stellar_model.lower() and stellar_metallicity != 0:
        logging.warning('\nPHOENIX models are only computed for solar metallicity stars. '
                      'Setting stellar_metallicity = 0.\n')
        stellar_metallicity = 0

    path = plc_data.exotethys()

    available_stellar_models = ['Atlas_2000', 'Phoenix_2012_13', 'Phoenix_2018', 'Phoenix_drift_2012',
                                'Stagger_2015', 'Stagger_2018']
    if not isinstance(stellar_model, str) or stellar_model not in available_stellar_models:
        raise PyLCInputError('Stellar_model {0} is not available. Available models: {1}. '
                             'Please consult EXOTETHYS documentation for more details'.format(
            stellar_model, ','.join(available_stellar_models)))

    method_map = {
        'claret': 'claret4',
        'quad': 'quadratic',
        'linear': 'linear'
    }
    if ldc_method not in method_map:
        raise PyLCInputError('Method {0} is not valid. Available methods: {1}'.format(
            ldc_method, ','.join(list(method_map.keys()))))

    run_id = ''.join([str(np.random.randint(0, 99)) for ff in range(10)])
    bandpass_file = os.path.join(path, 'ww_{0}.pass'.format(run_id))
    parameters_file = os.path.join(path, 'ww_{0}_parameters.txt'.format(run_id))
    output_file = os.path.join(path, 'ww_{0}_ldc.pickle'.format(run_id))

    sub_passband = plc_data.get_filter(filter_name, wlrange)
    np.savetxt(bandpass_file, sub_passband)

    r = open(os.path.join(path, 'ww_parameters.txt')).read()
    r = r.replace('{{output}}', path)
    r = r.replace('{{id}}', str(run_id))
    r = r.replace('{{model}}', stellar_model)
    r = r.replace('{{law}}', method_map[ldc_method])
    r = r.replace('{{teff}}', str(stellar_temperature))
    r = r.replace('{{logg}}', str(stellar_logg))
    r = r.replace('{{meta}}', str(stellar_metallicity))

    w = open(parameters_file, 'w')
    w.write(r)
    w.close()

    try:
        sail.ldc_calculate(parameters_file)

        results = open_dict(output_file)

        for ff in glob.glob(os.path.join(path, '*{0}*'.format(run_id))):
            os.remove(ff)

        ldcs = results['passbands']['ww_{0}.pass'.format(run_id)][method_map[ldc_method]]['coefficients']

        while len(ldcs) < 4:
            ldcs = np.append(ldcs, 0)

        return ldcs

    except pickle.UnpicklingError:
        plc_data.reset_exotethys_data()
        return exotethys(stellar_logg, stellar_temperature, stellar_metallicity, filter_name, wlrange, ldc_method,
                         stellar_model)


def _hjd_utc(target, jd):

    ra, dec = target
    ra *= np.pi / 180
    dec *= np.pi / 180

    ssb_ra, ssb_dec, ssb_d, ssb_dt = plc_data.heliocentre(jd)

    a = ssb_d / 60.0 / 24.0
    b = np.sin(dec) * np.sin(ssb_dec)
    c = np.cos(dec) * np.cos(ssb_dec) * np.cos(ra - ssb_ra)

    return jd - a * (b + c)


def _hjd_tdb(target, jd):

    ra, dec = target
    ra *= np.pi / 180
    dec *= np.pi / 180

    ssb_ra, ssb_dec, ssb_d, ssb_dt = plc_data.heliocentre(jd)

    a = ssb_d / 60.0 / 24.0
    b = np.sin(dec) * np.sin(ssb_dec)
    c = np.cos(dec) * np.cos(ssb_dec) * np.cos(ra - ssb_ra)

    return jd - a * (b + c) + ssb_dt / 60.0 / 60.0 / 24.0


def _bjd_utc(target, jd):

    ra, dec = target
    ra *= np.pi / 180
    dec *= np.pi / 180

    ssb_ra, ssb_dec, ssb_d, ssb_dt = plc_data.barycentre(jd)

    a = ssb_d / 60.0 / 24.0
    b = np.sin(dec) * np.sin(ssb_dec)
    c = np.cos(dec) * np.cos(ssb_dec) * np.cos(ra - ssb_ra)

    return jd - a * (b + c)


def _bjd_tdb(target, jd):

    ra, dec = target
    ra *= np.pi / 180
    dec *= np.pi / 180

    ssb_ra, ssb_dec, ssb_d, ssb_dt = plc_data.barycentre(jd)

    a = ssb_d / 60.0 / 24.0
    b = np.sin(dec) * np.sin(ssb_dec)
    c = np.cos(dec) * np.cos(ssb_dec) * np.cos(ra - ssb_ra)

    return jd - a * (b + c) + ssb_dt / 60.0 / 60.0 / 24.0


def convert_to_bjd_tdb(ra, dec, time_array, time_format):

    return_float = False
    if isinstance(time_array, float):
        time_array = np.array([time_array])
        return_float = True
    elif isinstance(time_array, np.ndarray):
        try:
            time_array = np.array(time_array, dtype=float)
        except:
            raise PyLCInputError('The time_array could not be converted to dtype=float')
    else:
        raise PyLCInputError('Not acceptable time_array type {0}. Please provide an np.ndarray '
                             'or a float object.'.format(type(time_array)))

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore",
                                message='Covariance of the parameters could not be estimated')

        if time_format in ['BJD_TDB', 'BJD_TT']:
            new_time_array = time_array
        elif time_format == 'JD_UTC':
            new_time_array = np.array([_bjd_tdb((ra, dec), ff) for ff in time_array])
        elif time_format == 'MJD_UTC':
            new_time_array = np.array([_bjd_tdb((ra, dec), ff + 2400000.5) for ff in time_array])
        elif time_format == 'BJD_UTC':
            new_time_array = np.array([_bjd_tdb((ra, dec), curve_fit(_bjd_utc, (ra, dec), [ff], p0=[ff])[0][0])
                             for ff in time_array])
        elif time_format in ['HJD_TDB', 'HJD_TT']:
            new_time_array = np.array([_bjd_tdb((ra, dec), curve_fit(_hjd_tdb, (ra, dec), [ff], p0=[ff])[0][0])
                             for ff in time_array])
        elif time_format == 'HJD_UTC':
            new_time_array = np.array([_bjd_tdb((ra, dec), curve_fit(_hjd_utc, (ra, dec), [ff], p0=[ff])[0][0])
                             for ff in time_array])
        else:
            raise PyLCInputError(
                'Not valid time format. Available formats: JD_UTC, MJD_UTC, BJD_UTC, BJD_TDB, BJD_TT, '
                'HJD_UTC, HJD_BJD, HJD_TT')

    if return_float:
        return new_time_array[0]
    else:
        return new_time_array


def convert_to_jd_utc(ra, dec, time_array, time_format):

    return_float = False
    if isinstance(time_array, float):
        time_array = np.array([time_array])
        return_float = True
    elif isinstance(time_array, np.ndarray):
        try:
            time_array = np.array(time_array, dtype=float)
        except:
            raise PyLCInputError('The time_array could not be converted to dtype=float')
    else:
        raise PyLCInputError('Not acceptable time_array type {0}. Please provide an np.ndarray '
                             'or a float object.'.format(type(time_array)))

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore",
                                message='Covariance of the parameters could not be estimated')

        if time_format in ['BJD_TDB', 'BJD_TT']:
            new_time_array = np.array([curve_fit(_bjd_tdb, (ra, dec), [ff], p0=[ff])[0][0] for ff in time_array])
        elif time_format == 'JD_UTC':
            new_time_array = time_array
        elif time_format == 'MJD_UTC':
            new_time_array = time_array + 2400000.5
        elif time_format == 'BJD_UTC':
            new_time_array = np.array([curve_fit(_bjd_utc, (ra, dec), [ff], p0=[ff])[0][0] for ff in time_array])
        elif time_format in ['HJD_TDB', 'HJD_TT']:
            new_time_array = np.array([curve_fit(_hjd_tdb, (ra, dec), [ff], p0=[ff])[0][0] for ff in time_array])
        elif time_format == 'HJD_UTC':
            new_time_array = np.array([curve_fit(_hjd_utc, (ra, dec), [ff], p0=[ff])[0][0] for ff in time_array])
        else:
            raise PyLCInputError(
                'Not valid time format. Available formats: JD_UTC, MJD_UTC, BJD_UTC, BJD_TDB, BJD_TT, '
                'HJD_UTC, HJD_BJD, HJD_TT')

    if return_float:
        return new_time_array[0]
    else:
        return new_time_array


def convert_to_relflux(flux_array, flux_unc_array, flux_format):

    if flux_format == 'mag':
        flux_unc_array = np.abs(flux_unc_array * (-0.921034 * np.exp(0.921034 * flux_array[0] - 0.921034 * np.array(flux_array))))
        flux_array = 10 ** ((flux_array[0] - flux_array) / 2.5)
    elif flux_format == 'flux':
        pass
    else:
        raise PyLCInputError('Not acceptable flux format {0}. Please choose between "flux" and '
                             '"mag".'.format(flux_format))

    return flux_array, flux_unc_array
