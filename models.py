__all__ = ['transit', 'eclipse', 'mcmc_transit']

import pylightcurve_tools
import numpy as np
import os
import pymc


class PYLCError(BaseException):
    pass


class PYLCMcmcError(PYLCError):
    pass


def transit(limb_darkening_coefficients, rp_over_rs,
            period, sma_over_rs, eccenticity, inclination, periastron, mid_time, time_array):

    position_vector = pylightcurve_tools.position_vector(period, sma_over_rs,
                                                         eccenticity, inclination, periastron, mid_time, time_array)

    projected_distance = np.where(
        position_vector[0] < 0, 1.0 + 5.0 * rp_over_rs,
        np.sqrt(position_vector[1] * position_vector[1] + position_vector[2] * position_vector[2]))

    return pylightcurve_tools.flux_drop(limb_darkening_coefficients, rp_over_rs, projected_distance)


def eclipse(fp_over_fs, rp_over_rs, period, sma_over_rs, eccenticity, inclination, periastron, mid_time, time_array):

    position_vector = pylightcurve_tools.position_vector(period, -sma_over_rs / rp_over_rs,
                                                         eccenticity, inclination, periastron, mid_time, time_array)

    projected_distance = np.where(
        position_vector[0] < 0, 1.0 + 5.0 * rp_over_rs,
        np.sqrt(position_vector[1] * position_vector[1] + position_vector[2] * position_vector[2]))

    return (1.0 + fp_over_fs * pylightcurve_tools.flux_drop((0, 0, 0, 0), 1 / rp_over_rs, projected_distance)) \
        / (1.0 + fp_over_fs)


def mcmc_transit(limb_darkening_coefficients, rp_over_rs,
                 period, sma_over_rs, eccenticity, inclination, periastron, mid_time,
                 data, fit_rp_over_rs, iterations, burn, directory, detrend_order=0,
                 fit_period=None, fit_sma_over_rs=None, fit_eccentricity=None,
                 fit_inclination=None, fit_periastron=None, fit_mid_time=None):

    if isinstance(data, basestring):
        data = np.loadtxt(data, unpack=True)

    if len(data) == 3:
        datax, datay, dataz = data
    else:
        datax, datay = data
        dataz = np.ones_like(datax) * np.sqrt(np.mean((datay[1:-1] - 0.5 * (datay[:-2] + datay[2:])) ** 2))

    if not os.path.isdir(directory):
        os.mkdir(directory)

    initial_directory = os.path.abspath('.')
    os.chdir(directory)
    os.system('rm *')

    if detrend_order not in [0, 1, 2]:
        raise PYLCMcmcError('Value for detrend_order must be 0, 1, or 2')

    detrend_zero_initial = np.mean(datay)
    detrend_zero_limits = [0.5 * np.mean(datay), 1.5 * np.mean(datay)]

    detrend_one_initial = 0
    detrend_two_initial = 0
    max_detrend_limits = 10 * (max(datay) - min(datay)) / (max(datax) - min(datax)) / np.mean(datay)

    if detrend_order == 0:
        detrend_one_limits = None
        detrend_two_limits = None
    elif detrend_order == 1:
        detrend_one_limits = [-max_detrend_limits, max_detrend_limits]
        detrend_two_limits = None
    else:
        detrend_one_limits = [-max_detrend_limits, max_detrend_limits]
        detrend_two_limits = [-max_detrend_limits, max_detrend_limits]

    names = ['N', 'L', 'Q', 'rp', 'P', 'a', 'e', 'i', 'w', 'mt']
    initial_values = [detrend_zero_initial, detrend_one_initial, detrend_two_initial,
                      rp_over_rs, period, sma_over_rs, eccenticity, inclination, periastron, mid_time]
    limits = [detrend_zero_limits, detrend_one_limits, detrend_two_limits,
              fit_rp_over_rs, fit_period, fit_sma_over_rs,
              fit_eccentricity, fit_inclination, fit_periastron, fit_mid_time]
    variables = []
    traces = []
    results = []
    errors = []

    for var in range(len(names)):
        if limits[var] is None:
            variables.append(initial_values[var])
        else:
            if initial_values[var] < limits[var][0] or initial_values[var] > limits[var][1]:
                raise PYLCMcmcError('Initial value for ' + names[var] + 'is outside the range of the prior')
            else:
                variables.append(pymc.Uniform(names[var], limits[var][0], limits[var][1], value=initial_values[var]))

    @pymc.deterministic
    def mcmc_f(model_variables=variables):
        d_time = datax - datax[0]
        return model_variables[0] * (1 + model_variables[1] * d_time + model_variables[2] * d_time * d_time) \
                                  * transit(limb_darkening_coefficients, *model_variables[3:], time_array=datax)

    y = pymc.Normal('y', mu=mcmc_f, tau=dataz ** (-2), observed=True, value=datay)
    mcmc = pymc.MCMC([variables, mcmc_f, y], db='pickle', dbname='model.pickle')
    mcmc.isample(iterations, burn=burn, verbose=1)

    for var in range(len(names)):
        if limits[var] is None:
            traces.append(np.array([initial_values[var]]))
        else:
            traces.append(np.array(mcmc.trace(names[var])[:]))

    for var in range(len(names)):
        if limits[var] is not None:
            try:
                dist = pylightcurve_tools.distribution(traces[var])
            except RuntimeError:
                print " gaussian could not be fitted for: " + names[var]
                dist = [0, np.mean(traces[var]), np.std(traces[var])]
            results.append(dist[1])
            errors.append(dist[2])

        else:
            results.append(initial_values[var])
            errors.append(0)

    def final_model(time_array):
        d_time = time_array - time_array[0]
        return results[0] * (1 + results[1] * d_time + results[2] * d_time * d_time) \
                          * transit(limb_darkening_coefficients, *results[3:], time_array=time_array)

    print 'Saving results...'
    pylightcurve_tools.save_results(names, results, errors)
    pylightcurve_tools.save_model(datax, datay, final_model)
    pylightcurve_tools.save_traces(names, traces)

    print 'Plotting...'
    pylightcurve_tools.plot_model(datax, datay, final_model)
    pylightcurve_tools.plot_traces(names, traces, results, errors)
    pylightcurve_tools.plot_correlations(names, traces, results, errors)

    os.system('rm model.pickle')
    os.chdir(initial_directory)

