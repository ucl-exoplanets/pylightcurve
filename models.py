__all__ = ['transit', 'eclipse', 'mcmc_transit']

import pylightcurve_tools
import flux
import position

import numpy as np
import os
import pymc


class PYLCError(BaseException):
    pass


class PYLCMcmcError(PYLCError):
    pass


def transit(limb_darkening_coefficients, rp_over_rs,
            period, sma_over_rs, eccentricity, inclination, periastron, mid_time, time_array,
            method='claret', precision=0):

    position_vector = position.position_vector(period, sma_over_rs,
                                               eccentricity, inclination, periastron, mid_time, time_array)

    projected_distance = np.where(
        position_vector[0] < 0, 1.0 + 5.0 * rp_over_rs,
        np.sqrt(position_vector[1] * position_vector[1] + position_vector[2] * position_vector[2]))

    return flux.flux_drop(limb_darkening_coefficients, rp_over_rs, projected_distance,
                          method=method, precision=precision)


def eclipse(fp_over_fs, rp_over_rs, period, sma_over_rs, eccentricity, inclination, periastron, mid_time, time_array):

    position_vector = position.position_vector(period, -sma_over_rs / rp_over_rs,
                                               eccentricity, inclination, periastron, mid_time, time_array)

    projected_distance = np.where(
        position_vector[0] < 0, 1.0 + 5.0 / rp_over_rs,
        np.sqrt(position_vector[1] * position_vector[1] + position_vector[2] * position_vector[2]))

    return (1.0 + fp_over_fs * flux.flux_drop((0, 0, 0, 0), 1 / rp_over_rs, projected_distance)) \
        / (1.0 + fp_over_fs)


def mcmc_transit(limb_darkening_coefficients, rp_over_rs,
                 period, sma_over_rs, eccentricity, inclination, periastron, mid_time,
                 data, fit_rp_over_rs, iterations, burn, directory, detrend_order=0,
                 fit_period=None, fit_sma_over_rs=None, fit_eccentricity=None,
                 fit_inclination=None, fit_periastron=None, fit_mid_time=None,
                 method='claret', precision=0, exp_time=0, time_factor=1):

    if method not in ['linear', 'quad', 'sqrt', 'claret']:
        raise PYLCMcmcError('method argument must be linear, quad, sqrt or claret')

    if len(limb_darkening_coefficients) != 4:
        raise PYLCMcmcError('limb_darkening_coefficients argument must be a 4-item array-like object')

    try:
        iterations = int(iterations)
    except:
        raise PYLCMcmcError('iterations argument must be an integer')

    try:
        burn = int(burn)
    except:
        raise PYLCMcmcError('burn argument must be an integer')

    if not isinstance(directory, basestring):
        raise PYLCMcmcError('directory argument must be a string')

    if detrend_order not in [0, 1, 2]:
        raise PYLCMcmcError('detrend_order argument must be 0, 1, or 2')

    try:
        exp_time *= (1.0 / 60 / 60 / 24)
    except:
        raise PYLCMcmcError('exp_time argument must be a float (in seconds)')

    try:
        time_factor = int(time_factor)
    except:
        raise PYLCMcmcError('time_factor argument must be an integer')

    if not os.path.isdir(directory):
        os.mkdir(directory)

    initial_directory = os.path.abspath('.')
    os.chdir(directory)

    if len(np.array(data[0]).shape) == 1:
        sets_n = 1
        data = [data]
    elif len(np.array(data[0]).shape) in [2, 3]:
        sets_n = np.array(data).shape[0]
    else:
        raise PYLCMcmcError('Improper input data format')

    datax = np.array([])
    datay = np.array([])
    dataz = np.array([])
    datai = np.array([])
    datadt = np.array([])

    data_mid_time = 0
    for set_n in range(sets_n):
        dataset = data[set_n]

        if len(dataset) == 3:
            datasetx, datasety, datasetz = dataset
        else:
            datasetx, datasety = dataset
            datasetz = np.ones_like(datasetx) * \
                np.std(datasety[1:-1] - 0.5 * (datasety[:-2] + datasety[2:]))

        if set_n == 0:
            data_mid_time = np.mean(datasetx)

        datax = np.append(datax, datasetx)
        datay = np.append(datay, datasety)
        dataz = np.append(dataz, datasetz)
        datai = np.append(datai, np.ones_like(datasetx) * set_n)
        datadt = np.append(datadt, datasetx - datasetx[0])

    datai = np.int_(datai)

    names = []
    initial = []
    limits = []
    variables = []
    traces = []
    results = []
    errors = []

    for set_n in range(sets_n):
        dataset = data[set_n]

        if len(dataset) == 3:
            datasetx, datasety, datasetz = dataset
        else:
            datasetx, datasety = dataset

        max_limit = 10 * (max(datasety) - min(datasety)) / (max(datasetx) - min(datasetx)) / np.mean(datasety)

        names.append('N' + str(set_n))
        initial.append(np.mean(datasety))
        limits.append([0.5 * np.mean(datay), 1.5 * np.mean(datay)])

        names.append('L' + str(set_n))
        initial.append(0)
        if detrend_order in [1, 2]:
            limits.append([-max_limit, max_limit])
        else:
            limits.append(None)

        names.append('Q' + str(set_n))
        initial.append(0)
        if detrend_order == 2:
            limits.append([-max_limit, max_limit])
        else:
            limits.append(None)

    time_shift = round((data_mid_time - mid_time) / period)
    mid_time += time_shift * period
    fit_mid_time = [fit_mid_time[0] + time_shift * period, fit_mid_time[1] + time_shift * period]

    names += ['rp', 'P', 'a', 'e', 'i', 'w', 'mt']
    initial += [rp_over_rs, period, sma_over_rs, eccentricity, inclination, periastron, mid_time]
    limits += [fit_rp_over_rs, fit_period, fit_sma_over_rs,
               fit_eccentricity, fit_inclination, fit_periastron, fit_mid_time]

    for var in range(len(names)):

        try:
            initial[var] = float(initial[var])
        except:
            raise PYLCMcmcError('Improper value for ' + names[var])

        if limits[var] is None:
            variables.append(initial[var])
        elif isinstance(limits[var], float):
            variables.append(pymc.Normal(names[var], initial[var], 1.0 / (limits[var] ** 2)))
            limits[var] = None
        else:
            try:
                if len(limits[var]) != 2:
                    raise PYLCMcmcError('Improper limits for ' + names[var])
            except:
                raise PYLCMcmcError('Improper limits for ' + names[var])

            if initial[var] < limits[var][0] or initial[var] > limits[var][1]:
                raise PYLCMcmcError('Initial value for ' + names[var] + ' is outside the range of the prior')
            else:
                variables.append(pymc.Uniform(names[var], limits[var][0], limits[var][1], value=initial[var]))

    if exp_time == 0:

        @pymc.deterministic
        def mcmc_f(model_variables=variables):

            detrend_zero = np.array([model_variables[3 * xx] for xx in range(sets_n)])
            detrend_zero = detrend_zero[datai]
            detrend_one = np.array([model_variables[3 * xx + 1] for xx in range(sets_n)])
            detrend_one = detrend_one[datai]
            detrend_two = np.array([model_variables[3 * xx + 2] for xx in range(sets_n)])
            detrend_two = detrend_two[datai]

            detrend = detrend_zero * (1 + detrend_one * datadt + detrend_two * datadt * datadt)
            transit_model = transit(limb_darkening_coefficients, *model_variables[3 * sets_n:],
                                    time_array=datax, method=method, precision=precision)

            return detrend * transit_model

    else:
        datax_hr = np.array([])
        for i in range(time_factor):
            datax_hr = np.append(datax_hr, datax - exp_time / 2.0 + (i + 0.5) * exp_time / time_factor)

        @pymc.deterministic
        def mcmc_f(model_variables=variables):

            detrend_zero = np.array([model_variables[3 * xx] for xx in range(sets_n)])
            detrend_zero = detrend_zero[datai]
            detrend_one = np.array([model_variables[3 * xx + 1] for xx in range(sets_n)])
            detrend_one = detrend_one[datai]
            detrend_two = np.array([model_variables[3 * xx + 2] for xx in range(sets_n)])
            detrend_two = detrend_two[datai]

            detrend = detrend_zero * (1 + detrend_one * datadt + detrend_two * datadt * datadt)
            transit_model = transit(limb_darkening_coefficients, *model_variables[3 * sets_n:],
                                    time_array=datax_hr, method=method, precision=precision)
            transit_model = np.mean(np.reshape(transit_model, (time_factor, len(datax))), 0)

            return detrend * transit_model

    y = pymc.Normal('y', mu=mcmc_f, tau=dataz ** (-2), observed=True, value=datay)
    mcmc = pymc.MCMC([variables, mcmc_f, y], db='pickle', dbname='model.pickle')
    mcmc.sample(iterations, burn=burn)

    for var in range(len(names)):
        if limits[var] is None:
            traces.append(np.array([initial[var]]))
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
            results.append(initial[var])
            errors.append(0)

    for set_n in range(sets_n):

        datasetx = datax[np.where(datai == set_n)]
        datasety = datay[np.where(datai == set_n)]

        def final_model(time_array):

            if exp_time == 0:
                time_array_hr = time_array
            else:
                time_array_hr = np.array([])
                for ii in range(time_factor):
                    time_array_hr = \
                        np.append(time_array_hr, time_array - exp_time / 2.0 + (ii + 0.5) * exp_time / time_factor)

            transit_model = transit(limb_darkening_coefficients, *results[3 * sets_n:],
                                    time_array=time_array_hr, method=method, precision=precision)
            transit_model = np.mean(np.reshape(transit_model, (time_factor, len(time_array))), 0)
            return transit_model

        def final_systematics_model(time_array):
            d_time = time_array - time_array[0]

            return results[3 * set_n] \
                * (1 + results[3 * set_n + 1] * d_time + results[3 * set_n + 2] * d_time * d_time)

        pylightcurve_tools.save_model(datasetx, datasety, final_model, final_systematics_model, set_n)
        pylightcurve_tools.plot_model(datasetx, datasety, names, initial, results, errors, final_model,
                                      final_systematics_model, set_n)

    print 'Saving results...'
    pylightcurve_tools.save_results(names, initial, results, errors, limb_darkening_coefficients)
    pylightcurve_tools.save_traces(names, traces)

    print 'Plotting...'
    pylightcurve_tools.plot_traces(names, traces, results, errors)
    pylightcurve_tools.plot_correlations(names, traces, results, errors)

    os.system('rm model.pickle')
    os.chdir(initial_directory)
