__all__ = ['transit', 'eclipse', 'mcmc_transit']

import pylightcurve_tools
import pylightcurve_quad
import numpy as np
import os
import pymc


class PYLCError(BaseException):
    pass


class PYLCMcmcError(PYLCError):
    pass


def transit(limb_darkening_coefficients, rp_over_rs,
            period, sma_over_rs, eccentricity, inclination, periastron, mid_time, time_array):

    position_vector = pylightcurve_tools.position_vector(period, sma_over_rs,
                                                         eccentricity, inclination, periastron, mid_time, time_array)

    projected_distance = np.where(
        position_vector[0] < 0, 1.0 + 5.0 * rp_over_rs,
        np.sqrt(position_vector[1] * position_vector[1] + position_vector[2] * position_vector[2]))

    return pylightcurve_tools.flux_drop(limb_darkening_coefficients, rp_over_rs, projected_distance)


def transit_quad(limb_darkening_coefficients, rp_over_rs,
                 period, sma_over_rs, eccentricity, inclination, periastron, mid_time, time_array):

    position_vector = pylightcurve_tools.position_vector(period, sma_over_rs,
                                                         eccentricity, inclination, periastron, mid_time, time_array)

    projected_distance = np.where(
        position_vector[0] < 0, 1.0 + 5.0 * rp_over_rs,
        np.sqrt(position_vector[1] * position_vector[1] + position_vector[2] * position_vector[2]))

    return pylightcurve_quad.flux_drop(limb_darkening_coefficients, rp_over_rs, projected_distance)


def eclipse(fp_over_fs, rp_over_rs, period, sma_over_rs, eccentricity, inclination, periastron, mid_time, time_array):

    position_vector = pylightcurve_tools.position_vector(period, -sma_over_rs / rp_over_rs,
                                                         eccentricity, inclination, periastron, mid_time, time_array)

    projected_distance = np.where(
        position_vector[0] < 0, 1.0 + 5.0 / rp_over_rs,
        np.sqrt(position_vector[1] * position_vector[1] + position_vector[2] * position_vector[2]))

    return (1.0 + fp_over_fs * pylightcurve_tools.flux_drop((0, 0, 0, 0), 1 / rp_over_rs, projected_distance)) \
        / (1.0 + fp_over_fs)


def mcmc_transit(limb_darkening_coefficients, rp_over_rs,
                 period, sma_over_rs, eccentricity, inclination, periastron, mid_time,
                 data, fit_rp_over_rs, iterations, burn, directory, detrend_order=0,
                 fit_period=None, fit_sma_over_rs=None, fit_eccentricity=None,
                 fit_inclination=None, fit_periastron=None, fit_mid_time=None):

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

    for set_n in range(sets_n):
        dataset = data[set_n]

        if len(dataset) == 3:
            datasetx, datasety, datasetz = dataset
        else:
            datasetx, datasety = dataset
            datasetz = np.ones_like(datasetx) * \
                np.sqrt(np.mean((datasety[1:-1] - 0.5 * (datasety[:-2] + datasety[2:])) ** 2))

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

    @pymc.deterministic
    def mcmc_f(model_variables=variables):

        detrend_zero = np.array([model_variables[3 * xx] for xx in range(sets_n)])
        detrend_zero = detrend_zero[datai]
        detrend_one = np.array([model_variables[3 * xx + 1] for xx in range(sets_n)])
        detrend_one = detrend_one[datai]
        detrend_two = np.array([model_variables[3 * xx + 2] for xx in range(sets_n)])
        detrend_two = detrend_two[datai]

        detrend = detrend_zero * (1 + detrend_one * datadt + detrend_two * datadt * datadt)

        return detrend * transit(limb_darkening_coefficients, *model_variables[3 * sets_n:], time_array=datax)

    y = pymc.Normal('y', mu=mcmc_f, tau=dataz ** (-2), observed=True, value=datay)
    mcmc = pymc.MCMC([variables, mcmc_f, y], db='pickle', dbname='model.pickle')
    mcmc.isample(iterations, burn=burn, verbose=1)

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
            d_time = time_array - time_array[0]
            return results[3 * set_n] \
                * (1 + results[3 * set_n + 1] * d_time + results[3 * set_n + 2] * d_time * d_time) \
                * transit(limb_darkening_coefficients, *results[3 * sets_n:], time_array=time_array)

        pylightcurve_tools.save_model(datasetx, datasety, final_model, set_n)
        pylightcurve_tools.plot_model(datasetx, datasety, final_model, set_n)

    print 'Saving results...'
    pylightcurve_tools.save_results(names, results, errors)
    pylightcurve_tools.save_traces(names, traces)

    print 'Plotting...'
    pylightcurve_tools.plot_traces(names, traces, results, errors)
    pylightcurve_tools.plot_correlations(names, traces, results, errors)

    os.system('rm model.pickle')
    os.chdir(initial_directory)


def mcmc_transit_quad(limb_darkening_coefficients, rp_over_rs,
                      period, sma_over_rs, eccentricity, inclination, periastron, mid_time,
                      data, fit_rp_over_rs, iterations, burn, directory, detrend_order=0,
                      fit_period=None, fit_sma_over_rs=None, fit_eccentricity=None,
                      fit_inclination=None, fit_periastron=None, fit_mid_time=None):

    if len(limb_darkening_coefficients) != 2:
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

    for set_n in range(sets_n):
        dataset = data[set_n]

        if len(dataset) == 3:
            datasetx, datasety, datasetz = dataset
        else:
            datasetx, datasety = dataset
            datasetz = np.ones_like(datasetx) * \
                np.sqrt(np.mean((datasety[1:-1] - 0.5 * (datasety[:-2] + datasety[2:])) ** 2))

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

    @pymc.deterministic
    def mcmc_f(model_variables=variables):

        detrend_zero = np.array([model_variables[3 * xx] for xx in range(sets_n)])
        detrend_zero = detrend_zero[datai]
        detrend_one = np.array([model_variables[3 * xx + 1] for xx in range(sets_n)])
        detrend_one = detrend_one[datai]
        detrend_two = np.array([model_variables[3 * xx + 2] for xx in range(sets_n)])
        detrend_two = detrend_two[datai]

        detrend = detrend_zero * (1 + detrend_one * datadt + detrend_two * datadt * datadt)

        return detrend * transit_quad(limb_darkening_coefficients, *model_variables[3 * sets_n:], time_array=datax)

    y = pymc.Normal('y', mu=mcmc_f, tau=dataz ** (-2), observed=True, value=datay)
    mcmc = pymc.MCMC([variables, mcmc_f, y], db='pickle', dbname='model.pickle')
    mcmc.isample(iterations, burn=burn, verbose=1)

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
            d_time = time_array - time_array[0]
            return results[3 * set_n] \
                * (1 + results[3 * set_n + 1] * d_time + results[3 * set_n + 2] * d_time * d_time) \
                * transit(limb_darkening_coefficients, *results[3 * sets_n:], time_array=time_array)

        pylightcurve_tools.save_model(datasetx, datasety, final_model, set_n)
        pylightcurve_tools.plot_model(datasetx, datasety, final_model, set_n)

    print 'Saving results...'
    pylightcurve_tools.save_results(names, results, errors)
    pylightcurve_tools.save_traces(names, traces)

    print 'Plotting...'
    pylightcurve_tools.plot_traces(names, traces, results, errors)
    pylightcurve_tools.plot_correlations(names, traces, results, errors)

    os.system('rm model.pickle')
    os.chdir(initial_directory)
