
__all__ = ['Fitting', 'EmceeFitting', 'curve_fit']

import emcee
import sys
import numpy as np
import warnings

from scipy.stats import shapiro
from scipy.optimize import curve_fit as scipy_curve_fit
from scipy.optimize import minimize

from .distributions import one_d_distribution
from .gaussian import gaussian
from .statistics import values_to_print, residual_statistics

from ..errors import *
from ..processes.counter import Counter
from ..processes.files import save_dict
from ..plots.plots_fitting import plot_mcmc_corner, plot_mcmc_traces, plot_mcmc_fitting


def curve_fit(*args, **kwargs):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message='Covariance of the parameters could not be estimated')
        return scipy_curve_fit(*args, **kwargs)


class Fitting:

    def __init__(self, input_data_x, input_data_y, input_data_y_unc,
                 model, initials, limits1, limits2,
                 logspace=None,
                 data_x_name='x', data_y_name='y', data_x_print_name='x', data_y_print_name='y',
                 walkers=None, iterations=None, burn_in=None, strech_prior=0.2, walkers_spread=0.01,
                 parameters_names=None, parameters_print_names=None,
                 counter=None,
                 optimise_initial_parameters=True,
                 optimise_initial_parameters_trials=4,
                 scale_uncertainties=False,
                 filter_outliers=False,
                 optimiser='emcee',
                 ):

        self.input_data_x = np.array(input_data_x, dtype=float)
        self.input_data_y = np.array(input_data_y, dtype=float)
        self.input_data_y_unc = np.array(input_data_y_unc, dtype=float)
        self.input_data_y_unc_backup = np.array(input_data_y_unc, dtype=float)

        self.data_x_name = str(data_x_name)
        self.data_y_name = str(data_y_name)
        self.data_x_print_name = str(data_x_print_name)
        self.data_y_print_name = str(data_y_print_name)

        self.model = model
        self.initials = np.array(initials, dtype=float)
        self.limits1 = np.array(limits1, dtype=float)
        self.limits2 = np.array(limits2, dtype=float)
        self.fitted_parameters_indices = np.where(~np.isnan(self.limits1 * self.limits2))[0]
        self.dimensions = len(self.fitted_parameters_indices)

        self.names = ['p{0}'.format(ff) for ff in range(len(initials))]
        if parameters_names is not None:
            self.names = parameters_names

        self.print_names = ['p{0}'.format(ff) for ff in range(len(initials))]
        if parameters_print_names is not None:
            self.print_names = parameters_print_names

        self.logspace = np.zeros_like(initials, dtype=bool)
        if logspace is not None:
            self.logspace = np.array(logspace, dtype=bool)

        for var in range(len(self.initials)):
            if np.isnan(self.limits1[var] * self.limits2[var]):
                self.logspace[var] = False

        self.internal_limits1 = []
        self.internal_limits2 = []
        self.internal_initials = []
        for var in self.fitted_parameters_indices:

            if self.initials[var] < self.limits1[var] or self.initials[var] > self.limits2[var]:
                raise PyLCInputError('Initial value for parameter {0} is outside the prior bounds.'.format(
                    self.names[var]))

            if self.logspace[var]:
                self.internal_limits1.append(np.log10(self.limits1[var]))
                self.internal_limits2.append(np.log10(self.limits2[var]))
                self.internal_initials.append((np.log10(self.initials[var]) - np.log10(self.limits1[var])) / (np.log10(self.limits2[var]) - np.log10(self.limits1[var])))
            else:
                self.internal_limits1.append(self.limits1[var])
                self.internal_limits2.append(self.limits2[var])
                self.internal_initials.append((self.initials[var] - self.limits1[var]) / (self.limits2[var] - self.limits1[var]))
        self.internal_limits1 = np.array(self.internal_limits1)
        self.internal_limits2 = np.array(self.internal_limits2)
        self.internal_initials = np.array(self.internal_initials)

        if optimiser not in ['emcee', 'curve_fit']:
            raise PyLCInputError('Optimiser {0} in not valid. Please choose between '
                             'emcee, scipy_minimize.'.format(optimiser))
        self.optimiser = optimiser

        if not walkers:
            self.walkers = 3 * self.dimensions
        else:
            self.walkers = int(walkers)

        if not iterations:
            self.iterations = 5000
        else:
            self.iterations = int(iterations)

        if burn_in is None:
            self.burn_in = int(self.iterations * 0.2)
        else:
            self.burn_in = int(burn_in)
            if self.burn_in >= self.iterations:
                raise PyLCInputError('burn_in must be lower than iterations.')

        self.strech_prior = float(strech_prior)
        self.walkers_spread = float(walkers_spread)
        self.walkers_initial_positions = None
        self.sampler = None
        self.progress = 0

        if counter:
            if not isinstance(counter, str):
                raise PyLCInputError('Counter should be a string.')
        self.counter_name = counter

        self.scale_uncertainties = bool(scale_uncertainties)
        self.filter_outliers = bool(filter_outliers)
        self.optimise_initial_parameters = bool(optimise_initial_parameters)
        self.optimise_initial_parameters_trials = int(optimise_initial_parameters_trials)

        self.results = {
            'original_data': {},
            'settings': {},
            'prefit': {},
            'input_series': {},
            'parameters': {},
            'parameters_final': [],
            'output_series': {},
            'statistics': {}}

        self.fitted_parameters = []

        self.mcmc_run_complete = False
        self.prefit_complete = False

        self.results['original_data'][self.data_x_name] = np.ones_like(self.input_data_x) * self.input_data_x
        self.results['original_data'][self.data_y_name] = np.ones_like(self.input_data_y) * self.input_data_y
        self.results['original_data'][self.data_y_name + '_unc'] = np.ones_like(self.input_data_y_unc) * self.input_data_y_unc
        self.results['settings']['walkers'] = self.walkers
        self.results['settings']['iterations'] = self.iterations
        self.results['settings']['burn_in'] = self.burn_in
        self.results['settings']['data_x_name'] = self.data_x_name
        self.results['settings']['data_x_print_name'] = self.data_x_print_name
        self.results['settings']['data_y_name'] = self.data_y_name
        self.results['settings']['data_y_print_name'] = self.data_y_print_name
        self.results['settings']['strech_prior'] = self.strech_prior
        self.results['settings']['walkers_spread'] = self.walkers_spread
        self.results['settings']['optimise_initial_parameters'] = self.optimise_initial_parameters
        self.results['settings']['scale_uncertainties'] = self.scale_uncertainties
        self.results['settings']['filter_outliers'] = self.filter_outliers
        self.results['settings']['optimiser'] = self.optimiser

    def _pass(self):
        pass

    def _internal_model(self, theta):
        parameters = np.ones_like(self.initials) * self.initials
        parameters[self.fitted_parameters_indices] = theta * (self.internal_limits2 - self.internal_limits1) + self.internal_limits1
        parameters[self.logspace] = 10**parameters[self.logspace]
        return self.model(self.input_data_x, *parameters)

    def _internal_model_curve_fit(self, x, *theta):
        theta = np.array(theta)
        parameters = np.ones_like(self.initials) * self.initials
        parameters[self.fitted_parameters_indices] = theta * (self.internal_limits2 - self.internal_limits1) + self.internal_limits1
        parameters[self.logspace] = 10**parameters[self.logspace]
        return self.model(self.input_data_x, *parameters)

    def _probability(self, theta):
        if np.prod((0 < theta) * (theta < 1)):
            chi = (self.input_data_y - self._internal_model(theta)) / self.input_data_y_unc
            return -0.5 * (np.sum(chi * chi) +
                           np.sum(np.log(2.0 * np.pi * (self.input_data_y_unc * self.input_data_y_unc))))
        else:
            return -np.inf

    def _prefit(self, verbose=False):

        if self.scale_uncertainties or self.filter_outliers or self.optimise_initial_parameters:

            nll = lambda *args: -self._probability(*args)
            soln_test = nll(self.internal_initials)

            test_initials = np.ones_like(self.internal_initials) * self.internal_initials
            optimisation_ok = False

            # import matplotlib.pyplot as plt
            # plt.figure()
            # plt.plot(self.input_data_x, self.input_data_y, 'ko')
            # plt.plot(self.input_data_x, self._internal_model(self.internal_initials), 'r-')
            # plt.show()

            for ii in range(self.optimise_initial_parameters_trials):

                if verbose:
                    print('Optimising initial parameters, attempt {0}: maximizing likelihood...'.format(ii+1))

                soln_i = minimize(nll, test_initials, method='Nelder-Mead')
                soln_test_i = nll(soln_i.x)
                optimisation_ok = soln_i.success

                if soln_i.success and soln_test_i <= soln_test:

                    if verbose:
                        print('Optimisation completed.')

                    self.internal_initials = soln_i.x

                    if self.filter_outliers:
                        norm_res = (self.input_data_y - self._internal_model(self.internal_initials)) / self.input_data_y_unc
                        outliers = len(np.where(np.abs(norm_res) >= 3 * np.std(norm_res))[0])

                        while outliers > 0:

                            if verbose:
                                print('Filtering outliers...'.format(ii+1))

                            flags = np.where(np.abs(norm_res) >= 3 * np.std(norm_res))[0]
                            self.input_data_y_unc[flags] = 1000000000

                            soln_i = minimize(nll, self.internal_initials, method='Nelder-Mead')
                            self.internal_initials = soln_i.x
                            norm_res = (self.input_data_y - self._internal_model(self.internal_initials)) / self.input_data_y_unc
                            # print(np.std(norm_res), np.median(np.abs(norm_res - np.median(norm_res))))

                            outliers = len(np.where(np.abs(norm_res) >= 3 * np.std(norm_res))[0])
                            # print(outliers)

                    optimisation_ok = True
                    break

                elif soln_test_i <= soln_test:
                    test_initials = soln_i.x
                    soln_test = nll(test_initials)
                else:
                    test_initials = self.internal_initials + np.random.normal(0, self.strech_prior/2.0, len(self.internal_initials))
                    test_initials = np.maximum(0, test_initials)
                    test_initials = np.minimum(1, test_initials)

            if not optimisation_ok:
                raise PyLCProcessError('Optimisation failed. You can try re-running, or increasing the strech_prior parameter, or the prior limits.')

        # import matplotlib.pyplot as plt
        # plt.figure()
        # plt.plot(self.input_data_x, self.input_data_y, 'ko')
        # plt.plot(self.input_data_x, self._internal_model(self.internal_initials), 'r-')
        # plt.show()

        self.results['prefit']['outliers_map'] = self.input_data_y_unc == 1000000000
        outliers = np.where(self.results['prefit']['outliers_map'])
        points_to_use = np.where(~self.results['prefit']['outliers_map'])
        self.results['prefit']['outliers'] = len(outliers[0])
        self.input_data_x = self.input_data_x[points_to_use]
        self.input_data_y = self.input_data_y[points_to_use]
        self.input_data_y_unc = self.input_data_y_unc[points_to_use]

        if self.scale_uncertainties:
            scale_factor = np.sqrt(np.nanmean(((self.input_data_y - self._internal_model(self.internal_initials)) ** 2) / (self.input_data_y_unc ** 2)))
            # import matplotlib.pyplot as plt
            # plt.figure()
            # plt.plot(self.input_data_x, self.input_data_y - self._internal_model(self.internal_initials), 'ko')
            # plt.plot(self.input_data_x, self.input_data_y, 'ko')
            # plt.plot(self.input_data_x, self._internal_model(self.internal_initials), 'r-')
            # plt.show()
        else:
            scale_factor = 1

        self.input_data_y_unc *= scale_factor
        self.results['prefit']['scale_factor'] = scale_factor

        if verbose:
            print('Data-points excluded:', self.results['prefit']['outliers'])
            print('Scaling uncertainties by:', scale_factor)
            print('Initial parameters:')

        self.results['prefit']['initials'] = []
        for var in range(len(self.names)):
            if np.isnan(self.limits1[var]):
                self.results['prefit']['initials'].append(self.initials[var])
            else:
                idx = np.where(self.fitted_parameters_indices == var)[0][0]

                self.initials[var] = (self.internal_initials * (self.internal_limits2 - self.internal_limits1) + self.internal_limits1)[idx]
                if self.logspace[var]:
                    self.initials[var] = 10 ** self.initials[var]

                if verbose:
                    print(self.names[var], ': ', self.initials[var])

                self.results['prefit']['initials'].append(self.initials[var])

        self.results['prefit']['initials'] = np.array(self.results['prefit']['initials'])
        self.results['input_series'][self.data_x_name] = self.input_data_x
        self.results['input_series'][self.data_y_name] = self.input_data_y
        self.results['input_series'][self.data_y_name + '_unc'] = self.input_data_y_unc

        self.prefit_complete = True

    def run(self, verbose=False):

        if not self.prefit_complete:
            self._prefit(verbose=verbose)

        # run sampler

        if self.optimiser == 'curve_fit':

            option = int(bool(self.counter_name))

            print(['Running curve_fit...', self.counter_name][option] + '...')

            popt, pcov = curve_fit(self._internal_model_curve_fit, self.input_data_x,
                                   self.input_data_y, sigma=self.input_data_y_unc, p0=self.internal_initials,
                                   maxfev=self.iterations
                                   )

            for var in range(len(self.names)):

                if np.isnan(self.limits1[var]):

                    variable = {'name': self.names[var], 'print_name': self.print_names[var],
                                'initial': None, 'min_allowed': None, 'max_allowed': None,
                                'trace': None, 'trace_bins': None, 'trace_counts': None,
                                'value': self.initials[var], 'm_error': None, 'p_error': None,
                                'print_value': self.initials[var], 'print_m_error': '-', 'print_p_error': '-'}

                else:

                    idx = np.where(self.fitted_parameters_indices == var)[0][0]

                    value = popt[idx]
                    min_value = value - np.sqrt(pcov[idx][idx])
                    max_value = value + np.sqrt(pcov[idx][idx])

                    value = value * (self.internal_limits2 - self.internal_limits1)[idx] + self.internal_limits1[idx]
                    min_value = min_value * (self.internal_limits2 - self.internal_limits1)[idx] + self.internal_limits1[idx]
                    max_value = max_value * (self.internal_limits2 - self.internal_limits1)[idx] + self.internal_limits1[idx]

                    if self.logspace[var]:
                        value = 10**value
                        min_value = 10**min_value
                        max_value = 10**max_value

                    m_error = value - min_value
                    p_error = max_value - value

                    print_value, print_m_error, print_p_error = values_to_print(value, m_error, p_error)

                    variable = {'name': self.names[var], 'print_name': self.print_names[var],
                                'initial': self.initials[var],
                                'min_allowed': self.limits1[var], 'max_allowed': self.limits2[var],
                                'trace': None, 'trace_bins': None, 'trace_counts': None,
                                'value': value, 'm_error': m_error, 'p_error': p_error,
                                'print_value': print_value, 'print_m_error': print_m_error, 'print_p_error': print_p_error}

                    self.fitted_parameters.append(self.names[var])

                self.results['parameters'][self.names[var]] = variable
                self.results['parameters_final'].append(variable['value'])

            self.results['statistics']['corr_matrix'] = pcov
            self.results['statistics']['corr_variables'] = ','.join(self.fitted_parameters)

            self._postfit()

        elif self.optimiser == 'emcee':

            sys.setrecursionlimit(self.iterations)

            self.sampler = emcee.EnsembleSampler(self.walkers, self.dimensions, self._probability)

            option = int(bool(self.counter_name))

            self.counter = Counter(['MCMC', self.counter_name][option], self.iterations - self.progress,
                                   show_every=[self.iterations, 0][option] + 10, increment=10)
            self._emcee_run()

    def _emcee_run(self):

        if self.progress == 0:
            while self.progress == 0:
                try:
                    self.walkers_initial_positions = np.random.uniform(
                        (self.internal_initials - 0.5 * self.walkers_spread)[:, None] * np.ones(self.walkers),
                        (self.internal_initials + 0.5 * self.walkers_spread)[:, None] * np.ones(self.walkers))
                    self.walkers_initial_positions = np.swapaxes(self.walkers_initial_positions,0, 1)
                    self.walkers_initial_positions = np.minimum(self.walkers_initial_positions, 1)
                    self.walkers_initial_positions = np.maximum(self.walkers_initial_positions, 0)

                    self.sampler.run_mcmc(self.walkers_initial_positions, 10)
                    self.progress += 10
                    self.counter.update()

                    self._emcee_run()

                except ValueError:
                    pass

        elif self.progress < self.iterations:

            self.sampler.run_mcmc(None, 10, skip_initial_state_check=True)
            self.progress += 10
            self.counter.update()

            self._emcee_run()

        else:

            mcmc_results = self.sampler.get_chain()

            trace_to_analyse = 0
            vars_check = 0
            for var in range(len(self.names)):

                if not np.isnan(self.limits1[var]):
                    trace = mcmc_results[self.burn_in:, :, np.where(self.fitted_parameters_indices == var)[0][0]]
                    trace = trace.flatten()

                    median = np.median(trace)
                    mad = np.sqrt(np.median((trace - median) ** 2))

                    trace_to_analyse += (trace > (median - 10 * mad)) * (trace < (median + 10 * mad))

                    vars_check += 1

            trace_to_analyse = np.where(trace_to_analyse == vars_check)

            for var in range(len(self.names)):

                if np.isnan(self.limits1[var]):

                    variable = {'name': self.names[var], 'print_name': self.print_names[var],
                                'initial': None, 'min_allowed': None, 'max_allowed': None,
                                'trace': None, 'trace_bins': None, 'trace_counts': None,
                                'value': self.initials[var], 'm_error': None, 'p_error': None,
                                'print_value': str(self.initials[var]), 'print_m_error': '-', 'print_p_error': '-'}

                else:
                    idx = np.where(self.fitted_parameters_indices == var)[0][0]

                    trace = mcmc_results[self.burn_in:, :, idx]
                    trace = trace.flatten()
                    trace = trace[trace_to_analyse]

                    bins, counts = one_d_distribution(trace)

                    min_value, value, max_value = np.quantile(trace, [0.16, 0.5, 0.84])

                    trace = trace * (self.internal_limits2 - self.internal_limits1)[idx] + self.internal_limits1[idx]
                    bins = bins * (self.internal_limits2 - self.internal_limits1)[idx] + self.internal_limits1[idx]
                    value = value * (self.internal_limits2 - self.internal_limits1)[idx] + self.internal_limits1[idx]
                    min_value = min_value * (self.internal_limits2 - self.internal_limits1)[idx] + self.internal_limits1[idx]
                    max_value = max_value * (self.internal_limits2 - self.internal_limits1)[idx] + self.internal_limits1[idx]

                    if self.logspace[var]:
                        trace = 10**trace
                        bins = 10**bins
                        value = 10**value
                        min_value = 10**min_value
                        max_value = 10**max_value

                    m_error = value - min_value
                    p_error = max_value - value

                    print_value, print_m_error, print_p_error = values_to_print(value, m_error, p_error)

                    variable = {'name': self.names[var], 'print_name': self.print_names[var],
                                'initial': self.initials[var],
                                'min_allowed': self.limits1[var], 'max_allowed': self.limits2[var],
                                'trace': trace, 'trace_bins': bins, 'trace_counts': counts,
                                'value': value, 'm_error': m_error, 'p_error': p_error,
                                'print_value': print_value, 'print_m_error': print_m_error, 'print_p_error': print_p_error}

                    self.fitted_parameters.append(self.names[var])

                self.results['parameters'][self.names[var]] = variable
                self.results['parameters_final'].append(variable['value'])

            to_correlate = []
            for parameter in self.fitted_parameters:
                to_correlate.append(self.results['parameters'][parameter]['trace'])
            correlation_matrix = np.corrcoef(to_correlate)
            self.results['statistics']['corr_matrix'] = correlation_matrix
            self.results['statistics']['corr_variables'] = ','.join(self.fitted_parameters)

            self._postfit()

    def _postfit(self):

        self.results['parameters_final'] = np.array(self.results['parameters_final'])

        self.results['output_series']['model'] = self.model(self.input_data_x, *self.results['parameters_final'])
        self.results['output_series']['residuals'] = self.input_data_y - self.results['output_series']['model']

        statistics = residual_statistics(self.input_data_x, self.input_data_y, self.input_data_y_unc,
                                          self.results['output_series']['model'], len(self.fitted_parameters))
        for statistic in statistics:
            self.results['statistics'][statistic] = statistics[statistic]

        self.mcmc_run_complete = True

    def save_all(self, export_file):

        if not self.mcmc_run_complete:
            raise PyLCProcessError('MCMC not completed')

        save_dict(self.results, export_file)

    def save_results(self, export_file):

        if not self.mcmc_run_complete:
            raise PyLCProcessError('MCMC not completed')

        cols = [
            ['# variable'],
            ['fix/fit'],
            ['value'],
            ['uncertainty'],
            ['initial'],
            ['min.allowed'],
            ['max.allowed']
        ]

        for i in self.names:

            cols[0].append(self.results['parameters'][i]['name'])
            if self.results['parameters'][i]['initial'] is None:
                cols[1].append('fix')
                cols[2].append(self.results['parameters'][i]['print_value'])
                cols[3].append(' ')
                cols[4].append(' ')
                cols[5].append(' ')
                cols[6].append(' ')
            else:
                cols[1].append('fit')
                cols[2].append(self.results['parameters'][i]['print_value'])
                cols[3].append(
                    '-{0} +{1}'.format(
                        self.results['parameters'][i]['print_m_error'],
                        self.results['parameters'][i]['print_p_error'])
                )
                cols[4].append(str(self.results['parameters'][i]['initial']))
                cols[5].append(str(self.results['parameters'][i]['min_allowed']))
                cols[6].append(str(self.results['parameters'][i]['max_allowed']))

        for col in cols:
            col_length = np.max([len(ff) for ff in col])
            for ff in range(len(col)):
                col[ff] = col[ff] + ' ' * (col_length - len(col[ff]))

        lines = []

        for row in range(len(cols[0])):
            lines.append('  '.join([col[row] for col in cols]))

        try:
            _ = self.results['prefit']
            lines.append('')
            lines.append('#Pre-fit:')
            lines.append('#Number of outliers removed: {0}'.format(self.results['prefit']['outliers']))
            lines.append('#Uncertainties scale factor: {0}'.format(self.results['prefit']['scale_factor']))
        except:
            pass

        lines.append('')
        lines.append('#Residuals:')
        lines.append('#Mean: {0}'.format(self.results['statistics']['res_mean']))
        lines.append('#STD: {0}'.format(self.results['statistics']['res_std']))
        lines.append('#RMS: {0}'.format(self.results['statistics']['res_rms']))
        lines.append('#Chi squared: {0}'.format(self.results['statistics']['res_chi_sqr']))
        lines.append('#Reduced chi squared: {0}'.format(self.results['statistics']['res_red_chi_sqr']))
        lines.append('#Max auto-correlation: {0}'.format(self.results['statistics']['res_max_autocorr']))
        lines.append('#Max auto-correlation flag: {0}'.format(self.results['statistics']['res_max_autocorr_flag']))
        lines.append('#Shapiro test: {0}'.format(self.results['statistics']['res_shapiro']))
        lines.append('#Shapiro test flag: {0}'.format(self.results['statistics']['res_shapiro_flag']))

        w = open(export_file, 'w')
        w.write('\n'.join(lines))
        w.close()

    def plot_fitting(self, export_file):

        plot_mcmc_fitting(self, export_file)

    def plot_corner(self, export_file):

        plot_mcmc_corner(self, export_file)

    def plot_traces(self, export_file):

        plot_mcmc_traces(self, export_file)



# for compatibility reasons
class EmceeFitting(Fitting):

    def __init__(self, *args, **kwargs):

        try:
            kwargs['walkers'] = args[7]
            kwargs['iterations'] = args[8]
            kwargs['burn_in'] = args[9]
            if args[8] > 10000:
                print('In version 4.1 of PyLightcurve, the total number of calculations is given by iterations x walkers.')
                kwargs['iterations'] /= kwargs['walkers']
                kwargs['burn_in'] /= kwargs['walkers']
        except:
            kwargs['walkers'] = None
            kwargs['iterations'] = None
            kwargs['burn_in'] = None

        print('Warning plc.EmceeFitting will not be supported in future versions.')
        print('Use plc.Fitting instead, and set optimiser="emcee".')

        Fitting.__init__(self, *args[:7], **kwargs)

        self.run_mcmc = self.run
