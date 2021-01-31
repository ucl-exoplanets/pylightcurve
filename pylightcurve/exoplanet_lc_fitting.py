from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from .exoplanet_orbit import *
from .exoplanet_lc import *
from .analysis_emcee_fitting import *


class TransitAndPolyFitting:

    def __init__(self, data,
                 method, limb_darkening_coefficients,
                 rp_over_rs, period, sma_over_rs, eccentricity, inclination, periastron, mid_time,
                 iterations, walkers, burn, precision=3,
                 exp_time=0, time_factor=1, fit_first_order=False, fit_second_order=False,
                 fit_rp_over_rs=False, fit_period=False, fit_sma_over_rs=False,
                 fit_eccentricity=False, fit_inclination=False, fit_periastron=False, fit_mid_time=False,
                 counter=True, counter_window=False):

        # TODO check input parameters

        self.data = data
        self.total_sets = len(self.data)

        self.method = method
        self.limb_darkening_coefficients = limb_darkening_coefficients
        self.fit_ld = False

        if isinstance(self.limb_darkening_coefficients, str):
            if self.limb_darkening_coefficients == 'fit':
                self.fit_ld = True
                if self.method == 'linear':
                    self.limb_darkening_coefficients = [0.5]
                elif self.method in ['quad', 'sqrt']:
                    self.limb_darkening_coefficients = [0.5, 0.5]
                elif self.method == 'claret':
                    self.limb_darkening_coefficients = [0.5, 0.5, 0.5, 0.5]

        if self.method == 'claret':
            self.total_ldcs = 4
        elif self.method in ['quad', 'sqrt']:
            self.total_ldcs = 2
        elif self.method == 'linear':
            self.total_ldcs = 1

        self.rp_over_rs = rp_over_rs
        self.period = period
        self.sma_over_rs = sma_over_rs
        self.eccentricity = eccentricity
        self.inclination = inclination
        self.periastron = periastron
        self.mid_time = mid_time

        self.precision = precision
        self.exp_time = exp_time / (60.0 * 60.0 * 24.0)
        self.time_factor = time_factor

        self.fit_first_order = fit_first_order
        self.fit_second_order = fit_second_order

        self.fit_rp_over_rs = fit_rp_over_rs
        self.fit_period = fit_period
        self.fit_sma_over_rs = fit_sma_over_rs
        self.fit_eccentricity = fit_eccentricity
        self.fit_inclination = fit_inclination
        self.fit_periastron = fit_periastron
        self.fit_mid_time = fit_mid_time

        self.iterations = iterations
        self.walkers = walkers
        self.burn = burn
        self.counter = counter
        self.counter_window = counter_window

        self.data_time = np.array([])
        self.data_flux = np.array([])
        self.data_flux_error = np.array([])
        self.data_set_number = np.array([])
        self.data_set_dt = np.array([])

        for set_number, set_arrays in enumerate(self.data):

            if set_number == 0:
                time_shift = round((np.mean(set_arrays[0]) - self.mid_time) / self.period)
                self.mid_time += time_shift * self.period
                if self.fit_mid_time:
                    self.fit_mid_time[0] += time_shift * self.period
                    self.fit_mid_time[1] += time_shift * self.period

            self.data_time = np.append(self.data_time, set_arrays[0])
            self.data_flux = np.append(self.data_flux, set_arrays[1])
            self.data_flux_error = np.append(self.data_flux_error, set_arrays[2])
            self.data_set_number = np.append(self.data_set_number, np.ones_like(set_arrays[0]) * set_number)
            self.data_set_dt = np.append(self.data_set_dt, set_arrays[0] - set_arrays[0][0])

        self.data_set_number = np.int_(self.data_set_number)

        self.names = []
        self.print_names = []
        self.limits1 = []
        self.limits2 = []
        self.initial = []

        for set_number, set_arrays in enumerate(self.data):

            max_limit = (10 * (max(set_arrays[1]) - min(set_arrays[1])) / (max(set_arrays[0]) - min(set_arrays[0])) /
                         np.mean(set_arrays[1]))

            self.names.append('N{0}'.format(str(set_number)))
            self.print_names.append('N_{0}'.format(str(set_number)))
            self.initial.append(np.mean(set_arrays[1]))
            self.limits1.append(0.5 * np.mean(set_arrays[1]))
            self.limits2.append(1.5 * np.mean(set_arrays[1]))

            self.names.append('L{0}'.format(str(set_number)))
            self.print_names.append('L_{0}'.format(str(set_number)))
            self.initial.append(0)
            if self.fit_first_order:
                self.limits1.append(-max_limit)
                self.limits2.append(max_limit)
            else:
                self.limits1.append(np.nan)
                self.limits2.append(np.nan)

            self.names.append('Q{0}'.format(str(set_number)))
            self.print_names.append('Q_{0}'.format(str(set_number)))
            self.initial.append(0)
            if self.fit_second_order:
                self.limits1.append(-max_limit)
                self.limits2.append(max_limit)
            else:
                self.limits1.append(np.nan)
                self.limits2.append(np.nan)

        self.names.append('ldc1')
        self.print_names.append('ldc_1')
        self.initial.append(self.limb_darkening_coefficients[0])
        if self.fit_ld:
            self.limits1.append(0.000001)
            self.limits2.append(0.999999)
        else:
            self.limits1.append(np.nan)
            self.limits2.append(np.nan)

        if self.method in ['claret', 'quad', 'sqrt']:

            self.names.append('ldc2')
            self.print_names.append('ldc_2')
            self.initial.append(self.limb_darkening_coefficients[1])
            if self.fit_ld:
                self.limits1.append(0.000001)
                self.limits2.append(0.999999)
            else:
                self.limits1.append(np.nan)
                self.limits2.append(np.nan)

        if self.method == 'claret':

            self.names.append('ldc3')
            self.print_names.append('ldc_3')
            self.initial.append(self.limb_darkening_coefficients[2])
            if self.fit_ld:
                self.limits1.append(0.000001)
                self.limits2.append(0.999999)
            else:
                self.limits1.append(np.nan)
                self.limits2.append(np.nan)

            self.names.append('ldc4')
            self.print_names.append('ldc_4')
            self.initial.append(self.limb_darkening_coefficients[3])
            if self.fit_ld:
                self.limits1.append(0.000001)
                self.limits2.append(0.999999)
            else:
                self.limits1.append(np.nan)
                self.limits2.append(np.nan)

        self.names += ['rp', 'P', 'a', 'e', 'i', 'w', 'mt']
        self.print_names += ['R_\mathrm{p}/R_*', 'P', 'a/R_*', 'e', 'i', '\omega', 'T_0']
        self.initial += [self.rp_over_rs, self.period, self.sma_over_rs, self.eccentricity,
                         self.inclination, self.periastron, self.mid_time]
        limits = self.limits1 + [self.fit_rp_over_rs, self.fit_period, self.fit_sma_over_rs, self.fit_eccentricity,
                                 self.fit_inclination, self.fit_periastron, self.fit_mid_time]

        for var in range(3 * self.total_sets + self.total_ldcs, len(self.names)):

            try:
                self.initial[var] = float(self.initial[var])
            except:
                raise PyLCInputError('Improper value for {0}'.format(self.names[var]))

            if limits[var] is False:
                self.limits1.append(np.nan)
                self.limits2.append(np.nan)

            elif limits[var] is None:
                self.limits1.append(np.nan)
                self.limits2.append(np.nan)

            else:
                try:
                    if len(np.array(limits[var])) != 2:
                        raise RuntimeError('Improper limits for {0}'.format(self.names[var]))
                except:
                    raise RuntimeError('Improper limits for {0}'.format(self.names[var]))

                if self.initial[var] < np.array(limits[var])[0] or self.initial[var] > np.array(limits[var])[1]:
                    raise RuntimeError('Initial value for {0} is outside the range of the prior.'.format(
                        self.names[var]))
                else:
                    self.limits1.append(np.array(limits[var])[0])
                    self.limits2.append(np.array(limits[var])[1])

        if self.exp_time == 0:

            self.data_time_hr = self.data_time

        else:
            self.data_time_hr = np.array([])
            for i in range(self.time_factor):
                self.data_time_hr = np.append(self.data_time_hr, self.data_time - self.exp_time / 2.0 +
                                              (i + 0.5) * self.exp_time / self.time_factor)

        self.fitting = EmceeFitting([self.data_flux, self.data_flux_error],
                                    self.full_model, self.initial, self.limits1, self.limits2,
                                    self.walkers, self.iterations, self.burn,
                                    names=self.names, print_names=self.print_names,
                                    counter=self.counter, counter_window=self.counter_window)

        self.results = 0
        self.mcmc_run_complete = False

    def detrend_model(self, *model_variables):

        detrend_zero = np.array([model_variables[3 * xx] for xx in range(self.total_sets)])
        detrend_zero = detrend_zero[self.data_set_number]
        detrend_one = np.array([model_variables[3 * xx + 1] for xx in range(self.total_sets)])
        detrend_one = detrend_one[self.data_set_number]
        detrend_two = np.array([model_variables[3 * xx + 2] for xx in range(self.total_sets)])
        detrend_two = detrend_two[self.data_set_number]

        return detrend_zero * (1 + detrend_one * self.data_set_dt +
                               detrend_two * self.data_set_dt * self.data_set_dt)

    def transit_model(self, *model_variables):

        limb_darkening_coefficients = list(model_variables)[3 * self.total_sets: 3 * self.total_sets + self.total_ldcs]

        rp_over_rs = list(model_variables)[3 * self.total_sets + self.total_ldcs]

        z_over_rs = transit_projected_distance(*model_variables[3 * self.total_sets + self.total_ldcs + 1:],
                                               time_array=self.data_time_hr)

        transit_hr = transit_flux_drop(self.method, limb_darkening_coefficients, rp_over_rs, z_over_rs,
                                       precision=self.precision)

        return np.mean(np.reshape(transit_hr, (self.time_factor, len(self.data_time))), 0)

    def full_model(self, *model_variables):

        return self.detrend_model(*model_variables) * self.transit_model(*model_variables)

    def run_mcmc(self):

        self.fitting.run_mcmc()
        self.results = self.fitting.results
        self.mcmc_run_complete = True

        self.results['input_series']['hjd'] = self.data_time

        period = self.results['parameters']['P']['value']
        mt = self.results['parameters']['mt']['value']
        self.results['output_series']['phase'] = \
            (self.data_time - mt) / period - np.round((self.data_time - mt) / period)

        self.results['detrended_input_series'] = {
            'hjd': self.results['input_series']['hjd'],
            'value': self.results['input_series']['value'] / self.detrend_model(*self.results['parameters_final']),
            'error': self.results['input_series']['error'] / self.detrend_model(*self.results['parameters_final'])}

        self.results['detrended_output_series'] = {
            'phase': self.results['output_series']['phase'],
            'model': self.results['output_series']['model'] / self.detrend_model(*self.results['parameters_final']),
            'residuals': (self.results['output_series']['residuals']
                          / self.detrend_model(*self.results['parameters_final']))}

        self.results['detrended_statistics'] = {ff: self.results['statistics'][ff] for ff in self.results['statistics']}
        self.results['detrended_statistics']['res_std'] = np.std(self.results['detrended_output_series']['residuals'])

    def save_all(self, export_file):

        pickle.dump(self.results, open(export_file, 'wb'))

    def save_results(self, export_file):

        self.fitting.save_results(export_file)

    def plot_corner(self, export_file):

        self.fitting.plot_corner(export_file)

    def plot_traces(self, export_file):

        self.fitting.plot_traces(export_file)

    def plot_models(self, export_file, target=None, data_dates=None):

        if target is None:
            target = ' '

        if data_dates is None:
            data_dates = map(str, ['set_{0}'.format(str(ff)) for ff in range(1, self.total_sets + 1)])

        for set_number in range(self.total_sets):

            fig = plt.figure(set_number + 1)
            fig.set_tight_layout(False)

            self.results = {ff: self.results[ff] for ff in self.results}

            period = self.results['parameters']['P']['value']
            mt = self.results['parameters']['mt']['value']
            mt += round((np.mean(self.data[set_number][0]) - mt) / period) * period

            prediction = (self.mid_time +
                          round((np.mean(self.data[set_number][0]) - self.mid_time) / self.period) * self.period)

            duration = transit_duration(self.rp_over_rs, self.period, self.sma_over_rs,
                                        self.inclination, self.eccentricity, self.periastron)

            ingress = prediction - duration / 2
            egress = prediction + duration / 2

            set_indices = np.where(self.data_set_number == set_number)

            plt.subplot2grid((4, 1), (0, 0), rowspan=3)

            plt.plot(self.results['output_series']['phase'][set_indices],
                     self.results['input_series']['value'][set_indices], 'ko', ms=2)
            plt.plot(self.results['output_series']['phase'][set_indices],
                     self.results['output_series']['model'][set_indices], 'r-')

            plt.ylim(min(self.results['output_series']['model'][set_indices])
                     - 5 * np.std(self.results['output_series']['residuals'][set_indices]),
                     max(self.results['output_series']['model'][set_indices])
                     + 5 * np.std(self.results['output_series']['residuals'][set_indices]))

            plt.yticks(plt.yticks()[0][1:])
            plt.ylabel(r'$\mathrm{relative} \ \mathrm{flux}$', fontsize=15)

            plt.ylim(min(self.results['output_series']['model'][set_indices])
                     - 5 * np.std(self.results['output_series']['residuals'][set_indices]),
                     max(self.results['output_series']['model'][set_indices])
                     + 5 * np.std(self.results['output_series']['residuals'][set_indices]))

            x_max = max(np.abs(self.results['output_series']['phase'][set_indices]) +
                        0.05 * (max(self.results['output_series']['phase'][set_indices]) -
                                min(self.results['output_series']['phase'][set_indices])))
            plt.xlim(-x_max, x_max)
            plt.tick_params(labelbottom='off')

            rpstr = '{0}{1}{2}{3}{4}{5}{6}{7}'.format(
                r'$R_\mathrm{p}/R_* = ', self.results['parameters']['rp']['print_value'], '_{-',
                self.results['parameters']['rp']['print_m_error'], '}', '^{+',
                self.results['parameters']['rp']['print_p_error'], '}$')
            mtstr = '{0}{1}{2}{3}{4}{5}{6}{7}'.format(
                r'$T_\mathrm{HJD} = ', self.results['parameters']['mt']['print_value'], '_{-',
                self.results['parameters']['mt']['print_m_error'], '}', '^{+',
                self.results['parameters']['mt']['print_p_error'], '}$')

            plt.text(plt.xlim()[0] + 0.5 * (plt.xlim()[-1] - plt.xlim()[0]),
                     plt.ylim()[0] + 0.07 * (plt.ylim()[-1] - plt.ylim()[0]),
                     '{0}\n{1}'.format(rpstr, mtstr), ha='center', va='center', fontsize=10)

            plt.axvline((ingress - mt) / period, 0.3, 1.0, ls='--', c='k', lw=0.75)
            plt.text((ingress - mt) / period, plt.ylim()[0] + 0.3 * (plt.ylim()[1] - plt.ylim()[0]),
                     '{0}{1}{2}{3}{4}'.format(r'$\mathrm{predicted}$', '\n', r'$\mathrm{ingress}$', '\n',
                                              r'$\mathrm{start}$'),
                     ha='right', va='top', fontsize=10)
            plt.axvline((egress - mt) / period, 0.3, 1.0, ls='--', c='k', lw=0.75)
            plt.text((egress - mt) / period, plt.ylim()[0] + 0.3 * (plt.ylim()[1] - plt.ylim()[0]),
                     '{0}{1}{2}{3}{4}'.format(r'$\mathrm{predicted}$', '\n', r'$\mathrm{egress}$', '\n',
                                              r'$\mathrm{end}$'),
                     ha='left', va='top', fontsize=10)

            plt.suptitle('{0}{1}{2}'.format(r'$\mathbf{', target, '}$'), fontsize=20)
            plt.text(plt.xlim()[1], plt.ylim()[1], '{0}{1}{2}'.format(r'$', data_dates[set_number], '$'),
                     fontsize=12, ha='right', va='bottom')

            plt.subplot(4, 1, 4)
            plt.cla()
            plt.plot(self.results['output_series']['phase'][set_indices],
                     self.results['output_series']['residuals'][set_indices], 'ko', ms=2)
            plt.plot(self.results['output_series']['phase'][set_indices],
                     np.zeros_like(self.results['output_series']['phase'][set_indices]), 'r-')

            plt.ylim(- 5 * np.std(self.results['output_series']['residuals'][set_indices]),
                     5 * np.std(self.results['output_series']['residuals'][set_indices]))

            plt.xlabel(r'$\mathrm{phase}$', fontsize=15)
            plt.ylabel(r'$\mathrm{residuals}$', fontsize=15)

            plt.xlim(-x_max, x_max)
            plt.text(plt.xlim()[0] + 0.02 * (plt.xlim()[-1] - plt.xlim()[0]),
                     plt.ylim()[0] + 0.07 * (plt.ylim()[-1] - plt.ylim()[0]),
                     r'$\mathrm{rms}_\mathrm{res} = %.1e$' %
                     np.std(self.results['output_series']['residuals'][set_indices]), fontsize=10)

            plt.subplots_adjust(left=0.15, right=0.975, bottom=0.12, top=0.9, hspace=0.0)

            plt.savefig(os.path.join(os.path.split(export_file)[0],
                                     'set_{0}_{1}'.format(str(set_number + 1), os.path.split(export_file)[1])),
                        transparent=True)
            plt.close('all')

    def plot_detrended_models(self, export_file, target=None, data_dates=None, return_plot=False):

        if target is None:
            target = ' '

        if data_dates is None:
            data_dates = map(str, ['set_{0}'.format(str(ff)) for ff in range(1, self.total_sets + 1)])

        for set_number in range(self.total_sets):

            fig = plt.figure(set_number + 1)
            fig.set_tight_layout(False)

            self.results = {ff: self.results[ff] for ff in self.results}

            period = self.results['parameters']['P']['value']
            mt = self.results['parameters']['mt']['value']
            mt += round((np.mean(self.data[set_number][0]) - mt) / period) * period

            prediction = (self.mid_time +
                          round((np.mean(self.data[set_number][0]) - self.mid_time) / self.period) * self.period)

            duration = transit_duration(self.rp_over_rs, self.period, self.sma_over_rs,
                                        self.inclination, self.eccentricity, self.periastron)

            ingress = prediction - duration / 2
            egress = prediction + duration / 2

            set_indices = np.where(self.data_set_number == set_number)

            plt.subplot2grid((4, 1), (0, 0), rowspan=3)

            plt.plot(self.results['detrended_output_series']['phase'][set_indices],
                     self.results['detrended_input_series']['value'][set_indices], 'ko', ms=2)
            plt.plot(self.results['detrended_output_series']['phase'][set_indices],
                     self.results['detrended_output_series']['model'][set_indices], 'r-')

            plt.ylim(min(self.results['detrended_output_series']['model'][set_indices])
                     - 5 * np.std(self.results['detrended_output_series']['residuals'][set_indices]),
                     max(self.results['detrended_output_series']['model'][set_indices])
                     + 5 * np.std(self.results['detrended_output_series']['residuals'][set_indices]))

            plt.yticks(plt.yticks()[0][1:])
            plt.ylabel(r'$\mathrm{relative} \ \mathrm{flux}$', fontsize=15)

            plt.ylim(min(self.results['detrended_output_series']['model'][set_indices])
                     - 5 * np.std(self.results['detrended_output_series']['residuals'][set_indices]),
                     max(self.results['detrended_output_series']['model'][set_indices])
                     + 5 * np.std(self.results['detrended_output_series']['residuals'][set_indices]))

            plt.ylim(-1.17647 * (- plt.ylim()[0] + 0.15 * plt.ylim()[1]), plt.ylim()[1])

            x_max = max(np.abs(self.results['detrended_output_series']['phase'][set_indices]) +
                        0.05 * (max(self.results['detrended_output_series']['phase'][set_indices]) -
                                min(self.results['detrended_output_series']['phase'][set_indices])))
            plt.xlim(-x_max, x_max)
            plt.tick_params(labelbottom='off')

            rpstr = '{0}{1}{2}{3}{4}{5}{6}{7}'.format(
                r'$R_\mathrm{p}/R_* = ', self.results['parameters']['rp']['print_value'], '_{-',
                self.results['parameters']['rp']['print_m_error'], '}', '^{+',
                self.results['parameters']['rp']['print_p_error'], '}$')
            mtstr = '{0}{1}{2}{3}{4}{5}{6}{7}'.format(
                r'$T_\mathrm{HJD} = ', self.results['parameters']['mt']['print_value'], '_{-',
                self.results['parameters']['mt']['print_m_error'], '}', '^{+',
                self.results['parameters']['mt']['print_p_error'], '}$')

            plt.text(plt.xlim()[0] + 0.5 * (plt.xlim()[-1] - plt.xlim()[0]),
                     plt.ylim()[0] + 0.07 * (plt.ylim()[-1] - plt.ylim()[0]),
                     '{0}{1}{2}'.format(rpstr, '\n', mtstr), ha='center', va='center', fontsize=10)

            plt.axvline((ingress - mt) / period, 0.3, 1.0, ls='--', c='k', lw=0.75)
            plt.text((ingress - mt) / period, plt.ylim()[0] + 0.3 * (plt.ylim()[1] - plt.ylim()[0]),
                     '{0}{1}{2}{3}{4}'.format(r'$\mathrm{predicted}$', '\n', r'$\mathrm{ingress}$', '\n',
                                              r'$\mathrm{start}$'),
                     ha='right', va='top', fontsize=10)
            plt.axvline((egress - mt) / period, 0.3, 1.0, ls='--', c='k', lw=0.75)
            plt.text((egress - mt) / period, plt.ylim()[0] + 0.3 * (plt.ylim()[1] - plt.ylim()[0]),
                     '{0}{1}{2}{3}{4}'.format(r'$\mathrm{predicted}$', '\n', r'$\mathrm{egress}$', '\n',
                                              r'$\mathrm{end}$'),
                     ha='left', va='top', fontsize=10)

            plt.suptitle('{0}{1}{2}'.format(r'$\mathbf{', target, '}$'), fontsize=20)
            plt.text(plt.xlim()[1], plt.ylim()[1], '{0}{1}{2}'.format(r'$', data_dates[set_number], '$'),
                     fontsize=12, ha='right', va='bottom')

            plt.subplot(4, 1, 4)
            plt.cla()
            plt.plot(self.results['detrended_output_series']['phase'][set_indices],
                     self.results['detrended_output_series']['residuals'][set_indices], 'ko', ms=2)
            plt.plot(self.results['detrended_output_series']['phase'][set_indices],
                     np.zeros_like(self.results['detrended_output_series']['phase'][set_indices]), 'r-')

            plt.ylim(- 5 * np.std(self.results['detrended_output_series']['residuals'][set_indices]),
                     5 * np.std(self.results['detrended_output_series']['residuals'][set_indices]))

            plt.xlabel(r'$\mathrm{phase}$', fontsize=15)
            plt.ylabel(r'$\mathrm{residuals}$', fontsize=15)

            plt.xlim(-x_max, x_max)

            plt.text(plt.xlim()[0] + 0.02 * (plt.xlim()[-1] - plt.xlim()[0]),
                     plt.ylim()[0] + 0.07 * (plt.ylim()[-1] - plt.ylim()[0]),
                     r'$\mathrm{rms}_\mathrm{res} = %.1e$' %
                     np.std(self.results['detrended_output_series']['residuals'][set_indices]),
                     fontsize=10)

            plt.subplots_adjust(left=0.15, right=0.975, bottom=0.12, top=0.9, hspace=0.0)

            plt.savefig(os.path.join(os.path.split(export_file)[0],
                                     'set_{0}_{1}'.format(str(set_number + 1), os.path.split(export_file)[1])),
                        transparent=True)
            if return_plot:
                return [plt.figure(ff + 1) for ff in range(self.total_sets)]
            else:
                plt.close('all')

    def save_models(self, export_file):

        for set_number in range(self.total_sets):

            self.results = {ff: self.results[ff] for ff in self.results}

            set_indices = np.where(self.data_time == self.data[set_number][0])

            np.savetxt(os.path.join(os.path.split(export_file)[0],
                                    'set_{0}_{1}'.format(str(set_number + 1), os.path.split(export_file)[1])),
                       np.swapaxes([self.results['input_series']['hjd'][set_indices],
                                    self.results['output_series']['phase'][set_indices],
                                    self.results['input_series']['value'][set_indices],
                                    self.results['input_series']['error'][set_indices],
                                    self.results['output_series']['model'][set_indices],
                                    self.results['output_series']['residuals'][set_indices]
                                    ], 0, 1))

    def save_detrended_models(self, export_file):

        for set_number in range(self.total_sets):

            self.results = {ff: self.results[ff] for ff in self.results}

            set_indices = np.where(self.data_time == self.data[set_number][0])

            np.savetxt(os.path.join(os.path.split(export_file)[0],
                                    'set_{0}_{1}'.format(str(set_number + 1), os.path.split(export_file)[1])),
                       np.swapaxes([self.results['detrended_input_series']['hjd'][set_indices],
                                    self.results['detrended_output_series']['phase'][set_indices],
                                    self.results['detrended_input_series']['value'][set_indices],
                                    self.results['detrended_input_series']['error'][set_indices],
                                    self.results['detrended_output_series']['model'][set_indices],
                                    self.results['detrended_output_series']['residuals'][set_indices]
                                    ], 0, 1))


class TransitAndHubbleFitting:

    def __init__(self, data,
                 apply_up_down_stream_correction,
                 exclude_initial_orbits, exclude_final_orbits, exclude_initial_orbit_points,
                 first_orbit_ramp, second_order_ramp, mid_orbit_ramps,
                 method, limb_darkening_coefficients,
                 rp_over_rs, period, sma_over_rs, eccentricity, inclination, periastron, mid_time,
                 iterations, walkers, burn, precision=3,
                 exp_time=0, time_factor=1,
                 fit_rp_over_rs=False, fit_period=False, fit_sma_over_rs=False,
                 fit_eccentricity=False, fit_inclination=False, fit_periastron=False, fit_mid_time=False,
                 counter=True, counter_window=False):

        # TODO check input parameters

        self.method = method
        self.limb_darkening_coefficients = limb_darkening_coefficients
        self.fit_ld = False

        if isinstance(self.limb_darkening_coefficients, str):
            if self.limb_darkening_coefficients == 'fit':
                self.fit_ld = True
                if self.method == 'linear':
                    self.limb_darkening_coefficients = [0.5]
                elif self.method in ['quad', 'sqrt']:
                    self.limb_darkening_coefficients = [0.5, 0.5]
                elif self.method == 'claret':
                    self.limb_darkening_coefficients = [0.5, 0.5, 0.5, 0.5]

        if self.method == 'claret':
            self.total_ldcs = 4
        elif self.method in ['quad', 'sqrt']:
            self.total_ldcs = 2
        elif self.method == 'linear':
            self.total_ldcs = 1

        self.rp_over_rs = rp_over_rs
        self.period = period
        self.sma_over_rs = sma_over_rs
        self.eccentricity = eccentricity
        self.inclination = inclination
        self.periastron = periastron
        self.mid_time = mid_time
        self.precision = precision
        self.time_factor = time_factor

        self.fit_second_order_ramp = second_order_ramp

        self.fit_rp_over_rs = fit_rp_over_rs
        self.fit_period = fit_period
        self.fit_sma_over_rs = fit_sma_over_rs
        self.fit_eccentricity = fit_eccentricity
        self.fit_inclination = fit_inclination
        self.fit_periastron = fit_periastron
        self.fit_mid_time = fit_mid_time

        self.data = {}
        self.sets = ['set_{0}'.format(str(ff).zfill(2)) for ff in range(len(data))]
        self.total_sets = len(data)
        self.data_set_number = np.array([])
        self.data_time = np.array([])
        data_flux = np.array([])
        data_flux_error = np.array([])

        for set_number, set_arrays in enumerate(data):

            new_set_arrays = [ff for ff in set_arrays]

            # up-stream / down-stream correction

            star_y_position_array = new_set_arrays[1]
            spectrum_direction_array = new_set_arrays[2]
            scan_length_array = new_set_arrays[3]

            if apply_up_down_stream_correction:
                test1 = star_y_position_array[0] - 507
                test2 = test1 + spectrum_direction_array[0] * scan_length_array[0]
                if test1 * test2 < 0:
                    apply_up_down_stream_correction = True
                else:
                    apply_up_down_stream_correction = False

            if apply_up_down_stream_correction:
                for scan_direction in [1.0, -1.0]:
                    fr = np.where(spectrum_direction_array == scan_direction)[0]
                    if len(fr) > 0:
                        zerofr = star_y_position_array.value[fr]
                        sigmfr = scan_length_array.value[fr]
                        begnfr = zerofr
                        fitfr = np.poly1d(np.polyfit(begnfr, sigmfr, 1))
                        for ii in range(4, len(new_set_arrays)):
                            new_set_arrays[ii][fr] = new_set_arrays[ii][fr] * fitfr(begnfr[0]) / fitfr(begnfr)
                            new_set_arrays[ii][fr] = new_set_arrays[ii][fr] * fitfr(begnfr[0]) / fitfr(begnfr)

            # exclude orbits / points

            heliocentric_julian_date_array = new_set_arrays[0]

            indices_to_remain = np.arange(len(heliocentric_julian_date_array))

            if exclude_initial_orbits > 0:
                htime = heliocentric_julian_date_array
                orbits = np.where(abs(htime - np.roll(htime, 1)) > 30.0 / 60.0 / 24.0)[0]
                indices_to_remain = indices_to_remain[orbits[exclude_initial_orbits]:]

            if exclude_final_orbits > 0:
                htime = heliocentric_julian_date_array[indices_to_remain]
                orbits = np.where(abs(htime - np.roll(htime, 1)) > 30.0 / 60.0 / 24.0)[0]
                indices_to_remain = indices_to_remain[:orbits[-exclude_final_orbits]]

            if exclude_initial_orbit_points > 0:
                htime = heliocentric_julian_date_array[indices_to_remain]
                orbits = np.where(abs(htime - np.roll(htime, 1)) > 30.0 / 60.0 / 24.0)[0]
                indices_to_remain = np.delete(indices_to_remain,
                                              np.concatenate(
                                                  [orbits + i for i in range(exclude_initial_orbit_points)]))

            new_set_arrays = [ff[indices_to_remain] for ff in new_set_arrays]

            # define hst orbital phases

            heliocentric_julian_date_array = new_set_arrays[0]

            if mid_orbit_ramps:
                htime = heliocentric_julian_date_array
                orbits = np.where(abs(htime - np.roll(htime, 1)) > 30.0 / 60.0 / 24.0)[0]
                dumps = np.where(abs(htime - np.roll(htime, 1)) > 5.0 / 60.0 / 24.0)[0]
                dphase = np.zeros(len(htime))
                for i in range(1, len(dumps)):
                    if dumps[i] not in orbits:
                        if i != len(dumps) - 1:
                            for j in range(dumps[i], dumps[i + 1]):
                                dphase[j] = 1
                        else:
                            for j in range(dumps[i], len(dphase)):
                                dphase[j] = 1
            else:
                htime = heliocentric_julian_date_array
                dphase = np.zeros(len(htime))

            if first_orbit_ramp:
                htime = heliocentric_julian_date_array
                if mid_orbit_ramps:
                    orbits = np.where(abs(htime - np.roll(htime, 1)) > 5.0 / 60.0 / 24.0)[0]
                else:
                    orbits = np.where(abs(htime - np.roll(htime, 1)) > 30.0 / 60.0 / 24.0)[0]
                orbits = htime[orbits]
                fphase = np.where(htime < orbits[1], 1, 0)

            else:
                htime = heliocentric_julian_date_array
                fphase = np.zeros(len(htime))

            htime = heliocentric_julian_date_array
            if mid_orbit_ramps:
                orbits = np.where(abs(htime - np.roll(htime, 1)) > 5.0 / 60.0 / 24.0)[0]
            else:
                orbits = np.where(abs(htime - np.roll(htime, 1)) > 30.0 / 60.0 / 24.0)[0]
            t0s = htime[orbits]
            ophase = []
            for pp in t0s:
                ppp = htime - pp
                ppp = np.where(ppp < 0, 1000, ppp)
                ophase.append(ppp)

            ophase = np.min(ophase, 0)

            # outliers filter

            lightcurves = [ff for ff in new_set_arrays[6::2]]
            ica = FastICA(n_components=len(lightcurves), max_iter=1000)
            components = ica.fit_transform(np.array(lightcurves).T).T

            indices_to_remain = []
            for i in components:
                indices_to_remain.append(
                    np.array(np.abs(i - np.median(i)) < 20 * np.median(np.abs(i - np.median(i)))))
            indices_to_remain = np.where(np.product(indices_to_remain, 0))[0]
            indices_to_remain = np.sort(np.unique(np.array(indices_to_remain)))

            new_set_arrays = [ff[indices_to_remain] for ff in new_set_arrays]
            ophase = ophase[indices_to_remain]
            dphase = dphase[indices_to_remain]
            fphase = fphase[indices_to_remain]

            # match forward and reverse scans

            spectrum_direction_array = new_set_arrays[2]
            flux_array = new_set_arrays[4]

            fr = np.where(spectrum_direction_array > 0)[0]
            if len(fr) != len(spectrum_direction_array):

                fr_out = np.where(spectrum_direction_array > 0)[0]
                rv_out = np.where(spectrum_direction_array < 0)[0]
                shift = np.mean(flux_array[fr_out]) / np.mean(flux_array[rv_out])
                for ii in range(4, len(new_set_arrays)):
                    new_set_arrays[ii][fr] = new_set_arrays[ii][fr] / shift

            if set_number == 0:
                time_shift = round((np.mean(new_set_arrays[0]) - self.mid_time) / self.period)
                self.mid_time += time_shift * self.period
                if self.fit_mid_time:
                    self.fit_mid_time[0] += time_shift * self.period
                    self.fit_mid_time[1] += time_shift * self.period

            data_flux = np.append(data_flux, new_set_arrays[4])
            data_flux_error = np.append(data_flux_error, new_set_arrays[5])
            self.data_time = np.append(self.data_time, new_set_arrays[0])
            self.data_set_number = np.int_(np.append(self.data_set_number,
                                                     np.ones_like(new_set_arrays[0]) * set_number))

            if exp_time == 0:
                hjd_hd = new_set_arrays[0]

            else:
                hjd_hd = np.array([])
                for i in range(self.time_factor):
                    hjd_hd = np.append(hjd_hd, new_set_arrays[0] - (exp_time / (60.0 * 60.0 * 24.0)) / 2.0 +
                                       (i + 0.5) * (exp_time / (60.0 * 60.0 * 24.0)) / self.time_factor)

            epoch = round((np.mean(new_set_arrays[0]) - self.mid_time) / self.period)

            self.data[self.sets[set_number]] = {'epoch': epoch,
                                                'hjd': new_set_arrays[0],
                                                'ophase': ophase,
                                                'dphase': dphase,
                                                'fphase': fphase,
                                                'scan': new_set_arrays[2],
                                                'hjd_hd': hjd_hd,
                                                'flux': new_set_arrays[4],
                                                'error': new_set_arrays[5],
                                                'pindices': []}

        names = []
        print_names = []
        limits1 = []
        limits2 = []
        initial = []

        for set_number, set_name in enumerate(self.sets):

            flux = self.data[set_name]['flux']
            scan = self.data[set_name]['scan']
            fphase = self.data[set_name]['fphase']
            dphase = self.data[set_name]['dphase']

            # forward scans normalisation factors
            names.append('n_w_for_{0}'.format(str(set_number)))
            print_names.append('n{}{}{}'.format('w', 'for', str(set_number)))
            initial.append(np.max(flux))
            if (scan < 0).all():
                limits1.append(np.nan)
                limits2.append(np.nan)
            else:
                limits1.append(np.max(flux) * 0.99)
                limits2.append(np.max(flux) * 1.01)
            self.data[set_name]['pindices'].append(len(names) - 1)

            print(np.median(flux), np.max(flux), np.max(flux) * 0.99, np.max(flux) * 1.01)

            # reverse scans normalisation factors
            names.append('n_w_rev_{0}'.format(str(set_number)))
            print_names.append('n{}{}{}'.format('w', 'rev', str(set_number)))
            initial.append(np.max(flux))
            if (scan > 0).all():
                limits1.append(np.nan)
                limits2.append(np.nan)
            else:
                limits1.append(np.max(flux) * 0.99)
                limits2.append(np.max(flux) * 1.01)
            self.data[set_name]['pindices'].append(len(names) - 1)

            # long term ramp - 1st order
            names.append('r_a1_{0}'.format(str(set_number)))
            print_names.append('r{}{}'.format('a1', str(set_number)))
            initial.append(0.001)
            limits1.append(-1.0)
            limits2.append(1.0)
            self.data[set_name]['pindices'].append(len(names) - 1)

            # long term ramp - 2nd order
            names.append('r_a2_{0}'.format(str(set_number)))
            print_names.append('r{}{}'.format('a2', str(set_number)))
            initial.append(0.0)
            if self.fit_second_order_ramp:
                limits1.append(-1.0)
                limits2.append(1.0)
            else:
                limits1.append(np.nan)
                limits2.append(np.nan)
            self.data[set_name]['pindices'].append(len(names) - 1)

            # sort term ramp - amplitude
            names.append('r_b1_{0}'.format(str(set_number)))
            print_names.append('r{}{}'.format('b1', str(set_number)))
            initial.append(0.001)
            limits1.append(-1.0)
            limits2.append(1.0)
            self.data[set_name]['pindices'].append(len(names) - 1)

            # sort term mid-orbit ramp - amplitude
            names.append('mor_b1_{0}'.format(str(set_number)))
            print_names.append('mor{}{}'.format('b1', str(set_number)))
            initial.append(0.001)
            if np.sum(dphase ** 2) == 0:
                limits1.append(np.nan)
                limits2.append(np.nan)
            else:
                limits1.append(-1.0)
                limits2.append(1.0)
            self.data[set_name]['pindices'].append(len(names) - 1)

            # sort term first-orbit ramp - amplitude
            names.append('for_b1_{0}'.format(str(set_number)))
            print_names.append('for{}{}'.format('b1', str(set_number)))
            initial.append(0.001)
            if np.sum(fphase ** 2) == 0:
                limits1.append(np.nan)
                limits2.append(np.nan)
            else:
                limits1.append(-1.0)
                limits2.append(1.0)
            self.data[set_name]['pindices'].append(len(names) - 1)

            # sort term ramp - decay
            names.append('r_b2_{0}'.format(str(set_number)))
            print_names.append('r{}{}'.format('b2', str(set_number)))
            initial.append(250.0)
            limits1.append(50.0)
            limits2.append(500.0)
            self.data[set_name]['pindices'].append(len(names) - 1)

            # sort term mid-orbit ramp - decay
            names.append('mor_b2_{0}'.format(str(set_number)))
            print_names.append('mor{}{}'.format('b2', str(set_number)))
            initial.append(250.0)
            if np.sum(dphase ** 2) == 0:
                limits1.append(np.nan)
                limits2.append(np.nan)
            else:
                limits1.append(50.0)
                limits2.append(500.0)
            self.data[set_name]['pindices'].append(len(names) - 1)

            # sort term first-orbit ramp - decay
            names.append('for_b2_{0}'.format(str(set_number)))
            print_names.append('for{}{}'.format('b2', str(set_number)))
            initial.append(150.0)
            if np.sum(fphase ** 2) == 0:
                limits1.append(np.nan)
                limits2.append(np.nan)
            else:
                limits1.append(50.0)
                limits2.append(500.0)
            self.data[set_name]['pindices'].append(len(names) - 1)

            # rp
            names.append('rp_{0}'.format(str(set_number)))
            print_names.append('Rp/R*{}'.format(str(set_number)))
            initial.append(self.rp_over_rs)
            limits1.append(self.fit_rp_over_rs[0])
            limits2.append(self.fit_rp_over_rs[1])
            self.data[set_name]['pindices'].append(len(names) - 1)

        self.len_systematics = int(len(names) / self.total_sets)

        names.append('ldc1')
        print_names.append('ldc_1')
        initial.append(self.limb_darkening_coefficients[0])
        if self.fit_ld:
            limits1.append(0.000001)
            limits2.append(0.999999)
        else:
            limits1.append(np.nan)
            limits2.append(np.nan)

        for set_name in self.sets:
            self.data[set_name]['pindices'].append(len(names) - 1)

        if self.method in ['claret', 'quad', 'sqrt']:

            names.append('ldc2')
            print_names.append('ldc_2')
            initial.append(self.limb_darkening_coefficients[1])
            if self.fit_ld:
                limits1.append(0.000001)
                limits2.append(0.999999)
            else:
                limits1.append(np.nan)
                limits2.append(np.nan)

            for set_name in self.sets:
                self.data[set_name]['pindices'].append(len(names) - 1)

        if self.method == 'claret':

            names.append('ldc3')
            print_names.append('ldc_3')
            initial.append(self.limb_darkening_coefficients[2])
            if self.fit_ld:
                limits1.append(0.000001)
                limits2.append(0.999999)
            else:
                limits1.append(np.nan)
                limits2.append(np.nan)

            for set_name in self.sets:
                self.data[set_name]['pindices'].append(len(names) - 1)

            names.append('ldc4')
            print_names.append('ldc_4')
            initial.append(self.limb_darkening_coefficients[3])
            if self.fit_ld:
                limits1.append(0.000001)
                limits2.append(0.999999)
            else:
                limits1.append(np.nan)
                limits2.append(np.nan)

            for set_name in self.sets:
                self.data[set_name]['pindices'].append(len(names) - 1)

        for pindex in range(len(names), len(names) + 6):
            for set_name in self.sets:
                self.data[set_name]['pindices'].append(pindex)

        for set_name in self.sets:
            self.data[set_name]['pindices'] = np.int_(self.data[set_name]['pindices'])

        names += ['P', 'a', 'e', 'i', 'w', 'mt']
        print_names += ['P', 'a/R_*', 'e', 'i', '\omega', 'T_0']

        initial += [self.period, self.sma_over_rs, self.eccentricity,
                    self.inclination, self.periastron, self.mid_time]

        limits = limits1 + [self.fit_period, self.fit_sma_over_rs, self.fit_eccentricity,
                            self.fit_inclination, self.fit_periastron, self.fit_mid_time]

        for var in range(self.len_systematics * self.total_sets + self.total_ldcs, len(names)):

            try:
                initial[var] = float(initial[var])
            except:
                raise RuntimeError('Improper value for {0}'.format(names[var]))

            if limits[var] is False:
                limits1.append(np.nan)
                limits2.append(np.nan)

            elif limits[var] is None:
                limits1.append(np.nan)
                limits2.append(np.nan)

            else:
                try:
                    if len(np.array(limits[var])) != 2:
                        raise RuntimeError('Improper limits for {0}'.format(names[var]))
                except:
                    raise RuntimeError('Improper limits for {0}'.format(names[var]))

                if initial[var] < np.array(limits[var])[0] or initial[var] > np.array(limits[var])[1]:
                    raise RuntimeError('Initial value for {0} is outside the range of the prior.'.format(
                        names[var]))
                else:
                    limits1.append(np.array(limits[var])[0])
                    limits2.append(np.array(limits[var])[1])

        self.fitting = EmceeFitting([data_flux, data_flux_error],
                                    self.full_model, initial, limits1, limits2,
                                    walkers, iterations, burn,
                                    names=names, print_names=print_names,
                                    counter=counter, counter_window=counter_window, strech_prior=1.0)

        self.results = {}
        self.mcmc_run_complete = False

    def detrend_model(self, *model_variables):

        model = []

        for set_name in self.sets:

            (model_norm_f, model_norm_r, model_r_a1, model_r_a2, model_r_b1, model_mor_b1, model_for_b1, model_r_b2,
             model_mor_b2, model_for_b2, model_rp) = np.array(
                model_variables)[self.data[set_name]['pindices']][:self.len_systematics]
            model_period = np.array(
                model_variables)[self.data[set_name]['pindices']][self.len_systematics + self.total_ldcs]
            model_mid_time = np.array(model_variables)[self.data[set_name]['pindices']][-1]

            model_ophase = self.data[set_name]['ophase']
            model_dphase = self.data[set_name]['dphase']
            model_fphase = self.data[set_name]['fphase']
            model_scan = self.data[set_name]['scan']
            model_hjd = self.data[set_name]['hjd']
            model_epoch = self.data[set_name]['epoch']
            model_vtime = model_hjd - (model_mid_time + model_epoch * model_period)

            normalisation = np.where(model_scan > 0, model_norm_f, model_norm_r)
            detrend1 = (1.0 - model_r_a1 * model_vtime + model_r_a2 * (model_vtime ** 2))
            ramp_ampl = np.where(model_dphase == 0, model_r_b1, model_mor_b1)
            ramp_ampl = np.where(model_fphase == 0, ramp_ampl, model_for_b1)
            ramp_decay = np.where(model_dphase == 0, model_r_b2, model_mor_b2)
            ramp_decay = np.where(model_fphase == 0, ramp_decay, model_for_b2)
            detrend2 = 1.0 - ramp_ampl * np.exp(- ramp_decay * model_ophase)

            model.append(normalisation * detrend1 * detrend2)

        return np.concatenate(model)

    def transit_model(self, *model_variables):

        model = []

        for set_name in self.sets:
            model_rp_over_rs = np.array(model_variables)[self.data[set_name]['pindices']][10]

            model_hjd = self.data[set_name]['hjd']
            model_hjd_hd = self.data[set_name]['hjd_hd']

            limb_darkening_coefficients = np.array(
                model_variables)[self.data[set_name]['pindices']][self.len_systematics:
                                                                  self.len_systematics + self.total_ldcs]

            z_over_rs = transit_projected_distance(*np.array(
                model_variables)[self.data[set_name]['pindices']][self.len_systematics + self.total_ldcs:],
                time_array=model_hjd_hd)

            transit_hr = transit_flux_drop(self.method, limb_darkening_coefficients, model_rp_over_rs, z_over_rs,
                                           precision=self.precision)

            model.append(np.mean(np.reshape(transit_hr, (self.time_factor, len(model_hjd))), 0))

        return np.concatenate(model)

    def full_model(self, *model_variables):

        return self.detrend_model(*model_variables) * self.transit_model(*model_variables)

    def run_mcmc(self):

        self.fitting.run_mcmc()
        self.results = self.fitting.results
        self.mcmc_run_complete = True

        self.results['input_series']['hjd'] = self.data_time

        period = self.results['parameters']['P']['value']
        mt = self.results['parameters']['mt']['value']
        self.results['output_series']['phase'] = \
            (self.data_time - mt) / period - np.round((self.data_time - mt) / period)

        self.results['detrended_input_series'] = {
            'hjd': self.results['input_series']['hjd'],
            'value': self.results['input_series']['value'] / self.detrend_model(*self.results['parameters_final']),
            'error': self.results['input_series']['error'] / self.detrend_model(*self.results['parameters_final'])}

        self.results['detrended_output_series'] = {
            'phase': self.results['output_series']['phase'],
            'model': self.results['output_series']['model'] / self.detrend_model(*self.results['parameters_final']),
            'residuals': (self.results['output_series']['residuals']
                          / self.detrend_model(*self.results['parameters_final']))}

        self.results['detrended_statistics'] = {ff: self.results['statistics'][ff] for ff in self.results['statistics']}
        self.results['detrended_statistics']['res_std'] = np.std(self.results['detrended_output_series']['residuals'])

    def save_all(self, export_file):

        pickle.dump(self.results, open(export_file, 'wb'))

    def save_results(self, export_file):

        self.fitting.save_results(export_file)

    def plot_corner(self, export_file):

        self.fitting.plot_corner(export_file)

    def plot_traces(self, export_file):

        self.fitting.plot_traces(export_file)

    def plot_models(self, export_file, target=None, data_dates=None):

        if target is None:
            target = ' '

        if data_dates is None:
            data_dates = map(str, ['set_{0}'.format(str(ff)) for ff in range(1, self.total_sets + 1)])

        for set_number in range(self.total_sets):

            set_indices = np.where(self.data_set_number == set_number)

            fig = plt.figure(set_number + 1)
            fig.set_tight_layout(False)

            self.results = {ff: self.results[ff] for ff in self.results}

            period = self.results['parameters']['P']['value']
            mt = self.results['parameters']['mt']['value']

            mt += round((np.mean(self.data_time[set_indices]) - mt) / period) * period

            prediction = (self.mid_time +
                          round((np.mean(self.data_time[set_indices]) - self.mid_time) / self.period) * self.period)

            duration = transit_duration(self.rp_over_rs, self.period, self.sma_over_rs,
                                        self.inclination, self.eccentricity, self.periastron)

            ingress = prediction - duration / 2
            egress = prediction + duration / 2

            plt.subplot2grid((4, 1), (0, 0), rowspan=3)

            plt.plot(self.results['output_series']['phase'][set_indices],
                     self.results['input_series']['value'][set_indices], 'ko', ms=2)
            plt.plot(self.results['output_series']['phase'][set_indices],
                     self.results['output_series']['model'][set_indices], 'r-')

            plt.ylim(min(self.results['output_series']['model'][set_indices])
                     - 5 * np.std(self.results['output_series']['residuals'][set_indices]),
                     max(self.results['output_series']['model'][set_indices])
                     + 5 * np.std(self.results['output_series']['residuals'][set_indices]))

            plt.yticks(plt.yticks()[0][1:])
            plt.ylabel(r'$\mathrm{relative} \ \mathrm{flux}$', fontsize=15)

            plt.ylim(min(self.results['output_series']['model'][set_indices])
                     - 5 * np.std(self.results['output_series']['residuals'][set_indices]),
                     max(self.results['output_series']['model'][set_indices])
                     + 5 * np.std(self.results['output_series']['residuals'][set_indices]))

            x_max = max(np.abs(self.results['output_series']['phase'][set_indices]) +
                        0.05 * (max(self.results['output_series']['phase'][set_indices]) -
                                min(self.results['output_series']['phase'][set_indices])))
            plt.xlim(-x_max, x_max)
            plt.tick_params(labelbottom='off')

            rpstr = '{0}{1}{2}{3}{4}{5}{6}{7}'.format(
                r'$R_\mathrm{p}/R_* = ', self.results['parameters']['rp_{0}'.format(str(set_number))]['print_value'],
                '_{-', self.results['parameters']['rp_{0}'.format(str(set_number))]['print_m_error'], '}', '^{+',
                self.results['parameters']['rp_{0}'.format(str(set_number))]['print_p_error'], '}$')
            mtstr = '{0}{1}{2}{3}{4}{5}{6}{7}'.format(
                r'$T_\mathrm{HJD} = ', self.results['parameters']['mt']['print_value'], '_{-',
                self.results['parameters']['mt']['print_m_error'], '}', '^{+',
                self.results['parameters']['mt']['print_p_error'], '}$')

            plt.text(plt.xlim()[0] + 0.5 * (plt.xlim()[-1] - plt.xlim()[0]),
                     plt.ylim()[0] + 0.07 * (plt.ylim()[-1] - plt.ylim()[0]),
                     '{0}\n{1}'.format(rpstr, mtstr), ha='center', va='center', fontsize=10)

            plt.axvline((ingress - mt) / period, 0.3, 1.0, ls='--', c='k', lw=0.75)
            plt.text((ingress - mt) / period, plt.ylim()[0] + 0.3 * (plt.ylim()[1] - plt.ylim()[0]),
                     '{0}{1}{2}{3}{4}'.format(r'$\mathrm{predicted}$', '\n', r'$\mathrm{ingress}$', '\n',
                                              r'$\mathrm{start}$'),
                     ha='right', va='top', fontsize=10)
            plt.axvline((egress - mt) / period, 0.3, 1.0, ls='--', c='k', lw=0.75)
            plt.text((egress - mt) / period, plt.ylim()[0] + 0.3 * (plt.ylim()[1] - plt.ylim()[0]),
                     '{0}{1}{2}{3}{4}'.format(r'$\mathrm{predicted}$', '\n', r'$\mathrm{egress}$', '\n',
                                              r'$\mathrm{end}$'),
                     ha='left', va='top', fontsize=10)

            plt.suptitle('{0}{1}{2}'.format(r'$\mathbf{', target, '}$'), fontsize=20)
            plt.text(plt.xlim()[1], plt.ylim()[1], '{0}{1}{2}'.format(r'$', data_dates[set_number], '$'),
                     fontsize=12, ha='right', va='bottom')

            plt.subplot(4, 1, 4)
            plt.cla()
            plt.plot(self.results['output_series']['phase'][set_indices],
                     self.results['output_series']['residuals'][set_indices], 'ko', ms=2)
            plt.plot(self.results['output_series']['phase'][set_indices],
                     np.zeros_like(self.results['output_series']['phase'][set_indices]), 'r-')

            plt.ylim(- 5 * np.std(self.results['output_series']['residuals'][set_indices]),
                     5 * np.std(self.results['output_series']['residuals'][set_indices]))

            plt.xlabel(r'$\mathrm{phase}$', fontsize=15)
            plt.ylabel(r'$\mathrm{residuals}$', fontsize=15)

            plt.xlim(-x_max, x_max)
            plt.text(plt.xlim()[0] + 0.02 * (plt.xlim()[-1] - plt.xlim()[0]),
                     plt.ylim()[0] + 0.07 * (plt.ylim()[-1] - plt.ylim()[0]),
                     r'$\mathrm{rms}_\mathrm{res} = %.1e$' %
                     np.std(self.results['output_series']['residuals'][set_indices]), fontsize=10)

            plt.subplots_adjust(left=0.15, right=0.975, bottom=0.12, top=0.9, hspace=0.0)

            plt.savefig(os.path.join(os.path.split(export_file)[0],
                                     'set_{0}_{1}'.format(str(set_number + 1), os.path.split(export_file)[1])),
                        transparent=True)
            plt.close('all')

    def plot_detrended_models(self, export_file, target=None, data_dates=None, return_plot=False):

        if target is None:
            target = ' '

        if data_dates is None:
            data_dates = map(str, ['set_{0}'.format(str(ff)) for ff in range(1, self.total_sets + 1)])

        for set_number in range(self.total_sets):

            set_indices = np.where(self.data_set_number == set_number)

            fig = plt.figure(set_number + 1)
            fig.set_tight_layout(False)

            self.results = {ff: self.results[ff] for ff in self.results}

            period = self.results['parameters']['P']['value']
            mt = self.results['parameters']['mt']['value']
            mt += round((np.mean(self.data_time[set_indices]) - mt) / period) * period

            prediction = (self.mid_time +
                          round((np.mean(self.data_time[set_indices]) - self.mid_time) / self.period) * self.period)

            duration = transit_duration(self.rp_over_rs, self.period, self.sma_over_rs,
                                        self.inclination, self.eccentricity, self.periastron)

            ingress = prediction - duration / 2
            egress = prediction + duration / 2

            plt.subplot2grid((4, 1), (0, 0), rowspan=3)

            plt.plot(self.results['detrended_output_series']['phase'][set_indices],
                     self.results['detrended_input_series']['value'][set_indices], 'ko', ms=2)
            plt.plot(self.results['detrended_output_series']['phase'][set_indices],
                     self.results['detrended_output_series']['model'][set_indices], 'r-')

            plt.ylim(min(self.results['detrended_output_series']['model'][set_indices])
                     - 5 * np.std(self.results['detrended_output_series']['residuals'][set_indices]),
                     max(self.results['detrended_output_series']['model'][set_indices])
                     + 5 * np.std(self.results['detrended_output_series']['residuals'][set_indices]))

            plt.yticks(plt.yticks()[0][1:])
            plt.ylabel(r'$\mathrm{relative} \ \mathrm{flux}$', fontsize=15)

            plt.ylim(min(self.results['detrended_output_series']['model'][set_indices])
                     - 5 * np.std(self.results['detrended_output_series']['residuals'][set_indices]),
                     max(self.results['detrended_output_series']['model'][set_indices])
                     + 5 * np.std(self.results['detrended_output_series']['residuals'][set_indices]))

            plt.ylim(-1.17647 * (- plt.ylim()[0] + 0.15 * plt.ylim()[1]), plt.ylim()[1])

            x_max = max(np.abs(self.results['detrended_output_series']['phase'][set_indices]) +
                        0.05 * (max(self.results['detrended_output_series']['phase'][set_indices]) -
                                min(self.results['detrended_output_series']['phase'][set_indices])))
            plt.xlim(-x_max, x_max)
            plt.tick_params(labelbottom='off')

            rpstr = '{0}{1}{2}{3}{4}{5}{6}{7}'.format(
                r'$R_\mathrm{p}/R_* = ', self.results['parameters']['rp_{0}'.format(str(set_number))]['print_value'],
                '_{-', self.results['parameters']['rp_{0}'.format(str(set_number))]['print_m_error'], '}', '^{+',
                self.results['parameters']['rp_{0}'.format(str(set_number))]['print_p_error'], '}$')
            mtstr = '{0}{1}{2}{3}{4}{5}{6}{7}'.format(
                r'$T_\mathrm{HJD} = ', self.results['parameters']['mt']['print_value'], '_{-',
                self.results['parameters']['mt']['print_m_error'], '}', '^{+',
                self.results['parameters']['mt']['print_p_error'], '}$')

            plt.text(plt.xlim()[0] + 0.5 * (plt.xlim()[-1] - plt.xlim()[0]),
                     plt.ylim()[0] + 0.07 * (plt.ylim()[-1] - plt.ylim()[0]),
                     '{0}{1}{2}'.format(rpstr, '\n', mtstr), ha='center', va='center', fontsize=10)

            plt.axvline((ingress - mt) / period, 0.3, 1.0, ls='--', c='k', lw=0.75)
            plt.text((ingress - mt) / period, plt.ylim()[0] + 0.3 * (plt.ylim()[1] - plt.ylim()[0]),
                     '{0}{1}{2}{3}{4}'.format(r'$\mathrm{predicted}$', '\n', r'$\mathrm{ingress}$', '\n',
                                              r'$\mathrm{start}$'),
                     ha='right', va='top', fontsize=10)
            plt.axvline((egress - mt) / period, 0.3, 1.0, ls='--', c='k', lw=0.75)
            plt.text((egress - mt) / period, plt.ylim()[0] + 0.3 * (plt.ylim()[1] - plt.ylim()[0]),
                     '{0}{1}{2}{3}{4}'.format(r'$\mathrm{predicted}$', '\n', r'$\mathrm{egress}$', '\n',
                                              r'$\mathrm{end}$'),
                     ha='left', va='top', fontsize=10)

            plt.suptitle('{0}{1}{2}'.format(r'$\mathbf{', target, '}$'), fontsize=20)
            plt.text(plt.xlim()[1], plt.ylim()[1], '{0}{1}{2}'.format(r'$', data_dates[set_number], '$'),
                     fontsize=12, ha='right', va='bottom')

            plt.subplot(4, 1, 4)
            plt.cla()
            plt.plot(self.results['detrended_output_series']['phase'][set_indices],
                     self.results['detrended_output_series']['residuals'][set_indices], 'ko', ms=2)
            plt.plot(self.results['detrended_output_series']['phase'][set_indices],
                     np.zeros_like(self.results['detrended_output_series']['phase'][set_indices]), 'r-')

            plt.ylim(- 5 * np.std(self.results['detrended_output_series']['residuals'][set_indices]),
                     5 * np.std(self.results['detrended_output_series']['residuals'][set_indices]))

            plt.xlabel(r'$\mathrm{phase}$', fontsize=15)
            plt.ylabel(r'$\mathrm{residuals}$', fontsize=15)

            plt.xlim(-x_max, x_max)

            plt.text(plt.xlim()[0] + 0.02 * (plt.xlim()[-1] - plt.xlim()[0]),
                     plt.ylim()[0] + 0.07 * (plt.ylim()[-1] - plt.ylim()[0]),
                     r'$\mathrm{rms}_\mathrm{res} = %.1e$' %
                     np.std(self.results['detrended_output_series']['residuals'][set_indices]),
                     fontsize=10)

            plt.subplots_adjust(left=0.15, right=0.975, bottom=0.12, top=0.9, hspace=0.0)

            plt.savefig(os.path.join(os.path.split(export_file)[0],
                                     'set_{0}_{1}'.format(str(set_number + 1), os.path.split(export_file)[1])),
                        transparent=True)
            if return_plot:
                return [plt.figure(ff + 1) for ff in range(self.total_sets)]
            else:
                plt.close('all')

    def save_models(self, export_file):

        for set_number in range(self.total_sets):

            self.results = {ff: self.results[ff] for ff in self.results}

            set_indices = np.where(self.data_set_number == set_number)

            np.savetxt(os.path.join(os.path.split(export_file)[0],
                                    'set_{0}_{1}'.format(str(set_number + 1), os.path.split(export_file)[1])),
                       np.swapaxes([self.results['input_series']['hjd'][set_indices],
                                    self.results['output_series']['phase'][set_indices],
                                    self.results['input_series']['value'][set_indices],
                                    self.results['input_series']['error'][set_indices],
                                    self.results['output_series']['model'][set_indices],
                                    self.results['output_series']['residuals'][set_indices]
                                    ], 0, 1))

    def save_detrended_models(self, export_file):

        for set_number in range(self.total_sets):

            self.results = {ff: self.results[ff] for ff in self.results}

            set_indices = np.where(self.data_set_number == set_number)

            np.savetxt(os.path.join(os.path.split(export_file)[0],
                                    'set_{0}_{1}'.format(str(set_number + 1), os.path.split(export_file)[1])),
                       np.swapaxes([self.results['detrended_input_series']['hjd'][set_indices],
                                    self.results['detrended_output_series']['phase'][set_indices],
                                    self.results['detrended_input_series']['value'][set_indices],
                                    self.results['detrended_input_series']['error'][set_indices],
                                    self.results['detrended_output_series']['model'][set_indices],
                                    self.results['detrended_output_series']['residuals'][set_indices]
                                    ], 0, 1))
