__all__ = ['TransitAndPolyFitting']


import os
import pickle

import numpy as np

import matplotlib.pyplot as plt

from exoplanet_orbit import *
from transit_flux_drop import *
from emcee_fitting import *
from transit_duration import *


class TransitAndPolyFitting():

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

            self.names.append('N' + str(set_number))
            self.print_names.append('N_' + str(set_number))
            self.initial.append(np.mean(set_arrays[1]))
            self.limits1.append(0.5 * np.mean(set_arrays[1]))
            self.limits2.append(1.5 * np.mean(set_arrays[1]))

            self.names.append('L' + str(set_number))
            self.print_names.append('L_' + str(set_number))
            self.initial.append(0)
            if self.fit_first_order:
                self.limits1.append(-max_limit)
                self.limits2.append(max_limit)
            else:
                self.limits1.append(np.nan)
                self.limits2.append(np.nan)

            self.names.append('Q' + str(set_number))
            self.print_names.append('Q_' + str(set_number))
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
                raise RuntimeError('Improper value for ' + self.names[var])

            if limits[var] is False:
                self.limits1.append(np.nan)
                self.limits2.append(np.nan)

            elif limits[var] is None:
                self.limits1.append(np.nan)
                self.limits2.append(np.nan)

            else:
                try:
                    if len(np.array(limits[var])) != 2:
                        raise RuntimeError('Improper limits for ' + self.names[var])
                except:
                    raise RuntimeError('Improper limits for ' + self.names[var])

                if self.initial[var] < np.array(limits[var])[0] or self.initial[var] > np.array(limits[var])[1]:
                    raise RuntimeError('Initial value for ' + self.names[var] + ' is outside the range of the prior')
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

        pickle.dump(self.results, open(export_file, 'w'))

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
            data_dates = map(str, ['set_' + str(ff) for ff in range(1, self.total_sets + 1)])

        for set_number in range(self.total_sets):

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

            x_max = max(self.results['output_series']['phase'][set_indices] +
                        0.05 * (max(self.results['output_series']['phase'][set_indices]) -
                                min(self.results['output_series']['phase'][set_indices])))
            plt.xlim(-x_max, x_max)
            plt.tick_params(labelbottom='off')

            rpstr = r'$R_\mathrm{p}/R_* = ' + self.results['parameters']['rp']['print_value'] + \
                    '_{-' + self.results['parameters']['rp']['print_m_error'] + '}' + \
                    '^{+' + self.results['parameters']['rp']['print_p_error'] + '}$'
            mtstr = r'$T_\mathrm{HJD} = ' + self.results['parameters']['mt']['print_value'] + \
                    '_{-' + self.results['parameters']['mt']['print_m_error'] + '}' + \
                    '^{+' + self.results['parameters']['mt']['print_p_error'] + '}$'

            plt.text(plt.xlim()[0] + 0.5 * (plt.xlim()[-1] - plt.xlim()[0]),
                     plt.ylim()[0] + 0.07 * (plt.ylim()[-1] - plt.ylim()[0]),
                     rpstr + '\n' + mtstr, ha='center', va='center', fontsize=10)

            plt.axvline((ingress - mt) / period, 0.3, 1.0, ls='--', c='k', lw=0.75)
            plt.text((ingress - mt) / period, plt.ylim()[0] + 0.3 * (plt.ylim()[1] - plt.ylim()[0]),
                     r'$\mathrm{predicted}$' + '\n' + r'$\mathrm{ingress}$' + '\n' + r'$\mathrm{start}$',
                     ha='right', va='top', fontsize=10)
            plt.axvline((egress - mt) / period, 0.3, 1.0, ls='--', c='k', lw=0.75)
            plt.text((egress - mt) / period, plt.ylim()[0] + 0.3 * (plt.ylim()[1] - plt.ylim()[0]),
                     r'$\mathrm{predicted}$' + '\n' + r'$\mathrm{egress}$' + '\n' + r'$\mathrm{end}$',
                     ha='left', va='top', fontsize=10)

            plt.suptitle(r'$\mathbf{' + target + '}$', fontsize=20)
            plt.text(plt.xlim()[1], plt.ylim()[1],
                     r'$' + data_dates[set_number] + '$', fontsize=12, ha='right', va='bottom')

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

            plt.subplots_adjust(hspace=0.0)

            plt.savefig(os.path.join(os.path.split(export_file)[0],
                                     'set_' + str(set_number + 1) + '_' + os.path.split(export_file)[1]),
                        bbox_inches='tight', transparent=True)
            plt.close('all')

    def plot_detrended_models(self, export_file, target=None, data_dates=None, return_plot=False):

        if target is None:
            target = ' '

        if data_dates is None:
            data_dates = map(str, ['set_' + str(ff) for ff in range(1, self.total_sets + 1)])

        for set_number in range(self.total_sets):

            plt.figure(set_number + 1)

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

            x_max = max(self.results['detrended_output_series']['phase'][set_indices] +
                        0.05 * (max(self.results['detrended_output_series']['phase'][set_indices]) -
                                min(self.results['detrended_output_series']['phase'][set_indices])))
            plt.xlim(-x_max, x_max)
            plt.tick_params(labelbottom='off')

            rpstr = r'$R_\mathrm{p}/R_* = ' + self.results['parameters']['rp']['print_value'] + \
                    '_{-' + self.results['parameters']['rp']['print_m_error'] + '}' + \
                    '^{+' + self.results['parameters']['rp']['print_p_error'] + '}$'
            mtstr = r'$T_\mathrm{HJD} = ' + self.results['parameters']['mt']['print_value'] + \
                    '_{-' + self.results['parameters']['mt']['print_m_error'] + '}' + \
                    '^{+' + self.results['parameters']['mt']['print_p_error'] + '}$'

            plt.text(plt.xlim()[0] + 0.5 * (plt.xlim()[-1] - plt.xlim()[0]),
                     plt.ylim()[0] + 0.07 * (plt.ylim()[-1] - plt.ylim()[0]),
                     rpstr + '\n' + mtstr, ha='center', va='center', fontsize=10)

            plt.axvline((ingress - mt) / period, 0.3, 1.0, ls='--', c='k', lw=0.75)
            plt.text((ingress - mt) / period, plt.ylim()[0] + 0.3 * (plt.ylim()[1] - plt.ylim()[0]),
                     r'$\mathrm{predicted}$' + '\n' + r'$\mathrm{ingress}$' + '\n' + r'$\mathrm{start}$',
                     ha='right', va='top', fontsize=10)
            plt.axvline((egress - mt) / period, 0.3, 1.0, ls='--', c='k', lw=0.75)
            plt.text((egress - mt) / period, plt.ylim()[0] + 0.3 * (plt.ylim()[1] - plt.ylim()[0]),
                     r'$\mathrm{predicted}$' + '\n' + r'$\mathrm{egress}$' + '\n' + r'$\mathrm{end}$',
                     ha='left', va='top', fontsize=10)

            plt.suptitle(r'$\mathbf{' + target + '}$', fontsize=20)
            plt.text(plt.xlim()[1], plt.ylim()[1],
                     r'$' + data_dates[set_number] + '$', fontsize=12, ha='right', va='bottom')

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

            plt.subplots_adjust(hspace=0.0)

            plt.savefig(os.path.join(os.path.split(export_file)[0],
                                     'set_' + str(set_number + 1) + '_' + os.path.split(export_file)[1]),
                        bbox_inches='tight', transparent=True)
            if return_plot:
                return [plt.figure(ff + 1) for ff in range(self.total_sets)]
            else:
                plt.close('all')

    def save_models(self, export_file):

        for set_number in range(self.total_sets):

            self.results = {ff: self.results[ff] for ff in self.results}

            set_indices = np.where(self.data_time == self.data[set_number][0])

            np.savetxt(os.path.join(os.path.split(export_file)[0],
                                    'set_' + str(set_number + 1) + '_' + os.path.split(export_file)[1]),
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
                                    'set_' + str(set_number + 1) + '_' + os.path.split(export_file)[1]),
                       np.swapaxes([self.results['detrended_input_series']['hjd'][set_indices],
                                    self.results['detrended_output_series']['phase'][set_indices],
                                    self.results['detrended_input_series']['value'][set_indices],
                                    self.results['detrended_input_series']['error'][set_indices],
                                    self.results['detrended_output_series']['model'][set_indices],
                                    self.results['detrended_output_series']['residuals'][set_indices]
                                    ], 0, 1))