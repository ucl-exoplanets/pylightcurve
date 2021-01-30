__all__ = ['Planet']

import os
import numpy as np

from pylightcurve.errors import *
from pylightcurve.models.exoplanet_lc import planet_orbit, planet_star_projected_distance, planet_phase, \
    transit, transit_integrated, \
    transit_duration, transit_depth, eclipse, eclipse_integrated, eclipse_mid_time, eclipse_duration, eclipse_depth,\
    exotethys, fp_over_fs, _get_filter
from pylightcurve.analysis.optimisation import EmceeFitting
from pylightcurve.processes.files import open_dict, save_dict
from pylightcurve.plots.plots_fitting import plot_transit_fitting_models
from pylightcurve.spacetime.angles import Degrees, _request_angle
from pylightcurve.spacetime.targets import FixedTarget


class Filter:

    def __init__(self, rp_over_rs, ldc1, ldc2, ldc3, ldc4, fp_over_fs):

        self.rp_over_rs = rp_over_rs
        self.fp_over_fs = fp_over_fs
        self.limb_darkening_coefficients = [ldc1, ldc2, ldc3, ldc4]


class Planet:

    def __init__(self, name, ra, dec, stellar_logg, stellar_temperature, stellar_metallicity,
                 rp_over_rs, period, sma_over_rs, eccentricity, inclination,
                 periastron, mid_time, mid_time_format, ww=0,
                 albedo=0.15, emissivity=1.0,
                 ldc_method='claret', ldc_stellar_model='phoenix'):

        if isinstance(ra, float):
            ra = Degrees(ra)
        else:
            _request_angle(ra)

        if isinstance(dec, float):
            dec = Degrees(dec)
        else:
            _request_angle(dec)

        if isinstance(inclination, float):
            pass
        else:
            _request_angle(inclination)
            inclination = inclination.deg()

        if isinstance(periastron, float):
            pass
        else:
            _request_angle(periastron)
            periastron = periastron.deg()

        self.name = name
        self.target = FixedTarget(ra, dec)

        self.stellar_logg = stellar_logg
        self.stellar_temperature = stellar_temperature
        self.stellar_metallicity = stellar_metallicity
        self.ldc_method = ldc_method
        self.ldc_stellar_model = ldc_stellar_model

        self.rp_over_rs = rp_over_rs
        self.period = period
        self.sma_over_rs = sma_over_rs
        self.eccentricity = eccentricity
        self.inclination = inclination
        self.periastron = periastron
        self.mid_time = self.target.convert_to_bjd_tdb(mid_time, mid_time_format)
        self.eclipse_mid_time = eclipse_mid_time(self.period, self.sma_over_rs, self.eccentricity, self.inclination,
                                                 self.periastron, self.mid_time)
        self.ww = ww
        self.albedo = albedo
        self.emissivity = emissivity

        self.filters = {}

        self.observations = {}

    # filter-dependent values

    def filter(self, filter_name):

        if filter_name in self.filters:
            return self.filters[filter_name]

        else:
            try:
                _get_filter(filter_name)
                ldc1, ldc2, ldc3, ldc4 = exotethys(self.stellar_logg, self.stellar_temperature, self.stellar_metallicity, filter_name,
                                                   method=self.ldc_method, stellar_model=self.ldc_stellar_model)
                fp = fp_over_fs(self.rp_over_rs, self.sma_over_rs, self.albedo, self.emissivity, self.stellar_temperature,
                                filter_name)
                print('Fp/Fs estimated using A={0}, e={1} for filter {2}.'.format(self.albedo, self.emissivity,
                                                                                  filter_name))

                self.filters[filter_name] = Filter(self.rp_over_rs, ldc1, ldc2, ldc3, ldc4, fp)

                return self.filters[filter_name]

            except PyLCInputError:
                raise PyLCInputError('Filter not available, you need to add {0} to the existing filters first.'.format
                                     (filter_name))

    def add_filter(self, filter_name, rp_over_rs, ldc1, ldc2, ldc3, ldc4, fp):

        self.filters[filter_name] = Filter(rp_over_rs, ldc1, ldc2, ldc3, ldc4, fp)

    # planet calculations

    def planet_orbit(self, time, time_format):
        time = self.target.convert_to_bjd_tdb(time, time_format)
        return planet_orbit(self.period, self.sma_over_rs, self.eccentricity, self.inclination, self.periastron,
                            self.mid_time, time, ww=self.ww)

    def planet_star_projected_distance(self, time, time_format):
        time = self.target.convert_to_bjd_tdb(time, time_format)
        return planet_star_projected_distance(self.period, self.sma_over_rs, self.eccentricity, self.inclination,
                                              self.periastron, self.mid_time, time)

    def planet_phase(self, time, time_format):
        time = np.array([self.target.convert_to_bjd_tdb(ff, time_format) for ff in time])
        return planet_phase(self.period, self.mid_time, time)

    def transit(self, time, time_format, filter_name, precision=3):
        filter_data = self.filter(filter_name)
        time = self.target.convert_to_bjd_tdb(time, time_format)

        return transit(filter_data.limb_darkening_coefficients,
                       filter_data.rp_over_rs, self.period, self.sma_over_rs, self.eccentricity,
                       self.inclination, self.periastron, self.mid_time, time,
                       method=self.ldc_method, precision=precision)

    def transit_integrated(self, time, time_format, exp_time, time_stamp, filter_name, max_sub_exp_time=10, precision=3):

        if time_stamp == 'start':
            time = np.array(time) + 0.5 * exp_time / (60.0 * 60.0 * 24.0)
        elif time_stamp == 'mid':
            time = np.array(time)
        elif time_stamp == 'end':
            time = np.array(time) - 0.5 * exp_time / (60.0 * 60.0 * 24.0)
        else:
            raise PyLCInputError(
                'Not acceptable time stamp {0}. Please choose between "mid", "start", "end".'.format(time_stamp))

        filter_data = self.filter(filter_name)
        time = self.target.convert_to_bjd_tdb(time, time_format)

        return transit_integrated(filter_data.limb_darkening_coefficients,
                                  filter_data.rp_over_rs, self.period, self.sma_over_rs, self.eccentricity,
                                  self.inclination, self.periastron, self.mid_time, time, exp_time=exp_time,
                                  max_sub_exp_time=max_sub_exp_time,
                                  method=self.ldc_method, precision=precision)

    def transit_duration(self, filter_name):
        filter_data = self.filter(filter_name)
        return transit_duration(filter_data.rp_over_rs,
                                self.period, self.sma_over_rs, self.eccentricity, self.inclination, self.periastron)

    def transit_depth(self, filter_name, precision=6):
        filter_data = self.filter(filter_name)
        return transit_depth(filter_data.limb_darkening_coefficients,
                             filter_data.rp_over_rs, self.period, self.sma_over_rs, self.eccentricity,
                             self.inclination, self.periastron,
                             method=self.ldc_method, precision=precision)

    def eclipse(self, time, time_format, filter_name, precision=3):
        filter_data = self.filter(filter_name)
        time = self.target.convert_to_bjd_tdb(time, time_format)

        return eclipse(filter_data.fp_over_fs, filter_data.rp_over_rs,
                       self.period, self.sma_over_rs, self.eccentricity,
                       self.inclination, self.periastron, self.eclipse_mid_time, time, precision=precision)

    def eclipse_integrated(self, time, time_format, exp_time, time_stamp, filter_name, max_sub_exp_time=10, precision=3):

        if time_stamp == 'start':
            time = np.array(time) + 0.5 * exp_time / (60.0 * 60.0 * 24.0)
        elif time_stamp == 'mid':
            time = np.array(time)
        elif time_stamp == 'end':
            time = np.array(time) - 0.5 * exp_time / (60.0 * 60.0 * 24.0)
        else:
            raise PyLCInputError(
                'Not acceptable time stamp {0}. Please choose between "mid", "start", "end".'.format(time_stamp))

        filter_data = self.filter(filter_name)
        time = self.target.convert_to_bjd_tdb(time, time_format)

        return eclipse_integrated(filter_data.fp_over_fs, filter_data.rp_over_rs,
                                  self.period, self.sma_over_rs, self.eccentricity,
                                  self.inclination, self.periastron, self.eclipse_mid_time, time, exp_time=exp_time,
                                  max_sub_exp_time=max_sub_exp_time, precision=precision)

    def eclipse_duration(self, filter_name):
        filter_data = self.filter(filter_name)
        return eclipse_duration(filter_data.rp_over_rs, self.period, self.sma_over_rs, self.eccentricity,
                                self.inclination, self.periastron)

    def eclipse_depth(self, filter_name, precision=6):
        filter_data = self.filter(filter_name)
        return eclipse_depth(filter_data.fp_over_fs, filter_data.rp_over_rs, self.period,
                             self.sma_over_rs, self.eccentricity, self.inclination, self.periastron,
                             precision=precision)

    # data fitting

    def add_observation(self, time, time_format, exp_time, time_stamp, flux, flux_unc, flux_format, filter_name):

        filter_data = self.filter(filter_name)

        original_data = {
            'time': time,
            'time_stamp': time_stamp,
            'time_format': time_format,
            'flux': flux,
            'flux_format': flux_format,
            'flux_unc': flux_unc,
            'exp_time': exp_time,
            'filter': filter_name
        }

        if time_stamp == 'start':
            time = np.array(time) + 0.5 * exp_time / (60.0 * 60.0 * 24.0)
        elif time_stamp == 'mid':
            time = np.array(time)
        elif time_stamp == 'end':
            time = np.array(time) - 0.5 * exp_time / (60.0 * 60.0 * 24.0)
        else:
            raise PyLCInputError(
                'Not acceptable time stamp {0}. Please choose between "mid", "start", "end".'.format(time_stamp))

        date = self.target.convert_to_jd(time[0], time_format)

        time = np.array([self.target.convert_to_bjd_tdb(ff, time_format) for ff in time])

        if flux_format == 'mag':
            flux_unc = np.abs(np.array(flux_unc) * (-0.921034 * np.exp(0.921034 * flux[0] - 0.921034 * np.array(flux))))
            flux = 10 ** ((flux[0] - np.array(flux)) / 2.5)
        elif flux_format == 'flux':
            flux = np.array(flux)
            flux_unc = np.array(flux_unc)
        else:
            raise PyLCInputError('Not acceptable flux format {0}. Please choose between "flux" and '
                                 '"mag".'.format(flux_format))

        obs_id = 0
        while obs_id in self.observations:
            obs_id += 1

        check_transit = (np.mean(time) - self.mid_time) / self.period
        check_transit = abs(check_transit - int(check_transit))

        check_eclipse = (np.mean(time) - self.eclipse_mid_time) / self.period
        check_eclipse = abs(check_eclipse - int(check_eclipse))

        if check_transit < check_eclipse:
            observation_type = 'transit'
            epoch = int(round((np.mean(time) - self.mid_time) / self.period, 0))
        else:
            observation_type = 'eclipse'
            epoch = int(round((np.mean(time) - self.eclipse_mid_time) / self.period, 0))

        self.observations[obs_id] = {
            'target': self.name,
            'time': time,
            'dtime': time - time[0],
            'flux': flux,
            'flux_unc': flux_unc,
            'exp_time': exp_time,
            'filter': filter_name,
            'epoch': epoch,
            'date': date,
            'observation_type': observation_type,
            'original_data': original_data
        }

    def transit_fitting(self, output_folder,

                        max_sub_exp_time=10, precision=3,

                        detrending_order=2,

                        iterations=130000, walkers=200, burn_in=30000,

                        fit_ldc1=False, fit_ldc2=False, fit_ldc3=False, fit_ldc4=False,

                        fit_rp_over_rs=True, fit_individual_rp_over_rs=True,

                        fit_sma_over_rs=False, fit_inclination=False,

                        fit_mid_time=True, fit_period=False,
                        fit_individual_times=True,

                        fit_ldc_limits=[0.0, 1.0],
                        fit_rp_over_rs_limits=[0.5, 2.0],
                        fit_sma_over_rs_limits=[0.5, 2.0],
                        fit_inclination_limits=[70.0, 90.0],
                        fit_mid_time_limits=[-0.2, 0.2],
                        fit_period_limits=[0.8, 1.2],

                        counter='MCMC'):

        if not os.path.isdir(output_folder):
            os.mkdir(output_folder)

        parameters_map = [[] for observation in self.observations]

        names = []
        print_names = []
        limits1 = []
        limits2 = []
        initial = []

        # de-trending parameters

        for observation_num, observation in enumerate(self.observations):

            max_limit = (10 * (max(self.observations[observation]['flux']) - min(self.observations[observation]['flux'])) /
                         (max(self.observations[observation]['flux']) - min(self.observations[observation]['flux'])) / np.mean(self.observations[observation]['flux']))

            names.append('N_{0}'.format(observation_num + 1))
            print_names.append('N_{0}'.format(observation_num + 1))
            initial.append(np.mean(self.observations[observation]['flux']))
            limits1.append(0.5 * np.mean(self.observations[observation]['flux']))
            limits2.append(1.5 * np.mean(self.observations[observation]['flux']))

            parameters_map[observation_num].append(len(names) - 1)

            names.append('L_{0}'.format(observation_num + 1))
            print_names.append('L_{0}'.format(observation_num + 1))
            initial.append(0)
            if detrending_order >= 1:
                limits1.append(-max_limit)
                limits2.append(max_limit)
            else:
                limits1.append(np.nan)
                limits2.append(np.nan)

            parameters_map[observation_num].append(len(names) - 1)

            names.append('Q_{0}'.format(observation_num + 1))
            print_names.append('Q_{0}'.format(observation_num + 1))
            initial.append(0)
            if detrending_order == 2:
                limits1.append(-5)
                limits2.append(5)
            else:
                limits1.append(np.nan)
                limits2.append(np.nan)

            parameters_map[observation_num].append(len(names) - 1)

        # limb-darkening and rp_over_rs parameters

        unique_filters = []
        for observation in self.observations:
            if self.observations[observation]['filter'] not in unique_filters:
                unique_filters.append(self.observations[observation]['filter'])

        for phot_filter in unique_filters:

            filter_data = self.filter(phot_filter)
            ldc1, ldc2, ldc3, ldc4 = filter_data.limb_darkening_coefficients
            rp_over_rs = filter_data.rp_over_rs

            names.append('LDC1_{0}'.format(phot_filter))
            print_names.append('LDC1_{0}'.format(phot_filter))
            initial.append(ldc1)
            if fit_ldc1:
                limits1.append(fit_ldc_limits[0])
                limits2.append(fit_ldc_limits[1])
            else:
                limits1.append(np.nan)
                limits2.append(np.nan)

            for observation_num, observation in enumerate(self.observations):
                if self.observations[observation]['filter'] == phot_filter:
                    parameters_map[observation_num].append(len(names) - 1)

            names.append('LDC2_{0}'.format(phot_filter))
            print_names.append('LDC2_{0}'.format(phot_filter))
            initial.append(ldc2)
            if fit_ldc2:
                limits1.append(fit_ldc_limits[0])
                limits2.append(fit_ldc_limits[1])
            else:
                limits1.append(np.nan)
                limits2.append(np.nan)

            for observation_num, observation in enumerate(self.observations):
                if self.observations[observation]['filter'] == phot_filter:
                    parameters_map[observation_num].append(len(names) - 1)

            names.append('LDC3_{0}'.format(phot_filter))
            print_names.append('LDC3_{0}'.format(phot_filter))
            initial.append(ldc3)
            if fit_ldc3:
                limits1.append(fit_ldc_limits[0])
                limits2.append(fit_ldc_limits[1])
            else:
                limits1.append(np.nan)
                limits2.append(np.nan)

            for observation_num, observation in enumerate(self.observations):
                if self.observations[observation]['filter'] == phot_filter:
                    parameters_map[observation_num].append(len(names) - 1)

            names.append('LDC4_{0}'.format(phot_filter))
            print_names.append('LDC4_{0}'.format(phot_filter))
            initial.append(ldc4)
            if fit_ldc4:
                limits1.append(fit_ldc_limits[0])
                limits2.append(fit_ldc_limits[1])
            else:
                limits1.append(np.nan)
                limits2.append(np.nan)

            for observation_num, observation in enumerate(self.observations):
                if self.observations[observation]['filter'] == phot_filter:
                    parameters_map[observation_num].append(len(names) - 1)

            if fit_individual_rp_over_rs:

                names.append('rp_{0}'.format(phot_filter))
                print_names.append('(R_\mathrm{p}/R_*)_{' + phot_filter + '}')
                initial.append(rp_over_rs)
                if fit_rp_over_rs:
                    limits1.append(rp_over_rs * fit_rp_over_rs_limits[0])
                    limits2.append(rp_over_rs * fit_rp_over_rs_limits[1])
                else:
                    limits1.append(np.nan)
                    limits2.append(np.nan)

                for observation_num, observation in enumerate(self.observations):
                    if self.observations[observation]['filter'] == phot_filter:
                        parameters_map[observation_num].append(len(names) - 1)

        if not fit_individual_rp_over_rs:

            names.append('rp')
            print_names.append('(R_\mathrm{p}/R_*)')
            initial.append(self.rp_over_rs)
            if fit_rp_over_rs:
                limits1.append(self.rp_over_rs * fit_rp_over_rs_limits[0])
                limits2.append(self.rp_over_rs * fit_rp_over_rs_limits[1])
            else:
                limits1.append(np.nan)
                limits2.append(np.nan)

            for observation_num, observation in enumerate(self.observations):
                parameters_map[observation_num].append(len(names) - 1)

        # orbital parameters

        names.append('P')
        print_names.append('P')
        initial.append(self.period)
        if fit_period:
            limits1.append(self.period * fit_period_limits[0])
            limits2.append(self.period * fit_period_limits[1])
        else:
            limits1.append(np.nan)
            limits2.append(np.nan)

        for observation_num, observation in enumerate(self.observations):
            parameters_map[observation_num].append(len(names) - 1)

        names.append('a')
        print_names.append('a/R_*')
        initial.append(self.sma_over_rs)
        if fit_sma_over_rs:
            limits1.append(self.sma_over_rs * fit_sma_over_rs_limits[0])
            limits2.append(self.sma_over_rs * fit_sma_over_rs_limits[1])
        else:
            limits1.append(np.nan)
            limits2.append(np.nan)

        for observation_num, observation in enumerate(self.observations):
            parameters_map[observation_num].append(len(names) - 1)

        names.append('e')
        print_names.append('e')
        initial.append(self.eccentricity)
        limits1.append(np.nan)
        limits2.append(np.nan)

        for observation_num, observation in enumerate(self.observations):
            parameters_map[observation_num].append(len(names) - 1)

        names.append('i')
        print_names.append('i')
        initial.append(self.inclination)
        if fit_inclination:
            limits1.append(fit_inclination_limits[0])
            limits2.append(fit_inclination_limits[1])
        else:
            limits1.append(np.nan)
            limits2.append(np.nan)

        for observation_num, observation in enumerate(self.observations):
            parameters_map[observation_num].append(len(names) - 1)

        names.append('w')
        print_names.append('\omega')
        initial.append(self.periastron)
        limits1.append(np.nan)
        limits2.append(np.nan)

        for observation_num, observation in enumerate(self.observations):
            parameters_map[observation_num].append(len(names) - 1)

        # time parameters

        test_epochs = []
        test_epochs_weights = []

        for observation in self.observations:
            test_epochs.append((self.observations[observation]['time'] - self.mid_time) / self.period)
            norm_errors = self.observations[observation]['flux_unc'] / self.observations[observation]['flux']
            test_epochs_weights.append(1 / (norm_errors * norm_errors))

        test_epochs = np.concatenate(test_epochs)
        test_epochs_weights = np.concatenate(test_epochs_weights)

        new_epoch = np.round(np.sum(test_epochs * test_epochs_weights) / np.sum(test_epochs_weights), 0)
        new_mid_time = self.mid_time + new_epoch * self.period

        for observation in self.observations:
            self.observations[observation]['epoch'] = int(round((np.mean(self.observations[observation]['time']) -new_mid_time) / self.period, 0))

        unique_epochs = []
        for observation in self.observations:
            if self.observations[observation]['epoch'] not in unique_epochs:
                unique_epochs.append(self.observations[observation]['epoch'])

        if not fit_individual_times:

            names.append('T_0')
            print_names.append('T_0')
            initial.append(new_mid_time)
            if fit_mid_time:
                limits1.append(new_mid_time + fit_mid_time_limits[0])
                limits2.append(new_mid_time + fit_mid_time_limits[1])
            else:
                limits1.append(np.nan)
                limits2.append(np.nan)

            for observation_num, observation in enumerate(self.observations):
                parameters_map[observation_num].append(len(names) - 1)

        else:

            if fit_period:
                raise PyLCInputError('Period and individual mid times cannot be fitted simultaneously.')

            for epoch in unique_epochs:

                names.append('T_mid_{0}'.format(epoch))
                print_names.append('T_\mathrm{mid_' + str(epoch) + '}')
                initial.append(new_mid_time + epoch * self.period)
                if fit_mid_time:
                    limits1.append(new_mid_time + epoch * self.period + fit_mid_time_limits[0])
                    limits2.append(new_mid_time + epoch * self.period + fit_mid_time_limits[1])
                else:
                    limits1.append(np.nan)
                    limits2.append(np.nan)

                for observation_num, observation in enumerate(self.observations):
                    if self.observations[observation]['epoch'] == epoch:
                        parameters_map[observation_num].append(len(names) - 1)

        model_time = np.array([])
        model_flux = np.array([])
        model_flux_unc = np.array([])
        for observation in self.observations:
            model_time = np.append(model_time, self.observations[observation]['time'])
            model_flux = np.append(model_flux, self.observations[observation]['flux'])
            model_flux_unc = np.append(model_flux_unc, self.observations[observation]['flux_unc'])

        parameters_map = np.array(parameters_map)

        def detrend_model(model_time, *model_variables):

            model = []

            for observation_num, observation in enumerate(self.observations):

                detrend_zero,  detrend_one, detrend_two, ldc1, ldc2, ldc3, ldc4, r, p, a, e, i, w, mt = \
                    np.array(model_variables)[parameters_map[observation_num]]

                deltatime = self.observations[observation]['dtime']

                model.append(detrend_zero * (1 + detrend_one * deltatime + detrend_two * deltatime * deltatime))

            return np.concatenate(model)

        def full_model(model_time, *model_variables):

            model_variables = np.array(model_variables)

            model = []

            for observation_num, observation in enumerate(self.observations):
                detrend_zero, detrend_one, detrend_two, ldc1, ldc2, ldc3, ldc4, r, p, a, e, i, w, mt = \
                    model_variables[parameters_map[observation_num]]

                deltatime = self.observations[observation]['dtime']

                model1 = detrend_zero * (1 + detrend_one * deltatime + detrend_two * deltatime * deltatime)

                model.append(model1 * transit_integrated(
                    [ldc1, ldc2, ldc3, ldc4], r, p, a, e, i, w, mt,
                    time_array=self.observations[observation]['time'],
                    exp_time=self.observations[observation]['exp_time'],
                    max_sub_exp_time=max_sub_exp_time,
                    method=self.ldc_method,
                    precision=precision
                ))

            return np.concatenate(model)

        fitting = EmceeFitting(model_time, model_flux, model_flux_unc,
                               full_model, initial, limits1, limits2, walkers,
                               iterations, burn_in,
                               data_x_name='time', data_y_name='flux', data_x_print_name='Time',
                               data_y_print_name='Relative Flux',
                               parameters_names=names, parameters_print_names=print_names,
                               counter=counter)

        fitting.run_mcmc()
        results = fitting.results

        # TODO
        # save observations separately and calculate some statistics individually

        results['settings'] = {}
        results['settings']['max_sub_exp_time'] = max_sub_exp_time
        results['settings']['precision'] = precision
        results['settings']['detrending_order'] = detrending_order
        results['settings']['iterations'] = iterations
        results['settings']['walkers'] = walkers
        results['settings']['burn_in'] = burn_in
        results['settings']['fit_ldc1'] = fit_ldc1
        results['settings']['fit_ldc2'] = fit_ldc2
        results['settings']['fit_ldc3'] = fit_ldc3
        results['settings']['fit_ldc4'] = fit_ldc4
        results['settings']['fit_rp_over_rs'] = fit_rp_over_rs
        results['settings']['force_same_rp_over_rs'] = fit_individual_rp_over_rs
        results['settings']['fit_sma_over_rs'] = fit_sma_over_rs
        results['settings']['fit_inclination'] = fit_inclination
        results['settings']['fit_mid_time'] = fit_mid_time
        results['settings']['fit_period'] = fit_period
        results['settings']['fit_individual_times'] = fit_individual_times
        results['settings']['fit_ldc_limits'] = fit_ldc_limits
        results['settings']['fit_rp_over_rs_limits'] = fit_rp_over_rs_limits
        results['settings']['fit_sma_over_rs_limits'] = fit_sma_over_rs_limits
        results['settings']['fit_inclination_limits'] = fit_inclination_limits
        results['settings']['fit_mid_time_limits'] = fit_mid_time_limits
        results['settings']['fit_period_limits'] = fit_period_limits
        results['settings']['filter_map'] = self.filters

        results['detrended_input_series'] = {
            'time': results['input_series']['time'],
            'flux': results['input_series']['flux'] / detrend_model(model_time, *results['parameters_final']),
            'flux_unc': results['input_series']['flux_unc'] / detrend_model(model_time, *results['parameters_final'])}

        results['detrended_output_series'] = {
            'model': results['output_series']['model'] / detrend_model(model_time, *results['parameters_final']),
            'residuals': (results['output_series']['residuals']
                          / detrend_model(model_time, *results['parameters_final']))}

        results['detrended_statistics'] = {
            'res_mean': np.mean(results['detrended_output_series']['residuals']),
            'res_std': np.std(results['detrended_output_series']['residuals']),
            'res_rms': np.sqrt(np.mean(results['detrended_output_series']['residuals'] ** 2))
        }

        observations_id_series = np.array([])
        for observation in self.observations:
            observations_id_series = np.append(observations_id_series,
                                               np.ones(len(self.observations[observation]['time'])) * observation)

        for observation in self.observations:

            id_series = np.where(observations_id_series == observation)

            self.observations[observation]['model'] = results['output_series']['model'][id_series]
            self.observations[observation]['residuals'] = results['output_series']['residuals'][id_series]

            res_autocorr = np.correlate(self.observations[observation]['residuals'],
                                        self.observations[observation]['residuals'], mode='full')
            res_autocorr = res_autocorr[res_autocorr.size // 2:] / res_autocorr[res_autocorr.size // 2:][0]

            self.observations[observation]['res_autocorr'] = res_autocorr
            self.observations[observation]['res_max_autocorr'] = np.max(res_autocorr[1:])
            self.observations[observation]['res_mean'] = np.mean(self.observations[observation]['residuals'])
            self.observations[observation]['res_std'] = np.std(self.observations[observation]['residuals'])
            self.observations[observation]['res_rms'] = np.sqrt(np.mean(self.observations[observation]['residuals'] ** 2))
            self.observations[observation]['res_chi_sqr'] = np.sum(
                (self.observations[observation]['residuals'] ** 2) / (self.observations[observation]['flux_unc'] ** 2))
            self.observations[observation]['res_red_chi_sqr'] = (
                    self.observations[observation]['res_chi_sqr'] / (
                        len(self.observations[observation]['flux_unc']) - len(results['statistics']['corr_variables'])))

            self.observations[observation]['detrended_flux'] = results['detrended_input_series']['flux'][id_series]
            self.observations[observation]['detrended_flux_unc'] = results['detrended_input_series']['flux_unc'][id_series]
            self.observations[observation]['detrended_model'] = results['detrended_output_series']['model'][id_series]
            self.observations[observation]['detrended_residuals'] = results['detrended_output_series']['residuals'][id_series]

            res_autocorr = np.correlate(self.observations[observation]['detrended_residuals'],
                                        self.observations[observation]['detrended_residuals'], mode='full')
            res_autocorr = res_autocorr[res_autocorr.size // 2:] / res_autocorr[res_autocorr.size // 2:][0]

            self.observations[observation]['detrended_res_autocorr'] = res_autocorr
            self.observations[observation]['detrended_res_max_autocorr'] = np.max(res_autocorr[1:])
            self.observations[observation]['detrended_res_mean'] = np.mean(self.observations[observation]['detrended_residuals'])
            self.observations[observation]['detrended_res_std'] = np.std(self.observations[observation]['detrended_residuals'])
            self.observations[observation]['detrended_res_rms'] = np.sqrt(
                np.mean(self.observations[observation]['detrended_residuals'] ** 2))
            self.observations[observation]['detrended_res_chi_sqr'] = np.sum(
                (self.observations[observation]['detrended_residuals'] ** 2) / (self.observations[observation]['detrended_flux_unc'] ** 2))
            self.observations[observation]['detrended_res_red_chi_sqr'] = (
                    self.observations[observation]['detrended_res_chi_sqr'] / (
                    len(self.observations[observation]['detrended_flux_unc']) - len(results['statistics']['corr_variables'])))

            w = open(os.path.join(output_folder, 'diagnostics_dataset_{0}.txt'.format(observation + 1)), 'w')
            w.write('\n#Residuals:\n')
            w.write('#Mean: {0}\n'.format(self.observations[observation]['res_mean']))
            w.write('#STD: {0}\n'.format(self.observations[observation]['res_std']))
            w.write('#RMS: {0}\n'.format(self.observations[observation]['res_rms']))
            w.write('#Max auto-correlation: {0}\n'.format(self.observations[observation]['res_max_autocorr']))
            w.write('#Chi squared: {0}\n'.format(self.observations[observation]['res_chi_sqr']))
            w.write('#Reduced chi squared: {0}\n'.format(self.observations[observation]['res_red_chi_sqr']))

            w.write('\n\n#Detrended Residuals:\n')
            w.write('#Mean: {0}\n'.format(self.observations[observation]['detrended_res_mean']))
            w.write('#STD: {0}\n'.format(self.observations[observation]['detrended_res_std']))
            w.write('#RMS: {0}\n'.format(self.observations[observation]['detrended_res_rms']))
            w.write('#Max auto-correlation: {0}\n'.format(self.observations[observation]['detrended_res_max_autocorr']))
            w.write('#Chi squared: {0}\n'.format(self.observations[observation]['detrended_res_chi_sqr']))
            w.write('#Reduced chi squared: {0}\n'.format(self.observations[observation]['detrended_res_red_chi_sqr']))

            w.close()

        results['observations'] = self.observations

        save_dict(results, os.path.join(output_folder, 'results.pickle'))

        fitting.save_results(os.path.join(output_folder, 'results.txt'))

        fitting.plot_corner(os.path.join(output_folder, 'correlations.pdf'))

        fitting.plot_traces(os.path.join(output_folder, 'traces.pdf'))

        plot_transit_fitting_models(results, os.path.join(output_folder, 'lightcurves.pdf'))

    def eclipse_fitting(self, output_folder,

                        max_sub_exp_time=10, precision=3,

                        detrending_order=2,

                        iterations=130000, walkers=200, burn_in=30000,

                        fit_fp_over_fs=True, fit_individual_fp_over_fs=False,
                        fit_rp_over_rs=False, fit_individual_rp_over_rs=False,

                        fit_sma_over_rs=False, fit_inclination=False,

                        fit_mid_time=False, fit_period=False,
                        fit_individual_times=False,

                        fit_rp_over_rs_limits=[0.5, 2.0],
                        fit_fp_over_fs_limits=[0.001, 1000.0],
                        fit_sma_over_rs_limits=[0.5, 2.0],
                        fit_inclination_limits=[70.0, 90.0],
                        fit_mid_time_limits=[-0.2, 0.2],
                        fit_period_limits=[0.8, 1.2],

                        counter='MCMC'):

        if not os.path.isdir(output_folder):
            os.mkdir(output_folder)

        parameters_map = [[] for observation in self.observations]

        names = []
        print_names = []
        limits1 = []
        limits2 = []
        initial = []

        # de-trending parameters

        for observation_num, observation in enumerate(self.observations):

            max_limit = (10 * (max(self.observations[observation]['flux']) - min(self.observations[observation]['flux'])) /
                         (max(self.observations[observation]['flux']) - min(self.observations[observation]['flux'])) / np.mean(self.observations[observation]['flux']))

            names.append('N_{0}'.format(observation_num + 1))
            print_names.append('N_{0}'.format(observation_num + 1))
            initial.append(np.mean(self.observations[observation]['flux']))
            limits1.append(0.5 * np.mean(self.observations[observation]['flux']))
            limits2.append(1.5 * np.mean(self.observations[observation]['flux']))

            parameters_map[observation_num].append(len(names) - 1)

            names.append('L_{0}'.format(observation_num + 1))
            print_names.append('L_{0}'.format(observation_num + 1))
            initial.append(0)
            if detrending_order >= 1:
                limits1.append(-max_limit)
                limits2.append(max_limit)
            else:
                limits1.append(np.nan)
                limits2.append(np.nan)

            parameters_map[observation_num].append(len(names) - 1)

            names.append('Q_{0}'.format(observation_num + 1))
            print_names.append('Q_{0}'.format(observation_num + 1))
            initial.append(0)
            if detrending_order == 2:
                limits1.append(-5)
                limits2.append(5)
            else:
                limits1.append(np.nan)
                limits2.append(np.nan)

            parameters_map[observation_num].append(len(names) - 1)

        # limb-darkening and rp_over_rs parameters

        unique_filters = []
        for observation in self.observations:
            if self.observations[observation]['filter'] not in unique_filters:
                unique_filters.append(self.observations[observation]['filter'])

        fp_over_fs = 0.00001
        for phot_filter in unique_filters:

            filter_data = self.filter(phot_filter)
            fp_over_fs = filter_data.fp_over_fs

            if fit_individual_fp_over_fs:

                names.append('fp_{0}'.format(phot_filter))
                print_names.append('(F_\mathrm{p}/F_*)_{' + phot_filter + '}')
                initial.append(fp_over_fs)
                if fit_fp_over_fs:
                    limits1.append(fp_over_fs * fit_fp_over_fs_limits[0])
                    limits2.append(fp_over_fs * fit_fp_over_fs_limits[1])
                else:
                    limits1.append(np.nan)
                    limits2.append(np.nan)

                for observation_num, observation in enumerate(self.observations):
                    if self.observations[observation]['filter'] == phot_filter:
                        parameters_map[observation_num].append(len(names) - 1)

        if not fit_individual_fp_over_fs:

            names.append('fp')
            print_names.append('(F_\mathrm{p}/F_*)')
            initial.append(fp_over_fs)
            if fit_fp_over_fs:
                limits1.append(fp_over_fs * fit_fp_over_fs_limits[0])
                limits2.append(fp_over_fs * fit_fp_over_fs_limits[1])
            else:
                limits1.append(np.nan)
                limits2.append(np.nan)

            for observation_num, observation in enumerate(self.observations):
                parameters_map[observation_num].append(len(names) - 1)

        rp_over_rs = self.rp_over_rs
        for phot_filter in unique_filters:

            filter_data = self.filter(phot_filter)
            rp_over_rs = filter_data.rp_over_rs

            if fit_individual_rp_over_rs:

                names.append('rp_{0}'.format(phot_filter))
                print_names.append('(R_\mathrm{p}/R_*)_{' + phot_filter + '}')
                initial.append(rp_over_rs)
                if fit_rp_over_rs:
                    limits1.append(rp_over_rs * fit_rp_over_rs_limits[0])
                    limits2.append(rp_over_rs * fit_rp_over_rs_limits[1])
                else:
                    limits1.append(np.nan)
                    limits2.append(np.nan)

                for observation_num, observation in enumerate(self.observations):
                    if self.observations[observation]['filter'] == phot_filter:
                        parameters_map[observation_num].append(len(names) - 1)

        if not fit_individual_rp_over_rs:

            names.append('rp')
            print_names.append('(R_\mathrm{p}/R_*)')
            initial.append(rp_over_rs)
            if fit_rp_over_rs:
                limits1.append(rp_over_rs * fit_rp_over_rs_limits[0])
                limits2.append(rp_over_rs * fit_rp_over_rs_limits[1])
            else:
                limits1.append(np.nan)
                limits2.append(np.nan)

            for observation_num, observation in enumerate(self.observations):
                parameters_map[observation_num].append(len(names) - 1)

        # orbital parameters

        names.append('P')
        print_names.append('P')
        initial.append(self.period)
        if fit_period:
            limits1.append(self.period * fit_period_limits[0])
            limits2.append(self.period * fit_period_limits[1])
        else:
            limits1.append(np.nan)
            limits2.append(np.nan)

        for observation_num, observation in enumerate(self.observations):
            parameters_map[observation_num].append(len(names) - 1)

        names.append('a')
        print_names.append('a/R_*')
        initial.append(self.sma_over_rs)
        if fit_sma_over_rs:
            limits1.append(self.sma_over_rs * fit_sma_over_rs_limits[0])
            limits2.append(self.sma_over_rs * fit_sma_over_rs_limits[1])
        else:
            limits1.append(np.nan)
            limits2.append(np.nan)

        for observation_num, observation in enumerate(self.observations):
            parameters_map[observation_num].append(len(names) - 1)

        names.append('e')
        print_names.append('e')
        initial.append(self.eccentricity)
        limits1.append(np.nan)
        limits2.append(np.nan)

        for observation_num, observation in enumerate(self.observations):
            parameters_map[observation_num].append(len(names) - 1)

        names.append('i')
        print_names.append('i')
        initial.append(self.inclination)
        if fit_inclination:
            limits1.append(fit_inclination_limits[0])
            limits2.append(fit_inclination_limits[1])
        else:
            limits1.append(np.nan)
            limits2.append(np.nan)

        for observation_num, observation in enumerate(self.observations):
            parameters_map[observation_num].append(len(names) - 1)

        names.append('w')
        print_names.append('\omega')
        initial.append(self.periastron)
        limits1.append(np.nan)
        limits2.append(np.nan)

        for observation_num, observation in enumerate(self.observations):
            parameters_map[observation_num].append(len(names) - 1)

        # time parameters

        unique_epochs = []
        for observation in self.observations:
            if self.observations[observation]['epoch'] not in unique_epochs:
                unique_epochs.append(self.observations[observation]['epoch'])

        self.eclipse_mid_time += min(unique_epochs) * self.period

        if not fit_individual_times:

            names.append('T_0')
            print_names.append('T_0')
            initial.append(self.eclipse_mid_time)
            if fit_mid_time:
                limits1.append(self.eclipse_mid_time + fit_mid_time_limits[0])
                limits2.append(self.eclipse_mid_time + fit_mid_time_limits[1])
            else:
                limits1.append(np.nan)
                limits2.append(np.nan)

            for observation_num, observation in enumerate(self.observations):
                parameters_map[observation_num].append(len(names) - 1)

        else:

            if fit_period:
                raise PyLCInputError('Period and individual mid times cannot be fitted simultaneously.')

            for epoch in unique_epochs:

                epoch -= min(unique_epochs)

                names.append('T_{0}'.format(epoch))
                print_names.append('T_{' + str(epoch) + '}')
                initial.append(self.eclipse_mid_time + epoch * self.period)
                if fit_mid_time:
                    limits1.append(self.eclipse_mid_time + epoch * self.period + fit_mid_time_limits[0])
                    limits2.append(self.eclipse_mid_time + epoch * self.period + fit_mid_time_limits[1])
                else:
                    limits1.append(np.nan)
                    limits2.append(np.nan)

                for observation_num, observation in enumerate(self.observations):
                    if self.observations[observation]['epoch'] - min(unique_epochs) == epoch:
                        parameters_map[observation_num].append(len(names) - 1)

        model_time = np.array([])
        model_flux = np.array([])
        model_flux_unc = np.array([])
        for observation in self.observations:
            model_time = np.append(model_time, self.observations[observation]['time'])
            model_flux = np.append(model_flux, self.observations[observation]['flux'])
            model_flux_unc = np.append(model_flux_unc, self.observations[observation]['flux_unc'])

        parameters_map = np.array(parameters_map)

        def detrend_model(model_time, *model_variables):

            model = []

            for observation_num, observation in enumerate(self.observations):

                detrend_zero,  detrend_one, detrend_two, f, r, p, a, e, i, w, mt = \
                    np.array(model_variables)[parameters_map[observation_num]]

                deltatime = self.observations[observation]['dtime']

                model.append(detrend_zero * (1 + detrend_one * deltatime + detrend_two * deltatime * deltatime))

            return np.concatenate(model)

        def full_model(model_time, *model_variables):

            model = []

            for observation_num, observation in enumerate(self.observations):
                detrend_zero, detrend_one, detrend_two, f, r, p, a, e, i, w, mt = \
                    np.array(model_variables)[parameters_map[observation_num]]

                deltatime = self.observations[observation]['dtime']

                model1 = detrend_zero * (1 + detrend_one * deltatime + detrend_two * deltatime * deltatime)

                model.append(model1 * eclipse_integrated(f, r, p, a, e, i, w, mt,
                                                         time_array=self.observations[observation]['time'],
                                                         exp_time=self.observations[observation]['exp_time'],
                                                         max_sub_exp_time=max_sub_exp_time,
                                                         precision=precision
                                                         ))

            return np.concatenate(model)

        fitting = EmceeFitting(model_time, model_flux, model_flux_unc,
                               full_model, initial, limits1, limits2, walkers,
                               iterations, burn_in,
                               data_x_name='time', data_y_name='flux', data_x_print_name='Time',
                               data_y_print_name='Relative Flux',
                               parameters_names=names, parameters_print_names=print_names,
                               counter=counter)

        fitting.run_mcmc()
        results = fitting.results

        # TODO
        # save observations separately and calculate some statistics individually

        results['settings'] = {}
        results['settings']['max_sub_exp_time'] = max_sub_exp_time
        results['settings']['precision'] = precision
        results['settings']['detrending_order'] = detrending_order
        results['settings']['iterations'] = iterations
        results['settings']['walkers'] = walkers
        results['settings']['burn_in'] = burn_in
        results['settings']['fit_fp_over_fs'] = fit_fp_over_fs
        results['settings']['force_same_fp_over_fs'] = fit_individual_fp_over_fs
        results['settings']['fit_rp_over_rs'] = fit_rp_over_rs
        results['settings']['force_same_rp_over_rs'] = fit_individual_rp_over_rs
        results['settings']['fit_sma_over_rs'] = fit_sma_over_rs
        results['settings']['fit_inclination'] = fit_inclination
        results['settings']['fit_mid_time'] = fit_mid_time
        results['settings']['fit_period'] = fit_period
        results['settings']['fit_individual_times'] = fit_individual_times
        results['settings']['fit_rp_over_rs_limits'] = fit_rp_over_rs_limits
        results['settings']['fit_rp_over_rs_limits'] = fit_fp_over_fs_limits
        results['settings']['fit_sma_over_rs_limits'] = fit_sma_over_rs_limits
        results['settings']['fit_inclination_limits'] = fit_inclination_limits
        results['settings']['fit_mid_time_limits'] = fit_mid_time_limits
        results['settings']['fit_period_limits'] = fit_period_limits
        results['settings']['filter_map'] = self.filters

        results['detrended_input_series'] = {
            'time': results['input_series']['time'],
            'flux': results['input_series']['flux'] / detrend_model(model_time, *results['parameters_final']),
            'flux_unc': results['input_series']['flux_unc'] / detrend_model(model_time, *results['parameters_final'])}

        results['detrended_output_series'] = {
            'model': results['output_series']['model'] / detrend_model(model_time, *results['parameters_final']),
            'residuals': (results['output_series']['residuals']
                          / detrend_model(model_time, *results['parameters_final']))}

        results['detrended_statistics'] = {
            'res_mean': np.mean(results['detrended_output_series']['residuals']),
            'res_std': np.std(results['detrended_output_series']['residuals']),
            'res_rms': np.sqrt(np.mean(results['detrended_output_series']['residuals'] ** 2))
        }

        observations_id_series = np.array([])
        for observation in self.observations:
            observations_id_series = np.append(observations_id_series,
                                               np.ones(len(self.observations[observation]['time'])) * observation)

        for observation in self.observations:

            id_series = np.where(observations_id_series == observation)

            self.observations[observation]['model'] = results['output_series']['model'][id_series]
            self.observations[observation]['residuals'] = results['output_series']['residuals'][id_series]

            res_autocorr = np.correlate(self.observations[observation]['residuals'],
                                        self.observations[observation]['residuals'], mode='full')
            res_autocorr = res_autocorr[res_autocorr.size // 2:] / res_autocorr[res_autocorr.size // 2:][0]

            self.observations[observation]['res_autocorr'] = res_autocorr
            self.observations[observation]['res_max_autocorr'] = np.max(res_autocorr[1:])
            self.observations[observation]['res_mean'] = np.mean(self.observations[observation]['residuals'])
            self.observations[observation]['res_std'] = np.std(self.observations[observation]['residuals'])
            self.observations[observation]['res_rms'] = np.sqrt(np.mean(self.observations[observation]['residuals'] ** 2))
            self.observations[observation]['res_chi_sqr'] = np.sum(
                (self.observations[observation]['residuals'] ** 2) / (self.observations[observation]['flux_unc'] ** 2))
            self.observations[observation]['res_red_chi_sqr'] = (
                    self.observations[observation]['res_chi_sqr'] / (
                        len(self.observations[observation]['flux_unc']) - len(results['statistics']['corr_variables'])))

            self.observations[observation]['detrended_flux'] = results['detrended_input_series']['flux'][id_series]
            self.observations[observation]['detrended_flux_unc'] = results['detrended_input_series']['flux_unc'][id_series]
            self.observations[observation]['detrended_model'] = results['detrended_output_series']['model'][id_series]
            self.observations[observation]['detrended_residuals'] = results['detrended_output_series']['residuals'][id_series]

            res_autocorr = np.correlate(self.observations[observation]['detrended_residuals'],
                                        self.observations[observation]['detrended_residuals'], mode='full')
            res_autocorr = res_autocorr[res_autocorr.size // 2:] / res_autocorr[res_autocorr.size // 2:][0]

            self.observations[observation]['detrended_res_autocorr'] = res_autocorr
            self.observations[observation]['detrended_res_max_autocorr'] = np.max(res_autocorr[1:])
            self.observations[observation]['detrended_res_mean'] = np.mean(self.observations[observation]['detrended_residuals'])
            self.observations[observation]['detrended_res_std'] = np.std(self.observations[observation]['detrended_residuals'])
            self.observations[observation]['detrended_res_rms'] = np.sqrt(
                np.mean(self.observations[observation]['detrended_residuals'] ** 2))
            self.observations[observation]['detrended_res_chi_sqr'] = np.sum(
                (self.observations[observation]['detrended_residuals'] ** 2) / (self.observations[observation]['detrended_flux_unc'] ** 2))
            self.observations[observation]['detrended_res_red_chi_sqr'] = (
                    self.observations[observation]['detrended_res_chi_sqr'] / (
                    len(self.observations[observation]['detrended_flux_unc']) - len(results['statistics']['corr_variables'])))

            w = open(os.path.join(output_folder, 'diagnostics_dataset_{0}.txt'.format(observation + 1)), 'w')
            w.write('\n#Residuals:\n')
            w.write('#Mean: {0}\n'.format(self.observations[observation]['res_mean']))
            w.write('#STD: {0}\n'.format(self.observations[observation]['res_std']))
            w.write('#RMS: {0}\n'.format(self.observations[observation]['res_rms']))
            w.write('#Max auto-correlation: {0}\n'.format(self.observations[observation]['res_max_autocorr']))
            w.write('#Chi squared: {0}\n'.format(self.observations[observation]['res_chi_sqr']))
            w.write('#Reduced chi squared: {0}\n'.format(self.observations[observation]['res_red_chi_sqr']))

            w.write('\n\n#Detrended Residuals:\n')
            w.write('#Mean: {0}\n'.format(self.observations[observation]['detrended_res_mean']))
            w.write('#STD: {0}\n'.format(self.observations[observation]['detrended_res_std']))
            w.write('#RMS: {0}\n'.format(self.observations[observation]['detrended_res_rms']))
            w.write('#Max auto-correlation: {0}\n'.format(self.observations[observation]['detrended_res_max_autocorr']))
            w.write('#Chi squared: {0}\n'.format(self.observations[observation]['detrended_res_chi_sqr']))
            w.write('#Reduced chi squared: {0}\n'.format(self.observations[observation]['detrended_res_red_chi_sqr']))

            w.close()

        results['observations'] = self.observations

        save_dict(results, os.path.join(output_folder, 'results.pickle'))

        fitting.save_results(os.path.join(output_folder, 'results.txt'))

        fitting.plot_corner(os.path.join(output_folder, 'correlations.pdf'))

        fitting.plot_traces(os.path.join(output_folder, 'traces.pdf'))

        plot_transit_fitting_models(results, os.path.join(output_folder, 'lightcurves.pdf'))
