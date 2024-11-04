
__all__ = ['Planet', 'get_planet', 'get_all_planets']

import os
import numpy as np
import shutil
import time as sys_time
import itertools
import exoclock


from inspect import signature
from astropy.time import Time as astrotime

from ..analysis.optimisation import Fitting
from ..analysis.statistics import residual_statistics

from ..errors import *
from ..models.exoplanet_lc import planet_orbit, planet_star_projected_distance, planet_phase, \
    transit, transit_integrated, \
    transit_duration, transit_depth, eclipse, eclipse_integrated, eclipse_mid_time, eclipse_duration, eclipse_depth, \
    fp_over_fs, transit_t12, exotethys, convert_to_bjd_tdb, convert_to_jd_utc, convert_to_relflux, airmass
from ..processes.files import save_dict, copy_dict, open_dict
from ..plots.plots_fitting import plot_transit_fitting_models


class Observation:

    def __init__(self,
                 time, time_format, exp_time, time_stamp,
                 flux, flux_unc, flux_format,
                 filter_name, wlrange, stellar_model,
                 auxiliary_data, observatory_latitude, observatory_longitude,
                 detrending_series, detrending_order,
                 trend_function, trend_parameters,
                 ):

        if time_stamp not in ['start', 'mid', 'end']:
            raise PyLCInputError(
                'Not acceptable time stamp {0}. Please choose between "mid", "start", "end".'.format(time_stamp))

        if trend_function:
            if not trend_parameters:
                raise PyLCInputError('You need to provide the detrending parameters information '
                                     '(initial guess, limits, names, names to print)')
            else:
                detrending_series = None
                detrending_order = None

        elif detrending_series:
            detrending_series_ok = True

            available_detrending_series = ['time'] + list(auxiliary_data.keys())
            if None not in [observatory_latitude, observatory_longitude]:
                available_detrending_series.append('airmass')

            if isinstance(detrending_series, str):
                if detrending_series not in available_detrending_series:
                    detrending_series_ok = False
                else:
                    detrending_series = [detrending_series]

            elif isinstance(detrending_series, list):
                for i in detrending_series:
                    if i not in available_detrending_series:
                        detrending_series_ok = False
            else:
                detrending_series_ok = False

            if not detrending_series_ok:
                raise PyLCInputError('Not acceptable detrending_series: {0}. '
                                     'Please provide a string or a list of strings. '
                                     'Available detrending_series: {1}'.format(
                    detrending_series, ','.join(available_detrending_series)))

            try:
                detrending_order = int(detrending_order)
                if detrending_order <= 0:
                    raise PyLCInputError('You need to detrending_order must be a positive integer.')
            except:
                raise PyLCInputError('You need to detrending_order must be a positive integer.')

        else:
            raise PyLCInputError('You need to indicate a detrending_series or a trend_function')

        self.dict = {
            'original_data': {
                'time': np.ones_like(time) * time,
                'time_format': time_format,
                'exp_time': exp_time,
                'time_stamp': time_stamp,
                'flux': np.ones_like(flux) * flux,
                'flux_unc': np.ones_like(flux_unc) * flux_unc,
                'flux_format': flux_format,
                'filter_name': filter_name,
                'wlrange': wlrange,
                'stellar_model': stellar_model,
                'observatory_longitude': observatory_longitude,
                'observatory_latitude': observatory_latitude,
                'auxiliary_data': auxiliary_data,
                'trend_function': trend_function,
                'trend_parameters': trend_parameters,
                'detrending_series': detrending_series,
                'detrending_order': detrending_order,
            }
        }

    def reset(self, planet):

        self.dict['input_series'] = {}
        self.dict['output_series'] = {}
        self.dict['detrended_series'] = {}
        self.dict['statistics'] = {}

        self.dict['data_conversion_info'] = {
            'notes': [],
            'ra': planet.ra,
            'dec': planet.dec,
        }
        self.dict['model_info'] = {
            'target': planet.name,
            'date': astrotime(
                convert_to_jd_utc(planet.ra, planet.dec, min(self.dict['original_data']['time']),
                                  self.dict['original_data']['time_format']),
                format='jd'
            ).datetime.isoformat().split('T')[0],
            'ldc_method': planet.method,
            'max_sub_exp_time': planet.max_sub_exp_time,
            'precision': planet.precision,
            'stellar_temperature': planet.stellar_temperature,
            'stellar_logg': planet.stellar_logg,
            'stellar_metallicity': planet.stellar_metallicity,
            'limb_darkening_coefficients': planet.exotethys(self.dict['original_data']['filter_name'],
                                                            self.dict['original_data']['wlrange'],
                                                            self.dict['original_data']['stellar_model']),
            'emissivity': planet.emissivity,
            'albedo': planet.albedo,
            'sma_over_rs': planet.sma_over_rs,
            'fp_over_fs': planet.fp_over_fs(self.dict['original_data']['filter_name'],
                                            self.dict['original_data']['wlrange'])
        }

        self.precision = planet.precision
        self.max_sub_exp_time = planet.max_sub_exp_time
        self.ldc_method = planet.method

        # convert time

        time_stamp_correction = {
            'start': 0.5,
            'mid': 0,
            'end': -0.5
        }[self.dict['original_data']['time_stamp']]
        time_stamp_correction = time_stamp_correction * self.dict['original_data']['exp_time'] / (60.0 * 60.0 * 24.0)

        time = np.array(self.dict['original_data']['time']) + time_stamp_correction
        self.dict['data_conversion_info']['time_stamp_correction'] = time_stamp_correction
        self.dict['data_conversion_info']['notes'].append('Time converted to mid-exposure time.')

        bjd_tdb_correction = planet.convert_to_bjd_tdb(time, self.dict['original_data']['time_format']) - time
        time = np.array(time) + bjd_tdb_correction
        self.dict['data_conversion_info']['bjd_tdb_correction'] = time_stamp_correction
        self.dict['data_conversion_info']['notes'].append('Time converted to BJD_TDB.')

        self.dict['input_series']['time'] = time

        # convert flux

        flux, flux_unc = planet.convert_to_relflux(self.dict['original_data']['flux'],
                                                   self.dict['original_data']['flux_unc'],
                                                   self.dict['original_data']['flux_format'])
        self.dict['input_series']['flux'], self.dict['input_series']['flux_unc'] = \
            flux, flux_unc
        self.dict['data_conversion_info']['notes'].append('Flux and flux_unc converted to relative flux.')

        #  auxiliary_data

        for auxiliary_timeseries in self.dict['original_data']['auxiliary_data']:
            self.dict['input_series'][auxiliary_timeseries] = \
                self.dict['original_data']['auxiliary_data'][auxiliary_timeseries]

        if (self.dict['original_data']['observatory_latitude'] is not None
                and self.dict['original_data']['observatory_longitude'] is not None):
            self.dict['input_series']['airmass'] = planet.airmass(time,
                                                                  self.dict['original_data']['observatory_latitude'],
                                                                  self.dict['original_data']['observatory_longitude']
                                                                  )
        # sort by time

        idx_sorted = sorted(np.arange(len(time)), key=lambda x: time[x])
        for input_series in self.dict['input_series']:
            self.dict['input_series'][input_series] = np.array(self.dict['input_series'][input_series])[idx_sorted]

        # ids
        if self.dict['original_data']['wlrange'] is None:
            self.filter_id = self.dict['original_data']['filter_name']
        else:
            self.filter_id = '{0}:{1}-{2}'.format(self.dict['original_data']['filter_name'],
                                                  self.dict['original_data']['wlrange'][0],
                                                  self.dict['original_data']['wlrange'][1])

        self.dict['model_info']['filter_id'] = self.filter_id

        transit_phase = (np.mean(time) - planet.mid_time) / planet.period
        transit_epoch = int(round(transit_phase, 0))

        eclipse_phase = (np.mean(time) - planet.eclipse_mid_time()) / planet.period
        eclipse_epoch = int(round(eclipse_phase, 0))

        if abs(transit_phase - transit_epoch) < abs(eclipse_phase - eclipse_epoch):
            self.observation_type = 'transit'
            self.epoch = transit_epoch
        else:
            self.observation_type = 'eclipse'
            self.epoch = eclipse_epoch

        self.dict['model_info']['observation_type'] = self.observation_type
        self.dict['model_info']['epoch'] = self.epoch

        self.time_factor = int(self.dict['original_data']['exp_time'] / self.max_sub_exp_time) + 1
        self.exp_time_h = self.dict['original_data']['exp_time'] / (60.0 * 60.0 * 24.0)

        self.series_length = len(time)
        self.time_hr = (time[:, None] + np.arange(
            -self.exp_time_h / 2 + self.exp_time_h / self.time_factor / 2,
            self.exp_time_h / 2,
            self.exp_time_h / self.time_factor
        )).flatten()

        self.parameters_map = np.array([], dtype=int)
        self.non_outliers_map = np.array(np.arange(self.series_length), dtype=int)

        # setup de-trending
        if not self.dict['original_data']['trend_function']:

            detrending_names = [ff for ff in self.dict['original_data']['detrending_series']]

            detrending_series = [self.dict['input_series'][ff] - np.min(self.dict['input_series'][ff]) for ff in detrending_names]

            detrending_names.append('???')
            detrending_series.append(np.ones_like(detrending_series[0]))

            detrending_names = np.array(detrending_names)
            detrending_series = np.array(detrending_series)

            combinations = list(itertools.combinations_with_replacement(np.arange(len(detrending_series)),
                                                                        self.dict['original_data']['detrending_order']))
            detrending_names = detrending_names[(combinations,)]
            detrending_series = np.prod(detrending_series[(combinations,)], 1)
            detrending_names = detrending_names[np.where(np.sum((detrending_series - 1) ** 2, 1) != 0)]
            detrending_series = detrending_series[np.where(np.sum((detrending_series - 1) ** 2, 1) != 0)]

            self.detrending = detrending_series
            self.number_of_detrending_parameters = 1 + len(detrending_series)
            self.names = ['n'] + ['_'.join(ff).replace('_???', '') for ff in detrending_names]
            self.print_names = ['n'] + ['-'.join(ff).replace('-???', '') for ff in detrending_names]
            self.limits1 = [0] + [-2 for ff in range(len(detrending_series))]
            self.limits2 = [2] + [2 for ff in range(len(detrending_series))]
            self.initial = [1] + [0 for ff in range(len(detrending_series))]
            self.logspace = [True] + [False for ff in range(len(detrending_series))]

        else:
            self.detrending = None

            self.number_of_detrending_parameters = len(str(signature(
                self.dict['original_data']['trend_function']))[1:-1].split(','))

            self.initial = [1] + [ff[0] for ff in self.dict['original_data']['trend_parameters']]
            self.limits1 = [0] + [ff[1] for ff in self.dict['original_data']['trend_parameters']]
            self.limits2 = [2] + [ff[2] for ff in self.dict['original_data']['trend_parameters']]
            self.names = ['n'] + [ff[3] for ff in self.dict['original_data']['trend_parameters']]
            self.print_names = ['n'] + [ff[4] for ff in self.dict['original_data']['trend_parameters']]
            self.logspace = [True] + [False for ff in self.dict['original_data']['trend_parameters']]

        # adjust normalisation factor
        df = max(np.median(flux) - np.min(flux), np.max(flux) - np.median(flux))
        self.initial[0] = np.median(flux)
        self.limits1[0] = np.min(flux) - 2 * df
        self.limits2[0] = np.max(flux) + 2 * df

    def add_to_parameters_map(self, idx):
        self.parameters_map = np.append(self.parameters_map, int(idx))

    def clear_outliers(self):

        if self.detrending is not None:
            self.detrending = self.detrending[:, self.non_outliers_map]

        self.dict['input_series'] = {ff:self.dict['input_series'][ff][self.non_outliers_map]
                                     for ff in self.dict['input_series']}

        self.series_length = len(self.dict['input_series']['time'])
        self.time_hr = (self.dict['input_series']['time'][:, None] + np.arange(
            -self.exp_time_h / 2 + self.exp_time_h / self.time_factor / 2, self.exp_time_h / 2, self.exp_time_h / self.time_factor
        )).flatten()

    def _signal_model(self, parameters):
        if self.dict['model_info']['observation_type'] == 'transit':
            ldc1, ldc2, ldc3, ldc4, r, p, a, e, i, w, mt = parameters
            return np.mean(np.reshape(
                transit([ldc1, ldc2, ldc3, ldc4], r, p, a, e, i, w, mt, self.time_hr,
                        method=self.ldc_method, precision=self.precision),
                (self.series_length, self.time_factor)), 1)

        elif self.dict['model_info']['observation_type'] == 'eclipse':
            d, r, p, a, e, i, w, mt = parameters
            return np.mean(np.reshape(
                eclipse(d, r, p, a, e, i, w, mt, self.time_hr, precision=self.precision),
                (self.series_length, self.time_factor)), 1)

    def _trend_model(self, parameters):
        if not self.dict['original_data']['trend_function']:
            return 1 + np.sum(self.detrending * np.array(parameters)[:, None], 0)
        else:
            return self.dict['original_data']['trend_function'](self.dict['input_series'], *parameters)

    def splitted_model(self, indices, parameters):

        parameters = parameters[self.parameters_map]
        trend_model = self._trend_model(parameters[1:self.number_of_detrending_parameters])
        signal_model = self._signal_model(parameters[self.number_of_detrending_parameters:])

        if indices is None:
            return (parameters[0] * trend_model), signal_model
        else:
            return (parameters[0] * trend_model)[indices], signal_model[np.int_(indices)]

    def full_model(self, indices, parameters):

        parameters = parameters[self.parameters_map]
        trend_model = self._trend_model(parameters[1:self.number_of_detrending_parameters])
        signal_model = self._signal_model(parameters[self.number_of_detrending_parameters:])

        if indices is None:
            return parameters[0] * trend_model * signal_model
        else:
            return (parameters[0] * trend_model * signal_model)[np.int_(indices)]


class Planet:

    def __init__(self, name, ra, dec, stellar_logg, stellar_temperature, stellar_metallicity,
                 rp_over_rs, period, sma_over_rs, eccentricity, inclination,
                 periastron, mid_time,
                 mid_time_format='BJD_TDB',
                 ldc_method='claret', max_sub_exp_time=10, precision=2,
                 asc_node=0,
                 albedo=0.15, emissivity=1.0):

        self.name = name
        self.ra = ra
        self.dec = dec
        if self.dec > 90 or self.dec < -90:
            raise PyLCInputError('Declination must be within -90, 90 degrees')
        self.stellar_logg = stellar_logg
        self.stellar_temperature = stellar_temperature
        self.stellar_metallicity = stellar_metallicity

        self.rp_over_rs = rp_over_rs
        self.period = period
        self.sma_over_rs = sma_over_rs
        self.eccentricity = eccentricity
        self.inclination = inclination
        self.periastron = periastron
        self.mid_time = convert_to_bjd_tdb(ra, dec, mid_time, mid_time_format)

        self.asc_node = asc_node
        self.albedo = albedo
        self.emissivity = emissivity

        self.method = ldc_method
        self.max_sub_exp_time = max_sub_exp_time
        self.precision = precision

        self.ldcs_memory = {}
        self.observations = []

    # conversions

    def convert_to_bjd_tdb(self, time_array, time_format):
        return convert_to_bjd_tdb(self.ra, self.dec, time_array, time_format)

    def convert_to_jd_utc(self, time_array, time_format):
        return convert_to_jd_utc(self.ra, self.dec, time_array, time_format)

    def convert_to_relflux(self, flux_array, flux_unc_array, flux_format):
        return convert_to_relflux(flux_array, flux_unc_array, flux_format)

    def airmass(self, time_array, observatory_latitude, observatory_longitude):
        return airmass(observatory_latitude, observatory_longitude, self.ra, self.dec, time_array)

    # filter-dependent values

    def exotethys(self, filter_name, wlrange=None, stellar_model='Phoenix_2018'):
        options = '{0}_{1}_{2}_{3}'.format(filter_name, wlrange, self.method, stellar_model)
        if options not in self.ldcs_memory:
            self.ldcs_memory[options] = exotethys(self.stellar_logg, self.stellar_temperature, self.stellar_metallicity,
                                                  filter_name, wlrange, self.method, stellar_model)
        return self.ldcs_memory[options]

    def add_custom_limb_darkening_coefficients(self, limb_darkening_coefficients, filter_name, wlrange, stellar_model):
        options = '{0}_{1}_{2}_{3}'.format(filter_name, wlrange, self.method, stellar_model)
        self.ldcs_memory[options] = limb_darkening_coefficients

    def fp_over_fs(self, filter_name, wlrange=None):
        return fp_over_fs(self.rp_over_rs, self.sma_over_rs, self.albedo, self.emissivity, self.stellar_temperature,
                          filter_name, wlrange)

    def planet_orbit(self, time):
        return planet_orbit(self.period, self.sma_over_rs, self.eccentricity, self.inclination, self.periastron,
                            self.mid_time, time, ww=self.asc_node)

    def planet_star_projected_distance(self, time):
        return planet_star_projected_distance(self.period, self.sma_over_rs, self.eccentricity, self.inclination,
                                              self.periastron, self.mid_time, time)

    def planet_phase(self, time):
        return planet_phase(self.period, self.mid_time, time)

    def transit(self, time, filter_name, wlrange=None, stellar_model='Phoenix_2018'):
        return transit(self.exotethys(filter_name, wlrange, stellar_model),
                       self.rp_over_rs, self.period, self.sma_over_rs, self.eccentricity,
                       self.inclination, self.periastron, self.mid_time, time,
                       self.method, self.precision)

    def transit_integrated(self, time, exp_time, filter_name, wlrange=None, stellar_model='Phoenix_2018'):
        return transit_integrated(self.exotethys(filter_name, wlrange, stellar_model),
                                  self.rp_over_rs, self.period, self.sma_over_rs, self.eccentricity,
                                  self.inclination, self.periastron, self.mid_time, time, exp_time,
                                  self.max_sub_exp_time, self.method, self.precision)

    def transit_duration(self):
        return transit_duration(self.rp_over_rs, self.period, self.sma_over_rs, self.eccentricity, self.inclination,
                                self.periastron)

    def transit_t12(self):
        return transit_t12(self.rp_over_rs, self.period, self.sma_over_rs, self.eccentricity, self.inclination,
                           self.periastron)

    def transit_depth(self, filter_name, wlrange=None, stellar_model='Phoenix_2018'):
        return transit_depth(self.exotethys(filter_name, wlrange, stellar_model),
                             self.rp_over_rs, self.period, self.sma_over_rs, self.eccentricity,
                             self.inclination, self.periastron,
                             self.method, self.precision)

    def eclipse(self, time, filter_name, wlrange=None):
        return eclipse(self.fp_over_fs(filter_name, wlrange), self.rp_over_rs,
                       self.period, self.sma_over_rs, self.eccentricity,
                       self.inclination, self.periastron, self.eclipse_mid_time(), time, self.precision)

    def eclipse_integrated(self, time, exp_time, filter_name, wlrange=None):
        return eclipse_integrated(self.fp_over_fs(filter_name, wlrange=wlrange), self.rp_over_rs,
                                  self.period, self.sma_over_rs, self.eccentricity,
                                  self.inclination, self.periastron, self.eclipse_mid_time(), time, exp_time,
                                  self.max_sub_exp_time, self.precision)

    def eclipse_duration(self):
        return eclipse_duration(self.rp_over_rs, self.period, self.sma_over_rs, self.eccentricity,
                                self.inclination, self.periastron)

    def eclipse_depth(self, filter_name, wlrange=None):
        return eclipse_depth(self.fp_over_fs(filter_name, wlrange), self.rp_over_rs, self.period,
                             self.sma_over_rs, self.eccentricity, self.inclination, self.periastron,
                             self.precision)

    def eclipse_mid_time(self):
        return eclipse_mid_time(self.period, self.sma_over_rs, self.eccentricity, self.inclination, self.periastron,
                                self.mid_time)

    # data fitting

    def add_observation(self, time, time_format, exp_time, time_stamp, flux, flux_unc, flux_format,
                        filter_name, wlrange=None, stellar_model='Phoenix_2018',
                        auxiliary_data={}, observatory_latitude=None, observatory_longitude=None,
                        detrending_series='time', detrending_order=2,
                        trend_function=None, trend_parameters=None,
                        **kwargs
                        ):

        self.observations.append(
            Observation(time, time_format, exp_time, time_stamp,
                        flux, flux_unc, flux_format,
                        filter_name, wlrange, stellar_model,
                        auxiliary_data, observatory_latitude, observatory_longitude,
                        detrending_series, detrending_order,
                        trend_function, trend_parameters))

    def add_observation_from_dict(self, dictionary):

        if isinstance(dictionary, str):
            dictionary = open_dict(dictionary)
        elif not isinstance(dictionary, dict):
            raise PyLCInputError('Dictionary should be a dict object or a path to a pickle file.')

        self.add_observation(**dictionary)

    def clear_observations(self):

        self.observations = []

    def transit_fitting(self, output_folder=None,

                        iterations=None, walkers=None, burn_in=None,

                        fit_ldc1=False, fit_ldc2=False, fit_ldc3=False, fit_ldc4=False,

                        fit_individual_rp_over_rs=True,
                        fit_sma_over_rs=False, fit_inclination=False,
                        fit_mid_time=True, fit_period=False,
                        fit_individual_times=True,

                        fit_ldc_limits=[-3.0, 3.0],
                        fit_rp_over_rs_limits=[0.1, 10],
                        fit_sma_over_rs_limits=[0.001, 1000],
                        fit_inclination_limits=[10.0, 90.0],
                        fit_mid_time_limits=[-0.2, 0.2],
                        fit_period_limits=[0.001, 1000],

                        counter='Transit fitting',
                        optimise_initial_parameters_trials=3,
                        scale_uncertainties=False,
                        filter_outliers=False,
                        optimiser='emcee',
                        strech_prior=0.2,
                        walkers_spread=0.01,
                        return_traces=True,
                        verbose=False,
                        ):

        for observation in self.observations:
            observation.reset(self)
            if observation.dict['model_info']['observation_type'] == 'eclipse':
                raise PyLCInputError('You need to add only transit observation to proceed.')

        if len(self.observations) == 0:
            raise PyLCInputError('You need to add at least one transit observation to proceed.')

        for observation_num, observation in enumerate(self.observations):
            observation.obs_id = 'obs{0}'.format(observation_num)
            observation.dict['model_info']['obs_id'] = observation.obs_id

        fit_rp_over_rs_limits = [self.rp_over_rs * fit_rp_over_rs_limits[0], self.rp_over_rs * fit_rp_over_rs_limits[1]]
        fit_sma_over_rs_limits = [self.sma_over_rs * fit_sma_over_rs_limits[0], self.sma_over_rs * fit_sma_over_rs_limits[1]]
        fit_period_limits = [self.period * fit_period_limits[0], self.period * fit_period_limits[1]]

        # separate filters and epochs

        unique_epochs = []
        for observation in self.observations:
            if observation.epoch not in unique_epochs:
                unique_epochs.append(observation.epoch)

        if len(unique_epochs) == 1:
            fit_individual_times = False
            if fit_period:
                raise PyLCInputError('Period cannot be fitted only on one epoch.')

        if fit_individual_times and fit_period:
            raise PyLCInputError('Period and individual mid times cannot be fitted simultaneously.')

        unique_filters = {}
        for observation in self.observations:
            if observation.filter_id not in unique_filters:
                unique_filters[observation.filter_id] = \
                    observation.dict['model_info']['limb_darkening_coefficients']

        if len(unique_filters) == 1:
            fit_individual_rp_over_rs = False

        # add parameters

        names = []
        print_names = []
        limits1 = []
        limits2 = []
        initial = []
        logspace = []

        # de-trending parameters

        for observation in self.observations:

            for coefficient in range(observation.number_of_detrending_parameters):

                if len(self.observations) > 1:
                    names.append('{0}::{1}'.format(observation.names[coefficient], observation.obs_id))
                    print_names.append('{0}::{1}'.format(observation.print_names[coefficient], observation.obs_id))
                else:
                    names.append('{0}'.format(observation.names[coefficient]))
                    print_names.append('{0}'.format(observation.print_names[coefficient]))

                initial.append(observation.initial[coefficient])
                limits1.append(observation.limits1[coefficient])
                limits2.append(observation.limits2[coefficient])
                logspace.append(observation.logspace[coefficient])

                observation.add_to_parameters_map(len(names) - 1)

        # limb-darkening

        for phot_filter in unique_filters:

            ldcs = unique_filters[phot_filter]
            fit_ldcs = [fit_ldc1, fit_ldc2, fit_ldc3, fit_ldc4]

            for ldc in range(4):

                if len(unique_filters) > 1:
                    names.append('a_{0}::{1}'.format(ldc + 1, phot_filter))
                    print_names.append(r'$a_{{{0}}}::{1}$'.format(ldc +1, phot_filter))
                else:
                    names.append('a_{0}'.format(ldc + 1))
                    print_names.append(r'$a_{{{0}}}$'.format(ldc +1))

                initial.append(ldcs[ldc])
                if fit_ldcs[ldc]:
                    limits1.append(ldcs[ldc] + fit_ldc_limits[0])
                    limits2.append(ldcs[ldc] + fit_ldc_limits[1])
                else:
                    limits1.append(np.nan)
                    limits2.append(np.nan)
                logspace.append(False)

                for observation in self.observations:
                    if observation.filter_id == phot_filter:
                        observation.add_to_parameters_map(len(names) - 1)

        # rp_over_rs

        if fit_individual_rp_over_rs:
            for phot_filter in unique_filters:
                names.append('rp_over_rs::{0}'.format(phot_filter))
                print_names.append(r'$R_\mathrm{{p}}/R_*$::{0}'.format(phot_filter))
                initial.append(self.rp_over_rs)
                limits1.append(fit_rp_over_rs_limits[0])
                limits2.append(fit_rp_over_rs_limits[1])
                logspace.append(False)
                for observation in self.observations:
                    if observation.filter_id == phot_filter:
                        observation.add_to_parameters_map(len(names) - 1)

            global_parameters_names = []
            global_parameters_print_names = []
            global_parameters_initial = []
            global_parameters_fit = []
            global_parameters_limits = []
            global_parameters_logspace = []

        else:
            global_parameters_names = ['rp_over_rs']
            global_parameters_print_names = [r'$R_\mathrm{p}/R_*$']
            global_parameters_initial = [self.rp_over_rs]
            global_parameters_fit = [True]
            global_parameters_limits = [fit_rp_over_rs_limits]
            global_parameters_logspace = [False]

        # orbital parameters

        global_parameters_names += ['period', 'sma_over_rs', 'eccentricity', 'inclination', 'periastron']
        global_parameters_print_names += [r'$P$', r'$a$', r'$e$', r'$i$', r'$\omega$']
        global_parameters_initial += [self.period, self.sma_over_rs, self.eccentricity, self.inclination, self.periastron]
        global_parameters_fit += [fit_period, fit_sma_over_rs, False, fit_inclination, False]
        global_parameters_limits += [fit_period_limits, fit_sma_over_rs_limits, False, fit_inclination_limits, False]
        global_parameters_logspace += [True, True, False, True, False]

        if not fit_mid_time or not fit_individual_times or len(unique_epochs) == 1:

            test_epochs = np.array([])
            test_epochs_weights = np.array([])

            for observation in self.observations:
                test_epochs = np.append(test_epochs,
                                        np.ones_like(observation.dict['input_series']['flux_unc']) * observation.epoch)
                norm_errors = observation.dict['input_series']['flux_unc'] / observation.dict['input_series']['flux']
                test_epochs_weights = np.append(test_epochs_weights, 1 / (norm_errors * norm_errors))

            new_epoch = np.round(np.sum(test_epochs * test_epochs_weights) / np.sum(test_epochs_weights), 0)
            new_mid_time = self.mid_time + new_epoch * self.period

            global_parameters_names += ['mid_time']
            global_parameters_print_names += [r'$T_\mathrm{mid}$']
            global_parameters_initial += [new_mid_time]
            global_parameters_fit += [fit_mid_time]
            global_parameters_limits += [[new_mid_time + fit_mid_time_limits[0], new_mid_time + fit_mid_time_limits[1]]]
            global_parameters_logspace += [False]

        for global_parameter in range(len(global_parameters_names)):
            names.append(global_parameters_names[global_parameter])
            print_names.append(global_parameters_print_names[global_parameter])
            initial.append(global_parameters_initial[global_parameter])
            if global_parameters_fit[global_parameter]:
                limits1.append(global_parameters_limits[global_parameter][0])
                limits2.append(global_parameters_limits[global_parameter][1])
            else:
                limits1.append(np.nan)
                limits2.append(np.nan)
            logspace.append(global_parameters_logspace[global_parameter])

            for observation in self.observations:
                observation.add_to_parameters_map(len(names) - 1)

        # individual time parameters
        if fit_mid_time and fit_individual_times and len(unique_epochs) > 1:

            for epoch in unique_epochs:

                names.append('mid_time::{0}'.format(epoch))
                print_names.append(r'$T_\mathrm{{mid}}$::{0}'.format(epoch))
                initial.append(self.mid_time + epoch * self.period)
                limits1.append(self.mid_time + epoch * self.period + fit_mid_time_limits[0])
                limits2.append(self.mid_time + epoch * self.period + fit_mid_time_limits[1])
                logspace.append(False)

                for observation in self.observations:
                    if observation.epoch == epoch:
                        observation.add_to_parameters_map(len(names) - 1)

        initial = np.array(initial)
        limits1 = np.array(limits1)
        limits2 = np.array(limits2)
        logspace = np.array(logspace)

        # single observation tests

        for observation in self.observations:

            def single_observation_full_model(indices, *model_variables):
                return observation.full_model(indices, np.array(model_variables))

            fitting = Fitting(observation.non_outliers_map,
                              observation.dict['input_series']['flux'],
                              observation.dict['input_series']['flux_unc'],
                              single_observation_full_model, initial, limits1, limits2,
                              logspace=logspace,
                              data_x_name='time', data_y_name='flux',
                              data_x_print_name='r$t_{BJD_{TDB}}$', data_y_print_name='Relative Flux',
                              parameters_names=names, parameters_print_names=print_names,
                              walkers=walkers, iterations=iterations, burn_in=burn_in,
                              counter=counter,
                              optimise_initial_parameters=True,
                              optimise_initial_parameters_trials=optimise_initial_parameters_trials,
                              scale_uncertainties=scale_uncertainties,
                              filter_outliers=filter_outliers,
                              optimiser='curve_fit',
                              strech_prior=strech_prior,
                              walkers_spread=walkers_spread,
                              )

            fitting._prefit(verbose=verbose)

            initial[observation.parameters_map] = fitting.results['prefit']['initials'][observation.parameters_map]
            observation.dict['input_series']['flux_unc'] *= fitting.results['prefit']['scale_factor']
            observation.dict['data_conversion_info']['scale_factor'] = fitting.results['prefit']['scale_factor']
            observation.dict['data_conversion_info']['notes'].append('Flux_unc multiplied by {0}'.format(fitting.results['prefit']['scale_factor']))
            observation.dict['data_conversion_info']['outliers'] = len(np.where(fitting.results['prefit']['outliers_map'])[0])
            observation.dict['data_conversion_info']['notes'].append('{0} outliers removed'.format(len(np.where(fitting.results['prefit']['outliers_map'])[0])))
            observation.dict['data_conversion_info']['non_outliers_map'] = np.where(~fitting.results['prefit']['outliers_map'])[0]
            observation.non_outliers_map = np.where(~fitting.results['prefit']['outliers_map'])[0]

            print()
            print('Observation: ', observation.obs_id)
            print('Filter: ', observation.filter_id)
            print('Epoch: ', observation.epoch)
            print('Data-points excluded: ', len(np.where(fitting.results['prefit']['outliers_map'])[0]))
            print('Scaling uncertainties by: ', fitting.results['prefit']['scale_factor'])

            observation.clear_outliers()

        model_time = np.array([])
        model_flux = np.array([])
        model_flux_unc = np.array([])
        for observation in self.observations:
            model_time = np.append(model_time, observation.dict['input_series']['time'])
            model_flux = np.append(model_flux, observation.dict['input_series']['flux'])
            model_flux_unc = np.append(model_flux_unc, observation.dict['input_series']['flux_unc'])

        def detrend_model(model_time, *model_variables):
            model = np.array([])
            for observation in self.observations:
                model = np.append(model, observation.splitted_model(None, np.array(model_variables))[0])
            return model

        def full_model(model_time, *model_variables):
            model = np.array([])
            for observation in self.observations:
                model = np.append(model, observation.full_model(None, np.array(model_variables)))
            return model

        fitting = Fitting(model_time, model_flux, model_flux_unc,
                          full_model, initial, limits1, limits2,
                          logspace=logspace,
                          data_x_name='time', data_y_name='flux',
                          data_x_print_name='r$t_{BJD_{TDB}}$', data_y_print_name='Relative Flux',
                          parameters_names=names, parameters_print_names=print_names,
                          walkers=walkers, iterations=iterations, burn_in=burn_in,
                          counter=counter,
                          optimise_initial_parameters=True,
                          scale_uncertainties=False,
                          filter_outliers=False,
                          optimiser=optimiser,
                          strech_prior=strech_prior,
                          walkers_spread=walkers_spread,
                          )

        print('\nOptimising initial parameters...')

        fitting.run(verbose=verbose)

        fitting.results['settings'] = {
            'output_folder': output_folder,
            'iterations': fitting.iterations,
            'walkers': fitting.walkers,
            'burn_in': fitting.burn_in,
            'strech_prior': strech_prior,
            'walkers_spread': walkers_spread,
            'optimise_initial_parameters': True,
            'scale_uncertainties': scale_uncertainties,
            'filter_outliers': filter_outliers,
            'optimiser': optimiser,
            'data_x_name': 'time',
            'data_x_print_name': r't_${BJD_{TDB}}$',
            'data_y_name': 'flux',
            'data_y_print_name': 'Relative Flux',
            'fit_ldc1': fit_ldc1,
            'fit_ldc2': fit_ldc2,
            'fit_ldc3': fit_ldc3,
            'fit_ldc4': fit_ldc4,
            'fit_individual_rp_over_rs': fit_individual_rp_over_rs,
            'fit_sma_over_rs': fit_sma_over_rs,
            'fit_inclination': fit_inclination,
            'fit_mid_time': fit_mid_time,
            'fit_period': fit_period,
            'fit_individual_times': fit_individual_times,
            'fit_ldc_limits': fit_ldc_limits,
            'fit_rp_over_rs_limits': fit_rp_over_rs_limits,
            'fit_sma_over_rs_limits': fit_sma_over_rs_limits,
            'fit_inclination_limits': fit_inclination_limits,
            'fit_mid_time_limits': fit_mid_time_limits,
            'fit_period_limits': fit_period_limits,
        }

        del fitting.results['prefit']
        del fitting.results['original_data']

        trend = detrend_model(model_time, *fitting.results['parameters_final'])

        fitting.results['output_series']['trend'] = trend

        fitting.results['detrended_series'] = {
            'time': model_time,
            'flux': fitting.results['input_series']['flux'] / trend,
            'flux_unc': fitting.results['input_series']['flux_unc'] / trend,
            'model': fitting.results['output_series']['model'] / trend,
            'residuals': fitting.results['output_series']['residuals'] / trend
        }

        fitting.results['detrended_statistics'] = residual_statistics(fitting.results['detrended_series']['time'],
                                                                      fitting.results['detrended_series']['flux'],
                                                                      fitting.results['detrended_series']['flux_unc'],
                                                                      fitting.results['detrended_series']['model'],
                                                                      len(fitting.fitted_parameters))

        for observation in self.observations:

            observation.dict['parameters'] = {}
            number_of_free_parameters = 0
            for parameter_index in observation.parameters_map:
                parameter_name = names[parameter_index]
                parameter_sub_name = parameter_name.split('::')[0]
                parameter_data = copy_dict(fitting.results['parameters'][parameter_name])
                parameter_data['trace'] = None
                parameter_data['name'] = parameter_sub_name
                parameter_data['print_name'] = parameter_data['print_name'].split(':')[0]
                observation.dict['parameters'][parameter_sub_name] = parameter_data
                if fitting.results['parameters'][parameter_name]['initial'] is not None:
                    number_of_free_parameters += 0

            observation.dict['output_series']['model'] = observation.full_model(None, fitting.results['parameters_final'])
            observation.dict['output_series']['trend'] = observation.splitted_model(None, fitting.results['parameters_final'])[0]
            observation.dict['output_series']['residuals'] = observation.dict['input_series']['flux'] - observation.dict['output_series']['model']
            observation.dict['detrended_series']['time'] = observation.dict['input_series']['time']
            observation.dict['detrended_series']['flux'] = observation.dict['input_series']['flux'] / observation.dict['output_series']['trend']
            observation.dict['detrended_series']['flux_unc'] = observation.dict['input_series']['flux_unc'] / observation.dict['output_series']['trend']
            observation.dict['detrended_series']['model'] = observation.dict['output_series']['model'] / observation.dict['output_series']['trend']
            observation.dict['detrended_series']['residuals'] = observation.dict['detrended_series']['flux'] - observation.dict['detrended_series']['model']

            observation.dict['statistics'] = residual_statistics(observation.dict['input_series']['time'],
                                                                 observation.dict['input_series']['flux'],
                                                                 observation.dict['input_series']['flux_unc'],
                                                                 observation.dict['output_series']['model'],
                                                                 number_of_free_parameters)

            observation.dict['detrended_statistics'] = residual_statistics(observation.dict['detrended_series']['time'],
                                                                           observation.dict['detrended_series']['flux'],
                                                                           observation.dict['detrended_series']['flux_unc'],
                                                                           observation.dict['detrended_series']['model'],
                                                                           number_of_free_parameters)

        results_copy = copy_dict(fitting.results)
        if not return_traces:
            for parameter in results_copy['parameters']:
                results_copy['parameters'][parameter]['trace'] = None

        if output_folder:

            if os.path.isdir(output_folder):
                shutil.rmtree(output_folder)
            os.mkdir(output_folder)

            save_dict(results_copy, os.path.join(output_folder, 'global_results.pickle'))

            fitting.save_results(os.path.join(output_folder, 'global_results.txt'))

            fitting.plot_corner(os.path.join(output_folder, 'global_correlations.pdf'))

            fitting.plot_traces(os.path.join(output_folder, 'global_traces.pdf'))

            for observation in self.observations:

                cols = [
                    ['# variable'],
                    ['fix/fit'],
                    ['value'],
                    ['uncertainty'],
                    ['initial'],
                    ['min.allowed'],
                    ['max.allowed']
                ]

                for parameter_index in observation.parameters_map:
                    parameter_name = names[parameter_index].split('::')[0]

                    cols[0].append(observation.dict['parameters'][parameter_name]['name'])
                    if observation.dict['parameters'][parameter_name]['initial'] is None:
                        cols[1].append('fix')
                        cols[2].append(observation.dict['parameters'][parameter_name]['print_value'])
                        cols[3].append('-- --')
                        cols[4].append('--')
                        cols[5].append('--')
                        cols[6].append('--')
                    else:
                        cols[1].append('fit')
                        cols[2].append(observation.dict['parameters'][parameter_name]['print_value'])
                        cols[3].append(
                            '-{0} +{1}'.format(
                                observation.dict['parameters'][parameter_name]['print_m_error'],
                                observation.dict['parameters'][parameter_name]['print_p_error'])
                        )
                        cols[4].append(str(observation.dict['parameters'][parameter_name]['initial']))
                        cols[5].append(str(observation.dict['parameters'][parameter_name]['min_allowed']))
                        cols[6].append(str(observation.dict['parameters'][parameter_name]['max_allowed']))

                for col in cols:
                    col_length = np.max([len(ff) for ff in col])
                    for ff in range(len(col)):
                        col[ff] = col[ff] + ' ' * (col_length - len(col[ff]))

                lines = []

                for row in range(len(cols[0])):
                    lines.append('  '.join([col[row] for col in cols]))

                lines.append('')
                lines.append('#Filter: {0}'.format(observation.filter_id))
                lines.append('#Epoch: {0}'.format(observation.epoch))
                lines.append('#Number of outliers removed: {0}'.format(observation.dict['data_conversion_info']['outliers']))
                lines.append('#Uncertainties scale factor: {0}'.format(observation.dict['data_conversion_info']['scale_factor']))

                lines.append('')
                lines.append('#Residuals:')
                lines.append('#Mean: {0}'.format(observation.dict['statistics']['res_mean']))
                lines.append('#STD: {0}'.format(observation.dict['statistics']['res_std']))
                lines.append('#RMS: {0}'.format(observation.dict['statistics']['res_rms']))
                lines.append('#Chi squared: {0}'.format(observation.dict['statistics']['res_chi_sqr']))
                lines.append('#Reduced chi squared: {0}'.format(observation.dict['statistics']['res_red_chi_sqr']))
                lines.append('#Max auto-correlation: {0}'.format(observation.dict['statistics']['res_max_autocorr']))
                lines.append('#Max auto-correlation flag: {0}'.format(observation.dict['statistics']['res_max_autocorr_flag']))
                lines.append('#Shapiro test: {0}'.format(observation.dict['statistics']['res_shapiro']))
                lines.append('#Shapiro test flag: {0}'.format(observation.dict['statistics']['res_shapiro_flag']))

                lines.append('')
                lines.append('#Detrended Residuals:')
                lines.append('#Mean: {0}'.format(observation.dict['detrended_statistics']['res_mean']))
                lines.append('#STD: {0}'.format(observation.dict['detrended_statistics']['res_std']))
                lines.append('#RMS: {0}'.format(observation.dict['detrended_statistics']['res_rms']))
                lines.append('#Chi squared: {0}'.format(observation.dict['detrended_statistics']['res_chi_sqr']))
                lines.append('#Reduced chi squared: {0}'.format(observation.dict['detrended_statistics']['res_red_chi_sqr']))
                lines.append('#Max auto-correlation: {0}'.format(observation.dict['detrended_statistics']['res_max_autocorr']))
                lines.append('#Max auto-correlation flag: {0}'.format(observation.dict['detrended_statistics']['res_max_autocorr_flag']))
                lines.append('#Shapiro test: {0}'.format(observation.dict['detrended_statistics']['res_shapiro']))
                lines.append('#Shapiro test flag: {0}'.format(observation.dict['detrended_statistics']['res_shapiro_flag']))

                w = open(os.path.join(output_folder, '{0}_results.txt'.format(observation.obs_id)), 'w')
                w.write('\n'.join(lines))
                w.close()

                observation_copy_dict = copy_dict(observation.dict)
                save_dict(observation_copy_dict, os.path.join(output_folder, '{0}_results.pickle'.format(observation.obs_id)))

                plot_transit_fitting_models(observation.dict, os.path.join(output_folder, '{0}_lightcurve.pdf'.format(observation.obs_id)))

        results_copy['observations'] = {}
        for observation in self.observations:
            results_copy['observations'][observation.obs_id] = observation.dict

        return results_copy

    def eclipse_fitting(self, output_folder=None,

                        iterations=None, walkers=None, burn_in=None,

                        fit_individual_fp_over_fs=True,
                        fit_rp_over_rs=False,
                        fit_individual_rp_over_rs=False,
                        fit_sma_over_rs=False, fit_inclination=False,
                        fit_mid_time=True, fit_period=False,
                        fit_individual_times=True,

                        fit_rp_over_rs_limits=[0.001, 1000],
                        fit_fp_over_fs_limits=[0.001, 1000],
                        fit_sma_over_rs_limits=[0.001, 1000],
                        fit_inclination_limits=[10.0, 90.0],
                        fit_mid_time_limits=[-0.2, 0.2],
                        fit_period_limits=[0.001, 1000],

                        counter='Eclipse fitting',
                        optimise_initial_parameters_trials=3,
                        scale_uncertainties=False,
                        filter_outliers=False,
                        optimiser='emcee',
                        strech_prior=0.2,
                        walkers_spread=0.01,
                        return_traces=True,
                        verbose=False,
                        ):

        for observation in self.observations:
            observation.reset(self)
            if observation.dict['model_info']['observation_type'] == 'transit':
                raise PyLCInputError('You need to add only eclipse observation to proceed.')

        if len(self.observations) == 0:
            raise PyLCInputError('You need to add at least one eclipse observation to proceed.')

        for observation_num, observation in enumerate(self.observations):
            observation.obs_id = 'obs{0}'.format(observation_num)
            observation.dict['model_info']['obs_id'] = observation.obs_id

        fit_rp_over_rs_limits = [self.rp_over_rs * fit_rp_over_rs_limits[0], self.rp_over_rs * fit_rp_over_rs_limits[1]]
        fit_sma_over_rs_limits = [self.sma_over_rs * fit_sma_over_rs_limits[0], self.sma_over_rs * fit_sma_over_rs_limits[1]]
        fit_period_limits = [self.period * fit_period_limits[0], self.period * fit_period_limits[1]]

        # separate filters and epochs

        unique_epochs = []
        for observation in self.observations:
            if observation.epoch not in unique_epochs:
                unique_epochs.append(observation.epoch)

        if len(unique_epochs) == 1:
            fit_individual_times = False
            if fit_period:
                raise PyLCInputError('Period cannot be fitted only on one epoch.')

        if fit_individual_times and fit_period:
            raise PyLCInputError('Period and individual mid times cannot be fitted simultaneously.')

        unique_filters = {}
        for observation in self.observations:
            if observation.filter_id not in unique_filters:
                unique_filters[observation.filter_id] = \
                    observation.dict['model_info']['fp_over_fs']

        if len(unique_filters) == 1:
            fit_individual_fp_over_fs = False
            fit_individual_rp_over_rs = False

        # add parameters

        names = []
        print_names = []
        limits1 = []
        limits2 = []
        initial = []
        logspace = []

        # de-trending parameters

        for observation in self.observations:

            for coefficient in range(observation.number_of_detrending_parameters):

                if len(self.observations) > 1:
                    names.append('{0}::{1}'.format(observation.names[coefficient], observation.obs_id))
                    print_names.append('{0}::{1}'.format(observation.print_names[coefficient], observation.obs_id))
                else:
                    names.append('{0}'.format(observation.names[coefficient]))
                    print_names.append('{0}'.format(observation.print_names[coefficient]))

                initial.append(observation.initial[coefficient])
                limits1.append(observation.limits1[coefficient])
                limits2.append(observation.limits2[coefficient])
                logspace.append(False)

                observation.add_to_parameters_map(len(names) - 1)

        # fp_over_fs

        if fit_individual_fp_over_fs:
            for phot_filter in unique_filters:
                names.append('fp_over_fs::{0}'.format(phot_filter))
                print_names.append(r'$F_\mathrm{{p}}/F_*$::{0}'.format(phot_filter))
                initial.append(unique_filters[phot_filter])
                limits1.append(unique_filters[phot_filter] * fit_fp_over_fs_limits[0])
                limits2.append(unique_filters[phot_filter] * fit_fp_over_fs_limits[1])
                logspace.append(True)
                for observation in self.observations:
                    if observation.filter_id == phot_filter:
                        observation.add_to_parameters_map(len(names) - 1)
        else:
            phot_filter = list(unique_filters.keys())[0]
            names.append('fp_over_fs')
            print_names.append(r'$F_\mathrm{{p}}/F_*$')
            initial.append(unique_filters[phot_filter])
            limits1.append(unique_filters[phot_filter] * fit_fp_over_fs_limits[0])
            limits2.append(unique_filters[phot_filter] * fit_fp_over_fs_limits[1])
            logspace.append(True)
            for observation in self.observations:
                observation.add_to_parameters_map(len(names) - 1)

        # rp_over_rs

        if fit_rp_over_rs and fit_individual_rp_over_rs:
            for phot_filter in unique_filters:
                names.append('rp_ovr_rs::{0}'.format(phot_filter))
                print_names.append(r'$R_\mathrm{{p}}/R_*$::{0}'.format(phot_filter))
                initial.append(self.rp_over_rs)
                limits1.append(fit_rp_over_rs_limits[0])
                limits2.append(fit_rp_over_rs_limits[1])
                logspace.append(True)
                for observation in self.observations:
                    if observation.filter_id == phot_filter:
                        observation.add_to_parameters_map(len(names) - 1)

            global_parameters_names = []
            global_parameters_print_names = []
            global_parameters_initial = []
            global_parameters_fit = []
            global_parameters_limits = []
            global_parameters_logspace = []

        else:
            global_parameters_names = ['rp_over_rs']
            global_parameters_print_names = [r'$R_\mathrm{p}/R_*$']
            global_parameters_initial = [self.rp_over_rs]
            global_parameters_fit = [fit_rp_over_rs]
            global_parameters_limits = [fit_rp_over_rs_limits]
            global_parameters_logspace = [True]

        # orbital parameters

        global_parameters_names += ['period', 'sma_over_rs', 'eccentricity', 'inclination', 'periastron']
        global_parameters_print_names += [r'$P$', r'$a$', r'$e$', r'$i$', r'$\omega$']
        global_parameters_initial += [self.period, self.sma_over_rs, self.eccentricity, self.inclination, self.periastron]
        global_parameters_fit += [fit_period, fit_sma_over_rs, False, fit_inclination, False]
        global_parameters_limits += [fit_period_limits, fit_sma_over_rs_limits, False, fit_inclination_limits, False]
        global_parameters_logspace += [True, True, False, True, False]

        if not fit_mid_time or not fit_individual_times or len(unique_epochs) == 1:

            test_epochs = np.array([])
            test_epochs_weights = np.array([])

            for observation in self.observations:
                test_epochs = np.append(test_epochs,
                                        np.ones_like(observation.dict['input_series']['flux_unc']) * observation.epoch)
                norm_errors = observation.dict['input_series']['flux_unc'] / observation.dict['input_series']['flux']
                test_epochs_weights = np.append(test_epochs_weights, 1 / (norm_errors * norm_errors))

            new_epoch = np.round(np.sum(test_epochs * test_epochs_weights) / np.sum(test_epochs_weights), 0)
            new_mid_time = self.eclipse_mid_time() + new_epoch * self.period

            global_parameters_names += ['eclipse_mid_time']
            global_parameters_print_names += [r'$T_\mathrm{ec}$']
            global_parameters_initial += [new_mid_time]
            global_parameters_fit += [fit_mid_time]
            global_parameters_limits += [[new_mid_time + fit_mid_time_limits[0], new_mid_time + fit_mid_time_limits[1]]]
            global_parameters_logspace += [False]

        for global_parameter in range(len(global_parameters_names)):
            names.append(global_parameters_names[global_parameter])
            print_names.append(global_parameters_print_names[global_parameter])
            initial.append(global_parameters_initial[global_parameter])
            if global_parameters_fit[global_parameter]:
                limits1.append(global_parameters_limits[global_parameter][0])
                limits2.append(global_parameters_limits[global_parameter][1])
            else:
                limits1.append(np.nan)
                limits2.append(np.nan)
            logspace.append(global_parameters_logspace[global_parameter])

            for observation in self.observations:
                observation.add_to_parameters_map(len(names) - 1)

        # individual time parameters
        eclipse_mid_time = self.eclipse_mid_time()
        if fit_mid_time and fit_individual_times and len(unique_epochs) > 1:

            for epoch in unique_epochs:

                names.append('eclipse_mid_time::{0}'.format(epoch))
                print_names.append(r'$T_\mathrm{{ec}}$::{0}'.format(epoch))
                initial.append(eclipse_mid_time + epoch * self.period)
                limits1.append(eclipse_mid_time + epoch * self.period + fit_mid_time_limits[0])
                limits2.append(eclipse_mid_time + epoch * self.period + fit_mid_time_limits[1])
                logspace.append(False)

                for observation in self.observations:
                    if observation.epoch == epoch:
                        observation.add_to_parameters_map(len(names) - 1)

        initial = np.array(initial)
        limits1 = np.array(limits1)
        limits2 = np.array(limits2)
        logspace = np.array(logspace, dtype=bool)

        # single observation tests

        for observation in self.observations:

            def single_observation_full_model(indices, *model_variables):
                return observation.full_model(indices, np.array(model_variables))

            fitting = Fitting(observation.non_outliers_map,
                              observation.dict['input_series']['flux'],
                              observation.dict['input_series']['flux_unc'],
                              single_observation_full_model, initial, limits1, limits2,
                              logspace=logspace,
                              data_x_name='time', data_y_name='flux',
                              data_x_print_name=r'$t_{BJD_{TDB}}$', data_y_print_name='Relative Flux',
                              parameters_names=names, parameters_print_names=print_names,
                              walkers=walkers, iterations=iterations, burn_in=burn_in,
                              counter=counter,
                              optimise_initial_parameters=True,
                              optimise_initial_parameters_trials=optimise_initial_parameters_trials,
                              scale_uncertainties=scale_uncertainties,
                              filter_outliers=filter_outliers,
                              optimiser='curve_fit',
                              strech_prior=strech_prior,
                              walkers_spread=walkers_spread,
                              )

            fitting._prefit(verbose=verbose)

            initial[observation.parameters_map] = fitting.results['prefit']['initials'][observation.parameters_map]
            observation.dict['input_series']['flux_unc'] *= fitting.results['prefit']['scale_factor']
            observation.dict['data_conversion_info']['scale_factor'] = fitting.results['prefit']['scale_factor']
            observation.dict['data_conversion_info']['notes'].append('Flux_unc multiplied by {0}'.format(fitting.results['prefit']['scale_factor']))
            observation.dict['data_conversion_info']['outliers'] = len(np.where(fitting.results['prefit']['outliers_map'])[0])
            observation.dict['data_conversion_info']['notes'].append('{0} outliers removed'.format(len(np.where(fitting.results['prefit']['outliers_map'])[0])))
            observation.dict['data_conversion_info']['non_outliers_map'] = np.where(~fitting.results['prefit']['outliers_map'])[0]
            observation.non_outliers_map = np.where(~fitting.results['prefit']['outliers_map'])[0]

            print()
            print('Observation: ', observation.obs_id)
            print('Filter: ', observation.filter_id)
            print('Epoch: ', observation.epoch)
            print('Data-points excluded: ', len(np.where(fitting.results['prefit']['outliers_map'])[0]))
            print('Scaling uncertainties by: ', fitting.results['prefit']['scale_factor'])

            observation.clear_outliers()

        model_time = np.array([])
        model_flux = np.array([])
        model_flux_unc = np.array([])
        for observation in self.observations:
            model_time = np.append(model_time, observation.dict['input_series']['time'])
            model_flux = np.append(model_flux, observation.dict['input_series']['flux'])
            model_flux_unc = np.append(model_flux_unc, observation.dict['input_series']['flux_unc'])

        def detrend_model(model_time, *model_variables):
            model = np.array([])
            for observation in self.observations:
                model = np.append(model, observation.splitted_model(None, np.array(model_variables))[0])
            return model

        def full_model(model_time, *model_variables):
            model = np.array([])
            for observation in self.observations:
                model = np.append(model, observation.full_model(None, np.array(model_variables)))
            return model

        fitting = Fitting(model_time, model_flux, model_flux_unc,
                          full_model, initial, limits1, limits2,
                          data_x_name='time', data_y_name='flux',
                          data_x_print_name=r'$t_{BJD_{TDB}}$', data_y_print_name='Relative Flux',
                          parameters_names=names, parameters_print_names=print_names,
                          walkers=walkers, iterations=iterations, burn_in=burn_in,
                          counter=counter,
                          optimise_initial_parameters=True,
                          scale_uncertainties=False,
                          filter_outliers=False,
                          optimiser=optimiser,
                          strech_prior=strech_prior,
                          walkers_spread=walkers_spread,
                          )

        print('\nOptimising initial parameters...')

        fitting.run(verbose=verbose)

        fitting.results['settings'] = {
            'output_folder': output_folder,
            'iterations': fitting.iterations,
            'walkers': fitting.walkers,
            'burn_in': fitting.burn_in,
            'strech_prior': strech_prior,
            'walkers_spread': walkers_spread,
            'optimise_initial_parameters': True,
            'scale_uncertainties': scale_uncertainties,
            'filter_outliers': filter_outliers,
            'optimiser': optimiser,
            'data_x_name': 'time',
            'data_x_print_name': r't_${BJD_{TDB}}$',
            'data_y_name': 'flux',
            'data_y_print_name': 'Relative Flux',
            'fit_individual_fp_over_fs': fit_individual_fp_over_fs,
            'fit_rp_over_rs': fit_rp_over_rs,
            'fit_individual_rp_over_rs': fit_individual_rp_over_rs,
            'fit_sma_over_rs': fit_sma_over_rs,
            'fit_inclination': fit_inclination,
            'fit_mid_time': fit_mid_time,
            'fit_period': fit_period,
            'fit_individual_times': fit_individual_times,
            'fit_fp_over_fs_limits': fit_fp_over_fs_limits,
            'fit_rp_over_rs_limits': fit_rp_over_rs_limits,
            'fit_sma_over_rs_limits': fit_sma_over_rs_limits,
            'fit_inclination_limits': fit_inclination_limits,
            'fit_mid_time_limits': fit_mid_time_limits,
            'fit_period_limits': fit_period_limits,
        }

        del fitting.results['prefit']
        del fitting.results['original_data']

        trend = detrend_model(model_time, *fitting.results['parameters_final'])

        fitting.results['output_series']['trend'] = trend

        fitting.results['detrended_series'] = {
            'time': model_time,
            'flux': fitting.results['input_series']['flux'] / trend,
            'flux_unc': fitting.results['input_series']['flux_unc'] / trend,
            'model': fitting.results['output_series']['model'] / trend,
            'residuals': fitting.results['output_series']['residuals'] / trend
        }

        fitting.results['detrended_statistics'] = residual_statistics(fitting.results['detrended_series']['time'],
                                                                      fitting.results['detrended_series']['flux'],
                                                                      fitting.results['detrended_series']['flux_unc'],
                                                                      fitting.results['detrended_series']['model'],
                                                                      len(fitting.fitted_parameters))

        for observation in self.observations:

            observation.dict['parameters'] = {}
            number_of_free_parameters = 0
            for parameter_index in observation.parameters_map:
                parameter_name = names[parameter_index]
                parameter_sub_name = parameter_name.split('::')[0]
                parameter_data = copy_dict(fitting.results['parameters'][parameter_name])
                parameter_data['trace'] = None
                parameter_data['name'] = parameter_sub_name
                parameter_data['print_name'] = parameter_data['print_name'].split(':')[0]
                observation.dict['parameters'][parameter_sub_name] = parameter_data
                if fitting.results['parameters'][parameter_name]['initial'] is not None:
                    number_of_free_parameters += 0

            observation.dict['output_series']['model'] = observation.full_model(None, fitting.results['parameters_final'])
            observation.dict['output_series']['trend'] = observation.splitted_model(None, fitting.results['parameters_final'])[0]
            observation.dict['output_series']['residuals'] = observation.dict['input_series']['flux'] - observation.dict['output_series']['model']
            observation.dict['detrended_series']['time'] = observation.dict['input_series']['time']
            observation.dict['detrended_series']['flux'] = observation.dict['input_series']['flux'] / observation.dict['output_series']['trend']
            observation.dict['detrended_series']['flux_unc'] = observation.dict['input_series']['flux_unc'] / observation.dict['output_series']['trend']
            observation.dict['detrended_series']['model'] = observation.dict['output_series']['model'] / observation.dict['output_series']['trend']
            observation.dict['detrended_series']['residuals'] = observation.dict['detrended_series']['flux'] - observation.dict['detrended_series']['model']

            observation.dict['statistics'] = residual_statistics(observation.dict['input_series']['time'],
                                                                 observation.dict['input_series']['flux'],
                                                                 observation.dict['input_series']['flux_unc'],
                                                                 observation.dict['output_series']['model'],
                                                                 number_of_free_parameters)

            observation.dict['detrended_statistics'] = residual_statistics(observation.dict['detrended_series']['time'],
                                                                           observation.dict['detrended_series']['flux'],
                                                                           observation.dict['detrended_series']['flux_unc'],
                                                                           observation.dict['detrended_series']['model'],
                                                                           number_of_free_parameters)

        results_copy = copy_dict(fitting.results)
        if not return_traces:
            for parameter in results_copy['parameters']:
                results_copy['parameters'][parameter]['trace'] = None

        if output_folder:

            if os.path.isdir(output_folder):
                shutil.rmtree(output_folder)
            os.mkdir(output_folder)

            save_dict(results_copy, os.path.join(output_folder, 'global_results.pickle'))

            fitting.save_results(os.path.join(output_folder, 'global_results.txt'))

            fitting.plot_corner(os.path.join(output_folder, 'global_correlations.pdf'))

            fitting.plot_traces(os.path.join(output_folder, 'global_traces.pdf'))

            for observation in self.observations:

                cols = [
                    ['# variable'],
                    ['fix/fit'],
                    ['value'],
                    ['uncertainty'],
                    ['initial'],
                    ['min.allowed'],
                    ['max.allowed']
                ]

                for parameter_index in observation.parameters_map:
                    parameter_name = names[parameter_index].split('::')[0]

                    cols[0].append(observation.dict['parameters'][parameter_name]['name'])
                    if observation.dict['parameters'][parameter_name]['initial'] is None:
                        cols[1].append('fix')
                        cols[2].append(observation.dict['parameters'][parameter_name]['print_value'])
                        cols[3].append('-- --')
                        cols[4].append('--')
                        cols[5].append('--')
                        cols[6].append('--')
                    else:
                        cols[1].append('fit')
                        cols[2].append(observation.dict['parameters'][parameter_name]['print_value'])
                        cols[3].append(
                            '-{0} +{1}'.format(
                                observation.dict['parameters'][parameter_name]['print_m_error'],
                                observation.dict['parameters'][parameter_name]['print_p_error'])
                        )
                        cols[4].append(str(observation.dict['parameters'][parameter_name]['initial']))
                        cols[5].append(str(observation.dict['parameters'][parameter_name]['min_allowed']))
                        cols[6].append(str(observation.dict['parameters'][parameter_name]['max_allowed']))

                for col in cols:
                    col_length = np.max([len(ff) for ff in col])
                    for ff in range(len(col)):
                        col[ff] = col[ff] + ' ' * (col_length - len(col[ff]))

                lines = []

                for row in range(len(cols[0])):
                    lines.append('  '.join([col[row] for col in cols]))

                lines.append('')
                lines.append('#Filter: {0}'.format(observation.filter_id))
                lines.append('#Epoch: {0}'.format(observation.epoch))
                lines.append('#Number of outliers removed: {0}'.format(observation.dict['data_conversion_info']['outliers']))
                lines.append('#Uncertainties scale factor: {0}'.format(observation.dict['data_conversion_info']['scale_factor']))

                lines.append('')
                lines.append('#Residuals:')
                lines.append('#Mean: {0}'.format(observation.dict['statistics']['res_mean']))
                lines.append('#STD: {0}'.format(observation.dict['statistics']['res_std']))
                lines.append('#RMS: {0}'.format(observation.dict['statistics']['res_rms']))
                lines.append('#Chi squared: {0}'.format(observation.dict['statistics']['res_chi_sqr']))
                lines.append('#Reduced chi squared: {0}'.format(observation.dict['statistics']['res_red_chi_sqr']))
                lines.append('#Max auto-correlation: {0}'.format(observation.dict['statistics']['res_max_autocorr']))
                lines.append('#Max auto-correlation flag: {0}'.format(observation.dict['statistics']['res_max_autocorr_flag']))
                lines.append('#Shapiro test: {0}'.format(observation.dict['statistics']['res_shapiro']))
                lines.append('#Shapiro test flag: {0}'.format(observation.dict['statistics']['res_shapiro_flag']))

                lines.append('')
                lines.append('#Detrended Residuals:')
                lines.append('#Mean: {0}'.format(observation.dict['detrended_statistics']['res_mean']))
                lines.append('#STD: {0}'.format(observation.dict['detrended_statistics']['res_std']))
                lines.append('#RMS: {0}'.format(observation.dict['detrended_statistics']['res_rms']))
                lines.append('#Chi squared: {0}'.format(observation.dict['detrended_statistics']['res_chi_sqr']))
                lines.append('#Reduced chi squared: {0}'.format(observation.dict['detrended_statistics']['res_red_chi_sqr']))
                lines.append('#Max auto-correlation: {0}'.format(observation.dict['detrended_statistics']['res_max_autocorr']))
                lines.append('#Max auto-correlation flag: {0}'.format(observation.dict['detrended_statistics']['res_max_autocorr_flag']))
                lines.append('#Shapiro test: {0}'.format(observation.dict['detrended_statistics']['res_shapiro']))
                lines.append('#Shapiro test flag: {0}'.format(observation.dict['detrended_statistics']['res_shapiro_flag']))

                w = open(os.path.join(output_folder, '{0}_results.txt'.format(observation.obs_id)), 'w')
                w.write('\n'.join(lines))
                w.close()

                observation_copy_dict = copy_dict(observation.dict)
                save_dict(observation_copy_dict, os.path.join(output_folder, '{0}_results.pickle'.format(observation.obs_id)))

                plot_transit_fitting_models(observation.dict, os.path.join(output_folder, '{0}_lightcurve.pdf'.format(observation.obs_id)))

        results_copy['observations'] = {}
        for observation in self.observations:
            results_copy['observations'][observation.obs_id] = observation.dict

        return results_copy

    def performance_report(self, array_size):

        duration = self.transit_duration()
        time_array = np.linspace(-duration/2 - 0.0000001, duration/2 + 0.0000001, array_size) + self.mid_time

        limb_darkening_coefficients = self.exotethys('COUSINS_R')

        ref = transit(limb_darkening_coefficients, self.rp_over_rs, self.period,
                           self.sma_over_rs, self.eccentricity, self.inclination,
                           self.periastron, self.mid_time,
                           time_array, precision='ref')

        print('Performance report')
        print('Planet: ', self.name)
        print('Array size: ', array_size)

        for precision in range(2, 11):
            plc_times = []
            for i in range(100000):
                t0 = sys_time.time()
                for j in range(10):
                    test = transit(limb_darkening_coefficients, self.rp_over_rs, self.period,
                                   self.sma_over_rs, self.eccentricity, self.inclination,
                                   self.periastron, self.mid_time,
                                   time_array, precision=precision)
                plc_times.append(100*(sys_time.time() - t0))
                if np.sum(plc_times) > 500:
                    break

            t1 = np.median(plc_times)
            err1 = np.max(np.abs(1000000*(test-ref)))

            print('precision = {0} : {1:.2e} ppm, {2:.2f} ms '.format(precision, err1, t1))


def get_planet(name):

    all_data = exoclock.get_planet(name)

    planet = Planet(
        all_data['name'],
        all_data['star']['ra_deg'],
        all_data['star']['dec_deg'],
        all_data['planet']['logg'],
        all_data['planet']['teff'],
        all_data['planet']['meta'],
        all_data['planet']['rp_over_rs'],
        all_data['planet']['ephem_period'],
        all_data['planet']['sma_over_rs'],
        all_data['planet']['eccentricity'],
        all_data['planet']['inclination'],
        all_data['planet']['periastron'],
        all_data['planet']['ephem_mid_time'],
    )

    planet.all_data = all_data

    return planet


def get_all_planets():
    return exoclock.get_all_planets()
