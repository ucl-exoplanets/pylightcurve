
__all__ = ['plc_data', 'all_filters', 'add_filter']

import os
import glob
import numpy as np
import pickle
import shutil
import datetime

from exotethys import ls_database
from scipy.interpolate import interp1d

from pylightcurve.processes.files import open_dict, save_dict, open_dict_online
from .errors import *
from .__databases_setup__ import _setup_database

from pylightcurve import __version__

databases_file = '__databases__.pickle'
package_name = 'pylightcurve4'
github_link = 'https://github.com/ucl-exoplanets/pylightcurve/raw/master/pylightcurve/__databases__.pickle?raw=true'


class PlcData:

    def __init__(self, _reset=False, _test=False):

        self.package_name = package_name
        self.version = '.'.join(__version__.split('.')[:2])

        self.build_in_databases_file_path = os.path.join(os.path.abspath(os.path.dirname(__file__)), databases_file)

        self.databases_directory_path = os.path.join(os.path.abspath(os.path.expanduser('~')),
                                                     '.{0}'.format(self.package_name))

        self.databases_file_path = os.path.join(self.databases_directory_path, databases_file)
        self.databases_file_path_new = os.path.join(self.databases_directory_path, databases_file + '_new')

        # initiate databases

        if not os.path.isdir(self.databases_directory_path):
            os.mkdir(self.databases_directory_path)

        # test = False
        # if test:
        #     shutil.copy(self.build_in_databases_file_path, self.databases_file_path)

        else:

            if not os.path.isfile(self.databases_file_path):
                shutil.copy(self.build_in_databases_file_path, self.databases_file_path)

            # check for updates in the databases (identified on github)

            test_online_db = open_dict_online(github_link)
            test_local_db = open_dict(self.databases_file_path)
            if test_online_db and test_online_db != test_local_db:
                save_dict(test_online_db, self.databases_file_path)

        # load databases

        self.databases = open_dict(self.databases_file_path)

        self.exotethys_loaded = _setup_database(self, 'exotethys')
        self.ephemerides_loaded = _setup_database(self, 'ephemerides')
        self.photometry_loaded = _setup_database(self, 'photometry')
        self.catalogues_loaded = _setup_database(self, 'catalogues')
        self.ut_loaded = _setup_database(self, 'ut')

        self.leap_seconds_data = None
        self.earth_rotation_data = None
        self.barycenter_data = None
        self.sun_data = None
        self.moon_data = None
        self.ecc_data = None

    def exotethys(self):
        return self.exotethys_loaded

    def ephemeris(self):
        return self.ephemerides_loaded

    def photometry(self):
        return self.photometry_loaded

    def catalogues(self):
        return self.catalogues_loaded

    def ut(self):
        return self.ut_loaded

    def all_filters(self):

        return [os.path.split(ff)[1].split('.')[0]
                for ff in glob.glob(os.path.join(self.photometry(), '*'))]

    def get_filter(self, filter_name, wlrange=None):

        if filter_name not in self.all_filters():
            raise PyLCInputError('\n\n{0} is not available. To see all the available bandpasses, type: '
                                 '\n    plc.all_filters()'
                                 '\nAlternatively you can define your own bandpass as follows: '
                                 '\n    plc.add_filter(filter_name, passband_file) '
                                 '\nThe passband file should be a txt file containing two columns '
                                 'separated by space or tab:'
                                 '\ncolumn 1: wavelength in A, '
                                 '\ncolumn 2: total throughput in electrons/photons.'.format(
                filter_name, ','.join(all_filters())))

        filter_data = np.loadtxt(os.path.join(self.photometry(), filter_name + '.pass'))

        passband_wlrange = [min(filter_data[:, 0]), max(filter_data[:, 0])]
        if not wlrange:
            wlrange = passband_wlrange
        else:
            wlrange = [min(wlrange), max(wlrange)]

        if wlrange[0] < passband_wlrange[0] or wlrange[1] > passband_wlrange[1]:
            raise PyLCInputError('Wavelength is not compatible with the {0} passband. '
                                 'Wavelength range requested: {1}, '
                                 'Passband wavelength range: {2}.'.format(filter_name, wlrange, passband_wlrange))

        return filter_data[np.where((filter_data[:, 0] <= wlrange[1]) * (filter_data[:, 0] >= wlrange[0]))]

    def add_filter(self, filter_name, passband_file):

        try:
            _ = np.loadtxt(passband_file)
        except:
            raise PyLCInputError('Wrong passband format or file path.')

        shutil.copy(passband_file, os.path.join(self.photometry(), filter_name + '.pass'))

    def reset_exotethys_data(self):
        for file in glob.glob(os.path.join(ls_database()[0], '*', '*')):
            try:
                _ = open_dict(file)
            except pickle.UnpicklingError:
                os.remove(file)

    def ecc(self):

        if not self.ecc_data:
            stars = open_dict(os.path.join(self.catalogues(), 'ecc_stars.pickle'))
            planets = open_dict(os.path.join(self.catalogues(), 'ecc_planets.pickle'))

            hosts = {stars[ff]['simbad_id']: stars[ff]['planets'] for ff in stars}

            def _flat_name(name):

                flat_name_list = [
                    [' ', ''],
                    ['-', ''],
                    ['cancri', 'cnc'],
                    ['hatp10', 'wasp11'],
                    ['wasp40', 'hatp27'],
                    ['wasp51', 'hatp30'],
                    ['wasp86', 'kelt12'],
                    ['kelt22', 'wasp173'],
                ]

                name = name.lower()

                for char in flat_name_list:
                    name = name.replace(char[0], char[1])

                return name

            flats = {_flat_name(ff): stars[ff]['simbad_id'] for ff in stars}
            for planet in planets:
                star = planets[planet]['star']
                flats[_flat_name(planet)] = star
                flats[_flat_name(star)] = star

            self.ecc_data = {'stars': stars, 'planets': planets, 'hosts': hosts, 'flats': flats}

        return self.ecc_data

    def barycentre(self, jd_utc):

        if not self.barycenter_data:
            self.barycenter_data = open_dict(os.path.join(self.ephemeris(), 'bjd_dict.pickle'))

        bjd_dict = self.barycenter_data[int(jd_utc)]

        ssb_t = bjd_dict['t']
        ssb_ra = bjd_dict['ra']
        ssb_dec = bjd_dict['dec']
        ssb_d = bjd_dict['d']
        ssb_dt = bjd_dict['dt']

        ssb_ra = interp1d(ssb_t, ssb_ra, kind='cubic')(jd_utc)
        ssb_dec = interp1d(ssb_t, ssb_dec, kind='cubic')(jd_utc)
        ssb_d = interp1d(ssb_t, ssb_d, kind='cubic')(jd_utc)
        ssb_dt = interp1d(ssb_t, ssb_dt, kind='cubic')(jd_utc)

        return ssb_ra, ssb_dec, ssb_d, ssb_dt

    def heliocentre(self, jd_utc):

        if not self.sun_data:
            self.sun_data = open_dict(os.path.join(self.ephemeris(), 'hjd_dict.pickle'))

        hjd_dict = self.sun_data[int(jd_utc)]

        ssb_t = hjd_dict['t']
        ssb_ra = hjd_dict['ra']
        ssb_dec = hjd_dict['dec']
        ssb_d = hjd_dict['d']
        ssb_dt = hjd_dict['dt']

        ssb_ra = interp1d(ssb_t, ssb_ra, kind='cubic')(jd_utc)
        ssb_dec = interp1d(ssb_t, ssb_dec, kind='cubic')(jd_utc)
        ssb_d = interp1d(ssb_t, ssb_d, kind='cubic')(jd_utc)
        ssb_dt = interp1d(ssb_t, ssb_dt, kind='cubic')(jd_utc)

        return ssb_ra, ssb_dec, ssb_d, ssb_dt

    def mooncentre(self, jd_utc):

        if not self.moon_data:
            self.moon_data = open_dict(os.path.join(self.ephemeris(), 'moon_dict.pickle'))

        moon_dict = self.moon_data[int(jd_utc)]

        ssb_t = moon_dict['t']
        ssb_ra = moon_dict['ra']
        ssb_dec = moon_dict['dec']
        ssb_d = moon_dict['d']
        ssb_dt = moon_dict['dt']

        ssb_ra = interp1d(ssb_t, ssb_ra, kind='cubic')(jd_utc)
        ssb_dec = interp1d(ssb_t, ssb_dec, kind='cubic')(jd_utc)
        ssb_d = interp1d(ssb_t, ssb_d, kind='cubic')(jd_utc)
        ssb_dt = interp1d(ssb_t, ssb_dt, kind='cubic')(jd_utc)

        return ssb_ra, ssb_dec, ssb_d, ssb_dt

    def leap_seconds(self, utc):

        if not self.leap_seconds_data:
            self.leap_seconds_data = [[datetime.datetime(int(line[3]), int(line[2]), int(line[1])), line[4]] for line in
                                      np.loadtxt(open(os.path.join(self.ut(), 'leap_seconds.txt')))]

        ls = self.leap_seconds_data[-1][1]
        for check in range(1, len(self.leap_seconds_data)):
            if utc < self.leap_seconds_data[check][0]:
                ls = self.leap_seconds_data[check - 1][1]
                break

        return ls

    def earth_rotation(self, jd_utc):

        if not self.earth_rotation_data:

            earth_rotation_data_x = []
            earth_rotation_data_y = []
            for ff in open(os.path.join(self.ut(), 'earth_rotation.txt')).readlines():
                if ff[58:68].replace(' ', '') != '':
                    earth_rotation_data_x.append(float(ff[7:12]) + 2400000.5)
                    earth_rotation_data_y.append(float(ff[58:68].replace(' ', '')))

            self.earth_rotation_data = interp1d(earth_rotation_data_x, earth_rotation_data_y, kind='cubic')

        return float(self.earth_rotation_data(jd_utc))


plc_data = PlcData()

def all_filters():

    xx = []
    for i in plc_data.all_filters():
        if 'hst_wfc3_g' not in i:
            xx.append(i)

    return(xx)

def add_filter(*args, **kwargs):
    plc_data.add_filter(*args, **kwargs)