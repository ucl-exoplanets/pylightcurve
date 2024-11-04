
__all__ = ['plc_data', 'all_filters', 'add_filter']

import os
import glob
import numpy as np
import pickle
import shutil

from exotethys import ls_database

from pylightcurve.processes.files import open_dict, save_dict, open_dict_online
from .errors import *
from .__databases_setup__ import _setup_database

databases_file = '__databases__.pickle'
package_name = 'pylightcurve4'
github_link = 'https://www.dropbox.com/scl/fi/cmpsnfosniz108l97nugv/__databases__.pickle?rlkey=lv6oh6dqgs2blyvs5gwezsoyf&dl=1'


class PlcData:

    def __init__(self, _reset=False, _test=False):

        self.package_name = package_name
        self.version = '.'.join(open(os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                                  '__version__.txt')).read().split('.')[:2])

        self.build_in_databases_file_path = os.path.join(os.path.abspath(os.path.dirname(__file__)), databases_file)

        self.databases_directory_path = os.path.join(os.path.abspath(os.path.expanduser('~')),
                                                     '.{0}'.format(self.package_name))

        self.databases_file_path = os.path.join(self.databases_directory_path, databases_file)
        self.databases_file_path_new = os.path.join(self.databases_directory_path, databases_file + '_new')

        # initiate databases

        if not os.path.isdir(self.databases_directory_path):
            os.mkdir(self.databases_directory_path)

        if _test or not os.path.isfile(self.databases_file_path):
            shutil.copy(self.build_in_databases_file_path, self.databases_file_path)

        # check for updates in the databases (identified on github)

        test_online_db = open_dict_online(github_link)
        test_local_db = open_dict(self.databases_file_path)
        if test_online_db and test_online_db != test_local_db:
            save_dict(test_online_db, self.databases_file_path)

        # load databases

        self.databases = open_dict(self.databases_file_path)

        self.exotethys_loaded = _setup_database(self, 'exotethys')
        self.photometry_loaded = _setup_database(self, 'photometry')

    def exotethys(self):
        return self.exotethys_loaded

    def photometry(self):
        return self.photometry_loaded

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

    @staticmethod
    def reset_exotethys_data():
        for file in glob.glob(os.path.join(ls_database()[0], '*', '*')):
            try:
                _ = open_dict(file)
            except pickle.UnpicklingError:
                os.remove(file)


plc_data = PlcData()


def all_filters():

    xx = []
    for i in plc_data.all_filters():
        if 'hst_wfc3_g' not in i:
            xx.append(i)

    return xx


def add_filter(*args, **kwargs):
    plc_data.add_filter(*args, **kwargs)
