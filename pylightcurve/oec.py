__all__ = ['oec_catalogue', 'find_oec_parameters']


import os
import time
import gzip
import urllib
import socket

import numpy as np

import thirdparty_exodata
import thirdparty_exodata.astroquantities as aq


def oec_catalogue():

    data_base_location = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'oec_data_base')

    data_base_url = 'https://github.com/OpenExoplanetCatalogue/oec_gzip/raw/master/systems.xml.gz'

    data_base_file_path = os.path.join(data_base_location, 'systems.xml')
    last_update_file_path = os.path.join(data_base_location, 'systems_last_update.txt')

    date = time.strftime('%y%m%d')
    update = False
    if not os.path.isfile(last_update_file_path):
        update = True
    elif not os.path.isfile(data_base_file_path):
        update = True
    elif int(open(last_update_file_path).readlines()[0]) < int(date):
        update = True

    if update:

        print 'Updating OEC...'

        try:
            socket.setdefaulttimeout(5)
            urllib.urlretrieve(data_base_url, data_base_file_path + '.gz')
            socket.setdefaulttimeout(30)

            w = open(data_base_file_path, 'w')
            for i in gzip.open(data_base_file_path + '.gz'):
                w.write(i)

            w.close()

            os.remove('{0}.gz'.format(data_base_file_path))

            w = open(last_update_file_path, 'w')
            w.write(date)
            w.close()

        except IOError:
            print 'Updating OEC failed.'
            pass

    return thirdparty_exodata.OECDatabase(data_base_file_path, stream=True)


def find_oec_parameters(target, catalogue=None):

    if catalogue is None:
        catalogue = oec_catalogue()

    planet = catalogue.searchPlanet(target)

    if isinstance(planet, list):
        planet = planet[0]

    name = planet.name

    stellar_logg = float(planet.star.calcLogg())

    stellar_temperature = float(planet.star.T)

    if np.isnan(planet.star.Z):
        stellar_metallicity = 0
    else:
        stellar_metallicity = planet.star.Z

    rp_over_rs = float(planet.R.rescale(aq.m) / planet.star.R.rescale(aq.m))

    planet_temperature = float(planet.calcTemperature())
    fp_over_fs = (rp_over_rs ** 2.0) * \
                 ((np.exp(6.626 * (10 ** (-5)) * 299792458.0 / (1.4 * 1.381 * stellar_temperature)) - 1.)
                  / (np.exp(6.626 * (10 ** (-5)) * 299792458.0 / (1.4 * 1.381 * planet_temperature)) - 1.))

    period = float(planet.P)

    sma_over_rs = float(planet.calcSMA().rescale(aq.m) / planet.star.R.rescale(aq.m))

    if np.isnan(planet.e):
        eccentricity = 0.0
    else:
        eccentricity = float(planet.e)

    if np.isnan(planet.i):
        inclination = 90.0
    else:
        inclination = float(planet.i)

    if np.isnan(planet.periastron):
        periastron = 0.0
    else:
        periastron = float(planet.periastron)

    mid_time = float(planet.transittime)

    return (name, stellar_logg, stellar_temperature, stellar_metallicity, rp_over_rs, fp_over_fs,
            period, sma_over_rs, eccentricity, inclination, periastron, mid_time)