from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

__all__ = ['oec_catalogue', 'find_oec_parameters', 'find_oec_stellar_parameters', 'find_next_transit',
           'find_current_phase']

import ephem
import numpy as np
import exodata
import seaborn as sns
import warnings
import exodata.astroquantities as aq

from .database_handling import *

sns.reset_orig()
warnings.filterwarnings('ignore', '\'second\' was found  to be \'60.0\', which is not in range [0,60). '
                                  'Treating as 0 sec, +1 min [astropy.coordinates.angle_utilities]')


def oec_catalogue():

        return exodata.OECDatabase(gzip.GzipFile(oec_database()), stream=True)


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

    # known mistakes in the catalogue

    if name == 'WASP-121 b':
        mid_time = 2456635.70832
    elif name == 'Qatar-1 b':
        mid_time = 2455518.4102

    return (name, stellar_logg, stellar_temperature, stellar_metallicity, rp_over_rs, fp_over_fs, period, sma_over_rs,
            eccentricity, inclination, periastron, mid_time)


def find_oec_stellar_parameters(target, catalogue=None):

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

    stellar_radius = float(planet.star.R)

    stellar_vmag = float(planet.star.vmag)

    return name, stellar_logg, stellar_temperature, stellar_metallicity, stellar_radius, stellar_vmag


def jd_to_hjd(julian_date, ra_target, dec_target):

    # calculate the RA and DEC of the sun for this time,
    # we need "- 2415020" to convert julian date to dublin julian date (used by ephem)

    sun = ephem.Sun()
    sun.compute(ephem.date(julian_date - 2415020))
    ra_sun, dec_sun = float(sun.ra), float(sun.dec)

    # get the RA and DEC of the target and convert degrees to radians

    ra_target *= np.pi / 180
    dec_target *= np.pi / 180

    # calculate the hjd correction (in days) for the given time and target

    hjd_correction = - ((149597870700.0 / ephem.c) *
                        (np.sin(dec_target) * np.sin(dec_sun) +
                         np.cos(dec_target) * np.cos(dec_sun) * np.cos(ra_target - ra_sun)) *
                        (1.0 / (24.0 * 60.0 * 60.0)))

    # return the heliocentric julian date

    return julian_date + hjd_correction


def find_next_transit(target, date, catalogue=None):

    if catalogue is None:
        catalogue = oec_catalogue()

    planet = catalogue.searchPlanet(target)

    (planet_name, stellar_logg, stellar_temperature, stellar_metallicity, rp_over_rs, fp_over_fs,
     period, sma_over_rs, eccentricity, inclination, periastron, mid_time) = find_oec_parameters(target, catalogue)

    date_hjd = jd_to_hjd(ephem.Date(date) + 2415020, planet.system.ra.deg, planet.system.dec.deg)

    next_date_hjd = mid_time + (int((date_hjd - mid_time) / period) + 1) * period

    next_date_dif = (next_date_hjd - date_hjd) * 24.0

    next_date = ephem.Date(ephem.Date(date) + (next_date_hjd - date_hjd))

    return planet_name, next_date_dif, '{0}/{1}/{2} {3}:{4}:{5:.0f}'.format(*next_date.tuple())


def find_current_phase(target, date, catalogue=None):

    if catalogue is None:
        catalogue = oec_catalogue()

    planet = catalogue.searchPlanet(target)

    (planet_name, stellar_logg, stellar_temperature, stellar_metallicity, rp_over_rs, fp_over_fs,
     period, sma_over_rs, eccentricity, inclination, periastron, mid_time) = find_oec_parameters(target, catalogue)

    date_hjd = jd_to_hjd(ephem.Date(date) + 2415020, planet.system.ra.deg, planet.system.dec.deg)

    return (date_hjd - mid_time) / period - int((date_hjd - mid_time) / period)