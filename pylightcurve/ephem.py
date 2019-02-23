from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from .oec import *


def ut_to_jd(ut_string):
    tt = astrotime(ut_string)
    return float(tt.jd)


def hjd(ra_target, dec_target, date, time_format):

    date = astrotime(date, format=time_format)
    ra_target *= np.pi / 180
    dec_target *= np.pi / 180

    julian_date = date.jd

    sun = astrosun(date)
    ra_sun, dec_sun = sun.ra.rad, sun.dec.rad

    a = sun.distance.m / 299792458.0
    b = np.sin(dec_target) * np.sin(dec_sun)
    c = np.cos(dec_target) * np.cos(dec_sun) * np.cos(ra_target - ra_sun)

    heliocentric_julian_date = julian_date - (a * (b + c)) / (24.0 * 60.0 * 60.0)

    return heliocentric_julian_date


def jd_to_hjd(ra_target, dec_target, julian_date):

    return hjd(ra_target, dec_target, julian_date, time_format='jd')


def mjd_to_hjd(ra_target, dec_target, modified_julian_date):

    return hjd(ra_target, dec_target, modified_julian_date, time_format='mjd')


def ut_to_hjd(ra_target, dec_target, ut_string):

    return hjd(ra_target, dec_target, ut_string, time_format='iso')


def ra_dec_string_to_deg(ra_dec_string):

    if len(ra_dec_string.split()) != 2:
        if len(ra_dec_string.split()) == 6:
            ra_dec_string = '{}:{}:{} {}:{}:{}'.format(*ra_dec_string.split())
        else:
            print('Wrong RA-DEC format.')
            return None, None

    if len(ra_dec_string.split()[0].split(':')) != 3:
        print('Wrong RA-DEC format.')
        return None, None

    if len(ra_dec_string.split()[0].split(':')[0]) != 2:
        print('Wrong RA-DEC format.')
        return None, None
    else:
        try:
            if int(ra_dec_string.split()[0].split(':')[0]) >= 24:
                print('Wrong RA-DEC format.')
                return None, None
        except ValueError:
            print('Wrong RA-DEC format.')
            return None, None

    if len(ra_dec_string.split()[0].split(':')[1]) != 2:
        print('Wrong RA-DEC format.')
        return None, None
    else:
        try:
            if int(ra_dec_string.split()[0].split(':')[1]) >= 60:
                print('Wrong RA-DEC format.')
                return None, None
        except ValueError:
            print('Wrong RA-DEC format.')
            return None, None

    if len(ra_dec_string.split()[0].split(':')[2]) < 2:
        print('Wrong RA-DEC format.')
        return None, None
    else:
        try:
            if float(ra_dec_string.split()[0].split(':')[2]) >= 60:
                print('Wrong RA-DEC format.')
                return None, None
        except ValueError:
            print('Wrong RA-DEC format.')
            return None, None

    if len(ra_dec_string.split()[1][1:].split(':')) != 3:
        print('Wrong RA-DEC format.')
        return None, None

    if len(ra_dec_string.split()[1].split(':')[0]) != 3:
        print('Wrong RA-DEC format.')
        return None, None
    else:
        try:
            if abs(int(ra_dec_string.split()[1].split(':')[0])) >= 90:
                print('Wrong RA-DEC format.')
                return None, None
        except ValueError:
            print('Wrong RA-DEC format.')
            return None, None

    if len(ra_dec_string.split()[1].split(':')[1]) != 2:
        print('Wrong RA-DEC format.')
        return None, None
    else:
        try:
            if int(ra_dec_string.split()[1].split(':')[1]) >= 60:
                print('Wrong RA-DEC format.')
                return None, None
        except ValueError:
            print('Wrong RA-DEC format.')
            return None, None

    if len(ra_dec_string.split()[1].split(':')[2]) < 2:
        print('Wrong RA-DEC format.')
        return None, None
    else:
        try:
            if float(ra_dec_string.split()[1].split(':')[2]) >= 60:
                print('Wrong RA-DEC format.')
                return None, None
        except ValueError:
            print('Wrong RA-DEC format.')
            return None, None

    try:
        ra_h, ra_m, ra_s = ra_dec_string.split()[0].split(':')
        ra = float(ra_h) / 1.0 + float(ra_m) / 60.0 + float(ra_s) / 3600.0
        ra *= 15.0
        dec_d, dec_m, dec_s = ra_dec_string.split()[1][1:].split(':')
        if ra_dec_string.split()[1][0] == '+':
            dec = float(dec_d) / 1.0 + float(dec_m) / 60.0 + float(dec_s) / 3600.0
        elif ra_dec_string.split()[1][0] == '-':
            dec = - float(dec_d) / 1.0 - float(dec_m) / 60.0 - float(dec_s) / 3600.0
        else:
            print('Wrong RA-DEC format.')
            return None, None

        return ra, dec

    except (ValueError, TypeError):
        print('Wrong RA-DEC format.')
        return None, None

#
# def find_next_transit(target, date, catalogue=None):
#
#     planet = find_target(target, catalogue)
#
#     (planet_name, stellar_logg, stellar_temperature, stellar_metallicity, rp_over_rs, fp_over_fs,
#      period, sma_over_rs, eccentricity, inclination, periastron, mid_time) = find_oec_parameters(target, catalogue)
#
#     date_hjd = jd_to_hjd(ephem.Date(date) + 2415020, planet.system.ra.deg, planet.system.dec.deg)
#
#     next_date_hjd = mid_time + (int((date_hjd - mid_time) / period) + 1) * period
#
#     next_date_dif = (next_date_hjd - date_hjd) * 24.0
#
#     next_date = ephem.Date(ephem.Date(date) + (next_date_hjd - date_hjd))
#
#     return planet_name, next_date_dif, '{0}/{1}/{2} {3}:{4}:{5:.0f}'.format(*next_date.tuple())
#
#
# def find_current_phase(target, julian_date, catalogue=None):
#
#     planet = find_target(target, catalogue)
#
#     (planet_name, stellar_logg, stellar_temperature, stellar_metallicity, rp_over_rs, fp_over_fs,
#      period, sma_over_rs, eccentricity, inclination, periastron, mid_time) = find_oec_parameters(target, catalogue)
#
#     if julian_date < mid_time:
#         mid_time -= int((mid_time - julian_date)/period + 10) * period
#
#     date_hjd = jd_to_hjd(julian_date, planet.system.ra.deg, planet.system.dec.deg)
#
#     return (date_hjd - mid_time) / period - int((date_hjd - mid_time) / period)
