from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from .catalogues_oec import *
from .tools_files import *


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


def from_hjd(ra_target, dec_target, heliocentric_julian_date, time_format):

    def func(x, jd):
        return jd_to_hjd(ra_target, dec_target, jd) - heliocentric_julian_date

    popt, pcov = curve_fit(func, [0], [0], p0=[heliocentric_julian_date])
    julian_date = astrotime(popt[0], format='jd')

    if time_format == 'jd':
        return julian_date.jd
    elif time_format == 'mjd':
        return julian_date.mjd
    elif time_format == 'iso':
        return julian_date.iso


def hjd_to_jd(ra_target, dec_target, heliocentric_julian_date):

    return from_hjd(ra_target, dec_target, heliocentric_julian_date, time_format='jd')


def hjd_to_mjd(ra_target, dec_target, heliocentric_julian_date):

    return from_hjd(ra_target, dec_target, heliocentric_julian_date, time_format='mjd')


def hjd_to_ut(ra_target, dec_target, heliocentric_julian_date):

    return from_hjd(ra_target, dec_target, heliocentric_julian_date, time_format='iso')


bjd_dict = open_dict(glob.glob(os.path.join(databases.ephemeris, 'bjd_dict.pickle''*'))[0])
hjd_dict = open_dict(glob.glob(os.path.join(databases.ephemeris, 'hjd_dict.pickle''*'))[0])


def bjd_tdb(ra_target, dec_target, date, time_format):

    date = astrotime(date, format=time_format)
    ra_target *= np.pi / 180
    dec_target *= np.pi / 180

    julian_date = date.jd

    ssb_t = bjd_dict[int(julian_date)]['t']
    ssb_ra = bjd_dict[int(julian_date)]['ra']
    ssb_dec = bjd_dict[int(julian_date)]['dec']
    ssb_d = bjd_dict[int(julian_date)]['d']
    ssb_dt = bjd_dict[int(julian_date)]['dt']

    ssb_ra = interp1d(ssb_t, ssb_ra, kind='cubic')(julian_date)
    ssb_dec = interp1d(ssb_t, ssb_dec, kind='cubic')(julian_date)
    ssb_d = interp1d(ssb_t, ssb_d, kind='cubic')(julian_date)
    ssb_dt = interp1d(ssb_t, ssb_dt, kind='cubic')(julian_date)

    a = ssb_d / 60 / 24
    b = np.sin(dec_target) * np.sin(ssb_dec)
    c = np.cos(dec_target) * np.cos(ssb_dec) * np.cos(ra_target - ssb_ra)

    return julian_date - a * (b + c) + ssb_dt / 60 / 60 / 24


def jd_utc_to_hjd_tdb(ra_target, dec_target, jd_utc):

    ra_target *= np.pi / 180
    dec_target *= np.pi / 180

    sun_t = hjd_dict[int(jd_utc)]['t']
    sun_ra = hjd_dict[int(jd_utc)]['ra']
    sun_dec = hjd_dict[int(jd_utc)]['dec']
    sun_d = hjd_dict[int(jd_utc)]['d']
    sun_dt = hjd_dict[int(jd_utc)]['dt']

    sun_ra = interp1d(sun_t, sun_ra, kind='cubic')(jd_utc)
    sun_dec = interp1d(sun_t, sun_dec, kind='cubic')(jd_utc)
    sun_d = interp1d(sun_t, sun_d, kind='cubic')(jd_utc)
    sun_dt = interp1d(sun_t, sun_dt, kind='cubic')(jd_utc)

    a = sun_d / 60 / 24
    b = np.sin(dec_target) * np.sin(sun_dec)
    c = np.cos(dec_target) * np.cos(sun_dec) * np.cos(ra_target - sun_ra)

    return jd_utc - a * (b + c) + sun_dt / 60 / 60 / 24


def jd_utc_to_hjd_utc(ra_target, dec_target, jd_utc):

    ra_target *= np.pi / 180
    dec_target *= np.pi / 180

    sun_t = hjd_dict[int(jd_utc)]['t']
    sun_ra = hjd_dict[int(jd_utc)]['ra']
    sun_dec = hjd_dict[int(jd_utc)]['dec']
    sun_d = hjd_dict[int(jd_utc)]['d']

    sun_ra = interp1d(sun_t, sun_ra, kind='cubic')(jd_utc)
    sun_dec = interp1d(sun_t, sun_dec, kind='cubic')(jd_utc)
    sun_d = interp1d(sun_t, sun_d, kind='cubic')(jd_utc)

    a = sun_d / 60 / 24
    b = np.sin(dec_target) * np.sin(sun_dec)
    c = np.cos(dec_target) * np.cos(sun_dec) * np.cos(ra_target - sun_ra)

    return jd_utc - a * (b + c)


def jd_utc_to_bjd_tdb(ra_target, dec_target, jd_utc):

    return bjd_tdb(ra_target, dec_target, jd_utc, time_format='jd')


def mjd_utc_to_bjd_tdb(ra_target, dec_target, mjd_utc):

    return bjd_tdb(ra_target, dec_target, mjd_utc, time_format='mjd')


def utc_to_bjd_tdb(ra_target, dec_target, utc):

    return bjd_tdb(ra_target, dec_target, utc, time_format='iso')


def hjd_utc_to_bjd_tdb(ra_target, dec_target, hjd_utc):

    def func(x, jd):
        return jd_utc_to_hjd_utc(ra_target, dec_target, jd) - hjd_utc

    popt, pcov = curve_fit(func, [0], [0], p0=[hjd_utc])

    return jd_utc_to_bjd_tdb(ra_target, dec_target, popt[0])


def hjd_tdb_to_bjd_tdb(ra_target, dec_target, hjd_tdb):

    def func(x, jd):
        return jd_utc_to_hjd_tdb(ra_target, dec_target, jd) - hjd_tdb

    popt, pcov = curve_fit(func, [0], [0], p0=[hjd_tdb])

    return jd_utc_to_bjd_tdb(ra_target, dec_target, popt[0])


def bjd_utc_to_bjd_tdb(ra_target, dec_target, bjd_utc):

    ssb_t = bjd_dict[int(bjd_utc)]['t']
    ssb_dt = bjd_dict[int(bjd_utc)]['dt']
    ssb_dt = interp1d(ssb_t, ssb_dt, kind='cubic')(bjd_utc)

    return bjd_utc + ssb_dt / 60 / 60 / 24


def bjd_tdb_to_jd_utc(ra_target, dec_target, bjd_tdb_in):

    def func(x, jd):
        return jd_utc_to_bjd_tdb(ra_target, dec_target, jd) - bjd_tdb_in

    popt, pcov = curve_fit(func, [0], [0], p0=[bjd_tdb_in])

    return popt[0]


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


def lat_long_string_to_deg(lat_long_string):

    if len(lat_long_string.split()) != 2:
        if len(lat_long_string.split()) == 6:
            lat_long_string = '{}:{}:{} {}:{}:{}'.format(*lat_long_string.split())
        else:
            print('Wrong RA-DEC format.')
            return None, None

    if len(lat_long_string.split()[0].split(':')) != 3:
        print('Wrong RA-DEC format.')
        return None, None

    if len(lat_long_string.split()[0].split(':')[0]) != 2:
        print('Wrong RA-DEC format.')
        return None, None
    else:
        try:
            if int(lat_long_string.split()[0].split(':')[0]) >= 24:
                print('Wrong RA-DEC format.')
                return None, None
        except ValueError:
            print('Wrong RA-DEC format.')
            return None, None

    if len(lat_long_string.split()[0].split(':')[1]) != 2:
        print('Wrong RA-DEC format.')
        return None, None
    else:
        try:
            if int(lat_long_string.split()[0].split(':')[1]) >= 60:
                print('Wrong RA-DEC format.')
                return None, None
        except ValueError:
            print('Wrong RA-DEC format.')
            return None, None

    if len(lat_long_string.split()[0].split(':')[2]) < 2:
        print('Wrong RA-DEC format.')
        return None, None
    else:
        try:
            if float(lat_long_string.split()[0].split(':')[2]) >= 60:
                print('Wrong RA-DEC format.')
                return None, None
        except ValueError:
            print('Wrong RA-DEC format.')
            return None, None

    if len(lat_long_string.split()[1][1:].split(':')) != 3:
        print('Wrong RA-DEC format.')
        return None, None

    if len(lat_long_string.split()[1].split(':')[0]) != 3:
        print('Wrong RA-DEC format.')
        return None, None
    else:
        try:
            if abs(int(lat_long_string.split()[1].split(':')[0])) >= 90:
                print('Wrong RA-DEC format.')
                return None, None
        except ValueError:
            print('Wrong RA-DEC format.')
            return None, None

    if len(lat_long_string.split()[1].split(':')[1]) != 2:
        print('Wrong RA-DEC format.')
        return None, None
    else:
        try:
            if int(lat_long_string.split()[1].split(':')[1]) >= 60:
                print('Wrong RA-DEC format.')
                return None, None
        except ValueError:
            print('Wrong RA-DEC format.')
            return None, None

    if len(lat_long_string.split()[1].split(':')[2]) < 2:
        print('Wrong RA-DEC format.')
        return None, None
    else:
        try:
            if float(lat_long_string.split()[1].split(':')[2]) >= 60:
                print('Wrong RA-DEC format.')
                return None, None
        except ValueError:
            print('Wrong RA-DEC format.')
            return None, None

    try:
        lat_d, lat_m, lat_s = lat_long_string.split()[0].split(':')
        lat = float(lat_d) / 1.0 + float(lat_m) / 60.0 + float(lat_s) / 3600.0
        long_d, long_m, long_s = lat_long_string.split()[1][1:].split(':')
        if lat_long_string.split()[1][0] == '+':
            long = float(long_d) / 1.0 + float(long_m) / 60.0 + float(long_s) / 3600.0
        elif lat_long_string.split()[1][0] == '-':
            long = - float(long_d) / 1.0 - float(long_m) / 60.0 - float(long_s) / 3600.0
        else:
            print('Wrong RA-DEC format.')
            return None, None

        return lat, long

    except (ValueError, TypeError):
        print('Wrong LAT-LONG format.')
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