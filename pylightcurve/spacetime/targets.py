
__all__ = ['FixedTarget']

import numpy as np

from pylightcurve.errors import *
from pylightcurve.__databases__ import plc_data
from pylightcurve.analysis.curve_fit import curve_fit
from pylightcurve.spacetime.angles import arccos, _request_angle


class _Target:

    def __init__(self, ra, dec):

        _request_angle(ra)
        _request_angle(dec)

        if 270 > dec.deg() > 90:
            raise PyLCInputError('Declination must within -90, 90 degrees')

        self.dec = dec
        self.ra = ra

        self.coord = '{0} {1}'.format(self.ra.hms(), self.dec.dms_coord())


class FixedTarget(_Target):

    def __init__(self, ra, dec):

        _Target.__init__(self, ra, dec)

    def convert_to_bjd_tdb(self, time, time_format):

        if isinstance(time, float):

            if time_format in ['BJD_TDB', 'BJD_TT']:
                return time
            elif time_format == 'JD_UTC':
                return self._bjd_tdb(None, time)
            elif time_format == 'MJD_UTC':
                return self._bjd_tdb(None, time + 2400000.5)
            elif time_format == 'BJD_UTC':
                return self._bjd_tdb(None, curve_fit(self._bjd_utc, [0], [time], p0=[time])[0][0])
            elif time_format in ['HJD_TDB', 'HJD_TT']:
                return self._bjd_tdb(None, curve_fit(self._hjd_tdb, [0], [time], p0=[time])[0][0])
            elif time_format == 'HJD_UTC':
                return self._bjd_tdb(None, curve_fit(self._hjd_utc, [0], [time], p0=[time])[0][0])
            else:
                raise PyLCInputError(
                    'Not valid time format. Available formats: JD_UTC, MJD_UTC, BJD_UTC, BJD_TDB, BJD_TT, '
                    'HJD_UTC, HJD_BJD, HJD_TT')

        else:
            try:
                time = np.array(time, dtype=float)

                if time_format in ['BJD_TDB', 'BJD_TT']:
                    return time
                elif time_format == 'JD_UTC':
                    return np.array([self._bjd_tdb(None, ff) for ff in time])
                elif time_format == 'MJD_UTC':
                    return np.array([self._bjd_tdb(None, ff + 2400000.5) for ff in time])
                elif time_format == 'BJD_UTC':
                    return np.array([self._bjd_tdb(None, curve_fit(self._bjd_utc, [0], [ff], p0=[ff])[0][0])
                                     for ff in time])
                elif time_format in ['HJD_TDB', 'HJD_TT']:
                    return np.array([self._bjd_tdb(None, curve_fit(self._hjd_tdb, [0], [ff], p0=[ff])[0][0])
                                     for ff in time])
                elif time_format == 'HJD_UTC':
                    return np.array([self._bjd_tdb(None, curve_fit(self._hjd_utc, [0], [ff], p0=[ff])[0][0])
                                     for ff in time])
                else:
                    raise PyLCInputError(
                        'Not valid time format. Available formats: JD_UTC, MJD_UTC, BJD_UTC, BJD_TDB, BJD_TT, '
                        'HJD_UTC, HJD_BJD, HJD_TT')

            except:
                raise PyLCInputError('Not valid input for time')

    def convert_to_jd(self, time, time_format):

        if isinstance(time, float):

            if time_format in ['BJD_TDB', 'BJD_TT']:
                return curve_fit(self._bjd_tdb, [0], [time], p0=[time])[0][0]
            elif time_format == 'JD_UTC':
                return time
            elif time_format == 'MJD_UTC':
                return time + 2400000.5
            elif time_format == 'BJD_UTC':
                return curve_fit(self._bjd_utc, [0], [time], p0=[time])[0][0]
            elif time_format in ['HJD_TDB', 'HJD_TT']:
                return curve_fit(self._hjd_tdb, [0], [time], p0=[time])[0][0]
            elif time_format == 'HJD_UTC':
                return curve_fit(self._hjd_utc, [0], [time], p0=[time])[0][0]
            else:
                raise PyLCInputError(
                    'Not valid time format. Available formats: JD_UTC, MJD_UTC, BJD_UTC, BJD_TDB, BJD_TT, '
                    'HJD_UTC, HJD_BJD, HJD_TT')

        else:
            try:
                time = np.array(time, dtype=float)

                if time_format in ['BJD_TDB', 'BJD_TT']:
                    return np.array([curve_fit(self._bjd_tdb, [0], [ff], p0=[time])[0][0] for ff in time])
                elif time_format == 'JD_UTC':
                    return time
                elif time_format == 'MJD_UTC':
                    return time + 2400000.5
                elif time_format == 'BJD_UTC':
                    return np.array([curve_fit(self._bjd_utc, [0], [ff], p0=[time])[0][0] for ff in time])
                elif time_format in ['HJD_TDB', 'HJD_TT']:
                    return np.array([curve_fit(self._hjd_tdb, [0], [ff], p0=[time])[0][0] for ff in time])
                elif time_format == 'HJD_UTC':
                    return np.array([curve_fit(self._hjd_utc, [0], [ff], p0=[time])[0][0] for ff in time])
                else:
                    raise PyLCInputError(
                        'Not valid time format. Available formats: JD_UTC, MJD_UTC, BJD_UTC, BJD_TDB, BJD_TT, '
                        'HJD_UTC, HJD_BJD, HJD_TT')

            except:
                raise PyLCInputError('Not valid input for time')

    def convert_to_mjd(self, time, time_format):

        return self.convert_to_jd(time, time_format) - 2400000.5

    def _hjd_utc(self, x, jd):

        ssb_ra, ssb_dec, ssb_d, ssb_dt = plc_data.heliocentre(jd)

        a = ssb_d / 60.0 / 24.0
        b = self.dec.sin() * np.sin(ssb_dec)
        c = self.dec.cos() * np.cos(ssb_dec) * np.cos(self.ra.rad() - ssb_ra)

        return jd - a * (b + c)

    def _hjd_tdb(self, x, jd):

        ssb_ra, ssb_dec, ssb_d, ssb_dt = plc_data.heliocentre(jd)

        a = ssb_d / 60.0 / 24.0
        b = self.dec.sin() * np.sin(ssb_dec)
        c = self.dec.cos() * np.cos(ssb_dec) * np.cos(self.ra.rad() - ssb_ra)

        return jd - a * (b + c) + ssb_dt / 60.0 / 60.0 / 24.0

    def _bjd_utc(self, x, jd):

        ssb_ra, ssb_dec, ssb_d, ssb_dt = plc_data.barycentre(jd)

        a = ssb_d / 60.0 / 24.0
        b = self.dec.sin() * np.sin(ssb_dec)
        c = self.dec.cos() * np.cos(ssb_dec) * np.cos(self.ra.rad() - ssb_ra)

        return jd - a * (b + c)

    def _bjd_tdb(self, x, jd):

        ssb_ra, ssb_dec, ssb_d, ssb_dt = plc_data.barycentre(jd)

        a = ssb_d / 60.0 / 24.0
        b = self.dec.sin() * np.sin(ssb_dec)
        c = self.dec.cos() * np.cos(ssb_dec) * np.cos(self.ra.rad() - ssb_ra)

        return jd - a * (b + c) + ssb_dt / 60.0 / 60.0 / 24.0

    def distance_on_sphere(self, other):

        _request_target(other)

        return arccos(self.dec.sin() * other.dec.sin() + self.dec.cos() * other.dec.cos() * (self.ra - other.ra).cos())

    def __str__(self):
        return 'plc.FixedTarget(RA(hms)/DEC(dms): {0})'.format(self.coord)

    def __repr__(self):
        return self.__str__()


def _is_target(item):
    if isinstance(item, FixedTarget):
        return True
    else:
        return False


def _request_target(item):
    if _is_target(item):
        pass
    else:
        raise PyLCInputError('A plc.Target object is required (plc.FixedTarget)')
