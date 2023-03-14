
__all__ = ['FixedTarget', 'Moon', 'Sun']

import numpy as np
from functools import lru_cache

from .angles import *
from .angles import _request_angle
from .moments import Moment, _request_time

from ..models.exoplanet_lc import convert_to_bjd_tdb, convert_to_jd_utc

from pylightcurve.errors import *
from pylightcurve.databases import plc_data


class _Target:

    def __init__(self, ra, dec, name=None, all_names=[], epoch=None):

        _request_angle(ra)
        _request_angle(dec)

        if 270 > dec.deg() > 90:
            raise PyLCInputError('Declination must be within -90, 90 degrees')

        self.ra, self.dec, self.name, self.all_names, self.epoch = ra, dec, name, all_names, epoch

    @lru_cache()
    def _get_class(self):
        return 'plc.{0}'.format(str(self.__class__).split('.')[-1][:-2])

    def reset(self, epoch):
        pass

    def __repr__(self):
        string = '{0}({1}'.format(self._get_class(), self.coord())
        if self.name:
            string += ', {0}'.format(self.name)
        return string + ')'

    @lru_cache()
    def coord(self):
        return '{0} {1}'.format(self.ra.hms(), self.dec.dms_coord())

    def distance_on_sphere(self, other):

        _request_target(other)

        return arccos(min(1, self.dec.sin() * other.dec.sin() +
                          self.dec.cos() * other.dec.cos() * (self.ra - other.ra).cos()))

    def convert_to_bjd_tdb(self, time_array, time_format):
        return convert_to_bjd_tdb(self.ra.deg(), self.dec.deg_coord(), time_array, time_format)

    def convert_to_jd_utc(self, time_array, time_format):
        return convert_to_jd_utc(self.ra.deg(), self.dec.deg_coord(), time_array, time_format)


class FixedTarget(_Target):

    def __init__(self, *args, **kwargs):

        _Target.__init__(self, *args, **kwargs)


class Sun(_Target):

    def __init__(self, epoch):

        _request_time(epoch)

        ssb_ra, ssb_dec, ssb_d, ssb_dt = plc_data.heliocentre(epoch.jd_utc())

        _Target.__init__(self, Rad(ssb_ra), Rad(ssb_dec), epoch=epoch)

    def reset(self, epoch):

        _request_time(epoch)

        if epoch.jd_utc() != self.epoch.jd_utc():

            ssb_ra, ssb_dec, ssb_d, ssb_dt = plc_data.heliocentre(epoch.jd_utc())

            _Target.__init__(self, Rad(ssb_ra), Rad(ssb_dec), epoch=epoch)

    def __repr__(self):
        return 'plc.Sun(RA(hms)/DEC(dms): {0} at {1} UTC)'.format(self.coord(), self.epoch.utc())


class Moon(_Target):

    def __init__(self, epoch):

        _request_time(epoch)

        ssb_ra, ssb_dec, ssb_d, ssb_dt = plc_data.mooncentre(epoch.jd_utc())

        _Target.__init__(self, Rad(ssb_ra), Rad(ssb_dec), epoch=epoch)

    def reset(self, epoch):

        _request_time(epoch)

        if epoch.jd_utc() != self.epoch.jd_utc():

            ssb_ra, ssb_dec, ssb_d, ssb_dt = plc_data.mooncentre(epoch.jd_utc())

            _Target.__init__(self, Rad(ssb_ra), Rad(ssb_dec), epoch=epoch)

    def illumination(self):

        sun_ra, sun_dec, sun_d, sun_dt = plc_data.heliocentre(self.epoch.jd_utc())
        moon_ra, moon_dec, moon_d, moon_dt = plc_data.mooncentre(self.epoch.jd_utc())

        theta = self.distance_on_sphere(FixedTarget(Rad(sun_ra), Rad(sun_dec)))

        sun_moon_d = np.sqrt(sun_d ** 2 + moon_d ** 2 - 2 * sun_d * moon_d * theta.cos())
        f = arcsin(theta.sin() * (sun_d / sun_moon_d))

        if theta.deg() < 90:
            return 0.5 * (1.0 - f.cos())
        else:
            return 0.5 * (1.0 + f.cos())

    def __repr__(self):
        return 'plc.Moon(RA(hms)/DEC(dms): {0} at {1} UTC)'.format(self.coord(), self.epoch.utc())


def _is_target(item):
    return isinstance(item, FixedTarget) or isinstance(item, Sun) or isinstance(item, Moon)


def _request_target(item):
    if _is_target(item):
        pass
    else:
        raise PyLCInputError('A plc.Target object is required')
