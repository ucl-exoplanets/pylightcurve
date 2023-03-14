__all__ = ['Observatory']

import datetime

import numpy as np
from functools import lru_cache

from scipy.interpolate import interp1d

from .angles import *
from .angles import _request_angle
from .moments import *
from .moments import _request_time
from .targets import *
from .targets import _request_target

from pylightcurve.errors import *


class _Horizon:

    def __init__(self, horizon):

        self.horizon_value = None
        self.horizon_function = None

        try:
            self.horizon_value = float(horizon)
        except:
            pass

        if self.horizon_value is None:

            try:

                if isinstance(horizon, str):
                    horizon_list = []
                    for horizon_line in horizon.split('\n'):
                        if horizon_line != '':
                            horizon_list.append([Degrees(horizon_line.split()[0]), Degrees(horizon_line.split()[1])])

                elif isinstance(horizon, list):
                    horizon_list = [[Degrees(ff[0]), Degrees(ff[1])] for ff in horizon]

                elif isinstance(horizon, np.ndarray):
                    horizon_list = [[Degrees(ff[0]), Degrees(ff[1])] for ff in horizon]

                else:
                    raise PyLCInputError('Not valid horizon format')

                horizon_list = [[ff[0].deg(), ff[1].deg()] for ff in horizon_list]
                horizon_list.append([360.0, horizon_list[0][1]])

                horizon_list = np.swapaxes(np.array(horizon_list, dtype=float), 0, 1)
                self.horizon_function = interp1d(horizon_list[0], horizon_list[1])

            except:
                raise PyLCInputError('Not valid horizon format')

    def horizon(self, azimuth):

        _request_angle(azimuth)

        if self.horizon_value is None:
            return Degrees(self.horizon_function(azimuth.deg()))

        else:
            return Degrees(self.horizon_value)


class Observatory:

    def __init__(self, latitude, longitude, time_zone=0, horizon=0):

        _request_angle(latitude)
        _request_angle(longitude)

        if 270 > latitude.deg() > 90:
            raise PyLCInputError('Latitude must within -90, 90 degrees')

        self.latitude = latitude
        self.longitude = longitude

        try:
            self.time_zone = float(time_zone)
        except:
            raise PyLCInputError('Not valid time zone')

        if self.time_zone > 12 or self.time_zone < -12:
            raise PyLCInputError('Time zone must be within -12 and 12 hours')

        self.horizon = _Horizon(horizon).horizon

    @lru_cache()
    def coord(self):
        return '{0} {1}'.format(self.latitude.dms_coord(), self.longitude.dms())

    def lt(self, time_obj):
        _request_time(time_obj)
        return time_obj.utc() + datetime.timedelta(hours=self.time_zone)

    def lera(self, time_obj):
        _request_time(time_obj)
        return time_obj.era() + self.longitude

    def target_azimuth_altitude(self, target, time_obj):
        _request_time(time_obj)
        _request_target(target)
        target.reset(time_obj)

        ha = Degrees(self.lera(time_obj).deg() - target.ra.deg())

        altitude = arcsin(np.clip(target.dec.sin() * self.latitude.sin()
                                         + target.dec.cos() * self.latitude.cos() * ha.cos(), -1,
                                         1))

        if ha.hours() < 12:
            azimuth = Rad(np.pi - np.arccos(np.clip((target.dec.sin() -
                                                     altitude.sin() * self.latitude.sin()) /
                                                    (altitude.cos() * self.latitude.cos()), -1, 1)))
        else:

            azimuth = Rad(np.pi + np.arccos(np.clip((target.dec.sin() -
                                                     altitude.sin() * self.latitude.sin()) /
                                                    (altitude.cos() * self.latitude.cos()), -1, 1)))

        return azimuth, altitude

    def airmass(self, target, time_obj):
        _request_time(time_obj)
        _request_target(target)

        azimuth, altitude = self.target_azimuth_altitude(target, time_obj)

        alt = altitude.deg_coord()

        return 1.0 / (np.cos((90 - alt) * np.pi / 180) + 0.50572 * ((6.07995 + alt) ** (-1.6364)))

    def is_target_visible(self, target, time_obj):
        _request_time(time_obj)
        _request_target(target)

        azimuth, altitude = self.target_azimuth_altitude(target, time_obj)
        horizon_altitude = self.horizon(azimuth)
        return altitude.deg_coord() > horizon_altitude.deg_coord()

    def target_horizon_crossings(self, target, time_obj_1, window):

        _request_time(time_obj_1)
        _request_target(target)

        def _target_over_horizon(m_target, time_obj):
            m_target.reset(time_obj)
            azimuth, altitude = self.target_azimuth_altitude(m_target, time_obj)
            horizon_altitude = self.horizon(azimuth)
            return altitude.deg_coord() - horizon_altitude.deg_coord()

        crossings = []

        jd1 = time_obj_1.jd_utc()
        jd2 = jd1 + window

        test_alt = _target_over_horizon(target, Moment(jd_utc=jd1))
        test_jd = jd1

        for jd in np.arange(jd1 + 0.007, jd2 + 0.0001, 0.007):

            new_test_alt = _target_over_horizon(target, Moment(jd_utc=jd))
            new_test_jd = jd

            if test_alt * new_test_alt < 0:

                if test_alt > 0:
                    crossings.append(
                        [Moment(jd_utc=test_jd + (0.007 * abs(test_alt / new_test_alt)) / (1 + abs(test_alt / new_test_alt))), 'set']
                    )
                else:
                    crossings.append(
                        [Moment(jd_utc=test_jd + (0.007 * abs(test_alt / new_test_alt)) / (1 + abs(test_alt / new_test_alt))), 'rise']
                    )

            else:
                pass

            test_alt = new_test_alt
            test_jd = new_test_jd

        return crossings

    def target_altitude_crossings(self, target, start, window, limit_altitude):

        _request_time(start)
        _request_target(target)
        _request_angle(limit_altitude)

        def _target_over_altitude(m_target, time_obj):
            m_target.reset(time_obj)
            azimuth, altitude = self.target_azimuth_altitude(m_target, time_obj)
            return altitude.deg_coord() - limit_altitude.deg_coord()

        crossings = []

        jd1 = start.jd_utc()
        jd2 = jd1 + window

        test_alt = _target_over_altitude(target, Moment(jd_utc=jd1))
        test_jd = jd1

        for jd in np.arange(jd1 + 0.007, jd2 + 0.0001, 0.007):

            new_test_alt = _target_over_altitude(target, Moment(jd_utc=jd))
            new_test_jd = jd

            if test_alt * new_test_alt < 0:

                if test_alt > 0:
                    crossings.append(
                        [Moment(jd_utc=test_jd + (0.007 * abs(test_alt / new_test_alt)) / (1 + abs(test_alt / new_test_alt))), 'set']
                    )
                else:
                    crossings.append(
                        [Moment(jd_utc=test_jd + (0.007 * abs(test_alt / new_test_alt)) / (1 + abs(test_alt / new_test_alt))), 'rise']
                    )

            else:
                pass

            test_alt = new_test_alt
            test_jd = new_test_jd

        return crossings

    def periodic_events_visibility(self, target, start, window, ref_epoch, ref_epoch_format, period, duration,
                                   target_min_altitude=Degrees(20), sun_max_altitude=Degrees(-18),
                                   max_moon_illumination=0.9, min_moon_distance=Degrees(30)):

        _request_time(start)
        _request_target(target)
        _request_angle(target_min_altitude)
        _request_angle(sun_max_altitude)
        _request_angle(min_moon_distance)

        t0 = target.convert_to_bjd_tdb(ref_epoch, time_format=ref_epoch_format)
        t1 = target.convert_to_bjd_tdb(start.jd_utc(), time_format='JD_UTC')
        t2 = t1 + window

        sun = Sun(now())

        events = {}

        t = (t0 - duration / 2.0) + int(1 + (t1 - (t0 - duration / 2.0)) / period) * period

        while t + duration < t2:
            e1 = Moment(jd_utc=target.convert_to_jd_utc(t, 'BJD_TDB'))
            if self.is_target_visible(target, e1):
                taz1, ta1 = self.target_azimuth_altitude(target, e1)
                if not target_min_altitude or ta1.deg_coord() > target_min_altitude.deg_coord():
                    sun.reset(e1)
                    sa1 = self.target_azimuth_altitude(sun, e1)[1]
                    if not sun_max_altitude or sa1.deg_coord() < sun_max_altitude.deg_coord():
                        e2 = e1 + datetime.timedelta(days=duration / 2.0)
                        if self.is_target_visible(target, e2):
                            taz2, ta2 = self.target_azimuth_altitude(target, e2)
                            if not target_min_altitude or ta2.deg_coord() > target_min_altitude.deg_coord():
                                sun.reset(e2)
                                sa2 = self.target_azimuth_altitude(sun, e2)[1]
                                if not sun_max_altitude or sa2.deg_coord() < sun_max_altitude.deg_coord():
                                    e3 = e1 + datetime.timedelta(days=duration)
                                    if self.is_target_visible(target, e3):
                                        taz3, ta3 = self.target_azimuth_altitude(target, e3)
                                        if not target_min_altitude or ta3.deg_coord() > target_min_altitude.deg_coord():
                                            sun.reset(e3)
                                            sa3 = self.target_azimuth_altitude(sun, e3)[1]
                                            if not sun_max_altitude or sa3.deg_coord() < sun_max_altitude.deg_coord():
                                                moon = Moon(e2)
                                                moon_illumination = moon.illumination()
                                                if not max_moon_illumination or moon_illumination < max_moon_illumination:
                                                    moon_distance = moon.distance_on_sphere(target)
                                                    if not min_moon_distance or moon_distance.deg() > min_moon_distance.deg():
                                                        events[len(list(events.keys())) + 1] = {
                                                            'event_start': {'time': e1,
                                                                            'target_altitude': ta1,
                                                                            'target_azimuth': taz1,
                                                                            'sun_altitude': sa1},
                                                            'event_middle': {'time': e2,
                                                                             'target_altitude': ta2,
                                                                             'target_azimuth': taz2,
                                                                             'sun_altitude': sa2},
                                                            'event_end': {'time': e3,
                                                                          'target_altitude': ta3,
                                                                          'target_azimuth': taz3,
                                                                          'sun_altitude': sa3},
                                                            'moon': {'illumination': moon_illumination,
                                                                     'distance': moon_distance.deg()}
                                                            }

            t = t + period

        return events

    def __repr__(self):
        return 'plc.Observatory(latitude(dms)/longitude(dms): {0}, time zone: {1} hours)'.format(
            self.coord(), self.time_zone)


def _is_observatory(item):
    return isinstance(item, Observatory)


def _request_observatory(item):
    if _is_observatory(item):
        pass
    else:
        raise PyLCInputError('A plc.Observatory object is required.')
