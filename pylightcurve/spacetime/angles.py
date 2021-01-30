
__all__ = ['Degrees', 'Hours', 'Rad', 'arcsin', 'arccos', 'arctan', 'pi']


import numpy as np
from functools import lru_cache

from pylightcurve.errors import *


def _break_seconds(seconds):

    seconds = float(seconds)
    seconds_int = int(seconds)
    seconds_dec = int(str(seconds).split('.')[1])
    a = int(seconds_int / 3600.0)
    seconds_int = int(round(seconds_int - a * 3600.0))
    b = int(seconds_int / 60.0)
    seconds_int = int(round(seconds_int - b * 60.0))
    c = round(float('{0}.{1}'.format(seconds_int, seconds_dec)), 6)
    return a, b, c


def _collapse_to_seconds(degrees_hours, minutes, seconds):

        try:
            degrees_hours = float(degrees_hours)
        except:

            if isinstance(degrees_hours, str):

                if minutes != 0 or seconds != 0:
                    raise PyLCInputError('Not valid angle format')
                else:
                    try:
                        degrees_hours, minutes, seconds = degrees_hours.replace(':', ' ').split()
                        if degrees_hours[0] == '-':
                            minutes = '-{0}'.format(minutes)
                            seconds = '-{0}'.format(seconds)
                    except:
                        raise PyLCInputError('Not valid angle format')

        try:

            degrees_hours = float(degrees_hours)
            minutes = float(minutes)
            seconds = float(seconds)

        except:
            raise PyLCInputError('Not valid angle format')

        if degrees_hours <= 0 and minutes <= 0 and seconds <= 0:
            pass
        elif degrees_hours >= 0 and minutes >= 0 and seconds >= 0:
            pass
        else:
            raise PyLCInputError('Not valid angle format. Hours, degrees, minutes and seconds should be either ALL '
                                 'positive or ALL negative.')

        return round(degrees_hours * 3600.0 + minutes * 60.0 + seconds, 6)


class _DMS:

    def __init__(self, arcseconds):
        self.d, self.m, self.s = _break_seconds(arcseconds)
        self.list = [self.d, self.m, self.s]

    def __str__(self):
        return '{0}:{1}:{2}{3}'.format(str(self.d).zfill(2), str(self.m).zfill(2), '0'*(self.s < 10), str(self.s))

    def __repr__(self):
        return self.__str__()


class _HMS:

    def __init__(self, seconds):
        self.h, self.m, self.s = _break_seconds(seconds)

        self.list = [self.h, self.m, self.s]

    def __str__(self):
        return '{0}:{1}:{2}{3}'.format(str(self.h).zfill(2), str(self.m).zfill(2), '0'*(self.s < 10), str(self.s))

    def __repr__(self):
        return self.__str__()


class _Angle:

    def __init__(self, arcseconds):

        arcseconds = round(arcseconds, 6)

        if arcseconds < 0:
            arcseconds += 1296000.0 * (int(abs(arcseconds) / 1296000.0) + 1.0)
        arcseconds -= 1296000.0 * int(arcseconds / 1296000.0)

        self.arcseconds = arcseconds

        self._definition = '{0} arcsec'.format(arcseconds)

    @lru_cache()
    def deg(self):
        return self.arcseconds / 3600.0

    @lru_cache()
    def deg_coord(self):
        if self.deg() <= 90:
            return self.deg()

        elif self.deg() <= 180:
            return 180 - self.deg()

        elif self.deg() <= 270:
            return - (self.deg() - 180.0)

        else:
            return - (360.0 - self.deg())

    @lru_cache()
    def hours(self):
        return self.arcseconds / 15.0 / 3600.0

    @lru_cache()
    def rad(self):
        return (self.arcseconds / 3600.0) * np.pi / 180.0

    @lru_cache()
    def dms(self):
        return _DMS(self.arcseconds)

    @lru_cache()
    def dms_coord(self):

        if self.deg() <= 90:
            sign = '+'
            dec_print = self

        elif self.deg() <= 180:
            sign = '+'
            dec_print = pi - self

        elif self.deg() <= 270:
            sign = '-'
            dec_print = self - pi

        else:
            sign = '-'
            dec_print = 2 * pi - self

        return '{0}{1}'.format(sign, dec_print.dms())

    @lru_cache()
    def hms(self):
        return _HMS(self.arcseconds / 15.0)

    @lru_cache()
    def sin(self):
        return np.sin(self.rad())

    @lru_cache()
    def cos(self):
        return np.cos(self.rad())

    @lru_cache()
    def tan(self):
        return np.tan(self.rad())

    @lru_cache()
    def _get_class(self):
        return 'plc.{0}'.format(str(self.__class__).split('.')[-1][:-2])

    def __add__(self, other):
        _request_angle(other)
        return Degrees(0, 0, self.arcseconds + other.arcseconds)

    def __sub__(self, other):
        _request_angle(other)
        return Degrees(0, 0, self.arcseconds - other.arcseconds)

    def __mul__(self, other):
        if isinstance(other, int) or isinstance(other, float):
            return Degrees(0, 0, self.arcseconds * other)
        else:
            raise PyLCError('Operation not supported between {0} and {1}.'.format(self._get_class(), type(other)))

    def __truediv__(self, other):
        if isinstance(other, int) or isinstance(other, float):
            return Degrees(0, 0, self.arcseconds / other)
        else:
            raise PyLCError('Operation not supported between {0} and {1}.'.format(self._get_class(), type(other)))

    def __rmul__(self, other):
        return self.__mul__(other)

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return '{0}({1} DMS, defined as {2})'.format(self._get_class(), self.dms(), self._definition)


class Degrees(_Angle):

    def __init__(self, degrees, arcminutes=0.0, arcseconds=0.0):

        total_arcseconds = _collapse_to_seconds(degrees, arcminutes, arcseconds)

        _Angle.__init__(self, total_arcseconds)

        if arcminutes == 0 and arcseconds == 0:
            self._definition = '{0} degrees'.format(degrees)
        else:
            self._definition = '{0} degrees, {1} arcminutes, {2} arcseconds'.format(degrees, arcminutes, arcseconds)


class Hours(_Angle):

    def __init__(self, hours, minutes=0.0, seconds=0.0):

        total_arcseconds = _collapse_to_seconds(hours, minutes, seconds) * 15.0

        _Angle.__init__(self, total_arcseconds)

        if minutes == 0 and seconds == 0:
            self._definition = '{0} hours'.format(hours)
        else:
            self._definition = '{0} hours, {1} minutes, {2} seconds'.format(hours, minutes, seconds)


class Rad(_Angle):

    def __init__(self, rad):

        try:
            arcseconds = float(rad) * 648000.0 / np.pi
        except:
            raise PyLCInputError('Not valid input for rad.')

        _Angle.__init__(self, arcseconds)

        self._definition = '{0} rad'.format(rad)


def _is_angle(obsject):
    if isinstance(obsject, Degrees) or isinstance(obsject, Hours) or isinstance(obsject, Rad):
        return True
    else:
        return False


def _request_angle(item):
    if _is_angle(item):
        pass
    else:
        raise PyLCInputError('An angle object is required (plc.Degrees, plc.Hours or plc.Rad)')


def arccos(number):
    return Rad(np.arccos(number))


def arcsin(number):
    return Rad(np.arcsin(number))


def arctan(number):
    return Rad(np.arctan(number))


pi = Degrees(180.0)

