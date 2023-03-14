
__all__ = ['Degrees', 'Hours', 'Rad', 'arcsin', 'arccos', 'arctan']


import numpy as np
from functools import lru_cache

from pylightcurve.errors import *


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
                        raise PyLCInputError('Not valid angle format: ', degrees_hours, minutes, seconds)

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

        return degrees_hours * 3600.0 + minutes * 60.0 + seconds


class _Angle:

    def __init__(self, arcseconds):

        arcseconds = arcseconds

        if arcseconds < 0:
            arcseconds += 1296000.0 * (int(abs(arcseconds) / 1296000.0) + 1.0)
        arcseconds -= 1296000.0 * int(arcseconds / 1296000.0)

        self._arcseconds = arcseconds

    def _break_seconds(self, seconds):

        seconds = float(seconds)
        seconds_int = int(seconds)
        seconds_dec = int(str(seconds).split('.')[1])
        a = int(seconds_int / 3600.0)
        seconds_int = int(round(seconds_int - a * 3600.0))
        b = int(seconds_int / 60.0)
        seconds_int = int(round(seconds_int - b * 60.0))
        c = round(float('{0}.{1}'.format(seconds_int, seconds_dec)), 6)

        return a, b, c

    @lru_cache()
    def deg(self):
        """

        :returns:
        """
        return self._arcseconds / 3600.0

    @lru_cache()
    def deg_coord(self):
        """

        :returns:
        """
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
        """

        :returns:
        """
        return self._arcseconds / 15.0 / 3600.0

    @lru_cache()
    def rad(self):
        """

        :returns:
        """
        return (self._arcseconds / 3600.0) * np.pi / 180.0

    @lru_cache()
    def dms(self):
        """

        :returns:
        """
        d, m, s = self._break_seconds(self._arcseconds)
        return '{0}:{1}:{2}{3}'.format(str(d).zfill(2), str(m).zfill(2), '0'*(s < 10), str(s))

    @lru_cache()
    def dms_coord(self):
        """

        :returns:
        """

        if self.deg() <= 90:
            sign = '+'
            dec_print = self

        elif self.deg() <= 180:
            sign = '+'
            dec_print = Degrees(180.0) - self

        elif self.deg() <= 270:
            sign = '-'
            dec_print = self - Degrees(180.0)

        else:
            sign = '-'
            dec_print = 2 * Degrees(180.0) - self

        return '{0}{1}'.format(sign, dec_print.dms())

    @lru_cache()
    def hms(self):
        """

        :returns:
        """
        h, m, s = self._break_seconds(self._arcseconds / 15.0)
        return '{0}:{1}:{2}{3}'.format(str(h).zfill(2), str(m).zfill(2), '0'*(s < 10), str(s))

    @lru_cache()
    def sin(self):
        """

        :returns:
        """
        return np.sin(self.rad())

    @lru_cache()
    def cos(self):
        return np.cos(self.rad())

    @lru_cache()
    def tan(self):
        """

        :returns:
        """
        return np.tan(self.rad())

    def __add__(self, other):
        _request_angle(other)
        return Degrees(0, 0, self._arcseconds + other._arcseconds)

    def __sub__(self, other):
        _request_angle(other)
        return Degrees(0, 0, self._arcseconds - other._arcseconds)

    def __mul__(self, other):
        if isinstance(other, int) or isinstance(other, float):
            return Degrees(0, 0, self._arcseconds * other)
        else:
            raise PyLCError('Operation not supported between {0} and {1}.'.format(self._get_class(), type(other)))

    def __truediv__(self, other):
        if isinstance(other, int) or isinstance(other, float):
            return Degrees(0, 0, self._arcseconds / other)
        else:
            raise PyLCError('Operation not supported between {0} and {1}.'.format(self._get_class(), type(other)))

    def __rmul__(self, other):
        return self.__mul__(other)

    @lru_cache()
    def _get_class(self):
        return 'plc.{0}'.format(str(self.__class__).split('.')[-1][:-2])

    def __repr__(self):

        if self._get_class() in ['plc.Angle', 'plc.Degrees']:
            return '{0}({1} DMS)'.format(self._get_class(), self.dms())
        elif self._get_class() == 'plc.Hours':
            return '{0}({1} HMS)'.format(self._get_class(), self.hms())
        elif self._get_class() == 'plc.Rad':
            return '{0}({1} RAD)'.format(self._get_class(), self.rad())


class Degrees(_Angle):

    def __init__(self, degrees, arcminutes=0.0, arcseconds=0.0):
        """

        :param degrees:
        :param arcminutes:
        :param arcseconds:
        """
        total_arcseconds = _collapse_to_seconds(degrees, arcminutes, arcseconds)

        _Angle.__init__(self, total_arcseconds)


class Hours(_Angle):

    def __init__(self, hours, minutes=0.0, seconds=0.0):
        """

        :param hours:
        :param minutes:
        :param seconds:
        """
        total_arcseconds = _collapse_to_seconds(hours, minutes, seconds) * 15.0

        _Angle.__init__(self, total_arcseconds)


class Rad(_Angle):

    def __init__(self, rad):
        """

        :param rad:
        """
        try:
            arcseconds = float(rad) * 648000.0 / np.pi
        except:
            raise PyLCInputError('Not valid input for rad.')

        _Angle.__init__(self, arcseconds)


def _request_angle(item):
    if isinstance(item, Degrees) or isinstance(item, Hours) or isinstance(item, Rad):
        pass
    else:
        raise PyLCInputError('An angle object is required (plc.Degrees, plc.Hours or plc.Rad)')


def _reformat_or_request_angle(item):
    try:
        return Degrees(float(item))
    except:
        _request_angle(item)
        return item


def arccos(number):
    return Rad(np.arccos(number))


def arcsin(number):
    return Rad(np.arcsin(number))


def arctan(number):
    return Rad(np.arctan(number))


