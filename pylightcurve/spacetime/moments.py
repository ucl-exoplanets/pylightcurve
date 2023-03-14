
__all__ = ['Moment', 'now']

import datetime

from functools import lru_cache

from .angles import *

from pylightcurve.errors import *
from pylightcurve.databases import plc_data


_date_format_error = ('Not valid date format, '
                      'acceptable formats are (seconds are rounded to the 6th digit):'
                      '\nYYYY-MM-DDTHH:MM:SS.SS'
                      '\nYYYY-MM-DDTHH:MM:SS.SSZ'
                      '\nYYYY-MM-DDTHH:MM:SS.SS+HH'
                      '\nYYYY-MM-DDTHH:MM:SS.SS+HHMM'
                      '\nYYYY-MM-DDTHH:MM:SS.SS+HH:MM'
                      '\nYYYY-MM-DD HH:MM:SS.SS'
                      '\nYYYY-MM-DD HH:MM:SS.SSZ'
                      '\nYYYY-MM-DD HH:MM:SS.SS+HH'
                      '\nYYYY-MM-DD HH:MM:SS.SS+HHMM'
                      '\nYYYY-MM-DD HH:MM:SS.SS+HH:MM'
                     )

class Moment:
    """
    .. |UTC| replace:: UTC (Coordinated Universal Time)
    .. _UTC: https://en.wikipedia.org/wiki/Coordinated_Universal_Time

    :param utc: |UTC|_.
    :type utc: ``datetime.datetime``, ``float``, or ``str``

    The :py:attr:`pylightcurve.Moment` class represents a specific moment in time, after 1970.
    It allows for conversions to different time systems and operations between different
    :py:attr:`pylightcurve.Moment` objects.

    **Definition**

    A :py:attr:`pylightcurve.Moment` cab be defined through:

    * a ``datetime.datetime`` object, representing |UTC|_, with no option for time-zone);
    * a ``float`` object, representing |UTC|_, exprested in ulian date (days since noon on January 1, 4713 BC);
    * an ``str`` object, following the `ISO 8601 <https://en.wikipedia.org/wiki/ISO_8601>`_  full date and time format:

      ``YYYY-MM-DDTHH:MM:SS.mmmmmmTIMEZONE``,

      where:

      * ``YYYY``: year, after 1970;
      * ``MM``: months, from 1-12;
      * ``DD``: days, from 1 to 31;
      * ``T``: a separator, could be a space as well;
      * ``HH``: hours, from 1 to 23;
      * ``MM``: minutes, from 0 to 59;
      * ``SS``: seconds from 0 to 59;
      * ``.mmmmmm`` (optional): decimal points of seconds
        (if more than 6 characters are given, they will be rounded);
      * ``TIMEZONE`` (optional): the time-zone that the rest of the string refers to;
        ``Z``, alone, represents UTC+00, while other acceptable formats are:

        ``±hh:mm``, ``±hhmm``, or ``±hh``,

        where:

        - ``hh``: hours, from 1 to 12;
        - ``mm``: minutes, from 0 to 59.

    .. code-block:: python

        >>> import pylightcurve as plc
        >>> import datetime

        >>> # all the definitions below will create the same object
        >>> plc.Moment('2022-01-01T12:00:00')
        >>> plc.Moment('2022-01-01T12:00:00+00')
        >>> plc.Moment('2022-01-01T12:00:00+0000')
        >>> plc.Moment('2022-01-01T12:00:00+00:00')
        >>> plc.Moment('2022-01-01T13:00:00+01:00')
        >>> plc.Moment('2022-01-01T12:00:00Z')
        >>> plc.Moment('2022-01-01 12:00:00')
        >>> plc.Moment('2022-01-01 12:00:00+00')
        >>> plc.Moment('2022-01-01 12:00:00+0000')
        >>> plc.Moment('2022-01-01 12:00:00+00:00')
        >>> plc.Moment('2022-01-01 13:00:00+01:00')
        >>> plc.Moment('2022-01-01 12:00:00Z')
        >>> plc.Moment(2459581.0)
        >>> plc.Moment(datetime.datetime(2022,1,1,12,0,0))
        plc.Moment(2022-01-01T12:00:00 UTC)

    **Operations**

    A ``datetime.timedelta`` object can be added to or subtracted from a :py:attr:`pylightcurve.Moment` object,
    returning a new :py:attr:`pylightcurve.Moment` object.

    .. code-block:: python

        >>> import pylightcurve as plc
        >>> import datetime

        >>> plc.Moment('2022-01-01T12:00:00') + datetime.timedelta(hours=1)
        plc.Moment(2022-01-01T13:00:00 UTC)
        >>> plc.Moment('2022-01-01T12:00:00') - datetime.timedelta(hours=1)
        plc.Moment(2022-01-01T11:00:00 UTC)

    Also, a :py:attr:`pylightcurve.Moment` object can be subtracted from another
    :py:attr:`pylightcurve.Moment` object, returning a ``datetime.timedelta`` object.

    .. code-block:: python

        >>> import pylightcurve as plc
        >>> import datetime

        >>> plc.Moment('2022-01-01T13:00:00') - plc.Moment('2022-01-01T12:00:00')
        datetime.timedelta(seconds=3600)
        >>> plc.Moment('2022-01-01T12:00:00') - plc.Moment('2022-01-01T13:00:00')
        datetime.timedelta(days=-1, seconds=82800)

    **Methods**

    """

    def __init__(self, utc=None, jd_utc=None):

        if utc is None:
            if not isinstance(jd_utc, float) + isinstance(jd_utc, int):
                raise PyLCInputError('You need to provide either utc as a datetime.datetime or string object '
                                     'OR jd_utc as float or int object')

            else:
                if jd_utc - 2440588 < 0:
                    raise PyLCInputError('plc.Moment is valid only for dates after 01/01/1972.')
                self._utc = datetime.datetime(1970, 1, 1, 12, 0, 0, 0) + datetime.timedelta(days=(jd_utc - 2440588))

        elif jd_utc is None:
            if not isinstance(utc, datetime.datetime) + isinstance(utc, str):
                raise PyLCInputError('You need to provide either utc as a datetime.datetime or string object '
                                     'OR jd_utc as float or int object')

            if isinstance(utc, datetime.datetime):

                self._utc = utc

            elif isinstance(utc, str):

                try:
                    utc = utc.replace('T', ' ')

                    year, month, day = utc.split(' ')[0].split('-')
                    hour = utc.split(' ')[1].split(':')[0]
                    minutes = utc.split(' ')[1].split(':')[1]
                    seconds = ':'.join(utc.split(' ')[1].split(':')[2:])

                    if '+' in seconds:
                        time_zone_sign = 1
                        seconds, time_zone = seconds.split('+')
                    elif '-' in seconds:
                        time_zone_sign = -1
                        seconds, time_zone = seconds.split('-')
                    elif seconds[-1] == 'Z':
                        seconds = seconds[:-1]
                        time_zone_sign = 1
                        time_zone = '00:00'
                    else:
                        time_zone_sign = 1
                        time_zone = '00:00'
                except:
                    raise PyLCInputError(_date_format_error)

                if ':' in time_zone:
                    time_zone_hours = time_zone.split(':')[0]
                    time_zone_minutes = time_zone.split(':')[1]
                else:
                    if len(time_zone) == 4:
                        time_zone_hours = time_zone[:2]
                        time_zone_minutes = time_zone[2:]
                    elif len(time_zone) == 2:
                        time_zone_hours = time_zone
                        time_zone_minutes = '00'
                    else:
                        raise PyLCInputError(_date_format_error)

                if '.' in seconds:
                    seconds, subseconds = seconds.split('.')
                else:
                    subseconds = '0'

                for element in [
                    year, month, day, hour, minutes, seconds, time_zone_hours, time_zone_minutes
                ]:
                    try:
                        if int(element) != float(element) or len(element) < 2:
                            raise PyLCInputError(_date_format_error)
                    except:
                        raise PyLCInputError(_date_format_error)

                year = int(year)
                month = int(month)
                day = int(day)
                hour = int(hour)
                minutes = int(minutes)
                seconds = int(seconds)
                time_zone_hours = time_zone_sign * int(time_zone_hours)
                time_zone_minutes = time_zone_sign * int(time_zone_minutes)
                microsec = subseconds.ljust(6, '0')
                if len(microsec) > 6:
                    microsec = int(microsec[:6]) + int(round(float('0.' + microsec[6:])))
                else:
                    microsec = int(microsec)

                self._utc = (datetime.datetime(year, month, day, hour, minutes, seconds, microsec) -
                             datetime.timedelta(hours=time_zone_hours) - datetime.timedelta(minutes=time_zone_minutes))

        else:
            raise PyLCInputError('You need to provide either utc as a datetime.datetime or string object '
                                 'OR jd_utc as float or int object. Not both!')

        if self._utc.year < 1972:
            raise PyLCInputError('plc.Moment is valid only for dates after 01/01/1972.')

    def _datetime_to_jd_2part(self,  datetime_object):
        # datetime does not support BC dates so the difference is calculated from 1970
        # and the remaining days are added on top.
        delta = datetime_object - datetime.datetime(1970, 1, 1, 12, 0, 0, 0) + datetime.timedelta(days=2440588)
        return delta.days, (delta.seconds + delta.microseconds / 1000000.0) / 60.0 / 60.0 / 24.0

    def _datetime_to_jd(self,  datetime_object):
        return sum(self._datetime_to_jd_2part(datetime_object))

    @lru_cache()
    def utc(self):
        """

        :return: |UTC|_ at this moment.
        :rtype: ``datetime.datetime``
        """
        return self._utc

    @lru_cache()
    def jd_utc(self):
        """Converts :meth:`utc` to Julian date (days since noon on January 1, 4713 BC).

        :return: Julian date based on |UTC|_.
        :rtype: ``float``
        """
        return self._datetime_to_jd(self.utc())

    @lru_cache()
    def ut1_utc_diff(self):
        """
        Calls :py:attr:`pylightcurve.plc_data.earth_rotation` to return the UT1-UTC seconds at this moment.
        The (UT1-UTC) seconds are the difference between |UT1|_ and |UTC|_. The UT1-UTC seconds are measured
        and reported by the `International Earth Rotation and Reference Systems database
        <https://www.iers.org/IERS/EN/DataProducts/EarthOrientationData/eop.html>`_ under the
        `latest version of the finals.all (IAU2000) file
        <https://datacenter.iers.org/data/latestVersion/finals.all.iau1980.txt>`_.

        .. |UT1| replace:: UT1 (Universal Time)
        .. _UT1: https://en.wikipedia.org/wiki/Universal_Time

        :return: (UT1-UTC) seconds at this moment.
        :rtype: ``float``
        """
        return plc_data.earth_rotation(self.jd_utc())

    @lru_cache()
    def ut1(self):
        """
        Calculates the |UT1|_ at this moment by adding the (UT1-UTC) seconds to the |UTC|_.
        The UT1-UTC seconds are measured and reported by the `International Earth Rotation and Reference Systems database
        <https://www.iers.org/IERS/EN/DataProducts/EarthOrientationData/eop.html>`_ under the
        `latest version of the finals.all (IAU2000) file
        <https://datacenter.iers.org/data/latestVersion/finals.all.iau1980.txt>`_.
        The :py:attr:`pylightcurve.plc_data.earth_rotation` function called to extimate

        UT1 = :meth:`utc` + :meth:`ut1_utc_diff`

        :return: |UT1|_ at this moment.
        :rtype: ``datetime.datetime``
        """
        return self.utc() + datetime.timedelta(seconds=self.ut1_utc_diff())

    @lru_cache()
    def jd_ut1(self):
        """Converts :meth:`ut1` to Julian date (days since noon on January 1, 4713 BC).

        :return: Julian date based on |UT1|_.
        :rtype: ``float``
        """
        return self._datetime_to_jd(self.ut1())

    @lru_cache()
    def leap_seconds(self):
        """
        Calls :py:attr:`pylightcurve.plc_data.leap_seconds` to return the leap seconds at this moment.
        The leap seconds are the difference between |UTC|_ and |TAI|_.

        .. |TAI| replace:: TAI (International Atomic Time)
        .. _TAI: https://en.wikipedia.org/wiki/International_Atomic_Time

        :return: Leap seconds at this moment.
        :rtype: ``float``
        """
        return plc_data.leap_seconds(self.utc())

    @lru_cache()
    def tai(self):
        """
        Calculates the |TAI|_ at this moment

        TAI = :meth:`utc` + :meth:`leap_seconds`

        .. |TAI| replace:: TAI (International Atomic Time)
        .. _TAI: https://en.wikipedia.org/wiki/International_Atomic_Time

        :return: |TAI|_ at this moment.
        :rtype: ``datetime.datetime``
        """
        return self.utc() + datetime.timedelta(seconds=self.leap_seconds())

    @lru_cache()
    def jd_tai(self):
        """Converts :meth:`tai` to Julian date (days since noon on January 1, 4713 BC).

        :return: Julian date based on |TAI|_.
        :rtype: ``float``
        """
        return self._datetime_to_jd(self.tai())

    @lru_cache()
    def tt(self):
        """
        Calculates the |TT|_ at this moment

        TT = :meth:`utc` + :meth:`leap_seconds` + 32.184s

        .. |TT| replace:: TT (Terrestrial Time)
        .. _TT: https://en.wikipedia.org/wiki/Terrestrial_Time

        :return: |TT|_ at this moment.
        :rtype: ``datetime.datetime``
        """
        return self.utc() + + datetime.timedelta(seconds=self.leap_seconds() + 32.184)

    @lru_cache()
    def jd_tt(self):
        """Converts :meth:`tt` to Julian date (days since noon on January 1, 4713 BC).

        :return: Julian date based on |TT|_.
        :rtype: ``float``
        """
        return self._datetime_to_jd(self.tt())

    @lru_cache()
    def era(self):
        """

        .. math::
            GMST (arcsec) = 1296000 \\times (1.00273781191135448 \\times T_u + 0.7790572732640)

                   + 0.014506
                   + 4612.156534 t
                   + 1.3915817 t^2
                   - 0.00000044 t^3

                   - 0.000029956 t^4
                   - 0.0000000368 t^5


        where:

        * :math:`T_u` = :meth:`jd_ut1` - 2451545, and
        * :math:`t` = (:meth:`jd_tt` - 2451545) / 36525;

        :return: Greenwich Mean Sidereal Time (GMST) at this moment.
        :rtype: ``datetime.datetime``
        """

        def mul(a,b):

            xx_int, xx_dec = divmod(a[1]*b[0], 1)
            yy_int, yy_dec = divmod(a[0]*b[1], 1)
            zz_int, zz_dec = divmod(xx_dec + yy_dec + a[1]*b[1], 1)

            return a[0]*b[0] + xx_int + yy_int + zz_int, zz_dec

        def sum(a,b):
            xx_int, xx_dec = divmod(a[1] + b[1], 1)
            return a[0] + b[0] + xx_int, xx_dec

        def sub(a,b):
            xx_int, xx_dec = divmod(a[1] - b[1], 1)
            return a[0] - b[0] + xx_int, xx_dec

        tu = sub(self._datetime_to_jd_2part(self.ut1()), (2451545, 0))

        era = sum(mul((1, 0.00273781191135448), tu), (0, 0.7790572732640))
        era = (0, era[1])
        era = mul(era, (1296000, 0))

        return Degrees(0, 0, era[0]) + Degrees(0, 0, era[1])

    def __add__(self, other):
        if isinstance(other, datetime.timedelta):
            return Moment(self.utc() + other)
        else:
            raise PyLCInputError('Only a datetime.timedelta object can be added to a {0} object.'.format(
                self._get_class()))

    def __sub__(self, other):
        if isinstance(other, datetime.timedelta):
            return Moment(self.utc() - other)
        elif isinstance(other, Moment):
            return self.utc() - other.utc()
        else:
            raise PyLCInputError('Only a datetime.timedelta or a plc.Moment object can be subtracted from a {0} object.'.format(
                self._get_class()))

    @lru_cache()
    def _get_class(self):
        return 'plc.{0}'.format(str(self.__class__).split('.')[-1][:-2])

    def __repr__(self):
        return '{0}({1} UTC)'.format(self._get_class(), self.utc().isoformat())


def _is_time(item):
    return isinstance(item, Moment)


def _request_time(item):
    if _is_time(item):
        pass
    else:
        raise PyLCInputError('A plc.Moment object is required.')


def now():
    return Moment(datetime.datetime.now(datetime.timezone.utc).isoformat())
