
import datetime
import pytest
import pylightcurve as plc
from astropy.time import Time

from pylightcurve.spacetime.moments import _request_time


def test_times():

    astropy_test_era = Time('2022-01-01T12:00:00.0').earth_rotation_angle(0).arcsec
    dtime = datetime.timedelta(hours=1)
    ref = plc.Moment('2022-01-01T11:30:00.0')
    now = plc.now()

    for test_time in [
        plc.Moment(datetime.datetime(2022, 1, 1, 12, 0, 0, microsecond=0)),
        plc.Moment('2022-01-01T12:00:00.0'),
        plc.Moment('2022-01-01T12:00:00.0000001'),
        plc.Moment('2022-01-01T12:00:00.0Z'),
        plc.Moment('2022-01-01T12:00:00.0+00'),
        plc.Moment('2022-01-01T12:00:00.0+00:00'),
        plc.Moment('2022-01-01T12:30:00.0+0030'),
        plc.Moment('2022-01-01T11:30:00.0-0030'),
        plc.Moment(jd_utc=2459581.0),
    ]:

        assert test_time.leap_seconds() == 37
        assert round(test_time.ut1_utc_diff(), 6) == -0.110451

        assert (test_time + dtime).utc() == datetime.datetime(2022, 1, 1, 13, 0, 0)
        assert (test_time - dtime).utc() == datetime.datetime(2022, 1, 1, 11, 0, 0)
        assert test_time - ref == datetime.timedelta(hours=0.5)

        assert test_time.utc() == datetime.datetime(2022, 1, 1, 12, 0, 0)
        assert test_time.tai() == datetime.datetime(2022, 1, 1, 12, 0, 37)
        assert test_time.tt() == datetime.datetime(2022, 1, 1, 12, 1, 9, 184000)
        assert test_time.ut1() == datetime.datetime(2022, 1, 1, 11, 59, 59, 889549)

        assert round(test_time.jd_utc(), 6) == 2459581.0
        assert round(test_time.jd_tai(), 6) == 2459581.000428
        assert round(test_time.jd_tt(), 6) == 2459581.000801
        assert round(test_time.jd_ut1(), 6) == 2459580.999999

        assert round(test_time.era()._arcseconds, 3) == round(astropy_test_era, 3)

        print(test_time.__repr__)
        _request_time(test_time)

    for utc in [
        10,
        'a',
        '2022-01-01T12:00:00.0+00100',
        '2022-01-01T1:00:00.0',
        '2022-01-01T1.1:00:00.0',
        '2022-01-01T1a:00:00.0',
        '1960-01-01T12:00:00',
    ]:
        with pytest.raises(plc.PyLCInputError):
            plc.Moment(utc)

    for utc in [
        '2022-01-01T121:00:00.0',
    ]:
        with pytest.raises(ValueError):
            plc.Moment(utc)

    for jd_utc in [
            10,
            'a',
        ]:
            with pytest.raises(plc.PyLCInputError):
                plc.Moment(jd_utc=jd_utc)

    with pytest.raises(plc.PyLCInputError):
        plc.Moment('2022-01-01T12:00:00.0', 2459581.0)

    with pytest.raises(plc.PyLCInputError):
        _request_time('a')

    with pytest.raises(plc.PyLCInputError):
        _ = test_time + 10

    with pytest.raises(plc.PyLCInputError):
        _ = test_time - 10
