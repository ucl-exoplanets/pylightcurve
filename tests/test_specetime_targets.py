
import pytest
import numpy as np
import pylightcurve as plc

from pylightcurve.spacetime.targets import _request_target


def test_targets():

    with pytest.raises(plc.PyLCInputError):
        plc.FixedTarget('a', plc.Degrees(23.5))

    with pytest.raises(plc.PyLCInputError):
        plc.FixedTarget(plc.Degrees(40.4), plc.Degrees(170))

    with pytest.raises(plc.PyLCInputError):
        _request_target('a')

    test_fixedstar = plc.simbad_search_by_name('XO-1')
    test_fixedstar.reset(plc.Moment('2022-03-20T23:00:00'))
    print(test_fixedstar.__repr__)
    _request_target(test_fixedstar)

    test_fixedstar = plc.FixedTarget(plc.Hours(3, 15), plc.Degrees(45.7))
    test_fixedstar.reset(plc.Moment('2022-03-20T23:00:00'))
    print(test_fixedstar.__repr__)
    _request_target(test_fixedstar)

    assert test_fixedstar.coord() == '03:15:00.0 +45:42:00.0'
    assert test_fixedstar.distance_on_sphere(test_fixedstar).deg() == 0
    assert round(test_fixedstar.convert_to_bjd_tdb(2458485.0, 'JD_UTC'), 7) == 2458485.0046340
    assert round(test_fixedstar.convert_to_bjd_tdb(58484.5, 'MJD_UTC'), 7) == 2458485.0046340
    assert round(test_fixedstar.convert_to_bjd_tdb(2458485.0046340, 'BJD_TDB'), 7) == 2458485.0046340
    assert round(test_fixedstar.convert_to_bjd_tdb(2458485.0038333, 'BJD_UTC'), 7) == 2458485.0046340
    assert round(test_fixedstar.convert_to_bjd_tdb(2458485.0046033, 'HJD_TDB'), 7) == 2458485.0046340
    assert round(test_fixedstar.convert_to_bjd_tdb(2458485.00380255, 'HJD_UTC'), 7) == 2458485.0046340
    assert round(test_fixedstar.convert_to_bjd_tdb(np.array([2458485.00000000]), 'JD_UTC')[0], 7) == 2458485.0046340
    assert round(test_fixedstar.convert_to_bjd_tdb(np.array([0058484.50000000]), 'MJD_UTC')[0], 7) == 2458485.0046340
    assert round(test_fixedstar.convert_to_bjd_tdb(np.array([2458485.00463400]), 'BJD_TDB')[0], 7) == 2458485.0046340
    assert round(test_fixedstar.convert_to_bjd_tdb(np.array([2458485.00383330]), 'BJD_UTC')[0], 7) == 2458485.0046340
    assert round(test_fixedstar.convert_to_bjd_tdb(np.array([2458485.00460330]), 'HJD_TDB')[0], 7) == 2458485.0046340
    assert round(test_fixedstar.convert_to_bjd_tdb(np.array([2458485.00380255]), 'HJD_UTC')[0], 7) == 2458485.0046340
    assert round(test_fixedstar.convert_to_jd_utc(2458485.0, 'JD_UTC'), 7) == 2458485.0
    assert round(test_fixedstar.convert_to_jd_utc(58484.5, 'MJD_UTC'), 7) == 2458485.0
    assert round(test_fixedstar.convert_to_jd_utc(2458485.0046340, 'BJD_TDB'), 7) == 2458485.0
    assert round(test_fixedstar.convert_to_jd_utc(2458485.0038333, 'BJD_UTC'), 7) == 2458485.0
    assert round(test_fixedstar.convert_to_jd_utc(2458485.0046033, 'HJD_TDB'), 7) == 2458485.0
    assert round(test_fixedstar.convert_to_jd_utc(2458485.00380255, 'HJD_UTC'), 7) == 2458485.0
    assert round(test_fixedstar.convert_to_jd_utc(np.array([2458485.00000000]), 'JD_UTC')[0], 7) == 2458485.0
    assert round(test_fixedstar.convert_to_jd_utc(np.array([0058484.50000000]), 'MJD_UTC')[0], 7) == 2458485.0
    assert round(test_fixedstar.convert_to_jd_utc(np.array([2458485.00463400]), 'BJD_TDB')[0], 7) == 2458485.0
    assert round(test_fixedstar.convert_to_jd_utc(np.array([2458485.00383330]), 'BJD_UTC')[0], 7) == 2458485.0
    assert round(test_fixedstar.convert_to_jd_utc(np.array([2458485.00460330]), 'HJD_TDB')[0], 7) == 2458485.0
    assert round(test_fixedstar.convert_to_jd_utc(np.array([2458485.00380255]), 'HJD_UTC')[0], 7) == 2458485.0

    with pytest.raises(plc.PyLCInputError):
        test_fixedstar.convert_to_bjd_tdb(2458485.0, 'aaa')

    with pytest.raises(plc.PyLCInputError):
        test_fixedstar.convert_to_bjd_tdb([2458485.0], 'aaa')

    with pytest.raises(plc.PyLCInputError):
        test_fixedstar.convert_to_bjd_tdb('a', 'JD_UTC')

    with pytest.raises(plc.PyLCInputError):
        test_fixedstar.convert_to_jd_utc(2458485.0, 'aaa')

    with pytest.raises(plc.PyLCInputError):
        test_fixedstar.convert_to_jd_utc([2458485.0], 'aaa')

    with pytest.raises(plc.PyLCInputError):
        test_fixedstar.convert_to_jd_utc('a', 'JD_UTC')

    sun = plc.Sun(plc.Moment('2022-03-20T23:00:00'))
    print(sun.__repr__)
    assert round(sun.ra.deg(), 1) == 0

    sun.reset(plc.Moment('2022-09-23T09:00:00'))
    print(sun.__repr__)
    assert round(sun.ra.deg(), 1) == 180

    moon = plc.Moon(plc.Moment('2022-03-20T23:00:00'))
    print(moon.__repr__)
    assert round(moon.ra.deg(), 1) == 212.3
    assert round(moon.illumination(), 1) == 0.9

    moon.reset(plc.Moment('2022-09-23T09:00:00'))
    print(moon.__repr__)
    assert round(moon.ra.deg(), 1) == 154.1
    assert round(moon.illumination(), 1) == 0.1
