import pytest
import pylightcurve as plc

from pylightcurve.spacetime.targets import _is_target, _request_target


def test_targets():

    test_fixedstar = plc.FixedTarget(plc.Hours(3, 15), plc.Degrees(-45.7))
    assert test_fixedstar.coord == '03:15:00.0 -45:42:00.0'

    ra, dec = plc.Hours('03 15 00'), plc.Degrees('45 42 00')
    test_fixedstar = plc.FixedTarget(ra, dec)
    assert test_fixedstar.coord == '03:15:00.0 +45:42:00.0'

    print(test_fixedstar.distance_on_sphere(test_fixedstar).deg())

    assert test_fixedstar.distance_on_sphere(test_fixedstar).deg() == 0

    _request_target(test_fixedstar)
    assert _is_target(test_fixedstar) is True
    print(test_fixedstar.__repr__)

    assert round(test_fixedstar.convert_to_bjd_tdb(2458485.0, 'JD_UTC'), 7) == 2458485.0046340
    assert round(test_fixedstar.convert_to_bjd_tdb(58484.5, 'MJD_UTC'), 7) == 2458485.0046340
    assert round(test_fixedstar.convert_to_bjd_tdb(2458485.0046340, 'BJD_TDB'), 7) == 2458485.0046340
    assert round(test_fixedstar.convert_to_bjd_tdb(2458485.0038333, 'BJD_UTC'), 7) == 2458485.0046340
    assert round(test_fixedstar.convert_to_bjd_tdb(2458485.0046033, 'HJD_TDB'), 7) == 2458485.0046340
    assert round(test_fixedstar.convert_to_bjd_tdb(2458485.00380255, 'HJD_UTC'), 7) == 2458485.0046340

    assert round(test_fixedstar.convert_to_bjd_tdb([2458485.0], 'JD_UTC')[0], 7) == 2458485.0046340
    assert round(test_fixedstar.convert_to_bjd_tdb([58484.5], 'MJD_UTC')[0], 7) == 2458485.0046340
    assert round(test_fixedstar.convert_to_bjd_tdb([2458485.0046340], 'BJD_TDB')[0], 7) == 2458485.0046340
    assert round(test_fixedstar.convert_to_bjd_tdb([2458485.0038333], 'BJD_UTC')[0], 7) == 2458485.0046340
    assert round(test_fixedstar.convert_to_bjd_tdb([2458485.0046033], 'HJD_TDB')[0], 7) == 2458485.0046340
    assert round(test_fixedstar.convert_to_bjd_tdb([2458485.00380255], 'HJD_UTC')[0], 7) == 2458485.0046340

    assert round(test_fixedstar.convert_to_jd(2458485.0, 'JD_UTC'), 7) == 2458485.0
    assert round(test_fixedstar.convert_to_jd(58484.5, 'MJD_UTC'), 7) == 2458485.0
    assert round(test_fixedstar.convert_to_jd(2458485.0046340, 'BJD_TDB'), 7) == 2458485.0
    assert round(test_fixedstar.convert_to_jd(2458485.0038333, 'BJD_UTC'), 7) == 2458485.0
    assert round(test_fixedstar.convert_to_jd(2458485.0046033, 'HJD_TDB'), 7) == 2458485.0
    assert round(test_fixedstar.convert_to_jd(2458485.00380255, 'HJD_UTC'), 7) == 2458485.0

    assert round(test_fixedstar.convert_to_jd([2458485.0], 'JD_UTC')[0], 7) == 2458485.0
    assert round(test_fixedstar.convert_to_jd([58484.5], 'MJD_UTC')[0], 7) == 2458485.0
    assert round(test_fixedstar.convert_to_jd([2458485.0046340], 'BJD_TDB')[0], 7) == 2458485.0
    assert round(test_fixedstar.convert_to_jd([2458485.0038333], 'BJD_UTC')[0], 7) == 2458485.0
    assert round(test_fixedstar.convert_to_jd([2458485.0046033], 'HJD_TDB')[0], 7) == 2458485.0
    assert round(test_fixedstar.convert_to_jd([2458485.00380255], 'HJD_UTC')[0], 7) == 2458485.0

    assert round(test_fixedstar.convert_to_mjd(2458485.0, 'JD_UTC'), 7) == 58484.5
    assert round(test_fixedstar.convert_to_mjd(58484.5, 'MJD_UTC'), 7) == 58484.5
    assert round(test_fixedstar.convert_to_mjd(2458485.0046340, 'BJD_TDB'), 7) == 58484.5
    assert round(test_fixedstar.convert_to_mjd(2458485.0038333, 'BJD_UTC'), 7) == 58484.5
    assert round(test_fixedstar.convert_to_mjd(2458485.0046033, 'HJD_TDB'), 7) == 58484.5
    assert round(test_fixedstar.convert_to_mjd(2458485.00380255, 'HJD_UTC'), 7) == 58484.5

    assert round(test_fixedstar.convert_to_mjd([2458485.0], 'JD_UTC')[0], 7) == 58484.5
    assert round(test_fixedstar.convert_to_mjd([58484.5], 'MJD_UTC')[0], 7) == 58484.5
    assert round(test_fixedstar.convert_to_mjd([2458485.0046340], 'BJD_TDB')[0], 7) == 58484.5
    assert round(test_fixedstar.convert_to_mjd([2458485.0038333], 'BJD_UTC')[0], 7) == 58484.5
    assert round(test_fixedstar.convert_to_mjd([2458485.0046033], 'HJD_TDB')[0], 7) == 58484.5
    assert round(test_fixedstar.convert_to_mjd([2458485.00380255], 'HJD_UTC')[0], 7) == 58484.5

    with pytest.raises(plc.PyLCInputError):
        test_fixedstar.convert_to_bjd_tdb(2458485.0, 'aaa')
    with pytest.raises(plc.PyLCInputError):
        test_fixedstar.convert_to_bjd_tdb([2458485.0], 'aaa')
    with pytest.raises(plc.PyLCInputError):
        test_fixedstar.convert_to_bjd_tdb('a', 'JD_UTC')

    with pytest.raises(plc.PyLCInputError):
        test_fixedstar.convert_to_jd(2458485.0, 'aaa')
    with pytest.raises(plc.PyLCInputError):
        test_fixedstar.convert_to_jd([2458485.0], 'aaa')
    with pytest.raises(plc.PyLCInputError):
        test_fixedstar.convert_to_jd('a', 'JD_UTC')

    with pytest.raises(plc.PyLCInputError):
        test_fixedstar.convert_to_mjd(2458485.0, 'aaa')
    with pytest.raises(plc.PyLCInputError):
        test_fixedstar.convert_to_mjd([2458485.0], 'aaa')
    with pytest.raises(plc.PyLCInputError):
        test_fixedstar.convert_to_mjd('a', 'JD_UTC')

    with pytest.raises(plc.PyLCInputError):
        plc.FixedTarget('a', plc.Degrees(23.5))
    with pytest.raises(plc.PyLCInputError):
        plc.FixedTarget(plc.Degrees(40.4), plc.Degrees(170))
    with pytest.raises(plc.PyLCInputError):
        plc.FixedTarget(plc.Degrees(40.4), 'a')

    with pytest.raises(plc.PyLCInputError):
        _request_target('a')

    assert _is_target('a') is False

