import numpy as np
import pytest
import pylightcurve as plc

from pylightcurve.spacetime.angles import _is_angle, _request_angle


def test_angles():

    for test_angle in [plc.Degrees(54.61),
                       plc.Degrees('54.61'),
                       plc.Degrees(54, 36, 36.0),
                       plc.Degrees('54:36:36.0'),
                       plc.Hours(54.61 / 15.0),
                       plc.Hours(54 / 15.0, 36.0 / 15.0, 36.0 / 15.0),
                       plc.Rad(54.61 * np.pi / 180.0)]:

        assert round(test_angle.arcseconds, 1) == 196596.0
        assert [round(ff, 1) for ff in test_angle.dms().list] == [54, 36, 36]
        assert round(test_angle.deg(), 2) == 54.61
        assert round(test_angle.deg_coord(), 2) == 54.61
        assert [round(ff, 1) for ff in test_angle.hms().list] == [3, 38, 26.4]
        assert round(test_angle.hours(), 4) == 3.6407
        assert round(test_angle.rad(), 10) == 0.9531243045
        assert round(test_angle.sin(), 10) == 0.8152288870
        assert round(test_angle.cos(), 10) == 0.5791388969
        assert round(test_angle.tan(), 10) == 1.4076569392
        assert round(plc.arcsin(test_angle.sin()).deg(), 2) == 54.61
        assert round(plc.arccos(test_angle.cos()).deg(), 2) == 54.61
        assert round(plc.arctan(test_angle.tan()).deg(), 2) == 54.61
        assert _is_angle(test_angle) is True
        _request_angle(test_angle)

    for test_angle in [plc.Degrees(-54.61),
                       plc.Degrees('-54.61'),
                       plc.Degrees(-54, -36, -36.0),
                       plc.Degrees('-54:36:36.0'),
                       plc.Hours(-54.61 / 15.0),
                       plc.Hours(-54 / 15.0, -36.0 / 15.0, -36.0 / 15.0),
                       plc.Rad(-54.61 * np.pi / 180.0)]:

        assert round(test_angle.arcseconds, 1) == 1099404.0
        assert [round(ff, 1) for ff in test_angle.dms().list] == [305, 23, 24.0]
        assert round(test_angle.deg(), 2) == 305.39
        assert round(test_angle.deg_coord(), 2) == -54.61
        assert [round(ff, 1) for ff in test_angle.hms().list] == [20, 21, 33.6]
        assert round(test_angle.hours(), 4) == 20.3593
        assert round(test_angle.rad(), 10) == 5.3300610027
        assert round(test_angle.sin(), 10) == -0.8152288870
        assert round(test_angle.cos(), 10) == 0.5791388969
        assert round(test_angle.tan(), 10) == -1.4076569392
        assert _is_angle(test_angle) is True
        _request_angle(test_angle)

    a1, a2 = plc.Degrees(55.0), plc.Degrees(5.0)
    assert (a1 + a2).deg() == 60.0
    assert (a1 - a2).deg() == 50.0
    assert (a1 * 2).deg() == 110.0
    assert (2 * a1).deg() == 110.0
    assert (a1 / 2).deg() == 27.5

    assert plc.Degrees(5.0).deg_coord() == 5.0
    assert plc.Degrees(5.0).dms_coord() == '+05:00:00.0'
    assert plc.Degrees(95.0).deg_coord() == 85.0
    assert plc.Degrees(95.0).dms_coord() == '+85:00:00.0'
    assert plc.Degrees(185.0).deg_coord() == -5.0
    assert plc.Degrees(185.0).dms_coord() == '-05:00:00.0'
    assert plc.Degrees(275.0).deg_coord() == -85.0
    assert plc.Degrees(275.0).dms_coord() == '-85:00:00.0'

    assert plc.pi.rad() == np.pi

    print(test_angle.__repr__, test_angle.hms().__repr__, test_angle.dms().__repr__)

    with pytest.raises(plc.PyLCInputError):
        plc.Degrees('a')
    with pytest.raises(plc.PyLCInputError):
        plc.Degrees('35:00:00', 10)
    with pytest.raises(plc.PyLCInputError):
        plc.Degrees(['a'])
    with pytest.raises(plc.PyLCInputError):
        plc.Degrees(35, -10)
    with pytest.raises(plc.PyLCInputError):
        plc.Rad('a')
    with pytest.raises(plc.PyLCError):
        a1 * a2
    with pytest.raises(plc.PyLCError):
        a1 / a2
    with pytest.raises(plc.PyLCError):
        _request_angle('a')

    assert _is_angle('a') is False
