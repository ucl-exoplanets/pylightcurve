
import datetime
import numpy as np
import pytest
import pylightcurve as plc

from pylightcurve.spacetime.observatories import _request_observatory, _Horizon


def test_observing():

    with pytest.raises(plc.PyLCInputError):
        _Horizon('xxx').horizon(plc.Degrees(10)).deg()

    with pytest.raises(plc.PyLCInputError):
        _Horizon({}).horizon(plc.Degrees(10)).deg()

    assert _Horizon(20).horizon(plc.Degrees(10)).deg() == 20
    assert _Horizon('0 20\n90 20\n180 20\n270 20\n').horizon(plc.Degrees(10)).deg() == 20
    assert _Horizon([[0,20],[90,20],[180,20],[270,20]]).horizon(plc.Degrees(10)).deg() == 20
    assert _Horizon(np.array([[0,20],[90,20],[180,20],[270,20]])).horizon(plc.Degrees(10)).deg() == 20

    with pytest.raises(plc.PyLCInputError):
            plc.Observatory(plc.Degrees(140), plc.Degrees(23), 2)

    with pytest.raises(plc.PyLCInputError):
        plc.Observatory(plc.Degrees(40), plc.Degrees(23), 'x')

    with pytest.raises(plc.PyLCInputError):
        plc.Observatory(plc.Degrees(40), plc.Degrees(23), 20)

    with pytest.raises(plc.PyLCError):
        _request_observatory('a')

    observatory = plc.Observatory(plc.Degrees(40), plc.Degrees(-23), -2)
    print(observatory.__repr__)
    _request_observatory(observatory)
    assert observatory.coord() == '+40:00:00.0 337:00:00.0'
    assert observatory.lt(plc.Moment(jd_utc=2459903)) == datetime.datetime(2022, 11, 19, 10, 0, 0)
    assert round(observatory.lera(plc.Moment(jd_utc=2459903)).deg(), 1) == 215.20

    observatory = plc.Observatory(plc.Degrees(40), plc.Degrees(23), 2)
    print(observatory.__repr__)
    _request_observatory(observatory)
    assert observatory.coord() == '+40:00:00.0 23:00:00.0'
    assert observatory.lt(plc.Moment(jd_utc=2459903)) == datetime.datetime(2022, 11, 19, 14, 0, 0)
    assert round(observatory.lera(plc.Moment(jd_utc=2459903)).deg(), 1) == 261.20

    moment = plc.Moment(jd_utc=2459903)
    target = plc.simbad_search_by_name('XO-1')
    assert round(observatory.target_azimuth_altitude(target,  moment)[0].deg()) == 62.0
    assert round(observatory.airmass(target,  moment), 1) == 1.1
    assert observatory.is_target_visible(target,  moment)
    assert len(observatory.target_horizon_crossings(target, moment, 1)) == 2
    assert len(observatory.target_altitude_crossings(target, moment, 1, plc.Degrees(20))) == 2

    moment = plc.Moment(jd_utc=2459903)
    target = plc.simbad_search_by_name('HD189733')
    assert round(observatory.target_azimuth_altitude(target,  moment)[0].deg()) == 286.0
    assert round(observatory.airmass(target,  moment), 1) == 1.3
    assert observatory.is_target_visible(target,  moment)
    assert len(observatory.target_horizon_crossings(target, moment, 1)) == 2
    assert len(observatory.target_altitude_crossings(target, moment, 1, plc.Degrees(20))) == 2

    moment = plc.Moment('2022-05-01T12:00:00')
    target = plc.simbad_search_by_name('XO-1')
    assert len(observatory.periodic_events_visibility(target, moment, 1, 2455787.553228, 'BJD_TDB', 3.94150468, 5/24)) == 1

