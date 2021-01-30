
import numpy as np
import pytest
import pylightcurve as plc
from pylightcurve.models.exoplanet_lc import _get_filter


def test_gaussian():

    data_x = np.arange(0, 10, 0.01)
    data_y = plc.gaussian(data_x, 100.0, 10.0, 5.0, 1.0)

    results = plc.fit_gaussian(data_x, data_y)

    assert len(results) == 2
    assert round(results[0][0]) == 100
    assert round(results[0][1]) == 10
    assert round(results[0][2]) == 5
    assert round(results[0][3]) == 1

    results = plc.fit_gaussian(data_x, data_y, positive=True)

    assert len(results) == 2
    assert round(results[0][0]) == 100
    assert round(results[0][1]) == 10
    assert round(results[0][2]) == 5
    assert round(results[0][3]) == 1

    results = plc.fit_gaussian(data_x, data_y, positive=True, sampled=True)

    assert len(results) == 2
    assert round(results[0][0]) == 10000
    assert round(results[0][1]) == 1000
    assert round(results[0][2]) == 5
    assert round(results[0][3]) == 1

    data_x = np.arange(0, 10, 0.01)
    data_y = np.arange(0, 10, 0.01)
    data_x, data_y = np.meshgrid(data_x, data_y)
    data_z = plc.two_d_gaussian(data_x, data_y, 100, 10, 5, 5, 1, 1, 0)

    results = plc.fit_two_d_gaussian(data_x, data_y, data_z)

    assert len(results) == 2
    assert round(results[0][0]) == 100
    assert round(results[0][1]) == 10
    assert round(results[0][2]) == 5
    assert round(results[0][3]) == 5
    assert round(results[0][4]) == 1
    assert round(results[0][5]) == 1

    results = plc.fit_two_d_gaussian(data_x, data_y, data_z, positive=True)

    assert len(results) == 2
    assert round(results[0][0]) == 100
    assert round(results[0][1]) == 10
    assert round(results[0][2]) == 5
    assert round(results[0][3]) == 5
    assert round(results[0][4]) == 1
    assert round(results[0][5]) == 1

    results = plc.fit_two_d_gaussian(data_x, data_y, data_z, symmetric=True)

    assert len(results) == 2
    assert round(results[0][0]) == 100
    assert round(results[0][1]) == 10
    assert round(results[0][2]) == 5
    assert round(results[0][3]) == 5
    assert round(results[0][4]) == 1

    results = plc.fit_two_d_gaussian(data_x, data_y, data_z, positive=True, symmetric=True)

    assert len(results) == 2
    assert round(results[0][0]) == 100
    assert round(results[0][1]) == 10
    assert round(results[0][2]) == 5
    assert round(results[0][3]) == 5
    assert round(results[0][4]) == 1

    results = plc.fit_two_d_gaussian(data_x, data_y, data_z, point_xy=(5, 5), sigma=None)

    assert len(results) == 2
    assert round(results[0][0]) == 100
    assert round(results[0][1]) == 10
    assert round(results[0][2]) == 5
    assert round(results[0][3]) == 5
    assert round(results[0][4]) == 1
    assert round(results[0][5]) == 1

    results = plc.fit_two_d_gaussian(data_x, data_y, data_z, sigma=1)

    assert len(results) == 2
    assert round(results[0][0]) == 100
    assert round(results[0][1]) == 10
    assert round(results[0][2]) == 5
    assert round(results[0][3]) == 5
    assert round(results[0][4]) == 1
    assert round(results[0][5]) == 1
