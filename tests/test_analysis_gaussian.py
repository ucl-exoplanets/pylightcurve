
import numpy as np
import pytest
import pylightcurve as plc


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

