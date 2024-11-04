
import numpy as np
import pylightcurve as plc


def test_gaussian():

    data_x = np.arange(0, 10, 0.1)
    data_y = plc.gaussian(data_x, 100.0, 10.0, 5.0, 1.0) + np.random.normal(0, 0.1, len(data_x))

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

    results = plc.fit_gaussian(data_x, data_y, positive=True, sampled=True, point_x=5, norm=1000, floor=100, sigma=1)

    assert len(results) == 2
    assert round(results[0][0]) in [999, 1000, 1001]
    assert round(results[0][1]) == 100
    assert round(results[0][2]) == 5
    assert round(results[0][3]) == 1

    data_x = np.arange(0, 10, 0.1)
    data_y = np.arange(0, 10, 0.1)
    data_x, data_y = np.meshgrid(data_x, data_y)
    data_z = plc.two_d_gaussian(data_x, data_y, 100, 10, 5, 5, 1, 1, 0) + np.random.normal(0, 0.1, data_x.shape)

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

    data_x = np.arange(0, 10, 0.1)
    data_y = np.arange(0, 10, 0.1)
    data_x, data_y = np.meshgrid(data_x, data_y)
    data_z = plc.two_d_gaussian(data_x, data_y, 100, 10, 5, 5, 3, 1, np.pi/4) + np.random.normal(0, 0.1, data_x.shape)
    errors = np.sqrt(data_z)

    results = plc.fit_two_d_gaussian(data_x, data_y, data_z, errors=errors,
                                     point_xy=(5, 5), norm=100, sigma=1, floor=10, theta=np.pi/4)

    assert len(results) == 2
    assert round(results[0][0]) == 100
    assert round(results[0][1]) == 10
    assert round(results[0][2]) == 5
    assert round(results[0][3]) == 5
    assert round(results[0][4]) == 3
    assert round(results[0][5]) == 1
    assert round(results[0][6]) == 1
