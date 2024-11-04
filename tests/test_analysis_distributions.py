
import pytest
import numpy as np
import pylightcurve as plc


def test_distribution_one_d():

    test_mean = 10.0
    test_std = 1.0
    distr = np.random.normal(test_mean, test_std, 100000)
    outlier_distr = np.append(distr, 1000000)
    narrow_distr = np.append(np.ones(100000) * test_mean, np.ones(10) * test_mean * 3)

    for a in [plc.one_d_distribution(distr),
              plc.one_d_distribution(distr, samples=10000),
              plc.one_d_distribution(outlier_distr, mad_filter=5)
              ]:
        assert round(a[0][np.argmax(a[1])]) in [9, 10, 11]
        assert round(5 * (a[0][1] - a[0][0])) == test_std

    for a in [
        plc.one_d_distribution(distr, confidence_interval=True),
        plc.one_d_distribution(distr, confidence_interval=0.68),
        plc.one_d_distribution(distr, gaussian_fit=True)
    ]:
        assert len(a) == 5
        assert round(a[0][np.argmax(a[1])]) in [9, 10, 11]
        assert round(5 * (a[0][1] - a[0][0])) == test_std
        assert round(a[2]) == test_mean
        assert round(a[3]) == test_std
        assert round(a[4]) == test_std

    b = plc.one_d_distribution(narrow_distr)
    assert round(b[0][np.argmax(b[1])]) == test_mean

    c = plc.one_d_distribution(distr, abs_step=1.6, min_value=5.0, max_value=15.0)
    assert round(c[0][np.argmax(c[1])], 1) == 10.6
    assert round(c[0][1]-c[0][0], 1) == 1.6
    assert round(c[0][0], 1) == 5.8
    assert round(c[0][-1], 1) == 15.4

    with pytest.raises(plc.PyLCProcessError):
        plc.one_d_distribution(distr, confidence_interval=0.68, gaussian_fit=True)


def test_distribution_two_d():

    test_mean_1 = 10.0
    test_std_1 = 1.0
    test_mean_2 = 10.0
    test_std_2 = 1.0

    distr_1 = np.random.normal(test_mean_1, test_std_1, 100000)
    distr_2 = np.random.normal(test_mean_2, test_std_2, 100000)
    outlier_distr_1 = np.append(distr_1, 100000)
    outlier_distr_2 = np.append(distr_2, 100000)
    narrow_distr_1 = np.append(np.ones(100000) * test_mean_1, np.ones(10) * test_mean_1 * 3)
    narrow_distr_2 = np.append(np.ones(100000) * test_mean_2, np.ones(10) * test_mean_2 * 3)

    for a in [plc.two_d_distribution(distr_1, distr_2),
              plc.two_d_distribution(distr_1, distr_2, samples=10000),
              plc.two_d_distribution(outlier_distr_1, outlier_distr_2, mad_filter=5)
              ]:
        assert round(a[0].flatten()[np.argmax(a[2])]) in [9, 10, 11]
        assert round(a[1].flatten()[np.argmax(a[2])]) in [9, 10, 11]
        assert round(5 * (a[0][0][1]-a[0][0][0])) == test_std_1
        assert round(5 * (a[1][1][0]-a[1][0][0])) == test_std_2

    b = plc.two_d_distribution(narrow_distr_1, narrow_distr_2)
    assert round(b[0].flatten()[np.argmax(b[2])]) == test_mean_1
    assert round(b[1].flatten()[np.argmax(b[2])]) == test_mean_2

    c = plc.two_d_distribution(distr_1, distr_2,
                               abs_step_x=1.6, min_value_x=5.0, max_value_x=15.0,
                               abs_step_y=1.6, min_value_y=5.0, max_value_y=15.0)

    assert round(c[0].flatten()[np.argmax(c[2])], 1) == 10.6
    assert round(c[1].flatten()[np.argmax(c[2])], 1) == 10.6
    assert round(c[0][0][1]-c[0][0][0], 1) == 1.6
    assert round(c[1][1][0]-c[1][0][0], 1) == 1.6
    assert round(c[0][0][0], 1) == 5.8
    assert round(c[1][0][0], 1) == 5.8
    assert round(c[0][-1][-1], 1) == 15.4
    assert round(c[1][-1][-1], 1) == 15.4
