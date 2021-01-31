
import pytest
import numpy as np
import pylightcurve as plc


def test_distribution_one_d():

    test_mean = 10.0
    test_std = 1.0
    distr = np.append(np.random.normal(test_mean, test_std, 100000), 1000000)

    a = plc.one_d_distribution(distr, mad_filter=5)
    assert round(a[0][np.argmax(a[1])]) == test_mean
    assert round(5 * (a[0][1]-a[0][0])) == test_std

    distr = np.random.normal(test_mean, test_std, 100000)

    a = plc.one_d_distribution(distr)
    assert round(a[0][np.argmax(a[1])]) == test_mean
    assert round(5 * (a[0][1]-a[0][0])) == test_std

    b = plc.one_d_distribution(distr, step=10)
    assert round(b[0][np.argmax(b[1])]) == test_mean
    assert round(10 * (b[0][1]-b[0][0])) == test_std

    c = plc.one_d_distribution(distr, abs_step=1.6, min_value=5.0, max_value=15.0)
    assert round(c[0][np.argmax(c[1])], 1) == 10.6
    assert round(c[0][1]-c[0][0], 1) == 1.6
    assert round(c[0][0], 1) == 5.8
    assert round(c[0][-1], 1) == 15.4

    c = plc.one_d_distribution(distr, abs_step=1.6, min_value=5.0)
    assert round(c[0][0], 1) == 5.8

    c = plc.one_d_distribution(distr, abs_step=1.6, max_value=15.0)
    assert round(c[0][-1], 1) < 16.6

    a = plc.one_d_distribution(distr, confidence_interval=0.68)
    assert len(a) == 5
    assert round(a[0][np.argmax(a[1])]) == test_mean
    assert round(5 * (a[0][1] - a[0][0])) == test_std
    assert round(a[2]) == test_mean
    assert round(a[3]) == test_std
    assert round(a[4]) == test_std
    a = plc.one_d_distribution(distr, gaussian_fit=True)
    assert len(a) == 5
    assert round(a[0][np.argmax(a[1])]) == test_mean
    assert round(5 * (a[0][1] - a[0][0])) == test_std
    assert round(a[2]) == test_mean
    assert round(a[3]) == test_std
    assert round(a[4]) == test_std

    with pytest.raises(plc.PyLCProcessError):
        plc.one_d_distribution(distr, confidence_interval=0.68, gaussian_fit=True)


def test_distribution_two_d():

    test_mean_1 = 10.0
    test_std_1 = 1.0
    test_mean_2 = 20.0
    test_std_2 = 2.0

    distr_1 = np.append(np.random.normal(test_mean_1, test_std_1, 100000), 100000)
    distr_2 = np.append(np.random.normal(test_mean_2, test_std_2, 100000), 100000)

    a = plc.two_d_distribution(distr_1, distr_2, mad_filter=5)
    assert abs(a[0].flatten()[np.argmax(a[2])] - test_mean_1) < test_std_1
    assert abs(a[1].flatten()[np.argmax(a[2])] - test_mean_2) < test_std_2
    assert round(5 * (a[0][0][1]-a[0][0][0]), 1) == 0.7
    assert round(5 * (a[1][1][0]-a[1][0][0]), 1) in [1.3, 1.4]

    distr_1 = np.random.normal(test_mean_1, test_std_1, 100000)
    distr_2 = np.random.normal(test_mean_2, test_std_2, 100000)

    a = plc.two_d_distribution(distr_1, distr_2)
    assert abs(a[0].flatten()[np.argmax(a[2])] - test_mean_1) < test_std_1
    assert abs(a[1].flatten()[np.argmax(a[2])] - test_mean_2) < test_std_2
    assert round(5 * (a[0][0][1]-a[0][0][0]), 1) == 0.7
    assert round(5 * (a[1][1][0]-a[1][0][0]), 1) in [1.3, 1.4]

    c = plc.two_d_distribution(distr_1, distr_2, step=1, min_value_x=5.0, min_value_y=15.0)
    assert round(np.min(c[0]), 1) >= 5.0
    assert round(np.min(c[1]), 1) >= 15.0

    c = plc.two_d_distribution(distr_1, distr_2, step=1, max_value_x=15.0, max_value_y=25.0)
    assert round(np.max(c[0]), 1) <= 15.5
    assert round(np.max(c[1]), 1) <= 25.7

    c = plc.two_d_distribution(distr_1, distr_2, step=1, min_value_x=5.0, max_value_x=15.0, min_value_y=15.0,
                               max_value_y=25.0)
    assert round(np.min(c[0]), 1) >= 5.0
    assert round(np.max(c[0]), 1) <= 15.6
    assert round(np.min(c[1]), 1) >= 15.0
    assert round(np.max(c[1]), 1) <= 25.7

    a = plc.two_d_distribution(distr_1, distr_2, step=10)
    assert abs(a[0].flatten()[np.argmax(a[2])] - test_mean_1) < test_std_1
    assert abs(a[1].flatten()[np.argmax(a[2])] - test_mean_2) < test_std_2
    assert round(10 * (a[0][0][1]-a[0][0][0]), 1) == 0.7
    assert round(10 * (a[1][1][0]-a[1][0][0]), 1) in [1.3, 1.4]

    a = plc.two_d_distribution(distr_1, distr_2, abs_step_x=0.1, abs_step_y=0.1)
    assert abs(a[0].flatten()[np.argmax(a[2])] - test_mean_1) < test_std_1
    assert abs(a[1].flatten()[np.argmax(a[2])] - test_mean_2) < test_std_2
    assert round(a[0][0][1]-a[0][0][0], 1) == 0.1
    assert round(a[1][1][0]-a[1][0][0], 1) == 0.1

    a = plc.two_d_distribution(distr_1, distr_2, confidence_interval=0.68)
    assert len(a) == 9
    assert abs(a[0].flatten()[np.argmax(a[2])] - test_mean_1) < test_std_1
    assert abs(a[1].flatten()[np.argmax(a[2])] - test_mean_2) < test_std_2
    assert round(5 * (a[0][0][1]-a[0][0][0]), 1) == 0.7
    assert round(5 * (a[1][1][0]-a[1][0][0]), 1) in [1.3, 1.4]
    assert round(a[3]) == test_mean_1
    assert round(a[4]) == test_std_1
    assert round(a[5]) == test_std_1
    assert round(a[6]) == test_mean_2
    assert round(a[7]) == test_std_2
    assert round(a[8]) == test_std_2

    a = plc.two_d_distribution(distr_1, distr_2, gaussian_fit=True)
    assert len(a) == 9
    assert abs(a[0].flatten()[np.argmax(a[2])] - test_mean_1) < test_std_1
    assert abs(a[1].flatten()[np.argmax(a[2])] - test_mean_2) < test_std_2
    assert round(5 * (a[0][0][1]-a[0][0][0]), 1) == 0.7
    assert round(5 * (a[1][1][0]-a[1][0][0]), 1) in [1.3, 1.4]
    assert round(a[3]) == test_mean_1
    assert round(a[4]) == test_std_1
    assert round(a[5]) == test_std_1
    assert round(a[6]) == test_mean_2
    assert round(a[7]) == test_std_2
    assert round(a[8]) == test_std_2

    with pytest.raises(plc.PyLCProcessError):
        plc.two_d_distribution(distr_1, distr_2, confidence_interval=0.68, gaussian_fit=True)
