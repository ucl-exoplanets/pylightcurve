
import pytest
import numpy as np
import pylightcurve as plc


def test_distribution_one_d():

    test_mean = 10.0
    test_std = 1.0

    distr = np.random.normal(test_mean, test_std, 100000)

    a = plc.one_d_distribution(distr)
    assert round(a[0][np.argmax(a[1])]) == test_mean
    assert round(5 * (a[0][1]-a[0][0])) == test_std

    b = plc.one_d_distribution(distr, step=10)
    assert round(b[0][np.argmax(b[1])]) == test_mean
    assert round(10 * (b[0][1]-b[0][0])) == test_std

