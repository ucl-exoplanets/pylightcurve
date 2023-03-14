
import numpy as np
import pylightcurve as plc


def test_stats():

    xx = np.arange(0, 10, 0.1)
    yy = 0.1 * xx + 5 + np.random.normal(0, 0.01, len(xx))
    zz = np.ones_like(xx) * 0.01
    model = 0.1 * xx + 5

    assert round(plc.waverage(yy - model, zz)[0], 2) == 0.0
    assert round(plc.mad(yy - model), 2) == 0.01
    assert round(plc.mean_std_from_median_mad(yy - model)[0], 2) <= 0.1
    assert round(plc.mean_std_from_median_mad(yy - model, samples=50)[0], 2) <= 0.1

    assert round(plc.correlation(xx, yy), 1) == 1

    residual_statistics = plc.residual_statistics(xx, yy, zz, model, 2)
    assert residual_statistics['res_max_autocorr'] < 1
    assert residual_statistics['res_shapiro'] < 1
    assert residual_statistics['res_mean'] < 1
    assert residual_statistics['res_std'] < 1
    assert residual_statistics['res_rms'] < 1
    assert not residual_statistics['res_max_autocorr_flag']
    assert not residual_statistics['res_shapiro_flag']

    assert (np.mean(plc.get_data_errorbars(yy - model))) < 1

    assert plc.values_to_print(3.1234567, 0.12456, 0.4345) == ('3.12', '0.12', '0.43')
    assert plc.values_to_print(3.1234567, 0.4345, 0.12456) == ('3.12', '0.43', '0.12')
    assert plc.values_to_print(3.1234567, 0.456, 0.657) == ('3.12', '0.46', '0.66')
    assert plc.values_to_print(3.1234567, 3.3, 2.65789) == ('3.1', '3.3', '2.7')
    assert plc.values_to_print(3.1234567, 3.0, 2.0) == ('3.1', '3.0', '2.0')
    assert plc.values_to_print(31.234567, 32.4, 22.1) == ('31.0', '32.0', '22.0')
    assert plc.values_to_print(31.234567, np.inf, 22.1) == ('31.23', 'NaN', 'NaN')