
__all__ = ['waverage', 'mad', 'mean_std_from_median_mad', 'correlation',
           'values_to_print', 'residual_statistics', 'get_data_errorbars']

import numpy as np
import warnings

import scipy


def waverage(data, uncertainties, axis=None):
    """Calculates the weighted average using the invert square of the uncertainties as weights.

    :param data: Measurements.
    :type data: numpy.ndarray
    :param uncertainties: Measurements uncertainties.
    :type uncertainties: numpy.ndarray
    :param axis: Axis to sum
    :type axis: int
    :return: Weighted average of the measurements
    :rtype: numpy.ndarray
    """
    weights = 1.0 / (uncertainties ** 2)

    weighted_average = np.sum(data * weights, axis=axis) / np.sum(weights, axis=axis)
    weighted_average_uncertainty = 1.0 / np.sqrt((np.sum(weights, axis=axis)))

    return weighted_average, weighted_average_uncertainty


def mad(data):

    data = 1.0 * data.flatten()

    return np.sqrt(np.median((data - np.median(data)) ** 2))


def mean_std_from_median_mad(data, samples=None):

    data = 1.0 * data.flatten()

    if samples:
        data = data[::int(len(data) / samples + 1)]

    median = np.median(data)
    dx = data - median
    mad = np.sqrt(np.median(dx * dx))

    return median, mad * 0.6745


def get_data_errorbars(data, window=3):
    sigma = np.ones_like(data)
    for i in range(len(sigma)):
        if i < window:
            sigma[i] = mean_std_from_median_mad(data[:i + window])[1]
        elif i < len(sigma) - window:
            sigma[i] = mean_std_from_median_mad(data[i - window:i + window])[1]
        else:
            sigma[i] = mean_std_from_median_mad(data[i - window:])[1]

    return sigma


def correlation(datax, datay):

    datax = np.array(datax)
    datay = np.array(datay)
    correlation_xy = (np.sum((datax - np.mean(datax)) * (datay - np.mean(datay))) /
                      ((len(datax) - 1) * np.std(datax) * np.std(datay)))

    return correlation_xy


# decimal points and rounding


def values_to_print(value, error_minus, error_plus):

    value = float(value)
    error_minus = float(error_minus)
    error_plus = float(error_plus)

    if np.isinf(error_minus) or np.isinf(error_plus):
        print_value = str(round(value, 2))
        print_m_error = 'NaN'
        print_p_error = 'NaN'

    else:

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            width = min(int("{:.1e}".format(error_minus).split('e')[1]), int("{:.1e}".format(error_plus).split('e')[1]))
            width -= 1
            width *= -1

        print_value = str(round(value, width))
        print_m_error = str(round(error_minus, width))
        print_p_error = str(round(error_plus, width))

    return print_value, print_m_error, print_p_error


def gaussian(x_array, model_norm, model_floor, model_mean, model_std):

    return model_floor + (model_norm *
                          np.exp(- 0.5 * (model_mean - x_array) * (model_mean - x_array) / (model_std * model_std)))


def residual_statistics(datax, datay, datay_unc, model, number_of_free_parameters):

    residuals = datay - model

    norm_residuals = residuals / datay_unc
    norm_residuals = np.swapaxes([datax, norm_residuals], 0, 1)
    norm_residuals = sorted(norm_residuals, key=lambda x: x[0])
    norm_residuals = np.swapaxes(norm_residuals, 0, 1)[1]

    res_autocorr = np.correlate(norm_residuals, norm_residuals, mode='full')
    res_autocorr = res_autocorr[res_autocorr.size // 2:]
    res_autocorr /= res_autocorr[0]

    log_len = np.log10(len(norm_residuals))

    res_autocorr_model_mean = np.poly1d([0.01337266, -0.14360713,  0.09546788, -0.3990635])
    res_autocorr_model_std = np.poly1d([0.00023335,  0.00171671, -0.03351896,  0.13821703])

    res_max_autocorr_flag_sigma = (np.log10(np.max(np.abs(res_autocorr[1:]))) - res_autocorr_model_mean(log_len)) / res_autocorr_model_std(log_len)

    res_shapiro = scipy.stats.shapiro(norm_residuals)

    res_shapiro_flag_sigma = scipy.stats.norm.ppf(1 - res_shapiro.pvalue)

    statistics = {
        'res_autocorr': res_autocorr,
        'res_max_autocorr': np.max(np.abs(res_autocorr[1:])),
        'res_max_autocorr_flag': res_max_autocorr_flag_sigma > 3,
        'res_max_autocorr_flag_sigma': res_max_autocorr_flag_sigma,
        'res_shapiro': res_shapiro.statistic,
        'res_shapiro_flag': res_shapiro_flag_sigma > 3,
        'res_shapiro_flag_sigma': res_shapiro_flag_sigma,
        'res_mean': np.mean(residuals),
        'res_std': np.std(residuals),
        'res_rms': np.sqrt(np.mean(residuals**2)),
        'res_chi_sqr': np.sum(norm_residuals ** 2),
        'res_red_chi_sqr': np.sum(norm_residuals ** 2) / (len(norm_residuals) - number_of_free_parameters)
    }

    return statistics
