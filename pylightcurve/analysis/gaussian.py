
__all__ = ['gaussian', 'fit_gaussian', 'two_d_gaussian', 'fit_two_d_gaussian']


import numpy as np
import warnings

from scipy.optimize import curve_fit

from .numerical_integration import sample_function
from .statistics import get_data_errorbars, mean_std_from_median_mad

# gaussian 1D

def gaussian(x_array, model_norm, model_floor, model_mean, model_std):

    return model_floor + (model_norm *
                          np.exp(- 0.5 * (model_mean - x_array) * (model_mean - x_array) / (model_std * model_std)))


def fit_gaussian(datax, datay, errors=None, positive=False, sampled=False, sampled_precision=3, maxfev=10000,
                 point_x=None, norm=None, floor=None, sigma=None,
                 filter_outliers=True,
                 filter_outliers_window=3,
                 filter_outliers_std_limit=3
                 ):

    datax = np.array(datax, dtype=float)
    datay = np.array(datay, dtype=float)

    if floor is None:
        initial_floor = np.median(datay)
    else:
        initial_floor = floor

    if norm is None:
        initial_norm = 3 * mean_std_from_median_mad(datay-initial_floor)[1]
    else:
        initial_norm = norm

    if sigma is None:
        initial_sigma = 1
    else:
        initial_sigma = sigma

    if point_x is None:
        initial_x_mean_arg = int(len(datax)/2)
        model = gaussian(datax, initial_norm, initial_floor, datax[initial_x_mean_arg], initial_sigma)
        dx = np.argmax(np.convolve(datay, model)) - np.argmax(np.convolve(model, model))
        initial_mean = datax[initial_x_mean_arg + dx]
    else:
        initial_mean = point_x

    if norm is None:
        test_x = np.clip(int(initial_mean), np.min(datax), np.max(datax))
        test_x = np.argmin(np.abs(datax - test_x))
        initial_norm = datay[test_x]

    if sigma is None:
        initial_sigma = max(min(datax[1:] - datax[:-1]),
                            np.abs(datax[np.argmin(np.abs((datay - initial_norm - initial_floor) / 2))] - initial_mean))

    initial_values = [initial_norm, initial_floor, initial_mean, initial_sigma]

    if positive:
        def gaussian_to_fit(x_array, model_norm, model_floor, model_mean, model_std):
            return gaussian(x_array, np.abs(model_norm), model_floor, model_mean, model_std)
    else:
        gaussian_to_fit = gaussian

    if sampled:
        splits = 0.5 * (datax[1:] + datax[:-1])
        datax1 = np.array([datax[0] - (splits[0] - datax[0])] + list(splits))
        datax2 = np.array(list(splits) + [datax[-1] + (datax[-1] - splits[-1])])

        datax = (datax1, datax2)
        gaussian_to_fit = sample_function(gaussian_to_fit, precision=sampled_precision)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if errors is not None:
            popt, pcov = curve_fit(gaussian_to_fit, datax, datay, p0=initial_values, sigma=errors, maxfev=maxfev)
        else:
            popt, pcov = curve_fit(gaussian_to_fit, datax, datay, p0=initial_values, maxfev=maxfev)

    if filter_outliers:

        if errors is not None:
            norm_res = (gaussian_to_fit(datax, *popt) - datay)/errors
        else:
            norm_res = (gaussian_to_fit(datax, *popt) - datay)/get_data_errorbars(datay, window=filter_outliers_window)

        non_outliers = np.where(np.abs(norm_res) < filter_outliers_std_limit * np.std(norm_res))

        if type(datax) == tuple:
            datax_clean = (datax[0][non_outliers], datax[1][non_outliers])
        else:
            datax_clean = datax[non_outliers]

        if errors is not None:
            popt, pcov = curve_fit(gaussian_to_fit, datax_clean, datay[non_outliers], p0=initial_values, sigma=errors[non_outliers], maxfev=maxfev)
        else:
            popt, pcov = curve_fit(gaussian_to_fit, datax_clean, datay[non_outliers], p0=initial_values, maxfev=maxfev)

    popt[3] = np.abs(popt[3])
    if positive:
        popt[0] = np.abs(popt[0])

    return popt, pcov


# gaussian 2D

def two_d_gaussian(x_array, y_array, model_norm, model_floor, model_x_mean, model_y_mean, model_x_sigma, model_y_sigma,
                   model_theta, extend=2):

    xt_array = x_array - model_x_mean
    yt_array = y_array - model_y_mean

    xt_array, yt_array = ((xt_array * np.cos(model_theta) + yt_array * np.sin(model_theta))/model_x_sigma,
                          (- xt_array * np.sin(model_theta) + yt_array * np.cos(model_theta))/model_y_sigma)

    return model_floor + model_norm * np.exp(-0.5*(yt_array * yt_array + xt_array**extend))


def fit_two_d_gaussian(datax, datay, dataz, errors=None, positive=False, symmetric=False,
                       point_xy=None, norm=None, floor=None, sigma=None, theta=None, maxfev=10000, extend=2):

    datax = np.array(datax, dtype=float)
    datay = np.array(datay, dtype=float)
    dataz = np.array(dataz, dtype=float)

    if floor is None:
        initial_floor = np.median(dataz)
    else:
        initial_floor = floor

    if norm is None:
        initial_norm = 3 * mean_std_from_median_mad(dataz-initial_floor)[1]
    else:
        initial_norm = norm

    if sigma is None:
        initial_x_sigma = 1.0
        initial_y_sigma = 1.0
    else:
        initial_x_sigma = sigma
        initial_y_sigma = sigma

    if theta is None:
        initial_theta = 0
    else:
        initial_theta = theta

    if point_xy is None:
        test_data_x = np.sum(dataz, 0)
        test_floor = np.median(test_data_x)
        test_norm = np.max(test_data_x) - test_floor

        initial_x_mean_arg = int(len(datax[0]) / 2)
        model = gaussian(datax[0], test_norm, test_floor, datax[0][initial_x_mean_arg], initial_x_sigma)
        dx = np.argmax(np.convolve(test_data_x, model)) - np.argmax(np.convolve(model, model))
        initial_x_mean = datax[0][initial_x_mean_arg + dx]

        test_data_y = np.sum(dataz, 1)
        test_floor = np.median(test_data_y)
        test_norm = np.max(test_data_y) - test_floor

        initial_y_mean_arg = int(len(datay[:, 0]) / 2)
        model = gaussian(datay[:, 0], test_norm, test_floor, datay[:, 0][initial_y_mean_arg], initial_y_sigma)
        dy = np.argmax(np.convolve(np.sum(dataz, 1), model)) - np.argmax(np.convolve(model, model))
        initial_y_mean = datay[:, 0][initial_y_mean_arg + dy]

    else:
        initial_x_mean, initial_y_mean = np.int_(point_xy)

    if norm is None:
        test_x = np.clip(int(initial_x_mean), np.min(datax), np.max(datax))
        test_y = np.clip(int(initial_y_mean), np.min(datay), np.max(datay))
        test_x = np.argmin(np.abs(datax[0] - test_x))
        test_y = np.argmin(np.abs(datay[:, 0] - test_y))
        initial_norm = dataz[test_y][test_x]

    if sigma is None:
        initial_x_sigma = np.abs(datax[0][np.argmin(np.abs(np.sum(dataz, 0) - np.max(np.sum(dataz, 1)) / 2))] -
                                 initial_x_mean)

        initial_y_sigma = np.abs(datay[:, 0][np.argmin(np.abs(np.sum(dataz, 1) - np.max(np.sum(dataz, 0)) / 2))] -
                                 initial_y_mean)

    if positive and symmetric:

        def gaussian_to_fit(xy_array, model_norm, model_floor, model_x_mean, model_y_mean, model_sigma):
            return two_d_gaussian(xy_array[0], xy_array[1], model_norm, model_floor, model_x_mean,
                                  model_y_mean, model_sigma, model_sigma, 0, extend=extend)

        initial_values = [initial_norm, initial_floor, initial_x_mean, initial_y_mean, initial_x_sigma]
        bounds_1 = [0, -np.inf, -np.inf, -np.inf, 0]
        bounds_2 = [np.inf, np.inf, np.inf, np.inf, np.inf]

    elif positive:

        def gaussian_to_fit(xy_array, model_norm, model_floor, model_x_mean, model_y_mean, model_x_sigma, model_y_sigma,
                            model_theta):
            return two_d_gaussian(xy_array[0], xy_array[1], model_norm, model_floor, model_x_mean,
                                  model_y_mean, model_x_sigma, model_y_sigma, model_theta, extend=extend)

        initial_values = [initial_norm, initial_floor, initial_x_mean, initial_y_mean,
                          initial_x_sigma, initial_y_sigma, initial_theta]
        bounds_1 = [0, -np.inf, -np.inf, -np.inf, 0, 0, -np.pi/2]
        bounds_2 = [np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.pi/2]

    elif symmetric:

        def gaussian_to_fit(xy_array, model_norm, model_floor, model_x_mean, model_y_mean, model_sigma):
            return two_d_gaussian(xy_array[0], xy_array[1], model_norm, model_floor, model_x_mean,
                                  model_y_mean, model_sigma, model_sigma, 0, extend=extend)

        initial_values = [initial_norm, initial_floor, initial_x_mean, initial_y_mean, initial_x_sigma]
        bounds_1 = [-np.inf, -np.inf, -np.inf, -np.inf, 0]
        bounds_2 = [np.inf, np.inf, np.inf, np.inf, np.inf]

    else:

        def gaussian_to_fit(xy_array, model_norm, model_floor, model_x_mean, model_y_mean, model_x_sigma, model_y_sigma,
                            model_theta):
            return two_d_gaussian(xy_array[0], xy_array[1], model_norm, model_floor, model_x_mean, model_y_mean,
                                  model_x_sigma, model_y_sigma, model_theta, extend=extend)

        initial_values = [initial_norm, initial_floor, initial_x_mean, initial_y_mean,
                          initial_x_sigma, initial_y_sigma, initial_theta]
        bounds_1 = [-np.inf, -np.inf, -np.inf, -np.inf, 0, 0, -np.pi/2]
        bounds_2 = [np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.pi/2]

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if errors is not None:
            popt, pcov = curve_fit(gaussian_to_fit, (datax.flatten(), datay.flatten()), dataz.flatten(),
                                   p0=initial_values, sigma=errors.flatten(), maxfev=maxfev,
                                   bounds=(np.array(bounds_1), np.array(bounds_2)))
        else:
            popt, pcov = curve_fit(gaussian_to_fit, (datax.flatten(), datay.flatten()), dataz.flatten(),
                                   p0=initial_values, maxfev=maxfev,
                                   bounds=(np.array(bounds_1), np.array(bounds_2)))

    if not symmetric:

        if popt[5] > popt[4]:
            popt[4], popt[5] = popt[5], popt[4]
            pcov[4][4], pcov[5][5] = pcov[5][5], pcov[4][4]
            popt[6] = -popt[6]

    return popt, pcov
