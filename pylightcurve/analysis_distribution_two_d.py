from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from .analysis_gauss_numerical_integration import *
from .analysis_distribution_one_d import *


def two_d_distribution(data_array, step=None, min_value=None, max_value=None, confidence_interval=None):

    try:
        data_array = data_array.flatten()
    except AttributeError:
        raise RuntimeError('Invalid input data')

    if not step:
        step = np.sqrt(np.median((data_array - np.median(data_array)) ** 2)) / 5.0
    else:
        try:
            step = float(step)
        except ValueError:
            raise RuntimeError('Invalid step value')

    if not min_value:
        min_value = min(data_array)
    else:
        try:
            min_value = float(min_value)
        except ValueError:
            raise RuntimeError('Invalid minimum value')

    if not max_value:
        max_value = max(data_array)
    else:
        try:
            min_value = float(min_value)
        except ValueError:
            raise RuntimeError('Invalid minimum value')

    bins_number = int((max_value - min_value) / step) + 1

    bins = np.array(min_value + step / 2. + np.arange(bins_number) * step)

    counts = np.bincount(np.int_((data_array - min_value) / step))
    counts = np.insert(counts, len(counts), np.zeros(bins_number - len(counts)))

    if confidence_interval is None:

        return bins, counts

    else:

        try:
            confidence_interval = float(confidence_interval)
        except ValueError:
            raise RuntimeError('Invalid confidence interval value')

        distrx = bins
        bin_width = bins[1] - bins[0]
        distr = counts

        # corresponds to the 1-sigma level probability

        pleft = 0.0
        centroid = np.argmax(distr)
        exp_val = distrx[centroid]

        total_probability_left = np.sum(bin_width * distr[:centroid]) * confidence_interval
        total_probability_right = np.sum(bin_width * distr[centroid:]) * confidence_interval

        num = centroid
        leftci = 0
        while pleft <= total_probability_left:
            if num == centroid:
                pleft += (bin_width / 2.0) * distr[num]
            else:
                pleft += bin_width * distr[num]
            leftci = distrx[num]
            num -= 1
            if num < 0:
                print('ERROR : confidence level can not be reached from left')
                break
        pright = 0.0
        num = centroid
        rightci = 0
        while pright <= total_probability_right:
            if num == centroid:
                pright += (bin_width / 2.0) * distr[num]
            else:
                pright += bin_width * distr[num]
            rightci = distrx[num]
            num += 1
            if num > len(distr) - 1:
                print('ERROR : confidence level can not be reached from right')
                break

        error_plus, error_minus = rightci - exp_val, exp_val - leftci

        try:
            digit1 = abs(int(np.log10(error_minus))) + 1
        except OverflowError:
            digit1 = 3
        try:
            digit2 = abs(int(np.log10(error_plus))) + 1
        except OverflowError:
            digit2 = 3
        try:
            done1 = 1 // int(error_minus * (10 ** digit1))
        except ZeroDivisionError:
            done1 = 0
        try:
            done2 = 1 // int(error_plus * (10 ** digit2))
        except ZeroDivisionError:
            done2 = 0

        if error_minus > 1. and error_plus > 1.:
            digit1, digit2, done1, done2 = 0, 0, 0, 0

        width = max(digit1 + done1, digit2 + done2)

        value = exp_val
        m_error = error_minus
        p_error = error_plus

        print_value = '{0:.{width}f}'.format(round(exp_val, width), width=width)
        print_m_error = '{0:.{width}f}'.format(round(error_minus, width), width=width)
        print_p_error = '{0:.{width}f}'.format(round(error_plus, width), width=width)

        return bins, counts, value, m_error, p_error, print_value, print_m_error, print_p_error


def two_d_gaussian(xy_array, model_norm, model_floor, model_x_mean, model_y_mean, model_x_std, model_y_std,
                   model_theta):

    x_array, y_array = xy_array

    a = (np.cos(model_theta) ** 2) / (2 * model_x_std ** 2) + (np.sin(model_theta) ** 2) / (2 * model_y_std ** 2)
    b = -(np.sin(2 * model_theta)) / (4 * model_x_std ** 2) + (np.sin(2 * model_theta)) / (4 * model_y_std ** 2)
    c = (np.sin(model_theta) ** 2) / (2 * model_x_std ** 2) + (np.cos(model_theta) ** 2) / (2 * model_y_std ** 2)

    return (model_floor + model_norm * np.exp(- (a * ((x_array - model_x_mean) ** 2)
                                              + 2.0 * b * (x_array - model_x_mean) * (y_array - model_y_mean)
                                              + c * ((y_array - model_y_mean) ** 2)))).flatten()


def two_d_gaussian_symmetric(xy_array, model_norm, model_floor, model_x_mean, model_y_mean, model_x_std):

    x_array, y_array = xy_array

    return (model_floor + model_norm * np.exp(-(((model_x_mean - x_array) ** 2) / (model_x_std ** 2)
                                              + ((model_y_mean - y_array) ** 2) / (model_x_std ** 2)) / 2)).flatten()


def two_d_gaussian_positive(xy_array, model_norm, model_floor, model_x_mean, model_y_mean, model_x_std, model_y_std,
                            model_theta):

    return two_d_gaussian(xy_array, np.abs(model_norm), np.abs(model_floor), model_x_mean, model_y_mean, model_x_std,
                          model_y_std, model_theta)


def two_d_gaussian_symmetric_positive(xy_array, model_norm, model_floor, model_x_mean, model_y_mean, model_x_std):

    return two_d_gaussian_symmetric(xy_array, np.abs(model_norm), np.abs(model_floor), model_x_mean, model_y_mean,
                                    model_x_std)


def fit_two_d_gaussian(dataxy, positive=False, symmetric=False, point_xy=None, window=None, sigma=None, maxfev=10000):

    if not point_xy:
        model = gaussian(np.arange(len(dataxy[0])), np.max(dataxy) - np.median(dataxy), np.median(dataxy), len(dataxy[0]) / 2, 1.0)
        dx = np.argmax(np.convolve(np.sum(dataxy, 0), model)) - np.argmax(np.convolve(model, model))
        model = gaussian(np.arange(len(dataxy)), np.max(dataxy) - np.median(dataxy), np.median(dataxy), len(dataxy) / 2, 1.0)
        dy = np.argmax(np.convolve(np.sum(dataxy, 1), model)) - np.argmax(np.convolve(model, model))

        x00 = int(dx + len(dataxy) / 2)
        y00 = int(dy + len(dataxy) / 2)
    else:
        x00, y00 = np.int_(point_xy)

    if window:
        window = np.int(window)
        min_y = int(max(y00 - window, 0))
        max_y = int(min(y00 + window, len(dataxy)))
        min_x = int(max(x00 - window, 0))
        max_x = int(min(x00 + window, len(dataxy[0])))
        dataxy = dataxy[min_y:max_y, min_x:max_x]
    else:
        min_y = 0
        max_y = len(dataxy)
        min_x = 0
        max_x = len(dataxy[0])

    if not sigma:
        sigma = np.abs(np.arange(len(dataxy[0]))[np.argmin(np.abs(np.sum(dataxy, 0)
                                                                  - np.max(np.sum(dataxy, 0)) / 2))] - x00)
    norm = np.max(dataxy) - np.median(dataxy)
    floor = np.median(dataxy)

    if symmetric:
        if positive:
            s_two_d_gaussian = two_d_gaussian_symmetric_positive
        else:
            s_two_d_gaussian = two_d_gaussian_symmetric
        pp0 = [norm, floor, x00, y00, sigma]
    else:
        if positive:
            s_two_d_gaussian = two_d_gaussian_positive
        else:
            s_two_d_gaussian = two_d_gaussian
        pp0 = [norm, floor, x00, y00,  sigma, sigma, 0]

    x = np.arange(min_x, max_x)
    y = np.arange(min_y, max_y)
    x, y = np.meshgrid(x, y)
    popt, pcov = curve_fit(s_two_d_gaussian, (x.flatten() + 0.5, y.flatten() + 0.5), dataxy.flatten(), p0=pp0,
                           maxfev=maxfev)

    return popt, pcov
