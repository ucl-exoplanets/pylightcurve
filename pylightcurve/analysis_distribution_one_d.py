from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from .analysis_gauss_numerical_integration import *
from .analysis_optimisation import *


def one_d_distribution(data_array, step=None, min_value=None, max_value=None, confidence_interval=None):

    try:
        data_array = np.array(data_array)
        data_array = data_array.flatten()
    except AttributeError:
        raise RuntimeError('Invalid input data')

    if not step:
        step = np.sqrt(np.median((data_array - np.median(data_array)) ** 2)) / 5.0
    else:
        try:
            step = np.sqrt(np.median((data_array - np.median(data_array)) ** 2)) / step
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
    counts = np.insert(counts, len(counts), np.zeros(int(bins_number) - len(counts)))

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

        value, p_error, m_error = exp_val, rightci - exp_val, exp_val - leftci

        print_value, print_m_error, print_p_error = values_to_print(exp_val, m_error, p_error)

        return bins, counts, value, m_error, p_error, print_value, print_m_error, print_p_error


def gaussian(x_array, model_norm, model_floor, model_mean, model_std):

    return model_floor + (model_norm *
                          np.exp(- 0.5 * (model_mean - x_array) * (model_mean - x_array) / (model_std * model_std)))


def gaussian_positive(x_array, model_norm, model_floor, model_mean, model_std):

    return gaussian(x_array, np.abs(model_norm), np.abs(model_floor), model_mean, model_std)


def fit_gaussian(datax, datay, positive=False, plot=False, maxfev=10000):

    if positive:
        sgaussian = gaussian_positive
    else:
        sgaussian = gaussian

    mean = datax[np.argmax(datay)]
    sigma = np.abs(datax[np.argmin(np.abs(datay - np.max(datay) / 2))] - mean)
    norm = (np.max(datay) - np.min(datay)) / gaussian(mean, mean, sigma, 1.0, 0.0)
    floor = np.min(datay)

    popt, pcov = curve_fit(sgaussian, datax, datay, p0=[norm, floor, mean, sigma], maxfev=maxfev)

    if plot:
        plt.plot(datax, datay, 'bo')
        plotx = np.arange(min(datax), max(datax), (max(datax) - min(datax)) / len(datax) / 100.0)
        plt.plot(plotx, gaussian(plotx, *popt), 'r-')

    return popt, pcov


def fit_sampled_gaussian(datax1, datax2, datay, positive=False, precision=3, plot=False):

    if positive:
        sgaussian = sample_function(gaussian, precision=precision)
    else:
        sgaussian = sample_function(gaussian_positive, precision=precision)

    datax = np.array(0.5 * (datax1 + datax2))
    binx = np.median(datax2 - datax1)

    mean = datax[np.argmax(datay)]
    sigma = np.abs(datax[np.argmin(np.abs(datay - np.max(datay) / 2))] - mean)
    norm = np.max(datay) - np.min(datay)
    floor = np.min(datay)

    popt, pcov = curve_fit(sgaussian, (datax1, datax2), datay, p0=[norm, floor, mean, sigma])

    if plot:
        plt.plot(datax, datay, 'bo')
        plotx = np.arange(min(datax), max(datax), (max(datax) - min(datax)) / len(datax) / 100.0)
        plt.plot(plotx, sgaussian((plotx - binx / 2.0, plotx + binx / 2.0), *popt), 'r-')

    return popt, pcov


def fit_one_d_distribution_gaussian(datax, step=None, min_value=None, plot=False):

    datax, datay = one_d_distribution(datax, step=step, min_value=min_value)

    popt, pcov = fit_gaussian(datax, datay, positive=True, plot=plot)

    return datax, datay, popt, pcov


def fit_one_d_distribution_sampled_gaussian(datax, step=None, min_value=None, precision=3, plot=False):

    datax1, datax2, datay = one_d_distribution(datax, step=step, min_value=min_value)

    popt, pcov = fit_sampled_gaussian(datax1, datax2, datay, positive=True, precision=precision, plot=plot)

    return datax1, datax2, datay, popt, pcov


def lorentzian(x_array, mean, gamma, norm, floor):

    return ((1.0 / abs(gamma) / np.pi) *
            (gamma * gamma / ((x_array - mean) * (x_array - mean) + gamma * gamma))) * norm + floor


def fit_lorentzian(datax, datay, plot=False):

    mean = datax[np.argmax(datay)]
    sigma = np.abs(datax[np.argmin(np.abs(datay - np.max(datay) / 2))] - mean)
    norm = (np.max(datay) - np.min(datay)) / lorentzian(mean, mean, sigma, 1.0, 0.0)
    floor = np.min(datay)

    popt, pcov = curve_fit(lorentzian, datax, datay, p0=[mean, sigma, norm, floor])

    if plot:
        plt.plot(datax, datay, 'bo')
        plotx = np.arange(min(datax), max(datax), (max(datax) - min(datax)) / len(datax) / 100.0)
        plt.plot(plotx, lorentzian(plotx, *popt), 'r-')

    return popt, pcov


def fit_sampled_lorentzian(datax1, datax2, datay, precision=3, plot=False):

    slorentzian = sample_function(lorentzian, precision=precision)
    datax = np.array(0.5 * (datax1 + datax2))
    binx = np.median(datax2 - datax1)

    mean = datax[np.argmax(datay)]
    sigma = np.abs(datax[np.argmin(np.abs(datay - np.max(datay) / 2))] - mean)
    norm = (np.max(datay) - np.min(datay)) / lorentzian(mean, mean, sigma, 1.0, 0.0)
    floor = np.min(datay)

    popt, pcov = curve_fit(slorentzian, (datax1, datax2), datay, p0=[mean, sigma, norm, floor])

    if plot:
        plt.plot(datax, datay, 'bo')
        plotx = np.arange(min(datax), max(datax), (max(datax) - min(datax)) / len(datax) / 100.0)
        plt.plot(plotx, slorentzian((plotx - binx / 2.0, plotx + binx / 2.0), *popt), 'r-')

    return popt, pcov


def fit_one_d_distribution_lorentzian(datax, xstep=None, xmin=None, precision=3, plot=False):

    datax1, datax2, datay = one_d_distribution(datax, xstep=xstep, xmin=xmin, plot=plot)

    popt, pcov = fit_sampled_lorentzian(datax1, datax2, datay, precision=precision, plot=plot)

    return datax1, datax2, datay, popt, pcov
