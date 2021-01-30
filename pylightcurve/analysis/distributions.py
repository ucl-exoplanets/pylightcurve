
__all__ = ['one_d_distribution', 'two_d_distribution']

import numpy as np

from pylightcurve.errors import *
from pylightcurve.analysis.gaussian import fit_gaussian, fit_two_d_gaussian


def one_d_distribution(datax, step=None, abs_step=None, min_value=None, max_value=None, mad_filter=None,
                       confidence_interval=None, gaussian_fit=False):

    datax = 1.0 * datax.flatten()

    if min_value is None and max_value is None:
        if mad_filter:
            data_median = np.median(datax)
            data_mad = np.sqrt(np.median((datax - data_median) ** 2))
            min_value = data_median - mad_filter * data_mad
            max_value = data_median + mad_filter * data_mad
        else:
            min_value = np.min(datax)
            max_value = np.max(datax)
    else:
        if min_value is None:
            min_value = np.min(datax)
        if max_value is None:
            max_value = np.max(datax)

    datax = datax[np.where(datax < max_value)]
    datax = datax[np.where(datax > min_value)]

    if abs_step is not None:
        step = abs_step
    else:
        if step is None:
            step = 5.0
        step = np.sqrt(np.median((datax - np.median(datax)) ** 2)) / step

    bins_number = int((max_value - min_value) / step) + 1
    bins = np.array(min_value + step / 2. + np.arange(bins_number) * step)

    counts = np.bincount(np.int_((datax - min_value) / step))
    counts = np.insert(counts, len(counts), np.zeros(int(bins_number) - len(counts)))

    if confidence_interval is None and not gaussian_fit:

        return bins, counts

    elif confidence_interval is not None and not gaussian_fit:

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

        return bins, counts, value, m_error, p_error

    elif confidence_interval is None and gaussian_fit:

        popt, pcov = fit_gaussian(bins, counts, positive=True, sampled=True)

        value, p_error, m_error = popt[2], popt[3], popt[3]

        return bins, counts, value, m_error, p_error

    else:

        raise PyLCProcessError('For the uncertainties in your distribution you should use either a confidence '
                               'interval or a gaussian fit, not both.')


def two_d_distribution(datax, datay, step=None, abs_step_x=None, abs_step_y=None, min_value_x=None, min_value_y=None,
                       max_value_x=None, max_value_y=None, mad_filter=None, confidence_interval=None,
                       gaussian_fit=False):

    datax = datax.flatten()
    datay = datay.flatten()

    if min_value_x is None and max_value_x is None:
        if mad_filter:
            data_median = np.median(datax)
            data_mad = np.sqrt(np.median((datax - data_median) ** 2))
            min_value_x = data_median - mad_filter * data_mad
            max_value_x = data_median + mad_filter * data_mad
        else:
            min_value_x = np.min(datax)
            max_value_x = np.max(datax)
    else:
        if min_value_x is None:
            min_value_x = np.min(datax)
        if max_value_x is None:
            max_value_x = np.max(datax)

    if min_value_y is None and max_value_y is None:
        if mad_filter:
            data_median = np.median(datay)
            data_mad = np.sqrt(np.median((datay - data_median) ** 2))
            min_value_y = data_median - mad_filter * data_mad
            max_value_y = data_median + mad_filter * data_mad
        else:
            min_value_y = np.min(datay)
            max_value_y = np.max(datay)
    else:
        if min_value_y is None:
            min_value_y = np.min(datay)
        if max_value_y is None:
            max_value_y = np.max(datay)

    values_to_keep = np.where((datax < max_value_x) * (datax > min_value_x) *
                              (datay < max_value_y) * (datay > min_value_y))

    datax = datax[values_to_keep]
    datay = datay[values_to_keep]

    if abs_step_x is not None:
        step_x = abs_step_x
    else:
        if step is None:
            step_x = 5.0
        else:
            step_x = step
        step_x = np.sqrt(np.median((datax - np.median(datax)) ** 2)) / step_x

    bins_number_x = int((max_value_x - min_value_x) / step_x) + 1
    bins_x = np.array(min_value_x + step_x / 2. + np.arange(bins_number_x) * step_x)

    if abs_step_y is not None:
        step_y = abs_step_y
    else:
        if step is None:
            step_y = 5.0
        else:
            step_y = step
        step_y = np.sqrt(np.median((datay - np.median(datay)) ** 2)) / step_y

    bins_number_y = int((max_value_y - min_value_y) / step_y) + 1
    bins_y = np.array(min_value_y + step_y / 2. + np.arange(bins_number_y) * step_y)

    bins_number_xy = bins_number_x * bins_number_y
    bins_xy = np.meshgrid(bins_x, bins_y)
    counts_xy = np.bincount((np.int_((datay - min_value_y) / step_y)) * bins_number_x +
                            (np.int_((datax - min_value_x) / step_x)))
    counts_xy = np.insert(counts_xy, len(counts_xy), np.zeros(bins_number_xy - len(counts_xy)))
    counts_xy = np.reshape(counts_xy, (bins_number_y, bins_number_x))

    if confidence_interval is None and not gaussian_fit:

        return bins_xy[0], bins_xy[1], counts_xy

    elif confidence_interval is not None and not gaussian_fit:

        value_x, m_error_x, p_error_x = one_d_distribution(
            datax, step=step, abs_step=abs_step_x, min_value=min_value_x, max_value=max_value_x,
            confidence_interval=confidence_interval)[2:]

        value_y, m_error_y, p_error_y = one_d_distribution(
            datay, step=step, abs_step=abs_step_y, min_value=min_value_y, max_value=max_value_y,
            confidence_interval=confidence_interval)[2:]

        return bins_xy[0], bins_xy[1], counts_xy, value_x, m_error_x, p_error_x, value_y, m_error_y, p_error_y

    elif confidence_interval is None and gaussian_fit:

        popt, pcov = fit_two_d_gaussian(bins_xy[0], bins_xy[1], counts_xy, positive=True)

        value_x, p_error_x, m_error_x = popt[2], popt[4], popt[4]
        value_y, p_error_y, m_error_y = popt[3], popt[5], popt[5]

        return bins_xy[0], bins_xy[1], counts_xy, value_x, m_error_x, p_error_x, value_y, m_error_y, p_error_y

    else:

        raise PyLCProcessError('For the uncertainties in your distribution you should use either a confidence '
                               'interval or a gaussian fit, not both.')
