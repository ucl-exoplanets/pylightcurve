
__all__ = ['one_d_distribution', 'two_d_distribution']

import numpy as np

from .gaussian import fit_gaussian

from ..errors import *


def one_d_distribution(datax, step=5.0, abs_step=None,
                       min_value=None, max_value=None, mad_filter=None,
                       confidence_interval=None, gaussian_fit=False,
                       samples=None):

    datax = datax.flatten()
    int_test = np.sum((datax-np.int_(datax))**2) == 0

    if samples:
        datax = datax[::int(len(datax)/samples + 1)]

    if min_value is None:
        min_value = np.min(datax)
    else:
        datax = datax[np.where(datax >= min_value)]

    if max_value is None:
        max_value = np.max(datax)
    else:
        datax = datax[np.where(datax <= max_value)]

    if mad_filter:
        data_median = np.median(datax)
        data_mad = np.sqrt(np.median((datax - data_median) ** 2))
        min_value = data_median - mad_filter * data_mad
        max_value = data_median + mad_filter * data_mad
        datax = datax[np.where(datax <= max_value)]
        datax = datax[np.where(datax >= min_value)]

    if abs_step is None:
        new_step = np.sqrt(np.median((datax - np.median(datax)) ** 2)) / step
        if new_step == 0:
            new_step = np.sqrt(np.mean((datax - np.mean(datax)) ** 2)) / step
        step = new_step
    else:
        step = abs_step

    if int_test:
        step = max(1, round(step, 0))

    bins_number = int((max_value - min_value) / step) + 1
    bins = np.array(min_value + step / 2. + np.arange(bins_number) * step)

    counts = np.bincount(np.int_((datax - min_value) / step))
    counts = np.insert(counts, len(counts), np.zeros(int(bins_number) - len(counts)))

    if not confidence_interval and not gaussian_fit:

        return bins, counts

    if confidence_interval and not gaussian_fit:

        if confidence_interval is True:
            confidence_interval = 0.6827

        min_value, value, max_value = np.quantile(datax, [0.5 - confidence_interval / 2.0,
                                                          0.5,
                                                          0.5 + confidence_interval / 2.0])
        m_error = value - min_value
        p_error = max_value - value

        return bins, counts, value, m_error, p_error

    elif not confidence_interval and gaussian_fit:

        popt, pcov = fit_gaussian(bins, counts, errors=np.where(counts>0, np.sqrt(counts), 1),
                                  positive=True, sampled=True, floor=0)

        value, p_error, m_error = popt[2], popt[3], popt[3]

        return bins, counts, value, m_error, p_error

    else:

        raise PyLCProcessError('For the uncertainties in your distribution you should use either a confidence '
                               'interval or a gaussian fit, not both.')


def two_d_distribution(datax, datay, step=5.0, abs_step_x=None, abs_step_y=None, min_value_x=None, min_value_y=None,
                       max_value_x=None, max_value_y=None, mad_filter=None, samples=None):

    datax = datax.flatten()
    datay = datay.flatten()
    int_test = (np.sum((datax-np.int_(datax))**2) == 0) + (np.sum((datax-np.int_(datax))**2) == 0)

    if samples:
        datax = datax[::int(len(datax)/samples + 1)]
        datay = datay[::int(len(datay)/samples + 1)]

    if min_value_x is None:
        min_value_x = np.min(datax)
    else:
        values_to_keep = np.where(datax >= min_value_x)
        datax = datax[values_to_keep]
        datay = datay[values_to_keep]

    if max_value_x is None:
        max_value_x = np.max(datax)
    else:
        values_to_keep = np.where(datax <= max_value_x)
        datax = datax[values_to_keep]
        datay = datay[values_to_keep]

    if min_value_y is None:
        min_value_y = np.min(datay)
    else:
        values_to_keep = np.where(datay >= min_value_y)
        datax = datax[values_to_keep]
        datay = datay[values_to_keep]

    if max_value_y is None:
        max_value_y = np.max(datay)
    else:
        values_to_keep = np.where(datay <= max_value_y)
        datax = datax[values_to_keep]
        datay = datay[values_to_keep]

    if mad_filter:
        data_median = np.median(datax)
        data_mad = np.sqrt(np.median((datax - data_median) ** 2))
        min_value_x = data_median - mad_filter * data_mad
        max_value_x = data_median + mad_filter * data_mad

        data_median = np.median(datay)
        data_mad = np.sqrt(np.median((datay - data_median) ** 2))
        min_value_y = data_median - mad_filter * data_mad
        max_value_y = data_median + mad_filter * data_mad

        values_to_keep = np.where((datax <= max_value_x) * (datax >= min_value_x) *
                                  (datay <= max_value_y) * (datay >= min_value_y))

        datax = datax[values_to_keep]
        datay = datay[values_to_keep]

    if abs_step_x is None:
        new_step_x = np.sqrt(np.median((datax - np.median(datax)) ** 2)) / step
        if new_step_x == 0:
            new_step_x = np.sqrt(np.mean((datax - np.mean(datax)) ** 2)) / step
        step_x = new_step_x
    else:
        step_x = abs_step_x

    if int_test:
        step_x = max(1, round(step_x, 0))

    bins_number_x = int((max_value_x - min_value_x) / step_x) + 1
    bins_x = np.array(min_value_x + step_x / 2. + np.arange(bins_number_x) * step_x)

    if abs_step_y is None:
        new_step_y = np.sqrt(np.median((datay - np.median(datay)) ** 2)) / step
        if new_step_y == 0:
            new_step_y = np.sqrt(np.mean((datay - np.mean(datay)) ** 2)) / step
        step_y = new_step_y
    else:
        step_y = abs_step_y

    if int_test:
        step_y = max(1, round(step_y, 0))

    bins_number_y = int((max_value_y - min_value_y) / step_y) + 1
    bins_y = np.array(min_value_y + step_y / 2. + np.arange(bins_number_y) * step_y)

    bins_number_xy = bins_number_x * bins_number_y
    bins_xy = np.meshgrid(bins_x, bins_y)
    counts_xy = np.bincount((np.int_((datay - min_value_y) / step_y)) * bins_number_x +
                            (np.int_((datax - min_value_x) / step_x)))
    counts_xy = np.insert(counts_xy, len(counts_xy), np.zeros(bins_number_xy - len(counts_xy)))
    counts_xy = np.reshape(counts_xy, (bins_number_y, bins_number_x))

    return bins_xy[0], bins_xy[1], counts_xy
