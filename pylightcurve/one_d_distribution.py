__all__ = ['one_d_distribution']


import numpy as np


def one_d_distribution(data_array, step=None, min_value=None, max_value=None, confidence_interval=None):

    try:
        data_array = data_array.flatten()
    except AttributeError:
        raise RuntimeError('Invalid input data')

    if step is None:
        step = np.sqrt(np.median((data_array - np.median(data_array)) ** 2)) / 5.0
    else:
        try:
            step = float(step)
        except ValueError:
            raise RuntimeError('Invalid step value')

    if min_value is None:
        min_value = min(data_array)
    else:
        try:
            min_value = float(min_value)
        except ValueError:
            raise RuntimeError('Invalid minimum value')

    if max_value is None:
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
                print "ERROR : confidence level can not be reached from left"
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
                print "ERROR : confidence level can not be reached from right"
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
            done1 = 1 / int(error_minus * (10 ** digit1))
        except ZeroDivisionError:
            done1 = 0
        try:
            done2 = 1 / int(error_plus * (10 ** digit2))
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
        print_p_error = '{0:.{width}f}'.format(round(error_minus, width), width=width)

        return bins, counts, value, m_error, p_error, print_value, print_m_error, print_p_error
