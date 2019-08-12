from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from ._1databases import *


def correlation(x, y):
    n = len(x)
    mx = np.mean(x)
    sx = np.std(x)
    my = np.mean(y)
    sy = np.std(y)
    return np.round(np.sum((x - mx) * (y - my)) / ((n - 1) * sx * sy), 2)


def waverage(x, x_err):
    x = np.array(x)
    x_err = np.array(x_err)

    wx = 1 / (x_err ** 2)

    wav = np.sum(x * wx, 0) / np.sum(wx, 0)
    wav_err = 1 / np.sqrt((np.sum(wx, 0)))
    return wav, wav_err


def fit_line(datax, datay):

    mx = np.mean(datax)
    my = np.mean(datay)

    ssxx = np.sum((datax - mx) ** 2)
    ssyy = np.sum((datay - my) ** 2)
    ssxy = np.sum((datax - mx) * (datay - my))

    bb = ssxy / ssxx
    aa = my - bb * mx

    n = len(datax)
    sss = np.sqrt((ssyy - bb * ssxy) / (n - 2))

    aerr = sss * np.sqrt(1.0 / n + (mx ** 2) / ssxx)
    berr = sss / np.sqrt(ssxx)

    return aa, bb, aerr, berr


def curve_fit_line(datax, datay, datay_err):

    def fline(x, a, b):
        return a * x+ b

    p0 = np.polyfit(datax, datay, 1)

    popt, pcov = curve_fit(fline, datax, datay, sigma=datay_err, p0=p0)

    return popt[0], popt[1], np.sqrt(pcov[0][0]), np.sqrt(pcov[1][1])


def box(x_arr, x0, aa, bb, cc):
    return aa * np.exp(-np.power(np.power(x0 - x_arr, 2) / (2 * (bb ** 2)), cc))


def fit_box(datax, datay):

    minlim = datax[np.argmax(datay[5:] - datay[:-5])]
    maxlim = datax[np.argmin(datay[5:] - datay[:-5])]

    for expand in range(10, 30):
        datax = np.append(datax, max(datax) + expand)
        datay = np.append(datay, 0)
        datax = np.append(datax, min(datax) - expand)
        datay = np.append(datay, 0)

    popt, pcov = curve_fit(box, datax, datay,
                           p0=[0.5 * (maxlim + minlim), max(datay), 0.5 * (maxlim - minlim), 1.0])

    center = popt[0]

    center_err = np.sqrt(pcov[0][0])

    fwhm = popt[2] * 2.0 * np.sqrt(2.0 * (np.log(2) ** (1.0 / popt[3])))

    s, c = popt[2], popt[3]
    ss, cs = np.sqrt(pcov[2][2]), np.sqrt(pcov[3][3])
    fwhm_err = np.sqrt(8.0 * (ss ** 2) * (np.log(2.0) ** (1.0 / c))
                       + (2.0 * (cs ** 2) * (s ** 2) * (np.log(2.0) ** (1.0 / c))
                          * (np.log(np.log(2.0)) ** 2)) / (c ** 4))

    return center, center_err, fwhm, fwhm_err, popt


def values_to_print(value, error_minus, error_plus):

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

    print_value = '{0:.{width}f}'.format(round(value, width), width=width)
    print_m_error = '{0:.{width}f}'.format(round(error_minus, width), width=width)
    print_p_error = '{0:.{width}f}'.format(round(error_plus, width), width=width)

    return print_value, print_m_error, print_p_error


def find_demicals(number):

    xx = 25

    while round(number, xx) == round(number, xx - 1):
        xx -= 1

    return xx
