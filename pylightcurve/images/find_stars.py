
__all__ = ['find_single_star']

import numpy as np
import warnings

from pylightcurve.analysis.gaussian import fit_two_d_gaussian
from pylightcurve.analysis.distributions import one_d_distribution


def find_single_star(data_array, predicted_x, predicted_y, mean=None, std=None, burn_limit=65000, star_std=2,
                     std_limit=5.0):
    star = None

    if 0 < predicted_x < len(data_array[0]) and 0 < predicted_x < len(data_array):

        if mean is None or std is None:
            fit_mean, fit_std = one_d_distribution(data_array, gaussian_fit=True, mad_filter=5)[2:4]

            if not mean:
                mean = fit_mean

            if not std:
                std = fit_std

        centroids = find_centroids(data_array, predicted_x - 5 * star_std, predicted_x + 5 * star_std,
                                   predicted_y - 5 * star_std, predicted_y + 5 * star_std, mean, std, burn_limit, star_std,
                                   std_limit)

        centroids = sorted(centroids, key=lambda x: np.sqrt((x[0] - predicted_x) ** 2 + (x[1] - predicted_y) ** 2))

        for centroid in centroids:
            star = _star_from_centroid(data_array, centroid[0], centroid[1], mean, std, burn_limit, star_std, std_limit)
            if star:
                star = [star[0][2], star[0][3], star[0][0], star[0][1], star[0][4], star[0][5], centroid[0], centroid[1]]
                break

    return star


def _star_from_centroid(data_array, centroid_x, centroid_y, mean, std, burn_limit, star_std, std_limit):

    star = None
    try:
        search_window = int(round(10 * star_std))
        y_min = int(max(int(centroid_y) - search_window, 0))
        y_max = int(min(int(centroid_y) + search_window, len(data_array) - 1))
        x_min = int(max(int(centroid_x) - search_window, 0))
        x_max = int(min(int(centroid_x) + search_window, len(data_array[0]) - 1))

        datax, datay = np.meshgrid(np.arange(x_min, x_max + 1) + 0.5,
                                   np.arange(y_min, y_max + 1) + 0.5)

        dataz = data_array[y_min: y_max + 1, x_min: x_max + 1]
        popt, pcov = fit_two_d_gaussian(datax, datay, dataz, positive=True, point_xy=(centroid_x, centroid_y),
                                        sigma=star_std, maxfev=1000)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if popt[0] > std_limit * std and popt[0] + popt[1] < burn_limit:
                if np.sqrt(pcov[0][0]) != np.inf:
                    if popt[0] > std_limit * np.sqrt(pcov[0][0]):
                        star = (popt, pcov)
                else:
                    star = (popt, pcov)
    except:
        pass

    return star


def find_centroids(data_array, x_low, x_upper, y_low, y_upper, mean, std, burn_limit, star_std, std_limit):

    x_upper = int(min(x_upper, len(data_array[0])))
    y_upper = int(min(y_upper, len(data_array)))
    x_low = int(max(0, x_low))
    y_low = int(max(0, y_low))

    data_array = np.full_like(data_array[y_low:y_upper + 1, x_low:x_upper + 1],
                              data_array[y_low:y_upper + 1, x_low:x_upper + 1])

    test = []

    for i in range(-star_std, star_std + 1):
        for j in range(-star_std, star_std + 1):
            rolled = np.roll(np.roll(data_array, i, 0), j, 1)
            test.append(rolled)

    median_test = np.median(test, 0)
    max_test = np.max(test, 0)
    del test
    stars = np.where((data_array < burn_limit) & (data_array > mean + std_limit * std) & (max_test == data_array)
                     & (median_test > mean + 2 * std))
    del data_array

    stars = [stars[1] + x_low, stars[0] + y_low]
    stars = np.swapaxes(stars, 0, 1)

    return stars
