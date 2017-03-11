__all__ = ['transit', 'transit_integrated']


import numpy as np

from transit_flux_drop import *
from exoplanet_orbit import *


def transit(method, limb_darkening_coefficients, rp_over_rs,
            period, sma_over_rs, eccentricity, inclination, periastron, mid_time, time_array, precision=3):

    position_vector = exoplanet_orbit(period, sma_over_rs, eccentricity, inclination, periastron, mid_time, time_array)

    projected_distance = np.where(
        position_vector[0] < 0, 1.0 + 5.0 * rp_over_rs,
        np.sqrt(position_vector[1] * position_vector[1] + position_vector[2] * position_vector[2]))

    return transit_flux_drop(method, limb_darkening_coefficients, rp_over_rs, projected_distance, precision=precision)


def transit_integrated(method, limb_darkening_coefficients, rp_over_rs,
                       period, sma_over_rs, eccentricity, inclination, periastron, mid_time, time_array,
                       exp_time, time_factor, precision=3):

    exp_time /= (60.0 * 60.0 * 24.0)

    time_factor = int(time_factor)

    time_array_hr = (time_array[:, None] + np.arange(-exp_time / 2 + exp_time / time_factor / 2, exp_time / 2,
                                                     exp_time / time_factor)).flatten()

    position_vector = exoplanet_orbit(period, sma_over_rs, eccentricity,
                                      inclination, periastron, mid_time, time_array_hr)

    projected_distance = np.where(
        position_vector[0] < 0, 1.0 + 5.0 * rp_over_rs,
        np.sqrt(position_vector[1] * position_vector[1] + position_vector[2] * position_vector[2]))

    return np.mean(np.reshape(transit_flux_drop(method, limb_darkening_coefficients,
                                                rp_over_rs, projected_distance, precision=precision),
                              (len(time_array), time_factor)), 1)

