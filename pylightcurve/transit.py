from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from .transit_flux_drop import *
from .exoplanet_orbit import *


def transit(method, limb_darkening_coefficients, rp_over_rs, period, sma_over_rs, eccentricity, inclination, periastron,
            mid_time, time_array, precision=3):

    position_vector = exoplanet_orbit(period, sma_over_rs, eccentricity, inclination, periastron, mid_time, time_array)

    projected_distance = np.where(
        position_vector[0] < 0, 1.0 + 5.0 * rp_over_rs,
        np.sqrt(position_vector[1] * position_vector[1] + position_vector[2] * position_vector[2]))

    return transit_flux_drop(method, limb_darkening_coefficients, rp_over_rs, projected_distance, precision=precision)


def eclipse(fp_over_fs, rp_over_rs, period, sma_over_rs, eccentricity, inclination, periastron, mid_time, time_array,
            precision=3):

    position_vector = exoplanet_orbit(period, -sma_over_rs / rp_over_rs, eccentricity, inclination, periastron,
                                      mid_time, time_array)

    projected_distance = np.where(
        position_vector[0] < 0, 1.0 + 5.0 / rp_over_rs,
        np.sqrt(position_vector[1] * position_vector[1] + position_vector[2] * position_vector[2]))

    return (1.0 + fp_over_fs * transit_flux_drop('claret', [0, 0, 0, 0], 1 / rp_over_rs, projected_distance,
                                                 precision=precision)) / (1.0 + fp_over_fs)


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


def transit_and_hst():
    import numpy as np
    from astropy import units as u

    def VisitPlanner(exptime, num_orbits,
                     time_per_orbit=54 * u.min, hst_period=95 * u.min,
                     exp_overhead=1 * u.min):

        exp_per_dump = 10000000
        time_buffer_dump = 5.8 * u.min  # IH 23 pg 209

        exp_times = []
        orbit_start_index = []
        buffer_dump_index = []
        for orbit_n in xrange(num_orbits):
            if orbit_n == 0:
                guide_star_aq = 6 * u.min
            else:
                guide_star_aq = 5 * u.min

            # record the expnum of each orbit start
            orbit_start_index.append(len(exp_times))
            start_time = hst_period * orbit_n

            visit_time = start_time + guide_star_aq
            visit_end_time = start_time + time_per_orbit

            exp_n = 0  # For buffer dumps - with mixed exp types we should really
            #  track headers and size
            while visit_time < visit_end_time:

                # you cant convert a list of quantities to an array so we have to
                #  either know the length to preset one or
                # use floats in a list and convert after.
                exp_times.append(visit_time.to(u.min).value)  # start of exposure
                visit_time += (exptime + exp_overhead)

                exp_n += 1
                if exp_n > exp_per_dump:
                    visit_time += time_buffer_dump
                    exp_n = 0

                    buffer_dump_index.append(len(exp_times))

        returnDict = {
            'exp_times': np.array(exp_times) * u.min,  # start_times?
            'NSAMP': NSAMP,
            'SAMPSEQ': SAMPSEQ,
            'SUBARRAY': SUBARRAY,
            'num_exp': len(exp_times),
            'exptime': exptime,
            'num_orbits': num_orbits,
            'exp_overhead': exp_overhead,
            'time_per_orbit': time_per_orbit,
            'hst_period': hst_period,
            'buffer_dump_index': buffer_dump_index,
            'orbit_start_index': orbit_start_index,
        }

        return returnDict
