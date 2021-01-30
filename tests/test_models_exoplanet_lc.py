
import numpy as np
import pytest
import pylightcurve as plc
from pylightcurve.models.exoplanet_lc import _get_filter
from pylightcurve.__databases__ import plc_data


def test_exoplanet_lc():

    stellar_logg = 4.36
    stellar_temperature = 6065.0
    stellar_metallicity = 0.0
    rp_over_rs = 0.12086
    period = 3.5247486
    sma_over_rs = 8.76
    eccentricity = 0.0
    inclination = 86.71
    periastron = 0.0
    mid_time = 2452826.62928
    mid_time_format = 'BJD_TDB'
    time_array = np.arange(mid_time - 0.1, mid_time + 0.1, 0.001)
    eclipse_mid_time = plc.eclipse_mid_time(period, sma_over_rs, eccentricity, inclination, periastron, mid_time)
    eclipse_time_array = np.arange(eclipse_mid_time - 0.1, eclipse_mid_time + 0.1, 0.001)

    planet_orbit = plc.planet_orbit(period, sma_over_rs, eccentricity, inclination, periastron, mid_time, time_array)

    assert len(planet_orbit) == 3
    assert len(planet_orbit[0]) == len(time_array)
    assert len(planet_orbit[1]) == len(time_array)
    assert len(planet_orbit[2]) == len(time_array)
    assert min(planet_orbit[0]) > 8.6
    assert min(planet_orbit[1]) > -1.6
    assert min(planet_orbit[2]) > -0.6
    assert max(planet_orbit[0]) < 8.8
    assert max(planet_orbit[1]) < 1.6
    assert max(planet_orbit[2]) < -0.4

    planet_star_projected_distance = plc.planet_star_projected_distance(
        period, sma_over_rs, eccentricity, inclination, periastron, mid_time, time_array)

    assert min(planet_star_projected_distance) > 0.5
    assert max(planet_star_projected_distance) < 1.7

    planet_phase = plc.planet_phase(period, mid_time, time_array)

    assert min(planet_phase) > -0.1
    assert max(planet_phase) < 0.1

    for photometric_filter in plc_data.all_filters():

        for method in ['claret', 'sqrt', 'quad', 'linear']:

            for stellar_model in ['phoenix', 'atlas']:

                limb_darkening_coefficients = plc.exotethys(stellar_logg, stellar_temperature,
                                                            stellar_metallicity, photometric_filter,
                                                            method=method, stellar_model=stellar_model)

                assert len(limb_darkening_coefficients) == 4

                transit = plc.transit(limb_darkening_coefficients, rp_over_rs, period, sma_over_rs, eccentricity, inclination,
                                      periastron, mid_time, time_array, method=method)

                assert min(transit) > 0.98
                assert max(transit) == 1.0

                transit_integrated = plc.transit_integrated(limb_darkening_coefficients, rp_over_rs, period, sma_over_rs,
                                                            eccentricity, inclination, periastron, mid_time, time_array,
                                                            30, max_sub_exp_time=10, method=method, precision=3)

                assert min(transit_integrated) > 0.98
                assert max(transit_integrated) == 1.0

                assert round(plc.transit_duration(rp_over_rs, period, sma_over_rs, eccentricity, inclination,
                                                  periastron) * 24, 2) == 3.09

                assert round(plc.transit_depth(limb_darkening_coefficients, rp_over_rs, period, sma_over_rs, eccentricity,
                                               inclination, periastron), 2) <= 0.02

        fp_over_fs = plc.fp_over_fs(rp_over_rs, sma_over_rs, 0.2, 1.0, stellar_temperature, photometric_filter)

        eclipse = plc.eclipse(fp_over_fs, rp_over_rs, period, sma_over_rs, eccentricity, inclination, periastron,
                              eclipse_mid_time, eclipse_time_array)

        assert min(eclipse) > 0.99
        assert max(eclipse) == 1.0

        eclipse_integrated = plc.eclipse_integrated(fp_over_fs, rp_over_rs, period, sma_over_rs, eccentricity,
                                                    inclination, periastron, eclipse_mid_time, eclipse_time_array,
                                                    30, max_sub_exp_time=10)

        assert min(eclipse_integrated) > 0.998
        assert max(eclipse_integrated) == 1.0

        assert round(plc.eclipse_duration(rp_over_rs, period, sma_over_rs, eccentricity, inclination,
                                          periastron) * 24, 2) == 3.09

        assert round(plc.eclipse_depth(fp_over_fs, rp_over_rs, period, sma_over_rs, eccentricity,
                                       inclination, periastron), 5) <= 0.002

        planet_orbit = plc.planet_orbit(period, sma_over_rs, 0.1, inclination, 80.0, mid_time,
                                        time_array)

        assert len(planet_orbit) == 3
        assert len(planet_orbit[0]) == len(time_array)
        assert len(planet_orbit[1]) == len(time_array)
        assert len(planet_orbit[2]) == len(time_array)
        assert min(planet_orbit[0]) > 7.6
        assert min(planet_orbit[1]) > -1.8
        assert min(planet_orbit[2]) > -0.6
        assert max(planet_orbit[0]) < 8.8
        assert max(planet_orbit[1]) < 1.8
        assert max(planet_orbit[2]) < -0.4

        planet_orbit = plc.planet_orbit(period, sma_over_rs, 0.1, inclination, 130.0, mid_time,
                                        time_array)

        assert len(planet_orbit) == 3
        assert len(planet_orbit[0]) == len(time_array)
        assert len(planet_orbit[1]) == len(time_array)
        assert len(planet_orbit[2]) == len(time_array)
        assert min(planet_orbit[0]) > 7.6
        assert min(planet_orbit[1]) > -1.8
        assert min(planet_orbit[2]) > -0.6
        assert max(planet_orbit[0]) < 8.8
        assert max(planet_orbit[1]) < 1.8
        assert max(planet_orbit[2]) < -0.4

    with pytest.raises(plc.PyLCError):
        _get_filter('aa')

    with pytest.raises(plc.PyLCError):
        plc.exotethys(stellar_logg, stellar_temperature, stellar_metallicity, photometric_filter,
                      method='bla', stellar_model='phoenix')

    with pytest.raises(plc.PyLCError):
        plc.exotethys(stellar_logg, stellar_temperature, stellar_metallicity, photometric_filter,
                      method='claret', stellar_model='bla')

    with pytest.raises(plc.PyLCProcessError):
        plc.exotethys(1.0, 100000, 30, 'JOHNSON_V', method='claret', stellar_model='phoenix')

    assert len(plc.exotethys(stellar_logg, stellar_temperature, 0.1, 'JOHNSON_V', method='claret',
                             stellar_model='phoenix')) == 4
