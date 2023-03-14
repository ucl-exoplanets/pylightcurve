
import os
import numpy as np
import pytest
import pylightcurve as plc


def test_exoplanet_lc():

    ra = plc.Hours(3, 15).deg()
    dec = plc.Degrees(45.7).deg_coord()
    stellar_logg = 4.36
    stellar_temperature = 6065.0
    stellar_metallicity = 0.1
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

    assert round(plc.transit_duration(rp_over_rs, period, sma_over_rs, eccentricity, inclination,
                                      periastron) * 24, 2) == 3.09

    assert round(plc.transit_t12(rp_over_rs, period, sma_over_rs, eccentricity, inclination,
                                 periastron) * 24, 2) == 0.43

    assert round(plc.eclipse_duration(rp_over_rs, period, sma_over_rs, eccentricity, inclination,
                                      periastron) * 24, 2) == 3.09

    planet_orbit = plc.planet_orbit(period, sma_over_rs, 0.1, inclination, 80.0, mid_time,
                                    time_array, ww=90)
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

    __location__ = os.path.abspath(os.path.dirname(__file__))

    with pytest.raises(plc.PyLCError):
        plc.plc_data.add_filter('custom_filter', os.path.join(__location__, 'custom_f.txt'))

    plc.plc_data.add_filter('custom_filter', os.path.join(__location__, 'custom_filter.txt'))

    photometric_filter = 'COUSINS_R'

    for method in ['claret', 'quad', 'linear']:

        for stellar_model in ['Phoenix_2018', 'Atlas_2000']:

            limb_darkening_coefficients = plc.exotethys(stellar_logg, stellar_temperature,
                                                        stellar_metallicity, photometric_filter,
                                                        ldc_method=method, stellar_model=stellar_model)

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

    assert round(plc.eclipse_depth(fp_over_fs, rp_over_rs, period, sma_over_rs, eccentricity,
                                   inclination, periastron), 5) <= 0.002

    with pytest.raises(plc.PyLCError):
        plc.exotethys(stellar_logg, stellar_temperature, stellar_metallicity, 'bla',
                      ldc_method='claret', stellar_model=stellar_model)

    with pytest.raises(plc.PyLCError):
        plc.exotethys(stellar_logg, stellar_temperature, stellar_metallicity, photometric_filter,
                      ldc_method='bla', stellar_model=stellar_model)

    with pytest.raises(plc.PyLCError):
        plc.exotethys(stellar_logg, stellar_temperature, stellar_metallicity, photometric_filter,
                      ldc_method='claret', stellar_model='bla')

    with pytest.raises(plc.PyLCError):
        plc.exotethys(stellar_logg, stellar_temperature, stellar_metallicity, photometric_filter,
                      ldc_method='claret', stellar_model=stellar_model, wlrange=[5, 10])

    with pytest.raises(plc.PyLCError):
        plc.fp_over_fs(rp_over_rs, sma_over_rs, 0.2, 1.0, stellar_temperature, 'bla')

    with pytest.raises(plc.PyLCError):
        plc.fp_over_fs(rp_over_rs, sma_over_rs, 0.2, 1.0, stellar_temperature, photometric_filter, wlrange=[5, 10])

    assert len(plc.exotethys(stellar_logg, stellar_temperature, 0.1, 'JOHNSON_V', ldc_method='claret',
                             stellar_model=stellar_model)) == 4

    assert round(plc.convert_to_bjd_tdb(ra, dec, 2458485.0, 'JD_UTC'), 7) == 2458485.0046340
    assert round(plc.convert_to_bjd_tdb(ra, dec, 58484.5, 'MJD_UTC'), 7) == 2458485.0046340
    assert round(plc.convert_to_bjd_tdb(ra, dec, 2458485.0046340, 'BJD_TDB'), 7) == 2458485.0046340
    assert round(plc.convert_to_bjd_tdb(ra, dec, 2458485.0038333, 'BJD_UTC'), 7) == 2458485.0046340
    assert round(plc.convert_to_bjd_tdb(ra, dec, 2458485.0046033, 'HJD_TDB'), 7) == 2458485.0046340
    assert round(plc.convert_to_bjd_tdb(ra, dec, 2458485.00380255, 'HJD_UTC'), 7) == 2458485.0046340
    assert round(plc.convert_to_bjd_tdb(ra, dec, np.array([2458485.00000000]), 'JD_UTC')[0], 7) == 2458485.0046340
    assert round(plc.convert_to_bjd_tdb(ra, dec, np.array([0058484.50000000]), 'MJD_UTC')[0], 7) == 2458485.0046340
    assert round(plc.convert_to_bjd_tdb(ra, dec, np.array([2458485.00463400]), 'BJD_TDB')[0], 7) == 2458485.0046340
    assert round(plc.convert_to_bjd_tdb(ra, dec, np.array([2458485.00383330]), 'BJD_UTC')[0], 7) == 2458485.0046340
    assert round(plc.convert_to_bjd_tdb(ra, dec, np.array([2458485.00460330]), 'HJD_TDB')[0], 7) == 2458485.0046340
    assert round(plc.convert_to_bjd_tdb(ra, dec, np.array([2458485.00380255]), 'HJD_UTC')[0], 7) == 2458485.0046340
    assert round(plc.convert_to_jd_utc(ra, dec, 2458485.0, 'JD_UTC'), 7) == 2458485.0
    assert round(plc.convert_to_jd_utc(ra, dec, 58484.5, 'MJD_UTC'), 7) == 2458485.0
    assert round(plc.convert_to_jd_utc(ra, dec, 2458485.0046340, 'BJD_TDB'), 7) == 2458485.0
    assert round(plc.convert_to_jd_utc(ra, dec, 2458485.0038333, 'BJD_UTC'), 7) == 2458485.0
    assert round(plc.convert_to_jd_utc(ra, dec, 2458485.0046033, 'HJD_TDB'), 7) == 2458485.0
    assert round(plc.convert_to_jd_utc(ra, dec, 2458485.00380255, 'HJD_UTC'), 7) == 2458485.0
    assert round(plc.convert_to_jd_utc(ra, dec, np.array([2458485.00000000]), 'JD_UTC')[0], 7) == 2458485.0
    assert round(plc.convert_to_jd_utc(ra, dec, np.array([0058484.50000000]), 'MJD_UTC')[0], 7) == 2458485.0
    assert round(plc.convert_to_jd_utc(ra, dec, np.array([2458485.00463400]), 'BJD_TDB')[0], 7) == 2458485.0
    assert round(plc.convert_to_jd_utc(ra, dec, np.array([2458485.00383330]), 'BJD_UTC')[0], 7) == 2458485.0
    assert round(plc.convert_to_jd_utc(ra, dec, np.array([2458485.00460330]), 'HJD_TDB')[0], 7) == 2458485.0
    assert round(plc.convert_to_jd_utc(ra, dec, np.array([2458485.00380255]), 'HJD_UTC')[0], 7) == 2458485.0

    with pytest.raises(plc.PyLCInputError):
        plc.convert_to_bjd_tdb(ra, dec, 2458485.0, 'aaa')

    with pytest.raises(plc.PyLCInputError):
        plc.convert_to_bjd_tdb(ra, dec, [2458485.0], 'aaa')

    with pytest.raises(plc.PyLCInputError):
        plc.convert_to_bjd_tdb(ra, dec, 'a', 'JD_UTC')

    with pytest.raises(plc.PyLCInputError):
        plc.convert_to_bjd_tdb(ra, dec, np.array(['a']), 'JD_UTC')

    with pytest.raises(plc.PyLCInputError):
        plc.convert_to_jd_utc(ra, dec, 2458485.0, 'aaa')

    with pytest.raises(plc.PyLCInputError):
        plc.convert_to_jd_utc(ra, dec, [2458485.0], 'aaa')

    with pytest.raises(plc.PyLCInputError):
        plc.convert_to_jd_utc(ra, dec, 'a', 'JD_UTC')

    with pytest.raises(plc.PyLCInputError):
        plc.convert_to_jd_utc(ra, dec, np.array(['a']), 'JD_UTC')


