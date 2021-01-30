
import os
import numpy as np
import pytest
import shutil
import pylightcurve as plc
from pylightcurve.__databases__ import plc_data


def test_exoplanet():

    ra = 330.795
    dec = 18.884
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

    planet = plc.Planet('HD209458b', plc.Degrees(ra), plc.Degrees(dec),
                        stellar_logg, stellar_temperature, stellar_metallicity,
                        rp_over_rs, period, sma_over_rs, eccentricity, plc.Degrees(inclination),
                        plc.Degrees(periastron), mid_time, mid_time_format)

    planet = plc.Planet('HD209458b', ra, dec,
                        stellar_logg, stellar_temperature, stellar_metallicity,
                        rp_over_rs, period, sma_over_rs, eccentricity, inclination,
                        periastron, mid_time, mid_time_format)

    eclipse_mid_time = planet.eclipse_mid_time
    eclipse_time_array = np.arange(eclipse_mid_time - 0.1, eclipse_mid_time + 0.1, 0.001)

    planet_orbit = planet.planet_orbit(time_array, 'BJD_TDB')

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

    planet_star_projected_distance = planet.planet_star_projected_distance(time_array, 'BJD_TDB')

    assert min(planet_star_projected_distance) > 0.5
    assert max(planet_star_projected_distance) < 1.7

    planet_phase = planet.planet_phase(time_array, 'BJD_TDB')

    assert min(planet_phase) > -0.1
    assert max(planet_phase) < 0.1

    for photometric_filter in plc_data.all_filters():

        print(photometric_filter)

        transit = planet.transit(time_array, 'BJD_TDB', photometric_filter)

        assert min(transit) > 0.98
        assert max(transit) == 1.0

        transit_integrated = planet.transit_integrated(time_array, 'BJD_TDB', 30, 'mid', photometric_filter)
        transit_integrated = planet.transit_integrated(time_array, 'BJD_TDB', 30, 'start', photometric_filter)
        transit_integrated = planet.transit_integrated(time_array, 'BJD_TDB', 30, 'end', photometric_filter)

        with pytest.raises(plc.PyLCInputError):
            transit_integrated = planet.transit_integrated(time_array, 'BJD_TDB', 30, 'bla', photometric_filter)

        assert min(transit_integrated) > 0.98
        assert max(transit_integrated) == 1.0

        assert round(planet.transit_duration(photometric_filter) * 24, 2) == 3.09

        assert round(planet.transit_depth(photometric_filter), 2) <= 0.02

        eclipse = planet.eclipse(eclipse_time_array, 'BJD_TDB', photometric_filter)

        assert min(eclipse) > 0.99
        assert max(eclipse) == 1.0

        eclipse_integrated = planet.eclipse_integrated(eclipse_time_array, 'BJD_TDB', 30, 'mid', photometric_filter)
        eclipse_integrated = planet.eclipse_integrated(eclipse_time_array, 'BJD_TDB', 30, 'start', photometric_filter)
        eclipse_integrated = planet.eclipse_integrated(eclipse_time_array, 'BJD_TDB', 30, 'end', photometric_filter)

        with pytest.raises(plc.PyLCInputError):
            eclipse_integrated = planet.eclipse_integrated(eclipse_time_array, 'BJD_TDB', 30, 'bla', photometric_filter)

        assert min(eclipse_integrated) > 0.998
        assert max(eclipse_integrated) == 1.0

        assert round(planet.eclipse_duration(photometric_filter) * 24, 2) == 3.09

        assert round(planet.eclipse_depth(photometric_filter), 5) <= 0.002

    with pytest.raises(plc.PyLCInputError):
        planet.filter('test')

    planet.add_filter('test', 0.1, 1, 1, 1, 1, 0.001)

    __location__ = os.path.abspath(os.path.dirname(__file__))

    planet.observations = {}

    data_time = time_array
    data_flux = planet.transit_integrated(time_array, 'BJD_TDB', 60, 'mid', 'COUSINS_R')

    planet.add_observation(data_time, 'BJD_TDB', 60, 'start',
                           data_flux + np.random.normal(0, 0.001, len(data_time)),
                           np.ones_like(data_flux) * 0.001, 'flux', 'COUSINS_R')

    planet.add_observation(data_time, 'BJD_TDB', 60, 'end',
                           data_flux + np.random.normal(0, 0.001, len(data_time)),
                           np.ones_like(data_flux) * 0.001, 'flux', 'COUSINS_R')

    planet.add_observation(data_time, 'BJD_TDB', 60, 'end',
                           data_flux + np.random.normal(0, 0.001, len(data_time)),
                           np.ones_like(data_flux) * 0.001, 'mag', 'COUSINS_R')

    with pytest.raises(plc.PyLCInputError):
        planet.add_observation(data_time, 'BJD_TDB', 60, 'bla',
                               data_flux + np.random.normal(0, 0.001, len(data_time)),
                               np.ones_like(data_flux) * 0.001, 'flux', 'COUSINS_R')

    with pytest.raises(plc.PyLCInputError):
        planet.add_observation(data_time, 'BJD_TDB', 60, 'mid',
                               data_flux + np.random.normal(0, 0.001, len(data_time)),
                               np.ones_like(data_flux) * 0.001, 'bla', 'COUSINS_R')

    planet.observations = {}

    data_time = time_array
    data_flux = planet.transit_integrated(time_array, 'BJD_TDB', 60, 'mid', 'COUSINS_R')

    planet.add_observation(data_time, 'BJD_TDB', 60, 'mid',
        data_flux + np.random.normal(0, 0.001, len(data_time)),
        np.ones_like(data_flux) * 0.001, 'flux', 'COUSINS_R')

    try:
        shutil.rmtree(os.path.join(__location__, 'test_transit'))
    except:
        pass

    planet.transit_fitting(os.path.join(__location__, 'test_transit'),
                           iterations=1300, walkers=200, burn_in=300,
                           detrending_order=2,
                           fit_rp_over_rs=True, fit_individual_rp_over_rs=True,
                           fit_mid_time=True, fit_individual_times=True)

    planet.transit_fitting(os.path.join(__location__, 'test_transit'),
                           iterations=1300, walkers=200, burn_in=300,
                           detrending_order=2,
                           fit_rp_over_rs=False, fit_individual_rp_over_rs=True,
                           fit_mid_time=False, fit_individual_times=True)

    planet.transit_fitting(os.path.join(__location__, 'test_transit'),
                           iterations=1300, walkers=200, burn_in=300,
                           detrending_order=1,
                           fit_rp_over_rs=True, fit_individual_rp_over_rs=False,
                           fit_mid_time=True, fit_individual_times=False)

    planet.transit_fitting(os.path.join(__location__, 'test_transit'),
                           iterations=1300, walkers=200, burn_in=300,
                           detrending_order=1,
                           fit_rp_over_rs=False, fit_individual_rp_over_rs=False,
                           fit_mid_time=False, fit_individual_times=False)

    planet.transit_fitting(os.path.join(__location__, 'test_transit'),
                           iterations=1300, walkers=200, burn_in=300,
                           detrending_order=0,
                           fit_ldc1=True, fit_ldc2=True, fit_ldc3=True, fit_ldc4=True,
                           fit_sma_over_rs=True, fit_inclination=True,
                           fit_mid_time=True, fit_individual_times=True
                           )

    planet.transit_fitting(os.path.join(__location__, 'test_transit'),
                           iterations=1300, walkers=200, burn_in=300,
                           detrending_order=0,
                           fit_period=True, fit_mid_time=True, fit_individual_times=False
                           )

    with pytest.raises(plc.PyLCInputError):
        planet.transit_fitting(os.path.join(__location__, 'test_transit'),
                               iterations=1300, walkers=200, burn_in=300,
                               detrending_order=0,
                               fit_period=True, fit_mid_time=True, fit_individual_times=True
                               )

    planet.observations = {}

    data_time = eclipse_time_array
    data_flux = planet.eclipse_integrated(eclipse_time_array, 'BJD_TDB', 60, 'mid', 'COUSINS_R')

    planet.add_observation(data_time, 'BJD_TDB', 60, 'mid',
                           data_flux + np.random.normal(0, 0.00001, len(data_time)),
                           np.ones_like(data_flux) * 0.00001, 'flux', 'COUSINS_R')

    try:
        shutil.rmtree(os.path.join(__location__, 'test_eclipse'))
    except:
        pass

    planet.eclipse_fitting(os.path.join(__location__, 'test_eclipse'),
                           iterations=1300, walkers=200, burn_in=300,
                           detrending_order=2,
                           fit_rp_over_rs=True, fit_individual_rp_over_rs=True,
                           fit_fp_over_fs=True, fit_individual_fp_over_fs=True,
                           fit_mid_time=True, fit_individual_times=True)

    planet.eclipse_fitting(os.path.join(__location__, 'test_eclipse'),
                           iterations=1300, walkers=200, burn_in=300,
                           detrending_order=2,
                           fit_rp_over_rs=False, fit_individual_rp_over_rs=True,
                           fit_fp_over_fs=False, fit_individual_fp_over_fs=True,
                           fit_mid_time=False, fit_individual_times=True)

    planet.eclipse_fitting(os.path.join(__location__, 'test_eclipse'),
                           iterations=1300, walkers=200, burn_in=300,
                           detrending_order=1,
                           fit_rp_over_rs=True, fit_individual_rp_over_rs=False,
                           fit_fp_over_fs=True, fit_individual_fp_over_fs=False,
                           fit_mid_time=True, fit_individual_times=False)

    planet.eclipse_fitting(os.path.join(__location__, 'test_eclipse'),
                           iterations=1300, walkers=200, burn_in=300,
                           detrending_order=1,
                           fit_rp_over_rs=False, fit_individual_rp_over_rs=False,
                           fit_fp_over_fs=False, fit_individual_fp_over_fs=False,
                           fit_mid_time=False, fit_individual_times=False)

    planet.eclipse_fitting(os.path.join(__location__, 'test_eclipse'),
                           iterations=1300, walkers=200, burn_in=300,
                           detrending_order=0,
                           fit_sma_over_rs=True, fit_inclination=True,
                           fit_mid_time=True, fit_individual_times=True
                           )

    planet.eclipse_fitting(os.path.join(__location__, 'test_eclipse'),
                           iterations=1300, walkers=200, burn_in=300,
                           detrending_order=0,
                           fit_period=True, fit_mid_time=True, fit_individual_times=False
                           )

    with pytest.raises(plc.PyLCInputError):
        planet.eclipse_fitting(os.path.join(__location__, 'test_eclipse'),
                               iterations=1300, walkers=200, burn_in=300,
                               detrending_order=0,
                               fit_period=True, fit_mid_time=True, fit_individual_times=True
                               )
