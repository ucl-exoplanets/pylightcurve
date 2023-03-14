
import os
import numpy as np
import pytest
import pylightcurve as plc


def test_exoplanet():

    __location__ = os.path.abspath(os.path.dirname(__file__))
    plc.plc_data.add_filter('custom_filter', os.path.join(__location__, 'custom_filter.txt'))

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
    transit_time_array = np.arange(mid_time - 0.1, mid_time + 0.1, 0.001)

    planet = plc.Planet('HD209458b', ra, dec,
                        stellar_logg, stellar_temperature, stellar_metallicity,
                        rp_over_rs, period, sma_over_rs, eccentricity, inclination,
                        periastron, mid_time)

    eclipse_mid_time = planet.eclipse_mid_time
    eclipse_time_array = np.arange(eclipse_mid_time - 0.1, eclipse_mid_time + 0.1, 0.001)

    planet_orbit = planet.planet_orbit(transit_time_array)

    assert len(planet_orbit) == 3
    assert len(planet_orbit[0]) == len(transit_time_array)
    assert len(planet_orbit[1]) == len(transit_time_array)
    assert len(planet_orbit[2]) == len(transit_time_array)
    assert min(planet_orbit[0]) > 8.6
    assert min(planet_orbit[1]) > -1.6
    assert min(planet_orbit[2]) > -0.6
    assert max(planet_orbit[0]) < 8.8
    assert max(planet_orbit[1]) < 1.6
    assert max(planet_orbit[2]) < -0.4

    planet_star_projected_distance = planet.planet_star_projected_distance(transit_time_array)

    assert min(planet_star_projected_distance) > 0.5
    assert max(planet_star_projected_distance) < 1.7

    planet_phase = planet.planet_phase(transit_time_array)

    assert min(planet_phase) > -0.1
    assert max(planet_phase) < 0.1

    assert round(planet.transit_duration() * 24, 2) == 3.09
    assert round(planet.transit_t12() * 24, 2) == 0.43
    assert round(planet.eclipse_duration() * 24, 2) == 3.09

    transit = planet.transit(transit_time_array, 'custom_filter', wlrange=[7000, 8000])

    assert min(transit) < 1.0 and min(transit) > 0.9
    assert max(transit) == 1.0

    eclipse = planet.eclipse(eclipse_time_array, 'custom_filter', wlrange=[7000, 8000])

    assert min(eclipse) < 1.0 and min(eclipse) > 0.9
    assert max(eclipse) == 1.0

    for photometric_filter in plc.plc_data.all_filters():

        if 'hst_wfc3_g141_' not in photometric_filter and 'hst_wfc3_g102_' not in photometric_filter:

            if 'irac' in photometric_filter:
                stellar_model = 'Atlas_2000'
            elif 'jwst_niriss' in photometric_filter:
                stellar_model = 'Stagger_2018'
            else:
                stellar_model = 'Phoenix_2018'

            transit = planet.transit(transit_time_array, photometric_filter, stellar_model=stellar_model)

            assert min(transit) > 0.98
            assert max(transit) == 1.0

            transit_integrated = planet.transit_integrated(transit_time_array, 30,
                                                           photometric_filter, stellar_model=stellar_model)

            assert min(transit_integrated) > 0.98
            assert max(transit_integrated) == 1.0

            assert round(planet.transit_depth(photometric_filter, stellar_model=stellar_model), 2) <= 0.02

            eclipse = planet.eclipse(eclipse_time_array, photometric_filter)

            assert min(eclipse) > 0.99
            assert max(eclipse) == 1.0

            eclipse_integrated = planet.eclipse_integrated(eclipse_time_array, 30, photometric_filter)

            assert min(eclipse_integrated) > 0.998
            assert max(eclipse_integrated) == 1.0

            assert round(planet.eclipse_depth(photometric_filter), 5) <= 0.002
            break

    planet.add_custom_limb_darkening_coefficients([0.1, 0.1, 0.1, 0.1], 'filter_name', None, 'custom')

    # test adding observations

    obs0_time = transit_time_array
    obs0_flux = 100000 * (planet.transit_integrated(obs0_time, 60, 'COUSINS_R') +
                          np.random.normal(0, 0.001, len(obs0_time)))
    obs0_flux_unc = 100000 * np.ones_like(obs0_flux) * 0.001

    planet.add_observation(obs0_time, 'BJD_TDB', 60, 'mid', obs0_flux, obs0_flux_unc, 'flux', 'COUSINS_R')
    planet.add_observation_from_dict({'time': obs0_time, 'time_format': 'BJD_TDB', 'exp_time': 60,
                                      'time_stamp': 'start', 'flux': obs0_flux, 'flux_unc': obs0_flux_unc,
                                      'flux_format': 'flux', 'filter_name': 'COUSINS_R'
                                      })
    planet.add_observation(obs0_time, 'BJD_TDB', 60, 'mid', obs0_flux, obs0_flux_unc, 'flux', 'COUSINS_R',
                           auxiliary_data={'test': obs0_time + 10},
                           observatory_latitude=40.0, observatory_longitude=23.0,
                           detrending_series=['time', 'test']
                           )

    for observation in planet.observations:
        observation.reset(planet)

    with pytest.raises(plc.PyLCInputError):
        planet.add_observation(obs0_time, 'BJD_TDB', 60, 'bla', obs0_flux, obs0_flux_unc, 'flux', 'COUSINS_R')

    with pytest.raises(plc.PyLCInputError):
        planet.add_observation(obs0_time, 'BJD_TDB', 60, 'mid', obs0_flux, obs0_flux_unc, 'flux', 'COUSINS_R',
                               detrending_series='a')

    with pytest.raises(plc.PyLCInputError):
        planet.add_observation(obs0_time, 'BJD_TDB', 60, 'mid', obs0_flux, obs0_flux_unc, 'flux', 'COUSINS_R',
                               detrending_series=['a', 'b'])

    with pytest.raises(plc.PyLCInputError):
        planet.add_observation(obs0_time, 'BJD_TDB', 60, 'mid', obs0_flux, obs0_flux_unc, 'flux', 'COUSINS_R',
                               detrending_series=None)

    with pytest.raises(plc.PyLCInputError):
        planet.add_observation(obs0_time, 'BJD_TDB', 60, 'mid', obs0_flux, obs0_flux_unc, 'flux', 'COUSINS_R',
                               detrending_series=5)

    with pytest.raises(plc.PyLCInputError):
        planet.add_observation(obs0_time, 'BJD_TDB', 60, 'mid', obs0_flux, obs0_flux_unc, 'flux', 'COUSINS_R',
                               detrending_series='time', detrending_order=-3)

    with pytest.raises(plc.PyLCInputError):
        planet.add_observation_from_dict(10)

    # test fitting transit observations

    obs0_time = transit_time_array + planet.period
    obs0_flux = 10000 * (planet.transit_integrated(obs0_time, 60, 'COUSINS_I', wlrange=[7500, 8500]) +
                         np.random.normal(0, 0.001, len(obs0_time)))
    obs0_flux_unc = 10000 * np.ones_like(obs0_flux) * 0.001

    obs1_time = transit_time_array
    obs1_flux = 100000 * (planet.transit_integrated(obs1_time, 60, 'COUSINS_R') +
                          np.random.normal(0, 0.001, len(obs1_time)))
    obs1_flux_unc = 100000 * np.ones_like(obs1_flux) * 0.001

    planet.clear_observations()
    planet.add_observation(obs0_time, 'BJD_TDB', 60, 'mid', obs0_flux, obs0_flux_unc, 'flux', 'COUSINS_I',
                           wlrange=[7500, 8500])
    planet.add_observation(obs1_time, 'BJD_TDB', 60, 'mid', obs1_flux, obs1_flux_unc, 'flux', 'COUSINS_R')

    planet.transit_fitting(
        os.path.join(__location__, 'test_transit'),
        iterations=200, burn_in=10,
        scale_uncertainties=True,
        filter_outliers=True,
        fit_individual_rp_over_rs=True,
    )
    planet.transit_fitting(
        iterations=200, burn_in=10,
        fit_individual_rp_over_rs=True,
        fit_individual_times=False,
        fit_mid_time=True, fit_period=True,
    )

    def trend_function(auxiliary_data, c1):
        return 1.0 + c1 * auxiliary_data['time']

    trend_parameters = [[0, -1, 1, 'ls', '$l_s$']]
    planet.clear_observations()
    planet.add_observation(obs1_time, 'BJD_TDB', 60, 'mid', obs1_flux, obs1_flux_unc, 'flux', 'COUSINS_R',
                           trend_function=trend_function, trend_parameters=trend_parameters)
    planet.transit_fitting(
        iterations=200, burn_in=10,
        fit_mid_time=False, fit_individual_times=False,
        fit_ldc1=True, fit_ldc2=True, fit_ldc3=True, fit_ldc4=True,
        fit_sma_over_rs=True, fit_inclination=True,
        return_traces=False, optimiser='curve_fit',
    )

    with pytest.raises(plc.PyLCInputError):
        planet.clear_observations()
        planet.add_observation(obs0_time, 'BJD_TDB', 60, 'mid', obs0_flux, obs0_flux_unc, 'flux', 'COUSINS_I',
                               wlrange=[7500, 8500])
        planet.add_observation(obs1_time, 'BJD_TDB', 60, 'mid', obs1_flux, obs1_flux_unc, 'flux', 'COUSINS_R')
        planet.transit_fitting(fit_period=True, fit_mid_time=True, fit_individual_times=True)

    with pytest.raises(plc.PyLCInputError):
        planet.clear_observations()
        planet.add_observation(obs1_time, 'BJD_TDB', 60, 'mid', obs1_flux, obs1_flux_unc, 'flux', 'COUSINS_R')
        planet.transit_fitting(fit_period=True)

    with pytest.raises(plc.PyLCInputError):
        planet.clear_observations()
        planet.add_observation(obs1_time, 'BJD_TDB', 60, 'mid', obs1_flux, obs1_flux_unc, 'flux', 'COUSINS_R')
        planet.eclipse_fitting()

    with pytest.raises(plc.PyLCInputError):
        planet.clear_observations()
        planet.transit_fitting()

    # test fitting eclipse observations

    obs0_time = eclipse_time_array + planet.period
    obs0_flux = 10000 * (planet.eclipse_integrated(obs0_time, 60, 'COUSINS_I', wlrange=[7500, 8500]) +
                         np.random.normal(0, 0.000001, len(obs0_time)))
    obs0_flux_unc = 10000 * np.ones_like(obs0_flux) * 0.000001

    obs1_time = eclipse_time_array
    obs1_flux = 100000 * (planet.eclipse_integrated(obs1_time, 60, 'COUSINS_R') +
                          np.random.normal(0, 0.000001, len(obs1_time)))
    obs1_flux_unc = 100000 * np.ones_like(obs1_flux) * 0.000001

    planet.clear_observations()
    planet.add_observation(obs0_time, 'BJD_TDB', 60, 'start', obs0_flux, obs0_flux_unc, 'flux', 'COUSINS_I',
                           wlrange=[7500, 8500])
    planet.add_observation(obs1_time, 'BJD_TDB', 60, 'mid', obs1_flux, obs1_flux_unc, 'flux', 'COUSINS_R')

    planet.eclipse_fitting(
        os.path.join(__location__, 'test_eclipse'),
        iterations=200, burn_in=10,
        fit_individual_fp_over_fs=True,
    )
    planet.eclipse_fitting(
        iterations=200, burn_in=10,
        fit_individual_fp_over_fs=True,
        fit_rp_over_rs=True,
        fit_individual_rp_over_rs=True,
        fit_individual_times=False,
        fit_mid_time=True, fit_period=True,
        optimise_initial_parameters_trials=10,
    )

    planet.clear_observations()
    planet.add_observation(obs1_time, 'BJD_TDB', 60, 'mid', obs1_flux, obs1_flux_unc, 'flux', 'COUSINS_R')
    planet.eclipse_fitting(
        iterations=200, burn_in=10,
        fit_individual_times=False, fit_mid_time=False,
        fit_rp_over_rs=True, fit_sma_over_rs=True, fit_inclination=True,
        return_traces=False, optimiser='curve_fit',
    )

    with pytest.raises(plc.PyLCInputError):
        planet.clear_observations()
        planet.add_observation(obs0_time, 'BJD_TDB', 60, 'start', obs0_flux, obs0_flux_unc, 'flux', 'COUSINS_I',
                               wlrange=[7500, 8500])
        planet.add_observation(obs1_time, 'BJD_TDB', 60, 'mid', obs1_flux, obs1_flux_unc, 'flux', 'COUSINS_R')
        planet.eclipse_fitting(fit_period=True, fit_mid_time=True, fit_individual_times=True)

    with pytest.raises(plc.PyLCInputError):
        planet.clear_observations()
        planet.add_observation(obs1_time, 'BJD_TDB', 60, 'mid', obs1_flux, obs1_flux_unc, 'flux', 'COUSINS_R')
        planet.eclipse_fitting(fit_period=True)

    with pytest.raises(plc.PyLCInputError):
        planet.clear_observations()
        planet.add_observation(obs1_time, 'BJD_TDB', 60, 'mid', obs1_flux, obs1_flux_unc, 'flux', 'COUSINS_R')
        planet.transit_fitting()

    with pytest.raises(plc.PyLCInputError):
        planet.clear_observations()
        planet.eclipse_fitting()

    planet.performance_report(100)
