from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import pylightcurve as plc
plt = plc.plt
np = plc.np

# The main functions of pyligthcurve are:
# plc.find_oec_parameters
# plc.clablimb
# plc.transit_projected_distance
# plc.transit_flux_drop
# plc.transit
# plc.mcmc_transit


# plc.find_oec_parameters(target)
#
#     returns the stellar and transit parametrs
#     (planet oec name, stellar logg, stellar effective temperature, stellar metallicity,
#      relative rafius, relative bolometric emission, period, relative semi-major axis,
#      eccentricity, inclination, periastron, transit mid-time)
#     given the name of the planet
#
#     for examples, for HD209458b we can have:

(planet, logg, effective_temperature, metallicity, rp_over_rs, fp_over_fs,
 period, sma_over_rs, eccentricity, inclination, periastron, mid_time) = plc.find_oec_parameters('hd209458b')

print(planet)
print(logg)
print(effective_temperature)
print(metallicity)
print(rp_over_rs)
print(fp_over_fs)
print(period)
print(sma_over_rs)
print(eccentricity)
print(inclination)
print(periastron)
print(mid_time)


# plc.clablimb(method, logg, effective_temperature, metallicity, photometric_filter, stellar_model='ATLAS')
# 	
#     returns the limb darkening coefficients 
#     (4 coefficients non-linear low from Claret 2000)
#     given the stellar parameters and the photometric filter used for the observation
#
#     method                 ['claret',
#                             'quad' (not available yet), 'sqrt' (not available yet), 'linear' (not available yet)]
#     logg                   [cgs units]
#     effective_temperature  [K]  
#     metallicity            [dex Fe/H]
#     photometric_filter     ['u', 'v', 'b', 'y', 'U', 'B', 'V', 'R', 'I', 'J', 'H', 'K']
#     stellar_model          ['ATLAS' (default), 'PHOENIX' (limited cases)]
# 
#     for examples, for an observation of HD209458b in the optical we can have:

photometric_filter = 'V'

limb_darkening_coefficients = plc.clablimb('claret', logg, effective_temperature, metallicity, photometric_filter)

print(limb_darkening_coefficients)


# plc.transit_projected_distance(period, sma_over_rs, eccentricity, inclination, periastron, mid_time, time_array)
#
#     returns the position vector of the planet in a coordinate system with
#     the parent star at (x,y,z) = (0,0,0), the observer at (x,y,z) = (+inf,0,0)
#     and the z-axis perpendicular to the plane of reference
#     given the orbital parameters of the planet
#     note: when the planet of further than the star, the values are negative
#
#     period                 [days]       float
#     sma_over_rs            [no units]   float       : semi-major axis over the stellar radius
#     eccentricity           [no units]   float
#     inclination            [degrees]    float
#     periastron             [degrees]    float
#     mid_time               [days]       float
#     time_array             [days]       numpy array
#
#     for examples, for a transit observation of HD209458b
#     from 2 hours before the mid-transit to 2 hours after the mid-transit
#     and one exposure every minute, we can have:

time_array = np.arange(mid_time - 200.0 / 24.0, mid_time + 200.0 / 24.0, 1.0 / 60 / 24)

z_over_rs = plc.transit_projected_distance(period, sma_over_rs, eccentricity, inclination, periastron,
                                           mid_time, time_array)

plt.plot(time_array, z_over_rs, 'ko', ms=3)
plt.axhline(1, color='k', ls='--')
plt.text(0.5 * (plt.xlim()[1] + plt.xlim()[0]), 0.99, 'transit', ha='center', va='top')
plt.xlabel('time (days)')
plt.ylabel('projected distance (stellar radii)')
plt.show()


# plc.transit_flux_drop(method, limb_darkening_coefficients, rp_over_rs, z_over_rs, precision=3)
# 
#     returns the observed stellar flux as a function of time
#     i.e. the transit light-curve, inputs:
#
#     method                      ['claret', 'quad', 'sqrt', 'linear']
#     limb_darkening_coefficients [no units]   array-like
#     z_over_rs                   [no units]   float              : drojected distance over the stellar radius
#     rp_over_rs                  [no units]   float              : planetary radius over the stellar radius
#     precision                   [0 - 6, default is 3]           : precision level for the numerical integration
#
#     for examples, for a transit observation of HD209458b in the optical
#     from 2 hours before the mid-transit to 2 hours after the midtransit
#     and one exposure every minute, we can have:

transit_light_curve = plc.transit_flux_drop(
    # method='linear',
    # limb_darkening_coefficients = [0.60],
    # method='quad',
    # limb_darkening_coefficients = [0.33, 0.36],
    # method='sqrt',
    # limb_darkening_coefficients = [0.04, 0.82],
    method='claret',
    limb_darkening_coefficients=limb_darkening_coefficients,
    z_over_rs=z_over_rs,
    rp_over_rs=rp_over_rs,
    precision=3,
    # precision=6
)

plt.plot(time_array, transit_light_curve, 'ko', ms=3)
plt.ylim(plt.ylim()[0], 1.001)
plt.xlabel('time (days)')
plt.ylabel('observed flux (%)')
plt.show()


# create some simulated data

time_array = np.arange(mid_time + 10.0 * period - 0.11, mid_time + 10.0 * period + 0.11, 2.0 / 60.0 / 24.0)
flux_array = plc.transit_integrated('claret', limb_darkening_coefficients, rp_over_rs, period, sma_over_rs, eccentricity,
                                    inclination, periastron, mid_time, time_array, 120, 120, precision=6)
systematics_array = 1.2 * (1 + 0.013 * (time_array - time_array[0]) + 0.03 * ((time_array - time_array[0]) ** 2))
error_array = np.random.normal(0, 0.002, len(time_array))

time0 = time_array
flux0 = flux_array * systematics_array + error_array
error0 = np.ones_like(error_array) * np.std(error_array)

time_array = np.arange(mid_time + 25.0 * period - 0.13, mid_time + 25.0 * period + 0.13, 2.0 / 60.0 / 24.0)
flux_array = plc.transit_integrated('claret', limb_darkening_coefficients, rp_over_rs, period, sma_over_rs, eccentricity,
                                    inclination, periastron, mid_time, time_array, 120, 120, precision=6)
systematics_array = 3.6 * (1 - 0.02 * (time_array - time_array[0]) + 0.05 * ((time_array - time_array[0]) ** 2))
error_array = np.random.normal(0, 0.005, len(time_array))

time1 = time_array
flux1 = flux_array * systematics_array + error_array
error1 = np.ones_like(error_array) * np.std(error_array)

time_array = np.arange(mid_time + 31.0 * period - 0.115, mid_time + 31.0 * period + 0.115, 2.0 / 60.0 / 24.0)
flux_array = plc.transit_integrated('claret', limb_darkening_coefficients, rp_over_rs, period, sma_over_rs, eccentricity,
                                    inclination, periastron, mid_time, time_array, 120, 120, precision=6)
systematics_array = 0.75 * (1 - 0.01 * (time_array - time_array[0]) + 0.0001 * ((time_array - time_array[0]) ** 2))
error_array = np.random.normal(0, 0.0009, len(time_array))

time2 = time_array
flux2 = flux_array * systematics_array + error_array
error2 = np.ones_like(error_array) * np.std(error_array)

plt.subplot(2, 2, 1)
plt.plot(time0, flux0, 'ko', ms=3)
plt.subplot(2, 2, 2)
plt.plot(time1, flux1, 'ko', ms=3)
plt.subplot(2, 2, 3)
plt.plot(time2, flux2, 'ko', ms=3)
plt.show()

fitting = plc.TransitAndPolyFitting(
    data=[[time0, flux0, error0], [time1, flux1, error1], [time2, flux2, error2]],
    # method='linear',
    # limb_darkening_coefficients=[0.60],
    # method='quad',
    # limb_darkening_coefficients=[0.33, 0.36],
    # method='sqrt',
    # limb_darkening_coefficients=[2.2223, -1.4052],
    method='claret',
    limb_darkening_coefficients=limb_darkening_coefficients,
    # method='linear',
    # limb_darkening_coefficients='fit',
    # method='quad',
    # limb_darkening_coefficients='fit',
    # method='sqrt',
    # limb_darkening_coefficients='fit',
    # method='claret',
    # limb_darkening_coefficients='fit',
    rp_over_rs=rp_over_rs,
    period=period,
    sma_over_rs=sma_over_rs,
    eccentricity=eccentricity,
    inclination=inclination,
    periastron=periastron,
    mid_time=mid_time,
    iterations=150000,
    walkers=50,
    burn=50000,
    time_factor=2,
    exp_time=120,
    precision=3,
    # precision=6,
    fit_first_order=True,
    # fit_first_order=False,
    fit_second_order=True,
    # fit_second_order=False,
    fit_rp_over_rs=[rp_over_rs / 2.0, rp_over_rs * 2.0],
    # fit_rp_over_rs=False,
    fit_period=[period / 2.0, period * 2.0],
    # fit_period=False,
    fit_sma_over_rs=[sma_over_rs / 2, sma_over_rs * 2.0],
    # fit_sma_over_rs=False,
    fit_inclination=[70, 90],
    # fit_inclination=False,
    fit_mid_time=[mid_time - 0.1, mid_time + 0.1],
    # fit_mid_time=False
)

fitting.run_mcmc()

fitting.save_all('simulation_data_base.pickle')
fitting.save_results('simulation_results.txt')
fitting.plot_corner('simulation_correlations.pdf')
fitting.plot_traces('simulation_traces.pdf')
fitting.plot_models('simulation_full_models.pdf')
fitting.plot_detrended_models('simulation_detrended_models.pdf')
