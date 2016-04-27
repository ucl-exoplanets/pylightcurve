import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pylightcurve as plc


# The main functions of pyligthcurve are:
# plc.limb_darkening
# plc.position_vector
# plc.transit_linear 
# plc.transit_quad 
# plc.transit_sqrt
# plc.transit_claret
# plc.mcmc_transit





# plc.limb_darkening(metallicity, effective_temperature, logg, photometric_filter)
# 	
#     returns the limb darkening coefficients 
#     (4 coefficients non-linear low from Claret 2000)
#     given the stellar parameters and the photometric filter used for the observation
# 
#     metallicity            [dex Fe/H]
#     effective_temperature  [K] 
#     logg                   [cgs units] 
#     photometric_filter     ['u', 'v', 'b', 'y', 'U', 'B', 'V', 'R', 'I', 'J', 'H', 'K']
# 
#     for example, for an observation of HD209458b in the optical we can have:

metallicity           = 0.0
effective_temperature = 6092
logg                  = 4.28
photometric_filter    = 'V'

limb_darkening_coefficients = plc.limb_darkening(metallicity, effective_temperature, 
                                                 logg, photometric_filter)
print limb_darkening_coefficients

#     however, for the functions that require limb darkening coefficients, 
#     the user can also provide a 4-item array-like object,
#     with coefficients caplculated independently

limb_darkening_coefficients = (0.608402, -0.206180, 0.262286, -0.133088)
print limb_darkening_coefficients





# plc.position_vector(period, sma_over_rs, eccentricity, inclination, 
#                     periastron, mid_time, time_array)
# 
#     returns the position vector of the planet in a coordinate system with
#     the parent star at (x,y,z) = (0,0,0), the observer at (x,y,z) = (+inf,0,0) 
#     and the z-axis perpendicular to the plane of reference
#     given the orbital parameters of the planet
# 
#     period                 [days]       float
#     sma_over_rs            [no units]   float       : semi-major axis over the stellar radius
#     eccentricity           [no units]   float
#     inclination            [degrees]    float
#     periastron             [degrees]    float
#     mid_time               [days]       float
#     time_array             [days]       numpy array
# 
#     for example, for HD209458b we can have:

period       = 3.52474859
sma_over_rs  = 8.76
eccentricity = 0.0
inclination  = 86.71
periastron   = 0.0
mid_time     = 2452826.625521
time_array   = np.arange(2452826.625521 - 2, 2452826.625521 + 2, 0.0005)

position_vector = plc.position_vector(period, sma_over_rs, eccentricity, inclination, 
                                      periastron, mid_time, time_array)
plc.pylightcurve_tools.plot_trajectory(position_vector)





# plc.transit_*(limb_darkening_coefficients, rp_over_rs, period, sma_over_rs, eccentricity, 
#               inclination, periastron, mid_time, time_array)
# 
#     is a combination of the above functions, as it uses the 
#     orbital parameters of the planet to calculate the projected distance 
#     as a function of time and return the observed stellar flux as a function of time
#     i.e. the transit light-curve, inputs:
# 
#     limb_darkening_coefficients [no units]   4-item array-like
#     rp_over_rs                  [no units]   float              : planetary radius over the stellar radius
#     period                      [days]       float
#     sma_over_rs                 [no units]   float              : semi-major axis over the stellar radius
#     eccentricity                [no units]   float
#     inclination                 [degrees]    float
#     periastron                  [degrees]    float
#     mid_time                    [days]       float
#     time_array                  [days]       numpy array
# 
#     for example, for HD209458b we can have:

limb_darkening_coefficients = (0.608402, -0.206180, 0.262286, -0.133088)
rp_over_rs                  = 0.120859422
period                      = 3.52474859
sma_over_rs                 = 8.76
eccentricity                = 0.0
inclination                 = 86.71
periastron                  = 0.0
mid_time                    = 2452826.625521
time_array                  = np.arange(2452826.625521 - 0.11, 2452826.625521 + 0.11, 1.0/60.0/60/24)

transit_light_curve = plc.transit((1.1584), rp_over_rs, period,
                                  sma_over_rs, eccentricity, inclination, periastron, 
                                  mid_time , time_array, method='linear')
plt.plot(time_array, transit_light_curve, c='b', lw=2)

transit_light_curve = plc.transit((1.9571, -0.9669), rp_over_rs, period,
                                  sma_over_rs, eccentricity, inclination, periastron, 
                                  mid_time , time_array, method='quad')
plt.plot(time_array, transit_light_curve, c='r', lw=2)

transit_light_curve = plc.transit((2.2223, -1.4052), rp_over_rs, period,
                                  sma_over_rs, eccentricity, inclination, periastron, 
                                  mid_time , time_array, method='sqrt')
plt.plot(time_array, transit_light_curve, c='g', lw=2)

transit_light_curve = plc.transit((-0.1725, 0.5323, -0.5748, 1.1889), rp_over_rs, period,
                                  sma_over_rs, eccentricity, inclination, periastron, 
                                  mid_time , time_array, method='claret')
plt.plot(time_array, transit_light_curve, c='c', lw=2)


plt.ylim(plt.ylim()[0], 1.001)
plt.xlabel('time [days]')
plt.ylabel('observed flux [%]')
plt.show()


time0, flux0 = np.loadtxt('simulation_0.txt', unpack=True)
time1, flux1 = np.loadtxt('simulation_1.txt', unpack=True)
time2, flux2 = np.loadtxt('simulation_2.txt', unpack=True)

plc.mcmc_transit(method = 'claret',
                 limb_darkening_coefficients = (0.608402, -0.206180, 0.262286, -0.133088),
                 rp_over_rs                  = 0.120859422,
                 period                      = 3.52474859,
                 sma_over_rs                 = 8.76,
                 eccentricity                = 0.0,
                 inclination                 = 86.71,
                 periastron                  = 0.0,
                 mid_time                    = 2452826.625521,
                 data                        = [[time0, flux0], [time1, flux1], [time2, flux2]],
                 fit_rp_over_rs              = [0.1, 0.2],
                 iterations                  = 230000,
                 burn                        = 30000,
                 directory                   = 'simulation_results',
                 detrend_order               =	2,
                 fit_period                  = [3,4],
                 fit_sma_over_rs             = [3,15],
                 fit_eccentricity            = None,
                 fit_inclination             = [70,90],
                 fit_periastron              = None,
                 fit_mid_time                = [2452826.625521 - 0.1,2452826.625521 + 0.1]
                )
