# PyLightcurve

<img src="https://github.com/ucl-exoplanets/pylightcurve/blob/master/logo.jpg" width="20%">

A python package for modeling and analysing transit light-curves.

* Easy search for parameters of current exoplanets.
* Calculation of limb darkening coefficients.
* Calculation of exoplanetary orbits.
* Calculation of transit models.
* Flexible fitting of transit light-curves.

This module makes use of:

* [exodata](https://github.com/ryanvarley/ExoData), [Varley (2016)](http://www.sciencedirect.com/science/article/pii/S0010465516301254)
* [emcee](https://github.com/dfm/emcee), [Foreman-Mackey et al. (2013)](http://iopscience.iop.org/article/10.1086/670067)


## Installation

For the latest stable version 2.3.2, open a terminal and type `pip install pylightcurve`.

For the new (under development) version 3.0.0, download this repo, unzip and type `python setup.py install`.


## Usage

The code in the examples below can be found in the example/example.py file in this repo.

	>>> import pylightcurve as plc
	>>> import matplotlib.pyplot as plt
	>>> import numpy as np
	
	>>> plt.ion()


##### plc.find_oec_parameters(target)

Returns the following stellar and transit parameters: planet oec name, logarithmic stellar surface gravity, stellar 
effective temperature, stellar metallicity, planetary radius relative to teh stellar radius, planetary bolometric 
emission relative o the stellar bolometric emision, orbital period, orbital semi-major axis relative to the stellar 
radius, orbital eccentricity, orbital inclination, orbital argument of periastron, transit mid-time.

Note: The database is automatically updated on a daily basis if internet access is available.

- target  
Name of the planet (str). 

For example, we can find the parameters of HD209458b:
	
	>>> (planet, stellar_logg, stellar_temperature, stellar_metallicity, rp_over_rs, fp_over_fs, 
	     period, sma_over_rs, eccentricity, inclination, periastron, mid_time) = plc.find_oec_parameters('hd209458b')

	>>> print (planet, stellar_logg, stellar_temperature, stellar_metallicity, rp_over_rs, fp_over_fs, 
	           period, sma_over_rs, eccentricity, inclination, periastron, mid_time)
	('HD 209458 b', 4.375254713815686, 6075.0, 0.02, 0.12035170971037652, 5.1956599618667065e-05, 3.52474859, 
	 8.8593557009493, 0.0004, 86.59, 0.0, 2451370.048)


##### plc.clablimb(method, stellar_logg, stellar_temperature, stellar_metallicity, photometric_filter, stellar_model='ATLAS')

Returns a list of limb darkening coefficients.

- method  
Limb darkening law (str, 'claret' is the only one currently supported).

- stellar_logg  
Logarithmic stellar surface gravity (float, in cm/s/s).

- stellar_temperature  
Stellar effective temperature (float, in Kelvin).

- stellar_metallicity  
Stellar metallicity (float, dex Fe/H).

- photometric_filter  
Photometric band of the observation (str, available filters: 'B', 'C', 'H', 'I', 'J', 'K', 'Kp', 'R', 'S1', 'S2', 
'S3', 'S4', 'U', 'V', 'b', 'g,', 'i,', 'r,', 'u', 'u,', 'v', 'y', 'z,').

For example, we can calculate the limb darkening coefficients for the claret law for HD209458b in the optical band:

	>>> limb_darkening_coefficients = plc.clablimb('claret', stellar_logg, stellar_temperature, 
	                                               stellar_metallicity, 'V')

	>>> print limb_darkening_coefficients
	[ 0.38606363  0.58637444 -0.19471546 -0.00559748]


##### plc.exoplanet_orbit(period, sma_over_rs, eccentricity, inclination, periastron, mid_time, time_array)

Returns the position vector of the planet in a coordinate system with the parent star at (x,y,z) = (0,0,0), the 
observer at (x,y,z) = (+inf,0,0) and the z-axis perpendicular to the plane of reference.

- period  
Orbital period (float, in days).

- sma_over_rs  
Orbital semi-major axis relative to the stellar radius (float, no units).

- eccentricity  
Orbital eccentricity (float, no units).

- inclination  
Orbital inclination (float, in degrees).

- periastron  
Orbital argument of periastron (float, in degrees).

- mid_time  
Transit mid-time (float, in days).

- time_array  
A time sequence (numpy array, in days).

For example, we can calculate the position vector of HD209458b from 2 hours before the mid-transit to 2 hours after 
the mid-transit with a frequency of 1 point per minute:

    >>> time_array = np.arange(mid_time - 2.0 / 24.0, mid_time + 2.0 / 24.0, 1.0 / 60.0 / 24.0)

    >>> (position_x, position_y, position_z) = plc.exoplanet_orbit(period, sma_over_rs, eccentricity, inclination, 
                                                                   periastron, mid_time, time_array)
                                        
    >>> plt.subplot(3,1,1)
    >>> plt.plot(time_array, position_x, 'ko', ms=3)
    >>> plt.ylabel('x (R star)')
    >>> plt.subplot(3,1,2)
    >>> plt.plot(time_array, position_y, 'ko', ms=3)
    >>> plt.ylabel('y (R star)')
    >>> plt.subplot(3,1,3)
    >>> plt.plot(time_array, position_z, 'ko', ms=3)
    >>> plt.ylabel('z (R star)')
    >>> plt.xlabel('time (days)')


##### plc.transit_projected_distance(period, sma_over_rs, eccentricity, inclination, periastron, mid_time, time_array, precision=3)

Returns the projected distance between the planet and its parent star. When the planet is further than the star, 
the values returned are negative.

- period  
Orbital period (float, in days).

- sma_over_rs  
Orbital semi-major axis relative to the stellar radius (float, no units).

- eccentricity  
Orbital eccentricity (float, no units).

- inclination  
Orbital inclination (float, in degrees).

- periastron  
Orbital argument of periastron (float, in degrees).

- mid_time  
Transit mid-time (float, in days).

- time_array  
A time sequence (numpy array, in days).

For example, we can calculate the projected distance of HD209458b from its host star from 2 hours before the 
mid-transit to 2 hours after the mid-transit with a frequency of 1 point per minute:
    
    >>> z_over_rs = plc.transit_projected_distance(period, sma_over_rs, eccentricity, inclination, periastron,
                                                   mid_time, time_array)

    >>> plt.plot(time_array, z_over_rs, 'ko', ms=3)
    >>> plt.axhline(1, color='k', ls='--')
    >>> plt.text(0.5 * (plt.xlim()[1] + plt.xlim()[0]), 0.99, 'transit', ha='center', va='top')
    >>> plt.xlabel('time (days)')
    >>> plt.ylabel('projected distance (R star)')


##### plc.transit_flux_drop(method, limb_darkening_coefficients, rp_over_rs, z_over_rs, precision=3)

Returns the observed stellar flux as a function of time - i.e. the transit light-curve.

- method  
Limb darkening law (str, available methods: 'claret', 'quad', 'sqrt' or 'linear').

- limb_darkening_coefficients  
A list containing the limb darkening coefficients. The list should contain 1 element if the method used is the 
'linear', 2 if the method used is the 'quad' or teh 'sqrt', and 4 if the method used is the 'claret'.

- rp_over_rs  
Planetary radius relative to the stellar radius (float, no units)

- z_over_rs  
Projected distance between the planet and its parent star relative to the stellar radius (numpy array, no units).

- precision  
The level of the numerical precision for the calculation (int, 0 to 6, default value is 3).

For example, we can calculate the transit light-curve of HD209458b from 2 hours before the mid-transit to 2 hours 
after the mid-transit with a frequency of 1 point per minute:

    >>> flux_drop = plc.transit_flux_drop('claret', limb_darkening_coefficients, rp_over_rs, z_over_rs)

    >>> plt.plot(time_array, flux_drop, 'ko', ms=3)
    >>> plt.ylim(plt.ylim()[0], 1.001)
    >>> plt.xlabel('time (days)')
    >>> plt.ylabel('observed flux (%)')


##### plc.transit(method, limb_darkening_coefficients, rp_over_rs, period, sma_over_rs, eccentricity, inclination, periastron, mid_time, time_array, precision=3)

Returns the transit light-curve, directly from the orbital parameters.

- method  
Limb darkening law (str, available methods: 'claret', 'quad', 'sqrt' or 'linear').

- limb_darkening_coefficients  
A list containing the limb darkening coefficients. The list should contain 1 element if the method used is the 
'linear', 2 if the method used is the 'quad' or teh 'sqrt', and 4 if the method used is the 'claret'.

- rp_over_rs  
Planetary radius relative to the stellar radius (float, no units)

- period  
Orbital period (float, in days).

- sma_over_rs  
Orbital semi-major axis relative to the stellar radius (float, no units).

- eccentricity  
Orbital eccentricity (float, no units).

- inclination  
Orbital inclination (float, in degrees).

- periastron  
Orbital argument of periastron (float, in degrees).

- mid_time  
Transit mid-time (float, in days).

- time_array  
A time sequence (numpy array, in days).

- precision  
The level of the numerical precision for the calculation (int, 0 to 6, default value is 3).

For example, we can calculate the transit light-curve of HD209458b from 2 hours before the mid-transit to 2 hours after 
the mid-transit with a frequency of 1 point per minute:

    >>> transit_light_curve = plc.transit('claret', limb_darkening_coefficients, rp_over_rs, period, sma_over_rs, 
                                          eccentricity, inclination, periastron, mid_time, time_array)

    >>> plt.plot(time_array, transit_light_curve, 'ko', ms=3)
    >>> plt.ylim(plt.ylim()[0], 1.001)
    >>> plt.xlabel('time (days)')
    >>> plt.ylabel('observed flux (%)')


##### plc.transit_integrated(method, limb_darkening_coefficients, rp_over_rs, period, sma_over_rs, eccentricity, inclination, periastron, mid_time, time_array, exp_time, time_factor, precision=3)

Returns the exposure-integrated transit light-curve, directly from the orbital parameters.

- method  
Limb darkening law (str, available methods: 'claret', 'quad', 'sqrt' or 'linear').

- limb_darkening_coefficients  
A list containing the limb darkening coefficients. The list should contain 1 element if the method used is the 
'linear', 2 if the method used is the 'quad' or teh 'sqrt', and 4 if the method used is the 'claret'.

- rp_over_rs  
Planetary radius relative to the stellar radius (float, no units)

- period  
Orbital period (float, in days).

- sma_over_rs  
Orbital semi-major axis relative to the stellar radius (float, no units).

- eccentricity  
Orbital eccentricity (float, no units).

- inclination  
Orbital inclination (float, in degrees).

- periastron  
Orbital argument of periastron (float, in degrees).

- mid_time  
Transit mid-time (float, in days).

- time_array  
A time sequence (numpy array, in days).

- exp_time  
Exposure time (float, in seconds).

- time_factor  
Number of sub-exposures to be calculated per exposure (int, no units).

- precision  
The level of the numerical precision for the calculation (int, 0 to 6, default value is 3).

For example, we can calculate the transit light-curve of HD209458b from 2 hours before the mid-transit to 2 hours after 
the mid-transit with a frequency of 1 point per minute, assuming an exposure time of 30 seconds which is divided into 
10 sub-exposures:

    >>> transit_light_curve = plc.transit_integrated('claret', limb_darkening_coefficients, rp_over_rs, period, 
                                                     sma_over_rs, eccentricity, inclination, periastron, mid_time, 
                                                     time_array, 30, 10)

    >>> plt.plot(time_array, transit_light_curve, 'ko', ms=3)
    >>> plt.ylim(plt.ylim()[0], 1.001)
    >>> plt.xlabel('time (days)')
    >>> plt.ylabel('observed flux (%)')


##### plc.TransitAndPolyFitting(data, method, limb_darkening_coefficients, rp_over_rs, period, sma_over_rs, eccentricity, inclination, periastron, mid_time, iterations, walkers, burn, precision=3, exp_time=0, time_factor=1, fit_first_order=False, fit_second_order=False, fit_rp_over_rs=False, fit_period=False, fit_sma_over_rs=False, fit_eccentricity=False, fit_inclination=False, fit_periastron=False, fit_mid_time=False, counter=True, counter_window=False):

Offers a range of options for fitting observed transit light-curves, simultaneously with a second-order polynomial 
de-trending function.

- data  
A list containing the input data sets. Each element in the list is a list of 3 arrays, 
representing the time (in Heliocentric Julian Date), the stellar flux and the uncertainty in the stellar flux
example: `data=[[time_0, flux_0, error_0], [time_1, flux_1, error_1], [time_2, flux_2, error_2]]`

- method  
Limb darkening law (str, available methods: 'claret', 'quad', 'sqrt' or 'linear').

- limb_darkening_coefficients  
A list containing the limb darkening coefficients. The list should contain 1 element if the method used is the 
'linear', 2 if the method used is the 'quad' or teh 'sqrt', and 4 if the method used is the 'claret'.
To fit for the limb darkening coefficients set `limb_darkening_coefficients='fit'`.

- rp_over_rs  
Initial value for the planetary radius relative to the stellar radius (float, no units)

- period  
Initial value for the orbital period (float, in days).

- sma_over_rs  
Initial value for the orbital semi-major axis relative to the stellar radius (float, no units).

- eccentricity  
Initial value for the orbital eccentricity (float, no units).

- inclination  
Initial value for the orbital inclination (float, in degrees).

- periastron  
Initial value for the orbital argument of periastron (float, in degrees).

- mid_time  
Initial value for the transit mid-time (float, in days).

- time_array  
A time sequence (numpy array, in days).

- iterations  
Number of total mcmc iterations (int, no units).

- walkers  
Number of walkers, as defined in the emcee package (int, no units).

- burn  
Number of iterations to be excluded from the beginning of the chains (int, no units).

- precision  
The level of the numerical precision for the calculation (int, 0 to 6, default value is 3).

- exp_time  
Exposure time (float, in seconds, default value is 0).

- time_factor  
Number of sub-exposures to be calculated per exposure (int, no units, default value is 1).

- fit_first_order  
Flag for including a first order time-dependent de-trending factor (bool, default value is False).

- fit_second_order  
Flag for including a second order time-dependent de-trending factor (bool, default value is False).

- fit_rp_over_rs  
A 2-element list containing the lower and upper limits for fitting the planetary radius relative to the stellar radius.
To avoid fitting set `fit_rp_over_rs=False`. Default value is False.

- fit_period  
A 2-element list containing the lower and upper limits for fitting the orbital period.
To avoid fitting set `fit_rp_over_rs=False`. Default value is False.

- fit_sma_over_rs  
A 2-element list containing the lower and upper limits for fitting the orbital semi-major axis relative to the stellar radius.
To avoid fitting set `fit_rp_over_rs=False`. Default value is False.

- fit_eccentricity  
A 2-element list containing the lower and upper limits for fitting the orbital eccentricity.
To avoid fitting set `fit_rp_over_rs=False`. Default value is False.

- fit_inclination  
A 2-element list containing the lower and upper limits for fitting the orbital inclination.
To avoid fitting set `fit_rp_over_rs=False`. Default value is False.

- fit_periastron  
A 2-element list containing the lower and upper limits for fitting the orbital argument of periastron.
To avoid fitting set `fit_rp_over_rs=False`. Default value is False.

- fit_mid_time  
A 2-element list containing the lower and upper limits for fitting the the transit mid-time.
To avoid fitting set `fit_rp_over_rs=False`. Default value is False.

- counter  
Flag for printing a counter of the completed iterations (bool, default value is True).

- counter_window=False  
Flag for showing a counter of the completed iterations in an additional Tk window (bool, default value is False).

##### plc.TransitAndPolyFitting methods:

###### .run_mcmc()

Sets up and runs the mcmc.

###### .save_all(export_file)

Saves all the mcmc results (including the chains) in the form of a pickle file.

- export_file  
File to be created (str).

###### .save_results(export_file)

Saves the final values and uncertainties of the fitted parameters in the form of a txt file.

- export_file  
File to be created (str).

###### .plot_corner(export_file)

Plots the correlations between the fitted parameters.

- export_file  
File to be created (str).

###### .plot_traces(export_file)

Plots the mcmc chains of the fitted parameters.

- export_file  
File to be created (str).

###### .plot_models(export_file)

Plots the original data and the full model fitted. A prefix is added to indicate the different data sets (set_1, set_2, etc.).

- export_file  
File to be created (str).

###### .plot_detrended_models(export_file)

Plots the data corrected by the de-trending function and the transit model fitted. A prefix is added to indicate the 
different data sets (set_1, set_2, etc.).

- export_file  
File to be created (str).

In the following example we will create 3 simulated observations of HD209458b, with an exposure time of 2 minutes and 
additional second-order time-dependent systematics, and fit them using the **plc.TransitAndPolyFitting** class. To 
avoid an extremely slow process, we will use a time factor of 2 for the fitting, while we will use a time factor of 
120 to create the simulated observations:

    >>> time_array = np.arange(mid_time + 10.0 * period - 0.11, mid_time + 10.0 * period + 0.11, 2.0 / 60.0 / 24.0)
    >>> flux_array = plc.transit_integrated('claret', limb_darkening_coefficients, rp_over_rs, period, sma_over_rs, 
                                            eccentricity, inclination, periastron, mid_time, time_array, 120, 120, 
                                            precision=6)
    >>> systematics_array = 1.2 * (1 + 0.013 * (time_array - time_array[0]) + 
                                   0.03 * ((time_array - time_array[0]) ** 2))
    >>> error_array = np.random.normal(0, 0.002, len(time_array))

    >>> time0 = time_array
    >>> flux0 = flux_array * systematics_array + error_array
    >>> error0 = np.ones_like(error_array) * np.std(error_array)

    >>> time_array = np.arange(mid_time + 25.0 * period - 0.13, mid_time + 25.0 * period + 0.13, 2.0 / 60.0 / 24.0)
    >>> flux_array = plc.transit_integrated('claret', limb_darkening_coefficients, rp_over_rs, period, sma_over_rs, 
                                            eccentricity, inclination, periastron, mid_time, time_array, 120, 120, 
                                            precision=6)
    >>> systematics_array = 3.6 * (1 - 0.02 * (time_array - time_array[0]) + 
                                   0.05 * ((time_array - time_array[0]) ** 2))
    >>> error_array = np.random.normal(0, 0.005, len(time_array))

    >>> time1 = time_array
    >>> flux1 = flux_array * systematics_array + error_array
    >>> error1 = np.ones_like(error_array) * np.std(error_array)

    >>> time_array = np.arange(mid_time + 31.0 * period - 0.115, mid_time + 31.0 * period + 0.115, 2.0 / 60.0 / 24.0)
    >>> flux_array = plc.transit_integrated('claret', limb_darkening_coefficients, rp_over_rs, period, sma_over_rs, 
                                            eccentricity, inclination, periastron, mid_time, time_array, 120, 120, 
                                            precision=6)
    >>> systematics_array = 0.75 * (1 - 0.01 * (time_array - time_array[0]) 
                                    + 0.0001 * ((time_array - time_array[0]) ** 2))
    >>> error_array = np.random.normal(0, 0.0009, len(time_array))

    >>> time2 = time_array
    >>> flux2 = flux_array * systematics_array + error_array
    >>> error2 = np.ones_like(error_array) * np.std(error_array)

    >>> plt.close('all')
    >>> plt.ioff()

    >>> fitting = plc.TransitAndPolyFitting(
            data=[[time0, flux0, error0], [time1, flux1, error1], [time2, flux2, error2]],
            method='claret',
            limb_darkening_coefficients=limb_darkening_coefficients,
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
            precision=3,
            time_factor=2,
            exp_time=120,
            fit_first_order=True,
            fit_second_order=True,
            fit_rp_over_rs=[rp_over_rs / 2.0, rp_over_rs * 2.0],
            fit_period=[period / 2.0, period * 2.0],
            fit_sma_over_rs=[sma_over_rs / 2, sma_over_rs * 2.0],
            fit_inclination=[70, 90],
            fit_mid_time=[mid_time - 0.1, mid_time + 0.1])

    >>> fitting.run_mcmc()

    >>> fitting.save_all('simulation_data_base.pickle')
    >>> fitting.save_results('simulation_results.txt')
    >>> fitting.plot_corner('simulation_correlations.pdf')
    >>> fitting.plot_traces('simulation_traces.pdf')
    >>> fitting.plot_models('simulation_full_models.pdf')
    >>> fitting.plot_detrended_models('simulation_detrended_models.pdf')


## Licence

MIT License

Copyright (c) 2016-2019 Angelos Tsiaras, Konstantinos Karpouzas and Ryan Varley

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
