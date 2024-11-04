[![Build Python 3.8](https://github.com/ucl-exoplanets/pylightcurve/actions/workflows/build-python-3.8.yml/badge.svg)](https://github.com/ucl-exoplanets/pylightcurve/actions/workflows/build-python-3.8.yml)
[![Build Python 3.9](https://github.com/ucl-exoplanets/pylightcurve/actions/workflows/build-python-3.9.yml/badge.svg)](https://github.com/ucl-exoplanets/pylightcurve/actions/workflows/build-python-3.9.yml)
[![Build Python 3.10](https://github.com/ucl-exoplanets/pylightcurve/actions/workflows/build-python-3.10.yml/badge.svg)](https://github.com/ucl-exoplanets/pylightcurve/actions/workflows/build-python-3.10.yml)
[![Build Python 3.11](https://github.com/ucl-exoplanets/pylightcurve/actions/workflows/build-python-3.11.yml/badge.svg)](https://github.com/ucl-exoplanets/pylightcurve/actions/workflows/build-python-3.11.yml)


[![codecov](https://codecov.io/gh/ucl-exoplanets/pylightcurve/branch/master/graph/badge.svg?)](https://codecov.io/gh/ucl-exoplanets/pylightcurve)

[![Downloads](https://pepy.tech/badge/pylightcurve)](https://pepy.tech/project/pylightcurve)

# PyLightcurve 

A python package for analysing exoplanet light-curves.

<img src="https://github.com/ucl-exoplanets/pylightcurve/blob/master/logo.png" width="25%">

In PyLightcurve you will find tools for:

* Calculation of limb darkening coefficients.
* Calculation of exoplanetary orbits.
* Calculation of exoplanet transit and eclipse properties.
* Exoplanet transit and eclipse modeling.
* Flexible fitting of multi-epoch and multi-colour exoplanet transit and eclipse 
light-curves.
* Transformations between different angle and timing systems.

Developed by [Angelos Tsiaras](www.angelostsiaras.com)

PyLightcurve makes use of the following packages:

* [Emcee](https://github.com/dfm/emcee), [Foreman-Mackey et al. 2013](http://iopscience.iop.org/article/10.1086/670067)
* [ExoTETHyS](https://github.com/ucl-exoplanets/ExoTETHyS), [Morello et al. 2020](https://iopscience.iop.org/article/10.3847/1538-3881/ab63dc)
* [Matplotlib](https://matplotlib.org), [Hunter 2007](https://ieeexplore.ieee.org/document/4160265)
* [Numpy](https://numpy.org), [Oliphant 2006](https://archive.org/details/NumPyBook)
* [SciPy](https://www.scipy.org), [Virtanen et al. 2020](https://www.nature.com/articles/s41592-019-0686-2)
* [Astropy](https://www.astropy.org), [Astropy Collaboration 2013](https://www.aanda.org/articles/aa/abs/2013/10/aa22068-13/aa22068-13.html)

and of the following catalogues:

* Exoplanet Characterisation Catalogue, developed as part of the [ExoClock Project](www.exoclock.space), [Kokori et al. 2020](https://ui.adsabs.harvard.edu/abs/2020arXiv201207478K/abstract)

If you are using PyLightcurve for your research please cite all the above references individually and also:

[Tsiaras et al. 2016](https://ui.adsabs.harvard.edu/abs/2016ApJ...832..202T/abstract)


# Updates in version 4.0

#### Deprecation notice

PyLightcurve 4.0 no longer supports the use of the Open Exoplanet Catalogue (OEC), due to the large number 
of mistakes in the catalogue and the absence of parameters updates. OEC has been replaced by the 
Exoplanet Characterisation Catalogue, dedicated to transiting exoplanets. This catalogue does not contain all the 
exoplanets at the moment but will be continuously updated with more planets.

#### Multi-epoch and multi-color data fitting

Version 4.0 brings a new easy way of simultaneously fitting light curves from different 
sources (different flux and time formats, different exposure times), different filters
(multi-colour) and different epochs (multi-epoch). 

#### Eclipse model and eclipse data fitting

We now have the option of modeling eclipse data (finally!) and also tools to calculate the 
eclipse mid-time, eclipse depth and eclipse duration, based of the orbital characteristics 
of the planet. 

#### Limb-darkening coefficients

In this version, the default limb-darkening coefficients are calculated using the new 
package [ExoTETHyS](https://github.com/ucl-exoplanets/ExoTETHyS). This allows the 
calculation of the limb-darkening coefficients even for stars cooler than 3500 K. 

#### Angles and times

PyLighrcurve 4.0 includes flexible objects for converting different expressions of 
angles (e.g. '+47:12:34.05' to degrees) and of timing systems (e.g. HJD_UTC to BJD_TDB).

# Installation

Install PyLightcurve through pip:

```bash
pip install pylightcurve
```

... or download this repo, cd in it and use the setup.py file:
 
```bash
git clone https://github.com/ucl-exoplanets/pylightcurve
cd pylightcurve
python setup.py install
```

If you are looking for the previous version of PyLightcurve (v.3), check the relevant branch:

https://github.com/ucl-exoplanets/pylightcurve/tree/pylightcurve-3 


# Usage

## Exoplanets - using the plc.Planet object

In Pylightcurve 4.0.0 we can access all the old and new calculations related to 
exoplanets through the newly introduced plc.Planet object. These calculation aare:
- limb-darkening coefficients
- transit model
- transit data fitting
- transit duration
- transit depth
- planet-to-star flux ratio
- eclipse time
- eclipse model
- eclipse data fitting
- eclipse duration
- eclipse depth
- exoplanet orbital position
- planet-star projected distance

We start by importing Pylightcurve, together with Numpy and Matplotib.

```python
import pylightcurve as plc
import matplotlib.pyplot as plt
import numpy as np
```

We can then define a plc.Planet object for our favorite exoplanet, HD 209458 b.

```python
planet = plc.Planet(
    name='HD209458b', 
    
    ra = 330.795,                  # float values are assumed to be in degrees,
                                   # alternatively, you can provide a plc.Hours or plc.Degrees object
                                   # here it would be plc.Hours('22:03:10.7729')
    
    dec = 18.884,                  # float values are assumed to be in degrees,
                                   # alternatively, you can provide a plc.Hours or plc.Degrees object
                                   # here it would be plc.Degrees('+18:53:03.548')
    
    stellar_logg = 4.36,           # float, in log(cm/s^2)
    
    stellar_temperature = 6065.0,  # float, in Kelvin
    
    stellar_metallicity = 0.0,     # float, in dex(Fe/H) or dex(M/H)
    
    rp_over_rs = 0.12086,          # float, no units
    
    period = 3.5247486,            # float, in days
    
    sma_over_rs = 8.76,            # float, no units
    
    eccentricity = 0.0,            # float, no units
    
    inclination = 86.71,           # float values are assumed to be in degrees,
                                   # alternatively, you can provide a plc.Hours or plc.Degrees object
                                   # here it would be plc.Degrees(86.71)           
    
    periastron = 0.0,              # float values are assumed to be in degrees,
                                   # alternatively, you can provide a plc.Hours or plc.Degrees object
                                   # here it would be plc.Degrees(0.0)
    
    mid_time = 2452826.62928,      # float, in days
    
    mid_time_format = 'BJD_TDB',   # str, available formats are JD_UTC, MJD_UTC, HJD_UTC, HJD_TDB, BJD_UTC, BJD_TDB    
    
    ldc_method = 'claret',         # str, default = claret, the other methods are: linear, quad, sqrt
    
    ldc_stellar_model = 'phoenix', # str, default = phoenix, the other model is atlas
    
    albedo = 0.15,                 # float, default = 0.15, no units 
    
    emissivity = 1.0,              # float, default = 1.0, no units
)
```

We can quickly create a plc.Planet object based on catalague data as follows: 

```python
planet = plc.get_planet('hd209458b')
```

At the moment the data provided are base on the Exoplanet Characterisation Catalogue (ECC) developed as part of the 
ExoClock Project. The catalogue contains 370 objects and will gradually expand in future releases.

To retrieve a list of all the available planet names we can type:

```python
all_planets = plc.get_all_planets()
```


### Filters

There is a number of parameters for every planet that do depend on the observing filter
and are necessary to calculate the correct transit/eclipse models at different 
wavelengths and to allow the simultaneous analysis of multi-wavelength light-curves:
- the planet-to-star radius ratio (rp_over_rs)
- the limb-darkening coefficients (limb_darkening_coefficients)
- the planet-to-star flux ratio (fp_over_fs)

By default, a plc.Planet object contains the above parameters for a number of standard 
filters:
- clear
- luminance
- JOHNSON_U 
- JOHNSON_B 
- JOHNSON_V 
- COUSINS_R 
- COUSINS_I
- 2mass_j 
- 2mass_h
- 2mass_ks
- sdss_u 
- sdss_g 
- sdss_r 
- sdss_i
- sdss_z
- Kepler 
- TESS
- irac1 (available only for the atlas model, if phoenix is chosen, it will change 
    automatically to atlas, with more restrictions on the minimum temperature 
        available, 3500K)
- irac2 (available only for the atlas model, if phoenix is chosen, it will change 
    automatically to atlas, with more restrictions on the minimum temperature 
        available, 3500K)
- irac3 (available only for the atlas model, if phoenix is chosen, it will change 
    automatically to atlas, with more restrictions on the minimum temperature 
        available, 3500K)
- irac4 (available only for the atlas model, if phoenix is chosen, it will change 
    automatically to atlas, with more restrictions on the minimum temperature 
        available, 3500K)

In these default calculations:
- the rp_over_rs is equal to the value defined when creating the plc.Planet object,
- the limb-darkening coefficients are calculated using the filter response curves, together 
    with the stellar parameters, the ldc_method and the the ldc_stellar_model defined when 
    creating the plc.Planet object (more in the [ExoTETHyS](https://github.com/ucl-exoplanets/ExoTETHyS) package),
- the planet-to-star flux ratio (reflected + emmitted) is calculated using the albedo and 
    the emissivity defined when creating the plc.Planet object, assuming the star are 
    emitting as black bodies.

The default calculations can be accessed by typing:
```python
limb_darkening_coefficients = planet.filter('COUSINS_R').limb_darkening_coefficients
fp_over_fs = planet.filter('COUSINS_R').fp_over_fs
rp_over_rs = planet.filter('COUSINS_R').rp_over_rs # no difference from planet.rp_over_rs if we have not defined our own filter
```

Of course we can define additional filters (that can be accessed in the same way) or alter 
the parameters in the existing filters by typing:

```python
planet.add_filter('my_filter', rp_over_rs, ldc1, ldc2, ldc3, ldc4, fp_over_fs)
```

### Fitting Transit / Eclipse light-curves

For the purpose of this excersise, we will create some simulated data. In real life, you will 
provide your own light-curves with their characteristics. 

To model-fit observations, we first need to add them to the plc.Planet object. We should 
be careful to add either only transit data or only eclipse data.

```python

# first observation
time = np.arange(planet.mid_time - 0.1, planet.mid_time + 0.1, 0.001)

transit = planet.transit_integrated(time, time_format='BJD_TDB', exp_time=120, time_stamp = 'mid', filter_name='COUSINS_R', max_sub_exp_time=1)
systematics = 1.2 * (1 + 0.013 * (time - time[0]) + 0.03 * ((time - time[0]) ** 2))
error = np.random.normal(0, 0.002, len(time))
flux = transit * systematics + error

flux_unc = np.ones_like(error) * np.std(error)

planet.add_observation(
    time = time,                # the time vector of our observation
                                # np.array of float values 
    
    time_format = 'BJD_TDB',    # format in which our time vector is expressed
                                # str, available formats are: JD_UTC, MJD_UTC, HJD_UTC, HJD_TDB, BJD_UTC, BJD_TDB 

    exp_time = 120,             # exposure time of our time vector
                                # float, in seconds
        
    time_stamp = 'mid',         # exposure time stamp for our time vector (do the numbers refer to the exposure start, the mid-exposure, or the exposure end?)
                                # str, available stamps are: start, mid, end 
    
    flux = flux,                # the flux vector of our observation
                                # np.array of float values, 
    
    flux_unc = flux_unc,        # the flux-uncertainty vector of our observation
                                # np.array of float values, 
    
    flux_format = 'flux',       # format in which our flux and flux-uncertainty vectors are expressed
                                # str, available formats are: flux, mag
    
    filter_name = 'COUSINS_R'   # filter used for this observation 
                                # str, available filters are: all the default filters and those added manually by us
)

# second observation
time = np.arange(planet.mid_time - 0.05, planet.mid_time + 0.15, 0.001)

transit = planet.transit_integrated(time, time_format='HJD_UTC', exp_time=30, time_stamp = 'mid', filter_name='TESS', max_sub_exp_time=1)
systematics = 3.6 * (1 - 0.02 * (time - time[0]) + 0.05 * ((time - time[0]) ** 2))
error = np.random.normal(0, 0.0005, len(time))
flux = transit * systematics + error

flux_unc = np.ones_like(error) * np.std(error)

planet.add_observation(
    time = time,                                         
    time_format = 'HJD_UTC',                                                              
    exp_time = 30,
    time_stamp = 'mid',                                  
    flux = flux,                                   
    flux_unc = flux_unc,                         
    flux_format = 'flux',                       
    filter_name = 'TESS'                   
)
```

Once added, we can fit all the observations simultaneously using the following command:

```python
planet.transit_fitting(output_folder)
```

where the output_folder is the path where we want to save the results. There is a number of options 
when fitting the data:

```python
detrending_order            # default:2, instance:float, accepted values: 0, 1, 2
                            # indicates the order of a polynomial that will be fitted together with the observation 
                            # for de-trending purposes (every observation will be de-trended by a different polynomial

iterations                  # default:130000, instance:float
                            # indicates the number of MCMC iterations

walkers                     # default:200, instance:float
                            # indicates the number of MCMC wakers

burn_in                     # default:30000, instance:float
                            # indicates the number of MCMC burn-in

fit_rp_over_rs              # default:True, instance:bool
                            # indicates whether to fit for the rp_over_rs or not

fit_individual_rp_over_rs   # default:True, instance:bool
                            # indicates whether to fit different value for the rp_over_rs over different filters, 
                            # or not

fit_sma_over_rs             # default:False, instance:bool 
                            # indicates whether to fit for the sma_over_rs or not

fit_inclination             # default:False, instance:bool
                            # indicates whether to fit for the inclination or not

fit_mid_time                # default:True, instance:bool
                            # indicates whether to fit for the transit mid-time or not

fit_individual_times        # default:True, instance:bool
                            # indicates whether to fit different value for the transit mid-time over different epochs, 
                            # or not

fit_period                  # default:False, instance:bool
                            # indicates whether to fit for the period or not (we cannot activate both fit_period and 
                            # fit_individual_times)

fit_ldc1                    # default:False, instance:bool 
                            # indicates whether to fit for first limb-darkening coefficient or not

fit_ldc2                    # default:False, instance:bool 
                            # indicates whether to fit for second limb-darkening coefficient or not

fit_ldc3                    # default:False, instance:bool  
                            # indicates whether to fit for third limb-darkening coefficient or not

fit_ldc4                    # default:False, instance:bool 
                            # indicates whether to fit for forth limb-darkening coefficient or not

fit_rp_over_rs_limits       # default:[0.5, 2.0], instance:list (length=2)
                            # indicates the prior limits for the rp_over_rs, as a factor - by default the limits are 
                            # from half the initial rp_over_rs value to twice the initial rp_over_rs value

fit_sma_over_rs_limits      # default:[0.5, 2.0], instance:list (length=2)
                            # indicates the prior limits for the sma_over_rs, as a factor - by default the limits are 
                            # from half the initial sma_over_rs value to twice the initial sma_over_rs value

fit_inclination_limits      # default:[70.0, 90.0], instance:list (length=2)
                            # indicates the prior limits for the inclination, as a value - by default the limits are 
                            # from 70 to 90 degrees
                            
fit_mid_time_limits         # default:[-0.2, 0.2], instance:list (length=2)
                            # indicates the prior limits for the transit mid_time, as a difference - by default the 
                            # limits are from 0.2 days before the initial transit mid_time to 0.2 days after the initial 
                            # transit mid-time

fit_period_limits           # default:[0.8, 1.2], instance:list (length=2)
                            # indicates the prior limits for the period, as a factor - by default the limits are from 
                            # 0.8 times the initial period value to 1.2 times the initial period value

fit_ldc_limits              # default:[0.0, 1.0], instance:list (length=2)
                            # indicates the prior limits for the limb-darkening coefficients, as a value - by default 
                            # the limits are from 0 to 1

max_sub_exp_time            # default:10, instance:float
                            # maximum sub-exposure to be used when calculaing the exposure-integrated models

precision                   # default:3, instance:float
                            # numerical precision to be used when calculating the models
```

For eclipse observations, we need to type:

```python
planet.eclipse_fitting(output_folder)
```

and the fitting options are very similar:


```python

detrending_order            # default:2, instance:float, accepted values: 0, 1, 2
                            # indicates the order of a polynomial that will be fitted together with the observation 
                            # for de-trending purposes (every observation will be de-trended by a different polynomial

iterations                  # default:130000, instance:float
                            # indicates the number of MCMC iterations

walkers                     # default:200, instance:float
                            # indicates the number of MCMC wakers

burn_in                     # default:30000, instance:float
                            # indicates the number of MCMC burn-in

fit_fp_over_fs              # default:True, instance:bool
                            # indicates whether to fit for the fp_over_fs or not

fit_individual_fp_over_fs   # default:False, instance:bool
                            # indicates whether to fit for different values of fp_over_fs over different filters or not

fit_rp_over_rs              # default:False, instance:bool
                            # indicates whether to fit for the rp_over_rs or not

fit_individual_rp_over_rs   # default:False, instance:bool
                            # indicates whether to fit for different values of rp_over_rs over different filters or not
                            
fit_sma_over_rs             # default:False, instance:bool 
                            # indicates whether to fit for the sma_over_rs or not

fit_inclination             # default:False, instance:bool
                            # indicates whether to fit for the inclination or not

fit_mid_time                # default:False, instance:bool
                            # indicates whether to fit for the eclipse mid-time or not

fit_individual_times        # default:False, instance:bool
                            # indicates whether to fit different value for the transit mid-time over different epochs, 
                            # or not

fit_period                  # default:False, instance:bool
                            # indicates whether to fit for the period or not (we cannot activate both fit_period and 
                            # fit_individual_times)

fit_fp_over_fs_limits       # default:[0.001, 1000.0], instance:list (length=2)
                            # indicates the prior limits for the fp_over_fs, as a factor - by default the limits are 
                            # 0.1% of the initial fp_over_fs value to 1000 times the initial fp_over_fs value

fit_rp_over_rs_limits       # default:[0.5, 2.0], instance:list (length=2)
                            # indicates the prior limits for the rp_over_rs, as a factor - by default the limits are 
                            # from half the initial rp_over_rs value to twice the initial rp_over_rs value

fit_sma_over_rs_limits      # default:[0.5, 2.0], instance:list (length=2)
                            # indicates the prior limits for the sma_over_rs, as a factor - by default the limits are 
                            # from half the initial sma_over_rs value to twice the initial sma_over_rs value

fit_inclination_limits      # default:[70.0, 90.0], instance:list (length=2)
                            # indicates the prior limits for the inclination, as a value - by default the limits are 
                            # from 70 to 90 degrees
                            
fit_mid_time_limits         # default:[-0.2, 0.2], instance:list (length=2)
                            # indicates the prior limits for the eclipse mid_time, as a difference - by default the 
                            # limits are  from 0.2 days before the initial eclipse mid_time to 0.2 days after the 
                            # initial eclipse mid-time

fit_period_limits           # default:[0.8, 1.2], instance:list (length=2)
                            # indicates the prior limits for the period, as a factor - by default the limits are from 
                            # 0.8 times the initial period value to 1.2 times the initial period value

max_sub_exp_time            # default:10, instance:float
                            # maximum sub-exposure to be used when calculaing the exposure-integrated models

precision                   # default:3, instance:float
                            # numerical precision to be used when calculating the models
```



### Models of Transits / Eclipses light-curves

Getting the forward model of a transit/eclipse light-curve is as easy as it used to be.

For a specific time series (let's call it time_array with size N), defined in 
any of the acceptable time formats (JD_UTC, MJD_UTC, HJD_UTC, HJD_TDB, BJD_UTC, BJD_TDB) 
we can calculate the transit/eclipse model or the integrated transit/eclipse model for a 
given exposure time and time stamp (start, mid, end, representing the tart the middle and 
the end of the exposure, respectively) as follows:

```python
time_array = np.arange(planet.mid_time - 0.1, planet.mid_time + 0.1, 0.001)
transit_model = planet.transit(time_array, time_format='BJD_TDB', filter_name='COUSINS_R')
transit_model_integrated = planet.transit_integrated(time_array, time_format='BJD_TDB', 
                                                     exp_time=120, time_stamp='mid',
                                                     filter_name='COUSINS_R')

time_array = np.arange(planet.eclipse_mid_time - 0.1, planet.eclipse_mid_time + 0.1, 0.001)
eclipse_model = planet.eclipse(time_array, time_format='BJD_TDB', filter_name='COUSINS_R')
eclipse_model_integrated = planet.eclipse_integrated(time_array, time_format='BJD_TDB',
                                                     exp_time=120, time_stamp='mid',
                                                     filter_name='COUSINS_R')
```

For the integrated models, the calculation is done on a sub-exposure basis. The code calculates
the model values for a number of shorter exposures, which are then averaged. The length of the 
sub-exposures can be controlled through the max_sub_exp_time argument, for which the default 
value is 10 (in seconds). The lower the number, the more precise the integrated model will be, 
but the computation time will be increased. For all exoplanets, a max_sub_exp_time=10 should be
sufficient.

### Other Transit / Eclipse calculations

Based on the parameters defined when creating the plc.Planet object, the eclipse mid-time
(in the BJD_TDB format) is calculated automatically and we can access it by typing:

```python
eclipse_mid_time = planet.eclipse_mid_time
```

Also, we can calculate the transit/eclipse depth and duration for different filters:

```python
transit_duration = planet.transit_duration('COUSINS_R')
transit_depth = planet.transit_depth('COUSINS_R')
eclipse_duration = planet.eclipse_duration('COUSINS_R')
eclipse_depth = planet.eclipse_depth('COUSINS_R')
```

### Orbital calculations

For a specific time series (let's call it time_array with size N), we can calculate 
other series, related to the orbit of the planet. These are:

1. the 3D position of the planet as a function of time. The planet.planet_orbit function 
will return a 3xN array where the elements corresponds to the x,y,z coordinates of the 
planet. The coordinate system assumes that the star is at (x,y,z)=(0,0,0), the Earth is 
at (x,y,z) = (+inf,0,0) and that the planet is at its periastron during the mid-transit 
time when periastron = 90 degrees. All coordinates are in units of stellar radii.

2. the projected distance between the planet and the star as a function of time. The 
planet.planet_star_projected_distance function will return an array of size N. All 
distances are in units of stellar radii.

3. the orbital phase of the planet. The planet.planet_phase function will return an array 
of size N.

```python
time_array = np.arange(planet.mid_time - 0.1, planet.mid_time + 0.1, 0.001)

x, y, z = planet.planet_orbit(time_array, 'BJD_TDB')
projected_distance = planet.planet_star_projected_distance(time_array, 'BJD_TDB')
planet_phase = planet.planet_phase(time_array, 'BJD_TDB')
```

## Data analysis toolkit

### Angles and times

In this version we can easily convert angles to different formats. For example:

```python
ra = plc.Hours('22:03:10.7729')

ra_in_degrees = ra.deg()
ra_in_hours = ra.hours()
ra_in_rad = ra.rad()
ra_in_dms = ra.dms()
ra_in_dms_coordinate = ra.dms_coord() # this will give the angle between -90 and 90 degrees
ra_in_hms = ra.hms()
ra_in_degrees_coordinate = ra.deg_coord() # this will give the angle between -90 and 90 degrees
```

To convert between different time systems we need forst to define aa plc.FixtedTarget 
object:

```python
ra = plc.Hours('22:03:10.7729')
dec = plc.Degrees('+18:53:03.548')
target = plc.FixedTarget(ra, dec)

mid_time_in_hjd_utc = 2458485.00380255
mid_time_in_bjd_tdb = target.convert_to_bjd_tdb(mid_time_in_hjd_utc, 'HJD_UTC')
# available formats: JD_UTC, MJD_UTC, HJD_UTC, HJD_TDB, BJD_UTC, BJD_TDB
```


## Updates
4.0.1
- Fixed packaging and test issues.

4.0.2 
- Database updates checked in memory.
- Fixed database updating loop.

4.0.3
- Fix np.float bug.

4.0.4
- Fixed packaging and test issues.
- Fixed latex strings in plots.
- Fixed matplotlib.cm.get_cmap bug.

## Licence

MIT License

Copyright (c) 2016-present Angelos Tsiaras, and collaborators

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

Collaborators:

Mario Morvan

Konstantinos Karpouzas 

Ryan Varley
