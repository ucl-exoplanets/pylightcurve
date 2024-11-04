

[![codecov](https://codecov.io/gh/ucl-exoplanets/pylightcurve/branch/master/graph/badge.svg?)](https://codecov.io/gh/ucl-exoplanets/pylightcurve)

[![Downloads](https://pepy.tech/badge/pylightcurve)](https://pepy.tech/project/pylightcurve)

# PyLightcurve 

A python package for analysing exoplanet light-curves.

<img src="https://github.com/ucl-exoplanets/pylightcurve/blob/master/logo.png" width="25%">

In PyLightcurve you will find tools for:

* Exoplanet transit and eclipse modeling.
* Flexible fitting of multi-epoch and multi-colour exoplanet transit and eclipse
      light-curves.
* Calculation of exoplanet transit and eclipse properties.
* Calculation of exoplanetary orbits.
* Transformations between different angle and timing systems.

Developed by [Angelos Tsiaras](https://www.angelostsiaras.com)

PyLightcurve makes use of the following packages:

* [Emcee](https://github.com/dfm/emcee), [Foreman-Mackey et al. 2013](http://iopscience.iop.org/article/10.1086/670067)
* [ExoTETHyS](https://github.com/ucl-exoplanets/ExoTETHyS), [Morello et al. 2020](https://iopscience.iop.org/article/10.3847/1538-3881/ab63dc)
* [Matplotlib](https://matplotlib.org), [Hunter 2007](https://ieeexplore.ieee.org/document/4160265)
* [Numpy](https://numpy.org), [Oliphant 2006](https://archive.org/details/NumPyBook)
* [SciPy](https://www.scipy.org), [Virtanen et al. 2020](https://www.nature.com/articles/s41592-019-0686-2)
* [Astropy](https://www.astropy.org), [Astropy Collaboration 2013](https://www.aanda.org/articles/aa/abs/2013/10/aa22068-13/aa22068-13.html)

and of the following catalogues:

* Exoplanet Characterisation Catalogue, developed as part of the [ExoClock Project](https://www.exoclock.space), [Kokori et al. 2020](https://ui.adsabs.harvard.edu/abs/2020arXiv201207478K/abstract)

If you are using PyLightcurve for your research please cite all the above references individually and also:

[Tsiaras et al. 2016](https://ui.adsabs.harvard.edu/abs/2016ApJ...832..202T/abstract)


# Installation

Install the dev brunch of PyLightcurve through pip:

by downloading this repo, cd in it and using the setup.py file:
 
```bash
python setup.py install
```


# Usage

Please check the "notebooks" directory for usage examples and explanations.

- Check notebook "1_the_planet_class" for easy access to planet parameters, 
properties and calculations.

- Check the notebooks under the "notebooks/2_detrending_examples" directory 
  to see how to model your light-curves.

- Check the notebook "2_core_calculations" for higher efficiency 
  (suggested for developers).



## History

### v4.1

#### Changes in usage:
- New way of accessing the LDCs and the Fp/Fs in the planet class
```python
# v4.0
limb_darkening_coefficients = planet.filter('COUSINS_R').limb_darkening_coefficients
fp_over_fs = planet.filter('COUSINS_R').fp_over_fs
rp_over_rs = planet.filter('COUSINS_R').rp_over_rs

# v4.1
limb_darkening_coefficients = planet.exotethys('COUSINS_R', method='claret', 
                                               wlrange=None, stellar_model='Phoenix_2018')
# available methods: claret, power2, quad, linear
# available stellar models: Atlas_2000, Phoenix_2012_13, Stagger_2015, Stagger_2018
#                           Phoenix_2018, Phoenix_drift_2012
fp_over_fs = planet.fp_over_fs('COUSINS_R', wlrange=None)
rp_over_rs = planet.rp_over_rs
```

- Planet eclipse time is not automatically calulated
```python
# v4.0
eclipse_mid_time = planet.eclipse_mid_time

# v4.1
eclipse_mid_time = planet.eclipse_mid_time()
```



- New way of adding custom filters or limb darkening coefficients in the planet class
```python
# v4.0
planet.add_filter('my_filter', rp_over_rs, ldc1, ldc2, ldc3, ldc4, fp_over_fs)
# The passband file should be a txt file containing two columns separated by space or tab:
# column 1: wavelength in A
# column 2: total throughput in electrons/photons

# v4.1 - add a filter to the PyLightcurve database and use it as the default filters
plc.add_filter('my_filter', 'path_to_passband.txt')
# v4.1 - add custom limb darkening coefficients to planet 
planet.add_custom_limb_darkening_coefficients([ldc1, ldc2, ldc3, ldc4], filter_name, wlrange=None, stellar_model='my_model')
```

- Stellar model is no longer defined when initialising a Planet object (no ```ldc_stellar_model``` argument), it is defined when adding an observation
  
- New way of defining the iterations in MCMC
```python
# v4.0 -  total model evaluations = iterations
iterations, walkers = 15000, 3

# v4.1 - total model evaluations = iterations x walkers
iterations, walkers = 5000, 3
```
- Time conversions are now available through the Planet class
```python
# v4.0
time_in_bjd_utc = planet.target.convert_to_bjd_tdb(time_in_hjd_utc, 'HJD_UTC')


# v4.1 
time_in_bjd_utc = planet.convert_to_bjd_tdb(time_in_hjd_utc, 'HJD_UTC')
# or
time_in_bjd_utc = plc.convert_to_bjd_tdb(ra_in_degrees, dec_in_degrees, time_in_hjd_utc, 'HJD_UTC')
# where ra, dec in degrees
# available formats: JD_UTC, MJD_UTC, HJD_UTC, HJD_TDB, BJD_UTC, BJD_TDB
```
- PyLightcurve now accepts angles as floats only (degrees). No angle transformation through 
  ```pylightcurve``` are available, please use the new ```exoclock``` package.
```python
# v4.0
ra = plc.Hours('22:03:10.7729')

ra_in_degrees = ra.deg()
ra_in_hours = ra.hours()
ra_in_rad = ra.rad()
ra_in_dms = ra.dms()
ra_in_dms_coordinate = ra.dms_coord() # this will give the angle between -90 and 90 degrees
ra_in_hms = ra.hms()
ra_in_degrees_coordinate = ra.deg_coord() # this will give the angle between -90 and 90 degrees



# v4.1 
import exoclock
ra = exoclock.Hours('22:03:10.7729')

ra_in_degrees = ra.deg()
ra_in_hours = ra.hours()
ra_in_rad = ra.rad()
ra_in_dms = ra.dms()
ra_in_dms_coordinate = ra.dms_coord() # this will give the angle between -90 and 90 degrees
ra_in_hms = ra.hms()
ra_in_degrees_coordinate = ra.deg_coord() # this will give the angle between -90 and 90 degrees
```


#### New features:
- Custom detrending
- Automatic scaling of uncertainties, outliers rejection and initial parameter optimisation
- Wrapper for ExoTETHyS and sub-bandpass LDCs
- More precise integration
- Airmass detrending is now available
- Flux/Mag conversions are now available through the Planet class

### v4.0 - [latest (v4.0.3)](https://github.com/ucl-exoplanets/pylightcurve/releases/tag/v4.0.3)

- PyLightcurve 4.0 no longer supports the use of the Open Exoplanet Catalogue (OEC), due to the large number
of mistakes in the catalogue and the absence of parameters updates. OEC has been replaced by the
Exoplanet Characterisation Catalogue, dedicated to transiting exoplanets. This catalogue does not contain all the
exoplanets at the moment but will be continuously updated with more planets.

- Version 4.0 brings a new easy way of simultaneously fitting light curves from different
sources (different flux and time formats, different exposure times), different filters
(multi-colour) and different epochs (multi-epoch).
  
- We now have the option of modeling eclipse data (finally!) and also tools to calculate the
eclipse mid-time, eclipse depth and eclipse duration, based of the orbital characteristics
of the planet.
  
- In this version, the default limb-darkening coefficients are calculated using the new
package [ExoTETHyS](https://github.com/ucl-exoplanets/ExoTETHyS). This allows the
calculation of the limb-darkening coefficients even for stars cooler than 3500 K.

- PyLighrcurve 4.0 includes flexible objects for converting different expressions of
angles (e.g. '+47:12:34.05' to degrees) and of timing systems (e.g. HJD_UTC to BJD_TDB).

### 4.0.1
- Fixed packaging and test issues.

### 4.0.2 
- Database updates checked in memory.
- Fixed database updating loop.

### 4.0.3
- Fix np.float bug.



## Licence

MIT License

Copyright (c) 2016-2023 Angelos Tsiaras, and collaborators

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

Mario Morvan, Arianna Saba, Konstantinos Karpouzas, Ryan Varley