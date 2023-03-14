![Python package](https://github.com/ucl-exoplanets/pylightcurve/workflows/Python%20package/badge.svg?branch=master)
![](https://travis-ci.com/ucl-exoplanets/pylightcurve.svg?&branch=master) 
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


# Usage

Please check the "notebooks" directory for usage examples and explanations.

- Check notebook "1_the_planet_class" for easy access to planet parameters, 
properties and calculations.

- Check the notebooks under the "notebooks/2_detrending_examples" directory 
  to see how to model your lightcurves.

- Check the notebook "2_core_calculations" for higher efficiency 
  (suggestedfor developers).



## History

### v4.1

#### Changes in usage:
- New way of acessing the LDCs and the Fp/Fs in the planet class
- New way of adding custom filters in the planet class
- Stellar model is no longer defined when initialising a Plannet object, it is defined when adding an observation
- New way of defining the iterations in MCMC

#### New features:
- Custom detrending
- Automatic scaling of uncertainties, outliers rejection and initial parameter optimisation
- Wrapper for ExoTehys and sub-bandpass LDCs
- More precise integration
- Airmass detrendig is now available
- Time conversions are now available through the Planet class 
- Flux/Mag conversions are now available through the Planet class

### v4.0 - [latest (v4.0.2)](https://github.com/ucl-exoplanets/pylightcurve/releases/tag/v4.0.2)

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