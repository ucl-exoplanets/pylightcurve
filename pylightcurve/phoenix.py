from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

__all__ = ['find_phoenix_spectrum']

import numpy as np
from astropy.io import fits as pf

from .database_handling import *


def find_phoenix_spectrum(stellar_logg, stellar_temperature, stellar_radius, stellar_vmag=None):

    all_files = glob.glob(os.path.join(phoenix_database(), 'lte*'))

    temperatures = np.unique([int(os.path.split(ff)[-1][4:8]) for ff in all_files])

    temperature = temperatures[np.argmin((temperatures - stellar_temperature) ** 2)]

    temperature_files = glob.glob(os.path.join(phoenix_database(), 'lte{0}-*'.format(str(int(temperature)).zfill(5))))

    loggs = np.unique([float(os.path.split(ff)[-1][9:12]) for ff in temperature_files])

    logg = loggs[np.argmin((loggs - stellar_logg) ** 2)]

    final_file = glob.glob(os.path.join(phoenix_database(), 'lte{0}-{1}*'.format(str(int(temperature)).zfill(5),
                                                                                 float(logg))))[0]

    flux = pf.open(final_file)[0].data / (10 ** 8)
    wavelength = pf.open(glob.glob(os.path.join(phoenix_database(), 'WAVE*'))[0])[0].data

    if isinstance(stellar_vmag, float) or isinstance(stellar_vmag, int):

        r_sun = 69550800000.0

        r = stellar_radius * r_sun

        vfilter = np.loadtxt(glob.glob(os.path.join(phoenix_database(), 'Bessel_V*'))[0], unpack=True)
        vfilter_flux = (flux * np.interp(wavelength, (vfilter[0] * 10)[::-1],
                                         (np.maximum(0, vfilter[1]) / 100.0)[::-1]))

        test_wl = np.arange(wavelength[0], wavelength[-1], 0.1)
        test_flux = np.interp(test_wl[:-1], wavelength, vfilter_flux) * (test_wl[1:] - test_wl[:-1]) * 4 * np.pi * r * r

        mabs = 4.82 - 2.5 * np.log10(np.sum(test_flux) / (5.09964777265 * (10 ** 32)))
        mapp = 11.08

        distance = (10 ** ((mapp - mabs + 5) / 5)) * 3.0857 * (10 ** 18)

        return wavelength, flux * (r * r / distance / distance)

    else:
        return wavelength, flux

    # units: wavelength: A, flux: erg/s/cm^2/A

