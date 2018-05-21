from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

__all__ = ['find_phoenix_spectrum']

import numpy as np
from scipy.interpolate import interp1d
from astropy.io import fits as pf

from .database_handling import *


def find_phoenix_spectrum(stellar_logg, stellar_temperature, stellar_radius, stellar_vmag=None):

    all_files = glob.glob(os.path.join(phoenix_database(), 'lte*'))

    temperatures = np.unique([int(os.path.split(ff)[-1][4:8]) for ff in all_files])

    stellar_temperature = max(stellar_temperature, 2300)
    stellar_temperature = min(stellar_temperature, 12000)

    temperature = temperatures[np.argmin((temperatures - stellar_temperature) ** 2)]

    if int(stellar_temperature) not in temperatures:
        if temperature < stellar_temperature:
            temperature1 = temperature
            temperature2 = temperature + 100
        else:
            temperature1 = temperature - 100
            temperature2 = temperature

        temperature1_files = glob.glob(
            os.path.join(phoenix_database(), 'lte{0}-*'.format(str(int(temperature1)).zfill(5))))

        loggs = np.unique([float(os.path.split(ff)[-1][9:12]) for ff in temperature1_files])

        logg = loggs[np.argmin((loggs - stellar_logg) ** 2)]

        final_file1 = glob.glob(os.path.join(phoenix_database(), 'lte{0}-{1}*'.format(str(int(temperature1)).zfill(5),
                                                                                      float(logg))))[0]

        temperature2_files = glob.glob(
            os.path.join(phoenix_database(), 'lte{0}-*'.format(str(int(temperature2)).zfill(5))))

        loggs = np.unique([float(os.path.split(ff)[-1][9:12]) for ff in temperature2_files])

        logg = loggs[np.argmin((loggs - stellar_logg) ** 2)]

        final_file2 = glob.glob(os.path.join(phoenix_database(), 'lte{0}-{1}*'.format(str(int(temperature2)).zfill(5),
                                                                                      float(logg))))[0]

        flux1 = pf.open(final_file1)[0].data / (10 ** 8)
        flux2 = pf.open(final_file2)[0].data / (10 ** 8)
        flux = flux1 + (flux2 - flux1) * (stellar_temperature - temperature1) / 100
        wavelength = pf.open(glob.glob(os.path.join(phoenix_database(), 'WAVE*'))[0])[0].data

    else:

        temperature_files = glob.glob(os.path.join(phoenix_database(), 'lte{0}-*'.format(str(int(temperature)).zfill(5))))

        loggs = np.unique([float(os.path.split(ff)[-1][9:12]) for ff in temperature_files])

        logg = loggs[np.argmin((loggs - stellar_logg) ** 2)]

        final_file = glob.glob(os.path.join(phoenix_database(), 'lte{0}-{1}*'.format(str(int(temperature)).zfill(5),
                                                                                     float(logg))))[0]

        flux = pf.open(final_file)[0].data / (10 ** 8)
        wavelength = pf.open(glob.glob(os.path.join(phoenix_database(), 'WAVE*'))[0])[0].data

    if isinstance(stellar_vmag, float) or isinstance(stellar_vmag, int):

        vfilter = np.loadtxt(glob.glob(os.path.join(phoenix_database(), 'Bessel_V*'))[0], unpack=True)
        vfilter_flux = (flux * np.interp(wavelength, (vfilter[0] * 10)[::-1],
                                         (np.maximum(0, vfilter[1]) / 100.0)[::-1]))

        testx = [wavelength[0]]
        testy = [0]
        for i in range(int(len(vfilter_flux)/10000)):
            yy = vfilter_flux[1+i*10000:1+(i+1)*10000]
            xx = wavelength[1+i*10000:1+(i+1)*10000]
            mm = np.argmax(yy)
            testx.append(xx[mm])
            testy.append(yy[mm])
        testx.append(wavelength[-1])
        testy.append(0)
        smooth = interp1d(np.array(testx), np.array(testy), kind='cubic')

        test_wl = np.sum(wavelength * vfilter_flux) / np.sum(vfilter_flux)

        factor = (3.735 * (10 ** (-9)) * (10 ** (stellar_vmag / (-2.5))) /
                  smooth(wavelength)[np.argmin(np.abs(wavelength - test_wl))])

        return wavelength / 10000, flux * factor

    else:
        return wavelength / 10000, flux

    # units: wavelength: micron, flux: erg/s/cm^2/A

