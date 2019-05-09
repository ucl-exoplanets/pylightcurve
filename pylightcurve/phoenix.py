from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from .files import *


def get_flux_integral(fits, lambda1, lambda2):

    fits = get_fits_arrays(os.path.join(databases.phoenix(), fits), 1, ['Wavelength', 'Flux'])
    wavelength_array = fits['Wavelength']
    flux_array = fits['Flux']

    binsedge = 0.5 * (wavelength_array[:-1] + wavelength_array[1:])
    binsedge1 = np.append(wavelength_array[0] - (binsedge[0] - wavelength_array[0]), binsedge)
    binsedge2 = np.append(binsedge, wavelength_array[-1] + (wavelength_array[-1] - binsedge[-1]))
    binswidth = binsedge2 - binsedge1
    arg1 = np.where((lambda1 >= binsedge1) * (lambda1 < binsedge2))[0][0]
    arg2 = np.where((lambda2 > binsedge1) * (lambda2 <= binsedge2))[0][0]

    flux = np.sum(binswidth[arg1 + 1: arg2] * wavelength_array[arg1 + 1: arg2])
    flux += flux_array[arg1] * (binsedge2[arg1]-lambda1)
    flux += flux_array[arg2] * (lambda2 - binsedge1[arg2])

    return flux


def get_spectrum(fits):

    fits = get_fits_arrays(os.path.join(databases.phoenix(), fits), 1, ['Wavelength', 'Flux'])
    wavelength_array = fits['Wavelength']
    flux_array = fits['Flux']

    return wavelength_array, flux_array


def get_flux(stellar_logg, stellar_temperature, lambda1, lambda2):

    # lambda in um

    grid = open_dict(os.path.join(databases.phoenix(), 'grid.pickle'))

    logg = list(grid.keys())[np.argmin((np.array(list(grid.keys())) - stellar_logg) ** 2)]

    temperatures = list(grid[logg])
    temperatures.sort()

    temperature = temperatures[np.argmin((np.array(temperatures) - stellar_temperature) ** 2)]

    if temperature == stellar_temperature and logg == stellar_logg:

        return get_flux_integral(grid[logg][temperature], lambda1, lambda2)

    else:

        loggs = list(grid.keys())
        loggs.sort()

        logg = loggs[np.argmin((np.array(loggs) - stellar_logg) ** 2)]
        logg1 = loggs[loggs.index(logg) - 1]
        logg2 = loggs[loggs.index(logg) + 1]

        # logg

        temperatures = list(grid[logg])
        temperatures.sort()

        temperature = temperatures[np.argmin((np.array(temperatures) - stellar_temperature) ** 2)]
        temperature1 = temperatures[temperatures.index(temperature) - 1]
        temperature2 = temperatures[temperatures.index(temperature) + 1]

        flux = get_flux_integral(grid[logg][temperature], lambda1, lambda2)
        flux1 = get_flux_integral(grid[logg][temperature1], lambda1, lambda2)
        flux2 = get_flux_integral(grid[logg][temperature2], lambda1, lambda2)

        interp = interp1d([temperature1, temperature, temperature2], [flux1, flux, flux2])

        flux_logg = float(interp(stellar_temperature))

        # logg1

        temperatures = list(grid[logg1])
        temperatures.sort()

        temperature = temperatures[np.argmin((np.array(temperatures) - stellar_temperature) ** 2)]
        temperature1 = temperatures[temperatures.index(temperature) - 1]
        temperature2 = temperatures[temperatures.index(temperature) + 1]

        flux = get_flux_integral(grid[logg1][temperature], lambda1, lambda2)
        flux1 = get_flux_integral(grid[logg1][temperature1], lambda1, lambda2)
        flux2 = get_flux_integral(grid[logg1][temperature2], lambda1, lambda2)

        interp = interp1d([temperature1, temperature, temperature2], [flux1, flux, flux2])

        flux_logg1 = float(interp(stellar_temperature))

        # logg2

        temperatures = list(grid[logg2])
        temperatures.sort()

        temperature = temperatures[np.argmin((np.array(temperatures) - stellar_temperature) ** 2)]
        temperature1 = temperatures[temperatures.index(temperature) - 1]
        temperature2 = temperatures[temperatures.index(temperature) + 1]

        flux = get_flux_integral(grid[logg2][temperature], lambda1, lambda2)
        flux1 = get_flux_integral(grid[logg2][temperature1], lambda1, lambda2)
        flux2 = get_flux_integral(grid[logg2][temperature2], lambda1, lambda2)

        interp = interp1d([temperature1, temperature, temperature2], [flux1, flux, flux2])

        flux_logg2 = float(interp(stellar_temperature))

        # final

        interp = interp1d([logg1, logg, logg2], [flux_logg1, flux_logg, flux_logg2])

        return float(interp(stellar_logg))


def find_phoenix_spectrum(stellar_logg, stellar_temperature, stellar_radius,
                          stellar_umag=None, stellar_bmag=None, stellar_vmag=None,
                          stellar_rmag=None, stellar_imag=None,
                          stellar_jmag=None, stellar_hmag=None, stellar_kmag=None,
                          stellar_lmag=None):

    grid = open_dict(os.path.join(databases.phoenix(), 'grid.pickle'))

    logg = list(grid.keys())[np.argmin((np.array(list(grid.keys())) - stellar_logg) ** 2)]

    temperatures = list(grid[logg])
    temperatures.sort()

    temperature = temperatures[np.argmin((np.array(temperatures) - stellar_temperature) ** 2)]

    if temperature == stellar_temperature and logg == stellar_logg:

        wavelength1, flux = get_spectrum(grid[logg][temperature])

    else:

        loggs = list(grid.keys())
        loggs.sort()

        logg = loggs[np.argmin((np.array(loggs) - stellar_logg) ** 2)]
        logg1 = loggs[loggs.index(logg) - 1]
        logg2 = loggs[loggs.index(logg) + 1]

        # logg1

        temperatures = list(grid[logg1])
        temperatures.sort()

        temperature = temperatures[np.argmin((np.array(temperatures) - stellar_temperature) ** 2)]
        temperature1 = temperatures[temperatures.index(temperature) - 1]
        temperature2 = temperatures[temperatures.index(temperature) + 1]

        wavelength1, flux1 = get_spectrum(grid[logg][temperature1])
        wavelength2, flux2 = get_spectrum(grid[logg][temperature2])

        flux_logg1 = flux1 + (flux2 - flux1) * (stellar_temperature - temperature1) / (temperature2 - temperature1)

        # logg2

        temperature = temperatures[np.argmin((np.array(temperatures) - stellar_temperature) ** 2)]
        temperature1 = temperatures[temperatures.index(temperature) - 1]
        temperature2 = temperatures[temperatures.index(temperature) + 1]

        wavelength1, flux1 = get_spectrum(grid[logg][temperature1])
        wavelength2, flux2 = get_spectrum(grid[logg][temperature2])

        flux_logg2 = flux1 + (flux2 - flux1) * (stellar_temperature - temperature1) / (temperature2 - temperature1)

        # final

        flux = flux_logg1 + (flux_logg2 - flux_logg1) * (stellar_logg - logg1) / (logg2 - logg1)

    wavelength = wavelength1 * 10000
    flux = flux / 10

    zmag = False
    stellar_mag = False
    phot_filter = False

    if isinstance(stellar_umag, float) or isinstance(stellar_umag, int):
        zmag = 4.175e-9
        phot_filter = 'Bessel_U.txt'
        stellar_mag = stellar_umag

    if isinstance(stellar_bmag, float) or isinstance(stellar_bmag, int):
        zmag = 6.32e-9
        phot_filter = 'Bessel_B.txt'
        stellar_mag = stellar_bmag

    if isinstance(stellar_vmag, float) or isinstance(stellar_vmag, int):
        zmag = 3.631e-9
        phot_filter = 'Bessel_V.txt'
        stellar_mag = stellar_vmag

    if isinstance(stellar_rmag, float) or isinstance(stellar_rmag, int):
        zmag = 2.177e-9
        phot_filter = 'Bessel_R.txt'
        stellar_mag = stellar_rmag

    if isinstance(stellar_imag, float) or isinstance(stellar_imag, int):
        zmag = 1.126e-9
        phot_filter = 'Bessel_I.txt'
        stellar_mag = stellar_imag

    elif isinstance(stellar_jmag, float) or isinstance(stellar_jmag, int):
        zmag = 3.147e-10
        phot_filter = 'Bessel_J.txt'
        stellar_mag = stellar_jmag

    elif isinstance(stellar_hmag, float) or isinstance(stellar_hmag, int):
        zmag = 1.138e-10
        phot_filter = 'Bessel_H.txt'
        stellar_mag = stellar_hmag

    elif isinstance(stellar_kmag, float) or isinstance(stellar_kmag, int):
        zmag = 9.961e-11
        phot_filter = 'Bessel_K.txt'
        stellar_mag = stellar_kmag

    elif isinstance(stellar_lmag, float) or isinstance(stellar_lmag, int):

        zmag = 7.08e-12
        phot_filter = 'Bessel_L.txt'
        stellar_mag = stellar_lmag

    if zmag:

        phot_filter = np.loadtxt(os.path.join(databases.phoenix(), phot_filter), unpack=True)
        band = interp1d(phot_filter[0]*10, np.maximum(0, phot_filter[1]), kind='cubic')
        test = np.where((wavelength > min(phot_filter[0]*10)) * (wavelength < max(phot_filter[0]*10)))
        wavelength_test = wavelength[test]
        flux_test = flux[test]

        f_flux = np.sum(((flux_test * band(wavelength_test) * wavelength_test)[:-1] *
                         (wavelength_test[1:] - wavelength_test[:-1]))) / np.sum(
            (band(wavelength_test) * wavelength_test)[:-1] * (wavelength_test[1:] - wavelength_test[:-1]))

        stellar_radius *= 695700

        od2 = zmag / (10 ** (stellar_mag / 2.5)) / f_flux / stellar_radius / stellar_radius

        factor = stellar_radius * stellar_radius * od2

        return wavelength / 10000, flux * factor

    else:
        return wavelength / 10000, flux

    # units: wavelength: micron, flux: erg/s/cm^2/A


def get_flux_all_temp(stellar_logg, lambda1, lambda2):

    grid = open_dict(os.path.join(databases.phoenix(), 'grid.pickle'))

    logg = list(grid.keys())[np.argmin((np.array(list(grid.keys())) - stellar_logg) ** 2)]

    temperatures = list(grid[logg])
    temperatures.sort()

    all_temperatures = np.arange(temperatures[0], temperatures[-1] + 0.1, 1.0)

    tmps = []
    fluxes = []

    for temperature in temperatures:

        try:
            fluxes.append(get_flux_integral(grid[logg][temperature], lambda1, lambda2))
            tmps.append(temperature)
        except:
            pass

    interp = interp1d(temperatures, fluxes)

    return all_temperatures, np.array(interp(all_temperatures))


def get_vmag(stellar_logg, stellar_temperature, stellar_radius, distance):

    vband = np.loadtxt(glob.glob(os.path.join(databases.phoenix(), 'Bessel_V*'))[0], unpack=True)
    vband[0] = vband[0] / 1000
    vband[1] = np.maximum(0, vband[1]) / 100

    ffcnc = find_phoenix_spectrum(stellar_logg, stellar_temperature, stellar_radius)
    test = np.where((ffcnc[0] > min(vband[0])) * (ffcnc[0] < max(vband[0])))
    ffcnc = [ffcnc[0][test], ffcnc[1][test]]

    vband = interp1d(vband[0], vband[1], kind='cubic')

    v_flux = np.sum((ffcnc[1] * vband(ffcnc[0]) * ffcnc[0])[:-1] * (ffcnc[0][1:] - ffcnc[0][:-1])) / np.sum(
        (vband(ffcnc[0]) * ffcnc[0])[:-1] * (ffcnc[0][1:] - ffcnc[0][:-1]))

    stellar_radius *= 695700
    distance *= 3.086e+13

    return 2.5 * np.log10(3.6337641128286194e-09/(v_flux * (stellar_radius * stellar_radius / distance / distance)))


def find_phoenix_spectrum_plc2(stellar_logg, stellar_temperature, stellar_radius, stellar_vmag=None):

    if not stellar_radius:
        pass

    all_files = glob.glob(os.path.join(databases.phoenix(), 'lte*'))

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
            os.path.join(databases.phoenix(), 'lte{0}-*'.format(str(int(temperature1)).zfill(5))))

        loggs = np.unique([float(os.path.split(ff)[-1][9:12]) for ff in temperature1_files])

        logg = loggs[np.argmin((loggs - stellar_logg) ** 2)]

        final_file1 = glob.glob(os.path.join(databases.phoenix(), 'lte{0}-{1}*'.format(str(int(temperature1)).zfill(5),
                                                                                       float(logg))))[0]

        temperature2_files = glob.glob(
            os.path.join(databases.phoenix(), 'lte{0}-*'.format(str(int(temperature2)).zfill(5))))

        loggs = np.unique([float(os.path.split(ff)[-1][9:12]) for ff in temperature2_files])

        logg = loggs[np.argmin((loggs - stellar_logg) ** 2)]

        final_file2 = glob.glob(os.path.join(databases.phoenix(), 'lte{0}-{1}*'.format(str(int(temperature2)).zfill(5),
                                                                                       float(logg))))[0]

        flux1 = pf.open(final_file1)[0].data / (10 ** 8)
        flux2 = pf.open(final_file2)[0].data / (10 ** 8)
        flux = flux1 + (flux2 - flux1) * (stellar_temperature - temperature1) / 100
        wavelength = pf.open(glob.glob(os.path.join(databases.phoenix(), 'WAVE*'))[0])[0].data

    else:

        temperature_files = glob.glob(os.path.join(databases.phoenix(),
                                                   'lte{0}-*'.format(str(int(temperature)).zfill(5))))

        loggs = np.unique([float(os.path.split(ff)[-1][9:12]) for ff in temperature_files])

        logg = loggs[np.argmin((loggs - stellar_logg) ** 2)]

        final_file = glob.glob(os.path.join(databases.phoenix(), 'lte{0}-{1}*'.format(str(int(temperature)).zfill(5),
                                                                                      float(logg))))[0]

        flux = pf.open(final_file)[0].data / (10 ** 8)
        wavelength = pf.open(glob.glob(os.path.join(databases.phoenix(), 'WAVE*'))[0])[0].data

    if isinstance(stellar_vmag, float) or isinstance(stellar_vmag, int):

        vfilter = np.loadtxt(glob.glob(os.path.join(databases.phoenix(), 'Bessel_V*'))[0], unpack=True)
        vband = interp1d(vfilter[0]*10, np.maximum(0, vfilter[1]) / 100.0, kind='cubic')

        test = np.where((wavelength > min(vfilter[0]*10)) * (wavelength < max(vfilter[0]*10)))
        wavelength_test = wavelength[test]
        flux_test = flux[test]

        v_flux = np.sum(((flux_test * vband(wavelength_test) * wavelength_test)[:-1] *
                         (wavelength_test[1:] - wavelength_test[:-1]))) / np.sum(
            (vband(wavelength_test) * wavelength_test)[:-1] * (wavelength_test[1:] - wavelength_test[:-1]))

        stellar_radius *= 695700

        od2 = 3.6337641128286194e-09 / (10 ** (stellar_vmag / 2.5)) / v_flux / stellar_radius / stellar_radius

        factor = stellar_radius * stellar_radius * od2

        return wavelength / 10000, flux * factor

    else:
        return wavelength / 10000, flux

    # units: wavelength: micron, flux: erg/s/cm^2/A
