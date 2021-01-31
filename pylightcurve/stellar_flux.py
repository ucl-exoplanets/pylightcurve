from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from .tools_files import *


def get_vega_flux(mag):

    zmag = None
    phot_filter = None
    if mag == 'U':
        zmag = 4.175e-9
        phot_filter = 'Bessel_U.txt'

    if mag == 'B':
        zmag = 6.32e-9
        phot_filter = 'Bessel_B.txt'

    if mag == 'V':
        zmag = 3.631e-9
        phot_filter = 'Bessel_V.txt'

    if mag == 'R':
        zmag = 2.177e-9
        phot_filter = 'Bessel_R.txt'

    if mag == 'I':
        zmag = 1.126e-9
        phot_filter = 'Bessel_I.txt'

    if mag == 'J':
        zmag = 3.147e-10
        phot_filter = 'Bessel_J.txt'

    if mag == 'H':
        zmag = 1.138e-10
        phot_filter = 'Bessel_H.txt'

    if mag == 'K':
        zmag = 9.961e-11
        phot_filter = 'Bessel_K.txt'

    if mag == 'L':
        zmag = 7.08e-12
        phot_filter = 'Bessel_L.txt'

    if mag == 'G':
        zmag = 2.4773273088492937e-09
        phot_filter = 'G.txt'

    wavelength = pf.open('/Users/angelos/Desktop/alpha_lyr_004.fits')[1].data['WAVELENGTH']
    flux = pf.open('/Users/angelos/Desktop/alpha_lyr_004.fits')[1].data['FLUX']

    phot_filter = np.loadtxt(os.path.join(databases.phoenix(), phot_filter), unpack=True)
    phot_filter[0] = phot_filter[0] * 10
    band = interp1d(phot_filter[0], np.maximum(0, phot_filter[1]), kind='cubic')
    test = np.where((wavelength > min(phot_filter[0])) * (wavelength < max(phot_filter[0])))
    wavelength_test = wavelength[test]
    flux_test = flux[test]

    flux = np.sum(((flux_test * band(wavelength_test) * wavelength_test)[:-1] *
                     (wavelength_test[1:] - wavelength_test[:-1]))) / np.sum(
        (band(wavelength_test) * wavelength_test)[:-1] * (wavelength_test[1:] - wavelength_test[:-1]))

    return flux


def open_spectrum_file(fits):

    fits = get_fits_arrays(os.path.join(databases.phoenix(), fits), 1, ['Wavelength', 'Flux'])
    wavelength_array = fits['Wavelength']
    flux_array = fits['Flux']

    return wavelength_array * 10000, flux_array / 10
    # units: wavelength: A, flux: erg/s/cm^2/A


def get_spectrum(stellar_logg, stellar_temperature, stellar_radius,
                 stellar_umag=None, stellar_bmag=None, stellar_vmag=None,
                 stellar_rmag=None, stellar_imag=None,
                 stellar_jmag=None, stellar_hmag=None, stellar_kmag=None,
                 stellar_lmag=None, stellar_gmag=None):

    grid = open_dict(os.path.join(databases.phoenix(), 'grid.pickle'))

    logg = list(grid.keys())[np.argmin((np.array(list(grid.keys())) - stellar_logg) ** 2)]

    temperatures = list(grid[logg])
    temperatures.sort()

    temperature = temperatures[np.argmin((np.array(temperatures) - stellar_temperature) ** 2)]

    if temperature == stellar_temperature and logg == stellar_logg:

        wavelength1, flux = open_spectrum_file(grid[logg][temperature])

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

        wavelength1, flux1 = open_spectrum_file(grid[logg][temperature1])
        wavelength2, flux2 = open_spectrum_file(grid[logg][temperature2])

        flux_logg1 = flux1 + (flux2 - flux1) * (stellar_temperature - temperature1) / (temperature2 - temperature1)

        # logg2

        temperature = temperatures[np.argmin((np.array(temperatures) - stellar_temperature) ** 2)]
        temperature1 = temperatures[temperatures.index(temperature) - 1]
        temperature2 = temperatures[temperatures.index(temperature) + 1]

        wavelength1, flux1 = open_spectrum_file(grid[logg][temperature1])
        wavelength2, flux2 = open_spectrum_file(grid[logg][temperature2])

        flux_logg2 = flux1 + (flux2 - flux1) * (stellar_temperature - temperature1) / (temperature2 - temperature1)

        # final

        flux = flux_logg1 + (flux_logg2 - flux_logg1) * (stellar_logg - logg1) / (logg2 - logg1)

    wavelength = wavelength1

    zmag = False
    stellar_mag = False
    phot_filter = False

    # Bessel 1990, PASP 102, 1181

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

    if isinstance(stellar_jmag, float) or isinstance(stellar_jmag, int):
        zmag = 3.147e-10
        phot_filter = 'Bessel_J.txt'
        stellar_mag = stellar_jmag

    if isinstance(stellar_hmag, float) or isinstance(stellar_hmag, int):
        zmag = 1.138e-10
        phot_filter = 'Bessel_H.txt'
        stellar_mag = stellar_hmag

    if isinstance(stellar_kmag, float) or isinstance(stellar_kmag, int):
        zmag = 9.961e-11
        phot_filter = 'Bessel_K.txt'
        stellar_mag = stellar_kmag

    if isinstance(stellar_lmag, float) or isinstance(stellar_lmag, int):
        zmag = 7.08e-12
        phot_filter = 'Bessel_L.txt'
        stellar_mag = stellar_lmag

    if isinstance(stellar_gmag, float) or isinstance(stellar_gmag, int):
        zmag = 2.4773273088492937e-09
        phot_filter = 'G.txt'
        stellar_mag = stellar_gmag

    if zmag:

        phot_filter = np.loadtxt(os.path.join(databases.phoenix(), phot_filter), unpack=True)
        phot_filter[0] = phot_filter[0] * 10
        band = interp1d(phot_filter[0], np.maximum(0, phot_filter[1]), kind='cubic')
        test = np.where((wavelength > min(phot_filter[0])) * (wavelength < max(phot_filter[0])))
        wavelength_test = wavelength[test]
        flux_test = flux[test]

        f_flux = np.sum(((flux_test * band(wavelength_test) * wavelength_test)[:-1] *
                         (wavelength_test[1:] - wavelength_test[:-1]))) / np.sum(
            (band(wavelength_test) * wavelength_test)[:-1] * (wavelength_test[1:] - wavelength_test[:-1]))

        # stellar_radius *= 695700

        factor = zmag / (10 ** (stellar_mag / 2.5)) / f_flux

        return wavelength, flux * factor

    else:
        return wavelength, flux


def get_flux(stellar_logg, stellar_temperature, stellar_radius, lambda1, lambda2,
             stellar_umag=None, stellar_bmag=None, stellar_vmag=None,
             stellar_rmag=None, stellar_imag=None,
             stellar_jmag=None, stellar_hmag=None, stellar_kmag=None,
             stellar_lmag=None, stellar_gmag=None):

    wavelength_array, flux_array = get_spectrum(stellar_logg, stellar_temperature, stellar_radius,
                                                stellar_umag, stellar_bmag, stellar_vmag, stellar_rmag,
                                                stellar_imag, stellar_jmag, stellar_hmag, stellar_kmag,
                                                stellar_lmag, stellar_gmag)

    binsedge = 0.5 * (wavelength_array[:-1] + wavelength_array[1:])
    binsedge1 = np.append(wavelength_array[0] - (binsedge[0] - wavelength_array[0]), binsedge)
    binsedge2 = np.append(binsedge, wavelength_array[-1] + (wavelength_array[-1] - binsedge[-1]))
    binswidth = binsedge2 - binsedge1
    arg1 = np.where((lambda1 >= binsedge1) * (lambda1 < binsedge2))[0][0]
    arg2 = np.where((lambda2 > binsedge1) * (lambda2 <= binsedge2))[0][0]

    flux = np.sum(binswidth[arg1 + 1: arg2] * flux_array[arg1 + 1: arg2])
    flux += flux_array[arg1] * (binsedge2[arg1]-lambda1)
    flux += flux_array[arg2] * (lambda2 - binsedge1[arg2])

    return flux


def get_mag(stellar_logg, stellar_temperature, stellar_radius, distance, mag):

    # Bessel 1990, PASP 102, 1181
    zmag = None
    phot_filter = None
    if mag == 'U':
        zmag = 4.175e-9
        phot_filter = 'Bessel_U.txt'

    if mag == 'B':
        zmag = 6.32e-9
        phot_filter = 'Bessel_B.txt'

    if mag == 'V':
        zmag = 3.631e-9
        phot_filter = 'Bessel_V.txt'

    if mag == 'R':
        zmag = 2.177e-9
        phot_filter = 'Bessel_R.txt'

    if mag == 'I':
        zmag = 1.126e-9
        phot_filter = 'Bessel_I.txt'

    if mag == 'J':
        zmag = 3.147e-10
        phot_filter = 'Bessel_J.txt'

    if mag == 'H':
        zmag = 1.138e-10
        phot_filter = 'Bessel_H.txt'

    if mag == 'K':
        zmag = 9.961e-11
        phot_filter = 'Bessel_K.txt'

    if mag == 'L':
        zmag = 7.08e-12
        phot_filter = 'Bessel_L.txt'

    if mag == 'G':
        zmag = 2.4773273088492937e-09
        phot_filter = 'G.txt'

    if not zmag:
        raise PyLCInputError('Filter not available.')
    else:
        band = np.loadtxt(os.path.join(databases.phoenix(), phot_filter), unpack=True)
        band[0] = band[0] * 10

        ffcnc = get_spectrum(stellar_logg, stellar_temperature, stellar_radius)
        test = np.where((ffcnc[0] > min(band[0])) * (ffcnc[0] < max(band[0])))
        ffcnc = [ffcnc[0][test], ffcnc[1][test]]

        band = interp1d(band[0], band[1], kind='cubic')

        flux = np.sum((ffcnc[1] * band(ffcnc[0]) * ffcnc[0])[:-1] * (ffcnc[0][1:] - ffcnc[0][:-1])) / np.sum(
            (band(ffcnc[0]) * ffcnc[0])[:-1] * (ffcnc[0][1:] - ffcnc[0][:-1]))

        stellar_radius *= 695700
        distance *= 3.08567758149137e+13

        return 2.5 * np.log10(zmag/(flux * (stellar_radius * stellar_radius / distance / distance)))


def g_to_r_mag(stellar_logg, stellar_temperature, stellar_radius, stellar_mag):

    wavelength, flux = get_spectrum(stellar_logg, stellar_temperature, stellar_radius, stellar_gmag=stellar_mag)
    phot_filter = 'G.txt'

    phot_filter = np.loadtxt(os.path.join(databases.phoenix(), phot_filter), unpack=True)
    phot_filter[0] = phot_filter[0] * 10
    band = interp1d(phot_filter[0], np.maximum(0, phot_filter[1]), kind='cubic')
    test = np.where((wavelength > min(phot_filter[0])) * (wavelength < max(phot_filter[0])))
    wavelength_test = wavelength[test]
    flux_test = flux[test]

    flux = np.sum(((flux_test * band(wavelength_test) * wavelength_test)[:-1] *
                     (wavelength_test[1:] - wavelength_test[:-1]))) / np.sum(
        (band(wavelength_test) * wavelength_test)[:-1] * (wavelength_test[1:] - wavelength_test[:-1]))

    return 2.5 * np.log10(2.2090695086397724e-09/flux)


def g_to_v_mag(stellar_logg, stellar_temperature, stellar_radius, stellar_mag):

    wavelength, flux = get_spectrum(stellar_logg, stellar_temperature, stellar_radius, stellar_gmag=stellar_mag)
    phot_filter = 'G.txt'

    phot_filter = np.loadtxt(os.path.join(databases.phoenix(), phot_filter), unpack=True)
    phot_filter[0] = phot_filter[0] * 10
    band = interp1d(phot_filter[0], np.maximum(0, phot_filter[1]), kind='cubic')
    test = np.where((wavelength > min(phot_filter[0])) * (wavelength < max(phot_filter[0])))
    wavelength_test = wavelength[test]
    flux_test = flux[test]

    flux = np.sum(((flux_test * band(wavelength_test) * wavelength_test)[:-1] *
                     (wavelength_test[1:] - wavelength_test[:-1]))) / np.sum(
        (band(wavelength_test) * wavelength_test)[:-1] * (wavelength_test[1:] - wavelength_test[:-1]))

    return 2.5 * np.log10(3.519416551428086e-09/flux)


def clablimb(method, stellar_logg, stellar_temperature, stellar_metallicity, photometric_filter, stellar_model='ATLAS'):

    if method != 'claret':
        raise ImportError('Limb darkening model not currently supported.')

    data_file = glob.glob(os.path.join(databases.clablimb,
                                       '*_{0}_{1}.txt'.format(stellar_model, photometric_filter)))[0]
    data = np.loadtxt(data_file, usecols=[0, 1, 2, 3, 4, 5, 6, 7], unpack=True)

    xx = np.unique(data[0])
    yy = np.unique(data[1])
    zz = np.unique(data[2])

    yin = float(stellar_temperature)
    zin = float(stellar_metallicity)
    xin = float(stellar_logg)

    if xin in xx:
        xmin, xmax = xin, xin
    else:
        test = np.argmin(np.abs(np.ones_like(xx) * xin - xx))
        if xx[test] < xin:
            xmin = xx[test]
            xmax = xx[min(len(xx) - 1, test + 1)]
        else:
            xmin = xx[max(0, test - 1)]
            xmax = xx[test]
    if yin in yy:
        ymin, ymax = yin, yin
    else:
        test = np.argmin(np.abs(np.ones_like(yy) * yin - yy))
        if yy[test] < yin:
            ymin = yy[test]
            ymax = yy[min(len(yy) - 1, test + 1)]
        else:
            ymin = yy[max(0, test - 1)]
            ymax = yy[test]
    if zin in zz:
        zmin, zmax = zin, zin
    else:
        test = np.argmin(np.abs(np.ones_like(zz) * zin - zz))
        if zz[test] < zin:
            zmin = zz[test]
            zmax = zz[min(len(zz) - 1, test + 1)]
        else:
            zmin = zz[max(0, test - 1)]
            zmax = zz[test]

    def tri_linear(x, y, z, x0, x1, y0, y1, z0, z1, v000, v100, v010, v001, v101, v011, v110, v111):
        c0 = v000
        c1 = v100 - v000
        c2 = v010 - v000
        c3 = v001 - v000
        c4 = v110 - v010 - v100 + v000
        c5 = v011 - v001 - v010 + v000
        c6 = v101 - v001 - v100 + v000
        c7 = v111 - v011 - v101 - v110 + v100 + v001 + v010 - v000
        if x == x0 == x1:
            dx = 0
        else:
            dx = (x - x0) / (x1 - x0)
        if y == y0 == y1:
            dy = 0
        else:
            dy = (y - y0) / (y1 - y0)
        if z == z0 == z1:
            dz = 0
        else:
            dz = (z - z0) / (z1 - z0)
        return c0 + c1 * dx + c2 * dy + c3 * dz + c4 * dx * dy + c5 * dy * dz + c6 * dz * dx + c7 * dx * dy * dz

    final_coefficients = []

    for index in [4, 5, 6, 7]:

        vv000 = data[index][np.where((data[0] == xmin) & (data[1] == ymin) & (data[2] == zmin))][0]
        vv100 = data[index][np.where((data[0] == xmax) & (data[1] == ymin) & (data[2] == zmin))][0]
        vv010 = data[index][np.where((data[0] == xmin) & (data[1] == ymax) & (data[2] == zmin))][0]
        vv001 = data[index][np.where((data[0] == xmin) & (data[1] == ymin) & (data[2] == zmax))][0]
        vv101 = data[index][np.where((data[0] == xmax) & (data[1] == ymin) & (data[2] == zmax))][0]
        vv011 = data[index][np.where((data[0] == xmin) & (data[1] == ymax) & (data[2] == zmax))][0]
        vv110 = data[index][np.where((data[0] == xmax) & (data[1] == ymax) & (data[2] == zmin))][0]
        vv111 = data[index][np.where((data[0] == xmax) & (data[1] == ymax) & (data[2] == zmax))][0]

        res = tri_linear(xin, yin, zin, xmin, xmax, ymin, ymax, zmin, zmax,
                         vv000, vv100, vv010, vv001, vv101, vv011, vv110, vv111)

        final_coefficients.append(res)

    return np.array(final_coefficients)
