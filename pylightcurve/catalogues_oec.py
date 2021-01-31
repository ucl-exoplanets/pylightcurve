from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from ._1databases import *


def oec_catalogue():

        oec_file = glob.glob(os.path.join(databases.oec, '*'))[0]

        return exodata.OECDatabase(gzip.GzipFile(oec_file), stream=True)


def find_target(target, catalogue=None):

    if catalogue is None:
        catalogue = oec_catalogue()

    if target == 'Kepler-16 (AB) b':
        planet = catalogue.searchPlanet('Kepler-16 b')
    else:
        planet = catalogue.searchPlanet(target)

    if target == 'XO-1 b':
        planet = planet[1]
    elif isinstance(planet, list):
        planet = planet[0]

    return planet


def find_oec_parameters(target, catalogue=None, binary_star=0):

    planet = find_target(target, catalogue)

    name = planet.name

    try:
        star = planet.star
    except exodata.astroclasses.HierarchyError:
        star = planet.binary.stars[binary_star]

    # known mistakes in the catalogue

    if name in ['HD 3167 b', 'HD 3167 c']:
        star.R = 0.828 * aq.R_s

    stellar_radius = star.R

    stellar_logg = float(star.calcLogg())

    stellar_temperature = float(star.T)

    if np.isnan(star.Z):
        stellar_metallicity = 0
    else:
        stellar_metallicity = star.Z

    rp_over_rs = float(planet.R.rescale(aq.m) / stellar_radius.rescale(aq.m))

    planet_temperature = float(planet.calcTemperature())
    fp_over_fs = (rp_over_rs ** 2.0) * \
                 ((np.exp(6.626 * (10 ** (-5)) * 299792458.0 / (1.4 * 1.381 * stellar_temperature)) - 1.)
                  / (np.exp(6.626 * (10 ** (-5)) * 299792458.0 / (1.4 * 1.381 * planet_temperature)) - 1.))

    period = float(planet.P)

    if np.isnan(planet.a):
        try:
            sma_over_rs = float(planet.calcSMA().rescale(aq.m) / stellar_radius.rescale(aq.m))
        except:
            sma_over_rs = None
    else:
        sma_over_rs = float(planet.a.rescale(aq.m) / stellar_radius.rescale(aq.m))

    if np.isnan(planet.e):
        eccentricity = 0.0
    else:
        eccentricity = float(planet.e)

    if np.isnan(planet.i):
        inclination = 90.0
    else:
        inclination = float(planet.i)

    if np.isnan(planet.periastron):
        periastron = 0.0
    else:
        periastron = float(planet.periastron)

    if np.isnan(planet.transittime):
        mid_time = None
    else:
        mid_time = float(planet.transittime)

    # known mistakes in the catalogue

    if name == 'WASP-121 b':
        mid_time = 2456635.70832
    elif name == 'Qatar-1 b':
        mid_time = 2455518.4102
    elif name == 'HD 209458 b':
        mid_time = 2452826.628521
    elif name == 'KELT-11 b':
        mid_time = 2457061.9098
    elif name == 'WASP-107 b':
        mid_time = 2456514.4106
    elif name == 'WASP-62 b':
        period = 4.411953
        mid_time = 2455855.39195
    elif name == 'HD 3167 b':
        mid_time = 2457394.37450
    elif name == 'HD 3167 c':
        mid_time = 2457394.9788
    elif name == 'WASP-127 b':
        mid_time = 2457248.74126
    elif name == 'HIP 41378 b':
        mid_time = 2457152.2844
    elif name == 'HIP 41378 c':
        mid_time = 2457163.1659
    elif name == 'HIP 41378 d':
        mid_time = 2457166.2629
    elif name == 'HIP 41378 e':
        mid_time = 2457142.01656
    elif name == 'HIP 41378 f':
        mid_time = 2457186.91451
    elif name == 'TrES-2':
        mid_time = 2454955.762517
    elif name == 'HD 97658 b':
        mid_time = 2456665.46415
    elif name == 'Kepler-9 b':
        mid_time += 0.075 * period
    elif name == 'Kepler-9 c':
        period = 39.070
        mid_time = 2456293.581329965 + 0.15
    elif name == 'Kepler-9 d':
        mid_time -= 0.03 * period
    elif name == 'Kepler-11 e':
        period = 31.99590
        mid_time = 2454987.1590
    elif name == 'Kepler-11 d':
        period = 22.68719
        mid_time = 2454981.4550
    elif name == 'KOI-314 c':
        mid_time = 2455558.1404
    elif name == 'HD 219134 b':
        mid_time = 2457463.82884
        period = 3.092926
    elif name == 'HD 106315 b':
        mid_time = 2457605.6521000
    elif name == 'HD 106315 c':
        mid_time = 2457611.1310000
    elif name == 'Kepler-16 (AB) b':
        mid_time = 2455397.522
    elif name == 'CoRoT-1 b':
        mid_time = 2454524.62324
    elif name == 'HAT-P-7 b':
        mid_time = 2454700.81308
    elif name == 'Kepler-1625 b':
        mid_time = 2456043.9587
    elif name == 'WASP-79 b':
        mid_time = 2455545.2356

    return (name, stellar_logg, stellar_temperature, stellar_metallicity, rp_over_rs, fp_over_fs, period, sma_over_rs,
            eccentricity, inclination, periastron, mid_time)


def find_oec_coordinates(target, catalogue=None, output='deg'):

    planet = find_target(target, catalogue)

    if output == 'deg':
        return planet.system.ra.deg, planet.system.dec.deg
    elif output == 'str':
        ra = '{0}:{1}:{2}'.format(str(int(planet.system.ra.hms[0])).zfill(2),
                                  str(int(planet.system.ra.hms[1])).zfill(2),
                                  str(round(planet.system.ra.hms[2],1)).zfill(4))
        if planet.system.dec.deg > 0:
            dec = '+{0}:{1}:{2}'.format(str(int(planet.system.dec.dms[0])).zfill(2),
                                        str(int(planet.system.dec.dms[1])).zfill(2),
                                        str(round(planet.system.dec.dms[2], 1)).zfill(4))
        else:
            dec = '-{0}:{1}:{2}'.format(str(int(abs(planet.system.dec.dms[0]))).zfill(2),
                                        str(int(abs(planet.system.dec.dms[1]))).zfill(2),
                                        str(round(abs(planet.system.dec.dms[2]), 1)).zfill(4))

        return '{0} {1}'.format(ra, dec)


def find_oec_stellar_parameters(target, catalogue=None, binary_star=0):

    planet = find_target(target, catalogue)

    name = planet.name

    try:
        star = planet.star
    except exodata.astroclasses.HierarchyError:
        star = planet.binary.stars[binary_star]

    # known mistakes in the catalogue

    if name in ['HD 3167 b', 'HD 3167 c']:
        star.R = 0.828 * aq.R_s
    try:
        star.R = float(star.R) * aq.R_s
        star.M = float(star.M) * aq.M_s
        stellar_logg = float(star.calcLogg())
    except:
        stellar_logg = None

    stellar_temperature = float(star.T)

    if np.isnan(star.Z):
        stellar_metallicity = 0
    else:
        stellar_metallicity = star.Z

    stellar_radius = float(planet.star.R)

    try:
        magU = float(planet.star.magU)
    except:
        magU = None

    try:
        magB = float(planet.star.magB)
    except:
        magB = None

    try:
        magV = float(planet.star.magV)
    except:
        magV = None

    magR = None
    magI = None

    try:
        magJ = float(planet.star.magJ)
    except:
        magJ = None

    try:
        magH = float(planet.star.magH)
    except:
        magH = None

    try:
        magK = float(planet.star.magK)
    except:
        magK = None

    try:
        magL = float(planet.star.magL)
    except:
        magL = None

    stellar_maglist = [magU, magB, magV, magR, magI, magJ, magH, magK, magL]

    if np.array([mag is None for mag in stellar_maglist]).all():
        stellar_maglist = None

    return name, stellar_logg, stellar_temperature, stellar_metallicity, stellar_radius, stellar_maglist
