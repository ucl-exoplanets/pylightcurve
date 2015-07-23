import numpy as np
import matplotlib.pyplot as plt
import os
import glob

import ephem

import models
import mcmc
import animation

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
pi = np.pi


class PyLCError(BaseException):
    pass


class PyLCValueError(PyLCError):
    pass


class PyLCFilterError(PyLCError):
    pass


def ldcoeff(metall, teff, logg, phot_filter):
    """ Looks up the non quadtractic limb darkening coefficients in the Claret table

    :param metall: stellar metallicity (Fe/H)
    :param teff: Effective temperature of the star (Kelvin)
    :param logg: log(g) of the star
    :param phot_filter: which filter to retreive the coefficents for out of u, v, b, y, U, B, V, R, I, J, H, K

    :return: The 4 non quadratic limb darkening coefficients
    :rtype: (a1, a2 ,a3, a4)

    :raises PyLC_FilterError: If invalid filter is given
    """
    filterlist = (('u', 'v', 'b', 'y', 'U', 'B', 'V', 'R', 'I', 'J', 'H', 'K'),
                  (4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15))

    if phot_filter not in filterlist[0]:
        raise PyLCFilterError("Invalid filter, got {} must be in {}".format(phot_filter, filterlist[0]))

    # This could probably all be cleaned up by importing to a pandas dataframe
    phot_filter = filterlist[1][filterlist[0].index(phot_filter)]
    tables, mett = np.loadtxt(glob.glob(__location__ + '/*claretinfo*')[0], usecols=(0, 4), unpack=True)
    table = str(int(tables[np.argmin(abs(metall - mett))]))
    table_file = glob.glob(__location__ + '/*/TABLE' + table)[0]
    logglist, tefflist = np.loadtxt(table_file, usecols=(1, 2), unpack=True, skiprows=5)
    teff0 = tefflist[np.argmin(abs(teff - tefflist))]
    logg0 = logglist[np.argmin(abs(logg - logglist))]
    ld_coeffs = []
    for i in open(table_file).readlines()[5:]:
        coef = float(i.split()[phot_filter])
        logg = float(i.split()[1])
        teff = float(i.split()[2])
        if logg == logg0 and teff == teff0:
            ld_coeffs.append(coef)
    return tuple(ld_coeffs)


class Planet:
    def __init__(self):
        self.rp_rs = np.nan
        self.period = np.nan
        self.a_rs = np.nan
        self.eccentricity = np.nan
        self.inclination = np.nan
        self.omega = np.nan
        self.Omega = np.nan
        self.mid_transit = np.nan
    
    def __repr__(self):
        """Present Planet parameters."""
        return ('Planet parameters: \n'
                '   rp_rs              = %g \n'
                '   period      [days] = %g \n'
                '   a_rs               = %g \n'
                '   eccenticity        = %g \n'
                '   inclination [deg]  = %g \n'
                '   omega       [deg]  = %g \n'
                '   Omega       [deg]  = %g \n'
                '   mid_transit [jd]   = %g \n'
                % (self.rp_rs, self.period,
                    self.a_rs, self.eccentricity, self.inclination,
                    self.omega, self.Omega, self.mid_transit))

    def set_example(self):
        self.rp_rs = 0.15
        self.period = 2.2
        self.a_rs = 9.0
        self.eccentricity = 0.0
        self.inclination = 87.5
        self.omega = 0.0
        self.Omega = 0.0
        self.mid_transit = 100.0

    def test(self):
        if np.isnan(self.rp_rs):
            raise PyLCValueError("rp_rs is not set")
        elif np.isnan(self.period):
            raise PyLCValueError("period is not set")
        elif np.isnan(self.a_rs):
            raise PyLCValueError("a_rs is not set")
        elif np.isnan(self.eccentricity):
            raise PyLCValueError("eccentricity is not set")
        elif np.isnan(self.inclination):
            raise PyLCValueError("inclination is not set")
        elif np.isnan(self.omega):
            raise PyLCValueError("omega is not set")
        elif np.isnan(self.Omega):
            raise PyLCValueError("Omega is not set")
        elif np.isnan(self.mid_transit):
            raise PyLCValueError("mid_transit is not set")

    def next_transit(self, timezone, number):
        self.test()
        ww = self.omega * np.pi / 180
        ii = self.inclination * np.pi / 180
        ee = self.eccentricity
        aa = self.a_rs
        ro_pt = (1 - ee ** 2) / (1 + ee * np.sin(ww))
        b_pt = aa * ro_pt * np.cos(ii)
        s_ps = 1.0 + self.rp_rs
        df = np.arcsin(np.sqrt((s_ps ** 2 - b_pt ** 2) / ((aa ** 2) * (ro_pt ** 2) - b_pt ** 2)))
        duration = (self.period * (ro_pt ** 2)) / (np.pi * np.sqrt(1 - ee ** 2)) * df
        now = float(ephem.now())
        t0 = self.mid_transit - 2415020.0
        next_transit = t0 + (int((now - t0) / self.period) + 1) * self.period
        for i in range(number):
            transit_start = ephem.date(next_transit + i * self.period - duration + timezone / 24.0)
            transit_midle = ephem.date(next_transit + i * self.period + timezone / 24.0)
            transit_close = ephem.date(next_transit + i * self.period + duration + timezone / 24.0)
            print transit_start, transit_midle, transit_close


class Star:
    def __init__(self, coefficients=None, properties=None):
        coeficients_set = False
        if isinstance(coefficients, tuple):
            if len(coefficients) == 4:
                self.metallicity, self.temperature, self.logg, self.filter = (np.nan, np.nan, np.nan, 'np.nan')
                self.ld_1 = coefficients[0]
                self.ld_2 = coefficients[1]
                self.ld_3 = coefficients[2]
                self.ld_4 = coefficients[3]
                coeficients_set = True
        if isinstance(properties, tuple):
            if len(properties) == 4:
                self.metallicity, self.temperature, self.logg, self.filter = properties
                coefficients = ldcoeff(self.metallicity, self.temperature, self.logg, self.filter)
                self.ld_1 = coefficients[0]
                self.ld_2 = coefficients[1]
                self.ld_3 = coefficients[2]
                self.ld_4 = coefficients[3]
                coeficients_set = True
        if not coeficients_set:
            raise PyLCValueError("\n The coefficients parameter should be a tuple with four elements "
                                 "(the four limb darkening coefisients) or"
                                 "\n the properties parameter should be a tuple with four elements "
                                 "(metallicity, teff, logg, photometric filter)")

    def __repr__(self):
        """Present Star parameters."""
        return ('Stellar parameters: \n'
                '   temperature [K]    = %g \n'
                '   metallicity [Fe/H] = %g \n'
                '   logg        [cgs]  = %g \n'
                '   filter             = %s \n'
                '   ld_1               = %g \n'
                '   ld_2               = %g \n'
                '   ld_3               = %g \n'
                '   ld_4               = %g \n'
                % (self.temperature, self.metallicity, self.logg, self.filter,
                   self.ld_1, self.ld_2, self.ld_3, self.ld_4))

    def set_example(self):
        self.temperature = 6590
        self.metallicity = 0.01
        self.logg = 4.1
        coefficients = ldcoeff(self.metallicity, self.temperature, self.logg, 'V')
        self.ld_1 = coefficients[0]
        self.ld_2 = coefficients[1]
        self.ld_3 = coefficients[2]
        self.ld_4 = coefficients[3]

    def transit_lightcurve(self, planet, time_seq, plot=False, save=False, file_name='Lightcurve'):
        planet.test()
        lightcurve = models.transit((self.ld_1, self.ld_2, self.ld_3, self.ld_4),
                                    planet.rp_rs, planet.period, planet.a_rs, planet.eccentricity,
                                    planet.inclination, planet.omega, planet.Omega, planet.mid_transit,
                                    time_seq)
        if plot:
            plt.plot((time_seq - planet.mid_transit) / planet.period, lightcurve, 'k-', lw=2)
            plt.xlabel(r'$ phase $')
            plt.ylabel(r'$ relative \, flux $')
            plt.ylim((plt.ylim()[0], 1.002))
            if save:
                plt.savefig(file_name, dpi=200)
            else:
                plt.show()
            plt.close()
        return lightcurve
    
    def animate_transit_lightcurve(self, planet, time_seq, save=False, file_name='Animation'):
        animation.animation((self.ld_1, self.ld_2, self.ld_3, self.ld_4),
                            planet.rp_rs, planet.period, planet.a_rs, planet.eccentricity,
                            planet.inclination, planet.omega, planet.Omega, planet.mid_transit,
                            time_seq, save, file_name)


class Data:
    def __init__(self, data):
        if isinstance(data, basestring):
            self.time, self.flux = np.loadtxt(data, usecols=[0, 1], unpack=True)
        else:
            self.time, self.flux = data

    def __repr__(self):
        return '<Observation object from {0} to {1}>'.format(self.time[0], self.time[-1])

    def compare(self, star, planet):
        model = star.transit_lightcurve(planet, self.time)
        plt.plot((self.time - planet.mid_transit) / planet.period, self.flux, 'bo')
        plt.plot((self.time - planet.mid_transit) / planet.period, model, 'r-')
        plt.xlabel(r'$ phase $')
        plt.ylabel(r'$ relative \, flux $')
        plt.show()

    def transit_fitting(self, star, planet, iterations, burn, binning, rpvar, avar, ivar, fit_limb_darkening=False):
        if fit_limb_darkening:
            mcmc.transit_ld((self.time, self.flux), iterations, burn, binning,
                            (star.ld_1, star.ld_2, star.ld_3, star.ld_4),
                            planet.rp_rs, rpvar, planet.period, planet.a_rs, avar, planet.eccentricity,
                            planet.inclination, ivar, planet.omega, planet.Omega, planet.mid_transit)
        else:
            mcmc.transit((self.time, self.flux), iterations, burn, binning,
                         (star.ld_1, star.ld_2, star.ld_3, star.ld_4),
                         planet.rp_rs, rpvar, planet.period, planet.a_rs, avar, planet.eccentricity,
                         planet.inclination, ivar, planet.omega, planet.Omega, planet.mid_transit)