import numpy as np
import matplotlib.pyplot as plt

import ephem

import models
import mcmc
import animation
import tasks


class PyLCError(BaseException):
    pass


class PyLCValueError(PyLCError):
    pass


class PyLCFilterError(PyLCError):
    pass


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
        if np.isnan(self.period):
            raise PyLCValueError("period is not set")
        if np.isnan(self.a_rs):
            raise PyLCValueError("a_rs is not set")
        if np.isnan(self.eccentricity):
            raise PyLCValueError("eccentricity is not set")
        if np.isnan(self.inclination):
            raise PyLCValueError("inclination is not set")
        if np.isnan(self.omega):
            raise PyLCValueError("omega is not set")
        if np.isnan(self.Omega):
            raise PyLCValueError("Omega is not set")
        if np.isnan(self.mid_transit):
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
    def __init__(self):
        self.temperature = np.nan
        self.metallicity = np.nan
        self.logg = np.nan
        self.ld_1 = np.nan
        self.ld_2 = np.nan
        self.ld_3 = np.nan
        self.ld_4 = np.nan
    
    def __repr__(self):
        """Present Star parameters."""
        return ('Stellar parameters: \n'
                '   temperature [K]    = %g \n'
                '   metallicity [Fe/H] = %g \n'
                '   logg        [cgs]  = %g \n'
                '   ld_1               = %g \n'
                '   ld_2               = %g \n'
                '   ld_3               = %g \n'
                '   ld_4               = %g \n'
                % (self.temperature, self.metallicity, self.logg,
                   self.ld_1, self.ld_2, self.ld_3, self.ld_4))

    def set_example(self):
        self.temperature = 6590
        self.metallicity = 0.01
        self.logg = 4.1

    def set_limb_darkening(self, coefficients):
        if isinstance(coefficients, str):
            if np.isnan(self.temperature):
                raise PyLCValueError("temperature is not set")
            if np.isnan(self.metallicity):
                raise PyLCValueError("metallicity is not set")
            if np.isnan(self.logg):
                raise PyLCValueError("logg is not set")
            filterlist = ['u', 'v', 'b', 'y', 'U', 'B', 'V', 'R', 'I', 'J', 'H', 'K']
            if coefficients not in filterlist:
                raise PyLCFilterError("Invalid filter, got {} must be in \n {}".format(coefficients, filterlist))
            else:
                coefficients = tasks.ldcoeff(self.metallicity, self.temperature, self.logg, coefficients)
                self.ld_1 = coefficients[0]
                self.ld_2 = coefficients[1]
                self.ld_3 = coefficients[2]
                self.ld_4 = coefficients[3]
        else:
            if len(coefficients) != 4:
                raise PyLCFilterError("Should give 4 coefficients")
            else:
                self.ld_1 = coefficients[0]
                self.ld_2 = coefficients[1]
                self.ld_3 = coefficients[2]
                self.ld_4 = coefficients[3]

    def test(self):
        if np.isnan(self.temperature):
            raise PyLCValueError("temperature is not set")
        if np.isnan(self.metallicity):
            raise PyLCValueError("metallicity is not set")
        if np.isnan(self.logg):
            raise PyLCValueError("logg is not set")
        if np.isnan(self.ld_1):
            raise PyLCValueError("limb darkening coefficients are not set")

    def transit_lightcurve(self, planet, time_seq, plot=False, save=False, file_name='Lightcurve'):
        self.test()
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