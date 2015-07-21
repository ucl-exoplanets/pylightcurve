import numpy as np
import matplotlib.pyplot as plt

import ephem

import models
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
                '   mid_transit [days] = %g \n'
                % (self.rp_rs, self.period,
                    self.a_rs, self.eccentricity, self.inclination,
                    self.omega, self.Omega, self.mid_transit)
                )

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

    def next_transit(self, number):
        self.test()
        ww = self.omega * np.pi / 180
        ii = self.inclination * np.pi / 180

        ro_pt = (1 - e ** 2)/(1+e*np.sin(ww))
        b_pt = a*ro_pt*np.cos(ii)
        s_ps = 1.0 + RpRs
        df = np.arcsin(np.sqrt((s_ps**2-b_pt**2)/((a**2)*(ro_pt**2)-b_pt**2)))
        aprox = (P*(ro_pt**2))/(np.pi*np.sqrt(1-e**2))*df
        now = float(ephem.now())
        t0 = self.mid_transit - 2415020.0
        next = t0 + (int((now - t0) / self.period) + 1) * self.period
        for i in range(number):
            print ephem.date(next + i * self.period)


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
                   self.ld_1, self.ld_2, self.ld_3, self.ld_4)
                )

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

    def transit_lightcurve(self, planet, time_seq):
        self.test()
        planet.test()
        return models.transit((self.ld_1, self.ld_2, self.ld_3, self.ld_4),
                              planet.rp_rs, planet.period, planet.a_rs, planet.eccentricity,
                              planet.inclination, planet.omega, planet.Omega, planet.mid_transit,
                              time_seq
                              )
    
    def plot_transit_lightcurve(self, planet, time_seq, save=False, file_name='Lightcurve'):
        lightcurve = models.transit((self.ld_1, self.ld_2, self.ld_3, self.ld_4),
                                    planet.rp_rs, planet.period, planet.a_rs, planet.eccentricity,
                                    planet.inclination, planet.omega, planet.Omega, planet.mid_transit,
                                    time_seq
                                    )
        plt.plot((time_seq - planet.mid_transit) / planet.period, lightcurve, 'k-', lw=2)
        plt.xlabel(r'$ phase $')
        plt.ylabel(r'$ relative \, flux $')
        plt.ylim((plt.ylim()[0], 1.002))
        if save:
            plt.savefig(file_name, dpi=200)
        else:
            plt.show()
    
    def animate_transit_lightcurve(self, planet, time_seq, save=False, file_name='Animation'):
        animation.animation((self.ld_1, self.ld_2, self.ld_3, self.ld_4),
                            planet.rp_rs, planet.period, planet.a_rs, planet.eccentricity,
                            planet.inclination, planet.omega, planet.Omega, planet.mid_transit,
                            time_seq, save, file_name
                            )







##########################################################################################
class Observation:
    def __init__(self,data):
        if isinstance(data, basestring):
            self.time, self.flux = np.loadtxt(data,usecols=[0,1],unpack=True)
        else:
            self.time, self.flux = data
    
    def __repr__(self):
        return '<Observation object from {0} to {1}>'.format(self.time[0],self.time[-1])
    
    def Compare(self, Star, Planet, save=False):
        model = Star.Lightcurve(Planet, self.time)
        plt.plot((self.time-Planet.mid_transit)/Planet.period, self.flux, 'bo')
        plt.plot((self.time-Planet.mid_transit)/Planet.period, model, 'r-')
        plt.xlabel(r'$ phase $')
        plt.ylabel(r'$ relative \, flux $')
        if save:
            plt.savefig('Compare.png',dpi=200)
        plt.show()
    
    def MCMC(self, Star, Planet, Iterations, Burn, RPvar, Avar, Ivar):
        fcmodel.model_fit((self.time, self.flux), Iterations, Burn,
                    (Star.ld_1, Star.ld_2, Star.ld_3, Star.ld_4),
                    Planet.rp_rs, RPvar, Planet.period, Planet.a_rs, Avar, Planet.eccentricity, 
                    Planet.inclination, Ivar, Planet.omega, Planet.Omega, Planet.mid_transit)
##########################################################################################


