"""
Contains code for simulating observations and calculating signal to noise of
various targets.

I have changed the logic of equations in this module. Each equation is now a
class in which you can give it n-1 of the parameters and then ask for the
 remaining one via the variable.

Equations are named based on their main purpose or common name. I.e the equation
scale height is H = \\frac{k T_eff}{\mu g} even though we could use it to
calculate g, given the other parameters.

Equations are designed to be user friendly and accurate, not fast. This means
if you are using this as part of a large simulation and use exodata to generate
initial parameter you'll be fine, if however you are running a certain equation
millions of times it may be worth looking at for optimisation, without all the
 checking and overhead we do here.

**Abbreviations used in this module**

Where abbreviations are ambiguous, i.e R could be the radius of anything, we use
subscript p for planet and s for star

* R_p - Planetary Radius
* M_p - Planetary Mass
* M_s - Stellar Mass
* R_s - Stellar Radius
* H - Scale height of the Planets Atmosphere
* i - orbital inclination
* e - orbit eccentricity
* T_eff_s - Effective Temperature of the Star
* A - Albedo
* mu - mean molecular weight
"""

from __future__ import division
from numpy import sqrt, arcsin, sin, cos, log10, nan

import numpy as np
import os
import sys
from pkg_resources import resource_stream, resource_filename
import math

import quantities.constants as const
from . import astroquantities as aq
from . import params


pi = const.pi
sigma = const.Stefan_Boltzmann_constant
G = const.Newtonian_constant_of_gravitation

_rootdir = os.path.dirname(__file__)


class ExoDataEqn(object):

    def __init__(self):
        self.vars = (None,)

    def __repr__(self):
        vs = ['{}={}'.format(v, eval('self._{}'.format(v)), self)
              for v in self.vars if v is not None]
        return '{}({})'.format(self.__class__.__name__, ', '.join(vs))  # skip final ', '


class ScaleHeight(ExoDataEqn):

    def __init__(self, T_eff=None, mu=None, g=None, H=None):
        """ Uses the scale height equation to calculate a parameter given the
        others.

        .. math::
            H = \\frac{k T_eff}{\mu g}

        Where H is the scale height of the planets atmosphere, :math:`T_{eff}`
        is the planetary effective temperature,
        :math:`\mu` is the mean molecular weight of the planetary atmosphere and
        g is the planets surface gravity.

        :param T_eff: Effective temperature of the planet
        :param mu: mean molecular weight for the atmosphere
        :param g: surface gravity of the planet
        :return: H (scale Height) of the atmosphere
        """

        ExoDataEqn.__init__(self)

        self._T_eff = T_eff
        self._mu = mu
        self._g = g
        self._H = H

        self.vars = ('H', 'T_eff', 'mu', 'g')  # list of input variables

        if (T_eff, mu, g, H).count(None) > 1:
            raise EqnInputError("You must give all parameters bar one")

    @property
    def H(self):

        H = self._H
        mu = self._mu
        g = self._g
        T_eff = self._T_eff

        if H is None:
            H = (const.k * T_eff) / (mu * g)

        return H.rescale(aq.m)

    @property
    def T_eff(self):

        H = self._H
        mu = self._mu
        g = self._g
        T_eff = self._T_eff

        if T_eff is None:
            T_eff = (H * mu * g)/const.k

        return T_eff.rescale(aq.K)

    @property
    def mu(self):

        H = self._H
        mu = self._mu
        g = self._g
        T_eff = self._T_eff

        if mu is None:
            mu = (const.k * T_eff) / (H * g)

        return mu.rescale(aq.atomic_mass_unit)

    @property
    def g(self):

        H = self._H
        mu = self._mu
        g = self._g
        T_eff = self._T_eff

        if g is None:
            g = (const.k * T_eff) / (H * mu)

        return g.rescale(aq.m / aq.s**2)

class MeanPlanetTemp(ExoDataEqn):

    def __init__(self, A, T_s, R_s, a, epsilon=0.7, T_p=None):
        """ Calculate the equilibrium planet temperature

        .. math::

        assumes epsilon = 0.7 by default http://arxiv.org/pdf/1111.1455v2.pdf
        """
        ExoDataEqn.__init__(self)

        self._A = A
        self._T_s = T_s
        self._T_p = T_p
        self._R_s = R_s
        self._a = a
        self._epsilon = epsilon

        self.vars = ('A', 'T_s', 'T_p', 'R_s', 'a', 'epsilon')  # list of input variables

        if (A, T_s, R_s, a, epsilon, T_p).count(None) > 1:
            raise EqnInputError("You must give all parameters bar one")

    @property
    def T_p(self):

        T_p = self._T_p
        A = self._A
        T_s = self._T_s
        R_s = self._R_s
        a = self._a
        epsilon = self._epsilon

        if T_p is None:
            T_p = T_s * ((1-A)/epsilon)**(1/4) * sqrt(R_s/(2*a))

        return T_p.rescale(aq.degK)

    @property
    def A(self):

        T_p = self._T_p
        A = self._A
        T_s = self._T_s
        R_s = self._R_s
        a = self._a
        epsilon = self._epsilon

        if A is None:
            A = 1 - epsilon * (T_p / (T_s * sqrt(R_s/(2*a))))**4

        return A.rescale(aq.dimensionless)

    @property
    def T_s(self):

        T_p = self._T_p
        A = self._A
        T_s = self._T_s
        R_s = self._R_s
        a = self._a
        epsilon = self._epsilon

        if T_s is None:
            T_s = T_p / (((1-A)/epsilon)**(1/4) * sqrt(R_s/(2*a)))

        return T_s.rescale(aq.K)

    @property
    def R_s(self):

        T_p = self._T_p
        A = self._A
        T_s = self._T_s
        R_s = self._R_s
        a = self._a
        epsilon = self._epsilon

        if R_s is None:
           R_s = 2*a*(T_p / (T_s * ((1-A)/epsilon)**(1/4)))**2

        return R_s.rescale(aq.R_s)

    @property
    def a(self):

        T_p = self._T_p
        A = self._A
        T_s = self._T_s
        R_s = self._R_s
        a = self._a
        epsilon = self._epsilon

        if a is None:
            a = R_s/(2*(T_p / (T_s * ((1-A)/epsilon)**(1/4)))**2)

        return a.rescale(aq.au)

    @property
    def epsilon(self):

        T_p = self._T_p
        A = self._A
        T_s = self._T_s
        R_s = self._R_s
        a = self._a
        epsilon = self._epsilon

        if epsilon is None:
            epsilon = (1 - A) / (T_p / (T_s * sqrt(R_s / (2 * a)))) ** 4

        return epsilon.rescale(aq.dimensionless)


class StellarLuminosity(ExoDataEqn):

    def __init__(self, R=None, T=None, L=None):
        """ Calculate stellar luminosity

        .. math::
            L_\star = 4\pi R^2_\star \sigma T^4_\star

        Where :math:`L_\star` is the Stellar luminosity, :math:`R_\star` stellar
        radius, :math:`\sigma` Stefan-Boltzman constant, :math:`T_\star` the
        temperature of the star

        :param R: stellar radius
        :param T: effective surface temperature of the star
        :param L: Stellar luminosity
        """

        ExoDataEqn.__init__(self)

        self._L = L
        self._R = R
        self._T = T

        self.vars = ('R', 'T', 'L')  # list of input variables

        if (R, T, L).count(None) > 1:
            raise EqnInputError("You must give all parameters bar one")

    @property
    def L(self):

        L = self._L
        R = self._R
        T = self._T

        if L is None:
            L = 4 * pi * R**2 * sigma * T**4

        return L.rescale(aq.W)

    @property
    def R(self):

        L = self._L
        R = self._R
        T = self._T

        if R is None:
            R = sqrt(L / (4 * pi * sigma * T**4))

        return R.rescale(aq.R_s)

    @property
    def T(self):

        L = self._L
        R = self._R
        T = self._T

        if T is None:
            T = (L / (4 * pi * sigma * R**2))**0.25

        return T.rescale(aq.K)


class KeplersThirdLaw(ExoDataEqn):

    def __init__(self, a=None, M_s=None, P=None, M_p=0. * aq.M_j):
        """Calculates the period of the orbit using the stellar mass, planet
        mass and sma using keplers Third Law

        .. math::
            P = \sqrt{\frac{4\pi^2a^3}{G \left(M_\star + M_p \right)}}
        """

        ExoDataEqn.__init__(self)

        self._a = a
        self._M_s = M_s
        self._P = P
        self._M_p = M_p

        if (a, M_s, P).count(None) > 1:
            raise EqnInputError("You must give all parameters bar one")

    @property
    def P(self):

        a = self._a
        M_s = self._M_s
        M_p = self._M_p
        P = self._P

        if P is None:
            P = 2 * pi * sqrt(a**3 / (G * (M_s + M_p)))

        return P.rescale(aq.day)

    @property
    def a(self):

        a = self._a
        M_s = self._M_s
        M_p = self._M_p
        P = self._P

        try:
            if a is None:
                a = ((P**2 * G*(M_s + M_p))/(4*pi**2))**(1./3)
            return a.rescale(aq.au)
        except ValueError:
            return np.nan

    @property
    def M_s(self):

        a = self._a
        M_s = self._M_s
        M_p = self._M_p
        P = self._P

        if M_s is None:
            M_s  = ((4*pi**2 * a**3)/(G*P**2)) - M_p

        return M_s.rescale(aq.M_s)

    @property
    def M_p(self):

        a = self._a
        M_s = self._M_s
        M_p = self._M_p
        P = self._P

        if M_p is None:
            M_p = ((4*pi**2 * a**3)/(G*P**2)) - M_s

        return M_p.rescale(aq.M_j)

class SurfaceGravity(ExoDataEqn):

    def __init__(self, M=None, R=None, g=None):
        """ Calculates the surface acceleration due to gravity on the planet

        .. math::
            g_p = \\frac{GM_p}{R_p^2}

        where :math:`g_p` is the acceleration ude to gravity on the planet surface, G is the gravitational constant,
        :math:`M_p` planetary mass, :math:`R_p` radius of the planet

        :param params: dict containing 'M_p', 'R_p' with units

        :return: g - acceleration due to gravity * m / s**2
        """

        ExoDataEqn.__init__(self)

        self._M = M
        self._R = R
        self._g = g

        if (M, R, g).count(None) > 1:
            raise EqnInputError("You must give all parameters bar one")

    @property
    def g(self):

        M = self._M
        R = self._R
        g = self._g

        if g is None:
            g = (G * M)/(R**2)

        return g.rescale(aq.m / aq.s**2)

    @property
    def M(self):

        M = self._M
        R = self._R
        g = self._g

        if M is None:
            M = (g * R**2) / G

        return M.rescale(aq.M_j)

    @property
    def R(self):

        M = self._M
        R = self._R
        g = self._g

        if R is None:
            R = sqrt(M * G / g)

        return R.rescale(aq.R_j)

class Logg(ExoDataEqn):

    def __init__(self, M=None, R=None, logg=None):
        """ Calculates the surface acceleration due to gravity on the planet as logg, the base 10 logarithm of g in cgs
        units. This function uses :py:func:`surfaceGravity` and then rescales it to cgs and takes the log
        """

        ExoDataEqn.__init__(self)

        self._M = M
        self._R = R
        self._logg = logg

        if (M, R, logg).count(None) > 1:
            raise EqnInputError("You must give all parameters bar one")

    @property
    def logg(self):

        M = self._M
        R = self._R
        logg = self._logg

        if logg is None:
            g = SurfaceGravity(M, R).g
            logg = log10(float(g.rescale(aq.cm / aq.s**2)))  # the float wrapper is needed to remove dimensionality

        return logg

    @property
    def M(self):

        M = self._M
        R = self._R
        logg = self._logg

        if M is None:
            g = (10**logg) * aq.cm / aq.s**2
            M = SurfaceGravity(None, R, g).M

        return M.rescale(aq.M_j)

    @property
    def R(self):

        M = self._M
        R = self._R
        logg = self._logg

        if R is None:
            g = (10**logg) * aq.cm / aq.s**2
            R = SurfaceGravity(M, None, g).R

        return R.rescale(aq.R_j)

class TransitDepth(ExoDataEqn):

    def __init__(self, R_s=None, R_p=None, depth=None):

        ExoDataEqn.__init__(self)

        self._R_s = R_s
        self._R_p = R_p
        self._depth = depth

        if (R_s, R_p, depth).count(None) > 1:
            raise EqnInputError("You must give all parameters bar one")

    @property
    def depth(self):

        R_s = self._R_s
        R_p = self._R_p
        depth = self._depth

        if depth is None:
            depth = (R_p / R_s)**2

        return depth.rescale(aq.dimensionless)

    @property
    def R_s(self):

        R_s = self._R_s
        R_p = self._R_p
        depth = self._depth

        if R_s is None:
            R_s = R_p / sqrt(depth)

        return R_s

    @property
    def R_p(self):

        R_s = self._R_s
        R_p = self._R_p
        depth = self._depth

        if R_p is None:
            R_p = R_s * sqrt(depth)

        return R_p.rescale(aq.R_j)


class Density(ExoDataEqn):

    def __init__(self, M=None, R=None, density=None):
        """ Calculates the density in g/cm**3
        :param R: radius
        :param M: mass
        :return:
        """

        ExoDataEqn.__init__(self)

        self._M = M
        self._R = R
        self._density = density

        if (M, R, density).count(None) > 1:
            raise EqnInputError("You must give all parameters bar one")

    @property
    def density(self):

        M = self._M
        R = self._R
        density = self._density

        if density is None:
            volume = 4. / 3 * pi * R**3
            density = (M/volume)

        return density.rescale(aq.g / aq.cm**3)

    @property
    def M(self):

        M = self._M
        R = self._R
        density = self._density

        if M is None:
            volume = 4. / 3 * pi * R**3
            M = density * volume

        return M.rescale(aq.M_j)

    @property
    def R(self):

        M = self._M
        R = self._R
        density = self._density

        if R is None:
            R = (M/(density * 4. / 3 * pi))**(1./3)

        return R.rescale(aq.R_j)


class TransitDuration(ExoDataEqn):

    def __init__(self, P=None, a=None, Rp=None, Rs=None, i=None, e=None, w=None):
        """ Eccentric transit duration equation from Kipping (2011)

        Currently only calculates TD given the other values. There is little demand for the rearangements and they
        are non trivial

        .. math::
            T_{14} = \frac{P}{\pi} \frac{\varrho_{P,T}^2}{\sqrt{1-e^2}}
            \arcsin{\left( \sqrt{\frac{S^2_{P*} - b^2_{P,T}}{(a_P/R_\star)^2\varrho^2_{P,T} -b^2_{P,T}}}\right)}

            a_R = (a / R_\star)

            \varrho_{P,T}(f_p) = \frac{1-e^2_P}{1+e_P \sin(\omega)}

            b_{P,T} = (a_P / R_\star)\varrho_{P,T}\cos{i}

        :param P:
        :param a:
        :param Rp:
        :param Rs:
        :param i:
        :param e:
        :param w:
        """

        ExoDataEqn.__init__(self)

        self.P = P
        self.a = a
        self.Rp = Rp
        self.Rs = Rs
        self.i = i
        self.e = e
        Td = None

        if w == 0:
            self.w = 0 * aq.rad
        else:
            self.w = w

        if (P, a, Rp, Rs, i, e, w, Td).count(None) > 1:
            raise EqnInputError("You must give all parameters bar one")

    @property
    def Td(self):

        a = (self.a / self.Rs).rescale(aq.dimensionless)
        i = self.i.rescale(aq.rad)
        w = self.w.rescale(aq.rad)
        P = self.P
        RpRs = (self.Rp / self.Rs).rescale(aq.dimensionless)
        e = self.e

        ro_pt = (1-e**2)/(1+e*np.sin(w))
        b_pt = a*ro_pt*np.cos(i)
        s_ps = 1.0 + RpRs
        df = np.arcsin(np.sqrt((s_ps**2-b_pt**2)/((a**2)*(ro_pt**2)-b_pt**2)))

        duration = (P*(ro_pt**2))/(np.pi*np.sqrt(1-e**2))*df

        return duration.rescale(aq.min)


class ImpactParameter(ExoDataEqn):

    def __init__(self, a=None, R_s=None, i=None, b=None):
        """ projected distance between the planet and star centers during mid transit
        .. math::
            b \equiv \frac{a}{R_*} \cos{i}
        (Seager & Mallen-Ornelas 2003).
        """

        ExoDataEqn.__init__(self)

        self._a = a
        self._R_s = R_s
        self._i = i
        self._b = b

        if (a, R_s, i, b).count(None) > 1:
            raise EqnInputError("You must give all parameters bar one")

    @property
    def b(self):

        a = self._a
        R_s = self._R_s
        i = self._i

        b = (a / R_s) * cos(i.rescale(aq.rad))

        return b.rescale(aq.dimensionless)

    @property
    def a(self):

        b = self._b
        R_s = self._R_s
        i = self._i

        a = (b * R_s) / cos(i.rescale(aq.rad))

        return a.rescale(aq.au)

    @property
    def R_s(self):

        a = self._a
        b = self._b
        i = self._i

        R_s = (a / b) * cos(i.rescale(aq.rad))

        return R_s.rescale(aq.R_s)

    @property
    def i(self):

        a = self._a
        b = self._b
        R_s = self._R_s

        i = np.arccos((R_s * b / a).rescale(aq.dimensionless))

        return i.rescale(aq.deg)


def ratioTerminatorToStar(H_p, R_p, R_s):  # TODO add into planet class
    """ Calculates the ratio of the terminator to the star assuming 5 scale heights large. If you dont know all of the
    input try :py:func:`calcRatioTerminatorToStar`

    .. math::
        \Delta F = \\frac{10 H R_p + 25 H^2}{R_\star^2}

    Where :math:`\Delta F` is the ration of the terminator to the star, H scale height planet atmosphere,
    :math:`R_p` radius of the planet, :math:`R_s` radius of the star

    :param H_p:
    :param R_p:
    :param R_s:
    :return: ratio of the terminator to the star
    """
    # TODO let user overwrite scale heights and auto function to select on planet type?

    deltaF = ((10 * H_p * R_p) + (25 * H_p**2)) / (R_s**2)
    return deltaF.simplified


def SNRPlanet(SNRStar, starPlanetFlux, Nobs, pixPerbin, NVisits=1):
    """ Calculate the Signal to Noise Ratio of the planet atmosphere

    .. math::
        \\text{SNR}_\\text{planet} = \\text{SNR}_\\text{star} \\times \Delta F \\times \sqrt{N_\\text{obs}}
        \\times \sqrt{N_\\text{pixPerbin}} \\times \sqrt{N_\\text{visits}}

    Where :math:`\\text{SNR}_\star` SNR of the star detection, :math:`\Delta F` ratio of the terminator to the star,
    :math:`N_\\text{obs}` number of exposures per visit, :math:`N_\\text{pixPerBin}` number of pixels per wavelength bin
    , :math:`N_\\text{visits}` number of visits

    :return:
    """

    SNRplanet = SNRStar * starPlanetFlux * sqrt(Nobs) * sqrt(pixPerbin) * sqrt(NVisits)

    return SNRplanet


def calcRatioTerminatorToStar(params):  # TODO update with new format
    """ Calculates the ratio of the Terminator to the Star using :py:func:`ratioTerminatorToStar` but calculates all
    intermediary params on route.

    This is a convenience function which will largely be used over the others as the intermediary params arent often given
    in detection papers.

    :param params: dict containing ...
    :return: ratio of the terminator to the star
    """

    # check for optional params and calc if missing

    if 'H_p' not in params:
        if 'T_eff_p' not in params:
            if 'L_s' not in params:
                params['L_s'] = starLuminosity(params)
            params['T_eff_p'] = meanPlanetTemp(params)

        if 'g_p' not in params:
            params['g_p'] = surfaceGravity(params)

        params['H_p'] = scaleHeight(params)

    params['delta_F_p_s'] = ratioTerminatorToStar(params)

    return params['delta_F_p_s']


def transitDurationCircular(P, R_s, R_p, a, i):
    """ Estimation of the primary transit time. Assumes a circular orbit.

    .. math::
        T_\\text{dur} = \\frac{P}{\pi}\sin^{-1} \left[\\frac{R_\star}{a}\\frac{\sqrt{(1+k)^2 + b^2}}{\sin{a}} \\right]

    Where :math:`T_\\text{dur}` transit duration, P orbital period, :math:`R_\star` radius of the star,
    a is the semi-major axis, k is :math:`\\frac{R_p}{R_s}`, b is :math:`\frac{a}{R_*} \cos{i}` (Seager & Mallen-Ornelas 2003)

    :param i: orbital inclination
    :return:
    """

    if i is nan:
        i = 90 * aq.deg

    i = i.rescale(aq.rad)
    k = R_p / R_s  # lit reference for eclipsing binaries
    b = (a * cos(i)) / R_s

    duration = (P / pi) * arcsin(((R_s * sqrt((1 + k) ** 2 - b ** 2)) / (a * sin(i))).simplified)

    return duration.rescale(aq.min)


def estimateStellarTemperature(M_s):
    """ Estimates stellar temperature using the main sequence relationship T ~ 5800*M^0.65 (Cox 2000)??
    """
    # TODO improve with more x and k values from Cox 2000
    try:
        temp = (5800 * aq.K * float(M_s.rescale(aq.M_s) ** 0.65)).rescale(aq.K)
    except AttributeError:
        temp = np.nan
    return temp


def estimateDistance(m, M, Av=0.0):
    """ estimate the distance to star based on the absolute magnitude, apparent magnitude and the
    absorbtion / extinction

    :param m: apparent magnitude
    :param M: absolute magnitude
    :param Av: absorbtion / extinction

    :return: d (distance to object) in parsecs
    """
    try:
        m = float(m)  # basic value checking as there is no units
        M = float(M)
        Av = float(Av)
    except TypeError:
        return np.nan

    d = 10 ** ((m - M + 5 - Av) / 5)

    if math.isnan(d):
        return np.nan
    else:
        return d * aq.pc


def _createAbsMagEstimationDict():
    """ loads magnitude_estimation.dat which is from
    http://xoomer.virgilio.it/hrtrace/Sk.htm on 24/01/2014 and based on
    Schmid-Kaler (1982)

    creates a dict in the form [Classletter][ClassNumber][List of values for
    each L Class]
    """
    magnitude_estimation_filepath = resource_filename(__name__, 'data/magnitude_estimation.dat')
    raw_table = np.loadtxt(magnitude_estimation_filepath, '|S5')

    absMagDict = {'O': {}, 'B': {}, 'A': {}, 'F': {}, 'G': {}, 'K': {}, 'M': {}}
    for row in raw_table:
        if sys.hexversion >= 0x03000000:
            starClass = row[0].decode("utf-8")  # otherwise we get byte ints or b' caused by 2to3
            absMagDict[starClass[0]][int(starClass[1])] = [float(x) for x in row[1:]]
        else:
            # dict of spectral type = {abs mag for each luminosity class}
            absMagDict[row[0][0]][int(row[0][1])] = [float(x) for x in row[1:]]

    # manually typed from table headers - used to match columns with the L class (header)
    LClassRef = {'V': 0, 'IV': 1, 'III': 2, 'II': 3, 'Ib': 4, 'Iab': 5, 'Ia': 6, 'Ia0': 7}

    return absMagDict, LClassRef

absMagDict, LClassRef = _createAbsMagEstimationDict()


def estimateAbsoluteMagnitude(spectralType):
    """ Uses the spectral type to lookup an approximate absolute magnitude for the star.
    """

    from .astroclasses import SpectralType

    specType = SpectralType(spectralType)

    if specType.classLetter == '':
        return np.nan
    elif specType.classNumber == '':
        specType.classNumber = 5  # aproximation using mid magnitude value

    if specType.lumType == '':
        specType.lumType = 'V'  # assume main sequence

    LNum = LClassRef[specType.lumType]
    classNum = specType.classNumber
    classLet = specType.classLetter

    try:
        return absMagDict[classLet][classNum][LNum]
    # value not in table. Assume the number isn't there (Key p2.7, Ind p3+)
    except (KeyError, IndexError):
        try:
            classLookup = absMagDict[classLet]
            values = np.array(list(classLookup.values()))[:, LNum]  # only select the right L Type
            return np.interp(classNum, list(classLookup.keys()), values)
        except (KeyError, ValueError):
            return np.nan  # class not covered in table


def _createMagConversionDict():
    """ loads magnitude_conversion.dat which is table A% 1995ApJS..101..117K
    """
    magnitude_conversion_filepath = resource_stream(__name__, 'data/magnitude_conversion.dat')
    raw_table = np.loadtxt(magnitude_conversion_filepath, '|S5')

    magDict = {}
    for row in raw_table:
        if sys.hexversion >= 0x03000000:
            starClass = row[1].decode("utf-8")  # otherwise we get byte ints or b' caused by 2to3
            tableData = [x.decode("utf-8") for x in row[3:]]
        else:
            starClass = row[1]
            tableData = row[3:]
        magDict[starClass] = tableData

    return magDict

magDict = _createMagConversionDict()


def magKtoMagV(*args, **kwargs):
    """ Converts K magnitude to V magnitude."""
    raise DeprecationWarning(
        "This function has is been phased out in favour of the"
        "astroclasses.Magnitude class which can convert between many magnitudes"
    )


class EqnInputError(params.ExoDataError):
    pass
