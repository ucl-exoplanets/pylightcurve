""" Contains structural classes ie binary, star, planet etc which mimic the xml
structure with objects
"""
import sys
import math
from pkg_resources import resource_stream
import logging

import numpy as np
import astropy.coordinates
import astropy.units as u

from . import equations as eq
from . import astroquantities as aq
from . import assumptions as assum
from . import flags
from . import params as ed_params

logger = logging.getLogger('')


class _BaseObject(object):

    def __init__(self, params=None):

        self.children = []
        self.parent = False
        self.classType = 'BaseObject'
        self.flags = flags.Flags()

        self.params = {}
        if params is not None:
            self._updateParams(params)  # TODO value validator?

    def _addChild(self, child):

        self.children.append(child)

    def _updateParams(self, params):
        """ This method updates parameters allowing for any validation / unit additions in the near future
        """

        self.params.update(params)

    def _getParentClass(self, startClass, parentClass):
        """ gets the parent class by calling successive parent classes with .parent until parentclass is matched.
        """
        try:
            if not startClass:  # reached system with no hits
                raise AttributeError
        except AttributeError:  # i.e calling binary on an object without one
                raise HierarchyError('This object ({0}) has no {1} as a parent object'.format(self.name, parentClass))

        if startClass.classType == parentClass:
            return startClass
        else:
            return self._getParentClass(startClass.parent, parentClass)

    @property
    def name(self):  # TODO variable for altnames
        try:
            return self.params['name']
        except KeyError:
            try:
                return self.parent.name
            except AttributeError:
                return 'Un-named ' + self.classType
        except AttributeError:
            return 'Un-named ' + self.classType

    def __repr__(self):
        return '{0}({1!r})'.format(self.classType, self.name)

    def getParam(self, paramKey):
        """ Fetches a parameter from the params dictionary. If it's not there it will return NaN. This allows the use
        of list comprehensions over the entire planet set without KeyErrors.

        NaN was used as unlike False and None, NaN < 1 and NaN > 1 are both False
        """

        try:
            return self.params[paramKey]
        except KeyError:
            return np.NaN

    def __eq__(self, other):
        """ check the parameter dictionaries for both clases are the same (and both are of the same class)
        """

        if type(self) == type(other):
            return self.params == other.params
        else:
            return False

    @property
    def system(self):
        return self._getParentClass(self.parent, 'System')


class System(_BaseObject):

    def __init__(self, *args, **kwargs):
        _BaseObject.__init__(self, *args, **kwargs)
        self.classType = 'System'

    @property
    def ra(self):
        return self.getParam('rightascension')

    @ra.setter
    def ra(self, ra):
        self.params['rightascension'] = ra

    @property
    def dec(self):
        return self.getParam('declination')

    @dec.setter
    def dec(self, dec):
        self.params['declination'] = dec

    @property
    def d(self):
        return self.getParam('distance')

    @d.setter
    def d(self, d):
        d = d.rescale(aq.pc)
        self.params['distance'] = d

    @property
    def stars(self):
        return self.children  # TODO child could be a binary or planet

    @property
    def epoch(self):
        return self.getParam('epoch')

    @epoch.setter
    def epoch(self, epoch):
        self.params['epoch'] = epoch


class PlanetAndBinaryCommon(_BaseObject):
    def __init__(self, *args, **kwargs):
        _BaseObject.__init__(self, *args, **kwargs)
        self.classType = 'PlanetAndBinaryCommon'

    @property
    def i(self):
        return self.getParam('inclination')

    @i.setter
    def i(self, i):
        i = i.rescale(aq.deg)
        self.params['inclination'] = i

    @property
    def e(self):
        return self.getParam('eccentricity')

    @e.setter
    def e(self, e):
        self.params['eccentricity'] = e

    @property
    def P(self):
        period = self.getParam('period')
        if period is not np.nan:
            return period
        elif ed_params.estimateMissingValues:
            self.flags.addFlag('Calculated Period')
            return self.calcPeriod()
        else:
            return np.nan

    @P.setter
    def P(self, P):
        P = P.rescale(aq.day)
        self.params['period'] = P

    def calcPeriod(self):
        raise NotImplementedError('Only implemented for Binary and Planet child classes')

    @property
    def a(self):
        sma = self.getParam('semimajoraxis')
        if sma is np.nan and ed_params.estimateMissingValues:
            if self.getParam('period') is not np.nan:
                sma = self.calcSMA()  # calc using period
                self.flags.addFlag('Calculated SMA')
                return sma
            else:
                return np.nan
        else:
            return sma

    @a.setter
    def a(self, a):
        a = a.rescale(aq.au)
        self.params['semimajoraxis'] = a

    def calcSMA(self):
        raise NotImplementedError('Only implemented for Binary and Planet child classes')

    @property
    def transittime(self):
        return self.getParam('transittime')

    @transittime.setter
    def transittime(self, transittime):
        self.params['transittime'] = transittime

    @property
    def periastron(self):
        peri = self.getParam('periastron')
        if math.isnan(peri) and self.e == 0:
            peri = 0 * aq.deg
        return peri

    @periastron.setter
    def periastron(self, periastron):
        self.params['periastron'] = periastron

    @property
    def longitude(self):
        return self.getParam('longitude')

    @longitude.setter
    def longitude(self, longitude):
        self.params['longitude'] = longitude

    @property
    def ascendingnode(self):
        return self.getParam('ascendingnode')

    @ascendingnode.setter
    def ascendingnode(self, ascendingnode):
        self.params['ascendingnode'] = ascendingnode

    @property
    def separation(self):
        return self.getParam('separation')

    @separation.setter
    def seperation(self, seperation):
        self.params['seperation'] = seperation


class StarAndBinaryCommon(_BaseObject):
    def __init__(self, *args, **kwargs):
        _BaseObject.__init__(self, *args, **kwargs)
        self.classType = 'StarAndBinaryCommon'

    @property
    def magU(self):
        return self.getParam('magU')

    @magU.setter
    def magU(self, mag):
        self.params['magU'] = mag

    @property
    def magB(self):
        return self.getParam('magB')

    @magB.setter
    def magB(self, mag):
        self.params['magB'] = mag

    @property
    def magH(self):
        return self.getParam('magH')

    @magH.setter
    def magH(self, mag):
        self.params['magH'] = mag

    @property
    def magI(self):
        return self.getParam('magI')

    @magI.setter
    def magI(self, mag):
        self.params['magI'] = mag

    @property
    def magJ(self):
        return self.getParam('magJ')

    @magJ.setter
    def magJ(self, mag):
        self.params['magJ'] = mag

    @property
    def magK(self):
        return self.getParam('magK')

    @magK.setter
    def magK(self, mag):
        self.params['magK'] = mag

    @property
    def magV(self):
        return self.getParam('magV')

    @magV.setter
    def magV(self, mag):
        self.params['magV'] = mag

    @property
    def magL(self):
        return self.getParam('magL')

    @magL.setter
    def magL(self, mag):
        self.params['magL'] = mag

    @property
    def magM(self):
        return self.getParam('magM')

    @magM.setter
    def magM(self, mag):
        self.params['magM'] = mag

    @property
    def magN(self):
        return self.getParam('magN')

    @magN.setter
    def magN(self, mag):
        self.params['magN'] = mag


class StarAndPlanetCommon(_BaseObject):
    def __init__(self, *args, **kwargs):
        _BaseObject.__init__(self, *args, **kwargs)
        self.classType = 'StarAndPlanetCommon'

    @property
    def age(self):
        return self.getParam('age')

    @age.setter
    def age(self, age):
        age = age.rescale(aq.Gyear)
        self.params['age'] = age

    @property  # allows stars and planets to access system values by propagating up
    def ra(self):
        return self.parent.ra

    @ra.setter
    def ra(self, ra):
        self.parent.ra = ra

    @property
    def dec(self):
        return self.parent.dec

    @dec.setter
    def dec(self, dec):
        self.parent.dec = dec

    @property
    def d(self):
        return self.parent.d

    @d.setter
    def d(self, d):
        self.parent.d = d

    @property
    def R(self):
        return self.getParam('radius')

    @R.setter
    def R(self, R):
        self.params['radius'] = R

    @property
    def T(self):
        """ Looks for the temperature in the catalogue, if absent it calculates it using calcTemperature()

        :return: planet temperature
        """
        paramTemp = self.getParam('temperature')

        if not paramTemp is np.nan:
            return paramTemp
        elif ed_params.estimateMissingValues:
            self.flags.addFlag('Calculated Temperature')
            return self.calcTemperature()
        else:
            return np.nan

    @T.setter
    def T(self, T):
        T = T.rescale(aq.K)
        self.params['temperature'] = T

    @property
    def M(self):
        return self.getParam('mass')

    @M.setter
    def M(self, M):
        M = M.rescale(aq.M_j)
        self.params['mass'] = M

    def calcTemperature(self):
        raise NotImplementedError('Only implemented for Stars and Planet child classes')

    @property
    def binary(self):
        return self._getParentClass(self, 'Binary')

    def calcSurfaceGravity(self):

        return eq.SurfaceGravity(self.M, self.R).g

    def calcLogg(self):

        return eq.Logg(self.M, self.R).logg

    def calcDensity(self):

        if self.M is np.nan or self.R is np.nan:
            return np.nan
        else:
            return eq.Density(self.M, self.R).density


class Binary(PlanetAndBinaryCommon, StarAndBinaryCommon):  # TODO add binary methods and variables, remove unused one from starcommon

    def __init__(self, *args, **kwargs):
        StarAndBinaryCommon.__init__(self, *args, **kwargs)
        PlanetAndBinaryCommon.__init__(self, *args, **kwargs)
        self.classType = 'Binary'

    @property
    def stars(self):
        return self.children

    @property
    def d(self):
        return self.parent.d

    def calcPeriod(self):
        raise NotImplementedError  # TODO

    def calcSMA(self):
        raise NotImplementedError  # TODO


class Star(StarAndPlanetCommon, StarAndBinaryCommon):

    def __init__(self, *args, **kwargs):
        StarAndPlanetCommon.__init__(self, *args, **kwargs)
        self.classType = 'Star'

    @property
    def d(self):
        """ Note this should work from child parents as .d propergates, calculates using the star estimation method
        estimateDistance and estimateAbsoluteMagnitude
        """
        # TODO this will only work from a star or below. good thing?
        d = self.parent.d
        if ed_params.estimateMissingValues:
            if d is np.nan:
                d = self.estimateDistance()
                if d is not np.nan:
                    self.flags.addFlag('Estimated Distance')
            return d
        else:
            return np.nan

    def calcLuminosity(self):

        return eq.StellarLuminosity(self.R, self.T).L

    def calcTemperature(self):
        """ uses equations.starTemperature to estimate temperature based on main sequence relationship
        """
        return eq.estimateStellarTemperature(self.M)

    def _get_or_convert_magnitude(self, mag_letter):
        """ Takes input of the magnitude letter and ouputs the magnitude fetched from the catalogue or a converted value
        :return:
        """
        allowed_mags = "UBVJIHKLMN"
        catalogue_mags = 'BVIJHK'

        if mag_letter not in allowed_mags or not len(mag_letter) == 1:
            raise ValueError("Magnitude letter must be a single letter in {0}".format(allowed_mags))

        mag_str = 'mag'+mag_letter
        mag_val = self.getParam(mag_str)

        if isNanOrNone(mag_val) and ed_params.estimateMissingValues:  # then we need to estimate it!
            # old style dict comprehension for python 2.6
            mag_dict = dict(('mag'+letter, self.getParam('mag'+letter)) for letter in catalogue_mags)
            mag_class = Magnitude(self.spectralType, **mag_dict)
            try:
                mag_conversion = mag_class.convert(mag_letter)
                # logger.debug('Star Class: Conversion to {0} successful, got {1}'.format(mag_str, mag_conversion))
                self.flags.addFlag('Estimated mag{0}'.format(mag_letter))
                return mag_conversion
            except ValueError as e:  # cant convert
                logger.exception(e)
                # logger.debug('Cant convert to {0}'.format(mag_letter))
                return np.nan
        else:
            # logger.debug('returning {0}={1} from catalogue'.format(mag_str, mag_val))
            return mag_val

    @property
    def magU(self):
        return self._get_or_convert_magnitude('U')

    @property
    def magB(self):
        return self._get_or_convert_magnitude('B')

    @property
    def magV(self):
        return self._get_or_convert_magnitude('V')

    @property
    def magJ(self):
        return self._get_or_convert_magnitude('J')

    @property
    def magI(self):
        return self._get_or_convert_magnitude('I')

    @property
    def magH(self):
        return self._get_or_convert_magnitude('H')

    @property
    def magK(self):
        return self._get_or_convert_magnitude('K')

    @property
    def magL(self):
        return self._get_or_convert_magnitude('L')

    @property
    def magM(self):
        return self._get_or_convert_magnitude('M')

    @property
    def magN(self):
        return self._get_or_convert_magnitude('N')

    @property
    def Z(self):
        return self.getParam('metallicity')

    @Z.setter
    def Z(self, Z):
        self.params['metallicity'] = Z

    @property
    def spectralType(self):
        return self.getParam('spectraltype')

    @spectralType.setter
    def spectralType(self, spectraltype):
        self.params['spectraltype'] = spectraltype

    @property
    def planets(self):
        return self.children

    def getLimbdarkeningCoeff(self, wavelength=1.22):  # TODO replace with pylightcurve
        """ Looks up quadratic limb darkening parameter from the star based on T, logg and metalicity.

        :param wavelength: microns
        :type wavelength: float

        :return: limb darkening coefficients 1 and 2
        """
        # TODO check this returns correct value - im not certain
        # The intervals of values in the tables
        tempind = [ 3500., 3750., 4000., 4250., 4500., 4750., 5000., 5250., 5500., 5750., 6000., 6250.,
                 6500., 6750., 7000., 7250., 7500., 7750., 8000., 8250., 8500., 8750., 9000., 9250.,
                 9500., 9750., 10000., 10250., 10500., 10750., 11000., 11250., 11500., 11750., 12000., 12250.,
                 12500., 12750., 13000., 14000., 15000., 16000., 17000., 19000., 20000., 21000., 22000., 23000.,
                 24000., 25000., 26000., 27000., 28000., 29000., 30000., 31000., 32000., 33000., 34000., 35000.,
                 36000., 37000., 38000., 39000., 40000., 41000., 42000., 43000., 44000., 45000., 46000., 47000.,
                 48000., 49000., 50000.]
        lggind = [0., 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5.]
        mhind = [-5., -4.5, -4., -3.5, -3., -2.5, -2., -1.5, -1., -0.5, -0.3, -0.2, -0.1, 0., 0.1, 0.2, 0.3, 0.5, 1.]

        # Choose the values in the table nearest our parameters
        tempselect = _findNearest(tempind, float(self.T))
        lgselect = _findNearest(lggind, float(self.calcLogg()))
        mhselect = _findNearest(mhind, float(self.Z))

        quadratic_filepath = resource_stream(__name__, 'data/quadratic.dat')
        coeffTable = np.loadtxt(quadratic_filepath)

        foundValues = False
        for i in range(len(coeffTable)):
            if coeffTable[i, 2] == lgselect and coeffTable[i, 3] == tempselect and coeffTable[i, 4] == mhselect:
                if coeffTable[i, 0] == 1:
                    u1array = coeffTable[i, 8:]  # Limb darkening parameter u1 for each wl in waveind
                    u2array = coeffTable[i+1, 8:]
                    foundValues = True
                    break

        if not foundValues:
            raise ValueError('No limb darkening values could be found')  # TODO replace with better exception

        waveind = [0.365, 0.445, 0.551, 0.658, 0.806, 1.22, 1.63, 2.19, 3.45]  # Wavelengths available in table

        # Interpolates the value at wavelength from values in the table (waveind)
        u1AtWavelength = np.interp(wavelength, waveind, u1array, left=0, right=0)
        u2AtWavelength = np.interp(wavelength, waveind, u2array, left=0, right=0)

        return u1AtWavelength, u2AtWavelength

    def estimateAbsoluteMagnitude(self):
        return eq.estimateAbsoluteMagnitude(self.spectralType)

    def estimateDistance(self):
        # TODO handle other mags than V
        if self.magV is not np.nan:
            return eq.estimateDistance(self.magV, self.estimateAbsoluteMagnitude())
        else:
            return np.nan


class Planet(StarAndPlanetCommon, PlanetAndBinaryCommon):

    def __init__(self, *args, **kwargs):
        StarAndPlanetCommon.__init__(self, *args, **kwargs)
        PlanetAndBinaryCommon.__init__(self, *args, **kwargs)
        self.classType = 'Planet'

    @property
    def isTransiting(self):
        """ Checks the the istransiting tag to see if the planet transits. Note that this only works as of catalogue
        version  ee12343381ae4106fd2db908e25ffc537a2ee98c (11th March 2014) where the istransiting tag was implemented
        """
        try:
            isTransiting = self.params['istransiting']
        except KeyError:
            return False

        if isTransiting == '1':
            return True
        else:
            return False

    def calcTransitDuration(self, circular=False):
        """ Estimation of the primary transit time assuming a circular orbit (see :py:func:`equations.transitDuration`)
        """

        try:
            if circular:
                return eq.transitDurationCircular(self.P, self.star.R, self.R, self.a, self.i)
            else:
                return eq.TransitDuration(self.P, self.a, self.R, self.star.R, self.i, self.e, self.periastron).Td
        except (ValueError,
                AttributeError,  # caused by trying to rescale nan i.e. missing i value
                HierarchyError):  # i.e. planets that dont orbit stars
            return np.nan

    def calcScaleHeight(self):
        raise NotImplementedError
        # return eq.scaleHeight(self.T, , self.g)  # TODO mu based on assumptions

    def calcTransitDepth(self):
        return eq.TransitDepth(self.star.R, self.R).depth

    def type(self):
        return assum.planetType(self.T, self.M, self.R)

    def massType(self):
        return assum.planetMassType(self.M)

    def radiusType(self):
        return assum.planetRadiusType(self.R)

    def tempType(self):
        return assum.planetTempType(self.T)

    @property
    def mu(self):  # TODO make getter look in params first calc if not

        molweight = self.getParam('molweight')

        if molweight is np.nan:  # Use assumptions
            if self.M is not np.nan:
                return assum.planetMu(self.massType())
            elif self.R is not np.nan:
                return assum.planetMu(self.radiusType())
            else:
                return np.nan
        else:
            return molweight

    @mu.setter
    def mu(self, mu):
        mu = mu.rescale(aq.atomic_mass_unit)
        self.params['moleight'] = mu

    @property
    def albedo(self):
        albedo = self.getParam('albedo')
        if albedo is not np.nan:
            return albedo
        elif self.getParam('temperature') is not np.nan:
            planetClass = self.tempType()
        elif self.M is not np.nan:
            planetClass = self.massType()
        elif self.R is not np.nan:
            planetClass = self.radiusType()
        else:
            return np.nan

        return assum.planetAlbedo(planetClass)

    @albedo.setter
    def albedo(self, albedo):
        albedo = albedo
        self.params['albedo'] = albedo

    def calcTemperature(self):
        """ Calculates the temperature using which uses equations.MeanPlanetTemp, albedo assumption and potentially
        equations.starTemperature.

        issues
        - you cant get the albedo assumption without temp but you need it to calculate the temp.
        """
        try:
            return eq.MeanPlanetTemp(self.albedo, self.star.T, self.star.R, self.a).T_p
        except (ValueError, HierarchyError):  # ie missing value (.a) returning nan
            return np.nan

    def estimateMass(self):

        density = assum.planetDensity(self.radiusType())

        return eq.Density(None, self.R, density).M

    def calcSMA(self):
        """ Calculates the semi-major axis from Keplers Third Law
        """
        try:
            return eq.KeplersThirdLaw(None, self.star.M, self.P).a
        except HierarchyError:
            return np.nan

    def calcSMAfromT(self, epsilon=0.7):
        """ Calculates the semi-major axis based on planet temperature
        """

        return eq.MeanPlanetTemp(self.albedo(), self.star.T, self.star.R, epsilon, self.T).a

    def calcPeriod(self):
        """ calculates period using a and stellar mass
        """

        return eq.KeplersThirdLaw(self.a, self.star.M).P

    @property
    def discoveryMethod(self):
        return self.getParam('discoverymethod')

    @discoveryMethod.setter
    def discoveryMethod(self, discoverymethod):
        self.params['discoverymethod'] = discoverymethod

    @property
    def discoveryYear(self):
        try:
            return int(self.getParam('discoveryyear'))
        except ValueError:  # np.nan
            return self.getParam('discoveryyear')

    @discoveryYear.setter
    def discoveryYear(self, discoveryYear):
        self.params['discoveryyear'] = discoveryYear

    @property
    def lastUpdate(self):
        return self.getParam('lastupdate')

    @property
    def description(self):
        return self.getParam('description')

    @property
    def star(self):
        return self._getParentClass(self.parent, 'Star')


class Parameters(object):  # TODO would this subclassing dict be more preferable?
    """ A class to hold parameter dictionaries, the input can be validated, units added and handling of multi valued
    fields. In future this may be better as a child of dict.
    """

    def __init__(self):

        self.params = {
            'altnames': [],
            'list': [],
        }

        self._defaultUnits = {  # this holds quantities with no current or projected ambiguity about their unit
            'age': aq.Gyear,
            'distance': aq.pc,  # TODO more specific unit handling here or in classes?
            'magB': 1,
            'magH': 1,
            'magI': 1,
            'magJ': 1,
            'magK': 1,
            'magV': 1,
            'temperature': aq.K,
        }

        self.rejectTags = ('system', 'binary', 'star', 'planet', 'moon')  # These are handled in their own classes

    def addParam(self, key, value, attrib=None):
        """ Checks the key dosnt already exist, adds alternate names to a seperate list

        Future
            - format input and add units
            - logging
        """

        if key in self.rejectTags:
            return False  # TODO Replace with exception

        # Temporary code to handle the seperation tag than can occur several times with different units.
        # TODO code a full multi unit solution (github issue #1)
        if key == 'separation':
            if attrib is None:
                return False  # reject seperations without a unit
            try:
                if not attrib['unit'] == 'AU':
                    return False  # reject for now
            except KeyError:  # a seperation attribute exists but not one for units
                return False

        if key in self.params:  # if already exists

            if key == 'name':
                try:  # if flagged as a primary or popular name use this one, an option should be made to use either
                    if attrib['type'] == 'pri':  # first names or popular names.
                        oldname = self.params['name']
                        self.params['altnames'].append(oldname)
                        self.params['name'] = value
                    else:
                        self.params['altnames'].append(value)
                except (KeyError, TypeError):  # KeyError = no type key in attrib dict, TypeError = not a dict
                    self.params['altnames'].append(value)
            elif key == 'list':
                self.params['list'].append(value)
            else:
                try:
                    name = self.params['name']
                except KeyError:
                    name = 'Unnamed'
                print('rejected duplicate {0}: {1} in {2}'.format(key, value, name))  # TODO: log rejected value
                return False  # TODO Replace with exception

        else:  # If the key doesn't already exist and isn't rejected

            # Some tags have no value but a upperlimit in the attributes
            if value is None and attrib is not None:
                try:
                    value = attrib['upperlimit']
                except KeyError:
                    try:
                        value = attrib['lowerlimit']
                    except KeyError:
                        return False

            if key == 'rightascension':
                value = _ra_string_to_unit(value)
            elif key == 'declination':
                value = _dec_string_to_unit(value)
            elif key in self._defaultUnits:
                try:
                    value = float(value) * self._defaultUnits[key]
                except:
                    print('caught an error with {0} - {1}'.format(key, value))
            self.params[key] = value


class BinaryParameters(Parameters):

    def __init__(self):

        Parameters.__init__(self)

        self._defaultUnits.update({
            'separation': aq.au,  # TODO there is actually 2 different measurements (other is arcsec)
            'periastron': aq.deg,
        })


class StarParameters(Parameters):

    def __init__(self):

        Parameters.__init__(self)

        self._defaultUnits.update({
            'mass': aq.M_s,
            'metallicity': 1,
            'radius': aq.R_s,
        })


class PlanetParameters(Parameters):

    def __init__(self):

        Parameters.__init__(self)

        self._defaultUnits.update({
            'discoveryyear': 1,
            'mass': aq.M_j,
            'radius': aq.R_j,
            'inclination': aq.deg,
            'eccentricity': 1,
            'periastron': aq.deg,
            'period': aq.day,
            'semimajoraxis': aq.au,
            'transittime': aq.JD,  # TODO specific JD, MJF etc
            'molweight': aq.atomic_mass_unit,
            'separation': aq.au,  # TODO there is actually 2 different measurements (other is arcsec)
        })


def _findNearest(arr, value):
    """ Finds the value in arr that value is closest to
    """
    arr = np.array(arr)
    # find nearest value in array
    idx = (abs(arr-value)).argmin()
    return arr[idx]


class SpectralType(object):
    """ Takes input of a spectral type as a string and interprets it into the luminosity class and stellar type.

    .. usage :
        self.lumType = Luminosity Class
        self.classLetter = Stellar Class (ie O B A etc)
        self.classNumber = Stellar Class number
        self.specClass = ie A8V will be A8
        self.specType = ie A*V will be A8V (default for calling the class)
        self.original = the original string

    This class ignores spaces, only considers the first class if given multiple options (ie K0/K1V, GIV/V, F8-G0)
    ignores non-typical star classes (ie ) and ignores extra statements like G8 V+
    """

    def __init__(self, classString):
        self.original = classString
        self.lumType = ''
        self.classLetter = ''
        self.classNumber = ''

        self._parseSpecType(classString)

    @property
    def specClass(self):
        """ Spectral class ie A8V is A8 """
        return self.classLetter + str(self.classNumber)

    @property
    def roundedSpecClass(self):
        """ Spectral class with rounded class number ie A8.5V is A9 """
        try:
            classnumber = str(int(np.around(self.classNumber)))
        except TypeError:
            classnumber = str(self.classNumber)

        return self.classLetter + classnumber

    @property
    def specType(self):
        """ Spectral class ie A8V is A8V """
        return self.classLetter + str(self.classNumber) + self.lumType

    @property
    def roundedSpecType(self):
        """ Spectral class with rounded class number ie A8.5V is A9V """

        return self.roundedSpecClass + self.lumType

    def __repr__(self):
        return self.specType

    def _parseSpecType(self, classString):
        """ This class attempts to parse the spectral type. It should probably use more advanced matching use regex
        """

        try:
            classString = str(classString)
        except UnicodeEncodeError:
            # This is for the benefit of 1RXS1609 which currently has the spectral type K7\pm 1V
            # TODO add unicode support and handling for this case / ammend the target
            return False

        # some initial cases
        if classString == '' or classString == 'nan':
            return False

        possNumbers = range(10)
        possLType = ('III', 'II', 'Iab', 'Ia0', 'Ia', 'Ib', 'IV', 'V')  # in order of unique matches

        # remove spaces, remove slashes
        classString = classString.replace(' ', '')

        classString = classString.replace('-', '/')
        classString = classString.replace('\\', '/')
        classString = classString.split('/')[0]  # TODO we do not consider slashed classes yet (intemediates)

        # check first 3 chars for spectral types
        stellarClass = classString[:3]
        if stellarClass in _possSpectralClasses:
            self.classLetter = stellarClass
        elif stellarClass[:2] in _possSpectralClasses:  # needed because A5V wouldnt match before
            self.classLetter = stellarClass[:2]
        elif stellarClass[0] in _possSpectralClasses:
            self.classLetter = stellarClass[0]
        else:
            return False  # assume a non standard class and fail

        # get number
        try:
            numIndex = len(self.classLetter)
            classNum = int(classString[numIndex])
            if classNum in possNumbers:
                self.classNumber = int(classNum)  # don't consider decimals here, done at the type check
                typeString = classString[numIndex+1:]
            else:
                return False  # invalid number received
        except IndexError:  # reached the end of the string
            return True
        except ValueError:  # i.e its a letter - fail # TODO multi letter checking
            typeString = classString[1:]

        if typeString == '':  # ie there is no more information as in 'A8'
            return True

        # Now check for a decimal and handle those cases
        if typeString[0] == '.':
            # handle decimal cases, we check each number in turn, add them as strings and then convert to float and add
            # to original number
            decimalNumbers = '.'
            for number in typeString[1:]:
                try:
                    if int(number) in possNumbers:
                        decimalNumbers += number
                    else:
                        print('Something went wrong in decimal checking') # TODO replace with logging
                        return False # somethings gone wrong
                except ValueError:
                    break  # recevied a non-number (probably L class)
            #  add decimal to classNum
            try:
                self.classNumber += float(decimalNumbers)
            except ValueError: # probably trying to convert '.' to a float
                pass
            typeString = typeString[len(decimalNumbers):]
            if len(typeString) is 0:
                return True

        # Handle luminosity class
        for possL in possLType:  # match each possible case in turn (in order of uniqueness)
            Lcase = typeString[:len(possL)]  # match from front with length to minimise matching say IV in '<3 CIV'
            if possL == Lcase:
                self.lumType = possL
                return True

        if not self.classNumber == '':
            return True
        else:  # if there no number asumme we have a name ie 'Catac. var.'
            self.classLetter = ''
            self.classNumber = ''
            self.lumType = ''
            return False

_ExampleSystemCount = 1  # Used by example.py - put here to enable global

#               main sequence
_possSingleLetterClasses = ('O', 'B', 'A', 'F', 'G', 'K', 'M',
               'L', 'T', 'Y',  # dwarfs
               'C', 'S',
               'W',  # Wolf-Rayet
               'P', 'Q',  # Non-stellar spectral types
)
# skipped carbon stars with dashes ie C-R
_possMultiLetterClasses = ('WNE', 'WNL', 'WCE', 'WCL', 'WO', 'WR', 'WN', 'WC',  # Wolf-Rayet stars, WN/C skipped
                          'MS', 'MC',  # intermediary carbon-related classes
                          'DAB', 'DAO', 'DAZ', 'DBZ',  # Extended white dwarf spectral types
                          'DAV', 'DBV', 'DCV',  # Variable star designations, GW Vir (DOV and PNNV) skipped
                          'DA', 'DB', 'DO', 'DQ', 'DZ', 'DC', 'DX',  # white dwarf spectral types
                          )

_possSpectralClasses = _possMultiLetterClasses + _possSingleLetterClasses  # multi first


class Magnitude(object):
    """ Holds measured magnitudes and can convert between them given a spectral class.
    """

    def __init__(self, spectral_type, magU=None, magB=None, magV=None, magI=None, magJ=None, magH=None, magK=None, magL=None,
                 magM=None, magN=None):

        if isinstance(spectral_type, SpectralType):
            self.spectral_type = spectral_type
        else:
            self.spectral_type = SpectralType(spectral_type)

        self.magU = magU
        self.magB = magB
        self.magV = magV
        self.magI = magI
        self.magJ = magJ
        self.magH = magH
        self.magK = magK
        self.magL = magL
        self.magM = magM
        self.magN = magN

        # For magDict, these should probably be grouped together
        self.column_for_V_conversion = {
        #   mag  column,  sign (most are V-Mag (+1), some are Mag-V (-1))
            'U': (2,  -1),
            'B': (3,  -1),
            'J': (8, +1),
            'H': (9, +1),
            'K': (10, +1),
            'L': (11, +1),
            'M': (12, +1),
            'N': (13, +1),
        }

    def convert(self, to_mag, from_mag=None):
        """ Converts magnitudes using UBVRIJHKLMNQ photometry in Taurus-Auriga (Kenyon+ 1995)
         ReadMe+ftp1995ApJS..101..117K Colors for main-sequence stars

         If from_mag isn't specified the program will cycle through provided magnitudes and choose one. Note that all
         magnitudes are first converted to V, and then to the requested magnitude.

        :param to_mag: magnitude to convert to
        :param from_mag: magnitude to convert from
        :return:
        """
        allowed_mags = "UBVJIHKLMN"

        if from_mag:
            if to_mag == 'V':  # If V mag is requested (1/3) - from mag specified
                return self._convert_to_from('V', from_mag)
            if from_mag == 'V':
                magV = self.magV
            else:
                magV = self._convert_to_from('V', from_mag)

            return self._convert_to_from(to_mag, 'V', magV)

        # if we can convert from any magnitude, try V first
        elif not isNanOrNone(self.magV):
            if to_mag == 'V':  # If V mag is requested (2/3) - no need to convert
                return self.magV
            else:
                return self._convert_to_from(to_mag, 'V', self.magV)
        else:  # Otherwise lets try all other magnitudes in turn
            order = "UBJHKLMN"  # V is the intermediate step from the others, done by default if possible
            for mag_letter in order:
                try:
                    magV = self._convert_to_from('V', mag_letter)
                    if to_mag == 'V':  # If V mag is requested (3/3) - try all other mags to convert
                        logging.debug('Converted to magV from {0} got {1}'.format(mag_letter, magV))
                        return magV
                    else:
                        mag_val = self._convert_to_from(to_mag, 'V', magV)
                        logging.debug('Converted to mag{0} from {1} got {2}'.format(to_mag, mag_letter, mag_val))
                        return mag_val
                except ValueError:
                    continue  # this conversion may not be possible, try another

            raise ValueError('Could not convert from any provided magnitudes')

    def _convert_to_from(self, to_mag, from_mag, fromVMag=None):
        """ Converts from or to V mag using the conversion tables

        :param to_mag: uppercase magnitude letter i.e. 'V' or 'K'
        :param from_mag: uppercase magnitude letter i.e. 'V' or 'K'
        :param fromVMag: MagV if from_mag is 'V'

        :return:  estimated magnitude for to_mag from from_mag
        """
        lumtype = self.spectral_type.lumType

        # rounds decimal types, TODO perhaps we should interpolate?
        specClass = self.spectral_type.roundedSpecClass

        if not specClass:  # TODO investigate implications of this
            raise ValueError('Can not convert when no spectral class is given')

        if lumtype not in ('V', ''):
            raise ValueError("Can only convert for main sequence stars. Got {0} type".format(lumtype))

        if to_mag == 'V':
            col, sign = self.column_for_V_conversion[from_mag]

            try:  # TODO replace with pandas table
                offset = float(magDict[specClass][col])
            except KeyError:
                raise ValueError('No data available to convert those magnitudes for that spectral type')

            if math.isnan(offset):
                raise ValueError('No data available to convert those magnitudes for that spectral type')
            else:
                from_mag_val = self.__dict__['mag'+from_mag]  # safer than eval
                if isNanOrNone(from_mag_val):
                    # logger.debug('2 '+from_mag)
                    raise ValueError('You cannot convert from a magnitude you have not specified in class')
                return from_mag_val + (offset*sign)
        elif from_mag == 'V':
            if fromVMag is None:
                # trying to second guess here could mess up a K->B calulation by using the intermediate measured V. While
                # this would probably be preferable it is not was was asked and therefore could give unexpected results
                raise ValueError('Must give fromVMag, even if it is self.magV')

            col, sign = self.column_for_V_conversion[to_mag]
            try:
                offset = float(magDict[specClass][col])
            except KeyError:
                raise ValueError('No data available to convert those magnitudes for that spectral type')

            if math.isnan(offset):
                raise ValueError('No data available to convert those magnitudes for that spectral type')
            else:
                return fromVMag + (offset*sign*-1)  # -1 as we are now converting the other way
        else:
            raise ValueError('Can only convert from and to V magnitude. Use .convert() instead')


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


def isNanOrNone(val):
    """ Tests if val is float('nan') or None using math.isnan and is None. Needed as isnan fails if a non float is given.
    :param val:
    :return:
    """

    if val is None:
        return True
    else:
        try:
            return math.isnan(val)
        except TypeError:  # not a float
            return False


def _ra_string_to_unit(ra_string):

    ra_split = ra_string.split(' ')
    hour, min, sec = ra_split
    ra_astropy_format = '{}h{}m{}s'.format(hour, min, sec)

    ra_unit = astropy.coordinates.Longitude(ra_astropy_format, unit=u.deg)

    return ra_unit


def _dec_string_to_unit(dec_string):

    deg_split = dec_string.split(' ')
    deg, arcmin, arcsec = deg_split
    deg_astropy_format = '{}d{}m{}s'.format(deg, arcmin, arcsec)

    dec_unit = astropy.coordinates.Latitude(deg_astropy_format, unit=u.deg)

    return dec_unit


class HierarchyError(ed_params.ExoDataError):
    pass

