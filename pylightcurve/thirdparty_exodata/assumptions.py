""" This module handles any assumptions the rest of the package uses.

Im not sure how to best implement assumptions that are fully flexible to the user. This could be how the planet type
is defined or other characteristics. Im using a simple 'OK' method here, overwriting assumptions would require loading
this module and defining or changing the variables.

If anyone has a good solution to this issue, please create it in a fork or email me!
"""
import numpy as np

from . import astroquantities as aq

# TODO open an issue about this module for community discussion

planetAssumptions = {
    # Contains all planet assumptions, the key refers to the type of assumption and the format of the value can vary
    # based on the type of assumption

    'massType':
        [
            # Planet types are defined by their Mass, using inf as an absolute upper limit
            # the format of the tuples are (mass upperlimit, name). The current format dosn't allow overlaps and must be
            # in order. If you append a value run .sort() after.
            (10 * aq.M_e, 'Super-Earth'),
            (20 * aq.M_e, 'Neptune'),
            (float('inf'), 'Jupiter')
        ],

    'radiusType':  # TODO interface with rest of module
        [
            (3 * aq.R_e, 'Super-Earth'),
            (6 * aq.R_e, 'Neptune'),
            (float('inf'), 'Jupiter')
        ],

    'tempType':
        [
            (800 * aq.K, 'Cold'),
            (2000 * aq.K, 'Warm'),
            (float('inf'), 'Hot'),
        ],

    'mu':  # depends on masstype so takes the masstype as the key, and mu as the value
        {
            'Super-Earth': 18 * aq.atomic_mass_unit,  # TODO these should be more inherently linked to masstype
            'Neptune': 2.3 * aq.atomic_mass_unit,
            'Jupiter': 2 * aq.atomic_mass_unit
        },

    'albedo':  # depends on temperature so it takes tempType as the key and the albedo as the value
        {
            'Cold': 0.3,
            'Warm': 0.3,
            'Hot': 0.1,
            'Super-Earth': 0.3,
            'Neptune': 0.1,
            'Jupiter': 0.1,
        },

    'density':
    {
        'Super-Earth': 4 * aq.g / aq.cm**3,
        'Neptune': 1.6 * aq.g / aq.cm**3,
        'Jupiter': 1.3 * aq.g / aq.cm**3
    }
}


def planetMassType(mass):
    """ Returns the planet masstype given the mass and using planetAssumptions['massType']
    """

    if mass is np.nan:
        return None

    for massLimit, massType in planetAssumptions['massType']:

        if mass < massLimit:
            return massType


def planetRadiusType(radius):
    """ Returns the planet radiustype given the mass and using planetAssumptions['radiusType']
    """

    if radius is np.nan:
        return None

    for radiusLimit, radiusType in planetAssumptions['radiusType']:

        if radius < radiusLimit:
            return radiusType


def planetTempType(temperature):
    """ Returns the planet masstype given the temperature and using planetAssumptions['tempType']
    """

    for tempLimit, tempType in planetAssumptions['tempType']:

        if temperature < tempLimit:
            return tempType


def planetType(temperature, mass, radius):
    """ Returns the planet type as 'temperatureType massType'
    """

    if mass is not np.nan:
        sizeType = planetMassType(mass)
    elif radius is not np.nan:
        sizeType = planetRadiusType(radius)
    else:
        return None

    return '{0} {1}'.format(planetTempType(temperature), sizeType)


def planetMu(sizeType):
    return planetAssumptions['mu'][sizeType]


def planetAlbedo(tempType):
    return planetAssumptions['albedo'][tempType]


def planetDensity(radiusType):
     return planetAssumptions['density'][radiusType]