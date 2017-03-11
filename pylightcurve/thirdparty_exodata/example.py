""" Code to contain example and other empty classes mainly for tests
"""

# I am making a make system so that unit tests in EASE can load fake targets that will stay constant for the tests
# Due to its wider applicability im putting these here

from .astroclasses import Planet, Star, Binary, System, Parameters, PlanetParameters, StarParameters, BinaryParameters
from . import astroclasses as ac


def genExampleSystem():
    systemPar = Parameters()
    systemPar.addParam('name', 'Example System ' + str(ac._ExampleSystemCount))
    systemPar.addParam('distance', 58)
    systemPar.addParam('declination', '+04 05 06')
    systemPar.addParam('rightascension', '01 02 03')

    exampleSystem = System(systemPar.params)
    exampleSystem.flags.addFlag('Fake')

    ac._ExampleSystemCount += 1

    return exampleSystem


def genExampleBinary():
    binaryPar = BinaryParameters()
    binaryPar.addParam('name', 'Example Binary {0}AB'.format(ac._ExampleSystemCount))
    binaryPar.addParam('semimajoraxis', 10)  # TODO add realistic value
    binaryPar.addParam('period', 10)  # TODO add realistic value
    # TODO add the rest of binary parameters

    exampleBinary = Binary(binaryPar.params)
    exampleBinary.flags.addFlag('Fake')

    # generate other star
    exampleStar2 = genExampleStar('B', False)
    exampleStar2.parent = exampleBinary
    exampleBinary._addChild(exampleStar2)

    exampleSystem = genExampleSystem()
    exampleSystem._addChild(exampleBinary)
    exampleBinary.parent = exampleSystem

    return exampleBinary


def genExampleStar(binaryLetter='', heirarchy=True):
    """ generates example star, if binaryLetter is true creates a parent binary object, if heirarchy is true will create a
    system and link everything up
    """

    starPar = StarParameters()
    starPar.addParam('age', '7.6')
    starPar.addParam('magB', '9.8')
    starPar.addParam('magH', '7.4')
    starPar.addParam('magI', '7.6')
    starPar.addParam('magJ', '7.5')
    starPar.addParam('magK', '7.3')
    starPar.addParam('magV', '9.0')
    starPar.addParam('mass', '0.98')
    starPar.addParam('metallicity', '0.43')
    starPar.addParam('name', 'Example Star {0}{1}'.format(ac._ExampleSystemCount, binaryLetter))
    starPar.addParam('name', 'HD {0}{1}'.format(ac._ExampleSystemCount, binaryLetter))
    starPar.addParam('radius', '0.95')
    starPar.addParam('spectraltype', 'G5')
    starPar.addParam('temperature', '5370')

    exampleStar = Star(starPar.params)
    exampleStar.flags.addFlag('Fake')

    if heirarchy:
        if binaryLetter:
            exampleBinary = genExampleBinary()
            exampleBinary._addChild(exampleStar)
            exampleStar.parent = exampleBinary
        else:
            exampleSystem = genExampleSystem()
            exampleSystem._addChild(exampleStar)
            exampleStar.parent = exampleSystem

    return exampleStar


def genExamplePlanet(binaryLetter=''):
    """ Creates a fake planet with some defaults
    :param `binaryLetter`: host star is part of a binary with letter binaryletter
    :return:
    """

    planetPar = PlanetParameters()
    planetPar.addParam('discoverymethod', 'transit')
    planetPar.addParam('discoveryyear', '2001')
    planetPar.addParam('eccentricity', '0.09')
    planetPar.addParam('inclination', '89.2')
    planetPar.addParam('lastupdate', '12/12/08')
    planetPar.addParam('mass', '3.9')
    planetPar.addParam('name', 'Example Star {0}{1} b'.format(ac._ExampleSystemCount, binaryLetter))
    planetPar.addParam('period', '111.2')
    planetPar.addParam('radius', '0.92')
    planetPar.addParam('semimajoraxis', '0.449')
    planetPar.addParam('temperature', '339.6')
    planetPar.addParam('transittime', '2454876.344')
    planetPar.addParam('separation', '330', {'unit': 'AU'})

    examplePlanet = Planet(planetPar.params)
    examplePlanet.flags.addFlag('Fake')

    exampleStar = genExampleStar(binaryLetter=binaryLetter)
    exampleStar._addChild(examplePlanet)
    examplePlanet.parent = exampleStar

    return examplePlanet


examplePlanet = genExamplePlanet()
exampleStar = examplePlanet.parent
exampleSystem = exampleStar.parent
