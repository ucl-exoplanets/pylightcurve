""" Handles database classes including search functions
"""

import re
import xml.etree.ElementTree as ET
import glob
import os.path
# import io
# import gzip
# import requests

from .astroclasses import System, Binary, Star, Planet, Parameters, BinaryParameters, StarParameters, PlanetParameters

compactString = lambda string: string.replace(' ', '').replace('-', '').lower()


class OECDatabase(object):
    """ This Class Handles the OEC database including search functions.
    """

    def __init__(self, databaseLocation, stream=False):
        """ Holds the Open Exoplanet Catalogue database in python

        :param databaseLocation: file path to the Open Exoplanet Catalogue systems folder ie
            ~/git/open-exoplanet-catalogue-atmospheres/systems/
            get the catalogue from https://github.com/hannorein/open_exoplanet_catalogue
            OR the stream object (used by load_db_from_url)
        :param stream: if true treats the databaseLocation as a stream object
        """

        self._loadDatabase(databaseLocation, stream)
        self._planetSearchDict = self._generatePlanetSearchDict()

        self.systemDict = dict((system.name, system) for system in self.systems)
        self.binaryDict = dict((binary.name, binary) for binary in self.binaries)
        self.starDict = dict((star.name, star) for star in self.stars)
        self.planetDict = dict((planet.name, planet) for planet in self.planets)

    def __repr__(self):
        return 'OECDatabase({} Systems, {} Binaries, {} Stars, {} Planets)'.format(len(self.systems), len(self.binaries),
                                                                                   len(self.stars), len(self.planets))

    def searchPlanet(self, name):
        """ Searches the database for a planet. Input can be complete ie GJ1214b, alternate name variations or even
        just 1214.

        :param name: the name of the planet to search
        :return: dictionary of results as planetname -> planet object
        """

        searchName = compactString(name)
        returnDict = {}

        for altname, planetObj in self._planetSearchDict.iteritems():
            if re.search(searchName, altname):
                returnDict[planetObj.name] = planetObj

        if returnDict:
            if len(returnDict) == 1:
                return returnDict.values()[0]
            else:
                return returnDict.values()

        else:
            return False

    @property
    def transitingPlanets(self):
        """ Returns a list of transiting planet objects
        """

        transitingPlanets = []

        for planet in self.planets:
            try:
                if planet.isTransiting:
                    transitingPlanets.append(planet)
            except KeyError:  # No 'discoverymethod' tag - this also filters Solar System planets
                pass

        return transitingPlanets

    def _generatePlanetSearchDict(self):
        """ Generates a search dictionary for planets by taking all names and 'flattening' them to the most compact form
        (lowercase, no spaces and dashes)
        """

        planetNameDict = {}
        for planet in self.planets:

            name = planet.name
            altnames = planet.params['altnames']
            altnames.append(name)  # as we also want the default name to be searchable

            for altname in altnames:
                reducedname = compactString(altname)
                planetNameDict[reducedname] = planet

        return planetNameDict

    def _loadDatabase(self, databaseLocation, stream=False):
        """ Loads the database from a given file path in the class

        :param databaseLocation: the location on disk or the stream object
        :param stream: if true treats the databaseLocation as a stream object
        """

        # Initialise Database
        self.systems = []
        self.binaries = []
        self.stars = []
        self.planets = []

        if stream:
            tree = ET.parse(databaseLocation)
            for system in tree.findall(".//system"):
                self._loadSystem(system)
        else:
            databaseXML = glob.glob(os.path.join(databaseLocation, '*.xml'))
            if not len(databaseXML):
                raise LoadDataBaseError('could not find the database xml files. Have you given the correct location '
                                        'to the open exoplanet catalogues /systems folder?')

            for filename in databaseXML:
                try:
                    with open(filename, 'r') as f:
                        tree = ET.parse(f)
                except ET.ParseError as e:  # this is sometimes raised rather than the root.tag system check
                    raise LoadDataBaseError(e)

                root = tree.getroot()

                # Process the system
                if not root.tag == 'system':
                    raise LoadDataBaseError('file {0} does not contain a valid system - could be an error with your version'
                                            ' of the catalogue'.format(filename))

                self._loadSystem(root)

    def _loadSystem(self, root):
        systemParams = Parameters()
        for systemXML in root:

            tag = systemXML.tag
            text = systemXML.text
            attrib = systemXML.attrib

            systemParams.addParam(tag, text, attrib)

        system = System(systemParams.params)
        self.systems.append(system)  # Add system to the index

        self._loadBinarys(root, system)
        self._loadStars(root, system)

    def _loadBinarys(self, parentXML, parent):

        binarysXML = parentXML.findall("binary")

        for binaryXML in binarysXML:
            binaryParams = BinaryParameters()

            for value in binaryXML:

                tag = value.tag
                text = value.text
                attrib = value.attrib

                binaryParams.addParam(tag, text, attrib)

            binary = Binary(binaryParams.params)
            binary.parent = parent

            parent._addChild(binary)  # Add star to the system

            self._loadBinarys(binaryXML, binary)
            self._loadStars(binaryXML, binary)
            self._loadPlanets(binaryXML, binary)  # Load planets

            self.binaries.append(binary)  # Add star to the index

    def _loadStars(self, parentXML, parent):

        starsXML = parentXML.findall("star")

        for starXML in starsXML:
            starParams = StarParameters()

            for value in starXML:

                tag = value.tag
                text = value.text
                attrib = value.attrib

                starParams.addParam(tag, text, attrib)

            star = Star(starParams.params)
            star.parent = parent

            parent._addChild(star)  # Add star to the system

            self._loadPlanets(starXML, star)  # Load planets

            self.stars.append(star)  # Add star to the index

    def _loadPlanets(self, parentXML, parent):

        planetsXML = parentXML.findall("planet")

        for planetXML in planetsXML:

            planetParams = PlanetParameters()

            for value in planetXML:

                tag = value.tag
                text = value.text
                attrib = value. attrib

                planetParams.addParam(tag, text, attrib)

            planet = Planet(planetParams.params)
            planet.parent = parent

            parent._addChild(planet)  # Add planet to the star
            self.planets.append(planet)  # Add planet to the index


# def load_db_from_url(url="https://github.com/OpenExoplanetCatalogue/oec_gzip/raw/master/systems.xml.gz"):
#     """ Loads the database from a gzipped version of the system folder, by default the one located in the oec_gzip repo
#     in the OpenExoplanetCatalogue GitHub group.
#
#     The database is loaded from the url in memory
#
#     :param url: url to load (must be gzipped version of systems folder)
#     :return: OECDatabase objected initialised with latest OEC Version
#     """
#
#     catalogue = gzip.GzipFile(fileobj=io.BytesIO(requests.get(url).content))
#     database = OECDatabase(catalogue, stream=True)
#
#     return database


class LoadDataBaseError(IOError):
    pass