import os

__version__ = open(os.path.join(os.path.abspath(os.path.dirname(__file__)), '__version__.txt')).read()
__author__ = open(os.path.join(os.path.abspath(os.path.dirname(__file__)), '__author__.txt')).read()
__author_email__ = open(os.path.join(os.path.abspath(os.path.dirname(__file__)), '__author_email__.txt')).read()
__description__ = open(os.path.join(os.path.abspath(os.path.dirname(__file__)), '__description__.txt')).read()
__url__ = open(os.path.join(os.path.abspath(os.path.dirname(__file__)), '__url__.txt')).read()

from pylightcurve.analysis.distributions import *
from pylightcurve.analysis.numerical_integration import *
from pylightcurve.analysis.optimisation import *
from pylightcurve.analysis.gaussian import *

from pylightcurve.catalogues.catalogues import *

from pylightcurve.images.find_stars import *

from pylightcurve.models.exoplanet import *
from pylightcurve.models.exoplanet_lc import *

from pylightcurve.processes.counter import *
from pylightcurve.processes.files import *

from pylightcurve.spacetime.angles import *
from pylightcurve.spacetime.targets import *

from pylightcurve.errors import *
