
import os
import exoclock

__version__ = open(os.path.join(os.path.abspath(os.path.dirname(__file__)), '__version__.txt')).read()
__author__ = open(os.path.join(os.path.abspath(os.path.dirname(__file__)), '__author__.txt')).read()
__author_email__ = open(os.path.join(os.path.abspath(os.path.dirname(__file__)), '__author_email__.txt')).read()
__description__ = open(os.path.join(os.path.abspath(os.path.dirname(__file__)), '__description__.txt')).read()
__url__ = open(os.path.join(os.path.abspath(os.path.dirname(__file__)), '__url__.txt')).read()


from .analysis.distributions import *
from .analysis.numerical_integration import *
from .analysis.gaussian import *
from .analysis.statistics import *
from .analysis.optimisation import *

from .models.exoplanet_lc import *
from .models.exoplanet import *

from .images.find_stars import *

from .processes.counter import *
from .processes.files import *
from .plots.plots_fitting import *

from .errors import *
from .databases import *

# for compatibility with iraclis v1.5
FixedTarget = exoclock.FixedTarget
Degrees = exoclock.Degrees

def locate_planet(ra, dec, radius):
    return get_planet(
        exoclock.locate_planet(exoclock.Degrees(ra), exoclock.Degrees(dec), exoclock.Degrees(radius))['name']
        )