from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from ._0errors import *
from ._0imports import *
from .tools_databases import *


class Databases:

    def __init__(self):

        self.clablimb = Database('clablimb', 'pylightcurve', __file__, date_to_update='191023').path
        self.ephemeris = Database('ephemeris', 'pylightcurve', __file__, date_to_update='200303').path
        self.oec = Database('oec', 'pylightcurve', __file__, date_to_update='daily').path

    def phoenix(self):
        return Database('phoenixplc3', 'pylightcurve', __file__, date_to_update='181213', ask_size='3.5GB').path


databases = Databases()
