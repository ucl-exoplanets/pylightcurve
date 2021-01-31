from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os
import sys

if sys.version_info[0] > 2:
    from urllib.request import urlretrieve
else:
    from urllib import urlretrieve
    input = raw_input

import matplotlib
if os.environ.get('DISPLAY', '') == '':
    print('no display found. Using non-interactive Agg backend')
    matplotlib.use('Agg')
else:
    try:
        matplotlib.use('TkAgg')
    except ImportError:
        print('matplotlib.pyplot has been already imported. Tk features will not be supported')
        pass

import os
import sys
import copy
import glob
import gzip
import time
import emcee
import numpy as np
import scipy
import pickle
import shutil
import socket
import exodata
import exodata.astroquantities as aq
import datetime
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from scipy.interpolate import griddata
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy import interpolate
from astropy.io import fits as pf
from astropy.io import ascii
from astropy.time import Time as astrotime
from astropy.coordinates import get_sun as astrosun
from sklearn.decomposition import FastICA, PCA
from astropy import units as astrounits

from astroquery.simbad import Simbad

if int(matplotlib.__version__[0]) < 3:
    import seaborn
    seaborn.reset_orig()
