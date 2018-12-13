from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import warnings
warnings.filterwarnings("ignore",
                        message='Matplotlib is building the font cache using fc-list. This may take a moment.')
warnings.filterwarnings("ignore",
                        message='The installed version of numexpr 2.4.4 is not supported in pandas and will be not be used')
warnings.filterwarnings("ignore",
                        message='\'second\' was found  to be \'60.0\', which is not in range [0,60). '
                                'Treating as 0 sec, +1 min [astropy.coordinates.angle_utilities]')

import matplotlib
matplotlib.use('TkAgg')



import os
import sys
import glob
import gzip
import time
import emcee
import ephem
import numpy as np
import scipy
import docopt
import pickle
import shutil
import socket
import exodata
import exodata.astroquantities as aq
import seaborn
import datetime
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from scipy.interpolate import griddata
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy import interpolate
from astropy.io import fits as pf
from sklearn.decomposition import FastICA, PCA
from matplotlib import rc

if sys.version_info[0] > 2:
    from urllib.request import urlretrieve
else:
    from urllib import urlretrieve
    input = raw_input

seaborn.reset_orig()

