
__all__ = ['open_dict', 'save_dict', 'copy_dict',
           'open_yaml', 'save_yaml', 'copy_yaml',
           'open_fits', 'save_fits', 'copy_fits',
           'fits_sci', 'fits_err', 'fits_sci_err',
           'download'
           ]

import os
import ssl
import copy
import yaml
import numpy as np
import pickle
import urllib

from urllib.request import urlretrieve
from astropy.io import fits as pf


def open_dict(path):
    try:
        return pickle.load(open(path, 'rb'))
    except UnicodeDecodeError:
        return pickle.load(open(path, 'rb'), encoding='latin-1')


def save_dict(dictionary, path):
    internal_copy = copy_dict(dictionary)
    pickle.dump(dictionary, open(path, 'wb'), protocol=2)
    del internal_copy


def copy_dict(dictionary):
    return copy.deepcopy(dictionary)


def open_yaml(path):
    return yaml.load(open(path, 'r'), Loader=yaml.SafeLoader)


def save_yaml(dictionary, path):
    yaml.dump(dictionary, open(path, 'w'), default_flow_style=False)


def copy_yaml(dictionary):
    return copy.deepcopy(dictionary)


def open_fits(path):
    with pf.open(path, memmap=False) as hdulist:
        internal_copy = copy_fits(hdulist)
    internal_copy.verify('fix')
    return internal_copy


def save_fits(hdulist, path):
    internal_copy = copy_fits(hdulist)
    try:
        internal_copy.writeto(path, output_verify='fix')
    except pf.VerifyError as e:
        hdu = int(str(e.args)[4:-4].split('\\n')[1].replace('HDU ', '').replace(':', ''))
        card = int(str(e.args)[4:-4].split('\\n')[2].replace('Card ', '').replace(':', ''))
        del internal_copy[hdu].header[card]
        internal_copy.writeto(path, output_verify='fix')
    
    del internal_copy


def copy_fits(hdulist):
    return pf.HDUList([ff.copy() for ff in hdulist])


def fits_sci(hdulist):
    return np.where(np.array([ff.name for ff in hdulist]) == 'SCI')[0]


def fits_err(hdulist):
    return np.where(np.array([ff.name for ff in hdulist]) == 'ERR')[0]


def fits_sci_err(hdulist):
    return np.swapaxes([fits_sci(hdulist), fits_err(hdulist)], 0, 1)


ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE


def download(link, destination, filename=None):

    if not filename:
        filename = os.path.split(destination)[1]

    print('    Downloading {0}...'.format(filename))
    try:
        urlretrieve(link, destination)
        print('    Done!')
        return True
    except:
        try:
            with urllib.request.urlopen(link, context=ctx) as u, \
                    open(destination, 'wb') as f:
                f.write(u.read())
            print('    Done!')
            return True
        except:
            print('    Could not download {0}'.format(filename))
            return False
