
__all__ = ['open_dict', 'save_dict', 'copy_dict',
           'open_yaml', 'save_yaml', 'copy_yaml',
           'open_fits', 'save_fits', 'copy_fits',
           'fits_sci', 'fits_err', 'fits_sci_err',
           'download', 'zip_files', 'zip_directory',
           'open_dict_online'
           ]

import os
import ssl
import copy
import yaml
import numpy as np
import pickle
import urllib
import zipfile

from astropy.io import fits as pf
from urllib.request import urlretrieve

from .counter import Counter


def open_dict(path):

    class Dummy(object):
        def __init__(self, *argv, **kwargs):
            pass

    _ = Dummy(5)

    class Unpickler(pickle._Unpickler):
        def find_class(self, module, name):
            try:
                return super().find_class(module, name)
            except Exception as e:
                print(e)
                return Dummy

    with open(path, 'rb') as f:
        unpickler = Unpickler(f)
        return unpickler.load()


def save_dict(dictionary, path):
    internal_copy = copy_dict(dictionary)
    pickle.dump(internal_copy, open(path, 'wb'), protocol=2)
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
        try:
            internal_copy.verify('fix')
        except pf.VerifyError as e:
            hdu = int(str(e.args)[4:-4].split('\\n')[1].replace('HDU ', '').replace(':', ''))
            card = int(str(e.args)[4:-4].split('\\n')[2].replace('Card ', '').replace(':', ''))
            del internal_copy[hdu].header[card]
            internal_copy.verify('fix')
        for i in internal_copy:
            while i.header[-1] == '':
                i.header = i.header[:-1]
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


def download(link, destination, filename=None, verbose=True):

    if not filename:
        filename = os.path.split(destination)[1]

    if verbose:
        print('    Downloading {0}...'.format(filename))
        print('           from {0} '.format(link))
    try:
        with urllib.request.urlopen(link, context=ctx) as u, \
                open(destination, 'wb') as f:
            f.write(u.read())
        if verbose:
            print('    Done!')
        return True
    except Exception as e:
        print('ERROR: {0}\n    Could not download {1}', e, filename)
        return False


def open_dict_online(link):

    try:
        return pickle.load(urllib.request.urlopen(link, context=ctx))
    except Exception as e:
        print('Could not open dictionary {0}: {1}'.format(link, e))
        return False


def zip_files(list_of_files, output):

    zip = zipfile.ZipFile(output, 'w', zipfile.ZIP_DEFLATED)
    counter = Counter('Zipping', len(list_of_files))
    for file in list_of_files:
        zip.write(file)
        counter.update('Adding file: {0}'.format(file))

    zip.close()


def zip_directory(directory, output):

    list_of_files = []

    for root, directories, files in os.walk(directory):
        for filename in files:
            filepath = os.path.join(root, filename)
            list_of_files.append(filepath)

    zip_files(list_of_files, output)
