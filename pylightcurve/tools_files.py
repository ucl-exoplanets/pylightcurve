from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from ._1databases import *


def replace_in_file(file_path, old_string, new_string):
    r = open(str(file_path)).read()
    r = r.replace(str(old_string), str(new_string))
    w = open(str(file_path), 'w')
    w.write(r)
    w.close()


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


def open_fits(path):

    with pf.open(path) as hdulist:
        internal_copy = copy_fits(hdulist)

    return internal_copy


def get_fits_arrays(path, index, arrays):

    if isinstance(arrays, str):
        arrays = [arrays]

    out = {}

    with pf.open(path) as fits:
        for array in arrays:
            out[array] = np.ones_like(fits[index].data[array]) * fits[index].data[array]

    return out


def save_fits(hdulist, path):

    internal_copy = copy_fits(hdulist)
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
