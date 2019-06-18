from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from ._1databases import *


def central_crop(original_array, destination_fits):
    crop1 = len(original_array) / 2 - len(destination_fits[1].data) / 2
    crop2 = len(original_array) / 2 + len(destination_fits[1].data) / 2

    return original_array[crop1:crop2, crop1:crop2]


def cartesian_to_polar(x_position, y_position, x_ref_position, y_ref_position):
    x_position, y_position = float(x_position), float(y_position)
    x_ref_position, y_ref_position = float(x_ref_position), float(y_ref_position)

    radius = np.sqrt((x_position - x_ref_position) ** 2 + (y_position - y_ref_position) ** 2)

    if (x_position - x_ref_position) > 0:
        if (y_position - y_ref_position) >= 0:
            angle = np.arctan((y_position - y_ref_position) / (x_position - x_ref_position))
        else:
            angle = 2.0 * np.pi + np.arctan((y_position - y_ref_position) / (x_position - x_ref_position))
    elif (x_position - x_ref_position) < 0:
        angle = np.arctan((y_position - y_ref_position) / (x_position - x_ref_position)) + np.pi
    else:
        if (y_position - y_ref_position) >= 0:
            angle = np.pi / 2
        else:
            angle = -np.pi / 2

    return radius, angle
