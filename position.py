__all__ = ['position_vector']

import numpy as np

pi = np.pi


class PYLCError(BaseException):
    pass


class PYLCOptimiseError(PYLCError):
    pass
    

def position_vector(period, sma_over_rs, eccentricity, inclination, periastron, mid_time, time_array, ww=0):
    if np.isnan(periastron):
        periastron = 0.
    inclination = inclination * pi / 180.0
    periastron = periastron * pi / 180.0
    ww = ww * pi / 180.0

    if eccentricity == 0 and ww == 0:
        vv = 2 * pi * (time_array - mid_time) / period
        bb = sma_over_rs * np.cos(vv)
        return [bb * np.sin(inclination), sma_over_rs * np.sin(vv), - bb * np.cos(inclination)]

    if periastron < pi / 2:
        aa = 1.0 * pi / 2 - periastron
    else:
        aa = 5.0 * pi / 2 - periastron
    bb = 2 * np.arctan(np.sqrt((1 - eccentricity) / (1 + eccentricity)) * np.tan(aa / 2))
    if bb < 0:
        bb += 2 * pi
    mid_time = float(mid_time) - (period / 2.0 / pi) * (bb - eccentricity * np.sin(bb))
    m = (time_array - mid_time - np.int_((time_array - mid_time) / period) * period) * 2.0 * pi / period
    u0 = m
    stop = False
    u1 = 0
    for ii in xrange(10000):  # setting a limit of 1k iterations - arbitrary limit
        u1 = u0 - (u0 - eccentricity * np.sin(u0) - m) / (1 - eccentricity * np.cos(u0))
        stop = (np.abs(u1 - u0) < 10 ** (-7)).all()
        if stop:
            break
        else:
            u0 = u1
    if not stop:
        raise PYLCOptimiseError("Failed to find a solution in 10000 loops")
    vv = 2 * np.arctan(np.sqrt((1 + eccentricity) / (1 - eccentricity)) * np.tan(u1 / 2))
    #
    rr = sma_over_rs * (1 - (eccentricity ** 2)) / (np.ones_like(vv) + eccentricity * np.cos(vv))
    aa = np.cos(vv + periastron)
    bb = np.sin(vv + periastron)
    x = rr * bb * np.sin(inclination)
    y = rr * (-aa * np.cos(ww) + bb * np.sin(ww) * np.cos(inclination))
    z = rr * (-aa * np.sin(ww) - bb * np.cos(ww) * np.cos(inclination))
    return [x, y, z]


def projected_distance(period, sma_over_rs, eccentricity, inclination, periastron, mid_time, time_array, ww=0):
    xv, yv, zv = position_vector(
        period, sma_over_rs, eccentricity, inclination, periastron, mid_time, time_array, ww=ww)
    return np.sqrt(yv * yv + zv * zv)
