__all__ = ['ldcoeff', 'transit', 'eclipse']

from pylightcurve_tools import *


def transit(ldcoeffs, rprs, p, a, e, i, w, t0, tt, ww=0):
    if np.isnan(w):
        w = 0.
    xyz = position(p, a, e, i * np.pi / 180, w * np.pi / 180, ww * np.pi / 180, t0, tt)
    return single_model(ldcoeffs, rprs, xyz, tt)


def eclipse(fpfs, rprs, p, a, e, i, w, t0, tt, ww=0):
    if np.isnan(w):
        w = 0.
    xyz = position(p, -a / rprs, e, i * np.pi / 180, w * np.pi / 180, ww * np.pi / 180, t0, tt)
    return (1.0 + fpfs * single_model((0, 0, 0, 0), 1 / rprs, xyz, tt)) / (1.0 + fpfs)

