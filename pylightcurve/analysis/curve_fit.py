
__all__ = ['curve_fit']

import warnings

from scipy.optimize import curve_fit as scipy_curve_fit


def curve_fit(*args, **kwargs):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message='Covariance of the parameters could not be estimated')
        return scipy_curve_fit(*args, **kwargs)
