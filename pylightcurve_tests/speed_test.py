import numpy as np
import timeit
import pyfits as pf
from scipy.optimize import curve_fit

import pylightcurve



# planet
RpRs = 0.120859422
P = 3.52474859
a = 8.76
e = 0.0
i = 86.71
w = 0.0
T = 100.0

# star
ldcoeffs = (+6.08395E-01, -2.06189E-01, +2.62324E-01, -1.33106E-01)

# time
tt = np.arange(99.9, 100.1, 0.0005)




iterations = 10000
benchtime = timeit.timeit("pylightcurve.transit_claret( ldcoeffs, RpRs, P, a, e, i, w, T, tt )",
                          setup="from __main__ import *", number=iterations)
print "\n{:.5}s total time\n".format(benchtime)
print "\n{:.5}ms per lightcurve ( from {} iterations )\n".format(benchtime / iterations * 1000, iterations)

try:
    from line_profiler import LineProfiler

    def do_profile(follow=[]):
        def inner(func):
            def profiled_func(*args, **kwargs):
                try:
                    profiler = LineProfiler()
                    profiler.add_function(func)
                    for f in follow:
                        profiler.add_function(f)
                    profiler.enable_by_count()
                    return func(*args, **kwargs)
                finally:
                    profiler.print_stats()
            return profiled_func
        return inner

except ImportError:
    def do_profile(follow=[]):
        "Helpful if you accidentally leave in production!"
        def inner(func):
            def nothing(*args, **kwargs):
                return func(*args, **kwargs)
            return nothing
        return inner


@do_profile(follow=[pylightcurve.flux_claret.flux_drop])
def expensive_function():
    for i in range(1000):
        pylightcurve.transit_claret(ldcoeffs, RpRs, P, a, e, i, w, T, tt)

result = expensive_function()

