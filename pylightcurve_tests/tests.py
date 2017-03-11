import numpy as np
import timeit

import pylightcurve as plc


# comparison with some analytically calculated values

print 1 - plc.transit_flux_drop('claret',
                                (0.6023, -0.5110, 0.4655, -0.1752), 0.15, np.array([0.4]))[0] - 0.023672010166127
print 1 - plc.transit_flux_drop('claret',
                                (0.6023, -0.5110, 0.4655, -0.1752), 0.15, np.array([0.15]))[0] - 0.023930512384295
print 1 - plc.transit_flux_drop('claret',
                                (0.6023, -0.5110, 0.4655, -0.1752), 0.15, np.array([0.0]))[0] - 0.023969120193732
print 1 - plc.transit_flux_drop('claret',
                                (0.6023, -0.5110, 0.4655, -0.1752), 0.15, np.array([0.05]))[0] - 0.023964876981085
print 1 - plc.transit_flux_drop('claret',
                                (0.6023, -0.5110, 0.4655, -0.1752), 0.01, np.array([0.0]))[0] - 0.0001066124712015837
print 1 - plc.transit_flux_drop('claret',
                                (0.6023, -0.5110, 0.4655, -0.1752), 0.5, np.array([0.5]))[0] - 0.255855066194966
print 1 - plc.transit_flux_drop('claret',
                                (0.6023, -0.5110, 0.4655, -0.1752), 0.5, np.array([0.4]))[0] - 0.259457198918689
print 1 - plc.transit_flux_drop('claret',
                                (0.0, 0.0, 0.0, 0.9), 0.5, np.array([0.4]))[0] - 0.337954545454545
print 1 - plc.transit_flux_drop('claret',
                                (0.0, 0.0, 0.0, 0.9), 0.5, np.array([0.6]))[0] - 0.255363492624757


# planet
RpRs = 0.120859422
P = 3.52474859
a = 8.76
e = 0.0
i = 86.71
w = 0.0
T = 100.0

# star
method='claret'
ldcoeffs = (+6.08395E-01, -2.06189E-01, +2.62324E-01, -1.33106E-01)

# time
tt = np.arange(99.9, 100.1, 0.0005)




iterations = 10000
benchtime = timeit.timeit("plc.transit(method, ldcoeffs, RpRs, P, a, e, i, w, T, tt )",
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


@do_profile(follow=[plc.transit_flux_drop])
def expensive_function():
    for i in range(1000):
        plc.transit(method, ldcoeffs, RpRs, P, a, e, i, w, T, tt)

result = expensive_function()

