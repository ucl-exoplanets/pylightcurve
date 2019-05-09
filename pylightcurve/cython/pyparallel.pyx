import numpy as np
cimport numpy as np
import cython


cdef extern from "stdlib.h":
    void *malloc(size_t size)
    void free(void *ptr)

cdef extern from "pyparallel_menu.h":
    cdef double *transit_flux_drop_claret_cython(double *flux, double a1, double a2, double a3, double a4, double rprs, double *z,
                                                 int size)

@cython.boundscheck(False)
@cython.wraparound(False)
def cc_transit_flux_drop_claret(np.ndarray[double, ndim=1, mode="c"] flux, a1, a2, a3, a4, rprs, np.ndarray z, N):

    transit_flux_drop_claret_cython((<double*> flux.data), a1, a2, a3, a4, rprs, (<double *> z.data), N)
    return None
