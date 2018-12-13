from __future__ import absolute_import
from __future__ import division
from __future__ import print_function


class PyLCError(BaseException):
    pass


class PyLCLibraryError(PyLCError):
    pass


class PyLCFileError(PyLCError):
    pass


class PyLCProcessError(PyLCError):
    pass


class PyLCInputError(PyLCError):
    pass
