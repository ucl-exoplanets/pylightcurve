
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
