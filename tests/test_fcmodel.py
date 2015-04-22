import unittest
from .. import fcmodel


class Test_ldcoeff(unittest.TestCase):

    def test_ldcoeff_raises_PyLCFilterError_on_invalid_filter(self):

        with self.assertRaises(fcmodel.PyLCFilterError):
            fcmodel.ldcoeff(0.01, 6590, 4.1, 'a')