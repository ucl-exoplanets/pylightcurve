import unittest
from .. import fcmodel


class Test_ldcoeff(unittest.TestCase):

    def test_ldcoeff_raises_PyLC_FilterError_on_invalid_filter(self):

        with self.assertRaises(fcmodel.PyLC_FilterError):
            fcmodel.ldcoeff(0.01, 6590, 4.1, 'a')