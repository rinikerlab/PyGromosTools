import unittest
from pygromos.gromos.gromosPP import GromosPP


class test_gromosPP_class(unittest.TestCase):
    def test_construct_with_None(self):
        gromPP = GromosPP(gromosPP_bin_dir=None)
        assert isinstance(gromPP, GromosPP)
