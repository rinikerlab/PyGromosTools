import unittest, numpy as np
from pygromos.gromos.gromosXX import GromosXX


class test_gromosXX_class(unittest.TestCase):
    def test_construct_with_None(self):
        test_gromosXX_class = GromosXX(gromosXX_bin_dir=None)