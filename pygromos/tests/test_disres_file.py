import unittest, os
from pygromos.files.topology import top

in_path = os.path.dirname(__file__)+"/testfiles/disres/in_disres.dat"
out_path = os.path.dirname(__file__)+"/testfiles/disres/out_disres.dat"

class test_disres(unittest.TestCase):
    def test_parsing_test_file(self):
        disres = top.disres(in_path)
        #print(disres)
        disres.write(out_path)
        return 0
