import unittest, os, tempfile

from pygromos.files.topology.disres import Disres
from pygromos.utils import bash

from pygromos.tests.in_testfiles import in_test_file_path
root_in =in_test_file_path+"/disres"
in_path = root_in+"/in_disres.dat"

from pygromos.tests.test_files import out_test_root_dir
root_out = tempfile.mkdtemp(dir=out_test_root_dir, prefix="disres_")
out_path = root_out+"/out_disres.dat"

class test_disres(unittest.TestCase):
    def test_parsing_test_file(self):
        disres = Disres(in_path)
        #print(disres)
        disres.write(out_path)
        return 0
