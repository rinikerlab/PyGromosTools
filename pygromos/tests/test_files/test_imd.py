import unittest, os, tempfile
from pygromos.files.simulation_parameters import imd
from pygromos.utils import bash


from pygromos.tests.in_testfiles import in_test_file_path
root_in = in_test_file_path+"/imd"
in_path = root_in+"/in_imd_REEDS1.imd"

from pygromos.tests.test_files import out_test_root_dir
root_out = tempfile.mkdtemp(dir=out_test_root_dir, prefix="imd_")
out_path = root_out+"/out_imd_REEDS1.imd"


class test_imd(unittest.TestCase):
    def test_parsing_test_file(self):
        imd_file = imd.Imd(in_path) #"./testfiles/BRD4_alternat_test_file.imd")
        return 0

    def test_to_string(self):
        imd_file = imd.Imd(in_path)
        #print(imd_file)
        return 0


    def test_edit_REEDS(self):
        imd_file = imd.Imd(in_path)
        svals = "1.0 1.0 1.0 1.0".split()
        EIR = 0.0 # #
        EIR_VECTOR =  [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6]   #write dev cases!
        #EIR_MATRIX = [[0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
         #             [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
         #             [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
         #             [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6]]
        EIR = EIR_VECTOR
        imd_file.edit_REEDS(SVALS=svals, EIR=EIR)
        return 0


    def test_write_out(self):
        imd_file = imd.Imd(in_path)
        imd_file.TITLE.content = "NEW TEST!"
        imd_file.write(out_path)
        return 0
