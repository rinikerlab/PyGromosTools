import unittest, os
from pygromos.files import imd

out_folder = os.path.dirname(__file__)+"/testfiles/gromos_files"
in_path = os.path.dirname(__file__)+"/testfiles/imd/in_imd_REEDS1.imd"
out_path = os.path.dirname(__file__)+"/testfiles/imd/out_imd_REEDS1.imd"

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
