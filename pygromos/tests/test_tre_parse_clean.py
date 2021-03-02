import unittest

from pygromos.files.energy import tre as ene
#from pygromos.files import _general_gromos_file
import os
path= os.path.dirname(__file__) +"/testfiles/tre/in_tre1.tre"
outpath=os.path.dirname(__file__) +"/testfiles/tre/out_tre1.tre"
class test_disres(unittest.TestCase):

    def test_IO(self):
        tre_file = ene.Tre(path)

        tre_file.clean_timestep_numbering()
        #print(tre_file.get_block_names())
        tre_file.write(outpath)
