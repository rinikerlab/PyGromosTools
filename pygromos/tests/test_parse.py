import os, unittest
from pygromos.files._basics import parser
from pygromos.files import coord, imd, topology

path= os.path.dirname(__file__) +"/testfiles/imd/in_imd_REEDS1.imd"
outpath=os.path.dirname(__file__) +"/testfiles/imd/out_imd_REEDS1.imd"

class test_parser(unittest.TestCase):

    def test_imd(self):
        path = os.path.dirname(__file__) + "/testfiles/imd/in_imd_REEDS1.imd"
        outpath = os.path.dirname(__file__) + "/testfiles/imd/out_imd_REEDS1.imd"
        test= parser.read_imd(path)




class test_repdat(unittest.TestCase):
    root_dir = os.path.dirname(__file__)
    in_path = root_dir + "/testfiles/trc/in_test.trc"
    out_file = root_dir + "/testfiles/trc/out.trc"

    def test_constructing_trc_file(self):
        trc = coord.Trc(self.in_path)
       #print(trc)
       #print(vars(trc))
        return 0

    def test_dt_func(self):
        trc = coord.Trc(self.in_path)
        dt=trc.get_dt()
       #print(dt)

    def test_dstep_func(self):
        trc = coord.Trc(self.in_path)
        dstep=trc.get_dsteps()
       #print(dstep)

    def test_clean_timesteps_func(self):
        trc = coord.Trc(self.in_path)
        trc.clean_timestep_numbering()

    def test_write_fnc(self):
        trc = coord.Trc(self.in_path)
        trc.write(self.out_file)
