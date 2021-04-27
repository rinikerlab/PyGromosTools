import os
import tempfile
import unittest

from pygromos.files.trajectory import _general_trajectory as gt
from pygromos.files.trajectory import trc, tre, trg
from pygromos.files.trajectory.tre_field_libs import ene_fields

from pygromos.tests.in_testfiles import in_test_file_path
from pygromos.tests.test_files import out_test_root_dir
root_out = tempfile.mkdtemp(dir=out_test_root_dir, prefix="trajs_")

class traj_standard_tests(unittest.TestCase):
    class_name: gt._General_Trajectory = gt._General_Trajectory
    in_file_path = in_file_h5_path = None
    outpath: str

    # Constructors
    def test_constructor_empty(self):
        t = self.class_name(input_value=None)
        print(t)

    def test_constructor_trg_file_path(self):
        t = self.class_name(input_value=self.in_file_path)
        print(t)
        print(t.database.columns)
        print(t.database.head())

    def test_constructor_trg_h5_file_path(self):
        t = self.class_name(input_value=self.in_file_path)
        print(t)

    def test_write(self):
        if (hasattr(self, 't')):
            getattr(self, 't').write(self.outpath)
        else:
            pass

    def setUp(self) -> None:
        self.t1 = self.class_name(input_value=self.in_file_path)

    def test_add(self):
        if("traj_standard_tests" == self.__class__.__name__):
            return 0
        # addition
        tre_3 = self.t1 + self.t1
        tre_3 += self.t1
        print(self.t1, tre_3)
        print(tre_3)


class test_trc(traj_standard_tests):
    class_name = trc.Trc
    in_file_path = in_test_file_path+ "/trc/in_test.trc"
    in_file_h5_path =  in_test_file_path+ "/trc/in_test.trc.h5"
    outpath = root_out + "/out_trc1.trc.h5"


class test_tre(traj_standard_tests):
    class_name = tre.Tre
    in_file_path =  in_test_file_path+ "/tre/in_tre1.tre"
    in_file2_path =  in_test_file_path+ "/tre/in_tre2.tre"
    in_file_eds_path = in_test_file_path+ "/tre/in_eds.tre"
    in_file_lam_path = in_test_file_path+ "/tre/in_lam.tre"
    in_file_h5_path = in_test_file_path+ "/tre/in_tre1.tre.h5"
    outpath = root_out + "/out_tre1.tre.h5"

    def test_get_totals(self):
        t = self.class_name(input_value=self.in_file_path)
        t.tre_block_name_table = ene_fields.gromos_2015_tre_block_names_table
        tots_ene = t.get_totals()
        print(tots_ene)
        pass

    def test_get_eds(self):
        t = self.class_name(input_value=self.in_file_eds_path)
        eds_ene = t.get_eds()
        print(eds_ene.shape)

        self.assertEqual(eds_ene.numstates[0], 9, msg="Number of num_states should be two")
        self.assertEqual(eds_ene.shape, (10,37), msg="The traj should have 500 timesteps and 23 fields for precalclam")
        pass

    def test_get_lam(self):
        t = self.class_name(input_value=self.in_file_lam_path)
        lam_ene = t.get_precalclam()
        print(lam_ene, lam_ene.shape)
        self.assertEqual(lam_ene.nr_lambdas[0], 2, msg="Number of lambdas should be two")
        self.assertEqual(lam_ene.shape, (500,25), msg="The traj should have 500 timesteps and 23 fields for precalclam")


class test_trg(traj_standard_tests):
    class_name = trg.Trg
    in_file_path =  in_test_file_path+ "/trg/test.trg"
    in_file_h5_path =  in_test_file_path+ "/trg/test.trg.h5"
    outpath =  root_out + "/out_tre1.tre.h5"

    # Test Get functions

    def test_get_totals(self):
        totals = self.t1.get_totals()
        print(totals)

    def test_get_lambdas(self):
        lambdas = self.t1.get_lambdas()
        print(lambdas)

    def test_get_precalclam(self):
        t = self.class_name(input_value=self.in_file_path)
        precalclam = t.get_precalclam()
        print(precalclam)
