import tempfile
import unittest
from numpy import testing
from pygromos.files.trajectory import _general_trajectory as gt
from pygromos.files.trajectory import trc, tre, trg
from pygromos.files.coord.cnf import Cnf
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
        t = self.class_name(input_value=None, auto_save=False)
        assert isinstance(t, self.class_name)
        # print(t)

    def test_constructor_trg_file_path(self):
        t = self.class_name(input_value=self.in_file_path, auto_save=False)
        assert isinstance(t, self.class_name)
        # print(t)
        # print(t.database.columns)
        # print(t.database.head())

    def test_constructor_trg_h5_file_path(self):
        t = self.class_name(input_value=self.in_file_path, auto_save=False)
        assert isinstance(t, self.class_name)
        # print(t)

    def test_write(self):
        if hasattr(self, "t"):
            getattr(self, "t").write(self.outpath)
        else:
            pass

    def setUp(self) -> None:
        self.t1 = self.class_name(input_value=self.in_file_path, auto_save=False)

    def test_add(self):
        if "traj_standard_tests" == self.__class__.__name__:
            return 0
        # addition
        tre_3 = self.t1 + self.t1
        tre_3 += self.t1
        # print(self.t1, tre_3)
        # print(tre_3)


class test_trc(unittest.TestCase):
    class_name = trc.Trc
    help_class = Cnf(in_test_file_path + "/trc/in_test.cnf")
    in_file_path = in_test_file_path + "/trc/in_test.trc"
    in_file_path_h5 = in_test_file_path + "/trc/in_test_trc.h5"
    in_file_w_genbox_path = in_test_file_path + "/trc/in_test_genbox.trc"
    in_file_w_genbox_cnf_path = in_test_file_path + "/trc/in_test_genbox.cnf"

    outpath = root_out + "/out_trc1.h5"
    trc_outpath = root_out + "/out_.trc.gz"

    # Constructors
    def test_constructor_empty(self):
        t = self.class_name()
        assert isinstance(t, self.class_name)
        # print(t)

    def test_constructor_trc_file_path(self):
        t = self.class_name(traj_path=self.in_file_path, in_cnf=self.help_class)
        assert isinstance(t, self.class_name)
        # print(t)

    def test_constructor_trc_file_noTop_path(self):
        t = self.class_name(traj_path=self.in_file_path)
        assert isinstance(t, self.class_name)
        # print(t)

    def test_constructor_trc_h5_file_path(self):
        t = self.class_name(traj_path=self.in_file_path_h5)
        assert isinstance(t, self.class_name)
        # print(t)

    def test_write(self):
        t = self.class_name(traj_path=self.in_file_path, in_cnf=self.help_class)
        t.save(self.outpath)

    def test_trc_with_boxes_traj(self):
        c = Cnf(self.in_file_w_genbox_cnf_path)
        t_origin = self.class_name(traj_path=self.in_file_w_genbox_path, in_cnf=self.in_file_w_genbox_cnf_path)

        # CNF was the last frame:
        testing.assert_allclose(actual=t_origin._unitcell_lengths[-1], desired=c.GENBOX.length)

        # these should not be equal! as the box changes in NPT over time
        assert t_origin._unitcell_lengths[0][0] != c.GENBOX.length[0]
        assert t_origin._unitcell_lengths[0][1] != c.GENBOX.length[1]
        assert t_origin._unitcell_lengths[0][2] != c.GENBOX.length[2]

    def test_to_trc_file(self):
        # Read in trc
        t_origin = self.class_name(traj_path=self.in_file_path, in_cnf=self.help_class)

        # take a subset of frames
        t = t_origin[[5, 7, 12]]

        # write trc
        t.write_trc(self.trc_outpath)

        # read in new trc
        t_new = self.class_name(traj_path=self.trc_outpath, in_cnf=self.help_class)

        # test if new trc coordinates have correct shapes
        assert t_new.xyz.shape[0] == 3
        assert t_new.xyz.shape[1] == t_origin.xyz.shape[1]
        assert t_new.xyz.shape[2] == t_origin.xyz.shape[2]

    def test_to_conf(self):
        t = self.class_name(traj_path=self.in_file_path, in_cnf=self.help_class)

        # TEST DUMMY
        conf_60 = t.to_cnf(60)
        conf_60_None = t[60].to_cnf()
        assert conf_60 == conf_60_None

        # TEST with base
        conf_60 = t.to_cnf(60, base_cnf=self.help_class)
        conf_60_None = t[60].to_cnf(base_cnf=self.help_class)
        assert conf_60 == conf_60_None


class test_tre(traj_standard_tests):
    class_name = tre.Tre
    in_file_path = in_test_file_path + "/tre/in_tre1.tre"
    in_file2_path = in_test_file_path + "/tre/in_tre2.tre"
    in_file_eds_path = in_test_file_path + "/tre/in_eds.tre"
    in_file_lam_path = in_test_file_path + "/tre/in_lam.tre"
    in_file_h5_path = in_test_file_path + "/tre/in_tre1.tre.h5"
    outpath = root_out + "/out_tre1.tre.h5"

    def test_get_totals(self):
        t = self.class_name(input_value=self.in_file_path, auto_save=False)
        t.tre_block_name_table = ene_fields.gromos_2015_tre_block_names_table
        tots_ene = t.get_totals()
        print(tots_ene)
        pass

    def test_get_eds(self):
        t = self.class_name(input_value=self.in_file_eds_path, auto_save=False)
        eds_ene = t.get_eds()
        print(eds_ene.shape)

        self.assertEqual(eds_ene.numstates[0], 9, msg="Number of num_states should be two")
        self.assertEqual(eds_ene.shape, (10, 37), msg="The traj should have 500 timesteps and 23 fields for precalclam")
        pass

    def test_get_lam(self):
        t = self.class_name(input_value=self.in_file_lam_path, auto_save=False)
        lam_ene = t.get_precalclam()
        print(lam_ene, lam_ene.shape)
        self.assertEqual(lam_ene.nr_lambdas[0], 2, msg="Number of lambdas should be two")
        self.assertEqual(lam_ene.shape, (5, 25), msg="The traj should have 5 timesteps and 23 fields for precalclam")


class test_trg(traj_standard_tests):
    class_name = trg.Trg
    in_file_path = in_test_file_path + "/trg/test.trg"
    in_file_h5_path = in_test_file_path + "/trg/test.trg.h5"
    outpath = root_out + "/out_tre1.tre.h5"

    # Test Get functions

    def test_get_totals(self):
        totals = self.t1.get_totals()
        print(totals)

    def test_get_lambdas(self):
        lambdas = self.t1.get_lambdas()
        print(lambdas)

    def test_get_precalclam(self):
        t = self.class_name(input_value=self.in_file_path, auto_save=False)
        precalclam = t.get_precalclam()
        print(precalclam)
