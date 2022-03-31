import tempfile
import unittest
from pygromos.files.otherfiles import repdat

from pygromos.tests.in_testfiles import in_test_file_path
from pygromos.tests.test_files import out_test_root_dir

in_dir = in_test_file_path + "/repdat"
in_path = in_dir + "/in_REEDS_repdat2_short.dat"
in_path2 = in_dir + "/2_ligs_4E96_4J3I_sopt1_2_repdat.dat"


root_out = tempfile.mkdtemp(dir=out_test_root_dir, prefix="repdat_")
out_path = root_out + "/out_REEDS_repdat2_short.dat"
out_path_clean = root_out + "/out_REEDS_repdat2_cleaned_cat.dat"
out_plot1 = root_out + "/out_REEDS_repdat2_transitions.png"
out_plot2 = root_out + "/out_REEDS_repdat2_transitions2.png"


class test_repdat(unittest.TestCase):  # general_file_tests): #Todo: make copy and deepcopy + pickable
    __test__ = True
    class_type = repdat.Repdat
    in_file_path = in_path
    root_out = root_out

    def test_parsing_test_file(self):
        repdat_file = self.class_type(self.in_file_path)
        assert isinstance(repdat_file, self.class_type)
        return 0

    def test_parsing_new_file(self):
        repdat_file = self.class_type(in_path2)
        assert isinstance(repdat_file, self.class_type)

    def test_write_out(self):
        repdat_file = self.class_type(self.in_file_path)
        repdat_file.write(out_path=out_path)

    def test_cleaning_concat(self):
        repdat_file = self.class_type(self.in_file_path)
        tmp_size1 = repdat_file.DATA.shape[0]
        tmp_run1 = int(repdat_file.DATA.run[tmp_size1 - 1])

        repdat_file2 = repdat.Repdat(in_path)

        # concat and clean
        repdat_file.append(repdat_file2)

        tmp_size2 = repdat_file.DATA.shape[0]
        tmp_run2 = int(repdat_file.DATA.run[tmp_size2 - 1])

        repdat_file.write(out_path=out_path_clean)

        # print("\nrepdats:\t size:\t", tmp_size1, "\truns:\t", tmp_run1)
        # print("repdats:\t size:\t", tmp_size2, "\truns:\t", tmp_run2)
        # print(repdat_file.DATA[:][:20])
        # print(repdat_file.DATA[:][-20:])

        # check_concat
        self.assertEqual(2 * tmp_size1, tmp_size2)
        # check_renumbering
        self.assertEqual(2 * tmp_run1, tmp_run2)

    def test_get_transitions(self):
        repdat_file = self.class_type(self.in_file_path)
        transitions = repdat_file.get_replica_traces()
        assert transitions is not None

    def test_test_transitions(self):
        repdat_file = self.class_type(self.in_file_path)
        transitions = repdat_file.get_replica_traces()
        assert transitions is not None
