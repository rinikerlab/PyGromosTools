import tempfile
from pygromos.files import coord as cnf
from pygromos.tests.test_files.general_file_functions import general_file_tests

from pygromos.tests.in_testfiles import in_test_file_path
from pygromos.tests.test_files import out_test_root_dir

in_file_path = in_test_file_path + "/cnf/in_cnf1.cnf"
in_file_renum_path = in_test_file_path + "/cnf/in_renumber_ligs.cnf"


root_out = tempfile.mkdtemp(dir=out_test_root_dir, prefix="cnf_")
out_path = root_out + "/out_cnf1.cnf"
out_red_path = root_out + "/out_cnf1_reduced.cnf"


class test_cnf(general_file_tests):
    __test__ = True

    in_file_path = in_test_file_path + "/cnf/in_cnf1.cnf"
    root_out = root_out
    class_type = cnf.Cnf

    def test_parse(self):
        cnf_file = self.class_type(self.in_file_path)
        self.assertEqual(type(cnf_file), cnf.Cnf)

    def test_delete_res(self):
        cnf_file = cnf.Cnf(in_file_path)
        cnf_file.delete_residue(resID=4, resName="4E96")
        cnf_file.write(out_red_path)

    def test_clean_posiResNumsByName(self):
        cnf_file = cnf.Cnf(in_file_renum_path, clean_resiNumbers_by_Name=True)
        cnf_file.clean_posiResNums()
        residues = cnf_file.get_residues()

        # print(residues)
        self.assertEqual(2, len(residues), "Found  more/less than two residues!")
        self.assertEqual(2, list(residues["3U5L"].keys())[0], "Could not recognzize expected resNum for 3U5l = 2")
        self.assertEqual(1, list(residues["3MXF"].keys())[0], "Could not recognzize expected resNum for 3MXF = 1")
        self.assertEqual(31, residues["3U5L"][2], "Could not recognzize expected atom ammount for 3U5l = 2")
        self.assertEqual(35, residues["3MXF"][1], "Could not recognzize expected atom ammount for 3MXF = 1")

    def test_write_out(self):
        cnf_file = cnf.Cnf(in_file_path)
        cnf_file.write(out_path)

    def test_visualize(self):
        cnf_file = cnf.Cnf(in_file_path)
        view = cnf_file.view
        assert view is not None
