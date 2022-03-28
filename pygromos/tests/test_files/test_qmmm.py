import tempfile
from pygromos.files.simulation_parameters import imd
from pygromos.files.qmmm import qmmm
from pygromos.tests.test_files.general_file_functions import general_file_tests


from pygromos.tests.in_testfiles import in_test_file_path
from pygromos.tests.test_files import out_test_root_dir

root_in = in_test_file_path + "/qmmm"
in_path_imd = root_in + "/md.imd"
in_path_qmmm = root_in + "/menthol-methanol-dmf.qmmm"

root_out = tempfile.mkdtemp(dir=out_test_root_dir, prefix="qmmm_")
out_path_imd = root_out + "/out_qmmm.imd"
out_path_qmmm = root_out + "/out_qmmm.qmmm"


class test_qmmm_imd(general_file_tests):
    __test__ = True
    class_type = imd.Imd
    in_file_path = in_path_imd
    root_out = root_out

    def test_parsing_test_file(self):
        imd_file = self.class_type(self.in_file_path)
        assert isinstance(imd_file, self.class_type)
        return 0

    def test_to_string(self):
        imd_file = self.class_type(self.in_file_path)
        print(imd_file)
        return 0

    def test_write_out(self):
        imd_file = self.class_type(self.in_file_path)
        imd_file.TITLE.content = "NEW TEST!"
        imd_file.write(out_path_imd)
        return 0


class test_qmmm(general_file_tests):
    __test__ = True
    class_type = qmmm.QMMM
    in_file_path = in_path_qmmm
    root_out = root_out

    def test_parsing_test_file(self):
        qmmm_file = self.class_type(self.in_file_path)
        assert isinstance(qmmm_file, self.class_type)
        return 0

    def test_to_string(self):
        qmmm_file = self.class_type(self.in_file_path)
        print(qmmm_file)
        return 0

    def test_write_out(self):
        qmmm_file = self.class_type(self.in_file_path)
        qmmm_file.TITLE.content = "NEW TEST!"
        qmmm_file.write(out_path_qmmm)
        return 0
