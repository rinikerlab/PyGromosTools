import tempfile
from pygromos.files.simulation_parameters import imd
from pygromos.tests.test_files.general_file_functions import general_file_tests


from pygromos.tests.in_testfiles import in_test_file_path
from pygromos.tests.test_files import out_test_root_dir

root_in = in_test_file_path + "/imd"
in_path = root_in + "/in_REEDS1.imd"

root_out = tempfile.mkdtemp(dir=out_test_root_dir, prefix="imd_")
out_path = root_out + "/out_imd_REEDS1.imd"


class test_imd(general_file_tests):
    __test__ = True
    class_type = imd.Imd
    in_file_path = in_path
    root_out = root_out

    def test_parsing_test_file(self):
        imd_file = self.class_type(self.in_file_path)
        assert isinstance(imd_file, self.class_type)
        return 0

    def test_to_string(self):
        imd_file = self.class_type(self.in_file_path)
        print(imd_file)
        return 0

    def test_edit_REEDS(self):
        imd_file = self.class_type(self.in_file_path)
        svals = "1.0 1.0 1.0 1.0".split()
        EIR = 0.0  # #
        EIR_VECTOR = [
            0.0,
            0.1,
            0.2,
            0.3,
            0.4,
        ]  # write dev cases!
        EIR = EIR_VECTOR
        imd_file.edit_REEDS(SVALS=svals, EIR=EIR)
        return 0

    def test_write_out(self):
        imd_file = self.class_type(self.in_file_path)
        imd_file.TITLE.content = "NEW TEST!"
        imd_file.write(out_path)
        return 0
