import tempfile
from pygromos.data import pdb_lib
from pygromos.files.otherfiles import residue_library
from pygromos.tests.test_files.general_file_functions import general_file_tests


# from pygromos.tests.in_testfiles import in_test_file_path

from pygromos.tests.test_files import out_test_root_dir

root_out = tempfile.mkdtemp(dir=out_test_root_dir, prefix="resnLib_")
outpath = root_out + "/out_lib.resnlib"


class test_resnlib(general_file_tests):
    __test__ = True
    class_type = residue_library.residue_library
    in_file_path = pdb_lib
    root_out = root_out

    def test_IO(self):
        resnLib = self.class_type(self.in_file_path)
        print(resnLib)
        rb = resnLib.write(outpath)
        print(rb)
