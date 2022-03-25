import tempfile
from pygromos.files.topology import disres
from pygromos.tests.test_files.general_file_functions import general_file_tests

from pygromos.tests.in_testfiles import in_test_file_path
from pygromos.tests.test_files import out_test_root_dir

root_in = in_test_file_path + "/top"
in_path = root_in + "/disres5.disres"

root_out = tempfile.mkdtemp(dir=out_test_root_dir, prefix="top_")
out_path = root_out + "/out_imd_REEDS1.imd"


class test_disres(general_file_tests):
    __test__ = True
    class_type = disres.Disres
    in_file_path = in_path
    root_out = root_out

    def test_parsing_test_file(self):
        imd_file = self.class_type(self.in_file_path)
        assert isinstance(imd_file, self.class_type)
        return 0
