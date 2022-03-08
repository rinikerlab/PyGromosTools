import tempfile

from pygromos.files.topology.disres import Disres
from pygromos.tests.test_files.general_file_functions import general_file_tests

from pygromos.tests.in_testfiles import in_test_file_path
from pygromos.tests.test_files import out_test_root_dir

root_in = in_test_file_path + "/disres"
in_path = root_in + "/in_disres.dat"

root_out = tempfile.mkdtemp(dir=out_test_root_dir, prefix="disres_")
out_path = root_out + "/out_disres.dat"


class test_disres(general_file_tests):
    __test__ = True
    class_type = Disres
    in_file_path = in_path
    root_out = root_out

    def test_parsing_test_file(self):
        disres = self.class_type(self.in_file_path)
        disres.write(out_path)
        return 0
