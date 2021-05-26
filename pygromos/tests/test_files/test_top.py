import unittest, os, tempfile
from pygromos.files.topology import top
from pygromos.tests.test_files.general_file_functions import general_file_tests

from pygromos.utils import bash


from pygromos.tests.in_testfiles import in_test_file_path
root_in = in_test_file_path+"/top"
in_path = root_in+"/test.top"

from pygromos.tests.test_files import out_test_root_dir
root_out = tempfile.mkdtemp(dir=out_test_root_dir, prefix="top_")
out_path = root_out+"/out_imd_REEDS1.imd"


class test_top(general_file_tests):
    __test__ = True
    class_type = top.Top
    in_file_path = in_path
    root_out = root_out

    def test_parsing_test_file(self):
        imd_file = self.class_type(self.in_file_path)
        return 0
