import tempfile
from pygromos.files.topology import top
from pygromos.tests.test_files.general_file_functions import general_file_tests


from pygromos.tests.in_testfiles import in_test_file_path
from pygromos.tests.test_files import out_test_root_dir

root_in = in_test_file_path + "/top"
in_path = root_in + "/test.top"
simple_path = root_in + "/simpleTest.top"


root_out = tempfile.mkdtemp(dir=out_test_root_dir, prefix="top_")
out_path = root_out + "/out_imd_REEDS1.imd"


class test_top(general_file_tests):
    __test__ = True
    class_type = top.Top
    in_file_path = in_path
    root_out = root_out

    def test_parsing_test_file(self):
        imd_file = self.class_type(self.in_file_path)
        assert isinstance(imd_file, self.class_type)
        return 0


class test_top_simple(general_file_tests):
    __test__ = True
    class_type = top.Top
    in_file_path = simple_path
    root_out = root_out

    def test_parsing_test_file(self):
        top_file = self.class_type(self.in_file_path)
        assert isinstance(top_file, self.class_type)
        return 0

    def test_additon(self):
        top_file = self.class_type(self.in_file_path)
        top_file += top_file
        top_file += top_file
        return 0

    def test_multiplication(self):
        top_file = self.class_type(self.in_file_path)
        top_file *= 3
        return 0

    def test_add_eq_mul(self):
        top_file = self.class_type(self.in_file_path)
        top_mul = top_file * 2
        top_add = top_file + top_file
        assert top_mul == top_add
        return 0
