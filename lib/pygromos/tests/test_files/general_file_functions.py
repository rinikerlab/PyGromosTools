import pickle
import copy
import unittest

from pygromos.files.topology.ptp import Pertubation_topology


class general_file_tests(unittest.TestCase):
    __test__ = False

    def test_pickle(self):
        obj_file = self.class_type(self.in_file_path)

        out_path = self.root_out + "/pickleTest.obj"
        pickle.dump(obj=obj_file, file=open(out_path, "wb"))

        obj_loaded = pickle.load(open(out_path, "rb"))

        print(obj_loaded)

    def test_copy(self):
        obj_file = self.class_type(self.in_file_path)
        obj_copy = copy.copy(obj_file)

        print(obj_copy)

    def test_deepcopy(self):
        obj_file = self.class_type(self.in_file_path)
        obj_copy = copy.deepcopy(obj_file)

        print(obj_copy)

    def test_equality(self):
        obj_file_1 = self.class_type(self.in_file_path)
        obj_file_2 = self.class_type(self.in_file_path)
        self.assertEqual(obj_file_1, obj_file_2)

    def test_equalAfterCopy(self):
        obj_file = self.class_type(self.in_file_path)
        obj_copy = copy.deepcopy(obj_file)

        if isinstance(obj_file, Pertubation_topology):
            pass
        else:
            self.assertEqual(obj_file, obj_copy)
