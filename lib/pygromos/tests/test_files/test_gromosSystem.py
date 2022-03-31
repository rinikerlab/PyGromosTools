import unittest
import tempfile
from pygromos.utils import bash

from pygromos.files.gromos_system.gromos_system import Gromos_System

from pygromos.tests.in_testfiles import in_test_file_path
from pygromos.data.simulation_parameters_templates import template_md

from pygromos.tests.test_files import out_test_root_dir

tmp_test_dir = tempfile.mkdtemp(dir=out_test_root_dir, prefix="gromSystem_")


class test_gromos_system(unittest.TestCase):
    __test__ = True
    class_type = Gromos_System
    verbose = True

    input_cnf_path = in_test_file_path + "/small_system/6J29.cnf"
    input_top_path = in_test_file_path + "/small_system/6J29.top"

    def test_construct_empty(self):
        subSys = self.class_type(work_folder=tmp_test_dir, system_name="Testing1")
        print(subSys)

    def test_construct_files(self):
        subSys_files = self.class_type(
            work_folder=tmp_test_dir,
            system_name="Testing1",
            in_cnf_path=self.input_cnf_path,
            in_top_path=self.input_top_path,
        )
        assert isinstance(subSys_files, self.class_type)

        subSys_files2 = self.class_type(
            work_folder=tmp_test_dir,
            system_name="Testing2",
            in_cnf_path=self.input_cnf_path,
            in_top_path=self.input_top_path,
            in_imd_path=template_md,
        )

        print(subSys_files2)

    def setUp(self) -> None:
        self.subSys_files = self.class_type(
            work_folder=tmp_test_dir,
            system_name="Testing",
            in_cnf_path=self.input_cnf_path,
            in_top_path=self.input_top_path,
            in_imd_path=template_md,
        )

    def test_rebase(self):
        new_base = bash.make_folder(tmp_test_dir + "/rebase")
        self.subSys_files.work_folder = new_base
        self.subSys_files.rebase_files()

    def test_write(self):
        subSys_files2 = self.class_type(
            work_folder=tmp_test_dir,
            system_name="Testing2",
            in_cnf_path=self.input_cnf_path,
            in_top_path=self.input_top_path,
            in_imd_path=template_md,
        )
        subSys_files2.save(tmp_test_dir + "/out_gromSystem.obj")

    def test_load(self):
        subSys_files2 = self.class_type(
            work_folder=tmp_test_dir,
            system_name="Testing2",
            in_cnf_path=self.input_cnf_path,
            in_top_path=self.input_top_path,
            in_imd_path=template_md,
        )
        load_path = subSys_files2.save(tmp_test_dir + "/out_gromSystem2.obj")

        grom_sys = self.class_type.load(load_path)
        print(grom_sys)
