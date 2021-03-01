import unittest, os, tempfile, numpy as np

from pygromos.simulation_runner import simulation_building_blocks
from pygromos.hpc_queuing.submission_systems.Submission_Systems import DUMMY

from pygromos.tests.in_testfiles import in_test_file_path
from pygromos.files.gromos_system.gromos_system import Gromos_System

from pygromos.tests.test_simulation_blocks import out_test_root_dir

class test_simulation_blocks(unittest.TestCase):
    submissionSystem = DUMMY()
    sim_block = simulation_building_blocks.emin
    input_cnf_path = in_test_file_path+"/small_system/6J29.cnf"
    input_top_path = in_test_file_path+"/small_system/6J29.top"

    def setUp(self) -> None:
        self.tmp_test_dir = out_test_root_dir
        self.gromSystem = Gromos_System(work_folder=self.tmp_test_dir, system_name=str(__name__),
                                        in_cnf_path=self.input_cnf_path, in_top_path=self.input_top_path)

        print(self.gromSystem)

    def test_smulation(self):
        from pygromos.data.simulation_parameters_templates import template_md
        self.gromSystem.imd = template_md
        simulation_building_blocks.simulation(in_gromos_system=self.gromSystem, project_dir=self.tmp_test_dir,
                      submission_system=self.submissionSystem)

    def test_emin(self):
        simulation_building_blocks.emin(in_gromos_system=self.gromSystem, project_dir=self.tmp_test_dir,
                      submission_system=self.submissionSystem)

    def test_sd(self):
        simulation_building_blocks.sd(in_gromos_system=self.gromSystem, project_dir=self.tmp_test_dir,
                      submission_system=self.submissionSystem)

    def test_md(self):
        simulation_building_blocks.md(in_gromos_system=self.gromSystem, project_dir=self.tmp_test_dir,
                      submission_system=self.submissionSystem)

    def test_lam_window(self):
        simulation_building_blocks._TI_lam_step(in_gromos_system=self.gromSystem, project_dir=self.tmp_test_dir,
                      submission_system=self.submissionSystem)

    def test_ti_sampling(self):
        from pygromos.data.simulation_parameters_templates import template_md
        self.gromSystem.imd = template_md
        simulation_building_blocks.TI_sampling(in_gromos_system=self.gromSystem, project_dir=self.tmp_test_dir,
                    lambda_values = np.arange(0, 1.1, 0.1), subSystem = self.submissionSystem,
                    n_simulation_repetitions = 3, n_equilibrations = 1)