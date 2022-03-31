import unittest
import numpy as np

import pygromos.simulations.modules.general_simulation_modules
import pygromos.simulations.modules.ti_modules
from pygromos.simulations.modules import preset_simulation_modules
from pygromos.simulations.hpc_queuing.submission_systems import DUMMY

from pygromos.tests.in_testfiles import in_test_file_path
from pygromos.files.gromos_system.gromos_system import Gromos_System

from pygromos.tests.test_simulation_blocks import out_test_root_dir


class test_simulation_blocks(unittest.TestCase):
    verbose = False
    submissionSystem = DUMMY(verbose=False)
    sim_block = preset_simulation_modules.emin
    input_cnf_path = in_test_file_path + "/small_system/6J29.cnf"
    input_top_path = in_test_file_path + "/small_system/6J29.top"

    def setUp(self) -> None:
        self.tmp_test_dir = out_test_root_dir
        self.gromSystem = Gromos_System(
            work_folder=self.tmp_test_dir,
            system_name=str(__name__),
            in_cnf_path=self.input_cnf_path,
            in_top_path=self.input_top_path,
            verbose=self.verbose,
        )

        print(self.gromSystem)

    def test_smulation(self):
        from pygromos.data.simulation_parameters_templates import template_md

        self.gromSystem.imd = template_md
        pygromos.simulations.modules.general_simulation_modules.simulation(
            in_gromos_simulation_system=self.gromSystem,
            override_project_dir=self.tmp_test_dir,
            submission_system=self.submissionSystem,
            verbose=self.verbose,
        )

    def test_emin(self):
        preset_simulation_modules.emin(
            in_gromos_system=self.gromSystem,
            override_project_dir=self.tmp_test_dir,
            submission_system=self.submissionSystem,
            verbose=self.verbose,
        )

    def test_sd(self):
        preset_simulation_modules.sd(
            in_gromos_system=self.gromSystem,
            override_project_dir=self.tmp_test_dir,
            submission_system=self.submissionSystem,
            verbose=self.verbose,
        )

    def test_md(self):
        preset_simulation_modules.md(
            in_gromos_system=self.gromSystem,
            override_project_dir=self.tmp_test_dir,
            submission_system=self.submissionSystem,
            verbose=self.verbose,
        )

    def test_lam_window(self):
        from pygromos.data import simulation_parameters_templates as imd_templates

        pygromos.simulations.modules.ti_modules._TI_lam_step(
            in_gromos_system=self.gromSystem,
            project_dir=self.tmp_test_dir,
            in_imd_path=imd_templates.template_md,
            submission_system=self.submissionSystem,
            verbose=self.verbose,
        )

    def test_ti_sampling(self):
        from pygromos.data.simulation_parameters_templates import template_md

        self.gromSystem.imd = template_md
        pygromos.simulations.modules.ti_modules.TI_sampling(
            in_gromos_system=self.gromSystem,
            project_dir=self.tmp_test_dir,
            lambda_values=np.arange(0, 1.1, 0.1),
            subSystem=self.submissionSystem,
            n_productions=3,
            n_equilibrations=1,
            verbose=self.verbose,
        )
