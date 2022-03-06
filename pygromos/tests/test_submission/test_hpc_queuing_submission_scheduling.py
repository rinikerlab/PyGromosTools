import unittest
import tempfile

from pygromos.simulations.hpc_queuing.job_scheduling.schedulers import simulation_scheduler

from pygromos.data.simulation_parameters_templates import template_md
from pygromos.data.topology_templates import blank_topo_template
from pygromos.simulations.hpc_queuing.submission_systems import DUMMY
from pygromos.files.gromos_system.gromos_system import Gromos_System

from pygromos.tests.in_testfiles import in_test_file_path
from pygromos.tests.test_files import out_test_root_dir


class test_MD_scheduler(unittest.TestCase):
    submissionSystem = DUMMY

    def setUp(self) -> None:
        self.tmp_test_dir = tempfile.mkdtemp(dir=out_test_root_dir, prefix="scheduling_Dummy_")

    def test_do(self):
        in_cnf = in_test_file_path + "/cnf/in_cnf1.cnf"
        out_dir_path = self.tmp_test_dir

        in_simSystem = Gromos_System(
            system_name="test_do",
            work_folder=out_dir_path,
            in_top_path=blank_topo_template,
            in_cnf_path=in_cnf,
            in_imd_path=template_md,
            in_gromosXX_bin_dir=None,
            in_gromosPP_bin_dir=None,
            verbose=False,
        )
        submission_system = self.submissionSystem()

        simulation_scheduler.do(
            in_simSystem=in_simSystem,
            out_dir_path=out_dir_path,
            submission_system=submission_system,
            simulation_run_num=2,
            verbose=False,
        )
