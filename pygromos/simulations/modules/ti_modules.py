"""
    Free Energy
"""

from collections import OrderedDict
from typing import List
import numpy as np

from pygromos.files.blocks.imd_blocks import PERTURBATION, PRECALCLAM
from pygromos.files.gromos_system import Gromos_System
from pygromos.simulations.hpc_queuing.job_scheduling.workers.analysis_workers import simulation_analysis
from pygromos.simulations.hpc_queuing.submission_systems._submission_system import _SubmissionSystem
from pygromos.simulations.hpc_queuing.submission_systems.local import LOCAL
from pygromos.simulations.modules.general_simulation_modules import simulation
from pygromos.utils import bash


def TI_sampling(in_gromos_system: Gromos_System, project_dir: str, step_name="lambda_sampling",
                lambda_values: List[float] = np.arange(0, 1.1, 0.1), subSystem: _SubmissionSystem = LOCAL(),
                n_simulation_repetitions: int = 3, n_equilibrations: int = 1):

    work_dir = bash.make_folder(project_dir + "/" + step_name)
    in_gromos_system.save(work_dir)
    sys_path = in_gromos_system.save(work_dir+"/"+step_name+"_input.obj") #Ugly workaround for deepcopy
    general_system_name_prefix = in_gromos_system.name

    lam_systems = []
    for lam in lambda_values:
        lam = np.round(lam, 2)

        lam_system = Gromos_System.load(sys_path)
        lam_system.name = general_system_name_prefix+"_l_" + str(lam)
        lam_system.work_folder = work_dir

        # IMD
        ## Pertubation
        ###Pertubation of system.
        pert_block = PERTURBATION(NTG=1, NRDGL=0, RLAM=lam, DLAMT=0,
                                  ALPHC=0.5, ALPHLJ=0.5, NLAM=2, NSCALE=0)
        lam_system.imd.add_block(block=pert_block)

        ###Calculate additional lambda points for end-states
        if(not hasattr(lam_system, "PRECALCLAM")):
            precalc_lam_block = PRECALCLAM(NRLAM=2, MINLAM=0, MAXLAM=1)
            lam_system.imd.add_block(block=precalc_lam_block)

        # Submit
        out_gromos_system, jobID = _TI_lam_step(in_gromos_system=lam_system, project_dir=work_dir,
                                                step_name=lam_system.name, submission_system=subSystem,
                                                in_imd_path=None,
                                                simulation_runs=n_simulation_repetitions,
                                                equilibration_runs=n_equilibrations)

        out_gromos_system.save(out_gromos_system.work_folder + "/sd_out_system.obj")
        lam_systems.append(out_gromos_system)

    return lam_system, jobID


def _TI_lam_step(in_gromos_system: Gromos_System, project_dir: str, step_name: str = "lam", in_imd_path=None,
                 submission_system: _SubmissionSystem = LOCAL(), simulation_runs: int = 1, equilibration_runs: int = 0,
                 previous_simulation_run: int = None, analysis_script: callable = simulation_analysis.do) ->(Gromos_System, int):
    template_control_dict = OrderedDict({
        "concat": {"do": True,
                   "sub": {
                       "cp_cnf": True,
                       "repair_NaN": True,
                       "cp_omd": True,
                       "cat_trc": True,
                       "cat_tre": True,
                       "cat_trg": True,
                       "convert_trcs": False,
                   }
                   },
        "simulation_folder": {
            "do": False,
            "sub": {
                "tar": False,
                "remove": False
            }
        }
    })
    return simulation(in_gromos_simulation_system=in_gromos_system, override_project_dir=project_dir, previous_simulation_run=previous_simulation_run,
                      step_name=step_name, in_imd_path=in_imd_path, submission_system=submission_system,
                      simulation_runs=simulation_runs, equilibration_runs=equilibration_runs, analysis_control_dict=template_control_dict,
                      analysis_script=analysis_script)