"""
    Free Energy
"""

from collections import OrderedDict
import numpy as np
from copy import deepcopy

from pygromos.files.blocks.imd_blocks import PERTURBATION, PRECALCLAM
from pygromos.files.gromos_system import Gromos_System
from pygromos.files.coord.cnf import Cnf
from pygromos.simulations.hpc_queuing.job_scheduling.workers.analysis_workers import simulation_analysis
from pygromos.simulations.hpc_queuing.submission_systems._submission_system import _SubmissionSystem
from pygromos.simulations.hpc_queuing.submission_systems.local import LOCAL
from pygromos.simulations.modules.general_simulation_modules import simulation
from pygromos.utils import bash
from pygromos.utils.typing import List


def TI_sampling(
    in_gromos_system: Gromos_System,
    project_dir: str = None,
    step_name="lambda_sampling",
    lambda_values: List[float] = np.arange(0, 1.1, 0.1),
    subSystem: _SubmissionSystem = LOCAL(),
    n_productions: int = 3,
    n_equilibrations: int = 1,
    randomize: bool = False,
    dual_cnf: List[str] = None,
    verbose: bool = True,
):
    """
    This function will automatically submit N independent (different lambda)
    MD simulations with a lambda dependent potential energy.

    Parameters
    ----------
        in_gromos_system: Gromos_System
            input gromos system
        project_dir: str
            directory in which simulation input files are found
        step_name: str
            subdirectory of project_dir, in which we will write the output
            important: allows to run multiple random seeds with a different "step_name"
        lambda_values: List [float]
            List of lambda values for each independent simulation
        subSystem: _SubmissionSystem
            where will the calculation run
        n_productions: int
            number of chunks each independent simulation is broken down into
        n_equilibrations: int
            number of chunks of equilibration preceding the production for each independent simulation
        randomize: bool
            Choose a random number for the initial random seed (same for all lambda values)
        dual_cnf: List [str], optional
            If provided, should be the path to two conformations (matching end states A and B) which
            can be used as initial conformations
            Simulations with a lambda value between 0 and 0.5 will use the first as starting conformation

    Returns
    --------
        lam_system: Gromos_System
            Gromos system of the simulation submitted last
    """

    if project_dir is None:
        project_dir = in_gromos_system.work_folder

    work_dir = bash.make_folder(project_dir + "/" + step_name)

    # Select a random seed here so all lambda windows have the same
    if randomize:
        in_gromos_system.imd.randomize_seed()

    in_gromos_system.save(work_dir)
    general_system_name_prefix = in_gromos_system.name

    lam_systems = []
    for lam in lambda_values:
        lam = np.round(lam, 4)

        lam_system = deepcopy(in_gromos_system)
        lam_system.name = general_system_name_prefix + "_l_" + str(lam)
        lam_system.work_folder = work_dir

        # Choose different conformations depending on lambda point.
        # e.g. this allows to use RE-EDS SSM conformations.
        # dual_cnf[i] contains the path of the cnf to use
        if dual_cnf is not None and lam <= 0.5:
            lam_system.cnf = Cnf(dual_cnf[0])
        elif dual_cnf is not None:
            lam_system.cnf = Cnf(dual_cnf[1])

        # IMD
        # Pertubation
        # Pertubation of system.
        pert_block = PERTURBATION(NTG=1, NRDGL=0, RLAM=lam, DLAMT=0, ALPHC=0.5, ALPHLJ=0.5, NLAM=2, NSCALE=0)
        lam_system.imd.add_block(block=pert_block)

        # Calculate additional lambda points
        if not hasattr(lam_system, "PRECALCLAM"):
            # Note: This assumes uniformely distributed lambda values
            precalc_lam_block = PRECALCLAM(NRLAM=len(lambda_values), MINLAM=0, MAXLAM=1)
            lam_system.imd.add_block(block=precalc_lam_block)

        # Submit
        out_gromos_system = _TI_lam_step(
            in_gromos_system=lam_system,
            project_dir=work_dir,
            step_name=lam_system.name,
            submission_system=subSystem,
            in_imd_path=None,
            simulation_runs=n_productions,
            equilibration_runs=n_equilibrations,
            verbose=verbose,
        )

        out_gromos_system.save(out_gromos_system.work_folder + "/sd_out_system.obj")
        lam_systems.append(out_gromos_system)

    return lam_system


def _TI_lam_step(
    in_gromos_system: Gromos_System,
    project_dir: str,
    step_name: str = "lam",
    in_imd_path=None,
    submission_system: _SubmissionSystem = LOCAL(),
    simulation_runs: int = 1,
    equilibration_runs: int = 0,
    previous_simulation_run: int = None,
    analysis_script: callable = simulation_analysis.do,
    verbose: bool = True,
) -> Gromos_System:
    template_control_dict = OrderedDict(
        {
            "concat": {
                "do": True,
                "sub": {
                    "cp_cnf": True,
                    "repair_NaN": True,
                    "cp_omd": True,
                    "cat_trc": True,
                    "cat_tre": True,
                    "cat_trg": True,
                    "convert_trcs": False,
                },
            },
            "simulation_folder": {"do": False, "sub": {"tar": False, "remove": False}},
        }
    )
    return simulation(
        in_gromos_simulation_system=in_gromos_system,
        override_project_dir=project_dir,
        previous_simulation_run=previous_simulation_run,
        step_name=step_name,
        in_imd_path=in_imd_path,
        submission_system=submission_system,
        simulation_runs=simulation_runs,
        equilibration_runs=equilibration_runs,
        analysis_control_dict=template_control_dict,
        analysis_script=analysis_script,
        verbose=verbose,
    )
