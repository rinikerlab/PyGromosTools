import numpy as np
from typing import Tuple
from pygromos.files.gromos_system import Gromos_System

from pygromos.data.simulation_parameters_templates import template_emin, template_md, template_sd
from pygromos.simulations.modules.general_simulation_modules import simulation

from pygromos.simulations.hpc_queuing.job_scheduling.workers.analysis_workers import simulation_analysis
from pygromos.simulations.hpc_queuing.submission_systems._submission_system import _SubmissionSystem
from pygromos.simulations.hpc_queuing.submission_systems.local import LOCAL

"""
    Simulations
"""


def emin(
    in_gromos_system: Gromos_System,
    step_name: str = "emin",
    override_project_dir: str = None,
    in_imd_path=None,
    submission_system: _SubmissionSystem = LOCAL(),
    simulation_runs: int = 1,
    equilibration_runs: int = 0,
    previous_simulation_run: int = None,
    _template_imd_path: str = template_emin,
    initialize_first_run=False,
    analysis_script: callable = simulation_analysis.do,
    verbose: bool = True,
) -> Tuple[Gromos_System, int]:
    template_emin_control_dict = simulation_analysis.template_control_dict
    template_emin_control_dict["concat"]["cat_trc"] = False
    template_emin_control_dict["concat"]["cat_tre"] = False
    template_emin_control_dict["concat"]["cat_trg"] = False

    if hasattr(in_gromos_system.imd, "WRITETRAJ"):
        if in_gromos_system.imd.WRITETRAJ.NTWX > 0:
            template_emin_control_dict["concat"]["cat_trc"] = False
        if in_gromos_system.imd.WRITETRAJ.NTWE > 0:
            template_emin_control_dict["concat"]["cat_tre"] = False
        if in_gromos_system.imd.WRITETRAJ.NTWG > 0:
            template_emin_control_dict["concat"]["cat_trg"] = False

    return simulation(
        in_gromos_simulation_system=in_gromos_system,
        override_project_dir=override_project_dir,
        previous_simulation_run=previous_simulation_run,
        step_name=step_name,
        in_imd_path=in_imd_path,
        submission_system=submission_system,
        initialize_first_run=initialize_first_run,
        simulation_runs=simulation_runs,
        equilibration_runs=equilibration_runs,
        analysis_control_dict=template_emin_control_dict,
        analysis_script=analysis_script,
        _template_imd_path=_template_imd_path,
        verbose=verbose,
    )


def md(
    in_gromos_system: Gromos_System,
    step_name: str = "md",
    override_project_dir: str = None,
    in_imd_path=None,
    submission_system: _SubmissionSystem = LOCAL(),
    simulation_runs: int = 1,
    equilibration_runs: int = 0,
    initialize_first_run=False,
    reinitialize_every_run=False,
    previous_simulation_run: int = None,
    _template_imd_path: str = template_md,
    analysis_script: callable = simulation_analysis.do,
    verbose: bool = True,
) -> Tuple[Gromos_System, int]:
    return simulation(
        in_gromos_simulation_system=in_gromos_system,
        override_project_dir=override_project_dir,
        previous_simulation_run=previous_simulation_run,
        step_name=step_name,
        in_imd_path=in_imd_path,
        submission_system=submission_system,
        initialize_first_run=initialize_first_run,
        simulation_runs=simulation_runs,
        equilibration_runs=equilibration_runs,
        analysis_script=analysis_script,
        _template_imd_path=_template_imd_path,
        verbose=verbose,
    )


def sd(
    in_gromos_system: Gromos_System,
    step_name: str = "sd",
    override_project_dir: str = None,
    in_imd_path=None,
    submission_system: _SubmissionSystem = LOCAL(),
    simulation_runs: int = 1,
    equilibration_runs: int = 0,
    initialize_first_run=False,
    reinitialize_every_run=False,
    previous_simulation_run: int = None,
    _template_imd_path: str = template_sd,
    analysis_script: callable = simulation_analysis.do,
    verbose: bool = True,
) -> Tuple[Gromos_System, int]:
    return simulation(
        in_gromos_simulation_system=in_gromos_system,
        override_project_dir=override_project_dir,
        previous_simulation_run=previous_simulation_run,
        step_name=step_name,
        in_imd_path=in_imd_path,
        submission_system=submission_system,
        initialize_first_run=initialize_first_run,
        simulation_runs=simulation_runs,
        equilibration_runs=equilibration_runs,
        analysis_script=analysis_script,
        _template_imd_path=_template_imd_path,
        verbose=verbose,
    )


def thermalisation(
    in_gromos_system: Gromos_System,
    temperatures=np.linspace(60, 300, 4),
    step_name: str = "eq_therm",
    override_project_dir: str = None,
    in_imd_path=None,
    submission_system: _SubmissionSystem = LOCAL(),
    simulation_runs: int = 1,
    equilibration_runs: int = 0,
    previous_simulation_run: int = None,
    _template_imd_path: str = template_sd,
    analysis_script: callable = simulation_analysis.do,
    verbose: bool = True,
) -> Tuple[Gromos_System, int]:

    for runID, temperature in enumerate(temperatures):
        print("run", runID, "T: ", temperature)

        # adapt temperature
        in_gromos_system.imd.MULTIBATH.TEMP0 = [temperature for x in range(in_gromos_system.imd.MULTIBATH.NBATHS)]

        # turn off the posres for the last run.
        if runID + 1 == len(temperatures):
            in_gromos_system.imd.POSITIONRES.NTPOR = 0
            in_gromos_system.imd.POSITIONRES.CPOR = 0

            # Last run
            return simulation(
                in_gromos_simulation_system=in_gromos_system,
                override_project_dir=override_project_dir,
                previous_simulation_run=previous_simulation_run,
                step_name=step_name,
                in_imd_path=in_imd_path,
                submission_system=submission_system,
                simulation_runs=simulation_runs,
                equilibration_runs=equilibration_runs,
                analysis_script=analysis_script,
                _template_imd_path=_template_imd_path,
                verbose=verbose,
            )

        else:
            simulation(
                in_gromos_simulation_system=in_gromos_system,
                override_project_dir=override_project_dir,
                previous_simulation_run=previous_simulation_run,
                step_name=step_name,
                in_imd_path=in_imd_path,
                submission_system=submission_system,
                simulation_runs=simulation_runs,
                equilibration_runs=equilibration_runs,
                analysis_script=analysis_script,
                _template_imd_path=_template_imd_path,
                verbose=verbose,
            )
