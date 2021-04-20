from typing import Tuple
from pygromos.files.gromos_system import Gromos_System

from pygromos.data.simulation_parameters_templates import template_emin, template_md, template_sd
from pygromos.simulations.modules.general_simulation_modules import simulation

from pygromos.hpc_queuing.job_scheduling.workers.analysis_workers import simulation_analysis
from pygromos.hpc_queuing.submission_systems.Submission_Systems import LOCAL, _SubmissionSystem

"""
    Simulations
"""
def emin(in_gromos_system: Gromos_System, project_dir: str, step_name: str = "emin", in_imd_path=template_emin,
         submission_system: _SubmissionSystem = LOCAL(), simulation_runs: int = 1, equilibration_runs: int = 0,
         previous_simulation_run: int = None, analysis_script: callable = simulation_analysis.do) -> Tuple[Gromos_System, int]:

    template_emin_control_dict = simulation_analysis.template_control_dict
    template_emin_control_dict['concat']['cat_trc'] = False
    template_emin_control_dict['concat']['cat_tre'] = False
    template_emin_control_dict['concat']['cat_trg'] = False


    return simulation(in_gromos_system=in_gromos_system, project_dir=project_dir, previous_simulation_run=previous_simulation_run,
                      step_name=step_name, in_imd_path=in_imd_path, submission_system=submission_system,
                      simulation_runs=simulation_runs, equilibration_runs=equilibration_runs, analysis_control_dict = template_emin_control_dict,
                      analysis_script=analysis_script)

def md(in_gromos_system: Gromos_System, project_dir: str, step_name: str = "md", in_imd_path=template_md,
       submission_system: _SubmissionSystem = LOCAL(), simulation_runs: int = 1, equilibration_runs: int = 0,
       previous_simulation_run: int = None, analysis_script: callable = simulation_analysis.do) -> Tuple[Gromos_System, int]:
    return simulation(in_gromos_system=in_gromos_system, project_dir=project_dir, previous_simulation_run=previous_simulation_run,
                      step_name=step_name, in_imd_path=in_imd_path, submission_system=submission_system,
                      simulation_runs=simulation_runs, equilibration_runs=equilibration_runs,
                      analysis_script=analysis_script)

def sd(in_gromos_system: Gromos_System, project_dir: str, step_name: str = "sd", in_imd_path=template_sd,
       submission_system: _SubmissionSystem = LOCAL(), simulation_runs: int = 1, equilibration_runs: int = 0,
       previous_simulation_run: int = None, analysis_script: callable = simulation_analysis.do) -> Tuple[Gromos_System, int]:
    return simulation(in_gromos_system=in_gromos_system, project_dir=project_dir, previous_simulation_run=previous_simulation_run,
                      step_name=step_name, in_imd_path=in_imd_path, submission_system=submission_system,
                      simulation_runs=simulation_runs, equilibration_runs=equilibration_runs,
                      analysis_script=analysis_script)
