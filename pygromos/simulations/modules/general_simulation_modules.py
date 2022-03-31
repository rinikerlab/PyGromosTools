import os
import sys
import glob
import traceback
import time
import warnings
from copy import deepcopy
from collections import OrderedDict

from pygromos.files.coord import cnf
from pygromos.files.gromos_system import Gromos_System
from pygromos.files.trajectory.trc import Trc
from pygromos.files.trajectory.tre import Tre
from pygromos.files.trajectory.trg import Trg
from pygromos.simulations.hpc_queuing.job_scheduling.schedulers import simulation_scheduler
from pygromos.simulations.hpc_queuing.job_scheduling.workers.analysis_workers import simulation_analysis
from pygromos.simulations.hpc_queuing.submission_systems._submission_system import _SubmissionSystem
from pygromos.simulations.hpc_queuing.submission_systems.local import LOCAL
from pygromos.utils import bash, utils
from pygromos.utils.utils import spacer as spacer, time_wait_s_for_filesystem


def simulation(
    in_gromos_simulation_system: Gromos_System,
    override_project_dir: str = None,
    step_name: str = "sim",
    in_imd_path: str = None,
    submission_system: _SubmissionSystem = LOCAL(),
    simulation_runs: int = 1,
    equilibration_runs: int = 0,
    previous_simulation_run: int = None,
    force_simulation: bool = False,
    initialize_first_run=False,
    reinitialize_every_run=False,
    analysis_script: callable = simulation_analysis.do,
    analysis_control_dict: dict = None,
    _no_double_submit_check: bool = False,
    _work_dir: str = None,
    verbose: bool = True,
    verbose_lvl: int = 1,
    _template_imd_path: str = None,
) -> Gromos_System:
    """
        This function is a generic simulation block, that can be used to run and schedule simulations.

    Parameters
    ----------
    in_gromos_simulation_system : Gromos_System
        gromos system that contains the information of the files, that are required for the simulation run.
    override_project_dir : str
        parent project directory
    step_name :  str, optional
        name of the step - used as jobname and folder in the project dir.
    in_imd_path : str, optional
        if a path to an imd file is given here, it will be used!
    submission_system : _SubmissionSystem, optional
        the system, that should be used to submit a job
    simulation_runs : int, optional
        number of sequential simulation runs
    equilibration_runs : int, optional
        number of equilibration runs, that will not appear in the final traj
    previous_simulation_run : int, optional
        job-ID of the previous simulation
    force_simulation : bool, optional
        if simulation already exist, shall it be overwritten?
    analysis_script : callable, optional
        script that is used to analyse the job
    analysis_control_dict : dict, optional
        sub selection of steps, for the analysis script, that should be executed
    verbose : bool, optional
        "baeh, baeh, baeh"  - a wise black-nose sheep from ausserberg.
    verbose_lvl : int, optional
        how much baeh?
    _template_imd_path : str, optional
        default template imd

    Returns
    -------
    Gromos_System,
        returns a new gromos system which is containing the simulation info
    """
    # PREPERATIONS
    try:
        try:
            gromos_system = deepcopy(in_gromos_simulation_system)
            if previous_simulation_run is not None:
                gromos_system._last_jobID = previous_simulation_run

            # check if override dir is given and set project to correct location
            if override_project_dir is not None:
                init_work_folder = override_project_dir
                step_dir = override_project_dir + "/" + step_name
            else:
                init_work_folder = gromos_system.work_folder
                step_dir = gromos_system.work_folder + "/" + step_name

            bash.make_folder(step_dir)
            out_input_dir = step_dir + "/input"
            out_simulation_dir = step_dir + "/simulation"
            out_analysis_dir = step_dir + "/analysis"
            bash.make_folder(out_input_dir)
            os.chdir(step_dir)

            # Prepare gromos system:
            gromos_system.work_folder = out_input_dir
            gromos_system.name = step_name

            if in_imd_path is not None:
                gromos_system.imd = in_imd_path
            elif hasattr(gromos_system.imd, "TITLE"):
                pass
            elif _template_imd_path is not None:
                if verbose:
                    warnings.warn("Template_imd_path was used: " + _template_imd_path)
                gromos_system.imd = _template_imd_path
                gromos_system.prepare_for_simulation()
            else:
                raise ValueError(
                    "Could not find any .imd path (gromos system has no imd, in_imd_path not given and also no _template_imd_path!)"
                )

            out_analysis_cnf = out_analysis_dir + "/data/" + gromos_system.name + ".cnf"

            if verbose:
                print(spacer)
                print(step_name)
                print(spacer)

            # Write out, all non promised files
            gromos_system.rebase_files()

            # Write Out Ana Script
            in_analysis_control_dict = analysis_control_dict
            n_analysis_processors = 1  # Maybe in the future: 5 if(nmpi>5) else 1
        except Exception as err:
            raise Exception("Could not prepare the gromos System\n\t" + "\n\t".join(map(str, err.args)))
        # do
        if analysis_script is not None:
            analysis_control_dict = (
                simulation_analysis.template_control_dict
                if (in_analysis_control_dict is None)
                else in_analysis_control_dict
            )

            analysis_vars = OrderedDict(
                {
                    "in_simulation_dir": out_simulation_dir,
                    "sim_prefix": gromos_system.name,
                    "out_analysis_dir": out_analysis_dir,
                    "gromosPP_bin_dir": gromos_system.gromosPP._bin,
                    "control_dict": analysis_control_dict,
                    "n_processes": n_analysis_processors,
                    "verbose": verbose,
                    "verbose_lvl": verbose_lvl,
                }
            )
            try:
                in_analysis_script_path = utils.write_job_script(
                    out_script_path=step_dir + "/job_analysis.py",
                    target_function=analysis_script,
                    variable_dict=analysis_vars,
                )
            except Exception as err:
                raise Exception("Could not prepare the analysis script\n\t" + "\n\t".join(map(str, err.args)))
        else:
            in_analysis_script_path = None

        # Write Out schedulling Script
        # Build analysis_script
        MD_job_vars = OrderedDict(
            {
                "in_simSystem": gromos_system,
                "out_dir_path": out_simulation_dir,
                "simulation_run_num": simulation_runs,
                "equilibration_run_num": equilibration_runs,
                "submission_system": submission_system,
                "analysis_script_path": in_analysis_script_path,
                "initialize_first_run": initialize_first_run,
                "reinitialize_every_run": reinitialize_every_run,
                "previous_job_ID": gromos_system._last_jobID,
                "_no_double_submit_check": _no_double_submit_check,
                "_work_dir": _work_dir,
                "verbose": verbose,
                "verbose_lvl": verbose_lvl,
            }
        )

        try:
            utils.write_job_script(
                out_script_path=step_dir + "/schedule_MD_job.py",
                target_function=simulation_scheduler.do,
                variable_dict=MD_job_vars,
            )
        except Exception as err:
            raise Exception("Could not prepare the scheduling script\n\t" + "\n\t".join(map(str, err.args)))
    except Exception as err:
        traceback.print_exception(*sys.exc_info())
        raise Exception("Could not prepare the command block\n\t" + "\n\t".join(map(str, err.args)))

    # schedule
    try:
        if (os.path.exists(out_analysis_cnf) and os.path.exists(out_simulation_dir + ".tar")) and not force_simulation:
            if verbose:
                print(utils.spacer2 + "FOUND RESULT: " + out_analysis_cnf + "\n GOING TO SKIPT THIS SUBMISSION!")
            # warnings.warn("Skipping active submission, as result CNF was found: \n"+out_analysis_cnf)
            last_jobID = 0
        else:
            last_jobID = simulation_scheduler.do(**MD_job_vars)
    except Exception as err:
        traceback.print_exception(*sys.exc_info())
        raise Exception("Could not submit the commands\n\t" + "\n\t".join(map(str, err.args)))

    time.sleep(time_wait_s_for_filesystem)

    # Return the promise final system
    if in_analysis_script_path is None and os.path.exists(out_simulation_dir):
        cnfs = list(
            sorted(glob.glob(out_simulation_dir + "/*/*cnf"), key=lambda x: int(x.split("_")[-1].split(".")[0]))
        )
        last_cnf = cnfs[-1]
        gromos_system.cnf = cnf.Cnf(last_cnf)
    elif os.path.exists(out_analysis_cnf):
        gromos_system.cnf = cnf.Cnf(out_analysis_cnf)
    else:
        gromos_system.cnf = cnf.Cnf(in_value=None)
        gromos_system.cnf._future_file = True
        gromos_system.cnf.path = out_analysis_dir + "/data/" + gromos_system.name + ".cnf"

    # Return trajectories if available
    gromos_system._traj_files_path = {}
    if hasattr(gromos_system.imd, "WRITETRAJ") and gromos_system.imd.WRITETRAJ.NTWX > 0:
        final_trc_file = out_analysis_dir + "/data/" + gromos_system.name + ".trc"
        if os.path.exists(final_trc_file + ".h5"):
            gromos_system.trc = final_trc_file + ".h5"
        elif os.path.exists(final_trc_file):
            gromos_system.trc = final_trc_file
        else:
            gromos_system.trc = Trc(traj_path=None, in_cnf=None)
            gromos_system.trc._future_file = True
            gromos_system.trc.path = final_trc_file

    if hasattr(gromos_system.imd, "WRITETRAJ") and gromos_system.imd.WRITETRAJ.NTWE > 0:
        final_tre_file = out_analysis_dir + "/data/" + gromos_system.name + ".tre"
        if os.path.exists(final_tre_file + ".h5"):
            gromos_system.tre = final_tre_file + ".h5"
        elif os.path.exists(final_tre_file):
            gromos_system.tre = final_tre_file
        else:
            gromos_system.tre = Tre(input_value=None)
            gromos_system.tre._future_file = True
            gromos_system.tre.path = final_tre_file

    if hasattr(gromos_system.imd, "WRITETRAJ") and gromos_system.imd.WRITETRAJ.NTWG > 0:
        final_trg_file = out_analysis_dir + "/data/" + gromos_system.name + ".trg"
        if os.path.exists(final_trg_file + ".h5"):
            gromos_system.trg = Trg(input_value=final_trg_file + ".h5")
        elif os.path.exists(final_trg_file):
            gromos_system.trg = Trg(input_value=final_trg_file)
        else:
            gromos_system.trg = Trg(input_value=None)
            gromos_system.trg._future_file = True
            gromos_system.trg.path = final_trg_file

    # Set work_folder to initial work folder path to avoid nested simulations
    # However, do not update all file names (work_folder_no_update) to avoid
    # compromising the file system structure
    gromos_system.work_folder_no_update(init_work_folder)
    gromos_system._last_jobID = last_jobID
    os.chdir(init_work_folder)

    return gromos_system
