"""
Wrapper for Long simulations -  similar to repex_EDS_long_production_run.
This script schedules Simulations ons euler into the queue.

"""

import os
import sys
import traceback

from pygromos.files.gromos_system import Gromos_System
from pygromos.simulations.hpc_queuing.job_scheduling.schedulers.scheduler_functions import chain_submission
from pygromos.simulations.hpc_queuing.job_scheduling.workers.simulation_workers import (
    simulation_run_worker as workerScript,
)
from pygromos.simulations.hpc_queuing.submission_systems._submission_system import _SubmissionSystem
from pygromos.simulations.hpc_queuing.submission_systems.submission_job import Submission_job
from pygromos.simulations.hpc_queuing.submission_systems.lsf import LSF
from pygromos.utils import bash, utils

spacer = utils.spacer


def do(
    in_simSystem: Gromos_System,
    out_dir_path: str,
    simulation_run_num: int,
    equilibration_run_num: int = 0,
    initialize_first_run=False,
    reinitialize_every_run=False,
    analysis_script_path: str = None,
    submission_system: _SubmissionSystem = LSF(),
    previous_job_ID: int = None,
    _no_double_submit_check: bool = False,
    _work_dir: str = None,
    verbose: bool = True,
    verbose_lvl: int = 1,
) -> int:
    """
        This function schedules simulations starting from the gromos system.

    Parameters
    ----------
    in_simSystem : Gromos_System
        system that should be submitted with the provided imd file
    out_dir_path : str
        out directory  path for the simulation
    simulation_run_num : int
        number of simulations
    equilibration_run_num : int, optional
        number of the equilibraitons, by default 0
    initialize_first_run : bool, optional
        should the velocities be initialized in the first run?, by default False
    reinitialize_every_run : bool, optional
        DEAPPRECIATED! should always the velocities be initialized, by default False
    analysis_script_path : str, optional
        path to the analysis script, that should be used for this simulaiton approach, by default None
    submission_system : _SubmissionSystem, optional
        system, to be used to submit the jobs, by default LSF()
    previous_job_ID : int, optional
        previous job ID, by default None
    _no_double_submit_check : bool, optional
        don't check if job was already submit to queue (increases performance!), by default False
    _work_dir : str, optional
        directory, to write out the tmp files, by default None
    verbose : bool, optional
        Baeh Baeh, by default True
    verbose_lvl : int, optional
        amount of Baehs, by default 1

    Returns
    -------
    int
        the last job id, that was submitted.

    Raises
    ------
    IOError
        If error happens in preperation of simulation or in the submission

    """

    submission_system.block_double_submission = _no_double_submit_check
    job_verb = True if (verbose and verbose_lvl > 2) else False

    # prepare
    try:
        if verbose:
            print("Script: ", __file__)

        if verbose and verbose_lvl > 2:
            print(spacer + "Simulation PREPERATION\n" + spacer + "\n")

        # Outdir
        bash.make_folder(out_dir_path)  # final output_folder

        # workdir:
        if not isinstance(_work_dir, type(None)) and _work_dir != "None":
            if verbose and verbose_lvl > 2:
                print("\t -> Generating given workdir: " + _work_dir)
            bash.make_folder(_work_dir, "-p")
            os.chdir(_work_dir)
        else:
            if verbose and verbose_lvl > 2:
                print("\t -> Using on node workdir")

        # sim vars logs
        out_prefix = in_simSystem.name
        worker_script = workerScript.__file__

        # CHECK PATH DEPENDENCIES - all Files present?
        # needed variables
        check_path_dependencies_paths = [
            worker_script,
            out_dir_path,
        ]  # Coord file is used by repex in_imd_path prepared_im
        # variable paths
        if _work_dir is not None and _work_dir != "out_dir":
            check_path_dependencies_paths.append(_work_dir)

        if not in_simSystem.top._future_file:
            check_path_dependencies_paths.append(in_simSystem.top.path)
        if not in_simSystem.cnf._future_file:
            check_path_dependencies_paths.append(in_simSystem.cnf.path)
        if not in_simSystem.imd._future_file:
            check_path_dependencies_paths.append(in_simSystem.imd.path)
        if in_simSystem.ptp is not None:
            check_path_dependencies_paths.append(in_simSystem.ptp.path)
        if in_simSystem.disres is not None:
            check_path_dependencies_paths.append(in_simSystem.disres.path)
        if in_simSystem.posres is not None:
            check_path_dependencies_paths.append(in_simSystem.posres.path)
        if in_simSystem.refpos is not None:
            check_path_dependencies_paths.append(in_simSystem.refpos.path)
        if in_simSystem.qmmm is not None:
            check_path_dependencies_paths.append(in_simSystem.qmmm.path)

        bash.check_path_dependencies(check_path_dependencies_paths, verbose=job_verb)

    except Exception as err:
        print("#####################################################################################")
        print("\t\tERROR in Preperations")
        print("#####################################################################################")

        traceback.print_exception(*sys.exc_info())

        raise IOError("ERROR in Preperations to submission!") from err

    # RUN Job
    try:

        if verbose:
            print("\n" + spacer + "Simulation Setup:\n" + spacer)
        if verbose:
            print("steps_per_run: ", in_simSystem.imd.STEP.NSTLIM)
        if verbose:
            print("equis: ", equilibration_run_num)
        if verbose:
            print("simulation runs: ", simulation_run_num)

        # Submission
        # EQ
        eq_job_id = None
        if equilibration_run_num > 0:  # EQUILIBRATION
            tmp_outprefix = "eq_" + out_prefix
            tmp_jobname = in_simSystem.name + "_eq"
            eq_job_id, tmp_jobname, in_simSystem = chain_submission(
                simSystem=in_simSystem,
                out_dir_path=out_dir_path,
                out_prefix=tmp_outprefix,
                jobname=tmp_jobname,
                chain_job_repetitions=equilibration_run_num,
                worker_script=workerScript.__file__,
                job_submission_system=submission_system,
                start_run_index=1,
                prefix_command="",
                previous_job_ID=previous_job_ID,
                work_dir=_work_dir,
                initialize_first_run=initialize_first_run,
                reinitialize_every_run=reinitialize_every_run,
                verbose=job_verb,
            )

        # MD
        tmp_outprefix = out_prefix
        tmp_jobname = in_simSystem.name
        previous_job_ID = previous_job_ID if (eq_job_id is None) else eq_job_id
        previous_job_ID, tmp_jobname, in_simSystem = chain_submission(
            simSystem=in_simSystem,
            out_dir_path=out_dir_path,
            out_prefix=tmp_outprefix,
            jobname=tmp_jobname,
            chain_job_repetitions=equilibration_run_num + simulation_run_num,
            start_run_index=equilibration_run_num + 1,
            worker_script=workerScript.__file__,
            job_submission_system=submission_system,
            prefix_command="",
            previous_job_ID=previous_job_ID,
            work_dir=_work_dir,
            initialize_first_run=initialize_first_run,
            reinitialize_every_run=reinitialize_every_run,
            verbose=job_verb,
        )

        ana_previous_job_ID = previous_job_ID
        if analysis_script_path is not None:
            tmp_jobname = in_simSystem.name + "_ana"

            ana_log = os.path.dirname(analysis_script_path) + "/ana_out.log"
            if verbose:
                print(spacer + "\n submit final analysis part \n")
            if verbose:
                print(ana_log)
            if verbose:
                print(analysis_script_path)

            sub_job = Submission_job(
                command="python3 " + analysis_script_path,
                jobName=tmp_jobname,
                outLog=ana_log,
                queue_after_jobID=previous_job_ID,
            )
            ana_previous_job_ID = submission_system.submit_to_queue(sub_job)

            if verbose:
                print("ANA jobID: " + str(previous_job_ID))

    except Exception as err:
        print("#####################################################################################")
        print("\t\tERROR in Submission")
        print("#####################################################################################")

        traceback.print_exception(*sys.exc_info())
        raise IOError("ERROR in SUBMISSION!") from err

    # in_simSystem._future_promise = False #reset future promising if necessary
    return ana_previous_job_ID
