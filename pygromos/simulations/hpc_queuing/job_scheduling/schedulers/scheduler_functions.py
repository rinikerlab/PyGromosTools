import os
import glob
import pandas as pd

from pygromos.utils.typing import Tuple
from pygromos.files.coord.cnf import Cnf
from pygromos.files.gromos_system import Gromos_System
from pygromos.simulations.hpc_queuing.submission_systems._submission_system import _SubmissionSystem
from pygromos.simulations.hpc_queuing.submission_systems.submission_job import Submission_job
from pygromos.simulations.hpc_queuing.submission_systems.local import LOCAL
from pygromos.utils import bash
from pygromos.utils.utils import spacer3


def do_skip_job(
    tmp_out_cnf: str,
    simSystem: Gromos_System,
    tmp_jobname: str,
    job_submission_system: _SubmissionSystem,
    previous_job: int,
    verbose: bool = True,
    verbose_lvl: int = 1,
) -> Tuple[bool, int]:
    """
        This function returns true if the job was already submitted or already finished, as well as the job id.

    Parameters
    ----------
    tmp_out_cnf : str
        the possible final out cnf
    simSystem : Gromos_System
        the simulation system
    tmp_jobname : str
        possible job-name
    job_submission_system : _SubmissionSystem
        submission system
    previous_job : int
        job, that was submitted before this one.
    verbose : bool, optional
        more bla bla, by default True
    verbose_lvl : int, optional
        nicely define ammount of bla bla, by default 1

    Returns
    -------
    Tuple[bool, int]
        the bool tells, if the job was already submitted,
        the id is the possibly found job id of the job, if it is in the queue.

    Raises
    ------
    IOError
        if a tmp_out_cnf is present, we assume the simulation run was successfull. If it can not be parsed in, something else is a problem, but we should not continue our chain!
    ValueError
        if multiple jobs with same name were found in the queue, we assume, the search regex did not work.
    """

    # Check if job with same name is already in the queue!
    if (verbose) and verbose_lvl >= 2:
        print("Checking if jobs was already submitted or done")
    if job_submission_system.block_double_submission:  # can we find an job with this name in the queue?
        if (verbose) and verbose_lvl >= 2:
            print("Checking for jobs with name: " + tmp_jobname)
        queued_job_ids = job_submission_system.search_queue_for_jobname(job_name=tmp_jobname)

        if isinstance(queued_job_ids, (pd.DataFrame, pd.Series)):
            queued_job_ids = list(queued_job_ids.where(queued_job_ids.STAT.isin(["RUN", "PEND"])).dropna().JOBID)

        # check if already submitted
        if len(queued_job_ids) > 1:  # check if job is already submitted:
            if verbose:
                print(
                    "\t\t\tSKIP job Submission: "
                    + tmp_jobname
                    + " was already submitted to the queue! \n\t\t\t\tSKIP\n"
                    + "\t\t\tsubmitted IDS: "
                    + "\n\t\t\t".join(map(str, queued_job_ids))
                    + "\n"
                )
            try:
                simSystem.cnf = tmp_out_cnf
            except Exception as err:
                raise IOError("Tried to read found cnf file: " + tmp_out_cnf + ". FAILED! \n" + str(err.args))

            if len(queued_job_ids) == 1:
                previous_job = queued_job_ids[0]
                if (verbose) and verbose_lvl >= 2:
                    print("\nTRY to attach next job to ", previous_job, "\n")
            else:
                raise ValueError(
                    "\nthere are multiple jobs, that could be the precessor. " + " ".join(map(str, queued_job_ids))
                )
            if verbose:
                print(True)
            if verbose:
                print()
            return True, previous_job

    # Check if run was already finished:
    tmp_out_cnfs_regex = "_".join(tmp_out_cnf.split("_")[:-1]) + "*.cnf"
    if verbose:
        print("Checking for resulting files: " + tmp_out_cnfs_regex)

    if len(glob.glob(tmp_out_cnfs_regex)) > 0:  # was this job already run and finished?
        if verbose:
            print("\t\t NOT SUBMITTED!(inScript) as these Files were found: \n\t" + tmp_out_cnfs_regex)
        setattr(simSystem, "coord_seeds", tmp_out_cnf)  # set next coord Files

        if (verbose) and verbose_lvl >= 2:
            print(simSystem.cnf.path)
        if (verbose) and verbose_lvl >= 2:
            print(True)
        if (verbose) and verbose_lvl >= 2:
            print()
        return True, None
    if verbose:
        print(False)
    if verbose:
        print()
    return False, previous_job


def chain_submission(
    simSystem: Gromos_System,
    out_dir_path: str,
    out_prefix: str,
    chain_job_repetitions: int,
    worker_script: str,
    job_submission_system: _SubmissionSystem,
    jobname: str,
    run_analysis_script_every_x_runs: int = 0,
    in_analysis_script_path: str = "",
    start_run_index: int = 1,
    prefix_command: str = "",
    previous_job_ID: int = None,
    work_dir: str = None,
    initialize_first_run: bool = True,
    reinitialize_every_run: bool = False,
    verbose: bool = False,
    verbose_lvl: int = 1,
) -> Tuple[int, str, Gromos_System]:
    """
        this function submits a chain of simulation steps to the queuing system and does the file managment.


    Parameters
    ----------
    simSystem : Gromos_System
        simulation system
    out_dir_path : str
        out directory path
    out_prefix : str
        out prefix for simulation files
    chain_job_repetitions : int
        how often, should the simulation be repeated (in continuation)
    worker_script : str
        worker, that should be submitted. This script will be executed at each scheduled job.
    job_submission_system : _SubmissionSystem
        submission system, what type of submission?
    jobname : str
        name of the simulation job
    run_analysis_script_every_x_runs : int, optional
        run analysis in between - (careful will not be overwritten, make sure final analysis is correct.), by default 0
    in_analysis_script_path : str, optional
        analysis script for simulation, that should be applied (will at least be applied after the full simulation chain.), by default ""
    start_run_index : int, optional
        start index of the job chain., by default 1
    prefix_command : str, optional
        any bash prefix commands, before submitting?, by default ""
    previous_job_ID : int, optional
        ID of the prefious job, to be chained to. , by default None
    work_dir : str, optional
        dir to wich the work in progress will be written. if None a tmp-srcatch dir will be used with LSF!, by default None
    initialize_first_run : bool, optional
        should the velocities for the first run be initialized?, by default True
    reinitialize_every_run : bool, optional
        should in every run, the velocities be reinitialized?, by default False
    verbose : bool, optional
        more bla bla, by default False
    verbose_lvl : int, optional
        nicely define ammount of bla bla, by default 1

    Returns
    -------
    Tuple[int, str, Gromos_System]
        Tuple[previous_job_ID, tmp_jobname, simSystem]
        will return the last job_ID, the last tmp_jobname and the final gromosSystem.

    Raises
    ------
    ValueError
        if submission fails. This can habe various reasons, always check also the present files! (*omd etc.)
    """

    if verbose:
        print("\nChainSubmission - " + out_prefix + "\n" + "=" * 30 + "\n")
    if (verbose) and verbose_lvl >= 2:
        print("start_run_index " + str(start_run_index))
    if (verbose) and verbose_lvl >= 2:
        print("job reptitions " + str(chain_job_repetitions))

    if job_submission_system is not LOCAL:
        simSystem._future_promise = True

    ana_id = None
    job_submission_system.job_duration = job_submission_system.job_duration
    for runID in range(start_run_index, chain_job_repetitions + 1):

        if verbose:
            print("\n submit  " + jobname + "_" + str(runID) + "\n" + spacer3)

        tmp_outprefix = out_prefix + "_" + str(runID)
        tmp_jobname = jobname + "_" + str(runID)
        tmp_outdir = out_dir_path + "/" + tmp_outprefix
        tmp_out_cnf = tmp_outdir + "/" + tmp_outprefix + ".cnf"

        # Checks if run should be skipped!
        do_skip, previous_job_ID = do_skip_job(
            tmp_out_cnf=tmp_out_cnf,
            simSystem=simSystem,
            tmp_jobname=tmp_jobname,
            job_submission_system=job_submission_system,
            previous_job=previous_job_ID,
            verbose=verbose,
        )

        if not do_skip:
            bash.make_folder(tmp_outdir)

            # build COMMANDS:
            if len(prefix_command) > 1:
                prefix_command += " && "

            # We will write the arguments to the python script in a bash array
            # to make it simpler to read in our input files.

            md_args = "md_args=(\n"

            md_args += "-out_dir " + tmp_outdir + "\n"
            md_args += "-in_cnf_path " + simSystem.cnf.path + "\n"
            md_args += "-in_imd_path " + simSystem.imd.path + "\n"
            md_args += "-in_top_path " + simSystem.top.path + "\n"
            md_args += "-runID " + str(runID) + "\n"

            # OPTIONAL ARGS
            if simSystem.disres is not None:
                md_args += "-in_disres_path " + simSystem.disres.path + "\n"
            if simSystem.ptp is not None:
                md_args += "-in_perttopo_path " + simSystem.ptp.path + "\n"
            if simSystem.refpos is not None:
                md_args += "-in_refpos_path " + simSystem.refpos.path + "\n"
            if simSystem.qmmm is not None:
                md_args += "-in_qmmm_path " + simSystem.qmmm.path + " "
            if simSystem.posres is not None:
                md_args += "-in_posres_path " + simSystem.posres.path + "\n"

            md_args += "-nmpi " + str(job_submission_system.nmpi) + "\n"
            md_args += "-nomp " + str(job_submission_system.nomp) + "\n"
            md_args += "-initialize_first_run " + str(initialize_first_run) + "\n"
            md_args += "-reinitialize_every_run " + str(reinitialize_every_run) + "\n"
            md_args += "-gromosXX_bin_dir " + str(simSystem.gromosXX.bin) + "\n"
            md_args += "-gromosXX_check_binary_paths " + str(simSystem.gromosXX._check_binary_paths) + "\n"

            if work_dir is not None:
                md_args += "-work_dir " + str(work_dir) + "\n"

            if hasattr(simSystem.imd, "WRITETRAJ"):
                if simSystem.imd.WRITETRAJ.NTWX > 0:
                    md_args += "-out_trc " + str(True) + "\n"
                if simSystem.imd.WRITETRAJ.NTWE > 0:
                    md_args += "-out_tre " + str(True) + "\n"
                if simSystem.imd.WRITETRAJ.NTWV > 0:
                    md_args += "-out_trv " + str(True) + "\n"
                if simSystem.imd.WRITETRAJ.NTWF > 0:
                    md_args += "-out_trf " + str(True) + "\n"
                if simSystem.imd.WRITETRAJ.NTWG > 0:
                    md_args += "-out_trg " + str(True) + "\n"

            md_args += "-zip_trajectories " + str(job_submission_system.zip_trajectories) + "\n"

            md_args += ")\n"  # closing the bash array which stores all arguments.

            # add zip option here

            # MAIN commands
            md_script_command = prefix_command + "\n\n" + md_args + "\n"
            md_script_command += "python3 " + worker_script + '  "${md_args[@]}" \n'

            if verbose:
                print("PREVIOUS ID: ", previous_job_ID)
                if verbose_lvl >= 2:
                    print("COMMAND: ", md_script_command)

            # SCHEDULE THE COMMANDS
            try:
                if verbose:
                    print("\tSIMULATION")
                os.chdir(tmp_outdir)
                sub_job = Submission_job(
                    command=md_script_command,
                    jobName=tmp_jobname,
                    submit_from_dir=tmp_outdir,
                    queue_after_jobID=previous_job_ID,
                    outLog=tmp_outdir + "/" + out_prefix + "_md.out",
                    errLog=tmp_outdir + "/" + out_prefix + "_md.err",
                    sumbit_from_file=True,
                )
                previous_job_ID = job_submission_system.submit_to_queue(sub_job)
                if verbose:
                    print("SIMULATION ID: ", previous_job_ID)
            except ValueError as err:  # job already in the queue
                raise ValueError(
                    "ERROR during submission of main job " + str(tmp_jobname) + ":\n" + "\n".join(err.args)
                )

            # OPTIONAL schedule - analysis inbetween.
            if (
                runID > 1
                and run_analysis_script_every_x_runs != 0
                and runID % run_analysis_script_every_x_runs == 0
                and runID < chain_job_repetitions
            ):

                if (verbose) and verbose_lvl >= 2:
                    print("\tINBETWEEN ANALYSIS")
                sub_job = Submission_job(
                    command=in_analysis_script_path,
                    jobName=jobname + "_intermediate_ana_run_" + str(runID),
                    outLog=tmp_outdir + "/" + out_prefix + "_inbetweenAna.out",
                    errLog=tmp_outdir + "/" + out_prefix + "_inbetweenAna.err",
                    queue_after_jobID=previous_job_ID,
                )
                try:
                    ana_id = job_submission_system.submit_to_queue(sub_job)
                    if (verbose) and verbose_lvl >= 2:
                        print("\n")
                except ValueError as err:  # job already in the queue
                    print("ERROR during submission of analysis command of " + sub_job.jobName + ":\n")
                    print("\n".join(err.args))
        else:
            if (verbose) and verbose_lvl >= 2:
                print("Did not submit!")
        if (verbose) and verbose_lvl >= 2:
            print("\n")
        if (verbose) and verbose_lvl >= 2:
            print("job_postprocess ")
        prefix_command = ""

        # Resulting cnf is provided to use it in further approaches.
        simSystem.cnf = Cnf(tmp_out_cnf, _future_file=True)

    if ana_id is not None:
        previous_job_ID = ana_id

    return previous_job_ID, tmp_jobname, simSystem
