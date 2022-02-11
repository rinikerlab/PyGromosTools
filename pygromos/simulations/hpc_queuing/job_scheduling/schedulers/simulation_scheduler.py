"""
Wrapper for Long simulations -  similar to repex_EDS_long_production_run.
This script schedules Simulations ons euler into the queue.

"""

import os, sys, traceback

from pygromos.files.gromos_system import Gromos_System
from pygromos.simulations.hpc_queuing.job_scheduling.schedulers.scheduler_functions import chain_submission
from pygromos.simulations.hpc_queuing.job_scheduling.workers.simulation_workers import simulation_run_worker as workerScript
from pygromos.simulations.hpc_queuing.submission_systems._submission_system import _SubmissionSystem
from pygromos.simulations.hpc_queuing.submission_systems.submission_job import Submission_job
from pygromos.simulations.hpc_queuing.submission_systems.lsf import LSF
from pygromos.utils import bash, utils

spacer = utils.spacer

def do(in_simSystem: Gromos_System,
       out_dir_path: str,
       simulation_run_num: int,
       equilibration_run_num: int = 0,
       work_dir: str = None,
       analysis_script_path: str = None,
       submission_system:_SubmissionSystem = LSF(),
       previous_job_ID: int = None, no_double_submit:bool=True,
       verbose: bool = True, verbose_lvl:int=1):
    """

    Parameters
    ----------
    in_simSystem
    out_dir_path
    simulation_run_num
    equilibration_run_num
    gromos_bin_dir
    work_dir
    analysis_script_path
    submission_system
    previous_job_ID
    no_double_submit
    verbose
    verbose_lv:
        1 - slightly verbose
        2 - more info
        3 - let me tell you a story... let's start with Adam & Eve!....

    Returns
    -------

    """
    job_verb = True if (verbose and verbose_lvl > 2) else False

    # prepare
    try:
        if (verbose): print("Script: ", __file__)

        if (verbose and verbose_lvl > 2): print(spacer+"Simulation PREPERATION\n"+spacer+"\n")

        # Outdir
        bash.make_folder(out_dir_path)  # final output_folder

        # workdir:
        if (not isinstance(work_dir, type(None)) and work_dir != "None"):
            if (verbose and verbose_lvl>2): print("\t -> Generating given workdir: " + work_dir)
            bash.make_folder(work_dir, "-p")
            os.chdir(work_dir)
            prepared_imd = work_dir + "/" + os.path.basename(in_simSystem.imd.path)
        else:
            if (verbose and verbose_lvl>2): print("\t -> Using on node workdir")
            prepared_imd = out_dir_path + "/" + os.path.basename(in_simSystem.imd.path)

        # sim vars logs
        out_prefix = in_simSystem.name
        slave_script = workerScript.__file__

        # CHECK PATH DEPENDENCIES - all Files present?
        ##needed variables
        check_path_dependencies_paths = [slave_script, out_dir_path]  # Coord file is used by repex in_imd_path prepared_im
        ##variable paths
        if (not work_dir is None):
            check_path_dependencies_paths.append(work_dir)

        if(not in_simSystem.top._future_file):
            check_path_dependencies_paths.append(in_simSystem.top.path)
        if (not in_simSystem.cnf._future_file):
            check_path_dependencies_paths.append(in_simSystem.cnf.path)
        if (not in_simSystem.imd._future_file):
            check_path_dependencies_paths.append(in_simSystem.imd.path)
        if (not in_simSystem.ptp is  None):
            check_path_dependencies_paths.append(in_simSystem.ptp.path)
        if (not in_simSystem.disres is None):
            check_path_dependencies_paths.append(in_simSystem.disres.path)
        if (not in_simSystem.posres is None):
            check_path_dependencies_paths.append(in_simSystem.posres.path)
        if (not in_simSystem.refpos is None):
            check_path_dependencies_paths.append(in_simSystem.refpos.path)
        if (not in_simSystem.qmmm is None):
            check_path_dependencies_paths.append(in_simSystem.qmmm.path)

            #prepared_imd = bash.copy_file(in_simSystem.imd.path, prepared_imd) #Todo: Remove? @bschroed

        bash.check_path_dependencies(check_path_dependencies_paths, verbose=job_verb)

    except Exception as err:
        print("#####################################################################################")
        print("\t\tERROR in Preperations")
        print("#####################################################################################")

        traceback.print_exception(*sys.exc_info())

        raise IOError("ERROR in Preperations to submission!") from err

    # RUN Job
    try:

        if(verbose): print("\n"+spacer+"Simulation Setup:\n"+spacer)
        if (verbose): print("steps_per_run: ", in_simSystem.imd.STEP.NSTLIM)
        if (verbose): print("equis: ", equilibration_run_num)
        if (verbose): print("simulation runs: ", simulation_run_num)

        #Submission
        ## EQ
        eq_job_id=None
        if(equilibration_run_num>0):   #EQUILIBRATION
            tmp_outprefix = "eq_" + out_prefix
            tmp_jobname = in_simSystem.name + "_eq"
            eq_job_id, tmp_jobname, in_simSystem = chain_submission(simSystem=in_simSystem,
                             out_dir_path= out_dir_path, out_prefix=tmp_outprefix,
                             jobname=tmp_jobname,
                             chain_job_repetitions=equilibration_run_num, worker_script=workerScript.__file__,
                             job_submission_system=submission_system, start_run_index = 1,
                             prefix_command =  "", previous_job_ID = previous_job_ID, work_dir = work_dir,
                             initialize_first_run = True, reinitialize = False,
                             verbose = job_verb)

        ## MD
        tmp_outprefix = out_prefix
        tmp_jobname = in_simSystem.name
        previous_job_ID = previous_job_ID if(eq_job_id is None) else eq_job_id
        previous_job_ID, tmp_jobname, in_simSystem = chain_submission(simSystem=in_simSystem,
                         out_dir_path= out_dir_path, out_prefix=tmp_outprefix,
                         jobname=tmp_jobname,
                         chain_job_repetitions=equilibration_run_num + simulation_run_num, start_run_index = equilibration_run_num+1,
                         worker_script=workerScript.__file__,
                         job_submission_system=submission_system,
                         prefix_command =  "", previous_job_ID = previous_job_ID, work_dir = None,
                         initialize_first_run= False, reinitialize= False,
                         verbose = job_verb)

        ana_previous_job_ID = previous_job_ID
        if (not analysis_script_path is None):
            tmp_jobname = in_simSystem.name + "_ana"

            ana_log = os.path.dirname(analysis_script_path) + "/ana_out.log"
            if (verbose): print(spacer + "\n submit final analysis part \n")
            if (verbose): print(ana_log)
            if (verbose): print(analysis_script_path)

            sub_job = Submission_job(command="python3 "+analysis_script_path,
                                    jobName=tmp_jobname,
                                    outLog=ana_log,
                                    queue_after_jobID=previous_job_ID)
            ana_previous_job_ID = submission_system.submit_to_queue(sub_job)

            if (verbose): print("ANA jobID: " + str(previous_job_ID))

    except Exception as err:
        print("#####################################################################################")
        print("\t\tERROR in Submission")
        print("#####################################################################################")

        traceback.print_exception(*sys.exc_info())
        raise IOError("ERROR in SUBMISSION!") from err

    #in_simSystem._future_promise = False #reset future promising if necessary
    return ana_previous_job_ID
            