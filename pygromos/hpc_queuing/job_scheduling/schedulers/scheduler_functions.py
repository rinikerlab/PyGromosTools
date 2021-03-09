import glob
import os

from pygromos.files.gromos_system import Gromos_System
from pygromos.hpc_queuing.submission_systems.Submission_Systems import _SubmissionSystem, LOCAL
from pygromos.hpc_queuing.job_scheduling.workers.simulation_workers import clean_up_simulation_files
from pygromos.utils import bash
from pygromos.utils.utils import spacer, spacer3


def do_skip_job(tmp_out_cnf: str, simSystem: Gromos_System,
                tmp_jobname: str, job_submission_system: _SubmissionSystem, previous_job: int,
                do_not_doubly_submit_to_queue: bool = True, 
                verbose: bool = True):

    # Check if job with same name is already in the queue!
    if (verbose): print("Checking if jobs was already submitted or done")
    if (do_not_doubly_submit_to_queue):  # can we find an job with this name in the queue?
        if (verbose): print("Checking for jobs with name: " + tmp_jobname)
        queued_job_ids = job_submission_system.get_jobs_from_queue(job_text=tmp_jobname)

        ## check if already submitted
        if (len(queued_job_ids) > 0):  # check if job is already submitted:
            if (True): print(
                "\t\t\tSKIP job Submission: " + tmp_jobname + " was already submitted to the queue! \n\t\t\t\tSKIP\n"
                + "\t\t\tsubmitted IDS: " + "\n\t\t\t".join(map(str, queued_job_ids)) + "\n")
            try:
                simSystem.cnf = tmp_out_cnf
            except Exception as err:
                raise IOError("Tried to read found cnf file: "+tmp_out_cnf+". FAILED! \n"+str(err.args))
            if (len(queued_job_ids) == 1):
                previous_job = queued_job_ids[0]
                if (verbose): print("\nTRY to attach next job to ", previous_job, "\n")
            else:
                raise ValueError("\nthere are multiple jobs, that could be the precessor. " + " ".join(map(str, queued_job_ids)))
            if (verbose): print(True)
            if (verbose): print()
            return True, previous_job

    # Check if run was already finished:
    tmp_out_cnfs_regex = "_".join(tmp_out_cnf.split("_")[:-1]) + "*.cnf"
    if (verbose): print("Checking for resulting files: " + tmp_out_cnfs_regex)

    if (len(glob.glob(tmp_out_cnfs_regex)) > 0):  # was this job already run and finished?
        if (True): print(
            "\t\t NOT SUBMITTED!(inScript) as these Files were found: \n\t" + tmp_out_cnfs_regex)
        setattr(simSystem, "coord_seeds", tmp_out_cnf)  # set next coord Files
        prefix_command = ""
        if (verbose): print(simSystem.cnf.path)
        if (verbose): print(True)
        if (verbose): print()
        return True, None
    if (verbose): print(False)
    if(verbose): print()
    return False, previous_job


def chain_submission(simSystem:Gromos_System,
                     out_dir_path: str, out_prefix: str,
                     chain_job_repetitions: int, worker_script: str,
                     job_submission_system: _SubmissionSystem, jobname: str, nmpi: int=1, nomp:int=1,
                     job_queue_duration: str = "24:00",
                     run_analysis_script_every_x_runs: int = 0, in_analysis_script_path: str = "",
                     do_not_doubly_submit_to_queue: bool = True, start_run_index: int = 1,
                     prefix_command: str = "", previous_job_ID: int = None, work_dir: str = None,
                     initialize_first_run: bool = True, reinitialize: bool = False,
                     verbose: bool = False):
    """

    Parameters
    ----------
    simSystem
    out_dir_path
    out_prefix
    chain_job_repetitions
    worker_script
    job_submission_system
    jobname
    nmpi
    job_queue_duration
    run_analysis_script_every_x_runs
    in_analysis_script_path
    do_not_doubly_submit_to_queue
    start_run_index
    prefix_command
    previous_job_ID
    work_dir
    initialize_first_run
    reinitialize
        initialize_first_run must be False
    verbose

    Returns
    -------

    """
    print("\nChainSubmission - "+out_prefix+"\n"+"="*30+"\n")
    if(verbose): print("start_run_index " + str(start_run_index))
    if(verbose): print("job reptitions " + str(chain_job_repetitions))
    orig_prefix_command = prefix_command

    if(not job_submission_system is LOCAL):
        simSystem._future_promise = True
    job_submission_system.job_duration = job_queue_duration
    for runID in range(start_run_index, chain_job_repetitions + 1):

        print("\n submit  " + jobname + "_" + str(runID) + "\n"+spacer3)

        tmp_outprefix = out_prefix + "_" + str(runID)
        tmp_jobname = jobname + "_" + str(runID)
        tmp_outdir = out_dir_path + "/" + tmp_outprefix
        tmp_out_cnf = tmp_outdir + "/" + tmp_outprefix + ".cnf"

        # Checks if run should be skipped!
        do_skip, previous_job_ID = do_skip_job(tmp_out_cnf=tmp_out_cnf, simSystem=simSystem,
                                               do_not_doubly_submit_to_queue=do_not_doubly_submit_to_queue,
                                               tmp_jobname=tmp_jobname, job_submission_system=job_submission_system,
                                               previous_job=previous_job_ID, verbose=verbose)

        if (not do_skip):
            bash.make_folder(tmp_outdir)

            # build COMMANDS:
            prefix_command += " sleep 60s "

            # MAIN commands
            md_script_command = prefix_command + " && "
            md_script_command += "python " + worker_script + " "
            md_script_command += "-out_dir " + tmp_outdir + " "
            md_script_command += "-in_cnf_path " + simSystem.cnf.path + " "
            md_script_command += "-in_imd_path " + simSystem.imd.path + " "
            md_script_command += "-in_top_path " + simSystem.top.path + " "
            md_script_command += "-runID " + str(runID) + " "

            ## OPTIONAL ARGS
            if (not simSystem.disres is None):
                md_script_command += "-in_disres_path " + simSystem.disres.path + " "
            if (not simSystem.ptp is None):
                md_script_command += "-in_perttopo_path " + simSystem.ptp.path + " "
            if (not simSystem.refpos is None):
                md_script_command += "-in_refpos_path " + simSystem.refpos.path + " "
            if (not simSystem.posres is None):
                md_script_command += "-in_posres_path " + simSystem.posres.path + " "

            #if(out_trg)

            md_script_command += "-nmpi " + str(nmpi) + " "
            md_script_command += "-nomp " + str(nomp) + " "
            md_script_command += "-initialize_first_run "+str(initialize_first_run)+ " "
            md_script_command += "-gromosXX_bin_dir " + str(simSystem.gromosXX.bin) + " "
            if(not work_dir is None):
                md_script_command += "-work_dir " + str(work_dir) + " "

            if(hasattr(simSystem.imd, "WRITETRAJ")):
                if(simSystem.imd.WRITETRAJ.NTWX > 0):
                    md_script_command += "-out_trc "+str(True)+" "
                if(simSystem.imd.WRITETRAJ.NTWE > 0):
                    md_script_command += "-out_tre "+str(True)+" "
                if(simSystem.imd.WRITETRAJ.NTWV > 0):
                    md_script_command += "-out_trv "+str(True)+" "
                if(simSystem.imd.WRITETRAJ.NTWF > 0):
                    md_script_command += "-out_trf "+str(True)+" "
                if(simSystem.imd.WRITETRAJ.NTWG > 0):
                    md_script_command += "-out_trg "+str(True)+" "

            if(verbose): print("COMMAND: ", md_script_command)

            ## POST COMMAND
            clean_up_processes = nomp if (nomp > nmpi) else nmpi
            clean_up_command = "python " + str(clean_up_simulation_files.__file__) + "  -in_simulation_dir " + str(
                tmp_outdir) + " -n_processes " + str(clean_up_processes)

            print("PREVIOUS ID: ", previous_job_ID)

            # SCHEDULE THE COMMANDS
            try:
                if (verbose): print("\tSIMULATION")
                os.chdir(tmp_outdir)
                outLog = tmp_outdir + "/" + out_prefix + "_md.out"
                errLog = tmp_outdir + "/" + out_prefix + "_md.err"
                previous_job_ID = job_submission_system.submit_to_queue(command=md_script_command, jobName=tmp_jobname,
                                                                        submit_from_dir=tmp_outdir,
                                                                        queue_after_jobID=previous_job_ID,
                                                                        outLog=outLog, errLog=errLog, sumbit_from_file=True,
                                                                        nmpi=nmpi, end_mail=True, verbose=verbose)
                print("SIMULATION ID: ", previous_job_ID)
            except ValueError as err:  # job already in the queue
                raise ValueError("ERROR during submission of main job "+str(tmp_jobname)+":\n"+"\n".join(err.args))
            try:
                # schedule - simulation cleanup:
                ##this mainly tars files.
                if (verbose): print("\tCLEANING")
                outLog = tmp_outdir + "/" + out_prefix + "_cleanUp.out"
                errLog = tmp_outdir + "/" + out_prefix + "_cleanUp.err"
                clean_id = job_submission_system.submit_to_queue(command=clean_up_command,
                                                                 jobName=tmp_jobname + "_cleanUP",
                                                                 queue_after_jobID=previous_job_ID,
                                                                 outLog=outLog, errLog=errLog,
                                                                 nmpi=nmpi, verbose=verbose)
                print("CLEANING ID: ", previous_job_ID)

            except ValueError as err:  # job already in the queue
                raise ValueError("ERROR during submission of clean-up command of "+str(tmp_jobname)+"_cleanUP:\n"+"\n".join(err.args))

            # OPTIONAL schedule - analysis inbetween.
            if (runID > 1 and run_analysis_script_every_x_runs != 0 and
                    runID % run_analysis_script_every_x_runs == 0
                    and runID < chain_job_repetitions):

                if (verbose): print("\tINBETWEEN ANALYSIS")
                tmp_ana_jobname = jobname + "_intermediate_ana_run_" + str(runID)
                outLog = tmp_outdir + "/" + out_prefix + "_inbetweenAna.out"
                errLog = tmp_outdir + "/" + out_prefix + "_inbetweenAna.err"

                try:
                    ana_id = job_submission_system.submit_to_queue(command=in_analysis_script_path,
                                                                   jobName=tmp_ana_jobname,
                                                                   outLog=outLog, errLog=errLog,
                                                                   maxStorage=20000, queue_after_jobID=clean_id, nmpi=5,
                                                                   verbose=verbose)
                    if (verbose): print("\n")
                except ValueError as err:  # job already in the queue
                    print("ERROR during submission of analysis command of "+str(tmp_ana_jobname)+":\n")
                    print("\n".join(err.args))
        else:
            if(verbose): print("Did not submit!")
            clean_id = None
        if (verbose): print("\n")
        if (verbose): print("job_postprocess ")
        prefix_command = ""

        #Resulting cnf is provided to use it in further approaches.
        simSystem.cnf = tmp_out_cnf

    previous_job_ID = clean_id
    return previous_job_ID, tmp_jobname, simSystem
