---
date: '2022-04-01T13:55:18.499Z'
docname: _source/pygromos.simulations.hpc_queuing.job_scheduling.schedulers
images: {}
path: /source-pygromos-simulations-hpc-queuing-job-scheduling-schedulers
title: pygromos.simulations.hpc_queuing.job_scheduling.schedulers package
---

# pygromos.simulations.hpc_queuing.job_scheduling.schedulers package

## Submodules

## pygromos.simulations.hpc_queuing.job_scheduling.schedulers.scheduler_functions module


### pygromos.simulations.hpc_queuing.job_scheduling.schedulers.scheduler_functions.chain_submission(simSystem: [pygromos.files.gromos_system.gromos_system.Gromos_System](#pygromos.files.gromos_system.gromos_system.Gromos_System), out_dir_path: str, out_prefix: str, chain_job_repetitions: int, worker_script: str, job_submission_system: pygromos.simulations.hpc_queuing.submission_systems._submission_system._SubmissionSystem, jobname: str, run_analysis_script_every_x_runs: int = 0, in_analysis_script_path: str = '', start_run_index: int = 1, prefix_command: str = '', previous_job_ID: Optional[int] = None, work_dir: Optional[str] = None, initialize_first_run: bool = True, reinitialize_every_run: bool = False, verbose: bool = False, verbose_lvl: int = 1)
> this function submits a chain of simulation steps to the queuing system and does the file managment.


* **Parameters**

    
    * **simSystem** (*Gromos_System*) – simulation system


    * **out_dir_path** (*str*) – out directory path


    * **out_prefix** (*str*) – out prefix for simulation files


    * **chain_job_repetitions** (*int*) – how often, should the simulation be repeated (in continuation)


    * **worker_script** (*str*) – worker, that should be submitted. This script will be executed at each scheduled job.


    * **job_submission_system** (*_SubmissionSystem*) – submission system, what type of submission?


    * **jobname** (*str*) – name of the simulation job


    * **run_analysis_script_every_x_runs** (*int, optional*) – run analysis in between - (careful will not be overwritten, make sure final analysis is correct.), by default 0


    * **in_analysis_script_path** (*str, optional*) – analysis script for simulation, that should be applied (will at least be applied after the full simulation chain.), by default “”


    * **start_run_index** (*int, optional*) – start index of the job chain., by default 1


    * **prefix_command** (*str, optional*) – any bash prefix commands, before submitting?, by default “”


    * **previous_job_ID** (*int, optional*) – ID of the prefious job, to be chained to. , by default None


    * **work_dir** (*str, optional*) – dir to wich the work in progress will be written. if None a tmp-srcatch dir will be used with LSF!, by default None


    * **initialize_first_run** (*bool, optional*) – should the velocities for the first run be initialized?, by default True


    * **reinitialize_every_run** (*bool, optional*) – should in every run, the velocities be reinitialized?, by default False


    * **verbose** (*bool, optional*) – more bla bla, by default False


    * **verbose_lvl** (*int, optional*) – nicely define ammount of bla bla, by default 1



* **Returns**

    Tuple[previous_job_ID, tmp_jobname, simSystem]
    will return the last job_ID, the last tmp_jobname and the final gromosSystem.



* **Return type**

    Tuple[int, str, [Gromos_System](#pygromos.files.gromos_system.gromos_system.Gromos_System)]



* **Raises**

    **ValueError** – if submission fails. This can habe various reasons, always check also the present files! (

    ```
    *
    ```

    omd etc.)



### pygromos.simulations.hpc_queuing.job_scheduling.schedulers.scheduler_functions.do_skip_job(tmp_out_cnf: str, simSystem: [pygromos.files.gromos_system.gromos_system.Gromos_System](#pygromos.files.gromos_system.gromos_system.Gromos_System), tmp_jobname: str, job_submission_system: pygromos.simulations.hpc_queuing.submission_systems._submission_system._SubmissionSystem, previous_job: int, verbose: bool = True, verbose_lvl: int = 1)
> This function returns true if the job was already submitted or already finished, as well as the job id.


* **Parameters**

    
    * **tmp_out_cnf** (*str*) – the possible final out cnf


    * **simSystem** (*Gromos_System*) – the simulation system


    * **tmp_jobname** (*str*) – possible job-name


    * **job_submission_system** (*_SubmissionSystem*) – submission system


    * **previous_job** (*int*) – job, that was submitted before this one.


    * **verbose** (*bool, optional*) – more bla bla, by default True


    * **verbose_lvl** (*int, optional*) – nicely define ammount of bla bla, by default 1



* **Returns**

    the bool tells, if the job was already submitted,
    the id is the possibly found job id of the job, if it is in the queue.



* **Return type**

    Tuple[bool, int]



* **Raises**

    
    * **IOError** – if a tmp_out_cnf is present, we assume the simulation run was successfull. If it can not be parsed in, something else is a problem, but we should not continue our chain!


    * **ValueError** – if multiple jobs with same name were found in the queue, we assume, the search regex did not work.


## pygromos.simulations.hpc_queuing.job_scheduling.schedulers.simulation_scheduler module

Wrapper for Long simulations -  similar to repex_EDS_long_production_run.
This script schedules Simulations ons euler into the queue.


### pygromos.simulations.hpc_queuing.job_scheduling.schedulers.simulation_scheduler.do(in_simSystem: pygromos.files.gromos_system.gromos_system.Gromos_System, out_dir_path: str, simulation_run_num: int, equilibration_run_num: int = 0, initialize_first_run=False, reinitialize_every_run=False, analysis_script_path: typing.Optional[str] = None, submission_system: pygromos.simulations.hpc_queuing.submission_systems._submission_system._SubmissionSystem = <pygromos.simulations.hpc_queuing.submission_systems.lsf.LSF object>, previous_job_ID: typing.Optional[int] = None, _no_double_submit_check: bool = False, _work_dir: typing.Optional[str] = None, verbose: bool = True, verbose_lvl: int = 1)
> This function schedules simulations starting from the gromos system.


* **Parameters**

    
    * **in_simSystem** (*Gromos_System*) – system that should be submitted with the provided imd file


    * **out_dir_path** (*str*) – out directory  path for the simulation


    * **simulation_run_num** (*int*) – number of simulations


    * **equilibration_run_num** (*int, optional*) – number of the equilibraitons, by default 0


    * **initialize_first_run** (*bool, optional*) – should the velocities be initialized in the first run?, by default False


    * **reinitialize_every_run** (*bool, optional*) – DEAPPRECIATED! should always the velocities be initialized, by default False


    * **analysis_script_path** (*str, optional*) – path to the analysis script, that should be used for this simulaiton approach, by default None


    * **submission_system** (*_SubmissionSystem, optional*) – system, to be used to submit the jobs, by default LSF()


    * **previous_job_ID** (*int, optional*) – previous job ID, by default None


    * **_no_double_submit_check** (*bool, optional*) – don’t check if job was already submit to queue (increases performance!), by default False


    * **_work_dir** (*str, optional*) – directory, to write out the tmp files, by default None


    * **verbose** (*bool, optional*) – Baeh Baeh, by default True


    * **verbose_lvl** (*int, optional*) – amount of Baehs, by default 1



* **Returns**

    the last job id, that was submitted.



* **Return type**

    int



* **Raises**

    **IOError** – If error happens in preperation of simulation or in the submission


## Module contents
