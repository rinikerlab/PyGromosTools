---
date: '2022-04-01T13:55:18.499Z'
docname: _source/pygromos.simulations.hpc_queuing.submission_systems
images: {}
path: /source-pygromos-simulations-hpc-queuing-submission-systems
title: pygromos.simulations.hpc_queuing.submission_systems package
---

# pygromos.simulations.hpc_queuing.submission_systems package

## Submodules

## pygromos.simulations.hpc_queuing.submission_systems.dummy module


### _class_ pygromos.simulations.hpc_queuing.submission_systems.dummy.DUMMY(verbose: bool = False, nomp: int = 1, nmpi: int = 1, job_duration: str = '24:00', submission: bool = True, environment=None)
Bases: `pygromos.simulations.hpc_queuing.submission_systems._submission_system._SubmissionSystem`

This SubmissionSystem is for testing submission systems. It basically prints everything out.


#### get_jobs_from_queue(job_text: str, \*\*kwargs)
search_queue_for_jobname

> this jobs searches the job queue for a certain job id.


* **Parameters**

    
    * **job_text** (*str*)


    * **regex** (*bool, optional*) – if the string is a Regular Expression



* **Returns**

    output contains all ids of fitting jobs to the querry



* **Return type**

    List[int]



#### job_queue_list(_: pandas.core.frame.DataFram_ )

#### search_queue_for_jobname(job_name: str, \*\*kwargs)
> this jobs searches the job queue for a certain job id.


* **Parameters**

    
    * **job_name** (*str*)


    * **regex** (*bool, optional*) – if the string is a Regular Expression



* **Returns**

    the output of the queue containing the jobname



* **Return type**

    List[str]



#### submit_jobAarray_to_queue(sub_job: pygromos.simulations.hpc_queuing.submission_systems.submission_job.Submission_job)
this function is submitting
:Parameters: **sub_job** (*Submission_job*) – submission job parameters


* **Returns**

    if a job was submitted the jobID is returned else None.



* **Return type**

    int, None



#### submit_to_queue(sub_job: pygromos.simulations.hpc_queuing.submission_systems.submission_job.Submission_job)
This function submits a str command to the submission system.


* **Parameters**

    **sub_job** (*Submission_job*) – submission job parameters



* **Returns**

    if a job was submitted the jobID is returned else None.



* **Return type**

    int, None



#### verbose(_: boo_ )
## pygromos.simulations.hpc_queuing.submission_systems.local module


### _class_ pygromos.simulations.hpc_queuing.submission_systems.local.LOCAL(submission: bool = True, nomp: int = 1, nmpi: int = 1, job_duration: str = '24:00', verbose: bool = False, environment=None, zip_trajectories: bool = True)
Bases: `pygromos.simulations.hpc_queuing.submission_systems._submission_system._SubmissionSystem`

This class handles local submission without a queueing system


#### job_queue_list(_: pandas.core.frame.DataFram_ )

#### search_queue_for_jobid(job_id: int, \*\*kwargs)
> this jobs searches the job queue for a certain job id.
> DUMMY FUNCTION!


* **Parameters**

    **job_id** (*int*) – id of the job



* **Raises**

    **NotImplemented** – Needs to be implemented in subclasses



#### search_queue_for_jobname(job_name: str, \*\*kwargs)
> this jobs searches the job queue for a certain job name.
> DUMMY FUNCTION!


* **Parameters**

    
    * **job_name** (*str*)


    * **regex** (*bool, optional*) – if the string is a Regular Expression



* **Returns**

    the output of the queue containing the jobname



* **Return type**

    List[str]



#### submit_jobAarray_to_queue(sub_job: pygromos.simulations.hpc_queuing.submission_systems.submission_job.Submission_job)
submitt a local job array


#### submit_to_queue(sub_job: pygromos.simulations.hpc_queuing.submission_systems.submission_job.Submission_job)
submitt a local job


#### verbose(_: boo_ )
## pygromos.simulations.hpc_queuing.submission_systems.lsf module


### _class_ pygromos.simulations.hpc_queuing.submission_systems.lsf.LSF(submission: bool = True, nomp: int = 1, nmpi: int = 1, job_duration: str = '24:00', max_storage: float = 1000, verbose: bool = False, environment=None, block_double_submission: bool = True, bjobs_only_same_host: bool = False, chain_prefix: str = 'done', begin_mail: bool = False, end_mail: bool = False, zip_trajectories: bool = True)
Bases: `pygromos.simulations.hpc_queuing.submission_systems._submission_system._SubmissionSystem`

This class is a wrapper for the LSF queueing system by IBM, like it is used on Euler.


#### get_queued_jobs()
> This function updates the job-list of the queueing system in the class.


* **Returns**

    returns the job_queue as pandas dataFrame.



* **Return type**

    pd.DataFrame



#### kill_jobs(job_name: Optional[str] = None, regex: bool = False, job_ids: Optional[Union[int, List[int]]] = None)
> this function can be used to terminate or remove pending jobs from the queue.


* **Parameters**

    
    * **job_name** (*str*) – name of the job to be killed


    * **regex** (*bool*) – if true, all jobs matching job_name get killed!


    * **job_ids** (*Union[List[int], int]*) – job Ids to be killed



#### search_queue_for_jobid(job_id: int)
> this jobs searches the job queue for a certain job id.


* **Parameters**

    **job_id** (*int*) – id of the job



* **Raises**

    **NotImplementedError** – Needs to be implemented in subclasses



#### search_queue_for_jobname(job_name: str, regex: bool = False)
> this jobs searches the job queue for a certain job name.


* **Parameters**

    
    * **job_name** (*str*)


    * **regex** (*bool, optional*) – if the string is a Regular Expression



* **Returns**

    the output of the queue containing the jobname



* **Return type**

    List[str]



#### submit_jobAarray_to_queue(sub_job: pygromos.simulations.hpc_queuing.submission_systems.submission_job.Submission_job)
This functioncan be used for submission of a job array. The ammount of jobs is determined by  the difference:

    end_job-start_job

An array index variable is defined called ${JOBID} inside the command representing job x in the array.


* **Parameters**

    **sub_job** (*Submission_job*) – the job to be submitted



* **Returns**

    return job ID



* **Return type**

    int



#### submit_to_queue(sub_job: pygromos.simulations.hpc_queuing.submission_systems.submission_job.Submission_job)
> This function submits the given command to the LSF QUEUE


* **Parameters**

    
    * **submission_job** (*Submission_job*) – the job to be submitted


    * **——-**


## pygromos.simulations.hpc_queuing.submission_systems.submission_job module


### _class_ pygromos.simulations.hpc_queuing.submission_systems.submission_job.Submission_job(command: Optional[str] = None, jobName: Optional[str] = None, outLog: Optional[str] = None, errLog: Optional[str] = None, start_job: Optional[int] = None, end_job: Optional[int] = None, jobLim: Optional[int] = None, queue_after_jobID: Optional[int] = None, post_execution_command: Optional[str] = None, submit_from_dir: Optional[str] = None, sumbit_from_file: bool = True, jobGroup: Optional[str] = None, jobID=None)
Bases: `object`

Description:
This class stores parameters for the submission of jobs. It is used by the submission_systems:
- submission dummy
- submission local
- submission lsf

It should handle all the information required for a single job, while the submission_systems handles more permanent settings.

It should provide an easy way to modify jobs, even from high level modules (e.g. the simulation module like TI, EMIN, …).

Author: Marc Lehner


#### _property_ command(_: st_ )

#### _property_ end_job(_: in_ )

#### _property_ errLog(_: st_ )

#### _property_ jobGroup(_: st_ )

#### _property_ jobID(_: in_ )

#### _property_ jobLim(_: in_ )

#### _property_ jobName(_: st_ )

#### _property_ outLog(_: st_ )

#### _property_ post_execution_command(_: st_ )

#### _property_ queue_after_jobID(_: in_ )

#### _property_ start_job(_: in_ )

#### _property_ submit_from_dir(_: st_ )

#### _property_ sumbit_from_file(_: boo_ )
## Module contents


### pygromos.simulations.hpc_queuing.submission_systems.get_submission_system(testing: bool = False)
