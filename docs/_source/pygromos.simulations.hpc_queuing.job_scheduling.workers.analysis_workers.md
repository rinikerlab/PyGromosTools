---
date: '2022-04-01T13:55:18.499Z'
docname: _source/pygromos.simulations.hpc_queuing.job_scheduling.workers.analysis_workers
images: {}
path: /source-pygromos-simulations-hpc-queuing-job-scheduling-workers-analysis-workers
title: pygromos.simulations.hpc_queuing.job_scheduling.workers.analysis_workers package
---

# pygromos.simulations.hpc_queuing.job_scheduling.workers.analysis_workers package

## Submodules

## pygromos.simulations.hpc_queuing.job_scheduling.workers.analysis_workers.simulation_analysis module


### pygromos.simulations.hpc_queuing.job_scheduling.workers.analysis_workers.simulation_analysis.do(in_simulation_dir: str, out_analysis_dir: str, sim_prefix: str, n_processes: int = 1, control_dict: Optional[dict] = None, verbose: bool = True)
> This function is a analysis framework structure, that starts an analysis folder containing a data folder with all concatenated files, from which analysis can be started.


* **Parameters**

    
    * **in_simulation_dir** (*str*) – input simulation directory (with succesfully finished simulations)


    * **out_analysis_dir** (*str*) – output directory


    * **sim_prefix** (*str*) – prefix of the simulation == name of simulation


    * **n_processes** (*int, optional*) – WARNING: parallelization is currently not implemented!, by default 1


    * **control_dict** (*dict, optional*) – control structure, steering the executions, by default None


    * **verbose** (*bool, optional*) – bla bla, by default True



### pygromos.simulations.hpc_queuing.job_scheduling.workers.analysis_workers.simulation_analysis.gather_simulation_step_file_paths(in_folder: str, filePrefix: str = '', fileSuffixes: Union[str, List[str]] = ['.tre', '.tre.tar.gz'], verbose: bool = False, finalNumberingSort=False)
gather_replica_file_paths

> Finds all trajectory paths in a simulation folder.


* **Parameters**

    
    * **in_folder** (*str*) – folder, containing the files


    * **replicas** (*int*) – Number of replicas


    * **filePrefix** (*str, optional*) – str prefix the desired files


    * **fileSuffixes** (*str, optional*) – str suffix of the desired files


    * **verbose** (*bool*) – toggle verbosity



### pygromos.simulations.hpc_queuing.job_scheduling.workers.analysis_workers.simulation_analysis.project_concatenation(in_folder: str, out_folder: str, in_prefix: str, control_dict: Dict[str, bool], verbose: bool = False)
concatenation of the simulation data.


* **Parameters**

    
    * **in_folder** (*str*) – folder containing the simulation results


    * **out_folder** (*str*) – folder that should contain the concatenated out files


    * **in_prefix** (*str*) – prefix of the simulation files.


    * **control_dict** (*Dict[str, bool]*) – control of what should be executed


    * **verbose** (*bool, optional*) – baeeeeh baeeeeh, by default False



* **Returns**

    resulting cnf path.



* **Return type**

    str


## Module contents
