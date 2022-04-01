---
date: '2022-04-01T13:55:18.499Z'
docname: _source/pygromos.simulations.hpc_queuing.job_scheduling.workers.simulation_workers
images: {}
path: /source-pygromos-simulations-hpc-queuing-job-scheduling-workers-simulation-workers
title: pygromos.simulations.hpc_queuing.job_scheduling.workers.simulation_workers
  package
---

# pygromos.simulations.hpc_queuing.job_scheduling.workers.simulation_workers package

## Submodules

## pygromos.simulations.hpc_queuing.job_scheduling.workers.simulation_workers.clean_up_simulation_files module

After simulation cleanup
This module is thought to be used as a clean up script after each run, reducing the need of storage for RE-EDS simulations on Euler or any other cluster :)
it mainly removes not needed temporary files and compresses long trajectory files.
It should be hanged in after each simulation step.


### pygromos.simulations.hpc_queuing.job_scheduling.workers.simulation_workers.clean_up_simulation_files.do(in_simulation_dir: str, n_processes: int = 1, verbose: bool = True)
do clean_up_simulation_files


* **Parameters**

    
    * **in_simulation_dir** (*str*) – path to folder, which contains the raw gromos output


    * **n_processes** (*int, optional*) – how many processes should be used?


    * **verbose** (*bool, optional*) – loud and noisy?



* **Return type**

    None


## pygromos.simulations.hpc_queuing.job_scheduling.workers.simulation_workers.simulation_run_worker module


### pygromos.simulations.hpc_queuing.job_scheduling.workers.simulation_workers.simulation_run_worker.work(out_dir: str, in_cnf_path: str, in_imd_path: str, in_top_path: str, runID: int = 1, in_perttopo_path: Optional[str] = None, in_disres_path: Optional[str] = None, in_posres_path: Optional[str] = None, in_refpos_path: Optional[str] = None, in_qmmm_path: Optional[str] = None, out_trc: bool = False, out_tre: bool = False, out_trg: bool = False, out_trv: bool = False, out_trf: bool = False, out_trs: bool = False, nmpi: int = 1, nomp: int = 1, reinitialize_every_run: bool = False, initialize_first_run: bool = True, gromosXX_bin_dir: Optional[str] = None, gromosXX_check_binary_paths: bool = True, work_dir: Optional[str] = None, zip_trajectories: bool = True, \*\*kwargs)
Executed by repex_EDS_long_production_run as workers


* **Parameters**

    
    * **out_dir** (*str*) – final output dir


    * **in_cnf_path** (*str*) – input coordinates


    * **in_imd_path** (*str*) – input imd-parameter file


    * **in_top_path** (*str*) – input topology


    * **in_perttopo_path** (*str*) – input pertubation


    * **in_disres_path** (*str*) – input disres


    * **in_qmmm_path** (*str*) – input qmmm


    * **nmpi** (*int, optional*) – number of mpi cores (def.=1)


    * **nomp** (*int, optional*) – number of omp cores (def.= 1)


    * **out_trg** (*bool, optional*) – True if trg shall be written out.


    * **out_trv** (*bool, optional*) – True if trv shall be written out.


    * **out_trf** (*bool, optional*) – True if trf shall be written out.


    * **out_trs** (*bool, optional*) – True if trs shall be written out.


    * **gromosXX_bin_dir** (*str, optional*) – path to gromos binary directory


    * **work_dir** (*str, optional*) – work directory


    * **zip_trajectories** (*bool*) – determines whether trajectories are zipped



* **Returns**

    return number



* **Return type**

    int


## Module contents
