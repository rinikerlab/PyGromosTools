---
date: '2022-04-01T13:55:18.499Z'
docname: _source/pygromos.simulations.hpc_queuing.job_scheduling
images: {}
path: /source-pygromos-simulations-hpc-queuing-job-scheduling
title: pygromos.simulations.hpc_queuing.job_scheduling package
---

# pygromos.simulations.hpc_queuing.job_scheduling package

## Subpackages


* [pygromos.simulations.hpc_queuing.job_scheduling.schedulers package]()


    * [Submodules](#submodules)


    * [pygromos.simulations.hpc_queuing.job_scheduling.schedulers.scheduler_functions module](#module-pygromos.simulations.hpc_queuing.job_scheduling.schedulers.scheduler_functions)


    * [pygromos.simulations.hpc_queuing.job_scheduling.schedulers.simulation_scheduler module](#module-pygromos.simulations.hpc_queuing.job_scheduling.schedulers.simulation_scheduler)


    * [Module contents](#module-pygromos.simulations.hpc_queuing.job_scheduling.schedulers)


* [pygromos.simulations.hpc_queuing.job_scheduling.workers package]()


    * [Subpackages](#subpackages)


        * [pygromos.simulations.hpc_queuing.job_scheduling.workers.analysis_workers package]()


            * [Submodules](#submodules)


            * [pygromos.simulations.hpc_queuing.job_scheduling.workers.analysis_workers.simulation_analysis module](#module-pygromos.simulations.hpc_queuing.job_scheduling.workers.analysis_workers.simulation_analysis)


            * [Module contents](#module-pygromos.simulations.hpc_queuing.job_scheduling.workers.analysis_workers)


        * [pygromos.simulations.hpc_queuing.job_scheduling.workers.simulation_workers package]()


            * [Submodules](#submodules)


            * [pygromos.simulations.hpc_queuing.job_scheduling.workers.simulation_workers.clean_up_simulation_files module](#module-pygromos.simulations.hpc_queuing.job_scheduling.workers.simulation_workers.clean_up_simulation_files)


            * [pygromos.simulations.hpc_queuing.job_scheduling.workers.simulation_workers.simulation_run_worker module](#module-pygromos.simulations.hpc_queuing.job_scheduling.workers.simulation_workers.simulation_run_worker)


            * [Module contents](#module-pygromos.simulations.hpc_queuing.job_scheduling.workers.simulation_workers)


    * [Module contents](#module-pygromos.simulations.hpc_queuing.job_scheduling.workers)


## Submodules

## pygromos.simulations.hpc_queuing.job_scheduling.file_management module

This module is doing all the post simulation file juggeling needed for gromos. CURRENTLY OLD DON”T USE


### pygromos.simulations.hpc_queuing.job_scheduling.file_management._thread_worker_cat_trc(job: int, replicaID_range: List[int], trc_files: Dict[int, List[str]], out_prefix: str, topology_path: str, out_trcs: dict, dt: float, time: float = 0, verbose: bool = False, boundary_conditions: str = 'r cog', include_all: bool = False)
This thread worker_scripts concatenates all .trc files of one replica into one file.


* **Parameters**

    
    * **job** (*rank of this thread*)


    * **replicaID_range** (*x_range - list of all*)


    * **trc_files** (*Dict[int, List[str]]*) – Dictionary containing all replicas, with list of all trc files concerning one trc.


    * **out_prefix** (*str*) – output prefix


    * **verbose** (*bool*) – verbosity?



* **Return type**

    None



### pygromos.simulations.hpc_queuing.job_scheduling.file_management.compress_files(in_paths: List[str], n_processes: int = 1)
compress a list of files


* **Parameters**

    
    * **in_paths** (*List[str]*)


    * **n_processes** (*int*) – how many processes can be used in parallel?



* **Returns**

    outpaths



* **Return type**

    List[str]



### pygromos.simulations.hpc_queuing.job_scheduling.file_management.find_and_unarchive_tar_files(trc_files: List[str], verbose: bool = False)

### pygromos.simulations.hpc_queuing.job_scheduling.file_management.find_header(path: str)

### pygromos.simulations.hpc_queuing.job_scheduling.file_management.gather_simulation_file_paths(in_folder: str, filePrefix: str = '', fileSuffixes: Union[str, List[str]] = ['.tre', '.tre.tar.gz'], files_per_folder: int = 1, verbose: bool = False)

### pygromos.simulations.hpc_queuing.job_scheduling.file_management.gather_simulation_replica_file_paths(in_folder: str, replicas: int, filePrefix: str = '', fileSuffixes: Union[str, List[str]] = ['.tre', '.tre.tar.gz'], verbose: bool = False, finalNumberingSort=False)
gather_replica_file_paths

> Finds all trajectory paths in a simulation folder and sorts them by replica.


* **Parameters**

    
    * **in_folder** (*str*) – folder, containing the files


    * **replicas** (*int*) – Number of replicas


    * **filePrefix** (*str, optional*) – str prefix the desired files


    * **fileSuffixes** (*str, optional*) – str suffix of the desired files


    * **verbose** (*bool*) – toggle verbosity



### pygromos.simulations.hpc_queuing.job_scheduling.file_management.parse_csv_energy_trajectories(in_folder: str, ene_trajs_prefix: str, verbose: bool = False)
> searches a directory and loads energy eds csvs as pandas dataframes.


* **Parameters**

    
    * **in_folder** (*str*) – folder with energy_traj - csvs


    * **ene_trajs_prefix** (*str*) – prefix name


    * **verbose** (*bool*) – loud?



* **Returns**

    return a list with pandas data frames containing all energy infos.



* **Return type**

    List[pd.DataFrame]



### pygromos.simulations.hpc_queuing.job_scheduling.file_management.parse_csv_energy_trajectory(in_ene_traj_path: str, verbose: bool = False)
> parse_one ene_ana csv


* **Parameters**

    
    * **in_ene_traj_path** (*str*) – path to input file


    * **verbose** (*bool*) – loud?



* **Returns**

    return a pandas data frame containing all energies



* **Return type**

    pd.DataFrame



### pygromos.simulations.hpc_queuing.job_scheduling.file_management.project_concatenation(in_folder: str, in_topology_path: str, in_imd: str, num_replicas: int, control_dict: Dict[str, bool], out_folder: str, in_ene_ana_lib_path: str, out_file_prefix: str = 'test', fit_traj_to_mol: int = 1, starting_time: float = 0, include_water_in_trc=True, additional_properties: Union[Tuple[str], List[str]] = ('solvtemp2', 'totdisres'), n_processes: int = 1, gromosPP_bin_dir: Optional[str] = None, verbose: bool = False, nofinal=False, boundary_conditions: str = 'r cog')

### pygromos.simulations.hpc_queuing.job_scheduling.file_management.thread_worker_concat_repdat(job: int, repdat_file_out_path: str, repdat_file_paths: Union[str, List[str]], verbose: bool = False)

### pygromos.simulations.hpc_queuing.job_scheduling.file_management.thread_worker_isolate_energies(in_en_file_paths: str, out_folder: str, properties: List[str], replicas: List[int], in_ene_ana_lib: str, gromosPP_path: str, out_prefix: str = '', tre_prefix: str = '', time=None, dt=None, job: int = - 1, verbose=True)
isolate_properties_from_tre

    This func uses Ene Ana from gromos to isolate potentials from out_tre Files
    in in_folder generated by reeds.


* **Parameters**

    
    * **in_en_file_paths** (*str*) – path, in which the input tre_folders are situated.


    * **out_folder** (*str*) – output folder, where to write the energy .csvs


    * **properties** (*List[str]*) – potentials to isolate from the .out_tre Files


    * **replicas** (*int*) – number of replicas, that should be found


    * **in_ene_ana_lib** (*str*) – path to the ene_ana lib, encoding the out_tre Files


    * **gromosPP_path** (*str*) – path to the ene_ana lib, encoding the out_tre Files


    * **out_prefix** (*str, optional*)


    * **tre_prefix** (*str, optional*)


    * **verbose** (*bool, optional*)



* **Returns**

    return list of result Files.



* **Return type**

    List[str]


## pygromos.simulations.hpc_queuing.job_scheduling.module_functions module


### pygromos.simulations.hpc_queuing.job_scheduling.module_functions.write_job_script(out_script_path: str, target_function: callable, variable_dict: dict, python_cmd: str = 'python3', verbose: bool = False)
## Module contents
