---
date: '2022-04-01T13:55:18.499Z'
docname: _source/pygromos.simulations.modules
images: {}
path: /source-pygromos-simulations-modules
title: pygromos.simulations.modules package
---

# pygromos.simulations.modules package

## Submodules

## pygromos.simulations.modules.general_analysis_modules module

ideas:

    
    * Automatic Equilibration detection


    * property analysis


    * default analysis for emin ? eq?

Start with local submission system first. :)

## pygromos.simulations.modules.general_simulation_modules module


### pygromos.simulations.modules.general_simulation_modules.simulation(in_gromos_simulation_system: pygromos.files.gromos_system.gromos_system.Gromos_System, override_project_dir: typing.Optional[str] = None, step_name: str = 'sim', in_imd_path: typing.Optional[str] = None, submission_system: pygromos.simulations.hpc_queuing.submission_systems._submission_system._SubmissionSystem = <pygromos.simulations.hpc_queuing.submission_systems.local.LOCAL object>, simulation_runs: int = 1, equilibration_runs: int = 0, previous_simulation_run: typing.Optional[int] = None, force_simulation: bool = False, initialize_first_run=False, reinitialize_every_run=False, analysis_script: callable = <function do>, analysis_control_dict: typing.Optional[dict] = None, _no_double_submit_check: bool = False, _work_dir: typing.Optional[str] = None, verbose: bool = True, verbose_lvl: int = 1, _template_imd_path: typing.Optional[str] = None)
> This function is a generic simulation block, that can be used to run and schedule simulations.


* **Parameters**

    
    * **in_gromos_simulation_system** (*Gromos_System*) – gromos system that contains the information of the files, that are required for the simulation run.


    * **override_project_dir** (*str*) – parent project directory


    * **step_name** (*str, optional*) – name of the step - used as jobname and folder in the project dir.


    * **in_imd_path** (*str, optional*) – if a path to an imd file is given here, it will be used!


    * **submission_system** (*_SubmissionSystem, optional*) – the system, that should be used to submit a job


    * **simulation_runs** (*int, optional*) – number of sequential simulation runs


    * **equilibration_runs** (*int, optional*) – number of equilibration runs, that will not appear in the final traj


    * **previous_simulation_run** (*int, optional*) – job-ID of the previous simulation


    * **force_simulation** (*bool, optional*) – if simulation already exist, shall it be overwritten?


    * **analysis_script** (*callable, optional*) – script that is used to analyse the job


    * **analysis_control_dict** (*dict, optional*) – sub selection of steps, for the analysis script, that should be executed


    * **verbose** (*bool, optional*) – “baeh, baeh, baeh”  - a wise black-nose sheep from ausserberg.


    * **verbose_lvl** (*int, optional*) – how much baeh?


    * **_template_imd_path** (*str, optional*) – default template imd



* **Returns**

    returns a new gromos system which is containing the simulation info



* **Return type**

    [Gromos_System](#pygromos.files.gromos_system.gromos_system.Gromos_System),


## pygromos.simulations.modules.preset_simulation_modules module


### pygromos.simulations.modules.preset_simulation_modules.emin(in_gromos_system: pygromos.files.gromos_system.gromos_system.Gromos_System, step_name: str = 'emin', override_project_dir: typing.Optional[str] = None, in_imd_path=None, submission_system: pygromos.simulations.hpc_queuing.submission_systems._submission_system._SubmissionSystem = <pygromos.simulations.hpc_queuing.submission_systems.local.LOCAL object>, simulation_runs: int = 1, equilibration_runs: int = 0, previous_simulation_run: typing.Optional[int] = None, _template_imd_path: str = '/home/bschroed/Documents/projects/PyGromosTools/pygromos/data/simulation_parameters_templates/emin.imd', _no_double_submit_check: bool = False, initialize_first_run=False, analysis_script: callable = <function do>, verbose: bool = True)

### pygromos.simulations.modules.preset_simulation_modules.md(in_gromos_system: pygromos.files.gromos_system.gromos_system.Gromos_System, step_name: str = 'md', override_project_dir: typing.Optional[str] = None, in_imd_path=None, submission_system: pygromos.simulations.hpc_queuing.submission_systems._submission_system._SubmissionSystem = <pygromos.simulations.hpc_queuing.submission_systems.local.LOCAL object>, simulation_runs: int = 1, equilibration_runs: int = 0, initialize_first_run=False, previous_simulation_run: typing.Optional[int] = None, _template_imd_path: str = '/home/bschroed/Documents/projects/PyGromosTools/pygromos/data/simulation_parameters_templates/md.imd', _no_double_submit_check: bool = False, analysis_script: callable = <function do>, verbose: bool = True)

### pygromos.simulations.modules.preset_simulation_modules.sd(in_gromos_system: pygromos.files.gromos_system.gromos_system.Gromos_System, step_name: str = 'sd', override_project_dir: typing.Optional[str] = None, in_imd_path=None, submission_system: pygromos.simulations.hpc_queuing.submission_systems._submission_system._SubmissionSystem = <pygromos.simulations.hpc_queuing.submission_systems.local.LOCAL object>, simulation_runs: int = 1, equilibration_runs: int = 0, initialize_first_run=False, previous_simulation_run: typing.Optional[int] = None, _template_imd_path: str = '/home/bschroed/Documents/projects/PyGromosTools/pygromos/data/simulation_parameters_templates/vacuum_sd.imd', _no_double_submit_check: bool = False, analysis_script: callable = <function do>, verbose: bool = True)
## pygromos.simulations.modules.ti_modules module

Free Energy


### pygromos.simulations.modules.ti_modules.TI_sampling(in_gromos_system: pygromos.files.gromos_system.gromos_system.Gromos_System, project_dir: str, step_name='lambda_sampling', lambda_values: typing.List[float] = array([0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1. ]), subSystem: pygromos.simulations.hpc_queuing.submission_systems._submission_system._SubmissionSystem = <pygromos.simulations.hpc_queuing.submission_systems.local.LOCAL object>, n_productions: int = 3, n_equilibrations: int = 1, randomize: bool = False, dual_cnf: typing.Optional[typing.List[str]] = None, verbose: bool = True)
This function will automatically submit N independent (different lambda)
MD simulations with a lambda dependent potential energy.


* **Parameters**

    
    * **in_gromos_system** (*Gromos_System*) – input gromos system


    * **project_dir** (*str*) – directory in which simulation input files are found


    * **step_name** (*str*) – subdirectory of project_dir, in which we will write the output
    important: allows to run multiple random seeds with a different “step_name”


    * **lambda_values** (*List [float]*) – List of lambda values for each independent simulation


    * **subSystem** (*_SubmissionSystem*) – where will the calculation run


    * **n_productions** (*int*) – number of chunks each independent simulation is broken down into


    * **n_equilibrations** (*int*) – number of chunks of equilibration preceding the production for each independent simulation


    * **randomize** (*bool*) – Choose a random number for the initial random seed (same for all lambda values)


    * **dual_cnf** (*List [str], optional*) – If provided, should be the path to two conformations (matching end states A and B) which
    can be used as initial conformations
    Simulations with a lambda value between 0 and 0.5 will use the first as starting conformation



* **Returns**

    **lam_system** – Gromos system of the simulation submitted last



* **Return type**

    [Gromos_System](#pygromos.files.gromos_system.gromos_system.Gromos_System)


## Module contents
