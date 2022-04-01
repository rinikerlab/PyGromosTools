---
date: '2022-04-01T13:55:18.499Z'
docname: _source/pygromos.tests.test_simulation_blocks
images: {}
path: /source-pygromos-tests-test-simulation-blocks
title: pygromos.tests.test_simulation_blocks package
---

# pygromos.tests.test_simulation_blocks package

## Submodules

## pygromos.tests.test_simulation_blocks.test_simulation_runner_blocks module


### _class_ pygromos.tests.test_simulation_blocks.test_simulation_runner_blocks.test_simulation_blocks(methodName='runTest')
Bases: `unittest.case.TestCase`


#### input_cnf_path(_ = '/home/bschroed/Documents/projects/PyGromosTools/pygromos/tests/in_testfiles/small_system/6J29.cnf_ )

#### input_top_path(_ = '/home/bschroed/Documents/projects/PyGromosTools/pygromos/tests/in_testfiles/small_system/6J29.top_ )

#### setUp()
Hook method for setting up the test fixture before exercising it.


#### sim_block(step_name: str = 'emin', override_project_dir: typing.Optional[str] = None, in_imd_path=None, submission_system: pygromos.simulations.hpc_queuing.submission_systems._submission_system._SubmissionSystem = <pygromos.simulations.hpc_queuing.submission_systems.local.LOCAL object>, simulation_runs: int = 1, equilibration_runs: int = 0, previous_simulation_run: typing.Optional[int] = None, _template_imd_path: str = '/home/bschroed/Documents/projects/PyGromosTools/pygromos/data/simulation_parameters_templates/emin.imd', _no_double_submit_check: bool = False, initialize_first_run=False, analysis_script: callable = <function do>, verbose: bool = True)

#### submissionSystem(_ = <pygromos.simulations.hpc_queuing.submission_systems.dummy.DUMMY object_ )

#### test_emin()

#### test_lam_window()

#### test_md()

#### test_sd()

#### test_smulation()

#### test_ti_sampling()

#### verbose(_ = Fals_ )
## Module contents
