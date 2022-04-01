---
date: '2022-04-01T13:55:18.499Z'
docname: _source/pygromos.tests.test_submission
images: {}
path: /source-pygromos-tests-test-submission
title: pygromos.tests.test_submission package
---

# pygromos.tests.test_submission package

## Submodules

## pygromos.tests.test_submission.test_hpc_queuing_submission_job module


### _class_ pygromos.tests.test_submission.test_hpc_queuing_submission_job.test_queuing_system(methodName='runTest')
Bases: `unittest.case.TestCase`


#### file_class()
alias of [`pygromos.simulations.hpc_queuing.submission_systems.submission_job.Submission_job`](#pygromos.simulations.hpc_queuing.submission_systems.submission_job.Submission_job)


#### test_acces()

#### test_construct_full()

#### test_construct_min()
## pygromos.tests.test_submission.test_hpc_queuing_submission_scheduling module


### _class_ pygromos.tests.test_submission.test_hpc_queuing_submission_scheduling.test_MD_scheduler(methodName='runTest')
Bases: `unittest.case.TestCase`


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### submissionSystem()
alias of [`pygromos.simulations.hpc_queuing.submission_systems.dummy.DUMMY`](#pygromos.simulations.hpc_queuing.submission_systems.dummy.DUMMY)


#### test_do()
## pygromos.tests.test_submission.test_hpc_queuing_submission_systems module


### _class_ pygromos.tests.test_submission.test_hpc_queuing_submission_systems.test_DUMMY(methodName='runTest')
Bases: `pygromos.tests.test_submission.test_hpc_queuing_submission_systems.test_queuing_system`


#### file_class()
alias of [`pygromos.simulations.hpc_queuing.submission_systems.dummy.DUMMY`](#pygromos.simulations.hpc_queuing.submission_systems.dummy.DUMMY)


### _class_ pygromos.tests.test_submission.test_hpc_queuing_submission_systems.test_LOCAL(methodName='runTest')
Bases: `pygromos.tests.test_submission.test_hpc_queuing_submission_systems.test_queuing_system`


#### file_class()
alias of [`pygromos.simulations.hpc_queuing.submission_systems.local.LOCAL`](#pygromos.simulations.hpc_queuing.submission_systems.local.LOCAL)


#### submission(_ = Fals_ )

### _class_ pygromos.tests.test_submission.test_hpc_queuing_submission_systems.test_LSF(methodName='runTest')
Bases: `pygromos.tests.test_submission.test_hpc_queuing_submission_systems.test_queuing_system`


#### file_class()
alias of [`pygromos.simulations.hpc_queuing.submission_systems.lsf.LSF`](#pygromos.simulations.hpc_queuing.submission_systems.lsf.LSF)


#### submission(_ = Fals_ )

### _class_ pygromos.tests.test_submission.test_hpc_queuing_submission_systems.test_queuing_system(methodName='runTest')
Bases: `unittest.case.TestCase`


#### file_class()
alias of `pygromos.simulations.hpc_queuing.submission_systems._submission_system._SubmissionSystem`


#### submission(_ = Tru_ )

#### test_construct()

#### test_get_jobs_from_queue()

#### test_search_queue_for_jobname()

#### test_submit()

#### test_submit_jobAarray_to_queue()

#### test_submit_jobAarray_to_queue1_10()

#### verbose(_ = Fals_ )
## Module contents
