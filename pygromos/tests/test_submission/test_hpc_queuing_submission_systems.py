import unittest
from pygromos.hpc_queuing.submission_systems._submission_system import _SubmissionSystem
from pygromos.hpc_queuing.submission_systems.lsf import LSF
from pygromos.hpc_queuing.submission_systems.local import LOCAL
from pygromos.hpc_queuing.submission_systems.dummy import DUMMY


class test_queuing_system(unittest.TestCase):
    file_class = _SubmissionSystem
    verbose = True
    submission = True

    def test_construct(self):
        subSys = self.file_class()


class test_DUMMY(test_queuing_system):
    file_class = DUMMY

    def setUp(self) -> None:
        self.SubmissionSystem = self.file_class(verbose=self.verbose, submission=self.submission)

    def test_submit(self):
        command = "echo \" WUHAHAHA\""

        self.SubmissionSystem.submit_to_queue(command=command, verbose=self.verbose)

    def test_submit_jobAarray_to_queue(self):
        command = "echo \" WUHAHAHA\""
        start_job = 1
        end_job = 1

        self.SubmissionSystem.submit_jobAarray_to_queue(command=command, start_Job=start_job, end_job=end_job, verbose=self.verbose)

        start_job = 1
        end_job = 10

        self.SubmissionSystem.submit_jobAarray_to_queue(command=command, start_Job=start_job, end_job=end_job, verbose=self.verbose)

    def test_search_queue_for_jobname(self):
        search_job = "Test"

        self.SubmissionSystem.search_queue_for_jobname(job_name=search_job, verbose=self.verbose)

    def test_get_jobs_from_queue(self):
        get_jobs_with = "Test"

        self.SubmissionSystem.get_jobs_from_queue(job_text=get_jobs_with, verbose=self.verbose)


class test_LOCAL(test_queuing_system):
    file_class = LOCAL

    def setUp(self) -> None:
        self.SubmissionSystem = self.file_class(verbose=self.verbose, submission=False)

    def test_submit(self):
        command = "echo \" WUHAHAHA\""
        self.SubmissionSystem.submit_to_queue(command=command)

    def test_submit_jobAarray_to_queue(self):
        command = "echo \" WUHAHAHA\""
        start_job = 1
        end_job = 1

        self.SubmissionSystem.submit_jobAarray_to_queue(command=command, start_Job=start_job, end_job=end_job)

        start_job = 1
        end_job = 10

        self.SubmissionSystem.submit_jobAarray_to_queue(command=command, start_Job=start_job, end_job=end_job)

    def test_search_queue_for_jobname(self):
        search_job = "Test"

        self.SubmissionSystem.search_queue_for_jobname(job_name=search_job, verbose=self.verbose)

    def test_get_jobs_from_queue(self):
        get_jobs_with = "Test"

        self.SubmissionSystem.get_jobs_from_queue(job_text=get_jobs_with, verbose=self.verbose)


class test_LSF(test_queuing_system):
    file_class = LSF

    def setUp(self) -> None:
        self.SubmissionSystem = self.file_class(verbose=self.verbose, submission=False)

    def test_submit(self):
        command = "echo \" WUHAHAHA\""

        self.SubmissionSystem.submit_to_queue(command=command, jobName="TEST", sumbit_from_file=False,
                                              dummyTesting=True)

    def test_submit_jobAarray_to_queue(self):
        command = "echo \" WUHAHAHA\""
        start_job = 1
        end_job = 1
        self.SubmissionSystem.submit_jobAarray_to_queue(command=command, start_Job=start_job, end_job=end_job, jobName="TEST",)

        start_job = 1
        end_job = 10
        self.SubmissionSystem.submit_jobAarray_to_queue(command=command, start_Job=start_job, end_job=end_job, jobName="TEST",)

    def test_search_queue_for_jobname(self):
        search_job = "Test"
        self.SubmissionSystem.search_queue_for_jobname(job_name=search_job,dummyTesting=True)

    def test_get_jobs_from_queue(self):
        get_jobs_with = "Test"

        self.SubmissionSystem.get_jobs_from_queue(job_text=get_jobs_with, dummyTesting=True)
