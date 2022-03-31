import unittest
from pygromos.simulations.hpc_queuing.submission_systems._submission_system import _SubmissionSystem
from pygromos.simulations.hpc_queuing.submission_systems.submission_job import Submission_job

from pygromos.simulations.hpc_queuing.submission_systems.lsf import LSF
from pygromos.simulations.hpc_queuing.submission_systems.local import LOCAL
from pygromos.simulations.hpc_queuing.submission_systems.dummy import DUMMY


class test_queuing_system(unittest.TestCase):
    file_class = _SubmissionSystem
    verbose = False
    submission = True

    def test_construct(self):
        subSys = self.file_class(verbose=self.verbose, submission=self.submission)
        assert isinstance(subSys, self.file_class)

    def test_submit(self):
        sub_job = Submission_job(jobName="test_job", command='echo " WUHAHAHA"')
        subSys = self.file_class(verbose=self.verbose, submission=self.submission)
        subSys.submit_to_queue(sub_job=sub_job)

    def test_submit_jobAarray_to_queue(self):

        sub_job = Submission_job(jobName="test_job", command='echo " WUHAHAHA"', start_job=1, end_job=1)

        subSys = self.file_class(verbose=self.verbose, submission=self.submission)
        subSys.submit_jobAarray_to_queue(sub_job=sub_job)

    def test_submit_jobAarray_to_queue1_10(self):
        sub_job2 = Submission_job(jobName="test_job", command='echo " WUHAHAHA"', start_job=1, end_job=10)
        subSys = self.file_class(verbose=self.verbose, submission=self.submission)
        subSys.submit_jobAarray_to_queue(sub_job=sub_job2)

    def test_search_queue_for_jobname(self):
        if self.file_class == _SubmissionSystem:
            return
        search_job = "Test"
        subSys = self.file_class(verbose=self.verbose, submission=self.submission)
        subSys.search_queue_for_jobname(job_name=search_job)

    def test_get_jobs_from_queue(self):
        if self.file_class == _SubmissionSystem:
            return
        get_jobs_with = "Test"
        subSys = self.file_class(verbose=self.verbose, submission=self.submission)
        subSys.get_jobs_from_queue(job_text=get_jobs_with)


class test_DUMMY(test_queuing_system):
    file_class = DUMMY


class test_LOCAL(test_queuing_system):
    submission = False
    file_class = LOCAL


class test_LSF(test_queuing_system):
    submission = False
    file_class = LSF
