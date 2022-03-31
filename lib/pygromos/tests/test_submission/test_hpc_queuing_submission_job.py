import unittest
from pygromos.simulations.hpc_queuing.submission_systems.submission_job import Submission_job


class test_queuing_system(unittest.TestCase):
    file_class = Submission_job

    def test_construct_min(self):
        sub_job = self.file_class(command='echo " WUHAHAHA"')
        assert isinstance(sub_job, self.file_class)

    def test_construct_full(self):
        sub_job = self.file_class(
            command='echo " WUHAHAHA"',
            jobName="test_job",
            start_job=1,
            end_job=10,
            jobLim=20,
            outLog="/dev/null/out.log",
            errLog="/dev/null/err.log",
            queue_after_jobID=1,
            post_execution_command='echo " WUHAHAHA2"',
            submit_from_dir=".",
            sumbit_from_file=True,
        )
        assert isinstance(sub_job, self.file_class)

    def test_acces(self):
        sub_job = self.file_class(command='echo " WUHAHAHA"', queue_after_jobID=1)
        id = sub_job.queue_after_jobID
        assert isinstance(id, int)
        assert id == 1
        assert sub_job.command == 'echo " WUHAHAHA"'
