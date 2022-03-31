import os
import warnings
import pandas as pd
from pygromos.utils.typing import List

from pygromos.simulations.hpc_queuing.submission_systems._submission_system import _SubmissionSystem
from pygromos.simulations.hpc_queuing.submission_systems.submission_job import Submission_job

from pygromos.utils import bash


class LOCAL(_SubmissionSystem):
    """
    This class handles local submission without a queueing system
    """

    def __init__(
        self,
        submission: bool = True,
        nomp: int = 1,
        nmpi: int = 1,
        job_duration: str = "24:00",
        verbose: bool = False,
        environment=None,
        zip_trajectories: bool = True,
    ):
        super().__init__(
            verbose=verbose,
            nmpi=nmpi,
            nomp=nomp,
            job_duration=job_duration,
            submission=submission,
            environment=environment,
            zip_trajectories=zip_trajectories,
        )

    def submit_to_queue(self, sub_job: Submission_job) -> int:
        """
        submitt a local job
        """

        orig_dir = os.getcwd()

        if isinstance(sub_job.submit_from_dir, str) and os.path.isdir(sub_job.submit_from_dir):
            os.chdir(sub_job.submit_from_dir)
            command_file_path = sub_job.submit_from_dir + "/job_" + str(sub_job.jobName) + ".sh"
        else:
            command_file_path = "./job_" + str(sub_job.jobName) + ".sh"

        sub_job.command = sub_job.command.strip()  # remove trailing linebreaks

        if self._nomp >= 1:
            command = "export OMP_NUM_THREADS=" + str(self._nomp) + ";\n " + sub_job.command + ""
        else:
            command = sub_job.command

        if sub_job.sumbit_from_file:
            command_file = open(command_file_path, "w")
            command_file.write("#!/bin/bash\n")
            command_file.write(command.replace("&& ", ";\n") + ";\n")
            command_file.close()
            command = command_file_path
            bash.execute("chmod +x " + command_file_path, env=self.environment)

        # finalize string

        if self.verbose:
            print("Submission Command: \t", " ".join(command))
        if self.submission:
            try:
                process = bash.execute(command=command, catch_STD=True, env=self.environment)
                std_out_buff = map(str, process.stdout.readlines())
                std_out = "\t" + "\n\t".join(std_out_buff)

                # next sopt_job is queued with id:
                if self.verbose:
                    print("STDOUT: \n\t" + std_out + "\nEND\n")
                if os.path.exists(orig_dir):
                    os.chdir(orig_dir)

                return 0
            except ChildProcessError:
                try:
                    print(process)
                except ChildProcessError:
                    pass
                raise ChildProcessError("command failed: \n" + str(command))
        else:
            print("Did not submit: ", command)
            return -1

    def submit_jobAarray_to_queue(self, sub_job: Submission_job) -> int:
        """
        submitt a local job array
        """

        # generate submission_string:
        submission_string = ""

        if isinstance(sub_job.submit_from_dir, str) and os.path.isdir(sub_job.submit_from_dir):
            submission_string += "cd " + sub_job.submit_from_dir + " && "

        if self._nomp > 1:
            command = submission_string + " export OMP_NUM_THREADS=" + str(self._nomp) + " && " + sub_job.command
        else:
            command = submission_string + sub_job.command

        if self.verbose:
            print("Submission Command: \t", " ".join(submission_string))
        if self.submission:
            try:
                for jobID in range(sub_job.start_job, sub_job.end_job + 1):
                    std_out_buff = bash.execute(
                        command="export JOBID=" + str(jobID) + " && " + command, env=self.environment
                    )
                    std_out = "\n".join(std_out_buff.readlines())
                    if self.verbose:
                        print("sdtout : " + str(std_out))
                return 0
            except ChildProcessError:
                raise ChildProcessError("could not submit this command: \n" + submission_string)
        else:
            print("Did note submit: ", command)
            return -1

    def search_queue_for_jobname(self, job_name: str, **kwargs) -> List[str]:
        """search_queue_for_jobname

            this jobs searches the job queue for a certain job name.
            DUMMY FUNCTION!

        Parameters
        ----------
        job_name :  str
        regex:  bool, optional
            if the string is a Regular Expression
        Returns
        -------
        List[str]
            the output of the queue containing the jobname
        """
        if self.verbose:
            print("Searching job Name: ", job_name)
            warnings.warn("Queue search was called, but no queue present!")
        return []

    def search_queue_for_jobid(self, job_id: int, **kwargs) -> pd.DataFrame:
        """search_queue_for_jobid

            this jobs searches the job queue for a certain job id.
            DUMMY FUNCTION!

        Parameters
        ----------
        job_id :  int
            id of the job
        Raises
        -------
        NotImplemented
            Needs to be implemented in subclasses
        """
        if self.verbose:
            print("Searching job ID: ", job_id)
            warnings.warn("Queue search was called, but no queue present!")
        return []
