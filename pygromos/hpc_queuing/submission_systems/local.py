import os
import warnings
import pandas as pd
from typing import Union, List

from pygromos.hpc_queuing.submission_systems._submission_system import _SubmissionSystem
from pygromos.utils import bash


class LOCAL(_SubmissionSystem):
    """LSF
        This class is a wrapper for the LSF queueing system by IBM, like it is used on Euler.
    """

    def __init__(self, submission: bool = True, nomp: int = 1, nmpi: int = 1, job_duration: str = "24:00",
                 verbose: bool = False, enviroment=None):
        super().__init__(verbose=verbose, nmpi=nmpi, nomp=nomp, job_duration=job_duration, submission=submission, enviroment=enviroment)

    def submit_to_queue(self, command: str, jobName: str = "local", submit_from_dir: str = None,
                        sumbit_from_file: bool = False, **kwargs) -> Union[int, None]:
        """submit_to_queue

                This function submits a job to the submission queue.

        Parameters
        ----------
        command :   str
        submit_from_dir:    str, optional
            from which dir shall the jobs be submitted? (def. current)
        verbose:    bool, optional
            print out some messages

        Returns
        -------
        int
            return job ID

        Raises
        ------
        ValueError
            if job already submitted a Value Error is raised
        ChildProcessError
            if submission to queue via bash fails an Child Process Error is raised.
        """
        orig_dir = os.getcwd()

        if (isinstance(submit_from_dir, str) and os.path.isdir(submit_from_dir)):
            os.chdir(submit_from_dir)
            command_file_path = submit_from_dir + "/job_" + str(jobName) + ".sh"
        else:
            command_file_path = "./job_" + str(jobName) + ".sh"

        if (self.nomp > 1):
            command = "export OMP_NUM_THREADS=" + str(self.nomp) + ";\n " + command + ""
        else:
            command = command.strip()

        if (sumbit_from_file):
            command_file = open(command_file_path, "w")
            command_file.write("#!/bin/bash\n")
            command_file.write(command.replace("&& ", ";\n") + ";\n")
            command_file.close()
            command = command_file_path
            bash.execute("chmod +x " + command_file_path, env=self._enviroment)

        ##finalize string

        if (self.verbose): print("Submission Command: \t", " ".join(command))
        if (self.submission):
            try:
                process = bash.execute(command=command, catch_STD=True, env=self._enviroment)
                std_out_buff = map(str, process.stdout.readlines())
                std_out = "\t" + "\n\t".join(std_out_buff)

                # next sopt_job is queued with id:
                if self.verbose: print("STDOUT: \n\t" + std_out + "\nEND\n")
                if (os.path.exists(orig_dir)): os.chdir(orig_dir)

                return 0
            except:
                try:
                    print(process)
                except:
                    pass
                raise ChildProcessError("command failed: \n" +
                                        str(command))
        else:
            print("Did not submit: ", command)
            return None

    def submit_jobAarray_to_queue(self, command: str, start_Job: int, end_job: int, submit_from_dir: str = None,
                                  **kwargs) -> Union[int, None]:
        # job_properties:Job_properties=None, <- currently not usd
        """submit_to_queue
        NOT IMPLEMENTED YET!
                This function submits a job to the submission queue.

        Parameters
        ----------
        command :   str
        submit_from_dir:    str, optional
            from which dir shall the jobs be submitted? (def. current)
        verbose:    bool, optional
            print out some messages

        Returns
        -------
        int
            return job ID

        Raises
        ------
        ValueError
            if job already submitted a Value Error is raised
        ChildProcessError
            if submission to queue via bash fails an Child Process Error is raised.
        """
        # generate submission_string:
        submission_string = ""

        if (isinstance(submit_from_dir, str) and os.path.isdir(submit_from_dir)):
            submission_string += "cd " + submit_from_dir + " && "

        if (self.nomp > 1):
            command = submission_string + " export OMP_NUM_THREADS=" + str(self.nomp) + " && " + command
        else:
            command = submission_string + command

        if (self.verbose): print("Submission Command: \t", " ".join(submission_string))
        if (self.submission):
            try:
                for jobID in range(start_Job, end_job + 1):
                    std_out_buff = bash.execute(command="export JOBID=" + str(jobID) + " && " + command, env=self._enviroment)
                    std_out = "\n".join(std_out_buff.readlines())
                    if self.verbose: print("sdtout : " + str(std_out))
                return 0
            except:
                raise ChildProcessError("could not submit this command: \n" + submission_string)
        else:
            print("Did note submit: ", command)
            return None


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
        if (self.verbose):
            print("Searching job Name: ", job_name)
            warnings.warn("Queue search was called, but no queue present!")
        return []


    def search_queue_for_jobid(self, job_id: int, **kwargs)->pd.DataFrame:
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
        if (self.verbose):
            print("Searching job ID: ", job_id)
            warnings.warn("Queue search was called, but no queue present!")
        return []