import os
import re
import socket
import warnings
from typing import List, Union

from pygromos.utils import bash


def get_submission_system(testing: bool = False):
    if (testing):
        return DUMMY
    if ("eu" in socket.gethostname()):
        return LSF
    else:
        return LOCAL


class _SubmissionSystem:
    verbose: bool
    submission: bool

    _job_duration: str = "24:00"
    _nmpi: int
    _nomp: int
    _max_storage: float

    def __init__(self, submission: bool = False,
                 nmpi: int = 1, nomp: int = 1, max_storage: float = 1000, job_duration: str = "24:00",
                 verbose: bool = False, enviroment=None):
        self.verbose = verbose
        self.submission = submission

        self._job_duration = job_duration
        self._nmpi = nmpi
        self._nomp = nomp
        self._max_storage = max_storage
        self._enviroment = enviroment

    def submit_to_queue(self, **kargs) -> Union[int, None]:
        raise NotImplemented("Do is not implemented for: " + self.__class__.__name__)

    def get_script_generation_command(self, var_name: str = None, var_prefixes: str = "") -> str:
        name = self.__class__.__name__
        if (var_name is None):
            var_name = var_prefixes + name

        gen_cmd = "#Generate " + name + "\n"
        gen_cmd += "from " + self.__module__ + " import " + name + " as " + name + "_obj" + "\n"
        gen_cmd += var_name + " = " + name + "_obj(submission=" + str(self.submission) + ", verbose=" + str(
            self.verbose) + ")\n\n"
        return gen_cmd

    def get_jobs_from_queue(self, job_text: str, **kwargs) -> List[int]:
        """get_jobs_from_queue

            this function searches the job queue for a certain job id.

        Parameters
        ----------
        job_text :  str
            text part of the job_line
        regex:  bool, optional
            if the string is a Regular Expression
        Raises
        -------
        NotImplemented
            Needs to be implemented in subclasses
        """

        raise NotImplemented("Do is not implemented for: " + self.__class__.__name__)

    def search_queue_for_jobname(self, job_name: str, **kwargs) -> List[str]:
        """search_queue_for_jobname

            this jobs searches the job queue for a certain job id.

        Parameters
        ----------
        job_name :  str
            name of the job
        regex:  bool, optional
            if the string is a Regular Expression
        Raises
        -------
        NotImplemented
            Needs to be implemented in subclasses
        """

        raise NotImplemented("Do is not implemented for: " + self.__class__.__name__)

    @property
    def nmpi(self) -> int:
        return self._nmpi

    @nmpi.setter
    def nmpi(self, nmpi: int):
        self._nmpi = int(nmpi)

    @property
    def nomp(self) -> int:
        return self._nomp

    @nomp.setter
    def nomp(self, nomp: int):
        self._nomp = int(nomp)

    @property
    def job_duration(self) -> str:
        return self._job_duration

    @job_duration.setter
    def job_duration(self, job_duration: str):
        self._job_duration = str(job_duration)

    @property
    def max_storage(self) -> str:
        return self._max_storage

    @max_storage.setter
    def max_storage(self, max_storage: float):
        self._max_storage = float(max_storage)


class LSF(_SubmissionSystem):
    """LSF
        This class is a wrapper for the LSF queueing system by IBM, like it is used on Euler.
    """

    def __init__(self, submission: bool = True, nomp: int = 1, nmpi: int = 1, job_duration: str = "24:00",
                 verbose: bool = False, enviroment=None):
        super().__init__(verbose=verbose, nmpi=nmpi, nomp=nomp, job_duration=job_duration, submission=submission, enviroment=enviroment)

    def submit_to_queue(self, command: str,
                        jobName: str, outLog=None, errLog=None,
                        queue_after_jobID: int = None, do_not_doubly_submit_to_queue: bool = False,
                        force_queue_start_after: bool = False,
                        projectName: str = None, jobGroup: str = None, priority=None,
                        begin_mail: bool = False, end_mail: bool = False,
                        post_execution_command: str = None,
                        submit_from_dir: str = None, sumbit_from_file: bool = True,
                        dummyTesting: bool = False, verbose: bool = None) -> Union[int, None]:
        """
            This function submits the given command to the LSF QUEUE

        Parameters
        ----------
        command : str
            command to be executed
        jobName : str
            name of the job in the queue
        outLog: str, optional
            out std-out log path
        errLog: str, optional
            out std-err log path
        submit_from_dir
        queue_after_jobID: int, optional
            shall this job be queued after another one?
        force_queue_start_after: bool, optional
            shall this job start after another job, no matter the exit state?
        projectName :  NOT IMPLEMENTED AT THE MOMENT
        jobGroup :  NOT IMPLEMENTED AT THE MOMENT
        priority :  NOT IMPLEMENTED AT THE MOMENT
        begin_mail :    bool, optional
            send a mail when job starts
        end_mail :  bool, optional
            send a mail, when job is finished
        post_execution_command: str, optional
            command which will be executed after the execution of command
        do_not_doubly_submit_to_queue:  bool, optional
            if True: script checks the submission queue and looks for an identical job_name. than raises ValueError if it is already submitted. (default: True)
        verbose:    bool, optional
            WARNING! - Will be removed in Future! use attribute verbose or constructor! WARNING!
            print out some messages
        dummyTesting:   bool, optional
            WARNING! - Will be removed in Future! use attribute or constructor as submission! WARNING!
            do not submit the job to the queue.
        stupid_mode
        Returns
        -------

        """
        # job_properties:Job_properties=None, <- currently not usd
        orig_dir = os.getcwd()

        # required parsing will be removed in future:
        if (not verbose is None):
            self.verbose = verbose
        if (dummyTesting):
            self.submission = not dummyTesting

        # generate submission_string:
        submission_string = ""

        # QUEUE checking to not double submit
        if (do_not_doubly_submit_to_queue and self.submission):
            if (self.verbose): print('check queue')
            ids = self.get_jobs_from_queue(jobName)

            if (len(ids) > 0):
                if (self.verbose): print(
                    "\tSKIP - FOUND JOB: \t\t" + "\n\t\t".join(map(str, ids)) + "\n\t\t with jobname: " + jobName)
                return ids[0]

        if (isinstance(submit_from_dir, str) and os.path.isdir(submit_from_dir)):
            os.chdir(submit_from_dir)
            command_file_path = submit_from_dir + "/job_" + str(jobName) + ".sh"
        else:
            command_file_path = "./job_" + str(jobName) + ".sh"

        submission_string += "bsub "
        submission_string += " -J" + jobName + " "
        submission_string += " -W " + str(self.job_duration) + " "

        if (not isinstance(post_execution_command, type(None))):
            submission_string += "-Ep \"" + post_execution_command + "\" "

        if (not isinstance(outLog, str) and not isinstance(errLog, str)):
            outLog = jobName + ".out"
            submission_string += " -o " + outLog
        elif (isinstance(outLog, str)):
            submission_string += " -o " + outLog

        if (isinstance(errLog, str)):
            submission_string += " -e " + errLog

        nCPU = self.nmpi * self.nomp
        submission_string += " -n " + str(nCPU) + " "
        add_string = ""
        # add_string= "-R \"select[model==XeonGold_5118 || model==XeonGold_6150 || model==XeonE3_1585Lv5 || model==XeonE3_1284Lv4 || model==XeonE7_8867v3 || model == XeonGold_6140 || model==XeonGold_6150 ]\""
        if (isinstance(self.max_storage, int)):
            submission_string += " -R rusage[mem=" + str(self.max_storage) + "] "

        if (isinstance(queue_after_jobID, (int, str))):
            prefix = "done"
            if (force_queue_start_after):
                prefix = "ended"
            submission_string += " -w \"" + prefix + "(" + str(queue_after_jobID) + ")\" "

        if (begin_mail):
            submission_string += " -B "
        if (end_mail):
            submission_string += " -N "

        if (self.nomp > 1):
            command = "\"export OMP_NUM_THREADS=" + str(self.nomp) + ";\n " + command + "\""
        else:
            command = "\n " + command.strip() + ""

        if (sumbit_from_file):
            if (self.verbose): print("writing tmp-submission-file to: ", command_file_path)
            command_file = open(command_file_path, "w")
            command_file.write("#!/bin/bash\n")
            command_file.write(command + ";\n")
            command_file.close()
            command = command_file_path

            bash.execute("chmod +x " + command_file_path, env=self._enviroment)

        ##finalize string
        submission_string = list(map(lambda x: x.strip(), submission_string.split())) + [command]

        if (self.verbose): print("Submission Command: \t", " ".join(submission_string))
        if (self.submission):
            try:
                out_process = bash.execute(command=submission_string, catch_STD=True, env=self._enviroment)
                std_out = "\n".join(map(str, out_process.stdout.readlines()))

                # next sopt_job is queued with id:
                id_start = std_out.find("<")
                id_end = std_out.find(">")
                job_id = str(std_out[id_start + 1:id_end]).strip()
                if self.verbose: print("process returned id: " + str(job_id))
                if (job_id == "" and job_id.isalnum()):
                    raise ValueError("Did not get at job ID!")
            except:
                raise ChildProcessError("could not submit this command: \n" +
                                        str(submission_string))
        else:
            job_id = -1

        os.chdir(orig_dir)
        return int(job_id)

    def submit_jobAarray_to_queue(self, command: str, jobName: str,
                                  start_Job: int, end_job: int, jobLim: int = None,
                                  outLog=None, errLog=None, submit_from_dir: str = None,
                                  queue_after_jobID: int = None, force_queue_start_after: bool = False,
                                  do_not_doubly_submit_to_queue: bool = True,
                                  jobGroup: str = None,
                                  begin_mail: bool = False, end_mail: bool = False,
                                  verbose: bool = None, dummyTesting: bool = False, ) -> Union[int, None]:
        """
        This functioncan be used for submission of a job array. The ammount of jobs is determined by  the difference:
                    end_Job-start_Job
        An array index variable is defined called ${JOBID} inside the command representing job x in the array.

        Parameters
        ----------
        command : str
            command to be executed
        jobName : str
            name of the job in the queue
        start_Job: int
            starting job_id
        end_job: int
            ending job_id
        jobLim: int, optional
            limits the in parallel execute of job arrays.
        duration: str, optional
            this string defines the max job-run duration like HHH:MM (default: 04:00)
        outLog: str, optional
            out std-out log path
        errLog: str, optional
            out std-err log path
        submit_from_dir
        nmpi :  int, optional
            integer number of mpi cores (default: 1)
        nomp :  int, optional
            integer number of omp cores (default: 1)
        maxStorage : int, optional
            Max memory per core.
        queue_after_jobID: int, optional
            shall this job be queued after another one?
        force_queue_start_after: bool, optional
            shall this job start after another job, no matter the exit state?
        projectName :  NOT IMPLEMENTED AT THE MOMENT
        jobGroup :  NOT IMPLEMENTED AT THE MOMENT
        priority :  NOT IMPLEMENTED AT THE MOMENT
        begin_mail :    bool, optional
            send a mail when job starts
        end_mail :  bool, optional
            send a mail, when job is finished
        verbose:    bool, optional
            WARNING! - Will be removed in Future! use attribute verbose or constructor! WARNING!
            print out some messages
        dummyTesting:   bool, optional
            WARNING! - Will be removed in Future! use attribute or constructor as submission! WARNING!
            do not submit the job to the queue.
        do_not_doubly_submit_to_queue:  bool, optional - NOT IMPLEMENTED AT THE MOMENT
            if True: script checks the submission queue and looks for an identical job_name. than raises ValueError if it is already submitted. (default: True)


        Returns
        -------
         Union[int, None]
            return job ID

        Raises
        ------
        ValueError
            if job already submitted a Value Error is raised
        ChildProcessError
            if submission to queue via bash fails an Child Process Error is raised.
        """

        # required parsing will be removed in future:
        if (not verbose is None):
            self.verbose = verbose
        if (dummyTesting):
            self.submission = not dummyTesting

        # QUEUE checking to not double submit
        if (self.submission and do_not_doubly_submit_to_queue):
            if (self.verbose): print('check queue')
            ids = self.get_jobs_from_queue(jobName)

            if (len(ids) > 0):
                if (self.verbose): print(
                    "\tSKIP - FOUND JOB: \t\t" + "\n\t\t".join(map(str, ids)) + "\n\t\t with jobname: " + jobName)
                return ids[0]

        # generate submission_string:
        submission_string = ""
        if (isinstance(submit_from_dir, str) and os.path.isdir(submit_from_dir)):
            submission_string += "cd " + submit_from_dir + " && "

        if (jobLim is None):
            jobLim = end_job - start_Job

        jobName = str(jobName) + "[" + str(start_Job) + "-" + str(end_job) + "]%" + str(jobLim)

        submission_string += "bsub -J \" " + jobName + " \" -W \"" + str(self.job_duration) + "\" "

        if (isinstance(jobGroup, str)):
            submission_string += " -g " + jobGroup + " "

        if (not isinstance(outLog, str) and not isinstance(errLog, str)):
            outLog = jobName + ".out"
            submission_string += " -oo " + outLog
        elif (isinstance(outLog, str)):
            submission_string += " -oo " + outLog

        if (isinstance(errLog, str)):
            submission_string += " -eo " + errLog

        nCPU = self.nmpi * self.nomp
        submission_string += " -n " + str(nCPU) + " "

        if (isinstance(self.max_storage, int)):
            submission_string += " -R \"rusage[mem=" + str(self.max_storage) + "]\" "

        if (isinstance(queue_after_jobID, (int, str))):
            prefix = "\"done"
            if (force_queue_start_after):
                prefix = "\"ended"
            submission_string += " -w " + prefix + "(" + str(queue_after_jobID) + ")\" "

        if (begin_mail):
            submission_string += " -B "
        if (end_mail):
            submission_string += " -N "

        if (self.nomp > 1):
            command = " \" export OMP_NUM_THREADS=" + str(self.nomp) + " && " + command + "\""
        else:
            command = " \"" + command + "\""

        ##finalize string
        submission_string = list(map(lambda x: x.strip(), submission_string.split())) + [command]

        if (self.verbose): print("Submission Command: \t", " ".join(submission_string))
        if (self.submission):
            try:
                std_out_buff = bash.execute(command=submission_string, env=self._enviroment)
                std_out = "\n".join(std_out_buff.readlines())

                # next sopt_job is queued with id:
                id_start = std_out.find("<")
                id_end = std_out.find(">")
                job_id = str(std_out[id_start + 1:id_end]).strip()
                if self.verbose: print("process returned id: " + str(job_id))
                if (job_id == "" and job_id.isalnum()):
                    raise ValueError("Did not get at job ID!")
            except:
                raise ChildProcessError("could not submit this command: \n" + " ".join(submission_string))
        else:
            job_id = 0
        return int(job_id)

    def get_jobs_from_queue(self, job_text: str, regex: bool = False,
                            dummyTesting: bool = False, verbose: bool = False) -> List[int]:
        """search_queue_for_jobname

            this jobs searches the job queue for a certain job id.

        Parameters
        ----------
        job_text :  str
        regex:  bool, optional
            if the string is a Regular Expression
        Returns
        -------
        List[int]
            output contains all ids of fitting jobs to the querry
        """
        out_job_lines = self.search_queue_for_jobname(job_text, regex=regex,
                                                      dummyTesting=dummyTesting, verbose=verbose)
        if (verbose): print("Isolate job IDs")
        get_job_ids = list(map(int, filter(lambda x: x.isdigit(), map(lambda x: x.split(" ")[0], out_job_lines))))
        if (verbose): print("jobID: ", get_job_ids)
        return get_job_ids

    def search_queue_for_jobname(self, job_name: str, regex: bool = False,
                                 dummyTesting: bool = False, verbose: bool = False) -> List[str]:
        """search_queue_for_jobname

            this jobs searches the job queue for a certain job id.

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
        try:
            if (not regex):
                job_name = " " + job_name + " "
            if (verbose): print("getting job data")
            if (dummyTesting):
                out_job_lines = ["123 TEST", "456 TEST2"]
            else:
                out_process = bash.execute("bjobs -w", catch_STD=True, env=self._enviroment)
                stdout = list(map(str, out_process.stdout.readlines()))
                out_job_lines = list(filter(lambda line: re.search(job_name, line), stdout))

            if (verbose): print("got job_lines: ", out_job_lines)
        except TimeoutError:
            return []
        return out_job_lines


class LOCAL(_SubmissionSystem):
    """LSF
        This class is a wrapper for the LSF queueing system by IBM, like it is used on Euler.
    """

    def __init__(self, submission: bool = True, nomp: int = 1, nmpi: int = 1, job_duration: str = "24:00",
                 verbose: bool = False, enviroment=None):
        super().__init__(verbose=verbose, nmpi=nmpi, nomp=nomp, job_duration=job_duration, submission=submission, enviroment=enviroment)

    def submit_to_queue(self, command: str, jobName: str = "local", submit_from_dir: str = None,
                        sumbit_from_file: bool = False, **kwargs) -> Union[int, None]:
        # job_properties:Job_properties=None, <- currently not usd
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

    def get_jobs_from_queue(self, job_text: str, **kwargs) -> List[int]:
        """search_queue_for_jobname

            this jobs searches the job queue for a certain job id.

        Parameters
        ----------
        job_text :  str
        regex:  bool, optional
            if the string is a Regular Expression
        Returns
        -------
        List[int]
            output contains all ids of fitting jobs to the querry
        """
        if (self.verbose): print("Searching ID: ", job_text)
        warnings.warn("Queue search was called, but no queue present!")
        return []

    def search_queue_for_jobname(self, job_name: str, **kwargs) -> List[str]:
        """search_queue_for_jobname

            this jobs searches the job queue for a certain job id.

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
        if (self.verbose): print("Searching job Name: ", job_name)
        warnings.warn("Queue search was called, but no queue present!")
        return []


class DUMMY(_SubmissionSystem):
    """LSF
        This class is a wrapper for the LSF queueing system by IBM, like it is used on Euler.
    """

    def __init__(self, verbose: bool = False, nomp: int = 1, nmpi: int = 1, job_duration: str = "24:00",
                 submission: bool = True, enviroment=None):
        super().__init__(verbose=verbose, nmpi=nmpi, nomp=nomp, job_duration=job_duration, submission=submission, enviroment=enviroment)

    def submit_to_queue(self, command: str, **kwargs) -> Union[int, None]:
        if (self.submission):
            print("\n", command, "\n")
            return 0
        else:
            print('did not submit')
            return None

    def submit_jobAarray_to_queue(self, command: str, start_Job: int, end_job: int,
                                  **kwargs) -> Union[int, None]:
        if (self.submission):
            print()
            for jobID in range(start_Job, end_job + 1):
                print("Job " + str(jobID) + ":", command, "\n")
            print()
            return 0
        else:
            print('did not submit')
            return None

    def get_jobs_from_queue(self, job_text: str, **kwargs) -> List[int]:
        """search_queue_for_jobname

            this jobs searches the job queue for a certain job id.

        Parameters
        ----------
        job_text :  str
        regex:  bool, optional
            if the string is a Regular Expression
        Returns
        -------
        List[int]
            output contains all ids of fitting jobs to the querry
        """
        if (self.verbose): print("Retrieving jobs from list with: ", job_text)
        warnings.warn("Queue search was called, but no queue present!")
        return []

    def search_queue_for_jobname(self, job_name: str, **kwargs) -> List[str]:
        """search_queue_for_jobname

            this jobs searches the job queue for a certain job id.

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
        if (self.verbose): print("Searching job Name: ", job_name)
        warnings.warn("Queue search was called, but no queue present!")
        return []
