import os
from datetime import datetime
import pandas as pd
from pygromos.utils.typing import Union, List

from pygromos.simulations.hpc_queuing.submission_systems._submission_system import _SubmissionSystem
from pygromos.simulations.hpc_queuing.submission_systems.submission_job import Submission_job

from pygromos.utils import bash


class LSF(_SubmissionSystem):
    """LSF
    This class is a wrapper for the LSF queueing system by IBM, like it is used on Euler.
    """

    _dummy: bool = False
    _refresh_job_queue_list_all_s: int = 60  # update the job-queue list every x seconds
    _job_queue_time_stamp: datetime

    def __init__(
        self,
        submission: bool = True,
        nomp: int = 1,
        nmpi: int = 1,
        job_duration: str = "24:00",
        max_storage: float = 1000,
        verbose: bool = False,
        environment=None,
        block_double_submission: bool = True,
        bjobs_only_same_host: bool = False,
        chain_prefix: str = "done",
        begin_mail: bool = False,
        end_mail: bool = False,
        zip_trajectories: bool = True,
    ):
        # general settings for the submission system
        super().__init__(
            verbose=verbose,
            nmpi=nmpi,
            nomp=nomp,
            job_duration=job_duration,
            max_storage=max_storage,
            submission=submission,
            environment=environment,
            block_double_submission=block_double_submission,
            chain_prefix=chain_prefix,
            begin_mail=begin_mail,
            end_mail=end_mail,
            zip_trajectories=zip_trajectories,
        )
        # Only LSF specific settings:
        self.bjobs_only_same_host = bjobs_only_same_host

    def submit_to_queue(self, sub_job: Submission_job) -> int:
        """
            This function submits the given command to the LSF QUEUE

        Parameters
        ----------
        submission_job : Submission_job
            the job to be submitted
        -------

        """
        # job_properties:Job_properties=None, <- currently not usd
        orig_dir = os.getcwd()

        # generate submission_string:
        submission_string = ""

        # QUEUE checking to not double submit
        if self._block_double_submission and self._submission:
            if self.verbose:
                print("check queue")
            ids = list(self.search_queue_for_jobname(sub_job.jobName).index)

            if len(ids) > 0:
                if self.verbose:
                    print(
                        "\tSKIP - FOUND JOB: \t\t"
                        + "\n\t\t".join(map(str, ids))
                        + "\n\t\t with jobname: "
                        + sub_job.jobName
                    )
                return ids[0]

        if isinstance(sub_job.submit_from_dir, str) and os.path.isdir(sub_job.submit_from_dir):
            os.chdir(sub_job.submit_from_dir)
            command_file_path = sub_job.submit_from_dir + "/job_" + str(sub_job.jobName) + ".sh"
        else:
            command_file_path = "./job_" + str(sub_job.jobName) + ".sh"

        submission_string += "bsub "
        submission_string += " -J" + sub_job.jobName + " "
        submission_string += " -W " + str(self._job_duration) + " "

        if not isinstance(sub_job.post_execution_command, type(None)):
            submission_string += '-Ep "' + sub_job.post_execution_command + '" '

        if not isinstance(sub_job.outLog, str) and not isinstance(sub_job.errLog, str):
            outLog = sub_job.jobName + ".out"
            submission_string += " -o " + outLog
        elif isinstance(sub_job.outLog, str):
            submission_string += " -o " + sub_job.outLog

        if isinstance(sub_job.errLog, str):
            submission_string += " -e " + sub_job.errLog

        nCPU = self._nmpi * self._nomp
        submission_string += " -n " + str(nCPU) + " "
        # TODO: add GPU support
        # add_string = ""
        # add_string= "-R \"select[model==XeonGold_5118 || model==XeonGold_6150 || model==XeonE3_1585Lv5 || model==XeonE3_1284Lv4 || model==XeonE7_8867v3 || model == XeonGold_6140 || model==XeonGold_6150 ]\""
        if isinstance(self._max_storage, int):
            submission_string += " -R rusage[mem=" + str(self._max_storage) + "] "

        if isinstance(sub_job.queue_after_jobID, (int, str)) and (
            sub_job.queue_after_jobID != 0 or sub_job.queue_after_jobID != "0"
        ):
            submission_string += ' -w "' + self._chain_prefix + "(" + str(sub_job.queue_after_jobID) + ')" '

        if self._begin_mail:
            submission_string += " -B "
        if self._end_mail:
            submission_string += " -N "

        sub_job.command = sub_job.command.strip()  # remove trailing line breaks

        if self._nomp >= 1:
            command = "export OMP_NUM_THREADS=" + str(self._nomp) + ";\n " + sub_job.command + " "
        else:
            command = "\n " + sub_job.command + ""

        if sub_job.sumbit_from_file:
            if self.verbose:
                print("writing tmp-submission-file to: ", command_file_path)
            command_file = open(command_file_path, "w")
            command_file.write("#!/bin/bash\n")
            command_file.write(command + ";\n")
            command_file.close()
            command = command_file_path

            bash.execute("chmod +x " + command_file_path, env=self._environment)

        # finalize string
        submission_string = list(map(lambda x: x.strip(), submission_string.split())) + [command]

        if self.verbose:
            print("Submission Command: \t", " ".join(submission_string))
        if self._submission and not self._dummy:
            try:
                out_process = bash.execute(command=submission_string, catch_STD=True, env=self._environment)
                std_out = "\n".join(map(str, out_process.stdout.readlines()))

                # next sopt_job is queued with id:
                id_start = std_out.find("<")
                id_end = std_out.find(">")
                job_id = int(str(std_out[id_start + 1 : id_end]).strip())
                if self.verbose:
                    print("process returned id: " + str(job_id))
                if str(job_id) == "" and job_id.isalnum():
                    raise ValueError("Did not get at job ID!")
            except Exception as e:
                raise ChildProcessError("could not submit this command: \n" + str(submission_string) + "\n\n" + str(e))
        else:
            job_id = -1

        os.chdir(orig_dir)
        sub_job.jobID = job_id
        return job_id

    def submit_jobAarray_to_queue(self, sub_job: Submission_job) -> int:
        """
        This functioncan be used for submission of a job array. The ammount of jobs is determined by  the difference:
                    end_job-start_job
        An array index variable is defined called ${JOBID} inside the command representing job x in the array.

        Parameters
        ----------
        sub_job: Submission_job
            the job to be submitted

        Returns
        -------
         int
            return job ID

        """

        # QUEUE checking to not double submit
        if self._submission and self._block_double_submission:
            if self.verbose:
                print("check queue")
            ids = self.search_queue_for_jobname(sub_job.jobName)

            if len(ids) > 0:
                if self.verbose:
                    print(
                        "\tSKIP - FOUND JOB: \t\t"
                        + "\n\t\t".join(map(str, ids))
                        + "\n\t\t with jobname: "
                        + sub_job.jobName
                    )
                return ids[0]

        # generate submission_string:
        submission_string = ""
        if isinstance(sub_job.submit_from_dir, str) and os.path.isdir(sub_job.submit_from_dir):
            submission_string += "cd " + sub_job.submit_from_dir + " && "

        if sub_job.jobLim is None:
            jobLim = sub_job.end_job - sub_job.start_job

        jobName = str(sub_job.jobName) + "[" + str(sub_job.start_job) + "-" + str(sub_job.end_job) + "]%" + str(jobLim)

        submission_string += 'bsub -J " ' + jobName + ' " -W "' + str(self._job_duration) + '" '

        if isinstance(sub_job.jobGroup, str):
            submission_string += " -g " + sub_job.jobGroup + " "

        if not isinstance(sub_job.outLog, str) and not isinstance(sub_job.errLog, str):
            outLog = jobName + ".out"
            submission_string += " -oo " + outLog
        elif isinstance(sub_job.outLog, str):
            submission_string += " -oo " + sub_job.outLog

        if isinstance(sub_job.errLog, str):
            submission_string += " -eo " + sub_job.errLog

        nCPU = self._nmpi * self._nomp
        submission_string += " -n " + str(nCPU) + " "

        if isinstance(self.max_storage, int):
            submission_string += ' -R "rusage[mem=' + str(self._max_storage) + ']" '

        if isinstance(sub_job.queue_after_jobID, (int, str)):
            submission_string += " -w " + self._chain_prefix + "(" + str(sub_job.queue_after_jobID) + ')" '

        if self._begin_mail:
            submission_string += " -B "
        if self._end_mail:
            submission_string += " -N "

        if self._nomp > 1:
            command = " export OMP_NUM_THREADS=" + str(self._nomp) + " && " + sub_job.command + " "
        else:
            command = " " + sub_job.command + " "

        # finalize string
        submission_string = list(map(lambda x: x.strip(), submission_string.split())) + [command]

        if self.verbose:
            print("Submission Command: \t", " ".join(submission_string))
        if self._submission and not self._dummy:
            try:
                std_out_buff = bash.execute(command=submission_string, env=self._environment)
                std_out = "\n".join(std_out_buff.readlines())

                # next sopt_job is queued with id:
                id_start = std_out.find("<")
                id_end = std_out.find(">")
                job_id = str(std_out[id_start + 1 : id_end]).strip()
                if self.verbose:
                    print("process returned id: " + str(job_id))
                if job_id == "" and job_id.isalnum():
                    raise ValueError("Did not get at job ID!")
            except Exception as e:
                raise ChildProcessError(
                    "could not submit this command: \n" + " ".join(submission_string) + "\n\n" + str(e)
                )
        else:
            job_id = -1
        sub_job.jobID = job_id
        return int(job_id)

    """
        Job Queue Managment
    """

    def get_queued_jobs(self) -> pd.DataFrame:
        """
            This function updates the job-list of the queueing system in the class.

        Returns
        -------
        pd.DataFrame
            returns the job_queue as pandas dataFrame.
        """
        # Do we need an update of the job list?
        check_job_list = True
        if hasattr(self, "_job_queue_time_stamp"):
            last_update = datetime.now() - self._job_queue_time_stamp
            check_job_list = last_update.seconds > self._refresh_job_queue_list_all_s
        if not self._submission:  # shortcut to reduce queue calls!
            self._job_queue_list = pd.DataFrame(
                columns=["JOBID      USER    STAT  QUEUE      FROM_HOST   EXEC_HOST   JOB_NAME   SUBMIT_TIME".split()]
            )
            return self._job_queue_list
        if check_job_list:
            # try getting the lsf queue
            if not self._dummy:
                try:
                    # get all running and pending jobs
                    if self.bjobs_only_same_host:
                        out_process = bash.execute("bjobs -w | grep '$HOSTNAME|JOBID'", catch_STD=True)
                    else:
                        out_process = bash.execute("bjobs -w", catch_STD=True)
                    job_list_str = list(map(lambda x: x.decode("utf-8"), out_process.stdout.readlines()))

                    # get all finished jobs
                    if self.bjobs_only_same_host:
                        out_process = bash.execute("bjobs -wd | grep '$HOSTNAME|JOBID'", catch_STD=True)
                    else:
                        out_process = bash.execute("bjobs -wd", catch_STD=True)
                    job_list_finished_str = list(map(lambda x: x.decode("utf-8"), out_process.stdout.readlines()))
                    self._job_queue_time_stamp = datetime.now()
                except Exception as err:
                    raise Exception("Could not get job_list!\nerr:\n" + "\n".join(err.args))
            else:
                job_list_str = []
                job_list_finished_str = []

            # format information:
            jlist = list(map(lambda x: x.strip().split(), job_list_str))
            jlist_fin = list(map(lambda x: x.strip().split(), job_list_finished_str))
            if len(jlist) > 1:
                header = jlist[0]
                jobs = jlist[1:] + jlist_fin[1:]

                jobs_dict = {}
                for job in jobs:
                    jobID = int(job[0].split("[")[0])
                    user = job[1]
                    status = job[2]
                    queue = job[3]
                    from_host = job[4]
                    exec_host = job[5]
                    job_name = " ".join(job[6:-3])
                    submit_time = datetime.strptime(
                        str(datetime.now().year) + " " + " ".join(job[-3:]), "%Y %b %d %H:%M"
                    )
                    values = [jobID, user, status, queue, from_host, exec_host, job_name, submit_time]
                    jobs_dict.update({jobID: {key: value for key, value in zip(header, values)}})

                self._job_queue_list = pd.DataFrame(jobs_dict, index=None).T
            else:
                self._job_queue_list = pd.DataFrame(
                    columns=[
                        "JOBID      USER    STAT  QUEUE      FROM_HOST   EXEC_HOST   JOB_NAME   SUBMIT_TIME".split()
                    ]
                )
        else:
            if self.verbose:
                print("Skipping refresh of job list, as the last update is " + str(last_update) + "s ago")
        return self._job_queue_list

    def search_queue_for_jobid(self, job_id: int) -> pd.DataFrame:
        self.get_queued_jobs()
        return self._job_queue_list.where(self._job_queue_list.JOBID == job_id).dropna()

    def search_queue_for_jobname(self, job_name: str, regex: bool = False) -> pd.DataFrame:
        """search_queue_for_jobname

            this jobs searches the job queue for a certain job name.

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

        self.get_queued_jobs()
        if regex:
            return self._job_queue_list.where(self._job_queue_list.JOB_NAME.str.match(job_name)).dropna()
        else:
            return self._job_queue_list.where(self._job_queue_list.JOB_NAME == job_name).dropna()

    """
        kill jobs
    """

    def kill_jobs(self, job_name: str = None, regex: bool = False, job_ids: Union[List[int], int] = None):
        """
            this function can be used to terminate or remove pending jobs from the queue.
        Parameters
        ----------
        job_name : str
            name of the job to be killed
        regex : bool
            if true, all jobs matching job_name get killed!
        job_ids : Union[List[int], int]
            job Ids to be killed

        """

        if job_name is not None:
            job_ids = list(self.search_queue_for_jobname(job_name, regex=regex).index)
        elif job_ids is not None:
            if isinstance(job_ids, int):
                job_ids = [job_ids]
        else:
            raise ValueError("Please provide either job_name or job_ids!")

        if self.verbose:
            print("Stopping: " + ", ".join(map(str, job_ids)))
        try:
            bash.execute("bkill " + " ".join(map(str, job_ids)))
        except Exception as err:
            if any(["Job has already finished" in x for x in err.args]):
                print("Job has already finished")
            else:
                raise ChildProcessError("could not execute this command: \n" + str(err.args))
