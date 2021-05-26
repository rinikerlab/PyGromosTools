from typing import List, Union
import pandas as pd

class _SubmissionSystem:
    verbose: bool
    submission: bool

    _job_duration: str = "24:00"
    _nmpi: int
    _nomp: int
    _max_storage: float
    job_queue_list: pd.DataFrame #contains all jobs in the queue (from this users)

    def __init__(self, submission: bool = False,
                 nmpi: int = 1, nomp: int = 1, max_storage: float = 1000, job_duration: str = "24:00",
                 verbose: bool = False, enviroment=None):
        """
            Construct a submission system with required parameters.

        Parameters
        ----------
        submission : bool, optional
            if the job should be submitted? or only a dry run?
        nmpi : int, optional
            number of desired nmpi threads  (default: 1 == No MPI)
        nomp : int, optional
            number of desired omp threads (default: 1 == No OMP)
        max_storage : float, optional
            maximally required storage for a job in mb (default, 1000MB=1GB)
        job_duration : str, optional
            the duration of the job as str("HHH:MM")  (default: "24:00")
        verbose : bool, optional
            let me write you a book!  (default: False)
        enviroment: dict, optional
            here you can pass environment variables as dict{varname: value} (default: None)
        """
        self.verbose = verbose
        self.submission = submission

        self._job_duration = job_duration
        self._nmpi = nmpi
        self._nomp = nomp
        self._max_storage = max_storage
        self._enviroment = enviroment

    def submit_to_queue(self, **kargs) -> Union[int, None]:
        """submit_to_queue
            This function submits a str command to the submission system.

        Parameters
        ----------
        command : str
            command to be submitted
        kwargs

        Returns
        -------
        int, None
            if a job was submitted the jobID is returned else None.

        """
        raise NotImplemented("Do is not implemented for: " + self.__class__.__name__)

    def get_script_generation_command(self, var_name: str = None, var_prefixes: str = "") -> str:
        """
            This command is generating a script, that is used to submit a given command.

        Parameters
        ----------
        var_name
        var_prefixes

        Returns
        -------
        str
            command as a string

        """
        name = self.__class__.__name__
        if (var_name is None):
            var_name = var_prefixes + name

        gen_cmd = "#Generate " + name + "\n"
        gen_cmd += "from " + self.__module__ + " import " + name + " as " + name + "_obj" + "\n"
        gen_cmd += var_name + " = " + name + "_obj(submission=" + str(self.submission) + ", verbose=" + str(
            self.verbose) + ", nmpi="+str(self.nmpi)+", nomp="+str(self.nomp)+ ", job_duration=\""+str(self.job_duration)+"\")\n\n"
        return gen_cmd

    def search_queue_for_jobname(self, job_name: str, regex:bool=False, **kwargs)->pd.DataFrame:
        """get_jobs_from_queue

            this function searches the job queue for a certain job id.

        Parameters
        ----------
        job_name :  str
            text part of the job_line
        regex:  bool, optional
            if the string is a Regular Expression
        Raises
        -------
        NotImplemented
            Needs to be implemented in subclasses
        """

        raise NotImplemented("Do is not implemented for: " + self.__class__.__name__)

    def search_queue_for_jobid(self, job_id: int, **kwargs)->pd.DataFrame:
        """search_queue_for_jobid

            this jobs searches the job queue for a certain job id.

        Parameters
        ----------
        job_id :  int
            id of the job
        Raises
        -------
        NotImplemented
            Needs to be implemented in subclasses
        """

        raise NotImplemented("search_queue_for_jobID is not implemented for: " + self.__class__.__name__)

    def is_job_in_queue(self, job_name: str=None, job_id:int=None, _onlyRUNPEND:bool=True) -> bool:
        """
        checks wether a function is still in the lsf queue

        Parameters
        ----------
        job_name : str
            name of the job.
        job_id : int
            id of the job.
        verbose : bool, optional
            extra prints, by default False

        Returns
        -------
        bool
            is the job in the lsf queue?
        """
        if(not job_name is None):
            if(_onlyRUNPEND):
                queued_job_ids = self.search_queue_for_jobname(job_name=job_name)
                queued_job_ids = queued_job_ids.where(queued_job_ids.STAT.isin(["RUN", "PEND"])).dropna()
                return len(queued_job_ids) > 0
            else:
                return len(self.search_queue_for_jobname(job_name=job_name)) > 0
        elif(not job_id is None):
            if(_onlyRUNPEND):
                queued_job_ids = self.search_queue_for_jobid(job_id=job_id)
                queued_job_ids = queued_job_ids.where(queued_job_ids.STAT.isin(["RUN", "PEND"])).dropna()
                return len(queued_job_ids) > 0
            else:
                return len(self.search_queue_for_jobid(job_id=job_id)) >0
        else:
            raise ValueError("Please provide either the job_name or the job_id!")

    def kill_jobs(self, job_name:str=None, regex:bool=False, job_ids: Union[List[int], int]=None):
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
        raise NotImplemented("kill_jobs is not implemented for: " + self.__class__.__name__)


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