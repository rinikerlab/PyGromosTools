from typing import List, Union

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
        name = self.__class__.__name__
        if (var_name is None):
            var_name = var_prefixes + name

        gen_cmd = "#Generate " + name + "\n"
        gen_cmd += "from " + self.__module__ + " import " + name + " as " + name + "_obj" + "\n"
        gen_cmd += var_name + " = " + name + "_obj(submission=" + str(self.submission) + ", verbose=" + str(
            self.verbose) + ", nmpi="+str(self.nmpi)+", nomp="+str(self.nomp)+ ", job_duration=\""+str(self.job_duration)+"\")\n\n"
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