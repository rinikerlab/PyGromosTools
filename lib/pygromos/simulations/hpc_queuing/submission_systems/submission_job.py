class Submission_job:
    """
    Description:
    This class stores parameters for the submission of jobs. It is used by the submission_systems:
    - submission dummy
    - submission local
    - submission lsf

    It should handle all the information required for a single job, while the submission_systems handles more permanent settings.

    It should provide an easy way to modify jobs, even from high level modules (e.g. the simulation module like TI, EMIN, ...).

    Author: Marc Lehner
    """

    def __init__(
        self,
        command: str = None,
        jobName: str = None,
        outLog: str = None,
        errLog: str = None,
        start_job: int = None,
        end_job: int = None,
        jobLim: int = None,
        queue_after_jobID: int = None,
        post_execution_command: str = None,
        submit_from_dir: str = None,
        sumbit_from_file: bool = True,
        jobGroup: str = None,
        jobID=None,
    ) -> None:
        self._command = command
        self._jobName = jobName
        self._outLog = outLog
        self._errLog = errLog
        self._start_job = start_job
        self._end_job = end_job
        self._jobLim = jobLim
        self._queue_after_jobID = queue_after_jobID
        self._post_execution_command = post_execution_command
        self._submit_from_dir = submit_from_dir
        self._sumbit_from_file = sumbit_from_file
        self._jobGroup = jobGroup
        self._jobID = jobID

    @property
    def command(self) -> str:
        if self._command is None:
            raise ValueError("command not set")
        return self._command

    @command.setter
    def command(self, command: str) -> None:
        if isinstance(command, str):
            self._command = command
        else:
            raise ValueError("command must be a string")

    @property
    def jobName(self) -> str:
        if self._jobName is None:
            return "test"
        return self._jobName

    @property
    def outLog(self) -> str:
        return self._outLog

    @property
    def errLog(self) -> str:
        return self._errLog

    @property
    def start_job(self) -> int:
        return self._start_job

    @property
    def end_job(self) -> int:
        return self._end_job

    @property
    def jobLim(self) -> int:
        return self._jobLim

    @jobLim.setter
    def jobLim(self, jobLim: int) -> None:
        if isinstance(jobLim, int):
            self._jobLim = jobLim
        else:
            raise ValueError("jobLim must be an int")

    @property
    def queue_after_jobID(self) -> int:
        return self._queue_after_jobID

    @property
    def post_execution_command(self) -> str:
        return self._post_execution_command

    @property
    def submit_from_dir(self) -> str:
        return self._submit_from_dir

    @property
    def sumbit_from_file(self) -> bool:
        return self._sumbit_from_file

    @property
    def jobGroup(self) -> str:
        return self._jobGroup

    @property
    def jobID(self) -> int:
        return self._jobID

    @jobID.setter
    def jobID(self, jobID: int) -> None:
        if isinstance(jobID, int):
            self._jobID = jobID
        else:
            raise ValueError("jobID must be an int")
