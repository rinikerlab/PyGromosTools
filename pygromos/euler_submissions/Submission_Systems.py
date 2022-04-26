from pygromos.utils import bash
import subprocess, os, re, time
from typing import List


class _SubmissionSystem:

    def __init__(self):
        pass
    def submit_to_queue(self):
        raise ValueError("Do is not implemented for: "+self.__class__.__name__)

class LSF(_SubmissionSystem):
    """LSF
        This class is a wrapper for the LSF queueing system by IBM, like it is used on Euler.
    """
    def __init__(self):
        super().__init__()


    def submit_to_queue(self, command:str, jobName:str, duration:str="4:00",
             outLog=None, errLog=None, submit_from_dir:str=None,
             nmpi:int=1, nomp:int=1,
             maxStorage:int=1000,
             queue_after_jobID:int=None, force_queue_start_after:bool=False,
             projectName: str = None, jobGroup:str=None, priority=None,
             begin_mail:bool =False, end_mail:bool=False, post_execution_command:str=None,
             verbose:bool=True, noSubmission:bool=False, do_not_doubly_submit_to_queue:bool=False, stupid_mode=False, no_queueing=True, sumbit_from_file:bool=True) -> int:
        #job_properties:Job_properties=None, <- currently not usd
        """submit_to_queue

                This function submits a job to the submission queue.

        Parameters
        ----------
        command :   str
        jobName :   str
        duration :  str, optional
            time for the job as a string (e.g.: 4:00 - HH:MM)
        submit_from_dir:    str, optional
            from which dir shall the jobs be submitted? (def. current)
        outLog :    str, optional
            path to out.log location
        errLog :
            path to error.log location
        nmpi :  int, optional
            integer number of mpi cores (default: 1)
        nomp :  int, optional
            integer number of omp cores (default: 1) - WARNING DOES NOT WORK STABLE AT THE MOMENT!
        maxStorage : int, optional
            Max memory per core.
        projectName :  NOT IMPLEMENTED AT THE MOMENT
        jobGroup :  NOT IMPLEMENTED AT THE MOMENT
        priority :  NOT IMPLEMENTED AT THE MOMENT
        queue_after_jobID: int, optional
            shall this job be queued after another one?
        force_queue_start_after: bool, optional
            shall this job start after another job, no matter the exit state?
        begin_mail :    bool, optional
            send a mail when job starts
        end_mail :  bool, optional
            send a mail, when job is finished
        verbose:    bool, optional
            print out some messages
        noSubmission:   bool, optional
            do not submit the job to the queue.
        do_not_doubly_submit_to_queue:  bool, optional
            if True: script checks the submission queue and looks for an identical job_name. than raises ValueError if it is already submitted. (default: True)

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

        #generate submission_string:
        submission_string = ""

        #QUEUE checking to not double submit
        if(not stupid_mode and do_not_doubly_submit_to_queue):
            if(verbose): print('check queue')
            ids = self.get_jobs_from_queue(jobName)

            if(len(ids) > 0):
                if(verbose): print("\tSKIP - FOUND JOB: \t\t"+"\n\t\t".join(map(str, ids))+"\n\t\t with jobname: "+jobName)
                return ids[0]

        if(isinstance(submit_from_dir, str) and os.path.isdir(submit_from_dir)):
            os.chdir(submit_from_dir)
            command_file_path = submit_from_dir+"/job_"+str(jobName)+".sh"
        else:
            command_file_path = "./job_" + str(jobName) + ".sh"

        submission_string += "bsub "
        submission_string += " -J"+jobName+" "
        submission_string += " -W " + str(duration) + " "

        if(not isinstance(post_execution_command, type(None))):
            submission_string += "-Ep \""+post_execution_command+"\" "

        if(not isinstance(outLog, str) and not isinstance(errLog, str)):
            outLog=jobName+".out"
            submission_string+= " -o "+outLog
        elif(isinstance(outLog, str)):
            submission_string+= " -o "+outLog

        if(isinstance(errLog, str)):
            submission_string+= " -e "+errLog
        
        # adding the ptile
        submission_string += " -R span[ptile="+str(nmpi)+"] "

        nCPU=nmpi*nomp
        submission_string+= " -n "+str(nCPU)+" "
        add_string = ""
        #add_string= "-R \"select[model==XeonGold_5118 || model==XeonGold_6150 || model==XeonE3_1585Lv5 || model==XeonE3_1284Lv4 || model==XeonE7_8867v3 || model == XeonGold_6140 || model==XeonGold_6150 ]\""
        if(isinstance(maxStorage, int)):
            submission_string += " -R rusage[mem="+str(maxStorage)+"] "

        if(isinstance(queue_after_jobID, (int, str))):
            prefix = "done"
            if(force_queue_start_after):
                prefix = "exit"
            submission_string += " -w "+prefix+"("+str(queue_after_jobID)+") "

        if(begin_mail):
            submission_string += " -B "
        if(end_mail):
            submission_string += " -N "

        if (nomp > 1):
            command = "\"export OMP_NUM_THREADS=" + str(nomp) + ";\n " + command + "\""
        else:
            command = "\n " + command.strip() + ""

        if(sumbit_from_file):
            command_file = open(command_file_path, "w")
            command_file.write("#!/bin/bash\n")
            command_file.write(command+";\n")
            command_file.close()
            command = command_file_path

            bash.execute("chmod +x "+command_file_path)

        ##finalize string
        submission_string = list(map(lambda x: x.strip(), submission_string.split()))+[command]

        if(verbose): print("Submission Command: \t", " ".join(submission_string))
        if(not noSubmission and not stupid_mode):
            try:
                std_out_buff = bash.execute(command=submission_string)
                std_out = "\n".join(std_out_buff.readlines())

                # next sopt_job is queued with id:
                id_start = std_out.find("<")
                id_end = std_out.find(">")
                job_id = str(std_out[id_start + 1:id_end]).strip()
                if verbose: print("process returned id: " + str(job_id))
                if(job_id == "" and job_id.isalnum()):
                    raise ValueError("Did not get at job ID!")
            except:
                raise ChildProcessError("could not submit this command: \n" +
                str(submission_string))
        elif(not noSubmission and stupid_mode):
            os.system(submission_string)
            job_id = -1
        else:
            job_id = -1

        os.chdir(orig_dir)
        return int(job_id)



    def submit_jobAarray_to_queue(self, command:str, jobName:str, duration:str="4:00",
             outLog=None, errLog=None, submit_from_dir:str=None,
             nmpi:int=1, nomp:int=1,
             maxStorage:int=1000,
             queue_after_jobID:int=None, force_queue_start_after:bool=False,
             projectName: str = None, jobGroup:str=None, priority=None,
             begin_mail:bool =False, end_mail:bool=False,
             verbose:bool=True, noSubmission:bool=False, do_not_doubly_submit_to_queue:bool=True) -> int:
        #job_properties:Job_properties=None, <- currently not usd
        """submit_to_queue
        NOT IMPLEMENTED YET!
                This function submits a job to the submission queue.

        Parameters
        ----------
        command :   str
        jobName :   str
        duration :  str, optional
            time for the job as a string (e.g.: 4:00 - HH:MM)
        submit_from_dir:    str, optional
            from which dir shall the jobs be submitted? (def. current)
        outLog :    str, optional
            path to out.log location
        errLog :
            path to error.log location
        nmpi :  int, optional
            integer number of mpi cores (default: 1)
        nomp :  int, optional
            integer number of omp cores (default: 1) - WARNING DOES NOT WORK STABLE AT THE MOMENT!
        maxStorage : int, optional
            Max memory per core.
        projectName :  NOT IMPLEMENTED AT THE MOMENT
        jobGroup :  NOT IMPLEMENTED AT THE MOMENT
        priority :  NOT IMPLEMENTED AT THE MOMENT
        queue_after_jobID: int, optional
            shall this job be queued after another one?
        force_queue_start_after: bool, optional
            shall this job start after another job, no matter the exit state?
        begin_mail :    bool, optional
            send a mail when job starts
        end_mail :  bool, optional
            send a mail, when job is finished
        verbose:    bool, optional
            print out some messages
        noSubmission:   bool, optional
            do not submit the job to the queue.
        do_not_doubly_submit_to_queue:  bool, optional
            if True: script checks the submission queue and looks for an identical job_name. than raises ValueError if it is already submitted. (default: True)

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
        #generate submission_string:
        submission_string = ""

        if(isinstance(submit_from_dir, str) and os.path.isdir(submit_from_dir)):
            submission_string += "cd "+submit_from_dir+ " && "

        submission_string += "bsub -J \" "+jobName+" \" -W \""+str(duration)+"\" "

        if(isinstance(jobGroup, str)):
            submission_string += " -g "+jobGroup+" "

        if(not isinstance(outLog, str) and not isinstance(errLog, str)):
            outLog=jobName+".out"
            submission_string+= " -oo "+outLog
        elif(isinstance(outLog, str)):
            submission_string+= " -oo "+outLog

        if(isinstance(errLog, str)):
            submission_string+= " -eo "+errLog

        nCPU=nmpi*nomp
        submission_string+= " -n "+str(nCPU)+" "

        if(isinstance(maxStorage, int)):
            submission_string += " -R \"rusage[mem="+str(maxStorage)+"]\" "

        if(isinstance(queue_after_jobID, (int, str))):
            prefix = "\"done"
            if(force_queue_start_after):
                prefix = "\"exit"
            submission_string += " -w "+prefix+"("+str(queue_after_jobID)+")\" "

        if(begin_mail):
            submission_string += " -B "
        if(end_mail):
            submission_string += " -N "


        if(nomp>1):
            command = " \" export OMP_NUM_THREADS="+str(nomp)+" && "+command+"\""
        else:
            command = " \""+command+"\""

        ##finalize string
        ##finalize string
        submission_string = list(map(lambda x: x.strip(), submission_string.split())) + [command]

        if (verbose): print("Submission Command: \t", " ".join(submission_string))
        if (not noSubmission):
            try:
                std_out_buff = bash.execute(command=submission_string)
                std_out = "\n".join(std_out_buff.readlines())

                # next sopt_job is queued with id:
                id_start = std_out.find("<")
                id_end = std_out.find(">")
                job_id = str(std_out[id_start + 1:id_end]).strip()
                if verbose: print("process returned id: " + str(job_id))
                if (job_id == "" and job_id.isalnum()):
                    raise ValueError("Did not get at job ID!")
            except:
                raise ChildProcessError("could not submit this command: \n" + submission_string)
        else:
            job_id = None
        return int(job_id)

    def get_jobs_from_queue(self, job_name:str , regex:bool=False)->List[int]:
        """search_queue_for_jobname

            this jobs searches the job queue for a certain job id.

        Parameters
        ----------
        job_name :  str
        regex:  bool, optional
            if the string is a Regular Expression
        Returns
        -------
        List[int]
            output contains all ids of fitting jobs to the querry
        """
        out_job_lines = self.search_queue_for_jobname(job_name, regex=regex)
        get_job_ids = list(map(int, filter(lambda x: x.isdigit(), map(lambda x: x.split(" ")[0], out_job_lines))))
        return get_job_ids

    def search_queue_for_jobname(self, job_name:str, regex:bool=False )->List[str]:
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
            if(not regex):
                job_name = " "+job_name+" "
            out_queue_lines = bash.execute("bjobs -w", wait_fail=True)
            out_job_lines = list(filter(lambda line: re.search(job_name, line), out_queue_lines.read().split("\n")))
        except TimeoutError:
            return []
        return out_job_lines

class LOCAL(_SubmissionSystem):
    def __init__(self):
        super().__init__()

    def submit_to_queue(self, command):
        bash.execute(command=command)
