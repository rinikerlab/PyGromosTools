import os
import sys
import traceback
from collections import OrderedDict
from copy import deepcopy
from typing import Union

from pygromos.files.coord import cnf
from pygromos.files.gromos_system import Gromos_System
from pygromos.files.trajectory.trc import Trc
from pygromos.files.trajectory.tre import Tre
from pygromos.files.trajectory.trg import Trg
from pygromos.hpc_queuing.job_scheduling.schedulers import simulation_scheduler
from pygromos.hpc_queuing.job_scheduling.workers.analysis_workers import simulation_analysis
from pygromos.hpc_queuing.submission_systems._submission_system import _SubmissionSystem
from pygromos.hpc_queuing.submission_systems.local import LOCAL
from pygromos.utils import bash, utils
from pygromos.utils.utils import spacer as spacer


def simulation(in_gromos_system:Gromos_System, project_dir:str,
               step_name:str="sim", in_imd_path = None,
               submission_system:_SubmissionSystem=LOCAL(), simulation_runs:int=1, equilibration_runs:int = 0,
               previous_simulation_run:int=None, nmpi:int = 1,
               force_simulation:bool=False,
               analysis_script:callable = simulation_analysis.do, analysis_control_dict:dict = None,
               verbose:bool = True, verbose_lvl:int=1) -> Union[Gromos_System, int]:
    #PREPERATIONS
    try:
        try:
            step_dir = project_dir + "/" + step_name
            bash.make_folder(step_dir)

            out_input_dir = step_dir + "/input"
            out_simulation_dir = step_dir + "/simulation"
            out_analysis_dir = step_dir + "/analysis"
            bash.make_folder(out_input_dir)

            ##Prepare gromos system:
            in_gromos_system = deepcopy(in_gromos_system)
            in_gromos_system.work_folder = out_input_dir
            in_gromos_system.name = step_name
            if in_imd_path is None:
                in_gromos_system.adapt_imd()
            else:
                in_gromos_system.imd = in_imd_path
                in_gromos_system.adapt_imd()

            out_analysis_cnf = out_analysis_dir + "/data/" + in_gromos_system.name + ".cnf"

            if verbose:
                print(spacer)
                print(step_name)
                print(spacer)

            #Write out, all non promised files
            #TODO: REMOVE - print(in_gromos_system.all_file_paths)
            in_gromos_system._update_all_file_paths()
            in_gromos_system.write_files()
            #TODO: REMOVE - print(in_gromos_system.all_file_paths)

            #Write Out Ana Script
            in_analysis_control_dict = analysis_control_dict
            n_analysis_processors = 1 #Maybe in the future: 5 if(nmpi>5) else 1
        except Exception as err:
            raise Exception("Could not prepare the gromos System\n\t" + "\n\t".join(map(str, err.args)))
        # do
        analysis_control_dict = simulation_analysis.template_control_dict if (in_analysis_control_dict is None) else in_analysis_control_dict

        analysis_vars = OrderedDict({
            "in_simulation_dir": out_simulation_dir,
            "sim_prefix": in_gromos_system.name,
            "out_analysis_dir": out_analysis_dir,
            "gromosPP_bin_dir": in_gromos_system.gromosPP._bin,
            "control_dict": analysis_control_dict,
            "n_processes": n_analysis_processors,
            "verbose": verbose,
            "verbose_lvl": verbose_lvl
        })
        try:
            in_analysis_script_path = utils.write_job_script(out_script_path=step_dir + "/job_analysis.py",
                                                             target_function=analysis_script,
                                                             variable_dict=analysis_vars)
        except Exception as err:
            raise Exception("Could not prepare the analysis script\n\t" + "\n\t".join(map(str, err.args)))

        ##Write Out schedulling Script
        ###Build analysis_script
        MD_job_vars = OrderedDict({
            "in_simSystem": in_gromos_system,
            "out_dir_path": out_simulation_dir,
            "simulation_run_num": simulation_runs,
            "equilibration_run_num": equilibration_runs,
            "submission_system": submission_system,
            "analysis_script_path": in_analysis_script_path,
            "verbose": verbose,
            "verbose_lvl": verbose_lvl
        })
        try:
            in_scheduler_script_path = utils.write_job_script(out_script_path=step_dir + "/schedule_MD_job.py",
                                                             target_function=simulation_scheduler.do,
                                                             variable_dict=MD_job_vars)
        except Exception as err:
            raise Exception("Could not prepare the scheduling script\n\t" + "\n\t".join(map(str, err.args)))
    except Exception as err:
        traceback.print_exception(*sys.exc_info())
        raise Exception("Could not prepare the command block\n\t"+"\n\t".join(map(str, err.args)))

    ##schedule
    try:
        if((os.path.exists(out_analysis_cnf) and os.path.exists(out_simulation_dir+".tar")) and not force_simulation):
            if verbose:
                print(utils.spacer2+"FOUND RESULT: "+out_analysis_cnf+"\n GOING TO SKIPT THIS SUBMISSION!")
            #warnings.warn("Skipping active submission, as result CNF was found: \n"+out_analysis_cnf)
            last_jobID = None
        else:
            last_jobID = simulation_scheduler.do(in_simSystem=in_gromos_system, out_dir_path=out_simulation_dir,
                                                    simulation_run_num=simulation_runs, equilibration_run_num=equilibration_runs,
                                                    submission_system=submission_system, previous_job_ID=previous_simulation_run,
                                                    analysis_script_path=in_analysis_script_path, verbose=verbose, verbose_lvl=verbose_lvl)
    except Exception as err:
        traceback.print_exception(*sys.exc_info())
        raise Exception("Could not submit the commands\n\t"+"\n\t".join(map(str, err.args)))

    # Return the promise final system
    final_cnf_file = out_analysis_dir + "/data/" + in_gromos_system.name + ".cnf"
    if(os.path.exists(final_cnf_file)):
        in_gromos_system.cnf = cnf.Cnf(final_cnf_file)
    else:
        in_gromos_system.cnf = cnf.Cnf(in_value=None)
        in_gromos_system.cnf._future_file = True
        in_gromos_system.cnf.path = out_analysis_dir + "/data/" + in_gromos_system.name + ".cnf"

    # Return trajectories if available
    final_trc_file = out_analysis_dir + "/data/" + in_gromos_system.name + ".trc"
    if os.path.exists(final_trc_file+".h5"):
        in_gromos_system.trc = Trc(input_value=final_trc_file+".h5")
    elif os.path.exists(final_trc_file):
        in_gromos_system.trc = Trc(input_value=final_trc_file)

    final_tre_file = out_analysis_dir + "/data/" + in_gromos_system.name + ".tre"
    if os.path.exists(final_tre_file+".h5"):
        in_gromos_system.tre = Tre(input_value=final_tre_file+".h5")
    elif os.path.exists(final_tre_file):
        in_gromos_system.tre = Tre(input_value=final_tre_file)

    final_trg_file = out_analysis_dir + "/data/" + in_gromos_system.name + ".trg"
    if os.path.exists(final_trg_file+".h5"):
        in_gromos_system.trg = Trg(input_value=final_trg_file+".h5")
    elif os.path.exists(final_trg_file):
        in_gromos_system.trg = Trg(input_value=final_trg_file)

    in_gromos_system.work_folder = step_dir
    return in_gromos_system, last_jobID