import os, sys
import warnings
import traceback
from typing import List
from collections import OrderedDict
import numpy as np

from pygromos.files.gromos_system import Gromos_System
from pygromos.files.coord import cnf
from pygromos.files.blocks.imd_blocks import PERTURBATION, PRECALCLAM, WRITETRAJ

from pygromos.data.simulation_parameters_templates import template_emin, template_md, template_sd

from pygromos.utils import bash
from pygromos.utils import utils
from pygromos.utils.utils import spacer as spacer

from pygromos.hpc_queuing.job_scheduling.schedulers import simulation_scheduler
from pygromos.hpc_queuing.job_scheduling.workers.analysis_workers import simulation_analysis
from pygromos.hpc_queuing.submission_systems.Submission_Systems import LOCAL, _SubmissionSystem

"""
    Simulations
"""
def emin(in_gromos_system: Gromos_System, project_dir: str, step_name: str = "emin", in_imd_path=template_emin,
         submission_system: _SubmissionSystem = LOCAL(), simulation_runs: int = 1, equilibration_runs: int = 0,
         previous_simulation_run: int = None, analysis_script: callable = simulation_analysis.do) ->(Gromos_System, int):

    template_emin_control_dict = OrderedDict({
        "concat": {"do": True,
                   "sub": {
                       "cp_cnf": True,
                       "cat_trc": False,
                       "cat_tre": False,
                       "convert_trcs": False,
                   }
                   },
        "simulation_folder": {
            "do": True,
            "sub": {
                "tar": True,
                "remove": False
            }
        }
    })


    return simulation(in_gromos_system=in_gromos_system, project_dir=project_dir,  previous_simulation_run=previous_simulation_run,
                      step_name=step_name, in_imd_path=in_imd_path, submission_system=submission_system,
                      simulation_runs=simulation_runs, equilibration_runs=equilibration_runs, analysis_control_dict = template_emin_control_dict,
                      analysis_script=analysis_script)

def md(in_gromos_system: Gromos_System, project_dir: str, step_name: str = "md", in_imd_path=template_md,
       submission_system: _SubmissionSystem = LOCAL(), simulation_runs: int = 1, equilibration_runs: int = 0,
       previous_simulation_run: int = None, analysis_script: callable = simulation_analysis.do) ->(Gromos_System, int):
    return simulation(in_gromos_system=in_gromos_system, project_dir=project_dir,  previous_simulation_run=previous_simulation_run,
                      step_name=step_name, in_imd_path=in_imd_path, submission_system=submission_system,
                      simulation_runs=simulation_runs, equilibration_runs=equilibration_runs,
                      analysis_script=analysis_script)

def sd(in_gromos_system: Gromos_System, project_dir: str, step_name: str = "sd", in_imd_path=template_sd,
       submission_system: _SubmissionSystem = LOCAL(), simulation_runs: int = 1, equilibration_runs: int = 0,
       previous_simulation_run: int = None, analysis_script: callable = simulation_analysis.do) ->(Gromos_System, int):
    return simulation(in_gromos_system=in_gromos_system, project_dir=project_dir, previous_simulation_run=previous_simulation_run,
                      step_name=step_name, in_imd_path=in_imd_path, submission_system=submission_system,
                      simulation_runs=simulation_runs, equilibration_runs=equilibration_runs,
                      analysis_script=analysis_script)

"""
    Free Energy
"""
def TI_sampling(in_gromos_system: Gromos_System, project_dir: str, step_name="lambda_sampling",
                lambda_values: List[float] = np.arange(0, 1.1, 0.1), write_coordinates: int = 100,
                write_energies: int = 100, simulation_steps: int = 500, subSystem: _SubmissionSystem = LOCAL(),
                n_simulation_repetitions: int = 3, n_equilibrations: int = 1):

    work_dir = bash.make_folder(project_dir + "/" + step_name)
    in_gromos_system.save(work_dir)
    sys_path = in_gromos_system.save(work_dir+"/"+step_name+"_input.obj") #Ugly workaround for deepcopy
    general_system_name_prefix = in_gromos_system.name

    lam_systems = []
    for lam in lambda_values:
        lam = np.round(lam, 2)

        lam_system = Gromos_System.load(sys_path)
        lam_system.name = general_system_name_prefix+"_l_" + str(lam)
        lam_system.work_folder = work_dir

        # IMD
        ## Pertubation
        ###Pertubation of system.
        pert_block = PERTURBATION(NTG=1, NRDGL=0, RLAM=lam, DLAMT=0,
                                  ALPHC=0.5, ALPHLJ=0.5, NLAM=2, NSCALE=0)
        lam_system.imd.add_block(block=pert_block)

        ###Calculate additional lambda points for end-states
        if(not hasattr(lam_system, "PRECALCLAM")):
            precalc_lam_block = PRECALCLAM(NRLAM=2, MINLAM=0, MAXLAM=1)
            lam_system.imd.add_block(block=precalc_lam_block)

        ## write out trajs:
        write_traj = WRITETRAJ(NTWX=write_coordinates, NTWE=write_energies, NTWG=write_energies)
        lam_system.imd.add_block(block=write_traj)

        lam_system.imd.STEP.NSTLIM = simulation_steps

        # Submit
        out_gromos_system, jobID = _TI_lam_step(in_gromos_system=lam_system, project_dir=work_dir,
                                                step_name=lam_system.name, submission_system=subSystem,
                                                simulation_runs=n_simulation_repetitions,
                                                equilibration_runs=n_equilibrations)

        out_gromos_system.save(out_gromos_system.work_folder + "/sd_out_system.obj")
        lam_systems.append(out_gromos_system)

    return lam_system, jobID



def _TI_lam_step(in_gromos_system: Gromos_System, project_dir: str, step_name: str = "lam", in_imd_path=template_sd,
                 submission_system: _SubmissionSystem = LOCAL(), simulation_runs: int = 1, equilibration_runs: int = 0,
                 previous_simulation_run: int = None, analysis_script: callable = simulation_analysis.do) ->(Gromos_System, int):
    template_emin_control_dict = OrderedDict({
        "concat": {"do": True,
                   "sub": {
                       "cp_cnf": True,
                       "cat_trc": True,
                       "cat_tre": True,
                       "cat_trg": True,
                       "convert_trcs": False,
                   }
                   },
        "simulation_folder": {
            "do": False,
            "sub": {
                "tar": False,
                "remove": False
            }
        }
    })

    return simulation(in_gromos_system=in_gromos_system, project_dir=project_dir, previous_simulation_run=previous_simulation_run,
                      step_name=step_name, in_imd_path=in_imd_path, submission_system=submission_system,
                      simulation_runs=simulation_runs, equilibration_runs=equilibration_runs, analysis_control_dict=template_emin_control_dict,
                      analysis_script=analysis_script)



def simulation(in_gromos_system:Gromos_System, project_dir:str,
               step_name:str="sim", in_imd_path = None,
               submission_system:_SubmissionSystem=LOCAL(), simulation_runs:int=1, equilibration_runs:int = 0,
               previous_simulation_run:int=None, nmpi:int = 1,
               force_simulation:bool=False,
               analysis_script:callable = simulation_analysis.do, analysis_control_dict:dict = None,
               verbose:bool = False) ->(Gromos_System, int):
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
            in_gromos_system.work_folder = out_input_dir
            in_gromos_system.name = step_name
            if(not in_imd_path is None and in_gromos_system.imd.path is None):
                in_gromos_system.imd = in_imd_path
                in_gromos_system.adapt_imd()
            elif(not (in_gromos_system.cnf._future_file or in_gromos_system.imd._future_file)):
                in_gromos_system.adapt_imd()

            out_analysis_cnf = out_analysis_dir + "/data/" + in_gromos_system.name + ".cnf"

            print(spacer)
            print(step_name)
            print(spacer)

            #Write out, all non promised files
            #TODO: REMOVE - print(in_gromos_system.all_file_paths)
            in_gromos_system._update_all_file_paths()
            in_gromos_system.write_files()
            #TODO: REMOVE - print(in_gromos_system.all_file_paths)

            # %%
            ##Write Out Ana Script
            # IN args
            in_analysis_control_dict = analysis_control_dict
            n_analysis_processors = 1 #Maybe in the future: 5 if(nmpi>5) else 1
            verbose = verbose
        except Exception as err:
            raise Exception("Could not prepare the gromos System\n\t" + "\n\t".join(map(str, err.args)))
        # do
        analysis_control_dict = simulation_analysis.template_control_dict if (in_analysis_control_dict is None) else in_analysis_control_dict

        analysis_vars = OrderedDict({
            "in_simulation_dir": out_simulation_dir,
            "sim_prefix": in_gromos_system.name,
            "out_analysis_dir": out_analysis_dir,
            "gromosPP_bin_dir": in_gromos_system.gromosPP.bin,
            "control_dict": analysis_control_dict,
            "n_processes": n_analysis_processors,
            "verbose": verbose,
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
            "nmpi" : nmpi,
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

    # %%
    ##schedule
    try:
        if((os.path.exists(out_analysis_cnf) and os.path.exists(out_simulation_dir+".tar")) and not force_simulation):
            print(utils.spacer2+"FOUND RESULT: "+out_analysis_cnf+"\n GOING TO SKIPT THIS SUBMISSION!")
            #warnings.warn("Skipping active submission, as result CNF was found: \n"+out_analysis_cnf)
            last_jobID = None
        else:
            last_jobID = simulation_scheduler.do(in_simSystem=in_gromos_system, out_dir_path=out_simulation_dir,
                                                    simulation_run_num=simulation_runs, equilibration_run_num=equilibration_runs,
                                                    submission_system=submission_system, previous_job_ID=previous_simulation_run,
                                                    analysis_script_path=in_analysis_script_path, nmpi=nmpi)
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

    return in_gromos_system, last_jobID

