#!/usr/bin/env python3
import os
import sys
import time
import math

package_path = os.path.abspath(__file__+"/../../../../../..")
sys.path.append(package_path)

from pygromos.gromos import gromosXX as mdGromosXX
from pygromos.files.simulation_parameters import imd
from pygromos.files.coord import cnf

import glob

from pygromos.utils import bash as bash
from pygromos.utils.utils import spacer3 as spacer, dynamic_parser, time_wait_s_for_filesystem

import pygromos.simulations.hpc_queuing.job_scheduling.workers.simulation_workers.clean_up_simulation_files as zip_files

def work(out_dir : str, in_cnf_path : str, in_imd_path : str, in_top_path : str, runID:int=1,
         in_perttopo_path: str = None, in_disres_path: str= None, in_posres_path:str = None, in_refpos_path:str=None,
         in_qmmm_path:str=None, out_trc:bool=False, out_tre: bool=False,
         out_trg: bool = False, out_trv: bool = False, out_trf: bool = False, out_trs: bool = False,
         nmpi: int = 1, nomp: int = 1,
         reinitialize_every_run: bool = False, initialize_first_run:bool = True,
         gromosXX_bin_dir: str = None, work_dir: str = None, zip_trajectories: bool = True, **kwargs):
    """
    Executed by repex_EDS_long_production_run as workers

    Parameters
    ----------
    out_dir : str
         final output dir
    in_cnf_path : str
        input coordinates
    in_imd_path : str
        input imd-parameter file
    in_top_path : str
        input topology
    in_perttopo_path : str
        input pertubation
    in_disres_path : str
        input disres
    in_qmmm_path : str
        input qmmm
    nmpi : int, optional
        number of mpi cores (def.=1)
    nomp : int, optional
        number of omp cores (def.= 1)
    out_trg : bool, optional
        True if trg shall be written out.
    out_trv : bool, optional
        True if trv shall be written out.
    out_trf : bool, optional
        True if trf shall be written out.
    out_trs : bool, optional
        True if trs shall be written out.
    gromosXX_bin_dir : str, optional
         path to gromos binary directory
    work_dir : str, optional
         work directory
    zip_trajectories: bool
        determines whether trajectories are zipped

    Returns
    -------
    int
        return number

    """
    time.sleep(time_wait_s_for_filesystem)
    # WORKDIR SetUP
    if (work_dir is None or work_dir == "None") and "TMPDIR" in os.environ:
        work_dir = os.environ["TMPDIR"]
        print("using TMPDIR")
    else:
        print("Could not find TMPDIR!\n Switched to outdir for work")
        work_dir = out_dir
        if (not os.path.isdir(work_dir)): bash.make_folder(work_dir)
    print("workDIR: " + str(work_dir))
    
    # Check if the calculation is running on multiple nodes:
    if 'LSB_HOSTS' in os.environ:
        hosts = os.environ['LSB_HOSTS'].split()
    else:
        hosts = []
    multi_node = True if len(hosts) > 1 else False
    
    # run a euler script to create tmpdir on all nodes
    if multi_node: os.system('remote_tmpdir create')
    
    os.chdir(work_dir)
    
    #Prepare IMD file:
    tmp_prefix = os.path.basename(out_dir)
    tmp_imd_path = out_dir+"/"+tmp_prefix+".imd"
    imd_file = imd.Imd(in_imd_path)
    
    is_repex_run = True if (hasattr(imd_file, "REPLICA") and imd_file.REPLICA is not None) \ 
                        or (hasattr(imd_file, "REPLICA_EDS") and imd_file.REPLICA_EDS is not None) else False

    #Prepare CNF file(s):
    
    if is_repex_run:
        cnf_paths = glob.glob(in_cnf_path.replace("1.cnf", "*.cnf")) 
        cnf_file = cnf.Cnf(cnf_paths[0]) # used to make sure imd/cnf is in sync
    else: 
        cnf_file = cnf.Cnf(in_cnf_path)

    ##check init_block - if specified!
    ###What kind of simulation
    is_stochastic_dynamics_sim = False
    is_vacuum = False
    is_energymin_sim = False

    if(imd_file.BOUNDCOND.NTB == 0):
        is_vacuum = True

    if (hasattr(imd_file, 'STOCHDYN')):
        if (imd_file.STOCHDYN.NTSD):
            is_stochastic_dynamics_sim = True

    if(hasattr(imd_file, "ENERGYMIN")):
        if (imd_file.ENERGYMIN.NTEM > 0):
            is_energymin_sim = True
            
    #Adapt Initializations:
    if (reinitialize_every_run or (initialize_first_run and runID == 1)):
        imd_file.INITIALISE.NTIVEL = 1

        if(hasattr(imd_file, "CONSTRAINT")):
                 imd_file.INITIALISE.NTISHK = 0 if(imd_file.CONSTRAINT.NTC > 0) else 1
         
        if(hasattr(imd_file, "MULTIBATH")):
                 imd_file.INITIALISE.NTINHT = 0 if(imd_file.MULTIBATH.ALGORITHM <= 1) else 1
                  
        imd_file.INITIALISE.NTISHI = 0 if(hasattr(cnf_file, "LATTICESHIFT")) else 1

        imd_file.INITIALISE.NTIRTC = 0
        imd_file.INITIALISE.NTICOM = 0
        imd_file.INITIALISE.NTISTI = 0

        if (is_stochastic_dynamics_sim):
            imd_file.INITIALISE.NTISTI = 1

    elif(initialize_first_run and runID > 1):
        imd_file.INITIALISE.NTIVEL = 0
        imd_file.INITIALISE.NTISHI = 0 
        imd_file.INITIALISE.NTISHK = 0 
        imd_file.INITIALISE.NTINHT = 0 
        imd_file.INITIALISE.NTINHB = 0
        imd_file.INITIALISE.NTIRTC = 0
        imd_file.INITIALISE.NTICOM = 0
        imd_file.INITIALISE.NTISTI = 0

        if (is_stochastic_dynamics_sim or is_vacuum):
            imd_file.INITIALISE.NTISHI = 1

    ##Write out:
    tmp_imd_path = imd_file.write(tmp_imd_path)

    # RUN
    gromosXX = mdGromosXX.GromosXX(gromosXX_bin_dir=gromosXX_bin_dir)
    try:
        print(spacer + "\n start MD " + str(os.path.basename(tmp_imd_path)) + "\n")

        #TODO: This is a stupid workaround as Euler tends to place nans in the euler angles, that should not be there!
        if is_repex_run: # do the workaround for each cnf one by one
            for tmp_cnf_path in cnfs_paths:
                tmp_cnf = cnf.Cnf(tmp_cnf_path)
                if hasattr(tmp_cnf, "GENBOX") and any([math.isnan(x) for x in tmp_cnf.GENBOX.euler])):
                    tmp_cnf.GENBOX.euler = [0.0, 0.0, 0.0]
                    tmp_cnf.write(tmp_cnf_path)
    
        elif:(hasattr(cnf_file, "GENBOX") and any([math.isnan(x) for x in cnf_file.GENBOX.euler])):
            cnf_file.GENBOX.euler = [0.0, 0.0, 0.0]
            cnf_file.write(in_cnf_path)
        
        # Start the execution of the gromosXX binary
        try:
            if is_repex_run:
                omd_file_path = gromosXX.repex_run(in_topo_path=in_top_path, in_coord_path=in_cnf_path, in_imd_path=tmp_imd_path,
                                      in_pert_topo_path=in_perttopo_path, in_disres_path=in_disres_path,
                                      in_posresspec_path=in_posres_path, in_refpos_path=in_refpos_path,
                                      nmpi=nmpi, nomp=nomp,
                                      out_prefix=tmp_prefix,
                                      out_tre=out_tre, out_trc=out_trc,
                                      out_trg=out_trg, out_trs=out_trs, out_trf=out_trf, out_trv=out_trv,
                                      verbose=True)
            else:
                omd_file_path = gromosXX.md_run(in_topo_path=in_top_path, in_coord_path=in_cnf_path, in_imd_path=tmp_imd_path,
                                     in_pert_topo_path=in_perttopo_path, in_disres_path=in_disres_path,
                                     in_posresspec_path=in_posres_path, in_refpos_path=in_refpos_path,
                                     in_qmmm_path=in_qmmm_path, nmpi=nmpi, nomp=nomp,
                                     out_prefix=tmp_prefix,
                                     out_tre=out_tre, out_trc=out_trc,
                                     out_trg=out_trg, out_trs=out_trs, out_trf=out_trf, out_trv=out_trv,
                                     verbose=True)

            print("Waiting to find: ", omd_file_path.replace(".omd", ".cnf"))
            bash.wait_for_fileSystem(omd_file_path.replace(".omd", ".cnf"))

            md_failed = False
        except Exception as err:
            print("Failed! process returned: \n Err: \n" + "\n".join(err.args))
            md_failed = True
        
        # zip the files after the simulation.
        n_cpu_zip = nmpi if nmpi >= nomp else nomp
        if not multi_node and zip_trajectories:
            zip_files.do(in_simulation_dir=work_dir, n_processes=n_cpu_zip)    

        # Copy the files back to the proper directory when calc occured on scratch
        if (out_dir != work_dir):
            if not multi_node:
                bash.move_file(work_dir + "/*", out_dir)
            else: 
                for host in hosts:
                    command = 'ssh ' + host + '  \"mv ${TMPDIR}/* ' + out_dir + '\"'
                    os.system(command)
            os.system('remote_tmpdir delete') # Works for both multi or single node                
        
        # Note: If job is multi-node, it is simpler to zip things in out_dir after copying back
        if multi_node and zip_trajectories: 
            zip_files.do(in_simulation_dir=out_dir, n_processes=n_cpu_zip)

    except Exception as err:
        print("\nFailed during simulations: ", file=sys.stderr)
        print(type(err), file=sys.stderr)
        print(err.args, file=sys.stderr)
        exit(1)
    if(md_failed):
        print("\nFailed during simulations: \n Checkout: \n "+str(tmp_prefix)+".omd", file=sys.stderr)
        exit(1)
    exit(0)

if __name__ == "__main__":
    # INPUT JUGGELING
    args = dynamic_parser(work, title="Run MD-Worker")
    # WORK Command
    work(**vars(args))
