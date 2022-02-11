#!/usr/bin/env python3
import os
import sys
import time
package_path = os.path.abspath(__file__+"/../../../../../..")
# print(package_path)
sys.path.append(package_path)

from pygromos.gromos import gromosXX as mdGromosXX
from pygromos.files.simulation_parameters import imd
from pygromos.utils import bash as bash
from pygromos.utils.utils import spacer3 as spacer, dynamic_parser, time_wait_s_for_filesystem



def work(out_dir : str, in_cnf_path : str, in_imd_path : str, in_top_path : str, runID:int=1,
         in_perttopo_path: str = None, in_disres_path: str= None, in_posres_path:str = None, in_refpos_path:str=None,
         in_qmmm_path:str=None, out_trc:bool=False, out_tre: bool=False,
         out_trg: bool = False, out_trv: bool = False, out_trf: bool = False, out_trs: bool = False,
         nmpi: int = 1, nomp: int = 1,
         reinitialize: bool = False, initialize_first_run:bool = True,
         gromosXX_bin_dir: str = None, work_dir: str = None, **kwargs):
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
    Returns
    -------
    int
        return number

    """
    time.sleep(time_wait_s_for_filesystem)
    # WORKDIR SetUP
    if ((work_dir is None or work_dir == "None") and "TMPDIR" in os.environ):
        work_dir = os.environ["TMPDIR"]
        print("using TmpDir")
    elif (work_dir  is None or work_dir == "None"):
        print("Could not find TMPDIR!\n Switched to outdir for work")
        work_dir = out_dir

        if (not os.path.isdir(work_dir)):
            bash.make_folder(work_dir)
    print("workDIR: " + str(work_dir))
    os.chdir(work_dir)

    #Prepare IMD file:
    tmp_prefix = os.path.basename(out_dir)
    tmp_imd_path = out_dir+"/"+tmp_prefix+".imd"
    imd_file = imd.Imd(in_imd_path)
    ##check init_block

    ###StochDyn
    is_stochastic_dynamics_sim = False
    is_vacuum = False
    print(type(imd_file.BOUNDCOND.NTB))
    if(imd_file.BOUNDCOND.NTB == 0):
        is_vacuum = True

    if (hasattr(imd_file, 'STOCHDYN')):
        if (imd_file.STOCHDYN.NTSD):
            is_stochastic_dynamics_sim = True

    if (reinitialize or (initialize_first_run and runID == 1)):
        imd_file.INITIALISE.NTIVEL = 1
        imd_file.INITIALISE.NTISHK = 0
        imd_file.INITIALISE.NTINHT = 0
        imd_file.INITIALISE.NTINHB = 0
        imd_file.INITIALISE.NTISHI = 0
        imd_file.INITIALISE.NTIRTC = 0
        imd_file.INITIALISE.NTICOM = 0
        imd_file.INITIALISE.NTISTI = 0

        if (is_stochastic_dynamics_sim):
            imd_file.INITIALISE.NTISHI = 1
            imd_file.INITIALISE.NTISTI = 1
        else:
            imd_file.INITIALISE.NTISHK = 3
            imd_file.INITIALISE.NTISHI = 1

    else:
        imd_file.INITIALISE.NTIVEL = 0
        imd_file.INITIALISE.NTISHK = 0
        imd_file.INITIALISE.NTINHT = 0
        imd_file.INITIALISE.NTINHB = 0
        imd_file.INITIALISE.NTISHI = 0
        imd_file.INITIALISE.NTIRTC = 0
        imd_file.INITIALISE.NTICOM = 0
        imd_file.INITIALISE.NTISTI = 0

        if (is_stochastic_dynamics_sim or is_vacuum):
            imd_file.INITIALISE.NTISHI = 1

        print(imd_file.INITIALISE.NTISHI, imd_file.BOUNDCOND.NTB, is_vacuum)
    ##Write out:
    tmp_imd_path = imd_file.write(tmp_imd_path)

    # RUN
    gromosXX = mdGromosXX.GromosXX(gromosXX_bin_dir=gromosXX_bin_dir)
    try:
        print(spacer + "\n start MD " + str(os.path.basename(tmp_imd_path)) + "\n")

        #TODO: This is a stupid workaround as Euler tends to place nans in the euler angles, that should not be there!
        from pygromos.files.coord import cnf
        import math
        cnf_file = cnf.Cnf(in_cnf_path)
        if(hasattr(cnf_file, "GENBOX") and any([math.isnan(x) for x in cnf_file.GENBOX.euler])):
            cnf_file.GENBOX.euler = [0.0, 0.0, 0.0]
            cnf_file.write(in_cnf_path)


        out_prefix = out_dir+"/"+tmp_prefix

        try:
            omd_file_path = gromosXX.md_run(in_topo_path=in_top_path, in_coord_path=in_cnf_path, in_imd_path=tmp_imd_path,
                                     in_pert_topo_path=in_perttopo_path, in_disres_path=in_disres_path,
                                     in_posresspec_path=in_posres_path, in_refpos_path=in_refpos_path,
                                     in_qmmm_path=in_qmmm_path, nmpi=nmpi, nomp=nomp,
                                     out_prefix=out_prefix,
                                     out_tre=out_tre, out_trc=out_trc,
                                     out_trg=out_trg, out_trs=out_trs, out_trf=out_trf, out_trv=out_trv,
                                     verbose=True)

            print("Waiting to find: ", omd_file_path.replace(".omd", ".cnf"))
            bash.wait_for_fileSystem(omd_file_path.replace(".omd", ".cnf"))

            md_failed = False
        except Exception as err:
            print("Failed! process returned: \n Err: \n" + "\n".join(err.args))
            md_failed = True

        if (out_dir != work_dir):
            bash.move_file(work_dir + "/*", out_dir)

        # post simulation cleanup
        if (not (work_dir is None or work_dir == "None") and work_dir != out_dir and not "TMPDIR" in os.environ):
            bash.remove_folder(work_dir, verbose=True)


    except Exception as err:
        print("\nFailed during simulations: ", file=sys.stderr)
        print(type(err), file=sys.stderr)
        print(err.args, file=sys.stderr)
        exit(1)
    if(md_failed):
        print("\nFailed during simulations: \n Checkout: \n "+str(out_prefix)+".omd", file=sys.stderr)
        exit(1)
    exit(0)


if __name__ == "__main__":
    # INPUT JUGGELING
    args = dynamic_parser(work, title="Run MD-Worker")
    # WORK Command
    work(**vars(args))
