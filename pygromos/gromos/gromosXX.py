"""
FUNCTIONLIB:            wrapper for gromosXX
Description:
    This file contains python wrappers for the bash commandline of gromosXX

Author: Benjamin Schroeder
"""

import os
import datetime
import time
from typing import Union

from pygromos.utils import bash
from pygromos.utils.utils import time_wait_s_for_filesystem
from pygromos.gromos.utils import gromosTypeConverter

class _Gromos:
    """
    GromosXX

    This is the gromosXX baseclass. This should be inherited by a concrete class that might reimplement some new features, that are version dependent.

    Attributes:
    -----------
    bin :   str, optional
        This is the path to the folder containing the binaries of gromosXX. If None, the bash enviroment variables  will be used.
    """
    _bin: str = ""

    def __init__(self, gromosXX_bin_dir: str = None):
        """
        Constructing a gromosXX object.

        Parameters
        ----------
        gromosXX_bin_dir :   str, optional
            This is the path to the folder containing the binaries of gromosXX. If None, the bash enviroment variables  will be used.
        """
        #lazy me - doc text for functions:
        functions_text = "\n    Methods:\n    ---------\n" + "\n".join(["\t" + x for x in dir(self) if (not x.startswith("_") and callable(getattr(self, x)))])
        self.__doc__ = self.__doc__+functions_text

        if (isinstance(gromosXX_bin_dir, type(None)) or gromosXX_bin_dir == "None"):
            self._bin = ""
        else:
            self._bin = gromosXX_bin_dir + "/"

    def __str__(self):
        return self.__doc__
    def __repr__(self):
        return self.__str__()

    @property
    def bin(self)->Union[str, None]:
        if(not hasattr(self, "_bin") or self._bin == ""):
            return None
        else:
            return self._bin

    @gromosTypeConverter
    def md_run(self, in_topo_path: str, in_coord_path: str, in_imd_path: str, out_prefix: str,
                   in_pert_topo_path:str=None, in_disres_path:str=None, in_posresspec_path:str=None, in_refpos_path:str=None,
                   in_qmmm_path:str=None, nomp: int = 1, nmpi: int = 1,
                   out_trc: bool = False, out_tre: bool = False, out_trv: bool = False, out_trf: bool = False,
                   out_trs: bool = False, out_trg:bool = False,
                   verbose:bool = False) -> str:
        """
        This function is a wrapper for gromosXX md_mpi. You can directly execute the gromosXX md_mpi in a bash enviroment here.

        Warnings
        --------
            Hybrid jobs are possible, but very difficult to implement correctly to Euler and performance gain is questionable.
            If OMP should be used, I suggest the md_run - function.

        Parameters
        ----------
        in_topo_path :    str
                    This is the path to the input topology file (x.top)

        in_coord_path :    str
                    This is the path to the input coordinate file (x.cnf)

        in_imd_path :    str
                    This is the path to the input simulation parameter file (x.imd)

        out_prefix :    str
                    This prefix, define the output name.

        in_pert_topo_path : str, optional
                    This is the path to the pertubation file (x.ptp)

        in_disres_path :    str, optional
                    This is the path to the distance restraint file (x.dat)

        in_posresspec_path :    str, optional
                    This is the path to the position restraint file (x.pos)

        in_refpos_path :    str, optional
                    This is the path to the reference position file (x.rpf)

        nomp :  int, optional
                    How many omp cores shall be used? Prerequesite, gromos was compiled with -enableOMP

        nmpi : int, optional
                    How many mpi cores shall be used? Prerequesite, gromos was compiled with -enableMPI (and suggested to have -disableOMP)

        out_trc :   bool, optional
                    do you want to output a coordinate trajectory (x.trc) file? (needs also an output number in write block of imd!)

        out_tre :   bool, optional
                    do you want to output a coordinate trajectory (x.trc) file? (needs also an output number in write block of imd!)

        out_trs :   bool, optional
                    do you want to output a coordinate trajectory (x.trc) file? (needs also an output number in write block of imd!)

        out_trg :   bool, optional
                    do you want to output the free energy trajectory (x.trg) file? (needs also an output number in write block of imd!)

        Returns
        -------
        str
            returns the log_file_path of the simulation.

        Raises
        ------
        ChildProcessError
            If the execution of the simulation fails, this error is raised.

        See Also
        --------
        md_run, repex_mpi

        """

        #BUILD UP str. Command
        command = []
        if(nmpi > 1 and nomp > 1):
            raise ValueError("There are no Hybrid NMPI and NOMP jobs possible with gromos!")
        elif(nmpi > 1):
            command += ["mpirun -n " + str(nmpi*nomp)+" "] # --loadbalance  " --cpus-per-proc " +  + " "
            command += [self._bin + "md_mpi"]
        elif(nomp >= 1):
            command += ["export OMP_NUM_THREADS=" + str(nomp) + "  && "]
            command += [self._bin + "md"]
        else:
            command += [self._bin + "md"]

        command += ["@topo", str(in_topo_path)]
        command += ["@conf", str(in_coord_path)]
        command += ["@input", str(in_imd_path)]

        if not in_pert_topo_path is None:
            command += ["@pttopo", str(in_pert_topo_path)]

        if not in_disres_path is None:
            command += ["@distrest", str(in_disres_path)]

        if not in_posresspec_path is None:
            print("POSRES", in_posresspec_path)
            command += ["@posresspec", str(in_posresspec_path)]

        if not in_refpos_path is None:
            command += ["@refpos", str(in_refpos_path)]

        if not in_qmmm_path is None:
            command += ["@qmmm", str(in_qmmm_path)]

        if isinstance(out_prefix, str):
            command += ["@fin", str(out_prefix) + ".cnf"]
            log_file_path = out_prefix + ".omd"

            if out_trc:
                command += ["@trc", str(out_prefix) + ".trc"]
            if out_trv:
                command += ["@trv", str(out_prefix) + ".trv"]
            if out_trf:
                command += ["@trf", str(out_prefix) + ".trf"]
            if out_tre:
                command += ["@tre", str(out_prefix) + ".tre"]
            if out_trs:
                command += ["@trs", str(out_prefix) + ".trs"]
            if out_trg:
                command += ["@trg", str(out_prefix) + ".trg"]
        else:
            raise ValueError("Outprefix needs to be string got: "+type(out_prefix)+" - "+str(out_prefix))

        command_text = " ".join(command) + " > " + log_file_path + "\n"
        if verbose: print("COMMAND: ", command_text)

        start_time = datetime.datetime.now()
        process = bash.execute(command_text)
        md_run_return = process.poll()
        #bash.wait_for_fileSystem(out_prefix+".cnf")
        end_time = datetime.datetime.now()
        duration = end_time-start_time

        time.sleep(time_wait_s_for_filesystem)
        log_file = open(log_file_path, "a")
        log_file.write("\n\nREMARKS\n")

        failed = False
        if (md_run_return == 0):
            log_file.write("\tRUN:\tSUCESSFUL\n")
        else:
            log_file.write("\tRUN:\tFAILED\n")
            failed = True

        log_file.write("\tTIME:\n\tstart: " + str(start_time) + "\tend: " + str(end_time) + "\n")
        log_file.write("\tDuration:\t " + str(duration) + " d:hh:s.ms\n")
        log_file.write("END\n")
        log_file.close()

        if (md_run_return!=0):
            raise ChildProcessError("GromosXX MD Run Failed!\n\nLOG:" + "\n".join(open(log_file_path, "r").readlines()))

        return log_file_path

    @gromosTypeConverter
    def repex_run(self, in_topo_path: str, in_coord_path: str, in_imd_path: str, out_prefix: str,
                  in_pert_topo_path: str = None, in_disres_path: str = None, in_posresspec_path: bool = False, in_refpos_path: bool = False,
                  out_trc: bool = True, out_tre: bool = True, out_trs: bool = False, out_trg:bool = False,
                  nomp: int = 1, nmpi: int = 1, verbose: bool = True) -> str:
        """
        This function is a wrapper for gromosXX repex_mpi. You can directly execute the gromosXX repex_mpi in a bash enviroment here.

        Warnings
        --------
            Hybrid jobs are possible, but very difficult to implement correctly to Euler and performance gain is questionable.

        Parameters
        ----------
        in_topo_path :    str
                    This is the path to the input topology file (x.top)

        in_coord_path :    str
                    This is the path to the input coordinate file (x.cnf)

        in_imd_path :    str
                    This is the path to the input simulation parameter file (x.imd) - needs to contain a Replica Exchange or RE-EDS block.

        out_prefix :    str
                    This prefix, define the output name.

        in_pert_topo_path : str, optional
                    This is the path to the pertubation file (x.ptp)

        in_disres_path :    str, optional
                    This is the path to the distance restraint file (x.dat)

        in_posresspec_path :    str, optional
                    This is the path to the position restraint file (x.pos)

        in_refpos_path :    str, optional
                    This is the path to the reference position file (x.rpf)

        nomp :  int, optional
                    How many omp cores shall be used? Prerequesite, gromos was compiled with -enableOMP

        nmpi : int, optional
                    How many mpi cores shall be used? Prerequesite, gromos was compiled with -enableMPI (and suggested to have -disableOMP)

        out_trc :   bool, optional
                    do you want to output a coordinate trajectory (x.trc) file? (needs also an output number in write block of imd!)

        out_tre :   bool, optional
                    do you want to output a coordinate trajectory (x.trc) file? (needs also an output number in write block of imd!)

        out_trs :   bool, optional
                    do you want to output a coordinate trajectory (x.trc) file? (needs also an output number in write block of imd!)

        out_trg :   bool, optional
                    do you want to output the free energy trajectory (x.trg) file? (needs also an output number in write block of imd!)

        queueing_systems : NONE
            This var is not in use yet! - under development

        Returns
        -------
        str
            returns the log_file_path of the simulation.

        Raises
        ------
        ChildProcessError
            If the execution of the simulation fails, this error is raised.

        See Also
        --------
        md_run, md_mpi_run

        """

        if (nomp > 1):  # OMP Hybrid sopt_job?
            command = ["export OMP_NUM_THREADS=" + str(nomp) + " && mpirun -n " + str(nmpi) + " --loadbalance --cpus-per-proc " + str(nomp) + " "]
        else:
            command = ["mpirun", "-n ", str(nmpi)]

        command.append(self._bin + "repex_mpi")

        # input params check
        if in_topo_path:
            command += ["@topo", str(in_topo_path)]
        else:
            raise IOError("Did not get an input top file. Got: " + in_topo_path)

        if in_coord_path:
            command += ["@conf", str(in_coord_path)]
        else:
            raise IOError("Did not get an input coord file. Got: " + in_coord_path)

        if in_imd_path:
            command += ["@input", str(in_imd_path)]
        else:
            raise IOError("Did not get an input imd file. Got: " + in_imd_path)

        if in_pert_topo_path:
            command += ["@pttopo", str(in_pert_topo_path)]

        if in_disres_path:
            command += ["@distrest", str(in_disres_path)]

        if in_posresspec_path:
            command += ["@posresspec", str(in_posresspec_path)]

        if in_refpos_path:
            command += ["@refpos", str(in_refpos_path)]

        if out_prefix:
            command += ["@fin", str(out_prefix + ".cnf")]

            if out_trc:
                command += ["@trc", str(out_prefix + ".trc")]
            if out_trs:
                command += ["@trs", str(out_prefix + ".trs")]
            if out_tre:
                command += ["@tre", str(out_prefix + ".tre")]
            if out_trg:
                command += ["@trg", str(out_prefix + ".trg")]

            command += ["@repout", str(out_prefix + "_repout.dat")]
            command += ["@repdat", str(out_prefix + "_repdat.dat")]

        log_file_path = out_prefix + ".omd"
        command_text = " ".join(command) + " >> " + log_file_path + "\n"

        if verbose: print(command_text)

        os.system("echo \"\" >" + str(log_file_path))

        if verbose: start_time = time.ctime()
        if verbose: print("START: " + str(start_time))

        #bash.execute(command, verbose=True)
        md_run = os.system(command_text)
        time.sleep(time_wait_s_for_filesystem)
        if verbose: end_time = time.ctime()
        if verbose: print("END: " + str(end_time))
        print("MDRUN OUT: ", md_run)
        #if (md_run != 0):
        #    raise ChildProcessError("GromosXX REPEX Run Failed!\n \t return value: " + str(md_run) + "\n\tLOG:" + "\n\t\t".join(open(log_file_path, "r").readlines()))

        return log_file_path



class GromosXX(_Gromos):
    """
    GromosXX

    This is the class represents gromosXX.

    Attributes:
    -----------
    bin :   str, optional
        This is the path to the folder containing the binaries of gromosXX. If None, the bash enviroment variables  will be used.
    """

    def __init__(self, gromosXX_bin_dir: str = None):
        super().__init__(gromosXX_bin_dir=gromosXX_bin_dir)
