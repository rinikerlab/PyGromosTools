"""
FUNCTIONLIB:            wrapper for gromosXX
Description:
    This file contains python wrappers for the bash commandline of gromosXX

Author: Benjamin Schroeder
"""

import glob
import os
import datetime
import time
from pygromos.files.simulation_parameters.imd import Imd

from pygromos.utils import bash
from pygromos.utils.utils import time_wait_s_for_filesystem
from pygromos.gromos._gromosClass import _gromosClass


class _GromosXX(_gromosClass):
    """
    GromosXX

    This is the gromosXX baseclass. This should be inherited by a concrete class that might reimplement some new features, that are version dependent.

    Attributes:
    -----------
    bin :   str, optional
        This is the path to the folder containing the binaries of gromosXX. If None, the bash enviroment variables  will be used.
    """

    def __init__(
        self,
        gromosXX_bin_dir: str = None,
        _check_binary_paths: bool = True,
    ):
        """
        Constructing a gromosXX object.

        Parameters
        ----------
        gromosXX_bin_dir :   str, optional
            This is the path to the folder containing the binaries of gromosXX. If None, the bash enviroment variables  will be used.
        _dont_check_binary : bool, optional
            This flag removes the checks of the binary presence for this obj. This can make sense if system access is slow!, by default False - checks will be made
        """
        # lazy me - doc text for functions:
        functions_text = "\n    Methods:\n    ---------\n" + "\n".join(
            ["\t" + x for x in dir(self) if (not x.startswith("_") and callable(getattr(self, x)))]
        )
        self.__doc__ = self.__doc__ + functions_text

        super().__init__(
            in_bin_dir=gromosXX_bin_dir, _check_binary_paths=_check_binary_paths
        )  # initialises the binary checks

    def __str__(self):
        return self.__doc__

    def __repr__(self):
        return self.__str__()

    """
        GromosXX Programms
    """

    @_gromosClass._gromosTypeConverter
    def md_run(
        self,
        in_topo_path: str,
        in_coord_path: str,
        in_imd_path: str,
        out_prefix: str,
        in_pert_topo_path: str = None,
        in_disres_path: str = None,
        in_posresspec_path: str = None,
        in_refpos_path: str = None,
        in_qmmm_path: str = None,
        nomp: int = 1,
        nmpi: int = 1,
        out_trc: bool = False,
        out_tre: bool = False,
        out_trv: bool = False,
        out_trf: bool = False,
        out_trs: bool = False,
        out_trg: bool = False,
        verbose: bool = False,
        _binary_name: str = "md",
    ) -> str:
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

        # BUILD UP str. Command
        command = []
        if nmpi > 1 and nomp > 1:
            raise ValueError("There are no Hybrid NMPI and NOMP jobs possible with gromos!")
        elif nmpi > 1:
            command += ["mpirun -n " + str(nmpi * nomp) + " "]  # --loadbalance  " --cpus-per-proc " +  + " "
            command += [self._bin + _binary_name + "_mpi"]
        elif nomp >= 1:
            command += ["export OMP_NUM_THREADS=" + str(nomp) + "  && "]
            command += [self._bin + _binary_name]
        else:
            command += [self._bin + _binary_name]

        command += ["@topo", str(in_topo_path)]
        command += ["@conf", str(in_coord_path)]
        command += ["@input", str(in_imd_path)]

        if in_pert_topo_path is not None:
            command += ["@pttopo", str(in_pert_topo_path)]

        if in_disres_path is not None:
            command += ["@distrest", str(in_disres_path)]

        if in_posresspec_path is not None:
            print("POSRES", in_posresspec_path)
            command += ["@posresspec", str(in_posresspec_path)]

        if in_refpos_path is not None:
            command += ["@refpos", str(in_refpos_path)]

        if in_qmmm_path is not None:
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
            raise ValueError("Outprefix needs to be string got: " + type(out_prefix) + " - " + str(out_prefix))

        command_text = " ".join(command) + " > " + log_file_path + "\n"
        if verbose:
            print("COMMAND: ", command_text)

        start_time = datetime.datetime.now()
        process = bash.execute(command_text)
        md_run_return = process.poll()
        # bash.wait_for_fileSystem(out_prefix+".cnf")
        end_time = datetime.datetime.now()
        duration = end_time - start_time

        time.sleep(time_wait_s_for_filesystem)
        log_file = open(log_file_path, "a")
        log_file.write("\n\nREMARKS\n")

        failed = False
        if md_run_return == 0:
            log_file.write("\tRUN:\tSUCESSFUL\n")
        else:
            log_file.write("\tRUN:\tFAILED\n")

            omd_file_content = open(log_file_path, "r").read_lines()
            if len(omd_file_content) > 0:
                print("\t" + "\n\t".join(omd_file_content))
            else:
                print("\t None")
            failed = True

        log_file.write("\tTIME:\n\tstart: " + str(start_time) + "\tend: " + str(end_time) + "\n")
        log_file.write("\tDuration:\t " + str(duration) + " d:hh:s.ms\n")
        log_file.write("END\n")
        log_file.close()

        if failed:
            raise ChildProcessError("GromosXX MD Run Failed!\n\nLOG:" + "\n".join(open(log_file_path, "r").readlines()))

        return log_file_path

    @_gromosClass._gromosTypeConverter
    def repex_run(
        self,
        in_topo_path: str,
        in_coord_path: str,
        in_imd_path: str,
        out_prefix: str,
        in_pert_topo_path: str = None,
        in_disres_path: str = None,
        in_posresspec_path: bool = False,
        in_refpos_path: bool = False,
        out_trc: bool = True,
        out_tre: bool = True,
        out_trs: bool = False,
        out_trg: bool = False,
        out_trf: bool = False,
        out_trv: bool = False,
        nomp: int = 1,
        nmpi: int = 1,
        verbose: bool = True,
        _binary_name: str = "repex_mpi",
    ) -> str:
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

        out_trf :   bool, optional
                    do you want to output the free energy trajectory (x.trg) file? (needs also an output number in write block of imd!)

        out_trv :   bool, optional
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

        if nomp > 1:  # OMP Hybrid sopt_job?
            command = [
                "export OMP_NUM_THREADS="
                + str(nomp)
                + " && mpirun -n "
                + str(nmpi)
                + " --loadbalance --cpus-per-proc "
                + str(nomp)
                + " "
            ]
        else:
            command = ["mpirun", "-n ", str(nmpi)]

        command.append(self._bin + _binary_name)

        # input params check
        if in_topo_path:
            command += ["@topo", str(in_topo_path)]
        else:
            raise IOError("Did not get an input top file. Got: " + in_topo_path)

        if in_imd_path:
            command += ["@input", str(in_imd_path)]

            # Input cnf file depends if we have the CONT keyword or not
            # with CONT == 1, the convention is to give the name of the
            # file without the "_1" extension
            imd = Imd(in_imd_path)
            if hasattr(imd, "REPLICA") and imd.REPLICA is not None and imd.REPLICA.CONT:
                tmp_path = "/".join(os.path.abspath(in_coord_path).split("/")[:-1])
                in_coord_path = sorted(glob.glob(tmp_path + "/*.cnf"))[0]
                in_coord_path = in_coord_path.replace("_1.cnf", ".cnf")
            elif hasattr(imd, "REPLICA_EDS") and imd.REPLICA_EDS is not None and imd.REPLICA_EDS.CONT:
                tmp_path = "/".join(os.path.abspath(in_coord_path).split("/")[:-1])
                in_coord_path = sorted(glob.glob(tmp_path + "/*.cnf"))[0]
                in_coord_path = in_coord_path.replace("_1.cnf", ".cnf")
        else:
            raise IOError("Did not get an input imd file. Got: " + in_imd_path)

        if in_coord_path:
            command += ["@conf", str(in_coord_path)]
        else:
            raise IOError("Did not get an input coord file. Got: " + in_coord_path)

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
        if out_trf:
            command += ["@trf", str(out_prefix + ".trf")]
        if out_trv:
            command += ["@trv", str(out_prefix + ".trv")]

        command += ["@repout", str(out_prefix + "_repout.dat")]
        command += ["@repdat", str(out_prefix + "_repdat.dat")]

        log_file_path = out_prefix + ".omd"
        command_text = " ".join(command) + " >> " + log_file_path + "\n"

        if verbose:
            print(command_text)

        os.system('echo "" >' + str(log_file_path))

        if verbose:
            start_time = time.ctime()
        if verbose:
            print("START: " + str(start_time))

        # bash.execute(command, verbose=True)
        md_run = os.system(command_text)
        time.sleep(time_wait_s_for_filesystem)
        if verbose:
            end_time = time.ctime()
        if verbose:
            print("END: " + str(end_time))
        print("MDRUN OUT: ", md_run)
        # if (md_run != 0):
        #    raise ChildProcessError("GromosXX REPEX Run Failed!\n \t return value: " + str(md_run) + "\n\tLOG:" + "\n\t\t".join(open(log_file_path, "r").readlines()))

        return log_file_path


class GromosXX(_GromosXX):
    """
    GromosXX

    This is the class represents gromosXX.

    Attributes:
    -----------
    bin :   str, optional
        This is the path to the folder containing the binaries of gromosXX. If None, the bash enviroment variables  will be used.
    """

    def __init__(self, gromosXX_bin_dir: str = None, _check_binary_paths: bool = True):
        super().__init__(gromosXX_bin_dir=gromosXX_bin_dir, _check_binary_paths=_check_binary_paths)
