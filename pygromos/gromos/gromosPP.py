"""
FUNCTIONLIB:            wrapper for gromos++
Description:
    This file contains python wrappers for the bash commandline of gromos++
Author: Benjamin Schroeder
"""

import os
from typing import Union

import pandas as pd
from typing import List
from numbers import Number

from pygromos.utils import bash
from pygromos.data import pdb_lib
from pygromos.data.ff.Gromos2016H66 import ifp, mtb
from pygromos.gromos.gromosBashSyntaxParser import gromosBashSyntaxParser
from pygromos.gromos.utils import gromosTypeConverter



class _gromosPPbase:
    """
    GromosPP

    This is the gromosPP baseclass.
    This should be inherited by a concrete class that might reimplement some new features, that are version dependent.

    Attributes:
    -----------
    bin :   str, optional
        This is the path to the folder containing the binaries of gromosPP. If None, the bash enviroment variables  will be used.
    """

    _bin: str = ""
    _isValid:bool = False

    def __init__(self, gromosPP_bin_dir:str=None):
        """
        Constructing a gromosPP object.

        Parameters
        ----------
            bin :   str, optional
                This is the path to the folder containing the binaries of gromosXX. If None, the bash enviroment variables  will be used.
        """
        #lazy me - doc text for functions:
        functions_text = "\n    Methods:\n    ---------\n" + "\n".join(["\t\t" + x for x in dir(self) if (not x.startswith("_") and callable(getattr(self, x)))])
        self.__doc__ = self.__doc__+functions_text

        if(isinstance(gromosPP_bin_dir, str) and not "None" == gromosPP_bin_dir):
            self._bin = gromosPP_bin_dir + "/"
        else:
            self._bin = ""

        try:
            self.make_top(out_top_path=os.devnull, in_building_block_lib_path=mtb, in_parameter_lib_path=ifp, in_sequence="CH4")
            self._isValid = True
        except:
            self._isValid = False

    def __str__(self):
        return self.__doc__

    def __repr__(self):
        return self.__str__()

    @property
    def bin(self)->Union[str, None]:
        if(not hasattr(self, "_bin") or self._bin == "" ):
            return None
        else:
            return self._bin

    @gromosTypeConverter
    def pdb2gromos(self, in_pdb_path:str, in_top_path:str, out_cnf_path:str=None, in_lib_path:str=pdb_lib, _binary_name:str= "pdb2g96", verbose:bool = False)->str:
        """
        This is a wrapper for pdb2gromos. It executes the gromosPP binary.

        Parameters
        ----------
        in_pdb_path :    str
        in_top_path :    str
        in_lib_path :    str, optional
            The lib can be used as a look up table for residues or atoms etc. - see gromos manual
        out_cnf_path :    str, optional
            The out_cnf is the path for the output file. if not given, the pdb basename and dir is taken as default.
        _binary_name :  str, optional
        verbose : bool, optional

        Returns
        -------
        str
            Returns out_cnf path

        See Also
        --------
             For more information checkout the Gromos Manual
        """
        if(out_cnf_path==None):
            out_cnf_path = str(os.path.splitext(os.path.basename(in_pdb_path))[0]) + ".cnf"
        if(not out_cnf_path.endswith(".cnf")):
            out_cnf_path += ".cnf"


        command = self._bin + _binary_name + " @pdb " + in_pdb_path + " @topo " + in_top_path + " @out " + out_cnf_path
        if(in_lib_path!=None):
            command += " @lib " + in_lib_path
        command +="\n"

        if (verbose): print(command)
        p = bash.execute(command, catch_STD=True)
        if (verbose): print(p.stdout.readlines())

        return out_cnf_path

    @gromosTypeConverter
    def pdb2seq(self, in_pdb_path:str, out_path:str= "", pH:float=7.4, select:str = "ALL", gff:str = "54a7", add_head:str= "NH3+", add_tail:str= "COO-", _binary:str= "pdb2seq")->str:
        """
        This function is translating a pdb into a sequence file, that can be used to generate for example topologies.

        Parameters
        ----------
        in_pdb_path :   str

        out_path :   str, optional
            out_path for output File

        pH :    float, optional

        select :    str, optional

        gff :   str,optional
            which Gromos ForceField

        add_head :  str, optional
            protein N-term capping with this group.

        add_tail :  str,optional
            protein C-term capping with this group

        _binary :   str, optional

        Returns
        -------
        str
            out_path

        See Also
        --------
             For more information checkout the Gromos Manual
        """
        if(out_path== ""):
            out_path = os.path.dirname(in_pdb_path) + "/" + str(os.path.splitext(os.path.basename(in_pdb_path))[0]) + ".seq"
        command = self._bin + _binary + " @develop @pdb " + in_pdb_path + " @pH " + str(pH) + " @select " + select + " @gff " + gff
        if(add_head != ""):
            command += " @head "+add_head
        if(add_tail != ""):
            command += " @tail "+add_tail
        command += " > " + out_path + " \n"
        bash.execute(command)
        return out_path

    @gromosTypeConverter
    def make_top(self, out_top_path:str, in_building_block_lib_path: str, in_parameter_lib_path: str, in_sequence:str, in_solvent:str= "H2O", additional_options="\n", _binary_name:str= "make_top")->str:
        """
        This wrapper uses make_top to generate a topology file.

        Parameters
        ----------
        out_top_path :   str
        in_building_block_lib_path : str
        in_parameter_lib_path :  str
        in_sequence :   str
        in_solvent :    str, optional
        additional_options :    str, optional

        Returns
        -------
        str
            returns out_file_path

        See Also
        --------
             For more information checkout the Gromos Manual
        """

        if(os.path.exists(in_sequence)):
           seq_file = open(in_sequence, "r")
           in_sequence = "".join(seq_file.readlines())

        args = ["@build " + in_building_block_lib_path,
                "@param " + in_parameter_lib_path,
                "@seq " + in_sequence,
                "@solv " + in_solvent]
        #arg_path = os.path.dirname(out_top_path) + "/topargs.arg"
        #arg_file = open(arg_path, "w")
        #arg_file.write("\n".join(args))
        #arg_file.close()
        #"@f " + arg_path
        command = self._bin +_binary_name +" "+" ".join(args)
        bash.execute(command, catch_STD=out_top_path)
        return out_top_path

    @gromosTypeConverter
    def com_top(self, in_topo_paths:(str or List[str]), topo_multiplier:(int or List[int])=1, out_top_path:str= "combined_out.top",
                take_topology_params_of_file:int=1, take_solvent_parameters_of_file:int=1,
                _binary_name:str="com_top")->str: #Todo: also take lists as input ! bschroed
        """
        Combine multiple topologies
        Parameters
        ----------
        in_topo_paths : (str or List[str])
        out_top_path : str, optional
        take_topology_params_of_file :  int, optional
        take_solvent_parameters_of_file :   int, optional
        _binary_name :  str, optional

        Returns
        -------
        str
            output file_path

        See Also
        --------
        make_top
             For more information checkout the Gromos Manual
        """
        if(len(in_topo_paths)==1 and type(in_topo_paths[0])==list):
            in_topo_paths = in_topo_paths[0]
        if(not out_top_path.endswith(".top")):
            out_top_path+= ".top"

        if(isinstance(in_topo_paths, list) and isinstance(topo_multiplier, int)):
            topo_multiplier = [topo_multiplier for x in range(len(in_topo_paths))]

        topo_argument = gromosBashSyntaxParser.multiplyArgumentParser(in_topo_paths, topo_multiplier)

        command = self._bin + _binary_name + " @topo " + topo_argument + " @param " + str(take_topology_params_of_file) + " @solv " + str(take_solvent_parameters_of_file)
        bash.execute(command, catch_STD=out_top_path)
        return out_top_path

    @gromosTypeConverter
    def dfmult(self, in_endstate_file_paths: List[str], in_reference_state_file_path:str, out_file_path:str= "dfmult_temp.out",
               temperature:float=298, _binary_name:str= "dfmult", verbose:bool=False)->str:
        """
            This funciton wraps dfmult of the gromos suite, it is used to calculate the free energy of a EDS simulation.
                .. math:: \frac{\langle e^{V_{i}-V_{R} \rangle_{R}}{\langle e^{V_{i}-V_{R} \rangle_{R}}

        Parameters
        ----------
        in_endstate_file_paths :    List[str]
             potentials energy of single state (get with ene Ana)
        in_reference_state_file_path :  str
                     potential energy of reference State (get with ene Ana)
        out_file_path :   str, optional
        temperature :   float, optional
        _binary_name :  str, optional
        verbose :   bool, optional

        Returns
        -------
        str
            out_put file path

        See Also
        --------
             For more information checkout the Gromos Manual
        """

        #formulate command
        if(verbose): print(" ".join(in_endstate_file_paths))
        command = self._bin + _binary_name + " @stateR " + in_reference_state_file_path + " @temp " + str(temperature) + " @endstates " + " ".join(in_endstate_file_paths) + " > " + out_file_path + "\n"
        if(verbose): print(command)
        #do
        ret = bash.execute_os(command,  verbose=verbose)
        if(verbose): print(ret.readlines())
        bash.wait_for_fileSystem(out_file_path)

        return out_file_path

    @gromosTypeConverter
    def frameout(self, in_top_path:str, in_coord_path:str, periodic_boundary_condition:str,
                 out_file_path:str=None, out_file_format:str=None, single_file:bool=None,
                 gather:(int or str)=None, include:str= "SOLUTE",
                 reference_structure_path:str=None, atomsfit:str=None,
                 frames:int=None, time:int=None, dt:int=None, notimeblock:bool=None,
                 _binary_name:str= "frameout", verbose:bool=False)->str:
        """
            this wrapper wraps frameout.
                frameout is a tool, that can be used for example to convert gromos coordinate files or to recenter them or ....
                Optional parameters are not used in the command, when they are None!

        Parameters
        ----------
        in_top_path :   str
        in_coord_path : str
        periodic_boundary_condition :   str
        out_file_path : None or str, optional
        out_file_format :   None or str, optional
        single_file :   None or bool, optional
        gather :    None or int or str optional

        include :   None or str, optional
            ALL that also includes the solvent.
            SOLUTE this option filters the Solvent out of the output file. (default)
            alternative gromos selection synthax can be used.

        reference_structure_path :  None or str, optional
            This path should provide a refrence position for atom fitting.

        atomsfit :  None or str, optional
            This option can be used to fit all frames to one reference structure.
            The selection Syntax follow the Gromos manual (e.g. "1:a" - aligns all frames to the first molecule and all its atoms)
            requires reference_sturcture_path.

        frames :    None or int, optional
        time :  None or float, optional
        dt :    None or float, optional
        notimeblock :   None or bool, optional

        _binary_name :  str, optional
        verbose :   bool, optional

        Returns
        -------
        str
            out_file_path

        See Also
        --------
             For more information checkout the Gromos Manual
        """

        options = ""
        if(out_file_format != None):
            options += "@outformat " + str(out_file_format) + " "

        periodic_boundary_condition = "@pbc " + periodic_boundary_condition + " "
        if (gather != None):
            periodic_boundary_condition += " " + str(gather) + " "
        if(include!=None):
            options+="@include "+str(include)+" "
        if(single_file!=None):
            options += "@single "
        if(atomsfit!=None and reference_structure_path!=None):
            options += "@ref " + str(reference_structure_path) + " @atomsfit " + atomsfit + " "
        if(not isinstance(time, type(None))):
            options += "@time "+str(time)+" "
        if(not isinstance(time, type(None)) and not isinstance(dt, type(None))):
            options += " "+str(dt)+" "
        if(notimeblock):
            options += "@time "+str(0)

        if(frames !=None ):
            raise NotImplementedError("Chosen Options for frameout not implemented yet!")

        if(out_file_path!=None):
            orig =os.getcwd()+"/FRAME_00001."+out_file_format

        else:
            orig =os.getcwd()+"/FRAME* "
            if(not isinstance(out_file_format, type(None))):
                out_file_path = os.path.dirname(in_coord_path) + "/FRAME_00001."+out_file_format
            else:
                out_file_path = os.path.dirname(in_coord_path) + "/FRAME_00001."+in_coord_path.split(".")[-1]

        command = self._bin + _binary_name + " @topo " + str(in_top_path) + " @traj " + str(in_coord_path) + " " + periodic_boundary_condition + " " + options
        if verbose: print("gromosPP.frameout: command:\n" +command)

        #DO
        try:
            ret = bash.execute(command, catch_STD=True)
            if(verbose):    print("STDOUT: ", "\n".join(ret.stdout.readlines()))
        except Exception as err:
            print("gromosPP.frameout: could not exectue framout:\n"+str(err.args))
            raise Exception("gromosPP.frameout: could not exectue framout:\n"+str(err.args))


        bash.wait_for_fileSystem(check_paths=orig, regex_mode=True)


        #move to output
        try:
            bash.move_file(orig, out_file_path)
        except Exception as err:
            print("gromosPP.frameout: could not move new FRAME*.pdb file to out-parameter:\n"+str(err.args))
            raise Exception("gromosPP.frameout: could not move new FRAME*.pdb file to out-parameter:\n" + str(err.args))

        return out_file_path

    @gromosTypeConverter
    def ene_ana(self, in_ene_ana_library_path: str, in_en_file_paths: str, in_properties: str,
                out_energy_folder_path: str, out_files_prefix: str = None, out_files_suffix: str = None,
                in_topo: str = None, time: float = None, single_file: bool = False,return_outstat_also:bool=False,
                verbose: bool = False, _binary_name: str = "ene_ana", workdir: bool = False) -> (List[str] or str):

        """
            This is a wrapper for ene_ana.
                ene_ana is a tool to extract energy properties from a x.tre file.

        Parameters
        ----------
        in_ene_ana_library_path :    str
        in_en_file_paths :   str
        in_properties : str
            this str should name the properties, that shall be written out. For the variable names, check the ene_ana lib.
                (e.g.: "solvtemp1 e1 eR")
        out_energy_folder_path : str
            give the path to an directory, where the output should be stored in.
        out_files_prefix :  None or str, optional
        out_files_suffix :  None or str, optional
        in_topo :   None or str, optional
        time :  None or float, optional
        single_file: bool, optional
            if used a single csv is generated. the return value gets the file_path(str).
        verbose :   bool, optional
        _binary_name :  str, optional

        Returns
        -------
        List[str] or str
            if single_file used, return = str

        """

        # check input
        if (type(in_properties) == list):
            in_properties = " ".join(in_properties)
        elif (type(in_properties) == str):
            pass
        else:
            raise IOError(
                "gromosPP.ene_ana:  got an format for potentials, that is unknown! (Please give a list of strings or a final string)\n" + str(
                    in_properties))

        if (isinstance(in_en_file_paths, list)):
            in_en_file_paths = " ".join(in_en_file_paths)
        elif (isinstance(in_en_file_paths, str)):
            pass
        else:
            raise IOError(
                "gromosPP.ene_ana:  got an format for potentials, that is unknown! (Please give a list of strings or a final string)\n" + str(
                    in_en_file_paths))

        # manage nice prefix
        if (out_files_prefix):
            prefix = out_files_prefix
        else:
            prefix = ""

        # manage nice suffix
        if (out_files_suffix):
            if (not out_files_suffix.startswith("_") and not out_files_prefix.endswith("_")):
                suffix = "_" + out_files_suffix
            else:
                suffix = out_files_suffix
        else:
            suffix = ""

        additional_options = ""

        if (in_topo):
            additional_options += " @topo " + str(in_topo)
        if (not isinstance(time, type(None))):
            additional_options += " @time " + str(time)

        original_pos = os.getcwd()
        if (not workdir):
            os.chdir(out_energy_folder_path)

        if (in_en_file_paths.strip().endswith("trg") or in_en_file_paths.strip().endswith("trg.gz")):
            in_file_form = "@fr_files"
        else:
            in_file_form = "@en_files"

        ##ene_ana command and log deletion if ene_ana was successfull.:
        command = self._bin + _binary_name + " @library " + str(in_ene_ana_library_path) + " " + in_file_form + " " + str(
            in_en_file_paths) + " @prop " + in_properties + " " + additional_options

        if (verbose): print("Isolate_energies")
        try:
            out_fun  = bash.execute(command, verbose=verbose)
        except Exception as err:
            raise Exception("gromosPP.ene_ana: could not read out energies:\n" + "\n".join(err.args) + "\n command used: " + command)

        # Wait for file system.
        tmp_files = []
        for prop in in_properties.split():
            check_path = str(os.getcwd() + "/" + prop + ".dat")
            bash.wait_for_fileSystem(check_path, verbose=verbose)
            tmp_files.append(check_path)

        if (not single_file):
            try:
                # move results to outfolder and give them proper name
                result_files = []
                if (type(in_properties) == str):
                    in_properties = in_properties.split()

                for prop in in_properties:
                    orig_file = str(os.getcwd() + "/" + prop + ".dat")
                    target_file = str(out_energy_folder_path + "/" + prefix + prop + suffix + ".dat")

                    if (target_file == orig_file):
                        pass
                    else:
                        bash.move_file(in_file_path=orig_file, out_file_path=target_file)

                    result_files.append(target_file)
            except Exception as err:
                raise Exception("gromosPP.ene_ana: could not move and rename Files:\n" + "\n".join(err.args) + "\n before use: " + command)

            bash.wait_for_fileSystem(result_files)
        else:

            if (verbose): print("reading in tmp_files files:\t" + "\n\t".join(tmp_files))
            # energy_properties = [pd.read_csv(in_ene_traj_path, header=0, delim_whitespace=True) for in_ene_traj_path in tmp_files]

            # fix columns
            first = True
            for in_ene_traj_path in tmp_files:
                energy_property = pd.read_csv(in_ene_traj_path, header=0, delim_whitespace=True)
                new_cols = [x.replace("#", "").strip() for x in energy_property.columns if (not x == "#")] + [""]
                new_cols = [x if ("time" in x) else x for x in new_cols]
                energy_property.columns = new_cols
                energy_property.pop("")
                if (verbose): print(energy_property.columns)
                if (verbose): print(energy_property.shape)

                if (first):
                    first = False
                    concat_energy_traj = energy_property
                else:
                    energy_property.pop("time")
                    concat_energy_traj = pd.concat([concat_energy_traj, energy_property], axis=1)
                del energy_property

            tmp_pandas_out = out_energy_folder_path + "/" + prefix + suffix + ".dat"
            concat_energy_traj.to_csv(tmp_pandas_out, sep="\t", header=True, index=False)
            del concat_energy_traj

            #   remove old trajs
            for x in tmp_files:
                if (os.path.exists(x)):
                    bash.remove_file(x)
            result_files = tmp_pandas_out

        if (not workdir):
            os.chdir(original_pos)

        bash.wait_for_fileSystem(result_files)

        if(return_outstat_also):
            return result_files, out_fun
        else:
            return result_files

    @gromosTypeConverter
    def protonate(self, in_cnf_path:str, in_top_path:str, out_cnf_path:str,
                  tolerance:float=0.1, periodic_boundary_condition:str="v", gathering:str="cog", _binary_name="gch") -> str:
        """
                    This function protonates a coordinate file.

        Parameters
        ----------
        in_cnf_path :   str
        in_top_path :  str
        out_cnf_path :  str
        tolerance : float, optional
        periodic_boundary_condition : str, optional
        gathering : str, optional
        _binary_name :    str, optional

        Returns
        -------
            out_cnf_path

        """
        command = self._bin + _binary_name + " @topo " + in_top_path + " @pos " + in_cnf_path + " @tol " + str(tolerance) + "  " \
                                  "@pbc " + periodic_boundary_condition +" " + gathering

        bash.execute(command, catch_STD=out_cnf_path)
        
        return out_cnf_path

    @gromosTypeConverter
    def sim_box(self, in_top_path:str, in_cnf_path:str, in_solvent_cnf_file_path:str, out_cnf_path:str= "",
                periodic_boundary_condition: str = "r", gathering_method:str=None, minwall:float=0.8, threshold:float=None, rotate:str=None,
                boxsize:float=None, _binary_name="sim_box", verbose=False)->str:
        """
            This wrapper wraps sim_box programm of gromos. It can be used to set a system box and solvate the box.

        Parameters
        ----------
        in_top_path :   str
        in_cnf_path :   str
        in_solvent_cnf_file_path :   str
        periodic_boundary_condition :   str, optional
        out_cnf_path :  str, optional
        minwall :   float, optional
            minimal box wall distance to solute molecule.
        threshold :   float, optional
        rotate :   float, optional
        boxsize :   float, optional
        gathering_method :    str, optional
        _binary_name :  str, optional

        Returns
        -------
        str
            out_cnf_path
        """
        command_suffix = ""
        if(out_cnf_path== ""):
            out_cnf_path = os.path.dirname(in_cnf_path) + "/" + str(os.path.splitext(os.path.basename(in_cnf_path))[0]) + "_solvent.cnf"
        if(rotate!=None):
            command_suffix+= " @rotate "
        if(gathering_method!=None):
            command_suffix+= " @gather " + str(gathering_method)
        if(boxsize!=None):
            command_suffix+= " @boxsize "+str(boxsize)
        if(threshold!=None):
            command_suffix+= " @thresh "+str(threshold)
        if(minwall!=None):
            command_suffix+= " @minwall " + str(minwall)

        command= self._bin + _binary_name + " @topo " + in_top_path + " @pbc " + periodic_boundary_condition + " @pos " + in_cnf_path + " @solvent " + in_solvent_cnf_file_path + " " + command_suffix
        print(command)
        p = bash.execute(command, verbose=verbose, catch_STD=out_cnf_path)
        return out_cnf_path

    @gromosTypeConverter
    def ran_box(self, in_top_path:str, in_cnf_path:str, out_cnf_path:str= "",
                periodic_boundary_condition: str = "r", nmolecule:int = 1, dens:float = 1.0, threshold:float=None, layer:bool = False, 
                boxsize:float=None, fixfirst:bool = False, seed:float=None, _binary_name="ran_box", verbose=False, return_command_only=False)->str:

        command_suffix = ""
        if(out_cnf_path== ""):
            out_cnf_path = os.path.dirname(in_cnf_path) + "/" + str(os.path.splitext(os.path.basename(in_cnf_path))[0]) + "_ran-box.cnf"
        if(threshold!=None):
            command_suffix+= " @thresh " + str(threshold)
        if layer:
            command_suffix += " @layer "
        if(boxsize!=None):
            command_suffix+= " @boxsize "+ str(boxsize)
        if fixfirst:
            command_suffix += " @fixfirst "
        if(seed!=None):
            command_suffix+= " @seed " + str(seed)


        command= self._bin + _binary_name + " @topo " + in_top_path + " @pbc " + periodic_boundary_condition + " @pos " + in_cnf_path + " @nsm " + str(nmolecule) + " @dens " + str(dens) + " " + command_suffix + " > " + out_cnf_path + " \n"
        if not return_command_only:
            print(command)
            std_out = bash.execute(command, verbose=verbose)
            return out_cnf_path
        else:
            return command

    @gromosTypeConverter
    def build_box(self, in_top_path:str, in_cnf_path:str, out_cnf_path:str= "",
                periodic_boundary_condition: str = "r", nmolecule:int = 1, dens:float = 1.0, _binary_name="build_box", verbose=False, return_command_only=False)->str:

        if(out_cnf_path== ""):
            out_cnf_path = os.path.dirname(in_cnf_path) + "/" + str(os.path.splitext(os.path.basename(in_cnf_path))[0]) + "_ran-box.cnf"

        command= self._bin + _binary_name + " @topo " + in_top_path + " @pos " + in_cnf_path + " @nsm " + str(nmolecule) + " @dens " + str(dens) + " > " + out_cnf_path + " \n"
        if not return_command_only:
            print(command)
            std_out = bash.execute(command, verbose=verbose)
            return out_cnf_path
        else:
            return command

    @gromosTypeConverter
    def tser(self, in_trc_path:str, in_top_path:str, out_csv_path:str, property:str,
             periodic_boundary_condition:str= "r", time:float = None, solvent:str=None,
             normalise_distribution:bool=False, skip_first_n_frames:int=0, take_each_nth_frame:int=1,
             _binary_name:str= "tser")->str:
        """
                    Tser is a gromos programm, that can analyze trajectories.

        Parameters
        ----------
        in_trc_path :   str
        in_top_path :   str
        out_csv_path :   str
        property :   str
        periodic_boundary_condition :   str, optional
        time :   float, optional
        solvent :   str, optional
        normalise_distribution :   bool, optional
        skip_first_n_frames :   int, optional
        take_each_nth_frame :   int, optional
        _binary_name :   str, optional

        Warnings
        --------
            missing options: @nots, @dist

        Returns
        -------
        str
            out_csv_path
        """
        optional_string = ""
        if(take_each_nth_frame > 1):
            optional_string += " @stride "+str(take_each_nth_frame)+" "
        if(skip_first_n_frames>0):
            optional_string += " @skip "+str(skip_first_n_frames)+" "
        if(normalise_distribution):
            optional_string += " @norm "
        if(isinstance(time, type(None))):
            optional_string += " @time "+str(time)+" "
        if(isinstance(solvent, type(None))):
            optional_string += " @solv "+str(solvent)+" "

        command = self._bin + _binary_name + " @topo " + in_top_path + " @pbc " + periodic_boundary_condition + " @traj " + in_trc_path + " @prop \"" + str(property) + "\" " + optional_string + " > " + out_csv_path +" \n"
        bash.execute(command)
        return out_csv_path

    @gromosTypeConverter
    def red_top(self, in_top_path:str, atom_selection:str, out_top_path:str, _binary_name:str= "red_top")->str:
        """
            red_top is a gromos tool to reduce a gromos tool to a certain selection.
        Parameters
        ----------
        in_top_path :   str
        atom_selection :    str
        out_top_path :  str
        _binary_name :  str, optional

        Returns
        -------
        str
            out_top_path
        """

        command =  [self._bin + _binary_name, " @topo ", in_top_path, " @atoms", "'" + atom_selection + "'", "> " + out_top_path + " \n"]
        bash.execute(command)
        return out_top_path

    @gromosTypeConverter
    def prep_eds(self, in_top_paths:List[str], number_of_eds_states:int, param_top_index:int=1, solv_top_index:int=1,
                 out_file_path:str= "dev", _binary_name:str= "prep_eds", verbose:bool=False)->(str, str):
        """
            prepare eds topology.

        Parameters
        ----------
        in_top_paths : List[str]
        number_of_eds_states :   int
        param_top_index :   int, optional
        solv_top_index :   int, optional
        out_file_path :    str, optional
            output path without file ending. (prefix)
        _binary_name :  str, optional
        verbose :   bool, optional

        Returns
        -------
        tuple[str,str]
            out_top, out_ptp
        """
        out_file_path = os.path.splitext(out_file_path)[0]
        if(type(in_top_paths) == str):
            in_top_paths = [in_top_paths]

        if(len(in_top_paths) == 0):
            raise ValueError("no topos were passed to function. please provide at least two")

        command = self._bin + _binary_name + " @topo " + " ".join(in_top_paths) + " @numstat " + str(number_of_eds_states) + " @param " + str(param_top_index) + " @solv " + str(solv_top_index)
        if(verbose):
            print(command)
        ret = bash.execute(command)
        if(verbose):
            print(ret.readlines())

        bash.wait_for_fileSystem(["./pert_eds.ptp", "./com_eds.top"])

        out_ptp = bash.move_file("./pert_eds.ptp", out_file_path + ".ptp")
        out_top = bash.move_file("./com_eds.top", out_file_path + ".top")

        return out_top, out_ptp

    @gromosTypeConverter
    def prep_noe(self, in_top_path:str, in_noe_path:str, in_library_path:str, out_path:str,
                 dish:float=0.1, disc:float=0.153,
                 title:str="NOE", _binary_name:str= "prep_noe",
                 in_correction_path: str = None,
                 verbose:bool=False)->str:
        """

        Parameters
        ----------
        in_top_path: str
            molecular topology file
        in_noe_path: str
            NOE specification file
        in_library_path: str
            NOE specification library
        out_path: str
            path to the output file
        dish: float, optional
            carbon-hydrogen distance; default: 0.1 nm
        disc: float, optional
            carbon-carbon distance; default: 0.153 nm
        title: str, optional
            NOE title for output, default: "NOE"
        correction: str, optional
            Correction file -> <correction_file> [correction type]

        Returns
        -------
        str
            path to the output file

        NotImplemented
        -------
        parsetype: <1,2,3>
         Choices are:
	        1: Upper bound == first number
	        2: Upper bound == first + third number (most common, default)
	        3: Upper bound == first - second number (commonly the lower bound)
        action: <add> or <substraction> = add
        filter: discard nNOE's above a certain distance[nm] = 10000 nm
        factor: conversion factor ang to nm , = 10

        """

        additional_options = ""
        if(isinstance(in_correction_path, str)):
            additional_options += "@correction "+in_correction_path+" "

        command = self._bin + _binary_name + " @topo " + in_top_path + " @title " + title +" @noe " + in_noe_path +" @lib " + in_library_path + " " \
                 " @dish " + str(dish) +" @disc " + str(disc) +" " + additional_options +" &> " + out_path

        if (verbose): print(command)
        ret = bash.execute(command)
        if (verbose): print(ret.readlines())

        return out_path

    @gromosTypeConverter
    def rmsf(self, in_top_path:str, in_trcs:Union[str, List[str]], atom_selection:str, out_file_path:str, pbc: str= "r", _binary_name:str= "rmsf")->str:
        """
            this is a wrapper for gromosPP rmsf programm. (Root mean square fluctuation

        Parameters
        ----------
        in_top_path : str
            path to topology file
        in_trcs : Union[str, List[str]]
            Path OR paths to trc coordinate files
        atom_selection : str
            selection of atoms
        out_file_path :
            out path.
        pbc : str
            periodic boundary condition of trc files
        _binary_name: str
            binary name of gromos file

        Returns
        -------
        str
            outpath of the traj
        """

        if(isinstance(in_trcs, list)):
            in_trcs = " ".join(in_trcs)

        additional_options = [""]
        additional_options = " ".join(map(str, additional_options))

        command = " ".join([self._bin + _binary_name, " @topo ", in_top_path, " @atomsrmsf", "'" + atom_selection + "'", "@pbc", pbc, "@traj", in_trcs,  additional_options, " \n"])
        bash.execute(command, catch_STD=out_file_path)
        return out_file_path

    @gromosTypeConverter
    def rmsd(self, in_top_path:str, in_trcs:Union[str, List[str]], atom_selection:str, out_file_path:str, pbc: str= "r", _binary_name:str= "rmsd")->str:
        """

        Parameters
        ----------
        in_top_path
        in_trcs
        atom_selection
        out_file_path
        pbc
        _binary_name

        Returns
        -------

        """

        if(isinstance(in_trcs, list)):
            in_trcs = " ".join(in_trcs)

        additional_options = [""]
        additional_options = " ".join(map(str, additional_options))
        command = " ".join([self._bin + _binary_name, " @topo ", in_top_path, " @atomsrmsd", "'" + atom_selection + "'", "@pbc", pbc, "@traj", in_trcs,  additional_options, " \n"])
        bash.execute(command, catch_STD=out_file_path)
        return out_file_path

    @gromosTypeConverter
    def cog(self, in_top_path:str, in_trcs:Union[str, List[str]], out_file_path:str,
            atom_selection:str=None, outformat:str=None, cog_com:str=None,
            add_repl:str=None, solv:str=None, nthframe:str=None,
            pbc: str="r cog", _binary_name:str= "cog", verbose:bool=False)->str:
        """

        
        Parameters
        ----------
        in_top_path : str
            path to topo file
        in_trcs : str
            path to trc files
        out_file_path : str
            outpath
        atom_selection: str, optional
            atomspecifier(s) for which to calculate cog/com
        outformat : str, optional
            output coordinates format
        cog_com : str, optional
            calculate centre of geometry (cog) or mass (com); default: cog
        add_repl: str, optional
            add (add) the cog/com or replace (repl) the solutes; default: repl
        solv: str, optional
            include solvent in outcoordinates
        nthframe: str, optional
            write every nth frame (default: 1)
        _binary_name : str, otpional
            binary name of gromosPP programm (default: cog)

        Returns
        -------
        str
            output_path of the generated csv file
        """

        if (isinstance(in_trcs, list)):
            in_trcs = " ".join(in_trcs)

        additional_options = [""]
        if(not atom_selection is None):
            additional_options+= [" @atomspec ", "'" + atom_selection + "'"]
        if(not outformat is None):
            additional_options+= [" @outformat ", outformat ]
        if(not cog_com is None):
            additional_options+= [" @cog_com ", cog_com ]
        if(not add_repl is None):
            additional_options += [" @add_repl ", add_repl]
        if(not solv is None):
            additional_options += [" @solv ", solv]
        if(not nthframe is None):
            additional_options += [" @nthframe ", nthframe]

        additional_options = " ".join(map(str, additional_options))

        command = " ".join([self._bin + _binary_name, " @topo ", in_top_path, "@pbc", pbc, "@traj", in_trcs, additional_options, " \n"])

        if(verbose): print(command)
        bash.execute(command, catch_STD=out_file_path, verbose=verbose)

        return out_file_path

    @gromosTypeConverter
    def noe(self,  in_top_path:str, in_noe_path:str, in_traj_path:str, out_path:str,
                pbc:str="v", gathering:str="cog",
            _binary_name:str= "noe", verbose:bool=False) -> str:
        """

        Parameters
        ----------
        in_top_path: str
            topology path
        in_noe_path: str
            output path of prep_noe
        in_traj_path: str
            coordinate file
        out_path: str
            path to the output file
        pbc: str, optional
            default: v
            periodic boundary condition of the coordinates:
                v - vacuum
                r - rectangular box
        gathering: str, optional
            default: cog
            how the coordinates shall be gathered before the calculation:
                cog - center of geometry
                com - center of mass

        Returns
        -------
        str
            out_path

        NotImplemented
        ---------------
        time: float, float (time, dt)
        """
        command =  self._bin + _binary_name + " @topo " + in_top_path + " @traj "+in_traj_path+" @noe "+in_noe_path+" @pbc "+str(pbc)+" "+str(gathering)+"\n"

        if (verbose): print(command)
        p = bash.execute(command, catch_STD=out_path, verbose=verbose)

        return out_path

    @gromosTypeConverter
    def jval(self,  in_top_path:str, in_jval_path:str, in_traj_path:(str, List[str]), out_path:str,
              pbc:str="v", gathering:str="cog",
              timeseries:bool=False, rmsd:bool = False, time:float = None,
              _binary_name:str= "jval", verbose:bool=False):
        """

        Parameters
        ----------
        in_top_path: str
            topology path
        in_jval_path: str
            jval specification path
        in_traj_path: str
            coordinate file
        out_path: str
            path to the output file
        pbc: str, optional
            default: v
            periodic boundary condition of the coordinates:
                v - vacuum
                r - rectangular box
        gathering: str, optional
            default: cog
            how the coordinates shall be gathered before the calculation:
                cog - center of geometry
                com - center of mass
        timeseries: bool, optional

        rmsd: bool, optional

        time: (Number, str), optional


        Returns
        -------
        str
            out_path

        NotImplemented
        ---------------
        time: float, float (time, dt)
        """

        additional_options= ""
        if(timeseries):
            additional_options += " @timeseries "
        if(rmsd):
            additional_options += " @rmsd "
        if(isinstance(time, (Number, str))):
            additional_options += " @time "+str(time)+" "

        if(isinstance(in_traj_path, List)):
            in_traj_path = "\n".join(in_traj_path)

        command =  self._bin + _binary_name + " @topo " + in_top_path + " @traj "+in_traj_path+" @jval "+in_jval_path+" @pbc "+str(pbc)+" "+str(gathering)+" "+additional_options+" &> "+out_path

        if (verbose): print(command)
        ret = bash.execute(command)
        if (verbose): print(ret.readlines())

        return out_path

    @gromosTypeConverter
    def ion(self, in_top_path:str,in_cnf_path:str, out_cnf_path:str,
                periodic_boundary_condition:str="v",
            negative:list=None, positive:list=None,
            potential:float=0.8, mindist:float=0.8, _binary_name = "ion", verbose:bool=False
            ):
        optional_args = []
        if(not positive is None):
            opt = "@positive "+" ".join(map(str, positive))
            optional_args.append(opt)

        if(not negative is None):
            opt = "@negative "+" ".join(map(str, negative))
            optional_args.append(opt)

        command = self._bin + _binary_name + " @topo " + in_top_path + " @pos " + in_cnf_path + " @pbc " + periodic_boundary_condition+" "
        command+= "@potential "+str(potential)+" @mindist "+str(mindist)+" "+" ".join(optional_args)

        if(verbose): print(command)
        bash.execute(command, catch_STD=out_cnf_path)

        return out_cnf_path

    #To implement
    def _gr962pdb(self):
        raise Exception('not implemented yet!')

    @gromosTypeConverter
    def rgyr(self, out_rgyr_path:str, in_coord_path:str, in_top_path:str, atom_selection:str, periodic_boundary_condition:str="r cog", time:int=None, dt:int=None,  
    mass_weighted:bool=False, _binary_name:str= "rgyr")->str:
        """
        This wrapper uses rgyr to compute the radius of gyration for a given atom selection. 

        Parameters
        ----------
        out_rgyr_path: str
        in_coord_path: str
        in_top_path: str
        atom_selection: str
        periodic_boundary_condition: str
        time: int, optional
        dt: int, optional
	mass_weighted:bool, optional
	_binary_name: str, optional

        Returns
        -------
        str
            returns out_rgyr_path

        See Also
        --------
             For more information checkout the Gromos Manual
        """

        args = ["@topo " + in_top_path,
                "@pbc " + periodic_boundary_condition,
                "@atoms " + atom_selection,
                "@traj " + in_coord_path]

        if(not isinstance(time, type(None))):
            args += "@time "+str(time)+" "
        if(not isinstance(time, type(None)) and not isinstance(dt, type(None))):
            args += " "+str(dt)+" "
	
        command = self._bin +_binary_name+" "+" ".join(args)
        bash.execute(command, catch_STD=out_rgyr_path)
        return out_rgyr_path


class GromosPP(_gromosPPbase):
    """
    GromosPP

    This is the class represents gromosPP.

    Attributes:
    -----------
    bin :   str, optional
        This is the path to the folder containing the binaries of gromosPP If None, the bash enviroment variables  will be used.
    """

    def __init__(self, gromosPP_bin_dir: str =None):
        super().__init__(gromosPP_bin_dir=gromosPP_bin_dir)
