"""
File:            gromos++ system folder managment
Warnings: this class is WIP!

Description:
    This is a super class for the collection of all gromos files used inside a project
Author: Marc Lehner, Benjamin Ries
"""

# imports
import os
import copy
import inspect
import functools
import importlib
import warnings
from typing import Dict, Union, List, Callable

from pygromos.files._basics._general_gromos_file import _general_gromos_file
from pygromos.files.blocks import imd_blocks
from pygromos.files.coord.refpos import Reference_Position
from pygromos.files.coord.posres import Position_Restraints

from pygromos.files.coord import cnf
from pygromos.files.coord.cnf import Cnf
from pygromos.files.simulation_parameters.imd import Imd
from pygromos.files.topology.top import Top
from pygromos.files.topology.disres import Disres
from pygromos.files.topology.ptp import Pertubation_topology
from pygromos.files.gromos_system.ff.forcefield_system import forcefield_system

from pygromos.gromos import GromosXX, GromosPP
from pygromos.utils import bash, utils

import pickle, io

if (importlib.util.find_spec("rdkit") != None):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    has_rdkit = True
else:
    has_rdkit = False

if (importlib.util.find_spec("openff") != None):
    from openff.toolkit.topology import Molecule
    from pygromos.files.gromos_system.ff.openforcefield2gromos import openforcefield2gromos
    from pygromos.files.gromos_system.ff.serenityff.serenityff import serenityff

    has_openff = True
else:
    has_openff = False

skip = {"solute_info":cnf.solute_infos,
        "protein_info":cnf.protein_infos,
        "non_ligand_info": cnf.non_ligand_infos,
        "solvent_info": cnf.solvent_infos}


class Gromos_System():
    required_files = {"imd": Imd,
                      "top": Top,
                      "cnf": Cnf}
    optional_files = {"disres": Disres,
                      "ptp": Pertubation_topology,
                      "posres": Position_Restraints,
                      "refpos": Reference_Position}

    residue_list:Dict
    solute_info:cnf.solute_infos
    protein_info:cnf.protein_infos
    non_ligand_info:cnf.non_ligand_infos
    solvent_info:cnf.solvent_infos
    checkpoint_path:str
    _future_promise:bool #for interest if multiple jobs shall be chained.
    _future_promised_files:list

    _single_multibath:bool  = False
    _single_energy_group:bool = False

    _gromosPP_bin_dir : Union[None, str] = None
    _gromosXX_bin_dir : Union[None, str] = None
    _gromosPP : GromosPP
    _gromosXX : GromosXX

    def __init__(self, work_folder: str, system_name: str, in_smiles: str = None,
                 in_top_path: str = None, in_cnf_path: str = None, in_imd_path: str = None,
                 in_disres_path: str = None, in_ptp_path: str = None, in_posres_path:str = None, in_refpos_path:str=None,
                 in_gromosXX_bin_dir:str = None, in_gromosPP_bin_dir:str=None,
                 rdkitMol: Chem.rdchem.Mol = None, readIn=True, Forcefield:forcefield_system=forcefield_system(), 
                 auto_convert:bool=False, adapt_imd_automatically:bool=True, verbose:bool=True):
        """

        Parameters
        ----------
        work_folder
        system_name
        in_smiles
        in_top_path
        in_cnf_path
        in_imd_path
        in_disres_path
        in_ptp_path
        in_posres_path
        in_refpos_path
        in_gromosXX_bin_dir
        in_gromosPP_bin_dir
        rdkitMol
        readIn
        Forcefield
        auto_convert
        adapt_imd_automatically
        """
        self.hasData = False
        self._name = system_name
        self._work_folder = work_folder
        self.smiles = in_smiles
        self.Forcefield = Forcefield
        self.mol = Chem.Mol()
        self.checkpoint_path = None
        self.verbose = verbose

        #GromosFunctionality
        self._gromosPP_bin_dir = in_gromosPP_bin_dir
        self._gromosXX_bin_dir = in_gromosXX_bin_dir

        self._gromosPP = GromosPP(gromosPP_bin_dir=in_gromosPP_bin_dir)
        self._gromosXX = GromosXX(gromosXX_bin_dir=in_gromosXX_bin_dir)

        #add functions of gromosPP to system
        self.__bind_gromosPPFuncs()

        ## For HPC-Queueing
        self._future_promise= False
        self._future_promised_files = []

        if (in_smiles == None and rdkitMol == None) or readIn == False:
            if verbose: warnings.warn("No data provided to gromos_system\nmanual work needed")

        # import files:
        file_mapping = {"imd": in_imd_path,
                        "top": in_top_path, "ptp": in_ptp_path,
                        "cnf": in_cnf_path,
                        "disres": in_disres_path,
                        "posres": in_posres_path, 
                        "refpos": in_refpos_path,
                        }
        self.parse_attribute_files(file_mapping, readIn=readIn, verbose=verbose)

        ##System Information:
        if(not self._cnf._future_file):
            self.residue_list, self.solute_info, self.protein_info, self.non_ligand_info, self.solvent_info = self._cnf.get_system_information(
                not_ligand_residues=[])
        else:
            self.residue_list = None
            self.solute_info = None
            self.protein_info = None
            self.non_ligand_info = None

        # import molecule from smiles using rdkit
        if in_smiles:
            self.mol = Chem.MolFromSmiles(in_smiles)
            self.mol = Chem.AddHs(self.mol)
            AllChem.EmbedMolecule(self.mol)
            AllChem.UFFOptimizeMolecule(self.mol)
            self.hasData = True

        # import  molecule from RDKit
        if rdkitMol:
            self.mol = rdkitMol
            self.smiles = Chem.MolToSmiles(self.mol)
            self.hasData = True

        # automate all conversions for rdkit mols if possible
        if auto_convert:
            if self.hasData:
                self.auto_convert()
            else:
                raise Warning("auto_convert active but no data provided -> auto_convert NOT done!")

        if in_cnf_path is None and type(self.mol) == Chem.rdchem.Mol and self.mol.GetNumAtoms() >= 1:
            self.cnf = Cnf(in_value=self.mol)

        
        if(adapt_imd_automatically and not self._cnf._future_file and not  self.imd._future_file):
            self.adapt_imd()

        #misc
        self._all_files_key = list(map(lambda x: "_"+x, self.required_files.keys()))
        self._all_files_key.extend(list(map(lambda x: "_"+x, self.optional_files.keys())))
        self._all_files = copy.copy(self.required_files)
        self._all_files.update(copy.copy(self.optional_files))

    def __str__(self)->str:
        msg = "\n"
        msg += "GROMOS SYSTEM: " + self.name + "\n"
        msg += utils.spacer
        msg += "WORKDIR: " + self._work_folder + "\n"
        msg += "LAST CHECKPOINT: " + str(self.checkpoint_path) + "\n"
        msg += "\n"
        msg += "GromosXX_bin: " + str(self._gromosXX.bin) + "\n"
        msg += "GromosPP_bin: " + str(self._gromosPP.bin) + "\n"
        msg += "FILES: \n\t"+"\n\t".join([str(key)+": "+str(val) for key,val in self.all_file_paths.items()])+"\n"
        msg += "FUTURE PROMISE: "+str(self._future_promise)+"\n"
        if(hasattr(self, "solute_info")
            or hasattr(self, "protein_info")
            or hasattr(self, "non_ligand_info")
            or hasattr(self, "solvent_info")):
            msg += "SYSTEM: \n"
            if(hasattr(self, "protein_info")  and not self.protein_info is None and self.protein_info.number_of_residues > 0):
                if(hasattr(self, "protein_info")  and not self.protein_info is None):
                    #+" resIDs: "+str(self.protein_info.residues[0])+"-"+str(self.protein_info.residues[-1])
                    msg += "\tPROTEIN:\t"+str(self.protein_info.name)+"  nresidues: "+str(self.protein_info.number_of_residues)+" natoms: "+str(self.protein_info.number_of_atoms)+"\n"
                if(hasattr(self, "solute_info") and not self.solute_info is None):
                    msg += "\tLIGANDS:\t" + str(self.solute_info.names) + "  resID: " + str(self.solute_info.positions) + "  natoms: " + str(self.solute_info.number_of_atoms) + "\n"
                if (hasattr(self, "non_ligand_info")  and not self.non_ligand_info is None):
                    #+"  resID: "+str(self.non_ligand_info.positions)
                    msg += "\tNon-LIGANDS:\t"+str(self.non_ligand_info.names)+"  nmolecules: "+str(self.non_ligand_info.number)+"  natoms: "+str(self.non_ligand_info.number_of_atoms)+"\n"
                if (hasattr(self, "solvent_info")  and not self.solvent_info is None):
                    #" resIDs: "+str(self.solvent_info.positions[0])+"-"+str(self.solvent_info.positions[-1])+
                    msg += "\tSOLVENT:\t"+str(self.solvent_info.name)+"  nmolecules: "+str(self.solvent_info.number)+"  natoms: "+str(self.solvent_info.number_of_atoms)+"\n"
            else:
                if (hasattr(self, "solute_info") and not self.solute_info is None):
                    msg += "\tSolute:\t" + str(self.solute_info.names) + "  resID: " + str(
                        self.solute_info.positions) + "  natoms: " + str(
                        self.solute_info.number_of_atoms) + "\n"
                if (hasattr(self, "solvent_info") and not self.solvent_info is None):
                    # " resIDs: "+str(self.solvent_info.positions[0])+"-"+str(self.solvent_info.positions[-1])+
                    msg += "\tSOLVENT:\t" + str(self.solvent_info.name) + "  nmolecules: " + str(
                        self.solvent_info.number) + "  natoms: " + str(self.solvent_info.number_of_atoms) + "\n"

        msg +="\n\n"
        return msg

    def __repr__(self)-> str:
        return str(self)

    def __getstate__(self):
        """
        preperation for pickling:
        remove the non trivial pickling parts
        """
        attribute_dict = self.__dict__
        new_dict ={}
        for key in attribute_dict.keys():
            if (not isinstance(attribute_dict[key], Callable) and not key in skip):
                new_dict.update({key:attribute_dict[key]})
            elif(not attribute_dict[key] is None and key in skip):
                new_dict.update({key: attribute_dict[key]._asdict()})
            else:
                new_dict.update({key: None})

        return new_dict

    def __setstate__(self, state):
        self.__dict__ = state
        for key in skip:
            if(key in self.__dict__ ) and not self.__dict__[key] is None:
                setattr(self, key, skip[key](**self.__dict__[key]))

        #misc
        self._all_files_key = list(map(lambda x: "_"+x, self.required_files.keys()))
        self._all_files_key.extend(list(map(lambda x: "_"+x, self.optional_files.keys())))
        self._all_files = copy.copy(self.required_files)
        self._all_files.update(copy.copy(self.optional_files))

        self._gromosPP = GromosPP(self._gromosPP_bin_dir)
        self._gromosXX = GromosXX(self._gromosXX_bin_dir)

        self.__bind_gromosPPFuncs()

        #are promised files now present?
        self._check_promises()

    def __deepcopy__(self, memo):
        copy_obj = self.__class__(system_name="Test", work_folder=self.work_folder, readIn=False, verbose=False)
        copy_obj.verbose = self.verbose
        copy_obj.__setstate__(copy.deepcopy(self.__getstate__()))
        return copy_obj

    #def __copy__(self):



    def copy(self):
        return copy.deepcopy(self)

    """
        Properties
    """
    @property
    def work_folder(self)->str:
        return self._work_folder

    @work_folder.setter
    def work_folder(self, work_folder:str):
        self._work_folder = work_folder
        self._update_all_file_paths()

    @property
    def name(self)->str:
        return self._name

    @name.setter
    def name(self, system_name:str):
        self._name = system_name

    @property
    def all_files(self) -> Dict[str, _general_gromos_file]:
        """

        Returns
        -------

        """
        self._all_files_key = list(self.required_files.keys())
        self._all_files_key.extend(list(self.optional_files.keys()))
        self._all_files = {key: getattr(self, key) for key in self._all_files_key if
                           (hasattr(self, key) and (not getattr(self, key) is None))}
        return self._all_files

    @property
    def all_file_paths(self) -> Dict[str, str]:
        """

        Returns
        -------

        """
        self._all_files_key = list(self.required_files.keys())
        self._all_files_key.extend(list(self.optional_files.keys()))
        self._all_file_paths = {key: getattr(self, key).path for key in self._all_files_key if
                                (hasattr(self, key) and (not getattr(self, key) is None))}
        return self._all_file_paths

    @property
    def top(self)->Top:
        return self._top

    @top.setter
    def top(self, input_value:Union[str, Top]):
        if(isinstance(input_value, str)):
            if(os.path.exists(input_value)):
                self._top = Top(in_value=input_value)
            elif(self._future_promise):
                self._top = Top(in_value=input_value, _future_file=self._future_promise)
                self._future_promised_files.append("top")
                self._future_promise = True

            else:
                raise FileNotFoundError("Could not find file: "+str(input_value))
        elif(isinstance(input_value, Top)):
            self._top = input_value
        else:
            raise ValueError("Could not parse input type: " + str(type(input_value)) + " " + str(input_value))

    @property
    def cnf(self)->Cnf:
        return self._cnf

    @cnf.setter
    def cnf(self, input_value:Union[str, Cnf]):
        if(isinstance(input_value, str)):
            if(os.path.exists(input_value)):
                self._cnf = Cnf(in_value=input_value, _future_file=False)
                self.residue_list = self._cnf.get_residues()
                self.residue_list, self.solute_info, self.protein_info, self.non_ligand_info, self.solvent_info = self._cnf.get_system_information(not_ligand_residues=[])

            elif(self._future_promise):
                self._cnf = Cnf(in_value=input_value, _future_file=self._future_promise)
                self._future_promised_files.append("cnf")
                self._future_promise = True

            else:
                raise FileNotFoundError("Could not find file: "+str(input_value))
        elif(isinstance(input_value, Cnf)):
            self._cnf = input_value

        else:
            raise ValueError("Could not parse input type: " + str(type(input_value)) + " " + str(input_value))

    @property
    def imd(self)->Imd:
        return self._imd

    @imd.setter
    def imd(self, input_value:Union[str, Imd], adapt_imd_automatically:bool=True):
        if(isinstance(input_value, str)):
            if(os.path.exists(input_value) or self._future_promise):
                self._imd = Imd(in_value=input_value)
                if (adapt_imd_automatically
                    and not (self.cnf._future_file
                        and (self.residue_list is None
                             or self.solute_info is None))):
                    self.adapt_imd()
            elif(self._future_promise):
                self._imd = Imd(in_value=input_value, _future_file=self._future_promise)
                self._future_promised_files.append("imd")
                self._future_promise = True
            else:
                raise FileNotFoundError("Could not find file: "+str(input_value))
        elif(isinstance(input_value, Imd)):
            self._imd = input_value
        elif(input_value is None):
            self._imd = None
        else:
            raise ValueError("Could not parse input type: " + str(type(input_value)) + " " + str(input_value))

    @property
    def refpos(self)->Imd:
        return self._refpos

    @refpos.setter
    def refpos(self, input_value:Union[str, Reference_Position]):
        if(isinstance(input_value, str)):
            if(os.path.exists(input_value) or self._future_promise):
                self._refpos = Reference_Position(in_value=input_value)
            elif(self._future_promise):
                self._refpos = Reference_Position(in_value=input_value, _future_file=self._future_promise)
                self._future_promised_files.append("refpos")
            else:
                raise FileNotFoundError("Could not find file: "+str(input_value))
        elif(isinstance(input_value, Reference_Position)):
            self._refpos = input_value
        elif(input_value is None):
            self._refpos = None
        else:
            raise ValueError("Could not parse input type: " + str(type(input_value)) + " " + str(input_value))

    @property
    def posres(self)->Position_Restraints:
        return self._posres

    @posres.setter
    def posres(self, input_value:Union[str, Position_Restraints]):
        if(isinstance(input_value, str)):
            if(os.path.exists(input_value) or self._future_promise):
                self._posres = Position_Restraints(in_value=input_value)
            elif(self._future_promise):
                self._posres = Position_Restraints(in_value=input_value, _future_file=self._future_promise)
                self._future_promised_files.append("posres")
            else:
                raise FileNotFoundError("Could not find file: "+str(input_value))
        elif(isinstance(input_value, Position_Restraints)):
            self._posres = input_value
        elif(input_value is None):
            self._posres = None
        else:
            raise ValueError("Could not parse input type: " + str(type(input_value)) + " " + str(input_value))

    @property
    def disres(self)->Disres:
        return self._disres

    @disres.setter
    def disres(self, input_value:Union[str, Disres]):
        if(isinstance(input_value, str)):
            if(os.path.exists(input_value)):
                self._disres = Disres(in_value=input_value)
            else:
                raise FileNotFoundError("Could not find file: "+str(input_value))
        elif(isinstance(input_value, Disres)):
            self._disres = input_value
        else:
            raise ValueError("Could not parse input type: " + str(type(input_value)) + " " + str(input_value))

    @property
    def ptp(self)->Pertubation_topology:
        return self._ptp

    @ptp.setter
    def ptp(self, input_value:Union[str, Pertubation_topology]):
        if(isinstance(input_value, str)):
            if(os.path.exists(input_value)):
                self._ptp = Pertubation_topology(in_value=input_value)
            else:
                raise FileNotFoundError("Could not find file: "+str(input_value))
        elif(isinstance(input_value,  Pertubation_topology)):
            self._ptp  = input_value
        else:
            raise ValueError("Could not parse input type: " + str(type(input_value)) + " " + str(input_value))

    @property
    def gromosXX(self)->GromosXX:
        return self._gromosXX

    @property
    def gromosPP(self)->GromosPP:
        return self._gromosPP


    """
        Functions
    """

    """
    System generation
    
    """
    def get_script_generation_command(self, var_name:str=None, var_prefixes:str ="")->str:
        if(var_name is None):
            var_name=var_prefixes + self.__class__.__name__.lower()

        gen_cmd = "#Generate "+ self.__class__.__name__+ "\n"
        gen_cmd += "from " + self.__module__ + " import "+ self.__class__.__name__+" as "+ self.__class__.__name__+"_obj"+ "\n"
        gen_cmd += var_name+" = "+__class__.__name__+"_obj(work_folder=\""+self.work_folder+"\", system_name=\""+self.name+"\")\n"
        for arg, path in self.all_file_paths.items():
            gen_cmd+= str(var_name)+"."+str(arg)+" = \""+str(path)+"\"\n"
        gen_cmd+="\n"

        return gen_cmd

    def parse_attribute_files(self, file_mapping: Dict[str, str], readIn:bool=True, verbose:bool=False):
        """
            This function sets dynamically builds the output folder, the file objs of this class and checks their dependencies.

        Parameters
        ----------
        file_mapping: Dict[str, Union[str, None]]
            attribute name: input path

        Returns
        -------

        """
        # File/Folder Attribute managment
        check_file_paths = []

        ##Check if system folder is present
        if (not os.path.exists(self._work_folder)):
            bash.make_folder(self._work_folder)
        check_file_paths.append(self._work_folder)

        ## Check if file- paths are valid
        all_files = {key:val for key, val in self.required_files.items()}
        all_files.update(self.optional_files)
        [check_file_paths.append(x) for k, x in file_mapping.items() if (not x is None)]
        if (len(check_file_paths) > 0):
            bash.check_path_dependencies(check_file_paths, verbose=verbose)

        # SET Attribute-FILES
        for attribute_name, file_path in file_mapping.items():
            if(readIn and not file_path is None):
                if (verbose): print("Parsing File: ", attribute_name)
                obj = all_files[attribute_name](file_path)
            elif(attribute_name in self.required_files):
                if (verbose): print("Generate Empty: ", attribute_name)
                obj = all_files[attribute_name](None, _future_file = True)
                obj.path = file_path
                if(not file_path is None):
                    self._future_promise = True
                    self._future_promised_files.append(attribute_name)
            else:
                obj = None
            setattr(self, "_"+attribute_name, obj)

        #if read in of files - try to adapt the imd file necessaries

    def _update_all_file_paths(self):
        for file_name in self._all_files_key:
            if(hasattr(self,  file_name) and not getattr(self, file_name) is None):
                file_obj = getattr(self, file_name,)
                if(file_obj._future_file):
                    if(self.verbose or True): warnings.warn("Did not change file path as its only promised " + str(file_obj.path))
                else:
                    file_obj.path = self._work_folder + "/" + self.name + "." + getattr(self, file_name)._gromos_file_ending

    def rebase_files(self):
        if(not os.path.exists(self.work_folder) and os.path.exists(os.path.dirname(self.work_folder))):
            bash.make_folder(self.work_folder)

        self._update_all_file_paths()
        self.write_files()

    def _check_promises(self):
        #misc
        self._all_files_key = list(map(lambda x: "_"+x, self.required_files.keys()))
        self._all_files_key.extend(list(map(lambda x: "_"+x, self.optional_files.keys())))

        for promised_file_key in self._future_promised_files:
            promised_file = getattr(self, promised_file_key)
            if(os.path.exists(promised_file.path)):
                if (self.verbose):
                    print("READING FILE")
                setattr(self, "_" + promised_file_key, self._all_files[promised_file_key](promised_file.path))
                self._future_promised_files.remove(promised_file_key)
            else:
                warnings.warn("Promised file was not found: "+promised_file_key)
        if(len(self._future_promised_files) == 0):
            self._future_promise = False

    def auto_convert(self):
        #create topology
        if self.Forcefield.name == "2016H66" or self.Forcefield.name == "54A7":
            # set parameters for make_top
            out=self.work_folder+"/make_top.top"
            mtb_temp = self.Forcefield.mtb_path
            if hasattr(self.Forcefield, "mtb_orga_path"):
                mtb_temp += " " + self.Forcefield.mtb_orga_path
            ifp_temp = self.Forcefield.path
            if self.Forcefield.mol_name is None:
                name = self.rdkit2GromosName()
            else:
                name = self.Forcefield.mol_name
            # make top
            if self._gromosPP._isValid:
                self.gromosPP.make_top(out_top_path=out, in_building_block_lib_path=mtb_temp, in_parameter_lib_path=ifp_temp, in_sequence=name)
                self.top = Top(in_value=out)
            else:
                warnings.warn("could notfind a gromosPP version. Please provide a valid version for Gromos auto system generation")
             
        elif self.Forcefield.name == "smirnoff" or self.Forcefield.name == "off" or self.Forcefield.name == "openforcefield":
            if not has_openff:
                raise ImportError("Could not import smirnoff FF as openFF toolkit was missing! "
                                  "Please install the package for this feature!")
            else:
                self.top = openforcefield2gromos(Molecule.from_rdkit(self.mol), self.top,
                                                 forcefield=self.Forcefield).convert_return()
        elif self.Forcefield.name == "serenityff" or self.Forcefield.name == "sff":
            if not has_openff:
                raise ImportError("Could not import smirnoff FF as openFF toolkit was missing! "
                                  "Please install the package for this feature!")
            else:
                self.serenityff = serenityff(mol=self.mol, forcefield=self.Forcefield, top=self.top)
                self.serenityff.create_top(C12_input=self.Forcefield.C12_input, partial_charges=self.Forcefield.partial_charges)
                self.serenityff.top.make_ordered()
                self.top = self.serenityff.top

        else:
            raise ValueError("I don't know this forcefield: " + self.Forcefield.name)

    def rdkit2GromosName(self) -> str:
        raise NotImplementedError("please find your own way to get the Gromos Name for your molecule.")

    def adapt_imd(self, not_ligand_residues:List[str]=[]):
        #Get residues
        if(self.cnf._future_file and self.residue_list is None
            and self.solute_info is None
            and self.protein_info is None
            and self.non_ligand_info is None):
            raise ValueError("The residue_list, solute_info, protein_info adn non_ligand_info are required for automatic imd-adaptation.")
        else:
            if(not self.cnf._future_file):
                cres, lig, prot, nonLig, solvent = self.cnf.get_system_information(not_ligand_residues=not_ligand_residues)
                self.residue_list = cres
                self.solute_info = lig
                self.protein_info = prot
                self.non_ligand_info = nonLig
                self.solvent_info = solvent
        ##System
        self.imd.SYSTEM.NPM = 1
        self.imd.SYSTEM.NSM = len(self.residue_list ["SOLV"]) if("SOLV" in self.residue_list ) else 0

        #boundary condition:
        from pygromos.files.blocks.coord_blocks import Pbc
        ##vacuum if no box or GENBOX-pbc == Vacuum
        if (not self.cnf._future_file and (not hasattr(self.cnf, "GENBOX") or self.cnf.GENBOX.pbc == Pbc.vacuum)):
            self.imd.BOUNDCOND.NTB = Pbc.vacuum.value
        elif(not self.cnf._future_file and (not hasattr(self.cnf, "GENBOX") or self.cnf.GENBOX.pbc == Pbc.rectangular)):
            self.imd.BOUNDCOND.NTB = Pbc.rectangular.value
            self.imd.BOUNDCOND.NDFMIN = 3


        ##Energy and Force Group
        energy_groups = {}

        if(self._single_energy_group):
            energy_groups = { self.solute_info.number_of_atoms + self.protein_info.number_of_atoms + self.non_ligand_info.number_of_atoms + self.solvent_info.number_of_atoms:1}
        else:
            #solute
            if(self.solute_info.number_of_atoms>0):
                energy_groups.update({self.solute_info.positions[0]: self.solute_info.number_of_atoms})
            #protein
            if (self.protein_info.number_of_atoms > 0):
                energy_groups.update({self.protein_info.start_position: self.protein_info.number_of_atoms})
            #ligand
            if (self.non_ligand_info.number_of_atoms > 0):
                solv_add = self.non_ligand_info.number_of_atoms
            else:
                solv_add = 0
            #solvent
            if(self.solvent_info.number_of_atoms > 0):
                    max_key = max(energy_groups)+1 if(len(energy_groups)>0) else 1
                    energy_groups.update({max_key: self.solvent_info.number_of_atoms+solv_add})

            #sort all entries in list
            last_atom_count = 0
            sorted_energy_groups = {}
            for ind,key in enumerate(sorted(energy_groups)):
                value = energy_groups[key]+last_atom_count
                sorted_energy_groups.update({1+ind : value})
                last_atom_count=value

        #adapt energy groups in IMD with sorted list of energy groups created above
        self.imd.FORCE.adapt_energy_groups(energy_groups=sorted_energy_groups)

        ##Multibath:
        if (hasattr(self.imd, "MULTIBATH") and not getattr(self.imd, "MULTIBATH") is None):
            last_atoms_baths = {}
            if(self._single_multibath):
                sorted_last_atoms_baths = { self.solute_info.number_of_atoms + self.protein_info.number_of_atoms + self.non_ligand_info.number_of_atoms + self.solvent_info.number_of_atoms:1}
            else:
                if(self.solute_info.number_of_atoms>0):
                    last_atoms_baths.update({self.solute_info.positions[0]: self.solute_info.number_of_atoms})
                if (self.protein_info.number_of_atoms > 0):
                    last_atoms_baths.update({self.protein_info.start_position: self.protein_info.number_of_atoms})
                if (self.non_ligand_info.number_of_atoms > 0):
                    solv_add = self.non_ligand_info.number_of_atoms
                    #raise Exception("The imd adaptation for nonLigand res and multibath was not written yet!")
                    # TODO: Do something for non-Ligands better in future
                else:
                    solv_add = 0
                if(self.solvent_info.number_of_atoms > 0):
                    max_key = max(last_atoms_baths)+1 if(len(last_atoms_baths)>0) else 1
                    last_atoms_baths.update({max_key: self.solvent_info.number_of_atoms+solv_add})

                last_atom_count = 0
                sorted_last_atoms_baths = {}
                for ind,key in enumerate(sorted(last_atoms_baths)):
                    value = last_atoms_baths[key]+last_atom_count
                    sorted_last_atoms_baths.update({value:1+ind})
                    last_atom_count=value
            
            self.imd.MULTIBATH.adapt_multibath(last_atoms_bath=sorted_last_atoms_baths)

        ff_name = self.Forcefield.name
        if ff_name == "openforcefield" or ff_name == "smirnoff" or ff_name == "off":
            #adjust harmonic forces
            if (hasattr(self.imd, "COVALENTFORM") and not getattr(self.imd, "COVALENTFORM") is None):
                if self.verbose:
                    print("Please make sure to use harmonic forces for simulations with OpenForceField torsions")
            else:
                setattr(self.imd,"COVALENTFORM", imd_blocks.COVALENTFORM(NTBBH=1, NTBAH=1, NTBDN=0))
            #adjust amberscale for LJ forces
            if (hasattr(self.imd, "AMBER") and not getattr(self.imd, "AMBER") is None):
                if self.verbose:
                    print("Please make sure to use amber LJ forces for simulations with OpenForceField LJ parameters")
            else:
                setattr(self.imd,"AMBER", imd_blocks.AMBER(AMBER=1, AMBSCAL=1.2))

    def generate_posres(self, residues:list=[], keep_residues:bool=True, verbose:bool=False):
        self.posres = self.cnf.gen_possrespec(residues=residues, keep_residues=keep_residues, verbose=verbose)
        self.refpos = self.cnf.gen_refpos()


    """
        IO - Files:
    """
    def read_files(self, verbose:bool=False):
        for file_name in self.required_files:
            try:
                getattr(self, file_name).read_file()
            except:
                warnings.warn("did not find the required " + file_name + " found")

        for file_name in self.optional_files:
            try:
                getattr(self, file_name).read_file()
            except:
                if (verbose):
                    warnings.warn("did not find the optional " + file_name + " found")


    def write_files(self, cnf:bool=False, imd:bool=False,
                    top: bool = False, ptp:bool=False,
                    disres:bool=False, posres: bool = False, refpos: bool = False,
                    mol: bool = False,
                    all:bool=True,
                    verbose:bool=False):
        if(all):
            control = {k.replace("_", ""): all for k in self._all_files_key}
        else:
            control = {"top": top, "imd": imd, "cnf": cnf, "ptp":ptp, "disres":disres, "posres":posres, "refpos":refpos}

        for file_name in control:
            if (control[file_name] and hasattr(self, file_name) and (not getattr(self, file_name) is None)):
                file_obj = getattr(self, file_name)
                if(file_obj.path is None):
                    print("File "+str(file_name)+" is empty , can not be written!")
                elif(file_obj._future_file):
                    print("File " + str(file_name) +" is promised , can not be written: " + str(file_obj.path))
                else:
                    if(verbose): print("Write out: " + str(file_name) +"\t" + file_obj.path)
                    file_obj.write(file_obj.path)
            else:
                if(file_name in self.required_files):
                    warnings.warn("did not find the required " + file_name + " found")
                elif(verbose):
                    warnings.warn("did not find the required " + file_name + " found")

        if mol:
            print(Chem.MolToMolBlock(self.mol), file=open(self._work_folder + "/" + self.name + ".mol", 'w+'))

    def get_file_paths(self) -> Dict[str, str]:
        """
            get the paths of the files in a dict.
        Returns
        -------
        Dict[str, str]
            returns alle file paths, with attribute file name as key.
        """
        return {x: file_obj.path for x, file_obj in self.all_files.items()}

    """
        RDKIT Fun
    """

    def rdkitImport(self, inputMol: Chem.rdchem.Mol):
        self.mol = inputMol

    def rdkit2Gromos(self):
        # generate top
        self.top.add_block(blocktitle="TITLE", content=[
            "topology generated with PyGromos from RDKit molecule: " + str(Chem.MolToSmiles(self.mol)) + "    " + str(
                self.name)], verbose=True)
        self.top.add_block(blocktitle="PHYSICALCONSTANTS", content=[""], verbose=True)
        self.top.add_block(blocktitle="TOPVERSION", content=[["2.0"]])

    """
    Serialize
    """
    def save(self, path: Union[str, io.FileIO] = None, safe:bool=True) -> str:
        """
        This method stores the Class as binary obj to a given path or fileBuffer.
        """
        safe_skip = False
        if (isinstance(path, str)):
            if (os.path.exists(path) and safe):
                warnings.warn("FOUND ALREADY A FILE! SKIPPING!")
                safe_skip = True
            else:
                bufferdWriter = open(path, "wb")
        elif (isinstance(path, io.BufferedWriter)):
            bufferdWriter = path
            path = bufferdWriter.name
        else:
            raise IOError("Please give as parameter a path:str or a File Buffer. To " + str(self.__class__) + ".save")

        if(not safe_skip):
            pickle.dump(obj=self, file=bufferdWriter)
            bufferdWriter.close()
            self.checkpoint_path = path
        return path

    @classmethod
    def load(cls, path: Union[str, io.FileIO] = None) -> object:
        """
        This method stores the Class as binary obj to a given path or fileBuffer.
        """
        if (isinstance(path, str)):
            bufferedReader = open(path, "rb")
        elif (isinstance(path, io.BufferedReader)):
            bufferedReader = path
        else:
            raise IOError("Please give as parameter a path:str or a File Buffer.")

        obj = pickle.load(file=bufferedReader)

        bufferedReader.close()
        if(hasattr(obj, "cnf") and hasattr(obj.cnf, "POSITION")):
            obj.residue_list, obj.solute_info, obj.protein_info, obj.non_ligand_info, obj.solvent_info = obj._cnf.get_system_information()
        obj.checkpoint_path = path
        return obj


    """
    super privates - don't even read!
    """
    def __SystemConstructionAttributeFinder(self, func:callable) -> callable:
        """
            ** DECORATOR **

            This decorator trys to find input parameters of the function in the gromossystem and will automatically assign those to the function call!
            functional programming

        Parameters
        ----------
        func : callable

        Returns
        -------
        callable
            returns the wrapped function

        """

        @functools.wraps(func)
        def findGromosSystemAttributes(*args, **kwargs):
            #print(func.__name__, args, kwargs)

            # collect input parameters present in system/ replace them with
            tmp_files = []
            for k in inspect.signature(func).parameters:
                attr_key = k.replace("in_", "").replace("_path", "")
                #print(attr_key)
                if ("in" in k and "path" in k and attr_key in dir(self)):
                    grom_obj = getattr(self, attr_key)
                    if (grom_obj.path is None):
                        tmp_file_path = self.work_folder + "/tmp_file." + grom_obj._gromos_file_ending
                        grom_obj.write(tmp_file_path)
                        kwargs.update({k: tmp_file_path})
                        tmp_files.append(tmp_file_path)
                    else:
                        grom_obj.write(grom_obj.path)  # make sure filestatus is good :)
                        kwargs.update({k: grom_obj.path})

            # execute function
            r = func(*args, **kwargs)

            # remove tmp_files
            [bash.remove_file(p) for p in tmp_files]

            return r
        return findGromosSystemAttributes

    def __SystemConstructionUpdater(self, func:callable) -> callable:
        """
            ** DECORATOR **
            This decorator trys to find output parameters of the function in the gromossystem and will automatically update the state of those attributes!
            functional programming

        Parameters
        ----------
        func: callable
            the function to be wrapped

        Returns
        -------
        func

        """

        @functools.wraps(func)
        def updateGromosSystem(*args, **kwargs):
            #rint(func.__name__, args, kwargs)

            # collect out_paths
            update_dict = {}
            for k in inspect.signature(func).parameters:
                if ("out" in k and "path" in k):
                    attr_key = k.replace("out_", "").replace("_path", "")
                    kwargs.update({k: self.work_folder + "/tmp_file." + attr_key})
                    update_dict.update({k: attr_key})

            # execute function
            r = func(*args, **kwargs)

            # update attribute states and remove tmp files.
            for k in update_dict:
                setattr(self, update_dict[k], kwargs[k])
                getattr(self, update_dict[k]).path = None
                bash.remove_file(kwargs[k])

            return r
        return updateGromosSystem

    def __bind_gromosPPFuncs(self):
        if(not self._gromosPP is None):
            func = [k for k in dir(self._gromosPP) if (not k.startswith("_") and k != "bin")]
            v = {f: self.__SystemConstructionUpdater(self.__SystemConstructionAttributeFinder(getattr(self._gromosPP, f)))
                 for f in func}
            self.__dict__.update(v)