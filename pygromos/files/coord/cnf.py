import copy
import os

import __main__

from collections import namedtuple
from typing import Dict, List, Tuple, Union
from typing import TypeVar
import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation
from rdkit import Chem
from rdkit.Chem import AllChem

from pygromos.files._basics import parser
from pygromos.files._basics._general_gromos_file import _general_gromos_file
from pygromos.files.blocks import coord_blocks as blocks
from pygromos.files.blocks._general_blocks import TITLE
from pygromos.utils import amino_acids as aa
from pygromos.utils.utils import _cartesian_distance
from pygromos.files.trajectory.trc import Trc
from pygromos.files.blocks.coord_blocks import GENBOX

Reference_Position = TypeVar("refpos.Reference_Position")
Position_Restraints = TypeVar("posres.Position_Restraints")
# Constructs describing the system
solute_infos = namedtuple("solute_info", ["names", "number", "positions", "number_of_atoms"])
protein_infos = namedtuple("protein_info",
                           ["name", "residues", "number_of_residues", "position", "start_position", "end_position",
                            "number_of_atoms"])
non_ligand_infos = namedtuple("ligands_info", ["names", "number", "positions", "number_of_atoms"])
solvent_infos = namedtuple("solvent_info", ["name", "number", "atoms_per_residue", "positions", "number_of_atoms"])





# make infos pickle-able

class Cnf(_general_gromos_file):
    """
    This class is a representation of the gromos .cnf coordinate files. It
    allows reading, analysis and modifying of the coordinate files.

    is a child of general_gromos_file
    """

    # general
    _orig_file_path: str
    gromos_file_ending: str = "cnf"
    residues: Dict[str, Dict[str, int]]

    # Standard Gromos blocks
    TITLE: blocks.TITLE  # required
    POSITION: blocks.POSITION  # required

    TIMESTEP: blocks.TIMESTEP
    LATTICESHIFTS: blocks.LATTICESHIFTS
    VELOCITY: blocks.VELOCITY
    GENBOX: blocks.GENBOX
    PERTDATA: blocks.PERTDATA
    atom_ref_pos_block: blocks.REFPOSITION

    # private
    _block_order: List[str] = ["TITLE", "TIMESTEP", "POSITION", "LATTICESHIFTS", "VELOCITY", "REFPOSITION"]
    _required_blocks: List[str] = ["TITLE", "POSITION"]
    _main_block:str = "POSITION"

    def __init__(self, in_value: (str or dict or None or __class__),
                 clean_resiNumbers_by_Name=False,
                 verbose: bool = False, _future_file: bool = False):
        #import for rdkit molecules
        if type(in_value) == Chem.rdchem.Mol:
            super().__init__(in_value=None, _future_file=_future_file)
            self.createRDKITconf(mol=in_value)
        #general import
        else:
            super().__init__(in_value=in_value, _future_file=_future_file)
            if (hasattr(self, "POSITION")):
                if clean_resiNumbers_by_Name: self.clean_posiResNums()  # carefull! if two resis same name after an another than, here is  a problem.
                self.residues = self.get_residues(verbose=verbose)

    """
        Manage coordinates
    """
    def delete_block(self, blockName: str):
        self._blocksset_names.remove(blockName)
        setattr(self, blockName, None)

    def add_empty_box(self):
        self.add_block(block=GENBOX())

    def read_file(self) -> Dict[str, any]:
        """

        Parameters
        ----------
        path :

        Returns
        -------

        """
        return parser.read_cnf(self._orig_file_path)
    
    def write(self, out_path: str) -> str:
        # write out
        out_file = open(out_path, "w")
        out_file.write(self.__str__())
        out_file.close()
        self.path = os.path.abspath(out_path)
        return out_path

    """
        manipulate/analysis of coordinates
    """
    def add_residue_positions(self, coords: object):
        """This function adds all residues of an coords file to @DEVELOP

        This is a very crude functio at the moment! It only takes positions
        of a residue and merges them! if there are residues with the same name,
        this might lead to problems, as clean_posiResNumByName function is not
        sensitive for that! todo: make more robust! bschroed

        Parameters
        ----------
        coords :function_libs.gromos.files.coord.Cnf obj
            object - CNF object

        Returns
        -------

        """

        positions = coords.POSITION.content
        self.POSITION.content.extend(positions)
        self.clean_posiResNums()
        self.get_residues(verbose=True)

    def get_system_information(self, not_ligand_residues: List[str] = [],
                                     ligand_resn_prefix: (str or List[str]) = None,
                                     solvent_name: str = "SOLV") -> \
            (Dict[str, Dict[int, int]], namedtuple, namedtuple, namedtuple, namedtuple):
        """get_system_information
        This function utilizes a dictionary containing all residues and atom numbers (e.g. cnf.get_residues()) and modifies them such, that the result can be used to set up a standard REEDS gromos_simulation

        Parameters
        ----------
        residues : Dict[str, Dict[int,int]]
             input a cnf.residues:dict that shall be cleaned and return a reduced form for parameter file.
        not_ligand_residues :List[str]
            here all molecules, that are not considered as ligand or protein.
        ligand_resn_prefix :List[str]
            here all molecules, that are considered as ligand are listed.

        Returns
        -------
        Dict[str, Dict[int,int]]
            cleaned_residue dict
        NamedTuple
            ligands
        NamedTuple
            protein
        NamedTuple
            non_ligands
        """

        residues = self.get_residues()
        # Build Up new list
        ##Criterium, when a ligand is considered
        if (isinstance(ligand_resn_prefix, str)):
            ligand_resn_prefix = [ligand_resn_prefix]

        # print(ligand_resn_prefix)
        ligand_residue = lambda res: ((res != "SOLV" and (
                res not in aa.three_letter_aa_lib and res != "prot")) and not res in not_ligand_residues) or (
                                             type(ligand_resn_prefix) != type(None) and res in ligand_resn_prefix)
        # print(ligand_residue)

        ##get ligand parameters
        ###Avoid multi ligands with same resi name!
        ligand_names = [res for res in residues if ligand_residue(res)]
        if (any([len(residues[name]) > 1 for name in ligand_names])):
            multi_res_ligands = [name for name in ligand_names if (len(residues[name]) > 1)]
            clean_residues = {}
            for name in ligand_names:
                if (name in multi_res_ligands):
                    for resi in residues[name].keys():
                        new_name = (name[:-1] + str(resi))
                        clean_residues.update({new_name: {resi: sum(list([residues[name][resi]]))}})
                else:
                    clean_residues.update({name: {min(residues[name].keys()): sum(list(residues[name].values()))}})
        else:
            clean_residues = {name: {min(residues[name].keys()): sum(list(residues[name].values()))} for name in
                              ligand_names}  # ligands as resi

        ###Update
        ligand_names = [res for res in clean_residues if ligand_residue(res)]
        number_of_ligands_atoms = sum([sum(list(clean_residues[res].values())) for res in ligand_names])
        number_of_ligands = len(ligand_names)
        ligand_positions = [min(clean_residues[res]) for res in ligand_names]
        ligands = solute_infos(names=ligand_names, number=number_of_ligands, positions=ligand_positions,
                               number_of_atoms=number_of_ligands_atoms)

        ## get protein parameters if present
        protein_residues = {res: val for res, val in residues.items() if (res in aa.three_letter_aa_lib)}
        if (len(protein_residues) > 0):
            protein_name = "protein"
            number_of_protein_residues = sum([len(list(residues[res].keys())) for res in protein_residues])
            number_of_protein_atoms = sum([sum(list(residues[res].values())) for res in protein_residues])
            protein_start_position = min([min(val) for res, val in protein_residues.items()])
            protein_end_position = max([max(val) for res, val in protein_residues.items()])
            clean_residues.update(
                {protein_name: {len(clean_residues) + 1: number_of_protein_atoms}})  # protein Atoms as one residue
            protein = protein_infos(name=protein_name, residues=protein_residues,
                                    number_of_residues=number_of_protein_residues,
                                    position=protein_start_position, start_position=protein_start_position,
                                    end_position=protein_end_position,
                                    number_of_atoms=number_of_protein_atoms)

        else:
            protein = protein_infos(name="", residues=0, number_of_residues=0, position=0, number_of_atoms=0,
                                    start_position=-1, end_position=-1)

        ##get non_Ligand_residues, e.g.: Cofactors
        excluded_resi_names = [res for res in residues if (res in not_ligand_residues)]
        if (len(excluded_resi_names) > 0):
            number_of_non_ligands_atoms = sum([sum(list(residues[res].values())) for res in excluded_resi_names])
            number_of_non_ligands = len(ligand_names)
            non_ligand_positions = 0  # min([min(clean_residues[res]) for res in ligand_names]) #TODO: Fix
            clean_residues.update({name: {min(residues[name].keys()): sum(list(residues[name].values()))} for name in
                                   excluded_resi_names})  # protein Atoms as one residue
            non_ligands = non_ligand_infos(names=excluded_resi_names, number=number_of_non_ligands,
                                           positions=non_ligand_positions, number_of_atoms=number_of_non_ligands_atoms)
        else:
            non_ligands = non_ligand_infos(names=[], number=0, positions=0, number_of_atoms=0)
        ##get Solvent
        if (solvent_name in residues):
            clean_residues.update({solvent_name: residues[solvent_name]})  # solvent
            positions = list(residues[solvent_name].keys())
            all_atoms = sum(list(residues[solvent_name].values()))
            natoms_per_resi = int(residues[solvent_name][positions[0]])
            solvent = solvent_infos(name=solvent_name, number=len(positions), atoms_per_residue=natoms_per_resi, positions=positions, number_of_atoms=all_atoms)
        else:
            solvent = solvent_infos(name=[], number=0, atoms_per_residue=0, positions=[], number_of_atoms=0)

        return clean_residues, ligands, protein, non_ligands, solvent

    def rename_residue(self, new_resName: str, resID: int = False, resName: str = False, verbose=False) -> int:
        """ rename_residue

        this function is renaming residues from a cnf file with taking following _blocks into account:
            "POSITION", "VELOCITY"

        it additionally recounts all atomIds and residueIDs afterwards.
        you can provide a residue ID or a residue Name or both (than only exact match will be deleted).


        Parameters
        ----------
        new_resName :   str
            new Name of the residue
        resID : int, optional
            Id of the residue to be renamed
        resName :   str, optional
            Name of the residue to be renamed
        verbose :   bool, optional
            Text... lots of it!
        Returns
        -------
        int
            0 if succesfull
        """

        check_blocks = ["POSITION", "VELOCITY"]

        atom_nums = []
        if ((not resID or not resName) and verbose):
            print("WARNING giving only resID or resName can be ambigous!")
        if resID or resName:
            # deletable with resID
            for block in check_blocks[:2]:
                result_sub = []  # build up new data list

                if block in dir(self):
                    subblock = getattr(self, block)
                    for i, x in enumerate(subblock.content):  # go through block
                        if ((int(x.resID) == resID or not resID) and (
                                x.resName == resName or not resName)):  # found correct residue to be changed
                            x.resName = new_resName
                        result_sub.append(x)
                    setattr(subblock, "content", result_sub)
                    setattr(self, block, subblock)
                elif verbose:
                    print("Skip block: " + block + ", as was not found in file.\n")
        else:
            raise IOError("No residue number or resName given.")
        return 0

    def clean_posiResNums(self) -> None:
        """clean_posiResNums
            This function recount the Residue number with respect to residue name and residue number.

        Warnings
        --------
        only in "Position_BLOCK!@development

        Returns
        -------
        None
        """
        position_copy = self.POSITION
        pos = position_copy.content
        tmpN = ""
        tmpID = 0
        tmpOldID = pos[0].resID

        for p in pos:
            # print(p)
            # print(tmpN,tmpID)
            if p.resName == tmpN and p.resID == tmpOldID:  # same residue as before
                p.resID = tmpID
            elif (
                    p.resName == tmpN and p.resID != tmpOldID):  # same resiname but diff ID (double? - this is a problem!)
                tmpOldID = p.resID
                tmpID += 1
                p.resID = tmpID
            else:  # next name and residue id
                tmpID += 1
                tmpN = p.resName
                tmpOldID = p.resID
                p.resID = tmpID

        self.POSITION.content = pos

    def delete_residue(self, resID: int = False, resName: str = False, verbose=False) -> int:
        """delete_residue

        this function is deleting residues from a cnf file with taking following _blocks into account:
            "POSITION", "VELOCITY", "LATTICESHIFTS", "REFPOSITION"

        it additionally recounts all atomIds and residueIDs afterwards.
        you can provide a residue ID or a residue Name or both (than only exact match will be deleted).



        Parameters
        ----------
        resID : int
            Id of the residue to be deleted
        resName :   str
            Name of the residue to be deleted
        verbose :   bool
             Text... lots of it!

        Returns
        -------
        int
             0 if succesfull
        """

        check_blocks = ["POSITION", "VELOCITY", "LATTICESHIFTS", "REFPOSITION"]
        atom_nums = []

        if ((not resID or not resName) and verbose):
            print("WARNING giving only resID or resName can be ambigous!")
        if resID or resName:
            # deletable with resID
            for block in check_blocks[:2]:
                result_sub = []

                if block in dir(self):
                    subblock = getattr(self, block)
                    if (isinstance(subblock, type(None))):
                        check_blocks.remove(block)
                        continue

                    offset_atomID = 0  # recount atomids
                    res_off = 0  # recount residues
                    tmp_del_resID = -10
                    for i, x in enumerate(subblock.content):
                        if ((int(x.resID) == resID or not resID) and (x.resName == resName or not resName)):
                            atom_nums.append(x.atomID)  # append to list of being deleted
                            offset_atomID += 1
                            res_off = 1
                            # catch multiple to be deleted residues after each other in a residue offset : for recounting
                            if (tmp_del_resID != x.resID - 1 and tmp_del_resID != x.resID):
                                res_off += 1
                            else:
                                res_off = 1
                            tmp_del_resID = x.resID
                            continue
                        else:
                            x.atomID -= offset_atomID
                            if (x.resID - res_off <= 0):
                                res_off = -1 * res_off - 1
                                x.resID -= res_off
                            else:
                                x.resID -= res_off
                            result_sub.append(x)
                    setattr(subblock, "content", result_sub)
                    setattr(self, block, subblock)
                elif verbose:
                    print("Skip block: " + block + ", as was not found in file.\n")

            # delete by atomnum
            for block in check_blocks[2:]:
                if (hasattr(self, block)):
                    for x in self.__getattribute__(check_blocks[2]).content:
                        if (x.atomID in atom_nums):
                            list(self.__getattribute__(check_blocks[2]).content).remove(x)
                elif verbose:
                    print("Skip block: " + block + ", as was not found in file.\n")

                    for x in self.__getattribute__(check_blocks[2]).content:
                        if (x.atomID in atom_nums):
                            list(self.__getattribute__(check_blocks[2]).content).remove(x)


        else:
            raise IOError("No residue number or resName given.")

        return 0

    def delete_atom(self, resID: int = False, resName: str = False, atomID: int = False, atomType: str = False):

        check_blocks = ["POSITION", "VELOCITY", "LATTICESHIFTS", "REFPOSITION"]
        atom_nums = []

        if (not resID or not resName):
            print("WARNING giving only resID or resName can be ambigous!")
        if resID or resName:
            # deletable with resID
            for block in check_blocks[:2]:
                result_sub = []

                if block in dir(self):
                    subblock = getattr(self, block)
                    offset_atomID = 0
                    res_off = 0
                    for i, x in enumerate(subblock.content):
                        if ((x.resID == resID or not resID) and (x.resName == resName or not resName) and (
                                x.atomID == atomID or not atomID) and (x.atomType == atomType or not atomType)):
                            atom_nums.append(x.atomID)
                            offset_atomID += 1
                            res_off = 1
                            continue
                        else:
                            x.atomID -= offset_atomID
                            if (x.resID - res_off == 0 and x.resName != "Solv"):
                                res_off = 0
                            else:
                                x.resID -= res_off
                            result_sub.append(x)

                    setattr(subblock, "content", result_sub)
                    setattr(self, block, subblock)

            # deletable with atomnum
            for block in check_blocks[2:]:
                for x in self._blocks[check_blocks[2]].content:
                    if (x.atomID in atom_nums):
                        list(self._blocks[block].content).remove(x)

        else:
            raise IOError("No residue number or resName given.")

        return 0

    def count_residue_atoms(self, resID: int = False, resName: str = False, verbose=False) -> (int or dict):

        if (resName and resID):
            print("\t" + resName + "\t" + str(resID) + "\t" + str(self.residues[resName][resID]))
            return self.residues[resName][resID]
        if (resName and not resID):
            print(resName + "\n".join(
                ["\t" + str(x) + "\t" + str(self.residues[resName][x]) for x in self.residues[resName]]))

            return self.residues[resName]
        if (not resName and resID):
            result = []
            for x in self.residues:
                if (resID in self.residues[x]):
                    result.update({x: self.residues[x]})
                else:
                    continue
            if (verbose):
                # make string
                if verbose:
                    out_string = "#\tresName\tresID\tcounts\n"
                    # generate string
                    for resn in result:
                        out_string += "\t" + resn
                        for id in result[resn]:
                            out_string += "\t" + str(id) + "\t" + str(result[resn][id]) + "\n"
                    print(out_string)

            return result
        else:
            raise ValueError("Please verify at least resName or resId for residue atom count ")

    def get_residues(self, verbose: bool = False) -> Dict[str, Dict[str, int]]:
        """get_residues
            This function is getting all residues of the used cnf file.
            it gives back,

        Parameters
        ----------
        verbose :   bool, optional
            texty?

        Returns
        -------
        Dict[str, Dict[str, Any]]
            returns dict containing all residues, numbers and atoms.
        """

        if ("POSITION" in dir(self)):
            counts = {}
            subblock = self.POSITION.content
            # find Residues and Count them.
            for i, x in enumerate(subblock):

                ##TREAT SOLVENT
                # if (x.resName == "SOLV" or x.resName == "SOL"):
                #    if not (x.resName in counts):
                #        counts.update({x.resName: 1})
                #    else:
                #        counts[x.resName] += 1

                # TREAT Unknown Molecules
                if (x.resName in counts and x.resID in counts[x.resName]):
                    counts[x.resName][x.resID] += 1
                elif (x.resName in counts and not x.resID in counts[x.resName]):
                    counts[x.resName].update({x.resID: 1})
                else:
                    counts.update({x.resName: {x.resID: 1}})

            if verbose: print("COUNTS", counts)

            # make string
            if verbose:
                out_string = "#\tresName\tresID:counts\n"
                # generate string
                for resn in counts:
                    out_string += "\t" + resn
                    out_string += "\t \t" + str(counts[resn]) + "\n"
                print(out_string)
            return counts
        else:
            raise ValueError("NO POSITION block in cnf-Object: " + self.path)

    def get_atomP(self, atomID: int = None, atomType: str = None, resID: int = None, resName: str = False) -> list:

        if ("POSITION" in dir(self)):
            atoms = []
            for x in self.POSITION.content:
                if ((x.resID == resID or not resID) and (x.resName == resName or not resName) and (
                        x.atomID == atomID or not atomID) and (x.atomType == atomType or not atomType)):
                    atoms.append(x)

            if (len(atoms) > 0):
                return atoms
            else:
                raise ValueError("COULD not find Atom with ID: " + str(atomID) + " in CNF object.")
        else:
            raise ValueError("NO POSITION block in cnf-Object: " + self.path)

    def get_atoms_distance(self, atomI: (int or tuple) = None, atomJ: int = None,
                           atoms: (list or dict) = None) -> float or list:

        if ("POSITION" in dir(self)):
            # one
            if (atomI != None and atomJ != None) and atoms == None:
                ai = self.get_atomP(atomI)[0]
                aj = self.get_atomP(atomJ)[0]

                distance = _cartesian_distance(x1=ai.xp, x2=aj.xp, y1=ai.yp, y2=aj.yp, z1=ai.zp, z2=aj.zp)
                return distance

            # one atom coord set to all atoms of cnf
            elif atomI != None and atomJ == None and atoms == None and type(atomI) is tuple:
                if (len(atomI) != 3):
                    raise ValueError("atomI-Coordinate tuple have more or less than 3 dims in get_atoms_distance")
                distances = {}
                for x in self._blocks["POSITION"].content:
                    aj = x
                    distance = _cartesian_distance(x1=atomI[0], x2=aj.xp, y1=atomI[1], y2=aj.yp, z1=atomI[2], z2=aj.zp)
                    distances.update({x.atomID: distance})
                return distances

            # one coord set + list of atoms
            elif atomI != None and atomJ == None and atoms != None and type(atomI) is tuple:
                if (len(atomI) != 3):
                    raise ValueError("atomI-Coordinate tuple have more or less than 3 dims in get_atoms_distance")
                distances = {}
                for x in atoms:
                    aj = self.get_atomP(x)[0]
                    distance = _cartesian_distance(x1=atomI[0], x2=aj.xp, y1=atomI[1], y2=aj.yp, z1=atomI[2], z2=aj.zp)
                    distances.update({aj.atomID: distance})
                return distances

            elif atomI != None and atomJ == None and atoms != None and type(atomI) is int:
                ai = self.get_atomP(atomI)[0]
                distances = {}
                for x in atoms:
                    aj = self.get_atomP(x)[0]
                    distance = _cartesian_distance(x1=ai.xp, x2=aj.xp, y1=ai.yp, y2=aj.yp, z1=ai.zp, z2=aj.zp)
                    distances.update({x: distance})
                return distances

            elif atomI == None and atomJ == None and atoms != None and atoms is list:
                distances_all = []
                for x in atoms:
                    ai = self.get_atomP(x)[0]
                    distances = []
                    for y in atoms[atoms.index(x):]:
                        aj = self.get_atomP(y)[0]
                        distance = _cartesian_distance(x1=ai.xp, x2=aj.xp, y1=ai.yp, y2=aj.yp, z1=ai.zp, z2=aj.zp)
                        distances.append(distance)
                    distances_all.append(distances)
                return distances_all

            elif atomI != None and atomJ != None and atoms != None and atoms is dict:
                distances_all = {}
                for x in atoms:
                    ai = self.get_atomP(x)
                    distances = {}
                    for y in atoms[atoms.index(x):]:
                        aj = self.get_atomP(y)
                        distance = _cartesian_distance(x1=ai.xp, x2=aj.xp, y1=ai.yp, y2=aj.yp, z1=ai.zp, z2=aj.zp)
                        distances.update({y: distance})
                    distances_all.append({x: distances})
                return distances_all
            else:
                raise ValueError(
                    "Combination of given parameters for get_atoms_distance in cnf class is unknown!\n Please give: int int or int list or list")
        else:
            raise ValueError("NO POSITION block in cnf-Object: " + self.path)

    def center_of_geometry(self, selectedAtoms:list=None) -> list:
        """calculates the center of geometry for asingle molecule or the selected Atoms

        Returns
        -------
        list
            cog
        """
        if ("POSITION" in dir(self)):
            cogx = 0.0
            cogy = 0.0
            cogz = 0.0
            if selectedAtoms is None:
                iterator = self.POSITION.content
            else:
                iterator = []
                for i in selectedAtoms:
                    iterator.append(self.POSITION.content[i])
            for i in iterator:
                cogx += i.xp
                cogy += i.yp
                cogz += i.zp
            n = len(self.POSITION.content)
            return [cogx/n, cogy/n, cogz/n]
        else:
            raise ValueError("NO POSITION block in cnf-Object: " + self.path)

    def rotate(self, rotationCenter:np.array=None, selectedAtoms=None, alpha=0, beta=0, gamma=0):
        #define rotation center
        if rotationCenter is None:
            rotationCenter = np.array(self.center_of_geometry())
        
        #select atoms to rotate
        if selectedAtoms is None:
            iterator = range(len(self.POSITION.content))

        #define Rotation
        rotation = Rotation.from_euler('XYZ', [alpha, beta, gamma], degrees=True).as_matrix()

        #apply rotation
        for i in iterator:
            #get atom
            atom = self.POSITION.content[i]
            vector = np.array([atom.xp, atom.yp, atom.zp])
            #do rotation around rotation center
            new_vec = np.dot((vector-rotationCenter), rotation)+rotationCenter
            #rewrite atom
            atom.xp, atom.yp, atom.zp = new_vec
            self.POSITION.content[i] = atom

    def supress_atomPosition_singulrarities(self) -> None:
        """supress_atomPosition_singulrarities
        This function adds a very small deviation to the position of an atom, dependent on the atom number.
        This might be needed to avoid singularities in gromosXX.

        Returns
        -------
        None
        """

        if ("POSITION" in dir(self)):
            for ind, atom in enumerate(self.POSITION.content):
                atom.xp = atom.xp + 10 ** (-7) * ind
                atom.yp = atom.yp - 10 ** (-7) * ind
                atom.zp = atom.zp - 10 ** (-7) * ind



    """
        generate further file
    """
    def generate_position_restraints(self, out_path_prefix: str, residues: dict or list, verbose: bool = False) -> \
    Tuple[str, str]:
        """

        Parameters
        ----------
        out_path_prefix
        residues

        parameters - Not implemented yet

        verbose

        Returns
        -------

        """
        possrespec = self.write_possrespec(out_path=out_path_prefix + ".por", residues=residues, verbose=verbose)
        refpos = self.write_refpos(out_path=out_path_prefix + ".rpf")

        return possrespec, refpos

    def gen_possrespec(self, residues: Union[dict or list], keep_residues:bool=True, verbose: bool = False) -> Position_Restraints:
        """
                This function writes out a gromos file, containing a atom list. that is to be position restrained!
                Raises not Implemented error, if a input variant of the residues


        Parameters
        ----------
        residues : dict, list
            residues to be restrained (dict containing resname:[res ids] or list of resnames or residue IDs)
        keep_residues : bool, optional
            should the passed residues be kept or deleted?
        verbose : bool, optional
            loud and noisy?

        Returns
        -------
        str
            out_path
        Position_Restraints
            posrespec-file-obj
        """

        from pygromos.files.coord import posres

        delete_res = {}  # dict for the atoms, that should be deleted.
        cnf = copy.deepcopy(self)  # deepcopied cnf for output

        # Get all atoms not included to posres restrain and store them in delete_res
        try:
            if (type(residues) == dict):  # if adict was given - todo: not tested
                for res in cnf.residues:
                    if (res not in residues):
                        delete_res.update({res: [-1]})
                    else:
                        ids = [resi for resi in list(cnf.residues[res].keys()) if (resi not in residues[res])]
                        delete_res.update({res: ids})

            if (type(residues) == list):  # if a list of resIDS was given
                if (all([type(x) == str for x in residues])):  # for resn
                    for res in cnf.residues:
                        if (keep_residues and res not in residues):
                            delete_res.update({res: [-1]})
                        elif(not keep_residues and res in residues):
                            delete_res.update({res: [-1]})

                elif (all([type(x) == int for x in residues])):  # for resids
                    for res in cnf.residues:
                        res_ids = cnf.residues[res]
                        if (type(res_ids) == dict):
                            for resi in res_ids:
                                if (not resi in residues):
                                    if (res in delete_res):
                                        delete_res[res].append(res_ids)
                                    else:
                                        delete_res.update({res: [res_ids]})
                        else:
                            delete_res.update({res: [res_ids]})
            else:
                raise Exception("I will be catched and translated in the except >)")
        except Exception as err:
            raise NotImplementedError("Posres _file input arg combination : Not implemented! " + "\n".join(err.args))

        # Remove all not to constrain atoms:
        if verbose: print("delete residues: ", delete_res)
        for resn, resi_list in delete_res.items():
            if (type(resi_list) == dict):
                for resi in resi_list:
                    cnf.delete_residue(resName=resn, resID=resi)
            else:
                cnf.delete_residue(resName=resn)

        if verbose: print("remaining: ", cnf.get_residues())

        return posres.Position_Restraints(cnf)

    def write_possrespec(self, out_path: str, residues: dict or list, verbose: bool = False) -> str:
        """write_possrespec
                This function writes out a gromos file, containing a atom list. that is to be position restrained!
                Raises not Implemented error, if a input variant of the residues

        Parameters
        ----------
        out_path :  str
            path to the outputfile
        residues : dict, list
            residues to be restrained (dict containing resname:[res ids] or list of resnames or residue IDs)
        verbose :   bool, optional
            loud and noisy?

        Raises
        ------
         NotImplementedError

        Returns
        -------
        str
            out_path
        Position_Restraints
            posrespec-file-obj
        """
        posres_class = self.gen_possrespec(residues=residues, verbose=verbose)
        posres_class.write(out_path)
        return out_path

    def gen_refpos(self) -> Reference_Position:
        from pygromos.files.coord import refpos
        return refpos.Reference_Position(self)

    def write_refpos(self, out_path: str) -> str:
        refpos_class = self.gen_refpos()
        refpos_class.write(out_path)
        return out_path

    """
        convert file: 
    """
    def createRDKITconf(self, mol:Chem.rdchem.Mol):
        """creates a PyGromosTools CNF type from a rdkit molecule. If a conformation exists the first one will be used.

        Parameters
        ----------
        mol : Chem.rdchem.Mol
            Molecule, possibly with a conformation
        """
        inchi=Chem.MolToInchi(mol).split("/")
        if len(inchi) >= 2:
            name=inchi[1]
        else:
            name="XXX"
        self.__setattr__("TITLE", TITLE("\t"+ name +" created from RDKit"))

        #check if conformations exist else create a new one
        if mol.GetNumConformers() < 1:
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol)
        conf = mol.GetConformer(0)

        #fill a list with atomP types from RDKit data
        atomList=[]
        for i in range(mol.GetNumAtoms()):
            x = conf.GetAtomPosition(i).x
            y = conf.GetAtomPosition(i).y
            z = conf.GetAtomPosition(i).z
            atomType = mol.GetAtomWithIdx(i).GetSymbol()
            atomList.append(blocks.atomP(resID=1, resName=name, atomType=atomType, atomID=i+1, xp=x, yp=y, zp=z))

        # set POSITION attribute
        self.__setattr__("POSITION", blocks.POSITION(atomList))
        self.__setattr__("GENBOX", blocks.GENBOX(pbc=1, length=[4,4,4], angles=[90,90,90]))

    def get_pdb(self)->str:
        """
            translate cnf to pdb.

        Returns
        -------
        str
            pdb str.

        """

        # 2) CONSTUCT PDB BLOCKS
        # ref: https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
        pdb_format = "ATOM  {:>5d} {:<4}{:1}{:<4} {:1}{:>3d}{:1}   {:>7.3f} {:>7.3f} {:>7.3f} {:>5}{:>6}{:<3}{:>2} {:>2d}"

        dummy_occupancy = dummy_bfactor = dummy_charge = 0.0
        dummy_alt_location = dummy_chain = dummy_insertion_code = dummy_segment = ""

        frame_positions = []
        for ind, atom in enumerate(self.POSITION):
            frame_positions.append(
                pdb_format.format(atom.atomID, atom.atomType, dummy_alt_location, atom.resName, dummy_chain, int(
                    atom.resID), dummy_insertion_code, atom.xp * 10, atom.yp * 10, atom.zp * 10, dummy_occupancy,
                                  dummy_bfactor,
                                  dummy_segment, atom.atomType, int(dummy_charge)))  # times *10 because pdb is in A

        pdb_str = "TITLE " + str(self.TITLE).replace("END", "") + "\n"
        pdb_str += "\n".join(frame_positions) + '\n'
        pdb_str += '\nEND'

        #future:
        ##add bond block with top optionally

        return pdb_str

    def get_xyz(self)->str:
        """
            translate cnf to xyz

        Returns
        -------
        str
            in xyz format
        """
        xyz_str = str(len(self.POSITION)) + "\n"
        xyz_str += "# "+str(self.TITLE.content[0])
        xyz_str += "# exported wit PyGromosTools\n"
        xyz_format = "{:<3}\t{:> 3.9f}   {:> 3.9f}   {:> 3.9f}\n"

        for position in self.POSITION:
            xyz_line = xyz_format.format(position.atomType[0], position.xp * 10, position.yp * 10, position.zp * 10)
            xyz_str += xyz_line

        return xyz_str

    def write_xyz(self, out_path: str)->str:
        """
            This function converts the atom POS db of the traj into a xyz structure.

        Parameters
        ----------
        out_path : str
            path, were the file should be written to.

        Returns
        -------
        str
            outpath of the file

        """
        return self._write_to_file(out_path=out_path, content_str=self.get_xyz())

    def write_pdb(self, out_path: str)->str:
        """
            This function converts the atom POS db of the traj into a pdb traj.

        Parameters
        ----------
        out_path : str
            path, were the file should be written to.

        Returns
        -------
        str
            outpath of the file

        """
        return self._write_to_file(out_path=out_path, content_str=self.get_pdb())

    def cnf2trc(self) -> Trc:
        """This function converts a cnf to a trajectory with a single frame

        Returns
        -------
        Trc
            Trc with informations from Cnf(self)
        """
        
        #create empty Trc
        trc = Trc(input_value=None)

        #set normal blocks
        trc.TITLE = self.TITLE

        #create dict for pd DataFrame
        dict = {}

        if hasattr(self,"TIMESTEP"):
            dict["TIMESTEP_step"] = self.TIMESTEP.step
            dict["TIMESTEP_time"] = self.TIMESTEP.t
        else:
            dict["TIMESTEP_step"] = 0
            dict["TIMESTEP_time"] = 0
        
        if hasattr(self,"POSITION"):
            for pos in self.POSITION:
                dict["POS_"+str(pos.atomID)] = np.array([pos.xp,pos.yp,pos.zp]).astype(np.float)
        else:
            raise ValueError("Cnf has no POSITION data")

        #TODO: Maybe add othe blocks to translate?

        #build dataframe and return
        trc.database = pd.DataFrame([dict])
        return trc

    """
        visualize
    """

    def visualize(self):
        from pygromos.visualization.coordinates_visualization import show_cnf
        return  show_cnf(self)