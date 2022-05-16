"""
FUNCTIONLIB:            gromos++ input file functions
Description:
    in this lib, gromosXX input file mainpulating functions are gathered
Author: Kay Schaller & Benjamin Schroeder
"""
import numpy as np
from numbers import Number
from typing import List, Dict, NamedTuple, Iterable
from collections import namedtuple

from pygromos.files._basics import _general_gromos_file, parser
from pygromos.files.blocks import imd_blocks as blocks
from pygromos.utils import amino_acids as aa


ligand_infos = namedtuple("ligands_info", ["names", "number", "positions", "number_of_atoms"])
protein_infos = namedtuple("protein_info", ["name", "residues", "number_of_residues", "position", "start_position", "end_position", "number_of_atoms"])
non_ligand_infos = namedtuple("ligands_info", ["names", "number", "positions",  "number_of_atoms"])

def nice_s_vals(svals:Iterable, base10=False) ->list:
    """

    Parameters
    ----------
    svals :
    base10 :

    Returns
    -------

    """

    nicer_labels = []
    if(base10):
        for val in svals:
            if(float(np.log10(val)).is_integer() or val == min(svals)):
                nicer_labels.append(round(val, str(val).count("0")+3))
            else:
                nicer_labels.append("")
    else:
        for val in svals:
            # nicer_labels.append(round(float(val), str(val).count("0")+2)) # old line
            nicer_labels.append(round(float(val), 5))
    return nicer_labels

class Imd(_general_gromos_file._general_gromos_file):

    _orig_file_path:str
    path:str
    _required_blocks = ["TITLE", "SYSTEM", "STEP", "BOUNDCOND", "FORCE", "CONSTRAINT", "PAIRLIST", "NONBONDED"] # Not jet implemented

    #POSSIBLE GROMOS BLOCKS
    TITLE: blocks.TITLE

    STEP: blocks.STEP
    WRITETRAJ: blocks.WRITETRAJ
    PRINTOUT: blocks.PRINTOUT

    SYSTEM: blocks.SYSTEM
    FORCE: blocks.FORCE

    ENERGYMIN: blocks.ENERGYMIN = None
    EDS:blocks.EDS=None

    CONSTRAINT: blocks.CONSTRAINT = None
    BOUNDCOND: blocks.BOUNDCOND = None
    COMTRANSROT: blocks.COMTRANSROT = None

    PRESSURESCALE: blocks.PRESSURESCALE = None
    MULTIBATH: blocks.MULTIBATH = None

    PAIRLIST: blocks.PAIRLIST = None
    NONBONDED: blocks.NONBONDED = None

    POSITIONRES: blocks.POSITIONRES = None
    DISTANCERES: blocks.DISTANCERES = None
    INITIALISE: blocks.INITIALISE = None

    REPLICA_EDS: blocks.REPLICA_EDS = None
    #OLD_REPLICA_EDS: blocks.OLD_REPLICA_EDS = None
    #NEW_REPLICA_EDS: blocks.NEW_REPLICA_EDS = None

    def __init__(self, input:str):
        super().__init__(input=input)

    def __str__(self):
        text = ""
        text += self.__getattribute__("TITLE").block_to_string()
        for block in sorted(self.get_block_names()):
            if(block == "TITLE" or isinstance(block, type(None))):
                continue
            text += str(self.__getattribute__(block))
        return text

    def read_file(self):
        data = parser.read_imd(self._orig_file_path)
        for key, sub_content in data.items():
            try:
                self.add_block(blocktitle=key, content=sub_content)
            except Exception as err:
                try:
                    if (key == "REPLICA_EDS"):  # TODO: remove as soon as new block is established! or change to old >)
                        self.add_block(blocktitle="NEW_REPLICA_EDS", content=sub_content)
                        self.REPLICA_EDS = self.NEW_REPLICA_EDS
                        self.NEW_REPLICA_EDS = ""
                    else:
                        raise IOError("Could not read in imd " + key + " block!\n values: \n\t" + str(sub_content) + "\n\n" + "\n\t".join(err.args))
                except Exception as err:
                    try:
                        if(key == "REPLICA_EDS"): #TODO: remove as soon as new block is established! or change to old >)
                            self.add_block(blocktitle="OLD_REPLICA_EDS", content=sub_content)
                            self.REPLICA_EDS = self.OLD_REPLICA_EDS
                        else:
                            raise IOError("Could not read in imd "+key+" block!\n values: \n\t"+str(sub_content)+"\n\n"+"\n\t".join(err.args))
                    except Exception as err:
                        raise IOError("could not read in reeds_imd "+key+" block!\n values: \n\t"+str(sub_content)+"\n\n"+"\n\t".join(err.args))
        return {}

    def edit_EDS(self, NUMSTATES:int, S:float, EIR:list, EDS:int=1, ALPHLJ:float=0.0, ALPHC:float=0.0, FUNCTIONAL:int=1):
        if(hasattr(self, "eds_block")):
            self.EDS.ALPHLJ = ALPHLJ
            self.EDS.ALPHC = ALPHC
            self.EDS.FUNCTIONAL = FUNCTIONAL
            self.EDS.NUMSTATES = NUMSTATES
            self.EDS.S = S
            self.EDS.EIR = EIR

        else:
            if(type(EIR) == float or type(EIR) == str or type(EIR) == int):
                EIR = [float(EIR) for x in range(NUMSTATES)]

            gromos_name = "EDS"
            setattr(self, gromos_name, blocks.EDS(NUMSTATES, S, EIR, EDS, ALPHLJ, ALPHC, FUNCTIONAL))


    def randomize_seed(self):
        self.INITIALISE.IG = np.random.randint(low=0, high=999999)

    def edit_REEDS(self, REEDS:(bool or int)=None, NUMSTATES:int=None, SVALS: (Number, List[Number])=None, EIR:(Number or Iterable[Number])=None,
                   NRETRIAL:int=None, NREQUIL:int=None, CONT:(bool, int)=None, EDS_STAT_OUT:(bool, int)=None,
                   RETS:List[float]= None, RET:int=None, NATOM:int=None) :   #TODO: old params - to be REMOVED!

        # specific relations are rescued here
        reeds_block = self.REPLICA_EDS

        if(isinstance(reeds_block, blocks.NEW_REPLICA_EDS)):
            if(isinstance(REEDS, bool)):
                reeds_block.REEDS = REEDS

            if(isinstance(NUMSTATES, (Number, str))):
                reeds_block.NUMSTATES = NUMSTATES

            if(isinstance(SVALS, Iterable)): #edit SVALS
                SVALS = nice_s_vals(SVALS)

                reeds_block.RES = list(map(str, SVALS))
                reeds_block.NRES = len(SVALS)   #adjust number of Svals
                
                if(int(reeds_block.REEDS) == 1): #1D REEDS
                    reeds_block.NEOFF = len(SVALS)
                elif(int(reeds_block.REEDS) > 1): #2D REEDS
                    reeds_block.NEOFF = len(SVALS)*(reeds_block.NUMSTATES+1)
                else: #No REEDS
                    reeds_block.NEOFF = 0
                    
                if(not isinstance(EIR, (Number, Iterable))): #expand energy offsets to new s value ammount
                    if(any([len(row)!= len(SVALS) for row in reeds_block.EIR])):
                        newEIR = []
                        for EIR_row in reeds_block.EIR:
                            newEIRrow = [EIR_row[0] for r in range(len(SVALS))]
                            newEIR.append(newEIRrow)
                        reeds_block.EIR = newEIR

            if(isinstance(EIR, (Number , Iterable))):
                EIR_matrix = []
                #print(EIR)
                # single number
                if isinstance(EIR, Number):  # depends on SVALS and NRES
                    EIR_vector = [str(EIR) for x in range(reeds_block.NUMSTATES)]

                    for z in EIR_vector:
                        EIR_matrix.append([z for i in range(int(reeds_block.NRES))])

                 # vector
                elif (isinstance(EIR, Iterable) and len(EIR) == int(reeds_block.NUMSTATES) and all(
                        [isinstance(x, (Number, str)) for x in EIR])):
                    for z in EIR:
                        EIR_matrix.append([float(z) for i in range(int(reeds_block.NRES))])
                # matrix
                elif(isinstance(EIR, Iterable) and all([isinstance(x, Iterable) and all([isinstance(y, (Number, str)) for y in x]) for x in EIR])):
                        if (self.REPLICA_EDS.NRES == len(EIR)):
                            EIR = np.array(EIR).T
                        EIR_matrix = list(map(lambda x: list(map(float, x)), EIR))
                else:
                    raise Exception(
                        "not enough EIR-vals for making REEDS Block. Got " + str(len(EIR)) + " for " + str(
                            reeds_block.NRES) + " SVals\n Number of states: "+str(reeds_block.NUMSTATES))
                reeds_block.EIR = EIR_matrix
                
            if (isinstance(NRETRIAL, (str, int))):
                reeds_block.NRETRIAL = NRETRIAL
            if (isinstance(NREQUIL, (str, int))):
                reeds_block.NREQUIL = NREQUIL
            if (isinstance(CONT, (str, int, bool))):
                reeds_block.CONT = CONT
            if (isinstance(EDS_STAT_OUT, (str, int, bool))):
                reeds_block.EDS_STAT_OUT = EDS_STAT_OUT

     
        else:
            raise Exception("neither old nor new REPLICA_EDS block recognized\n")

    @staticmethod
    def clean_residue_list_for_imd(residues:Dict[str, Dict[int,int]], not_ligand_residues:List[str]=[], ligand_resn_prefix:(str or List[str])=None)-> (Dict[str, Dict[int,int]], NamedTuple, NamedTuple, NamedTuple):
        """clean_residue_list_for_imd
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
        
        #Build Up new list
        ##Criterium, when a ligand is considered
        if(isinstance(ligand_resn_prefix, str)):
            ligand_resn_prefix = [ligand_resn_prefix]

        #print(ligand_resn_prefix)
        ligand_residue = lambda res: ((res != "SOLV" and (res not in aa.three_letter_aa_lib and res != "prot")) and not res in not_ligand_residues ) or (type(ligand_resn_prefix) != type(None) and res in ligand_resn_prefix)

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
            clean_residues = {name: {min(residues[name].keys()): sum(list(residues[name].values()))} for name in ligand_names}  # ligands as resi

        ###Update
        ligand_names = [res for res in clean_residues if ligand_residue(res)]
        number_of_ligands_atoms = sum([sum(list(clean_residues[res].values())) for res in ligand_names])
        number_of_ligands = len(ligand_names)
        ligand_positions = [min(clean_residues[res]) for res in ligand_names]
        ligands = ligand_infos(names=ligand_names, number=number_of_ligands, positions=ligand_positions, number_of_atoms=number_of_ligands_atoms)

        ## get protein parameters if present
        protein_residues = {res: val for res, val in residues.items() if (res in aa.three_letter_aa_lib)}
        if (len(protein_residues) > 0):
            protein_name= "protein"
            number_of_protein_residues = sum([len(list(residues[res].keys())) for res in protein_residues])
            number_of_protein_atoms = sum([sum(list(residues[res].values())) for res in protein_residues])
            protein_start_position = min([min(val) for res, val in protein_residues.items()])
            protein_end_position = max([max(val) for res, val in protein_residues.items()])
            clean_residues.update({protein_name: {len(clean_residues) + 1: number_of_protein_atoms}})  # protein Atoms as one residue
            protein = protein_infos(name=protein_name, residues=protein_residues, number_of_residues= number_of_protein_residues,
                                    position=protein_start_position, start_position=protein_start_position, end_position=protein_end_position,
                                    number_of_atoms=number_of_protein_atoms)

        else:
            protein = protein_infos(name="", residues=0, number_of_residues= 0, position=0, number_of_atoms=0, start_position=-1, end_position=-1)
            
        ##get non_Ligand_residues, e.g.: Cofactors
        excluded_resi_names = [res for res in residues if (res in not_ligand_residues)]
        if (len(excluded_resi_names) > 0):
            number_of_non_ligands_atoms = sum([sum(list(residues[res].values())) for res in excluded_resi_names])
            number_of_non_ligands = len(ligand_names)
            non_ligand_positions = min([min(clean_residues[res]) for res in ligand_names])
            clean_residues.update({name: {min(residues[name].keys()): sum(list(residues[name].values()))} for name in excluded_resi_names})  # protein Atoms as one residue
            non_ligands= non_ligand_infos(names=excluded_resi_names, number=number_of_non_ligands, positions=non_ligand_positions, number_of_atoms=number_of_non_ligands_atoms)
        else:
            non_ligands= non_ligand_infos(names=[], number=0, positions=0, number_of_atoms=0)
        ##get Solvent
        if ("SOLV" in residues):
            clean_residues.update({"SOLV": residues["SOLV"]})  # solvent

        return clean_residues, ligands, protein, non_ligands
