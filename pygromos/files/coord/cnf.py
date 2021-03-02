import copy
from typing import Dict, List, Tuple, NamedTuple

from pygromos.files._basics import parser
from pygromos.files._basics._general_gromos_file import _general_gromos_file

from pygromos.files.blocks import coord_blocks as blocks
from pygromos.utils.utils import _cartesian_distance


class Cnf(_general_gromos_file):
    """
    This class is a representation of the gromos .cnf coordinate files. It
    allows reading, analysis and modifying of the coordinate files.

    is a child of general_gromos_file
    """

    #general
    _orig_file_path:str
    residues:Dict[str, Dict[str, int]]

    #Standard Gromos blocks
    TITLE: blocks.TITLE   #required
    POSITION: blocks.POSITION  #required

    TIMESTEP: blocks.TIMESTEP = None
    LATTICESHIFTS: blocks.LATTICESHIFTS = None
    VELOCITY: blocks.VELOCITY = None
    GENBOX: blocks.GENBOX = None
    atom_ref_pos_block: blocks.REFPOSITION = None

    #private
    _block_order: List[str] = ["TITLE", "TIMESTEP", "POSITION", "LATTICESHIFTS", "VELOCITY", "REFPOSITION"]
    _required_blocks = ["TITLE", "POSITION"]

    def __init__(self, input:(str or dict or None or __class__), verbose:bool=False, clean_resiNumbers_by_Name=False):
        super().__init__(input=input)

        if(hasattr(self, "POSITION")):
            if clean_resiNumbers_by_Name: self.clean_posiResNums()  #carefull! if two resis same name after an another than, here is  a problem.
            self.residues = self.get_residues(verbose=verbose)

    def read_file(self)->Dict[str, any]:
        """

        Parameters
        ----------
        path :

        Returns
        -------

        """
        return parser.read_cnf(self._orig_file_path)

    def add_residue_positions(self, coords:object):
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

    def clean_residue_list_for_imd(self, not_ligand_residues:List[str], ligand_resn_prefix:(str or List[str])=[])-> (Dict[str, Dict[int,int]], NamedTuple, NamedTuple, NamedTuple):
        """clean_residue_list_for_imd
               This function utilizes a dictionary containing all residues and atom numbers (e.g. cnf.get_residues()) and modifies them such, that the result can be used to set up a standard REEDS gromos_simulation

               Pattern: @Decorator

        Parameters
        ----------
        not_ligand_residues :    List[str]
            here all molecules, that are not considered as ligand or protein.
        ligand_resn_prefix :     List[str]
             here all molecules, that are considered as ligand are listed.

        Returns
        -------
        (Dict[str, Dict[int,int]], NamedTuple, NamedTuple, NamedTuple)
            cleaned_residue dict, ligands, protein, non_ligands
        """

        if(type(self.residues) == type(None)):
            raise Exception("Residue dict field was not initialized! please do so.")

        from pygromos.files import imd
        return imd.Imd.clean_residue_list_for_imd(self.residues, not_ligand_residues, ligand_resn_prefix)

    def clean_posiResNums(self)->None:
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
        tmpOldID=pos[0].resID

        for p in pos:
            #print(p)
            #print(tmpN,tmpID)
            if p.resName == tmpN and p.resID == tmpOldID:               #same residue as before
                p.resID = tmpID
            elif(p.resName == tmpN and p.resID != tmpOldID):            #same resiname but diff ID (double? - this is a problem!)
                tmpOldID = p.resID
                tmpID += 1
                p.resID = tmpID
            else:                                                       #next name and residue id
                tmpID +=1
                tmpN = p.resName
                tmpOldID = p.resID
                p.resID = tmpID

        self.POSITION.content = pos

    def supress_atomPosition_singulrarities(self)->None:
        """supress_atomPosition_singulrarities
        This function adds a very small deviation to the position of an atom, dependent on the atom number.
        This might be needed to avoid singularities in gromosXX.

        Returns
        -------
        None
        """

        if("POSITION" in dir(self)):
            for ind, atom in enumerate(self.POSITION.content):
                atom.xp = atom.xp + 10 ** (-7) * ind
                atom.yp = atom.yp - 10 ** (-7) * ind
                atom.zp = atom.zp - 10 ** (-7) * ind

    def delete_block(self, blockName:str):
        self._blocksset_names.remove(blockName)
        setattr(self, blockName, None)

    def rename_residue(self, new_resName:str, resID:int=False, resName:str=False, verbose=False) -> int:
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
                result_sub = [] #build up new data list

                if block in dir(self):
                    subblock = getattr(self, block)
                    for i, x in enumerate(subblock.content):    #go through block
                        if ((int(x.resID) == resID or not resID) and (x.resName == resName or not resName)):    #found correct residue to be changed
                            x.resName = new_resName
                        result_sub.append(x)
                    setattr(subblock, "content", result_sub)
                    setattr(self, block, subblock)
                elif verbose:
                    print("Skip block: " + block + ", as was not found in file.\n")
        else:
            raise IOError("No residue number or resName given.")
        return 0

    def delete_residue(self, resID:int=False, resName:str=False, verbose=False) -> int:
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

        if((not resID or not resName ) and verbose):
            print("WARNING giving only resID or resName can be ambigous!")
        if resID or resName:
            #deletable with resID
            for block in check_blocks[:2]:
                result_sub = []

                if block in dir(self):
                    subblock = getattr(self, block)
                    if(isinstance(subblock, type(None))):
                        check_blocks.remove(block)
                        continue

                    offset_atomID = 0       #recount atomids
                    res_off = 0             #recount residues
                    tmp_del_resID = -10
                    for i,x in enumerate(subblock.content):
                        if((int(x.resID) == resID or not resID) and (x.resName == resName or not resName)):
                            atom_nums.append(x.atomID)  #append to list of being deleted
                            offset_atomID+=1
                            res_off = 1
                            #catch multiple to be deleted residues after each other in a residue offset : for recounting
                            if(tmp_del_resID != x.resID-1 and tmp_del_resID != x.resID):
                                res_off +=1
                            else:
                                res_off = 1
                            tmp_del_resID = x.resID
                            continue
                        else:
                            x.atomID -= offset_atomID
                            if( x.resID - res_off <= 0):
                                res_off = -1*res_off-1
                                x.resID -= res_off
                            else:
                                x.resID -= res_off
                            result_sub.append(x)
                    setattr(subblock, "content", result_sub)
                    setattr(self, block, subblock)
                elif verbose:
                    print("Skip block: "+block+", as was not found in file.\n")

            # delete by atomnum
            for block in check_blocks[2:]:
                if(hasattr(self, block)):
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

    def count_residue_atoms(self, resID:int=False, resName:str =False, verbose=False)->(int or dict):

            if(resName and resID):
                print("\t"+resName+"\t"+str(resID)+"\t"+str(self.residues[resName][resID]))
                return self.residues[resName][resID]
            if(resName and not resID):
                print(resName+"\n".join(["\t"+str(x)+"\t"+str(self.residues[resName][x]) for x in self.residues[resName]]))

                return self.residues[resName]
            if(not resName and resID):
                result = []
                for x in self.residues:
                    if(resID in self.residues[x]):
                        result.update({x: self.residues[x]})
                    else:
                        continue
                if(verbose):
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

    def get_residues(self, verbose:bool=False)->Dict[str, Dict[str, int]]:
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
            #find Residues and Count them.
            for i, x in enumerate(subblock):

                #TREAT SOLVENT
                if (x.resName == "SOLV" or x.resName == "SOL"):
                    if not (x.resName in counts):
                        counts.update({x.resName: 1})
                    else:
                        counts[x.resName] += 1

                #TREAT Unknown Molecules
                elif (x.resName in counts and x.resID in counts[x.resName]):
                    counts[x.resName][x.resID] += 1
                elif (x.resName in counts and not x.resID in counts[x.resName]):
                    counts[x.resName].update({x.resID: 1})
                else:
                    counts.update({x.resName: {x.resID: 1}})

            if verbose: print("COUNTS", counts)

            #make string
            if verbose:
                out_string = "#\tresName\tresID:counts\n"
                #generate string
                for resn in counts:
                    out_string += "\t"+resn
                    out_string += "\t \t" + str(counts[resn]) + "\n"
                print(out_string)
            return counts
        else:
            raise ValueError("NO POSITION block in cnf-Object: " + self.file_path)

    def get_atomP(self, atomID:int=None, atomType:str=None, resID:int=None, resName:str=False)->list:

        if("POSITION" in dir(self)):
            atoms = []
            for x in self.POSITION.content:
                if ((x.resID == resID or not resID) and (x.resName == resName or not resName) and (x.atomID == atomID or not atomID) and (x.atomType == atomType or not atomType)):
                    atoms.append(x)

            if(len(atoms)>0):
                return atoms
            else:
                raise ValueError("COULD not find Atom with ID: "+str(atomID)+" in CNF object.")
        else:
            raise ValueError("NO POSITION block in cnf-Object: " + self.file_path)

    def get_atoms_distance(self, atomI:(int or tuple)=None, atomJ:int=None, atoms:(list or dict)=None)->float or list:

        if("POSITION" in dir(self)):
            #one
            if (atomI != None and atomJ != None) and atoms == None:
                ai = self.get_atomP(atomI)[0]
                aj = self.get_atomP(atomJ)[0]

                distance =_cartesian_distance(x1=ai.xp ,x2=aj.xp,y1=ai.yp, y2=aj.yp, z1=ai.zp,z2=aj.zp )
                return distance

            #one atom coord set to all atoms of cnf
            elif atomI != None and atomJ == None and atoms == None and type(atomI) is tuple:
                if(len(atomI) != 3):
                    raise ValueError("atomI-Coordinate tuple have more or less than 3 dims in get_atoms_distance")
                distances = {}
                for x in self._blocks["POSITION"].content:
                    aj = x
                    distance =_cartesian_distance(x1=atomI[0] ,x2=aj.xp,y1=atomI[1], y2=aj.yp, z1=atomI[2],z2=aj.zp )
                    distances.update({x.atomID: distance})
                return distances

            #one coord set + list of atoms
            elif atomI != None and atomJ == None and atoms != None and type(atomI) is tuple:
                if(len(atomI) != 3):
                    raise ValueError("atomI-Coordinate tuple have more or less than 3 dims in get_atoms_distance")
                distances = {}
                for x in atoms:
                    aj = self.get_atomP(x)[0]
                    distance =_cartesian_distance(x1=atomI[0] ,x2=aj.xp,y1=atomI[1], y2=aj.yp, z1=atomI[2],z2=aj.zp )
                    distances.update({aj.atomID: distance})
                return distances

            elif atomI != None and atomJ == None and atoms != None and type(atomI) is int:
                ai = self.get_atomP(atomI)[0]
                distances = {}
                for x in atoms:
                    aj = self.get_atomP(x)[0]
                    distance =_cartesian_distance(x1=ai.xp ,x2=aj.xp,y1=ai.yp, y2=aj.yp, z1=ai.zp,z2=aj.zp )
                    distances.update({x: distance})
                return distances

            elif atomI == None and atomJ == None and atoms != None and atoms is list:
                distances_all = []
                for x in atoms:
                    ai = self.get_atomP(x)[0]
                    distances = []
                    for y in atoms[atoms.index(x):]:
                        aj = self.get_atomP(y)[0]
                        distance = _cartesian_distance(x1=ai.xp ,x2=aj.xp,y1=ai.yp, y2=aj.yp, z1=ai.zp,z2=aj.zp )
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
                        distance =_cartesian_distance(x1=ai.xp ,x2=aj.xp,y1=ai.yp, y2=aj.yp, z1=ai.zp,z2=aj.zp )
                        distances.update({y: distance})
                    distances_all.append({x: distances})
                return distances_all
            else:
                raise ValueError("Combination of given parameters for get_atoms_distance in cnf class is unknown!\n Please give: int int or int list or list")
        else:
            raise ValueError("NO POSITION block in cnf-Object: " + self.file_path)

    def write(self, out_path:str)->0:
        #write out
        out_file = open(out_path, "w")
        out_file.write(self.__str__())
        out_file.close()
        return out_path

    def generate_position_restraints(self, out_path_prefix:str, residues: dict or list, verbose:bool=False) -> Tuple[str,str]:
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
        possrespec = self.write_possrespec(out_path=out_path_prefix+".por", residues=residues,  verbose=verbose)
        refpos = self.write_refpos(path= out_path_prefix+".rpf")

        return possrespec, refpos

    def write_possrespec(self, out_path:str, residues: dict or list, verbose:bool=False)->str:
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
        Cnf
            posrespec-file-obj
        """
        delete_res={}   #dict for the atoms, that should be deleted.
        cnf = copy.deepcopy(self)   #deepcopied cnf for output

        #Get all atoms not included to posres restrain and store them in delete_res
        try:
            if(type(residues) == dict): #if adict was given - todo: not tested
                for res in cnf.residues:
                    if(res not in residues):
                        delete_res.update({res:[-1]})
                    else:
                        ids = [resi for resi in list(cnf.residues[res].keys()) if(resi not in residues[res])]
                        delete_res.update({res: ids})

            if(type(residues) == list): #if a list of resIDS was given
                if(all([type(x) == str for x in residues])):    #for resn
                    for res in cnf.residues:
                        if(res not in residues):
                            delete_res.update({res:[-1]})

                elif(all([type(x) == int for x in residues])):  #for resids
                    for res in cnf.residues:
                        res_ids = cnf.residues[res]
                        if(type(res_ids) == dict):
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
            raise NotImplementedError("Posres _file input arg combination : Not implemented! "+"\n".join(err.args))

        #Remove all not to constrain atoms:
        if verbose: print("delete residues: ", delete_res)
        for resn,resi_list in delete_res.items():
            if(type(resi_list) == dict):
                for resi in resi_list:
                    cnf.delete_residue(resName=resn, resID=resi)
            else:
                cnf.delete_residue(resName=resn)

        if verbose: print("remaining: ", cnf.get_residues())

        #del _blocks of cnf, that are not needed is posres
        not_del_blocks = ["TITLE", "POSITION"]
        for del_block in cnf._blocksset_names:
            if (not del_block in not_del_blocks):
                cnf.delete_block(del_block)

        #write out
        out_file = open(out_path, "w")
        out_file.write(cnf.__str__().replace("POSITION", "POSRESSPEC"))
        out_file.close()
        return out_path, cnf

    def write_refpos(self, path: str) -> str:

        out_file = open(path, "w")

        cnf = copy.deepcopy(self)

        #del _blocks, not needed
        not_del_blocks = ["TITLE", "POSITION", "GENBOX"]
        for del_block in cnf._blocksset_names:
            if (not del_block in not_del_blocks):
                cnf.delete_block(del_block)

        out_file.write(cnf.__str__().replace("POSITION", "REFPOSITION"))
        out_file.close()
        return path, cnf


