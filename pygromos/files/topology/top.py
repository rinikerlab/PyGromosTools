"""
File:            gromos++ topo file functions
Warnings: this CLASS IS NOT IMPLEMENTED!
TODO:REWORK
Description:
    in this lib, gromos topo file mainpulating functions are gathered
Author: Marc Lehner, Benjamin Ries
"""

#imports
from copy import deepcopy
from typing import TypeVar
import warnings
import math

from pygromos.utils import bash as bash
from pygromos.files._basics import _general_gromos_file, parser
from pygromos.files.blocks import topology_blocks as blocks


TopType = TypeVar("Top")


#functions
def make_topolog(input_arg, build, param, seq, solve="H2O"):
    #define python command
    command="make_top "+input_arg+" "+param+" "+seq+" "+solve+" \n"

    #execute command
    try:
        bash.execute(command=command)
    except Exception as err:
        bash.increment_error_level(err_prefix="Could not make_topology due to: ", old_err=err)

    return command

def combine_topologies():
    raise Exception('not implemented yet!')
    return 0

def check_top():
    raise Exception('not implemented yet!')
    return 0

#file Classes
class Top(_general_gromos_file._general_gromos_file):
    gromos_file_ending:str = "top"

    def __init__(self, in_value:(str or dict or None or TopType), _future_file:bool=False):
        if type(in_value) is str:
            super().__init__(in_value=in_value, _future_file=_future_file)
        elif(in_value == None):
            self.path = ""
            self.block_names = {}
            super().__init__(in_value=None)

        elif(type(in_value) is __class__):
            raise Exception('not implemented yet!')
        else:
            raise Exception('not implemented yet!')

    def __add__(self, top:TopType)->TopType:
        return self._add_top(top=top)

    def _add_top(self, top:TopType, solvFrom1:bool=True, verbose:bool=False)->TopType:
        """
        combines two topologies. Parameters are taken from the initial topology. 
        But missing parameters from the second topology will be added.
        Can be used like com_top from Gromos++ 

        Parameters
        ----------
        top : TopType
            second topology to add to the first topology
        solvFrom1 : bool, optional
            should the solvent be taken from the first topology? (else second), by default True
        verbose : bool, optional
            extra print statements, by default True
        sanityCheck : bool, optional
            feature and compatibility check, by default True

        Returns
        -------
        TopType
            returns a topology made by combing two topologies
        """
        # create the return top
        retTop = deepcopy(self)
        # add solv
        if not solvFrom1:
            if verbose: print("taking solvent from second topology")
            retTop.SOLVENTATOM = top.SOLVENTATOM
            retTop.SOLVENTCONSTR = top.SOLVENTCONSTR

        #calculate the shift of atom types of the second topology and add new atomtypes
        atomTypeShift = {}
        for idx, atomT in enumerate(top.ATOMTYPENAME.content[1:]): #new atomtypes to find names for
            foundAtomType = False
            for mainIdx, mainAtomT in enumerate(retTop.ATOMTYPENAME.content[1:]): #AtomTypes in self to match against
                if atomT == mainAtomT:
                    foundAtomType = True
                    atomTypeShift.update({idx+1:mainIdx+1})
                    break
            if not foundAtomType:
                retTop.ATOMTYPENAME.content[0] = str(int(retTop.ATOMTYPENAME.content[0]) + 1)
                retTop.ATOMTYPENAME.content.append(atomT)

        if verbose: print("atomTypeShift: " + str(atomTypeShift))

        #add RESNAME
        for resname in top.RESNAME.content[1:]:
            retTop.add_new_resname(resname[0])

        #add SOLUTEATOM
        atnmShift = retTop.SOLUTEATOM.content[-1].ATNM #Number of atoms found in main top. Shift secondary top atoms accordingly
        if verbose: print("atom number shift: " + str(atnmShift))
        mresShift = retTop.SOLUTEATOM.content[-1].MRES #Number of molecules found in main top.
        if verbose: print("molecule number shift: " + str(mresShift))

        for atom in top.SOLUTEATOM.content:
            retTop.add_new_soluteatom(ATNM = atnmShift + atom.ATNM,
                                    MRES = mresShift + atom.MRES,
                                    PANM = atom.PANM,
                                    IAC = atomTypeShift[atom.IAC],
                                    MASS = atom.MASS,
                                    CG = atom.CG,
                                    CGC = atom.CGC,
                                    INE = [str(int(x)+atnmShift) for x in atom.INEvalues],
                                    INE14 = [str(int(x)+atnmShift) for x in atom.INE14values])

        # add bonds and bonds with H
        for bond in top.BOND.content:
            bondType = top.BONDSTRETCHTYPE.content[bond.ICB - 1]
            retTop.add_new_bond(k=bondType.CHB, 
                                b0=bondType.B0, 
                                atomI=bond.IB + atnmShift, 
                                atomJ=bond.JB + atnmShift)
        for bond in top.BONDH.content:
            bondType = top.BONDSTRETCHTYPE.content[bond.ICB - 1]
            retTop.add_new_bond(k=bondType.CHB, 
                                b0=bondType.B0, 
                                atomI=bond.IB + atnmShift, 
                                atomJ=bond.JB + atnmShift,
                                includesH=True)

        # add angles and angles with H
        for angle in top.BONDANGLE.content:
            angleType = top.BONDANGLEBENDTYPE.content[angle.ICT - 1]
            retTop.add_new_angle(k=angleType.CB, 
                                kh=angleType.CHB, 
                                b0=angleType.B0, 
                                atomI=angle.IT + atnmShift, 
                                atomJ=angle.JT + atnmShift,
                                atomK=angle.KT + atnmShift)
        for angle in top.BONDANGLEH.content:
            angleType = top.BONDANGLEBENDTYPE.content[angle.ICT - 1]
            retTop.add_new_angle(k=angleType.CB, 
                                kh=angleType.CHB, 
                                b0=angleType.B0, 
                                atomI=angle.IT + atnmShift, 
                                atomJ=angle.JT + atnmShift,
                                atomK=angle.KT + atnmShift, includesH=True)

        # add diheadrals and diheadrals with H
        for dihdrl in top.DIHEDRAL.content:
            dihdrlType = top.TORSDIHEDRALTYPE.content[dihdrl.ICP - 1]
            retTop.add_new_torsiondihedral(CP=dihdrlType.CP, 
                                        PD=dihdrlType.PD, 
                                        NP=dihdrlType.NP, 
                                        atomI=dihdrl.IP + atnmShift, 
                                        atomJ=dihdrl.JP + atnmShift, 
                                        atomK=dihdrl.KP + atnmShift, 
                                        atomL=dihdrl.LP + atnmShift)
        for dihdrl in top.DIHEDRALH.content:
            dihdrlType = top.TORSDIHEDRALTYPE.content[dihdrl.ICPH - 1]
            retTop.add_new_torsiondihedral(CP=dihdrlType.CP, 
                                        PD=dihdrlType.PD, 
                                        NP=dihdrlType.NP, 
                                        atomI=dihdrl.IPH + atnmShift, 
                                        atomJ=dihdrl.JPH + atnmShift, 
                                        atomK=dihdrl.KPH + atnmShift, 
                                        atomL=dihdrl.LPH + atnmShift,
                                        includesH=True)

        # add impdihedrals with and without H
        for dihdrl in top.IMPDIHEDRAL.content:
            dihdrlType = top.IMPDIHEDRALTYPE.content[dihdrl.ICQ - 1]
            retTop.add_new_impdihedral(CQ=dihdrlType.CQ, 
                                        Q0=dihdrlType.Q0,
                                        atomI=dihdrl.IQ + atnmShift, 
                                        atomJ=dihdrl.JQ + atnmShift, 
                                        atomK=dihdrl.KQ + atnmShift, 
                                        atomL=dihdrl.LQ + atnmShift)
        for dihdrl in top.IMPDIHEDRALH.content:
            dihdrlType = top.IMPDIHEDRALTYPE.content[dihdrl.ICQH - 1]
            retTop.add_new_torsiondihedral(CQ=dihdrlType.CQ, 
                                        Q0=dihdrlType.Q0, 
                                        atomI=dihdrl.IQH + atnmShift, 
                                        atomJ=dihdrl.JQH + atnmShift, 
                                        atomK=dihdrl.KQH + atnmShift, 
                                        atomL=dihdrl.LQH + atnmShift,
                                        includesH=True)

        # add SOLUTEMOLECULES
        retTop.SOLUTEMOLECULES.content[0][0] = str(int(retTop.SOLUTEMOLECULES.content[0][0])+int(top.SOLUTEMOLECULES.content[0][0]))
        for solmol in top.SOLUTEMOLECULES.content[1:]:
            retTop.SOLUTEMOLECULES.content.append([str(int(solmol[0]) + atnmShift)])

        # add TEMPERATUREGROUPS
        retTop.TEMPERATUREGROUPS.content[0][0] = str(int(retTop.TEMPERATUREGROUPS.content[0][0])+int(top.SOLUTEMOLECULES.content[0][0]))
        for solmol in top.TEMPERATUREGROUPS.content[1:]:
            retTop.TEMPERATUREGROUPS.content.append([str(int(solmol[0]) + atnmShift)])

        # add PRESSUREGROUPS
        retTop.PRESSUREGROUPS.content[0][0] = str(int(retTop.PRESSUREGROUPS.content[0][0])+int(top.SOLUTEMOLECULES.content[0][0]))
        for solmol in top.PRESSUREGROUPS.content[1:]:
            retTop.PRESSUREGROUPS.content.append([str(int(solmol[0]) + atnmShift)])

        return retTop

    def read_file(self):
        #Read blocks to string
        data = parser.read_general_gromos_file(self._orig_file_path)

        #translate the string subblocks
        blocks = {}
        for block_title in data:
            #print(block_title)
            self.add_block(blocktitle=block_title, content=data[block_title])
            blocks.update({block_title: self.__getattribute__(block_title)})
        return blocks

    def make_ordered(self, orderList:list=None):
        if orderList:
            self._block_order = orderList
        else:
            self._block_order = ["TITLE", "PHYSICALCONSTANTS","TOPVERSION","ATOMTYPENAME","RESNAME","SOLUTEATOM","BONDSTRETCHTYPE","BONDH","BOND","BONDANGLEBENDTYPE","BONDANGLEH","BONDANGLE","IMPDIHEDRALTYPE","IMPDIHEDRALH","IMPDIHEDRAL","TORSDIHEDRALTYPE","DIHEDRALH","DIHEDRAL","CROSSDIHEDRALH","CROSSDIHEDRAL","LJPARAMETERS","SOLUTEMOLECULES","TEMPERATUREGROUPS","PRESSUREGROUPS","LJEXCEPTIONS","SOLVENTATOM","SOLVENTCONSTR"]

    def add_new_atomtype(self, name:str, verbose=False):
        if not hasattr(self, "ATOMTYPENAME"):
            defaultContent=['0', 'Dummy']
            self.add_block(blocktitle="ATOMTYPENAME", content=defaultContent, verbose=verbose)
            self.ATOMTYPENAME.content.append([name])
            self.ATOMTYPENAME.content.remove(['Dummy'])
        else:
            if len(self.ATOMTYPENAME.content) < 1:
                self.ATOMTYPENAME.content.append(["0"])
            self.ATOMTYPENAME.content.append([name])
        self.ATOMTYPENAME.content[0][0] = str(int(self.ATOMTYPENAME.content[0][0])+1)

    def add_new_resname(self, name:str, verbose=False):
        if not hasattr(self, "RESNAME"):
            defaultContent=['0', 'Dummy']
            self.add_block(blocktitle="RESNAME", content=defaultContent, verbose=verbose)
            self.RESNAME.content.append([name])
            self.RESNAME.content.remove(['Dummy'])
        else:
            if len(self.RESNAME.content) < 1:
                self.RESNAME.content.append(["0"])
            self.RESNAME.content.append([name])
        self.RESNAME.content[0][0] = str(int(self.RESNAME.content[0][0])+1)

    def add_new_soluteatom(self, ATNM:int=0, MRES:int=0, PANM:str="", IAC:int=0, MASS:float=0, CG:float=0, CGC:int=0, INE:list=[], INE14:list=[], verbose=False):
        if not hasattr(self, "SOLUTEATOM"):
            self.add_block(blocktitle="SOLUTEATOM", content=[], verbose=verbose)
            self.SOLUTEATOM.NRP = 0
        # some auto set methods
        if ATNM == 0:
            ATNM = len(self.SOLUTEATOM.content) + 1
        if MRES == 0:
            if len(self.SOLUTEATOM.content) >= 1:
                MRES = self.SOLUTEATOM.content[-1].MRES + 1
            else:
                MRES = 1
        #create new entry
        entry = blocks.soluteatom_type(ATNM=ATNM, MRES=MRES, PANM=PANM, IAC=IAC, MASS=MASS, CG=CG, CGC=CGC, INE=len(INE), INEvalues=INE, INE14=len(INE14), INE14values=INE14)
        self.SOLUTEATOM.content.append(entry)
        self.SOLUTEATOM.NRP += 1


    def add_new_bond(self, k:float, b0:float, atomI:int, atomJ:int, includesH:bool = False, verbose=False):
        #check if all classes are ready, if not create
        if not hasattr(self, "BONDSTRETCHTYPE"):
            self.add_block(blocktitle="BONDSTRETCHTYPE", content=list(), verbose=verbose)
        if includesH:
            if not hasattr(self, "BONDH"):
                self.add_block(blocktitle="BONDH", content=list(), verbose=verbose)
        else:
            if not hasattr(self, "BOND"):
                self.add_block(blocktitle="BOND", content=list(), verbose=verbose)
        
        
        # find the bondstretchtype number or create new bondstretchtype
        # TODO: add quartic force (CB)
        bond_type_number = 0
        iterator = 1
        quartic = k/(2*(b0**2))
        newBondStretchType = blocks.bondstretchtype_type(CB=quartic, CHB=k, B0=b0)
        for bond_type in self.BONDSTRETCHTYPE.content:
            if bond_type.CHB == newBondStretchType.CHB and bond_type.B0 == newBondStretchType.B0:
                break
            else:
                iterator += 1
        bond_type_number = iterator
        if iterator > len(self.BONDSTRETCHTYPE.content):#bond type was not found -> add new bondtype
            self.BONDSTRETCHTYPE.content.append(newBondStretchType)
            self.BONDSTRETCHTYPE.NBTY += 1

        #create new bond TODO: maybe check if already exists. But I will asume smart users
        newBond = blocks.top_bond_type(IB=atomI, JB=atomJ, ICB=bond_type_number)

        #check if we are adding a bond to BOND or BONDH
        if includesH:
            self.BONDH.content.append(newBond)
            self.BONDH.NBONH += 1
        else:
            self.BOND.content.append(newBond)
            self.BOND.NBON += 1

    def add_new_angle(self, k:float, kh:float, b0:float, atomI:int, atomJ:int, atomK:int, includesH:bool = False, verbose=False):
        #check if all classes are ready, if not create
        if not hasattr(self, "BONDANGLEBENDTYPE"):
            self.add_block(blocktitle="BONDANGLEBENDTYPE", content=[], verbose=verbose)
        if includesH:
            if not hasattr(self, "BONDANGLEH"):
                self.add_block(blocktitle="BONDANGLEH", content=[], verbose=verbose)
        else:
            if not hasattr(self, "BONDANGLE"):
                self.add_block(blocktitle="BONDANGLE", content=[], verbose=verbose)
        
        # find the BONDANGLEBENDTYPE number or create new BONDANGLEBENDTYPE
        # TODO: add harmonic in the angle cosine force (CT)
        angle_type_number = 0
        iterator = 1
        for angle_type in self.BONDANGLEBENDTYPE.content:
            if angle_type.CB == k and angle_type.B0 == b0:
                break
            else:
                iterator += 1
        angle_type_number = iterator
        if iterator > len(self.BONDANGLEBENDTYPE.content):#angle type was not found -> add new bondtype
            newBONDANGLEBENDTYPE = blocks.bondstretchtype_type(CB=k, CHB=kh, B0=b0)
            self.BONDANGLEBENDTYPE.content.append(newBONDANGLEBENDTYPE)
            self.BONDANGLEBENDTYPE.NBTY += 1
        
        #create new angle TODO: maybe check if already exists. But I will asume smart users
        newAngle = blocks.bondangle_type(IT=atomI, JT=atomJ, KT=atomK, ICT=angle_type_number)
        #check if we are adding a bond to BONDANGLE or BONDANGLEH
        if includesH:
            self.BONDANGLEH.content.append(newAngle)
            self.BONDANGLEH.NTHEH += 1
        else:
            self.BONDANGLE.content.append(newAngle)
            self.BONDANGLE.NTHE += 1

    def add_new_torsiondihedral(self, CP:float, PD:float, NP:int, atomI:int, atomJ:int, atomK:int, atomL:int, includesH:bool = False, verbose=False):
        #check if all classes are ready, if not create
        if not hasattr(self, "TORSDIHEDRALTYPE"):
            self.add_block(blocktitle="TORSDIHEDRALTYPE", content=[], verbose=verbose)
        if includesH:
            if not hasattr(self, "DIHEDRALH"):
                self.add_block(blocktitle="DIHEDRALH", content=[], verbose=verbose)
        else:
            if not hasattr(self, "DIHEDRAL"):
                self.add_block(blocktitle="DIHEDRAL", content=[], verbose=verbose)
        
        # find the TORSDIHEDRALTYPE number or create new TORSDIHEDRALTYPE
        torsion_type_number = 0
        iterator = 1
        for torsion_type in self.TORSDIHEDRALTYPE.content:
            if torsion_type.CP == CP and torsion_type.PD == PD and torsion_type.NP == NP:
                break
            else:
                iterator += 1
        torsion_type_number = iterator #found the torsion
        if iterator > len(self.TORSDIHEDRALTYPE.content):#torsion type was not found -> add new bondtype
            newTORSDIHEDRALTYPE = blocks.torsdihedraltype_type(CP=CP, PD=PD, NP=NP)
            self.TORSDIHEDRALTYPE.content.append(newTORSDIHEDRALTYPE)
            self.TORSDIHEDRALTYPE.NPTY += 1
        
        #check if we are adding a bond to DIHEDRAL or DIHEDRALH
        if includesH:
            self.DIHEDRALH.content.append(blocks.dihedralh_type(IPH=atomI, JPH=atomJ, KPH=atomK, LPH=atomL, ICPH=torsion_type_number))
            self.DIHEDRALH.NPHIH += 1
        else:
            self.DIHEDRAL.content.append(blocks.top_dihedral_type(IP=atomI, JP=atomJ, KP=atomK, LP=atomL, ICP=torsion_type_number))
            self.DIHEDRAL.NPHI += 1

    
    def add_new_impdihedral_type(self, CQ:float, Q0:float, verbose=False):
        #check if all classes are ready, if not create
        if not hasattr(self, "IMPDIHEDRALTYPE"):
            self.add_block(blocktitle="IMPDIHEDRALTYPE", content=[], verbose=verbose)
        newIMPDIHEDRALTYPE = blocks.impdihedraltype_type(CQ=CQ, Q0=Q0)
        self.IMPDIHEDRALTYPE.content.append(newIMPDIHEDRALTYPE)
        self.IMPDIHEDRALTYPE.NQTY += 1


    def add_new_impdihedral(self, CQ:float, Q0:float, atomI:int, atomJ:int, atomK:int, atomL:int, includesH:bool = False, verbose=False):
        #check if all classes are ready, if not create
        if not hasattr(self, "IMPDIHEDRALTYPE"):
            self.add_block(blocktitle="IMPDIHEDRALTYPE", content=[], verbose=verbose)
        if includesH:
            if not hasattr(self, "IMPDIHEDRALH"):
                self.add_block(blocktitle="IMPDIHEDRALH", content=[], verbose=verbose)
        else:
            if not hasattr(self, "IMPDIHEDRAL"):
                self.add_block(blocktitle="IMPDIHEDRAL", content=[], verbose=verbose)
        
        # find the IMPDIHEDRALTYPE number or create new IMPDIHEDRALTYPE
        impdihedral_type_number = 1
        iterator = 1
        for imp_type in self.IMPDIHEDRALTYPE.content:
            if imp_type.CQ == CQ and imp_type.Q0 == Q0:
                break
            else:
                iterator += 1
            impdihedral_type_number = iterator #found the torsion
        if iterator > len(self.IMPDIHEDRALTYPE.content):#torsion type was not found -> add new bondtype
            self.add_new_impdihedral_type(CQ=CQ, Q0=Q0)
        
        #check if we are adding a bond to IMPDIHEDRALH or IMPDIHEDRALH
        if includesH:
            self.IMPDIHEDRALH.content.append(blocks.impdihedralh_type(IQH=atomI, JQH=atomJ, KQH=atomK, LQH=atomL, ICQH=impdihedral_type_number))
            self.IMPDIHEDRALH.NQHIH += 1
        else:
            self.IMPDIHEDRAL.content.append(blocks.impdihedral_type(IQ=atomI, JQ=atomJ, KQ=atomK, LQ=atomL, ICQ=impdihedral_type_number))
            self.IMPDIHEDRAL.NQHI += 1


    #TODO: add implementation
    def add_new_crossdihedral(self, verbose=False):
        raise NotImplementedError("Who needs this???? Could you plox implement it. UwU")

    def add_new_LJparameter(self, C6:float, C12:float, CS6:float=0, CS12:float=0, combination_rule:str="geometric", verbose=False, AddATOMTYPENAME:str=None, lowerBound:float=1e-100):
        if not hasattr(self, "LJPARAMETERS"):
            self.add_block(blocktitle="LJPARAMETERS", content=[], verbose=verbose)
            self.LJPARAMETERS.NRATT2 = 0
        #safety
        if C6 < lowerBound:
            C6 = lowerBound
        if C12 < lowerBound:
            C12 = lowerBound
        if CS6 < lowerBound:
            CS6 = lowerBound
        if CS12 < lowerBound:
            CS12 = lowerBound
        # add LJ parameter for all existing combinations
        num=0
        nratt=int((math.sqrt(8*self.LJPARAMETERS.NRATT2+1)-1)/2)
        for i in range(nratt):
            if combination_rule == "geometric":
                c6 = math.sqrt(float(C6 * self.LJPARAMETERS.content[num].C6))
                c12 = math.sqrt(float(C12 * self.LJPARAMETERS.content[num].C12))
                cs6 = math.sqrt(float(CS6 * self.LJPARAMETERS.content[num].CS6))
                cs12 = math.sqrt(float(CS12 * self.LJPARAMETERS.content[num].CS12))
            else:
                raise NotImplementedError("Error in add_new_LJparameter: desired combination rule not implemented")
            add = blocks.ljparameters_type(IAC=i+1, JAC=nratt+1, C6=c6, C12=c12, CS12=cs12, CS6=cs6)
            self.LJPARAMETERS.append(add)
            num += i+2
        #add new LJ paramter to self
        add = blocks.ljparameters_type(IAC=nratt+1, JAC=nratt+1, C6=C6, C12=C12, CS12=CS12, CS6=CS6)
        self.LJPARAMETERS.append(add)
        self.LJPARAMETERS.NRATT2 += nratt + 1
        
        if AddATOMTYPENAME != None:
            if not hasattr(self, "ATOMTYPENAME"):
                self.add_block(blocktitle="ATOMTYPENAME", content=[], verbose=verbose)
                self.LJPARAMETERS.NRATT = 0
            self.add_new_atomtype(AddATOMTYPENAME)
            if(int(self.ATOMTYPENAME.content[0][0]) != self.LJPARAMETERS.content[-1].IAC):
                raise IndexError("Missmatch between number of ATOMTYPNAMEs and LJPARAMETERS")


    def find_LJparameterNumber(self, C12:float, C6:float) -> int:
        if not hasattr(self, "LJPARAMETERS"):
            return 0
        elif self.LJPARAMETERS.NRATT2 < 1:
            return 0
        else:
            for lj in self.LJPARAMETERS.content:
                if C12 == lj.C12 and C6 == lj.C6:
                    return lj.IAC
            return 0 # LJ parameter not found


    def add_new_SOLUTEATOM(self, ATNM:int, MRES:int=1, PANM:str='_', IAC:int=1, MASS:float=1.0, CG:int=0, CGC:int=0, INE:list=[], INE14:list=[], verbose=None, C6:float=None, C12:float=None, CS6:float=0, CS12:float=0, IACname:str=None):
        if not hasattr(self, "SOLUTEATOM"):
             self.add_block(blocktitle="SOLUTEATOM", content=[], verbose=verbose)
             self.SOLUTEATOM.NRP = 0
        
        # Find IAC and (if needed) add a new LJ Parameter
        if C6 != None or C12 != None:           #need to find PANM and IAC
            if hasattr(self, "LJPARAMETERS"):
                IAC = self.find_LJparameterNumber(C6=C6, C12=C12)
                if IAC == 0: #IAC not found -> add new LJ parameter
                    self.add_new_LJparameter(C6=C6, C12=C12, CS6=CS6, CS12=CS12, verbose=verbose, AddATOMTYPENAME=IACname)
                    IAC = self.LJPARAMETERS.content[-1].IAC
            else:
                self.add_new_LJparameter(C6=C6, C12=C12, CS6=CS6, CS12=CS12, verbose=verbose, AddATOMTYPENAME=IACname)
                IAC = 1
        
        # IAC should be known at this point -> we can search for PANM if not known
        #if PANM == None:
        #    if not (IAC >=1):
        #        raise "You miss treated your IAC or created a different unexpected error"
        #    else:
        #        if not hasattr(self, "ATOMTYPENAME"):
        #            raise "How did you think we could find PANM if ATOMTYPENAME does not even exist"
        #        elif len(self.ATOMTYPENAME.content) <= IAC:
        #            raise "The desired IAC is not yet written into ATOMTYPENAME"
        #        else:
        #            PANM = self.ATOMTYPENAME.content[IAC][0]

        #TODO: Maybe add further automation for ATNM, MRES, MASS, CG, ...
        #Now all variables of the new SOLUTEATOM should be known
        newSoluteAtom = blocks.soluteatom_type(ATNM=ATNM, MRES=MRES, PANM=PANM, IAC=IAC, MASS=MASS, CG=CG, CGC=CGC, INE=len(INE), INEvalues=INE, INE14=len(INE14), INE14values=INE14)
        self.SOLUTEATOM.content.append(newSoluteAtom)
        self.SOLUTEATOM.NRP += 1

    def add_new_CONSTRAINT(self, IC:int, JC:int, ICC:float, verbose=False):
        """
        adds a CONSTRAINT entry to the topology

        Parameters
        ----------
        IC : int
            atom index I
        JC : int
            atom index J
        ICC : float
            constraint length 
        verbose : bool, optional
        """
        if not hasattr(self, "CONSTRAINT"):
            self.add_block(blocktitle="CONSTRAINT", content=[], verbose=verbose)
            self.CONSTRAINT.NCON = 0
        if not hasattr(self, "BONDSTRETCHTYPE"):
            self.add_block(blocktitle="BONDSTRETCHTYPE", content=list(), verbose=verbose)
        
        # find the bondstretchtype number or create new bondstretchtype
        bond_type_number = 0
        iterator = 1
        newBondStretchType = blocks.bondstretchtype_type(CB=1, CHB=1, B0=ICC)
        for bond_type in self.BONDSTRETCHTYPE.content:
            if bond_type.B0 == newBondStretchType.B0:
                break
            else:
                iterator += 1
        bond_type_number = iterator
        if iterator > len(self.BONDSTRETCHTYPE.content):#bond type was not found -> add new bondtype
            self.BONDSTRETCHTYPE.content.append(newBondStretchType)
            self.BONDSTRETCHTYPE.NBTY += 1
        self.CONSTRAINT.content.append(blocks.constraint_type(IC=IC, JC=JC, ICC=bond_type_number))
        self.CONSTRAINT.NCON += 1

    def get_mass(self) -> float:
        """
        Calculates the total mass of the solute molecule

        Returns
        -------
        float
            total mass in a.u.
        """
        mass = 0
        if hasattr(self, "SOLUTEATOM"):
            for i in self.SOLUTEATOM.content:
                mass += i.MASS
        return mass
        

        

