import re
from enum import Enum
from typing import Union, Iterable, List
import inspect
import math

from pygromos.files.blocks.topology_blocks import FORCEFIELD, MAKETOPVERSION, TITLE
from pygromos.files.blocks._general_blocks import _generic_gromos_block, _iterable_gromos_block, _generic_field


class mtb_atoms_field(_generic_field):
    def __init__(
        self, ATOM: int, ANM: str, IACM: int, MASS: float, CGMI: float, CGM: int, MAE: int, MSAE: list[float]
    ) -> None:
        self.ATOM = int(ATOM)
        self.ANM = ANM
        self.IACM = int(IACM)
        self.MASS = float(MASS)
        self.CGMI = float(CGMI)
        self.CGM = int(CGM)
        self.MAE = int(MAE)
        self.MSAE = [float(x) for x in MSAE]

    def to_string(self) -> str:
        return_str = ""
        return_str += self.fieldseperator + str(self.ATOM)
        return_str += self.fieldseperator + self.ANM
        return_str += self.fieldseperator + str(self.IACM)
        return_str += self.fieldseperator + str(self.MASS)
        return_str += self.fieldseperator + str(self.CGMI)
        return_str += self.fieldseperator + str(self.CGM)
        return_str += self.fieldseperator + str(self.MAE)
        for iter in self.MSAE:
            return_str += "\t" + str(iter).strip()
            lcounter += 1
            if (lcounter % 6) == 0 and len(self.MSAE) > 6:
                str_line += "\n\t\t\t\t\t\t\t\t\t\t"
        return_str += self.lineseperator
        return return_str


class mtb_bonds_field(_generic_field):
    def __init__(self, IB: int, JB: int, MCB: int) -> None:
        self.IB = int(IB)
        self.JB = int(JB)
        self.MCB = int(MCB)

    def to_string(self) -> str:
        return_str = ""
        return_str += self.fieldseperator + str(self.IB)
        return_str += self.fieldseperator + str(self.JB)
        return_str += self.fieldseperator + str(self.MCB)
        return_str += self.lineseperator
        return return_str


class mtb_angles_field(_generic_field):
    def __init__(self, IB: int, JB: int, KB: int, MCB: int) -> None:
        self.IB = int(IB)
        self.JB = int(JB)
        self.KB = int(KB)
        self.MCB = int(MCB)

    def to_string(self) -> str:
        return_str = ""
        return_str += self.fieldseperator + str(self.IB)
        return_str += self.fieldseperator + str(self.JB)
        return_str += self.fieldseperator + str(self.KB)
        return_str += self.fieldseperator + str(self.MCB)
        return_str += self.lineseperator
        return return_str


class mtb_dihedral_field(_generic_field):
    def __init__(self, IB: int, JB: int, KB: int, LB: int, MCB: int) -> None:
        self.IB = int(IB)
        self.JB = int(JB)
        self.KB = int(KB)
        self.LB = int(LB)
        self.MCB = int(MCB)

    def to_string(self) -> str:
        return_str = ""
        return_str += self.fieldseperator + str(self.IB)
        return_str += self.fieldseperator + str(self.JB)
        return_str += self.fieldseperator + str(self.KB)
        return_str += self.fieldseperator + str(self.LB)
        return_str += self.fieldseperator + str(self.MCB)
        return_str += self.lineseperator
        return return_str


class mtb_lj_exceptions_field(_generic_field):
    def __init__(self, iac: int, jac: int, mcb: int):
        self.iac = int(iac)
        self.jac = int(jac)
        self.mcb = int(mcb)

    def to_string(self) -> str:
        return_str = ""
        return_str += self.fieldseperator + str(self.iac)
        return_str += self.fieldseperator + str(self.jac)
        return_str += self.fieldseperator + str(self.mcb)
        return_str += self.lineseperator
        return return_str


class MTBUILDBLSOLUTE(_generic_gromos_block):
    FORCEFIELD: FORCEFIELD
    MAKETOPVERSION: MAKETOPVERSION
    atoms: List[mtb_atoms_field]

    def __init__(self, FORCEFIELD: FORCEFIELD = None, MAKETOPVERSION: MAKETOPVERSION = None, content=None):
        super().__init__(name=self.__class__.__name__, used=True, content=content)
        self.FORCEFIELD = FORCEFIELD
        self.MAKETOPVERSION = MAKETOPVERSION

    def read_content_from_str(self, content: str):
        # reset all storage
        self.atoms = []
        self.bonds = []

        # first line
        first_line = content[0].split()
        self.filename = first_line[1]
        self.residuecode = first_line[3]
        self.function = first_line[4]
        self.type = first_line[6]
        self.fullname = first_line[8]

        self.RNME = content[3].strip()

        self.NMAT = int(content[6].split()[0])
        self.NLIN = int(content[6].split()[1])

        # TODO: implement NLIN
        itr = 6 + 4
        nmat_found = 0

        # try to find all atoms in the ATOM subblock and hope the while does not go to far...
        while nmat_found < self.NMAT and itr < (len(content) - 10):
            if content[itr].startswith("#"):
                itr += 1
                continue
            else:
                dump1 = content[itr].strip().split()
                atom, anm, iacm, mass, cgm, icgm, mae = dump1[0:7]
                if 1 <= int(mae):
                    msae_values = [int(i) for i in dump1[7:]]
                    # keep reading in lines until we have all the data needed.
                    while int(mae) > len(msae_values):
                        itr += 1
                        try:
                            if any(i in content[itr] for i in ["\t\t\t\t\t", "                   "]):
                                msae_values.extend([int(i) for i in content[itr].pop(0).strip().split()])
                        except:
                            raise IOError("Problem reading MSAE for anm=" + str(anm) + " mae=" + str(mae))
                            break
                else:
                    msae_values = []
                self.atoms.append(mtb_atoms_field(atom, anm, iacm, mass, cgm, icgm, mae, msae_values))
            itr += 1
        # all atoms from ATOM subblock should be found and parsed at this point

        # TODO: parse trailing atoms subblock

        # read bonds
        while itr < len(content):
            if content[itr].startswith("#"):
                itr += 1
                continue
            else:
                try:
                    self.NB = int(content[itr].strip())
                    itr += 1
                    break
                except:
                    raise IOError("Problem reading number of bonds")
        bonds_found = 0
        while itr < len(content) and bonds_found < self.NB:
            if content[itr].startswith("#"):
                itr += 1
                continue
            else:
                dump1 = content[itr].strip().split()
                ib, jb, mcb = dump1[0:3]
                self.bonds.append(mtb_bonds_field(ib, jb, mcb))
                bonds_found += 1
                itr += 1

        # read angles
        while itr < len(content):
            if content[itr].startswith("#"):
                itr += 1
                continue
            else:
                try:
                    self.NBA = int(content[itr].strip())
                    itr += 1
                    break
                except:
                    raise IOError("Problem reading number of angles")
        angles_found = 0
        while itr < len(content) and angles_found < self.NBA:
            if content[itr].startswith("#"):
                itr += 1
                continue
            else:
                dump1 = content[itr].strip().split()
                ib, jb, kb, mcb = dump1[0:4]
                self.angles.append(mtb_angles_field(ib, jb, kb, mcb))
                angles_found += 1
                itr += 1

        # read improper dihedrals
        while itr < len(content):
            if content[itr].startswith("#"):
                itr += 1
                continue
            else:
                try:
                    self.NIDA = int(content[itr].strip())
                    itr += 1
                    break
                except:
                    raise IOError("Problem reading number of improper dihedrals")
        improper_dihedrals_found = 0
        while itr < len(content) and improper_dihedrals_found < self.NIDA:
            if content[itr].startswith("#"):
                itr += 1
                continue
            else:
                dump1 = content[itr].strip().split()
                ib, jb, kb, lb, mcb = dump1[0:5]
                self.improper_dihedrals.append(mtb_dihedral_field(ib, jb, kb, lb, mcb))
                improper_dihedrals_found += 1
                itr += 1

        # read dihedrals
        while itr < len(content):
            if content[itr].startswith("#"):
                itr += 1
                continue
            else:
                try:
                    self.NDA = int(content[itr].strip())
                    itr += 1
                    break
                except:
                    raise IOError("Problem reading number of dihedrals")
        dihedrals_found = 0
        while itr < len(content) and dihedrals_found < self.NDA:
            if content[itr].startswith("#"):
                itr += 1
                continue
            else:
                dump1 = content[itr].strip().split()
                ib, jb, kb, lb, mcb = dump1[0:5]
                self.dihedrals.append(mtb_dihedral_field(ib, jb, kb, lb, mcb))
                dihedrals_found += 1
                itr += 1

        # read LJ exceptions
        while itr < len(content):
            if content[itr].startswith("#"):
                itr += 1
                continue
            else:
                try:
                    self.NEX = int(content[itr].strip())
                    itr += 1
                    break
                except:
                    raise IOError("Problem reading number of LJ exceptions")
        lj_exceptions_found = 0
        while itr < len(content) and lj_exceptions_found < self.NEX:
            if content[itr].startswith("#"):
                itr += 1
                continue
            else:
                dump1 = content[itr].strip().split()
                iac, jac, mcb = dump1[0:3]
                self.lj_exceptions.append(mtb_lj_exceptions_field(iac, jac, mcb))
                lj_exceptions_found += 1
                itr += 1

    def block_to_string(self) -> str:
        result = "MTBUILDBLSOLUTE" + self.line_seperator
        result += (
            "#@BLOCKTYPE "
            + self.filename
            + " BLK "
            + self.residuecode
            + " "
            + self.function
            + " TYPE "
            + self.type
            + " NAME "
            + self.fullname
            + self.line_seperator
        )
        result += "# building block" + self.line_seperator

        result += "# RNME" + self.line_seperator
        result += self.RNME + self.line_seperator
        result += "# number of atoms, number of preceding exclusions" + self.line_seperator
        result += "# NMAT NLIN" + self.line_seperator
        result += str(self.NMAT) + self.field_seperator + str(self.NLIN) + self.line_seperator
        result += "# preceding exclusions" + self.line_seperator

        result += "#ATOM                               MAE MSAE" + self.line_seperator
        result += "# atoms" + self.line_seperator
        result += "#ATOM ANM  IACM MASS        CGMICGM MAE MSAE" + self.line_seperator
        for atom in self.atoms:
            result += atom.to_string()

        result += "# trailing atoms" + self.line_seperator
        result += "#ATOM ANM  IACM MASS        CGMICGM" + self.line_seperator
        # TODO: add trailing atoms

        result += "# bonds" + self.line_seperator
        result += "#NB" + self.line_seperator
        result += str(self.NB) + self.line_seperator
        result += "#  IB   JB  MCB" + self.line_seperator
        for bond in self.bonds:
            result += bond.to_string()

        result += "# angles" + self.line_seperator
        result += "#NBA" + self.line_seperator
        result += str(self.NBA) + self.line_seperator
        result += "#  IB   JB   KB  MCB" + self.line_seperator
        for angle in self.angles:
            result += angle.to_string()

        result += "# improper dihedrals" + self.line_seperator
        result += "#NIDA" + self.line_seperator
        result += str(self.NIDA) + self.line_seperator
        result += "#  IB   JB   KB   LB  MCB" + self.line_seperator
        for dihedral in self.improper_dihedrals:
            result += dihedral.to_string()

        result += "# dihedrals" + self.line_seperator
        result += "#NDA" + self.line_seperator
        result += str(self.NDA) + self.line_seperator
        result += "#  IB   JB   KB   LB  MCB" + self.line_seperator
        for dihedral in self.dihedrals:
            result += dihedral.to_string()

        result += "# LJ exceptions" + self.line_seperator
        result += "#NEX" + self.line_seperator
        result += str(self.NEX) + self.line_seperator
        result += "# IAC  JAC  MCB" + self.line_seperator
        for lj_exception in self.lj_exceptions:
            result += lj_exception.to_string()

        result += "END" + self.line_seperator
        return result
