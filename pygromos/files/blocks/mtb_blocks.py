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

        result += "END" + self.line_seperator
        return result
