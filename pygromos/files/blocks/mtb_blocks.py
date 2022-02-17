import re
from enum import Enum
from typing import Union, Iterable, List
import inspect
import math

from pygromos.files.blocks.topology_blocks import FORCEFIELD, MAKETOPVERSION, TITLE
from pygromos.files.blocks._general_blocks import _generic_gromos_block, _iterable_gromos_block, _generic_field

class MTBUILDBLSOLUTE(_generic_gromos_block):
    FORCEFIELD: FORCEFIELD
    MAKETOPVERSION: MAKETOPVERSION

    def __init__(self, FORCEFIELD:FORCEFIELD = None, MAKETOPVERSION:MAKETOPVERSION = None, content=None):
        super().__init__(name=self.__class__.__name__, used=True, content=content)
        self.FORCEFIELD = FORCEFIELD
        self.MAKETOPVERSION = MAKETOPVERSION

    def read_content_from_str(self, content: str):
        #first line
        first_line = content[0].split()
        self.filename = first_line[1]
        self.residuecode = first_line[3]
        self.function = first_line[4]
        self.type = first_line[6]
        self.fullname = first_line[8]

        self.RNME = content[3].strip()

        self.NMAT = int(content[6].split()[0])
        self.NLIN = int(content[6].split()[1])

        #TODO: implement NLIN
        itr = 6+4
        nmat_found = 0
        while nmat_found < self.NMAT and itr < (len(content)-10):
            if content[itr].startswith("#"):
                itr += 1
                continue
            else:
                atom, anm, iacm, mass, cgm, icgm, mae = content[itr].split()[7:]





    def block_to_string(self) -> str:
        result = "MTBUILDBLSOLUTE"+self.line_seperator
        result += "#@BLOCKTYPE "+self.filename+" BLK "+self.residuecode+" "+self.function+" TYPE "+self.type+" NAME "+self.fullname+self.line_seperator
        result += "# building block"+self.line_seperator
        result += "# RNME"+self.line_seperator
        result += self.RNME+self.line_seperator
        result += "# number of atoms, number of preceding exclusions"+self.line_seperator
        result += "# NMAT NLIN"+self.line_seperator
        result += str(self.NMAT)+self.field_seperator+str(self.NLIN)+self.line_seperator
        result += "# preceding exclusions"+self.line_seperator
        result += "#ATOM                               MAE MSAE"+self.line_seperator
        result += "# atoms"+self.line_seperator

        result += "END"+self.line_seperator
        return result