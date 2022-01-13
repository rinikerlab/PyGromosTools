from typing import Iterable, Union
from pygromos.files.blocks._general_blocks import TITLE, _generic_gromos_block, _generic_field
from pygromos.files.blocks.topology_blocks import _topology_table_block

# forward declarations
TITLE: TITLE = TITLE


class qmzone_field(_generic_field):
    def __init__(self, name: str, A: int, B: int, C: int):
        self.name = name
        self.A = A
        self.B = B
        self.C = C

    def to_string(self):
        str_line = "\t" + str(self.name) + "\t\t\t\t\t\t\t\t" + str(self.A) + "\t" + str(self.B) + "\t" + str(self.C) +"\n"
        return str_line

class QMZONE(_topology_table_block):
    NBON:int = 1
    table_header: Iterable[str] = ["IB", "JB", "ICB"]
    table_line_type = qmzone_field

    def __init__(self, content: Union[Iterable[qmzone_field], str]):
        super().__init__(content=content, FORCEFIELD=None, MAKETOPVERSION=None)

class QMUNIT(_generic_gromos_block):
    """System Block

        The system block defines the number of solute molecules and solvent molecules

    Attributes
    ----------
    NPM:    int
        Number of Solute Molecules
    NSM:    int
        Number of Solvent Molecules



    """
    name: str = "SYSTEM"

    # fields
    NPM: int  # number of Solute Molecules
    NSM: int  # number of Solvent Molecules

    _order = [[["NPM", "NSM"]]]

    def __init__(self, NPM:int=0, NSM:int=0, content=None):
        super().__init__(used=True, content=content)
        if content is None:
            self.NPM = int(NPM)
            self.NSM = int(NSM)

class XTBELEMENTS(_generic_gromos_block):
    """System Block

        The system block defines the number of solute molecules and solvent molecules

    Attributes
    ----------
    NPM:    int
        Number of Solute Molecules
    NSM:    int
        Number of Solvent Molecules



    """
    name: str = "SYSTEM"

    # fields
    NPM: int  # number of Solute Molecules
    NSM: int  # number of Solvent Molecules

    _order = [[["NPM", "NSM"]]]

    def __init__(self, NPM:int=0, NSM:int=0, content=None):
        super().__init__(used=True, content=content)
        if content is None:
            self.NPM = int(NPM)
            self.NSM = int(NSM)


class ORCAELEMENTS(_generic_gromos_block):
    """System Block

        The system block defines the number of solute molecules and solvent molecules

    Attributes
    ----------
    NPM:    int
        Number of Solute Molecules
    NSM:    int
        Number of Solvent Molecules



    """
    name: str = "SYSTEM"

    # fields
    NPM: int  # number of Solute Molecules
    NSM: int  # number of Solvent Molecules

    _order = [[["NPM", "NSM"]]]

    def __init__(self, NPM:int=0, NSM:int=0, content=None):
        super().__init__(used=True, content=content)
        if content is None:
            self.NPM = int(NPM)
            self.NSM = int(NSM)