from pygromos.files.blocks._general_blocks import TITLE
from pygromos.files.blocks._general_blocks import _generic_gromos_block

# forward declarations
TITLE: TITLE = TITLE

class QMZONE(_generic_gromos_block):
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