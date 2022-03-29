from pygromos.files.blocks._general_blocks import _generic_field, TITLE as generic_TITLE
from pygromos.files.blocks.topology_blocks import _topology_table_block, _topology_block
from pygromos.utils.typing import (
    Iterable,
    Union,
    QMUNIT_Type,
    MNDOELEMENTS_Type,
    TURBOMOLEELEMENTS_Type,
    DFTBELEMENTS_Type,
    MOPACELEMENTS_Type,
    ORCAELEMENTS_Type,
    XTBELEMENTS_Type,
)

# Note that while many classes in this file inherit from _topology_block or _topology_table_block
# there is no obvious connection other than the tabular structure of .qmmm and .top files

# forward declarations
TITLE: generic_TITLE = generic_TITLE


class qmzone_field(_generic_field):
    def __init__(self, QMEN: str, QMEI: int, QMEZ: int, QMEB: int):
        """
        One line of the QMZONE block

        Parameters
        ----------
        QMEN : str
            (QM element name), indicates the atom identifier for this position
        QMEI : int
            (QM element iterator), specifies an iterator over the atom positions
        QMEZ : int
            (QM element Z), specifies the nuclear charge Z of the atom position
        QMEB : int
            (QM element bond), specifies whether a bond can be broken or not (== 0)
        """
        self.QMEN = str(QMEN)
        self.QMEI = int(QMEI)
        self.QMEZ = int(QMEZ)
        self.QMEB = int(QMEB)

    def to_string(self):
        # the first 24 characters or so are ignored by Gromos
        # use spaces instead of tabs and format the resulting
        # table nicely using ljust and rjust
        str_line = (
            str(self.QMEN).ljust(27)
            + str(self.QMEI).rjust(6)
            + str(self.QMEZ).rjust(6)
            + str(self.QMEB).rjust(6)
            + "\n"
        )
        return str_line


class QMZONE(_topology_table_block):
    NBON: int = 1
    table_header: Iterable[str] = ["QMEN", "QMEI", "QMEZ", "QMEB"]
    table_line_type = qmzone_field

    def __init__(self, content: Union[Iterable[qmzone_field], str]):
        super().__init__(content=content, FORCEFIELD=None, MAKETOPVERSION=None)

    def block_to_string(self) -> str:
        result = self.name + "\n"  # QMZONE
        result += f"{self.comment_char} {self.table_header[0]}".ljust(27)
        for element in self.table_header[1:]:
            result += element.rjust(6)
        result += "\n"
        for element in self.content:
            result += element.to_string()
        result += "END\n"
        return result


class QMUNIT(_topology_block):

    table_header: Iterable[str] = ["QLGL", "QEGE", "QCGC", "QIGI"]

    """
        The QMUNIT block

        Parameters
        ----------
        QLGL : float
             QM length to Gromos length (e.g. Bohr to nm)
        QEGE : float
             QM energy to Gromos energy (e.g. Hartree to kJ / mol)
        QCGC : float
             Gromos charge to QM charge
        QIGI : float
            QM input units to Gromos input units (e.g. Angstrom to nm)
        """

    def __init__(
        self,
        content: Union[str, dict, QMUNIT_Type],
        QLGL: float = 0.052918,
        QEGE: float = 2625.50,
        QCGC: float = 1.0,
        QIGI: float = 0.1,
    ):
        super().__init__(FORCEFIELD=None, MAKETOPVERSION=None, content=content)
        if content is None:
            self.QLGL = float(QLGL)
            self.QEGE = float(QEGE)
            self.QCGC = float(QCGC)
            self.QIGI = float(QIGI)

    def read_content_from_str(self, content: str):
        super().read_content_from_str(content)
        try:
            self.QLGL = float(self.content[0][0])
            self.QEGE = float(self.content[0][1])
            self.QCGC = float(self.content[0][2])
            self.QIGI = float(self.content[0][3])
        except IOError as e:
            print("Error while reading QMUNIT block: " + str(e))

    def block_to_string(self) -> str:
        result = self.name + "\n"  # QMUNIT
        result += f"{self.comment_char} {self.table_header[0]}".ljust(17)
        for element in self.table_header[1:]:
            result += element.ljust(15)
        result += "\n"
        result += (
            str(self.QLGL)
            + self.field_seperator
            + str(self.QEGE)
            + self.field_seperator
            + str(self.QCGC)
            + self.field_seperator
            + str(self.QIGI)
            + self.line_seperator
        )
        result += "END\n"
        return result


# There should be only of these blocks in the .qmmm file


class MNDOELEMENTS(_topology_block):
    def __init__(self, content: Union[str, dict, MNDOELEMENTS_Type]):
        super().__init__(FORCEFIELD=None, MAKETOPVERSION=None, content=content)


class TURBOMOLEELEMENTS(_topology_block):
    def __init__(self, content: Union[str, dict, TURBOMOLEELEMENTS_Type]):
        super().__init__(FORCEFIELD=None, MAKETOPVERSION=None, content=content)


class DFTBELEMENTS(_topology_block):
    def __init__(self, content: Union[str, dict, DFTBELEMENTS_Type]):
        super().__init__(FORCEFIELD=None, MAKETOPVERSION=None, content=content)


class MOPACELEMENTS(_topology_block):
    def __init__(self, content: Union[str, dict, MOPACELEMENTS_Type]):
        super().__init__(FORCEFIELD=None, MAKETOPVERSION=None, content=content)


class ORCAELEMENTS(_topology_block):
    def __init__(self, content: Union[str, dict, ORCAELEMENTS_Type]):
        super().__init__(FORCEFIELD=None, MAKETOPVERSION=None, content=content)


class XTBELEMENTS(_topology_block):
    def __init__(self, content: Union[str, dict, XTBELEMENTS_Type]):
        super().__init__(FORCEFIELD=None, MAKETOPVERSION=None, content=content)
