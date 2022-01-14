from typing import Iterable, Union
from pygromos.files.blocks._general_blocks import TITLE, _generic_gromos_block, _generic_field
from pygromos.files.blocks.topology_blocks import _topology_table_block, _topology_block

# forward declarations
TITLE: TITLE = TITLE


class qmzone_field(_generic_field):
    def __init__(self, QMEN: str, QMEI: int, QMEZ: int, QMEB: int):
        # Write comment line
        # QMEN (QM element name), indicates the atom identifier for this position
        # QMEI (QM element iterator), specifies an iterator over the atom positions
        # QMEZ (QM element Z), specifies the nuclear charge Z of the atom position
        # QMEB (QM element bond), specifies whether a bond can be broken or not (== 0).
        self.QMEN = QMEN
        self.QMEI = QMEI
        self.QMEZ = QMEZ
        self.QMEB = QMEB

    def to_string(self):
        str_line = str(self.QMEN) + "\t\t\t\t\t" + str(self.QMEI) + "\t" + str(self.QMEZ) + "\t" + str(self.QMEB) +"\n"
        return str_line

class QMZONE(_topology_table_block):
    NBON:int = 1
    table_header: Iterable[str] = ["QMEN", "QMEI", "QMEZ", "QMEB"]
    table_line_type = qmzone_field

    def __init__(self, content: Union[Iterable[qmzone_field], str]):
        super().__init__(content=content, FORCEFIELD=None, MAKETOPVERSION=None)

    def block_to_string(self) -> str:
        result = self.name + "\n"
        result += super().block_to_string()
        result += "END\n"
        return result

class QMUNIT(_topology_block):
    def __init__(self, content:(str or dict or None or __class__)):
        super().__init__(FORCEFIELD=None, MAKETOPVERSION=None, content=content)

class XTBELEMENTS(_topology_block):
    def __init__(self, content:(str or dict or None or __class__)):
        super().__init__(FORCEFIELD=None, MAKETOPVERSION=None, content=content)

class ORCAELEMENTS(_topology_block):
    def __init__(self, content:(str or dict or None or __class__)):
        super().__init__(FORCEFIELD=None, MAKETOPVERSION=None, content=content)