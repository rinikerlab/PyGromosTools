from typing import Iterable, Union
from pygromos.files.blocks._general_blocks import TITLE, _generic_gromos_block, _generic_field
from pygromos.files.blocks.topology_blocks import _topology_table_block, _topology_block

# forward declarations
TITLE: TITLE = TITLE


class qmzone_field(_generic_field):
    def __init__(self, name: str, A: int, B: int, C: int):
        self.name = name
        self.A = A
        self.B = B
        self.C = C

    def to_string(self):
        str_line = str(self.name) + "\t\t\t\t\t" + str(self.A) + "\t" + str(self.B) + "\t" + str(self.C) +"\n"
        return str_line

class QMZONE(_topology_table_block):
    NBON:int = 1
    table_header: Iterable[str] = ["name", "A", "B", "C"]
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