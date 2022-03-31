from pygromos.files.coord.cnf import Cnf
from pygromos.files.blocks import coord_blocks as blocks
from pygromos.utils.typing import Dict, List, Union, Reference_Position_Type, Cnf_Type


class Reference_Position(Cnf):
    """
    This class is a representation of the gromos .cnf coordinate files. It
    allows reading, analysis and modifying of the coordinate files.

    is a child of general_gromos_file
    """

    _future_file: bool
    path: str

    # general
    _gromos_file_ending: str = "rpf"
    residues: Dict[str, Dict[str, int]]

    # Standard Gromos blocks
    TITLE: blocks.TITLE  # required
    REFPOSITION: blocks.REFPOSITION  # required
    LATTICESHIFTS: blocks.LATTICESHIFTS
    GENBOX: blocks.GENBOX
    # private
    _block_order: List[str] = ["TITLE", "REFPOSITION", "LATTICESHIFTS", "GENBOX"]
    _required_blocks: List[str] = ["TITLE", "REFPOSITION"]
    _main_block: str = "REFPOSITION"

    def __init__(
        self,
        in_value: Union[str, dict, Reference_Position_Type, Cnf_Type],
        verbose: bool = False,
        _future_file: bool = False,
    ):
        if isinstance(in_value, Cnf):
            for block in self._block_order:
                if hasattr(in_value, block):
                    setattr(self, block, getattr(in_value, block))
            self.REFPOSITION = blocks.REFPOSITION(in_value.POSITION.content)
            self.path = None
            self._future_file = False

        else:
            super().__init__(in_value=in_value, _future_file=_future_file)
