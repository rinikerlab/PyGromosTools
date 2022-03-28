from pygromos.files.coord.cnf import Cnf
from pygromos.files.blocks import coord_blocks as blocks
from pygromos.utils.typing import Dict, List, Union, Position_Restraints_Type, Cnf_Type


class Position_Restraints(Cnf):
    """
    This class is a representation of the gromos .cnf coordinate files. It
    allows reading, analysis and modifying of the coordinate files.

    is a child of general_gromos_file
    """

    _future_file: bool
    path: str

    # general
    _gromos_file_ending: str = "pos"
    residues: Dict[str, Dict[str, int]]

    # Standard Gromos blocks
    TITLE: blocks.TITLE  # required
    POSRESSPEC: blocks.POSRESSPEC  # required
    LATTICESHIFTS: blocks.LATTICESHIFTS

    # private
    _block_order: List[str] = [
        "TITLE",
        "POSRESSPEC",
    ]
    _required_blocks: List[str] = ["TITLE", "POSRESSPEC"]
    _main_block: str = "POSRESSPEC"

    def __init__(
        self,
        in_value: Union[str, dict, Position_Restraints_Type, Cnf_Type],
        clean_resiNumbers_by_Name: bool = False,
        verbose: bool = False,
        _future_file: bool = False,
    ):
        if isinstance(in_value, Cnf):
            for block in self._block_order:
                if hasattr(in_value, block):
                    setattr(self, block, getattr(in_value, block))
            self.POSRESSPEC = blocks.POSRESSPEC(in_value.POSITION.block_to_string().split("\n")[1:-2])
            self.path = None
            self._future_file = False
        else:
            super().__init__(in_value=in_value, _future_file=_future_file)
