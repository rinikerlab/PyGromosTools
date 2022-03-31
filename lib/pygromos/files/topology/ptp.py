from copy import deepcopy

from pygromos.files._basics import _general_gromos_file, parser
from pygromos.files.blocks import pertubation_blocks as blocks
from pygromos.utils.typing import Union


class Pertubation_topology(_general_gromos_file._general_gromos_file):

    _block_order = ["TITLE"]
    required_blocks = [
        "TITLE",
    ]
    TITLE: blocks.TITLE
    MPERATOM: blocks.MPERTATOM
    PERTATOMPARAM: blocks.PERTATOMPARAM
    PERTBONDSTRETCH: blocks.PERTBONDSTRETCH
    PERTBONDSTRETCHH: blocks.PERTBONDSTRETCHH
    PERTBONDANGLE: blocks.PERTBONDANGLE
    PERTBONDANGLEH: blocks.PERTBONDANGLEH
    PERTPROPERDIH: blocks.PERTPROPERDIH

    _gromos_file_ending: str = "ptp"

    def __init__(self, in_value: Union[str, dict] = None):
        super().__init__(in_value=in_value)

        # TODO: maybe somebody can make a better solution for this. This is a ugly fix to unify the structure of the blocks
        for block in sorted(self.get_block_names()):
            setattr(self, block, deepcopy(getattr(self, block)))

    def read_blocks(self):
        # parse file into dicts
        data = parser.read_ptp(self.path)

        for key in data:
            self.add_block(block=data[key])


class Ptp(Pertubation_topology):
    pass
