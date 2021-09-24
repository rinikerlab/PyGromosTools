import warnings

import pygromos.files.blocks.pertubation_blocks
from pygromos.files._basics import _general_gromos_file, parser
from pygromos.files.blocks import pertubation_blocks as blocks


class Pertubation_topology(_general_gromos_file._general_gromos_file):

    _block_order = ["TITLE"]
    required_blocks = ["TITLE", ]
    TITLE: blocks.TITLE
    MPERATOM: blocks.MPERTATOM
    PERTATOMPARAM: blocks.PERTATOMPARAM
    PERTBONDSTRETCH: blocks.PERTBONDSTRETCH
    PERTBONDSTRETCHH: blocks.PERTBONDSTRETCHH
    PERTBONDANGLE: blocks.PERTBONDANGLE
    PERTBONDANGLEH: blocks.PERTBONDANGLEH
    PERTPROPERDIH: blocks.PERTPROPERDIH

    gromos_file_ending:str = "ptp"

    def __init__(self, in_value:(str or dict)=None):
        super().__init__(in_value=in_value)


    def read_blocks(self):
        #parse file into dicts
        data = parser.read_ptp(self.path)

        for key in data:
            self.add_block(block=data[key])



class Ptp(Pertubation_topology):
    pass