import warnings

from pygromos.files._basics import _general_gromos_file, parser
from pygromos.files.blocks import topology_blocks as blocks


class Pertubation_topology(_general_gromos_file._general_gromos_file):

    _block_order = ["TITLE"]
    required_blocks = ["TITLE", ]
    TITLE: blocks.TITLE
    MPERATOM: blocks.MPERTATOM
    PERTATOMPARAM:blocks.PERTATOMPARAM

    gromos_file_ending:str = "ptp"

    def __init__(self, in_value:(str or dict)=None):
        super().__init__(in_value=in_value)


    def read_blocks(self):
        #parse file into dicts
        data = parser.read_ptp(self.path)

        for key in data:
            if("MPERTATOM" == key):
                data[key]["STATEATOMS"] = list(map(lambda x: blocks.atom_eds_pertubation_state(**x), data["MPERTATOM"]["STATEATOMS"]))
                self.add_block(blocktitle=key, content=data[key])

            elif ("PERTATOMPARAM" == key):
                data["PERTATOMPARAM"]["STATEATOMS"] = list(map(lambda x: blocks.atom_lam_pertubation_state(**x), data["PERTATOMPARAM"]["STATEATOMS"]))
                self.add_block(blocktitle=key, content=data[key])
            else:
                self.add_block(blocktitle=key, content=data[key])


class Ptp(Pertubation_topology):
    pass