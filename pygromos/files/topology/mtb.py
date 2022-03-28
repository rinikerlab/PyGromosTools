"""
File:
Description:

Author: Marc Lehner
"""

from pygromos.files._basics import _general_gromos_file
from pygromos.files.blocks import mtb_blocks as blocks
from pygromos.utils.typing import Dict, List, Union


class Mtb(_general_gromos_file._general_gromos_file):
    _gromos_file_ending: str = "mtb"
    mtb_solutes: Dict[str, blocks.MTBUILDBLSOLUTE]
    mtb_solvents: Dict[str, blocks.MTBUILDBLSOLVENT]
    mtb_ends: Dict[str, blocks.MTBUILDBLEND]
    all_res_names: List

    def __init__(self, in_value: Union[str, Dict], _future_file: bool = False):
        self.mtb_solutes = {}
        self.mtb_solvents = {}
        self.mtb_ends = {}
        self.all_res_names = []
        super().__init__(in_value, _future_file)

    def __str__(self):
        ret_str = super().__str__()
        for res_name in self.mtb_solutes:
            ret_str += str(self.mtb_solutes[res_name])
        for res_name in self.mtb_solvents:
            ret_str += str(self.mtb_solvents[res_name])
        for res_name in self.mtb_ends:
            ret_str += str(self.mtb_ends[res_name])
        return ret_str

    def read_file(self):
        # define some containers
        MTBUILDBLSOLUTE_list = []
        MTBUILDBLSOLVENT_list = []
        MTBUILDBLEND_list = []
        # Read blocks to string
        data = self.read_mtb_file(self._orig_file_path)

        # translate the string subblocks
        block_dict = {}
        for block_title, block_data in data:
            if block_title == "MTBUILDBLSOLUTE":
                mtb_block = blocks.MTBUILDBLSOLUTE(content=block_data)
                MTBUILDBLSOLUTE_list.append(mtb_block)
            elif block_title == "MTBUILDBLSOLVENT":
                mtb_block = blocks.MTBUILDBLSOLVENT(content=block_data)
                MTBUILDBLSOLVENT_list.append(mtb_block)
            elif block_title == "MTBUILDBLEND":
                mtb_block = blocks.MTBUILDBLEND(content=block_data)
                MTBUILDBLEND_list.append(mtb_block)
            else:
                self.add_block(blocktitle=block_title, content=block_data)
                block_dict.update({block_title: self.__getattribute__(block_title)})
        # convert mtb lists to dicts
        self.mtb_solutes = {b.RNME: b for b in MTBUILDBLSOLUTE_list}
        self.mtb_solvents = {b.RNMES: b for b in MTBUILDBLSOLVENT_list}
        self.mtb_ends = {b.RNME: b for b in MTBUILDBLEND_list}
        self.all_res_names = list(self.mtb_solutes.keys()) + list(self.mtb_solvents.keys()) + list(self.mtb_ends.keys())
        return block_dict

    def read_mtb_file(self, path: str) -> List[str]:
        data = []

        first_key = True
        key = ""
        block = []

        with open(path, "r") as infile:
            for line in infile:
                if first_key:
                    if line.startswith("#"):
                        continue
                    key = line.strip().upper()
                    first_key = False
                elif "END" in line:
                    data.append([key, block])
                    block = []
                    first_key = True
                else:
                    block.append(line)
        infile.close()
        return data
