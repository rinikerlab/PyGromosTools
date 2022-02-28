"""
File:
Description:

Author: Marc Lehner
"""

from typing import Dict, List
import warnings
from pygromos.files._basics import _general_gromos_file
from pygromos.files.blocks import mtb_blocks as blocks


class Mtb(_general_gromos_file._general_gromos_file):
    _gromos_file_ending: str = "mtb"
    MTBUILDBLSOLUTE_list: List[blocks.MTBUILDBLSOLUTE]
    MTBUILDBLSOLVENT_list: List[blocks.MTBUILDBLSOLVENT]

    def __init__(self, in_value: (str or dict or None), _future_file: bool = False):
        self.MTBUILDBLSOLUTE_list = []
        self.MTBUILDBLSOLVENT_list = []
        super().__init__(in_value, _future_file)

    def __str__(self):
        ret_str = super().__str__()
        for block in self.MTBUILDBLSOLUTE_list:
            ret_str += str(block)
        for block in self.MTBUILDBLSOLVENT_list:
            ret_str += str(block)
        return ret_str

    def read_file(self):
        # Read blocks to string
        data = self.read_mtb_file(self._orig_file_path)

        # translate the string subblocks
        block_dict = {}
        for block_title, block_data in data:
            if block_title == "MTBUILDBLSOLUTE":
                mtb_block = blocks.MTBUILDBLSOLUTE(content=block_data)
                self.MTBUILDBLSOLUTE_list.append(mtb_block)
            elif block_title == "MTBUILDBLSOLVENT":
                mtb_block = blocks.MTBUILDBLSOLVENT(content=block_data)
                self.MTBUILDBLSOLVENT_list.append(mtb_block)
            else:
                self.add_block(blocktitle=block_title, content=block_data)
                block_dict.update({block_title: self.__getattribute__(block_title)})
        return block_dict

    def read_mtb_file(self, path: str) -> List:
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
