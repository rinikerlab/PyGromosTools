"""
File:
Description:

Author:
"""

from pygromos.files._basics import _general_gromos_file, parser
from pygromos.utils.typing import Union, Dict


class Ifp(_general_gromos_file._general_gromos_file):
    _gromos_file_ending = "ifp"

    def __init__(self, in_value: Union[str, Dict]):
        super().__init__(in_value=in_value)

    def read_file(self):
        # Read blocks to string
        data = parser.read_general_gromos_file(self._orig_file_path)

        # translate the string subblocks
        blocks = {}
        for block_title in data:
            # print(block_title)
            # print("\t", data[block_title])
            self.add_block(blocktitle=block_title, content=data[block_title])
            blocks.update({block_title: self.__getattribute__(block_title)})
        return blocks

        # for key, sub_content in data.items():
        #    try:
        #        self.add_block(blocktitle=key, content=sub_content)
        #    except Exception as err:
        #        #Here new block updates can be caugth
        #        raise IOError("Could not read in imd " + key + " block!\n values: \n\t" + str(sub_content) + "\n\n" + "\n\t".join(err.args))
