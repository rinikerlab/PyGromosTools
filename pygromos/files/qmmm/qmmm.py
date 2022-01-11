import numpy as np
from numbers import Number
import copy, json
from typing import List, Dict, NamedTuple, Iterable
from collections import namedtuple
from copy import deepcopy

from pygromos.files._basics import _general_gromos_file, parser
from pygromos.files.blocks import qmmm_blocks as blocks



class QMMM(_general_gromos_file._general_gromos_file):
    _gromos_file_ending:str = "qmmm"

    _orig_file_path:str
    path:str
    _required_blocks = ["TITLE", "QMZONE", "QMUNIT"] # Not yet implemented; taken from "class Imd", other blocks are e.g. "XTBELEMENTS" and "ORCAELEMENTS"

    # POSSIBLE GROMOS BLOCKS
    TITLE:  blocks.TITLE
    QMZONE: blocks.QMZONE
    QMUNIT: blocks.QMUNIT

    XTBELEMENTS: blocks.XTBELEMENTS = None
    ORCAELEMENTS: blocks.ORCAELEMENTS = None

    def __init__(self, in_value:str, _future_file:bool=False):
        super().__init__(in_value=in_value, _future_file=_future_file)

        #TODO: maybe somebody can make a better solution for this. This is a ugly fix to unify the structure of the blocks
        for block in sorted(self.get_block_names()):
            setattr(self, block, deepcopy(getattr(self, block)))

    def __str__(self):
        text = ""
        if(hasattr(self, "TITLE")):
            text += self.__getattribute__("TITLE").block_to_string()
        for block in sorted(self.get_block_names()):
            if(block == "TITLE" or isinstance(block, type(None))):
                continue
            text += str(self.__getattribute__(block))
        return text

    def read_file(self):
        data = parser.read_imd(self._orig_file_path)
        for key, sub_content in data.items():
            try:
                self.add_block(blocktitle=key, content=sub_content)
            except Exception as err:
                raise Exception("Error while reading file: " + str(err))
        return {}
