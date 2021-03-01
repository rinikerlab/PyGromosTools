"""
File:          gromos residue library
    needed if top and pdb resns or atoms are not the same.

Author: Benjamin Ries
"""

#imports
from inspect import istraceback
import warnings
from typing import Dict, List, NamedTuple
import math

import os
from pygromos.utils import bash as bash
from pygromos.files._basics import _general_gromos_file, parser
from pygromos.files.blocks import miscBlocks as blocks, _general_blocks as titleBlock


class residue_library(_general_gromos_file._general_gromos_file):
    required_blocks = ["TITLE", "RESIDUENAMELIB", "ATOMNAMELIB"]
    RESIDUENAMELIB: blocks.RESIDUENAMELIB
    ATOMNAMELIB: blocks.ATOMNAMELIB

    def __init__(self, path:(str or dict)=None):

        self.blocksset = []
        if(type(path) is str):
            self.path = path
            self.read_resnlib(path)
        elif(path==None):
            print("Warning!: generated empty REsidue Lib obj!")
            self.TITLE = titleBlock.TITLE(content="New empyt resn_lib-file\n\tgenerated with PyGromosTools.\n")
            self.RESIDUENAMELIB = blocks.RESIDUENAMELIB({})
            self.ATOMNAMELIB = blocks.ATOMNAMELIB({})


        else:
            raise IOError("pertubation_topology class got "+str(type(path))+" as input. Unknown input type for disres.")

    def read_resnlib(self, path:str):
        data = parser.read_general_gromos_file(path)
        #add _blocks as attribute to objects
        for key, sub_content in data.items():
            self.add_block(blocktitle=key, content=sub_content)
