"""
File:          gromos residue library
    needed if top and pdb resns or atoms are not the same.

Author: Benjamin Ries
"""

#imports
from typing import Union
from pygromos.files.blocks import _all_blocks as blocks
from pygromos.files._basics import parser

from pygromos.files._basics import _general_gromos_file

from pygromos.data import pdb_lib

class residue_library(_general_gromos_file._general_gromos_file):
    required_blocks = ["TITLE", "RESIDUENAMELIB", "ATOMNAMELIB"]
    RESIDUENAMELIB: blocks.RESIDUENAMELIB
    ATOMNAMELIB: blocks.ATOMNAMELIB
    verbose : bool = False

    _gromos_file_ending = "res"

    def __init__(self, in_value:Union[str, dict]=pdb_lib):
        """
            This class represents a file that is used for the gromosPP program - pdb2g96
            it contains two blocks for residue naming and atom naming

        Parameters
        ----------
        in_value : Union[str, dict]
        """

        self.blocksset = []
        if(type(in_value) is str):
            self.path = in_value
            self.read_resnlib(in_value)
        elif(in_value==None):
            if(self.verbose): print("Warning!: generated empty REsidue Lib obj!")
            self.TITLE = blocks.TITLE(content="New empyt resn_lib-file\n\tgenerated with PyGromosTools.\n")
            self.RESIDUENAMELIB = blocks.RESIDUENAMELIB({})
            self.ATOMNAMELIB = blocks.ATOMNAMELIB({})
        else:
            raise IOError("pertubation_topology class got "+str(type(in_value))+" as input. Unknown input type for disres.")

    def read_resnlib(self, path:str):
        data = parser.read_general_gromos_file(path)
        #add _blocks as attribute to objects
        for key, sub_content in data.items():
            self.add_block(blocktitle=key, content=sub_content)
