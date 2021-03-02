"""
FUNCTIONLIB:            gromos baseclasses for files
Description:
    this file is giving base_classes
Author: Benjamin Schroeder
"""
import os
from typing import List, Dict
from pygromos.files.blocks import _all_blocks as blocks


##file class
class _general_gromos_file():
    """_general_gromos_file
        This class is the generic gromos class, that takes care of all common gromos file featuers.
    """
    _orig_file_path:str
    path:str
    _required_blocks = ["TITLE"]
    _blocksset_names:List[str]

    #private
    _blocks:Dict[str, blocks._generic_gromos_block]
    _block_order: List[str] = []

    def __init__(self, input:(str or dict or None or __class__)):
        if(isinstance(input, str)):
            if (not os.path.exists(input)):
                raise IOError("Could not find File: ", input)

            self.file_path = self._orig_file_path = input
            self.read_blocks()
        elif(isinstance(type(input), __class__)):
            raise NotImplementedError("This variant is not implemented")
        elif(input is None):
            self.file_path = None
            self._orig_file_path = None
            #print("Empty class")
        else:
            raise ValueError("The given type of input could not be translated in "+str(__class__)+".__init__")
    
    def read_blocks(self):
        self._blocks = self.read_file()
        self._blocksset_names = list(self._blocks.keys())
        [setattr(self, key, value) for key, value in self._blocks.items()]
        del self._blocks

    def __str__(self):
        ##first write out certain _blocks
        out_text=""
        for block in self._block_order:
            if(block in self.get_block_names() and not isinstance(getattr(self, block,), type(None))):
                out_text += getattr(self, block).block_to_string()

        ##write out rest of _blocks
        rest_blocks = [block for block in self.get_block_names() if (not block in self._block_order and not isinstance(getattr(self, block,), type(None)))]
        for block in rest_blocks:
            if(block is None or getattr(self, block) is None):
                continue
            else:
                out_text += getattr(self, block).block_to_string()

        return out_text

    def get_block_names(self)->List[str]:
        """get_block_names
                    This is a function, that returns all _blocks, contained in an imd file.

        Returns
        -------
        List[str]
            returns a list of all str block names contained in the gromos file. (they are also the attribute names of the obj)
        """

        # also a nice idea:  self._block_names = list(filter(lambda x: not x.startswith("_") and x.isupper(), vars(imd).keys()))list(self._blocks.keys())

        return list(filter(lambda x: not x.startswith("_") and x.isupper(), vars(self).keys()))

    def add_block(self, blocktitle:str=None, content:dict=None, block:blocks._generic_gromos_block=None, verbose:bool=False):
        """add_block
            This method adds a gromos block to a gromos file object.

        Parameters
        ----------
        blocktitle :    str, optional
            title of a block
        content :   str, optional
            block content
        block : blocks._generic_gromos_block, optional
            block class
        verbose :   bool, optional
            shall messages be printed?

        Returns
        -------
        None

        """
        #"todo: rework"
        if(block and not (blocktitle or content)):
            if verbose:
                print("import block from " + block.name)
            blocktitle = block.name
            content = {}
            setattr(self, blocktitle, block)  #modern way
            

        elif(blocktitle != None and content != None):
            #if blocktitle in self._block_names:
            if(isinstance(content, dict)):
                if blocktitle == "TITLE":   #TODO fIX IN PARSER
                    self.__setattr__(blocktitle, blocks.__getattribute__(blocktitle)(content))
                else:
                    try:
                        content = {k.split("(")[0]: v for k, v in content.items()}
                        content = {k.split(":")[0]: v for k, v in content.items()}
                        content = {k.split(" ")[0]: v for k, v in content.items()}
                        # Required for add_block
                        ##For nasty block seperation, as someone did not care about unique block names.... I'm looking at you vienna!
                        from pygromos.files.imd import Imd
                        from pygromos.files.topology.top import Top
                        from pygromos.files.blocks import imd_blocks, topology_blocks

                        if(issubclass(self.__class__, Imd)):
                            self.__setattr__(blocktitle, imd_blocks.__getattribute__(blocktitle)(**content))
                        elif(issubclass(self.__class__, Top)):
                            self.__setattr__(blocktitle, topology_blocks.__getattribute__(blocktitle)(**content))
                        else:
                            self.__setattr__(blocktitle, blocks.__getattribute__(blocktitle)(**content))
                    
                        if verbose:
                            print("++++++++++++++++++++++++++++++")
                            print("New Block: Adding " + blocktitle + " block")
                            print(content)
                    except:
                        msg ="Error while adding new value - can not resolve value names in \'" + blocktitle + "\' block!\n"
                        msg += "Content is " + str(tuple(content.keys()))+"\n"
                        msg += "Block knows " + str((blocks.__getattribute__(blocktitle).__init__.__code__.co_varnames)[1:]) + "\n"
                        raise IOError(msg)
                if verbose:
                    print("Block "+blocktitle+" added to gromos File object.")

            elif(isinstance(content, list)):
                block_class = blocks.__getattribute__(blocktitle)
                block = block_class(content)

                self.__setattr__(blocktitle, block)

                if verbose:
                    print("++++++++++++++++++++++++++++++")
                    print("New Block: Adding " + blocktitle + " block")
                    print(content)   

            else:
                raise Exception("Not implemented")


    def read_file(self)->Dict[str, any]:
        """read_file

            give back the content. WARNING DEAPRECEATED.
        Warnings
        --------
            DEAPRECEATED

        Returns
        -------
        Dict[str, any]
            key is the block name of the gromos file, any is the content of a block
        """
        print("Begin read in of file: " + self._orig_file_path)
        raise NotImplementedError("The read in parser function for general gromos files is not implemented yet!")

    def write(self, out_path:str)->str:
        """write
            writes a general Gromos File out

        Parameters
        ----------
        out_path :  str
            out path where the file should be.

        Returns
        -------
        str
            out_path
        """
        file = open(out_path, "w")
        file.write(str(self))
        file.close()
        self.path = out_path

        return out_path


class _general_gromos_trajectory_file():

    def __init__(self, input:str, every_step:int=1):
        if isinstance(input, str):
            if (not os.path.exists(input)):
                raise IOError("Could not find File: ", input)

            self._orig_file_path = input

            self._blocks = self.read_file(every_step)
            self._blocksset_names = list(self._blocks.keys())
            [setattr(self, key, value) for key, value in self._blocks.items()]
            del self._blocks
        else:
            raise ValueError("The given type of input could not be translated in "+str(__class__)+".__init__")

        def read_file(self, every_step: int = 1) -> Dict[str, any]:
            """read_file

                give back the content. WARNING DEAPRECEATED.
            Warnings
            --------
                DEAPRECEATED

            Returns
            -------
            Dict[str, any]
                key is the block name of the gromos file, any is the content of a block
            """
            print("Begin read in of file: " + self._orig_file_path)
            raise NotImplementedError("The read in parser function for general gromos Files is not implemented yet!")
