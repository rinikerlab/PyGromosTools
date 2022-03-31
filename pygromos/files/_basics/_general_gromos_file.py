"""
FUNCTIONLIB:            gromos baseclasses for files
Description:
    this file is giving base_classes
Author: Benjamin Schroeder
"""
import os
import copy
import inspect
import warnings

from pygromos.files.blocks import all_blocks
from pygromos.utils.typing import List, Dict, Callable, Union


# file class
class _general_gromos_file:
    """_general_gromos_file
    This class is the generic gromos class, that takes care of all common gromos file featuers.
    """

    _orig_file_path: str
    path: str
    _required_blocks = ["TITLE"]
    _blocksset_names: List[str]

    _gromos_file_ending: str
    # private
    _blocks: Dict[str, all_blocks._generic_gromos_block]
    _block_order: List[str] = []
    _future_file: bool

    def __init__(self, in_value: Union[str, dict, None], _future_file: bool = False):
        self._future_file = _future_file
        if isinstance(in_value, str):
            self.path = self._orig_file_path = in_value
            if self._future_file:
                pass
            elif os.path.exists(in_value):
                self.read_blocks()
            else:
                raise IOError("Could not find File: ", in_value)

        elif isinstance(type(in_value), __class__):
            raise NotImplementedError("This variant is not implemented")

        elif in_value is None:
            self.path = None
            self._orig_file_path = None
            # print("Empty class")

        else:
            raise ValueError("The given type of input could not be translated in " + str(__class__) + ".__init__")

    def read_blocks(self):
        self._blocks = self.read_file()
        self._blocksset_names = list(self._blocks.keys())
        [setattr(self, key, value) for key, value in self._blocks.items()]
        del self._blocks

    def __str__(self):
        # first write out certain _blocks
        out_text = ""
        for block in self._block_order:
            if block in self.get_block_names() and not isinstance(
                getattr(
                    self,
                    block,
                ),
                type(None),
            ):
                out_text += getattr(self, block).block_to_string()

        # write out rest of _blocks
        rest_blocks = [
            block
            for block in self.get_block_names()
            if (
                block not in self._block_order
                and not isinstance(
                    getattr(
                        self,
                        block,
                    ),
                    type(None),
                )
            )
        ]
        for block in rest_blocks:
            if block is None or getattr(self, block) is None:
                continue
            else:
                out_text += getattr(self, block).block_to_string()

        return out_text

    def __repr__(self):
        return str(self)

    def __getstate__(self):
        """
        preperation for pickling:
        remove the non trivial pickling parts
        """
        skip = ["_view"]
        attribute_dict = self.__dict__
        new_dict = {}
        for key in attribute_dict.keys():
            if not isinstance(attribute_dict[key], Callable) and key not in skip:
                new_dict.update({key: attribute_dict[key]})

        return new_dict

    def __setstate__(self, state):
        self.__dict__ = state

    def __deepcopy__(self, memo):
        # print(self.__class__.__name__)
        copy_obj = self.__class__(in_value=None)
        copy_obj.__setstate__(copy.deepcopy(self.__getstate__()))
        return copy_obj

    def __eq__(self, __o: object) -> bool:
        if self is not None and __o is not None:
            self_blocks = self.get_block_names()
            other_blocks = __o.get_block_names()
            if self_blocks != other_blocks:
                return False
            else:
                for block in self_blocks:
                    if not getattr(self, block) == getattr(__o, block):
                        return False
                return True
        elif self is None and __o is None:
            return True
        else:
            return False

    def get_block_names(self) -> List[str]:
        """get_block_names
                    This is a function, that returns all _blocks, contained in an imd file.

        Returns
        -------
        List[str]
            returns a list of all str block names contained in the gromos file. (they are also the attribute names of the obj)
        """

        # also a nice idea:  self._block_names = list(filter(lambda x: not x.startswith("_") and x.isupper(), vars(imd).keys()))list(self._blocks.keys())

        return list(filter(lambda x: not x.startswith("_") and x.isupper(), vars(self).keys()))

    def add_block(
        self,
        blocktitle: str = None,
        content: dict = None,
        block: all_blocks._generic_gromos_block = None,
        verbose: bool = False,
    ):
        """add_block
            This method adds a gromos block to a gromos file object.

        Parameters
        ----------
        blocktitle :    str, optional
            title of a block
        content :   str, optional
            block content
        block : all_blocks._generic_gromos_block, optional
            block class
        verbose :   bool, optional
            shall messages be printed?

        Returns
        -------
        None

        """
        if block and not (blocktitle or content):
            if verbose:
                print("import block from " + block.name)
            blocktitle = block.name
            setattr(self, blocktitle, block)  # modern way

        elif blocktitle is not None and content is not None:
            # if blocktitle in self._block_names:
            if blocktitle == "TITLE":
                self.__setattr__(blocktitle, all_blocks.__getattribute__(blocktitle)(content))
            elif isinstance(content, dict):
                try:
                    content = {k.split("(")[0]: v for k, v in content.items()}
                    content = {k.split(":")[0]: v for k, v in content.items()}
                    content = {k.split(" ")[0]: v for k, v in content.items()}
                    # Required for add_block
                    # For nasty block seperation, as someone did not care about unique block names.... I'm looking at you vienna!
                    from pygromos.files.simulation_parameters.imd import Imd
                    from pygromos.files.topology.top import Top
                    from pygromos.files.blocks import imd_blocks, topology_blocks

                    if issubclass(self.__class__, Imd):
                        self.kwCreateBlock(blocktitle, content, imd_blocks)
                    elif issubclass(self.__class__, Top):
                        self.kwCreateBlock(blocktitle, content, topology_blocks)
                    else:
                        self.kwCreateBlock(blocktitle, content, all_blocks.all_blocks)

                    if verbose:
                        print("++++++++++++++++++++++++++++++")
                        print("New Block: Adding " + blocktitle + " block")
                        print(content)
                except Exception as e:
                    msg = (
                        "Error while adding new value - can not resolve value names in '"
                        + blocktitle
                        + "' block!\n"
                        + str(e)
                    )
                    msg += "Content is " + str(tuple(content.keys())) + "\n"
                    msg += (
                        "Block knows "
                        + str((all_blocks.get_all_blocks()[blocktitle].__init__.__code__.co_varnames)[1:])
                        + "\n"
                    )
                    raise IOError(msg)
                if verbose:
                    print("Block " + blocktitle + " added to gromos File object.")

            elif isinstance(content, list):
                block_class = all_blocks.__getattribute__(blocktitle)
                block = block_class(content)

                self.__setattr__(blocktitle, block)

                if verbose:
                    print("++++++++++++++++++++++++++++++")
                    print("New Block: Adding " + blocktitle + " block")
                    print(content)

            else:
                raise Exception("Not implemented")

    def delete_block(self, blockName: str):
        """
        delete_block
            This method deletes a block from a gromos file object.

        Parameters
        ----------
        blockName :    str
            title of a block
        """
        if blockName in self._blocksset_names:
            self._blocksset_names.remove(blockName)
        if hasattr(self, blockName):
            delattr(self, blockName)

    def kwCreateBlock(self, blocktitle, content, blocks_lib):
        block_type = blocks_lib.__getattribute__(blocktitle)  # get the blocktype
        sig = inspect.signature(block_type.__init__)  # block init signature
        known_params = list(sig.parameters.keys())  # the params the function knows
        known_content = {k: v for k, v in content.items() if (k in known_params)}
        unknown_content = {k: v for k, v in content.items() if (k not in known_params)}

        # construct class
        self.__setattr__(blocktitle, block_type(**known_content))

        if len(unknown_content) > 0:  # add unkown parameters
            warnings.warn(
                "FOUND UNKOWN Parameters for "
                + str(block_type.__name__)
                + " adding these params: "
                + str(unknown_content)
            )
            [self.__getattribute__(blocktitle).__setattr__(k, v) for k, v in unknown_content.items()]

    def read_file(self) -> Dict[str, any]:
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

    def write(self, out_path: str) -> str:
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
        self._write_to_file(out_path=out_path, content_str=str(self))
        self.path = os.path.abspath(out_path)
        self._future_file = False

        return out_path

    def _write_to_file(self, out_path: str, content_str: str) -> str:
        """
            write to file
        Returns
        -------

        """
        # 1) OpenFile
        if isinstance(out_path, str):
            if os.path.exists(os.path.dirname(out_path)):
                out_file = open(out_path, "w")
            else:
                raise IOError("Could not find directory to write to: " + str(os.path.dirname(out_path)))
        else:
            raise ValueError("Did not understand the Value of out_path. Must be str.")

        # 3) Write File
        out_file.write(content_str)
        out_file.close()
        return out_path
