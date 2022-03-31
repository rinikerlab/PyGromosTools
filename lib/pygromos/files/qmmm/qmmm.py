from copy import deepcopy
import warnings

from pygromos.files._basics import _general_gromos_file, parser
from pygromos.files.blocks import qmmm_blocks as blocks
from pygromos.utils.typing import List


class QMMM(_general_gromos_file._general_gromos_file):
    _gromos_file_ending: str = "qmmm"

    _orig_file_path: str
    path: str
    _required_blocks = [
        "TITLE",
        "QMZONE",
        "QMUNIT",
    ]  # Not yet implemented; taken from "class Imd", other blocks are e.g. "XTBELEMENTS" and "ORCAELEMENTS"

    # POSSIBLE GROMOS BLOCKS
    TITLE: blocks.TITLE
    QMZONE: blocks.QMZONE
    QMUNIT: blocks.QMUNIT

    # One of those blocks should be available
    # Consult the Gromos documentation for more details
    # Also consult the corresponding QMMM block in imd files
    MNDOELEMENTS: blocks.MNDOELEMENTS
    TURBOMOLEELEMENTS: blocks.TURBOMOLEELEMENTS
    DFTBELEMENTS: blocks.DFTBELEMENTS
    MOPACELEMENTS: blocks.MOPACELEMENTS
    ORCAELEMENTS: blocks.ORCAELEMENTS
    XTBELEMENTS: blocks.XTBELEMENTS

    def __init__(self, in_value: str, _future_file: bool = False):
        super().__init__(in_value=in_value, _future_file=_future_file)

        # TODO: maybe somebody can make a better solution for this. This is a ugly fix to unify the structure of the blocks
        for block in sorted(self.get_block_names()):
            setattr(self, block, deepcopy(getattr(self, block)))

        # Perform some sanity checks
        # only if not a future file and if there is at least one block
        # (to avoid being run while deep_copy)
        if not _future_file and len(self.get_block_names()) > 0:
            self._health_check()

    def __str__(self):
        text = ""
        if hasattr(self, "TITLE"):
            text += self.__getattribute__("TITLE").block_to_string()
        for block in sorted(self.get_block_names()):
            if block == "TITLE" or isinstance(block, type(None)):
                continue
            text += str(self.__getattribute__(block))
        return text

    def read_file(self):
        # Read blocks to string
        data = parser.read_general_gromos_file(self._orig_file_path)

        # translate the string subblocks
        blocks = {}
        for block_title in data:
            self.add_block(blocktitle=block_title, content=data[block_title])
            blocks.update({block_title: self.__getattribute__(block_title)})
        return blocks

    def get_qm_engines(self) -> List[str]:
        """
        Returns the QM engine used

        Returns
        -------
        List[str]
            A list of strings (in case the user wanted for some reason specify more than one engine)
        """
        engine = [element.replace("ELEMENTS", "") for element in self.get_block_names() if element.endswith("ELEMENTS")]
        return engine

    def _health_check(self):
        """
        Runs tests on the integrity of the file and spits out warnings
        members = [attr for attr in dir(self) if not callable(getattr(self, attr)) and attr.endswith("ELEMENTS")]
        """
        members = [element for element in self.get_block_names() if element.endswith("ELEMENTS")]
        if len(members) > 1:
            warnings.warn(
                "Declaration of more than one (QMENGINE)ELEMENTS blocks in QM/MM specification file detected: "
                + ", ".join(members)
            )
        elif len(members) < 1:
            warnings.warn("No (QMENGINE)ELEMENTS block in QM/MM specification file detected.")
