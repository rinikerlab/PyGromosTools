from typing import Dict

from pygromos.files.blocks._general_blocks import _generic_gromos_block
from pygromos.files.blocks._general_blocks import TITLE, TIMESTEP, TRAJ

# forward declarations
TITLE:TITLE = TITLE
TIMESTEP:TIMESTEP = TIMESTEP
TRAJ:TRAJ = TRAJ

class ENEVERSION(_generic_gromos_block):
    """ene_ana_version
        This class gives the version of the in_ene_ana_lib
        It is usually contained in the in_ene_ana_lib used by ene_ana.
        There this Block is called ENVERSION
    """

    def __init__(self, content: str):
        super().__init__(used=True, name="ENEVERSION")
        self.content = content


class ENERTRJ(_generic_gromos_block):
    """ene_traj_structure
        This class gives the structure definition of an gromos Energy traj .out_tre file.
        It is usually contained in the in_ene_ana_lib used by ene_ana.
        There this Block is called ENETRJ

    """

    def __init__(self, comments: str, blocks: Dict):
        super().__init__(used=True, name="ENERTRJ")
        self.comments = comments
        self.blocks = blocks

    def block_to_string(self) -> str:
        str_blocks = ""
        for block_tit in self.blocks:
            # main blcok
            str_blocks += "block " + block_tit + self.line_seperator
            # resolve subblocks
            sub_blocks = self.blocks[block_tit]
            sub_block_str = ""
            if ("size" in sub_blocks):
                size_var = [self.field_seperator + "size " + size for size in sub_blocks["size"]]
                sub_block_str += self.line_seperator.join(size_var)
            if ("subblock" in sub_blocks):
                sub_blocks_str = [self.field_seperator + "subblock " + self.field_seperator.join(sb) for sb in sub_blocks["subblock"]]
                sub_block_str += self.line_seperator.join(sub_blocks_str)
            str_blocks += sub_block_str + self.line_seperator

        result = self.name + self.line_seperator
        result += self.comments + self.line_seperator
        result += str_blocks + self.line_seperator
        result += "END" + self.line_seperator
        return result

    def parse_block_string(self, block_str):
        """
                    Parse the block input from str.

        Parameters
        ----------
        block_str :

        Returns
        -------

        """

        comments = []
        block_dict = {}
        current_block = False
        lines = block_str.split("\n")
        for line in lines:
            if (self.name in line or "END" in line or "" in line):
                continue
            elif (line.strip().startswith("#")):
                comments.append(line)
            elif (line.strip().startswith("block")):
                current_block = line.replace("block", "").strip()
                block_dict.update({current_block: {}})
            elif (line.strip().startswith("subblock")):
                if (not "subblock" in block_dict[current_block]):
                    block_dict[current_block].update({"subblock": [line.replace("subblock", "").strip().split()]})
                else:
                    block_dict[current_block]["subblock"].append(line.replace("subblock", "").strip().split())
            elif (line.strip().startswith("size")):
                if (not "size" in block_dict[current_block]):
                    block_dict[current_block].update({"size": [line.replace("size", "").strip()]})
                else:
                    block_dict[current_block]["size"].append(line.replace("size", "").strip())
            else:
                raise ValueError("Ene_ana did not understand line in " + current_block + ": ", line)

        self.comments = comments
        self.blocks = block_dict


class FRENERTRJ(_generic_gromos_block):
    """ene_ana_free_traj
        This class gives the structure definition of an free Energy .trf file.
        It is usually contained in the in_ene_ana_lib used by ene_ana.
        There this Block is called FRENETRJ

    """

    def __init__(self, comments: str, blocks: Dict):
        super().__init__(used=True, name="FRENERTRJ")
        self.comments = comments
        self.blocks = blocks

    def block_to_string(self) -> str:
        str_blocks = ""
        for block_tit in self.blocks:
            # main blcok
            str_blocks += "block " + block_tit + self.line_seperator
            # resolve subblocks
            sub_blocks = self.blocks[block_tit]
            sub_block_str = ""
            if ("size" in sub_blocks):
                size_var = [self.field_seperator + "size " + size for size in sub_blocks["size"]]
                sub_block_str += self.line_seperator.join(size_var)
            if ("subblock" in sub_blocks):
                sub_blocks_str = [self.field_seperator + "subblock " + self.field_seperator.join(sb) for sb in sub_blocks["subblock"]]
                sub_block_str += self.line_seperator.join(sub_blocks_str)
            str_blocks += sub_block_str + self.line_seperator

        result = self.name + self.line_seperator
        result += self.comments + self.line_seperator
        result += str_blocks + self.line_seperator
        result += "END" + self.line_seperator
        return result


class ENEANAVARS(_generic_gromos_block):
    """ene_ana_variables
        This class gives the variable definitions for ene_ana analysing a gromos Energy traj .out_tre or Free Energy file .trf.
        It is usually contained in the in_ene_ana_lib used by ene_ana.
        There this Block is called VARIABLES
    """

    def __init__(self, comments: str, variables: Dict):
        super().__init__(used=True, name="ENEANAVARS")
        self.comments = comments
        self.variables = variables

    def block_to_string(self) -> str:
        str_blocks = ""
        for variable_name, variable_value in self.variables.items():
            str_blocks += variable_name + " = " + variable_value + self.line_seperator
        result = self.name + self.line_seperator
        result += self.comments + self.line_seperator
        result += str_blocks
        result += "END" + self.line_seperator
        return result

    pass