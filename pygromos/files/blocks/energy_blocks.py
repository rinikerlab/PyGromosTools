
from pygromos.files.blocks._general_blocks import TITLE, TIMESTEP, TRAJ

# forward declarations
TITLE:TITLE = TITLE
TIMESTEP:TIMESTEP = TIMESTEP
TRAJ:TRAJ = TRAJ

"""
This was an Idea!


class ENEVERSION(_generic_gromos_block):
    ""ene_ana_version
        This class gives the version of the in_ene_ana_lib
        It is usually contained in the in_ene_ana_lib used by ene_ana.
        There this Block is called ENVERSION
    ""

    def __init__(self, content: str):
        super().__init__(used=True, name="ENEVERSION")
        self.content = content

class ENEANAVARS(_generic_gromos_block):
    ""ene_ana_variables
        This class gives the variable definitions for ene_ana analysing a gromos Energy traj .out_tre or Free Energy file .trf.
        It is usually contained in the in_ene_ana_lib used by ene_ana.
        There this Block is called VARIABLES
    ""

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
    
"""