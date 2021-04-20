    """
    Python implementation of the Gromos++ program ran_box which is used to generate randomized configurations for liquids (and gases)


    Author: Marc Lehner
    """

from pygromos.files.topology.top import Top
from pygromos.files.coord.cnf import Cnf

def ran_box(self, in_top_path:str, 
                    in_cnf_path:str, out_cnf_path:str= "", 
                    periodic_boundary_condition: str = "r", 
                    nmolecule:int = 1, 
                    dens:float = 1.0, 
                    threshold:float=None, 
                    layer:bool = False, 
                    boxsize:float=None, 
                    fixfirst:bool = False, 
                    seed:float=None, 
                    _binary_name="ran_box", 
                    verbose=False, 
                    return_command_only=False)->str:

    top = Top(in_value=in_top_path)
    cnf = Cnf(in_value=in_cnf_path)

    mol_mass = top.get_mass()
    return ""