from enum import Enum
from typing import List, Tuple

from pygromos.gromos.pyGromosPP.generate_FETopologyApproach import generate_dual_topology_approach, generate_hybrid_topology_approach, generate_single_topology_approach


from pygromos.files.coord.cnf import Cnf
from pygromos.files.top.top import Top
from pygromos.files.top.ptp import Ptp
from pygromos.files.gromos_system.gromos_system import Gromos_System

    """
        UTILS
    """
class topology_approach(Enum):
    dual = 1
    hybrid = 2
    single = 3
class AlchemicalFE_System(Gromos_System):
    """ 
        The Alchemical Free Energy Systems, should help to quickly setup and run FE simulations.
    
    Todo: Consider how to insert the FF-Conversions for single topologies (e.g. amber2Gromos)

    """


    def __init__(work_folder: str, system_name: str, in_endstate_cnfs:List[Cnf], in_endstate_tops:List[Top],
                 in_imd_path: str = None, topology_approach:topology_approach=topology_approach.dual, eds_simulation:bool=False,
                 in_gromosXX_bin_dir:str = None, in_gromosPP_bin_dir:str=None,
                 readIn=True, Forcefield:forcefield_system=forcefield_system(),
                 auto_convert:bool=False, adapt_imd_automatically:bool=True, verbose:bool=False):

        self.single_endstate_coordinates = in_endstate_cnfs
        self.single_endstate_topologies = in_endstate_tops

        super().__init__(work_folder=work_folder, system_name=system_name, in_imd_path=in_imd_path, 
                         in_gromosXX_bin_dir=in_gromosXX_bin_dir, in_gromosPP_bin_dir=in_gromosPP_bin_dir,
                         readIn=readIn, Forcefield=Forcefield, auto_convert=auto_convert, adapt_imd_automatically=adapt_imd_automatically
                         verbose=verbose)

        self.eds_simulation = self._is_eds_simulation(eds_simulation)

        self.generate_alchemical_system(in_endstate_cnfs=in_endstate_cnfs, in_endstate_tops=in_endstate_tops,
                                        topology_approach=topology_approach, eds_simulation=eds_simulation)


    def _is_eds_simulation(self, eds_simulation:bool)->bool:
        if(hasattr(self, "imd") and eds_simulation is None):
            if(hasattr(self.imd, "EDS") and self.imd.EDS.EDS):
                eds_simulation = True
            elif(hasattr(self.imd, "REEDS") and self.imd.REEDS.REEDS>0):    #will work in future
                eds_simulation = True
            else:
                eds_simulation = False
        return eds_simulation
    
    @classmethod
    def generate_alchemical_system(cls, in_endstate_cnfs:List[Cnf], in_endstate_tops:List[Top], topology_approach:topology_approach, eds_simulation:bool=False):
        if(topology_approach == case topology_approach.dual):
                self.cnf, self.top, self.ptp, self.disres = generate_dual_topology_approach(endstate_cnfs=self.single_endstate_coordinates,
                                                                                            endstate_tops=self.single_endstate_topologies,
                                                                                            eds_simulation=self.eds_simulation,
                                                                                            generate_distance_restraints=True)
        elif(topology_approach == case topology_approach.hybrid):
                 self.cnf, self.top, self.ptp = generate_hybrid_topology_approach(endstate_cnfs=self.single_endstate_coordinates,
                                                                                            endstate_tops=self.single_endstate_topologies,
                                                                                            atom_mapping=None
                                                                                            eds_simulation=self.eds_simulation,)

        elif(topology_approach == case topology_approach.single):
                 self.cnf, self.top, self.ptp = generate_single_topology_approach(endstate_cnfs=self.single_endstate_coordinates,
                                                                                            endstate_tops=self.single_endstate_topologies,
                                                                                            atom_mapping=None
                                                                                            eds_simulation=self.eds_simulation,)
        else:
            raise IOError("No valid topology_approach was passed. please use:\n"+str(help(topology_approach)))
