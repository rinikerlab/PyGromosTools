"""_summary_

"""
    
from typing import List, Tuple

from pygromos.files.coord.cnf import Cnf
from pygromos.files.top.top import Top
from pygromos.files.top.ptp import Ptp
from pygromos.files.top.disres import Disres



def generate_dual_topology_approach(endstate_cnfs:List[Cnf], endstate_tops:List[Top], eds_simulation:bool=False, generate_distance_restraints:bool=True)->Tuple[Cnf, Top, Ptp, Disres]:
    """
        This function generates a Gromos System, according to the dual topology paradigm.
        In this process a coordinate, topology and pertubation file will be generated.
        Additionally a distance restraint file will be generated with RestraintMaker.
        
        The definition of the dual topology follows the description in:
            RestraintMaker: A Graph-Based Approach to Select Distance Restraints in Free-Energy Calculations with Dual Topology
            Benjamin Ries$^\text{\dag}$, Salomé Rieder$^\text{\dag}$, Clemens Rhiner, Philippe H. Hünenberger$^\text{*}$ and Sereina Riniker$^\text{*}$ (2022, JCAMD)

    Parameters
    ----------
    endstate_cnfs : List[Cnf]
        The cnf files of the single end-states (for protein-ligands FE -> [ligandA.cnf, ligandB.cnf, ...])
    endstate_tops : List[Top]
        The top files of the single end-states (for protein-ligands FE -> [ligandA.top, ligandB.top, ...])
    eds_simulation : bool, optional
        If you want do perform a eds simulation set True, as the ptp blocks vary between eds and lambda dependent calculations, by default False
    generate_distance_restraints : bool, optional
        If True, the dual topology approach will additionally contain distance restrained following the linked dual topology principal. 
        These restraints will be automatically attempted by RestraintMaker. (But double check!)
    
    Returns
    -------
    Tuple[Cnf, Top, Ptp, Disres]
        Returns the prepared Alchemical System files.
    """
    pass

def generate_hybrid_topology_approach(endstate_cnfs:List[Cnf], endstate_tops:List[Cnf], atom_mapping:List=None, eds_simulation:bool=False)->Tuple[Cnf, Top, Ptp]:
    """
        This function generates a Gromos System, according to the hybrid topology paradigm.
        In this process a coordinate, topology and pertubation file will be generated.
        A atomMapping for the different regions (single topology region/dual topology region) can be provided.
        If no mapping is provided, the function will automatically generate an MCS and estimate the single and dual topology regions.
        
        The definition of the hybrid topology follows the description in:
            RestraintMaker: A Graph-Based Approach to Select Distance Restraints in Free-Energy Calculations with Dual Topology
            Benjamin Ries$^\text{\dag}$, Salomé Rieder$^\text{\dag}$, Clemens Rhiner, Philippe H. Hünenberger$^\text{*}$ and Sereina Riniker$^\text{*}$ (2022, JCAMD)

    Parameters
    ----------
    endstate_cnfs : List[Cnf]
        The cnf files of the single end-states (for protein-ligands FE -> [ligandA.cnf, ligandB.cnf, ...])
    endstate_tops : List[Top]
        The top files of the single end-states (for protein-ligands FE -> [ligandA.top, ligandB.top, ...])
    atom_mapping : List, optional
        The mapping of the atoms for the different end-states, by default None
    eds_simulation : bool, optional
        If you want do perform a eds simulation set True, as the ptp blocks vary between eds and lambda dependent calculations, by default False

    Returns
    -------
    Tuple[Cnf, Top, Ptp]
        Returns the prepared Alchemical System files.
    """
    pass

def generate_single_topology_approach(endstate_cnfs:List[Cnf], endstate_tops:List[Top], atom_mapping:List=None, eds_simulation:bool=False)->Tuple[Cnf, Top, Ptp]:
    """
        This function generates a Gromos System, according to the single topology paradigm.
        In this process a coordinate, topology and pertubation file will be generated.
        An atomMapping for the single-topology regio can be provided. 
        If no mapping is provided, the function will automatically generate an MCS and estimate the molecule overlaps.

        The definition of the single topology follows the description in:
            RestraintMaker: A Graph-Based Approach to Select Distance Restraints in Free-Energy Calculations with Dual Topology
            Benjamin Ries$^\text{\dag}$, Salomé Rieder$^\text{\dag}$, Clemens Rhiner, Philippe H. Hünenberger$^\text{*}$ and Sereina Riniker$^\text{*}$ (2022, JCAMD)

    Parameters
    ----------
    endstate_cnfs : List[Cnf]
        The cnf files of the single end-states (for protein-ligands FE -> [ligandA.cnf, ligandB.cnf, ...])
    endstate_tops : List[Top]
        The top files of the single end-states (for protein-ligands FE -> [ligandA.top, ligandB.top, ...])
    atom_mapping : List, optional
        The mapping of the atoms for the different end-states, by default None
    eds_simulation : bool, optional
        If you want do perform a eds simulation set True, as the ptp blocks vary between eds and lambda dependent calculations, by default False


    Returns
    -------
    Tuple[Cnf, Top, Ptp]
        Returns the prepared Alchemical System files.
    """
    pass