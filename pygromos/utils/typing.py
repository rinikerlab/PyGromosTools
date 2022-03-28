from numbers import Number  # noqa: F401
from typing import TypeVar
from typing import Dict, List, Union, Tuple, Iterable, Callable  # noqa: F401

# Files
# Topology
Top_Type = TypeVar("Top_Type")

# CoordinateFiles
Cnf_Type = TypeVar("Cnf_Type")
Reference_Position_Type = TypeVar("Reference_Position_Type")
Position_Restraints_Type = TypeVar("Position_Restraints_Type")

# QM/MM - Files
QMUNIT_Type = TypeVar("QMUNIT_Type")
MNDOELEMENTS_Type = TypeVar("MNDOELEMENTS_Type")
TURBOMOLEELEMENTS_Type = TypeVar("TURBOMOLEELEMENTS_Type")
DFTBELEMENTS_Type = TypeVar("DFTBELEMENTS_Type")
MOPACELEMENTS_Type = TypeVar("MOPACELEMENTS_Type")
ORCAELEMENTS_Type = TypeVar("ORCAELEMENTS_Type")
XTBELEMENTS_Type = TypeVar("XTBELEMENTS_Type")


# Trajs
_General_Trajectory_Type = TypeVar("_General_Trajectory_Type")
Repdat_Type = TypeVar("Repdat_Type")
Tre_Type = TypeVar("Tre_Type")

_topology_table_block_Type = TypeVar("_topology_table_block_Type")
_iterable_topology_block_Type = TypeVar("_iterable_topology_block_Type")
PHYSICALCONSTANTS_Type = TypeVar("PHYSICALCONSTANTS_Type")
TOPVERSION_Type = TypeVar("TOPVERSION_Type")
ATOMTYPENAME_Type = TypeVar("ATOMTYPENAME_Type")
RESNAME_Type = TypeVar("RESNAME_Type")
TOPVERSION_Type = TypeVar("TOPVERSION_Type")
