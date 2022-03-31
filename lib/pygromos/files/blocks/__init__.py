import inspect

# from pygromos.files.blocks import _general_block
from pygromos.files.blocks import coord_blocks as coords
from pygromos.files.blocks import imd_blocks as imd
from pygromos.files.blocks import topology_blocks as topo
from pygromos.files.blocks import mtb_blocks as mtb
from pygromos.files.blocks import pertubation_blocks as ptp
from pygromos.files.blocks import qmmm_blocks as qmm
from pygromos.files.blocks import replica_exchange_blocks as repdat
from pygromos.files.blocks import miscBlocks as misc

from pygromos.files.blocks._general_blocks import _generic_gromos_block as _ggb
from pygromos.files.blocks._general_blocks import _iterable_gromos_block as _igb
from pygromos.files.blocks._general_blocks import TIMESTEP as generic_TIMESTEP
from pygromos.files.blocks._general_blocks import TITLE as generic_TITLE
from pygromos.files.blocks._general_blocks import TRAJ as generic_TRAJ

# forward declarations
TITLE: generic_TITLE = generic_TITLE
TIMESTEP: generic_TIMESTEP = generic_TIMESTEP
TRAJ: generic_TRAJ = generic_TRAJ
_generic_gromos_block: _ggb = _ggb
_iterable_gromos_block: _igb = _igb


class all_blocks_class:
    """
    This class is a very simplistic block library containing all present gromos blocks.
    """

    def __init__(self):
        self.all_block_types = [coords, imd, topo, mtb, ptp, qmm, repdat, misc]

        for block_types in self.all_block_types:
            [
                setattr(self, name, b)
                for name, b in inspect.getmembers(block_types)
                if (inspect.isclass(b) and issubclass(b, _generic_gromos_block))
            ]

    def get_all_blocks(self) -> dict:
        return vars(self)


all_blocks = all_blocks_class()
