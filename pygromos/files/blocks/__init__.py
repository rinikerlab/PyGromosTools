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

from pygromos.files.blocks._general_blocks import _generic_gromos_block, _iterable_gromos_block, TIMESTEP, TITLE, TRAJ

# forward declarations
TITLE: TITLE = TITLE
TIMESTEP: TIMESTEP = TIMESTEP
TRAJ: TRAJ = TRAJ
_generic_gromos_block: _generic_gromos_block = _generic_gromos_block
_iterable_gromos_block: _iterable_gromos_block = _iterable_gromos_block


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
