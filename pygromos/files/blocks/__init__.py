import inspect

from pygromos.files.blocks import coord_blocks as coords
from pygromos.files.blocks import imd_blocks as imd
from pygromos.files.blocks import topology_blocks as topo
from pygromos.files.blocks import mtb_blocks as mtb
from pygromos.files.blocks import pertubation_blocks as ptp
from pygromos.files.blocks import qmmm_blocks as qmm
from pygromos.files.blocks import replica_exchange_blocks as repdat

from pygromos.files.blocks._general_blocks import _generic_gromos_block, _iterable_gromos_block, TIMESTEP, TITLE, TRAJ

# forward declarations
TITLE: TITLE = TITLE
TIMESTEP: TIMESTEP = TIMESTEP
TRAJ: TRAJ = TRAJ

_generic_gromos_block: _generic_gromos_block = _generic_gromos_block
_iterable_gromos_block: _iterable_gromos_block = _iterable_gromos_block

all_block_types = [coords, imd, topo, mtb, ptp, qmm, repdat]
all_blocks = {b.__name__: b for b in [_generic_gromos_block, _iterable_gromos_block, TIMESTEP, TITLE, TRAJ]}


class all_blocks_class:
    def __init__(self):
        for block_types in all_block_types:
            [
                setattr(self, name, b)
                for name, b in inspect.getmembers(block_types)
                if (inspect.isclass(b) and issubclass(b, _generic_gromos_block))
            ]

    def get_all_blocks() -> dict:
        return vars(self)


all_blocks = all_blocks_class()
