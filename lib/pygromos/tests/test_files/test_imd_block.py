import unittest
from pygromos.files.blocks import imd_blocks


class test_imd_block(unittest.TestCase):
    def test_some_blocks(self):
        test = imd_blocks.REPLICA_EDS(NRES=2, NUMSTATES=5, RES=6, EIR=8, NRETRIAL=9, NREQUIL=10, CONT=11)
        assert isinstance(test, imd_blocks.REPLICA_EDS)

        test2 = imd_blocks.MULTIBATH(
            ALGORITHM=1, NBATHS=1, TEMP0=[3], TAU=[4], DOFSET=5, LAST=[6], COMBATH=[6], IRBATH=[7]
        )
        assert isinstance(test2, imd_blocks.MULTIBATH)

        test3 = imd_blocks.PRESSURESCALE(
            COUPLE=1, SCALE=2, COMP=1, TAUP=1, VIRIAL=3, SEMIANISOTROPIC=[4], PRES0=[[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        )
        assert isinstance(test3, imd_blocks.PRESSURESCALE)

        test4 = imd_blocks.RANDOMNUMBERS(NTRNG=0, NTGSL=0)
        assert isinstance(test4, imd_blocks.RANDOMNUMBERS)
