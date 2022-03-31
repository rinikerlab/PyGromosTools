"""
FUNCTIONLIB:            gromos++ input file functions
Description:
    in this lib, gromosXX input file mainpulating functions are gathered
Author: Kay Schaller & Benjamin Schroeder
"""
import numpy as np
import copy
import json

from pygromos.files._basics import _general_gromos_file, parser
from pygromos.files.blocks import imd_blocks as blocks
from pygromos.utils.utils import nice_s_vals
from pygromos.utils.typing import List, Iterable, Union, Number


class Imd(_general_gromos_file._general_gromos_file):
    _gromos_file_ending: str = "imd"

    _orig_file_path: str
    path: str
    _required_blocks = [
        "TITLE",
        "SYSTEM",
        "STEP",
        "BOUNDCOND",
        "FORCE",
        "CONSTRAINT",
        "PAIRLIST",
        "NONBONDED",
    ]  # Not yet implemented

    # POSSIBLE GROMOS BLOCKS
    TITLE: blocks.TITLE

    STEP: blocks.STEP
    WRITETRAJ: blocks.WRITETRAJ
    PRINTOUT: blocks.PRINTOUT

    SYSTEM: blocks.SYSTEM
    FORCE: blocks.FORCE

    ENERGYMIN: blocks.ENERGYMIN
    EDS: blocks.EDS

    CONSTRAINT: blocks.CONSTRAINT
    BOUNDCOND: blocks.BOUNDCOND
    COMTRANSROT: blocks.COMTRANSROT

    PRESSURESCALE: blocks.PRESSURESCALE
    MULTIBATH: blocks.MULTIBATH

    PAIRLIST: blocks.PAIRLIST
    NONBONDED: blocks.NONBONDED

    POSITIONRES: blocks.POSITIONRES
    DISTANCERES: blocks.DISTANCERES
    INITIALISE: blocks.INITIALISE

    REPLICA_EDS: blocks.REPLICA_EDS
    NEW_REPLICA_EDS: blocks.NEW_REPLICA_EDS

    REPLICA: blocks.REPLICA

    QMMM: blocks.QMMM

    def __init__(self, in_value: str, _future_file: bool = False):
        super().__init__(in_value=in_value, _future_file=_future_file)

        # TODO: maybe somebody can make a better solution for this. This is a ugly fix to unify the structure of the blocks
        for block in sorted(self.get_block_names()):
            setattr(self, block, copy.deepcopy(getattr(self, block)))

    def __str__(self):
        text = ""
        if hasattr(self, "TITLE"):
            text += self.__getattribute__("TITLE").block_to_string()
        for block in sorted(self.get_block_names()):
            if block == "TITLE" or isinstance(block, type(None)):
                continue
            text += str(self.__getattribute__(block))
        return text

    def read_file(self):
        data = parser.read_imd(self._orig_file_path)
        for key, sub_content in data.items():
            try:
                self.add_block(blocktitle=key, content=sub_content)
            except Exception as err:
                raise Exception("Error while reading file: " + str(err))
        return {}

    def edit_EDS(
        self,
        NUMSTATES: int,
        S: float,
        EIR: list,
        EDS: int = 1,
        ALPHLJ: float = 0.0,
        ALPHC: float = 0.0,
        FUNCTIONAL: int = 1,
    ):
        if hasattr(self, "eds_block"):
            self.EDS.ALPHLJ = ALPHLJ
            self.EDS.ALPHC = ALPHC
            self.EDS.FUNCTIONAL = FUNCTIONAL
            self.EDS.NUMSTATES = NUMSTATES
            self.EDS.S = S
            self.EDS.EIR = EIR

        else:
            print("Setting new EDS_block")
            if type(EIR) == float or type(EIR) == str or type(EIR) == int:
                EIR = [float(EIR) for x in range(NUMSTATES)]

            gromos_name = "EDS"
            setattr(self, gromos_name, blocks.EDS(NUMSTATES, S, EIR, EDS, ALPHLJ, ALPHC, FUNCTIONAL))

    def randomize_seed(self):
        self.INITIALISE.IG = np.random.randint(low=0, high=999999)

    def edit_REEDS(
        self,
        REEDS: int = None,
        NUMSTATES: int = None,
        SVALS: Union[Number, List[Number]] = None,
        EIR: Union[Number, Iterable[Number]] = None,
        NRETRIAL: int = None,
        NREQUIL: int = None,
        CONT: Union[bool, int] = None,
        EDS_STAT_OUT: Union[bool, int] = None,
    ):

        # specific relations are rescued here
        reeds_block = self.REPLICA_EDS
        print(type(reeds_block))

        if isinstance(REEDS, bool):
            reeds_block.REEDS = REEDS

        if isinstance(NUMSTATES, (Number, str)):
            reeds_block.NUMSTATES = NUMSTATES

        if isinstance(SVALS, Iterable):  # edit SVALS
            SVALS = nice_s_vals(SVALS)

            reeds_block.RES = list(map(str, SVALS))
            reeds_block.NRES = len(SVALS)  # adjust number of Svals

            if isinstance(EIR, (Number, Iterable)):  # expand energy offsets to new s value ammount
                # set new EIR with 3 different types of input (single number, vector or matrix)
                EIR_matrix = []

                # single number
                if isinstance(EIR, Number):  # depends on SVALS and NRES
                    EIR_vector = [EIR for x in range(reeds_block.NUMSTATES)]

                    for z in EIR_vector:
                        EIR_matrix.append([z for i in range(int(reeds_block.NRES))])

                # vector or matrix
                elif (
                    isinstance(EIR, Iterable)
                    and len(EIR) == int(reeds_block.NUMSTATES)
                    and all([isinstance(x, Number) for x in EIR])
                ):
                    for z in EIR:
                        EIR_matrix.append([z for i in range(int(reeds_block.NRES))])
                else:
                    raise Exception(
                        "not enough EIR-vals for making REEDS Block. Got "
                        + str(len(EIR))
                        + " for "
                        + str(reeds_block.NRES)
                        + " SVals\n Number of states: "
                        + str(reeds_block.NUMSTATES)
                    )
                reeds_block.EIR = EIR_matrix
            else:
                if any([len(row) != len(SVALS) for row in reeds_block.EIR]):
                    newEIR = []
                    for EIR_row in reeds_block.EIR:
                        newEIRrow = [EIR_row[0] for r in range(len(SVALS))]
                        newEIR.append(newEIRrow)
                    reeds_block.EIR = newEIR

        if isinstance(EIR, (Number, Iterable)):
            EIR_matrix = []
            # single number
            if isinstance(EIR, Number):  # depends on SVALS and NRES
                EIR_vector = [str(EIR) for x in range(reeds_block.NUMSTATES)]

                for z in EIR_vector:
                    EIR_matrix.append([z for i in range(int(reeds_block.NRES))])

            # vector
            elif (
                isinstance(EIR, Iterable)
                and len(EIR) == int(reeds_block.NUMSTATES)
                and all([isinstance(x, (Number, str)) for x in EIR])
            ):
                for z in EIR:
                    EIR_matrix.append([float(z) for i in range(int(reeds_block.NRES))])
            # matrix
            elif isinstance(EIR, Iterable) and all(
                [isinstance(x, Iterable) and all([isinstance(y, (Number, str)) for y in x]) for x in EIR]
            ):
                if len(self.REPLICA_EDS.NRES) == len(EIR):
                    EIR = np.array(EIR).T
                EIR_matrix = list(map(lambda x: list(map(float, x)), EIR))
            else:
                raise Exception(
                    "not enough EIR-vals for making REEDS Block. Got "
                    + str(len(EIR))
                    + " for "
                    + str(reeds_block.NRES)
                    + " SVals and states: "
                    + str(reeds_block.NUMSTATES)
                    + "\n"
                )
            reeds_block.EIR = EIR_matrix

        if isinstance(NRETRIAL, (str, int)):
            reeds_block.NRETRIAL = NRETRIAL
        if isinstance(NREQUIL, (str, int)):
            reeds_block.NREQUIL = NREQUIL
        if isinstance(CONT, (str, int, bool)):
            reeds_block.CONT = CONT
        if isinstance(EDS_STAT_OUT, (str, int, bool)):
            reeds_block.EDS_STAT_OUT = EDS_STAT_OUT

    def write_json(self, out_path: str):
        d = copy.deepcopy(vars(self))

        for v in vars(self):
            if isinstance(d[v], (blocks.TITLE, blocks._generic_imd_block, blocks._generic_gromos_block)) or issubclass(
                d[v].__class__, _general_gromos_file._general_gromos_file
            ):
                d.update({v: vars(d[v])})

        json.dump(d, open(out_path, "w"))
        return out_path
