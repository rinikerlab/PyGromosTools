"""
File:            Class for trc files in pandas
Description:
    The pandas trajectory TRC class offers a easy method to process GROMOS's .trc files in python
    The trc files are parsed into an easy to use pandas dataframe

Author: Marc Thierry Lehner

TODO: add stride option to all member functions
TODO: add support for periodic boundary condition
TODO: add support for rdkit mol selector
TODO: add support for rdkit conformers

"""

# imports
import os
import pandas as pd
import numpy as np
from typing import TypeVar, Union, Dict
from pandas.core.base import DataError
import pygromos.files.trajectory._general_trajectory as traj
from pygromos.files.coord.cnf import Cnf
from pygromos.analysis import coordinate_analysis as ca

from pygromos.files.blocks._general_blocks import TITLE

import nglview as nj
import mdtraj

TrcType = TypeVar("Trc")
CnfType = TypeVar("Cnf")


class Trc(mdtraj.Trajectory):
    # Attributes
    TITLE: TITLE

    def __init__(
        self,
        xyz=None,
        topology=None,
        time=None,
        unitcell_lengths=None,
        unitcell_angles=None,
        traj_path=None,
        in_cnf : [str,Cnf]=None,
    ):
        if not (traj_path is None and in_cnf is None):
            #Parse TRC
            if isinstance(traj_path, str):
                xyz, time, step = self.parse_trc_efficiently(traj_path)

            #Topology from Cnf
            if isinstance(in_cnf, str):
                in_cnf = Cnf(in_cnf)
            in_cnf.write_pdb("/localhome/kpaul/pygromosday/test.pdb")
            single = mdtraj.load_pdb("/localhome/kpaul/pygromosday/test.pdb")


            super().__init__(xyz=xyz, topology=single.topology, time=time)
            self._step = step
        else:
            super().__init__(xyz, topology, time, unitcell_lengths=unitcell_lengths, unitcell_angles=unitcell_angles)

    def parse_trc_efficiently(self, traj_path:str)->(np.array, np.array, np.array):
        self._block_map = self._generate_blockMap(in_trc_path=traj_path)

        #build block mapping
        rep_time = 1
        start = self._block_map['TIMESTEP']
        title = self._block_map['TITLE']
        end = start + self._block_map['POSITIONRED'] - 1
        timestep_block_length = sum([self._block_map[key] for key in self._block_map if (not key in ["TITLE", "commentLines"])])
        chunk = self._block_map['POSITIONRED'] - self._block_map['commentLines'] - 2 + 1

        ##block mapping logic
        rows = lambda x: not (
                (((x - title) % timestep_block_length > start) and ((x - title) % timestep_block_length < end)) or (
                x - title) % timestep_block_length == rep_time) if x > title else True

        #parsing
        data = []
        time = []
        step = []
        for b in pd.read_table(
                traj_path,
                delim_whitespace=True, skiprows=rows, names=['x', 'y', 'z'], chunksize=chunk, comment="#"):
            data.append(b.values[1:, :])
            time.append(b.values[0, 1])
            step.append(b.values[0, 0])

        #make np.Arrays
        xyz = np.array(data)
        time = np.array(time)
        step = np.array(step)

        return xyz, step,time

    @property
    def step(self)->np.array:
        return self._step

    rmsd = lambda self, x: mdtraj.rmsd(self, self, x)

    def _generate_blockMap(self, in_trc_path: str) -> Dict[str, int]:
        block_map = {}
        with open(in_trc_path, "r") as file_handle:
            inBlock = False
            inTitleBlock = False
            blockKey = None
            nLines = 0
            nCommentLines = 0
            titleStr = []

            while True:
                line = file_handle.readline().strip()

                if(line.startswith("#")):
                    nCommentLines+=1

                if "END" == line.strip():
                    block_map.update({blockKey: nLines})
                    if blockKey == TITLE.__name__:
                        inTitleBlock = False
                        self.TITLE = TITLE(content="\n".join(titleStr))
                    inBlock = False
                elif not inBlock:
                    blockKey = line.strip()
                    if blockKey in block_map:
                        break
                    elif blockKey == TITLE.__name__:
                        inTitleBlock = True
                    inBlock = True
                    nLines = 1
                elif inTitleBlock:
                    titleStr.append(line)
                nLines += 1

        if not hasattr(self, TITLE.__name__):
            self.TITLE = TITLE(content="Empty TITLE")

        block_map.update({"commentLines":nCommentLines})

        return block_map

    @property
    def view(self, re_create: bool = False) -> nj.NGLWidget:
        if not hasattr(self, '_view') or isinstance(self, "_view", nj.NGLWidget) or re_create:
            self._view = nj.show_mdtraj(self)
        return self._view

    def write_trc(self):
        raise NotImplemented("Not Implemented")