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
import tempfile
import mdtraj
import pandas as pd
import numpy as np
import nglview as nj
from copy import deepcopy
from typing import Dict, List, Tuple
from pygromos.utils import bash
from pygromos.files.blocks._general_blocks import TITLE
from pygromos.files.coord.cnf import Cnf
from pygromos.visualization.coordinates_visualization import visualize_system


class Trc(mdtraj.Trajectory):
    # Attributes
    TITLE: TITLE

    _future_file: bool
    path: str  # Todo: we need to set this variable.

    def __init__(
        self,
        xyz=None,
        topology=None,
        time=None,
        unitcell_lengths=None,
        unitcell_angles=None,
        traj_path=None,
        in_cnf: [str, Cnf] = None,
    ):
        print("traj_path", traj_path)
        print("in_cnf", in_cnf)
        self._future_file = False

        if traj_path is not None and (traj_path.endswith(".h5") or traj_path.endswith(".hf5")):
            trj = self.load(traj_path)
            self.__dict__.update(vars(trj))

        elif traj_path is not None and (traj_path.endswith(".trc") or traj_path.endswith(".trc.gz")):

            # Parse TRC
            compress = False
            if traj_path.endswith(".gz"):
                traj_path = bash.compress_gzip(in_path=traj_path, extract=True)
                compress = True

            if isinstance(traj_path, str):
                xyz, time, step = self.parse_trc_efficiently(traj_path)

            if compress:
                traj_path = bash.compress_gzip(in_path=traj_path)

            # Topology from Cnf
            if isinstance(in_cnf, str):
                in_cnf = Cnf(in_cnf)
            elif isinstance(in_cnf, Cnf) and hasattr(in_cnf, "POSITION"):
                pass
            else:
                in_cnf = self.get_dummy_cnf(xyz)

            # Topo tmp file
            tmpFile = tempfile.NamedTemporaryFile(suffix="_tmp.pdb")
            in_cnf.write_pdb(tmpFile.name)
            single = mdtraj.load_pdb(tmpFile.name)

            tmpFile.close()

            super().__init__(xyz=xyz, topology=single.topology, time=time)
            self._step = step
        elif not (xyz is None and topology is None):
            super().__init__(xyz, topology, time, unitcell_lengths=unitcell_lengths, unitcell_angles=unitcell_angles)

        else:
            self._unitcell_lengths = []
            self._unitcell_angles = []
            self._xyz = np.array([], ndmin=2)
            self._topology = None
            self._future_file = True

        self.path = traj_path

    def __copy__(self):
        attribs = {
            "xyz": deepcopy(self._xyz),
            "topology": deepcopy(self._topology),
            "path": deepcopy(self.path),
            "_future_file": self._future_file,
        }
        for additional_key in ["unitcell_angles", "unitcell_angles"]:
            if hasattr(self, additional_key) and getattr(self, additional_key) is not None:
                attribs.update({additional_key: deepcopy(getattr(self, additional_key))})

        return self.__class__(**attribs)

    def get_dummy_cnf(self, xyz) -> Cnf:
        from pygromos.files.blocks import coords

        new_Cnf = Cnf(None)
        new_Cnf.add_block(blocktitle="TITLE", content="THis is a dummy top depending on the first Frame.")
        new_Cnf.add_block(
            blocktitle="POSITION",
            content=[
                coords.atomP(resID=1, resName="DUM", atomType="C", atomID=i, xp=coord[0], yp=coord[1], zp=coord[2])
                for i, coord in enumerate(xyz[0])
            ],
        )
        return new_Cnf

    def parse_trc_efficiently(self, traj_path: str) -> (np.array, np.array, np.array):
        self._block_map = self._generate_blockMap(in_trc_path=traj_path)

        # build block mapping
        rep_time = 1
        start = self._block_map["TIMESTEP"]
        title = self._block_map["TITLE"]
        end = start + self._block_map["POSITIONRED"] - 1
        timestep_block_length = sum(
            [self._block_map[key] for key in self._block_map if (key not in ["TITLE", "commentLines"])]
        )
        chunk = self._block_map["POSITIONRED"] - self._block_map["commentLines"] - 2 + 1

        # block mapping logic
        rows = (
            lambda x: not (
                (((x - title) % timestep_block_length > start) and ((x - title) % timestep_block_length < end))
                or (x - title) % timestep_block_length == rep_time
            )
            if x > title
            else True
        )

        # parsing
        data = []
        time = []
        step = []
        for b in pd.read_table(
            traj_path, delim_whitespace=True, skiprows=rows, names=["x", "y", "z"], chunksize=chunk, comment="#"
        ):
            data.append(b.values[1:, :])
            time.append(b.values[0, 1])
            step.append(b.values[0, 0])

        # make np.Arrays
        xyz = np.array(data)
        time = np.array(time)
        step = np.array(step)

        return xyz, step, time

    @property
    def step(self) -> np.array:
        return self._step

    # Analysis of traj
    def rmsd(self, reference_frame: int = 0, reference: mdtraj.Trajectory = None) -> pd.DataFrame:
        if reference is None:
            reference = self
        time_scale = pd.Series(data=self.time, name="time")
        return pd.DataFrame({"rmsd": mdtraj.rmsd(self, reference, reference_frame)}, index=time_scale)

    def distances(
        self,
        atom_pairs: List[Tuple[int, int]],
        periodic: bool = True,
        opt: bool = True,
    ) -> pd.DataFrame:
        arr = mdtraj.compute_distances(self, atom_pairs=atom_pairs, periodic=periodic, opt=opt)
        time_scale = pd.Series(data=self.time, name="time")
        return pd.DataFrame(
            {str(key[0]) + "-" + str(key[1]): val for key, val in zip(atom_pairs, arr.T)}, index=time_scale
        )

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

                if line.startswith("#"):
                    nCommentLines += 1

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

        block_map.update({"commentLines": nCommentLines})

        return block_map

    # Visualizaton
    @property
    def view(self, re_create: bool = False) -> nj.NGLWidget:
        if not hasattr(self, "_view") or not isinstance(self._view, nj.NGLWidget) or re_create:
            self._view = nj.show_mdtraj(self)
            self.view.clear_representations()
            self.view.add_representation("licorice", selection="all")
        return self._view

    def recreate_view(self) -> nj.NGLWidget:
        self._view = visualize_system(traj=self)
        return self._view

    # io
    def write_trc(self, out_path: str) -> str:
        raise NotImplementedError("Not Implemented")

    def save(self, out_path: str) -> str:

        compress = False
        if out_path.endswith(".trc.gz"):
            compress = True
            out_path.replace(".gz", "")

        # write out
        if out_path.endswith(".trc"):
            out_path = self.write_trc(out_path)
            # compress if desired
            if compress:
                out_path = bash.compress_gzip(in_path=out_path)
            return out_path
        else:
            super().save(out_path)
            return out_path

    def write(self, out_path: str) -> str:
        return self.save(out_path)
