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
import gzip
import tempfile
import mdtraj
import pandas as pd
import numpy as np
import nglview as nj
from copy import deepcopy

from pygromos.files.blocks._general_blocks import TITLE
from pygromos.files.coord.cnf import Cnf
from pygromos.files.blocks.coord_blocks import POSITION
from pygromos.utils.typing import Dict, List, Tuple, Union
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
        in_cnf: Union[str, Cnf] = None,
        timestep_duration: float = 0.002,
        _future_file: bool = False,
    ):

        self._future_file = _future_file

        if traj_path is not None:
            if traj_path.endswith(".h5") or traj_path.endswith(".hf5"):
                trj = self.load(traj_path)
                self.__dict__.update(vars(trj))

            elif traj_path.endswith(".trc") or traj_path.endswith(".trc.gz"):

                parser = TrcParser()
                # Parse TRC
                if traj_path.endswith(".gz"):
                    parser.load_gzipped(traj_path)
                else:
                    parser.load_filename(traj_path)

                xyz = parser.positions

                # Topology from Cnf
                if isinstance(in_cnf, str):
                    in_cnf = Cnf(in_cnf)
                elif not (isinstance(in_cnf, Cnf) and hasattr(in_cnf, "POSITION")):
                    in_cnf = self.get_dummy_cnf(xyz)

                unitcell_angles = parser.unitcell_angles
                unitcell_lengths = parser.unitcell_lengths
                # get cnf boxDims
                if hasattr(in_cnf, "GENBOX") and len(unitcell_angles) == 0:
                    unitcell_angles = np.full((len(xyz), 3), in_cnf.GENBOX.angles)
                if hasattr(in_cnf, "GENBOX") and len(parser.unitcell_lengths) == 0:
                    unitcell_lengths = np.full((len(xyz), 3), in_cnf.GENBOX.length)

                # Topo tmp file
                tmpFile = tempfile.NamedTemporaryFile(suffix="_tmp.pdb")
                in_cnf.write_pdb(tmpFile.name)
                single = mdtraj.load_pdb(tmpFile.name)
                tmpFile.close()

                self._step = np.array(parser.step)
                self.TITLE = TITLE("\n".join(parser.title))
                super().__init__(
                    xyz=xyz,
                    topology=single.topology,
                    time=time or parser.time or None,
                    unitcell_lengths=unitcell_lengths if np.size(unitcell_lengths) else None,
                    unitcell_angles=unitcell_angles if np.size(unitcell_angles) else None,
                )
            else:
                raise ValueError("Unknown file extension for file: " + str(traj_path))

        elif not (xyz is None and topology is None):
            super().__init__(
                xyz=xyz,
                topology=topology,
                time=time,
                unitcell_lengths=unitcell_lengths,
                unitcell_angles=unitcell_angles,
            )

            self._step = np.array(np.round(self._time / timestep_duration), dtype=int)
            self.TITLE = TITLE(content=" Generic Title... to be changed by YOU!")

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
            "unitcell_lengths": deepcopy(self._unitcell_lengths),
            "unitcell_angles": deepcopy(self._unitcell_angles),
            "traj_path": deepcopy(self.path),
            "_future_file": self._future_file,
        }
        for additional_key in ["unitcell_angles", "unitcell_angles"]:
            if hasattr(self, additional_key) and getattr(self, additional_key) is not None:
                attribs.update({additional_key: deepcopy(getattr(self, additional_key))})

        cCls = self.__class__(**attribs)
        cCls._step = deepcopy(self._step)
        cCls.TITLE = deepcopy(self.TITLE)
        cCls._time = deepcopy(self._time)

        return cCls

    def __getitem__(self, key: int):

        t = self.slice(key)
        if hasattr(t, "_step"):
            t._step = deepcopy(np.array(self._step[key], ndmin=1))
        if hasattr(t, "_time"):
            t._time = deepcopy(np.array(self._time[key], ndmin=1))
        if hasattr(t, "TITLE"):
            t.TITLE = deepcopy(self.TITLE)

        return t

    def get_dummy_cnf(self, xyz: np.array) -> Cnf:
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
        """Computes the distance between two atoms using mdtraj. Hint: atoms are 0 indexed.

        Args:
            atom_pairs (List[Tuple[int, int]]): Indices of atoms.
            periodic (bool, optional): If periodic is True and the trajectory contains unitcell information, we will treat angles that cross periodic images using the minimum image convention. Defaults to True.
            opt (bool, optional): Use an optimized native library to calculate distances. Our optimized SSE angle calculation implementation is 10-20x faster than the (itself optimized) numpy implementation. Defaults to True.

        Returns:
            pd.DataFrame: Distances of the two atoms over the trajectory.
        """
        arr = mdtraj.compute_distances(self, atom_pairs=atom_pairs, periodic=periodic, opt=opt)
        time_scale = pd.Series(data=self.time, name="time")
        return pd.DataFrame(
            {str(key[0]) + "-" + str(key[1]): val for key, val in zip(atom_pairs, arr.T)}, index=time_scale
        )

    def angles(
        self, atom_pairs: List[Tuple[int, int, int]], degrees: bool = True, periodic: bool = True, opt: bool = True
    ) -> pd.DataFrame:
        """Computes the angle between two three using mdtraj. Hint: atoms are 0 indexed.

        Args:
            atom_pairs (List[Tuple[int, int, int]]): Indices of atoms.
            degrees (bool, optional): convert to degrees or return radians
            periodic (bool, optional): If periodic is True and the trajectory contains unitcell information, we will treat angles that cross periodic images using the minimum image convention. Defaults to True.
            opt (bool, optional): Use an optimized native library to calculate distances. Our optimized SSE angle calculation implementation is 10-20x faster than the (itself optimized) numpy implementation. Defaults to True.

        Returns:
            pd.DataFrame: Angles of the two atoms over the trajectory.
        """
        arr = mdtraj.compute_angles(self, angle_indices=atom_pairs, periodic=periodic, opt=opt)
        time_scale = pd.Series(data=self.time, name="time")
        df = pd.DataFrame(
            {str(key[0]) + "-" + str(key[1]) + "-" + str(key[2]): val for key, val in zip(atom_pairs, arr.T)},
            index=time_scale,
        )
        if degrees:
            df[df.columns[0:]] = df[df.columns[0:]].apply(lambda x: np.rad2deg(x))
        return df

    def dihedrals(
        self, atom_pairs: List[Tuple[int, int, int, int]], degrees: bool = True, periodic: bool = True, opt: bool = True
    ) -> pd.DataFrame:
        """Computes the dihedrals between four atoms using mdtraj. Hint: atoms are 0 indexed.

        Args:
            atom_pairs (List[Tuple[int, int, int, int]]): Indices of atoms.
            degrees (bool, optional): convert to degrees or return radians
            periodic (bool, optional): If periodic is True and the trajectory contains unitcell information, we will treat angles that cross periodic images using the minimum image convention. Defaults to True.
            opt (bool, optional): Use an optimized native library to calculate distances. Our optimized SSE angle calculation implementation is 10-20x faster than the (itself optimized) numpy implementation. Defaults to True.

        Returns:
            pd.DataFrame: Dihedrals of the two atoms over the trajectory.
        """
        arr = mdtraj.compute_dihedrals(self, indices=atom_pairs, periodic=periodic, opt=opt)
        time_scale = pd.Series(data=self.time, name="time")
        df = pd.DataFrame(
            {
                str(key[0]) + "-" + str(key[1]) + "-" + str(key[2]) + "-" + str(key[3]): val
                for key, val in zip(atom_pairs, arr.T)
            },
            index=time_scale,
        )
        if degrees:
            df[df.columns[0:]] = df[df.columns[0:]].apply(lambda x: np.rad2deg(x))
        return df


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

    _formatting = np.vectorize(np.format_float_positional)

    def write_trc(self, out_path: str) -> str:

        first_entry = self.generate_TITLE_entry()
        array_list = [first_entry]

        for i in range(len(self.time)):
            array_list.append(self.generate_entry_for_frame(i))

        array = np.concatenate(array_list, axis=0)
        array = np.array(array, dtype=str)
        array[array == "None"] = ""

        np.savetxt(out_path, array, fmt="%s", delimiter="\t")

    def generate_entry_for_frame(self, frame_id: int):
        length = 3  # TIMESTEP
        length += 2  # Positionred
        length += self.xyz[0].shape[0]

        # Add comments
        length += self.xyz[0].shape[0] // 10

        has_unitcells = self.unitcell_lengths is not None and len(self.unitcell_lengths) > 0
        if has_unitcells:
            length += 7

        array = np.empty((length, 4), dtype=object)

        # Add TIMESTEP
        array[0, 0] = "TIMESTEP"
        array[1, 1] = self.step[frame_id]
        array[1, 2] = self._formatting(self.time[frame_id], precision=9, unique=False, pad_left=2)
        array[2, 0] = "END"

        # Add POSITIONRED
        array[3, 0] = "POSITIONRED"
        start = 4

        for block in range(self.xyz[frame_id].shape[0] // 10):
            array[start + block * 10 : start + (block + 1) * 10, 1:] = self._formatting(
                self.xyz[frame_id][block * 10 : (block + 1) * 10, :], precision=9, unique=False, pad_left=2
            )
            array[start + (block + 1) * 10, 0] = "#"
            array[start + (block + 1) * 10, 1] = (block + 1) * 10
            start += 1

        last = self.xyz[frame_id].shape[0] % 10

        if has_unitcells:
            array_last = last + 7
        else:
            array_last = last

        array[-(array_last + 1) : -(array_last - last) - 1, 1:] = self._formatting(
            self.xyz[frame_id][-last:, :], precision=9, unique=False, pad_left=2
        )
        array[-(array_last - last) - 1, 0] = "END"

        if has_unitcells:
            array[-(array_last - last), 0] = "GENBOX"
            array[-(array_last - last) + 1, 1] = 1
            array[-(array_last - last) + 2, 1:] = self._formatting(
                self.unitcell_lengths[frame_id], precision=9, unique=False, pad_left=2
            )
            array[-(array_last - last) + 3, 1:] = self._formatting(
                self.unitcell_angles[frame_id], precision=9, unique=False, pad_left=2
            )
            array[-(array_last - last) + 4, 1:] = self._formatting(0, precision=9, unique=False, pad_left=2)
            array[-(array_last - last) + 5, 1:] = self._formatting(0, precision=9, unique=False, pad_left=2)
            array[-(array_last - last) + 6, 0] = "END"

        return array

    def generate_TITLE_entry(self):
        length = 3  # TITLE will only be one line

        # Add comments
        array = np.empty((length, 4), dtype=object)

        # Add Title
        array[0, 0] = "TITLE"

        titlestring = ""
        for t in self.TITLE.content:
            titlestring += t

        array[1, 0] = titlestring
        array[2, 0] = "END"

        return array

    def to_cnf(self, frame_id: int = None, base_cnf: Cnf = None):
        from pygromos.files.blocks import coords

        if frame_id is None:
            frame_id = 0

        content_str = (
            "THIS IS THE FRAME AT TIMESTEP: "
            + str(self.time[frame_id])
            + " OF THE TRAJECTORY WITH TITLE: \n"
            + "\n".join(self.TITLE.content)
        )

        if base_cnf is None:
            new_Cnf = self.get_dummy_cnf(self.xyz)
        else:
            if type(base_cnf) is Cnf:
                new_Cnf = base_cnf
            else:
                new_Cnf = Cnf(base_cnf)
        new_Cnf.TITLE = TITLE(content_str)

        new_Cnf.POSITION = POSITION(
            [
                coords.atomP(
                    resID=new_Cnf.POSITION.content[i].resID,
                    resName=new_Cnf.POSITION.content[i].resName,
                    atomType=new_Cnf.POSITION.content[i].atomType,
                    atomID=i,
                    xp=coord[0],
                    yp=coord[1],
                    zp=coord[2],
                )
                for i, coord in enumerate(self.xyz[frame_id])
            ]
        )

        if hasattr(new_Cnf, "GENBOX"):
            new_Cnf.GENBOX.length = list(self.unitcell_lengths[frame_id])
            new_Cnf.GENBOX.angles = list(self.unitcell_angles[frame_id])
            new_Cnf.GENBOX.euler = [0.0, 0.0, 0.0]
            new_Cnf.GENBOX.origin = [0.0, 0.0, 0.0]
        else:
            from pygromos.files.blocks.coord_blocks import GENBOX

            box_block = GENBOX(
                pbc=1, length=list(self.unitcell_lengths[frame_id]), angles=list(self.unitcell_angles[frame_id])
            )
            new_Cnf.add_block(block=box_block)

        if hasattr(new_Cnf, "TIMESTEP"):
            new_Cnf.TIMESTEP.t = self.time[frame_id]
            new_Cnf.TIMESTEP.step = self.step[frame_id]

        return new_Cnf

    def save(self, out_path: str) -> str:

        # write out
        if out_path.endswith(".trc") or out_path.endswith(".trc.gz"):
            out_path = self.write_trc(out_path)
            return out_path
        else:
            super().save(out_path)
            return out_path

    def write(self, out_path: str) -> str:
        return self.save(out_path)

    @classmethod
    def load(cls, in_path: str, in_cnf_path: str = None, timestep_duration: float = 0.002) -> any:

        if in_path.endswith(".trc") or in_path.endswith(".trc.gz"):
            o = cls(traj_path=in_path, in_cnf_path=in_cnf_path)
        else:
            so = super().load(in_path)
            o = cls(
                xyz=so.xyz,
                topology=so.topology,
                time=so.time,
                unitcell_lengths=so.unitcell_lengths,
                unitcell_angles=so.unitcell_angles,
                timestep_duration=timestep_duration,
            )

        return o


class TrcParser:
    def __init__(self):
        self.title = []
        self.step = []
        self.time = []
        self.positions = []
        self.unitcell_lengths = []
        self.unitcell_angles = []

    def load_filename(self, fname):
        with open(fname) as f:
            return self.load(f)

    def load_gzipped(self, fname):
        with gzip.open(fname, mode='rt') as f:
            return self.load(f)

    def load(self, fh):
        parsers = {
            'TITLE': self.parse_title,
            'TIMESTEP': self.parse_timestep,
            'POSITIONRED': self.parse_positions,
            'GENBOX': self.parse_genbox,
            '': lambda x: None  # just skip the line
        }
        for line in fh:
            blocktype = line.strip()
            parsers[blocktype](fh)

    def parse_positions(self, fh):
        positions = []
        for line in fh:
            if line[0] == '#':
                continue
            entries = line.split()
            if len(entries) == 1 and entries[0] == "END":
                break
            assert len(entries) == 3, "Wrong number of entries in a coordinate line: " + line
            positions.append((float(entries[0]), float(entries[1]), float(entries[2])))
        self.positions.append(positions)

    def parse_timestep(self, fh):
        entries = next(fh).split()
        if len(entries) != 2:
            raise ValueError("Invalid number of entries in timestep: " + len(entries))
        self.step.append(int(entries[0]))
        self.time.append(float(entries[1]))
        if next(fh).strip() != "END":
            raise ValueError("More than 1 line in TIMESTEP entry")

    def parse_genbox(self, fh):
        next(fh)
        lengths = tuple(float(l) for l in next(fh).split())
        angles = tuple(float(a) for a in next(fh).split())
        next(fh)
        next(fh)
        if next(fh).strip() != "END":
            raise ValueError("Wrong number of lines in GENBOX entry")
        self.unitcell_lengths.append(lengths)
        self.unitcell_angles.append(angles)

    def parse_title(self, fh):
        for line in fh:
            if line.strip() == "END":
                break
            self.title.append(line.strip())