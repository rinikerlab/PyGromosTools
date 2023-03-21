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
from pygromos.utils.typing import List, Tuple, Union
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

    def write_trc(self, out_path: str) -> None:
        writer = TrcWriter(self)
        if out_path.endswith(".gz"):
            writer.write_gzipped(out_path)
        else:
            writer.write_filename(out_path)

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
        "Write to out_path, using TrcWriter for .trc and .trc.gz, mdtraj otherwise."
        # write out
        if out_path.endswith(".trc") or out_path.endswith(".trc.gz"):
            self.write_trc(out_path)
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
        self.title: List[str] = []
        self.step: List[int] = []
        self.time: List[float] = []
        self.positions: List[List[Tuple]] = []
        self.unitcell_lengths: List[Tuple] = []
        self.unitcell_angles: List[Tuple] = []

    def load_filename(self, fname):
        "Open and read a trc file"
        with open(fname) as f:
            return self.load(f)

    def load_gzipped(self, fname):
        "Open and read a trc.gz file"
        with gzip.open(fname, mode="rt") as f:
            return self.load(f)

    def load(self, fh):
        "Read a trc file from an open file handle"
        parsers = {
            "TITLE": self.parse_title,
            "TIMESTEP": self.parse_timestep,
            "POSITIONRED": self.parse_positions,
            "GENBOX": self.parse_genbox,
            "": lambda _x: None,  # just skip the line
        }
        for line in fh:
            blocktype = line.strip()
            parsers[blocktype](fh)

    def parse_positions(self, fh):
        "Parse a positions block from an iterator yielding lines (i.e., an open file)"
        positions = []
        for line in fh:
            if line[0] == "#":
                continue
            entries = line.split()
            if len(entries) == 1 and entries[0] == "END":
                break
            assert len(entries) == 3, "Wrong number of entries in a coordinate line: " + line
            positions.append((float(entries[0]), float(entries[1]), float(entries[2])))
        self.positions.append(positions)

    def parse_timestep(self, fh):
        "Parse a timestep block from an iterator yielding lines (i.e., an open file)"
        entries = next(fh).split()
        if len(entries) != 2:
            raise ValueError("Invalid number of entries in timestep: " + len(entries))
        self.step.append(int(entries[0]))
        self.time.append(float(entries[1]))
        if next(fh).strip() != "END":
            raise ValueError("More than 1 line in TIMESTEP entry")

    def parse_genbox(self, fh):
        "Parse a genbox block from an iterator yielding lines (i.e., an open file)"
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
        "Parse a title block from an iterator yielding lines (i.e., an open file)"
        for line in fh:
            if line.strip() == "END":
                break
            self.title.append(line.strip())


class TrcWriter:
    """Class to write basic Trc trajectory info to a file.

    Supported blocks:
    * TITLE
    * POSITIONRED
    * TIMESTEP
    * GENBOX
    """
    def __init__(self, traj: Trc, float_format="{:>14.9f}"):
        self.traj = traj
        self.float_format = float_format
        self.three_float_fmt = " {} {} {}\n".format(float_format, float_format, float_format)
        self.int_and_float_fmt = " {:>14d} " + float_format + "\n"

    def write(self, fh):
        "Write all frames to an open file handle"
        self.write_title(fh)
        for frame_number in range(self.traj.n_frames):
            self.write_timestep(fh, frame_number)
            self.write_positions(fh, frame_number)
            if self.has_unitcells():
                self.write_unitcells(fh, frame_number)

    def write_timestep(self, fh, frame_number: int):
        "Write TIMESTEP for frame *frame_number* to an open file handle"
        fh.write("TIMESTEP\n")
        time = self.traj.time[frame_number]
        step = self.traj.step[frame_number]
        fh.write(self.int_and_float_fmt.format(step, time))
        fh.write("END\n")

    def write_positions(self, fh, frame_number: int):
        "Write POSITIONRED for frame *frame_number* to an open file handle"
        fh.write("POSITIONRED\n")
        positions = self.traj.xyz[frame_number].tolist()
        for atom_number in range(self.traj.n_atoms):
            fh.write(self.three_float_fmt.format(*positions[atom_number]))
            if (atom_number + 1) % 10 == 0:
                fh.write("#{:>10d}\n".format(atom_number + 1))
        fh.write("END\n")

    def write_title(self, fh):
        "Write TITLE entry to an open file handle"
        fh.write("TITLE\n")
        for line in self.traj.TITLE.content:
            fh.write("\t" + line + "\n")
        fh.write("END\n")

    def write_unitcells(self, fh, frame_number):
        "Write GENBOX entry for frame *frame_number* to an open file handle"
        fh.write("GENBOX\n")
        fh.write("    1\n")
        lengths = self.traj.unitcell_lengths[frame_number]
        angles = self.traj.unitcell_angles[frame_number]
        fh.write(self.three_float_fmt.format(*lengths))
        fh.write(self.three_float_fmt.format(*angles))
        fh.write(self.three_float_fmt.format(0.0, 0.0, 0.0))
        fh.write(self.three_float_fmt.format(0.0, 0.0, 0.0))
        fh.write("END\n")

    def write_filename(self, filename):
        "Open a file and write to it"
        with open(filename, 'wt') as f:
            self.write(f)

    def write_gzipped(self, filename):
        "Open a gzipped file and write to it"
        with gzip.open(filename, 'wt') as f:
            self.write(f)

    def has_unitcells(self) -> bool:
        return self.traj.unitcell_lengths is not None and np.size(self.traj.unitcell_lengths) > 0
