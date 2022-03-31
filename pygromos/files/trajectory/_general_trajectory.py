"""
FUNCTIONLIB:            trajectory template with pandas
Description:
    in this parent class for all trajectories the general parsing to a pandas database is handeled.
    This class is intended as quick data evaluation in python without gromos++ programs like ene_ana

    Specific block structures are handeled in the respectve classes (trc/tre/trf/...)

Author: Marc Thierry Lehner

Remark: WIP

###########################
ATTENTION: NOT FULLY TESTET
ATTENTION: WORK IN PROGRESS
###########################

TODO: MTL: implement, check, test

"""

# imports
import collections
import os
import re
import pandas
import numpy
import pathlib

from pygromos.files.trajectory.blocks import trajectory_blocks as blocks
from pygromos.utils import bash
from pygromos.utils.typing import _General_Trajectory_Type


class _General_Trajectory:

    # attribute annotation:
    database: pandas.DataFrame
    path: str
    _future_file: bool  # if code is executed async, this helps organizing.
    _gromos_file_ending: str

    def __init__(self, input_value: str, auto_save: bool = True, stride: int = 1, skip: int = 0):
        if input_value is None:
            self.TITLE = "Empty Trajectory"
            self.database = pandas.DataFrame({"": []})

        elif isinstance(input_value, str):
            if input_value.endswith(".gz"):
                tmp_input = bash.compress_gzip(input_value, extract=True)
                self._read_from_file(input_path=tmp_input, auto_save=auto_save, stride=stride, skip=skip)
                bash.compress_gzip(tmp_input)
            else:
                self._read_from_file(input_path=input_value, auto_save=auto_save, stride=stride, skip=skip)
            self.path = input_value

        elif isinstance(input_value, list) and all([x is str for x in input_value]):
            self._read_from_file(input_path=input_value[0], auto_save=False, stride=stride, skip=skip)

            for tmp_traj_file in input_value[1:]:
                self += __class__(tmp_traj_file, auto_save=False)

            if auto_save:
                auto_save_path = input_value[0] + ".h5"
                self.write(auto_save_path)

        elif isinstance(input_value.__class__, __class__) or issubclass(input_value.__class__, __class__):
            for attr in vars(input_value):
                setattr(self, attr, getattr(input_value, attr))
            self.database = input_value.database.copy(deep=True)
        else:
            raise IOError("Constructor not found")

    def __str__(self) -> str:
        if isinstance(self.TITLE, str):
            msg = "Trajectory: \n\t" + "\n\t".join(self.TITLE.split("\n")) + "\n"
        elif isinstance(self.TITLE, list):
            msg = "Trajectory: \n\t" + "\n\t".join(self.TITLE) + "\n"
        else:
            print(type(self.TITLE), self.TITLE)
            msg = "Trajectory: \n\t" + "\n\t".join(str(self.TITLE).split("\n")) + "\n"

        msg += "Type: \n\t" + str(self.__class__.__name__) + "\n"
        msg += "Frames: \t" + str(self.database.shape[0]) + "\t Columns:\t" + str(self.database.shape[1]) + "\n"
        msg += "\n"
        return msg

    def __repr__(self):
        return str(self)

    def __add__(self, traj: _General_Trajectory_Type):
        return self.add_traj(traj)

    def __copy__(self):
        traj = type(self)(input_value=None)
        for attr in vars(self):
            setattr(traj, attr, getattr(self, attr))
        traj.database = self.database.copy(deep=False)
        return traj

    def __deepcopy__(self, memo):
        traj = type(self)(input_value=None)
        for attr in vars(self):
            setattr(traj, attr, getattr(self, attr))
        traj.database = self.database.copy(deep=True)
        return traj

    def add_traj(
        self,
        traj: _General_Trajectory_Type,
        skip_new_0: bool = False,
        auto_detect_skip: bool = True,
        correct_time: bool = False,
    ):
        """Combine (Catenate) two trajectories to a longer trajectory. Important: A+B!=B+A

        Parameters
        ----------
        traj : trajectory of the same format and type as the base class
            the trajectory to be added to the base class

        skip_new_0 : bool
            if True, the first frame of the new trajectory will be skipped. This will override the auto_detect argument

        auto_detect_skip : bool
            compares the time of the first frame of the new trajectory to the time of the last frame of the base class and skips the frame if they are the same

        correct_time : bool
            if True, the time of the new trajectory will be corrected to continue after the last frame of the base class
            Not Implemented currently!

        Returns
        -------
        traj
            trajA + trajB
        """
        if type(traj) != type(self):
            raise Exception(
                "Different types of trajectories can not be added (concatenated).\nTry to add a "
                + str(type(self))
                + " class of the same type"
            )
        if traj.database.shape[1] != self.database.shape[1]:
            raise Warning(
                "trajectories database shapes do not match!\n Please check if this is expected\n"
                "first shape: " + str(self.database.shape) + "\tsecond shape: " + str(traj.database.shape) + "\n"
            )
        # get end data from first trajectory
        time_offset = float(self.database.time.iloc[-1])
        delta_time_self = self.get_time_step()

        # copy and modify second trajectory
        new_data = traj.database.copy(deep=True)

        if skip_new_0:
            new_data = new_data.iloc[1:]
        elif auto_detect_skip:
            new_frame = new_data.iloc[0].iloc[2:]
            old_frame = self.database.iloc[-1].iloc[2:]
            if (
                new_frame.equals(old_frame)
            ) and new_frame.keys() == old_frame.keys():  # check if the firstStep==lastStep without considering the time
                if all([numpy.allclose(new_frame[x], old_frame[x]) for x in new_frame.keys()]):
                    new_data = new_data.iloc[1:]
            elif correct_time:
                if delta_time_self == traj.get_time_step():
                    time_offset += delta_time_self

        if correct_time:
            # is the time access all right?
            # last_step_traj1 = int(self.database.time.iloc[-1])
            # first_step_traj2 = int(new_data.database.time.iloc[0])
            # consecutive_trajs = last_step_traj1+delta_time_self == first_step_traj2
            # correct the axis:
            #
            raise NotImplementedError("This is not implemented currently!")

        # create output trajectory (copy of first traj) and combine trajectories
        new_traj = self.__class__(input_value=self)
        new_traj.database = self.database.append(new_data, ignore_index=True)
        del new_data
        return new_traj

    def _read_from_file(self, input_path: str, auto_save: bool = True, stride: int = 1, skip: int = 0):
        if input_path.endswith(".h5"):
            self._read_db_from_hf5(input_path=input_path)
        elif re.search("\.tr.$", input_path):  # noqa: W605
            self._read_trajectory(input_path=input_path, auto_save=auto_save, stride=stride, skip=skip)
        else:
            raise IOError("Did not understand the file ending of given file path: ", input_path)

    def _read_db_from_hf5(self, input_path: str, title: str = "Read from hdf save \nContains only database\n"):
        self.TITLE = title
        self.database = pandas.read_hdf(path_or_buf=input_path, key=input_path.split(".")[-1])

    def _raw_read_trajectory(self, input_path: str, stride: int = 1, skip: int = 0) -> collections.defaultdict:
        # define temp storage
        header = {}
        data = []
        dataInTimestep = {}
        block = []
        blockname = ""

        # set contorl bool
        isInTimeStep = False
        isInBlock = False
        timeStepCounter = 0 - skip

        # start to read the trajectory
        with open(input_path, "r") as infile:
            for line in infile:
                if isInTimeStep:
                    if timeStepCounter >= 0 and timeStepCounter % stride == 0:
                        if isInBlock:
                            if not line.strip().startswith("END"):
                                block.append(line)
                            else:
                                isInBlock = False
                                dataInTimestep.update({blockname: block})
                        else:
                            if line.strip().startswith("#") or line.strip() == "":
                                continue
                            blockname = line.strip()
                            block = []
                            isInBlock = True
                            if blockname.startswith("TIMESTEP"):
                                data.append(dataInTimestep)
                                dataInTimestep = {}
                                timeStepCounter += 1
                    else:
                        if line.strip().startswith("TIMESTEP"):
                            timeStepCounter += 1
                else:
                    if isInBlock:
                        if not line.strip().startswith("END"):
                            block.append(line)
                        else:
                            isInBlock = False
                            header.update({blockname: block})
                    else:
                        if line.strip().startswith("#") or line.strip() == "":
                            continue
                        blockname = line.strip()
                        block = []
                        isInBlock = True
                        if blockname.startswith("TIMESTEP"):
                            isInTimeStep = True
            if isInTimeStep:
                # final time step finish since no smarter way to detect end of Timestep
                data.append(dataInTimestep)
            else:
                raise ValueError("No timestep found")
        del dataInTimestep, block, blockname
        return {"header": header, "body": data}

    def _read_trajectory(self, input_path: str, stride: int = 1, skip: int = 0, auto_save=True):
        if auto_save:
            # check if parsed file exists and is up to date
            if os.path.isfile(input_path + ".h5"):
                if pathlib.Path(input_path).stat().st_ctime < pathlib.Path(input_path + ".h5").stat().st_ctime:
                    self._read_db_from_hf5(
                        input_path=input_path + ".h5",
                        title="Reread from hdf save \nContains only database\nfor all other blocks please make a fresh import",
                    )

        if not os.path.exists(input_path):
            raise IOError("Could not find File: ", input_path)
        else:
            table = []
            data = self._raw_read_trajectory(input_path=input_path, stride=stride, skip=skip)
            header = data["header"]
            body = data["body"]
            for key, block in header.items():
                setattr(self, key, block)
            for time_step_entry in body:
                table_entry = {}
                for blocktitle, block in time_step_entry.items():
                    if not hasattr(blocks, blocktitle):
                        raise IOError("No trajectory block found named: " + blocktitle)
                    tmp_block = getattr(blocks, blocktitle)(block)
                    table_entry.update(tmp_block.to_dict())
                table.append(table_entry)
        db = pandas.DataFrame.from_dict(table)
        del table, table_entry, data, header, body
        self.database = db
        if auto_save:
            self.write(input_path + ".h5")

    def write(self, output_path: str) -> str:
        if not output_path.endswith(".h5"):
            output_path += ".h5"
        if not os.path.exists(os.path.dirname(output_path)):
            raise IOError(
                "Could not find target directory for outfile! target dir: " + str(os.path.dirname(output_path))
            )
        self.database.to_hdf(
            path_or_buf=output_path, key=output_path.split(".")[-1]
        )  # TODO: @Marc is the key arg here correct, or rather not using it?
        self.path = output_path
        return output_path

    def get_time_step(self):
        if len(self.database) == 0:
            return 0
        elif len(self.database) == 1:
            return float(self.database.time.iloc[0])
        else:
            return float(self.database.time.iloc[1]) - float(self.database.time.iloc[0])
