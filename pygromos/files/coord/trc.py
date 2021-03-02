"""
FUNCTIONLIB:            gromos++ coord files
Description:
    in this lib, gromos coordinate files are wrapped in python and functions are added to modify them.

Author: Benjamin Schroeder
"""

from typing import List, Dict, Iterable

from pygromos.files._basics import parser
from pygromos.files.blocks import coord_blocks as blocks
from pygromos.files._basics._general_gromos_file import _general_gromos_trajectory_file
from pygromos.files.coord.cnf import Cnf

#Util:
#FILES
class Trc(_general_gromos_trajectory_file):

    #general
    _orig_file_path:str

    #private:
    _block_order: List[str] = ["TITLE", "TRAJ"]

    #Standard Gromos blocks
    TRAJ: blocks.TRAJ
    TITLE: blocks.TITLE

    def __init__(self, input: (str, List[Dict]), title:str=None, every_step:int=1):
        if(isinstance(input, str)):
            super().__init__(input=input, every_step=every_step)
        elif(isinstance(input, List)):
            if(isinstance(title, type(None))):
                title = "THIS IS A GENERIC GROMOS TITLE, because you did not specify one! <3"

            self._orig_file_path = ""

            self.TITLE = blocks.TITLE(content=title)
            self.TRAJ = blocks.TRAJ(timestep_blocks=input)

        elif(isinstance(input, type(None))):
            pass
        else:
            raise IOError("Don't understand the input!")

    @classmethod
    def cnfs_to_trc(cls, cnfs:Iterable[str], title=None, start_time:float=0.0, dt:float=0.002 ):
        return cls(input=cls._cnf_to_trc_dict(cnfs, start_time=start_time, dt=dt), title=title)

    @classmethod
    def _cnf_to_trc_dict(cls, cnfs:Iterable[str], start_time:float=0.0, dt:float=0.002)->List[Dict]:

        time = start_time
        out_frames_list:List[Dict] = []
        for steps, cnf_path in enumerate(cnfs):
            cnf = Cnf(cnf_path)
            frame_traj = {'steps': steps, 'time': time, "POSITIONRED": []}
            block_str = ""
            for atom in sorted(cnf.POSITION, key=lambda at: at.atomID):
                if (atom.atomID % 10 == 0): block_str += "#\t" + str(atom.atomID) + "\n"
                block_str += "\t" + str(atom.xp) + "\t" + str(atom.yp) + "\t" + str(atom.zp) + "\n"
            frame_traj.update({"POSITIONRED": block_str})
            frame_traj.update({"GENBOX": str(cnf.GENBOX).replace("GENBOX\n", "").replace("END\n", "")}) #a bit hacky

            out_frames_list.append(frame_traj)
            time += dt

        return out_frames_list

    def concat(self, trcs:(object, List[object]), renumber_timesteps:bool=False, newTitle:str=None):
        if(not isinstance(trcs, List)):
            trcs = [trcs]

        for trc in trcs:
            self.TRAJ.extend(trc.TRAJ)

        if(renumber_timesteps):
            self.clean_timestep_numbering()

        if(isinstance(newTitle, str)):
            self.TITLE = blocks.TITLE(newTitle)

    def read_file(self, every_step:int=1)-> Dict[str, any]:
        """
            Blocks which were read in.
        Returns
        -------
        Dict[str, any]
            refactored blocks
        """

        tmp_blocks = parser.read_trc(self._orig_file_path, every_step=every_step)
        refactored_blocks = {}
        for content_key in tmp_blocks:
            block = tmp_blocks[content_key]
            if (content_key == "TITLE"):
                block = blocks.TITLE(block)
                content_key = "TITLE"
            if (content_key == "TIMESTEP"):
                content_key = "TRAJ"
                block = blocks.TRAJ(list(block.values()))
            refactored_blocks.update({content_key: block})

        return refactored_blocks

    def get_dt(self):
        if(len(self.TRAJ) >2):
            step1 = self.TRAJ[1]
            step2 = self.TRAJ[2]
            return step2["time"]-step1["time"]
        else:
            raise ValueError("This out_trc file has less than two timesteps! Can not do "+__name__+" for out_trc: "+str(self._orig_file_path))

    def get_dsteps(self):
        if(len(self.TRAJ) >2):
            step1 = self.TRAJ[1]
            step2 = self.TRAJ[2]
            return step2["steps"] - step1["steps"]
        else:
            raise ValueError("This out_trc file has less than two timesteps! Can not do " + __name__ + " for out_trc: " + str(self._orig_file_path))

    def clean_timestep_numbering(self, dsteps=False, dt=False):
        """
            Assumes dT is equal in each file!
        Parameters
        ----------
        dsteps :
        dt :

        Returns
        -------

        """

        if(not dsteps and not dt):
            dt = self.get_dt()
            dsteps = self.get_dsteps()
            for ind, frame in enumerate(self.TRAJ):
                new_step = (ind)*dsteps
                new_time = (ind)*dt
                frame["steps"] = new_step
                frame["time"] = new_time
        else:
            raise NotImplementedError(__name__+" function is not implemented for: "+__class__)

    def write(self, out_path)->str:
        file = open(out_path, "w")
        data = self.TRAJ

        #write:
        file.write(str(self.TITLE))

        ##TIMESTEPS
        prior_order = ["POSITIONRED"]
        last_order = ["GENBOX"]
        for frame in data:
            key = "TIMESTEP"
            file.write(key + "\n")
            file.write("\t\t"+str(frame["steps"])+"\t"+str(frame["time"])+"\n")
            file.write("END\n")

            for y in prior_order:
                key = y
                file.write(key + "\n")
                if(type(frame[key]) is list):
                    file.write("".join(frame[key]))
                else:
                    file.write(frame[key])
                file.write("END\n")
            for y in sorted(frame):
                if( y != "steps" and y!= "time" and not y in prior_order and not y in last_order):
                    key = y
                    file.write(key + "\n")
                    if(type(frame[key]) is list):
                        file.write("".join(frame[key]))
                    else:
                        file.write(frame[key])
                    file.write("END\n")
            for y in last_order:
                key = y
                file.write(key + "\n")
                if(type(frame[key]) is list):
                    file.write("".join(frame[key]))
                else:
                    file.write(frame[key])
                file.write("END\n")

        file.close()
        return out_path
