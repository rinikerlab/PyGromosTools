from typing import List, Dict

from pygromos.files._basics import parser, _general_gromos_file
from pygromos.files.blocks import energy_blocks as blocks


class Tre(_general_gromos_file._general_gromos_trajectory_file):  #

    # general
    _orig_file_path: str

    # private:
    _block_order: List[str] = ["TITLE", "TIMESTEP", "TRAJ"]

    # Standard Gromos blocks
    TRAJ: blocks.TRAJ
    TITLE: blocks.TITLE

    def __init__(self, input: (str or dict), every_step:int = 1):
        super().__init__(input=input, every_step=every_step)

    def read_file(self, every_step:int=1) -> Dict[str, any]:
        """
              Blocks which were read in. - very crude
          Returns
          -------
          Dict[str, any]
              refactored blocks
          """

        tmp_blocks = parser.read_tre(self._orig_file_path, every_step=every_step)
        refactored_blocks = {}
        for content_key in tmp_blocks:
            block = tmp_blocks[content_key]
            if (content_key == "TITLE"):
                block = blocks.TITLE(block)
                content_key = "TITLE"
            if (content_key == "ENEVERSION"):
                block = blocks.ENEVERSION(content=block)
            if (content_key == "TIMESTEP"):
                content_key = "TRAJ"
                block = blocks.TRAJ(list(block.values()))
            refactored_blocks.update({content_key: block})

        return refactored_blocks

    def append_tre(self, in_tre_obj):
        self.TRAJ.extend(in_tre_obj.TRAJ)
        self.clean_timestep_numbering()

    def get_dt(self):
        step1 = getattr(self, "TRAJ")[0]['time']
        step2 = getattr(self, "TRAJ")[1]['time']
        return step2 - step1

    def get_dsteps(self):
        step1 = getattr(self, "TRAJ")[0]['steps']
        step2 = getattr(self, "TRAJ")[1]['steps']
        return step2 - step1

    def concat(self, tres:(object, List[object]), renumber_timesteps:bool=False, newTitle:str=None):
        if(not isinstance(tres, List)):
            tres = [tres]

        for tre in tres:
            self.TRAJ.extend(tre.TRAJ)

        if(renumber_timesteps):
            self.clean_timestep_numbering()

        if(isinstance(newTitle, str)):
            self.TITLE = blocks.TITLE(newTitle)

    def clean_timestep_numbering(self, dsteps=False, dt=False):
        """clean_timestep_numbering
            Assumes dT is equal in each file!
        Parameters
        ----------
        dsteps :
        dt :

        Returns
        -------

        """

        if (not dsteps and not dt):
            dt = self.get_dt()
            dsteps = self.get_dsteps()
            traj = getattr(self, "TRAJ")
            timestep_list = []

            first = True
            zero_time = 0
            zero_step = 0
            step = 0
            for timestep in traj:
                if (first):
                    zero_step = timestep["steps"]
                    zero_time = timestep["time"]
                    first = False
                else:
                    new_step = step * dsteps + zero_step
                    timestep.update({"steps": new_step})
                    new_time = step * dt + zero_time
                    timestep.update({"time": new_time})
                step += 1
                timestep_list.append(timestep)
            setattr(traj, "timestep", timestep_list)
            setattr(self, "TRAJ", traj)
        else:
            raise NotImplementedError("Please give dsteps and dt to concat files. in clean_timestep_numbering")

    def write(self, out_path)->str:
        file = open(out_path, "w")

        # write:
        ##TITLE
        key = "TITLE"
        file.write(str(self.__getattribute__(key)))

        ##ENEVERS
        key = "ENEVERSION"
        file.write(str(self.__getattribute__(key)))

        ##TIMESTEPS
        timesteps = self.__getattribute__("TRAJ")
        for timestep in timesteps:
            key = "TIMESTEP"
            file.write(key + "\n")
            file.write("\t\t" + str(timestep["steps"]) + "\t" + str(timestep["time"]) + "\n")
            file.write("END\n")

            for y in sorted(timestep):
                if (y != "steps" and y != "time"):
                    key = y
                    file.write(key + "\n")
                    if (type(timestep[key]) is list):
                        file.write("".join(timestep[key]))
                    else:
                        file.write(timestep[key])
                    file.write("END\n")
            continue

        file.close()
        return out_path
