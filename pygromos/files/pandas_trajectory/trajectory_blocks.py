import numpy as np
import pygromos.files.pandas_trajectory.energy_trajectory_subblock as ene_sub_block

class _general_pandas_trajectory_block():
    def __init__(self, content):
        self.content = content

    def to_dict(self)->dict:
        return {"default_block": self.content}  #Only use for Debuging

class _general_pandas_tre_block(_general_pandas_trajectory_block):
    def __init__(self, content):
        super().__init__(content)
        # reset the static variables in the energy trajectory subsubblock
        ene_sub_block._general_pandas_energy_trajectory_subblock_numerated.num_energy_baths = 0
        ene_sub_block._general_pandas_energy_trajectory_subblock_numerated.num_force_groups = 0
        ene_sub_block._general_pandas_energy_trajectory_subblock_numerated.num_num_states = 0
        ene_sub_block._general_pandas_energy_trajectory_subblock_numerated.num_temp_groups = 0
        ene_sub_block._general_pandas_energy_trajectory_subblock_numerated.num_lambdas = 0

    def to_dict(self) -> dict:
        return_dict = {}
        subblock_content = {}
        subblock_name = "dummy"

        for line in self.content:
            if line.strip() == "":
                continue
            if line.startswith("# "):
                subblock_name = line.strip().replace("# ","")
                subblock_content.update({subblock_name:[]})
            else:
                subblock_content[subblock_name].append(line.strip())
                

        for subblocktitle, subblock in subblock_content.items():
            tmp_sub_block = getattr(ene_sub_block, subblocktitle)(subblock)
            return_dict.update(tmp_sub_block.to_dict())
        
        return return_dict

class TIMESTEP(_general_pandas_trajectory_block):
    def __init__(self, content):
        super().__init__(content)
        entries = self.content[0].strip().split()
        self.step = entries[0]
        self.time = entries[1]

    def to_dict(self) -> dict:
        return {"TIMESTEP_step":self.step, "TIMESTEP_time":self.time}

class POSITIONRED(_general_pandas_trajectory_block):
    def __init__(self, content):
        super().__init__(content)
        
    def to_dict(self) -> dict:
        return_dict = {}
        iterator = 1
        for line in self.content:
            if line.strip().startswith('#') or line.strip() == '':
                continue
            return_dict.update({"POS_"+str(iterator):np.array(line.strip().split()).astype(np.float)})
            iterator += 1
        return return_dict

class ENERGY03(_general_pandas_tre_block):
    def __init__(self, content):
        super().__init__(content)

    def to_dict(self) -> dict:
        return super().to_dict()


class VOLUMEPRESSURE03(_general_pandas_tre_block):
    def __init__(self, content):
        super().__init__(content)

    def to_dict(self) -> dict:
        return super().to_dict()
    
                    

