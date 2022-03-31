import numpy as np
from pygromos.utils.typing import List, Dict
from pygromos.files.trajectory.blocks import energy_trajectory_subblock as ene_sub_block


class _general_pandas_trajectory_block:
    def __init__(self, content: List):
        self.content = content

    def __eq__(self, __o: object) -> bool:
        return (self.content == __o.content) and (self.__class__ == __o.__class__)

    def to_dict(self) -> Dict:
        return {"default_block": self.content}  # Only use for Debuging


class _general_pandas_tre_block(_general_pandas_trajectory_block):
    def __init__(self, content: List):
        super().__init__(content)
        # reset the static variables in the energy trajectory subsubblock
        ene_sub_block._general_pandas_energy_trajectory_subblock_numerated.num_energy_baths = 0
        ene_sub_block._general_pandas_energy_trajectory_subblock_numerated.num_force_groups = 0
        ene_sub_block._general_pandas_energy_trajectory_subblock_numerated.num_num_states = 0
        ene_sub_block._general_pandas_energy_trajectory_subblock_numerated.num_temp_groups = 0
        ene_sub_block._general_pandas_energy_trajectory_subblock_numerated.num_lambdas = 0

    def to_dict(self) -> Dict:
        return_dict = {}
        subblock_content = {}
        subblock_name = "dummy"

        for line in self.content:
            if line.strip() == "":
                continue
            if line.startswith("# "):
                subblock_name = line.strip().replace("# ", "")
                subblock_content.update({subblock_name: []})
            else:
                subblock_content[subblock_name].append(line.strip())

        tmp_sub2 = None
        tmp_sub_block = None
        for subblocktitle, subblock in subblock_content.items():
            if "lambda" == subblocktitle:  # Exceptions for TRG
                tmp_sub_block = getattr(ene_sub_block, "lam")(
                    subblock
                )  # stupid block naming.... but they couldn't have known
            elif "ABdih" == subblocktitle:
                subblocktitle = "precalclam"
                tmp_sub2.extend(subblock)
                tmp_sub_block = getattr(ene_sub_block, subblocktitle)(tmp_sub2)
            elif "numstates" == subblocktitle:
                subblocktitle = "eds"
                tmp_sub_block = getattr(ene_sub_block, subblocktitle)(subblock)
            elif "precalclam" == subblocktitle or "eds" == subblocktitle:
                continue
            elif "nr_lambdas" == subblocktitle:
                tmp_sub2 = subblock
                continue
            else:
                tmp_sub_block = getattr(ene_sub_block, subblocktitle)(subblock)

            return_dict.update(tmp_sub_block.to_dict())

        return return_dict


class TIMESTEP(_general_pandas_trajectory_block):
    def __init__(self, content: List):
        super().__init__(content)
        entries = self.content[0].strip().split()
        self.step = int(entries[0])
        self.time = float(entries[1])

    def to_dict(self) -> Dict:
        return {"step": self.step, "time": self.time}


"""
TRC Blocks:
"""


class POSITIONRED(_general_pandas_trajectory_block):
    def __init__(self, content: List):
        super().__init__(content)

    def to_dict(self) -> Dict:
        return_dict = {}
        iterator = 1
        for line in self.content:
            if line.strip().startswith("#") or line.strip() == "":
                continue
            return_dict.update({"POS_" + str(iterator): np.array(line.strip().split()).astype(np.float)})
            iterator += 1
        return return_dict


class SHAKEFAILPOSITION(_general_pandas_trajectory_block):
    def __init__(self, content: List):
        super().__init__(content)

    def to_dict(self) -> Dict:
        return_dict = {}
        iterator = 1
        for line in self.content:
            if line.strip().startswith("#") or line.strip() == "":
                continue
            return_dict.update({"SHK_" + str(iterator): np.array(line.strip().split()).astype(np.float)})
            iterator += 1
        return return_dict


class SHAKEFAILPREVPOSITION(_general_pandas_trajectory_block):
    def __init__(self, content: List):
        super().__init__(content)

    def to_dict(self) -> Dict:
        return_dict = {}
        iterator = 1
        for line in self.content:
            if line.strip().startswith("#") or line.strip() == "":
                continue
            return_dict.update({"SHKP_" + str(iterator): np.array(line.strip().split()).astype(np.float)})
            iterator += 1
        return return_dict


class REFPOSITION(_general_pandas_trajectory_block):
    def __init__(self, content: List):
        super().__init__(content)

    def to_dict(self) -> Dict:
        return_dict = {}
        iterator = 1
        for line in self.content:
            if line.strip().startswith("#") or line.strip() == "":
                continue
            return_dict.update({"LAT_" + str(iterator): np.array(line.strip().split()).astype(np.float)})
            iterator += 1
        return return_dict


class GENBOX(_general_pandas_trajectory_block):
    def __init__(self, content: List):
        super().__init__(content)

    def to_dict(self) -> Dict:
        if len(self.content) < 4:
            raise IOError("GENBOX block incomplete in trajectory")
        return_dict = {}
        return_dict.update({"genbox": self.content[0]})
        return_dict.update({"length": np.array(self.content[1].strip().split()).astype(np.float)})
        return_dict.update({"angles": np.array(self.content[2].strip().split()).astype(np.float)})
        return_dict.update({"watEver": np.array(self.content[3].strip().split()).astype(np.float)})

        return return_dict


"""
TRE Blocks:
"""


class ENERGY03(_general_pandas_tre_block):
    def __init__(self, content: List):
        super().__init__(content)

    def to_dict(self) -> Dict:
        return super().to_dict()


class VOLUMEPRESSURE03(_general_pandas_tre_block):
    def __init__(self, content: List):
        super().__init__(content)

    def to_dict(self) -> Dict:
        return super().to_dict()


"""
TRG Blocks:
"""


class FREEENERDERIVS03(_general_pandas_tre_block):
    def __init__(self, content: List):
        super().__init__(content)

    def to_dict(self) -> Dict:
        return super().to_dict()
