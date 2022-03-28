import numpy as np
from pygromos.utils.typing import List, Dict


class _general_pandas_energy_trajectory_subblock:
    def __init__(self, content: List):
        self.content = []
        for i in content:
            if i.startswith("#"):  # cometimes there are inconsistent comment lines...(
                continue
            try:
                self.content.append(i.split())
            except ValueError:
                self.content.append(i)
        self.blockName = "default_block"
        self.list = []
        for i in self.content:
            for j in i:
                self.list.append(j)

    def to_dict(self) -> Dict:
        return {self.blockName: np.array(self.list, dtype=np.float64)}


class _general_pandas_energy_trajectory_subblock_numerated(_general_pandas_energy_trajectory_subblock):
    # first define some static variables for the stupid gromos format
    num_energy_baths = 0
    num_force_groups = 0
    num_nbforce_groups = 0
    num_num_states = 0
    num_temp_groups = 0
    num_lambdas = 0

    def __init__(self, content: List, subsubblock_number_code: str = "energy"):
        super().__init__(content)
        self.matrixList = []
        self.extra_line_for_num = 0

        if subsubblock_number_code == "energy":
            num_subsubblocks = _general_pandas_energy_trajectory_subblock_numerated.num_energy_baths
        elif subsubblock_number_code == "force":
            num_subsubblocks = _general_pandas_energy_trajectory_subblock_numerated.num_force_groups
        elif subsubblock_number_code == "nbforce":
            num_subsubblocks = _general_pandas_energy_trajectory_subblock_numerated.num_nbforce_groups
        elif subsubblock_number_code == "states":
            num_subsubblocks = _general_pandas_energy_trajectory_subblock_numerated.num_num_states
        elif subsubblock_number_code == "lambda":
            num_subsubblocks = _general_pandas_energy_trajectory_subblock_numerated.num_lambdas
        elif subsubblock_number_code == "temp":
            num_subsubblocks = _general_pandas_energy_trajectory_subblock_numerated.num_temp_groups
        else:
            raise KeyError("subblock number coder does not corresponde to variable")

        # get the number of baths
        if num_subsubblocks == 0:
            try:
                num_subsubblocks = int(self.content[0][0])
            except ValueError:
                print("no good value found for the number of subsubblocks")
            if num_subsubblocks <= 0:
                raise KeyError("invalid iterator with value: " + str(num_subsubblocks))

            # set value to static varibale for all futer values
            if subsubblock_number_code == "energy":
                _general_pandas_energy_trajectory_subblock_numerated.num_energy_baths = num_subsubblocks
            elif subsubblock_number_code == "force":
                _general_pandas_energy_trajectory_subblock_numerated.num_force_groups = num_subsubblocks
                _general_pandas_energy_trajectory_subblock_numerated.num_nbforce_groups = sum(
                    [x for x in range(1, num_subsubblocks + 1)]
                )
            elif subsubblock_number_code == "states":
                _general_pandas_energy_trajectory_subblock_numerated.num_num_states = num_subsubblocks
            elif subsubblock_number_code == "lambda":
                _general_pandas_energy_trajectory_subblock_numerated.num_lambdas = num_subsubblocks
            elif subsubblock_number_code == "temp":
                _general_pandas_energy_trajectory_subblock_numerated.num_temp_groups = num_subsubblocks

            # print("subblock number set to: " + str(num_subsubblocks) + "  for:  " + subsubblock_number_code)
            self.extra_line_for_num += 1
        # loop over all energy baths, import the data and format to numpy matrix
        for itr_line in range(num_subsubblocks):
            new_line = self.content[itr_line + self.extra_line_for_num]
            self.matrixList.append(new_line)

        self.matrix = np.array(self.matrixList, dtype=np.float64)

    def to_dict(self) -> Dict:
        return {self.blockName: self.matrix}


class totals(_general_pandas_energy_trajectory_subblock):
    def __init__(self, content: List):
        super().__init__(content)
        self.blockName = "totals"

    def to_dict(self) -> Dict:
        return super().to_dict()


class baths(_general_pandas_energy_trajectory_subblock_numerated):
    def __init__(self, content: List):
        super().__init__(content, subsubblock_number_code="energy")
        self.blockName = "baths"

    def to_dict(self) -> Dict:
        return super().to_dict()


class bonded(_general_pandas_energy_trajectory_subblock_numerated):
    def __init__(self, content: List):
        super().__init__(content, subsubblock_number_code="force")
        self.blockName = "bonded"

    def to_dict(self) -> Dict:
        return super().to_dict()


class nonbonded(_general_pandas_energy_trajectory_subblock_numerated):
    def __init__(self, content: List):
        super().__init__(content, subsubblock_number_code="nbforce")
        self.blockName = "nonbonded"

    def to_dict(self) -> Dict:
        return super().to_dict()


class special(_general_pandas_energy_trajectory_subblock_numerated):
    def __init__(self, content: List):
        super().__init__(content, subsubblock_number_code="force")
        self.blockName = "special"

    def to_dict(self) -> Dict:
        return super().to_dict()


class eds(_general_pandas_energy_trajectory_subblock):
    def __init__(self, content: List):
        super().__init__(content)
        self.blockName = "eds"

    def to_dict(self) -> Dict:
        return super().to_dict()


class precalclam(_general_pandas_energy_trajectory_subblock):
    def __init__(self, content: List):
        super().__init__(content)
        self.blockName = "precalclam"

    def to_dict(self) -> Dict:
        return super().to_dict()


class mass(_general_pandas_energy_trajectory_subblock):
    def __init__(self, content: List):
        super().__init__(content)
        self.blockName = "mass"

    def to_dict(self) -> Dict:
        return super().to_dict()


class temperature(_general_pandas_energy_trajectory_subblock_numerated):
    def __init__(self, content: List):
        super().__init__(content, subsubblock_number_code="temp")
        self.blockName = "temperature"

    def to_dict(self) -> Dict:
        return super().to_dict()


class volume(_general_pandas_energy_trajectory_subblock):
    def __init__(self, content: List):
        super().__init__(content)
        self.blockName = "volume"

    def to_dict(self) -> Dict:
        return super().to_dict()


class pressure(_general_pandas_energy_trajectory_subblock):
    def __init__(self, content: List):
        super().__init__(content)
        self.blockName = "pressure"

    def to_dict(self) -> Dict:
        return super().to_dict()


class xxx_copy_xxx(_general_pandas_energy_trajectory_subblock):
    def __init__(self, content: List):
        super().__init__(content)
        self.blockName = "xxx"

    def to_dict(self) -> Dict:
        return super().to_dict()


"""
TRG specials:
"""


class lam(_general_pandas_energy_trajectory_subblock):
    def __init__(self, content: List):
        super().__init__(content)
        self.blockName = "lambda"

    def to_dict(self) -> Dict:
        return super().to_dict()
