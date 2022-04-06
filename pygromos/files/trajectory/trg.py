"""
File:            Class for tre files in pandas
Description:
    The pandas trajectory TRE class offers a easy method to process GROMOS's .trg files in python
    The tre files are parsed into an easy to use pandas dataframe.

    This class should be a alternative for the data post processing with ene_ana in gromos++

Author:  Marc Thierry Lehner & Benjamin Ries

TODO: add stride option to all member functions
TODO: add support for periodic boundary condition

TODO: add ene_ana functions

"""

# imports
import pandas as pd
import numpy as np

import pygromos.files.trajectory._general_trajectory as traj


class gromos_2020_trg_block_names_table:
    totals_subblock_names = [
        "dHdl",
        "dKdl",
        "dVdl",
    ] + ["WIP" for x in range(40)]
    precalclam_subblock = [
        "nr_lambdas",
        "A_e_lj",
        "B_e_lj",
        "A_e_crf",
        "B_e_crf",
        "AB_kinetic",
        "AB_bond",
        "AB_angle",
        "AB_improper",
        "AB_disres",
        "AB_dihres",
        "AB_disfld",
    ]


class Trg(traj._General_Trajectory):
    _gromos_file_ending: str = "trg"

    def __init__(self, input_value: str or None, auto_save=True, stride: int = 1, skip: int = 0):
        super().__init__(input_value, auto_save=auto_save, stride=stride, skip=skip)
        self.block_name_table = gromos_2020_trg_block_names_table

    def get_totals(self) -> pd.DataFrame:
        self.totals = pd.DataFrame(
            data=np.stack(self.database["totals"].to_numpy()), columns=self.block_name_table.totals_subblock_names
        )
        return self.totals

    def get_lambdas(self) -> pd.DataFrame:
        self.totals = pd.DataFrame(data=np.stack(self.database["lambda"].to_numpy()), columns=["lambda"])
        return self.totals

    def get_precalclam(self) -> pd.DataFrame:

        groups = self.database["precalclam"].iloc[0][0]
        adapted_cols = [self.block_name_table.precalclam_subblock[0]]
        for i in range(1, int(groups) + 1):
            adapted_cols.extend(
                list(map(lambda x: str(x) + "_" + str(i), self.block_name_table.precalclam_subblock[1:]))
            )
        adapted_cols.extend(["A_dihedral", "B_dihedral"])

        self.totals = pd.DataFrame(data=np.stack(self.database["precalclam"].to_numpy()), columns=adapted_cols)
        return self.totals
