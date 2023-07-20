"""
FUNCTIONLIB:            gromos++ input file functions
Description:
    in this lib, gromosXX input file mainpulating functions are gathered
Author: Benjamin Schroeder
"""

import numpy as np
import pandas as pd
from pygromos.files._basics import _general_gromos_file
from pygromos.files.blocks import TITLE


class NOE(_general_gromos_file._general_gromos_file):

    _orig_file_path: str
    _required_blocks = ["TITLE", "AVERAGE_NOE", "NOE_VIOLATIONS", "RESTRAINT_LEGEND"]
    _gromos_file_ending = "noe"
    # POSSIBLE GROMOS BLOCKS
    TITLE: TITLE
    AVERAGE_NOE: pd.DataFrame
    NOE_VIOLATIONS: pd.DataFrame
    RESTRAINT_LEGEND: pd.DataFrame

    def __init__(self, in_value: str):
        super().__init__(in_value=in_value)

    def __str__(self):
        text = ""
        text += self.__getattribute__("TITLE").block_to_string()
        for block in sorted(self.get_block_names()):
            if block == "TITLE" or isinstance(block, type(None)):
                continue
            text += str(self.__getattribute__(block))
        return text

    def read_file(self):
        in_file = open(self._orig_file_path, "r")
        in_file_lines = in_file.readlines()

        # read_blocks
        known_blocks = ["TITLE", "NOE VIOLATIONS", "AVERAGE NOE"]
        in_block = False
        tmp_list = []
        state = 0
        noe_states = {state: {}}

        for line in in_file_lines:
            if any([line.startswith(block) for block in known_blocks]):
                in_block = line.strip()
                if in_block in noe_states[state].keys():
                    state = state + 1
                    noe_states.update({state: {}})
            elif in_block and line.startswith("END"):
                noe_states[state].update({in_block: tmp_list})
                tmp_list = []
                in_block = False
            elif in_block:
                tmp_list.append(line)

        # Blockify
        # TITLE
        title = noe_states[0]["TITLE"]

        # NOE - VIOLATIONS
        # read header info
        # NOE Restraint
        violations_header = [s.replace("#", "") for s in noe_states[0]["NOE VIOLATIONS"] if (s.startswith("#"))]

        t_legend = []
        for line in violations_header:
            if "\n" == line:
                break
            legend_entry = [int(line.split()[0])] + line.split()[1:5]
            if len(legend_entry) < 5:
                raise ValueError("Too few columns in line: " + line)
            t_legend.append(legend_entry)

        restraint_legend = pd.DataFrame(t_legend, columns=["Nr.", "resI", "atomI", "resJ", "atomJ"])

        # data Header
        column_header_nviols = [x for x in violations_header[-1].strip().split("  ") if (x != "")]

        header = [s.replace("#", "") for s in noe_states[0]["AVERAGE NOE"] if (s.startswith("#"))]
        column_header_avgnoe = [x for x in header[-1].strip().split("  ") if (x != "")]

        # DATA
        # Remove header Info & pandafy data
        NOE_violations = []
        average_NOE = []

        for state, noe_data in noe_states.items():
            noe_viol = noe_data["NOE VIOLATIONS"]
            av_noe = noe_data["AVERAGE NOE"]
            viol_data = [[float(x) for x in line.split()] for line in noe_viol if not line.startswith("#")]
            nvi = pd.DataFrame(viol_data, columns=column_header_nviols)
            nvi["state"] = state

            av_data = [[float(x) for x in line.split()] for line in av_noe if not line.startswith("#")]
            av_n = pd.DataFrame(av_data, columns=column_header_avgnoe)
            av_n["state"] = state

            NOE_violations.append(pd.DataFrame(nvi, columns=column_header_nviols + ["state"]))
            average_NOE.append(pd.DataFrame(av_n, columns=column_header_avgnoe + ["state"]))

        return {
            "TITLE": title,
            "AVERAGE_NOE": pd.concat(average_NOE, ignore_index=True),
            "NOE_VIOLATIONS": pd.concat(NOE_violations, ignore_index=True),
            "RESTRAINT_LEGEND": restraint_legend,
        }


class JVAL(_general_gromos_file._general_gromos_file):

    _orig_file_path: str

    # POSSIBLE GROMOS BLOCKS
    content: pd.DataFrame
    _data_header: dict = {
        "num": int,
        "mol": int,
        "resI": int,
        "resN": str,
        "atomNI": str,
        "atomNJ": str,
        "atomNK": str,
        "atomNL": str,
        "atomIdI": str,
        "atomIdJ": str,
        "atomIdK": str,
        "atomIdL": str,
        "A": float,
        "B": float,
        "C": float,
        "delta": float,
        "J_0": float,
        "phi_ave": float,
        "phi_rmsd": float,
        "J_ave": float,
        "J_rmsd": float,
        "|Jave-J0|": float,
    }

    def __init__(self, in_value: str):
        super().__init__(in_value=in_value)

        self.average_J_deviation = np.mean(self.content["J_ave"] - self.content["J_0"])
        self.average_abs_deviation = np.mean(np.abs(self.content["J_ave"] - self.content["J_0"]))
        self.root_mean_square_deviation = np.sqrt(np.mean(np.square(self.content["J_ave"] - self.content["J_0"])))
        self.root_mean_square_deviation_over_deviations = np.sqrt(np.mean(np.square(self.content["J_rmsd"])))

    def __str__(self):
        text = ""
        text += "#\taverage_J_deviation\t" + str(self.average_J_deviation) + "\n"
        text += "#\taverage_abs_deviation\t" + str(self.average_abs_deviation) + "\n"
        text += "#\troot_mean_square_deviation\t" + str(self.root_mean_square_deviation) + "\n"
        text += (
            "#\troot_mean_square_deviation_over_deviations\t"
            + str(self.root_mean_square_deviation_over_deviations)
            + "\n"
        )
        text += self.content.to_string() + "\n"
        return text

    def write(self, out_path: str) -> str:
        text = ""
        text += "#\taverage_J_deviation\t" + str(self.average_J_deviation) + "\n"
        text += "#\taverage_abs_deviation\t" + str(self.average_abs_deviation) + "\n"
        text += "#\troot_mean_square_deviation\t" + str(self.root_mean_square_deviation) + "\n"
        text += (
            "#\troot_mean_square_deviation_over_deviations\t"
            + str(self.root_mean_square_deviation_over_deviations)
            + "\n"
        )
        self.content.to_csv(out_path, header=text)
        return out_path

    def __getitem__(self, item: int):
        return self.content[item]

    def __len__(self):
        return len(self.content)

    def __iter__(self):
        return iter(self.content)

    def read_file(self):

        jval_file = open(self._orig_file_path, "r")
        jval_lines = jval_file.readlines()
        jval_file.close()

        in_block = False
        data_lines = []
        for line in jval_lines:
            if line.startswith("# num mol"):
                in_block = True
            elif line.startswith("#") or line == "\n":
                continue
            elif in_block:
                data_lines.append(line.strip().split())

        data = pd.DataFrame(data=data_lines, columns=self._data_header)
        data = data.astype(self._data_header)
        return {"content": data}
