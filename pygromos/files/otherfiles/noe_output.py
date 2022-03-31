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
        header = [s.replace("#", "") for s in noe_states[0]["NOE VIOLATIONS"] if (s.startswith("#"))]

        t_legend = []
        for i, line in enumerate(header):
            if "\n" == line:
                break
            else:
                t_legend.append(line)

        t_legend = [
            [int(x.strip().split(" ")[0])] + [y for y in x.strip().split(" ")[1:] if (y != "")] for x in t_legend
        ]
        for legend_entry in t_legend:
            if len(legend_entry) == 5:
                continue
            elif len(legend_entry) > 5:
                t_ind = t_legend.index(legend_entry)
                t_legend.remove(legend_entry)
                t_legend.insert(t_ind, legend_entry[:5])

            else:
                raise Exception("WHAT SHALL I DO WITH THAT? " + str(legend_entry))
        restraint_legend = pd.DataFrame(t_legend, columns=["Nr.", "resI", "atomI", "resJ", "atomJ"])

        # data Header
        column_header_nviols = [x for x in header[-1].strip().split("  ") if (x != "")]

        header = [s.replace("#", "") for s in noe_states[0]["AVERAGE NOE"] if (s.startswith("#"))]
        column_header_avgnoe = [x for x in header[-1].strip().split("  ") if (x != "")]

        # DATA
        # Remove header Info & pandafy data
        NOE_violations = None
        average_NOE = None

        for state, noe_data in noe_states.items():
            noe_viol = noe_data["NOE VIOLATIONS"]
            av_noe = noe_data["AVERAGE NOE"]
            no_header = list(
                map(lambda x: map(float, x.strip().split()), filter(lambda y: not y.startswith("#"), noe_viol))
            )
            nvi = pd.DataFrame(no_header, columns=column_header_nviols)
            nvi["state"] = state

            no_header = list(
                map(lambda x: map(float, x.strip().split()), filter(lambda y: not y.startswith("#"), av_noe))
            )
            av_n = pd.DataFrame(no_header, columns=column_header_avgnoe)
            av_n["state"] = state

            if not isinstance(NOE_violations, pd.DataFrame):
                NOE_violations = pd.DataFrame(nvi, columns=column_header_nviols + ["state"])
            else:
                NOE_violations = NOE_violations.append(nvi, ignore_index=True)

            if not isinstance(average_NOE, pd.DataFrame):
                average_NOE = pd.DataFrame(av_n, columns=column_header_avgnoe + ["state"])
            else:
                average_NOE = NOE_violations.append(av_n, ignore_index=True)

        return {
            "TITLE": title,
            "AVERAGE_NOE": average_NOE,
            "NOE_VIOLATIONS": NOE_violations,
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
