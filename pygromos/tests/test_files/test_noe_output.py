import unittest
import numpy as np
import pandas as pd
from pygromos.files.otherfiles import noe_output

from pygromos.tests.in_testfiles import in_test_file_path


class test_noe_output(unittest.TestCase):
    class_name: type = noe_output.NOE
    in_file_path = in_test_file_path + "/noe_output/NOE_IR.out"
    in_file_path_duplicated = in_test_file_path + "/noe_output/NOE_IR_2_FRAMES.out"

    def test_constructor_empty(self):
        t = self.class_name(None)
        assert isinstance(t, self.class_name)

    def test_constructor_noe_output_file(self):
        t = self.class_name(self.in_file_path)
        assert isinstance(t, self.class_name)

    def test_average_noe(self):
        t = self.class_name(self.in_file_path)
        content = t.read_file()
        avg_noe = content["AVERAGE_NOE"]
        assert avg_noe.shape == (179, 11)
        print(avg_noe.iloc[-1])
        assert np.all(avg_noe.iloc[0].values == (1, 0.321, 0.296, 0.275, 0.052, 0.074, 0.075, 0.010, 0.014, 0.013, 0))
        assert np.all(
            avg_noe.iloc[-1].values == (179, 0.469, 0.465, 0.462, 0.030, 0.031, 0.032, 0.001, 0.001, 0.001, 0)
        )

    def test_noe_violations(self):
        t = self.class_name(self.in_file_path)
        content = t.read_file()
        avg_noe = content["NOE_VIOLATIONS"]
        assert avg_noe.shape == (179, 6)
        print(avg_noe.iloc[-1])
        assert np.all(avg_noe.iloc[0].values == (1, 0.350, -0.029, -0.054, -0.075, 0))
        assert np.all(avg_noe.iloc[-1].values == (179, 0.720, -0.251, -0.255, -0.258, 0))

    def test_average_noe(self):
        t = self.class_name(self.in_file_path_duplicated)
        content = t.read_file()
        avg_noe = content["AVERAGE_NOE"]
        assert avg_noe.shape == (2 * 179, 11)
        assert avg_noe.iloc[-1].name == 2 * 179 - 1  # test that index is continous

    def test_noe_violations_2_frames(self):
        t = self.class_name(self.in_file_path_duplicated)
        content = t.read_file()
        avg_noe = content["NOE_VIOLATIONS"]
        assert avg_noe.shape == (2 * 179, 6)
        assert avg_noe.iloc[-1].name == 2 * 179 - 1

    def test_noe_restraint_legend(self):
        expected = []
        with open(self.in_file_path) as f:
            for _ in range(9):
                next(f)
            for _ in range(179):
                elements = next(f).split()[1:6]
                elements[0] = int(elements[0])
                expected.append(elements)
        t = self.class_name(self.in_file_path)
        content = t.read_file()
        output = content["RESTRAINT_LEGEND"]
        assert np.all(output == pd.DataFrame(expected, columns=["Nr.", "resI", "atomI", "resJ", "atomJ"]))
