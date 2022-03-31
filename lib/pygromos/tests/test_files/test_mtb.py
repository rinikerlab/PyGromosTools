import tempfile

from pygromos.files.topology.mtb import Mtb
from pygromos.data.ff.Gromos2016H66 import mtb, mtb_orga

from pygromos.data.ff.Gromos54A7 import mtb as mtb_g54a7
from pygromos.tests.test_files.general_file_functions import general_file_tests

from pygromos.tests.test_files import out_test_root_dir

root_out = tempfile.mkdtemp(dir=out_test_root_dir, prefix="mtb_")
out_path = root_out + "/out_mtb.mtb"


class test_mtb(general_file_tests):
    __test__ = True
    class_type = Mtb
    in_file_path = mtb
    root_out = root_out

    def test_parsing_test_file(self):
        mtb_file = self.class_type(self.in_file_path)
        assert isinstance(mtb_file, self.class_type)
        return 0

    def test_write(self):
        mtb_file = self.class_type(self.in_file_path)
        mtb_file.write(out_path)
        return 0


class test_mtb_g54a7(test_mtb):
    __test__ = True
    in_file_path = mtb_g54a7


class test_mtb_orga(test_mtb):
    __test__ = True
    in_file_path = mtb_orga

    def test_all_mtb_solutes_read(self):
        mtb_file = self.class_type(self.in_file_path)
        assert len(mtb_file.mtb_solutes) == 63  # 63 solutes in 2016H66 orga
        return 0

    def test_CHE_read(self):
        mtb_file = self.class_type(self.in_file_path)
        mol = mtb_file.mtb_solutes["CHE"]

        assert mol.RNME == "CHE"

        assert mol.NMAT == 6
        assert mol.NB == 6
        assert mol.NBA == 6
        assert mol.NIDA == 0
        assert mol.NDA == 6
        assert mol.NEX == 0

        atom = mol.atoms[0]
        assert atom.ATOM == 1
        assert atom.ANM == "C1"
        assert atom.IACM == 18
        assert atom.MASS == 4
        assert atom.CGMI == 0
        assert atom.CGM == 0
        assert atom.MAE == 4
        assert atom.MSAE == [2, 3, 5, 6]

        bond = mol.bonds[0]
        assert bond.IB == 1
        assert bond.JB == 2
        assert bond.MCB == 27

        angle = mol.angles[0]
        assert angle.IB == 1
        assert angle.JB == 2
        assert angle.KB == 3
        assert angle.MCB == 8

        dihedral = mol.dihedrals[0]
        assert dihedral.IB == 1
        assert dihedral.JB == 2
        assert dihedral.KB == 3
        assert dihedral.LB == 4
        assert dihedral.MCB == 34

        return 0
