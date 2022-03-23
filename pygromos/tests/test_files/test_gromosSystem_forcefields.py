import unittest
import tempfile
import importlib

from pygromos.files.gromos_system.gromos_system import Gromos_System
from pygromos.files.forcefield.gromos.gromosff import GromosFF
from pygromos.files.forcefield.openff.openff import OpenFF

from pygromos.tests.test_files import out_test_root_dir

tmp_test_dir = tempfile.mkdtemp(dir=out_test_root_dir, prefix="gromSystem_ff_")

if importlib.util.find_spec("openforcefield") is not None:
    has_openff = True
else:
    has_openff = False


class test_gromos_system_forcefields(unittest.TestCase):
    file_class = Gromos_System
    verbose = True

    smiles = "CO"
    ff = GromosFF()
    top_residue_list = ["MTL"]

    def test_construct_empty(self):
        grSys = self.file_class(work_folder=tmp_test_dir, system_name="Testing1", forcefield=self.ff)
        print(grSys)

    def test_construct_top_from_ff(self):
        grSys = self.file_class(
            work_folder=tmp_test_dir,
            system_name="Testing1",
            forcefield=self.ff,
            in_smiles=self.smiles,
            auto_convert=True,
            in_residue_list=self.top_residue_list,
        )
        print(grSys)


class test_gromos_system_54A7(test_gromos_system_forcefields):
    ff = GromosFF(name="54A7")
    top_residue_list = ["CH3OH"]


class test_gromos_system_2016H66(test_gromos_system_forcefields):
    ff = GromosFF(name="2016H66")
    top_residue_list = ["MTL"]


if has_openff:

    class test_gromos_system_openforcefield(test_gromos_system_forcefields):
        ff = OpenFF()


"""
#TODO: include as soon as serenityff is in production state
class test_gromos_system_serenityforcefield(test_gromos_system_forcefields):
    ff = forcefield_system(name="serenityff")
"""
