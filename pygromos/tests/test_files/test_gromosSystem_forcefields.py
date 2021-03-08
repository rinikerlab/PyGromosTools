import unittest, tempfile, importlib

from pygromos.files.gromos_system.gromos_system import Gromos_System
from pygromos.files.gromos_system.ff.forcefield_system import forcefield_system

from pygromos.tests.test_files import out_test_root_dir
tmp_test_dir = tempfile.mkdtemp(dir=out_test_root_dir, prefix="gromSystem_ff_")

if (importlib.util.find_spec("openforcefield") != None):
    has_openff = True
else:
    has_openff = False


class test_gromos_system_forcefields(unittest.TestCase):
    file_class = Gromos_System
    verbose = True

    smiles="CO"
    ff=forcefield_system()
    ff.mol_name = "MTL"
    

    def test_construct_empty(self):
        grSys = self.file_class(work_folder=tmp_test_dir, system_name="Testing1", Forcefield=self.ff)
        print(grSys)

    def test_construct_top_from_ff(self):
        grSys = self.file_class(work_folder=tmp_test_dir, system_name="Testing1", Forcefield=self.ff, in_smiles=self.smiles, auto_convert=True)
        print(grSys)

class test_gromos_system_54A7(test_gromos_system_forcefields):
    ff = forcefield_system(name="54A7")
    ff.mol_name = "MTL"

class test_gromos_system_2016H66(test_gromos_system_forcefields):
    ff = forcefield_system(name="2016H66")
    ff.mol_name = "MTL"

if has_openff:
    class test_gromos_system_openforcefield(test_gromos_system_forcefields):
        ff = forcefield_system(name="off")

"""
#TODO: include as soon as serenityff is in production state
class test_gromos_system_serenityforcefield(test_gromos_system_forcefields):
    ff = forcefield_system(name="serenityff")
"""
