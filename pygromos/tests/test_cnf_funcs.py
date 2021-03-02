import unittest, os
from pygromos.files import coord as cnf
in_file_path= os.path.dirname(__file__)+"/testfiles/cnf/in_cnf1.cnf"
in_file_renum_path = os.path.dirname(__file__)+"/testfiles/cnf/in_renumber_ligs.cnf"

out_path = os.path.dirname(__file__)+"/testfiles/cnf/out_cnf1.cnf"
out_red_path = os.path.dirname(__file__)+"/testfiles/cnf/out_cnf1_reduced.cnf"

class test_cnf(unittest.TestCase):
    def test_parse(self):
        cnf_file = cnf.Cnf(in_file_path)
        self.assertEqual(type(cnf_file), cnf.Cnf)

    def test_delete_res(self):
        cnf_file = cnf.Cnf(in_file_path)
        cnf_file.delete_residue(resID=4, resName="4E96")
        cnf_file.write(out_red_path)

    def test_clean_posiResNumsByName(self):
        cnf_file = cnf.Cnf(in_file_renum_path, clean_resiNumbers_by_Name=True)
        cnf_file.clean_posiResNums()
        residues = cnf_file.get_residues()

        #print(residues)
        self.assertEqual(2, len(residues), "Found  more/less than two residues!")
        self.assertEqual(2, list(residues["3U5L"].keys())[0], "Could not recognzize expected resNum for 3U5l = 2" )
        self.assertEqual(1, list(residues["3MXF"].keys())[0], "Could not recognzize expected resNum for 3MXF = 1" )
        self.assertEqual(31, residues["3U5L"][2], "Could not recognzize expected atom ammount for 3U5l = 2" )
        self.assertEqual(35, residues["3MXF"][1], "Could not recognzize expected atom ammount for 3MXF = 1" )

    def test_write_out(self):
        cnf_file = cnf.Cnf(in_file_path)
        cnf_file.write(out_path)
