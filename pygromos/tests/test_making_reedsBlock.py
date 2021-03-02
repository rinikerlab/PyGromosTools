import os
import unittest
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
        cnf_file = cnf.Cnf(in_file_renum_path)
        cnf_file.clean_posiResNums()
        residues = cnf_file.get_residues()

        self.assertEqual(2, len(residues), "Found  more/less than two residues!")
        self.assertEqual(2, list(residues["3U5L"].keys())[0], "Could not recognzize expected resNum for 3U5l = 2" )
        self.assertEqual(1, list(residues["3MXF"].keys())[0], "Could not recognzize expected resNum for 3MXF = 1" )
        self.assertEqual(31, residues["3U5L"][2], "Could not recognzize expected atom ammount for 3U5l = 2" )
        self.assertEqual(35, residues["3MXF"][1], "Could not recognzize expected atom ammount for 3MXF = 1" )

    def test_write_out(self):
        cnf_file = cnf.Cnf(in_file_path)
        cnf_file.write(out_path)

if __name__ == '__main__':
    unittest.main()


    from pygromos.files import imd

    in_path=os.path.dirname(__file__)+"/testfiles/imd/in_imd_REEDS1.imd"
    out_path=os.path.dirname(__file__)+"/testfiles/imd/out_imd_REEDS1.imd"
    svals = [ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.3162, 0.1, 0.0316
    , 0.0285, 0.0254, 0.0238, 0.02305, 0.0223, 0.022, 0.02185, 0.0217, 0.02155, 0.0214
    , 0.02125, 0.0211, 0.02095, 0.0208, 0.0207, 0.0206, 0.0205, 0.0204, 0.0203, 0.0202
    , 0.0201, 0.02, 0.0199, 0.0198, 0.0197, 0.0196, 0.0195, 0.0194, 0.0193, 0.0191
    , 0.0189, 0.0187, 0.0184, 0.0181, 0.0174, 0.0168, 0.0162, 0.015167, 0.014133, 0.0131
    , 0.012756, 0.012411, 0.012067, 0.011722, 0.011378, 0.011033, 0.010689, 0.010344, 0.01, 0.0066
    , 0.0032, 0.001, 0.0003]
    #svals = svals



    imd_file = imd.Imd(in_path)
    imd_file.edit_REEDS(SVALS=svals)
    imd_file.write(out_path)
    exit(0)
