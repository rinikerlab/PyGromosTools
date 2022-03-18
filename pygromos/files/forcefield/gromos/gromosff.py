from typing import List

# from pygromos.data import topology_templates
from pygromos.data.ff.Gromos2016H66 import ifp, mtb, mtb_orga
from pygromos.data.ff.Gromos54A7 import ifp as ifp_54a7, mtb as mtb_54a7

from pygromos.files.forcefield._generic_force_field import _generic_force_field
from pygromos.files.topology.top import Top
from pygromos.files.coord.cnf import Cnf
from pygromos.files.topology.ifp import Ifp
from pygromos.files.topology.mtb import Mtb

from pygromos.gromos.gromosPP import GromosPP


class GromosFF(_generic_force_field):
    def __init__(
        self, name: str = "gromos", path_to_files: List(str) = None, auto_import: bool = True, verbose: bool = False
    ):
        super().__init__(name, path_to_files, auto_import, verbose)

    def auto_import_ff(self):
        if self.name in ["2016H66", "gromos2016H66", "gromos"]:
            self.ifp = Ifp(ifp)
            self.mtb = Mtb(mtb)
            self.mtb_orga = Mtb(mtb_orga)
            self.in_building_block_lib_path = [mtb, mtb_orga]
            self.in_parameter_lib_path = [ifp]
        else:
            self.ifp = Ifp(ifp_54a7)
            self.mtb = Mtb(mtb_54a7)
            self.mtb_orga = None
            self.in_building_block_lib_path = [mtb_54a7]
            self.in_parameter_lib_path = [ifp_54a7]

    def create_top(self, mol: str, in_top: Top = None, **kwargs) -> Top:
        if "residue_list" in kwargs:
            residue_list = kwargs["residue_list"]
            work_dir = kwargs["work_dir"]
            GromosPP.make_top(
                out_top_path=work_dir + "/temp.top",
                in_building_block_lib_path=self.in_building_block_lib_path,
                in_parameter_lib_path=self.in_parameter_lib_path,
                in_sequence=residue_list,
            )
        else:
            self.creteTopoFromSmiles(mol, in_top, **kwargs)

    def create_cnf(self, mol: str, in_cnf: Cnf = None, **kwargs) -> Cnf:
        return Cnf(None)  # TODO: implement

    def creteTopoFromSmiles(self, mol: str, in_top: Top = None, **kwargs) -> Top:
        raise NotImplementedError("WIP")
