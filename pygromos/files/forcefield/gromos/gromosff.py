import warnings

# from pygromos.data import topology_templates
from pygromos.data.ff.Gromos2016H66 import ifp, mtb, mtb_orga
from pygromos.data.ff.Gromos54A7 import ifp as ifp_54a7, mtb as mtb_54a7

from pygromos.files.forcefield._generic_force_field import _generic_force_field
from pygromos.files.topology.top import Top
from pygromos.files.coord.cnf import Cnf
from pygromos.files.topology.ifp import Ifp
from pygromos.files.topology.mtb import Mtb
from pygromos.files.blocks.coord_blocks import atomP

from pygromos.gromos.gromosPP import GromosPP


class GromosFF(_generic_force_field):
    def __init__(
        self, name: str = "gromos", path_to_files: str = None, auto_import: bool = True, verbose: bool = False
    ):
        super().__init__(name, path_to_files, auto_import, verbose)
        if auto_import:
            self.auto_import_ff()
        self.gromosPP = GromosPP()

    def auto_import_ff(self):
        if self.name in ["2016H66", "gromos2016H66", "gromos"]:
            self.ifp = Ifp(ifp)
            self.mtb = Mtb(mtb)
            self.mtb_orga = Mtb(mtb_orga)
            self.in_building_block_lib_path = [mtb, mtb_orga]
            self.in_parameter_lib_path = ifp
        else:
            self.ifp = Ifp(ifp_54a7)
            self.mtb = Mtb(mtb_54a7)
            self.mtb_orga = None
            self.in_building_block_lib_path = [mtb_54a7]
            self.in_parameter_lib_path = ifp_54a7

    def create_top(self, mol: str, in_top: Top = None, **kwargs) -> Top:
        if "residue_list" in kwargs:
            if self.gromosPP._found_binary["make_top"]:
                residue_list = kwargs["residue_list"]
                work_dir = kwargs["work_folder"]
                self.gromosPP.make_top(
                    out_top_path=work_dir + "/temp.top",
                    in_building_block_lib_path=" ".join(self.in_building_block_lib_path),
                    in_parameter_lib_path=self.in_parameter_lib_path,
                    in_sequence=" ".join(residue_list),
                )
                self.top = Top(work_dir + "/temp.top")
                return self.top
            else:
                warnings.warn("Could not create topology, gromos binaries not found")
                return Top(None)
        else:
            return self.creteTopoFromSmiles(mol, in_top, **kwargs)

    def create_cnf(self, mol: str, in_cnf: Cnf = None, **kwargs) -> Cnf:
        # prepare Cnf with rdkit
        cnf = Cnf(in_value=mol)
        # if in_cnf is None:
        #     cnf = Cnf(in_value=mol)
        # elif isinstance(in_cnf, Cnf):
        #     cnf = in_cnf + Cnf(in_value=mol)
        # elif isinstance(in_cnf, str):
        #     cnf = Cnf(in_value=in_cnf) + Cnf(in_value=mol)
        # else:
        #     raise TypeError("in_cnf must be a Cnf object!")

        # reformat for gromos
        if self.gromosPP._check_all_binaries():
            try:
                new_pos = [
                    atomP(
                        xp=atom.xp,
                        yp=atom.yp,
                        zp=atom.zp,
                        resID=atom.resID,
                        atomType=atom.atomType + str(i + 1),
                        atomID=atom.atomID,
                        resName=kwargs["residue_list"][0],
                    )
                    for i, atom in enumerate(cnf.POSITION)
                ]
                cnf.POSITION = new_pos
                cnf.write_pdb(kwargs["work_folder"] + "/tmp.pdb")
                self.gromosPP.pdb2gromos(
                    in_top_path=self.top.path,
                    in_pdb_path=kwargs["work_folder"] + "/tmp.pdb",
                    out_cnf_path=kwargs["work_folder"] + "/tmp.cnf",
                )
                self.gromosPP.add_hydrogens(
                    in_cnf_path=kwargs["work_folder"] + "/tmp.cnf",
                    out_cnf_path=kwargs["work_folder"] + "/tmp_H.cnf",
                    in_top_path=self.top.path,
                )
                cnf = Cnf(in_value=kwargs["work_folder"] + "/tmp_H.cnf")
            except IOError:
                warnings.warn("Could not convert cnf from rdkit to gromos, will use rdkit cnf")
        return cnf

    def creteTopoFromSmiles(self, mol: str, in_top: Top = None, **kwargs) -> Top:
        raise NotImplementedError("WIP")
