from typing import List
from pygromos.data import topology_templates
from pygromos.files.forcefield._generic_force_field import _generic_force_field
from pygromos.files.topology.top import Top
from pygromos.files.topology.ifp import Ifp
from pygromos.files.topology.mtb import Mtb
from pygromos.data.ff.Gromos2016H66 import ifp, mtb, mtb_orga
from pygromos.data.ff.Gromos54A7 import ifp as ifp_54a7, mtb as mtb_54a7


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
        else:
            self.ifp = Ifp(ifp_54a7)
            self.mtb = Mtb(mtb_54a7)
            self.mtb_orga = None

    def create_top(self, mol: str, in_top: Top = None) -> Top:
        if in_top is None:
            in_top = Top(in_value=topology_templates.topology_template_dir + "/blank_template+spc.top")

        # check if mol is a residue or residue list
        known_res_names = self.mtb.all_res_names
        if self.mtb_orga is not None:
            known_res_names += self.mtb_orga.all_res_names
        if isinstance(mol, str):
            if mol in known_res_names:
                return self.convert_residue_to_top(in_top, mol)
            else:
                raise ValueError(
                    f"{mol} is not a valid residue name"
                )  # TODO: implement RDKit mol to residue name conversion
        elif isinstance(mol, list):
            if all([x in known_res_names for x in mol]):
                for res in mol:
                    in_top = self.convert_residue_to_top(in_top, res)
                return in_top
            else:
                raise ValueError(
                    f"{mol} is not a valid residue name"
                )  # TODO: implement RDKit mol to residue name conversion
        return in_top

    def convert_residue_to_top(self, in_top: Top, residue: str) -> Top:
        # find the mtb block for the resname. and store the block in data
        data = None
        if residue in self.mtb.mtb_solutes.keys():
            data = self.mtb.mtb_solutes[residue]
        elif residue in self.mtb.mtb_solvents.keys():
            data = self.mtb.mtb_solvents[residue]
        elif residue in self.mtb.mtb_ends.keys():
            data = self.mtb.mtb_ends[residue]
        elif self.mtb_orga is not None:
            if residue in self.mtb.mtb_solutes.keys():
                data = self.mtb.mtb_solutes[residue]
            elif residue in self.mtb.mtb_solvents.keys():
                data = self.mtb.mtb_solvents[residue]
            elif residue in self.mtb.mtb_ends.keys():
                data = self.mtb.mtb_ends[residue]
        else:
            raise ValueError(f"{residue} is not a valid residue name")

        # create the topology by going through all fields in data
        for atom in data.atoms:
            ATNM = 0
            MRES = 0
            PANM = ""
            IAC = 0
            MASS = 0
            CG = 0
            CGC = 0
            INE = []
            INE14 = []
            C6 = 0
            C12 = 0
            IACname = ""

            in_top.add_new_atom(
                ATNM=ATNM,
                MRES=MRES,
                PANM=PANM,
                IAC=IAC,
                MASS=MASS,
                CG=CG,
                CGC=CGC,
                INE=INE,
                INE14=INE14,
                C6=C6,
                C12=C12,
                IACname=IACname,
            )
        return in_top
