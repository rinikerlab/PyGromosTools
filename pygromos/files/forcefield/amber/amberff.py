import os
import shutil
from typing import List
from pygromos.files.coord.cnf import Cnf

from pygromos.files.forcefield._generic_force_field import _generic_force_field
from pygromos.files.topology.top import Top
from pygromos.files.forcefield.amber.amber2gromos import amber2gromos


class AmberFF(_generic_force_field):
    def __init__(
        self, name: str = "amber", path_to_files: List(str) = None, auto_import: bool = True, verbose: bool = False
    ):
        super().__init__(name, path_to_files=path_to_files, auto_import=auto_import, verbose=verbose)

    def auto_import_ff(self):
        # check path
        if self.path is not None:
            if isinstance(self.path, List) and len(self.path) > 0 and isinstance(self.path[0], str):
                self.amber_basedir = self.path[0]
        elif shutil.which("tleap") is not None:
            has_amber = True  # ambertools in path
            self.amber_basedir = os.path.abspath(os.path.dirname(shutil.which("tleap")) + "/../")
        else:
            has_amber = False
            raise ImportError(
                "Could not import GAFF FF as ambertools was missing! " "Please install the package for this feature!"
            )

        if self.verbose:
            print("Found amber: " + str(has_amber))

        self.amber_bindir = self.amber_basedir + "/bin"
        self.leaprc_files = [
            self.amber_basedir + "/dat/leap/cmd/leaprc.gaff",
            self.amber_basedir + "/dat/leap/cmd/leaprc.water.tip3p",
        ]
        self.frcmod_files = [self.amber_basedir + "/dat/leap/parm/frcmod.chcl3"]

        for leaprc in self.leaprc_files:
            if not os.path.isfile(leaprc):
                raise ImportError("could not find ff file " + leaprc)

        for frcmod in self.frcmod_files:
            if not os.path.isfile(frcmod):
                raise ImportError("could not find ff file " + frcmod)

    def create_top(self, mol: str, in_top: Top = None, in_mol2_file: str = None) -> Top:
        self.mol = mol
        self.in_mol2_file = in_mol2_file

        if self.amber is None:
            if in_mol2_file is None:
                self.create_mol2()
            self.amber = amber2gromos(
                in_mol2_file=self.in_mol2_file,
                mol=self.mol,
                forcefield=self.Forcefield,
                gromosPP=self.gromosPP,
                work_folder=self.work_folder,
            )
        if in_top is None:
            self.top = Top(self.amber.get_gromos_topology())
        elif isinstance(in_top, Top):
            self.top = in_top + Top(self.amber.get_gromos_topology())
        elif isinstance(in_top, str):
            self.top = Top(in_top) + Top(self.amber.get_gromos_topology())
        else:
            raise TypeError("in_top is of wrong type")

    def create_cnf(self, mol: str, in_cnf: Top = None) -> Cnf:
        if self.amber is None:
            self.create_mol2()
            self.amber = amber2gromos(
                in_mol2_file=self.in_mol2_file,
                mol=self.mol,
                forcefield=self.Forcefield,
                gromosPP=self.gromosPP,
                work_folder=self.work_folder,
            )
        if in_cnf is None:
            self.cnf = Cnf(self.amber.get_gromos_coordinate_file())
        elif isinstance(in_cnf, Cnf):
            self.cnf = in_cnf + Cnf(self.amber.get_gromos_coordinate_file())
        elif isinstance(in_cnf, str):
            self.cnf = Cnf(in_cnf) + Cnf(self.amber.get_gromos_coordinate_file())
        else:
            raise TypeError("in_cnf is of wrong type")

    def create_mol2(self, mol: str = None):
        pass
