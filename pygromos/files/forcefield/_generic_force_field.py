"""
File:            forcefield managment for PyGromosTools

Description:
    This class handles all possible forcefields gromos_system accepts.
    It implements a few methods to handle the forcefield files and store some general information about the forcefield.

    Attention: This is only a template class! It'should be used as a super class for all forcefields! But not on its own!


Author: Marc Lehner
"""

from pygromos.files.coord.cnf import Cnf
from pygromos.files.topology.top import Top


class _generic_force_field:
    def __init__(
        self, name: str = "generic", path_to_files: str = None, auto_import: bool = True, verbose: bool = False
    ):
        self.name = name
        self.path_to_files = path_to_files
        self.verbose = verbose

    def auto_import_ff(self, **kwargs):
        pass

    def create_top(self, mol: str, in_top: Top = None, **kwargs) -> Top:
        raise NotImplementedError("This is a template class! It'should be used as a super class for all forcefields!")

    def create_cnf(self, mol: str, in_cnf: Cnf = None, **kwargs) -> Cnf:
        raise NotImplementedError("This is a template class! It'should be used as a super class for all forcefields!")
