from pygromos.simulations.approaches.solvation_free_energy_calculation.solvation_free_energy import (
    Solvation_free_energy_calculation,
)
from pygromos.files.forcefield.openff.openff import OpenFF
from pygromos.simulations.hpc_queuing.submission_systems.dummy import DUMMY
import unittest


class test_sfe(unittest.TestCase):
    smiles = "c1ccccc1"
    workfolder = "/tmp/test_solvation_free_energy"
    number_of_atoms = 12

    sf = Solvation_free_energy_calculation(
        input_system=smiles,  # Gromos_System, SMILES (str) or rdkit Mol
        work_folder=workfolder,  # Folder to do calculations in
        system_name=smiles,
        # Name of the system (does not need to be smiles but convenient)
        forcefield=OpenFF(),  # Force field to use
        density=789,  # density of the liquid in kg/L
        num_molecules=512,  # number of molecules used for the calculation
        num_atoms=number_of_atoms,  # number of atoms in one molecule
        subsystem=DUMMY(),  # Subsystem to use for calculation local or lsf
        amberscaling=False,
    )  # Whether to use amberscaling (for openforcefield recommended)

    def test_constructor(self):
        print(self.sf)

    def test_min_imd(self):
        print(self.sf.create_liq_min_imd())

    def test_eq_imd(self):
        print(self.sf.create_liq_eq_imd())

    def test_ti_imd(self):
        print(self.sf.create_liq_ti_imd())

    def test_energy_groups(self):
        ti_imd = self.sf.create_liq_ti_imd()
        assert ti_imd.FORCE.NRE == [int(self.sf.num_atoms), int(self.sf.num_molecules * self.sf.num_atoms)]

    def test_ptp_file(self):
        ptp_system = self.sf.add_ptp_file(self.sf.groSys_liq)
        content_dict = [1, 1, 2, 1, 0]
        assert ptp_system.ptp.SCALEDINTERACTIONS.values == content_dict

    def test_box_generation(self):
        self.sf.create_liq()
