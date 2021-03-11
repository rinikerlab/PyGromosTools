"""
File:            automatic calculation of Hvap
Warnings: this class is WIP!

Description:
    For a given gromos_system (or smiles) the heat of vaporization is automaticaly calculated.

    Main elements:
        1) create single molecule topo and conformation
        2) run gas (single molecule) minimization
        3) run gas SD simulation (and equilibaration)
        4) generate multi molecule (liquid) topo and conformation
        5) run liquid minimization
        6) run liquid MD run (and equilibaration)
        7) calculate Hvap from gas and liquid trajectories

Author: Marc Lehner
"""

from copy import copy, deepcopy
from warnings import WarningMessage
from rdkit import Chem
import os
import warnings

from pygromos.files.gromos_system.gromos_system import Gromos_System
from pygromos.files.gromos_system.hvap_calculation import hvap_input_files

from pygromos.hpc_queuing.submission_systems.Submission_Systems import LOCAL as subSys
from pygromos.simulation_runner.simulation_building_blocks import simulation
from pygromos.hpc_queuing.job_scheduling.workers.analysis_workers import simulation_analysis

from pygromos.files.coord.cnf import Cnf
from pygromos.files.simulation_parameters.imd import Imd
from pygromos.files.topology.top import Top

from pygromos.gromos import GromosXX, GromosPP
from pygromos.utils import bash, utils

class Hvap_calculation():
    def __init__(self, input_system:Gromos_System or str or Chem.rdchem.Mol, work_folder:str, system_name:str="dummy") -> None:
        """For a given gromos_system (or smiles) the heat of vaporization is automaticaly calculated

        Parameters
        ----------
        input_system : Gromos_SystemorstrorChem.rdchem.Mol
            single molecule gromos_sytem or rdkit Molecule or SMILES
        """
        # system variables
        if type(input_system) is Gromos_System:
            self.groSys_gas = input_system
        elif type(input_system) is str:
            raise NotImplementedError("WIP use Gromos_System")
        elif type(input_system) is Chem.rdchem.Mol:
            raise NotImplementedError("WIP use Gromos_System")

        self.work_folder = work_folder
        self.system_name = system_name

        # create folders and structure
        try:
            os.mkdir(path=work_folder)
        except:
            warnings.warn("Folder does already exist")
        self.groSys_gas.work_folder = work_folder + "/" + system_name +"_gas"
        self.groSys_liq = deepcopy(self.groSys_gas)
        self.groSys_liq.work_folder = work_folder + "/" + system_name +"_liq"

        self.submissonSystem = subSys()

        # define template imd files (overwritte for specific systems) 
        # made for small molecule Hvap calculation
        self.imd_gas_min = Imd(hvap_input_files.imd_hvap_emin)
        self.imd_gas_eq = Imd(hvap_input_files.imd_hvap_gas_sd)
        self.imd_gas_sd = Imd(hvap_input_files.imd_hvap_gas_sd)
        self.imd_liq_min = Imd(hvap_input_files.imd_hvap_emin)
        self.imd_liq_eq = Imd(hvap_input_files.imd_hvap_liquid_eq)
        self.imd_liq_md = Imd(hvap_input_files.imd_hvap_liquid_md)

        # parameters for liquid simulation
        # used to multiply the single molecule system
        # made for small molecule Hvap calculation
        self.num_molecules = 512
        self.density = 1000

    def run(self) -> int:
        self.run_gas()
        self.run_liq()
        return self.calc_hvap()

    def run_gas(self):
        self.groSys_gas.rebase_files()

        #min
        simulation(in_gromos_system=self.groSys_gas, 
                    project_dir=self.groSys_gas.work_folder,
                    step_name="1_emin", 
                    in_imd_path=self.imd_gas_min,
                    submission_system=self.submissonSystem,
                    analysis_script=simulation_analysis.do)

        #eq
        simulation(in_gromos_system=self.groSys_gas, 
                    project_dir=self.groSys_gas.work_folder,
                    step_name="2_eq", 
                    in_imd_path=self.imd_gas_eq,
                    submission_system=self.submissonSystem,
                    analysis_script=simulation_analysis.do)

        #sd
        simulation(in_gromos_system=self.groSys_gas, 
                    project_dir=self.groSys_gas.work_folder,
                    step_name="3_sd", 
                    in_imd_path=self.imd_gas_sd,
                    submission_system=self.submissonSystem,
                    analysis_script=simulation_analysis.do)


    def run_liq(self):
        self.groSys_liq.rebase_files()

        #min
        simulation(in_gromos_system=self.groSys_liq, 
                    project_dir=self.groSys_liq.work_folder,
                    step_name="1_emin", 
                    in_imd_path=self.imd_liq_min,
                    submission_system=self.submissonSystem,
                    analysis_script=simulation_analysis.do)

        #eq
        simulation(in_gromos_system=self.groSys_liq, 
                    project_dir=self.groSys_liq.work_folder,
                    step_name="2_eq", 
                    in_imd_path=self.imd_liq_eq,
                    submission_system=self.submissonSystem,
                    analysis_script=simulation_analysis.do)

        #md
        simulation(in_gromos_system=self.groSys_liq, 
                    project_dir=self.groSys_liq.work_folder,
                    step_name="3_sd", 
                    in_imd_path=self.imd_liq_md,
                    submission_system=self.submissonSystem,
                    analysis_script=simulation_analysis.do)

    def calc_hvap(self) -> int:
        return 0
