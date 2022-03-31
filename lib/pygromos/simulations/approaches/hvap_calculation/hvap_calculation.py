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
import os
import time
import warnings
from copy import deepcopy

from rdkit import Chem

from pygromos.gromos.pyGromosPP.ran_box import ran_box
from pygromos.files.gromos_system.gromos_system import Gromos_System
from pygromos.simulations.approaches.hvap_calculation import hvap_input_files
from pygromos.files.forcefield._generic_force_field import _generic_force_field

from pygromos.simulations.hpc_queuing.submission_systems import get_submission_system, _submission_system
from pygromos.simulations.modules.general_simulation_modules import simulation
from pygromos.simulations.hpc_queuing.job_scheduling.workers.analysis_workers import simulation_analysis

from pygromos.files.simulation_parameters.imd import Imd
from pygromos.files.topology.top import Top
from pygromos.utils.utils import time_wait_s_for_filesystem

# automatically get Local or LSF submission system (depending on hostname)
subSystem = get_submission_system()


class Hvap_calculation:
    dens_modifier: float = 0.7
    submissonSystem_gas: _submission_system
    submissonSystem_liq: _submission_system

    def __init__(
        self,
        input_system: Gromos_System or str or Chem.rdchem.Mol,
        work_folder: str,
        system_name: str = "dummy",
        forcefield: _generic_force_field = _generic_force_field(),
        in_gromosXX_bin_dir: str = None,
        in_gromosPP_bin_dir: str = None,
        useGromosPlsPls: bool = True,
        verbose: bool = True,
    ) -> None:
        """For a given gromos_system (or smiles) the heat of vaporization is automaticaly calculated

        Parameters
        ----------
        input_system : Gromos_SystemorstrorChem.rdchem.Mol
            single molecule gromos_sytem or rdkit Molecule or SMILES
        """
        # system variables
        if type(input_system) is Gromos_System:
            self.groSys_gas = input_system
            if in_gromosXX_bin_dir is not None:
                self.groSys_gas.gromosXX = in_gromosXX_bin_dir
            if in_gromosPP_bin_dir is not None:
                self.groSys_gas.gromosPP = in_gromosPP_bin_dir
        elif (type(input_system) is str) or (type(input_system) is Chem.rdchem.Mol):
            self.groSys_gas = Gromos_System(
                work_folder=work_folder,
                system_name=system_name,
                in_smiles=input_system,
                in_gromosXX_bin_dir=in_gromosXX_bin_dir,
                in_gromosPP_bin_dir=in_gromosPP_bin_dir,
                forcefield=forcefield,
                in_imd_path=hvap_input_files.imd_hvap_gas_sd,
                verbose=verbose,
            )

        self.work_folder = work_folder
        self.system_name = system_name

        self.submissonSystem_gas = subSystem(job_duration="4:00")
        # give the liquid simulation as much cores as possible to speedup the long liquid simulation.
        # if you use a mpi version of gromos use nmpi=4 (or more)
        # the attribute can be overwritten by the user at runtime after the class is initiated
        self.submissonSystem_liq = subSystem(nomp=4, job_duration="24:00")

        # create folders and structure
        try:
            os.mkdir(path=work_folder)
        except FileExistsError:
            if verbose:
                warnings.warn("Folder does already exist")
            else:
                pass
        self.groSys_gas.work_folder = work_folder + "/" + system_name + "_gas"
        self.groSys_gas.rebase_files()
        self.groSys_liq = deepcopy(self.groSys_gas)
        self.groSys_liq.work_folder = work_folder + "/" + system_name + "_liq"
        self.groSys_liq.rebase_files()

        self.gromosXX = self.groSys_gas.gromosXX
        self.gromosPP = self.groSys_gas.gromosPP

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
        self.density = 700
        self.temperature = 298.15

        self.groSys_gas_final = None
        self.groSys_liq_final = None

        self.useGromosPlsPls = useGromosPlsPls

        self.verbose = verbose

    def run(self) -> int:
        self.create_liq()
        self.run_gas()
        self.run_liq()
        return self.calc_hvap()

    def create_liq(self):
        # create liq top
        if self.useGromosPlsPls:
            try:
                self.gromosPP.com_top(
                    self.groSys_gas.top.path,
                    topo_multiplier=self.num_molecules,
                    out_top_path=self.work_folder + "/temp.top",
                )
                tempTop = Top(in_value=self.work_folder + "/temp.top")
                tempTop.write(out_path=self.work_folder + "temp.top")
                time.sleep(time_wait_s_for_filesystem)  # wait for file to write and close
                self.groSys_liq.top = tempTop
            except Exception as e:
                self.groSys_liq.top = self.groSys_gas.top * self.num_molecules
                if self.verbose:
                    print(e)
        else:
            self.groSys_liq.top = self.groSys_gas.top * self.num_molecules

        # create liq cnf
        if self.useGromosPlsPls:
            self.gromosPP.ran_box(
                in_top_path=self.groSys_gas.top.path,
                in_cnf_path=self.groSys_gas.cnf.path,
                out_cnf_path=self.work_folder + "/temp.cnf",
                nmolecule=self.num_molecules,
                dens=self.dens_modifier * self.density,
                threshold=0.12,
                layer=True,
            )
        else:
            ran_box(
                in_top_path=self.groSys_gas.top.path,
                in_cnf_path=self.groSys_gas.cnf.path,
                out_cnf_path=self.work_folder + "/temp.cnf",
                nmolecule=self.num_molecules,
                dens=self.dens_modifier * self.density,
            )
        time.sleep(time_wait_s_for_filesystem)  # wait for file to write and close
        self.groSys_liq.cnf = self.work_folder + "/temp.cnf"

        # reset liq system
        self.groSys_liq.rebase_files()

    def run_gas(self):

        # min
        self.groSys_gas.imd = self.imd_gas_min
        self.groSys_gas.prepare_for_simulation()
        sys_emin_gas = simulation(
            in_gromos_simulation_system=self.groSys_gas,
            override_project_dir=self.groSys_gas.work_folder,
            step_name="1_emin",
            submission_system=self.submissonSystem_gas,
            analysis_script=simulation_analysis.do,
            verbose=self.verbose,
        )

        # eq
        sys_emin_gas.imd = self.imd_gas_eq
        sys_emin_gas.prepare_for_simulation()
        sys_eq_gas = simulation(
            in_gromos_simulation_system=sys_emin_gas,
            override_project_dir=self.groSys_gas.work_folder,
            step_name="2_eq",
            submission_system=self.submissonSystem_gas,
            analysis_script=simulation_analysis.do,
            verbose=self.verbose,
        )

        # sd
        sys_eq_gas.imd = self.imd_gas_eq
        sys_eq_gas.prepare_for_simulation()
        sys_sd_gas = simulation(
            in_gromos_simulation_system=sys_eq_gas,
            override_project_dir=self.groSys_gas.work_folder,
            step_name="3_sd",
            submission_system=self.submissonSystem_gas,
            analysis_script=simulation_analysis.do,
            verbose=self.verbose,
        )

        self.groSys_gas_final = sys_sd_gas

    def run_liq(self):

        # minsys_emin_liq, jobID
        self.groSys_liq.imd = self.imd_liq_min
        self.groSys_liq.prepare_for_simulation()
        sys_emin_liq = simulation(
            in_gromos_simulation_system=self.groSys_liq,
            override_project_dir=self.groSys_liq.work_folder,
            step_name="1_emin",
            submission_system=self.submissonSystem_liq,
            analysis_script=simulation_analysis.do,
            verbose=self.verbose,
        )

        # eq
        sys_emin_liq.imd = self.imd_liq_eq
        sys_emin_liq.prepare_for_simulation()
        sys_eq_liq = simulation(
            in_gromos_simulation_system=sys_emin_liq,
            override_project_dir=self.groSys_liq.work_folder,
            step_name="2_eq",
            submission_system=self.submissonSystem_liq,
            analysis_script=simulation_analysis.do,
            verbose=self.verbose,
        )

        # md
        sys_eq_liq.imd = self.imd_liq_md
        sys_eq_liq.prepare_for_simulation()
        sys_md_liq = simulation(
            in_gromos_simulation_system=sys_eq_liq,
            override_project_dir=self.groSys_liq.work_folder,
            step_name="3_sd",
            submission_system=self.submissonSystem_liq,
            analysis_script=simulation_analysis.do,
            verbose=self.verbose,
        )

        self.groSys_liq_final = sys_md_liq

    def calc_hvap(self) -> float:
        h_vap = self.groSys_liq_final.tre.get_Hvap(
            gas_traj=self.groSys_gas_final.tre, nMolecules=self.num_molecules, temperature=self.temperature
        )
        return h_vap
