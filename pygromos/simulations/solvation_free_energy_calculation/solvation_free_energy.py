"""
File: automatic calculation of Solvation Free Energies
Warnings: this class is WIP!

Description:
    For a given gromos_system (or smiles) the solvation Free Energy is automaticaly calculated.

    Main elements:
        1) create single molecule topo and conformation
        2) Generate a box filled with that molecule and according TOPO
        3) Energy-Minimize that system
        4) Equilibrate that system
        5) Run simulations at different lambda points (lambda = 0 one molecule is not present; lambda = 1 same
        molecule is present)
        6) Calculate solvation free energies by integrating over the different lambda points

Author: Paul Katzberger
"""

# Imports
from pygromos.files.gromos_system import Gromos_System
from pygromos.files.gromos_system.ff.forcefield_system import forcefield_system
from openff.toolkit.topology import Molecule
from pygromos.files.gromos_system.ff.openforcefield2gromos import openforcefield2gromos
from pygromos.files.topology.top import Top
from pygromos.files.coord import Cnf
from pygromos.files.gromos_system.ff.openforcefield2gromos import openforcefield2gromos
from pygromos.utils import bash
from rdkit import Chem
import os
import numpy as np
import warnings
from copy import deepcopy
from pygromos.simulations.hpc_queuing.submission_systems.local import LOCAL as subSys_local
from pygromos.simulations.hpc_queuing.submission_systems.lsf import LSF as subSys_lsf
from pygromos.simulations.modules.preset_simulation_modules import emin
from pygromos.data.simulation_parameters_templates import template_emin, template_md, template_sd
from pygromos.simulations.modules.preset_simulation_modules import sd
from pygromos.files.simulation_parameters.imd import Imd
from pygromos.files.blocks.imd_blocks import COVALENTFORM, COMTRANSROT, AMBER
from pygromos.gromos.pyGromosPP.ran_box import ran_box
from pygromos.simulations.modules.general_simulation_modules import simulation
from pygromos.simulations.hpc_queuing.job_scheduling.workers.analysis_workers import simulation_analysis
from pygromos.files.topology.ptp import Pertubation_topology
from pygromos.files.blocks.topology_blocks import pertubation_lam_state, atom_lam_pertubation_state, PERTATOMPARAM, \
    TITLE, SCALEDINTERACTIONS
from pygromos.files.blocks.imd_blocks import PRESSURESCALE, PERTURBATION, MULTIBATH
from pygromos.files.blocks.imd_blocks import RANDOMNUMBERS
from pygromos.analysis.error_estimate import ee
from pygromos.utils import bash
import scipy.integrate as integrate
import pandas as pd
from pygromos.files.trajectory.trg import Trg

class Solvation_free_energy_calculation:
    '''
    A class used to calculate Solvation Free Energies
    ...

    Attributes
    ----------
    input_system : Gromos_System or str or Chem.rdchem.Mol
        A prebuilt Gromos Syste an SMILES string or an RDKit molecule to use for the calculation
    work_folder : str
        Folder to work in
    system_name : str
        Name of the newly created Gromos System
    forcefield : forcefield_system
        forcefield_system to use for calculations
    gromosXX : str
        which GromosXX to use
    gromosPP : str
        which GromosPP to use
    useGromosPlsPls : bool
        Whether GromosPP should be used (Not implemented!)
    verbose : bool
        verbose flag
    num_molecules : int
        Number of molecules to build box with
    density : float
        Density of the used molecule
    num_atoms : int
        Number of atoms in the molecule
    subsystem : str
        Subsystem to use (either lsf for LSF or other for LOCAL)
    nmpi : int
        Number of MPI cores to use
    nomp : int
        Number of OMP cores to use
    amberscaling : bool
        whether to use amberscaling (not implemented on all GROMOS Branches)
    '''
    def __init__(self, input_system: Gromos_System or str or Chem.rdchem.Mol, work_folder: str,
                 system_name: str = "dummy", forcefield: forcefield_system = forcefield_system(name="54A7"),
                 gromosXX: str = None, gromosPP: str = None, useGromosPlsPls: bool = True,
                 verbose: bool = True, num_molecules: int = 512, density: float = 1000, num_atoms: int = 0,
                 subsystem: str = "lsf",nmpi = 6, nomp = 1, provided_topo : str = None, amberscaling = False) -> None:

        self.verbose = verbose
        self.amberscaling = amberscaling

        # system variables
        if type(input_system) is Gromos_System:
            self.groSys_liq = input_system
        elif (type(input_system) is str) or (type(input_system) is Chem.rdchem.Mol):
            self.groSys_liq = Gromos_System(system_name="ff2_" + input_system, work_folder=work_folder,
                                 in_smiles=input_system, auto_convert=True, Forcefield=forcefield,adapt_imd_automatically=False)

            if provided_topo:
                self.groSys_liq.top = Top(provided_topo)

            # test if openforcefield is used
            # Create Topology and convert to GROMOS
            elif forcefield.name == "openforcefield":
                molecule = Molecule.from_smiles(input_system)
                molecule.name = system_name[:4]
                # Change name to system_name for residue
                top = openforcefield2gromos(molecule).convert_return()

                # Setter not working, manual check needed
                assert isinstance(top,Top)
                self.groSys_liq._top = top
                self.groSys_liq.rebase_files()


        self.work_folder = work_folder
        self.system_name = system_name

        # create folders and structure
        try:
            os.mkdir(path=work_folder)
        except:
            if verbose:
                warnings.warn("Folder does already exist")
            else:
                pass


        self.groSys_liq.work_folder = work_folder + "/" + system_name + "_liq"
        self.groSys_liq.rebase_files()

        # Setup Submission system
        self._nmpi = nmpi
        self._nomp = nomp
        self._subsystem = subsystem

        #Create System
        self.create_new_submission_system()


        self.gromosXX = self.groSys_liq.gromosXX
        self.gromosPP = self.groSys_liq.gromosPP

        # parameters for liquid simulation
        self.num_molecules = num_molecules
        self.density = density
        self.temperature = 298.15
        self.num_atoms = num_atoms

        # read in defined imd from functions
        self.imd_liq_min = self.create_liq_min_imd()
        self.imd_liq_min.write(self.work_folder + system_name + "temp_min.imd")
        self.imd_liq_min.path = self.work_folder + system_name + "temp_min.imd"
        self.imd_liq_ti = self.create_liq_ti_imd()

        self.groSys_liq_final = None

        self.useGromosPlsPls = useGromosPlsPls


    def run(self):
        '''
        Run five step Process to calculate Self Solvation Free Energies
        1)  Create Topology and Box with Molecules
        2)  Energy minimize the system
        3)  Equilibrate the system
        4)  Setup Pertubation and sample different lambda points
        5)  Analyse the simulations at different lambda points and calculate Solvation Free Energies with Error
        Estimates
        Returns
        -------

        '''
        self.create_liq()
        emin_sys, jobID = self.minimize_liq(gromos_system=self.groSys_liq,prev_JobID=-1)
        eq_sys, jobID = self.eq_liq(gromos_system=emin_sys,prev_JobID=jobID)
        ti_sys, jobID = self.ti_liq(gromos_system=eq_sys,prev_JobID=jobID)
        self.calculate_solvation_free_energy(ti_sys,jobID)

    def create_liq(self):
        '''
        Create Coordinates and Topology for Liquid System
        '''

        # Test for fast add up
        if not np.log2(self.num_molecules) % 1:
            doubleing = int(np.log2(self.num_molecules))
            for i in range(doubleing):
                self.groSys_liq.top += self.groSys_liq.top
        else:
            new_top = Top(None)
            for i in range(self.num_molecules):
                new_top += self.groSys_liq.top
            self.groSys_liq.top = new_top


        # Create location for coordinates
        coord_dir = self.work_folder + "/coord/"
        try:
            os.mkdir(path=coord_dir)
        except:
            if self.verbose:
                warnings.warn("Folder does already exist")
            else:
                pass

        out_cnf_path = coord_dir + self.system_name + str(self.num_molecules) + ".cnf"

        # Create Box
        box_cnf_path = ran_box(in_top_path=self.groSys_liq.top.path,
                               in_cnf_path=self.groSys_liq.cnf.path,
                               out_cnf_path=out_cnf_path,
                               nmolecule=self.num_molecules,
                               dens=self.density)

        self.groSys_liq.cnf = Cnf(in_value=box_cnf_path)

        # reset liq system
        self.groSys_liq.rebase_files()

    def minimize_liq(self, gromos_system: Gromos_System, prev_JobID: int):
        '''
        Minimize Liquide
        Parameters
        ----------
        gromos_system : Gromos_System
            input Gromos System to minimize
        prev_JobID
            previous jobID to wait for
        Returns
        -------
        emin_sys : Gromos_System
            Energy minimized Gromos System
        jobID : int
            ID of submitted job for the next job to wait for (only relevant for LSF)
        '''
        gromos_system.adapt_imd_automatically = False
        gromos_system.imd = self.imd_liq_min
        emin_sys, jobID = simulation(in_gromos_system=gromos_system, project_dir=self.groSys_liq.work_folder,
                                     step_name=self.system_name + "_emin",
                                     submission_system=self.submissonSystem,
                                     analysis_script=simulation_analysis.do,
                                     verbose=self.verbose)
        return emin_sys, jobID

    def equilibration_gromos_system_step(self, gromos_system: Gromos_System, NTIVAL: int = 0, NTISHI: int = 0,
                                         TEMPI: float = 0, TEMPO: float = 0, prevID: int = -1, run_name: str = "",
                                         natoms: int = 0):
        '''
        Perform one step in equilibration routine
        Parameters
        ----------
        gromos_system : Gromos_System
            Gromos System to work on
        NTIVAL : int
            NTIVAL Parameter to set
        NTISHI : int
            NTISHI Parameter to set
        TEMPI : float
            TEMPI Parameter to set
        TEMPO : float
            TEMPO Parameter to set
        prevID : int
            ID of previous run (for LSF Submission Syste)
        run_name : str
            Name of the current run (note needs to be unique for LSF not to crash)
        natoms : int
            Number of Atoms in the molecule

        Returns
        -------
        md_sys : Gromos_System
            Equilibrated gromos system
        JobID : int
            ID of the submitted job for the next job to wait on
        '''

        # adapt imd
        eq_imd = Imd(template_md)
        eq_imd.SYSTEM.NSM = 0
        eq_imd.STEP.NSTLIM = 25000

        multibath = MULTIBATH(ALGORITHM=0, NBATHS=1, TEMP0=[TEMPO], TAU=[0.1], DOFSET=1,
                              LAST=[natoms], COMBATH=[1], IRBATH=[1])
        eq_imd.MULTIBATH = multibath

        eq_imd.PRESSURESCALE.COUPLE = 1

        eq_imd.FORCE.NEGR = 1
        eq_imd.FORCE.NRE = [natoms]

        eq_imd.NONBONDED.EPSRF = 4
        eq_imd.NONBONDED.NSHAPE = -1
        eq_imd.NONBONDED.NQEVAL = 100000

        # add COVALENTFORM Block
        covalentform = COVALENTFORM(NTBBH=0, NTBAH=0, NTBDN=0)
        eq_imd.add_block(block=covalentform)

        if self.amberscaling:
            # Add AMBER Block
            amber_block = AMBER(True, 1.2)
            eq_imd.add_block(block=amber_block)

        eq_imd.INITIALISE.IG = 3
        eq_imd.INITIALISE.NTIVEL = NTIVAL
        eq_imd.INITIALISE.NTISHI = NTISHI
        eq_imd.INITIALISE.TEMPI = TEMPI

        eq_imd.WRITETRAJ.NTWX = 500
        eq_imd.WRITETRAJ.NTWE = 500

        eq_imd.PRINTOUT.NTPR = 500
        eq_imd.randomize_seed()

        if prevID == -1:
            md_sys, JobID = simulation(in_gromos_system=gromos_system, project_dir=self.groSys_liq.work_folder,
                                       step_name="eq" + run_name, in_imd_path=eq_imd,
                                       submission_system=self.submissonSystem,
                                       analysis_script=simulation_analysis.do,
                                       verbose=self.verbose)
        else:
            md_sys, JobID = simulation(in_gromos_system=gromos_system, project_dir=self.groSys_liq.work_folder,
                                       step_name="eq" + run_name, in_imd_path=eq_imd,
                                       submission_system=self.submissonSystem,
                                       analysis_script=simulation_analysis.do,
                                       verbose=self.verbose, previous_simulation_run=prevID)

        return md_sys, JobID

    def eq_liq(self, gromos_system: Gromos_System, prev_JobID: int):
        '''
        Function to equilibrate the liquid following a 5 step scheme by sequentially calling
        self.equilibration_gromos_system_step() with different paramteres
        Parameters
        ----------
        gromos_system : Gromos_System
            Gromos System to use (normally energy minimized system)
        prev_JobID : int
            ID of previous run for LSF queue

        Returns
        -------

        '''
        # Start equilibration
        gromos_system.imd = Imd(template_md)

        # Equilibrate system and pass on for next step
        eq_sys1, JobID1 = self.equilibration_gromos_system_step(gromos_system=gromos_system, NTIVAL=1, NTISHI=1, TEMPI=60,
                                                                TEMPO=60,
                                                                prevID=prev_JobID, run_name=self.system_name + "_eq_1",
                                                                natoms=self.num_atoms * self.num_molecules)
        eq_sys2, JobID2 = self.equilibration_gromos_system_step(gromos_system=eq_sys1, NTIVAL=0, NTISHI=0, TEMPI=0,
                                                                TEMPO=120,
                                                                prevID=JobID1, run_name=self.system_name + "_eq_2",
                                                                natoms=self.num_atoms * self.num_molecules)
        eq_sys3, JobID3 = self.equilibration_gromos_system_step(gromos_system=eq_sys2, NTIVAL=0, NTISHI=0, TEMPI=0,
                                                                TEMPO=180,
                                                                prevID=JobID2, run_name=self.system_name + "_eq_3",
                                                                natoms=self.num_atoms * self.num_molecules)
        eq_sys4, JobID4 = self.equilibration_gromos_system_step(gromos_system=eq_sys3, NTIVAL=0, NTISHI=0, TEMPI=0,
                                                                TEMPO=240,
                                                                prevID=JobID3, run_name=self.system_name + "_eq_4",
                                                                natoms=self.num_atoms * self.num_molecules)
        eq_sys5, JobID5 = self.equilibration_gromos_system_step(gromos_system=eq_sys4, NTIVAL=0, NTISHI=0, TEMPI=0,
                                                                TEMPO=285,
                                                                prevID=JobID4, run_name=self.system_name + "_eq_5",
                                                                natoms=self.num_atoms * self.num_molecules)

        return eq_sys5, JobID5

    def ti_liq(self, gromos_system: Gromos_System, prev_JobID: int, n_points:int = 21):
        '''
        Perform TI of Liquid
        Parameters
        ----------
        gromos_system : Gromos_System
            Gromos system to start TI with
        prev_JobID : int
            job ID of previous job to wait for (in LSF queue)
        n_points : int
            number of lambda points to simulate (typically 21 or 41)

        Returns
        -------
        JobID : int
            ID of last submitted job (WARNING it could be that this job finishes before some others manual check needed)
        '''
        # First generate ptp file
        gromos_system = self.add_ptp_file(gromos_system)

        # Set up dictunaries
        systems = {}
        jobIDs = {}

        ti_dir = self.groSys_liq.work_folder + "/ti/"
        # create folders and structure
        try:
            os.mkdir(path=ti_dir)
        except:
            if self.verbose:
                warnings.warn("Folder does already exist")
            else:
                pass

        # Do Calculations
        for lambda_point in range(n_points):

            # calculate lambda value
            rlam = lambda_point / (n_points - 1)

            # Ensure that rlam is on the order of 3 (i.e. not 1.99999999999998)
            rlam = np.round(rlam, 3)

            iteration = 0
            # First run has longer equilibration time and different naming scheme
            if rlam == 0:
                system_name = self.system_name + "_L" + str(rlam) + "_" + str(iteration)
                c_ti_sys0_prep = self.set_up_ti_lambda_run(gromos_system=gromos_system, NSTLIM=175000, RLAM=rlam,
                                                        system_name=system_name)
                print(c_ti_sys0_prep.imd)
                c_ti_sys0_prep.adapt_imd_automatically = False
                c_ti_sys0, JobID = sd(in_gromos_system=c_ti_sys0_prep, project_dir=ti_dir + system_name,
                                      step_name=system_name, submission_system=self.submissonSystem,
                                      previous_simulation_run=prev_JobID)
            else:
                # Set system Name and previous system
                system_name = self.system_name + "_L" + str(rlam) + "_" + str(iteration)

                # Calculate previous Lambda step
                previous_rlam = np.round(rlam - 1 / (n_points - 1), 3)
                previous_system_name = self.system_name + "_L" + str(previous_rlam) + "_0"

                # Setup second step
                c_ti_sys0_prep = self.set_up_ti_lambda_run(gromos_system=systems[previous_system_name], NSTLIM=25000,
                                                           RLAM=rlam,
                                                           system_name=system_name)
                c_ti_sys0, JobID = sd(in_gromos_system=c_ti_sys0_prep, project_dir=ti_dir + system_name,
                                      step_name=system_name, submission_system=self.submissonSystem,
                                      previous_simulation_run=jobIDs[previous_system_name])

            # Put System into dictionary
            systems[system_name] = c_ti_sys0
            jobIDs[system_name] = JobID

            # Do simulation for actual lambda calc
            iteration = 1
            system_name = self.system_name + "_L" + str(rlam) + "_" + str(iteration)
            previous_system_name = self.system_name + "_L" + str(rlam) + "_" + str(iteration - 1)

            # Setup second step
            c_ti_sys1_prep = self.set_up_ti_lambda_run(gromos_system=systems[previous_system_name], NSTLIM=25000,
                                                    RLAM=rlam,
                                                  system_name=system_name)
            c_ti_sys1, JobID = sd(in_gromos_system=c_ti_sys1_prep, project_dir=ti_dir + system_name,
                                  step_name=system_name, submission_system=self.submissonSystem,
                                  previous_simulation_run=jobIDs[previous_system_name])

            # Put System into dictionary
            systems[system_name] = c_ti_sys1
            jobIDs[system_name] = JobID

            iteration = 2
            system_name = self.system_name + "_L" + str(rlam) + "_" + str(iteration)
            previous_system_name = self.system_name + "_L" + str(rlam) + "_" + str(iteration - 1)

            c_ti_sys2_prep = self.set_up_ti_lambda_run(gromos_system=systems[previous_system_name], NSTLIM=500000,
                                                   RLAM=rlam,
                                                  system_name=system_name)
            c_ti_sys2, JobID = sd(in_gromos_system=c_ti_sys2_prep, project_dir=ti_dir + system_name,
                                  step_name=system_name, submission_system=self.submissonSystem,
                                  previous_simulation_run=jobIDs[previous_system_name])

            # Put System into dictionary
            systems[system_name] = c_ti_sys2
            jobIDs[system_name] = JobID

        return JobID

    def set_up_ti_lambda_run(self, gromos_system: Gromos_System, NSTLIM: int = 0, RLAM: float = 0,
                             system_name: str = ""):
        '''
        Setup on TI run by setting the lambda value and steps
        Parameters
        ----------
        gromos_system : Gromos_System
            Gromos System to work with
        NSTLIM : int
            NSTLIM step parameter to use
        RLAM : float
            lambda value to use for simulation
        system_name : str
            Name of the system (must be unique)

        Returns
        -------
        gromos_system : Gromos_System
            Prepared Gromos System
        '''
        # Do not change imd automatically
        gromos_system.adapt_imd_automatically = False

        # Setup storage
        ti_dir = self.groSys_liq.work_folder + "/ti/"

        # Rename System and set storage path
        gromos_system.imd = self.imd_liq_ti
        gromos_system.name = system_name
        gromos_system.work_folder = ti_dir + system_name

        # Set number of steps and Lambda Parameter
        gromos_system.imd.STEP.NSTLIM = NSTLIM
        gromos_system.imd.PERTURBATION.RLAM = RLAM

        # Save files
        gromos_system.rebase_files()

        print('Setup imd',gromos_system.imd)
        return gromos_system

    def add_ptp_file(self, gromos_system: Gromos_System):
        '''
        Function to generate a perturbation Topology
        Parameters
        ----------
        gromos_system : Gromos_System
            Gromos system to add perturbation Topology to

        Returns
        -------
        gromos_system : Gromos_System
            Gromos system with added perturbation Topology
        '''

        # Setup lambda States
        pert_atoms = []

        # Set only atoms of the first molecule
        for atom_line in gromos_system.top.SOLUTEATOM[:self.num_atoms]:
            # Get atoms and the correct states
            states = {}
            phys_state = pertubation_lam_state(IAC=atom_line.IAC, MASS=atom_line.MASS, CHARGE=atom_line.CG)

            # Generate a second set of states (a bit hackish)
            states = {atom_line.MRES: phys_state,
                      atom_line.MRES + 1: phys_state}
            pert_atom = atom_lam_pertubation_state(atom_line.ATNM, RES=atom_line.MRES, NAME=atom_line.PANM,
                                                   STATES=states)
            pert_atoms.append(pert_atom)

        # Make Pertubation block
        pert_atom_block = PERTATOMPARAM(pert_atoms)

        # Generate ptp file
        gromos_system.ptp = Pertubation_topology(in_value=None)
        gromos_system.ptp.PERTATOMPARAM = pert_atom_block
        gromos_system.ptp.TITLE = TITLE("Automatic generated pertubation file. For Solvation Free Energy calculation")

        content_dict = [1, 1, 2, 1, 0]
        scaled = SCALEDINTERACTIONS(values=content_dict)
        gromos_system.ptp.add_block(block=scaled)

        return gromos_system

    def create_liq_min_imd(self):
        """
        Read in template_emin and make changes
        """
        emin_imd = Imd(template_emin)

        # adjust parameters
        emin_imd.SYSTEM.NPM = 1
        emin_imd.SYSTEM.NSM = 0

        # Adjust Energymin
        emin_imd.ENERGYMIN.NMIN = 5000

        # Adjust number of Steps
        emin_imd.STEP.NSTLIM = 20000

        # add COVALENTFORM Block
        covalentform = COVALENTFORM(NTBBH=0, NTBAH=0, NTBDN=0)
        emin_imd.add_block(block=covalentform)

        # Adjust force block
        emin_imd.FORCE.BONDS = 0
        emin_imd.FORCE.NEGR = 1
        emin_imd.FORCE.NRE = [self.num_molecules * self.num_atoms]

        # Adjust Pairlist
        emin_imd.PAIRLIST.ALGORITHM = 1
        emin_imd.PAIRLIST.RCUTP = 1.4
        emin_imd.PAIRLIST.RCUTL = 1.4
        emin_imd.PAIRLIST.SIZE = 0.7

        # Adjust Nonebonded
        emin_imd.NONBONDED.EPSRF = 4

        if self.amberscaling:
            print("Amber Scling is enabled")
            # Add AMBER Block
            amber_block = AMBER(True,1.2)
            emin_imd.add_block(block=amber_block)

        # adjust initialize
        emin_imd.INITIALISE.NTIVEL = 1
        emin_imd.INITIALISE.NTISHK = 3
        emin_imd.INITIALISE.NTICOM = 1
        emin_imd.INITIALISE.IG = 2
        emin_imd.INITIALISE.TEMPI = 298.15

        return emin_imd

    def create_liq_eq_imd(self):
        '''
        Not implemented done on the fly in the eq_liq function
        '''
        pass

    def create_liq_ti_imd(self):
        '''
        Create IMD to perform TI
        Returns
        -------
        ti_imd : Imd
            Processed Imd file
        '''

        # set new template
        ti_imd = Imd(template_sd)

        # Add correct random number generator
        randomnumbers = RANDOMNUMBERS(0, -1)

        ti_imd.randomize_seed()
        ti_imd.add_block(block=randomnumbers)

        # Change Friction Constant
        ti_imd.STOCHDYN.CFRIC = 10
        ti_imd.STOCHDYN.TEMPSD = 298.15

        # Add Pressurescale
        pressurescale = PRESSURESCALE(COUPLE=2, SCALE=1, COMP=0.0004575,
                                      TAUP=0.5, VIRIAL=2, SEMIANISOTROPIC=[1, 1, 1],
                                      PRES0=[[0.6102, 0, 0], [0, 0.6102, 0], [0, 0, 0.6102]])
        ti_imd.add_block(block=pressurescale)

        # Adjust forces
        ti_imd.FORCE.BONDS = 0
        ti_imd.FORCE.NEGR = 2
        ti_imd.FORCE.NRE = [int(self.num_atoms), int(self.num_molecules * self.num_atoms)]

        # Adjust covalentform
        ti_imd.COVALENTFORM.NTBBH = 0
        ti_imd.COVALENTFORM.NTBAH = 0
        ti_imd.COVALENTFORM.NTBDN = 0

        # Adjust Constraint
        ti_imd.CONSTRAINT.NTC = 3

        # Adjust PAIRLIST
        ti_imd.PAIRLIST.NSNB = 5
        ti_imd.PAIRLIST.RCUTP = 0.8
        ti_imd.PAIRLIST.RCUTL = 1.4
        ti_imd.PAIRLIST.SIZE = 0.4
        ti_imd.PAIRLIST.TYPE = 0

        # Adjust NONBONDED
        ti_imd.NONBONDED.RCRF = 1.4
        ti_imd.NONBONDED.EPSRF = 4
        ti_imd.NONBONDED.NSLFEXCL = 1
        ti_imd.NONBONDED.NSHAPE = -1
        ti_imd.NONBONDED.ASHAPE = 1.4
        ti_imd.NONBONDED.NA2CLC = 2
        ti_imd.NONBONDED.TOLA2 = "1e-10"
        ti_imd.NONBONDED.EPSLS = 0

        # Adjust WRITETRAJ
        ti_imd.WRITETRAJ.NTWX = 0
        ti_imd.WRITETRAJ.NTWSE = 0
        ti_imd.WRITETRAJ.NTWV = 0
        ti_imd.WRITETRAJ.NTWF = 0
        ti_imd.WRITETRAJ.NTWE = 250
        ti_imd.WRITETRAJ.NTWG = 250
        ti_imd.WRITETRAJ.NTWB = 0

        # Adjust BONDCOND
        ti_imd.BOUNDCOND.NTB = 1
        ti_imd.BOUNDCOND.NDFMIN = 3

        # Adjust Printout
        ti_imd.PRINTOUT.NTPR = 0

        if self.amberscaling:
            # Add AMBER Block
            amber_block = AMBER(True,1.2)
            ti_imd.add_block(block=amber_block)

        # Adjust Inintialize
        ti_imd.INITIALISE.NTIVEL = 0
        ti_imd.INITIALISE.NTISHI = 0
        ti_imd.INITIALISE.NTICOM = 1
        ti_imd.INITIALISE.TEMPI = 298.15


        # Add PERTURBATION Block
        pert_block = PERTURBATION(NTG=1, NRDGL=0, RLAM=0, DLAMT=0,
                                  ALPHLJ=0.5, ALPHC=0.5, NLAM=1, NSCALE=1)
        ti_imd.add_block(block=pert_block)
        print("TI imd",ti_imd)
        return ti_imd

    def calculate_solvation_free_energy(self,n_points : int = 21):
        '''
        Function to Calculate solvation free energy by integrating over lambda points
        Parameters
        ----------
        n_points : int
            Number of Lambda points used

        Returns
        -------
        solv_energy : float
            Solvation Free Energy
        error : float
            Error Estimate
        '''

        # Get dhdl information
        Metrics = pd.DataFrame()

        # Setup Storage
        lambdas = []
        averages = []
        rmsds = []
        ees = []

        # Set Lambda Points
        lambda_points = n_points

        for l in range(lambda_points):
            rlam = l / (lambda_points - 1)
            rlam = np.round(rlam, 3)

            # Get system Name
            iteration = 2
            system_name = self.system_name + "_L" + str(rlam) + "_" + str(iteration)
            path = self.groSys_liq.work_folder + "/ti/" + system_name + "/" + system_name + "/"

            # Extract files
            command = "tar" + " "
            command += "-xf" + " "
            command += path + "simulation.tar "
            command += "-C " + path
            bash.execute(command, verbose=True)

            # get trg info
            trg = Trg(path + "simulation/" + system_name + "_1/" + system_name + "_1.trg.gz")

            # Get dhdl
            dhdls = [trg.database.totals[i][2] for i in range(len(trg.database.totals))]

            # Remove simulation directory
            bash.remove_file(path + "simulation", recursive=True)

            # Get values
            lambdas.append(rlam)
            averages.append(np.mean(dhdls))
            rmsds.append(ee(dhdls).calculate_rmsd())
            ees.append(ee(dhdls).calculate_ee())

        Metrics["Lambda"] = lambdas
        Metrics["avg"] = averages
        Metrics["rmsd"] = rmsds
        Metrics["ee"] = ees

        # Integrate over lambda points for value and error
        solv_energy = integrate.simpson(Metrics["avg"], Metrics["Lambda"]) * -1
        error = integrate.simpson(Metrics["ee"], Metrics["Lambda"]) * -1

        if self.verbose:
            print(Metrics)
            print("Solvation Free Energy:",solv_energy,"pm",error)

        return solv_energy, error

    def create_new_submission_system(self):
        if self._subsystem == "lsf":
            self.submissonSystem = subSys_lsf(nmpi=self._nmpi, nomp=self._nomp, verbose=self.verbose)
        else:
            self.submissonSystem = subSys_local(nmpi=self._nmpi, nomp=self._nomp, verbose=self.verbose)

    @property
    def nmpi(self):
        return self._nmpi

    @nmpi.setter
    def nmpi(self,nmpi=1):
        self._nmpi = nmpi
        self.create_new_submission_system()

    @property
    def nomp(self):
        return self._nomp

    @nomp.setter
    def nomp(self, nomp=1):
        self._nomp = nomp
        self.create_new_submission_system()

    @property
    def subsystem(self):
        return self._subsystem

    @subsystem.setter
    def subsystem(self, subsystem="local"):
        self._subsystem = subsystem
        self.create_new_submission_system()