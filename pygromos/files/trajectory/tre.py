"""
File:            Class for tre files in pandas
Description:
    The pandas trajectory TRE class offers a easy method to process GROMOS's .tre files in python
    The tre files are parsed into an easy to use pandas dataframe.

    This class should be a alternative for the data post processing with ene_ana in gromos++

Author: Marc Thierry Lehner

TODO: add stride option to all member functions
TODO: add support for periodic boundary condition

TODO: add ene_ana functions

"""

#imports
import pandas as pd
import numpy as np

import pygromos.files.trajectory._general_trajectory as traj
from pygromos.files.trajectory.tre_field_libs.ene_fields import gromos_2020_tre_block_names_table
from pygromos.analysis import energy_analysis as ea

class Tre(traj._General_Trajectory):
    def __init__(self, input_value: str or None, auto_save=True, stride:int=1, skip:int=0, _ene_ana_names = gromos_2020_tre_block_names_table):
        super().__init__(input_value, auto_save=auto_save, stride=stride, skip=skip)
        self.tre_block_name_table = _ene_ana_names


    """--------------------------------------------------------------
    Basic getters for Subblocks and Subsubblocks
    -----------------------------------------------------------------
    """

    def get_totals(self) -> pd.DataFrame:
        #print(self.database["totals"][0].shape, self.database["totals"][0])
        self.totals = pd.DataFrame(data = np.stack(self.database["totals"].to_numpy()), columns=self.tre_block_name_table.totals_subblock_names)
        return self.totals

    def get_totals_total(self) -> pd.DataFrame:
        self.totals_total = self.database["totals"].apply(lambda x: x[0])
        return self.totals_total

    def get_totals_bonded(self) -> pd.DataFrame:
        self.totals_bonded = self.database["totals"].apply(lambda x: x[1])
        return self.totals_bonded

    def get_totals_nonbonded(self) -> pd.DataFrame:
        self.totals_nonbonded = self.database["totals"].apply(lambda x: x[2])
        return self.totals_nonbonded

    def get_eds(self)->pd.DataFrame:
        num_states = self.database["eds"][0][0]
        if(num_states>0):
            state_strings =[[column+"_"+str(state_row) for column in self.tre_block_name_table.eds_subblock_names_singleState] for state_row in range(1, 1+int(num_states))]
            self.tre_block_name_table.eds_subblock_names = ["numstates"]
            self.tre_block_name_table.eds_subblock_names.extend(list(np.concatenate(state_strings)))
        else:
            self.tre_block_name_table.eds_subblock_names = ["numstates"]

        self.eds = pd.DataFrame(data=np.stack(self.database["eds"].to_numpy()), columns=self.tre_block_name_table.eds_subblock_names)
        return self.eds

    def get_precalclam(self)->pd.DataFrame:
        num_lam = self.database["precalclam"][0][0]
        if(num_lam>0):
            state_strings =[[column+"_"+str(state_row) for column in self.tre_block_name_table.lam_subblock_names_singleLam] for state_row in range(1, 1+int(num_lam))]
            self.tre_block_name_table.lam_subblock_names = ["nr_lambdas"]
            self.tre_block_name_table.lam_subblock_names.extend(list(np.concatenate(state_strings)))
            self.tre_block_name_table.lam_subblock_names.extend(["A_dihedral", "B_dihedral"])
        else:
            self.tre_block_name_table.lam_subblock_names = ["nr_lambdas", "A_dihedral", "B_dihedral"]

        self.precalclam = pd.DataFrame(data=np.stack(self.database["precalclam"].to_numpy()), columns=self.tre_block_name_table.lam_subblock_names)
        return self.precalclam

    """--------------------------------------------------------------
    Additional predefined function for common analysis
    -----------------------------------------------------------------
    """

    def get_density(self) -> pd.DataFrame:
        """
        Calculate the density for every frame.
        (Uses mass and the first volume entry)

        Returns
        -------
        pd.DataFrame
            Dataframe with the densities for all time steps
        """
        return self.database[["mass","volume"]].apply(lambda x: ea.get_density(mass=x[0][0], volume=x[1][0]), axis=1)


    def get_temperature(self) -> pd.DataFrame:
        """
        Get the temperature in Kelvin for all temperature baths for every time step

        Returns
        -------
        pd.DataFrame
            pandas dataframe with all temperatures
        """
        return self.database["temperature"].apply(lambda x: x[:,0])


    def get_Hvap(self, gas_traj, nMolecules=1, temperature=None) -> float:
        gas_nonbonded_energy = 0
        if type(gas_traj) == type(self):
            gas_nonbonded_energy = gas_traj.get_totals_nonbonded().mean()
        elif type(gas_traj) == float:
            gas_nonbonded_energy = gas_traj
        else:
            raise TypeError("Did not understand the type of gas. Allowed are float (E_gas) or Tre (gas_trajectory)")
        liquid_nonbonded_energyself = self.get_totals_nonbonded().mean()

        # calculate heat of vaporization
        self.heat_vap = ea.get_Hvap(liq=liquid_nonbonded_energyself, gas=gas_nonbonded_energy, temperature=temperature, nMolecules=nMolecules)
        return self.heat_vap
