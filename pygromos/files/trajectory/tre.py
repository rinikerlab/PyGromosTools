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

class gromos_2020_tre_block_names_table():
    totals_subblock_names = ["totene","totkin","totpot","totcov","totbond","totangle","totimproper","totdihedral","totcrossdihedral","totnonbonded",
                             "totlj","totcrf","totls","totlspair","totlsreal","totlsk","totlsa","totlsself","totlssurf","totpolself","totspecial",
                             "totsasa","totsasavol","totconstraint","totdisres","totdisfieldres","totdihres","totposres","totjval","totxray","totle",
                             "totorder","totsymm","eds_vr,entropy","totqm","totbsleus","totrdc","wip1","wip2","wip3","wip4","wip5","wip6"]


class Tre(traj._General_Trajectory):
    def __init__(self, input_value: str or None, auto_save=True, stride:int=1, skip:int=0):
        super().__init__(input_value, auto_save=auto_save, stride=stride, skip=skip)
        self.tre_block_name_table = gromos_2020_tre_block_names_table


    """--------------------------------------------------------------
    Basic getters for Subblocks and Subsubblocks
    -----------------------------------------------------------------
    """

    def get_totals(self) -> pd.DataFrame:
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
        return self.database[["mass","volume"]].apply(lambda x: 1.66056 * x[0][0]/x[1][0], axis=1)

    def get_temperature(self) -> pd.DataFrame:
        """
        Get the temperature in Kelvin for all temperature baths for every time step

        Returns
        -------
        pd.DataFrame
            pandas dataframe with all temperatures
        """
        return self.database["temperature"].apply(lambda x: x[:,0])

    def get_Hvap(self, gas, nMolecules=1, temperature=None):
        #get gas nonbonded energy from multiple different gas arguments
        gas_nonbonded_energy = 0
        if type(gas) == type(self):
            gas_nonbonded_energy = gas.get_totals_nonbonded().mean()
        elif type(gas) == float:
            gas_nonbonded_energy = gas
        else:
            raise TypeError("Did not understand the type of gas. Allowed are float (E_gas) or Tre (gas_trajectory)")

        # get liquid nonbonded energy
        liquid_nonbonded_energyself = self.get_totals_nonbonded().mean()

        # calculate heat of vaporization
        rt_constant = 0.008314462618153239 * temperature # R in kilojoule_per_mole/kelvin * T
        self.heat_vap = gas_nonbonded_energy - liquid_nonbonded_energyself/nMolecules + rt_constant
        return self.heat_vap