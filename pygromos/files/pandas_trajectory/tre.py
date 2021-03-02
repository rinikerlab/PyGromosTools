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
import pygromos.files.pandas_trajectory._general_trajectory as traj

class Tre(traj._General_Trajectory):
    def __init__(self, input: str or None):
        super().__init__(input)


    def get_density(self) -> pd.DataFrame:
        """
        Calculate the density for every frame.
        (Uses mass and the first volume entry)

        Returns
        -------
        pd.DataFrame
            Dataframe with the densities for all time steps
        """
        return self.database[["mass","volume"]].apply(lambda x: x[0][0]/x[1][0], axis=1)

    def get_temperature(self) -> pd.DataFrame:
        """
        Get the temperature in Kelvin for all temperature baths for every time step

        Returns
        -------
        pd.DataFrame
            pandas dataframe with all temperatures
        """
        return self.database["temperature"].apply(lambda x: x[:,0])