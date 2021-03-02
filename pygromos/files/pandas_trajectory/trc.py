"""
File:            Class for trc files in pandas
Description:
    The pandas trajectory TRC class offers a easy method to process GROMOS's .trc files in python
    The trc files are parsed into an easy to use pandas dataframe

Author: Marc Thierry Lehner

TODO: add stride option to all member functions
TODO: add support for periodic boundary condition
TODO: add support for rdkit mol selector
TODO: add support for rdkit conformers

"""

#imports
import pandas as pd
import numpy as np
import pygromos.files.pandas_trajectory._general_trajectory as traj

class Trc(traj._General_Trajectory):
    def __init__(self, input: str or None):
        super().__init__(input)


    def get_atom_pair_distance_series(self, atomI:int, atomJ:int, periodicBoundary=False) -> pd.DataFrame:
        """
        creates a series with the euclidian distance for the atom pair i-j for every time step

        Parameters
        ----------
        atomI : int
        atomJ : int
        periodicBoundary : bool, optional
            WIP

        Returns
        -------
        pd.DataFrame
        """
        return self.database[["POS_"+str(atomI),"POS_"+str(atomJ)]].apply(lambda x: np.sqrt(np.sum((x[1]-x[0])**2)), axis=1)


    def get_atom_pair_distance_mean(self, atomI:int, atomJ:int) -> float:
        """
        creates the mean of the euclidian distances for every timestep between the atoms i and j

        Parameters
        ----------
        atomI : int
        atomJ : int

        Returns
        -------
        float
        """
        return self.get_atom_pair_distance_series(atomI=atomI, atomJ=atomJ).mean()


    def get_atom_movement_series(self, atomI:int, periodicBoundary=False) -> pd.DataFrame:
        """
        for every time step the difference to the last coordinates is calculated for an atom i
        Returns the euclidian difference to the last step as a numpy array

        Parameters
        ----------
        atomI : int
        periodicBoundary : bool, optional
            WIP

        Returns
        -------
        pd.DataFrame
        """
        return (self.database["POS_"+str(atomI)]-self.database["POS_"+str(atomI)].shift())[1:]
    

    def get_atom_movement_length_series(self, atomI:int, periodicBoundary=False) -> pd.DataFrame:
        """
        for every timestep the total euclidian difference to the last time step is calculate for atom i

        Parameters
        ----------
        atomI : int
        periodicBoundary : bool, optional
            WIP

        Returns
        -------
        pd.DataFrame
        """
        return self.get_atom_movement_series(atomI=atomI, periodicBoundary=periodicBoundary).apply(lambda x: np.sqrt(np.sum(x**2)))


    def get_atom_movement_length_mean(self, atomI:int, periodicBoundary=False) -> float:
        """
        The average euclidian movement between to consecutive timesteps is calculated for atom i

        Parameters
        ----------
        atomI : int
        periodicBoundary : bool, optional
            WIP

        Returns
        -------
        float
        """
        return self.get_atom_movement_length_series(atomI=atomI, periodicBoundary=periodicBoundary).mean()


    def get_atom_movement_length_total(self, atomI:int, periodicBoundary=False) -> float:
        """
        The total euclidian movement between to consecutive timesteps is calculated for atom i

        Parameters
        ----------
        atomI : int
        periodicBoundary : bool, optional
            WIP

        Returns
        -------
        float
        """
        return self.get_atom_movement_length_series(atomI=atomI, periodicBoundary=periodicBoundary).sum()


    def get_cog_movement_vector_series_for_atom_group(self, atoms:list, periodicBoundary=False) -> pd.DataFrame:
        """
        Calculate the movemnt vector between to consecutive time steps for the center of geometry defined by a list of atoms

        Parameters
        ----------
        atoms : list[int]
            the list of atoms defining the center of geometry
        periodicBoundary : bool, optional
            WIP

        Returns
        -------
        pd.DataFrame
            numpy array with the movement vector for every time step
        """
        #create a local copy which can be modified
        temp_database = self.database.copy()
        # create the difference to the last time step for every atom i in the list atoms
        for i in atoms:
            temp_database["mix"+str(i)] = temp_database["POS_"+str(i)]-temp_database["POS_"+str(i)].shift()
        #calculate vector to cog
        return temp_database[["mix"+str(i) for i in atoms]][1:].sum(axis=1).apply(lambda x: x/len(atoms))

    def get_cog_movement_total_series_for_atom_group(self, atoms:list, periodicBoundary=False) -> pd.DataFrame:
        """
        Calculate the total movemnt between to consecutive time steps for the center of geometry defined by a list of atoms

        Parameters
        ----------
        atoms : list[int]
            the list of atoms defining the center of geometry
        periodicBoundary : bool, optional
            WIP

        Returns
        -------
        pd.DataFrame
            total movement (float) for every time step
        """
        return self.get_cog_movement_vector_series_for_atom_group(atoms=atoms, periodicBoundary=periodicBoundary).apply(lambda x: np.sqrt(np.sum(x**2)))