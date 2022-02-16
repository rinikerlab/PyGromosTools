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
from typing import Tuple, Dict

import pygromos.files.trajectory._general_trajectory as traj
from pygromos.files.trajectory.tre_field_libs.ene_fields import gromos_2020_tre_block_names_table
from pygromos.analysis import energy_analysis as ea

class Tre(traj._General_Trajectory):
    _gromos_file_ending:str = "tre"
    contributions_names:Tuple[str] = ("VdW", "Culomb", "Special", "?")

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
        self.totals_total = self.database["totals"].apply(lambda x: x[self.tre_block_name_table.totals_subblock_names.index('totene')])
        return self.totals_total

    def get_totals_totkin(self) -> pd.DataFrame:
        self.totals_totkin = self.database["totals"].apply(lambda x: x[self.tre_block_name_table.totals_subblock_names.index('totkin')])
        return self.totals_totkin
    
    def get_totals_totpot(self) -> pd.DataFrame:
        self.totals_totpot= self.database["totals"].apply(lambda x: x[self.tre_block_name_table.totals_subblock_names.index('totpot')])
        return self.totals_totpot

    def get_totals_totcov(self) -> pd.DataFrame:
        self.totals_totcov = self.database["totals"].apply(lambda x: x[self.tre_block_name_table.totals_subblock_names.index('totcov')])
        return self.totals_totcov

    def get_totals_totbonded(self) -> pd.DataFrame:
        self.totals_totbonded = self.database["totals"].apply(lambda x: x[self.tre_block_name_table.totals_subblock_names.index('totbond')])
        return self.totals_totbonded

    def get_totals_totangle(self) -> pd.DataFrame:
        self.totals_totangle = self.database["totals"].apply(lambda x: x[self.tre_block_name_table.totals_subblock_names.index('totangle')])
        return self.totals_totangle

    def get_totals_totdihedral(self) -> pd.DataFrame:
        self.totals_totdihedral = self.database["totals"].apply(lambda x: x[self.tre_block_name_table.totals_subblock_names.index('totdihedral')])
        return self.totals_totdihedral

    def get_totals_totnonbonded(self) -> pd.DataFrame:
        self.totals_totnonbonded = self.database["totals"].apply(lambda x: x[self.tre_block_name_table.totals_subblock_names.index('totnonbonded')])
        return self.totals_totnonbonded

    def get_totals_totlj(self) -> pd.DataFrame:
        self.totals_totlj = self.database["totals"].apply(lambda x: x[self.tre_block_name_table.totals_subblock_names.index('totlj')])
        return self.totals_totlj

    def get_totals_totcrf(self) -> pd.DataFrame:
        self.totals_totcrf = self.database["totals"].apply(lambda x: x[self.tre_block_name_table.totals_subblock_names.index('totcrfs')])
        return self.totals_totcrf

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

    def get_nonbondedForceGroupContributions(self)->Dict[int, Dict[int, pd.DataFrame]]:
        """
            This function returns a nice formatted dictionary for the nonbonded Contributions according to the Force groups of the tre file.

        Returns
        -------
        Dict[int, Dict[int, pd.DataFrame]]
            The dictionary containing the nonbonded contributions of the single ForceGroups with each other.
            Dict[ForceGroupI, Dict[ForceGroupJ, NonbondedEnergyContribs]]

        Raises
        ------
        ValueError
            _description_
        """
        #Check contibution_dimensionalities
        nFFContributions = self.database.nonbonded[0].shape[1]
        if(nFFContributions!=len(self.contributions_names)):
            raise ValueError("The dimensionality of the NonbondedContributions is not corresponding to the expected dimensionality.\n expected: "+str(len(contributions_names))+" \t found: "+str(nFFContributions))
        
        #Get the number of force groups:
        nForceGroup = self._get_numberOfForceGroupsFromNonbondeds()
        
        #Generate dictionary for the different contributions
        t = 0
        forceGroupNonbondedContributions = {}
        for i in range(1, 1+nForceGroup):
            forceGroupNonbondedContributions[i]={}
            for j in range(1, 1+nForceGroup):
                forceGroupNonbondedContributions[i][j]=pd.DataFrame(list(energy_traj.database.nonbonded.apply(lambda x: x[i+j-2])), columns=self.contributions_names)
                t+=1
        
        return forceGroupNonbondedContributions
    
    def _get_numberOfForceGroupsFromNonbondeds(self)->int:
        """
            This function gets the number of Force groups in the simulation from the nonbonded block.
        
        Returns
        -------
        int
            number of ForceGroups used for this tre.
        """
        nForceGroupDim = self.database.nonbonded[0].shape[0]

        r=0
        for nForceGroup in range(1,nForceGroupDim+1):
            r+=nForceGroup
            if(r == nForceGroupDims):
                return nForceGroup
            elif(r>nForceGroupDims):
                raise Exception("That should not happen!")
            
        
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

        tmps = []
        for i, row in self.database.iterrows():
            row_d = {"TIMESTEP_time": row["TIMESTEP_time"], "TIMESTEP_step": row["TIMESTEP_step"]}
            row_d.update({"bath" + str(i + 1): temp for i, temp in enumerate(row["temperature"][:, 0])})
            tmps.append(row_d)

        return pd.DataFrame(tmps)


    def get_Hvap(self, gas_traj, nMolecules=1, temperature=None) -> float:
        gas_totpot_energy = 0
        if type(gas_traj) == type(self):
            gas_totpot_energy = gas_traj.get_totals_totpot().mean()
        elif type(gas_traj) == float:
            gas_totpot_energy = gas_traj
        else:
            raise TypeError("Did not understand the type of gas. Allowed are float (E_gas) or Tre (gas_trajectory)")
        liq_totpot_energy = self.get_totals_totpot().mean()

        #get temperature from liq trajectory if not given
        if temperature is None:
            temperature = float(self.get_temperature().mean()[0])

        # calculate heat of vaporization
        self.heat_vap = ea.get_Hvap(liq_totpot_energy=liq_totpot_energy, gas_totpot_energy=gas_totpot_energy, temperature=temperature, nMolecules=nMolecules)
        return self.heat_vap
