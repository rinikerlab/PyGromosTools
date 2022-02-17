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
    _contributions_nonbonded_names: Tuple[str] = ("Lennard-Jones", "Coulomb/RF", "lattice sum real",  "lattice sum reciproc")
    _contributions_bonded_names: Tuple[str] = ("bond", "angle", "improper", "dihedral", "crossdihedral")
    _contributions_temp_baths: Tuple[str] = ("kinetic total", "centre of mass", "internal/rotational")
    _contributions_special: Tuple[str]  = ("constraints", "pos. restraints", "dist. restraints", "disfield res", "dihe. restr.", "SASA", "SASA volume","jvalue","rdc","local elevation", "path integral", "angle restraint")

    def __init__(self, input_value: str or None, auto_save=True, stride:int=1, skip:int=0, _ene_ana_names = gromos_2020_tre_block_names_table):
        super().__init__(input_value, auto_save=auto_save, stride=stride, skip=skip)
        self.tre_block_name_table = _ene_ana_names


    """--------------------------------------------------------------
    Basic getters for Subblocks and Subsubblocks
    -----------------------------------------------------------------
    """

    def get_totals(self) -> pd.DataFrame:
        #print(self.database["totals"][0].shape, self.database["totals"][0])
        if(not hasattr(self, "totals")):
            self.totals = pd.DataFrame(data =np.stack(self.database["totals"].to_numpy()), index=self.database.time, columns=self.tre_block_name_table.totals_subblock_names)
        else:
            pass
        return self.totals

    def get_totene(self) -> pd.DataFrame:
        self.totene = self.get_totals()['totene']
        return self.totene

    def get_totkin(self) -> pd.DataFrame:
        self.totkin = self.get_totals()['totkin']
        return self.totkin

    def get_totpot(self) -> pd.DataFrame:
        self.totpot = self.get_totals()['totpot']
        return self.totpot

    def get_totcov(self) -> pd.DataFrame:
        self.totcov = self.get_totals()['totcov']
        return self.totcov

    def get_totbonded(self) -> pd.DataFrame:
        self.totbonded = self.get_totals()['totbond']
        return self.totbonded

    def get_totangle(self) -> pd.DataFrame:
        self.totangle = self.get_totals()['totangle']
        return self.totangle

    def get_totdihedral(self) -> pd.DataFrame:
        self.totdihedral = self.get_totals()['totdihedral']
        return self.totdihedral

    def get_totnonbonded(self) -> pd.DataFrame:
        self.totnonbonded = self.get_totals()['totnonbonded']
        return self.totnonbonded

    def get_totals_totlj(self) -> pd.DataFrame:
        self.totlj = self.get_totals()['totlj']
        return self.totlj

    def get_totals_totcrf(self) -> pd.DataFrame:
        self.totcrf = self.get_totals()['totcrf']
        return self.totcrf

    def get_eds(self)->pd.DataFrame:
        """
            Get EDS energies if present.

        Returns
        -------
        pd.DataFrame
            returns datafrae with columns for each endstate.
        """
        if(not hasattr(self, "eds")):
            num_states = self.database["eds"][0][0]
            if(num_states>0):
                state_strings =[[column+"_"+str(state_row) for column in self.tre_block_name_table.eds_subblock_names_singleState] for state_row in range(1, 1+int(num_states))]
                self.tre_block_name_table.eds_subblock_names = ["numstates"]
                self.tre_block_name_table.eds_subblock_names.extend(list(np.concatenate(state_strings)))
            else:
                self.tre_block_name_table.eds_subblock_names = ["numstates"]
            self.eds = pd.DataFrame(data=np.stack(self.database["eds"].to_numpy()), index=self.database.time, columns=self.tre_block_name_table.eds_subblock_names)
        else:
            pass
        
        return self.eds

    def get_precalclam(self)->pd.DataFrame:
        """
            Get the energies calculated for the different defined lambda values in a trajectory.

        Returns
        -------
        pd.DataFrame
            return the energies calculated for the different lambda values.
        """
        if(not hasattr(self, "precalclam")):
            num_lam = self.database["precalclam"][0][0]
            if(num_lam>0):
                state_strings =[[column+"_"+str(state_row) for column in self.tre_block_name_table.lam_subblock_names_singleLam] for state_row in range(1, 1+int(num_lam))]
                self.tre_block_name_table.lam_subblock_names = ["nr_lambdas"]
                self.tre_block_name_table.lam_subblock_names.extend(list(np.concatenate(state_strings)))
                self.tre_block_name_table.lam_subblock_names.extend(["A_dihedral", "B_dihedral"])
            else:
                self.tre_block_name_table.lam_subblock_names = ["nr_lambdas", "A_dihedral", "B_dihedral"]

            self.precalclam = pd.DataFrame(data=np.stack(self.database["precalclam"].to_numpy()), index=self.database.time, columns=self.tre_block_name_table.lam_subblock_names)
        else:
            pass
        
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
            returns Value error, if the dimensionality of the different contributions does not fit to the _nonbonded_contribution_names.
        """
        if(not hasattr(self, "forceGroupNonbondedContributions")):
            #Check contibution_dimensionalities
            nFFContributions = self.database.nonbonded[0].shape[1]
            if(nFFContributions!=len(self._contributions_nonbonded_names)):
                raise ValueError("The dimensionality of the NonbondedContributions is not corresponding to the expected dimensionality.\n expected: "+str(len(contributions_names))+" \t found: "+str(nFFContributions))
            
            #Get the number of force groups:
            nForceGroup = self._get_numberOfForceGroupsFromNonbondeds()
            
            #Generate dictionary for the different contributions
            t = 0
            forceGroupNonbondedContributions = {}
            for i in range(1, 1+nForceGroup):
                forceGroupNonbondedContributions[i]={}
                for j in range(1, 1+nForceGroup):
                    forceGroupNonbondedContributions[i][j]=pd.DataFrame(list(self.database.nonbonded.apply(lambda x: x[i+j-2])), columns=self.contributions_names)
                    t+=1
            self.forceGroupNonbondedContributions = forceGroupNonbondedContributions
        return self.forceGroupNonbondedContributions

    def get_bondedContributions(self)->Dict[int,pd.DataFrame]:
        if(not hasattr(self, "bondedContributions")):
            self.bondedContributions = pd.DataFrame(list(map(lambda x: list(x[0]), list(self.database.bonded))),
                        columns=self._contributions_bonded_names,
                        index=self.database.time)
        return self.bondedContributions

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
            if(r == nForceGroupDim):
                return nForceGroup
            elif(r>nForceGroupDim):
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
        return pd.Series(list(self.database[["mass","volume"]].apply(lambda x: ea.get_density(mass=x[0][0], volume=x[1][0]), axis=1)), 
                         index=self.database.time, name="density")



    def get_temperature(self) -> pd.DataFrame:  #CHECK THIS - tempcontrib:total             com               ir                scaling factor
        """
        Get the temperature in Kelvin for all temperature baths for every time step

        Returns
        -------
        pd.DataFrame
            pandas dataframe with all temperatures
        """

        tmps = []
        for i, row in self.database.iterrows():
            row_d = {"bath" + str(i + 1): temp for i, temp in enumerate(row["temperature"][:, 0])}
            tmps.append(row_d)

        return pd.DataFrame(tmps, index=self.database.time)

    def get_baths(self)->pd.DataFrame:  #CHECK THIS
        _contributions_temp_baths

    def get_special(self)->pd.DataFrame:    #CHECK THIS
        if(not hasattr(self, "specialContributions")):
            self.specialContributions = pd.DataFrame(list(map(lambda x: list(x[0]), list(self.database.special))),
                        columns=self._contributions_special,
                        index=self.database.time)
        return self.specialContributions
        

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
