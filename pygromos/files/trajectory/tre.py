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

# imports
import pandas as pd
import numpy as np

import pygromos.files.trajectory._general_trajectory as traj
from pygromos.files.trajectory.tre_field_libs.ene_fields import (
    gromos_2021_tre_block_names_table,
    gromos_tre_block_names_table,
)
from pygromos.utils.typing import Tuple, Dict, Tre_Type
from pygromos.analysis import energy_analysis as ea


class Tre(traj._General_Trajectory):
    """
    The Tre files are results from Gromos simulations, that store all the calculated energies and properties during the simulation.

    """

    _gromos_file_ending: str = "tre"
    _contributions_nonbonded_names: Tuple[str] = (
        "Lennard-Jones",
        "Coulomb/RF",
        "lattice sum real",
        "lattice sum reciproc",
    )
    _contributions_bonded_names: Tuple[str] = ("bond", "angle", "improper", "dihedral", "crossdihedral")
    _contributions_baths: Tuple[str] = ("kinetic total", "centre of mass", "internal/rotational")
    _contributions_special: Tuple[str] = (
        "constraints",
        "pos. restraints",
        "dist. restraints",
        "disfield res",
        "dihe. restr.",
        "SASA",
        "SASA volume",
        "jvalue",
        "rdc",
        "local elevation",
        "path integral",
        "angle restraint",
    )
    _contributions_temperature: Tuple[str] = ("total", "com", "ir", "scaling factor")

    def __init__(
        self,
        input_value: str,
        auto_save: bool = True,
        stride: int = 1,
        skip: int = 0,
        _ene_ana_names: gromos_tre_block_names_table = gromos_2021_tre_block_names_table,
    ):
        """
            Build a Gromos energy trajectory file (.tre)

        Parameters
        ----------
        input_value : str,None
            The input value can be None, or a string path to the .tre/.tre.gz file.
        auto_save : bool, optional
            automatically save the file, by default True
        stride : int, optional
            only read every x value, by default 1
        skip : int, optional
            skip the first x timesteps, by default 0
        _ene_ana_names : gromos_tre_block_names_table, optional
            get the field names after the provided standard., by default gromos_2020_tre_block_names_table
        """
        super().__init__(input_value, auto_save=auto_save, stride=stride, skip=skip)
        self.tre_block_name_table = _ene_ana_names

        if "time" in self.database.columns:
            self.time = self.database.time
        if "step" in self.database.columns:
            self.step = self.database.step

    """--------------------------------------------------------------
    Basic getters for Subblocks and Subsubblocks
    -----------------------------------------------------------------

        The following functions, return well formated values from the tre trajectories.


      ENERGY03 - fields:
    """

    def get_totals(self) -> pd.DataFrame:
        """
        get all totals of the system
        """
        # print(self.database["totals"][0].shape, self.database["totals"][0])
        if not hasattr(self, "totals"):
            totals_data = np.stack(self.database["totals"].to_numpy())
            self.totals = pd.DataFrame(
                data=totals_data,
                index=self.database.time,
                columns=self.tre_block_name_table.totals_subblock_names[: totals_data.shape[1]],
            )
        else:
            pass
        return self.totals

    def get_totene(self) -> pd.DataFrame:
        """
        get the total System energy / per time
        """
        self.totene = self.get_totals()["totene"]
        return self.totene

    def get_totkin(self) -> pd.DataFrame:
        """
        get the total kinetic Energy / per time
        """
        self.totkin = self.get_totals()["totkin"]
        return self.totkin

    def get_totpot(self) -> pd.DataFrame:
        """
        get the total potential Energy / per time
        """
        self.totpot = self.get_totals()["totpot"]
        return self.totpot

    def get_totcov(self) -> pd.DataFrame:
        """
        get the total covalent contribution/ per time
        """
        self.totcov = self.get_totals()["totcov"]
        return self.totcov

    def get_totbonded(self) -> pd.DataFrame:
        """
        get the total bonded contribution/ per time
        """
        self.totbonded = self.get_totals()["totbond"]
        return self.totbonded

    def get_totangle(self) -> pd.DataFrame:
        """
        get the total angle contribution/ per time
        """
        self.totangle = self.get_totals()["totangle"]
        return self.totangle

    def get_totdihedral(self) -> pd.DataFrame:
        """
        get the total dihedral contribution/ per time
        """
        self.totdihedral = self.get_totals()["totdihedral"]
        return self.totdihedral

    def get_totnonbonded(self) -> pd.DataFrame:
        """
        get the total nonbonded contribution/ per time
        """
        self.totnonbonded = self.get_totals()["totnonbonded"]
        return self.totnonbonded

    def get_totlj(self) -> pd.DataFrame:
        """
        get the total lennard jones contribution/ per time
        """
        self.totlj = self.get_totals()["totlj"]
        return self.totlj

    def get_totcrf(self) -> pd.DataFrame:
        """
        get the total columbic reactionfield contribution/ per time
        """
        self.totcrf = self.get_totals()["totcrf"]
        return self.totcrf

    def get_baths(self) -> pd.DataFrame:
        """
        extract data of the baths block
        """
        return self._set_data(attibute_name="baths", rows_name="baths", field_names=self._contributions_baths)

    def get_bondedContributions(self) -> Dict[int, pd.DataFrame]:
        """
        extract data of the bonded block
        """
        return self._set_data(
            attibute_name="bondedContributions", rows_name="bonded", field_names=self._contributions_bonded_names
        )

    def get_nonbondedContributions(self) -> Dict[int, Dict[int, pd.DataFrame]]:
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
        if not hasattr(self, "forceGroupNonbondedContributions"):
            # Check contibution_dimensionalities
            nFFContributions = self.database.nonbonded[0].shape[1]
            if nFFContributions != len(self._contributions_nonbonded_names):
                raise ValueError(
                    "The dimensionality of the NonbondedContributions is not corresponding to the expected dimensionality.\n expected: "
                    + str(len(self._contributions_nonbonded_names))
                    + " \t found: "
                    + str(nFFContributions)
                )

            # Get the number of force groups:
            nForceGroup = self._get_numberOfForceGroupsFromNonbondeds()

            # Generate dictionary for the different contributions
            t = 0
            forceGroupNonbondedContributions = {}
            for i in range(1, 1 + nForceGroup):
                forceGroupNonbondedContributions[i] = {}
                for j in range(1, 1 + nForceGroup):
                    forceGroupNonbondedContributions[i][j] = pd.DataFrame(
                        list(self.database.nonbonded.apply(lambda x: x[i + j - 2])),
                        columns=self._contributions_nonbonded_names,
                    )
                    t += 1
            self.forceGroupNonbondedContributions = forceGroupNonbondedContributions
        return self.forceGroupNonbondedContributions

    def get_specialContributions(self) -> Dict[int, pd.DataFrame]:  # CHECK THIS
        """
        extract data of the special block
        """
        return self._set_data(
            attibute_name="specialContributions", rows_name="special", field_names=self._contributions_special
        )

    def get_eds(self) -> pd.DataFrame:
        """
            Get EDS energies if present.

        Returns
        -------
        pd.DataFrame
            returns datafrae with columns for each endstate.
        """
        if not hasattr(self, "eds"):
            num_states = self.database["eds"][0][0]
            if num_states > 0:
                state_strings = [
                    [
                        column + "_" + str(state_row)
                        for column in self.tre_block_name_table.eds_subblock_names_singleState
                    ]
                    for state_row in range(1, 1 + int(num_states))
                ]
                self.tre_block_name_table.eds_subblock_names = ["numstates"]
                self.tre_block_name_table.eds_subblock_names.extend(list(np.concatenate(state_strings)))
            else:
                self.tre_block_name_table.eds_subblock_names = ["numstates"]
            self.eds = pd.DataFrame(
                data=np.stack(self.database["eds"].to_numpy()),
                index=self.database.time,
                columns=self.tre_block_name_table.eds_subblock_names,
            )
        else:
            pass

        return self.eds

    def get_precalclam(self) -> pd.DataFrame:
        """
            Get the energies calculated for the different defined lambda values in a trajectory.

        Returns
        -------
        pd.DataFrame
            return the energies calculated for the different lambda values.
        """
        if not hasattr(self, "precalclam"):
            num_lam = self.database["precalclam"][0][0]
            if num_lam > 0:
                state_strings = [
                    [column + "_" + str(state_row) for column in self.tre_block_name_table.lam_subblock_names_singleLam]
                    for state_row in range(1, 1 + int(num_lam))
                ]
                self.tre_block_name_table.lam_subblock_names = ["nr_lambdas"]
                self.tre_block_name_table.lam_subblock_names.extend(list(np.concatenate(state_strings)))
                self.tre_block_name_table.lam_subblock_names.extend(["A_dihedral", "B_dihedral"])
            else:
                self.tre_block_name_table.lam_subblock_names = ["nr_lambdas", "A_dihedral", "B_dihedral"]

            self.precalclam = pd.DataFrame(
                data=np.stack(self.database["precalclam"].to_numpy()),
                index=self.database.time,
                columns=self.tre_block_name_table.lam_subblock_names,
            )
        else:
            pass

        return self.precalclam

    """
      VOLUMEPRESSURE03 - fields:
    """

    def get_mass(self) -> pd.Series:
        """
            returns the systems mass per timestep

        Returns
        -------
        pd.Series
            series of mass per time
        """
        return pd.Series(map(float, self.database.mass), name="mass", index=self.database.time)

    def get_temperature_Info(self) -> Dict[int, pd.DataFrame]:
        """
            temperature baths

        Returns
        -------
        Dict[int,pd.DataFrame]
            returns the full info of the temperature baths per bath
        """
        return self._set_data(
            attibute_name="temperatureInfo", rows_name="temperature", field_names=self._contributions_temperature
        )

    def get_temperature(
        self,
    ) -> pd.DataFrame:  # CHECK THIS - tempcontrib:total             com               ir                scaling factor
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

    """
      UTILS
    """

    def _set_data(self, attibute_name: str, rows_name: str, field_names: Tuple[str]) -> pd.DataFrame:
        """_summary_
            This function extracts generially the information of a column per time

        Parameters
        ----------
        attibute_name : str
            name of the target attribute
        rows_name : str
            name of the block, that shall be extracted
        field_names : Tuple[str]
            name of the fields in each row

        Returns
        -------
        pd.DataFrame
            contains the extracted information
        """
        if not hasattr(self, attibute_name):
            nDimGroups = self.database[rows_name][0].shape[0]

            setattr(self, attibute_name, {})
            for i in range(nDimGroups):
                getattr(self, attibute_name).update(
                    {
                        i: pd.DataFrame(
                            list(map(lambda x: list(x[i]), list(self.database[rows_name]))),
                            columns=field_names,
                            index=self.database.time,
                        )
                    }
                )
        return getattr(self, attibute_name)

    def _get_numberOfForceGroupsFromNonbondeds(self) -> int:
        """
            This function gets the number of Force groups in the simulation from the nonbonded block.

        Returns
        -------
        int
            number of ForceGroups used for this tre.
        """
        nForceGroupDim = self.database.nonbonded[0].shape[0]

        r = 0
        for nForceGroup in range(1, nForceGroupDim + 1):
            r += nForceGroup
            if r == nForceGroupDim:
                return nForceGroup
            elif r > nForceGroupDim:
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
        return pd.Series(
            list(
                self.database[["mass", "volume"]].apply(lambda x: ea.get_density(mass=x[0][0], volume=x[1][0]), axis=1)
            ),
            index=self.database.time,
            name="density",
        )

    def get_Hvap(self, gas_traj: Tre_Type, nMolecules: int = 1, temperature: float = None) -> float:
        gas_totpot_energy = 0
        if type(gas_traj) == type(self):
            gas_totpot_energy = gas_traj.get_totpot().mean()
        elif type(gas_traj) == float:
            gas_totpot_energy = gas_traj
        else:
            raise TypeError("Did not understand the type of gas. Allowed are float (E_gas) or Tre (gas_trajectory)")
        liq_totpot_energy = self.get_totpot().mean()

        # get temperature from liq trajectory if not given
        if temperature is None:
            temperature = float(self.get_temperature().mean()[0])

        # calculate heat of vaporization
        self.heat_vap = ea.get_Hvap(
            liq_totpot_energy=liq_totpot_energy,
            gas_totpot_energy=gas_totpot_energy,
            temperature=temperature,
            nMolecules=nMolecules,
        )
        return self.heat_vap
