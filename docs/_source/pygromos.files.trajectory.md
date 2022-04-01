---
date: '2022-04-01T13:55:18.499Z'
docname: _source/pygromos.files.trajectory
images: {}
path: /source-pygromos-files-trajectory
title: pygromos.files.trajectory package
---

# pygromos.files.trajectory package

## Subpackages


* [pygromos.files.trajectory.blocks package]()


    * [Submodules](#submodules)


    * [pygromos.files.trajectory.blocks.energy_trajectory_subblock module](#module-pygromos.files.trajectory.blocks.energy_trajectory_subblock)


    * [pygromos.files.trajectory.blocks.trajectory_blocks module](#module-pygromos.files.trajectory.blocks.trajectory_blocks)


    * [Module contents](#module-pygromos.files.trajectory.blocks)


* [pygromos.files.trajectory.tre_field_libs package]()


    * [Submodules](#submodules)


    * [pygromos.files.trajectory.tre_field_libs.ene_fields module](#module-pygromos.files.trajectory.tre_field_libs.ene_fields)


    * [Module contents](#module-pygromos.files.trajectory.tre_field_libs)


## Submodules

## pygromos.files.trajectory.trc module

File:            Class for trc files in pandas
Description:

> The pandas trajectory TRC class offers a easy method to process GROMOS’s .trc files in python
> The trc files are parsed into an easy to use pandas dataframe

Author: Marc Thierry Lehner

TODO: add stride option to all member functions
TODO: add support for periodic boundary condition
TODO: add support for rdkit mol selector
TODO: add support for rdkit conformers


### _class_ pygromos.files.trajectory.trc.Trc()
Bases: `mdtraj.core.trajectory.Trajectory`


#### TITLE(_: [pygromos.files.blocks._general_blocks.TITLE](#pygromos.files.blocks.pertubation_blocks.TITLE_ )

#### distances(atom_pairs: List[Tuple[int, int]], periodic: bool = True, opt: bool = True)

#### generate_TITLE_entry()

#### generate_entry_for_frame(frame_id: int)

#### get_dummy_cnf(xyz: numpy.array)

#### _classmethod_ load(in_path: str, in_cnf_path: Optional[str] = None, timestep_duration: float = 0.002)
Load a trajectory from disk


* **Parameters**

    **filenames** (*{path-like, [path-like]}*) – Either a path or list of paths



* **Other Parameters**

    **As requested by the various load functions – it depends on the extension**



#### parse_trc_efficiently(traj_path: str)

#### path(_: st_ )

#### recreate_view()

#### rmsd(reference_frame: int = 0, reference: Optional[mdtraj.core.trajectory.Trajectory] = None)

#### save(out_path: str)
Save trajectory to disk, in a format determined by the filename extension


* **Parameters**

    **filename** (*path-like*) – filesystem path in which to save the trajectory. The extension will
    be parsed and will control the format.



* **Other Parameters**

    
    * **lossy** (*bool*) – For .h5 or .lh5, whether or not to use compression.


    * **no_models** (*bool*) – For .pdb. TODO: Document this?


    * **force_overwrite** (*bool*) – If filename already exists, overwrite it.



#### _property_ step(_: numpy.arra_ )

#### to_cnf(frame_id: Optional[int] = None, base_cnf: Optional[[pygromos.files.coord.cnf.Cnf](#pygromos.files.coord.cnf.Cnf)] = None)

#### _property_ view(_: nglview.widget.NGLWidge_ )

#### write(out_path: str)

#### write_trc(out_path: str)
## pygromos.files.trajectory.tre module

File:            Class for tre files in pandas
Description:

> The pandas trajectory TRE class offers a easy method to process GROMOS’s .tre files in python
> The tre files are parsed into an easy to use pandas dataframe.

> This class should be a alternative for the data post processing with ene_ana in gromos++

Author: Marc Thierry Lehner

TODO: add stride option to all member functions
TODO: add support for periodic boundary condition

TODO: add ene_ana functions


### _class_ pygromos.files.trajectory.tre.Tre(input_value: str, auto_save: bool = True, stride: int = 1, skip: int = 0, _ene_ana_names: pygromos.files.trajectory.tre_field_libs.ene_fields.gromos_tre_block_names_table = <class 'pygromos.files.trajectory.tre_field_libs.ene_fields.gromos_2021_tre_block_names_table'>)
Bases: `pygromos.files.trajectory._general_trajectory._General_Trajectory`

The Tre files are results from Gromos simulations, that store all the calculated energies and properties during the simulation.


#### \__init__(input_value: str, auto_save: bool = True, stride: int = 1, skip: int = 0, _ene_ana_names: pygromos.files.trajectory.tre_field_libs.ene_fields.gromos_tre_block_names_table = <class 'pygromos.files.trajectory.tre_field_libs.ene_fields.gromos_2021_tre_block_names_table'>)
> Build a Gromos energy trajectory file (.tre)


* **Parameters**

    
    * **input_value** (*str,None*) – The input value can be None, or a string path to the .tre/.tre.gz file.


    * **auto_save** (*bool, optional*) – automatically save the file, by default True


    * **stride** (*int, optional*) – only read every x value, by default 1


    * **skip** (*int, optional*) – skip the first x timesteps, by default 0


    * **_ene_ana_names** (*gromos_tre_block_names_table, optional*) – get the field names after the provided standard., by default gromos_2020_tre_block_names_table



#### _get_numberOfForceGroupsFromNonbondeds()
> This function gets the number of Force groups in the simulation from the nonbonded block.


* **Returns**

    number of ForceGroups used for this tre.



* **Return type**

    int



#### _set_data(attibute_name: str, rows_name: str, field_names: Tuple[str])
_summary_

    This function extracts generially the information of a column per time


* **Parameters**

    
    * **attibute_name** (*str*) – name of the target attribute


    * **rows_name** (*str*) – name of the block, that shall be extracted


    * **field_names** (*Tuple[str]*) – name of the fields in each row



* **Returns**

    contains the extracted information



* **Return type**

    pd.DataFrame



#### get_Hvap(gas_traj: pygromos.utils.typing.Tre_Type, nMolecules: int = 1, temperature: Optional[float] = None)

#### get_baths()
extract data of the baths block


#### get_bondedContributions()
extract data of the bonded block


#### get_density()
Calculate the density for every frame.
(Uses mass and the first volume entry)


* **Returns**

    Dataframe with the densities for all time steps



* **Return type**

    pd.DataFrame



#### get_eds()
> Get EDS energies if present.


* **Returns**

    returns datafrae with columns for each endstate.



* **Return type**

    pd.DataFrame



#### get_mass()
> returns the systems mass per timestep


* **Returns**

    series of mass per time



* **Return type**

    pd.Series



#### get_nonbondedContributions()
> This function returns a nice formatted dictionary for the nonbonded Contributions according to the Force groups of the tre file.


* **Returns**

    The dictionary containing the nonbonded contributions of the single ForceGroups with each other.
    Dict[ForceGroupI, Dict[ForceGroupJ, NonbondedEnergyContribs]]



* **Return type**

    Dict[int, Dict[int, pd.DataFrame]]



* **Raises**

    **ValueError** – returns Value error, if the dimensionality of the different contributions does not fit to the _nonbonded_contribution_names.



#### get_precalclam()
> Get the energies calculated for the different defined lambda values in a trajectory.


* **Returns**

    return the energies calculated for the different lambda values.



* **Return type**

    pd.DataFrame



#### get_specialContributions()
extract data of the special block


#### get_temperature()
Get the temperature in Kelvin for all temperature baths for every time step


* **Returns**

    pandas dataframe with all temperatures



* **Return type**

    pd.DataFrame



#### get_temperature_Info()
> temperature baths


* **Returns**

    returns the full info of the temperature baths per bath



* **Return type**

    Dict[int,pd.DataFrame]



#### get_totals()
get all totals of the system


#### get_totangle()
get the total angle contribution/ per time


#### get_totbonded()
get the total bonded contribution/ per time


#### get_totcov()
get the total covalent contribution/ per time


#### get_totcrf()
get the total columbic reactionfield contribution/ per time


#### get_totdihedral()
get the total dihedral contribution/ per time


#### get_totene()
get the total System energy / per time


#### get_totkin()
get the total kinetic Energy / per time


#### get_totlj()
get the total lennard jones contribution/ per time


#### get_totnonbonded()
get the total nonbonded contribution/ per time


#### get_totpot()
get the total potential Energy / per time

## pygromos.files.trajectory.trg module

File:            Class for tre files in pandas
Description:

> The pandas trajectory TRE class offers a easy method to process GROMOS’s .trg files in python
> The tre files are parsed into an easy to use pandas dataframe.

> This class should be a alternative for the data post processing with ene_ana in gromos++

Author:  Marc Thierry Lehner & Benjamin Ries

TODO: add stride option to all member functions
TODO: add support for periodic boundary condition

TODO: add ene_ana functions


### _class_ pygromos.files.trajectory.trg.Trg(input_value: str, auto_save=True, stride: int = 1, skip: int = 0)
Bases: `pygromos.files.trajectory._general_trajectory._General_Trajectory`


#### get_lambdas()

#### get_precalclam()

#### get_totals()

### _class_ pygromos.files.trajectory.trg.gromos_2020_trg_block_names_table()
Bases: `object`


#### precalclam_subblock(_ = ['nr_lambdas', 'A_e_lj', 'B_e_lj', 'A_e_crf', 'B_e_crf', 'AB_kinetic', 'AB_bond', 'AB_angle', 'AB_improper', 'AB_disres', 'AB_dihres', 'AB_disfld'_ )

#### totals_subblock_names(_ = ['dHdl', 'dKdl', 'dVdl', 'WIP', 'WIP', 'WIP', 'WIP', 'WIP', 'WIP', 'WIP', 'WIP', 'WIP', 'WIP', 'WIP', 'WIP', 'WIP', 'WIP', 'WIP', 'WIP', 'WIP', 'WIP', 'WIP', 'WIP', 'WIP', 'WIP', 'WIP', 'WIP', 'WIP', 'WIP', 'WIP', 'WIP', 'WIP', 'WIP', 'WIP', 'WIP', 'WIP', 'WIP', 'WIP', 'WIP', 'WIP', 'WIP', 'WIP', 'WIP'_ )
## Module contents
