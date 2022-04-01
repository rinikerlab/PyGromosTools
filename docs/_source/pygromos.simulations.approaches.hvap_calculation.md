---
date: '2022-04-01T13:55:18.499Z'
docname: _source/pygromos.simulations.approaches.hvap_calculation
images: {}
path: /source-pygromos-simulations-approaches-hvap-calculation
title: pygromos.simulations.approaches.hvap_calculation package
---

# pygromos.simulations.approaches.hvap_calculation package

## Subpackages


* [pygromos.simulations.approaches.hvap_calculation.hvap_input_files package]()


    * [Module contents](#module-pygromos.simulations.approaches.hvap_calculation.hvap_input_files)


## Submodules

## pygromos.simulations.approaches.hvap_calculation.hvap_calculation module

File:            automatic calculation of Hvap
Warnings: this class is WIP!

Description:

    For a given gromos_system (or smiles) the heat of vaporization is automaticaly calculated.

    Main elements:


        1. create single molecule topo and conformation


        2. run gas (single molecule) minimization


        3. run gas SD simulation (and equilibaration)


        4. generate multi molecule (liquid) topo and conformation


        5. run liquid minimization


        6. run liquid MD run (and equilibaration)


        7. calculate Hvap from gas and liquid trajectories

Author: Marc Lehner


### _class_ pygromos.simulations.approaches.hvap_calculation.hvap_calculation.Hvap_calculation(input_system: pygromos.files.gromos_system.gromos_system.Gromos_System, work_folder: str, system_name: str = 'dummy', forcefield: pygromos.files.forcefield._generic_force_field._generic_force_field = <pygromos.files.forcefield._generic_force_field._generic_force_field object>, in_gromosXX_bin_dir: typing.Optional[str] = None, in_gromosPP_bin_dir: typing.Optional[str] = None, useGromosPlsPls: bool = True, verbose: bool = True)
Bases: `object`


#### \__init__(input_system: pygromos.files.gromos_system.gromos_system.Gromos_System, work_folder: str, system_name: str = 'dummy', forcefield: pygromos.files.forcefield._generic_force_field._generic_force_field = <pygromos.files.forcefield._generic_force_field._generic_force_field object>, in_gromosXX_bin_dir: typing.Optional[str] = None, in_gromosPP_bin_dir: typing.Optional[str] = None, useGromosPlsPls: bool = True, verbose: bool = True)
For a given gromos_system (or smiles) the heat of vaporization is automaticaly calculated


* **Parameters**

    **input_system** (*Gromos_SystemorstrorChem.rdchem.Mol*) – single molecule gromos_sytem or rdkit Molecule or SMILES



#### calc_hvap()

#### create_liq()

#### dens_modifier(_: floa_ _ = 0._ )

#### run()

#### run_gas()

#### run_liq()

#### submissonSystem_gas(_: <module 'pygromos.simulations.hpc_queuing.submission_systems._submission_system' from '/home/bschroed/Documents/projects/PyGromosTools/pygromos/simulations/hpc_queuing/submission_systems/_submission_system.py'_ )

#### submissonSystem_liq(_: <module 'pygromos.simulations.hpc_queuing.submission_systems._submission_system' from '/home/bschroed/Documents/projects/PyGromosTools/pygromos/simulations/hpc_queuing/submission_systems/_submission_system.py'_ )
## Module contents
