---
date: '2022-04-01T13:55:18.499Z'
docname: _source/pygromos.files.simulation_parameters
images: {}
path: /source-pygromos-files-simulation-parameters
title: pygromos.files.simulation_parameters package
---

# pygromos.files.simulation_parameters package

## Submodules

## pygromos.files.simulation_parameters.imd module

FUNCTIONLIB:            gromos++ input file functions
Description:

> in this lib, gromosXX input file mainpulating functions are gathered

Author: Kay Schaller & Benjamin Schroeder


### _class_ pygromos.files.simulation_parameters.imd.Imd(in_value: str, _future_file: bool = False)
Bases: `pygromos.files._basics._general_gromos_file._general_gromos_file`


#### BOUNDCOND(_: [pygromos.files.blocks.imd_blocks.BOUNDCOND](#pygromos.files.blocks.imd_blocks.BOUNDCOND_ )

#### COMTRANSROT(_: [pygromos.files.blocks.imd_blocks.COMTRANSROT](#pygromos.files.blocks.imd_blocks.COMTRANSROT_ )

#### CONSTRAINT(_: [pygromos.files.blocks.imd_blocks.CONSTRAINT](#pygromos.files.blocks.imd_blocks.CONSTRAINT_ )

#### DISTANCERES(_: [pygromos.files.blocks.imd_blocks.DISTANCERES](#pygromos.files.blocks.imd_blocks.DISTANCERES_ )

#### EDS(_: [pygromos.files.blocks.imd_blocks.EDS](#pygromos.files.blocks.imd_blocks.EDS_ )

#### ENERGYMIN(_: [pygromos.files.blocks.imd_blocks.ENERGYMIN](#pygromos.files.blocks.imd_blocks.ENERGYMIN_ )

#### FORCE(_: [pygromos.files.blocks.imd_blocks.FORCE](#pygromos.files.blocks.imd_blocks.FORCE_ )

#### INITIALISE(_: [pygromos.files.blocks.imd_blocks.INITIALISE](#pygromos.files.blocks.imd_blocks.INITIALISE_ )

#### MULTIBATH(_: [pygromos.files.blocks.imd_blocks.MULTIBATH](#pygromos.files.blocks.imd_blocks.MULTIBATH_ )

#### NEW_REPLICA_EDS(_: [pygromos.files.blocks.imd_blocks.NEW_REPLICA_EDS](#pygromos.files.blocks.imd_blocks.NEW_REPLICA_EDS_ )

#### NONBONDED(_: [pygromos.files.blocks.imd_blocks.NONBONDED](#pygromos.files.blocks.imd_blocks.NONBONDED_ )

#### PAIRLIST(_: [pygromos.files.blocks.imd_blocks.PAIRLIST](#pygromos.files.blocks.imd_blocks.PAIRLIST_ )

#### POSITIONRES(_: [pygromos.files.blocks.imd_blocks.POSITIONRES](#pygromos.files.blocks.imd_blocks.POSITIONRES_ )

#### PRESSURESCALE(_: [pygromos.files.blocks.imd_blocks.PRESSURESCALE](#pygromos.files.blocks.imd_blocks.PRESSURESCALE_ )

#### PRINTOUT(_: [pygromos.files.blocks.imd_blocks.PRINTOUT](#pygromos.files.blocks.imd_blocks.PRINTOUT_ )

#### QMMM(_: [pygromos.files.blocks.imd_blocks.QMMM](#pygromos.files.blocks.imd_blocks.QMMM_ )

#### REPLICA(_: [pygromos.files.blocks.imd_blocks.REPLICA](#pygromos.files.blocks.imd_blocks.REPLICA_ )

#### REPLICA_EDS(_: [pygromos.files.blocks.imd_blocks.REPLICA_EDS](#pygromos.files.blocks.imd_blocks.REPLICA_EDS_ )

#### STEP(_: [pygromos.files.blocks.imd_blocks.STEP](#pygromos.files.blocks.imd_blocks.STEP_ )

#### SYSTEM(_: [pygromos.files.blocks.imd_blocks.SYSTEM](#pygromos.files.blocks.imd_blocks.SYSTEM_ )

#### TITLE(_: [pygromos.files.blocks._general_blocks.TITLE](#pygromos.files.blocks.pertubation_blocks.TITLE_ )

#### WRITETRAJ(_: [pygromos.files.blocks.imd_blocks.WRITETRAJ](#pygromos.files.blocks.imd_blocks.WRITETRAJ_ )

#### edit_EDS(NUMSTATES: int, S: float, EIR: list, EDS: int = 1, ALPHLJ: float = 0.0, ALPHC: float = 0.0, FUNCTIONAL: int = 1)

#### edit_REEDS(REEDS: Optional[int] = None, NUMSTATES: Optional[int] = None, SVALS: Optional[Union[numbers.Number, List[numbers.Number]]] = None, EIR: Optional[Union[numbers.Number, Iterable[numbers.Number]]] = None, NRETRIAL: Optional[int] = None, NREQUIL: Optional[int] = None, CONT: Optional[Union[bool, int]] = None, EDS_STAT_OUT: Optional[Union[bool, int]] = None)

#### path(_: st_ )

#### randomize_seed()

#### read_file()
> give back the content. WARNING DEAPRECEATED.

**WARNING**: DEAPRECEATED


* **Returns**

    key is the block name of the gromos file, any is the content of a block



* **Return type**

    Dict[str, any]



#### write_json(out_path: str)
## Module contents
