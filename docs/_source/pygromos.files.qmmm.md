---
date: '2022-04-01T13:55:18.499Z'
docname: _source/pygromos.files.qmmm
images: {}
path: /source-pygromos-files-qmmm
title: pygromos.files.qmmm package
---

# pygromos.files.qmmm package

## Submodules

## pygromos.files.qmmm.qmmm module


### _class_ pygromos.files.qmmm.qmmm.QMMM(in_value: str, _future_file: bool = False)
Bases: `pygromos.files._basics._general_gromos_file._general_gromos_file`


#### DFTBELEMENTS(_: [pygromos.files.blocks.qmmm_blocks.DFTBELEMENTS](#pygromos.files.blocks.qmmm_blocks.DFTBELEMENTS_ )

#### MNDOELEMENTS(_: [pygromos.files.blocks.qmmm_blocks.MNDOELEMENTS](#pygromos.files.blocks.qmmm_blocks.MNDOELEMENTS_ )

#### MOPACELEMENTS(_: [pygromos.files.blocks.qmmm_blocks.MOPACELEMENTS](#pygromos.files.blocks.qmmm_blocks.MOPACELEMENTS_ )

#### ORCAELEMENTS(_: [pygromos.files.blocks.qmmm_blocks.ORCAELEMENTS](#pygromos.files.blocks.qmmm_blocks.ORCAELEMENTS_ )

#### QMUNIT(_: [pygromos.files.blocks.qmmm_blocks.QMUNIT](#pygromos.files.blocks.qmmm_blocks.QMUNIT_ )

#### QMZONE(_: [pygromos.files.blocks.qmmm_blocks.QMZONE](#pygromos.files.blocks.qmmm_blocks.QMZONE_ )

#### TITLE(_: [pygromos.files.blocks._general_blocks.TITLE](#pygromos.files.blocks.pertubation_blocks.TITLE_ )

#### TURBOMOLEELEMENTS(_: [pygromos.files.blocks.qmmm_blocks.TURBOMOLEELEMENTS](#pygromos.files.blocks.qmmm_blocks.TURBOMOLEELEMENTS_ )

#### XTBELEMENTS(_: [pygromos.files.blocks.qmmm_blocks.XTBELEMENTS](#pygromos.files.blocks.qmmm_blocks.XTBELEMENTS_ )

#### _health_check()
Runs tests on the integrity of the file and spits out warnings
members = [attr for attr in dir(self) if not callable(getattr(self, attr)) and attr.endswith(“ELEMENTS”)]


#### get_qm_engines()
Returns the QM engine used


* **Returns**

    A list of strings (in case the user wanted for some reason specify more than one engine)



* **Return type**

    List[str]



#### path(_: st_ )

#### read_file()
> give back the content. WARNING DEAPRECEATED.

**WARNING**: DEAPRECEATED


* **Returns**

    key is the block name of the gromos file, any is the content of a block



* **Return type**

    Dict[str, any]


## Module contents
