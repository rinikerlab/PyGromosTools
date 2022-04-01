---
date: '2022-04-01T13:55:18.499Z'
docname: _source/pygromos.files.otherfiles
images: {}
path: /source-pygromos-files-otherfiles
title: pygromos.files.otherfiles package
---

# pygromos.files.otherfiles package

## Submodules

## pygromos.files.otherfiles.new_repdat module

FUNCTIONLIB:            Repdat File
Description:

> From a Replica-Exchange simulation a repdat file will be created, that gives insight on the replica exchanges of the simulation

Author: Benjamin Schroeder


### _class_ pygromos.files.otherfiles.new_repdat.Repdat(input_path: str)
Bases: `pandas.core.frame.DataFrame`

Replica exchange statistic file
This class is a representation for all transition information during a replica exchange run. it adds some useful functionality.


#### DATA(_: pandas.core.frame.DataFram_ )

#### SYSTEM(_: [pygromos.files.blocks.replica_exchange_blocks.repex_system](#pygromos.files.blocks.replica_exchange_blocks.repex_system_ )

#### \__init__(input_path: str)
Repdat_Constructor


* **Parameters**

    **input_path** (*str*) – path to gromos repdadat file



#### _caculate_transition_traces()
TODO: refactor code!
..autofunction: _caculate_transition_traces

> calculates the transition traces for all replicas from raw data and stores them in self.transition_traces.
> In the end you recieve the trace a replica coord system moved through s dist
> format: {replicaID: {[trials…], [position…], [PotE…]}}


* **Returns**

    None



* **Return type**

    None



#### _calculate_ndowns_nups_for_each_state(time_stride: int = - 1, min_state_potential_treshold: Optional[float] = None, verbose: bool = False)
> > calculates the visit counts for each replicaID position (Temperature or s_value).

> It splits into substates depending on the state potentials, to destinguish which state passed by.
> Up and Downs are counted from top to bottom.
> Additionally also a time dependend series is generated. This can be binned by the  argument time_window.
> If the min_state_potential_treshold is given, than a minimal state is also dependent on the other states, if they are below the threshold, the state is undefined.

> In the end you recieve the position state visit counts in a dict:
> format: {{replicaposition:{“tot_nup”:[], “tot_ndown”:[], “dt”:float, “dt_nup”:[], “dt_ndown”:[]}}


* **Parameters**

    
    * **time_stride** (*int, optional*) – determines the window bin size of flow trajectory for each replicaID. This there are total_transitions/time_window bins containing time_window many flow values. default -1 counts all frames


    * **min_state_potential_treshold** (*float, optional*) – a threshold, defining if a state is governing a system at a time point t


    * **verbose** (*bool, optional*)



* **Return type**

    None



#### _calculate_replica_roundtrips()
..autofunction: _calculate_replica_roundtrips

    This function is calculating the roundtrips over all replica positions for each replica.


* **Returns**

    None



* **Return type**

    None



#### _clean_replica_round_trips(replica_round_trips: Dict[int, int])
_clean_replica_round_trips - privat

> This function cleans up so that the minimal rountrip number in a roundtrip dict is 0


* **Parameters**

    **replica_round_trips** (*Dict[int:int]*) – a dictionary containing all replica roundtrip counts



* **Returns**

    **Dict[int** – a dictionary containing all replica roundtrip counts with lowest value 0



* **Return type**

    int]



#### added_property(_ = _ )

#### append(repdat: Union[List[pygromos.utils.typing.Repdat_Type], pygromos.utils.typing.Repdat_Type])
> This function concatenates two repdat files into the executing obj.


* **Parameters**

    **repdat** (*List[Repdat] or Repdat*) – one or multiple Repdat files.



* **Return type**

    None



#### clean_file_runs()
clean_file
DEAPRECEATED - ABOUT TO BE REMOVED!

> Updates the run numbers to be continous sequential. (for example needed for concatenation)


#### count_state_per_position(_: Dict[int, Dict[str, Union[List[int], int]]_ _ = Non_ )

#### get_replicaPosition_dependend_nup_ndown(time_window_size: int = - 1, potential_treshold: Optional[float] = None, recalculate: bool = False)
..autofunction: get_replicaPosition_dependend_nup_ndown

    This function is returning the replica position visit counts by all simulation state.


* **Parameters**

    
    * **time_window_size** (*int*) – how many timesteps shall be binned


    * **potential_treshold** (*float*) – if defined, and there is a time window, in which multiple states are below this threshold, the count is in an undefined state.


    * **recalculate** (*bool*) – shall the dict be recalculated?



* **Returns**

    returns a dict for all replica positions and the visit counts.



* **Return type**

    Dict[int, Dict[str, Union[List or float]]]



#### get_replicaPosition_dependend_nup_ndown_for_each_state(time_window_size: int = - 1, potential_treshold: Optional[float] = None, recalculate: bool = False)
..autofunction: get_replicaPosition_dependend_nup_ndown_for_each_state

    This function is returning the replica position visit counts by each simulation state, per state.


* **Parameters**

    
    * **time_window_size** (*int*) – how many timesteps shall be binned


    * **potential_treshold** (*float*) – if defined, and there is a time window, in which multiple states are below this threshold, the count is in an undefined state.


    * **recalculate** (*bool*) – shall the dict be recalculated?



* **Returns**

    returns a dict for all replica positions and their state visit counts.



* **Return type**

    Dict[int, Dict[str, Union[List or float]]]



#### get_replica_roundtrips(recalculate: bool = False)
..autofunction: get_replica_roundtrips

    This function is returning the count of rountrips (RT) for each replica.


* **Parameters**

    **recalculate** (*bool*) – shall the dict be recalculated?



* **Returns**

    returns a dict for all replica and their rountrip counts.



* **Return type**

    Dict[int,int]



#### get_replica_traces(recalculate: bool = False)
returns a replica_traces dictionary.
:Parameters: **recalculate** (*bool, optional*) – shall the dict be recalculated, if already present?


* **Returns**

    dictionary containing all individual replica_traces



* **Return type**

    Dict[int, Dict[str,List[float]]]



#### replica_round_trips(_: Dict[int, int_ _ = Non_ )

#### transition_traces(_: Dict[int, Dict[str, List[float]]_ _ = Non_ )

#### write(out_path: str)
..autofunction:

    Write out a repdat file to the outpath.


* **Parameters**

    **out_path** (*str*) – determines the output path



* **Returns**

    out_path


:rtype:str

## pygromos.files.otherfiles.noe_output module

FUNCTIONLIB:            gromos++ input file functions
Description:

> in this lib, gromosXX input file mainpulating functions are gathered

Author: Benjamin Schroeder


### _class_ pygromos.files.otherfiles.noe_output.JVAL(in_value: str)
Bases: `pygromos.files._basics._general_gromos_file._general_gromos_file`


#### content(_: pandas.core.frame.DataFram_ )

#### read_file()
> give back the content. WARNING DEAPRECEATED.

**WARNING**: DEAPRECEATED


* **Returns**

    key is the block name of the gromos file, any is the content of a block



* **Return type**

    Dict[str, any]



#### write(out_path: str)
writes a general Gromos File out


* **Parameters**

    **out_path** (*str*) – out path where the file should be.



* **Returns**

    out_path



* **Return type**

    str



### _class_ pygromos.files.otherfiles.noe_output.NOE(in_value: str)
Bases: `pygromos.files._basics._general_gromos_file._general_gromos_file`


#### AVERAGE_NOE(_: pandas.core.frame.DataFram_ )

#### NOE_VIOLATIONS(_: pandas.core.frame.DataFram_ )

#### RESTRAINT_LEGEND(_: pandas.core.frame.DataFram_ )

#### TITLE(_: [pygromos.files.blocks._general_blocks.TITLE](#pygromos.files.blocks.pertubation_blocks.TITLE_ )

#### read_file()
> give back the content. WARNING DEAPRECEATED.

**WARNING**: DEAPRECEATED


* **Returns**

    key is the block name of the gromos file, any is the content of a block



* **Return type**

    Dict[str, any]


## pygromos.files.otherfiles.repdat module

FUNCTIONLIB:            Repdat File
Description:

> From a Replica-Exchange simulation a repdat file will be created, that gives insight on the replica exchanges of the simulation

Author: Benjamin Schroeder


### _class_ pygromos.files.otherfiles.repdat.Repdat(input_path: str)
Bases: `pygromos.files._basics._general_gromos_file._general_gromos_file`

Replica exchange statistic file
This class is a representation for all transition information during a replica exchange run. it adds some useful functionality.


#### DATA(_: pandas.core.frame.DataFram_ )

#### SYSTEM(_: [pygromos.files.blocks.replica_exchange_blocks.repex_system](#pygromos.files.blocks.replica_exchange_blocks.repex_system_ )

#### \__init__(input_path: str)
Repdat_Constructor


* **Parameters**

    **input_path** (*str*) – path to gromos repdadat file



#### _caculate_transition_traces()
calculates the transition traces for all replicas from raw data and stores them in self.transition_traces.

    In the end you recieve the trace a replica coord system moved through s dist
    format: {replicaID: {[trials…], [position…], [PotE…]}}

TODO: refactor code!


* **Return type**

    None



#### _calculate_ndowns_nups_for_each_state(time_stride: int = - 1, min_state_potential_treshold: Optional[float] = None, verbose: bool = False)
> > calculates the visit counts for each replicaID position (Temperature or s_value).

> It splits into substates depending on the state potentials, to destinguish which state passed by.
> Up and Downs are counted from top to bottom.
> Additionally also a time dependend series is generated. This can be binned by the  argument time_window.
> If the min_state_potential_treshold is given, than a minimal state is also dependent on the other states, if they are below the threshold, the state is undefined.

> In the end you recieve the position state visit counts in a dict:
> format: {{replicaposition:{“tot_nup”:[], “tot_ndown”:[], “dt”:float, “dt_nup”:[], “dt_ndown”:[]}}


* **Parameters**

    
    * **time_stride** (*int, optional*) – determines the window bin size of flow trajectory for each replicaID. This there are total_transitions/time_window bins containing time_window many flow values. default -1 counts all frames


    * **min_state_potential_treshold** (*float, optional*) – a threshold, defining if a state is governing a system at a time point t


    * **verbose** (*bool, optional*)



* **Return type**

    None



#### _calculate_replica_roundtrips()
This function is calculating the roundtrips over all replica positions for each replica.


* **Return type**

    None



#### _clean_replica_round_trips(replica_round_trips: Dict[int, int])
_clean_replica_round_trips - privat

> This function cleans up so that the minimal rountrip number in a roundtrip dict is 0


* **Parameters**

    **replica_round_trips** (*Dict[int:int]*) – a dictionary containing all replica roundtrip counts



* **Returns**

    **Dict[int** – a dictionary containing all replica roundtrip counts with lowest value 0



* **Return type**

    int]



#### append(repdat: Union[List[pygromos.utils.typing.Repdat_Type], pygromos.utils.typing.Repdat_Type])
> This function concatenates two repdat files into the executing obj.


* **Parameters**

    **repdat** (*List[Repdat] or Repdat*) – one or multiple Repdat files.



* **Return type**

    None



#### clean_file_runs(starting_trial: float = 1)
clean_file

    Updates the run numbers to be continous sequential. (for example needed for concatenation)


* **Parameters**

    **starting_trial** (*int, optional*)



#### count_state_per_position(_: Dict[int, Dict[str, Union[List, float]]_ _ = Non_ )

#### get_replicaPosition_dependend_nup_ndown(time_window_size: int = - 1, potential_treshold: Optional[float] = None, recalculate: bool = False)
This function is returning the replica position visit counts by all simulation state.


* **Parameters**

    
    * **time_window_size** (*int,optional*) – how many timesteps shall be binned into one bin?


    * **potential_treshold** (*float*) – if defined, and there is a time window, in which multiple states are below this threshold, the count is in an undefined state.


    * **recalculate** (*bool*) – shall the dict be recalculated?



* **Returns**

    returns a dict for all replica positions and the visit counts.



* **Return type**

    Dict[int, Dict[str, Union[List or float]]]



#### get_replicaPosition_dependend_nup_ndown_for_each_state(time_window_size: int = - 1, potential_treshold: Optional[float] = None, recalculate: bool = False)
This function is returning the replica position visit counts by each simulation state, per state.


* **Parameters**

    
    * **time_window_size** (*int, optional*) – how many timesteps shall be binned into one bin?


    * **potential_treshold** (*float, optional*) – if defined, and there is a time window, in which multiple states are below this threshold, the count is in an undefined state.


    * **recalculate** (*bool, optional*) – shall the dict be recalculated?



* **Returns**

    returns a dict for all replica positions and their state visit counts.



* **Return type**

    Dict[int, Dict[str, Union[List or float]]]



#### get_replica_roundtrips(recalculate: bool = False)
This function is returning the count of rountrips (RT) for each replica.


* **Parameters**

    **recalculate** (*bool*) – shall the dict be recalculated?



* **Returns**

    returns a dict for all replica and their rountrip counts.



* **Return type**

    Dict[int, int]



#### get_replica_traces(recalculate: bool = False)
returns a replica_traces dictionary.
:Parameters: **recalculate** (*bool, optional*) – shall the dict be recalculated, if already present?


* **Returns**

    dictionary containing all individual replica_traces



* **Return type**

    Dict[int, Dict[str,List[float]]]



#### replica_round_trips(_: Dict[int, int_ _ = Non_ )

#### transition_traces(_: Dict[int, Dict[str, List[float]]_ _ = Non_ )

#### write(out_path: str)

* **Parameters**

    **out_path** (*str*) – determines the output path for repdat file



* **Returns**

    out_path



* **Return type**

    str


## pygromos.files.otherfiles.residue_library module

File:          gromos residue library

    needed if top and pdb resns or atoms are not the same.

Author: Benjamin Ries


### _class_ pygromos.files.otherfiles.residue_library.residue_library(in_value: Union[str, Dict] = '/home/bschroed/Documents/projects/PyGromosTools/pygromos/data/pdb2g96.lib')
Bases: `pygromos.files._basics._general_gromos_file._general_gromos_file`


#### ATOMNAMELIB(_: [pygromos.files.blocks.miscBlocks.ATOMNAMELIB](#pygromos.files.blocks.miscBlocks.ATOMNAMELIB_ )

#### RESIDUENAMELIB(_: [pygromos.files.blocks.miscBlocks.RESIDUENAMELIB](#pygromos.files.blocks.miscBlocks.RESIDUENAMELIB_ )

#### \__init__(in_value: Union[str, Dict] = '/home/bschroed/Documents/projects/PyGromosTools/pygromos/data/pdb2g96.lib')
> This class represents a file that is used for the gromosPP program - pdb2g96
> it contains two blocks for residue naming and atom naming


* **Parameters**

    **in_value** (*Union[str, dict]*)



#### read_resnlib(path: str)

#### required_blocks(_ = ['TITLE', 'RESIDUENAMELIB', 'ATOMNAMELIB'_ )

#### verbose(_: boo_ _ = Fals_ )
## Module contents
