---
date: '2022-04-01T13:55:18.499Z'
docname: _source/pygromos.gromos
images: {}
path: /source-pygromos-gromos
title: pygromos.gromos package
---

# pygromos.gromos package

## Subpackages


* [pygromos.gromos.pyGromosPP package]()


    * [Submodules](#submodules)


    * [pygromos.gromos.pyGromosPP.com_top module](#module-pygromos.gromos.pyGromosPP.com_top)


    * [pygromos.gromos.pyGromosPP.ran_box module](#module-pygromos.gromos.pyGromosPP.ran_box)


    * [Module contents](#module-pygromos.gromos.pyGromosPP)


## Submodules

## pygromos.gromos.compile_gromos module


### pygromos.gromos.compile_gromos._configure_gromosPP_autotools(build_dir: str, binary_dir: Optional[str] = None, with_omp: bool = False, with_debug: bool = False, verbose: bool = True, _timing_dict: dict = {})
> Setting the configurations for the compiling gromosPP process. (uses autotools)


* **Parameters**

    
    * **build_dir** (*str*) – directory, that should be used for building


    * **binary_dir** (*str, optional*) – directory in which the binaries should be written to, by default None


    * **with_omp** (*bool, optional*) – should gromosPP be compiled with omp, by default False


    * **with_debug** (*bool, optional*) – set gromos debug flag, by default False


    * **verbose** (*bool, optional*) – compiling is fun, I can tell you more!, by default True


    * **_timing_dict** (*dict, optional*) – structure for storing timings of process, by default {}



### pygromos.gromos.compile_gromos._configure_gromosXX_autotools(build_dir: str, binary_dir: Optional[str] = None, with_cuda_dir: Optional[str] = None, with_omp: bool = False, with_mpi: bool = False, with_debug: bool = False, verbose: bool = True, _timing_dict: dict = {})
> Setting the configurations for the compiling gromosXX process. (uses autotools)


* **Parameters**

    
    * **build_dir** (*str*) – directory, that should be used for building


    * **binary_dir** (*str, optional*) – directory in which the binaries should be written to, by default None


    * **with_cuda_dir** (*str, optional*) – use the following cuda path and activate cuda support, by default None


    * **with_omp** (*bool, optional*) – should gromosXX be compiled with omp - Warning dont combine wiht mpi!, by default False


    * **with_mpi** (*bool, optional*) – should gromosXX be compiled wiht mpi - Warning dont combine with omp!, by default False


    * **with_debug** (*bool, optional*) – set gromos debug flag, by default False


    * **verbose** (*bool, optional*) – compiling is fun, I can tell you more!, by default True


    * **_timing_dict** (*dict, optional*) – structure for storing timings of process, by default {}



* **Raises**

    **ValueError** – if wrong value was passed



### pygromos.gromos.compile_gromos._make_clean(build_dir: str, nCore: int = 1, verbose: bool = True, _timing_dict: dict = {}, _timing_prefix: str = '')
> This function triggers make clean and removes the build_dir.


* **Parameters**

    
    * **build_dir** (*str*) – directory prepared for make.


    * **nCore** (*int, optional*) – how many cores for each make command?, by default 1


    * **verbose** (*bool, optional*) – cleaning, cleaning, every thing gets so shiny and I can sing!, by default True


    * **_timing_dict** (*dict, optional*) – structure for storing timings of process, by default {}


    * **_timing_prefix** (*str, optional*) – prefix for timing keys, by default “”



### pygromos.gromos.compile_gromos._make_compile(build_dir: str, nCore: int = 1, verbose: bool = True, _timing_dict: dict = {}, _timing_prefix: str = '')
> This function triggers make and make install in the build_dir.


* **Parameters**

    
    * **build_dir** (*str*) – directory prepared for make.


    * **nCore** (*int, optional*) – how many cores for each make command?, by default 1


    * **verbose** (*bool, optional*) – make, make things and talk loudly about it! , by default True


    * **_timing_dict** (*dict, optional*) – structure for storing timings of process, by default {}


    * **_timing_prefix** (*str, optional*) – prefix for timing keys, by default “”



### pygromos.gromos.compile_gromos.default_install(nCores: int = 1, _timing_dict: dict = {}, verbose: bool = False)

### pygromos.gromos.compile_gromos.install_gromos(root_dir: Optional[str] = None, nCore: int = 3, gromosXX_with_mpi: bool = False, gromosXX_with_omp: bool = False, gromosXX_with_cuda: Optional[str] = None, gromosPP_with_omp=False, gromosPP_with_debug: bool = False, gromosXX_with_debug: bool = False, do_compile: bool = True, do_clean: bool = True, recompile: bool = False, recompile_from_scratch: bool = False, _do_gromosPP: bool = True, _do_gromosXX: bool = True, _timing_dict: dict = {}, verbose: bool = True)
> Install the gromos simulation packages. As a helper: to get the dependencies right, install and activate the provided pygromos environments!


* **Parameters**

    
    * **root_dir** (*str*) – this dir contains the gromosXX and the gromosPP git repositories


    * **nCore** (*int, optional*) – how many cores should be used to compile?, by default 1


    * **gromosXX_with_cuda_dir** (*str, optional*) – use the following cuda path and activate cuda support, by default None


    * **gromosXX_with_omp** (*bool, optional*) – should gromosXX be compiled with omp - Warning dont combine wiht mpi!, by default False


    * **gromosXX_with_mpi** (*bool, optional*) – should gromosXX be compiled wiht mpi - Warning dont combine with omp!, by default False


    * **gromosPP_with_omp** (*bool, optional*) – should gromosPP be compiled with omp, by default False


    * **gromosPP_with_debug** (*bool, optional*) – set gromosPP debug flag, by default False


    * **gromosXX_with_debug** (*bool, optional*) – set gromosXX debug flag, by default False


    * **do_compile** (*bool, optional*) – compile the programms, by default True


    * **do_clean** (*bool, optional*) – clean up and remove the programs (can be used together with do_compile to get a clean start), by default False


    * **recompile** (*bool, optional*) – recompile the programs with make, by default False


    * **recompile_from_scratch** (*bool, optional*) – recompile with configure and make, by default False


    * **_do_gromosPP** (*bool, optional*) – install the gromosPP program, by default True


    * **_do_gromosXX** (*bool, optional*) – install the gromosXX program, by default True


    * **_timing_dict** (*dict, optional*) – structure for storing timings of process, by default {}


    * **verbose** (*bool, optional*) – compiling is fun, I can tell you more!, by default True


## pygromos.gromos.gromosBashSyntaxParser module


### _class_ pygromos.gromos.gromosBashSyntaxParser.gromosBashSyntaxParser()
Bases: `object`

Helper class to parse general gromos bash syntax

all methods should be static


#### _static_ atomSliceParser()

#### _static_ moleculeSliceParser()

#### _static_ multiplyArgumentParser(args: Union[str, List[str]], multiplier: Union[int, List[int]] = 1)
Parser for multiplier syntax to gromos scripts

example: com_top @topo 1:Protein 2:Na 2:Cl 100:SPC …..


* **Parameters**

    
    * **args** (*str or list(str)*) – The actual argument which should be multiplied (ex. top)


    * **multiplier** (*int or list(int)*) – the multiplier for each argument provided in args


## pygromos.gromos.gromosPP module

FUNCTIONLIB:            wrapper for gromos++
Description:

> This file contains python wrappers for the bash commandline of gromos++

Author: Benjamin Schroeder


### _class_ pygromos.gromos.gromosPP.GromosPP(gromosPP_bin_dir: Optional[str] = None, _check_binary_paths: bool = True, verbose: bool = False)
Bases: `pygromos.gromos.gromosPP._gromosPPbase`

This is the class represents gromosPP.

bin

    This is the path to the folder containing the binaries of gromosPP If None, the bash enviroment variables  will be used.


### _class_ pygromos.gromos.gromosPP._gromosPPbase(gromosPP_bin_dir: Optional[str] = None, _check_binary_paths: bool = True, verbose: bool = False)
Bases: `pygromos.gromos._gromosClass._gromosClass`

GromosPP

This is the gromosPP baseclass.
This should be inherited by a concrete class that might reimplement some new features, that are version dependent.

bin

    This is the path to the folder containing the binaries of gromosPP. If None, the bash enviroment variables  will be used.


#### \__init__(gromosPP_bin_dir: Optional[str] = None, _check_binary_paths: bool = True, verbose: bool = False)
Constructing a gromosPP object.


* **Parameters**

    
    * **bin** (*Union[str, None], optional*) – This is the path to the folder containing the binaries of gromosXX. If None, the bash enviroment variables  will be used.


    * **_dont_check_binary** (*bool, optional*) – This flag removes the checks of the binary presence for this obj. This can make sense if system access is slow!, by default False - checks will be made



#### add_hydrogens(in_cnf_path: str, in_top_path: str, out_cnf_path: str, tolerance: float = 0.1, periodic_boundary_condition: str = 'v', gathering: str = 'cog', _binary_name: str = 'gch')
> This function protonates a coordinate file.


* **Parameters**

    
    * **in_cnf_path** (*str*)


    * **in_top_path** (*str*)


    * **out_cnf_path** (*str*)


    * **tolerance** (*float, optional*)


    * **periodic_boundary_condition** (*str, optional*)


    * **gathering** (*str, optional*)


    * **_binary_name** (*str, optional*)



* **Return type**

    out_cnf_path



#### amber2gromos(ambertop: str, solvent: str, ljscaling: float = 2.0, chargegroups: str = None, atomic_chargegroups: bool = False, out_path: str = 'topology.top', _binary_name: str = 'amber2gromos', verbose: bool = False)

* **Parameters**

    
    * **ambertop** (*str*) – path to AMBER molecular topology file


    * **solvent** (*str*) – path to GROMOS topology file with solvent


    * **ljscaling** (*float, optional*) – scaling factor for LJ parameters (default: 2.0)


    * **atomic_chargegroups** (*bool, optional*) – assign each atom to its own chargegroup (default: False)


    * **chargegroups** (*str, optional*) – path to chargegroup file


    * **_binary_name** (*str, otpional*) – binary name of gromosPP programm (default: amber2gromos)



* **Returns**

    converted GROMOS topology



* **Return type**

    str



#### build_box(in_top_path: str, in_cnf_path: str, out_cnf_path: str = '', periodic_boundary_condition: str = 'r', nmolecule: int = 1, dens: float = 1.0, _binary_name: str = 'build_box', verbose: bool = False, return_command_only: bool = False)

#### cog(in_top_path: str, in_trcs: Union[str, List[str]], out_file_path: str, atom_selection: str = None, outformat: str = None, cog_com: str = None, add_repl: str = None, solv: str = None, nthframe: str = None, pbc: str = 'r cog', _binary_name: str = 'cog', verbose: bool = False)

* **Parameters**

    
    * **in_top_path** (*str*) – path to topo file


    * **in_trcs** (*str*) – path to trc files


    * **out_file_path** (*str*) – outpath


    * **atom_selection** (*str, optional*) – atomspecifier(s) for which to calculate cog/com


    * **outformat** (*str, optional*) – output coordinates format


    * **cog_com** (*str, optional*) – calculate centre of geometry (cog) or mass (com); default: cog


    * **add_repl** (*str, optional*) – add (add) the cog/com or replace (repl) the solutes; default: repl


    * **solv** (*str, optional*) – include solvent in outcoordinates


    * **nthframe** (*str, optional*) – write every nth frame (default: 1)


    * **_binary_name** (*str, otpional*) – binary name of gromosPP programm (default: cog)



* **Returns**

    output_path of the generated csv file



* **Return type**

    str



#### com_top(in_topo_paths: Union[str, List[str]], topo_multiplier: Union[int, List[int]] = 1, out_top_path: str = 'combined_out.top', take_topology_params_of_file: int = 1, take_solvent_parameters_of_file: int = 1, _binary_name: str = 'com_top')
Combine multiple topologies
:Parameters: \* **in_topo_paths** (*(str or List[str])*)

> 
> * **out_top_path** (*str, optional*)


> * **take_topology_params_of_file** (*int, optional*)


> * **take_solvent_parameters_of_file** (*int, optional*)


> * **_binary_name** (*str, optional*)


* **Returns**

    output file_path



* **Return type**

    str



#### dfmult(in_endstate_file_paths: List[str], in_reference_state_file_path: str, out_file_path: str = 'dfmult_temp.out', temperature: float = 298, _binary_name: str = 'dfmult', verbose: bool = False)
> This funciton wraps dfmult of the gromos suite, it is used to calculate the free energy of a EDS simulation.


* **Parameters**

    
    * **in_endstate_file_paths** (*List[str]*) – potentials energy of single state (get with ene Ana)


    * **in_reference_state_file_path** (*str*) – potential energy of reference State (get with ene Ana)


    * **out_file_path** (*str, optional*)


    * **temperature** (*float, optional*)


    * **_binary_name** (*str, optional*)


    * **verbose** (*bool, optional*)



* **Returns**

    out_put file path



* **Return type**

    str



#### ene_ana(in_ene_ana_library_path: str, in_en_file_paths: str, in_properties: str, out_energy_folder_path: str, out_files_prefix: str = None, out_files_suffix: str = None, in_topo: str = None, time: float = None, single_file: bool = False, return_outstat_also: bool = False, verbose: bool = False, _binary_name: str = 'ene_ana', workdir: bool = False)
> This is a wrapper for ene_ana.

>     ene_ana is a tool to extract energy properties from a x.tre file.


* **Parameters**

    
    * **in_ene_ana_library_path** (*str*)


    * **in_en_file_paths** (*str*)


    * **in_properties** (*str*) –

    this str should name the properties, that shall be written out. For the variable names, check the ene_ana lib.

        (e.g.: “solvtemp1 e1 eR”)


    * **out_energy_folder_path** (*str*) – give the path to an directory, where the output should be stored in.


    * **out_files_prefix** (*None or str, optional*)


    * **out_files_suffix** (*None or str, optional*)


    * **in_topo** (*None or str, optional*)


    * **time** (*None or float, optional*)


    * **single_file** (*bool, optional*) – if used a single csv is generated. the return value gets the file_path(str).


    * **verbose** (*bool, optional*)


    * **_binary_name** (*str, optional*)



* **Returns**

    if single_file used, return = str



* **Return type**

    List[str] or str



#### filter(out_filter_path: str, in_coord_path: str, in_top_path: str, atom_selection: Optional[str] = None, periodic_boundary_condition: str = 'r cog', cutoff: Optional[float] = None, pairlist: Optional[str] = None, select: str = '1:a', reject: Optional[str] = None, time: Optional[int] = None, dt: Optional[int] = None, outformat: Optional[str] = None, _binary_name: str = 'filter')
This wrapper uses filter to reduce a given trajectory to a selection of atoms. By default, only the first residue (1:a) is retained.


* **Parameters**

    
    * **out_filter_path** (*str*)


    * **in_coord_path** (*str*)


    * **in_top_path** (*str*)


    * **atom_selection** (*str, optional*)


    * **periodic_boundary_condition** (*str, optional*)


    * **cutoff** (*float, optional*)


    * **pairlist** (*str, optional*)


    * **select** (*str, optional*)


    * **reject** (*str, optional*)


    * **time** (*int, optional*)


    * **dt** (*int, optional*)


    * **outformat** (*str, optional*)


    * **_binary_name** (*str, optional*)



* **Returns**

    returns out_filter_path



* **Return type**

    str



#### frameout(in_top_path: str, in_coord_path: str, periodic_boundary_condition: str, out_file_path: str = None, out_file_format: str = None, single_file: bool = None, gather: Union[int, str] = None, include: str = 'SOLUTE', reference_structure_path: str = None, atomsfit: str = None, frames: int = None, time: int = None, dt: int = None, notimeblock: bool = None, _binary_name: str = 'frameout', verbose: bool = False)
> this wrapper wraps frameout.

>     frameout is a tool, that can be used for example to convert gromos coordinate files or to recenter them or ….
>     Optional parameters are not used in the command, when they are None!


* **Parameters**

    
    * **in_top_path** (*str*)


    * **in_coord_path** (*str*)


    * **periodic_boundary_condition** (*str*)


    * **out_file_path** (*None or str, optional*)


    * **out_file_format** (*None or str, optional*)


    * **single_file** (*None or bool, optional*)


    * **gather** (*None or int or str optional*)


    * **include** (*None or str, optional*) – ALL that also includes the solvent.
    SOLUTE this option filters the Solvent out of the output file. (default)
    alternative gromos selection synthax can be used.


    * **reference_structure_path** (*None or str, optional*) – This path should provide a refrence position for atom fitting.


    * **atomsfit** (*None or str, optional*) – This option can be used to fit all frames to one reference structure.
    The selection Syntax follow the Gromos manual (e.g. “1:a” - aligns all frames to the first molecule and all its atoms)
    requires reference_sturcture_path.


    * **frames** (*None or int, optional*)


    * **time** (*None or float, optional*)


    * **dt** (*None or float, optional*)


    * **notimeblock** (*None or bool, optional*)


    * **_binary_name** (*str, optional*)


    * **verbose** (*bool, optional*)



* **Returns**

    out_file_path



* **Return type**

    str



#### gch(in_cnf_path: str, in_top_path: str, out_cnf_path: str, tolerance: float = 0.1, periodic_boundary_condition: str = 'v', gathering: str = 'cog', _binary_name: str = 'gch')
> This function adds reasonable hydrogenpositions a coordinate file.


* **Parameters**

    
    * **in_cnf_path** (*str*)


    * **in_top_path** (*str*)


    * **out_cnf_path** (*str*)


    * **tolerance** (*float, optional*)


    * **periodic_boundary_condition** (*str, optional*)


    * **gathering** (*str, optional*)


    * **_binary_name** (*str, optional*)



* **Return type**

    out_cnf_path



#### ion(in_top_path: str, in_cnf_path: str, out_cnf_path: str, periodic_boundary_condition: str = 'v', negative: list = None, positive: list = None, potential: float = 0.8, mindist: float = 0.8, random_seed: int = None, exclude: str = None, _binary_name: str = 'ion', verbose: bool = False)
When simulating a charged solute in solution, one may wish to include
counter-ions in the molecular system in order to obtain a neutral system, or
a system with a specific ionic strength. The program ion can replace solvent
molecules by atomic ions by placing the
ion at the position of the first atom of a solvent molecule. Substitution of
solvent molecules by positive or negative ions can be performed by selecting
the solvent positions with the lowest or highest Coulomb potential, respectively,
or by random selection. In order to prevent two ions being placed too
close together, a sphere around each inserted ion can be specified from which
no solvent molecules will be substituted by additional ions. In addition, the user can
specify specific water molecules that should not be considered for
replacement.


* **Parameters**

    
    * **in_top_path** (*str*) – the path to the input topology file (.top)


    * **in_cnf_path** (*str*) – the path to the input coordinate file (.cnf), to which the ions shall be added


    * **out_cnf_path** (*str*) – the path to the resulting coordinate (.cnf) file


    * **periodic_boundary_condition** (*str, optional*) – describes the boundary condition of the given system in the cnf. (r - rectangle, v - vacuum, ). a gathering method can be optionally added with a whitespace seperation., by default “v”


    * **negative** (*list, optional*) – the first element of the list is the number of ions and the second element of the list is the type of ion, optionally a third element can be passed giving the residue name, by default None


    * **positive** (*list, optional*) – the first element of the list is the number of ions and the second element of the list is the type of ion, optionally a third element can be passed giving the residue name, by default None


    * **potential** (*float, optional*) – cutoff for potential calculation[nm], by default 0.8


    * **mindist** (*float, optional*) – minimum distance between ions[nm], by default 0.8


    * **random_seed** (*int, optional*) – provide the used random seed, by default None


    * **exclude** (*str, optional*) – if you want to exclude solvent molecules, define a gromos selection here, by default None


    * **_binary_name** (*str, optional*) – the program name, by default “ion”


    * **verbose** (*bool, optional*) – stay a while and listen, by default False



* **Returns**

    returns the resulting cnf-file path



* **Return type**

    str



#### jval(in_top_path: str, in_jval_path: str, in_traj_path: Union[str, List[str]], out_path: str, pbc: str = 'v', gathering: str = 'cog', timeseries: bool = False, rmsd: bool = False, time: float = None, _binary_name: str = 'jval', verbose: bool = False)

* **Parameters**

    
    * **in_top_path** (*str*) – topology path


    * **in_jval_path** (*str*) – jval specification path


    * **in_traj_path** (*str*) – coordinate file


    * **out_path** (*str*) – path to the output file


    * **pbc** (*str, optional*) – default: v
    periodic boundary condition of the coordinates:

    > v - vacuum
    > r - rectangular box


    * **gathering** (*str, optional*) – default: cog
    how the coordinates shall be gathered before the calculation:

    > cog - center of geometry
    > com - center of mass


    * **timeseries** (*bool, optional*)


    * **rmsd** (*bool, optional*)


    * **time** (*(Number, str), optional*)



* **Returns**

    
    * *str* – out_path


    * *NotImplemented*


    * *—————*


    * **time** (*float, float (time, dt)*)




#### make_top(out_top_path: str, in_building_block_lib_path: str, in_parameter_lib_path: str, in_sequence: str, in_solvent: str = 'H2O', _binary_name: str = 'make_top')
This wrapper uses make_top to generate a topology file.


* **Parameters**

    
    * **out_top_path** (*str*)


    * **in_building_block_lib_path** (*str*)


    * **in_parameter_lib_path** (*str*)


    * **in_sequence** (*str*)


    * **in_solvent** (*str, optional*)


    * **additional_options** (*str, optional*)



* **Returns**

    returns out_file_path



* **Return type**

    str



#### noe(in_top_path: str, in_noe_path: str, in_traj_path: str, out_path: str, pbc: str = 'v', gathering: str = 'cog', _binary_name: str = 'noe', verbose: bool = False)

* **Parameters**

    
    * **in_top_path** (*str*) – topology path


    * **in_noe_path** (*str*) – output path of prep_noe


    * **in_traj_path** (*str*) – coordinate file


    * **out_path** (*str*) – path to the output file


    * **pbc** (*str, optional*) – default: v
    periodic boundary condition of the coordinates:

    > v - vacuum
    > r - rectangular box


    * **gathering** (*str, optional*) – default: cog
    how the coordinates shall be gathered before the calculation:

    > cog - center of geometry
    > com - center of mass



* **Returns**

    
    * *str* – out_path


    * *NotImplemented*


    * *—————*


    * **time** (*float, float (time, dt)*)




#### pdb2gromos(in_pdb_path: str, in_top_path: str, out_cnf_path: str = None, in_lib_path: str = '/home/bschroed/Documents/projects/PyGromosTools/pygromos/data/pdb2g96.lib', _binary_name: str = 'pdb2g96', verbose: bool = False)
This is a wrapper for pdb2gromos. It executes the gromosPP binary.


* **Parameters**

    
    * **in_pdb_path** (*str*)


    * **in_top_path** (*str*)


    * **in_lib_path** (*str, optional*) – The lib can be used as a look up table for residues or atoms etc. - see gromos manual


    * **out_cnf_path** (*str, optional*) – The out_cnf is the path for the output file. if not given, the pdb basename and dir is taken as default.


    * **_binary_name** (*str, optional*)


    * **verbose** (*bool, optional*)



* **Returns**

    Returns out_cnf path



* **Return type**

    str



#### pdb2seq(in_pdb_path: str, out_path: str = '', pH: float = 7.4, select: str = 'ALL', gff: str = '54a7', add_head: str = 'NH3+', add_tail: str = 'COO-', _binary_name: str = 'pdb2seq')
This function is translating a pdb into a sequence file, that can be used to generate for example topologies.


* **Parameters**

    
    * **in_pdb_path** (*str*)


    * **out_path** (*str, optional*) – out_path for output File


    * **pH** (*float, optional*)


    * **select** (*str, optional*)


    * **gff** (*str,optional*) – which Gromos ForceField


    * **add_head** (*str, optional*) – protein N-term capping with this group.


    * **add_tail** (*str,optional*) – protein C-term capping with this group


    * **_binary** (*str, optional*)



* **Returns**

    out_path



* **Return type**

    str



#### prep_eds(in_top_paths: List[str], number_of_eds_states: int, param_top_index: int = 1, solv_top_index: int = 1, out_file_path: str = 'dev', _binary_name: str = 'prep_eds', verbose: bool = False)
> prepare eds topology.


* **Parameters**

    
    * **in_top_paths** (*List[str]*)


    * **number_of_eds_states** (*int*)


    * **param_top_index** (*int, optional*)


    * **solv_top_index** (*int, optional*)


    * **out_file_path** (*str, optional*) – output path without file ending. (prefix)


    * **_binary_name** (*str, optional*)


    * **verbose** (*bool, optional*)



* **Returns**

    out_top, out_ptp



* **Return type**

    tuple[str,str]



#### prep_noe(in_top_path: str, in_noe_path: str, in_library_path: str, out_path: str, dish: float = 0.1, disc: float = 0.153, title: str = 'NOE', _binary_name: str = 'prep_noe', in_correction_path: str = None, verbose: bool = False)

* **Parameters**

    
    * **in_top_path** (*str*) – molecular topology file


    * **in_noe_path** (*str*) – NOE specification file


    * **in_library_path** (*str*) – NOE specification library


    * **out_path** (*str*) – path to the output file


    * **dish** (*float, optional*) – carbon-hydrogen distance; default: 0.1 nm


    * **disc** (*float, optional*) – carbon-carbon distance; default: 0.153 nm


    * **title** (*str, optional*) – NOE title for output, default: “NOE”


    * **correction** (*str, optional*) – Correction file -> <correction_file> [correction type]



* **Returns**

    
    * *str* – path to the output file


    * *NotImplemented*


    * *——-*


    * **parsetype** (*<1,2,3>*) –

    Choices are:

        1: Upper bound == first number
        2: Upper bound == first + third number (most common, default)
        3: Upper bound == first - second number (commonly the lower bound)


    * **action** (*<add> or <substraction> = add*)


    * **filter** (*discard nNOE’s above a certain distance[nm] = 10000 nm*)


    * **factor** (*conversion factor ang to nm , = 10*)




#### ran_box(in_top_path: str, in_cnf_path: str, out_cnf_path: str = '', periodic_boundary_condition: str = 'r', nmolecule: int = 1, dens: float = 1.0, threshold: float = None, layer: bool = False, boxsize: float = None, fixfirst: bool = False, seed: float = None, _binary_name: str = 'ran_box', verbose: bool = False, return_command_only: bool = False)

#### red_top(in_top_path: str, atom_selection: str, out_top_path: str, _binary_name: str = 'red_top')
> red_top is a gromos tool to reduce a gromos tool to a certain selection.


* **Parameters**

    
    * **in_top_path** (*str*)


    * **atom_selection** (*str*)


    * **out_top_path** (*str*)


    * **_binary_name** (*str, optional*)



* **Returns**

    out_top_path



* **Return type**

    str



#### rgyr(out_rgyr_path: str, in_coord_path: str, in_top_path: str, atom_selection: str = '1:a', periodic_boundary_condition: str = 'r cog', time: Optional[int] = None, dt: Optional[int] = None, mass_weighted: bool = False, _binary_name: str = 'rgyr')
This wrapper uses rgyr to compute the radius of gyration for a given atom selection.


* **Parameters**

    
    * **out_rgyr_path** (*str*)


    * **in_coord_path** (*str*)


    * **in_top_path** (*str*)


    * **atom_selection** (*str, optional*)


    * **periodic_boundary_condition** (*str, optional*)


    * **time** (*int, optional*)


    * **dt** (*int, optional*) – mass_weighted:bool, optional
    _binary_name: str, optional



* **Returns**

    returns out_rgyr_path



* **Return type**

    str



#### rmsd(in_top_path: str, in_trcs: Union[str, List[str]], atom_selection: str, out_file_path: str, pbc: str = 'r', _binary_name: str = 'rmsd')

* **Parameters**

    
    * **in_top_path**


    * **in_trcs**


    * **atom_selection**


    * **out_file_path**


    * **pbc**


    * **_binary_name**



#### rmsf(in_top_path: str, in_trcs: Union[str, List[str]], atom_selection: str, out_file_path: str, pbc: str = 'r', _binary_name: str = 'rmsf')
> this is a wrapper for gromosPP rmsf programm. (Root mean square fluctuation


* **Parameters**

    
    * **in_top_path** (*str*) – path to topology file


    * **in_trcs** (*Union[str, List[str]]*) – Path OR paths to trc coordinate files


    * **atom_selection** (*str*) – selection of atoms


    * **out_file_path** – out path.


    * **pbc** (*str*) – periodic boundary condition of trc files


    * **_binary_name** (*str*) – binary name of gromos file



* **Returns**

    outpath of the traj



* **Return type**

    str



#### sasa(out_sasa_path: str, in_coord_path: str, in_top_path: str, atom_selection: str = '1:a', sasa_atoms: str = '1:a', probe: str = '4 1.4', periodic_boundary_condition: str = 'r cog', zslice: Optional[float] = None, time: Optional[int] = None, dt: Optional[int] = None, verbose: bool = False, _binary_name: str = 'sasa')
This wrapper uses sasa to compute the solvent accessible surface area (SASA) for a given atom selection. By default,
this is done for the first residue (1:a) with parameters for water (IAC type: 4, radius: 0.4 nm)


* **Parameters**

    
    * **out_sasa_path** (*str*)


    * **in_coord_path** (*str*)


    * **in_top_path** (*str*)


    * **atom_selection** (*str, optional*)


    * **sasa_atoms** (*str, optional*)


    * **probe** (*str, optional*)


    * **periodic_boundary_condition** (*str, optional*)


    * **zslice** (*float, optional*)


    * **time** (*int, optional*)


    * **dt** (*int, optional*)


    * **verbose** (*bool, optional*)


    * **_binary_name** (*str, optional*)



* **Returns**

    returns out_sasa_path



* **Return type**

    str



#### sim_box(in_top_path: str, in_cnf_path: str, in_solvent_cnf_file_path: str, out_cnf_path: str = '', periodic_boundary_condition: str = 'r', gathering_method: str = None, minwall: float = 0.8, threshold: float = None, rotate: str = None, boxsize: bool = False, _binary_name: str = 'sim_box', verbose: bool = False)
When simulating a molecule in solution or in a crystal containing solvent
molecules, the atomic coordinates of the solvent molecules are to be
generated, if they are not available from experiment. Program sim_box can
solvate a solute in a pre-equilibrated box of solvent molecules. The file
specifying the solvent configuration should contain a BOX block with the
dimensions corresponding to the pre-equilibrated density. The solvent
topology is read from the solvent block in the specified topology.


* **Parameters**

    
    * **in_top_path** (*str*) – the path to the input topology file (.top)


    * **in_cnf_path** (*str*) – the path to the input coordinate file (.cnf), which shall be solvated


    * **in_solvent_cnf_file_path** (*str*) – the path to the input coordinate file of the solvent  (.cnf), that shall be used to solvate (checkout pygromos.data.solvent_coordinates for templates)


    * **out_cnf_path** (*str, optional*) – the path to the resulting coordinate (.cnf) file, by default “”


    * **periodic_boundary_condition** (*str, optional*) – describes the boundary condition of the given system in the cnf. (r - rectangle, v - vacuum, ), by default “r”


    * **gathering_method** (*str, optional*) – the gathering method to be used, by default None


    * **minwall** (*float, optional*) – minimum solute to wall distance, by default 0.8


    * **threshold** (*float, optional*) – minimum solvent-solute distance, by default None ->  0.23 nm


    * **rotate** (*str, optional*) – rotate solute: biggest axis along z, second along y, by default None


    * **boxsize** (*bool, optional*) – use boxsize specified in solute coordinate file, by default False


    * **_binary_name** (*str, optional*) – name of the binary, by default “sim_box”


    * **verbose** (*bool, optional*) – stay a while and listen!, by default False



* **Returns**

    return the path to the resulting cnf path.



* **Return type**

    str



#### tser(in_trc_path: str, in_top_path: str, out_csv_path: str, property: str, periodic_boundary_condition: str = 'r', time: float = None, solvent: str = None, normalise_distribution: bool = False, skip_first_n_frames: int = 0, take_each_nth_frame: int = 1, _binary_name: str = 'tser')
> Tser is a gromos programm, that can analyze trajectories.


* **Parameters**

    
    * **in_trc_path** (*str*)


    * **in_top_path** (*str*)


    * **out_csv_path** (*str*)


    * **property** (*str*)


    * **periodic_boundary_condition** (*str, optional*)


    * **time** (*float, optional*)


    * **solvent** (*str, optional*)


    * **normalise_distribution** (*bool, optional*)


    * **skip_first_n_frames** (*int, optional*)


    * **take_each_nth_frame** (*int, optional*)


    * **_binary_name** (*str, optional*)


**WARNING**: missing options: @nots, @dist


* **Returns**

    out_csv_path



* **Return type**

    str


## pygromos.gromos.gromosXX module

FUNCTIONLIB:            wrapper for gromosXX
Description:

> This file contains python wrappers for the bash commandline of gromosXX

Author: Benjamin Schroeder


### _class_ pygromos.gromos.gromosXX.GromosXX(gromosXX_bin_dir: Optional[str] = None, _check_binary_paths: bool = True)
Bases: `pygromos.gromos.gromosXX._GromosXX`

This is the class represents gromosXX.

bin

    This is the path to the folder containing the binaries of gromosXX. If None, the bash enviroment variables  will be used.


### _class_ pygromos.gromos.gromosXX._GromosXX(gromosXX_bin_dir: Optional[str] = None, _check_binary_paths: bool = True)
Bases: `pygromos.gromos._gromosClass._gromosClass`

GromosXX

This is the gromosXX baseclass. This should be inherited by a concrete class that might reimplement some new features, that are version dependent.

bin

    This is the path to the folder containing the binaries of gromosXX. If None, the bash enviroment variables  will be used.


#### \__init__(gromosXX_bin_dir: Optional[str] = None, _check_binary_paths: bool = True)
Constructing a gromosXX object.


* **Parameters**

    
    * **gromosXX_bin_dir** (*str, optional*) – This is the path to the folder containing the binaries of gromosXX. If None, the bash enviroment variables  will be used.


    * **_dont_check_binary** (*bool, optional*) – This flag removes the checks of the binary presence for this obj. This can make sense if system access is slow!, by default False - checks will be made



#### md_run(in_topo_path: str, in_coord_path: str, in_imd_path: str, out_prefix: str, in_pert_topo_path: str = None, in_disres_path: str = None, in_posresspec_path: str = None, in_refpos_path: str = None, in_qmmm_path: str = None, nomp: int = 1, nmpi: int = 1, out_trc: bool = False, out_tre: bool = False, out_trv: bool = False, out_trf: bool = False, out_trs: bool = False, out_trg: bool = False, verbose: bool = False, _binary_name: str = 'md')
This function is a wrapper for gromosXX md_mpi. You can directly execute the gromosXX md_mpi in a bash enviroment here.

**WARNING**: Hybrid jobs are possible, but very difficult to implement correctly to Euler and performance gain is questionable.
If OMP should be used, I suggest the md_run - function.


* **Parameters**

    
    * **in_topo_path** (*str*) – This is the path to the input topology file (x.top)


    * **in_coord_path** (*str*) – This is the path to the input coordinate file (x.cnf)


    * **in_imd_path** (*str*) – This is the path to the input simulation parameter file (x.imd)


    * **out_prefix** (*str*) – This prefix, define the output name.


    * **in_pert_topo_path** (*str, optional*) – This is the path to the pertubation file (x.ptp)


    * **in_disres_path** (*str, optional*) – This is the path to the distance restraint file (x.dat)


    * **in_posresspec_path** (*str, optional*) – This is the path to the position restraint file (x.pos)


    * **in_refpos_path** (*str, optional*) – This is the path to the reference position file (x.rpf)


    * **nomp** (*int, optional*) – How many omp cores shall be used? Prerequesite, gromos was compiled with -enableOMP


    * **nmpi** (*int, optional*) – How many mpi cores shall be used? Prerequesite, gromos was compiled with -enableMPI (and suggested to have -disableOMP)


    * **out_trc** (*bool, optional*) – do you want to output a coordinate trajectory (x.trc) file? (needs also an output number in write block of imd!)


    * **out_tre** (*bool, optional*) – do you want to output a coordinate trajectory (x.trc) file? (needs also an output number in write block of imd!)


    * **out_trs** (*bool, optional*) – do you want to output a coordinate trajectory (x.trc) file? (needs also an output number in write block of imd!)


    * **out_trg** (*bool, optional*) – do you want to output the free energy trajectory (x.trg) file? (needs also an output number in write block of imd!)



* **Returns**

    returns the log_file_path of the simulation.



* **Return type**

    str



* **Raises**

    **ChildProcessError** – If the execution of the simulation fails, this error is raised.



#### repex_run(in_topo_path: str, in_coord_path: str, in_imd_path: str, out_prefix: str, in_pert_topo_path: str = None, in_disres_path: str = None, in_posresspec_path: bool = False, in_refpos_path: bool = False, out_trc: bool = True, out_tre: bool = True, out_trs: bool = False, out_trg: bool = False, out_trf: bool = False, out_trv: bool = False, nomp: int = 1, nmpi: int = 1, verbose: bool = True, _binary_name: str = 'repex_mpi')
This function is a wrapper for gromosXX repex_mpi. You can directly execute the gromosXX repex_mpi in a bash enviroment here.

**WARNING**: Hybrid jobs are possible, but very difficult to implement correctly to Euler and performance gain is questionable.


* **Parameters**

    
    * **in_topo_path** (*str*) – This is the path to the input topology file (x.top)


    * **in_coord_path** (*str*) – This is the path to the input coordinate file (x.cnf)


    * **in_imd_path** (*str*) – This is the path to the input simulation parameter file (x.imd) - needs to contain a Replica Exchange or RE-EDS block.


    * **out_prefix** (*str*) – This prefix, define the output name.


    * **in_pert_topo_path** (*str, optional*) – This is the path to the pertubation file (x.ptp)


    * **in_disres_path** (*str, optional*) – This is the path to the distance restraint file (x.dat)


    * **in_posresspec_path** (*str, optional*) – This is the path to the position restraint file (x.pos)


    * **in_refpos_path** (*str, optional*) – This is the path to the reference position file (x.rpf)


    * **nomp** (*int, optional*) – How many omp cores shall be used? Prerequesite, gromos was compiled with -enableOMP


    * **nmpi** (*int, optional*) – How many mpi cores shall be used? Prerequesite, gromos was compiled with -enableMPI (and suggested to have -disableOMP)


    * **out_trc** (*bool, optional*) – do you want to output a coordinate trajectory (x.trc) file? (needs also an output number in write block of imd!)


    * **out_tre** (*bool, optional*) – do you want to output a coordinate trajectory (x.trc) file? (needs also an output number in write block of imd!)


    * **out_trs** (*bool, optional*) – do you want to output a coordinate trajectory (x.trc) file? (needs also an output number in write block of imd!)


    * **out_trg** (*bool, optional*) – do you want to output the free energy trajectory (x.trg) file? (needs also an output number in write block of imd!)


    * **out_trf** (*bool, optional*) – do you want to output the free energy trajectory (x.trg) file? (needs also an output number in write block of imd!)


    * **out_trv** (*bool, optional*) – do you want to output the free energy trajectory (x.trg) file? (needs also an output number in write block of imd!)


    * **queueing_systems** (*NONE*) – This var is not in use yet! - under development



* **Returns**

    returns the log_file_path of the simulation.



* **Return type**

    str



* **Raises**

    **ChildProcessError** – If the execution of the simulation fails, this error is raised.


## Module contents
