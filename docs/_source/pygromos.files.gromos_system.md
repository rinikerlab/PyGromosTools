---
date: '2022-04-01T13:55:18.499Z'
docname: _source/pygromos.files.gromos_system
images: {}
path: /source-pygromos-files-gromos-system
title: pygromos.files.gromos_system package
---

# pygromos.files.gromos_system package

## Submodules

## pygromos.files.gromos_system.gromos_system module

File:            gromos_system.py
Warnings:        this class is WIP!

Description:

This is a super class for the collection of all Gromos files used inside a project.
Bundle files with the system’s topology, coordinates, input parameters, etc. and
start your simulations from here.

Author: Marc Lehner, Benjamin Ries, Felix Pultar
test


### _class_ pygromos.files.gromos_system.gromos_system.Gromos_System(work_folder: str, system_name: str, rdkitMol: typing.Optional[rdkit.Chem.rdchem.Mol] = None, in_mol2_file: typing.Optional[str] = None, readIn: bool = True, forcefield: pygromos.files.forcefield._generic_force_field._generic_force_field = <pygromos.files.forcefield._generic_force_field._generic_force_field object>, auto_convert: bool = False, adapt_imd_automatically: bool = True, verbose: bool = False, in_smiles: typing.Optional[str] = None, in_residue_list: typing.Optional[typing.List] = None, in_top_path: typing.Optional[str] = None, in_cnf_path: typing.Optional[str] = None, in_imd_path: typing.Optional[str] = None, in_disres_path: typing.Optional[str] = None, in_ptp_path: typing.Optional[str] = None, in_posres_path: typing.Optional[str] = None, in_refpos_path: typing.Optional[str] = None, in_qmmm_path: typing.Optional[str] = None, in_gromosXX_bin_dir: typing.Optional[str] = None, in_gromosPP_bin_dir: typing.Optional[str] = None)
Bases: `object`


#### \__init__(work_folder: str, system_name: str, rdkitMol: typing.Optional[rdkit.Chem.rdchem.Mol] = None, in_mol2_file: typing.Optional[str] = None, readIn: bool = True, forcefield: pygromos.files.forcefield._generic_force_field._generic_force_field = <pygromos.files.forcefield._generic_force_field._generic_force_field object>, auto_convert: bool = False, adapt_imd_automatically: bool = True, verbose: bool = False, in_smiles: typing.Optional[str] = None, in_residue_list: typing.Optional[typing.List] = None, in_top_path: typing.Optional[str] = None, in_cnf_path: typing.Optional[str] = None, in_imd_path: typing.Optional[str] = None, in_disres_path: typing.Optional[str] = None, in_ptp_path: typing.Optional[str] = None, in_posres_path: typing.Optional[str] = None, in_refpos_path: typing.Optional[str] = None, in_qmmm_path: typing.Optional[str] = None, in_gromosXX_bin_dir: typing.Optional[str] = None, in_gromosPP_bin_dir: typing.Optional[str] = None)
> The Gromos_System class is the central unit of PyGromosTools for files and states.
> With this class all files can be read-in or the files can be automatically generated from smiles.
> Additionally to that can all gromos++ functions be used from the Gromos System, so system generation can be easily accomplished.

> if you want to remove all binary checks, do the following:
> >>> from pygromos.files.gromos_system.gromos_system import Gromos_System
> >>> Gromos_System._gromos_noBinary_checks = True
> >>> sys = Gromos_System()


* **Parameters**

    
    * **work_folder** (*str*) – This gives the initial working folder for the system.


    * **system_name** (*str*) – the name of the system, also used as file prefix


    * **in_smiles** (*str, optional*) – Molecule input SMILES for file generation, by default None


    * **in_top_path** (*str, optional*) – input Gromos topology path (.top), by default None


    * **in_cnf_path** (*str, optional*) – input Gromos coordinate path (.cnf), by default None


    * **in_imd_path** (*str, optional*) – input Gromos simulation parameters path (.imd), by default None


    * **in_disres_path** (*str, optional*) – input Gromos distance restraint path (.disres), by default None


    * **in_ptp_path** (*str, optional*) – input pertubation file for free energy calculations (.ptp), by default None


    * **in_posres_path** (*str, optional*) – input position restraints file (.por), by default None


    * **in_qmmm_path** (*str, optional*) – qmmm parameter file (.qmmm), by default None


    * **in_refpos_path** (*str, optional*) – input reference position file (.rpf), by default None


    * **in_gromosXX_bin_dir** (*str, optional*) – path to the binary dir of GromosXX, by default None -> uses the set binaries in the PATH variable


    * **in_gromosPP_bin_dir** (*str, optional*) – path to the binary dir of GromosPP, by default None -> uses the set binaries in the PATH variable


    * **rdkitMol** (*Chem.rdchem.Mol, optional*) – input rdkit Molecule, by default None


    * **in_mol2_file** (*str, optional*) – path to input mol2 file, by default None


    * **readIn** (*bool, optional*) – readIn all provided files?, by default True


    * **Forcefield** (*forcefield_system, optional*) – input PyGromos - forcefield Class , by default forcefield_system()


    * **auto_convert** (*bool, optional*) – automatically convert rdkit MOL and smiles to gromos files, by default False


    * **adapt_imd_automatically** (*bool, optional*) – adjust the input imd file to the GromosSystem, by default True


    * **verbose** (*bool, optional*) – Stay a while and listen!, by default False



* **Raises**

    **Warning** – Rises warning if files are not present.



#### adapt_imd(not_ligand_residues: List[str] = [])

#### _property_ all_file_paths(_: Dict[str, str_ )

#### _property_ all_files(_: Dict[str, pygromos.files._basics._general_gromos_file._general_gromos_file_ )

#### auto_convert()

#### _property_ cnf(_: [pygromos.files.coord.cnf.Cnf](#pygromos.files.coord.cnf.Cnf_ )

#### copy()

#### _property_ disres(_: [pygromos.files.topology.disres.Disres](#pygromos.files.topology.disres.Disres_ )

#### generate_posres(residues: typing.List = [<class 'int'>], keep_residues: bool = True, verbose: bool = False)

#### get_file_paths()
> get the paths of the files in a dict.


* **Returns**

    returns alle file paths, with attribute file name as key.



* **Return type**

    Dict[str, str]



#### get_script_generation_command(var_name: Optional[str] = None, var_prefixes: str = '')

#### _property_ gromosPP(_: [pygromos.gromos.gromosPP.GromosPP](#pygromos.gromos.gromosPP.GromosPP_ )

#### _property_ gromosPP_bin_dir(_: st_ )

#### _property_ gromosXX(_: [pygromos.gromos.gromosXX.GromosXX](#pygromos.gromos.gromosXX.GromosXX_ )

#### _property_ gromosXX_bin_dir(_: st_ )

#### _property_ imd(_: [pygromos.files.simulation_parameters.imd.Imd](#pygromos.files.simulation_parameters.imd.Imd_ )

#### _classmethod_ load(path: Optional[Union[str, _io.FileIO]] = None)
This method stores the Class as binary obj to a given path or fileBuffer.


#### _property_ name(_: st_ )

#### non_ligand_info(_: pygromos.files.coord.cnf.ligands_inf_ )

#### optional_files(_ = {'disres': <class 'pygromos.files.topology.disres.Disres'>, 'posres': <class 'pygromos.files.coord.posres.Position_Restraints'>, 'ptp': <class 'pygromos.files.topology.ptp.Pertubation_topology'>, 'qmmm': <class 'pygromos.files.qmmm.qmmm.QMMM'>, 'refpos': <class 'pygromos.files.coord.refpos.Reference_Position'>_ )

#### parse_attribute_files(file_mapping: Dict[str, str], readIn: bool = True, verbose: bool = False)
> This function sets dynamically builds the output folder, the file objs of this class and checks their dependencies.


* **Parameters**

    **file_mapping** (*Dict[str, Union[str, None]]*) – attribute name: input path



#### _property_ posres(_: [pygromos.files.coord.posres.Position_Restraints](#pygromos.files.coord.posres.Position_Restraints_ )

#### prepare_for_simulation(not_ligand_residues: List[str] = [])

#### protein_info(_: pygromos.files.coord.cnf.protein_inf_ )

#### _property_ ptp(_: [pygromos.files.topology.ptp.Pertubation_topology](#pygromos.files.topology.ptp.Pertubation_topology_ )

#### _property_ qmmm(_: [pygromos.files.qmmm.qmmm.QMMM](#pygromos.files.qmmm.qmmm.QMMM_ )

#### rdkit2Gromos()

#### rdkit2GromosName()

#### rdkitImport(inputMol: rdkit.Chem.rdchem.Mol)

#### read_files(verbose: bool = False)

#### rebase_files()

#### _property_ refpos(_: [pygromos.files.coord.refpos.Reference_Position](#pygromos.files.coord.refpos.Reference_Position_ )

#### required_files(_ = {'cnf': <class 'pygromos.files.coord.cnf.Cnf'>, 'imd': <class 'pygromos.files.simulation_parameters.imd.Imd'>, 'top': <class 'pygromos.files.topology.top.Top'>_ )

#### residue_list(_: Dic_ )

#### save(path: Optional[Union[str, _io.FileIO]] = None, safe: bool = True)
This method stores the Class as binary obj to a given path or fileBuffer.


#### solute_info(_: pygromos.files.coord.cnf.solute_inf_ )

#### solvent_info(_: pygromos.files.coord.cnf.solvent_inf_ )

#### _property_ top(_: [pygromos.files.topology.top.Top](#pygromos.files.topology.top.Top_ )

#### traj_files(_ = {'trc': <class 'pygromos.files.trajectory.trc.Trc'>, 'tre': <class 'pygromos.files.trajectory.tre.Tre'>_ )

#### _property_ trc(_: [pygromos.files.trajectory.trc.Trc](#pygromos.files.trajectory.trc.Trc_ )

#### _property_ tre(_: [pygromos.files.trajectory.tre.Tre](#pygromos.files.trajectory.tre.Tre_ )

#### _property_ work_folder(_: st_ )

#### work_folder_no_update(work_folder: str)

#### write_files(cnf: bool = False, imd: bool = False, top: bool = False, ptp: bool = False, disres: bool = False, posres: bool = False, refpos: bool = False, qmmm: bool = False, mol: bool = False, all: bool = True, verbose: bool = False)
## Module contents
