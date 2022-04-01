---
date: '2022-04-01T13:55:18.499Z'
docname: _source/pygromos.files.coord
images: {}
path: /source-pygromos-files-coord
title: pygromos.files.coord package
---

# pygromos.files.coord package

## Submodules

## pygromos.files.coord.cnf module


### _class_ pygromos.files.coord.cnf.Cnf(in_value: Optional[Union[str, dict]], clean_resiNumbers_by_Name: bool = False, verbose: bool = False, _future_file: bool = False)
Bases: `pygromos.files._basics._general_gromos_file._general_gromos_file`

This class is a representation of the gromos .cnf coordinate files. It
allows reading, analysis and modifying of the coordinate files.

is a child of general_gromos_file


#### GENBOX(_: [pygromos.files.blocks.coord_blocks.GENBOX](#pygromos.files.blocks.coord_blocks.GENBOX_ )

#### LATTICESHIFTS(_: [pygromos.files.blocks.coord_blocks.LATTICESHIFTS](#pygromos.files.blocks.coord_blocks.LATTICESHIFTS_ )

#### PERTDATA(_: [pygromos.files.blocks.coord_blocks.PERTDATA](#pygromos.files.blocks.coord_blocks.PERTDATA_ )

#### POSITION(_: [pygromos.files.blocks.coord_blocks.POSITION](#pygromos.files.blocks.coord_blocks.POSITION_ )

#### TIMESTEP(_: pygromos.files.blocks._general_blocks.TIMESTE_ )

#### TITLE(_: [pygromos.files.blocks._general_blocks.TITLE](#pygromos.files.blocks.pertubation_blocks.TITLE_ )

#### VELOCITY(_: [pygromos.files.blocks.coord_blocks.VELOCITY](#pygromos.files.blocks.coord_blocks.VELOCITY_ )

#### add_empty_box()

#### add_residue_positions(coords: object)
This function adds all residues of an coords file to @DEVELOP

This is a very crude functio at the moment! It only takes positions
of a residue and merges them! if there are residues with the same name,
this might lead to problems, as clean_posiresnumbyname function is not
sensitive for that! todo: make more robust! bschroed


* **Parameters**

    **coords** (*function_libs.gromos.files.coord.Cnf obj*) – object - CNF object



#### atom_ref_pos_block(_: [pygromos.files.blocks.coord_blocks.REFPOSITION](#pygromos.files.blocks.coord_blocks.REFPOSITION_ )

#### center_of_geometry(selectedAtoms: Optional[List[int]] = None)
calculates the center of geometry for asingle molecule or the selected Atoms


* **Returns**

    cog



* **Return type**

    list



#### clean_posiResNums()
This function recount the Residue number with respect to residue name and residue number.

**WARNING**: only in “[Position_BLOCK!@development](mailto:Position_BLOCK!@development)


* **Return type**

    None



#### count_residue_atoms(resID: int = False, resName: str = False, verbose=False)

#### createRDKITconf(mol: rdkit.Chem.rdchem.Mol, conversionFactor: float = 0.1)
creates a PyGromosTools CNF type from a rdkit molecule. If a conformation exists the first one will be used.


* **Parameters**

    
    * **mol** (*Chem.rdchem.Mol*) – Molecule, possibly with a conformation


    * **conversionFactor** (*float*) – the factor used to convert length from rdkit to Gromos
    (default: angstrom -> nano meter = 0.1)



#### delete_atom(resID: int = False, resName: str = False, atomID: int = False, atomType: str = False)

#### delete_residue(resID: int = False, resName: str = False, verbose=False)
this function is deleting residues from a cnf file with taking following _blocks into account:

    “POSITION”, “VELOCITY”, “LATTICESHIFTS”, “REFPOSITION”

it additionally recounts all atomIds and residueIDs afterwards.
you can provide a residue ID or a residue Name or both (than only exact match will be deleted).


* **Parameters**

    
    * **resID** (*int*) – Id of the residue to be deleted


    * **resName** (*str*) – Name of the residue to be deleted


    * **verbose** (*bool*) – Text… lots of it!



* **Returns**

    0 if succesfull



* **Return type**

    int



#### gen_possrespec(residues: Union[Dict[str, List[int]], List[int]], keep_residues: bool = True, verbose: bool = False)
> This function writes out a gromos file, containing a atom list. that is to be position restrained!
> Raises not Implemented error, if a input variant of the residues


* **Parameters**

    
    * **residues** (*dict, list*) – residues to be restrained (dict containing resname:[res ids] or list of resnames or residue IDs)


    * **keep_residues** (*bool, optional*) – should the passed residues be kept or deleted?


    * **verbose** (*bool, optional*) – loud and noisy?



* **Returns**

    
    * *str* – out_path


    * *Position_Restraints* – posrespec-file-obj




#### gen_refpos()

#### generate_position_restraints(out_path_prefix: str, residues: Union[Dict[str, List[int]], List[int]], verbose: bool = False)
> This function generates position restraints for the selected residues.


* **Parameters**

    
    * **out_path_prefix** (*str*) – target path prefix for the out files.


    * **residues** (*dict or list*) – residues to be restrained (dict containing resname:[res ids] or list of resnames or residue IDs)


    * **verbose** (*bool, optional*) – Loud and noisy, by default False



* **Returns**

    return the two resulting paths to the generated files.



* **Return type**

    Tuple[str, str]



#### get_atomP(atomID: Optional[int] = None, atomType: Optional[str] = None, resID: Optional[int] = None, resName: str = False)

#### get_atom_coordinates()
> This function returns a np.array containing all system coordinates.


* **Returns**

    dims are atoms[x,y,z]



* **Return type**

    np.array



#### get_atoms_distance(atomI: Optional[Union[int, List[int]]] = None, atomJ: Optional[int] = None, atoms: Optional[Union[List[int], Dict[int, int]]] = None)

#### get_density(mass: float = 1, autoCalcMass: bool = False)
> This function calculates the density of the cnf.


* **Returns**

    density of the cnf.



* **Return type**

    float



#### get_last_atomID()
get_last atom

    A very simple convenience function that returns the last atom


* **Returns**

    Returns the last atom of the system.



* **Return type**

    int



#### get_mass()
> This function calculates the mass of the cnf.


* **Returns**

    mass of the cnf.



* **Return type**

    float



#### get_mdtraj()

#### get_pdb(rdkit_ready: bool = False, connectivity_top=None)
> translate cnf to pdb.


* **Parameters**

    
    * **rdkit_ready** (*bool, optional*) – str output was tested with RDKIT (default: False)


    * **connectivity_top** (*top.Top, optional*) – if the pygromos top class is provided (containing a BOND block), then the pdb gets a connect block.



* **Returns**

    pdb str.



* **Return type**

    str



#### get_residues(verbose: bool = False)
This function is getting all residues of the used cnf file.

    it gives back,


* **Parameters**

    **verbose** (*bool, optional*) – texty?



* **Returns**

    returns dict containing all residues, numbers and atoms.



* **Return type**

    Dict[str, Dict[str, Any]]



#### get_system_information(not_ligand_residues: List[str] = [], ligand_resn_prefix: Optional[Union[str, List[str]]] = None, solvent_name: str = 'SOLV')
This function utilizes a dictionary containing all residues and atom numbers (e.g. cnf.get_residues()) and modifies them such, that the result can be used to set up a standard REEDS gromos_simulation


* **Parameters**

    
    * **residues** (*Dict[str, Dict[int,int]]*) – input a cnf.residues:dict that shall be cleaned and return a reduced form for parameter file.


    * **not_ligand_residues** (*List[str]*) – here all molecules, that are not considered as ligand or protein.


    * **ligand_resn_prefix** (*List[str]*) – here all molecules, that are considered as ligand are listed.



* **Returns**

    
    * *Dict[str, Dict[int,int]]* – cleaned_residue dict


    * *NamedTuple* – ligands


    * *NamedTuple* – protein


    * *NamedTuple* – non_ligands




#### get_volume()
> This function calculates the volume of the cnf.


* **Returns**

    volume of the cnf.



* **Return type**

    float



#### get_xyz()
> translate cnf to xyz


* **Returns**

    in xyz format



* **Return type**

    str



#### read_file()

* **Parameters**

    **path**



#### recenter_pbc()
This function is shifting the coordinates such that the solute molecule is centered, if it is placed in the corners.
However, this function might break down with more complicated situations. Careful the function is not fail safe ;)


#### recreate_view()

#### rename_residue(new_resName: str, resID: int = False, resName: str = False, verbose: bool = False)
this function is renaming residues from a cnf file with taking following _blocks into account:

    “POSITION”, “VELOCITY”

it additionally recounts all atomIds and residueIDs afterwards.
you can provide a residue ID or a residue Name or both (than only exact match will be deleted).


* **Parameters**

    
    * **new_resName** (*str*) – new Name of the residue


    * **resID** (*int, optional*) – Id of the residue to be renamed


    * **resName** (*str, optional*) – Name of the residue to be renamed


    * **verbose** (*bool, optional*) – Text… lots of it!



* **Returns**

    0 if succesfull



* **Return type**

    int



#### residues(_: Dict[str, Dict[str, int]_ )

#### rotate(rotationCenter: Optional[numpy.array] = None, selectedAtoms=None, alpha: float = 0, beta: float = 0, gamma: float = 0)

#### shift_periodic_boundary()
This function is shifting the coordinates such that the solute molecule is centered, if it is placed in the corners.
However, this function might break down with more complicated situations. Careful the function is not fail safe ;)


#### supress_atomPosition_singulrarities()
This function adds a very small deviation to the position of an atom, dependent on the atom number.
This might be needed to avoid singularities in gromosXX.


* **Return type**

    None



#### _property_ view(_: nglview.widget.NGLWidge_ )

#### write_pdb(out_path: str)
> This function converts the atom POS db of the traj into a pdb traj.


* **Parameters**

    **out_path** (*str*) – path, were the file should be written to.



* **Returns**

    outpath of the file



* **Return type**

    str



#### write_possrespec(out_path: str, residues: dict, verbose: bool = False)
This function writes out a gromos file, containing a atom list. that is to be position restrained!

    Raises not Implemented error, if a input variant of the residues


* **Parameters**

    
    * **out_path** (*str*) – path to the outputfile


    * **residues** (*dict, list*) – residues to be restrained (dict containing resname:[res ids] or list of resnames or residue IDs)


    * **verbose** (*bool, optional*) – loud and noisy?



* **Raises**

    **NotImplementedError** – 



* **Returns**

    
    * *str* – out_path


    * *Position_Restraints* – posrespec-file-obj




#### write_refpos(out_path: str)

#### write_xyz(out_path: str)
> This function converts the atom POS db of the traj into a xyz structure.


* **Parameters**

    **out_path** (*str*) – path, were the file should be written to.



* **Returns**

    outpath of the file



* **Return type**

    str



### pygromos.files.coord.cnf.non_ligand_infos()
alias of `pygromos.files.coord.cnf.ligands_info`


### pygromos.files.coord.cnf.protein_infos()
alias of `pygromos.files.coord.cnf.protein_info`


### pygromos.files.coord.cnf.solute_infos()
alias of `pygromos.files.coord.cnf.solute_info`


### pygromos.files.coord.cnf.solvent_infos()
alias of `pygromos.files.coord.cnf.solvent_info`

## pygromos.files.coord.posres module


### _class_ pygromos.files.coord.posres.Position_Restraints(in_value: Union[str, dict, pygromos.utils.typing.Position_Restraints_Type, pygromos.utils.typing.Cnf_Type], clean_resiNumbers_by_Name: bool = False, verbose: bool = False, _future_file: bool = False)
Bases: `pygromos.files.coord.cnf.Cnf`

This class is a representation of the gromos .cnf coordinate files. It
allows reading, analysis and modifying of the coordinate files.

is a child of general_gromos_file


#### LATTICESHIFTS(_: [pygromos.files.blocks.coord_blocks.LATTICESHIFTS](#pygromos.files.blocks.coord_blocks.LATTICESHIFTS_ )

#### POSRESSPEC(_: [pygromos.files.blocks.coord_blocks.POSRESSPEC](#pygromos.files.blocks.coord_blocks.POSRESSPEC_ )

#### TITLE(_: [pygromos.files.blocks._general_blocks.TITLE](#pygromos.files.blocks.pertubation_blocks.TITLE_ )

#### path(_: st_ )

#### residues(_: Dict[str, Dict[str, int]_ )
## pygromos.files.coord.refpos module


### _class_ pygromos.files.coord.refpos.Reference_Position(in_value: Union[str, dict, pygromos.utils.typing.Reference_Position_Type, pygromos.utils.typing.Cnf_Type], verbose: bool = False, _future_file: bool = False)
Bases: `pygromos.files.coord.cnf.Cnf`

This class is a representation of the gromos .cnf coordinate files. It
allows reading, analysis and modifying of the coordinate files.

is a child of general_gromos_file


#### GENBOX(_: [pygromos.files.blocks.coord_blocks.GENBOX](#pygromos.files.blocks.coord_blocks.GENBOX_ )

#### LATTICESHIFTS(_: [pygromos.files.blocks.coord_blocks.LATTICESHIFTS](#pygromos.files.blocks.coord_blocks.LATTICESHIFTS_ )

#### REFPOSITION(_: [pygromos.files.blocks.coord_blocks.REFPOSITION](#pygromos.files.blocks.coord_blocks.REFPOSITION_ )

#### TITLE(_: [pygromos.files.blocks._general_blocks.TITLE](#pygromos.files.blocks.pertubation_blocks.TITLE_ )

#### path(_: st_ )

#### residues(_: Dict[str, Dict[str, int]_ )
## Module contents
