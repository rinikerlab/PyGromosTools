---
date: '2022-04-01T13:55:18.499Z'
docname: _source/pygromos.files.forcefield.amber
images: {}
path: /source-pygromos-files-forcefield-amber
title: pygromos.files.forcefield.amber package
---

# pygromos.files.forcefield.amber package

## Submodules

## pygromos.files.forcefield.amber.amberff module


### _class_ pygromos.files.forcefield.amber.amberff.AmberFF(name: str = 'amber', path_to_files: Optional[str] = None, auto_import: bool = True, verbose: bool = False)
Bases: `pygromos.files.forcefield._generic_force_field._generic_force_field`


#### auto_import_ff()

#### clean(_ = Fals_ )

#### create_cnf(mol: str, in_cnf: Optional[[pygromos.files.topology.top.Top](#pygromos.files.topology.top.Top)] = None, work_folder: Optional[str] = None, in_mol2_file: Optional[str] = None, gromosPP: Optional[[pygromos.gromos.gromosPP.GromosPP](#pygromos.gromos.gromosPP.GromosPP)] = None)

#### create_mol2(mol: Optional[str] = None)

#### create_top(mol: str, in_top: [pygromos.files.topology.top.Top](#pygromos.files.topology.top.Top), in_mol2_file: Optional[str] = None, work_folder: Optional[str] = None, gromosPP: Optional[[pygromos.gromos.gromosPP.GromosPP](#pygromos.gromos.gromosPP.GromosPP)] = None)

#### solvate(_ = Fals_ )

#### solventbox(_ = 'TIP3PBOX_ )

### _class_ pygromos.files.forcefield.amber.amberff.amber2gromos(in_mol2_file: str, mol: rdkit.Chem.rdchem.Mol, gromosPP: [pygromos.gromos.gromosPP.GromosPP](#pygromos.gromos.gromosPP.GromosPP), forcefield: pygromos.files.forcefield._generic_force_field._generic_force_field, work_folder: str = '.', solvate: bool = False, solventbox: Optional[str] = None, clean: bool = False)
Bases: `object`


#### \__init__(in_mol2_file: str, mol: rdkit.Chem.rdchem.Mol, gromosPP: [pygromos.gromos.gromosPP.GromosPP](#pygromos.gromos.gromosPP.GromosPP), forcefield: pygromos.files.forcefield._generic_force_field._generic_force_field, work_folder: str = '.', solvate: bool = False, solventbox: Optional[str] = None, clean: bool = False)
uses the ambertools programs antechamber, parmchk, and tleap together with
the GROMOS++ program amber2gromos to generate a GROMOS topology and coordinate file
for a given molecule


* **Parameters**

    
    * **in_mol2_file** (*str*) – mol2 file for molecule to be parameterized


    * **mol** (*Chem.rdchem.Mol*) – rdkit molecule of the molecule in in_mol2_file


    * **gromosPP** (*GromosPP*)


    * **forcefield** (*forcefield_system*)


    * **work_folder** (*str, optional*) – where to generate the topology + cnf (default: “.”)


    * **solvate** (*bool, optional*) – should the topology be solvated? (default: False)


    * **solventbox** (*str, optional*) – what solvent should be used for solvation? e.g. TIP3PBOX or CHCL3BOX (default: TIP3PBOX)


    * **clean** (*bool, optional*) – should temporary ambertool files be removed after parameterization? (default: False)



#### amber2gromos()
converts an AMBER (GAFF) parameter and crd file to a GROMOS top and cnf file


#### antechamber()
executes the ambertools program antechamber


#### cleanup()
removes temporary parmchk, antechamber, and tleap directories


#### get_gromos_coordinate_file()

* **Returns**

    returns the path to the converted GROMOS coordinate file



* **Return type**

    str



#### get_gromos_topology()

* **Returns**

    returns the path to the converted GROMOS topology



* **Return type**

    str



#### parmchk()
executes the ambertools program parmchk


#### tleap()
executes the ambertools program tleap

## Module contents
