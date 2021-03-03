Welcome to PyGromosTools
==============================
[//]: # (Badges)
[![CI](https://github.com/SchroederB/PyGromosTools/actions/workflows/CI.yaml/badge.svg)](https://github.com/SchroederB/PyGromosTools/actions/workflows/CI.yaml)
[![Documentation](https://img.shields.io/badge/Documentation-here-white.svg)](https://schroederb.github.io/PyGromosTools/)
General
-------------
   The aim of the module is to bring Gromos to the Python3 World!
   This repository should make it easier to work with gromos in python and should enable you to write cleaner, more reliable and adapteble code.

   General informations about functions can be found in our wimi and usage example for many general functions and theire relations are shown in jupyter notebooks in the examples in the example folder.

Content

-------------

* Gromos wrappers
  * GromosXX wrapper: for simulation execution
  * GromosPP wrapper: for gromosPP-tool useage

* File handing of all gromos file types for automated creation/modifications/analysis :
  * coordinatte files CNF:
    * read and analysie CNF files
    * generate CNF files from RDKit
    * generate CNF files from SDF

    ```python
    cnf = Cnf(input_value="file_name")
    print(cnf.GENBOX)
    ```

  * topology files:
    * create topologies from a forcefield
      * Gromos 2016H66 / 54A7
      * OpenForceField
      * SerenityForceField
    * modify topologies
      * add new atoms
      * modify force parameters

    ```python
    top = Top(input_value="file_path")
    top.add_new_SOLUTEATOM(ATNM=42)
    print(top)
    ```

  * simulation parameter files IMD
    * a wide option of templates provided
    * modify IMD files to fit your simulation

    ```pythons
    imd = Imd(input_value="file_path")
    imd.INITIALISE.TEMPI = 137
    print(imd)
    ```

  * trajectories (tre, trc, trg, ...)
    * analyse trajectories with pandas dataframes
    * standart analysis like RSMD, RDF, ... for trc
    * auto saving of results for later use as hdf5
    * ene_ana like tools for tre
    * easy to add costume analysis tools

    ```python
    trc = Trc(inpout_value="file_path")
    print(trc.rmsd().mean())
    ```

  * replica exchange files:
        repdat.dat
  * classes for single blocks of each of these files.

* Automatin and file managment system `gromos_system`
  * offers clean file managment for simulations
  * offers a high level of automation
  * equiped with simulation queuing system
  * includes many forcefields

  ```python
  ff=forcefield_system(name="off")
  gsys = Gromos_System(work_folder="dir", in_smiles="c1ccccc1", auto_convert=True, Forcefield=ff)
  print(gsys)
  ```

* Other Utilities:
  * Automated queueing and submission for:
    * Local calculation
    * Cluster calculation (LSF system)
    * Dummy
  * Bash wrappers for GROMOS
  * amino acid library

General Informations
-------------

### Specifications

 * Python >=3.7:
 * requires: numpy, scipy, pandas, rdkit

 * optional: openforcefield for OpenForceField and Serenityff functions

### SETUP

see INSTALLATION.md file for more informations

### Copyright

Copyright (c) 2020, Benjamin Ries, Marc Lehner, Salome Rieder

### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.3.

