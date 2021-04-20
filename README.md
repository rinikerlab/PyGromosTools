![](.img/PyGromosToolsBanner.png)

Welcome to PyGromosTools
==============================
[//]: # (Badges)
[![CI](https://github.com/rinikerlab/PyGromosTools/actions/workflows/CI.yaml/badge.svg)](https://github.com/rinikerlab/PyGromosTools/actions/workflows/CI.yaml)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/rinikerlab/PyGromosTools.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/rinikerlab/PyGromosTools/context:python)
[![Documentation](https://img.shields.io/badge/Documentation-here-white.svg)](https://rinikerlab.github.io/PyGromosTools/)

General
-------------
   The aim of the module is to bring GROMOS to the Python3 World!
   This repository should make it easier to work with GROMOS in Python and should enable the user to write cleaner, more reliable and adaptable code.

   General information about functions can be found in our wimi. Usage examples for many general functions and their relations are shown in Jupyter notebooks in the examples in the example folder.

Content

-------------

* GROMOS wrappers
  * GromosXX wrapper: for simulation execution
  * GromosPP wrapper: for GROMOS++ program usage

* File handling of all GROMOS file types for automated creation/modification/analysis :
  * coordinate files CNF:
    * read and analyse CNF files
    * generate CNF files from RDKit
    * generate CNF files from SDF

    ```python
    cnf = Cnf(input_value="file_name")
    print(cnf.GENBOX)
    ```

  * topology files:
    * create topologies from a forcefield
      * GROMOS 2016H66 / 54A7
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
    * analyse trajectories with Pandas data frames
    * standard analysis like RSMD, RDF, ... for trc
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

* Automation and file management system `gromos_system`
  * offers clean file management for simulations
  * offers a high level of automation
  * equiped with simulation queuing system
  * includes many force fields

  ```python
  ff=forcefield_system(name="off")
  gsys = Gromos_System(work_folder="dir", in_smiles="c1ccccc1", auto_convert=True, Forcefield=ff)
  print(gsys)
  ```

* Other utilities:
  * Automated queueing and submission for:
    * Local calculation
    * Cluster calculation (LSF system)
    * Dummy
  * Bash wrappers for GROMOS
  * Amino acid library

General Information
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

