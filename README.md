Welcome to PyGromosTools
==============================
[//]: # (Badges)


General
-------------
   The aim of the module is to bring Gromos to the Python3 World!
   This repository should make it easier to work with gromos in python and should enable you to write cleaner code.
   Checkout the examples in the example folder.

Content
-------------
* Gromos
    * GromosXX wrapper: for simulation execution
    * GromosPP wrapper: for gromosPP-tool useage
* File handing of following gromos files:
    * coordinatte files:
        .cnf, .trc
    * energy files:
        .tre
    * parameter files:
        .imd
    * replica exchange files:
        repdat.dat
    * topology   
        top, mtb, ifp -> under development     
    * classes for single blocks of each of these files.

* Other Utilities:
    * Automated and system creation (from RDKit or topos)
    * Automated queueing and submission for:
        * Local calculation
        * Cluster calculation (LSF system)
        * Dummy 
    * Bash wrappers for GROMOS
    * amino acid library
    
        
## Specifications
 * Python >=3.7:
 * requires: numpy, scipy, pandas, rdkit

 * optional: openforcefield for off and serenityff functions


## SETUP

### Installing the package - _UNDER Developmten_

    cd pygromos
    python setup.py install

### Package Development
For using this repository and developing it, clone it into a directory on your machine and add the path to the repo to your python path.

    PYTHONPATH=${PYTHONPATH}:/path/to/pygromos/containint/folder/pygromos
    
If you are using Anaconda, you might need to use to this instead:
   
    conda develop -n <EnvironmentName> /path/to/pygromos/containint/folder/pygromos

Please if your writing code for this repository, first develop it on an own branch.

     git branch <MyBranch>    #generate your branch
     git checkout <MyBranch>  #switch to your branch
     git merge master   #for adding new features from master to your branch
     
If you implemented in your branch features, that you would like to share, just issue a merge/pull request with the master branch on gitlab.
If you find a bug or have an feature request, please raise an Issue.

P.s.: I can recommend Pycharm or VisualStuudioCode from dstar, for exploring the repository.

### Copyright

Copyright (c) 2020, Benjamin Ries, Salome Rieder, Marc Lehner 


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.3.

