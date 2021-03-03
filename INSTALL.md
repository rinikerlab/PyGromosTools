# Installation guide

## Specifications

Make sure to have the following packages installed:
 * Python >=3.7:
 * requires: numpy, scipy, pandas, rdkit
 * optional: openforcefield for OpenForceField and Serenityff functions

```bash
#with pip installer:
pip install <package_name> 
#or if you use conda:
conda install <package_name>
```

## SETUP

### Installing the package

```bash
    cd pygromos
    python setup.py install
```

### Package Development

For using this repository and developing it, clone it into a directory on your machine and add the path to the repo to your python path.

```bash
PYTHONPATH=${PYTHONPATH}:/path/to/pygromos/containint/folder/pygromos
```

If you are using Anaconda, you might need to use to this instead:

```bash
conda develop -n <EnvironmentName> /path/to/pygromos/containint/folder/pygromos
```

Please if your writing code for this repository, first develop it on an own branch.

```bash
     git branch <MyBranch>    #generate your branch
     git checkout <MyBranch>  #switch to your branch
     git merge master   #for adding new features from master to your branch
```

If you implemented in your branch features, that you would like to share, just issue a merge/pull request with the master branch on gitlab.

If you find a bug or have an feature request, please raise an Issue.

P.s.: I can recommend Pycharm or VisualStuudioCode, for exploring the repository.
