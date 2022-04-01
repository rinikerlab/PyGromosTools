# Installation guide

## Environment Setup

To use PyGromosTools, you need an environment setup, like specified in the conda_env.yml.
You simply can build such an environment with conda from the provided yml enviroment file shipped with the package:
```bash
conda env create -f conda_env.yml
```
This will build an envrionment named pygromos.

**Note**: if you want to develop PyGromosTools, checkout the dev/conda_envs folder, there you find environments, that contain all packages for constructing PyGromosTools.

## PygromosTools SETUP
### Installing the package
This is recommended if you want to use PyGromosTools, but not develop it:

```bash
    cd PyGromosTools
    python setup.py install
```
Make sure, that you have the GROMOS binaries around, as currently the binaries can not be shipped with the package, as it is not open-source.

### Package Development
For using this repository and developing it, clone it into a directory on your machine.

Set your python envs path to the following:
```bash
conda develop -n <EnvironmentName> /path/to/pygromos/containint/folder/pygromos
```

If you write code for this repository, first develop it on a seperated branch or a fork (see GitHub fork).

```bash
     git branch <MyBranch>    #generate your branch
     git checkout <MyBranch>  #switch to your branch
```

If you implemented in your branch features, that you would like to share, just issue a merge/pull request with the master branch on GitHub. We are looking forward to it!

**IMPORTANT:**
If you decide to do a merge/pull request, please make sure to follow the coding style guidelines in styleguide.md and make sure, that your code is well documented and passes pre-commit and the tests.

If you find a bug or have an feature request, please raise an Issue.
