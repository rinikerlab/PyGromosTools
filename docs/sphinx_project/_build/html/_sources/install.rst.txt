
Installation guide
==================

Environment Setup
-----------------

To use PyGromosTools, you need an environment setup, like specified in the conda_env.yml.
You simply can build such an environment with conda from the provided yml enviroment file shipped with the package:

.. code-block:: bash

   conda env create -f conda_env.yml

This will build an envrionment named pygromos.

**Note**\ : if you want to develop PyGromosTools, checkout the dev/conda_envs folder, there you find environments, that contain all packages for constructing PyGromosTools.

PyGromosTools Setup
-------------------

**Note**\ : If you want to use PyGromosTools with the GROMOS package (gromos.net), you need to generate the binaries for this code seperately and then can use them with PyGromosTools.

Using the Package
^^^^^^^^^^^^^^^^^

This is recommended if you want to use PyGromosTools, but not develop it. First activate the correct python environment, then:

.. code-block:: bash

       cd PyGromosTools
       python setup.py install

Make sure, that you have the GROMOS binaries around, as currently the binaries can not be shipped with the package, as it is not open-source.

Developing the Package
^^^^^^^^^^^^^^^^^^^^^^

For using this repository and developing it, clone it into a directory on your machine.

Build your python environment, activate it and add the path to the repository root to your environment:

.. code-block:: bash

   conda develop -n <EnvironmentName> /path/to/pygromos/containint/folder/pygromos

If you write code for this repository, first develop it on a seperated branch or a fork (see GitHub fork).

.. code-block:: bash

        git branch <MyBranch>    #generate your branch
        git checkout <MyBranch>  #switch to your branch

If you implemented in your branch features, that you would like to share, just issue a merge/pull request with the master branch on GitHub. We are looking forward to it!

**IMPORTANT:**
If you decide to do a merge/pull request, please make sure to follow the coding style guidelines in styleguide.md and make sure, that your code is well documented and passes pre-commit and the tests.
