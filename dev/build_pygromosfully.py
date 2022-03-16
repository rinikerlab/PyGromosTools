
import os, sys

abs_file = os.path.abspath(__file__)
package_path = os.path.dirname(os.path.dirname(abs_file))
connda_env_path = os.path.dirname(abs_file)+"/conda_envs/dev_env_withGromos.yaml"
env_name = 'pygromosWithGrom'

#Conda commands
conda_install_env = "conda env create -f "+connda_env_path
conda_setDevelop = "conda develop "+package_path
conda_activation = "conda activate "+env_name

os.system(conda_install_env)
os.system(conda_setDevelop)
os.system(conda_activation)

#Compile gromos
import pygromos
from pygromos.gromos.compile_gromos import default_install

default_install()
