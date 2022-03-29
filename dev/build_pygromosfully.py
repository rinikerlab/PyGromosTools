import os
from collections import OrderedDict
from datetime import datetime

timings = OrderedDict({})

abs_file = os.path.abspath(__file__)
package_path = os.path.dirname(os.path.dirname(abs_file))
connda_env_path = os.path.dirname(abs_file) + "/conda_envs/dev_env_withGromos.yaml"
env_name = "pygromosWithGrom"

# SetUp Conda Environment
# Conda commands
timings["conda_env_build_start"] = datetime.now()
conda_install_env = "conda env create -f " + connda_env_path
conda_setDevelop = "conda develop -n " + env_name + " " + package_path

print("Start Conda Env Build: \n\n")
os.system(conda_install_env)
os.system(conda_setDevelop)

timings["conda_env_build_end"] = datetime.now()
conda_duration = timings["conda_env_build_end"] - timings["conda_env_build_start"]
print("\n\n Finished Conda Env Build: \n\n")
print("CONDA BUILD - Duration:", conda_duration)


# Compile gromos
print("Start Gromos Build: \n\n")
from pygromos.gromos import compile_gromos  # noqa: E402

timings["gromos_build_start"] = datetime.now()
os.system("conda run -v --no-capture-output --live-stream -n " + env_name + " python " + compile_gromos.__file__)
timings["gromos_build_end"] = datetime.now()
grom_duration = timings["gromos_build_end"] - timings["gromos_build_start"]
print("\n\n Finished Gromos Build: \n\n")
print("Gromos BUILD:", grom_duration)

# TIMINGS Printout
print("\n" + ">" * 10 + " TIMINGS:")
for key, val in timings.items():
    print(key, val)

# get duration
print("\n\n" + ">" * 10 + " DURATION:")
durations = OrderedDict({})
keys = list(timings.keys())
for key in keys:
    prefix = "_".join(key.split("_")[:-1])
    if "start" in key:
        print(key)
        print(keys)
        end_key = list(filter(lambda x: x == (prefix + "_end"), keys))[0]
    else:
        pass
    durations[prefix] = timings[end_key] - timings[key]
    print(prefix, str(timings[end_key] - timings[key]))
