import os
from pygromos.utils import bash

root_dir = os.path.dirname(__file__)    
gromosXX_path = root_dir+"/gromosXX/gromosXX"
gromosPP_path = root_dir+"/gromosPlsPls/gromos++"
gromosXX_build_path = gromosXX_path+"/build"
gromosPP_build_path = gromosPP_path+"/build"

nCores=1
make_clean = False
configure_options = {}
recompile = False
recompile_from_scratch = False

do_gromosXX=False
do_gromosPP=True
print(root_dir)


if(make_clean):
    if(do_gromosPP):    bash.remove_file(gromosPP_build_path, recursive=True)
    if(do_gromosXX): bash.remove_file(gromosXX_build_path, recursive=True)     
else:
    #gromosXX
    if(do_gromosXX):
        if(not os.path.exists(gromosXX_build_path) and not recompile_from_scratch):
            os.chdir(gromosXX_path)
            print(os.getcwd())
            # Initial config
            os.system("./Config.sh")

            # Configure
            bash.make_folder(gromosXX_build_path)
            os.chdir(gromosXX_build_path)
            print("\n\n CONFIGURE \n\n")
            os.system("../configure")

        if(not os.path.exists(gromosXX_build_path+"/bin") or recompile or recompile_from_scratch):
            # Compile
            print("\n\n MAKE \n\n")
            x = os.system("make -j"+str(nCores))
            print("OUT: ", x)
            # Create Binaries
            print("\n\n INSTALL \n\n")
            os.system("make -j"+str(nCores)+" install")
            
    #gromosPP
    if(do_gromosPP):
        if(not os.path.exists(gromosPP_build_path) and not recompile_from_scratch):
            os.chdir(gromosPP_path)
            print(os.getcwd())
            # Initial config
            os.system("./Config.sh")

            # Configure
            bash.make_folder(gromosPP_build_path)
            print(os.getcwd())
            os.chdir(gromosPP_build_path)

            print("\n\n CONFIGURE \n\n")
            os.system("../configure")

        if(not os.path.exists(gromosPP_build_path+"/bin") or recompile or recompile_from_scratch):
            # Compile
            print("\n\n MAKE \n\n")
            os.chdir(gromosPP_build_path)
            print(os.getcwd())
            x = os.system("make -j"+str(nCores))
            print("OUT: ", x)

            # Create Binaries
            print("\n\n INSTALL \n\n")
            os.system("make -j"+str(nCores)+" install")

