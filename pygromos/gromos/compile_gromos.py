import os
from pygromos.utils import bash 
from pygromos.utils.utils import spacer, spacer2

def install_gromos(root_dir:str, 
                   gromosXX_with_mpi:bool=False, gromosXX_with_omp:bool=False, gromosXX_with_cuda:bool=False, 
                   gromosPP_with_omp=False, 
                   do_compile:bool=True, do_clean:bool=False, recompile:bool=False, recompile_from_scratch:bool=False, 
                   _do_gromosPP:bool=True, _do_gromosXX:bool=True,
                   with_debug:bool=False, silent:bool=False, nCore:int=1):
    if(do_clean):
        if(do_gromosPP):
            print(spacer+"\n CLEAN GROMOSPP \n"+spacer)
            gromosPP_path = root_dir+"/gromosPlsPls/gromos++"
            gromosPP_build_path = gromosPP_path+"/build"
            if(os.path.exists(gromosPP_build_path)):
                bash.remove_file(gromosPP_build_path, recursive=True)
            
        if(do_gromosXX): 
            print(spacer+"\n CLEAN GROMOSXX \n"+spacer)
            gromosXX_path = root_dir+"/gromosXX/gromosXX"
            gromosXX_build_path = gromosXX_path+"/build"
            if(os.path.exists(gromosXX_build_path)):
                bash.remove_file(gromosXX_build_path, recursive=True)     
    
    if(do_compile):
        binary_dir = root_dir
        #gromosXX
        if(_do_gromosXX):
            print(spacer+"\n BUILD GROMOSXX \n"+spacer)
            gromosXX_path = root_dir+"/gromosXX/gromosXX"
            gromosXX_build_path = gromosXX_path+"/build"

            if(not os.path.exists(gromosXX_build_path) and recompile_from_scratch):
                _configure_gromosXX_autotools(build_dir=gromosXX_build_path,  with_mpi=gromosXX_with_mpi, with_omp=gromosXX_with_omp, with_cuda=gromosXX_with_cuda)

            if(not os.path.exists(gromosXX_build_path+"/bin") or recompile or recompile_from_scratch):
                _make_compile(build_dir=gromosXX_build_path, nCore=nCore)

        #gromosPP
        if(_do_gromosPP):
            print(spacer+"\n BUILD GROMOSPP \n"+spacer)
            gromosPP_path = root_dir+"/gromosPlsPls/gromos++"
            gromosPP_build_path = gromosPP_path+"/build"
            
            if(not os.path.exists(gromosPP_build_path) and recompile_from_scratch):
                _configure_gromosPP_autotools(build_dir=gromosPP_build_path,binary_dir=binary_dir, with_omp=gromosPP_with_omp, with_debug=with_debug, silent=silent)

            if(not os.path.exists(gromosPP_build_path+"/bin") or recompile or recompile_from_scratch):
                _make_compile(build_dir=gromosPP_build_path, nCore=nCore)


def _configure_gromosPP_autotools(build_dir:str, binary_dir:str=None, with_omp:bool=False, with_debug: bool = False, silent:bool=False):
    configure_options = {}

    root_dir = os.path.dirname(build_dir)
    os.chdir(root_dir)
    print(os.getcwd())

    # Initial config
    print(spacer2+"\t\t> INIT \n"+spacer2)
    bash.execute("./Config.sh")
    print()

    # Configure
    bash.make_folder(build_dir)
    os.chdir(build_dir)

    print(spacer2+"\t\t> CONFIGURE \n"+spacer2)
    log_file = build_dir+"/configure.log"
    print("log_file: ", log_file)
    
    options={}
    if(binary_dir is not None):
        options.update({"--bindir": binary_dir})
        
    flags = []
    if(with_omp):
        flags.append("--enable-openmp")
    if(with_debug):
        flags.append("--enable-debug")
    if(silent):
        flags.append("--enable-silent-rule")
        
    cmd = "../configure "+" ".join([key+"="+val for key, val in options.items()])+" "+" ".join([flag for flag in flags])
    print("command: ", cmd)
    
    bash.execute(cmd, catch_STD=log_file)
    print()


def _configure_gromosXX_autotools(build_dir:str, binary_dir:str=None, with_cuda_dir:str=None, with_omp:bool=False, with_mpi:bool=False, with_cuda:bool=False, with_debug:bool=False):
    configure_options = {}

    root_dir = os.path.dirname(build_dir)
    os.chdir(root_dir)
    print(os.getcwd())

    # Initial config
    print(spacer2+"\t\t> INIT \n"+spacer2)
    bash.execute("./Config.sh")
    print()

    # Configure
    bash.make_folder(build_dir)
    os.chdir(build_dir)

    print(spacer2+"\t\t> CONFIGURE \n"+spacer2)
    log_file = build_dir+"/configure.log"
    print("log_file: ", log_file)
    
    options={}
    if(binary_dir is not None):
        options.update({"--bindir": binary_dir})
    if(with_cuda_dir is not None):
        options.update({"--with-cuda": with_cuda_dir})
        
                
    flags = []
    if(with_omp and with_mpi):
        raise ValueError("Can not use with_omp and with_mpi at the same time!")
    if(with_omp):
        flags.append("--enable-openmp")
    if(with_mpi):
        flags.append("--enable-mpi")
    if(with_cuda):
        raise NotImplementedError("sorry the with_cuda was not tested or implemented yet")
        flags.append("--with-cuda") #not the correct flag!
    if(with_debug):
        flags.append("--enable-debug")   
        
    cmd = "../configure "+" ".join([key+"="+val for key, val in options.items()])+" "+" ".join([flag for flag in flags])
    print("command: ", cmd)

    bash.execute(cmd, catch_STD=log_file)
    print()
    
def _make_compile(build_dir:str, nCore:int = 1): 
    # Compile      
    print(spacer2+"\t\t> MAKE \n"+spacer2)
    os.chdir(build_dir)
    
    log_file = build_dir+"/make.log"
    print("log_file: ", log_file)

    cmd = "make -j"+str(nCore)
    print("command: ", cmd)
    
    bash.execute(cmd, catch_STD=log_file)
    print()

    # Create Binaries
    print(spacer2+"\t\t> INSTALL \n"+spacer2)
    
    log_file = build_dir+"/make.log"
    print("log_file: ", log_file)

    cmd = "make -j"+str(nCore)+" install"
    print("command: ", cmd)

    bash.execute(cmd, catch_STD=build_dir+"/makeInstall.log")
    print()


if __name__ == "__main__":
    
    root_dir = os.path.dirname(__file__)    
    _do_gromosXX=False
    _do_gromosPP=True

    make_clean = False
    recompile = False
    recompile_from_scratch = True
    nCores=1


    
    install_gromos(root_dir=root_dir, 
                   _do_gromosXX = _do_gromosXX,
                   _do_gromosPP=_do_gromosPP, 
                   make_clean = True
                   recompile=recompile, 
                   recompile_from_scratch=recompile_from_scratch,
                   nCore=nCores)