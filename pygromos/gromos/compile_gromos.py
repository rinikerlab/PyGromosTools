import os
from datetime import datetime
from pygromos.utils import bash
from pygromos.utils.utils import spacer, spacer2


def default_install(_timing_dict: dict):
    root_dir = os.path.dirname(__file__)
    nCores = 1
    install_gromos(
        root_dir=root_dir,
        do_clean=True,
        gromosXX_with_mpi=False,
        gromosPP_with_omp=False,
        recompile_from_scratch=True,
        nCore=nCores,
        _timing_dict=_timing_dict,
    )


def install_gromos(
    root_dir: str,
    gromosXX_with_mpi: bool = False,
    gromosXX_with_omp: bool = False,
    gromosXX_with_cuda: str = None,
    gromosPP_with_omp=False,
    do_compile: bool = True,
    do_clean: bool = False,
    recompile: bool = False,
    recompile_from_scratch: bool = False,
    _do_gromosPP: bool = True,
    _do_gromosXX: bool = True,
    with_debug: bool = False,
    nCore: int = 1,
    _timing_dict: dict = None,
):
    """

        Install the gromos simulation packages. As a helper: to get the dependencies right, install and activate the provided pygromos environments!

    Parameters
    ----------
    root_dir : str
        this dir contains the gromosXX and the gromosPP git repositories
    gromosXX_with_cuda_dir : str, optional
        use the following cuda path and activate cuda support, by default None
    gromosXX_with_omp : bool, optional
        should gromosXX be compiled with omp - Warning dont combine wiht mpi!, by default False
    gromosXX_with_mpi : bool, optional
        should gromosXX be compiled wiht mpi - Warning dont combine with omp!, by default False
    gromosPP_with_omp : bool, optional
        should gromosPP be compiled with omp, by default False
    do_compile : bool, optional
        compile the programms, by default True
    do_clean : bool, optional
        clean up and remove the programs (can be used together with do_compile to get a clean start), by default False
    recompile : bool, optional
        recompile the programs with make, by default False
    recompile_from_scratch : bool, optional
        recompile with configure and make, by default False
    _do_gromosPP : bool, optional
        install the gromosPP program, by default True
    _do_gromosXX : bool, optional
        install the gromosXX program, by default True
    with_debug : bool, optional
        set gromos debug flag, by default False
    nCore : int, optional
        how many cores should be used to compile?, by default 1
    """
    if do_clean:
        if _do_gromosPP:
            print(spacer + "\n CLEAN GROMOSPP \n" + spacer)
            gromosPP_path = root_dir + "/gromosPlsPls/gromos++"
            gromosPP_build_path = gromosPP_path + "/build"

            if os.path.exists(gromosPP_build_path):
                if _timing_dict is not None:
                    _timing_dict["gromosPP_clean_start"] = datetime.now()
                os.chdir(gromosPP_build_path)
                try:
                    cmd = "make -j" + str(nCore) + " clean"
                    print("command: ", cmd)
                    bash.execute(cmd)
                except Exception:
                    pass
                bash.remove_file(gromosPP_build_path, recursive=True)
                if _timing_dict is not None:
                    _timing_dict["gromosPP_clean_end"] = datetime.now()

        if _do_gromosXX:
            print(spacer + "\n CLEAN GROMOSXX \n" + spacer)
            gromosXX_path = root_dir + "/gromosXX/gromosXX"
            gromosXX_build_path = gromosXX_path + "/build"

            if os.path.exists(gromosXX_build_path):
                if _timing_dict is not None:
                    _timing_dict["gromosXX_clean_start"] = datetime.now()

                os.chdir(gromosXX_build_path)
                try:
                    cmd = "make -j" + str(nCore) + " clean"
                    print("command: ", cmd)
                    bash.execute(cmd)
                except Exception:
                    pass
                bash.remove_file(gromosXX_build_path, recursive=True)
                if _timing_dict is not None:
                    _timing_dict["gromosXX_clean_end"] = datetime.now()

    if do_compile:
        binary_dir = root_dir + "/bin"
        if _do_gromosXX:
            print(spacer + "\n BUILD GROMOSXX \n" + spacer)
            gromosXX_path = root_dir + "/gromosXX/gromosXX"
            gromosXX_build_path = gromosXX_path + "/build"

            if not os.path.exists(gromosXX_build_path) or recompile_from_scratch:
                if _timing_dict is not None:
                    _timing_dict["gromosXX_conf_start"] = datetime.now()
                _configure_gromosXX_autotools(
                    build_dir=gromosXX_build_path,
                    binary_dir=binary_dir,
                    with_mpi=gromosXX_with_mpi,
                    with_omp=gromosXX_with_omp,
                    with_cuda_dir=gromosXX_with_cuda,
                )
                if _timing_dict is not None:
                    _timing_dict["gromosXX_conf_end"] = datetime.now()

            if not os.path.exists(gromosXX_build_path + "/bin") or recompile or recompile_from_scratch:
                if _timing_dict is not None:
                    _timing_dict["gromosXX_make_start"] = datetime.now()
                _make_compile(build_dir=gromosXX_build_path, nCore=nCore)
                if _timing_dict is not None:
                    _timing_dict["gromosXX_make_end"] = datetime.now()

        if _do_gromosPP:
            print(spacer + "\n BUILD GROMOSPP \n" + spacer)
            gromosPP_path = root_dir + "/gromosPlsPls/gromos++"
            gromosPP_build_path = gromosPP_path + "/build"

            if not os.path.exists(gromosPP_build_path) or recompile_from_scratch:
                if _timing_dict is not None:
                    _timing_dict["gromosPP_conf_start"] = datetime.now()
                _configure_gromosPP_autotools(
                    build_dir=gromosPP_build_path,
                    binary_dir=binary_dir,
                    with_omp=gromosPP_with_omp,
                    with_debug=with_debug,
                )
                if _timing_dict is not None:
                    _timing_dict["gromosPP_conf_end"] = datetime.now()

            if not os.path.exists(gromosPP_build_path + "/bin") or recompile or recompile_from_scratch:
                if _timing_dict is not None:
                    _timing_dict["gromosPP_make_start"] = datetime.now()
                _make_compile(build_dir=gromosPP_build_path, nCore=nCore)
                if _timing_dict is not None:
                    _timing_dict["gromosPP_make_end"] = datetime.now()


def _configure_gromosPP_autotools(
    build_dir: str, binary_dir: str = None, with_omp: bool = False, with_debug: bool = False
):
    """
        Setting the configurations for the compiling gromosPP process. (uses autotools)


    Parameters
    ----------
    build_dir : str
        directory, that should be used for building
    binary_dir : str, optional
        directory in which the binaries should be written to, by default None
    with_omp : bool, optional
        should gromosPP be compiled with omp, by default False
    with_debug : bool, optional
        set gromos debug flag, by default False

    """

    root_dir = os.path.dirname(build_dir)
    os.chdir(root_dir)
    print(os.getcwd())

    # Initial config
    print(spacer2 + "\t\t> INIT \n" + spacer2)
    bash.execute("./Config.sh")
    print()

    # Configure
    bash.make_folder(build_dir)
    os.chdir(build_dir)

    print(spacer2 + "\t\t> CONFIGURE \n" + spacer2)
    log_file = build_dir + "/configure.log"
    print("log_file: ", log_file)

    options = {}
    if binary_dir is not None:
        options.update({"--bindir": binary_dir})

    flags = []
    if with_omp:
        flags.append("--enable-openmp")
    if with_debug:
        flags.append("--enable-debug")

    cmd = (
        "../configure "
        + " ".join([key + "=" + val for key, val in options.items()])
        + " "
        + " ".join([flag for flag in flags])
    )
    print("command: ", cmd)

    bash.execute(cmd, catch_STD=log_file)
    print()


def _configure_gromosXX_autotools(
    build_dir: str,
    binary_dir: str = None,
    with_cuda_dir: str = None,
    with_omp: bool = False,
    with_mpi: bool = False,
    with_debug: bool = False,
):
    """
        Setting the configurations for the compiling gromosXX process. (uses autotools)

    Parameters
    ----------
    build_dir : str
        directory, that should be used for building
    binary_dir : str, optional
        directory in which the binaries should be written to, by default None
    with_cuda_dir : str, optional
        use the following cuda path and activate cuda support, by default None
    with_omp : bool, optional
        should gromosXX be compiled with omp - Warning dont combine wiht mpi!, by default False
    with_mpi : bool, optional
        should gromosXX be compiled wiht mpi - Warning dont combine with omp!, by default False
    with_debug : bool, optional
        set gromos debug flag, by default False

    Raises
    ------
    ValueError
        if wrong value was passed
    """

    root_dir = os.path.dirname(build_dir)
    os.chdir(root_dir)
    print(os.getcwd())

    # Initial config
    print(spacer2 + "\t\t> INIT \n" + spacer2)
    bash.execute("./Config.sh")
    print()

    # Configure
    bash.make_folder(build_dir)
    os.chdir(build_dir)

    print(spacer2 + "\t\t> CONFIGURE \n" + spacer2)
    log_file = build_dir + "/configure.log"
    print("log_file: ", log_file)

    options = {}
    if binary_dir is not None:
        options.update({"--bindir": binary_dir})
    if with_cuda_dir is not None:
        options.update({"--with-cuda": with_cuda_dir})

    flags = []
    if with_omp and with_mpi:
        raise ValueError("Can not use with_omp and with_mpi at the same time!")
    if with_omp:
        flags.append("--enable-openmp")
    if with_mpi:
        flags.append("--enable-mpi")
    if with_debug:
        flags.append("--enable-debug")

    cmd = (
        "../configure "
        + " ".join([key + "=" + val for key, val in options.items()])
        + " "
        + " ".join([flag for flag in flags])
    )
    print("command: ", cmd)

    bash.execute(cmd, catch_STD=log_file)
    print()


def _make_compile(build_dir: str, nCore: int = 1):
    """
        This function triggers make and make install in the build_dir.

    Parameters
    ----------
    build_dir : str
        directory prepared for make.
    nCore : int, optional
        how many cores for each make command?, by default 1
    """
    # Compile
    print(spacer2 + "\t\t> MAKE \n" + spacer2)
    os.chdir(build_dir)

    log_file = build_dir + "/make.log"
    print("log_file: ", log_file)

    cmd = "make -j" + str(nCore)
    print("command: ", cmd)

    bash.execute(cmd, catch_STD=log_file)
    print()

    # Create Binaries
    print(spacer2 + "\t\t> INSTALL \n" + spacer2)

    log_file = build_dir + "/makeInstall.log"
    print("log_file: ", log_file)

    cmd = "make -j" + str(nCore) + " install"
    print("command: ", cmd)

    bash.execute(cmd, catch_STD=log_file)
    print()


if __name__ == "__main__":

    root_dir = os.path.dirname(__file__)
    _do_gromosXX = True
    _do_gromosPP = False

    make_clean = True
    recompile = False
    recompile_from_scratch = False
    nCores = 3

    """
    install_gromos(root_dir=root_dir,
                   _do_gromosXX = _do_gromosXX,
                   _do_gromosPP=_do_gromosPP,
                   do_clean = make_clean,
                   recompile=recompile,
                   recompile_from_scratch=recompile_from_scratch,
                   nCore=nCores)
    """
    default_install()
