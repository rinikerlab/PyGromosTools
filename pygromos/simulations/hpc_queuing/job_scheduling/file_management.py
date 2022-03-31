"""
    This module is doing all the post simulation file juggeling needed for gromos. CURRENTLY OLD DON"T USE
"""
import os
import glob
import tempfile
import warnings
import pandas as pd
import multiprocessing as mult
from subprocess import SubprocessError

from pygromos.files.simulation_parameters import imd
from pygromos.files.otherfiles import repdat
from pygromos.files.trajectory import tre
from pygromos.gromos import gromosPP
from pygromos.utils import bash
from pygromos.utils.typing import List, Dict, Union, Tuple, Cnf_Type

"""
    PARALLEL WORKER - These functions are required for parallelized code Execution
"""


def _thread_worker_cat_trc(
    job: int,
    replicaID_range: List[int],
    trc_files: Dict[int, List[str]],
    out_prefix: str,
    topology_path: str,
    out_trcs: dict,
    dt: float,
    time: float = 0,
    verbose: bool = False,
    boundary_conditions: str = "r cog",
    include_all: bool = False,
):
    """_thread_worker_cat_trc
        This thread worker_scripts concatenates all .trc files of one replica into one file.

    Parameters
    ----------
    job :   rank of this thread
    replicaID_range :   x_range - list of all
    trc_files :     Dict[int, List[str]]
        Dictionary containing all replicas, with list of all trc files concerning one trc.
    out_prefix : str
        output prefix
    verbose : bool
        verbosity?

    Returns
    -------
    None
    """

    gromPP = gromosPP.GromosPP()
    start_dir = os.getcwd()
    if verbose:
        print("JOB " + str(job) + ": range " + str(list(replicaID_range)))
    for replicaID in replicaID_range:
        out_path = out_prefix + str(replicaID) + ".trc"
        compress_out_path = out_path + ".gz"

        out_trcs.update({replicaID: compress_out_path})

        if os.path.exists(compress_out_path):  # found perfect compressed trc file:)
            warnings.warn("Skipped generating file as I found: " + compress_out_path)
            if os.path.exists(out_path):
                bash.remove_file(out_path)
            continue
        elif os.path.exists(out_path):  # did not find compressed file. will compress
            warnings.warn("Skipped generating file as I found: " + out_path)
            continue
        else:  # concat files
            if verbose:
                print("JOB " + str(job) + ": " + "write out " + out_path + "\n")
            out_dir = os.path.dirname(out_path)
            tmp_dir = bash.make_folder(out_dir + "/TMP_replica_" + str(replicaID), additional_option="-p")
            os.chdir(tmp_dir)
            if include_all:
                out_path = gromPP.frameout(
                    in_top_path=topology_path,
                    in_coord_path=" ".join(trc_files[replicaID]),
                    periodic_boundary_condition=boundary_conditions,
                    single_file=True,
                    out_file_format="trc",
                    out_file_path=out_path,
                    time=time,
                    dt=dt,
                    include="ALL",
                )
            else:
                out_path = gromPP.frameout(
                    in_top_path=topology_path,
                    in_coord_path=" ".join(trc_files[replicaID]),
                    periodic_boundary_condition=boundary_conditions,
                    single_file=True,
                    out_file_format="trc",
                    out_file_path=out_path,
                    time=time,
                    dt=dt,
                )
            os.chdir(start_dir)
            bash.wait_for_fileSystem(out_path)
            bash.remove_folder(tmp_dir)
            if verbose:
                print("JOB " + str(job) + ": " + "write out " + out_path + "\t DONE\n")

        if verbose:
            print("JOB " + str(job) + ": " + "compress " + compress_out_path + "\n")
        compressed_trc = bash.compress_gzip(out_path, out_path=compress_out_path)

        if verbose:
            print("JOB " + str(job) + ": " + "compress " + compressed_trc + "\t DONE\n")


def _thread_worker_cat_tre(
    job: int,
    replicaID_range: List[int],
    tre_files: Dict[int, List[str]],
    out_prefix: str,
    out_tres: dict,
    verbose: bool = False,
):
    if verbose:
        print("JOB " + str(job) + ": range " + str(list(replicaID_range)))

    for replicaID in replicaID_range:
        use_tre_file_paths, unarchived_files = find_and_unarchive_tar_files(tre_files[replicaID], verbose=verbose)
        if verbose:
            print("FILES: ", use_tre_file_paths)
        if verbose:
            print("Archs:", unarchived_files)

        out_path = out_prefix + str(replicaID) + ".tre"
        compressed_tre = out_path + ".gz"
        if os.path.exists(compressed_tre):
            warnings.warn("Skipped generating .tre.gz file as I found: " + out_path)
        else:
            if os.path.exists(out_path):
                warnings.warn("Skipped generating .tre file as I found: " + out_path + "\n\t Continue Compressing.")
            else:
                tre_file = tre.Tre(use_tre_file_paths[0])
                if verbose:
                    print("JOB " + str(job) + ": parsing " + os.path.basename(use_tre_file_paths[0]))
                if len(use_tre_file_paths) > 1:
                    for tre_file_path in use_tre_file_paths[1:]:
                        if verbose:
                            print("JOB " + str(job) + ": append " + os.path.basename(tre_file_path))
                        tre_file += tre.Tre(tre_file_path)
                if verbose:
                    print("JOB " + str(job) + ": write out " + os.path.basename(out_path))
                tre_file.write(out_path)
                bash.wait_for_fileSystem(out_path)
                if verbose:
                    print("JOB " + str(job) + ": done " + os.path.basename(out_path))
                del tre_file

                if verbose:
                    print("JOB " + str(job) + ":  compress " + os.path.basename(out_path))
                compressed_tre = bash.compress_gzip(out_path, out_path=compressed_tre)
                if verbose:
                    print("JOB " + str(job) + ": " + "compress " + compressed_tre + "\t DONE\n")

        # file cleaning:
        for file_path in use_tre_file_paths:
            bash.compress_gzip(file_path)
        out_tres.update({replicaID: compressed_tre})


def thread_worker_concat_repdat(
    job: int, repdat_file_out_path: str, repdat_file_paths: Union[str, List[str]], verbose: bool = False
) -> str:
    if os.path.exists(repdat_file_out_path):
        warnings.warn("Skipped repdat creation as already existed!: " + repdat_file_out_path)
    else:
        if verbose:
            print("JOB " + str(job) + ": Found repdats: " + str(repdat_file_paths))  # verbose_repdat
        if isinstance(repdat_file_paths, str):
            repdat_file_paths = [repdat_file_paths]

        if verbose:
            print("JOB " + str(job) + ": Concatenate repdats: \tSTART")  # verbose_repdat
        repdat_file = repdat.Repdat(repdat_file_paths.pop(0))  # repdat Class
        for repdat_path in repdat_file_paths:
            if verbose:
                print("JOB " + str(job) + ": concat:\t", repdat_path)
            tmp_repdat_file = repdat.Repdat(repdat_path)
            repdat_file.append(tmp_repdat_file)
            del tmp_repdat_file

        if verbose:
            print("JOB " + str(job) + ": Concatenate repdats: \tDONE")  # verbose_repdat
        if verbose:
            print("JOB " + str(job) + ": Write out repdat: " + str(repdat_file_out_path))  # verbose_repdat
        repdat_file.write(repdat_file_out_path)
        del repdat_file


def _thread_worker_cnfs(
    job: int,
    out_cnfs: List[str],
    in_cnfs: List[Cnf_Type],
    replica_range: List[int],
    out_folder: str,
    verbose: bool = False,
):
    if verbose:
        print("JOB: " + str(job) + " copy to " + out_folder)
    for replicaID in replica_range:
        out_cnfs.update(
            {
                replicaID: bash.copy_file(
                    in_cnfs[replicaID][-1], out_folder + "/" + os.path.basename(in_cnfs[replicaID][-1])
                )
            }
        )


def _thread_worker_conv_trc(
    job: int,
    replica_range: List[int],
    trc_files: List[str],
    in_topology_path: str,
    gromos_path: str,
    out_traj: dict,
    fit_traj_to_mol: int = 1,
    boundary_conditions: str = "r",
    verbose: bool = False,
):
    if verbose:
        print("JOB: " + str(job) + " RANGE\t" + str(replica_range))
    gromPP = gromosPP.GromosPP(gromos_path)
    first = True
    import mdtraj

    start_dir = os.getcwd()
    for replicaID in replica_range:
        trc_path = trc_files[replicaID]
        if first:
            temp = tempfile.mkdtemp(suffix="_job" + str(job), prefix="convert_", dir=os.path.dirname(trc_path))
            os.chdir(temp)
            first = False

        unarchived = False
        use_trc = trc_path

        if verbose:
            print("using trc:", use_trc)
        out_top_path = use_trc.replace(".trc", "_last.pdb").replace(".gz", "")
        out_traj_path = use_trc.replace(".trc", ".dcd").replace(".gz", "")

        out_traj.update({replicaID: {"top": out_top_path, "traj": out_traj_path}})

        if os.path.exists(out_top_path) and os.path.exists(out_traj_path):
            warnings.warn("Found already the traj and its top!:" + out_traj_path)
            continue

        if verbose:
            print(
                "JOB "
                + str(job)
                + ": Converting trc_path to -> "
                + use_trc.replace("trc", "pdb").replace(".tar.gz", "")
            )
        pdb = gromPP.frameout(
            in_top_path=in_topology_path,
            in_coord_path=use_trc,
            periodic_boundary_condition=boundary_conditions,
            gather="cog",
            out_file_format="pdb",
            atomsfit=str(fit_traj_to_mol) + ":a",
            single_file=True,
            out_file_path=use_trc.replace("trc", "pdb").replace("tar", "").replace(".gz", ""),
        )

        if verbose:
            print("JOB " + str(job) + ": loading pdb : " + pdb)
        traj = mdtraj.load_pdb(pdb)
        if verbose:
            print("JOB " + str(job) + ": write out data to " + use_trc.replace(".trc", ".dcd/.pdb"))
        traj.save(out_traj_path)
        traj[0].save(out_top_path)

        del traj
        if verbose:
            print("Done writing out")
        bash.remove_file(pdb)

        if unarchived:
            if verbose:
                print("Clean unarchiving")
            bash.remove_file(use_trc)
    bash.remove_folder(temp)
    os.chdir(start_dir)


def thread_worker_isolate_energies(
    in_en_file_paths: str,
    out_folder: str,
    properties: List[str],
    replicas: List[int],
    in_ene_ana_lib: str,
    gromosPP_path: str,
    out_prefix: str = "",
    tre_prefix: str = "",
    time=None,
    dt=None,
    job: int = -1,
    verbose=True,
) -> List[str]:
    """isolate_properties_from_tre
        This func uses Ene Ana from gromos to isolate potentials from out_tre Files
        in in_folder generated by reeds.

    Parameters
    ----------
    in_en_file_paths : str
        path, in which the input tre_folders are situated.
    out_folder :    str
         output folder, where to write the energy .csvs
    properties :    List[str]
        potentials to isolate from the .out_tre Files
    replicas :  int
        number of replicas, that should be found
    in_ene_ana_lib :    str
         path to the ene_ana lib, encoding the out_tre Files
    gromosPP_path : str
        path to the ene_ana lib, encoding the out_tre Files
    out_prefix :    str, optional
    tre_prefix :    str, optional
    verbose :   bool, optional

    Returns
    -------
    List[str]
        return list of result Files.

    """

    gromos = gromosPP.GromosPP(gromosPP_path)
    result_files = []
    temp = tempfile.mkdtemp(suffix="_job" + str(job), prefix="ene_ana_", dir=out_folder)

    start_dir = os.getcwd()
    os.chdir(temp)

    # get the potentials from each replica.tre*
    if verbose:
        print("JOB" + str(job) + ": working with job: " + str(list(replicas)))
    for replicaID in replicas:
        in_en_file_path = in_en_file_paths[replicaID]

        out_suffix = "energies_s" + str(replicaID)
        out_path = out_folder + "/" + out_prefix + "_" + out_suffix + ".dat"

        if verbose:
            print("CHECKING: " + out_path)
        if not os.path.exists(out_path):
            if verbose:
                print("JOB" + str(job) + ":\t" + str(replicaID))
            if verbose:
                print("JOB" + str(job) + ":", in_en_file_path)
            tmp_out = gromos.ene_ana(
                in_ene_ana_library_path=in_ene_ana_lib,
                in_en_file_paths=in_en_file_path,
                out_energy_folder_path=out_folder,
                out_files_suffix=out_suffix,
                out_files_prefix=out_prefix,
                time=str(time) + " " + str(dt),
                in_properties=properties,
                verbose=verbose,
                single_file=True,
                workdir=True,
            )
            result_files.append(tmp_out)
            bash.remove_file(temp + "/*")  # remove logs if succesfull

    os.chdir(start_dir)
    bash.remove_folder(temp)
    if verbose:
        print("JOB" + str(job) + ": DONE")
    return result_files


def _thread_worker_delete(job: int, file_paths: List[str], verbose: bool = False) -> int:
    for file_path in file_paths:
        if verbose:
            print("JOB" + str(job) + " - rm: " + file_path + "")
        bash.remove_file(file_path)
    return 0


def _thread_worker_compress(job: int, in_file_paths: List[str], verbose: bool = False) -> int:
    for file_path in in_file_paths:
        if verbose:
            print("JOB" + str(job) + " - gz: " + file_path)
        bash.compress_gzip(in_path=file_path, verbose=verbose)
    return 0


"""
    FILE FINDER
"""


def find_and_unarchive_tar_files(trc_files: List[str], verbose: bool = False):
    # archive handling
    archived_files = list(filter(lambda x: (".tar" in x or ".gz" in x or ".tar.gz" in x), trc_files))
    not_archived_files = list(filter(lambda x: not ("tar" in x or ".gz" in x or ".tar.gz" in x), trc_files))
    unarchived_files = []
    if verbose:
        print("archives: ", archived_files)
    if verbose:
        print("narchives: ", not_archived_files)

    # untar files:
    for tared_file in archived_files:
        if len(not_archived_files) == 0 or not any([noAfile in tared_file for noAfile in not_archived_files]):
            try:
                # print("Ungzip ->\t", tared_file)
                out_path = bash.compress_gzip(
                    in_path=tared_file, out_path=tared_file.replace(".tar", "").replace(".gz", ""), extract=True
                )
            except SubprocessError:
                # print("Failed gzip, trying tar")
                out_path = bash.extract_tar(
                    in_path=tared_file,
                    out_path=tared_file.replace(".tar", "").replace(".gz", ""),
                    gunzip_compression=True,
                )

            # fix for stupid taring!    #todo: remove part
            if any(["cluster" == xfile for xfile in os.listdir(os.path.dirname(tared_file))]) and not os.path.exists(
                out_path
            ):
                nfound = True
                for cpath, tdir, files in os.walk(os.path.dirname(tared_file) + "/cluster"):
                    if os.path.basename(tared_file).replace(".tar", "").replace(".gz", "") in files:
                        if verbose:
                            print(
                                "FOUND PATH: ",
                                cpath + "/" + os.path.basename(tared_file).replace(".tar", "").replace(".gz", ""),
                            )
                        wrong_path = cpath + "/" + os.path.basename(tared_file).replace(".tar", "").replace(".gz", "")
                        out_file = bash.move_file(wrong_path, tared_file.replace(".tar", "").replace(".gz", ""))
                        unarchived_files.append(out_file)
                        nfound = False
                        break
                if nfound:
                    raise IOError("could not find untarred file!")
            else:
                unarchived_files.append(out_path)

                # raise Exception("this tar needs special treatment as it is in cluster/yadda/yadda/fun.trc")
        else:
            if verbose:
                print([noAfile for noAfile in not_archived_files if (noAfile in tared_file)])
    new_files = not_archived_files
    new_files.extend(unarchived_files)

    use_tre_file_paths = sorted(new_files, key=lambda x: int(x.split("_")[-1].split(".")[0]))
    return use_tre_file_paths, unarchived_files


def gather_simulation_replica_file_paths(
    in_folder: str,
    replicas: int,
    filePrefix: str = "",
    fileSuffixes: Union[str, List[str]] = [".tre", ".tre.tar.gz"],
    verbose: bool = False,
    finalNumberingSort=False,
) -> Dict[int, List[str]]:
    """gather_replica_file_paths

        Finds all trajectory paths in a simulation folder and sorts them by replica.


    Parameters
    ----------
    in_folder : str
        folder, containing the files
    replicas :  int
        Number of replicas
    filePrefix :    str, optional
        str prefix the desired files
    fileSuffixes :    str, optional
        str suffix of the desired files
    verbose :   bool
        toggle verbosity

    Returns
    -------

    """

    if isinstance(fileSuffixes, str):
        fileSuffixes = [fileSuffixes]

    # browse folders
    # init statewise dic
    files = {}
    for replica in range(1, replicas + 1):
        files.update({replica: []})

    if verbose:
        print("SEARCH PATTERN: " + filePrefix + " + * +" + str(fileSuffixes))

    for dirname, dirnames, filenames in os.walk(in_folder):
        if str(dirname[-1]).isdigit() and os.path.basename(dirname).startswith("eq"):
            continue
        # check actual in_dir for fle pattern
        tmp_files = [
            file for file in filenames if (filePrefix in file and any([suffix in file for suffix in fileSuffixes]))
        ]

        if len(tmp_files) != 0:
            for x in range(1, replicas + 1):
                grom_file_prefix = sorted(
                    [
                        dirname + "/" + file
                        for file in tmp_files
                        if (any(["_" + str(x) + suffix in file for suffix in fileSuffixes]))
                    ]
                )
                files[x] += grom_file_prefix
        if verbose:
            print("walking to in_dir: ", os.path.basename(dirname), "found: ", len(tmp_files))

    if not finalNumberingSort:
        # final_cleanup
        for x in files:
            files[x].sort(key=lambda x: int(x.split("_")[-2]))

    if verbose:
        print("\nfoundFiles:\n")
        for x in sorted(files):
            print("\n" + str(x))
            print("\t" + "\t".join([y + "\n" for y in files[x]]))

    if len(files[1]) == 0:
        raise ValueError("could not find any file with the prefix: " + filePrefix + " in folder : \n" + in_folder)

    return files


def gather_simulation_file_paths(
    in_folder: str,
    filePrefix: str = "",
    fileSuffixes: Union[str, List[str]] = [".tre", ".tre.tar.gz"],
    files_per_folder: int = 1,
    verbose: bool = False,
) -> List[str]:
    files = []
    if isinstance(fileSuffixes, str):
        fileSuffixes = [fileSuffixes]

    if verbose:
        print("SEARCH PATTERN: " + filePrefix + " + * +" + str(fileSuffixes))

    for dirname, dirnames, filenames in os.walk(in_folder):
        if str(dirname[-1]).isdigit() and os.path.basename(dirname).startswith("eq"):
            continue
        # check actual in_dir for fle pattern
        tmp_files = [
            file for file in filenames if (filePrefix in file and any([suffix in file for suffix in fileSuffixes]))
        ]

        if len(tmp_files) == files_per_folder:
            files.extend(list(map(lambda x: dirname + "/" + x, tmp_files)))

        if verbose:
            print("walking to in_dir: ", os.path.basename(dirname), "found: ", len(tmp_files))

    try:
        keys = [[int(y) for y in x.split("_") if (y.isdecimal())][-1] for x in files]
        sorted_files = list(map(lambda y: y[1], sorted(zip(keys, files), key=lambda x: x[0])))
    except SubprocessError:
        warnings.warn("Files are not all enumerated! no file sorting.")
        sorted_files = files

    if verbose:
        print("\nfoundFiles:\n")
        print("\t" + "\n\t".join(sorted_files))

    if len(sorted_files) == 0:
        raise ValueError("could not find any file with the prefix: " + filePrefix + " in folder : \n" + in_folder)

    return sorted_files


"""
    ENERGY FILE FUNCTIONS
"""


def find_header(path: str) -> int:
    comment_lines = -1
    with open(path, "r") as file:
        for line in file.readlines():
            if line.strip().startswith("#"):
                comment_lines += 1
                continue
            else:
                break
        if comment_lines == -1:
            return 0
    return comment_lines


def parse_csv_energy_trajectory(in_ene_traj_path: str, verbose: bool = False) -> pd.DataFrame:
    """
        parse_one ene_ana csv

    Parameters
    ----------
    in_ene_traj_path : str
        path to input file

    verbose :   bool
        loud?

    Returns
    -------
    pd.DataFrame
        return a pandas data frame containing all energies
    """
    if verbose:
        print("deal with: ", in_ene_traj_path)
    ene_traj = pd.read_csv(in_ene_traj_path, header=find_header(in_ene_traj_path), delim_whitespace=True)
    ene_traj.columns = [x.replace("#", "").strip() for x in ene_traj.columns]
    setattr(ene_traj, "in_path", in_ene_traj_path)
    return ene_traj


def parse_csv_energy_trajectories(in_folder: str, ene_trajs_prefix: str, verbose: bool = False) -> List[pd.DataFrame]:
    """
        searches a directory and loads energy eds csvs as pandas dataframes.

    Parameters
    ----------
    in_folder : str
        folder with energy_traj - csvs
    ene_trajs_prefix :  str
        prefix name
    verbose :   bool
        loud?

    Returns
    -------
    List[pd.DataFrame]
        return a list with pandas data frames containing all energy infos.
    """

    if verbose:
        print("SEARCH: ", in_folder + "/" + ene_trajs_prefix + "*.dat")
    in_ene_traj_paths = sorted(
        glob.glob(in_folder + "/" + ene_trajs_prefix + "*.dat"),
        key=lambda x: int(x.split("_")[-1].split(".")[0].replace("s", "")),
    )
    ene_trajs: List[pd.DataFrame] = []
    if verbose:
        print("FOUND: ", "\n".join(in_ene_traj_paths))

    for in_ene_traj_path in in_ene_traj_paths:
        ene_traj = parse_csv_energy_trajectory(in_ene_traj_path, verbose=verbose)
        if verbose:
            print("csv columns: \t", ene_traj.columns)
        # note: the previous version created problems for filenames which contained periods
        # setattr(ene_traj, "s", ((ene_traj.in_path.split("."))[0]).split("_")[-1])
        # setattr(ene_traj, "replicaID", int(((ene_traj.in_path.split("."))[0]).split("_")[-1].replace("s", "")))
        setattr(ene_traj, "s", ((ene_traj.in_path.split("."))[-2]).split("_")[-1])
        setattr(ene_traj, "replicaID", int(((ene_traj.in_path.split("."))[-2]).split("_")[-1].replace("s", "")))
        ene_trajs.append(ene_traj)

    if len(ene_trajs) == 0:
        raise ValueError("could not find any energy_trajectory in: ", in_folder + "/" + ene_trajs_prefix + "*.dat")

    ene_trajs = list(sorted(ene_trajs, key=lambda x: int(x.s.replace("s", ""))))
    return ene_trajs


"""
    concatenation wrapper
"""


def project_concatenation(
    in_folder: str,
    in_topology_path: str,
    in_imd: str,
    num_replicas: int,
    control_dict: Dict[str, bool],
    out_folder: str,
    in_ene_ana_lib_path: str,
    out_file_prefix: str = "test",
    fit_traj_to_mol: int = 1,
    starting_time: float = 0,
    include_water_in_trc=True,
    additional_properties: Union[Tuple[str], List[str]] = ("solvtemp2", "totdisres"),
    n_processes: int = 1,
    gromosPP_bin_dir: str = None,
    verbose: bool = False,
    nofinal=False,
    boundary_conditions: str = "r cog",
) -> dict:
    if verbose:
        print("reading imd file: " + in_imd)

    imd_file = imd.Imd(in_imd)
    dt = float(imd_file.STEP.DT)
    dt_trc = int(imd_file.WRITETRAJ.NTWX) * dt
    dt_tre = int(imd_file.WRITETRAJ.NTWE) * dt

    tmp_dir = out_folder + "/tmp_file"
    if os.path.isdir(tmp_dir):
        bash.make_folder(tmp_dir)

    out_cnfs = out_tres = out_trcs = out_dcd = out_repdat = None
    p_conv = p_cnf = p_repdat = p_trc = p_tre = p_ene_ana = False
    submitted_trc_job = submitted_tre_job = False

    if n_processes > 1:
        p = mult.Pool(n_processes)
        manager = mult.Manager()

    if control_dict["cp_cnf"]:
        if verbose:
            print("\tStart cnfs")
        # find all cnf files in this project
        sim_dir_cnfs = gather_simulation_replica_file_paths(
            in_folder, num_replicas, filePrefix="", fileSuffixes=".cnf", verbose=verbose, finalNumberingSort=nofinal
        )

        # do parallel
        if n_processes > 1:
            out_cnfs = manager.dict()
            distribute = [
                (n, out_cnfs, sim_dir_cnfs, range(n, len(sim_dir_cnfs) + 1, n_processes), out_folder, verbose)
                for n in range(1, n_processes + 1)
            ]
            # _async
            p_cnf = p.starmap(_thread_worker_cnfs, distribute)
        else:
            out_cnfs = {}
            _thread_worker_cnfs(
                job=-1,
                out_cnfs=out_cnfs,
                in_cnfs=sim_dir_cnfs,
                replica_range=list(sim_dir_cnfs.keys()),
                out_folder=out_folder,
                verbose=verbose,
            )
        if verbose:
            print("Out cnfs: ", out_cnfs)

    if control_dict["cat_trc"]:
        print("\tStart Trc Cat")
        # find all trc files in this project
        trc_files = gather_simulation_replica_file_paths(
            in_folder,
            num_replicas,
            filePrefix="",
            fileSuffixes=[".trc", ".trc.gz", ".trc.tar.gz"],
            verbose=verbose,
            finalNumberingSort=nofinal,
        )

        out_prefix = out_folder + "/" + out_file_prefix + "_"

        # concat all files to a single .trc
        if n_processes > 1:  # parallel
            submitted_trc_job = True
            if verbose:
                print("going parallel: n_processes - " + str(n_processes))
            out_trcs = manager.dict()
            distributed_jobs = [
                (
                    n,
                    range(n, len(trc_files) + 1, n_processes),
                    trc_files,
                    out_prefix,
                    in_topology_path,
                    out_trcs,
                    dt_trc,
                    starting_time,
                    verbose,
                    include_water_in_trc,
                )
                for n in range(1, n_processes + 1)
            ]
            # _async
            p_trc = p.starmap(_thread_worker_cat_trc, distributed_jobs)
            p.close()
            p.join()

        else:
            out_trcs = {}
            _thread_worker_cat_trc(
                job=-1,
                topology_path=in_topology_path,
                replicaID_range=list(trc_files.keys()),
                trc_files=trc_files,
                out_prefix=out_prefix,
                dt=dt_trc,
                time=starting_time,
                out_trcs=out_trcs,
                verbose=verbose,
                boundary_conditions=boundary_conditions,
                include_all=include_water_in_trc,
            )

    if control_dict["cat_tre"]:
        print("\tStart Tre Cat")

        # find all trc files in this project
        tre_files = gather_simulation_replica_file_paths(
            in_folder,
            num_replicas,
            filePrefix="",
            fileSuffixes=[".tre", ".tre.tar.gz"],
            verbose=verbose,
            finalNumberingSort=nofinal,
        )

        out_prefix = out_folder + "/" + out_file_prefix + "_"
        # concat all files to a single .trc
        if n_processes > 1:
            if verbose:
                print("going parallel: n_processes - " + str(n_processes), " for ", len(tre_files))
            submitted_tre_job = True
            out_tres = manager.dict()
            distributed_jobs = [
                (n, range(n, len(tre_files) + 1, n_processes), tre_files, out_prefix, out_tres, verbose)
                for n in range(1, n_processes + 1)
            ]
            p = mult.Pool(n_processes)
            p_tre = p.starmap(_thread_worker_cat_tre, distributed_jobs)
            p.close()
            p.join()
        else:
            out_tres = {}
            _thread_worker_cat_tre(
                job=-1,
                replicaID_range=tre_files.keys(),
                tre_files=tre_files,
                out_prefix=out_prefix,
                out_tres=out_tres,
                verbose=verbose,
            )

    if control_dict["ene_ana"]:
        print("\tStart ene ana")

        # wait for async job creating the trcs.
        if submitted_tre_job:
            p_tre.wait()

        # gather potentials
        properties = list(additional_properties)
        # find all trc files in this project
        tre_files = gather_simulation_replica_file_paths(
            in_folder,
            num_replicas,
            filePrefix="",
            fileSuffixes=[".tre", ".tre.gz"],
            verbose=verbose,
            finalNumberingSort=nofinal,
        )  # ".tre.tar.gz"
        # isolate potentials
        if verbose:
            print("Isolate ene_ana:")
        if n_processes > 1:
            p = mult.Pool(n_processes)
            distribute_jobs = [
                (
                    out_folder,
                    out_folder,
                    properties,
                    list(tre_files.keys())[n::n_processes],
                    in_ene_ana_lib_path,
                    gromosPP_bin_dir,
                    out_file_prefix,
                    "",
                    dt_tre,
                    n,
                    verbose,
                )
                for n in range(n_processes)
            ]
            p_ene_ana = p.starmap_async(thread_worker_isolate_energies, distribute_jobs)
        else:
            thread_worker_isolate_energies(
                in_en_file_paths=tre_files,
                out_folder=out_folder,
                properties=properties,
                out_prefix=out_file_prefix,
                in_ene_ana_lib=in_ene_ana_lib_path,
                gromosPP_path=gromosPP_bin_dir,
                dt=dt_tre,
                replicas=list(tre_files.keys()),
                verbose=verbose,
            )

    if control_dict["convert_trcs"]:
        print("\tStart Trc Conversion")
        # wait for async job creating the trcs.
        if submitted_trc_job:
            p_trc.wait()

        # get files:
        final_trc_files = list(
            sorted(glob.glob(out_folder + "/*.trc*"), key=lambda x: int(x.split("_")[-1].split(".")[0]))
        )

        if n_processes > 1:
            out_dcd = manager.dict()
            distributed_jobs = [
                (
                    n,
                    range(n, num_replicas, n_processes),
                    final_trc_files,
                    in_topology_path,
                    gromosPP_bin_dir,
                    out_dcd,
                    fit_traj_to_mol,
                    verbose,
                )
                for n in range(n_processes)
            ]
            p_conv = p.starmap_async(_thread_worker_conv_trc, distributed_jobs)
        else:
            out_dcd = {}
            _thread_worker_conv_trc(
                job=-1,
                replica_range=range(num_replicas),
                trc_files=final_trc_files,
                in_topology_path=in_topology_path,
                gromos_path=gromosPP_bin_dir,
                out_traj=out_dcd,
                fit_traj_to_mol=1,
                verbose=verbose,
                boundary_conditions=boundary_conditions,
            )

    if n_processes > 1:
        # wait for the jobs to finish
        if (
            (not p_conv or p_conv.wait())
            and (not p_cnf or p_cnf.wait())
            and (not p_repdat or p_repdat.wait())
            and (not p_trc or p_trc.wait())
            and (not p_tre or p_tre.wait())
            and (not p_ene_ana or p_ene_ana.wait())
        ):
            raise ChildProcessError("A process failed! ")

        p.close()
        p.join()

        out_dict = {
            "out_folder": out_folder,
            "cnfs": dict(out_cnfs),
            "repdat": out_repdat,
            "tres": dict(out_tres),
            "trcs": dict(out_trcs),
            "dcds": dict(out_dcd),
        }
        manager.shutdown()

    else:
        out_dict = {
            "out_folder": out_folder,
            "cnfs": out_cnfs,
            "repdat": out_repdat,
            "tres": out_tres,
            "trcs": out_trcs,
            "dcds": out_dcd,
        }
    if verbose:
        print("all jobs finished")
    return out_dict


"""
    COMPRESS FUNCTION
"""


def compress_files(in_paths: List[str], n_processes: int = 1) -> List[str]:
    """compress a list of files

    Parameters
    ----------
    in_paths :  List[str]

    n_processes: int
        how many processes can be used in parallel?

    Returns
    -------
    List[str]
        outpaths
    """

    if type(in_paths) == str:
        in_paths = [in_paths]
    out_paths = []

    # check:
    for path in in_paths:
        if os.path.exists(path):
            archive_path = path + ".gz"
            out_paths.append(archive_path)
        else:
            warnings.warn("File Path: " + path + " was not found!")

    # do:
    print("Gen Gzips:")
    if n_processes == 1:
        for path in in_paths:
            bash.compress_gzip(in_path=path, out_path=path + ".gz")
    else:  # do parallel
        p = mult.Pool(n_processes)
        distribute = [(job, in_paths[job::n_processes], True) for job in range(n_processes)]
        p.starmap(_thread_worker_compress, distribute)
        p.close()
        p.join()
    return out_paths
