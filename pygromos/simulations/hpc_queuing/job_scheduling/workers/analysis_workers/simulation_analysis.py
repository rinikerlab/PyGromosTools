#!/usr/bin/env python3
import os, sys

package_path = os.path.abspath(__file__ + "/../../../../../..")
# print(package_path)
sys.path.append(package_path)

import warnings
import shutil
import math
from typing import Dict, Union, List
from collections import OrderedDict
from pygromos.gromos import gromosPP
from pygromos.files.coord import cnf
from pygromos.files.trajectory import trc, tre, trg
from pygromos.utils import bash

template_control_dict = OrderedDict(
    {
        "concat": {
            "do": True,
            "sub": {
                "cp_cnf": True,
                "repair_NaN": True,
                "cp_omd": True,
                "cat_trc": True,
                "cat_tre": True,
                "cat_trg": False,
                "convert_trcs": False,
            },
        },
        "simulation_folder": {"do": True, "sub": {"tar": True, "remove": False}},
    }
)


def do(
    in_simulation_dir: str,
    out_analysis_dir: str,
    sim_prefix: str,  # in_system:Gromos_System,
    gromosPP_bin_dir: str = None,
    n_processes: int = 1,
    control_dict: dict = None,
    verbose: bool = True,
):
    """

    Parameters
    ----------
    in_simulation_dir
    in_system
    out_analysis_dir
    gromosPP_bin_dir
    n_processes
    control_dict
    verbose

    Returns
    -------

    """
    data = {}

    gromos = gromosPP.GromosPP(gromosPP_bin_dir)

    if not os.path.exists(out_analysis_dir) and not os.path.isdir(out_analysis_dir):
        bash.make_folder(out_analysis_dir)
    if not isinstance(control_dict, dict):
        warnings.warn("Recieved a non dict control-dict! Using the template. \n\tGot:" + str(control_dict))
        control_dict = template_control_dict

    #   ORGANIZE FILES
    out_data_dir = out_analysis_dir + "/data"
    if not os.path.isdir(out_data_dir):
        os.mkdir(out_data_dir)

    if control_dict["concat"]["do"]:
        project_concatenation(
            in_folder=in_simulation_dir,
            in_prefix=sim_prefix,
            out_folder=out_data_dir,
            control_dict=control_dict["concat"]["sub"],
            verbose=verbose,
        )
        # in_simSystem=in_system,
    # Other Analysis parts

    print(control_dict)
    if control_dict["simulation_folder"]["do"]:
        sub_dict = control_dict["simulation_folder"]["sub"]
        if os.path.exists(in_simulation_dir):
            if sub_dict["tar"]:
                out_tar_dir = bash.compress_tar(in_path=in_simulation_dir)
                bash.remove_file(in_simulation_dir, additional_options="-r")

            if sub_dict["remove"]:
                bash.remove_file(in_simulation_dir, additional_options="-r")
        else:
            warnings.warn("Simulation dir was not present. Skipped Compression.\nGiven Path: " + in_simulation_dir)

    return 0


def project_concatenation(
    in_folder: str,
    out_folder: str,
    in_prefix: str,  # in_simSystem: Gromos_System,
    control_dict: Dict[str, bool],
    verbose: bool = False,
) -> dict:
    # in_simSystem.work_folder = out_folder
    out_prefix = out_folder + "/" + in_prefix  # in_simSystem.name
    out_cnf = None
    if control_dict["cp_cnf"]:

        out_cnf = out_prefix + ".cnf"
        if os.path.exists(out_prefix + ".cnf"):
            warnings.warn(
                "Skipping concatenation of tres, as there already exists one!\n "
                "If you want to generate a new trc. Please delete the old one!\n\t Path:" + out_prefix + ".cnf"
            )
        else:
            if verbose:
                print("\tStart cnfs")
            # find all cnf files in this project
            sim_dir_cnfs = gather_simulation_step_file_paths(
                in_folder, filePrefix=in_prefix, fileSuffixes=".cnf", verbose=verbose  # in_simSystem.name,
            )

            import glob

            print(in_folder + "/" + in_prefix + "/*.cnf")
            cnf_path = glob.glob(in_folder + "/" + in_prefix + "*/*.cnf")[-1]
            sim_cnf = cnf.Cnf(cnf_path)
            if control_dict["repair_NaN"]:
                if any([math.isnan(x) for x in sim_cnf.GENBOX.euler]):
                    sim_cnf.GENBOX.euler = [0.0, 0.0, 0.0]
            out_cnf = sim_cnf.write(out_prefix + ".cnf")

            # in_simSystem.cnf = sim_dir_cnfs[-1]
            # in_simSystem.cnf.path = out_prefix+".cnf"
            # in_simSystem.write_files(cnf=True, all=False)

    if control_dict["cp_omd"]:

        out_omd = out_prefix + ".omd"
        if os.path.exists(out_prefix + ".omd"):
            warnings.warn(
                "Skipping concatenation of tres, as there already exists one!\n "
                "If you want to generate a new trc. Please delete the old one!\n\t Path:" + out_prefix + ".omd"
            )
        else:
            if verbose:
                print("\tStart omd")
            # find all omd files in this project
            sim_dir_omd = gather_simulation_step_file_paths(
                in_folder, filePrefix=in_prefix, fileSuffixes=".omd", verbose=verbose  # in_simSystem.name,
            )
            for omd in sim_dir_omd:
                shutil.copy2(omd, out_prefix + omd.split("/")[-1])

    if control_dict["cat_trc"]:
        # find all trc files in this project
        out_traj_path = out_prefix + ".trc.h5"
        if os.path.exists(out_traj_path):
            warnings.warn(
                "Skipping concatenation of trcs, as there already exists one!\n "
                "If you want to generate a new trc. Please delete the old one!\n\t Path:" + out_traj_path
            )
        else:
            if verbose:
                print("\tStart Trc Cat")
            traj_files = gather_simulation_step_file_paths(
                in_folder, filePrefix=in_prefix, fileSuffixes=[".trc", ".trc.gz", ".trc.tar.gz"], verbose=verbose
            )
            if len(traj_files) > 0:
                out_trc_file = trc.Trc(traj_path=traj_files[0], in_cnf=out_cnf)
                if len(traj_files) > 0:
                    for tmp_trc_file in traj_files[1:]:
                        tmp_trc = trc.Trc(tmp_trc_file, in_cnf=out_cnf)
                        out_trc_file += tmp_trc
                out_trc_file.write(output_path=out_traj_path)
    if control_dict["cat_tre"]:
        out_traj_path = out_prefix + ".tre.h5"
        if os.path.exists(out_traj_path):
            warnings.warn(
                "Skipping concatenation of tres, as there already exists one!\n "
                "If you want to generate a new trc. Please delete the old one!\n\t Path:" + out_traj_path
            )
        else:
            if verbose:
                print("\tStart Tre Cat")
            traj_files = gather_simulation_step_file_paths(
                in_folder, filePrefix=in_prefix, fileSuffixes=[".tre", ".tre.gz", ".tre.tar.gz"], verbose=verbose
            )
            if len(traj_files) > 0:
                out_trc_file = tre.Tre(traj_files[0], auto_save=False)
                if len(traj_files) > 0:
                    for tmp_trc_file in traj_files[1:]:
                        tmp_trc = tre.Tre(tmp_trc_file, auto_save=False)
                        out_trc_file += tmp_trc
                out_trc_file.write(output_path=out_traj_path)

    if control_dict["cat_trg"]:
        out_traj_path = out_prefix + ".trg.h5"
        if os.path.exists(out_traj_path):
            warnings.warn(
                "Skipping concatenation of tres, as there already exists one!\n "
                "If you want to generate a new trc. Please delete the old one!\n\t Path:" + out_traj_path
            )
        else:
            if verbose:
                print("\tStart Tre Cat")
            traj_files = gather_simulation_step_file_paths(
                in_folder, filePrefix=in_prefix, fileSuffixes=[".trg", ".trg.gz", ".trg.tar.gz"], verbose=verbose
            )
            if len(traj_files) > 0:
                out_trg_file = trg.Trg(traj_files[0], auto_save=False)
                if len(traj_files) > 0:
                    for tmp_trg_file in traj_files[1:]:
                        tmp_trg = trg.Trg(tmp_trg_file, auto_save=False)
                        out_trg_file += tmp_trg
                out_trg_file.write(output_path=out_traj_path)

    """

    if (control_dict["convert_trcs"]):
        print("\tStart Trc Conversion")
        # wait for async job creating the trcs.
        if (submitted_trc_job):
            p_trc.wait()

        # get files:
        final_trc_files = list(
            sorted(glob.glob(out_folder + "/*.trc*"), key=lambda x: int(x.split("_")[-1].split(".")[0])))

        if (n_processes > 1):
            out_dcd = manager.dict()
            distributed_jobs = [
                (n, range(n, num_replicas, n_processes), final_trc_files, in_topology_path, gromosPP_bin_dir,
                 out_dcd,
                 fit_traj_to_mol, verbose) for
                n in range(n_processes)]
            p_conv = p.starmap_async(_thread_worker_conv_trc, distributed_jobs)
        else:
            out_dcd = {}
            _thread_worker_conv_trc(job=-1, replica_range=range(num_replicas), trc_files=final_trc_files,
                                    in_topology_path=in_topology_path,
                                    gromos_path=gromosPP_bin_dir, out_traj=out_dcd, fit_traj_to_mol=1,
                                    verbose=verbose,
                                    boundary_conditions=boundary_conditions)
    """
    if verbose:
        print("all jobs finished")
    return out_cnf  # in_simSystem


def gather_simulation_step_file_paths(
    in_folder: str,
    filePrefix: str = "",
    fileSuffixes: Union[str, List[str]] = [".tre", ".tre.tar.gz"],
    verbose: bool = False,
    finalNumberingSort=False,
) -> Dict[int, List[str]]:
    """gather_replica_file_paths

        Finds all trajectory paths in a simulation folder.


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
    ##init statewise dic
    files = []
    if verbose:
        print("SEARCH PATTERN: " + filePrefix + " + * +" + str(fileSuffixes))

    for dirname, dirnames, filenames in os.walk(in_folder):
        if str(dirname[-1]).isdigit() and os.path.basename(dirname).startswith("eq"):
            continue
        # check actual in_dir for file pattern
        tmp_files = [
            file
            for file in filenames
            if (filePrefix in file and any([file.endswith(suffix) for suffix in fileSuffixes]))
        ]

        if len(tmp_files) != 0:
            grom_file_prefix = sorted(
                [dirname + "/" + file for file in tmp_files if (any([suffix in file for suffix in fileSuffixes]))]
            )
            files += grom_file_prefix
        if verbose:
            print("walking to in_dir: ", os.path.basename(dirname), "found: ", len(tmp_files))
    files = sorted(files, key=lambda x: int(x.split("_")[-1].split(".")[0]))

    if verbose:
        print("\nfoundFiles:\n")
        print("\t" + "\t".join([y + "\n" for y in files]))

    if len(files) == 0:
        warnings.warn(
            "could not find any file with the prefix: "
            + str(filePrefix)
            + " * "
            + str(fileSuffixes)
            + " in folder : \n"
            + in_folder
        )
    return files
