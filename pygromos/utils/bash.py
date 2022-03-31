"""bash-wrapper

Description
-----------
   This lib contains bash-wrappers for general use in a workflow, with nice error messages and other features.
   Additionally there are some convenience functions for error managment and dependencies check

:author: Benjamin Schroeder
"""

import io
import os
import glob
import time
import warnings
import shutil

import subprocess as sub
from pygromos.utils.utils import time_wait_s_for_filesystem
from pygromos.utils.typing import Union, List, Dict

#################################
#   General functions:


def wait_for_fileSystem(
    check_paths: Union[str, List[str]],
    regex_mode: bool = False,
    max_waiting_iterations: int = 1000,
    verbose: bool = False,
) -> bool:
    """
    This function can be used to circumvent lsf lag times.

    Parameters
    ----------
    check_path: str - Path to file to check if existant
    max_waiting_iterations: int - maximal iteration Time
    verbose: bool - tell me if found

    Returns
    -------
    True
        on success
    """
    if isinstance(check_paths, str):
        check_paths = [check_paths]

    for check_path in check_paths:
        it = 0
        waiting = True
        while waiting and it < max_waiting_iterations:
            if regex_mode:
                waiting = len(glob.glob(check_path)) <= 0
            else:
                waiting = not os.path.exists(check_path)
            time.sleep(time_wait_s_for_filesystem)
            it += 1
        if waiting:
            raise IOError("Could not find file: " + check_path)
        elif verbose:
            print("File Check FOUND: \t", check_path)

    return True


def check_path_dependencies(
    check_required_paths: Union[Dict[any, str], List[str]],
    check_warn_paths: Union[str, List[str]] = [],
    verbose: bool = True,
) -> str:
    """check_path_dependencies

        checks a list of dependencies if each path is present or not. throws an
        IOError, if an Path is not existing.

    Parameters
    ----------
    check_required_paths :    Union[t.Dict[any, str], List[str]]
        if path does not exist, an error wil be raised
    check_warn_paths:   Union[t.Dict[any, str], List[str]]
        if path does not exist, a warning will be written out.
    verbose :   bool, optional

    Raises
    ------
    IOERROR
        if a file is not existing

    Returns
    -------
    Union[t.Dict[any, str], List[str]]
        in_dependencies

    """

    found_error = False
    missing = []
    if verbose and type(check_required_paths) is list:
        print("\n\n==================\n\tCHECK dependencies\n")
        print("\n".join(list(map(lambda s: "Check " + str(s), check_required_paths))))
    elif verbose and type(check_required_paths) is dict:
        print("\nCHECK dependencies")
        print(
            "\n".join(list(map(lambda s: "Check " + str(s), [check_required_paths[x] for x in check_required_paths])))
        )

    # ERROR
    if type(check_required_paths) is dict:
        for x in check_required_paths:
            if "*" in x or "?" in x:
                if verbose:
                    print("Skipping regex")
                continue
            if verbose:
                print(x)
            if not isinstance(check_warn_paths[x], str) or (
                isinstance(check_required_paths[x], str) and not os.path.exists(check_required_paths[x])
            ):
                found_error = True
                missing.append(x)
    elif type(check_required_paths) is list:
        for x in check_required_paths:
            if "*" in x or "?" in x:
                if verbose:
                    print("Skipping regex")
                continue
            if verbose:
                print(x)
            if not isinstance(x, str) or (isinstance(x, str) and not os.path.exists(x)):
                found_error = True
                missing.append(x)

    # WARN
    if type(check_warn_paths) is dict:
        for x in check_warn_paths:
            if verbose:
                print(x)
            if not isinstance(check_required_paths[x], str) or (
                isinstance(check_required_paths[x], str) and not os.path.exists(check_required_paths[x])
            ):
                warnings.warn("\tDid not find: " + str(x) + " with path: " + check_required_paths[x])
    elif type(check_warn_paths) is list:
        for x in check_warn_paths:
            if verbose:
                print(x)
            if not isinstance(x, str) or (isinstance(x, str) and not os.path.exists(x)):
                warnings.warn("\tDid not find: " + str(x))

    if found_error:
        print("\n==================\nAUTSCH\n==================\n")
        missing_str = "\n\t".join(map(str, missing))
        raise IOError("COULD NOT FIND all DEPENDENCY!\n\t Could not find path to: \n\t" + str(missing_str), "\n\n")
    elif verbose:
        print("All dependencies are correct!", "\n\n")
    return 0


#################################
#   bash wrapper:


def extract_tar(
    in_path: str,
    out_path: str = None,
    gunzip_compression: bool = False,
    remove_tar_afterwards: bool = False,
    verbose: bool = False,
) -> str:
    """extract_tar

        this wrapper helps you to unpack a tar file via the OS
        -- since I always used to forget the command...
    Parameters
    ----------
    in_path :   str
    out_path :  str
        optional output, to which the untared in is moved.
    gunzip_compression : bool, optional
        is the tar file in the gunzip compression format? (-z added)
    remove_tar_afterwards : bool, optional
    verbose:    bool, optional

    Raises
    -------
    OSError
        if the bas command fails an OSError will be raised.

    Returns
    -------
    str
        out_path
    FileIO
        process return log.
    """

    option = "-"
    if gunzip_compression:
        option += "xzf"
    else:
        option += "xf"

    command = "tar " + option + " " + in_path

    if isinstance(out_path, str) and in_path.replace(".tar", "").replace(".gz", "") != out_path:
        command += " && mv  " + in_path.replace(".tar", "").replace(".gz", "") + " " + out_path
    else:
        out_path = in_path.replace(".tar", "").replace(".gz", "")

    if out_path == in_path:
        raise Exception(
            "Outpath is not allowed to be equal with in_path!\n In_path: " + in_path + "\n Out_path:" + out_path + "\n"
        )

    if verbose:
        print("cmd: ", command)
    orig_path = os.getcwd()
    os.chdir(os.path.dirname(in_path))
    ret = execute(command, verbose=verbose)
    os.chdir(orig_path)

    if verbose:
        print("\n".join(ret.readlines()))
    if remove_tar_afterwards:
        if os.path.exists(out_path):
            remove_file(in_path)

    wait_for_fileSystem(out_path)

    return out_path


def compress_tar(
    in_path: str,
    out_path: str = None,
    gunzip_compression: bool = False,
    remove_in_file_afterwards: bool = False,
    remove_in_dir_afterwards: bool = False,
    verbose: bool = False,
) -> str:
    """compress_tar

        compress a file or directory with tar via the OS

    Parameters
    ----------
    in_path :   str
    out_path :  str
    gunzip_compression :    bool, optional
    remove_in_file_afterwards : bool, optional
    verbose:    bool, optional

    Raises
    -------
    OSError
        if the bas command fails an OSError will be raised.

    Returns
    -------
    str
        out_path
    FileIO
        process return log.
    """
    print(out_path)
    if out_path is None:
        out_path = in_path
    if not out_path.endswith(".tar.gz") and gunzip_compression:
        out_path += ".tar.gz"
    elif not out_path.endswith(".tar") and not gunzip_compression:
        out_path += ".tar"

    option = "-c"
    if gunzip_compression:
        option += "z"
    option += "f"

    command = "tar " + option + " " + out_path + " " + os.path.basename(in_path) + " "

    # command
    orig_path = os.getcwd()
    os.chdir(os.path.dirname(in_path))
    execute(command, verbose=verbose)
    os.chdir(orig_path)

    wait_for_fileSystem(out_path, verbose=verbose)

    if remove_in_dir_afterwards:
        remove_file(in_path, recursive=True)
    elif remove_in_file_afterwards:
        remove_file(in_path)

    return out_path


def compress_gzip(
    in_path: str, out_path: str = None, extract: bool = False, force: bool = True, verbose: bool = False
) -> str:
    """compress_gzip

        compress a file or directory with tar via the OS

    Parameters
    ----------
    in_path :   str
    out_path :  str
    gunzip_compression :    bool, optional
    remove_in_file_afterwards : bool, optional
    verbose:    bool, optional

    Raises
    -------
    OSError
        if the bas command fails an OSError will be raised.

    Returns
    -------
    str
        out_path
    FileIO
        process return log.
    """

    if isinstance(out_path, type(None)) and not extract:
        out_path = in_path + ".gz"
    elif isinstance(out_path, type(None)) and extract:
        out_path = in_path.replace(".gz", "").replace(".tar", "")

    option = ""

    if extract:
        option += " -d "
    if force:
        option += " -f "

    command = "gzip " + option + " " + in_path + " "

    # command
    execute(command, verbose=verbose)
    if in_path + ".gz" != out_path and not extract:
        wait_for_fileSystem(in_path + ".gz", verbose=verbose)
        out_path = move_file(in_path + ".gz", out_file_path=out_path)
        wait_for_fileSystem(out_path, verbose=verbose)
    elif in_path != out_path + ".gz" and extract:
        wait_for_fileSystem(in_path.replace(".gz", ""), verbose=verbose)
        out_path = move_file(in_path, out_file_path=out_path)
        wait_for_fileSystem(out_path, verbose=verbose)

    return out_path


def copy_file(
    in_file_path: str, out_file_path: str, copy_a_directory: bool = False, additional_option: str = ""
) -> str:
    """copy_file

        copy files via the OS
        TODO: refactor to use execute command!

    Parameters
    ----------
    in_file_path :  str
    out_file_path : str
    copy_a_directory: bool, optional
    additional_option : str, optional

    Raises
    -------
    OSError
        if the bas command fails an OSError will be raised.

    Returns
    -------
    str
        out_file_path

    """

    if copy_a_directory:
        additional_option += " -r "

    copy_files = "cp " + str(additional_option) + " " + str(in_file_path) + " " + str(out_file_path) + "\n"

    if os.system(copy_files):
        raise OSError(
            "could not copy:\n "
            + str(in_file_path)
            + "\n \t to \n"
            + str(out_file_path)
            + "\n \t options: "
            + str(additional_option)
        )

    wait_for_fileSystem(out_file_path)

    return out_file_path


def move_file(in_file_path: str, out_file_path: str, additional_options: str = "", verbose: bool = False) -> str:
    """move_file

        copy files via the OS
        TODO: refactor to use execute command!

    Parameters
    ----------
    in_file_path :  str
    out_file_path : str
    additional_options : str, optional

    Raises
    -------
    OSError
        if the bash command fails an OSError will be raised.

    Returns
    -------
    str
        out_file_path

    """
    regex = False

    if "*" in in_file_path or "?" in in_file_path:
        regex = True
    elif os.path.isfile(in_file_path) and os.path.isdir(out_file_path):
        out_file_path = os.path.dirname(out_file_path) + "/" + os.path.basename(in_file_path)
    command = "mv " + additional_options + " " + in_file_path + " " + out_file_path + "\n"

    try:
        execute(command=command, verbose=verbose)
    except Exception as err:
        raise ValueError("BASH could not move file! " + "\n".join(map(str, err.args)))

    if not regex:
        wait_for_fileSystem(out_file_path)
    return out_file_path


def concatenate_text_files(in_file_paths: List[str], out_file_path: str, verbose: bool = False):
    """concatenate_text_files

        concatenated multiple text files via the OS
        TODO: refactor to use execute command!

    Parameters
    ----------
    in_file_paths :  List[str]
    out_file_path : str
    verbose: bool, optional

    Raises
    -------
    OSError
        if the bash command fails an OSError will be raised.
    IOError
        if the the directory for the output file cann  not be found or not all in_paths exist.

    Returns
    -------
    None

    """
    if isinstance(in_file_paths, str):
        in_file_paths = [in_file_paths]

    if isinstance(in_file_paths, str):
        in_file_paths = [in_file_paths]

    command = "cat " + " ".join(in_file_paths) + " > " + out_file_path + " \n"

    if verbose:
        print("CONCAT: " + command)
    if os.path.exists(os.path.dirname(out_file_path)):
        if all([os.path.exists(in_path) for in_path in in_file_paths]):
            if os.system(command):
                raise OSError("could not concate files:\n " + str(" ".join(in_file_paths)) + "\n to\n" + out_file_path)
        else:
            raise IOError(
                "could not find all in_paths!:\n "
                + "\n".join([in_path for in_path in in_file_paths if (not os.path.exists(in_path))])
            )
    else:
        raise IOError("could not find folder for:\n " + os.path.dirname(out_file_path))

    wait_for_fileSystem(out_file_path)
    return out_file_path


def replace_text_in_text_file(
    in_file_path: str, find_pattern: str, replace_pattern: str, out_file_path: str = None
) -> str:
    """replace_text_in_text_file

        this file replaces a regex pattern in a text file

    Parameters
    ----------
    in_file_path :   str
    find_pattern :  str
    replace_pattern :   str
    out_file_path :  str, optional
        if no out_file_path given, the text is replaced in place.

    Raises
    -------
    OSError
        If the bash command failed.

    Returns
    -------
    str
        out_file
    """

    command = ""
    if isinstance(out_file_path, type(None)):
        command = "sed s/" + find_pattern + "/" + replace_pattern + "/g " + in_file_path + " > " + out_file_path + " \n"
    else:
        command = "sed s/" + find_pattern + "/" + replace_pattern + "/g " + in_file_path + " \n"
        out_file_path = in_file_path

    if os.system(command):
        raise OSError("could not replace text in file:\n " + str(" ".join(in_file_path)) + "\n to\n" + out_file_path)

    wait_for_fileSystem(out_file_path)
    return out_file_path


def make_folder(in_directory_path: str, additional_option: str = "", verbose: bool = False) -> str:
    """make_folder

        make a new folder via the os.

    Parameters
    ----------
    in_directory_path :    str
    additional_option : str, optional
    verbose :   bool, optional

    Raises
    -------
    OSError
        If the bash command failed.

    Warnings
    ---------
        if folder already existed, a warning will be printed.

    Returns
    -------
    str
        directory_path
    """
    mk_folder = "mkdir " + additional_option + " " + str(in_directory_path) + "\n"
    if not os.path.isdir(in_directory_path):
        if os.system(mk_folder):
            raise OSError("could not make folder:\n " + str(in_directory_path))

    elif verbose:
        warnings.warn("WARNING:\n\t Did not build already existing folder: " + in_directory_path, category=UserWarning)

    wait_for_fileSystem(in_directory_path)
    return in_directory_path


def remove_file(in_file_path: str, recursive: bool = False, additional_options: str = "") -> None:
    """remove_file

        delete a file via the OS.
        TODO: refactor to use execute command!

    Parameters
    ----------
    in_file_path :  str
    recursive: bool, optional
        remove recursivly! CAREFULL is DANGEROUS!
    additional_options : str, optional

    Raises
    -------
    OSError
        if the bash command fails an OSError will be raised.

    Returns
    -------
    None

    """
    if recursive:
        additional_options += " -r "
    rm_command = "rm " + str(in_file_path) + " " + str(additional_options)
    if os.path.exists(in_file_path):
        if os.system(rm_command):
            raise OSError(
                "could not delete file/folder: " + str(in_file_path) + "\n options: " + str(additional_options)
            )


def remove_folder(in_directory_path: str, additional_options: str = "", verbose: bool = False) -> None:
    """remove_folder

        make a new folder via the os.

    Parameters
    ----------
    in_directory_path :    str
    additional_options : str, optional
    verbose :   bool, optional

    Raises
    -------
    OSError
        If the bash command failed.

    Warnings
    ---------
        if folder did not exist, a warning will be printed.

    Returns
    -------
        None
    """

    remove_folder = "rmdir " + additional_options + " " + str(in_directory_path) + "\n"
    if os.path.isdir(in_directory_path):
        if os.system(remove_folder):
            raise OSError("could not remove folder:\n " + str(in_directory_path))

    elif verbose:
        warnings.warn("Warning! Did not remove non existing folder: " + in_directory_path, category=UserWarning)


def save_make_folder(in_directory_path: str, additional_options: str = "") -> str:
    """save_make_folder
        This function checks if a folder already exists and adds an integer to it, so that there will be a guaranteed unique folder.


    Parameters
    ----------
    in_directory_path : str
    additional_options :    str, optional

    Raises
    -------
    OSError
        If the bash command failed.

    Returns
    -------
    str
        in_directory_path
    """

    offset = 1
    dir_versions = list(filter(lambda x: ".tar" not in x or ".gz" not in x, sorted(glob.glob(in_directory_path + "*"))))
    print("dirversions:", dir_versions)
    if len(dir_versions) > 0:
        last_dir = dir_versions[len(dir_versions) - 1]
        print("last:", last_dir)

        suffix = str(last_dir.replace(in_directory_path + "_", ""))
        print(suffix)
        if suffix != "" and suffix.isalnum():
            offset = int(suffix) + 1
        in_directory_path += "_" + str(offset)

    mk_folder = "mkdir " + additional_options + " " + str(in_directory_path) + "\n"
    if os.system(mk_folder):
        raise Exception("could not make folder:\n " + str(in_directory_path))

    return in_directory_path


def link_folder(in_directory_path: str, out_link_path: str, additional_options: str = ""):
    """link_folder

        linking a folder via the OS

    Parameters
    ----------
    in_directory_path : str
    out_link_path : str
    additional_options :    str, optional

    Raises
    -------
    IOError
        If already a file with the link exists.

    OSError
        If the bash command failed.

    Returns
    -------
    str
        out_link_path

    """

    link_folders = "ln -s " + additional_options + " " + in_directory_path + " " + out_link_path + "\n"
    if not os.path.exists(out_link_path):
        if os.system(link_folders):
            raise OSError(
                "could not link:\n "
                + str(in_directory_path)
                + "\n \t to \n"
                + str(out_link_path)
                + "\n \t options: "
                + str(additional_options)
            )
    else:
        raise IOError("A link with the give path already exists!: \n path:\t" + out_link_path)
    return out_link_path


def execute_os(command: Union[str, List[str]], verbose: bool = False) -> io.FileIO:
    """execute

        DEAPRECIATED

        sadly I can only use os.popen, as otherwise euler refuses.
        this command executes a command on the os-layer. (e.g.: on linux a bash command)
        therefore it uses os.popen

    Parameters
    ----------
    command :   str, Lists[str]
    verbose: bool, optional

    Raises
    ------
    OSError
        if bash command fails an OSError is raised

    Returns
    -------
    io.FileIO
        returns the execution log of the process.

    """

    class dummyProcess:
        def __init__(self, stdout, stderr, ret):
            self.stdout = stdout
            self.stderr = stderr
            self.poll = lambda x: int(ret)

    if type(command) == list:
        command = " ".join(command)

    if verbose:
        print(command + "\n")

    try:
        ret = os.popen(command)
        print(ret)
    except Exception as err:
        raise OSError(
            "could not execute bash command:\n  error: "
            + "\n\t".join(err.args)
            + "\n\t"
            + str(command)
            + "\n\tCommand returned: \t"
            + str(ret.read())
        )

    if verbose:
        print("\t" + "\n\t".join(ret.readlines()) + "\n")

    return dummyProcess(stdout=ret, stderr=[], ret=0)


def execute_subprocess(
    command: Union[str, List[str]], catch_STD: Union[bool, str] = False, env: dict = None, verbose: bool = False
) -> sub.CompletedProcess:
    """execute_subprocess
        This command starts a subprocess, that is executing the str command in bash.

    Parameters
    ----------
    command : str
        bash command
    catch_STD :
        if bool: catch the output and past it into the command line
        if str: catch output and write it into this file
    env: dict
        environment
    verbose : bool
        verbosity

    Returns
    -------
    CompletedProcess
        return the executed process obj. (from subprocess)
    """

    if isinstance(command, list):
        command = " ".join(command)
    if verbose:
        print("\texecute command: \n\t" + command)

    kwargs = {}
    if isinstance(catch_STD, bool):
        kwargs.update({"stdout": sub.PIPE})
    elif isinstance(catch_STD, str):
        kwargs.update({"stdout": open(catch_STD, "w")})

    if env is None:
        env = os.environ.copy()

    p = sub.Popen(args=command, shell=True, stderr=sub.PIPE, env=env, **kwargs)

    # print(p, vars(p))
    try:
        p.wait(120)  # Wait for process to finish
    except sub.TimeoutExpired:
        warnings.warn("TIME OUT WITH: " + str(command))
        print("Continue Waiting: ")
        p.wait()  # Wait for process to finish
    p.terminate()  # Make sure its terminated
    r = p.poll()
    if r:  # Did an Error occure?
        msg = "SubProcess Failed due to returncode: " + str(r) + "\n COMMAND: \n\t" + str(command)
        msg += "\nSTDOUT:\n\t"
        msg += "NONE" if (p.stdout is None) else "\n\t".join(map(str, p.stdout.readlines()))
        msg += "\nSTDERR:\n\t"
        msg += "NONE" if (p.stdout is None) else "\n\t".join(map(str, p.stderr.readlines()))
        raise ChildProcessError(msg)
    if verbose:
        print("RETURN: ", r)

    return p


def execute_old(
    command: Union[str, List[str]],
    verbose: bool = False,
    ignore_return_code: bool = False,
    wait_fail=False,
    out_cnf_path: str = None,
) -> io.FileIO:
    """execute
        this command executes a command on the os-layer. (e.g.: on linux a bash command)

    Parameters
    ----------
    command :   str, Lists[str]
    verbose: bool, optional

    Raises
    ------
    OSError
        if bash command fails an OSError is raised

    Returns
    -------
    io.FileIO
        returns the execution log of the process.

    """

    if isinstance(command, str):
        command = command.split()

    # TODO: maybe pass path directly?
    # This block overwrites the pipe of the sub process
    std_out = sub.PIPE
    std_err = sub.PIPE

    if verbose:
        print("COMMAND: " + " ".join(command) + "\n")

    try:
        p = sub.Popen(command, stdout=std_out, stderr=std_err)
        # testing this!:
        try:
            p.wait(timeout=15)
        except sub.TimeoutExpired:
            warnings.warn("Wait threw error!:( for cmd: " + str(command) + " If it is a long command... ok :)")

            if wait_fail:
                p.kill()
                raise TimeoutError("Wait failed for proces: cmd:" + str(command))

        try:
            while p.poll() is None:
                time.sleep(time_wait_s_for_filesystem)
        except KeyboardInterrupt:
            p.kill()

        if verbose:
            print("start writing")
        ret_stdout = io.StringIO("".join(map(lambda x: x.decode("utf-8"), p.stdout.readlines())))
        ret_stderr = "\n".join(map(lambda x: x.decode("utf-8"), p.stderr.readlines()))
    except TimeoutError as err:
        raise err
    except Exception as err:
        raise ChildProcessError(
            "Process failed executing Bash command:\n Command:\n\t"
            + str(command)
            + "\n error: "
            + "\n\t".join(map(str, err.args))
            + "\n"
        )

    # bash command failed?
    if p.returncode > 0 and not ignore_return_code:
        raise OSError(
            "Bash command return code was Non-Zero!\nRETURN CODE: "
            + str(p.returncode)
            + "\nERR:\n"
            + ret_stderr
            + "\nstdOut:\n\t"
            + str("\n\t".join(ret_stdout.readlines()))
            + "\nCommand: "
            + str(command)
        )

    if verbose:
        print("STDOUT:\n \t" + "\n\t".join(ret_stdout.readlines()) + "\n\n")
    if verbose:
        print("STDERR:\n \t" + ret_stderr + "\n\n")
    if verbose:
        print("RETURN CODE: \t" + str(p.returncode) + "\n")

    p.terminate()
    del p
    return ret_stdout


def execute(
    command: Union[str, List[str]], verbose: bool = False, catch_STD: Union[bool, str] = False, env: dict = None
):
    return execute_subprocess(command=command, verbose=verbose, catch_STD=catch_STD, env=env)


def which(command: str) -> str:
    """Finds the full path of a command.

    Args:
        command (str): Name of the command.

    Returns:
        str: Full path of the command.
    """

    return shutil.which(command)


def command_exists(command: str) -> bool:
    """Does the command exists in the current system / path?

    Args:
        command (str): Name of the command.

    Returns:
        bool: Does command exists or not.
    """

    path = which(command)

    if path is None:
        return False
    return True


def path_exists(path: str) -> bool:
    """Does the provided path exists? Gives no information on whether
    it is a directory or a file.

    Args:
        path (str): Path of the file or directory.

    Returns:
        bool: Does path exists or not.
    """
    return os.path.exists(path)


def is_directory(path: str) -> bool:
    """Is the provided path a directory.

    Args:
        path (str): Path of the file or directory.

    Returns:
        bool: Is a directory or not.
    """
    return os.path.isdir(path)


def is_file(path: str) -> bool:
    """Is the provided path a file.

    Args:
        path (str): Path of the file or directory.

    Returns:
        bool: Is a file or not.
    """
    return os.path.isfile(path)


def directory_exists(path: str) -> bool:
    """Tests whether the provided path is valid and
    also is a directory. Returns false if either condition
    is not fullfilled.

    Args:
        path (str): Path of the file or directory.

    Returns:
        bool: Is a directory with a valid path or not.
    """
    if is_directory(path) and path_exists(path):
        return True
    return False


def file_exists(path: str) -> bool:
    """Tests whether the provided path is valid and
    also is a file. Returns false if either condition
    is not fullfilled.

    Args:
        path (str): Path of the file or directory.

    Returns:
        bool: Is a file with a valid path or not
    """

    if is_file(path) and path_exists(path):
        return True
    return False
