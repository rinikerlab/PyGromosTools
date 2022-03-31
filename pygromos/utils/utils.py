"""utils

    This module should contain small usefull functions.

"""
import os
import math
import argparse
import inspect

from pygromos.utils.typing import Dict, List

spacer0 = "!" * 90 + "\n"
spacer = "#" * 80 + "\n"
spacer2 = "=" * 60 + "\n"
spacer3 = "_" * 40 + "\n"


time_wait_s_for_filesystem = 0


def _cartesian_distance(x1: float, x2: float, y1: float, y2: float, z1: float, z2: float) -> float:
    return math.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2)


"""
File and submission
"""


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ("yes", "true", "t", "y", "1"):
        return True
    elif v.lower() in ("no", "false", "f", "n", "0"):
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value expected.")


def dynamic_parser(func: callable, title: str):
    """
        This function builds dynamically a parser obj for any function, that has parameters with annotated types.
        Result is beeing able to parse any function dynamically via bash.

    Parameters
    ----------
    func : callable
        the function that should be parsed
    title : str
        title for the parser

    Returns
    -------
    args
        The parsed arguments

    Raises
    ------
    IOError
        error if a parsing arg is unknown
    """
    parser = argparse.ArgumentParser(description=title)
    args = inspect.getfullargspec(func)
    total_defaults = len(args.defaults)
    total_args = len(args.args)
    total_required = total_args - total_defaults

    parser.description = func.__name__ + " - " + func.__doc__.split("Parameters")[0]
    for argument, argument_type in args.annotations.items():
        index = args.args.index(argument)
        required = True if (index < total_required) else False
        default = None if (required) else args.defaults[index - total_required]
        if argument_type is bool:
            argument_type = str2bool
        parser.add_argument("-" + argument, type=argument_type, required=required, default=default)

    args, unkown_args = parser.parse_known_args()
    if len(unkown_args) > 0:
        raise IOError(__name__ + " got unexpected argument(s) for parser:\n" + str(unkown_args))
    return args


def dict_to_nice_string(control_dict: Dict) -> str:
    """
        Converts a dictionary of options (like template_control_dict)
          to a more human readable format. Which can then be printed to a text file,
          which can be manually modified before submiting analysis jobs.

    Parameters
    ----------
    control_dict : Dict
        analysis control dictonary

    Returns
    -------
    str
        nice formatting of the control dictionary for printing.

    """
    script_text = "control_dict = {\n"
    for key, value in control_dict.items():
        script_text += '\t"' + key + '": '
        first = False
        if type(value) == dict:
            if "do" in value:  # do should always be first in this list
                script_text += '{"do":' + str(value["do"]) + ","
                if len(value) > 1:
                    script_text += "\n"
                first = True
            for key2, value2 in value.items():  # alternative keys

                # prefix
                if first:
                    prefix = " "
                    first = False
                else:
                    prefix = "\t\t"

                # key_val
                if key2 == "do":
                    continue
                elif type(value2) == dict:
                    script_text += prefix + '"' + str(key2) + '": ' + _inline_dict(value2, "\t\t\t") + ",\n"
                else:
                    script_text += prefix + '"' + str(key2) + '": ' + str(value2) + ","
            script_text += prefix + " },\n"
        else:
            script_text += str(value) + ",\n"
    script_text += "}\n"
    return script_text


def _inline_dict(in_dict: Dict, prefix: str = "\t") -> str:
    """
        translate dictionary to one code line. can be used for meta-scripting

    Parameters
    ----------
    in_dict: Dict
        analysis control dict
    prefix : str, optional
        prfix symbol to dict write out.

    Returns
    -------
    str
        code line.

    """
    msg = "{\n"
    for key, value in in_dict.items():
        if type(value) == dict:
            msg += prefix + '"' + str(key) + '": ' + _inline_dict(in_dict=value, prefix=prefix + "\t") + ","
        else:
            msg += prefix + '"' + str(key) + '": ' + str(value) + ",\n"
    return msg + prefix + "}"


def write_job_script(
    out_script_path: str,
    target_function: callable,
    variable_dict: dict,
    python_cmd: str = "python3",
    verbose: bool = False,
) -> str:
    """
        this function writes submission commands into a file. The command will be started from a bash env into python.


    Parameters
    ----------
    out_script_path: str
        path of the output script.
    target_function : callable
        the function, that shall be submitted
    variable_dict : dict
        variables for this function
    python_cmd : str, optional
        which python command shall be supplied
    verbose : bool, optional
        c'est la vie

    Returns
    -------
    str
        returns an out script path.

    Raises
    ------
    IOERROR
        if outpath is not possible
    ValueError
        if required variable from the var-dict for the function is missing
    """

    if not os.path.exists(os.path.dirname(out_script_path)):
        raise IOError(
            "Could not find path of dir, that should contain the schedule script!\n\t Got Path: " + out_script_path
        )
    import pygromos

    # Build str:
    s = inspect.signature(target_function)  # to lazy for numpydoc
    import_string = ""
    import_string += "#IMPORTS\n"
    import_string += "import sys\n"
    import_string += 'sys.path.append("' + os.path.abspath(os.path.dirname(pygromos.__file__) + "/..") + '")\n'
    import_string += "from " + str(target_function.__module__) + " import " + target_function.__name__
    vars_string = "#VARIABLES: \n"
    cmd_options = ""

    from pygromos.files.gromos_system.gromos_system import Gromos_System
    from pygromos.simulations.hpc_queuing.submission_systems._submission_system import _SubmissionSystem

    missed_keys = []
    for key in s.parameters:
        if key in variable_dict:
            value = variable_dict[key]
            if isinstance(value, Gromos_System) or issubclass(
                value.__class__, _SubmissionSystem
            ):  # this is a nasty way! ... tends to fail!
                sys = value
                vars_string += sys.get_script_generation_command(var_name=key)
            elif isinstance(value, Dict):
                vars_string += dict_to_nice_string(value)
            elif isinstance(value, List):
                vars_string += key + "= [ " + ", ".join(map(str, value)) + "]\n"
            elif isinstance(value, str):
                vars_string += key + ' = "' + str(value) + '"\n'
            else:
                vars_string += key + " = " + str(value) + "\n"
            cmd_options += key + "=" + key + ", "
        elif s.parameters[key].default == inspect._empty:
            missed_keys.append(key)

    if len(missed_keys) > 0:
        raise ValueError(
            "Found some variables missing in variable dict,that are required!\n\t" + "\n\t".join(missed_keys)
        )

    cmd_string = "\n#DO\n"
    cmd_string += target_function.__name__ + "(" + cmd_options + ")"
    cmd_string += "\nexit(0)\n\n"

    script_text = (
        "#!/usr/bin/env " + python_cmd + "\n\n" + import_string + "\n\n" + vars_string + "\n\n" + cmd_string + "\n"
    )
    if verbose:
        print(script_text)

    # write out file
    out_script_file = open(out_script_path, "w")
    out_script_file.write(script_text)
    out_script_file.close()

    os.system("chmod +x " + out_script_path)

    return out_script_path


def nice_s_vals(svals: List[float]) -> List[float]:
    """
        This helper function formats s-vals for RE-EDS to nice readable values.
        Is/was used in RE-EDS applications. Main functionality is rounding the different s-vals to their significance digits.


    Parameters
    ----------
    svals : List[float]
        smoothing parameters of a reeds approach

    Returns
    -------
        List[float]
            nicely rounded s-values, such that only significant digits for each number is around.

    """

    nicer_labels = []
    for val in svals:
        nicer_labels.append(round(float(val), str(val).count("0") + 2))
    return nicer_labels
