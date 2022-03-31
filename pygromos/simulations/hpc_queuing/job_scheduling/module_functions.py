import os
from inspect import signature, _empty
from pygromos.utils.typing import List


def write_job_script(
    out_script_path: str,
    target_function: callable,
    variable_dict: dict,
    python_cmd: str = "python3",
    verbose: bool = False,
) -> str:
    if not os.path.exists(os.path.dirname(out_script_path)):
        raise IOError(
            "Could not find path of dir, that should contain the schedule script!\n\t Got Path: " + out_script_path
        )

    # Build str:
    s = signature(target_function)
    import_string = "#IMPORTS\n"
    import_string += "from " + str(target_function.__module__) + " import " + target_function.__name__
    vars_string = "#VARIABLES: \n"
    cmd_options = ""

    missed_keys = []
    for key in s.parameters:
        if key in variable_dict:
            value = variable_dict[key]
            if key == "in_simSystem":  # this is a nasty way! ... tends to fail!
                sys = value
                vars_string += sys.get_script_generation_command(var_name=key, var_prefixes="system")
            # elif(isinstance(value, Dict)):
            #    if(key  == "control_dict"):
            #        if(no_reeds_control_dict):
            #            vars_string += reeds_analysis.dict_to_nice_string(value)
            #        else:
            #            vars_string += reeds_analysis.dict_to_nice_string(reeds_analysis.check_script_control(value))
            #    else:
            #        vars_string +=  reeds_analysis.dict_to_nice_string(value)
            elif isinstance(value, List):
                vars_string += key + "= [ " + ", ".join(map(str, value)) + "]\n"
            elif isinstance(value, str):
                vars_string += key + ' = "' + str(value) + '"\n'
            else:
                vars_string += key + " = " + str(value) + "\n"
            cmd_options += key + "=" + key + ", "
        elif s.parameters[key].default == _empty:
            missed_keys.append(key)

    if len(missed_keys) > 0:
        raise ValueError(
            "Found some variables missing in variable dict,that are required!\n\t" + "\n\t".join(missed_keys)
        )

    cmd_string = "\n#DO\n"
    cmd_string += target_function.__name__ + "(" + cmd_options + ")"

    script_text = (
        "#!/usr/bin/env " + python_cmd + "\n\n" + import_string + "\n\n" + vars_string + "\n\n" + cmd_string + "\n"
    )
    if verbose:
        print(script_text)

    # write out file
    out_script_file = open(out_script_path, "w")
    out_script_file.write(script_text)
    out_script_file.close()

    return out_script_path
