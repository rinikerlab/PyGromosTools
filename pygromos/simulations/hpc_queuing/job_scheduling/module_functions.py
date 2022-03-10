import os
from inspect import signature, _empty
from typing import List
from pygromos.files import imd
from pygromos.files.coord import cnf as cnf_cls
from pygromos.utils import amino_acids as aa


def adapt_imd_template(
    in_template_imd_path: str,
    out_imd_dir: str,
    cnf: cnf_cls.Cnf,
    exclude_residues: list = [],
    simulation_steps: int = 100000,
    solvent_keyword: str = "SOLV",
):
    """
    .. autofunction: modify_imds - for generateOptimizedStates
        This function prepares the imd. protocol file for a gromos simulation.
    :return:    imd_template_path, lig_nums
    :rtype: str, list
    """
    orig_residues = cnf.get_residues()

    ignore_residues = (
        lambda res: res != solvent_keyword
        and res not in aa.three_letter_aa_lib
        and res != "prot"
        and (res not in exclude_residues)
    )
    protein_residues = {res: val for res, val in orig_residues.items() if (res in aa.three_letter_aa_lib)}
    n_atoms_prot = sum(
        [sum(list(orig_residues[res].values())) for res in orig_residues if res in aa.three_letter_aa_lib]
    )  # if protein is present!
    prot_position = min([min(val) for res, val in protein_residues.items()]) if (n_atoms_prot > 0) else 0

    # adapt resis
    residues = {res: val for res, val in orig_residues.items() if (res not in aa.three_letter_aa_lib)}

    # get ligand parameters
    lig_atoms = sum([sum(list(residues[res].values())) for res in residues if ignore_residues(res)])
    lig_num = sum([1 for res in residues if ignore_residues(res)])
    lig_names = [res for res in residues if ignore_residues(res)]
    lig_position = min([min(residues[res]) for res in residues if ignore_residues(res)])

    if n_atoms_prot > 1:
        print("protein_position:", prot_position)
        residues.update({"prot": {prot_position: n_atoms_prot}})

    if solvent_keyword not in residues:
        residues.update({solvent_keyword: 0})

    imd_file = imd.Imd(in_template_imd_path)
    imd_file.SYSTEM.NSM = int(round(residues["SOLV"] / 3)) if ("SOLV" in residues) else 0
    imd_file.FORCE.adapt_energy_groups(residues)
    imd_file.STEP.NSTLIM = simulation_steps
    imd_file.STEP.NSTLIM = simulation_steps

    # hack for TIP3P explicit solvent
    if len(exclude_residues) > 0 and "prot" not in residues:
        solvent_bath = (
            lig_atoms + sum([sum(list(residues[x].values())) for x in exclude_residues]) + residues[solvent_keyword]
        )
        temp_baths = {lig_atoms: 1, solvent_bath: 2}

    # Temperature baths
    elif "prot" in residues:
        temp_baths = {}
        if len(exclude_residues) < 0:  # TODO: DIRTY HACK: in PNMT is Cofactor at end of the file.
            solvent_bath = lig_atoms + n_atoms_prot + residues[solvent_keyword]
        else:
            solvent_bath = (
                lig_atoms
                + n_atoms_prot
                + sum([sum(list(residues[x].values())) for x in exclude_residues])
                + residues[solvent_keyword]
            )

        if lig_position < prot_position:
            temp_baths = {lig_atoms: 1, (lig_atoms + n_atoms_prot): 2, solvent_bath: 3}
        else:
            temp_baths = {n_atoms_prot: 1, (lig_atoms + n_atoms_prot): 2, solvent_bath: 3}
    else:
        temp_baths = (
            {lig_atoms: 1, (lig_atoms + residues[solvent_keyword]): 2}
            if (solvent_keyword in residues)
            else {lig_atoms: 1}
        )

    if not isinstance(imd_file.MULTIBATH, type(None)):
        imd_file.MULTIBATH.adapt_multibath(last_atoms_bath=temp_baths)

    # edit EDS part
    imd_template_path = out_imd_dir + "/opt_structs_" + "_".join(lig_names)
    s_values = [1.0 for x in range(lig_num)]
    for state, s_values in zip(range(1, lig_num + 1), s_values):
        imd_file.edit_EDS(
            NUMSTATES=lig_num, S=s_values, EIR=[500 if (x == state) else -500 for x in range(1, lig_num + 1)]
        )
        imd_file.write(imd_template_path + "_" + str(state) + ".imd")

    return imd_template_path, s_values, lig_num


def write_job_script(
    out_script_path: str,
    target_function: callable,
    variable_dict: dict,
    python_cmd: str = "python3",
    verbose: bool = False,
    no_reeds_control_dict: bool = False,
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


def build_MD_analysis_script(
    in_simulation_name: str,
    in_folder: str,
    in_topology_path: str,
    out_dir: str,
    out_script_path: str,
    in_ene_ana_lib: str,
    gromosPP_path: str,
):
    """
        BAD FUNCTION! HISTORIC RELICT
    Parameters
    ----------
    in_simulation_name
    in_folder
    in_topology_path
    out_dir
    out_script_path
    in_ene_ana_lib
    gromosPP_path

    Returns
    -------

    """
    import_text = "#!/usr/bin/env python3 \n\n"
    import_text += "#################################\n"
    import_text += "# This is a small automatic generated ana_script wrapper.\n"
    import_text += "#################################\n"
    import_text += "from reeds.function_libs.jobScheduling_scripts import MD_simulation_analysis as ana\n\n"

    dependencies_text = "# INPUT DIRS/PATHS\n"
    dependencies_text += "##binary or file paths\n"
    dependencies_text += 'in_ene_ana_lib = "' + in_ene_ana_lib + '"\n'
    dependencies_text += 'gromosPP_path = "' + gromosPP_path + '"\n'
    dependencies_text += "\n"
    dependencies_text += "##dependencies\n"
    dependencies_text += 'in_simulation_name = "' + in_simulation_name + '"\n'
    dependencies_text += 'in_folder = "' + in_folder + '"\n'
    dependencies_text += 'in_topology_path = "' + in_topology_path + '"\n'
    dependencies_text += 'out_dir = "' + out_dir + '"\n'
    dependencies_text += "\n"
    dependencies_text += "##controls\n"
    dependencies_text += "control_dict = {\n"
    dependencies_text += '"concat": {"do": True,\n'
    dependencies_text += '\t"sub": {\n'
    dependencies_text += '\t\t"tre": True,\n'
    dependencies_text += '\t\t"cnf": True,\n'
    dependencies_text += '\t\t"trc": True,\n'
    dependencies_text += '\t\t"convert_trc_xtc": True\n'
    dependencies_text += "\t\t}\n"
    dependencies_text += "\t},\n"
    dependencies_text += '"ene_ana":{"do":True}\n'
    dependencies_text += "}\n"
    dependencies_text += "trc_take_every_step = 1\n"
    dependencies_text += "\n\n"

    command_text = (
        "ana.do(in_simulation_name=in_simulation_name, in_folder=in_folder, "
        "in_topology_path=in_topology_path, out_dir=out_dir,"
        "in_ene_ana_lib=in_ene_ana_lib, gromosPP_path=gromosPP_path)\n\n"
    )

    script_text = import_text + dependencies_text + command_text

    ana_file = open(out_script_path, "w")
    ana_file.write(script_text)
    ana_file.close()

    return out_script_path
