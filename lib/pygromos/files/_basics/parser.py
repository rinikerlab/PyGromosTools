"""
FUNCTIONLIB: parser - gromos
Description:

Author: Benjamin Schroeder
"""
import re
import warnings
import pandas as pd
from itertools import chain, takewhile

from pygromos.files import blocks
from pygromos.files.blocks import topology_blocks as tb
from pygromos.files.blocks import pertubation_blocks as pb
from pygromos.files.blocks._general_blocks import _generic_gromos_block
from pygromos.utils.typing import List, Dict


# translation dicts
imd_field_translation_dict = {
    "FORCE": {
        "BONDS": ["bonds"],
        "ANGLES": ["angles"],
        "IMPROPER": ["imp.", "improper"],
        "DIHEDRAL": ["dihe", "dihedral"],
        "ELECTROSTATIC": ["charge", "electrostatic"],
        "VDW": ["nonbonded", "vdW"],
    },
    "PAIRLIST": {"ALGORITHM": ["algorithm"]},
}


# private functions
def _gather_bracket_key(keys: List[str]) -> List[str]:
    """_gather_bracket_key
        private
        this function is needed for bracket fields for exapmle in IMDS NRE[1...N]   makes these field recognizable
        after this function split can be used.

    Parameters
    ----------
    keys :  List[str]
        field names

    Returns
    -------
    List[str]
        field names without whitespaces
    """
    gathered_keys = []
    collecting = False
    collect = []
    tmp_key = False

    try:
        for x in keys:  # gather brackets from comments into one subkey Eir(, ... ,..xs)]
            if "(" in x and (")" not in x) and not collecting:  # if a bracket is not closed in a key
                if x.startswith("("):  # if ( was seperated from root term
                    tmp_key = gathered_keys.pop()
                collecting = True
                collect.append(x)
            elif collecting and ")" in x:  # bracket comment finished
                collecting = False
                collect.append(x)

                key_bracket = " ".join(collect)
                if tmp_key:
                    key_bracket = "".join(tmp_key + key_bracket)
                gathered_keys.append(key_bracket)
                collect = []
            elif collecting:
                collect.append(x)
            else:
                gathered_keys.append(x)
        return gathered_keys
    except IndexError as err:
        print('Index Error - While gathering bracket, unknown key: "' + str(x) + '" (' + str(err) + ")")


def _update_line_information(subkeys: List[str], values: List[any], sub_content: Dict[str, any]):
    """_update_line_information
        private
        Subfunction of read_gromos_top - read a gromos block - read update_line_information:


    Parameters
    ----------
    subkeys :  List[str]
     are the field names inside a block
    values :    List[any]
        values of these fields

    sub_content :   Dict[str, any]
        gets updated by key:value - scheme

    Returns
    -------
    None
    """
    # call by ref
    if len(subkeys) == len(values) and len(subkeys) > 1:
        sub_content.update({key: values[ind] for ind, key in enumerate(subkeys)})
    elif len(subkeys) == 1 and len(values) == 1:
        sub_content.update({subkeys[0]: values[0]})
    elif len(values) == 1:
        sub_content.update({" ".join(subkeys): values[0]})
    elif len(values) > 0:
        sub_content.update({" ".join(subkeys): values})


def _read_gromos_subblock(block: List[str]) -> Dict[str, any]:
    """Subfunction of read_gromos_imd
        private
        read a gromos block

    Parameters
    ----------
    block : List[str]
        block to be parsed

    Returns
    -------
    Dict[str, any]
        parsed content block
    """
    sub_content = {}
    subkeys = []
    values = []
    first_key = True
    first_value = True
    extend_first = True
    temp = []

    for line in block:
        if line.strip().startswith("#"):
            temp = line.lstrip("#").strip()
            if not temp.isdigit():
                # print("Is NOT number: ",temp,"<")
                if not first_key:
                    _update_line_information(subkeys=subkeys, values=values, sub_content=sub_content)
                subkeys = line.replace("#", "").split()
                values = []
                first_value = True
                first_key = False
                extend_first = True

        elif not first_key:  # no fields?
            tmp_values = line.split()
            if len(tmp_values) > 0:  # collect all values
                if first_value:
                    values = tmp_values
                    first_value = False
                elif extend_first:
                    values = [values, tmp_values]
                    extend_first = False
                else:
                    values.append(tmp_values)
        else:
            values.append(line)
            first_value = False

    if not first_key:  # block is finished, do last update of values
        _update_line_information(subkeys=subkeys, values=values, sub_content=sub_content)
    elif not first_value and first_key:  # no fields? than give block value here
        sub_content = values

    return sub_content


def _read_gromos_block(lines: List[str]) -> Dict[str, str]:
    """_read_gromos_block
        private
            seperate gromos plocks into a dictionary


    Parameters
    ----------
    lines : List[str]
        sepeartes an input file into gromos blocks.
    Returns
    -------
    Dict[str, str]
        contains all blocks and the block contents
    """
    first_key = True
    blocks = {}
    key = ""
    block_lines = []
    for line in lines:
        if first_key and not line.startswith("#"):
            key = line.strip()
            first_key = False

        elif line.strip().startswith("END"):
            blocks.update({key: block_lines})
            first_key = True
            block_lines = []
        else:
            block_lines.append(line)
    return blocks


# FUNCTION - parser
def read_general_gromos_file(path: str) -> Dict:
    data = {}

    first_key = True
    key = ""
    block = []

    with open(path, "r") as infile:
        for line in infile:
            if first_key:
                if line.startswith("#"):
                    continue
                key = line.strip().upper()
                first_key = False
            elif "END" in line:
                data.update({key: block})
                block = []
                first_key = True
            else:
                block.append(line)
    infile.close()
    return data


# TOP
def read_disres(in_path: str) -> Dict:
    """read_disres
        This function can read distance restraints

    Parameters
    ----------
    in_path :   input path

    Returns
    -------
    Dict
    """
    # specific parser for subblocks of disres files
    def _read_disres_subblock(blocks):
        result_data = {}
        for block in blocks:

            if block == "TITLE":
                block_obj = getattr(tb, block)(content=blocks[block])
                result_data.update({block: block_obj})

            elif block == "DISTANCERESSPEC":
                block_obj = getattr(tb, block)(content=blocks[block])
                result_data.update({block: block_obj})
            else:
                raise IOError(
                    "DISRES parser does not know block: " + str(block) + "\n with content: " + "\n".join(blocks[block])
                )
        return result_data

    with open(in_path, "r") as infile:
        lines = infile.readlines()
        # parse the coarse gromos_block structure
        blocks = _read_gromos_block(lines)
        # parse data of the _blocks
        data = _read_disres_subblock(blocks)
    return data


def read_ptp(in_path: str) -> Dict:
    """read_disres
        This function can read distance restraints

    Parameters
    ----------
    in_path :   input path

    Returns
    -------
    Dict
    """
    # specific parser for subblocks of disres files
    def _read_ptp_subblock(blocks):
        result_data = {}
        for block in blocks:
            if hasattr(pb, block):
                atom_block = getattr(pb, block)(content=blocks[block])
                result_data.update({block: atom_block})
            else:
                raise IOError(
                    "PTP parser does not know block: " + str(block) + "\n with content: " + "\n".join(blocks[block])
                )
        return result_data

    with open(in_path, "r") as infile:
        lines = infile.readlines()
        # parse the coarse gromos_block structure
        blocks = _read_gromos_block(lines)
        # parse data of the _blocks
        data = _read_ptp_subblock(blocks)
    return data


# COORDS
def read_cnf(in_file_path: str, verbose: bool = False) -> Dict[str, str]:
    """read_cnf
        This function is reading in cnf files from gromos and translates them to a dict if possible.
            atomP = namedtuple("AtomP", "resID resName atomType atomID xp yp zp")
            atomV = namedtuple("AtomV", "resID resName atomType atomID xv yv zv")
            lattice_shift = namedtuple("Latticeshifts", "atomID x y z")
            box = namedtuple("GenBox", "pbc length angles euler origin")

    Parameters
    ----------
    in_file_path :  str
        input path to cnf file.
    verbose :

    Returns
    -------

    """

    known_blocks = {
        "TITLE",
        "TIMESTEP",
        "POSITION",
        "LATTICESHIFTS",
        "VELOCITY",
        "POSRESSPEC",
        "REFPOSITION",
        "GENBOX",
        "PERTDATA",
        "STOCHINT",
    }
    file = open(in_file_path, "r")

    data = {}
    first_key = False
    block = "none"
    subblock = []
    # Translate string
    for line in file.readlines():
        if not first_key and "\#" not in line and len(line.split()) == 1:  # noqa: W605
            first_key = True
            block = line.strip().split()[0]
        elif "END" == line.strip():
            if block in known_blocks:

                atom_block = getattr(blocks.coords, block)(content=subblock)
                data.update({block: atom_block})
            else:
                print("WARNING!: Don't know block: " + block + ". read it simple.")
                data.update({block: _generic_gromos_block(name=block, used=True, content=subblock)})
            first_key = False
            block = "none"
            subblock = []
        elif block == "none":
            if line.startswith("#"):
                continue
            else:
                raise IOError(
                    "coord file has inconsistent gromos Block structure!\n This line is out of block: \n" + line
                )
        else:
            if "#" not in line:
                subblock.append(line)
    return data


# REEDS
def read_repdat(path: str, Vj_header=False) -> (Dict, List[blocks.repdat.replica_stat]):
    """
    Careful old description!
    Parameters
    ----------
    path :
    Vj_header :

    Returns
    -------
        result_dict{"system":..., "header":..., "data":....}
    """

    file = open(path, "r")

    system_options = {}
    data_header = []
    data_values = []
    data_values_start = False

    first = True
    new_repdat = False
    for line in file.readlines():
        if first and "#======================" in line:
            # HACK! TODO: make this hackless... currently due to the multiple gromos verions needed.
            new_repdat = True
            break

        first = False
        if line.startswith("#"):
            if Vj_header:
                data_header = list(
                    map(
                        lambda x: x.split("(")[0],
                        filter(lambda x: not str(x).startswith("Vj"), line.replace("#", "").split()),
                    )
                )
            else:
                data_header = list(line.replace("#", ""))
            data_values_start = True
        elif data_values_start:
            data_values.append(line.split())
        else:
            if line.strip().startswith("E"):
                fields = line.split()
                fields = [fields[0], fields[1:]]
            else:
                fields = line.strip().split("\t")
            if len(fields) > 1:
                key = (
                    fields[0]
                    .replace(":", "")
                    .strip()
                    .replace("Number of ", "")
                    .replace(" ", "_")
                    .replace(",_num", "_")
                    .replace("=", "")
                    .replace("_(RE-EDS)", "")
                )
                value = fields[1]
                system_options.update({key: value})

    if new_repdat:  # TODO: NEW PARSER for actual gromos version!
        file.close()

        eir = {}
        system_options_dict = {}
        eir_ind = 1
        with open(path, "r") as repdat:
            header = takewhile(lambda s: s.startswith("#"), repdat)
            for line in header:
                fields = line.replace("(RE-EDS)", "").replace("#", "").split()
                if len(fields) > 1:
                    if fields[0] == "T":
                        system_options_dict.update({"T": float(fields[1])})
                    elif fields[0] == "s":
                        system_options_dict.update({"s": list(map(float, fields[1:]))})
                    elif "eir" in fields[0] or re.search("E[0-9]*R\(s\)", fields[0]):  # noqa: W605
                        if "eir" in fields[0]:
                            eir.update({eir_ind: list(map(float, fields[7:]))})
                        else:
                            eir.update({eir_ind: list(map(float, fields[1:]))})
                        eir_ind += 1

            system_options_dict.update({"state_eir": eir})
        system_options = blocks.repdat.repex_system(
            s=system_options_dict["s"], T=system_options_dict["T"], state_eir=system_options_dict["state_eir"]
        )

        # read in data
        df = pd.read_csv(path, comment="#", delim_whitespace=True)

        # To do: messy juggling to match old fun
        df.rename(columns={"exch": "s"}, inplace=True)

        # orig  column order- to be maintained, just replace V columns by state_potential dict
        data_header = df.columns
        V_cols = [col for col in df.columns if (col.startswith("V"))]
        rest_cols = [col for col in df.columns if (col not in V_cols)]

        state_potential = list(df[V_cols].T.to_dict().values())
        dict_list = list(df[rest_cols].T.to_dict().values())
        [dict_vals.update({"state_potentials": state_V}) for state_V, dict_vals in zip(state_potential, dict_list)]

        new_header = [x for x in data_header if (not x.startswith("V"))] + ["state_potentials"]

        df = pd.DataFrame(dict_list)
        df = df[new_header]
        df = df.astype({"ID": "int32", "partner": "int32", "run": "int32", "s": "int32"})
        return system_options, df

    else:
        # cleanup system_options
        new_dict = {
            i: system_options[i]
            for i in system_options
            if not (i.startswith("eir") or i.strip().startswith("E")) and (i == "s" or i == "T")
        }
        # print(new_dict.keys())
        if "T" in new_dict:
            new_dict.update({"T": float(new_dict["T"])})
        if "s" in new_dict:
            new_dict.update({"s": list(map(float, new_dict["s"].split()))})
        if "lambda" in new_dict:
            new_dict.update({"s": list(map(float, new_dict["lambda"].split()))})

        eir = {}
        eir_state = 1
        for i in system_options:
            if i.startswith("eir"):
                eir.update({eir_state: list(map(float, system_options[i].split()))})
                eir_state += 1

            elif i.strip().startswith("E" + str(eir_state) + "R"):
                eir.update({eir_state: list(map(float, system_options[i]))})
                eir_state += 1

        new_dict.update({"state_eir": eir})
        system_options = blocks.repdat.repex_system(s=new_dict["s"], T=new_dict["T"], state_eir=new_dict["state_eir"])

        # data
        # clean up data
        head_num = len(data_header)
        values = []
        for line in data_values:
            pots = {data_header[i]: float(line[i]) for i in range(head_num) if (data_header[i].startswith("V"))}
            value = {data_header[i]: float(line[i]) for i in range(head_num) if not (data_header[i].startswith("V"))}
            value.update({"state_potentials": pots})
            values.append(value)
        new_header = [x for x in data_header if (not x.startswith("V"))] + ["state_potentials"]
        df = pd.DataFrame(values)
        df = df[new_header]
        return system_options, df


# IMD
def read_imd(in_file_path: str) -> Dict:
    """imd parser

    Parameters
    ----------
    in_file_path :  str
        input file path

    Returns
    -------
    Dict
        dict - containing all _blocks as keys and possible fields as subkeys
    """

    def translate_subcontent_keys(key: str, sub_content: dict, translation_dict: dict) -> dict:
        if key in translation_dict:
            # print(key)
            # print(list(sub_content))
            # print(translation_dict[key])
            translated_sub_content = {}
            sub_trans_dict = translation_dict[key]
            for sub_key in sub_content:  # assign new standardized keys fot the block fields
                tmp_sub_key = None
                for x in sub_trans_dict:
                    if sub_key in sub_trans_dict[x]:
                        tmp_sub_key = x
                        break
                if tmp_sub_key is not None:
                    # print("replace: "+sub_key+" with "+tmp_sub_key)
                    translated_sub_content.update({tmp_sub_key: sub_content[sub_key]})
                else:
                    # print("skipping: "+sub_key)
                    translated_sub_content.update({sub_key: sub_content[sub_key]})
            return translated_sub_content
        else:
            return sub_content

    def read_gromos_imd_block(block) -> dict:
        """read_gromos_imd_block
            Subfunction of read_gromos_imd - read a gromos block:

        Parameters
        ----------
        block :

        Returns
        -------
        Dict
            dict containing data
        """

        sub_content = {}
        subkeys = []
        values = []
        first_key = True
        first_value = True
        extend_first = True

        for line in block:
            if line.strip().startswith("#"):  # comment line defining field => keys for fields
                if not first_key and len(values) > 0:
                    update_line_information(subkeys=subkeys, values=values, sub_content=sub_content)
                subkeys = line.replace("#", "").split()
                values = []
                first_value = True
                first_key = False
                extend_first = True

            elif not first_key:  # no fields?
                tmp_values = line.split()
                if len(tmp_values) > 0:  # collect all values
                    if first_value:
                        values = tmp_values
                        first_value = False
                    elif extend_first:
                        values = [values, tmp_values]
                        extend_first = False
                    else:
                        values.append(tmp_values)
                        # print("TEST "+str(tmp_values))
            else:
                values.append(line)
                first_value = False

        if not first_key:  # block is finished, do last update of values
            update_line_information(subkeys=subkeys, values=values, sub_content=sub_content)
        elif not first_value and first_key:  # no fields? than give block value here
            sub_content = values

        return sub_content

    def update_line_information(subkeys, values, sub_content: dict):
        """update_line_information
            Subfunction of read_gromos_imd - read a gromos block - read update_line_information:

        Parameters
        ----------
        subkeys :
        values :
        sub_content :
            is updated during process
        Returns
        -------

        """

        # avoid brackets screwing the read in for keys collapses brackets to preceeding key
        subkeys = _gather_bracket_key(subkeys)
        subkeys = [str(subkey) for subkey in subkeys]

        # call by ref check how comment keys fit to data
        try:
            if len(subkeys) == len(values) and len(subkeys) > 1:  # more keys and equally ammount of values
                # Hack MULTIBATH block
                if "TEMP0" in subkeys[0]:
                    if isinstance(values[0], list):
                        values = list(chain.from_iterable(values))
                    sub_content.update({subkeys[0]: values[::2]})
                    sub_content.update({subkeys[1]: values[1::2]})
                # Hack MULTIBATH block
                elif "LAST" in subkeys[0]:
                    if isinstance(values[0], list):
                        values = list(chain.from_iterable(values))
                    sub_content.update({subkeys[0]: values[::3]})
                    sub_content.update({subkeys[1]: values[1::3]})
                    sub_content.update({subkeys[2]: values[2::3]})

                else:
                    sub_content.update({key: values[ind] for ind, key in enumerate(subkeys)})

            elif len(subkeys) == 1 and len(values) == 1:  # 1 keys and 1 value todo: really needed??? bschroed
                sub_content.update({subkeys[0]: values[0]})
            elif len(subkeys) == 1 and len(values) > 1:  # 1 key and more values
                sub_content.update({subkeys[0]: values})
            elif len(values) == 1:  # more keys and 1 value
                sub_content.update({" ".join(subkeys): values[0]})
            elif len(values) > 0:  # unrecognizable LINE! - more values than keys and more than one key
                # if("..." in subkeys):
                #    make_keys=len(values)-2
                #    tmp_subkeys=[subkeys[0]]+[subkeys[1].replace("1", "").replace("(", "").replace(")","")+"_"+str(i) for i in range(make_keys) ]+ [subkeys[len(subkeys)-1].replace("(", "_").replace(")","")]
                #    sub_content.update({co: values[ind] for ind,co in enumerate(tmp_subkeys)})

                if len(values) % len(subkeys) == 0:  # possible matrix as content?
                    n_vals = int(len(values) / len(subkeys))
                    # Hack MULTIBATH block
                    if "LAST" in subkeys[0]:
                        if isinstance(values[0], list):
                            values = list(chain.from_iterable(values))
                        sub_content.update({subkeys[0]: values[0::3]})
                        sub_content.update({subkeys[1]: values[1::3]})
                        sub_content.update({subkeys[2]: values[2::3]})
                    else:
                        for ind, subkey in enumerate(subkeys):
                            start = ind * n_vals
                            end = start + n_vals
                            sub_content.update({subkey: values[start:end]})
                # Hack PRESSURE block - SEMIANISOTROPIC has three parameters per parametername
                elif str(subkeys[0]) == "SEMIANISOTROPIC":
                    sub_content.update({subkeys[0]: values})
                # Hack FORCE block - energy groups do not have a parameter name per parameter
                elif str(subkeys[0]) == "NEGR":
                    if not ((len(values) - 1) == int((values[0]))):
                        print(
                            "Warning - In FORCE block: Number of energy groups does not match number of last atoms given!"
                        )
                    sub_content.update({subkeys[0]: values[0], "NRE": values[1:]})
                # Hack INNERLOOP block - four parameter names but only 3 parameters (NDEVG not implemented)
                elif str(subkeys[0]) == "NTILM":
                    sub_content.update({subkeys[0]: values[0]})
                    sub_content.update({subkeys[1]: values[1]})
                    sub_content.update({subkeys[2]: values[2]})
                else:
                    # Hack MULTIBATH block
                    if "TEMP0" in subkeys[0]:
                        if isinstance(values[0], list):
                            values = list(chain.from_iterable(values))
                        sub_content.update({subkeys[0]: values[::2]})
                        sub_content.update({subkeys[1]: values[1::2]})
                    elif "LAST" in subkeys[0]:
                        if isinstance(values[0], list):
                            values = list(chain.from_iterable(values))
                        sub_content.update({subkeys[0]: values[::3]})
                        sub_content.update({subkeys[1]: values[1::3]})
                        sub_content.update({subkeys[2]: values[2::3]})

                    else:
                        print("COULD not seperate:")
                        print(subkeys)
                        print(values)

                        sub_content.update({" ".join(subkeys): values})
        except IndexError:
            print("Errors while reading imd-file.")
        # TypeError:
        # if subkeys is None:
        #    print("Empty comment line found! Line ignored.")
        # else:
        #    print("Unknown TypeError.")

    data = {}

    first_key = True
    key = ""
    block = []

    with open(in_file_path, "r") as infile:
        for line in infile:
            if first_key:
                if line.startswith("#"):
                    continue
                key = line.strip().upper()
                first_key = False
            elif "END" in line:
                # print(key)
                # print(block)
                sub_content = read_gromos_imd_block(block=block)
                # print(sub_content)
                sub_content = translate_subcontent_keys(
                    key=key, sub_content=sub_content, translation_dict=imd_field_translation_dict
                )
                # print(sub_content)
                data.update({key: sub_content})
                block = []
                first_key = True
            else:
                block.append(line)
    infile.close()
    return data


# SIMPLIFIED Pasers- CRUDE - not necessarily correct
def read_simple_trx(in_file_path: str, every_step: int = 1, skip_steps: int = 0, verbose: bool = True) -> Dict:
    """
        Needs output for checknig.
        Simple tre_ read - reads only _blocks! does not seperate fields.
        rather ugly! but not much time!
    Parameters
    ----------
    in_file_path : str

    Returns
    -------
    Dict
        dict - containing all _blocks as keys (Title, Eneversion, num_trials {1:{}, 2:{} ....})
    """

    data = {"TIMESTEP": {}}

    key = ""
    block = []
    subblock = {}
    timestep = 0

    first_key = True
    issubblock = False
    firstTitleBlock = True
    firstEneVersionBlock = True
    skip = False
    with open(in_file_path, "r") as infile:
        current_step = 0
        if skip_steps > 0:
            step_filter = (
                lambda step_number: (step_number == 0 or step_number % every_step == 0) and step_number > skip_steps
            )
        else:
            step_filter = lambda step_number: step_number == 0 or step_number % every_step == 0

        for line in infile:
            if first_key:  # BlockHeader
                # Keep only one Title
                if "TITLE" in line and not firstTitleBlock:
                    skip = True

                # Keep only one eneversion
                elif "ENEVERSION" in line and not firstEneVersionBlock:
                    skip = True

                elif "TIMESTEP" in line and step_filter(current_step):
                    if timestep != 0:
                        data["TIMESTEP"][timestep].update(subblock)
                        subblock = {}

                    timestep += 1
                    vals = next(infile).split()
                    data["TIMESTEP"].update({timestep: {"steps": int(vals[0]), "time": float(vals[1])}})
                    skip = True

                elif "TITLE" in line:
                    key = line.strip()
                    firstTitleBlock = False

                elif "ENEVERSION" in line:
                    key = line.strip()
                    firstEneVersionBlock = False

                else:
                    issubblock = True
                    key = line.strip()
                first_key = False

            # BlockEnds
            elif "END" in line and issubblock:
                subblock.update({key: block})
                block = []
                issubblock = False
                first_key = True
                continue

            elif "END" in line and skip:
                skip = False
                first_key = True
                continue

            elif "END" in line:
                data.update({key: block})
                block = []
                first_key = True
                if key == "TIMESTEP":
                    current_step += 1
                continue

            # Contentlines
            elif skip:
                continue

            else:
                block.append(line)
                continue

    if not first_key:
        raise IOError("The Final gromos Block was not terminated!")

    # lastBlock
    if timestep != 0:
        data["TIMESTEP"][timestep].update(subblock)
    infile.close()
    return data


def read_tre(in_file_path: str, every_step: int = 1, verbose: bool = True) -> Dict:
    """
        Needs output for checknig.
        Simple tre_ read - reads only _blocks! does not seperate fields.
        rather ugly! but not much time!
    Parameters
    ----------
    in_file_path : str

    Returns
    -------
    Dict
        dict - containing all _blocks as keys (Title, Eneversion, num_trials {1:{}, 2:{} ....})
    """

    data = {"TIMESTEP": {}}

    timestep = 0
    key = ""
    block = []
    frameSubBlocks = {}

    # control bools
    isBlockTitle = True
    isFrameSubBlock = False
    isFirstTitleBlock = True
    isFirstEneVersionBlock = True
    skipBlockLines = False
    skipFrame = False
    doKeepFrameData = lambda step_number: step_number == 0 or step_number % every_step == 0

    with open(in_file_path, "r") as infile:
        for line in infile:
            if isBlockTitle:  # BlockHeader
                line = line.strip()
                if "TITLE" in line:
                    if isFirstTitleBlock:  # only keep one title block
                        key = line.strip()
                        isFirstTitleBlock = False
                    else:
                        skipBlockLines = True
                elif "ENEVERSION" in line:
                    if isFirstEneVersionBlock:  # only keep one title block
                        key = line.strip()
                        isFirstEneVersionBlock = False
                    else:
                        skipBlockLines = True
                elif "TIMESTEP" in line:
                    if doKeepFrameData(timestep):  # store data starts with
                        if timestep != 0:
                            data["TIMESTEP"].update(frameSubBlocks)
                        timestep += 1
                        vals = next(infile).split()
                        frameSubBlocks = {timestep: {"steps": int(vals[0]), "time": float(vals[1])}}
                        skipFrame = False
                    else:
                        timestep += 1
                        skipFrame = True
                    skipBlockLines = True
                elif "ENERGY03" in line or "VOLUMEPRESSURE03" in line:  # Known TRC subblocks
                    isFrameSubBlock = True
                    key = line.strip()
                else:
                    raise IOError("Found a not known tre block: " + line.strip())
                isBlockTitle = False
                continue
            elif "END" in line:  # BlockEnds
                if skipBlockLines or skipFrame:
                    skipBlockLines = False
                elif isFrameSubBlock:
                    frameSubBlocks[timestep].update({key: block})
                    block = []
                    isFrameSubBlock = False
                else:
                    data.update({key: block})
                    block = []
                isBlockTitle = True
                continue
            # Contentlines
            elif skipBlockLines or skipFrame:
                continue
            else:
                block.append(line)
                continue

        if not isBlockTitle:
            raise IOError("The Final gromos Block was not terminated by an END!")

        # lastBlock
        if timestep != 0 and len(frameSubBlocks) > 0:
            data["TIMESTEP"].update(frameSubBlocks)

    return data


def read_trc(in_file_path: str, every_step: int = 1, verbose: bool = True) -> Dict:
    """
        Needs output for checknig.
        Simple tre_ read - reads only _blocks! does not seperate fields.
        rather ugly! but not much time!
    Parameters
    ----------
    in_file_path : str

    Returns
    -------
    Dict
        dict - containing all _blocks as keys (Title, Eneversion, num_trials {1:{}, 2:{} ....})
    """

    data = {"TIMESTEP": {}}

    timestep = 0
    key = ""
    block = []
    frameSubBlocks = {}

    isBlockTitle = True
    isFrameSubBlock = False
    isFirstTitleBlock = True
    skipBlockLines = False
    skipFrame = False
    doKeepFrameData = lambda step_number: step_number == 0 or step_number % every_step == 0

    with open(in_file_path, "r") as infile:
        for line in infile:
            if isBlockTitle:  # BlockHeader
                if "TITLE" in line:
                    if isFirstTitleBlock:  # only keep one title block
                        key = line.strip()
                        isFirstTitleBlock = False
                    else:
                        skipBlockLines = True
                elif "TIMESTEP" in line:
                    if doKeepFrameData(timestep):  # store data starts with
                        if timestep != 0:
                            data["TIMESTEP"].update(frameSubBlocks)
                        timestep += 1
                        vals = next(infile).split()
                        frameSubBlocks = {timestep: {"steps": int(vals[0]), "time": float(vals[1])}}
                        skipFrame = False
                    else:
                        timestep += 1
                        skipFrame = True
                    skipBlockLines = True
                elif "POSITIONRED" in line or "GENBOX" in line:  # Known TRC subblocks
                    isFrameSubBlock = True
                    key = line.strip()
                else:
                    raise IOError("Found a not known trc block: " + line.strip())
                isBlockTitle = False
                continue

            elif "END" in line:  # BlockEnds
                if skipBlockLines or skipFrame:
                    skipBlockLines = False
                elif isFrameSubBlock:
                    frameSubBlocks[timestep].update({key: block})
                    block = []
                    isFrameSubBlock = False
                else:
                    data.update({key: block})
                    block = []
                isBlockTitle = True
            # Contentlines
            elif skipBlockLines or skipFrame:
                continue
            else:
                block.append(line)
                continue

        if not isBlockTitle:
            frame_keys = list(frameSubBlocks.keys())
            print(frame_keys)
            warnings.warn("total number of frames: " + str(len(frame_keys)))
            del frameSubBlocks[frame_keys[-1]]
            # raise IOError("The Final gromos Block was not terminated by an END!")

        # lastBlock
        if timestep != 0 and len(frameSubBlocks) > 0:
            data["TIMESTEP"].update(frameSubBlocks)

    return data


def read_gromos_csv(fileName, sep=" "):
    with open(fileName, "r") as infile:
        header = next(infile).replace("#", "").split()

        data_dict = {x: [] for x in header}
        for line in infile:
            words = list(filter(None, line.strip().split(sep)))
            for index, x in enumerate(header):
                try:
                    data_dict[x].append(float(words[index]))
                except IndexError:
                    print("WARNING: Line in %s contains too few columns: %s" % (fileName, line))
    return data_dict


def read_ene_ana_lib(in_path: str):

    with open(in_path, "r") as infile:
        lines = infile.readlines()
        block_lines = _read_gromos_block(lines)
        block_set = {}
        for block_title in block_lines:
            # print(block_title)
            if block_title in ["TITLE", "ENEVERSION"]:
                block_set.update({block_title: {"content": block_lines[block_title]}})

            elif block_title in ["ENERTRJ", "FRENERTRJ"]:
                comments = []
                block_dict = {}
                current_block = False
                for line in block_lines[block_title]:
                    if line.strip().startswith("#"):
                        comments.append(line)
                    elif line.strip().startswith("block"):
                        current_block = line.replace("block", "").strip()
                        block_dict.update({current_block: {}})
                    elif line.strip().startswith("subblock"):
                        if "subblock" not in block_dict[current_block]:
                            block_dict[current_block].update(
                                {"subblock": [line.replace("subblock", "").strip().split()]}
                            )
                        else:
                            block_dict[current_block]["subblock"].append(line.replace("subblock", "").strip().split())
                    elif line.strip().startswith("size"):
                        if "size" not in block_dict[current_block]:
                            block_dict[current_block].update({"size": [line.replace("size", "").strip()]})
                        else:
                            block_dict[current_block]["size"].append(line.replace("size", "").strip())
                    else:
                        raise ValueError("Ene_ana did not understand line in " + current_block + ": ", line)
                block_set.update({block_title: {"comments": "".join(comments), "_blocks": block_dict}})

            elif block_title == "VARIABLES":
                comments = []
                block_dict = {}
                for line in block_lines[block_title]:
                    if line.strip().startswith("#"):
                        comments.append(line)
                    else:
                        var_name, var_value = line.split("=")
                        block_dict.update({var_name.strip(): var_value.strip()})

                block_set.update({block_title: {"comments": "".join(comments), "variables": block_dict}})
            else:
                raise ValueError(
                    "FOUND A UKNOWN BLOCK! - PLEASE IMPLEMENT IT\n\t",
                    block_title,
                    "\n\t",
                    "\t".join(block_lines[block_title]),
                )
    return block_set
