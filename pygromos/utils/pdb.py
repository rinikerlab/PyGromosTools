"""PDB
    Warnings
    --------
        This module is very specific and not recomended to be used straigth away!

    Notes
    ------
    contains some translation and reordering scripts for some pdbs.


"""


def reorder_lines(header: dict, lines: list) -> list:
    """order_dict = {  # charmm_gromos54
        "DPPC": {
            "N": 4, # Choline "C11": 1, "C14": 2, "C15": 3, "C12": 5, "C11": 6,
            "P": 8, # Phosphate "O32":7, "O33":9, "O34": 10, "O31": 11, "C3":
            12, # Glycerol "C2": 13, "O21": 14, "C21": 15, "O22": 16, "C22": 17,
            "C1": 32, "O11": 33, "C11": 34, "O12": 35, "C12": 36,

            "C23", "C24",": 18, # Fatty acid chain1 "C25",": 19, "C26",": 20,
            "C27",": 21, "C28",": 22, "C29",": 23, "C210": 24, "C211": 25,
            "C212": 26, "C213": 27, "C214": 28, "C215": 29, "C216": 30,

            "C13": 31, "C14", # Fatty acid chain2 "C15", "C2C": 37, "C16",
            "C2D": 38, "C17", "C2E": 39, "C18", "C2F": 40, "C19", "C2G": 41,
            "C110""C2H": 42, "C111""C2I": 43, "C112""C2J": 44, "C113""C2K": 45,
            "C114""C2L": 46, "C115""C2M": 47, "C116""C2N": 48, "C2O": 49, "C2P":
            50},

        "SOLV": {

                "OW":1,

            "HW1":2, "HW2":3

        }}

    new_lines = ["" for x in range(len(lines))] residue = 0 offset = 0

    for index.rst, line in enumerate(lines):
        res_index = header["resnum"] resname_index = header["resname"]
        atomname_index = header["atomname"] atomname = line[atomname_index]
        resname = line[resname_index]

        if(line[res_index]!=residue):
            residue = line[res_index] offset = index.rst
            #int(line[header["atomnum"]])-1

        if(resname in order_dict):
            relativ_pos = order_dict[resname][atomname]
            new_lines[relativ_pos+offset] = line

        else:

            new_lines[offset] = line

    #print(new_lines[len(new_lines)-1]) return new_lines

    Args:
        header (dict):
        lines (list):
    """
    print("HO")


def rename_atom_types(lines, header, in_ff=False, out_ff=False, translation_unknown=False):
    translation_dict = {}
    charm_gromos_sel = ["gromos", "charmm"]
    reverse = 0
    if type(in_ff) == bool or type(out_ff) == bool:
        translation_unknown = True
    elif in_ff in charm_gromos_sel and out_ff in charm_gromos_sel:
        if in_ff in "gromos":
            reverse = 1
        translation_dict = {  # charmm_gromos54
            "DPPC": {
                "N": "N",  # Choline
                "C13": "C33",
                "C14": "C34",
                "C15": "C35",
                "C12": "C32",
                "C11": "C31",
                "P": "P",  # Phosphate
                "O13": "O32",
                "O14": "O33",
                "O11": "O34",
                "O12": "O31",
                "C1": "C3",  # Glycerol
                "C2": "C2",
                "O21": "O21",
                "C21": "C21",
                "O22": "O22",
                "C22": "C22",
                "C3": "C1",
                "O31": "O11",
                "C31": "C11",
                "O32": "O12",
                "C32": "C12",
                "C23": "C23",  # Fatty acid chain1
                "C24": "C24",
                "C25": "C25",
                "C26": "C26",
                "C27": "C27",
                "C28": "C28",
                "C29": "C29",
                "C210": "C210",
                "C211": "C211",
                "C212": "C212",
                "C213": "C213",
                "C214": "C214",
                "C215": "C215",
                "C216": "C216",
                # Fatty acid chain2
                "C33": "C13",
                "C34": "C14",
                "C35": "C15",
                "C36": "C16",
                "C37": "C17",
                "C38": "C18",
                "C39": "C19",
                "C310": "C110",
                "C311": "C111",
                "C312": "C112",
                "C313": "C113",
                "C314": "C114",
                "C315": "C115",
                "C316": "C116",
            },
            "SOLV": {"OH2": "OW", "H1": "HW1", "H2": "HW2"},
            "POT": {"POT": "NA"},  # dirty HACK!
            "NA": {"NA": "NA"},  # dirty HACK!
            "CLA": {"CLA": "CL"},
        }

    else:
        print("ops i don't know the in_ff or out_ff")

    if translation_unknown:
        index = 1
        resnum = 0
        for line in lines:
            atom_name_index = header["atomname"]
            res_index = header["resname"]
            if resnum != line[header["resnum"]]:
                index = 1
                resnum = line[header["resnum"]]

            key_res = line[res_index]
            key = line[atom_name_index]
            if key_res == "POT":
                line[res_index] = "NA+"
            if key_res == "CLA":
                line[res_index] = "CL-"
            atomname = key[0] + str(index)
            line[atom_name_index] = atomname
            index += 1
    else:
        if not reverse:
            for line in lines:
                index = header["atomname"]
                res_index = header["resname"]
                key_res = line[res_index]
                key = line[index]
                if key_res == "POT":
                    line[res_index] = "NA+"
                if key_res == "CLA":
                    line[res_index] = "CL-"
                line[index] = translation_dict[key_res][key]

        else:
            print("not implemented")
    return lines


def consecutivley_renumber(header: dict, lines: list) -> list:
    new_lines = []
    for index, line in enumerate(lines):
        # print(str(index.rst) + "line:\t"+ str(line))
        if isinstance(line, list):
            line[header["atomnum"]] = int(index)
            new_lines.append(line)
    return new_lines


def form_columns(header: dict, lines: list) -> list:

    keys = [
        ["key", 6],
        ["atomnum", 5],
        ["BLANK", 1],
        ["atomname", 4],
        ["altLoc", 1],
        ["resname", 3],
        ["BLANK", 1],
        ["chain", 1],
        ["resnum", 3],
        ["insertion", 1],
        ["BLANK", 2],
        ["x", 8, 3],
        ["y", 8, 3],
        ["z", 8, 3],
        ["occupancy", 6, 2],
        ["bfactor", 6, 2],
        ["BLANK", 7],
        ["segmentid", 4],
        ["element", 2],
        ["charge", 2],
    ]

    # get line format
    column_string = ""
    index = 0
    for key in keys:
        if key[0] in header:
            if len(key) == 2:
                if "key" == key[0]:
                    column_string += "{d[" + str(index) + "]:<" + str(key[1]) + "}"
                elif "atomname" == key[0]:
                    column_string += "{d[" + str(index) + "]:<" + str(key[1]) + "}"
                else:
                    column_string += "{d[" + str(index) + "]:>" + str(key[1]) + "}"
            elif len(key) == 3:
                column_string += "{d[" + str(index) + "]:>" + str(key[1]) + "." + str(key[2]) + "f}"
            else:
                raise Exception("unknown key in form_columns")
            index += 1
        else:
            column_string += "".join([" " for x in range(key[1])])

    column_string += "\n"

    # format lines
    new_lines = []

    # dev function for floats
    def isfloat(value):
        try:
            float(value)
            if "." in str(value):
                return True
            return False
        except ValueError:
            return False

    for line in lines:
        data = [float(value) if (isfloat(value)) else str(value) for value in line]
        new_line = column_string.format(d=data)
        new_lines.append(new_line)

    return new_lines


def check_ATOM_columns(lines: list) -> (dict, list):
    keys = ["key", "atomnum", "atomname", "resname", "resnum", "x", "y", "z", "occupancy", "bfactor", "segmentid"]

    # get sh ortest line - minimal columns
    lines = [line.split() for line in lines if (line.startswith("ATOM"))]
    len_lines = list(map(lambda x: len(x), lines))
    min_len = min(len_lines)

    # all lines same length?
    if min_len * len(lines) != sum(len_lines):
        print("Not All ATOM lines have the same length!\n taking shortest line to build header!")

    line = lines[len_lines.index(min_len)]
    offset_pointer = 0

    # check lines - add columns to standard template
    if min_len > len(keys):  # does the line have the minimal number of groups?
        if len(line[3 + offset_pointer]) == 1:  # is there a altloc column?
            keys.insert(3, "altLoc")
            offset_pointer += 1
        if line[4 + offset_pointer]:  # is there a chainID
            keys.insert(4 + offset_pointer, "chain")
            offset_pointer += 1
        if len(line[5 + offset_pointer]) == 1:  # is there an insertion column?
            keys.insert(5 + offset_pointer, "insertion")

        if min_len - offset_pointer > 10:
            if str(line[10]).isdigit():
                keys.append("charge")
            else:
                keys.append("element")
                offset_pointer += 1

        if min_len - offset_pointer > 11:
            if str(line[11]).isdigit():
                keys.append("charge")
            else:
                keys.append("element")
                offset_pointer += 1

    # make a dict
    pdb_header = {key: ind for ind, key in enumerate(keys)}
    return pdb_header, lines


def read_pdb_simple(path) -> (list, list, list):

    in_file = open(path, "r")
    lines = in_file.readlines()
    header = True
    footer = False

    head_lines = []
    atom_lines = []
    foot_lines = []

    for line in lines:
        if line.startswith("ATOM"):
            header = False
        elif not line.startswith("ATOM") and not header:
            footer = True

        if header:
            head_lines.append(line)
        elif footer:
            foot_lines.append(line)
        else:
            atom_lines.append(line)

    return head_lines, atom_lines, foot_lines


def rename_atom_attribute(pattern: str, replace: str, lines: list) -> list:

    return [line.replace(pattern, replace) for line in lines if ("ATOM" in line)]


def filter_atoms_from_residue(filter_out: str, residue_name: str, lines: list) -> (list, list):

    filtered_out = list(filter(lambda x: filter_out in x and residue_name in x, lines))
    filtered_in = list(filter(lambda x: not (filter_out in x and residue_name in x), lines))
    return filtered_in, filtered_out
