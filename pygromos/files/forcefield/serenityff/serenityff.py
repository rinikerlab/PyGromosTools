import collections
from simtk import unit
from rdkit import Chem

from pygromos.data import topology_templates
from pygromos.files.forcefield.openff.openff import OpenFF
from pygromos.files.forcefield._generic_force_field import _generic_force_field
from pygromos.files.forcefield.serenityff.serenityff_data import serenityff_C12, serenityff_C6
from pygromos.files.topology.top import Top


class SerenityFF(_generic_force_field):
    def __init__(
        self, name: str = "serenity", path_to_files: str = None, auto_import: bool = True, verbose: bool = False
    ):
        super().__init__(name, path_to_files, auto_import, verbose)
        if auto_import:
            self.auto_import_ff()
        self.off = OpenFF(name, path_to_files, auto_import, verbose)

    def auto_import_ff(self):
        self.top = Top(in_value=topology_templates.topology_template_dir + "/blank_template+spc.top")
        self.develop = False
        self.C12_input = {}
        self.partial_charges = collections.defaultdict(float)
        self.serenityFFelements = ["H", "C", "N", "O", "F", "S", "Br", "I"]
        self.C6_pattern = collections.defaultdict(list)
        self.C12_pattern = collections.defaultdict(list)

    def create_top(
        self,
        mol: str,
        in_top: Top = None,
        C12_input={"H": 0.0, "C": 0.0},
        partial_charges=collections.defaultdict(float),
    ):
        self.top = in_top
        self.off.gromosTop = in_top
        self.off._init_mol_for_convert(mol)
        self.off.convertResname()
        self.off.convertBonds()
        self.off.convertAngles()
        self.off.convertTosions()
        self.off.convertImproper()
        self.off.convert_other_stuff()
        self.create_serenityff_nonBonded(C12_input=C12_input, partial_charges=partial_charges)
        self.top = self.off.gromosTop
        return self.top

    def read_pattern(self, C6orC12: str = "C6"):
        for element in self.serenityFFelements:
            if C6orC12 == "C6":
                folder = serenityff_C6
            else:
                folder = serenityff_C12
            try:
                infile = open(folder + element + ".dat", "r")
                for line in infile:
                    content = line.strip().split()
                    if C6orC12 == "C6":
                        self.C6_pattern[element].append(content)
                    else:
                        self.C12_pattern[element].append(content)
                infile.close()
            except IOError:
                raise Exception("WIP")

    def _pattern_matching_for_one_element(self, element: str = "H") -> dict():
        # TODO: add C12 support

        return_dict = collections.defaultdict(list)
        for pattern in reversed(self.C6_pattern[element]):
            # create pattern and init some variables
            mol_pattern = Chem.MolFromSmarts(pattern[0])
            idx = 0
            idx_in_rdkmol = 0

            # find the atom for which the pattern was made
            for atom in mol_pattern.GetAtoms():
                if atom.GetAtomMapNum() == 1:
                    idx = atom.GetIdx()

            # get all matches

            matches = self.mol.GetSubstructMatches(mol_pattern, uniquify=False)
            if len(matches) >= 1:
                for match in matches:
                    idx_in_rdkmol = match[idx]
                    return_dict[idx_in_rdkmol] = [element + str(pattern[1]), pattern[2]]
        return return_dict

    def get_LJ_parameters(self) -> dict:
        return_dict = collections.defaultdict(list)
        if self.develop:
            # get all elements in self.mol
            contained_elements_set = set()
            for atom in self.mol.GetAtoms():
                element = atom.GetSymbol()
                contained_elements_set.add(element)

            # pattern matching for all elements contained in self.mol
            for element in contained_elements_set:
                return_dict.update(self._pattern_matching_for_one_element(element=element))
        else:
            raise NotImplementedError("WIP")
        return return_dict

    def create_serenityff_nonBonded(
        self, C12_input={"H": 0.0, "C": 0.0}, partial_charges=collections.defaultdict(float)
    ):
        if self.develop:
            self.read_pattern()
        c6dict = self.get_LJ_parameters()
        if self.develop:
            self.off.createVdWexclusionList()
            moleculeItr = 1
            for molecule in self.off.molecule_force_list:
                panm_dict = collections.defaultdict(int)
                for key in molecule["vdW"]:
                    force = molecule["vdW"][key]
                    ATNM = int(key[0]) + 1
                    MRES = moleculeItr
                    element_symbol = self.mol.GetAtomWithIdx(int(key[0])).GetSymbol()
                    panm_dict[element_symbol] += 1
                    PANM = element_symbol + str(panm_dict[element_symbol])
                    IAC = 0  # will not be used if we use automatic
                    MASS = self.off.openFFmolecule.atoms[int(key[0])].mass.value_in_unit(unit.dalton)
                    CG = 0
                    if self.develop:
                        CG = partial_charges[int(key[0])]
                    if ATNM == self.mol.GetNumAtoms():
                        CGC = 1
                    else:
                        CGC = 0
                    if str(key[0]) in self.off.exclusionList13:
                        openFFexList13 = list(self.off.exclusionList13[str(key[0])])
                        INE = [int(x) + 1 for x in openFFexList13]
                    else:
                        INE = list()
                    if str(key[0]) in self.off.exclusionList14:
                        openFFexList14 = list(self.off.exclusionList14[str(key[0])])
                        INE14 = [int(x) + 1 for x in openFFexList14]
                    else:
                        INE14 = list()
                    epsilon = float(force.epsilon.value_in_unit(unit.kilojoule_per_mole))
                    rmin = 2 * force.rmin_half.value_in_unit(unit.nanometer)
                    C6 = float(c6dict[key[0]][1])  # ** 2
                    CS6 = 0.5 * C6
                    C12 = epsilon * (rmin**12)
                    if self.develop:
                        C12 = C12_input[(c6dict[key[0]][0])]
                    CS12 = 0.5 * C12
                    IACname = c6dict[key[0]][0]
                    self.off.gromosTop.add_new_atom(
                        ATNM=ATNM,
                        MRES=MRES,
                        PANM=PANM,
                        IAC=IAC,
                        MASS=MASS,
                        CG=CG,
                        CGC=CGC,
                        INE=INE,
                        INE14=INE14,
                        C6=C6,
                        C12=C12,
                        CS6=CS6,
                        CS12=CS12,
                        IACname=IACname,
                    )
                moleculeItr += 1
        else:
            raise NotImplementedError("WIP")
