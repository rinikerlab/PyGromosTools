"""
File:            serenityff system folder managment
Warnings: this class is WIP!
TODO:REWORK
Description:
    This is a super class for the collection of all files used inside a serenityff project
Author: Marc Lehner
"""

# general imports
import collections, importlib
from simtk import unit

# import math

# rdkit imports
from rdkit import Chem

# pygromos imports
from pygromos.files.topology.top import Top
from pygromos.files.gromos_system.ff.serenityff.serenityff_data import serenityff_C12, serenityff_C6
from pygromos.files.gromos_system.ff.forcefield_system import forcefield_system

if importlib.util.find_spec("openff") == None:
    raise ImportError("SerenityFF is not enabled without openFF toolkit package! Please install openFF toolkit.")
else:
    from openff.toolkit.topology import Molecule
    from pygromos.files.gromos_system.ff.openforcefield2gromos import openforcefield2gromos


class serenityff:
    def __init__(
        self,
        mol: Chem.rdchem.Mol,
        forcefield: forcefield_system or str = None,
        top: Top = None,
        mol_name=None,
        develop=False,
    ):
        self.serenityFFelements = ["H", "C", "N", "O", "F", "S", "Br", "I"]
        self.C6_pattern = collections.defaultdict(list)
        self.C12_pattern = collections.defaultdict(list)
        self.mol = mol
        self.offmol = Molecule.from_rdkit(self.mol)

        if isinstance(forcefield, forcefield_system):
            self.top = forcefield.top
            if top != None:
                self.top = top
            self.offmol.name = forcefield.mol_name
            self.develop = forcefield.develop
            self.off = forcefield.off
            self.off2g = openforcefield2gromos(openFFmolecule=self.offmol, gromosTop=self.top, forcefield=forcefield)
        else:
            self.top = top
            if mol_name != None:
                self.offmol.name = mol_name
            self.off2g = openforcefield2gromos(openFFmolecule=self.offmol, gromosTop=self.top, forcefield=forcefield)
            self.develop = develop

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
            except:
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
            self.off2g.createVdWexclusionList()
            moleculeItr = 1
            for molecule in self.off2g.molecule_force_list:
                panm_dict = collections.defaultdict(int)
                for key in molecule["vdW"]:
                    force = molecule["vdW"][key]
                    ATNM = int(key[0]) + 1
                    MRES = moleculeItr
                    element_symbol = self.mol.GetAtomWithIdx(int(key[0])).GetSymbol()
                    panm_dict[element_symbol] += 1
                    PANM = element_symbol + str(panm_dict[element_symbol])
                    IAC = 0  # will not be used if we use automatic
                    MASS = self.off2g.openFFmolecule.atoms[int(key[0])].mass.value_in_unit(unit.dalton)
                    CG = 0
                    if self.develop:
                        CG = partial_charges[int(key[0])]
                    if ATNM == self.mol.GetNumAtoms():
                        CGC = 1
                    else:
                        CGC = 0
                    if str(key[0]) in self.off2g.exclusionList13:
                        openFFexList13 = list(self.off2g.exclusionList13[str(key[0])])
                        INE = [int(x) + 1 for x in openFFexList13]
                    else:
                        INE = list()
                    if str(key[0]) in self.off2g.exclusionList14:
                        openFFexList14 = list(self.off2g.exclusionList14[str(key[0])])
                        INE14 = [int(x) + 1 for x in openFFexList14]
                    else:
                        INE14 = list()
                    epsilon = float(force.epsilon.value_in_unit(unit.kilojoule_per_mole))
                    rmin = 2 * force.rmin_half.value_in_unit(unit.nanometer)
                    C6 = (float(c6dict[key[0]][1])) ** 2
                    CS6 = 0.5 * C6
                    C12 = epsilon * (rmin**12)
                    if self.develop:
                        C12 = C12_input[(c6dict[key[0]][0])]
                    CS12 = 0.5 * C12
                    IACname = c6dict[key[0]][0]
                    self.off2g.gromosTop.add_new_atom(
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

    def create_top(self, C12_input={"H": 0.0, "C": 0.0}, partial_charges=collections.defaultdict(float)):
        self.off2g.convertResname()
        self.off2g.convertBonds()
        self.off2g.convertAngles()
        self.off2g.convertTosions()
        self.off2g.convertImproper()
        self.off2g.convert_other_stuff()
        self.create_serenityff_nonBonded(C12_input=C12_input, partial_charges=partial_charges)
        self.top = self.off2g.gromosTop
