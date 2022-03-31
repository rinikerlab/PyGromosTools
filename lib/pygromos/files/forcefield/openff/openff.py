import os
import glob
import importlib
import collections

from rdkit import Chem
from simtk import unit as u


if importlib.util.find_spec("openff") is None:
    raise ImportError(
        "openforcefield2gromos is not enabled without openFF toolkit package! Please install openFF toolkit."
    )
else:
    from openff.toolkit.topology import Molecule, Topology
    from openff.toolkit.typing.engines import smirnoff

from pygromos.files.coord.cnf import Cnf
from pygromos.files.forcefield._generic_force_field import _generic_force_field
from pygromos.files.topology.top import Top
from pygromos.data.ff import data_ff_SMIRNOFF
from pygromos.data import topology_templates


class OpenFF(_generic_force_field):
    def __init__(
        self, name: str = "openff", path_to_files: str = None, auto_import: bool = True, verbose: bool = False
    ):
        self.atomic_number_dict = collections.defaultdict(str)
        super().__init__(name, path_to_files=path_to_files, auto_import=auto_import, verbose=verbose)
        if auto_import:
            self.auto_import_ff()
        self.gromosTop = None

    def auto_import_ff(self):
        if self.path_to_files is not None:
            try:
                self.off = smirnoff.ForceField(self.path_to_files)
            except ImportError:
                raise ImportError("Could not import a OpenForceField from path: " + str(self.path_to_files))
        else:
            filelist = glob.glob(data_ff_SMIRNOFF + "/*.offxml")
            filelist.sort()
            filelist.reverse()
            for f in filelist:
                try:
                    self.off = smirnoff.ForceField(f)
                    self.path_to_files = f
                    break
                except ImportError:
                    pass
        print("Found off: " + str(self.path_to_files))

        # set atomic_number_dict
        self.atomic_number_dict[1] = "H"
        self.atomic_number_dict[2] = "He"
        self.atomic_number_dict[3] = "Li"
        self.atomic_number_dict[4] = "Be"
        self.atomic_number_dict[5] = "B"
        self.atomic_number_dict[6] = "C"
        self.atomic_number_dict[7] = "N"
        self.atomic_number_dict[8] = "O"
        self.atomic_number_dict[9] = "F"
        self.atomic_number_dict[10] = "Ne"
        self.atomic_number_dict[11] = "Na"
        self.atomic_number_dict[12] = "Mg"
        self.atomic_number_dict[13] = "Al"
        self.atomic_number_dict[14] = "Si"
        self.atomic_number_dict[15] = "P"
        self.atomic_number_dict[16] = "S"
        self.atomic_number_dict[17] = "Cl"
        self.atomic_number_dict[18] = "Ar"
        self.atomic_number_dict[19] = "K"
        self.atomic_number_dict[20] = "Ca"
        self.atomic_number_dict[35] = "Br"
        self.atomic_number_dict[53] = "I"

    def create_cnf(self, mol: str, in_cnf: Cnf = None, **kwargs) -> Cnf:
        return Cnf(in_value=mol)

    def create_top(
        self,
        mol: str,
        in_top: Top = None,
        **kwargs,
    ) -> Top:
        # prepare topology
        if in_top is not None:
            self.gromosTop = in_top
        else:
            self.gromosTop = Top(in_value=topology_templates.blank_topo_template)
            self.gromosTop._orig_file_path = os.getcwd()

        # create molecule
        self._init_mol_for_convert(mol=mol)

        # convert molecule
        self.convert()
        return self.gromosTop

    def _init_mol_for_convert(self, mol: str = None):
        if mol is None:
            raise ValueError("No molecule given!")
        elif isinstance(mol, Molecule):
            self.openFFmolecule = mol
        elif isinstance(mol, Chem.rdchem.Mol):
            self.openFFmolecule = Molecule.from_rdkit(mol)
        elif isinstance(mol, str):
            self.openFFmolecule = Molecule.from_smiles(mol)
        else:
            raise ValueError("mol is not a supported molecule type!")

        self.openFFTop = Topology.from_molecules(self.openFFmolecule)

        # create list of all forces
        self.molecule_force_list = []
        self.molecule_force_list = self.off.label_molecules(self.openFFTop)
        self.openmm_system = self.off.create_openmm_system(self.openFFTop)

        # 1-3 / 1-4 exclusion lists
        self.exclusionList13 = dict()
        self.exclusionList14 = dict()

    def convert(self):
        # print OpenFF warning in Title
        titleString = ""
        titleString += "\n\tname: " + self.openFFmolecule.name + "\t hill_formula: " + self.openFFmolecule.hill_formula
        titleString += (
            "\n\t"
            + 40 * "-"
            + "\n\t| created from OpenForceField topology |\n\t| use Amber Block for OpenFF topology! |\n\t"
            + 40 * "-"
            + "\n"
        )
        if hasattr(self.gromosTop, "TITLE"):
            self.gromosTop.TITLE.content += [titleString]
        else:
            self.gromosTop.add_block(blocktitle="TITLE", content=[titleString])
        # Do all the conversions
        self.convertResname()
        self.convertBonds()
        self.convertAngles()
        self.convertTosions()
        self.convertImproper()
        self.convertVdW()
        self.convert_other_stuff()

    def convertResname(self):
        if len(self.openFFmolecule.name) >= 1:
            self.gromosTop.add_new_resname(self.openFFmolecule.name)
        else:
            self.gromosTop.add_new_resname(self.openFFmolecule.hill_formula)

    def convertBonds(self):
        for molecule in self.molecule_force_list:
            for key in molecule["Bonds"]:
                force = molecule["Bonds"][key]
                # hQ = topology.atom(force[0]).atomic_number == 1 or topology.atom(force[1]).atomic_number == 1
                # hQ = not all([self.openFFTop.atom(x).atomic_number != 1 for x in key])  # noqa: F841
                atomI = key[0] + 1
                atomJ = key[1] + 1
                k = force.k.value_in_unit(u.kilojoule / (u.mole * u.nanometer**2))
                b0 = force.length.value_in_unit(u.nanometer)
                self.gromosTop.add_new_bond(k=k, b0=b0, atomI=atomI, atomJ=atomJ, includesH=False)  # hQ
        if not hasattr(self.gromosTop, "BONDSTRETCHTYPE"):
            self.gromosTop.add_block(blocktitle="BONDSTRETCHTYPE", content=[])
        if not hasattr(self.gromosTop, "BONDH"):
            self.gromosTop.add_block(blocktitle="BONDH", content=[])
        if not hasattr(self.gromosTop, "BOND"):
            self.gromosTop.add_block(blocktitle="BOND", content=[])

    def convertAngles(self):
        for molecule in self.molecule_force_list:
            for key in molecule["Angles"]:
                force = molecule["Angles"][key]
                # hQ = not all([self.openFFTop.atom(x).atomic_number != 1 for x in key])  # noqa: F841
                atomI = key[0] + 1
                atomJ = key[1] + 1
                atomK = key[2] + 1
                k = 0  # TODO: proper conversion to quartic
                kh = force.k.value_in_unit(u.kilojoule / (u.mole * u.degree**2))
                b0 = force.angle.value_in_unit(u.degree)
                self.gromosTop.add_new_angle(
                    k=k, kh=kh, b0=b0, atomI=atomI, atomJ=atomJ, atomK=atomK, includesH=False, convertToQuartic=True
                )  # hQ
        if not hasattr(self.gromosTop, "BONDANGLEBENDTYPE"):
            self.gromosTop.add_block(blocktitle="BONDANGLEBENDTYPE", content=[])
        if not hasattr(self.gromosTop, "BONDANGLEH"):
            self.gromosTop.add_block(blocktitle="BONDANGLEH", content=[])
        if not hasattr(self.gromosTop, "BONDANGLE"):
            self.gromosTop.add_block(blocktitle="BONDANGLE", content=[])

    def convertTosions(self):
        for molecule in self.molecule_force_list:
            for key in molecule["ProperTorsions"]:
                force = molecule["ProperTorsions"][key]
                # hQ = not all([self.openFFTop.atom(x).atomic_number != 1 for x in key])  # noqa: F841
                atomI = key[0] + 1
                atomJ = key[1] + 1
                atomK = key[2] + 1
                atomL = key[3] + 1
                k_list = force.k
                phase_list = force.phase
                per_list = force.periodicity
                for t in range(len(k_list)):
                    CP = k_list[t].value_in_unit(u.kilojoule_per_mole)
                    PD = phase_list[t].value_in_unit(u.degree)
                    NP = per_list[t]
                    # convert negativ CP by phase shifting
                    if CP < 0:
                        CP = abs(CP)
                        PD += 180
                    self.gromosTop.add_new_torsiondihedral(
                        CP=CP, PD=PD, NP=NP, atomI=atomI, atomJ=atomJ, atomK=atomK, atomL=atomL, includesH=False
                    )  # hQ
        if not hasattr(self.gromosTop, "TORSDIHEDRALTYPE"):
            self.gromosTop.add_block(blocktitle="TORSDIHEDRALTYPE", content=[])
        if not hasattr(self.gromosTop, "DIHEDRALH"):
            self.gromosTop.add_block(blocktitle="DIHEDRALH", content=[])
        if not hasattr(self.gromosTop, "DIHEDRAL"):
            self.gromosTop.add_block(blocktitle="DIHEDRAL", content=[])

    def convertImproper(self):
        for molecule in self.molecule_force_list:
            for key in molecule["ImproperTorsions"]:
                force = molecule["ImproperTorsions"][key]
                # hQ = not all([self.openFFTop.atom(x).atomic_number != 1 for x in key])  # noqa: F841
                atomI = key[0] + 1
                atomJ = key[1] + 1
                atomK = key[2] + 1
                atomL = key[3] + 1
                k_list = force.k
                phase_list = force.phase
                per_list = force.periodicity
                for t in range(len(k_list)):
                    CP = k_list[t].value_in_unit(u.kilojoule / u.mole)
                    PD = phase_list[t].value_in_unit(u.degree)
                    NP = per_list[t]
                    self.gromosTop.add_new_torsiondihedral(
                        CP=CP, PD=PD, NP=NP, atomI=atomI, atomJ=atomJ, atomK=atomK, atomL=atomL, includesH=False
                    )  # hQ
        if not hasattr(self.gromosTop, "IMPDIHEDRALTYPE"):
            self.gromosTop.add_block(blocktitle="IMPDIHEDRALTYPE", content=[])
        if not hasattr(self.gromosTop, "IMPDIHEDRALH"):
            self.gromosTop.add_block(blocktitle="IMPDIHEDRALH", content=[])
        if not hasattr(self.gromosTop, "IMPDIHEDRAL"):
            self.gromosTop.add_block(blocktitle="IMPDIHEDRAL", content=[])

    def createVdWexclusionList(self):
        bondDict = dict()
        ex13 = dict()
        ex14 = dict()
        # create a list of all bonds
        for molecule in self.molecule_force_list:
            for key in molecule["Bonds"]:
                if not str(key[0]) in bondDict.keys():
                    bondDict[str(key[0])] = {key[1]}
                bondDict[str(key[0])].add(key[1])
                if not str(key[1]) in bondDict.keys():
                    bondDict[str(key[1])] = {key[0]}
                bondDict[str(key[1])].add(key[0])
        # use bond dict to createexclusion lists
        for lvl1 in bondDict:
            ex13[lvl1] = bondDict[lvl1].copy()
            ex14[lvl1] = bondDict[lvl1].copy()  # just init 1-3 values will be removed later
            for lvl2 in bondDict[lvl1]:
                ex13[lvl1].add(lvl2)
                for lvl3 in bondDict[str(lvl2)]:
                    ex13[lvl1].add(lvl3)
                    for lvl4 in bondDict[str(lvl3)]:
                        ex14[lvl1].add(lvl4)
        # remove 1-3 from 1-4
        for key in ex14:
            for i in ex13[key]:
                ex14[key].discard(i)
        # remove all smaller entries
        for key in ex13:
            for i in range(0, int(key) + 1):
                ex13[key].discard(i)
        for key in ex14:
            for i in range(0, int(key) + 1):
                ex14[key].discard(i)

        # return
        self.exclusionList13 = ex13
        self.exclusionList14 = ex14

    def convertVdW(self):
        self.createVdWexclusionList()
        moleculeItr = 1
        prev_atom_counter = 0
        for molecule in self.molecule_force_list:
            panm_dict = collections.defaultdict(int)
            tot_len = len(molecule["vdW"])
            for key in molecule["vdW"]:
                force = molecule["vdW"][key]
                ATNM = int(key[0]) + 1 + prev_atom_counter
                MRES = moleculeItr
                # get element sympol:
                atomic_number = self.openFFmolecule.atoms[int(key[0])].atomic_number
                element_symbol = self.atomic_number_dict[atomic_number]
                panm_dict[element_symbol] += 1
                PANM = element_symbol + str(panm_dict[element_symbol])
                IAC = 0
                MASS = self.openFFmolecule.atoms[int(key[0])].mass.value_in_unit(u.dalton)
                CG = (
                    self.openmm_system.getForce(1)
                    .getParticleParameters(int(key[0]))[0]
                    .value_in_unit(u.elementary_charge)
                )
                CGC = 1 if (int(key[0]) + 1 == tot_len) else 0
                if str(key[0]) in self.exclusionList13:
                    openFFexList13 = list(self.exclusionList13[str(key[0])])
                    INE = [int(x) + 1 for x in openFFexList13]
                else:
                    INE = list()
                if str(key[0]) in self.exclusionList14:
                    openFFexList14 = list(self.exclusionList14[str(key[0])])
                    INE14 = [int(x) + 1 for x in openFFexList14]
                else:
                    INE14 = list()
                epsilon = float(force.epsilon.value_in_unit(u.kilojoule_per_mole))
                rmin = 2 * force.rmin_half.value_in_unit(u.nanometer)
                C6 = 2 * epsilon * (rmin**6)
                C12 = epsilon * (rmin**12)
                CS6 = 0.5 * C6  # factor 0.5 for 1-4 interaction. Standart in GROMOS and OpenFF
                CS12 = 0.5 * C12  # factor 0.5 for 1-4 interaction. Standart in GROMOS and OpenFF
                IACname = force.id
                self.gromosTop.add_new_atom(
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
            prev_atom_counter += tot_len

    def convert_other_stuff(self):
        if not hasattr(self.gromosTop, "SOLUTEMOLECULES"):
            self.gromosTop.add_block(blocktitle="SOLUTEMOLECULES", content=["1", str(self.openFFmolecule.n_atoms)])
        else:
            self.gromosTop.SOLUTEMOLECULES.content = [["1"], [str(self.openFFmolecule.n_atoms)]]
        if not hasattr(self.gromosTop, "TEMPERATUREGROUPS"):
            self.gromosTop.add_block(blocktitle="TEMPERATUREGROUPS", content=["1", str(self.openFFmolecule.n_atoms)])
        else:
            self.gromosTop.TEMPERATUREGROUPS.content = [["1"], [str(self.openFFmolecule.n_atoms)]]
        if not hasattr(self.gromosTop, "PRESSUREGROUPS"):
            self.gromosTop.add_block(blocktitle="PRESSUREGROUPS", content=["1", str(self.openFFmolecule.n_atoms)])
        else:
            self.gromosTop.PRESSUREGROUPS.content = [["1"], [str(self.openFFmolecule.n_atoms)]]
        if not hasattr(self.gromosTop, "LJEXCEPTIONS"):
            self.gromosTop.add_block(blocktitle="LJEXCEPTIONS", content=["0", ""])
        if not hasattr(self.gromosTop, "SOLVENTATOM"):
            self.gromosTop.add_block(blocktitle="SOLVENTATOM", content=["0", ""])
        if not hasattr(self.gromosTop, "SOLVENTCONSTR"):
            self.gromosTop.add_block(blocktitle="SOLVENTCONSTR", content=["0", ""])
        if not hasattr(self.gromosTop, "TOPVERSION"):
            self.gromosTop.add_block(blocktitle="TOPVERSION", content=["2.0"])
        if not hasattr(self.gromosTop, "PHYSICALCONSTANTS"):
            self.gromosTop.add_block(blocktitle="PHYSICALCONSTANTS", content=[""])
