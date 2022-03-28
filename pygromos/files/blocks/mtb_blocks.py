from pygromos.files.blocks.topology_blocks import TITLE as generic_TITLE
from pygromos.files.blocks.topology_blocks import FORCEFIELD as generic_FORCEFIELD
from pygromos.files.blocks.topology_blocks import MAKETOPVERSION as generic_MAKETOPVERSION
from pygromos.files.blocks._general_blocks import _generic_gromos_block, _generic_field
from pygromos.utils.typing import List


TITLE: generic_TITLE = generic_TITLE
FORCEFIELD: generic_FORCEFIELD = generic_FORCEFIELD
MAKETOPVERSION: generic_MAKETOPVERSION = generic_MAKETOPVERSION


class mtb_blocks(_generic_gromos_block):
    def __init__(self, name: str = None, used: bool = None, content: str = None):
        super().__init__(name, used, content)
        self.field_separator = "  "
        self.line_separator = " \n"
        self.field_continue_next_line = "\n\t\t\t\t\t\t\t\t\t\t"


class mtb_fields(_generic_field):
    fieldseperator = "  "
    lineseperator = " \n"
    field_continue_next_line = "\n\t\t\t\t\t\t\t\t\t\t"


class mtb_atoms_field(mtb_fields):
    def __init__(
        self, ATOM: int, ANM: str, IACM: int, MASS: int, CGMI: float, CGM: int, MAE: int, MSAE: List[int]
    ) -> None:
        self.ATOM = int(ATOM)
        self.ANM = ANM
        self.IACM = int(IACM)
        self.MASS = int(MASS)
        self.CGMI = float(CGMI)
        self.CGM = int(CGM)
        self.MAE = int(MAE)
        self.MSAE = [int(x) for x in MSAE]

    def to_string(self) -> str:
        return_str = ""
        return_str += str(self.ATOM).rjust(6) + " "
        return_str += self.ANM.rjust(6) + " "
        return_str += str(self.IACM).rjust(6) + " "
        return_str += str(self.MASS).rjust(6) + " "
        # return_str += self.fieldseperator + str(self.CGMI)
        return_str += "{:.5f}".format(self.CGMI).rjust(10) + " "
        return_str += str(self.CGM).rjust(4) + " "
        return_str += str(self.MAE).rjust(4) + " "
        lcounter = 0
        temp_MSAE = len(self.MSAE)
        for iter in self.MSAE:
            return_str += self.fieldseperator + str(iter).strip()
            lcounter += 1
            if (lcounter % 6) == 0 and temp_MSAE > 6:
                return_str += self.lineseperator
                return_str += 49 * " "
                temp_MSAE -= 6
        return_str += self.lineseperator
        return return_str


class mtb_preceding_exclusions_field(mtb_fields):
    def __init__(self, ATOM: int, MAE: int, MSAE: List[int]) -> None:
        self.ATOM = int(ATOM)
        self.MAE = int(MAE)
        self.MSAE = [int(x) for x in MSAE]

    def to_string(self) -> str:
        return_str = ""
        return_str += self.fieldseperator + str(self.ATOM)
        return_str += self.fieldseperator + str(self.MAE)
        lcounter = 0
        temp_MSAE = len(self.MSAE)
        for iter in self.MSAE:
            return_str += self.fieldseperator + str(iter).strip()
            lcounter += 1
            if (lcounter % 6) == 0 and temp_MSAE > 6:
                return_str += self.field_continue_next_line
                temp_MSAE -= 6
        return_str += self.lineseperator
        return return_str


class mtb_trailing_atoms_field(mtb_fields):
    def __init__(self, ATOM: int, ANM: str, IACM: int, MASS: int, CGMI: float, CGM: int) -> None:
        self.ATOM = int(ATOM)
        self.ANM = ANM
        self.IACM = int(IACM)
        self.MASS = int(MASS)
        self.CGMI = float(CGMI)
        self.CGM = int(CGM)

    def to_string(self) -> str:
        return_str = ""
        return_str += self.fieldseperator + str(self.ATOM)
        return_str += self.fieldseperator + self.ANM
        return_str += self.fieldseperator + str(self.IACM)
        return_str += self.fieldseperator + str(self.MASS)
        # return_str += self.fieldseperator + str(self.CGMI)
        return_str += self.fieldseperator + "{:.5f}".format(self.CGMI)
        return_str += self.fieldseperator + str(self.CGM)
        return_str += self.lineseperator
        return return_str


class mtb_atoms_solvent_field(mtb_fields):
    def __init__(self, ATOM: int, ANM: str, IACM: int, MASS: int, CG: float) -> None:
        self.ATOM = int(ATOM)
        self.ANM = ANM
        self.IACM = int(IACM)
        self.MASS = int(MASS)
        self.CG = float(CG)

    def to_string(self) -> str:
        return_str = ""
        return_str += self.fieldseperator + str(self.ATOM)
        return_str += self.fieldseperator + self.ANM
        return_str += self.fieldseperator + str(self.IACM)
        return_str += self.fieldseperator + str(self.MASS)
        return_str += self.fieldseperator + "{:.5f}".format(self.CG)
        # return_str += self.fieldseperator + str(self.CG)
        return_str += self.lineseperator
        return return_str


class mtb_bonds_field(mtb_fields):
    def __init__(self, IB: int, JB: int, MCB: int) -> None:
        self.IB = int(IB)
        self.JB = int(JB)
        self.MCB = int(MCB)

    def to_string(self) -> str:
        return_str = ""
        return_str += self.fieldseperator + str(self.IB)
        return_str += self.fieldseperator + str(self.JB)
        return_str += self.fieldseperator + str(self.MCB)
        return_str += self.lineseperator
        return return_str


class mtb_angles_field(mtb_fields):
    def __init__(self, IB: int, JB: int, KB: int, MCB: int) -> None:
        self.IB = int(IB)
        self.JB = int(JB)
        self.KB = int(KB)
        self.MCB = int(MCB)

    def to_string(self) -> str:
        return_str = ""
        return_str += self.fieldseperator + str(self.IB)
        return_str += self.fieldseperator + str(self.JB)
        return_str += self.fieldseperator + str(self.KB)
        return_str += self.fieldseperator + str(self.MCB)
        return_str += self.lineseperator
        return return_str


class mtb_dihedral_field(mtb_fields):
    def __init__(self, IB: int, JB: int, KB: int, LB: int, MCB: int) -> None:
        self.IB = int(IB)
        self.JB = int(JB)
        self.KB = int(KB)
        self.LB = int(LB)
        self.MCB = int(MCB)

    def to_string(self) -> str:
        return_str = ""
        return_str += self.fieldseperator + str(self.IB)
        return_str += self.fieldseperator + str(self.JB)
        return_str += self.fieldseperator + str(self.KB)
        return_str += self.fieldseperator + str(self.LB)
        return_str += self.fieldseperator + str(self.MCB)
        return_str += self.lineseperator
        return return_str


class mtb_lj_exceptions_field(mtb_fields):
    def __init__(self, iac: int, jac: int, mcb: int):
        self.iac = int(iac)
        self.jac = int(jac)
        self.mcb = int(mcb)

    def to_string(self) -> str:
        return_str = ""
        return_str += self.fieldseperator + str(self.iac)
        return_str += self.fieldseperator + str(self.jac)
        return_str += self.fieldseperator + str(self.mcb)
        return_str += self.lineseperator
        return return_str


class mtb_constraints_field(mtb_fields):
    def __init__(self, IB: int, JB: int, LENGTH: float) -> None:
        self.IB = int(IB)
        self.JB = int(JB)
        self.LENGTH = float(LENGTH)

    def to_string(self) -> str:
        return_str = ""
        return_str += self.fieldseperator + str(self.IB)
        return_str += self.fieldseperator + str(self.JB)
        return_str += self.fieldseperator + str(self.LENGTH)
        return_str += self.lineseperator
        return return_str


class MTBUILDBLSOLUTE(mtb_blocks):
    FORCEFIELD: FORCEFIELD
    MAKETOPVERSION: MAKETOPVERSION

    atoms: List[mtb_atoms_field]
    preceding_exclusions: List[mtb_preceding_exclusions_field]
    trailing_atoms: List[mtb_trailing_atoms_field]
    bonds: List[mtb_bonds_field]
    angles: List[mtb_angles_field]
    improper_dihedrals: List[mtb_dihedral_field]
    dihedrals: List[mtb_dihedral_field]
    lj_exceptions: List[mtb_lj_exceptions_field]

    skip_lines: int = 2

    def __init__(self, FORCEFIELD: FORCEFIELD = None, MAKETOPVERSION: MAKETOPVERSION = None, content=None):
        super().__init__(name=self.__class__.__name__, used=True, content=content)
        self.FORCEFIELD = FORCEFIELD
        self.MAKETOPVERSION = MAKETOPVERSION

    def read_content_from_str(self, content: str):
        # reset all storage
        self.atoms = []
        self.trailing_atoms = []
        self.preceding_exclusions = []
        self.bonds = []
        self.angles = []
        self.improper_dihedrals = []
        self.dihedrals = []
        self.lj_exceptions = []

        # first line, check if information on solute type is available
        if "@BLOCKTYPE" in content[0]:
            first_line = content[0].split()
            self.filename = first_line[1]
            self.residuecode = first_line[3]
            self.function = first_line[4]
            self.type = first_line[6]
            self.fullname = first_line[8]
        else:
            self.filename = None
            self.residuecode = None
            self.function = None
            self.type = None
            self.fullname = None

        itr = 1
        while itr < len(content):
            if content[itr].startswith("#"):
                itr += 1
                continue
            else:
                break

        self.RNME = content[itr].strip()

        itr += 1
        while itr < len(content):
            if content[itr].startswith("#"):
                itr += 1
                continue
            else:
                break

        self.NMAT = int(content[itr].split()[0])
        self.NLIN = int(content[itr].split()[1])
        itr += 1

        if self.NLIN != 0:
            nlin_found = 0
            while nlin_found < self.NLIN and itr < len(content):
                if content[itr].startswith("#"):
                    itr += 1
                    continue
                else:
                    dump1 = content[itr].strip().split()
                    if len(dump1) >= 2:
                        atom, mae = dump1[0:2]
                    else:
                        raise Exception(
                            "Error in ATOM block: \n"
                            + content[itr - 1]
                            + "\n"
                            + content[itr]
                            + "\n"
                            + content[itr + 1]
                            + "\n"
                            + content[itr + 2]
                            + "\n"
                        )
                    if 1 <= int(mae):
                        msae_values = [int(i) for i in dump1[2:]]
                        # keep reading in lines until we have all the data needed.
                        while int(mae) > len(msae_values):
                            itr += 1
                            try:
                                if any(i in content[itr] for i in ["\t\t\t\t\t", "                   "]):
                                    msae_values.extend([int(i) for i in content[itr].strip().split()])
                            except IOError:
                                raise IOError("Problem reading MSAE for anm=" + str(atom) + " mae=" + str(mae))
                    else:
                        msae_values = []
                    self.preceding_exclusions.append(mtb_preceding_exclusions_field(atom, mae, msae_values))
                    itr += 1
                    nlin_found += 1

            nmat_found = 0
            # try to find all atoms in the ATOM subblock and hope the while does not go to far...
            while nmat_found < (self.NMAT - self.NLIN) and itr < (len(content) - 10):
                if content[itr].startswith("#"):
                    itr += 1
                    continue
                else:
                    dump1 = content[itr].strip().split()
                    if len(dump1) >= 7:
                        atom, anm, iacm, mass, cgm, icgm, mae = dump1[0:7]
                    else:
                        raise Exception(
                            "Error in ATOM block: \n"
                            + content[itr - 1]
                            + "\n"
                            + content[itr]
                            + "\n"
                            + content[itr + 1]
                            + "\n"
                            + content[itr + 2]
                            + "\n"
                        )
                    if 1 <= int(mae):
                        msae_values = [int(i) for i in dump1[7:]]
                        # keep reading in lines until we have all the data needed.
                        while int(mae) > len(msae_values):
                            itr += 1
                            try:
                                if any(i in content[itr] for i in ["\t\t\t\t\t", "                   "]):
                                    msae_values.extend([int(i) for i in content[itr].strip().split()])
                            except IOError:
                                raise IOError("Problem reading MSAE for anm=" + str(anm) + " mae=" + str(mae))
                    else:
                        msae_values = []
                    self.atoms.append(mtb_atoms_field(atom, anm, iacm, mass, cgm, icgm, mae, msae_values))
                    itr += 1
                    nmat_found += 1
            if self.NMAT == 0:
                itr += self.skip_lines

            # trailing atoms
            trailing_atoms_found = 0
            while trailing_atoms_found < self.NLIN and itr < (len(content) - 10):
                if content[itr].startswith("#"):
                    itr += 1
                    continue
                else:
                    dump1 = content[itr].strip().split()
                    if len(dump1) >= 6:
                        atom, anm, iacm, mass, cgm, icgm = dump1[0:6]
                    else:
                        raise Exception(
                            "Error in ATOM block: \n"
                            + content[itr - 1]
                            + "\n"
                            + content[itr]
                            + "\n"
                            + content[itr + 1]
                            + "\n"
                            + content[itr + 2]
                            + "\n"
                        )
                    self.trailing_atoms.append(mtb_trailing_atoms_field(atom, anm, iacm, mass, cgm, icgm))
                    itr += 1
                    trailing_atoms_found += 1

        else:
            nmat_found = 0
            # try to find all atoms in the ATOM subblock and hope the while does not go to far...
            while nmat_found < self.NMAT and itr < (len(content) - 10):
                if content[itr].startswith("#"):
                    itr += 1
                    continue
                else:
                    dump1 = content[itr].strip().split()
                    if len(dump1) >= 7:
                        atom, anm, iacm, mass, cgm, icgm, mae = dump1[0:7]
                    else:
                        raise Exception(
                            "Error in ATOM block: \n"
                            + content[itr - 1]
                            + "\n"
                            + content[itr]
                            + "\n"
                            + content[itr + 1]
                            + "\n"
                            + content[itr + 2]
                            + "\n"
                        )
                    if 1 <= int(mae):
                        msae_values = [int(i) for i in dump1[7:]]
                        # keep reading in lines until we have all the data needed.
                        while int(mae) > len(msae_values):
                            itr += 1
                            try:
                                if any(i in content[itr] for i in ["\t\t\t\t\t", "                   "]):
                                    msae_values.extend([int(i) for i in content[itr].strip().split()])
                            except IOError:
                                raise IOError("Problem reading MSAE for anm=" + str(anm) + " mae=" + str(mae))
                    else:
                        msae_values = []
                    self.atoms.append(mtb_atoms_field(atom, anm, iacm, mass, cgm, icgm, mae, msae_values))
                    itr += 1
                    nmat_found += 1
            if self.NMAT == 0:
                itr += self.skip_lines
        # all atoms from ATOM subblock should be found and parsed at this point

        # read bonds
        while itr < len(content):
            if content[itr].startswith("#"):
                itr += 1
                continue
            else:
                try:
                    self.NB = int(content[itr].strip())
                    itr += 1
                    break
                except IOError:
                    raise IOError("Problem reading number of bonds")
        bonds_found = 0
        while itr < len(content) and bonds_found < self.NB:
            if content[itr].startswith("#"):
                itr += 1
                continue
            else:
                dump1 = content[itr].strip().split()
                ib, jb, mcb = dump1[0:3]
                self.bonds.append(mtb_bonds_field(ib, jb, mcb))
                bonds_found += 1
                itr += 1
        if self.NB == 0:
            itr += self.skip_lines

        # read angles
        while itr < len(content):
            if content[itr].startswith("#"):
                itr += 1
                continue
            else:
                try:
                    self.NBA = int(content[itr].strip())
                    itr += 1
                    break
                except IOError:
                    raise IOError("Problem reading number of angles")
        angles_found = 0
        while itr < len(content) and angles_found < self.NBA:
            if content[itr].startswith("#"):
                itr += 1
                continue
            else:
                dump1 = content[itr].strip().split()
                ib, jb, kb, mcb = dump1[0:4]
                self.angles.append(mtb_angles_field(ib, jb, kb, mcb))
                angles_found += 1
                itr += 1
        if self.NBA == 0:
            itr += self.skip_lines

        # read improper dihedrals
        while itr < len(content):
            if content[itr].startswith("#"):
                itr += 1
                continue
            else:
                try:
                    self.NIDA = int(content[itr].strip())
                    itr += 1
                    break
                except IOError:
                    raise IOError("Problem reading number of improper dihedrals")
        improper_dihedrals_found = 0
        while itr < len(content) and improper_dihedrals_found < self.NIDA:
            if content[itr].startswith("#"):
                itr += 1
                continue
            else:
                dump1 = content[itr].strip().split()
                ib, jb, kb, lb, mcb = dump1[0:5]
                self.improper_dihedrals.append(mtb_dihedral_field(ib, jb, kb, lb, mcb))
                improper_dihedrals_found += 1
                itr += 1
        if self.NIDA == 0:
            itr += self.skip_lines

        # read dihedrals
        while itr < len(content):
            if content[itr].startswith("#"):
                itr += 1
                continue
            else:
                try:
                    self.NDA = int(content[itr].strip())
                    itr += 1
                    break
                except IOError:
                    raise IOError("Problem reading number of dihedrals")
        dihedrals_found = 0
        while itr < len(content) and dihedrals_found < self.NDA:
            if content[itr].startswith("#"):
                itr += 1
                continue
            else:
                dump1 = content[itr].strip().split()
                ib, jb, kb, lb, mcb = dump1[0:5]
                self.dihedrals.append(mtb_dihedral_field(ib, jb, kb, lb, mcb))
                dihedrals_found += 1
                itr += 1
        if self.NDA == 0:
            itr += self.skip_lines

        # read LJ exceptions
        while itr < len(content):
            if content[itr].startswith("#"):
                itr += 1
                continue
            else:
                try:
                    self.NEX = int(content[itr].strip())
                    itr += 1
                    break
                except IOError:
                    raise IOError("Problem reading number of LJ exceptions")
        lj_exceptions_found = 0
        while itr < len(content) and lj_exceptions_found < self.NEX:
            if content[itr].startswith("#"):
                itr += 1
                continue
            else:
                dump1 = content[itr].strip().split()
                iac, jac, mcb = dump1[0:3]
                self.lj_exceptions.append(mtb_lj_exceptions_field(iac, jac, mcb))
                lj_exceptions_found += 1
                itr += 1

    def block_to_string(self) -> str:
        result = "MTBUILDBLSOLUTE" + self.line_separator
        if not [x for x in (self.filename, self.residuecode, self.function, self.type, self.fullname) if x is None]:
            result += (
                "# @BLOCKTYPE "
                + self.filename
                + " BLK "
                + self.residuecode
                + " "
                + self.function
                + " TYPE "
                + self.type
                + " NAME "
                + self.fullname
                + self.line_separator
            )
            result += "# building block" + self.line_separator

        atom_s = "ATOM".rjust(4)
        anm_s = "ANM".rjust(6)
        iacm_s = "IACM".rjust(6)
        mass_s = "MASS".rjust(6)
        cgm_s = "CGM".ljust(7)
        icgm_s = "ICGM".rjust(4)
        mae_s = "MAE".rjust(4)
        msae_s = "MSAE"

        result += "# RNME" + self.line_separator
        result += self.RNME + self.line_separator
        result += "# number of atoms, number of preceding exclusions" + self.line_separator
        result += "# NMAT NLIN" + self.line_separator
        result += str(self.NMAT) + self.field_separator + str(self.NLIN) + self.line_separator
        result += "# preceding exclusions" + self.line_separator
        result += "# ATOM                               MAE MSAE" + self.line_separator
        for pre in self.preceding_exclusions:
            result += pre.to_string()
        result += "# atoms" + self.line_separator
        result += f"# {atom_s} {anm_s} {iacm_s} {mass_s}    {cgm_s} {icgm_s} {mae_s}   {msae_s}" + self.line_separator
        for atom in self.atoms:
            result += atom.to_string()

        result += "# trailing atoms" + self.line_separator
        result += "# ATOM ANM  IACM MASS        CGM ICGM" + self.line_separator
        for trailing_atom in self.trailing_atoms:
            result += trailing_atom.to_string()

        result += "# bonds" + self.line_separator
        result += "# NB" + self.line_separator
        result += str(self.NB) + self.line_separator
        result += "#  IB   JB  MCB" + self.line_separator
        for bond in self.bonds:
            result += bond.to_string()

        result += "# angles" + self.line_separator
        result += "# NBA" + self.line_separator
        result += str(self.NBA) + self.line_separator
        result += "#  IB   JB   KB  MCB" + self.line_separator
        for angle in self.angles:
            result += angle.to_string()

        result += "# improper dihedrals" + self.line_separator
        result += "# NIDA" + self.line_separator
        result += str(self.NIDA) + self.line_separator
        result += "#  IB   JB   KB   LB  MCB" + self.line_separator
        for dihedral in self.improper_dihedrals:
            result += dihedral.to_string()

        result += "# dihedrals" + self.line_separator
        result += "# NDA" + self.line_separator
        result += str(self.NDA) + self.line_separator
        result += "#  IB   JB   KB   LB  MCB" + self.line_separator
        for dihedral in self.dihedrals:
            result += dihedral.to_string()

        result += "# LJ exceptions" + self.line_separator
        result += "# NEX" + self.line_separator
        result += str(self.NEX) + self.line_separator
        result += "# IAC  JAC  MCB" + self.line_separator
        for lj_exception in self.lj_exceptions:
            result += lj_exception.to_string()

        result += "# @FREELINE" + self.line_separator
        result += "END" + self.line_separator
        return result


class LINKEXCLUSIONS(mtb_blocks):
    FORCEFIELD: FORCEFIELD
    MAKETOPVERSION: MAKETOPVERSION
    NRNE: int

    def __init__(self, FORCEFIELD: FORCEFIELD = None, MAKETOPVERSION: MAKETOPVERSION = None, content=None):
        self.NRNE = 0
        super().__init__(name=self.__class__.__name__, used=True, content=content)
        self.FORCEFIELD = FORCEFIELD
        self.MAKETOPVERSION = MAKETOPVERSION

    def read_content_from_str(self, content: str):
        while content[0].startswith("#"):
            content = content[1:]
        self.NRNE = int(content[0].strip())

    def block_to_string(self) -> str:
        result = "LINKEXCLUSIONS" + self.line_separator
        result += "# nearest neighbour exclusions when linking" + self.line_separator
        result += "# NRNE" + self.line_separator
        result += str(self.NRNE) + self.line_separator
        result += "# @FREELINE" + self.line_separator
        result += "END" + self.line_separator
        return result


class MTBUILDBLSOLVENT(mtb_blocks):
    FORCEFIELD: FORCEFIELD
    MAKETOPVERSION: MAKETOPVERSION

    atoms: List[mtb_atoms_field]
    constraints: List[mtb_constraints_field]

    def __init__(self, FORCEFIELD: FORCEFIELD = None, MAKETOPVERSION: MAKETOPVERSION = None, content=None):
        super().__init__(name=self.__class__.__name__, used=True, content=content)
        self.FORCEFIELD = FORCEFIELD
        self.MAKETOPVERSION = MAKETOPVERSION

    def read_content_from_str(self, content: str):

        self.atoms = []
        self.constraints = []

        # first line, check if information on solute type is available
        if "@BLOCKTYPE" in content[0]:
            first_line = content[0].split()
            self.filename = first_line[1]
            self.residuecode = first_line[3]
            self.function = first_line[4]
            self.type = first_line[6]
            self.fullname = first_line[8]
        else:
            self.filename = None
            self.residuecode = None
            self.function = None
            self.type = None
            self.fullname = None

        while content[0].startswith("#"):
            content = content[1:]
        self.RNMES = content[0].strip()
        content = content[1:]

        while content[0].startswith("#"):
            content = content[1:]
        self.number_of_atoms = int(content[0].strip())
        content = content[1:]

        found_atoms = 0
        while found_atoms < self.number_of_atoms:
            if content[0].startswith("#"):
                content = content[1:]
            else:
                atom, anm, iac, mass, cg = content[0].split()
                self.atoms.append(mtb_atoms_solvent_field(atom, anm, iac, mass, cg))
                found_atoms += 1
                content = content[1:]

        while content[0].startswith("#"):
            content = content[1:]
        self.number_of_constraints = int(content[0].strip())
        content = content[1:]

        found_constraints = 0
        while found_constraints < self.number_of_constraints:
            if content[0].startswith("#"):
                content = content[1:]
            else:
                ib, jb, mcb = content[0].split()
                self.constraints.append(mtb_constraints_field(ib, jb, mcb))
                found_constraints += 1
                content = content[1:]

    def block_to_string(self) -> str:
        result = "MTBUILDBLSOLVENT" + self.line_separator
        if not [x for x in (self.filename, self.residuecode, self.function, self.type, self.fullname) if x is None]:
            result += (
                "# @BLOCKTYPE "
                + self.filename
                + " BLK "
                + self.residuecode
                + " "
                + self.function
                + " TYPE "
                + self.type
                + " NAME "
                + self.fullname
                + self.line_separator
            )
        result += "# solvent name" + self.line_separator
        result += "# RNMES" + self.line_separator
        result += str(self.RNMES) + self.line_separator

        result += "# number of atoms" + self.line_separator
        result += str(self.number_of_atoms) + self.line_separator

        result += "# atoms" + self.line_separator
        result += "# ATOM ANM  IACM MASS        CG" + self.line_separator
        for atom in self.atoms:
            result += atom.to_string()

        result += "# number of constraints" + self.line_separator
        result += "# constraints" + self.line_separator
        result += str(self.number_of_constraints) + self.line_separator
        result += "#  IB   JB  LENGTH" + self.line_separator
        for constraint in self.constraints:
            result += constraint.to_string()

        result += "# @FREELINE" + self.line_separator
        result += "END" + self.line_separator
        return result


class MTBUILDBLEND(mtb_blocks):
    FORCEFIELD: FORCEFIELD
    MAKETOPVERSION: MAKETOPVERSION

    atoms: List[mtb_atoms_field]
    replacing_atoms: List[mtb_trailing_atoms_field]
    bonds: List[mtb_bonds_field]
    angles: List[mtb_angles_field]
    improper_dihedrals: List[mtb_dihedral_field]
    dihedrals: List[mtb_dihedral_field]
    lj_exceptions: List[mtb_lj_exceptions_field]

    skip_lines: int = 2

    def __init__(self, FORCEFIELD: FORCEFIELD = None, MAKETOPVERSION: MAKETOPVERSION = None, content=None):
        super().__init__(name=self.__class__.__name__, used=True, content=content)
        self.FORCEFIELD = FORCEFIELD
        self.MAKETOPVERSION = MAKETOPVERSION

    def read_content_from_str(self, content: str):
        # reset all storage
        self.atoms = []
        self.replacing_atoms = []
        self.bonds = []
        self.angles = []
        self.improper_dihedrals = []
        self.dihedrals = []
        self.lj_exceptions = []

        # first line, check if information on solute type is available
        if "@BLOCKTYPE" in content[0]:
            first_line = content[0].split()
            self.filename = first_line[1]
            self.residuecode = first_line[3]
            self.function = first_line[4]
            self.type = first_line[6]
            self.fullname = first_line[8]
        else:
            self.filename = None
            self.residuecode = None
            self.function = None
            self.type = None
            self.fullname = None

        itr = 1
        while itr < len(content):
            if content[itr].startswith("#"):
                itr += 1
                continue
            else:
                break

        self.RNME = content[itr].strip()

        itr += 1
        while itr < len(content):
            if content[itr].startswith("#"):
                itr += 1
                continue
            else:
                break

        self.NMAT = int(content[itr].split()[0])
        self.NREP = int(content[itr].split()[1])
        itr += 1

        if self.NREP >= 0:
            nmat_found = 0
            # try to find all atoms in the ATOM subblock and hope the while does not go to far...
            while nmat_found < (self.NMAT - self.NREP) and itr < (len(content) - 10):
                if content[itr].startswith("#"):
                    itr += 1
                    continue
                else:
                    dump1 = content[itr].strip().split()
                    if len(dump1) >= 7:
                        atom, anm, iacm, mass, cgm, icgm, mae = dump1[0:7]
                    else:
                        raise Exception(
                            "Error in ATOM block: \n"
                            + content[itr - 1]
                            + "\n"
                            + content[itr]
                            + "\n"
                            + content[itr + 1]
                            + "\n"
                            + content[itr + 2]
                            + "\n"
                        )
                    if 1 <= int(mae):
                        msae_values = [int(i) for i in dump1[7:]]
                        # keep reading in lines until we have all the data needed.
                        while int(mae) > len(msae_values):
                            itr += 1
                            try:
                                if any(i in content[itr] for i in ["\t\t\t\t\t", "                   "]):
                                    msae_values.extend([int(i) for i in content[itr].strip().split()])
                            except IOError:
                                raise IOError("Problem reading MSAE for anm=" + str(anm) + " mae=" + str(mae))
                    else:
                        msae_values = []
                    self.atoms.append(mtb_atoms_field(atom, anm, iacm, mass, cgm, icgm, mae, msae_values))
                    itr += 1
                    nmat_found += 1
            if self.NMAT == 0:
                itr += self.skip_lines

            # replacing atoms
            replacing_atoms_found = 0
            while replacing_atoms_found < self.NREP and itr < (len(content) - 10):
                if content[itr].startswith("#"):
                    itr += 1
                    continue
                else:
                    dump1 = content[itr].strip().split()
                    if len(dump1) >= 6:
                        atom, anm, iacm, mass, cgm, icgm = dump1[0:6]
                    else:
                        raise Exception(
                            "Error in ATOM block: \n"
                            + content[itr - 1]
                            + "\n"
                            + content[itr]
                            + "\n"
                            + content[itr + 1]
                            + "\n"
                            + content[itr + 2]
                            + "\n"
                        )
                    self.replacing_atoms.append(mtb_trailing_atoms_field(atom, anm, iacm, mass, cgm, icgm))
                    itr += 1
                    replacing_atoms_found += 1
        # no NREP or negative NREP
        else:
            nmat_found = 0
            # try to find all atoms in the ATOM subblock and hope the while does not go to far...
            while nmat_found < self.NMAT and itr < (len(content) - 10):
                if content[itr].startswith("#"):
                    itr += 1
                    continue
                else:
                    dump1 = content[itr].strip().split()
                    if len(dump1) >= 7:
                        atom, anm, iacm, mass, cgm, icgm, mae = dump1[0:7]
                    else:
                        raise Exception(
                            "Error in ATOM block: \n"
                            + content[itr - 1]
                            + "\n"
                            + content[itr]
                            + "\n"
                            + content[itr + 1]
                            + "\n"
                            + content[itr + 2]
                            + "\n"
                        )
                    if 1 <= int(mae):
                        msae_values = [int(i) for i in dump1[7:]]
                        # keep reading in lines until we have all the data needed.
                        while int(mae) > len(msae_values):
                            itr += 1
                            try:
                                if any(i in content[itr] for i in ["\t\t\t\t\t", "                   "]):
                                    msae_values.extend([int(i) for i in content[itr].strip().split()])
                            except IOError:
                                raise IOError("Problem reading MSAE for anm=" + str(anm) + " mae=" + str(mae))
                    else:
                        msae_values = []
                    self.atoms.append(mtb_atoms_field(atom, anm, iacm, mass, cgm, icgm, mae, msae_values))
                    itr += 1
                    nmat_found += 1
            if self.NMAT == 0:
                itr += self.skip_lines
        # all atoms from ATOM subblock should be found and parsed at this point

        # read bonds
        while itr < len(content):
            if content[itr].startswith("#"):
                itr += 1
                continue
            else:
                try:
                    self.NB = int(content[itr].strip())
                    itr += 1
                    break
                except IOError:
                    raise IOError("Problem reading number of bonds")
        bonds_found = 0
        while itr < len(content) and bonds_found < self.NB:
            if content[itr].startswith("#"):
                itr += 1
                continue
            else:
                dump1 = content[itr].strip().split()
                ib, jb, mcb = dump1[0:3]
                self.bonds.append(mtb_bonds_field(ib, jb, mcb))
                bonds_found += 1
                itr += 1
        if self.NB == 0:
            itr += self.skip_lines

        # read angles
        while itr < len(content):
            if content[itr].startswith("#"):
                itr += 1
                continue
            else:
                try:
                    self.NBA = int(content[itr].strip())
                    itr += 1
                    break
                except IOError:
                    raise IOError("Problem reading number of angles")
        angles_found = 0
        while itr < len(content) and angles_found < self.NBA:
            if content[itr].startswith("#"):
                itr += 1
                continue
            else:
                dump1 = content[itr].strip().split()
                ib, jb, kb, mcb = dump1[0:4]
                self.angles.append(mtb_angles_field(ib, jb, kb, mcb))
                angles_found += 1
                itr += 1
        if self.NBA == 0:
            itr += self.skip_lines

        # read improper dihedrals
        while itr < len(content):
            if content[itr].startswith("#"):
                itr += 1
                continue
            else:
                try:
                    self.NIDA = int(content[itr].strip())
                    itr += 1
                    break
                except IOError:
                    raise IOError("Problem reading number of improper dihedrals")
        improper_dihedrals_found = 0
        while itr < len(content) and improper_dihedrals_found < self.NIDA:
            if content[itr].startswith("#"):
                itr += 1
                continue
            else:
                dump1 = content[itr].strip().split()
                ib, jb, kb, lb, mcb = dump1[0:5]
                self.improper_dihedrals.append(mtb_dihedral_field(ib, jb, kb, lb, mcb))
                improper_dihedrals_found += 1
                itr += 1
        if self.NIDA == 0:
            itr += self.skip_lines

        # read dihedrals
        while itr < len(content):
            if content[itr].startswith("#"):
                itr += 1
                continue
            else:
                try:
                    self.NDA = int(content[itr].strip())
                    itr += 1
                    break
                except IOError:
                    raise IOError("Problem reading number of dihedrals")
        dihedrals_found = 0
        while itr < len(content) and dihedrals_found < self.NDA:
            if content[itr].startswith("#"):
                itr += 1
                continue
            else:
                dump1 = content[itr].strip().split()
                ib, jb, kb, lb, mcb = dump1[0:5]
                self.dihedrals.append(mtb_dihedral_field(ib, jb, kb, lb, mcb))
                dihedrals_found += 1
                itr += 1
        if self.NDA == 0:
            itr += self.skip_lines

        # read LJ exceptions
        while itr < len(content):
            if content[itr].startswith("#"):
                itr += 1
                continue
            else:
                try:
                    self.NEX = int(content[itr].strip())
                    itr += 1
                    break
                except IOError:
                    raise IOError("Problem reading number of LJ exceptions")
        lj_exceptions_found = 0
        while itr < len(content) and lj_exceptions_found < self.NEX:
            if content[itr].startswith("#"):
                itr += 1
                continue
            else:
                dump1 = content[itr].strip().split()
                iac, jac, mcb = dump1[0:3]
                self.lj_exceptions.append(mtb_lj_exceptions_field(iac, jac, mcb))
                lj_exceptions_found += 1
                itr += 1

    def block_to_string(self) -> str:
        result = "MTBUILDBLSOLUTE" + self.line_separator
        if not [x for x in (self.filename, self.residuecode, self.function, self.type, self.fullname) if x is None]:
            result += (
                "# @BLOCKTYPE "
                + self.filename
                + " BLK "
                + self.residuecode
                + " "
                + self.function
                + " TYPE "
                + self.type
                + " NAME "
                + self.fullname
                + self.line_separator
            )
            result += "# building block" + self.line_separator

        result += "# RNME" + self.line_separator
        result += self.RNME + self.line_separator
        result += "# number of atoms, number of preceding exclusions" + self.line_separator
        result += "# NMAT NREP" + self.line_separator
        result += str(self.NMAT) + self.field_separator + str(self.NREP) + self.line_separator
        result += "# ATOM                               MAE MSAE" + self.line_separator
        result += "# atoms" + self.line_separator
        result += "# ATOM ANM  IACM MASS        CGMICGM MAE MSAE" + self.line_separator
        for atom in self.atoms:
            result += atom.to_string()

        result += "# replacing atoms" + self.line_separator
        result += "# ATOM ANM  IACM MASS        CGMICGM" + self.line_separator
        for atom in self.replacing_atoms:
            result += atom.to_string()

        result += "# bonds" + self.line_separator
        result += "# NB" + self.line_separator
        result += str(self.NB) + self.line_separator
        result += "#  IB   JB  MCB" + self.line_separator
        for bond in self.bonds:
            result += bond.to_string()

        result += "# angles" + self.line_separator
        result += "# NBA" + self.line_separator
        result += str(self.NBA) + self.line_separator
        result += "#  IB   JB   KB  MCB" + self.line_separator
        for angle in self.angles:
            result += angle.to_string()

        result += "# improper dihedrals" + self.line_separator
        result += "# NIDA" + self.line_separator
        result += str(self.NIDA) + self.line_separator
        result += "#  IB   JB   KB   LB  MCB" + self.line_separator
        for dihedral in self.improper_dihedrals:
            result += dihedral.to_string()

        result += "# dihedrals" + self.line_separator
        result += "# NDA" + self.line_separator
        result += str(self.NDA) + self.line_separator
        result += "#  IB   JB   KB   LB  MCB" + self.line_separator
        for dihedral in self.dihedrals:
            result += dihedral.to_string()

        result += "# LJ exceptions" + self.line_separator
        result += "# NEX" + self.line_separator
        result += str(self.NEX) + self.line_separator
        result += "# IAC  JAC  MCB" + self.line_separator
        for lj_exception in self.lj_exceptions:
            result += lj_exception.to_string()

        result += "# @FREELINE" + self.line_separator
        result += "END" + self.line_separator
        return result
