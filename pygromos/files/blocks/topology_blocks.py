import re
import inspect
import numpy as np
from enum import Enum
from pygromos.files.blocks._general_blocks import TITLE as generic_TITLE
from pygromos.files.blocks._general_blocks import _generic_gromos_block, _iterable_gromos_block, _generic_field
from pygromos.utils.typing import (
    Union,
    Iterable,
    List,
    Tuple,
    Dict,
    Number,
    _topology_table_block_Type,
    _iterable_topology_block_Type,
    PHYSICALCONSTANTS_Type,
    TOPVERSION_Type,
    ATOMTYPENAME_Type,
    RESNAME_Type,
)

from pygromos.files.blocks.pertubation_blocks import pertubation_lam_state_nonbonded as pertubation_lam_state

"""
   ENUMS
"""


class distant_Restraint_Type(Enum):
    half_harmonic_repulsive = -1
    full_harmonic_distance_restraint = 0
    half_harmonic_attractive = 1


class geometric_code(Enum):
    real_atom = 0
    virtual_H_atom_aliphaticC = 1
    virtual_H_atom_aromaticC = 2
    virtual_H_atoms_goc_aliph = 3
    virtual_H_atom_aliphaticC_strange = 4
    pseudo_H_atom_goc_H_atoms_CH3 = 5
    pseudo_H_atoms_goc_of_two_CH3 = 6
    pseudo_H_atoms_goc_of_three_CH3 = 7
    virtual_atoms_cog = -1
    virtual_atoms_com = -2


"""
   FIELDS
"""


class atom_pair_distanceRes(_generic_field):
    def __init__(
        self,
        i1: int,
        j1: int,
        k1: int,
        l1: int,
        type1: Union[geometric_code, int],
        i2: int,
        j2: int,
        k2: int,
        l2: int,
        type2: Union[geometric_code, int],
        r0: float,
        w0: float,
        rah: Union[distant_Restraint_Type, int],
        comment: str = "",
    ):
        """

        Parameters
        ----------
        i1 :    int
            id of atom i of first molecule
        j1 : int
            id of atom j of first molecule
        k1 : int
            id of atom k of first molecule
        l1 : int
             id of atom l of first molecule
        type1 : int
            geometric restraintype of first molecule
        i2 : int
            id of atom i of second molecule
        j2 : int
            id of atom j of second molecule
        k2 : int
            id of atom k of second molecule
        l2 : int
            id of atom l of second molecule
        type2 : int
            geometric restraintype of second molecule
        r0 :    float
            radius_0 of restraint
        w0 :    float
            weighting of restraint
        rah : int
            restraint_type
        comment :   str, optional
            comment for this restraint

        """

        try:
            self.atom1i = int(i1)
            self.atom1j = int(j1)
            self.atom1k = int(k1)
            self.atom1l = int(l1)

            if type(type1) is geometric_code:
                self.atom1ic = type1
            elif type(type1) is int or (type(type1) is str and str(type1).isdigit()):
                self.atom1ic = geometric_code(int(type1))
            else:
                raise ValueError("geometric index.rst atom1ic in atom_pair_distanceRes unknown\n" + str(type1))

            self.atom2i = int(i2)
            self.atom2j = int(j2)
            self.atom2k = int(k2)
            self.atom2l = int(l2)

            if type(type2) is geometric_code:
                self.atom2ic = type2
            elif type(type2) is int or (type(type2) is str and str(type2).isdigit()):
                self.atom2ic = geometric_code(int(type2))
            else:
                raise ValueError("geometric index.rst atom2ic in atom_pair_distanceRes unknown\n" + str(type2))

            self.radius_0 = float(r0)
            self.weight = float(w0)
            if type(rah) is int or (type(rah) is str and str(rah).isdigit()):
                self.disResType = distant_Restraint_Type(int(rah))
            elif type(rah) is distant_Restraint_Type:
                self.disResType = rah
            else:
                raise ValueError("DisresType in atom_pair_distanceRes unknown\n" + str(rah))
            self.comment = comment
        except IOError:
            raise IOError("COULD NOT convert a parameter for distancerestraint field into correct form!")

    def to_string(self):
        if len(self.comment) > 0 and not self.comment.endswith("\n"):
            self.comment += "\n"
        return (
            self.comment
            + "{:>5} {:>5} {:>5} {:>5} {:>5}    {:>5} {:>5} {:>5} {:>5} {:>5}  {:10.5f} {:10.5f} {:>3}\n".format(
                self.atom1i,
                self.atom1j,
                self.atom1k,
                self.atom1l,
                self.atom1ic.value,
                self.atom2i,
                self.atom2j,
                self.atom2k,
                self.atom2l,
                self.atom2ic.value,
                self.radius_0,
                self.weight,
                self.disResType.value,
            )
        )


class atom_mass_type(_generic_field):
    def __init__(self, N: int, ATMAS: float, ATMASN: str, comment: str = ""):
        self.N = N
        self.ATMAS = ATMAS
        self.ATMASN = ATMASN
        self.comment = comment

    def to_string(self):
        return self.comment + "{:<8} {:<3.4f} {:<5}\n".format(self.N, self.ATMAS, self.ATMASN)


"""
    BONDED TERM BLOCKS
"""


class bond_type(_generic_field):
    def __init__(
        self,
        ICB: int,
        CB: float,
        HB: float,
        B0: float,
        atomI: Union[str, Iterable[str]],
        atomJ: Union[str, Iterable[str]],
        specialNumber: int,
    ):
        """
               GROMOS bond-stretching parameters for one possible bond

        Parameters
        ----------
        ICB: int
            Bond type code
        CB: float
            Quartic bond-stretch force constant
        HB: float
            Harmonic bond-stretch force constant
        B0: float
            Ideal bond length
        atomI: str, Iterable[str]
            possible atom I useages
        atomJ: str, Iterable[str]
            possible atom J useages
        specialNumber: int
            No Idea @Todo: refactor later!
        """
        self.ICB = ICB
        self.CB = CB
        self.HB = HB
        self.B0 = B0
        self.atomI = atomI
        self.atomJ = atomJ
        self.specialNumber = specialNumber

    def to_string(self):
        str_line = self.comment + "\t{:<3} {:<1.7f}   {:<1.7f}    {:<1.7f}\n".format(
            self.ICB, self.CB, self.HB, self.B0
        )

        atomI = self.atomI
        atomJ = self.atomJ
        if isinstance(atomI, list):
            atomI = ",".join(atomI)
        if isinstance(atomJ, list):
            atomJ = ",".join(atomJ)

        str_line += "#\t{:<3} - {:<3}    {:<5}\n".format(atomI, atomJ, self.specialNumber)
        return str_line


class angle_type(_generic_field):
    def __init__(
        self,
        ICT: int,
        CT: float,
        CHT: float,
        T0: float,
        atomI: Union[str, Iterable[str]],
        atomJ: Union[str, Iterable[str]],
        atomK: Union[str, Iterable[str]],
        specialNumber: int,
    ):
        """
               GROMOS bond-stretching parameters for one possible bond

        Parameters
        ----------
        ICT: int
            Bond-angle type code
        CT: float
            Non-harmonic force constant
        CHT: float
            Harmonic force constant
        T0: float
            Ideal bond angle
        atomI: str, Iterable[str]
            possible atom I useages
        atomJ: str, Iterable[str]
            possible atom J useages
        atomK: str, Iterable[str]
            possible atom K useages
        specialNumber: int
            No Idea @Todo: refactor later!
        """
        self.ICT = ICT
        self.CT = CT
        self.CHT = CHT
        self.T0 = T0
        self.atomI = atomI
        self.atomJ = atomJ
        self.atomK = atomK
        self.specialNumber = specialNumber

    def to_string(self):
        str_line = self.comment + "\t{:<3} {:<1.7f}   {:<1.7f}    {:<1.7f}\n".format(
            self.ICT, self.CT, self.CHT, self.T0
        )

        atomI = self.atomI
        atomJ = self.atomJ
        atomK = self.atomK
        if isinstance(atomI, list):
            atomI = ",".join(atomI)
        if isinstance(atomJ, list):
            atomJ = ",".join(atomJ)
        if isinstance(atomK, list):
            atomK = ",".join(atomK)

        str_line += "#\t{:<3} - {:<3} - {:<3}     {:<5}\n".format(atomI, atomJ, atomK, self.specialNumber)
        return str_line


class dihedral_type(_generic_field):
    def __init__(
        self,
        ICP: int,
        CP: float,
        PD: float,
        NP: int,
        atomI: str,
        atomJ: str,
        atomK: str,
        atomL: str,
        special_number: float,
        concrete_example: str = "",
    ):
        """
            GROMOS improper (harmonic) dihedral angle parameters

        Parameters
        ----------
        ICQ: int
            Dihedral-angle type code
        CQ: float
            Force constant
        Q0: float
            Phase shift
        NP: float
            Multiplicity
        atomI: str
            Examples for atomI
        atomJ: str
            Examples for atomJ
        atomK: str
            Examples for atomK
        atomL: str
            Examples for atomL
        special_number: int
            No Idea @Todo: refactor later!
        concrete_example: str, optional
            giving a concrete example
        """
        self.ICP = ICP
        self.CP = CP
        self.PD = PD
        self.NP = NP
        self.atomI = atomI
        self.atomJ = atomJ
        self.atomK = atomK
        self.atomL = atomL
        self.concrete_example = concrete_example
        self.special_number = special_number

    def to_string(self):
        str_line = self.comment + "\t{:<3}   {:<3.3f}    {:<3.3f}   {:<3}\n".format(self.ICP, self.CP, self.PD, self.NP)

        atomI = self.atomI
        atomJ = self.atomJ
        atomK = self.atomK
        atomL = self.atomL

        if isinstance(atomI, list):
            atomI = ",".join(atomI)
        if isinstance(atomJ, list):
            atomJ = ",".join(atomJ)
        if isinstance(atomK, list):
            atomK = ",".join(atomK)
        if isinstance(atomL, list):
            atomL = ",".join(atomL)

        str_line += "#\t{:<3} - {:<3} - {:<3} - {:<3}    {:<5}\n".format(
            atomI, atomJ, atomK, atomL, self.special_number
        )
        str_line += "#\t{:20}\n".format(self.concrete_example.replace("#", ""))

        return str_line


class improper_dihedral_type(_generic_field):
    def __init__(self, ICQ: int, CQ: float, Q0: float, group_type: str, special_number: int):
        """
            GROMOS improper (harmonic) dihedral angle parameters

        Parameters
        ----------
        ICQ: int
            Improper dihedral-angle type code
        CQ: float
            Force constant
        Q0: float
            Ideal improper dihedral angle
        group_type: str
            Example Useage
        special_number: int
            No Idea @Todo: refactor later!
        """
        self.ICQ = ICQ
        self.CQ = CQ
        self.Q0 = Q0
        self.group_type = group_type
        self.special_number = special_number

    def to_string(self):
        str_line = self.comment + "\t{:<3} {:<1.7f}    {:<1.7f}\n".format(self.ICQ, self.CQ, self.Q0)
        str_line += "#\t{:<3}    {:<5}\n".format(self.group_type, self.special_number)
        return str_line


class single_atom_lj_pair_type(_generic_field):
    def __init__(
        self,
        IAC: int,
        TYPE: str,
        C6: float,
        C12_1: float,
        C12_2: float,
        C12_3: float,
        CS6: float,
        CS12: float,
        LJ14PAIR: Iterable[float],
    ):
        """

        Parameters
        ----------
        IAC: int
            vander wals type index.rst
        TYPE: str
            atom type
        C6: float
            square-root of C6
        C12_2: float
            square-root of C6
        C12_3: float
            square-root of C6
        CS6: float
            No Idea @Todo: refactor later!
        CS12: float
            No Idea @Todo: refactor later!
        LJ14PAIR: Iterable[float]
            No Idea @Todo: refactor later!
        """
        self.IAC = IAC
        self.TYPE = TYPE
        self.C6 = C6
        self.C12_1 = C12_1
        self.C12_2 = C12_2
        self.C12_3 = C12_3
        self.CS6 = CS6
        self.CS12 = CS12
        self.LJ14PAIR = LJ14PAIR

    def to_string(self):
        str_line = self.comment + "#\tIAC TYPE C6 C12_1 C12_2 C12_3\n"
        str_line += "\t{:<3} {:<3} {:<1.7f}   {:<1.7f} {:<1.7f}    {:<1.7f}\n".format(
            self.IAC, self.TYPE, self.C6, self.C12_1, self.C12_2, self.C12_3
        )
        str_line += "#\t CS6 CS12\n"
        str_line += "\t{:<3f} {:<3f}\n".format(self.CS6, self.CS12)
        str_line += "#\tLJPAIRS\n"
        str_line += "\t" + "\n\t".join([" ".join(map(str, x)) for x in self.LJ14PAIR]) + "\n"
        str_line += "#---\n"

        return str_line


class mixed_atom_lj_pair_type(_generic_field):
    def __init__(self, IACI: int, IACJ: int, C6: float, C12_1: float, C12_2: float, C12_3: float):
        self.IACI = IACI
        self.IACJ = IACJ
        self.C6 = C6
        self.C12_1 = C12_1
        self.C12_2 = C12_2
        self.C12_3 = C12_3

    def to_string(self):
        str_line = self.comment + "\t{:<3} {:<3} {:<1.7f}   {:<1.7f} {:<1.7f}    {:<1.7f}\n".format(
            self.IACI, self.IACI, self.C6, self.C12_1, self.C12_2, self.C12_3
        )
        return str_line


class special_atom_lj_pair_type(_generic_field):
    def __init__(self, c):
        """

        Parameters
        ----------
        IAC: int
            vander wals type index.rst
        TYPE: str
            atom type
        C6: float
            square-root of C6
        C12_2: float
            square-root of C6
        C12_3: float
            square-root of C6
        CS6: float
            No Idea @Todo: refactor later!
        CS12: float
            No Idea @Todo: refactor later!
        LJ14PAIR: Iterable[float]
            No Idea @Todo: refactor later!
        """
        self.c

    def to_string(self):
        str_line = self.comment + "\t{:<3}\n".format(
            self.c,
        )
        return str_line


"""
    NON-BONDED TERM BLOCKS
"""

"""
BLOCKS
"""
# forward declarations
TITLE: generic_TITLE = generic_TITLE


# general Topo Blocks


class FORCEFIELD(_generic_gromos_block):
    NAME: str

    def __init__(self, NAME: str = None, content: List[str] = None):
        if content is None:
            super().__init__(name=self.__class__.__name__, used=True)
            self.NAME = NAME[0].strip()
        else:
            super().__init__(name=self.__class__.__name__, used=True, content=content)
            self.NAME = "\n".join(content)

    def block_to_string(self) -> str:
        result = self.name + self.line_seperator
        if hasattr(self, "NAME") and type(self.NAME) is str:
            result += self.NAME + self.line_seperator
        result += "END" + self.line_seperator
        return result


class MAKETOPVERSION(_generic_gromos_block):
    VERSION: str

    def __init__(self, VERSION: str = None, content: List[str] = None):
        if content is None:
            super().__init__(name=self.__class__.__name__, used=True)
            self.VERSION = VERSION[0].strip()
        else:
            super().__init__(name=self.__class__.__name__, used=True, content=content)
            self.VERSION = "\n".join(content)

    def block_to_string(self) -> str:
        result = self.name + self.line_seperator
        result += self.VERSION + self.line_seperator
        result += "END" + self.line_seperator
        return result


class _topology_block(_generic_gromos_block):
    FORCEFIELD: FORCEFIELD
    MAKETOPVERSION: MAKETOPVERSION

    def __init__(self, FORCEFIELD, MAKETOPVERSION, content: List[str] = None):
        super().__init__(name=self.__class__.__name__, used=True, content=content)
        self.FORCEFIELD = FORCEFIELD
        self.MAKETOPVERSION = MAKETOPVERSION


class _iterable_topology_block(_iterable_gromos_block):
    FORCEFIELD: FORCEFIELD
    MAKETOPVERSION: MAKETOPVERSION

    def __init__(self, FORCEFIELD: FORCEFIELD = None, MAKETOPVERSION: MAKETOPVERSION = None, content: List[str] = None):
        super().__init__(self.__class__.__name__, used=True, content=content)
        self.FORCEFIELD = FORCEFIELD
        self.MAKETOPVERSION = MAKETOPVERSION

    def __deepcopy__(self, memo):
        # return block as string, split by line and cut block title and END
        newContent = self.block_to_string().split(self.line_seperator)[1:-2]
        block = type(self)(content=newContent)
        return block


"""
TOPOLOGY BLOCKS
"""
"""
Restraints Blocks
"""


class DISTANCERESSPEC(_generic_gromos_block):
    def __init__(
        self,
        KDISH: int = None,
        KDISC: int = None,
        RESTRAINTHEADER: list = None,
        RESTRAINTS: list = None,
        content: List[str] = None,
    ):
        """

        Parameters
        ----------
        KDISH :
        KDISC :
        RESTRAINTHEADER :
        RESTRAINTS :
        """

        if content is None:
            content = [
                "# KDISH, KDISC\n",
                str(KDISH) + "\t" + str(KDISC),
                "\t".join(RESTRAINTHEADER),
            ] + list(map(str, RESTRAINTS))
        super().__init__(used=True, name="DISTANCERESSPEC", content=content)

    def read_content_from_str(self, content: List[str]):
        # readout KDISH or KDISC
        # keys = content[0].replace("#", "").strip().split()
        KDISH, KDISC = content[1].split()

        # read list header:
        line_header = content[2].replace("#", "").split()
        # unify keys:
        key_dict = {"i": 1, "j": 1, "k": 1, "l": 1, "type": 1}
        renamed_header = []
        for x in line_header:
            if x in key_dict:
                renamed_header.append(x + str(key_dict[x]))
                key_dict[x] += 1
            else:
                renamed_header.append(x)
        RESTRAINTHEADER = renamed_header

        # read restraints
        RESTRAINTS = []
        for line in content[3:]:
            if not line.startswith("#") and len(line.split()) == len(RESTRAINTHEADER):
                values = line.split()
                RESTRAINTS_dict = {key: values[RESTRAINTHEADER.index(key)] for key in RESTRAINTHEADER}
                RESTRAINTS.append(atom_pair_distanceRes(**RESTRAINTS_dict))
            elif line.startswith("#"):
                continue
            else:
                print("WARNING! could not Read in :" + line)
                continue

        self.KDISH = KDISH
        self.KDISC = KDISC
        self.RESTRAINTHEADER = RESTRAINTHEADER
        self.RESTRAINTS = RESTRAINTS

    def block_to_string(self) -> str:
        result = self.name + self.line_seperator
        result += "# KDISH" + self.field_seperator + "KDISC" + self.line_seperator
        result += self.field_seperator + str(self.KDISH) + self.field_seperator + str(self.KDISC) + self.line_seperator
        result += "#{:>4} {:>5} {:>5} {:>5} {:>5}    {:>5} {:>5} {:>5} {:>5} {:>5}  {:>10} {:>10} {:>3}\n".format(
            *self.RESTRAINTHEADER
        )
        for x in self.RESTRAINTS:
            result += x.to_string()
        result += "END\n"
        return result


"""
FORCEFIELD BLOCKS
"""

"""
IFP - TYPE Blocks
"""


class MASSATOMTYPECODE(_iterable_topology_block):
    NRMATY: int
    NMATY: int
    table_header: Iterable[str] = ["N", "ATMAS", "ATMASN"]

    def __init__(
        self,
        content: Union[Iterable[atom_mass_type], Iterable[str]],
        FORCEFIELD: FORCEFIELD = None,
        MAKETOPVERSION: MAKETOPVERSION = None,
        NRMATY: int = None,
        NMATY: int = None,
    ):
        super().__init__(FORCEFIELD=FORCEFIELD, MAKETOPVERSION=MAKETOPVERSION, content=content)

        self._content = []
        if isinstance(content, list) and all([isinstance(x, str) for x in content]):
            self.read_content_from_str(content)
        elif isinstance(content, Iterable) and all([isinstance(x, atom_mass_type) for x in content]):
            self.content = content
        else:
            raise IOError("I don't understand the type of content: " + str(type(content)))

        if NRMATY is None:
            self.NRMATY = len(self.content)
        elif isinstance(NRMATY, int):
            if NRMATY == len(self.content):  # CHECK FOR POSSIBLE ERROR
                self.NRMATY = NRMATY
            else:
                raise ValueError("IN MASSATOMTYPECODE NRMATY is not equal to the ammount of MASSATOMTYPES.")
        else:
            raise IOError("I don't understand the type of NRMATY: " + str(type(NRMATY)))

        if NMATY is None:
            self.NMATY = max([x.N for x in self.content])
        elif isinstance(NMATY, int):
            if NMATY == max([x.N for x in self.content]):  # CHECK FOR POSSIBLE ERROR
                self.NMATY = NMATY
            else:
                raise ValueError("IN MASSATOMTYPECODE NMATY is not the maximal Mass atom type code.")
        else:
            raise IOError("I don't understand the type of NMATY: " + str(type(NMATY)))

    def read_content_from_str(self, content: Union[str, List[str]]):

        if isinstance(content, str):
            lines = content.split("\n")
        else:
            lines = content

        if "NRMATY" in lines[0] and "NMATY" in lines[0]:
            NRMATY, NMATY = list(map(int, lines[1].strip().split()))

        table_start = 0
        for line in lines:
            table_start += 1
            if all([field in line for field in self.table_header]):
                break

        if table_start == len(lines):
            raise ValueError("Could not find the TABLE start in " + self.name)

        for field in lines[table_start:]:
            mass_atom_type_code, mass, mass_atom_name = field.strip().split()
            self.content.append(atom_mass_type(int(mass_atom_type_code), float(mass), str(mass_atom_name)))

        # print(self.content)

    def block_to_string(self) -> str:
        result = self.name + "\n"
        result += "#" + self.field_seperator + "NRMATY" + self.field_seperator + "NMATY" + self.line_seperator
        result += self.field_seperator + str(self.NRMATY) + self.field_seperator + str(self.NMATY) + self.line_seperator
        result += "# TABLE CONTENT: " + self.line_seperator
        result += "#" + self.field_seperator + self.field_seperator.join(self.table_header) + self.line_seperator
        for x in self.content:
            result += x.to_string()
        result += "END\n"

        return result


class BONDSTRETCHTYPECODE(_iterable_topology_block):
    NRBTY: int
    NBTY: int
    table_header: Iterable[str] = ["ICB(H)[N]", "CB[N]", "HB[N]", "B0[N]"]

    def __init__(
        self,
        content: Union[Iterable[bond_type], str],
        FORCEFIELD: FORCEFIELD = None,
        MAKETOPVERSION: MAKETOPVERSION = None,
        NRBTY: int = None,
        NBTY: int = None,
    ):
        """
                       GROMOS bond-stretching parameters

        Parameters
        ----------
        content: Union[Iterable[atom_mass_type], str]
        FORCEFIELD : FORCEFIELD
        MAKETOPVERSION : MAKETOPVERSION
        NRBTY : int, optional
            Number of bond types
        NBTY : int, optional
            Number of maximal bond index.rst
        """
        super().__init__(FORCEFIELD=FORCEFIELD, MAKETOPVERSION=MAKETOPVERSION, content=content)

        if isinstance(content, list) and all([isinstance(x, str) for x in content]):
            self.read_content_from_str(content)
        elif isinstance(content, Iterable) and all([isinstance(x, atom_mass_type) for x in content]):
            self.content = content
        else:
            raise IOError("I don't understand the type of content: " + str(type(content)))

        if NRBTY is None:
            self.NRBTY = len(self.content)
        elif isinstance(NRBTY, int):
            if NRBTY == len(self.content):  # CHECK FOR POSSIBLE ERROR
                self.NRBTY = NRBTY
            else:
                raise ValueError("IN MASSATOMTYPECODE NRMATY is not equal to the ammount of MASSATOMTYPES.")
        else:
            raise IOError("I don't understand the type of NRMATY: " + str(type(NRBTY)))

        if NBTY is None:
            self.NBTY = max([x.ICB for x in self.content])
        elif isinstance(NBTY, int):
            if NBTY == max([x.ICB for x in self.content]):  # CHECK FOR POSSIBLE ERROR
                self.NBTY = NBTY
            else:
                raise ValueError("IN MASSATOMTYPECODE NMATY is not the maximal Mass atom type code.")
        else:
            raise IOError("I don't understand the type of NMATY: " + str(type(NBTY)))

    def read_content_from_str(self, content: str):
        if isinstance(content, str):
            lines = content.split("\n")
        else:
            lines = content

        table_start = 0
        for line in lines:
            table_start += 1
            if all([field in line for field in self.table_header]):
                break

        if table_start == len(lines):
            raise ValueError("Could not find the TABLE start in " + self.name)

        table_lines = lines[table_start:]
        for field in table_lines:
            if not field.strip().startswith("#") and not len(field.strip()) == 0:
                ICB, CB, HB, B0 = field.strip().split()

                # nasty next line formatting of gromos
                comment_index = table_lines.index(field) + 1
                clean = re.sub("\(.*?\)", "", table_lines[comment_index])  # noqa # flake8: noqa
                split_line = [x for x in clean.replace("#", "").replace("-", "").strip().split("  ") if (len(x) > 0)]
                split_line = [x.split(",") if ("," in x) else x for x in split_line]
                if len(split_line) == 3:
                    atomI, atomJ, specialNumber = split_line
                elif len(split_line) == 2:
                    atomI, atomJ = split_line
                    specialNumber = 0

                # generate line
                params = bond_type(
                    ICB=int(ICB),
                    CB=float(CB),
                    HB=float(HB),
                    B0=float(B0),
                    atomI=atomI,
                    atomJ=atomJ,
                    specialNumber=float(specialNumber),
                )

                self.content.append(params)

    def block_to_string(self) -> str:
        result = self.name + "\n"
        result += "#" + self.field_seperator + "NRBTY" + self.field_seperator + "NBTY" + self.line_seperator
        result += "# TABLE CONTENT: " + self.line_seperator
        result += self.field_seperator + str(self.NRBTY) + self.field_seperator + str(self.NBTY) + self.line_seperator
        result += "#" + self.field_seperator + self.field_seperator.join(self.table_header) + self.line_seperator
        for x in self.content:
            result += x.to_string()
        result += "END\n"
        return result


class BONDANGLEBENDTYPECODE(_iterable_topology_block):
    NRTTY: int
    NTTY: int
    table_header: Iterable[str] = ["ICT(H)[N]", "CT[N]", "CHT[N]", "(T0[N])"]

    def __init__(
        self,
        content: Union[Iterable[angle_type], str],
        FORCEFIELD: FORCEFIELD = None,
        MAKETOPVERSION: MAKETOPVERSION = None,
        NRTTY: int = None,
        NTTY: int = None,
    ):
        """
                       GROMOS bond-stretching parameters

        Parameters
        ----------
        content: Union[Iterable[atom_mass_type], str]
        FORCEFIELD : FORCEFIELD
        MAKETOPVERSION : MAKETOPVERSION
        NRBTY : int, optional
            Number of bond types
        NBTY : int, optional
            Number of maximal bond index.rst
        """
        super().__init__(FORCEFIELD=FORCEFIELD, MAKETOPVERSION=MAKETOPVERSION, content=content)

        if isinstance(content, list) and all([isinstance(x, str) for x in content]):
            self.read_content_from_str(content)
        elif isinstance(content, Iterable) and all([isinstance(x, atom_mass_type) for x in content]):
            self.content = content
        else:
            raise IOError("I don't understand the type of content: " + str(type(content)))

        if NRTTY is None:
            self.NRTTY = len(self.content)
        elif isinstance(NRTTY, int):
            if NRTTY == len(self.content):  # CHECK FOR POSSIBLE ERROR
                self.NRBTY = NRTTY
            else:
                raise ValueError("IN MASSATOMTYPECODE NRMATY is not equal to the ammount of MASSATOMTYPES.")
        else:
            raise IOError("I don't understand the type of NRMATY: " + str(type(NRTTY)))

        if NTTY is None:
            self.NTTY = max([x.ICT for x in self.content])
        elif isinstance(NTTY, int):
            if NTTY == max([x.ICT for x in self.content]):  # CHECK FOR POSSIBLE ERROR
                self.NTTY = NTTY
            else:
                raise ValueError("IN MASSATOMTYPECODE NMATY is not the maximal Mass atom type code.")
        else:
            raise IOError("I don't understand the type of NMATY: " + str(type(NTTY)))

    def read_content_from_str(self, content: str):
        if isinstance(content, str):
            lines = content.split("\n")
        else:
            lines = content

        table_start = 0
        for line in lines:
            table_start += 1
            if all([field in line for field in self.table_header]):
                break

        if table_start == len(lines):
            raise ValueError("Could not find the TABLE start in " + self.name)

        table_lines = lines[table_start:]
        for field in table_lines:
            if not field.strip().startswith("#") and not len(field.strip()) == 0:
                ICT, CT, CHT, T0 = field.strip().split()

                # nasty next line formatting of gromos
                comment_index = table_lines.index(field) + 1
                clean = re.sub("\(.*?\)", "", table_lines[comment_index])  # noqa # flake8: noqa
                split_line = [x for x in clean.replace("#", "").replace("-", "").strip().split("  ") if (len(x) > 0)]
                split_line = [x.split(",") if ("," in x) else x for x in split_line]

                if len(split_line) == 4:
                    atomI, atomJ, atomK, specialNumber = split_line
                elif len(split_line) == 3:
                    atomI, atomJ, atomK = split_line
                    specialNumber = 0

                params = angle_type(
                    ICT=int(ICT),
                    CT=float(CT),
                    CHT=float(CHT),
                    T0=float(T0),
                    atomI=atomI,
                    atomJ=atomJ,
                    atomK=atomK,
                    specialNumber=int(specialNumber),
                )
                self.content.append(params)

    def block_to_string(self) -> str:
        result = self.name + "\n"
        result += "#" + self.field_seperator + "NRTTY" + self.field_seperator + "NTTY" + self.line_seperator
        result += self.field_seperator + str(self.NRTTY) + self.field_seperator + str(self.NTTY) + self.line_seperator
        result += "# TABLE CONTENT: " + self.line_seperator
        result += "#" + self.field_seperator + self.field_seperator.join(self.table_header) + self.line_seperator
        for x in self.content:
            result += x.to_string()
        result += "END\n"
        return result


class TORSDIHEDRALTYPECODE(_iterable_topology_block):
    NRPTY: int
    NPTY: int

    table_header: Iterable[str] = ["ICP(H)[N]", "CP[N]", "PD", "NP"]

    def __init__(
        self,
        content: Union[Iterable[dihedral_type], str],
        FORCEFIELD: FORCEFIELD = None,
        MAKETOPVERSION: MAKETOPVERSION = None,
        NRPTY: int = None,
        NPTY: int = None,
    ):
        """
                       GROMOS (trigonometric) dihedral torsional angle parameters

        Parameters
        ----------
        content: Union[Iterable[atom_mass_type], str]
        FORCEFIELD : FORCEFIELD
        MAKETOPVERSION : MAKETOPVERSION
        NRPTY : int, optional
            Number of dihedral-angle types
        NPTY : int, optional
            Number of maximal dihedral-angle index.rst
        """
        super().__init__(FORCEFIELD=FORCEFIELD, MAKETOPVERSION=MAKETOPVERSION, content=content)

        if isinstance(content, list) and all([isinstance(x, str) for x in content]):
            self.read_content_from_str(content)
        elif isinstance(content, Iterable) and all([isinstance(x, atom_mass_type) for x in content]):
            self.content = content
        else:
            raise IOError("I don't understand the type of content: " + str(type(content)))

        if NRPTY is None:
            self.NRPTY = len(self.content)
        elif isinstance(NRPTY, int):
            if NRPTY == len(self.content):  # CHECK FOR POSSIBLE ERROR
                self.NRPTY = NRPTY
            else:
                raise ValueError("IN MASSATOMTYPECODE NRMATY is not equal to the ammount of MASSATOMTYPES.")
        else:
            raise IOError("I don't understand the type of NRPTY: " + str(type(NRPTY)))

        if NPTY is None:
            self.NPTY = max([x.ICP for x in self.content])
        elif isinstance(NPTY, int):
            if NPTY == max([x.ICP for x in self.content]):  # CHECK FOR POSSIBLE ERROR
                self.NPTY = NPTY
            else:
                raise ValueError("IN MASSATOMTYPECODE NPTY is not the maximal Mass atom type code.")
        else:
            raise IOError("I don't understand the type of NMATY: " + str(type(NPTY)))

    def read_content_from_str(self, content: str):
        if isinstance(content, str):
            lines = content.split("\n")
        else:
            lines = content

        table_start = 0
        for line in lines:
            table_start += 1
            if all([field in line for field in self.table_header]):
                break

        if table_start == len(lines):
            raise ValueError("Could not find the TABLE start in " + self.name)

        table_lines = lines[table_start:]
        for field in table_lines:
            if not field.strip().startswith("#") and not len(field.strip()) == 0:
                ICP, CP, PD, NP = field.strip().split()

                # nasty next line formatting of gromos
                comment_index = table_lines.index(field) + 1
                clean = re.sub("\(.*?\)", "", table_lines[comment_index])  # noqa # flake8: noqa
                split_line = [x for x in clean.replace("#", "").strip().split(" ") if (len(x) > 0)]
                # print("splits", split_line)

                special_number = float(split_line[-1]) if (split_line[-1].isdigit()) else 0
                atom_dihedral_string = ["X" if (len(x) == 0) else x for x in split_line[0].split("-")]
                atom_dihedral_string = [
                    ["X" if (len(y) == 0) else y for y in x.split(",")] if ("," in x) else x
                    for x in atom_dihedral_string
                ]
                atom_dihedral_string = [x for x in split_line[0].split("-") if (len(x) > 0 and not x == "-")]

                # print("Atoms-String: ", atom_dihedral_string)

                if len(atom_dihedral_string) == 4:
                    atomI, atomJ, atomK, atomL = atom_dihedral_string
                elif len(atom_dihedral_string) == 2:
                    atomI, atomJ, atomK, atomL = ["X"] + atom_dihedral_string + ["X"]
                elif len(atom_dihedral_string) == 1:
                    atomI, atomJ, atomK, atomL = ["X"] + atom_dihedral_string + ["X"] + ["X"]
                else:
                    # print(atom_dihedral_string)
                    atomI, atomJ, atomK, atomL = atom_dihedral_string + ["X"]

                concrete_example = (
                    table_lines[comment_index + 1] if (table_lines[comment_index + 1].startswith("#")) else ""
                )
                # print(concrete_example)
                params = dihedral_type(
                    ICP=int(ICP),
                    CP=float(CP),
                    PD=float(PD),
                    NP=float(NP),
                    atomI=atomI,
                    atomJ=atomJ,
                    atomK=atomK,
                    atomL=atomL,
                    special_number=int(special_number),
                    concrete_example=concrete_example,
                )
                self.content.append(params)

    def block_to_string(self) -> str:
        result = self.name + "\n"
        result += "#" + self.field_seperator + "NRPTY" + self.field_seperator + "NPTY" + self.line_seperator
        result += self.field_seperator + str(self.NRPTY) + self.field_seperator + str(self.NPTY) + self.line_seperator
        result += "# TABLE CONTENT: " + self.line_seperator
        result += "#" + self.field_seperator + self.field_seperator.join(self.table_header) + self.line_seperator
        for x in self.content:
            result += x.to_string()
        result += "END\n"
        return result


class IMPDIHEDRALTYPECODE(_iterable_topology_block):
    NRQTY: int
    NQTY: int
    table_header: Iterable[str] = ["ICQ", "CQ", "Q0"]

    def __init__(
        self,
        content: Union[Iterable[improper_dihedral_type], str],
        FORCEFIELD: FORCEFIELD = None,
        MAKETOPVERSION: MAKETOPVERSION = None,
        NRQTY: int = None,
        NQTY: int = None,
    ):
        """
                       GROMOS improper dihedral type parameters

        Parameters
        ----------
        content: Union[Iterable[atom_mass_type], str]
        FORCEFIELD : FORCEFIELD
        MAKETOPVERSION : MAKETOPVERSION
        NRQTY : int, optional
            Number of improperDihedrals types
        NQTY : int, optional
            Number of maximal improperDihedrals index.rst
        """
        super().__init__(FORCEFIELD=FORCEFIELD, MAKETOPVERSION=MAKETOPVERSION, content=content)

        if isinstance(content, list) and all([isinstance(x, str) for x in content]):
            self.read_content_from_str(content)
        elif isinstance(content, Iterable) and all([isinstance(x, atom_mass_type) for x in content]):
            self.content = content
        else:
            raise IOError("I don't understand the type of content: " + str(type(content)))

        # print(self.content)

        if NRQTY is None:
            self.NRQTY = len(self.content)
        elif isinstance(NRQTY, int):
            if NRQTY == len(self.content):  # CHECK FOR POSSIBLE ERROR
                self.NRQTY = NRQTY
            else:
                raise ValueError("IN MASSATOMTYPECODE NRQTY is not equal to the ammount of MASSATOMTYPES.")
        else:
            raise IOError("I don't understand the type of NRQTY: " + str(type(NRQTY)))

        if NQTY is None:
            self.NQTY = max([x.ICQ for x in self.content])
        elif isinstance(NQTY, int):
            if NQTY == max([x.ICQ for x in self.content]):  # CHECK FOR POSSIBLE ERROR
                self.NQTY = NQTY
            else:
                raise ValueError("IN MASSATOMTYPECODE NQTY is not the maximal Mass atom type code.")
        else:
            raise IOError("I don't understand the type of NQTY: " + str(type(NQTY)))

    def read_content_from_str(self, content: str):
        if isinstance(content, str):
            lines = content.split("\n")
        else:
            lines = content

        table_start = 0
        for line in lines:
            table_start += 1
            if all([field in line for field in self.table_header]):
                break

        if table_start == len(lines):
            raise ValueError("Could not find the TABLE start in " + self.name)

        table_lines = lines[table_start:]
        for field in table_lines:
            if not field.strip().startswith("#") and not len(field.strip()) == 0:
                ICQ, CQ, Q0 = field.strip().split()

                # nasty next line formatting of gromos
                comment_index = table_lines.index(field) + 1
                clean = re.sub("\(.*?\)", "", table_lines[comment_index])  # noqa # flake8 ignore
                split_line = [x for x in clean.replace("#", "").replace("-", "").strip().split(" ") if (len(x) > 0)]
                split_line = [x.split(",") if ("," in x) else x for x in split_line]

                if len(split_line) == 3:
                    groupType, special_number = " ".join(split_line[:2]), split_line[2]
                elif len(split_line) == 2:
                    groupType = split_line
                    special_number = 0

                params = improper_dihedral_type(
                    ICQ=int(ICQ), CQ=float(CQ), Q0=float(Q0), group_type=groupType, special_number=int(special_number)
                )
                self.content.append(params)

    def block_to_string(self) -> str:
        result = self.name + "\n"
        result += "#" + self.field_seperator + "NRQTY" + self.field_seperator + "NQTY" + self.line_seperator
        result += self.field_seperator + str(self.NRQTY) + self.field_seperator + str(self.NQTY) + self.line_seperator
        result += "# TABLE CONTENT: " + self.line_seperator
        result += "#" + self.field_seperator + self.field_seperator.join(self.table_header) + self.line_seperator
        for x in self.content:
            result += x.to_string()
        result += "END\n"
        return result


class SINGLEATOMLJPAIR(_iterable_topology_block):
    NRATT: int
    table_header = ["IAC", "TYPE", "C6", "C12(1)", "C12(2)", "C12(3)"]

    def __init__(
        self,
        content: List[str],
        NRATT: int = None,
        FORCEFIELD: FORCEFIELD = None,
        MAKETOPVERSION: MAKETOPVERSION = None,
    ):
        """

        Parameters
        ----------
        content
        NRATT
        FORCEFIELD
        MAKETOPVERSION
        """
        super().__init__(FORCEFIELD=FORCEFIELD, MAKETOPVERSION=MAKETOPVERSION, content=content)

        if isinstance(content, list) and all([isinstance(x, str) for x in content]):
            self.read_content_from_str(content)
        elif isinstance(content, Iterable) and all([isinstance(x, atom_mass_type) for x in content]):
            self.content = content
        else:
            raise IOError("I don't understand the type of content: " + str(type(content)))

        if NRATT is None:
            self.NRATT = len(self.content)
        elif isinstance(NRATT, int):
            if NRATT == len(self.content):  # CHECK FOR POSSIBLE ERROR
                self.NRATT = NRATT
            else:
                raise ValueError("IN number of terms NRATT is not equal to the ammount of terms.")
        else:
            raise IOError("I don't understand the type of NRATT: " + str(type(NRATT)))

    def read_content_from_str(self, content: str):
        if isinstance(content, str):
            lines = content.split("\n")
        else:
            lines = content

        table_start = 0
        for line in lines:
            table_start += 1
            if all([field in line for field in self.table_header]):
                break

        if table_start == len(lines):
            raise ValueError("Could not find the TABLE start in " + self.name)

        table_lines = lines[table_start:]
        subblock_end = -1
        end_string = "#--"
        for subblock_start, field in enumerate(table_lines):
            if not field.strip().startswith("#") and not len(field.strip()) == 0 and subblock_start >= subblock_end:
                # print("START: ", field)
                IAC, TYPE, C6, C12_1, C12_2, C12_3 = field.strip().split()
                IAC, TYPE, C6, C12_1, C12_2, C12_3 = [
                    int(IAC),
                    str(TYPE),
                    float(C6),
                    float(C12_1),
                    float(C12_2),
                    float(C12_3),
                ]
                LJ14PAIR = []

                subblock_lines = 0
                secondCS = False
                for subblock_line in table_lines[subblock_start + 1 :]:
                    subblock_lines += 1
                    # print(subblock_line)
                    if end_string in subblock_line:
                        subblock_end = subblock_start + subblock_lines
                        break
                    elif not subblock_line.startswith("#") and not secondCS:
                        CS6, CS12 = list(map(float, subblock_line.strip().split()))
                        secondCS = True
                    elif not subblock_line.startswith("#") and subblock_lines > 2:
                        LJ14PAIR.append(list(map(int, subblock_line.split())))

                params = single_atom_lj_pair_type(
                    IAC=IAC,
                    TYPE=TYPE,
                    C6=C6,
                    C12_1=C12_1,
                    C12_2=C12_2,
                    C12_3=C12_3,
                    CS6=CS6,
                    CS12=CS12,
                    LJ14PAIR=LJ14PAIR,
                )
                # print(params)
                self.content.append(params)

    def block_to_string(self) -> str:
        result = self.name + "\n"
        result += "#" + self.field_seperator + "NRATT"
        result += self.field_seperator + str(self.NRATT) + self.line_seperator
        result += "# TABLE CONTENT: " + self.line_seperator
        result += "#" + self.field_seperator + self.field_seperator.join(self.table_header) + self.line_seperator
        for x in self.content:
            result += x.to_string()
        result += "END\n"
        return result


class MIXEDATOMLJPAIR(_iterable_topology_block):
    NRMTT: int
    table_header: Iterable[str] = ["IACI", "IACJ", "C6", "C12_1", "C12_2", "C12_3"]
    """
        GROMOS 43A1 normal van der Waals parameters for mixed atom type pairs (I,J)
    """

    def __init__(
        self,
        content: List[str],
        NRMTT: int = None,
        FORCEFIELD: FORCEFIELD = None,
        MAKETOPVERSION: MAKETOPVERSION = None,
    ):
        """

        Parameters
        ----------
        content
        NRATT
        FORCEFIELD
        MAKETOPVERSION
        """
        super().__init__(FORCEFIELD=FORCEFIELD, MAKETOPVERSION=MAKETOPVERSION, content=content)

        if isinstance(content, list) and all([isinstance(x, str) for x in content]):
            self.read_content_from_str(content)
        elif isinstance(content, Iterable) and all([isinstance(x, atom_mass_type) for x in content]):
            self.content = content
        else:
            raise IOError("I don't understand the type of content: " + str(type(content)))

        if NRMTT is None:
            self.NRMTT = len(self.content)
        elif isinstance(NRMTT, int):
            if NRMTT == len(self.content):  # CHECK FOR POSSIBLE ERROR
                self.NRMTT = NRMTT
            else:
                raise ValueError("IN number of terms NRMTT is not equal to the ammount of terms.")
        else:
            raise IOError("I don't understand the type of NRMTT: " + str(type(NRMTT)))

    def read_content_from_str(self, content: str):
        if isinstance(content, str):
            lines = content.split("\n")
        else:
            lines = content
        table_start = 0

        for line in lines:
            table_start += 1
            if all([field in line for field in self.table_header]):
                break

        if table_start == len(lines):  # as there is no table header in the file
            table_start = 0

        table_lines = lines[table_start:]
        for field in table_lines:
            if not field.startswith("#"):
                IACI, IACJ, C6, C12_1, C12_2, C12_3 = field.strip().split()
                self.content.append(
                    mixed_atom_lj_pair_type(
                        IACI=int(IACI),
                        IACJ=int(IACJ),
                        C6=float(C6),
                        C12_1=float(C12_1),
                        C12_2=float(C12_2),
                        C12_3=float(C12_3),
                    )
                )

    def block_to_string(self) -> str:
        result = self.name + "\n"
        result += "#" + self.field_seperator + "NRMTT" + self.line_seperator
        result += self.field_seperator + str(self.NRMTT) + self.line_seperator
        result += "# TABLE CONTENT: " + self.line_seperator
        result += "#" + self.field_seperator + self.field_seperator.join(self.table_header) + self.line_seperator
        for x in self.content:
            result += x.to_string()
        result += "END\n"
        return result


class SPECATOMLJPAIR(_iterable_topology_block):
    NRST: int
    table_header: Iterable[str] = ["???"]

    def __init__(
        self,
        content: List[str],
        NRST: int = None,
        FORCEFIELD: FORCEFIELD = None,
        MAKETOPVERSION: MAKETOPVERSION = None,
    ):
        """

        Parameters
        ----------
        content
        NRATT
        FORCEFIELD
        MAKETOPVERSION
        """
        super().__init__(FORCEFIELD=FORCEFIELD, MAKETOPVERSION=MAKETOPVERSION, content=content)

        if isinstance(content, list) and all([isinstance(x, str) for x in content]):
            self.read_content_from_str(content)
        elif isinstance(content, Iterable) and all([isinstance(x, atom_mass_type) for x in content]):
            self.content = content
        else:
            raise IOError("I don't understand the type of content: " + str(type(content)))

        if NRST is None:
            self.NRST = len(self.content)
        elif isinstance(NRST, int):
            if NRST == len(self.content):  # CHECK FOR POSSIBLE ERROR
                self.NRST = NRST
            else:
                raise ValueError("IN number of terms NRST is not equal to the ammount of terms.")
        else:
            raise IOError("I don't understand the type of NRST: " + str(type(NRST)))

    def read_content_from_str(self, content: str):
        # TODO: implement
        # warnings.warn("SPECIAL LJ BLOCK IS NOT IMPLEMENTED!")
        return []


###################################################################
# Top() class blocks
###################################################################


class soluteatom_type(_generic_field):
    def __init__(
        self,
        ATNM: int,
        MRES: int,
        PANM: str,
        IAC: int,
        MASS: float,
        CG: float,
        CGC: int,
        INE: int,
        INEvalues: List[int],
        INE14: int,
        INE14values: List[int],
    ):
        """soluteatom_type

        Parameters
        ----------
        ATNM : int
            [description]
        MRES : int
            [description]
        PANM : str
            [description]
        IAC : int
            [description]
        MASS : float
            [description]
        CG : float
            [description]
        CGC : int
            [description]
        INE : int
            [description]
        INEvalues : [type]
            [description]
        INE14 : int
            [description]
        INE14values : [type]
            [description]
        """
        self.ATNM = ATNM
        self.MRES = MRES
        self.PANM = PANM
        self.IAC = IAC
        self.MASS = MASS
        self.CG = CG
        self.CGC = CGC
        self.INE = INE
        self.INEvalues = INEvalues
        self.INE14 = INE14
        self.INE14values = INE14values

    def to_string(self):
        str_line = (
            "\t"
            + str(self.ATNM)
            + "\t"
            + str(self.MRES)
            + "\t"
            + str(self.PANM)
            + "\t"
            + str(self.IAC)
            + "\t"
            + str(self.MASS)
            + "\t"
            + str(self.CG)
            + "\t"
            + str(self.CGC)
            + "\t"
            + str(self.INE)
        )
        lcounter = 0
        temp_INE = len(self.INEvalues)
        for iter in self.INEvalues:
            str_line += "\t" + str(iter).strip()
            lcounter += 1
            if (lcounter % 6) == 0 and temp_INE > 6:
                str_line += "\n\t\t\t\t\t\t\t\t\t\t"
                temp_INE -= 6
        str_line += "\n\t\t\t\t\t\t\t\t\t\t" + str(self.INE14)
        for iter in self.INE14values:
            str_line += "\t" + str(iter)
        str_line += "\n"
        return str_line


class bondstretchtype_type(_generic_field):
    def __init__(self, CB: float, CHB: float, B0: float):
        """
               GROMOS bondstretchtype for a single pair

        Parameters
        ----------
        CB: float
            quartic force constant
        CHB: float
            harmonic force constant
        B0: float
            bond length at minimum energy
        """
        self.CB = CB
        self.CHB = CHB
        self.B0 = B0

    def to_string(self):
        str_line = (
            "\t" + "{:.5e}".format(self.CB) + "\t" + "{:.5e}".format(self.CHB) + "\t" + "{:.5e}".format(self.B0) + "\n"
        )
        return str_line


class bondanglebendtype_type(_generic_field):
    def __init__(self, CT: float, CHT: float, T0: float):
        """
               GROMOS bondanglebendtype for a single angle

        Parameters
        ----------
        CT: float
            quartic force constant
        CHT: float
            harmonic force constant
        T0: float
            angle at minimum energy
        """
        self.CT = CT
        self.CHT = CHT
        self.T0 = T0

    def to_string(self):
        str_line = (
            "\t" + "{:.5e}".format(self.CT) + "\t" + "{:.5e}".format(self.CHT) + "\t" + "{:.5e}".format(self.T0) + "\n"
        )
        return str_line


class top_bond_type(_generic_field):
    def __init__(self, IB: int, JB: int, ICB: int):
        """
               GROMOS bond definition for a single pair

        Parameters
        ----------
        IB: int
            Atom number i
        JB: int
            Atom number j
        ICB: int
            Bond type code
        """
        self.IB = IB
        self.JB = JB
        self.ICB = ICB

    def to_string(self):
        str_line = "\t" + str(self.IB) + "\t" + str(self.JB) + "\t" + str(self.ICB) + "\n"
        return str_line


class bondangle_type(_generic_field):
    def __init__(self, IT: int, JT: int, KT: int, ICT: int):
        """
               GROMOS bondangleype for a single pair

        Parameters
        ----------
        IT: int
            atom number in angle
        JT: int
            atom number in angle
        KT: int
            atom number in angle
        ICT: int
            bond angle type code
        """
        self.IT = IT
        self.JT = JT
        self.KT = KT
        self.ICT = ICT

    def to_string(self):
        str_line = "\t" + str(self.IT) + "\t" + str(self.JT) + "\t" + str(self.KT) + "\t" + str(self.ICT) + "\n"
        return str_line


class impdihedraltype_type(_generic_field):
    def __init__(self, CQ: float, Q0: float):
        """
               GROMOS impdihedraltype for a single pair

        Parameters
        ----------
        CQ: float
            force constant of improper dihedral per degrees square
        Q0: float
            improper dihedral angle at minimum energy in degrees
        """
        self.CQ = CQ
        self.Q0 = Q0

    def to_string(self):
        str_line = "\t" + "{:.5e}".format(self.CQ) + "\t" + "{:.5e}".format(self.Q0) + "\n"
        return str_line


class impdihedralh_type(_generic_field):
    def __init__(self, IQH: int, JQH: int, KQH: int, LQH: int, ICQH: int):
        """
               GROMOS impdihedralH for a single pair

        Parameters
        ----------
        IQH:int
        JQH:int
        KQH:int
        LQH:int
            IQH,JQH,KQH,LQH: atom sequence numbers of atoms forming an improper dihedral
        ICQH:int
            improper dihedral type code

        """
        self.IQH = IQH
        self.JQH = JQH
        self.KQH = KQH
        self.LQH = LQH
        self.ICQH = ICQH

    def to_string(self):
        str_line = (
            "\t"
            + str(self.IQH)
            + "\t"
            + str(self.JQH)
            + "\t"
            + str(self.KQH)
            + "\t"
            + str(self.LQH)
            + "\t"
            + str(self.ICQH)
            + "\n"
        )
        return str_line


class impdihedral_type(_generic_field):
    def __init__(self, IQ: int, JQ: int, KQ: int, LQ: int, ICQ: int):
        """
               GROMOS impdihedral for a single pair

        Parameters
        ----------
        IQ:int
        JQ:int
        KQ:int
        LQ:int
            IQ,JQ,KQ,LQ: atom sequence numbers of atoms forming an improper dihedral
        ICQ:int
            improper dihedral type code

        """
        self.IQ = IQ
        self.JQ = JQ
        self.KQ = KQ
        self.LQ = LQ
        self.ICQ = ICQ

    def to_string(self):
        str_line = (
            "\t"
            + str(self.IQ)
            + "\t"
            + str(self.JQ)
            + "\t"
            + str(self.KQ)
            + "\t"
            + str(self.LQ)
            + "\t"
            + str(self.ICQ)
            + "\n"
        )
        return str_line


class torsdihedraltype_type(_generic_field):
    def __init__(self, CP: float, PD: float, NP: int):
        """
               GROMOS dihedraltype for a single pair

        Parameters
        ----------
        CP:float
            force constant
        PD:float
            phase-shift angle
        NP:int
            multiplicity
        """
        self.CP = CP
        self.PD = PD
        self.NP = NP

    def to_string(self):
        str_line = "\t" + "{:.5e}".format(self.CP) + "\t" + "{:.5e}".format(self.PD) + "\t" + str(self.NP) + "\n"
        return str_line


class dihedralh_type(_generic_field):
    def __init__(self, IPH: int, JPH: int, KPH: int, LPH: int, ICPH: int):
        """
               GROMOS dihedral for a single pair

        Parameters
        ----------
        """
        self.IPH = IPH
        self.JPH = JPH
        self.KPH = KPH
        self.LPH = LPH
        self.ICPH = ICPH

    def to_string(self):
        str_line = (
            "\t"
            + str(self.IPH)
            + "\t"
            + str(self.JPH)
            + "\t"
            + str(self.KPH)
            + "\t"
            + str(self.LPH)
            + "\t"
            + str(self.ICPH)
            + "\n"
        )
        return str_line


class top_dihedral_type(_generic_field):
    def __init__(self, IP: int, JP: int, KP: int, LP: int, ICP: int):
        """
               GROMOS dihedral for a single pair

        Parameters
        ----------
        """
        self.IP = IP
        self.JP = JP
        self.KP = KP
        self.LP = LP
        self.ICP = ICP

    def to_string(self):
        str_line = (
            "\t"
            + str(self.IP)
            + "\t"
            + str(self.JP)
            + "\t"
            + str(self.KP)
            + "\t"
            + str(self.LP)
            + "\t"
            + str(self.ICP)
            + "\n"
        )
        return str_line


class crossgihedralh_type(_generic_field):
    def __init__(self, APH: int, BPH: int, CPH: int, DPH: int, EPH: int, FPH: int, GPH: int, HPH: int, ICCH: int):
        """
        GROMOS Cross Dihedral type for H

        Parameters
        ----------
        APH : int
            number of atoms forming a dihedral
        BPH : int
            number of atoms forming a dihedral
        CPH : int
            number of atoms forming a dihedral
        DPH : int
            number of atoms forming a dihedral
        EPH : int
            number of atoms forming a dihedral
        FPH : int
            number of atoms forming a dihedral
        GPH : int
            number of atoms forming a dihedral
        HPH : int
            number of atoms forming a dihedral
        ICCH : int
            dihedral type code
        """

        self.APH = APH
        self.BPH = BPH
        self.CPH = CPH
        self.DPH = DPH
        self.EPH = EPH
        self.FPH = FPH
        self.GPH = GPH
        self.HPH = HPH
        self.ICCH = ICCH

    def to_string(self):
        str_line = (
            "\t"
            + str(self.APH)
            + "\t"
            + str(self.BPH)
            + "\t"
            + str(self.CPH)
            + "\t"
            + str(self.DPH)
            + "\t"
            + str(self.EPH)
            + "\t"
            + str(self.FPH)
            + "\t"
            + str(self.GPH)
            + "\t"
            + str(self.HPH)
            + "\t"
            + str(self.ICCH)
            + "\n"
        )
        return str_line


class crossgihedral_type(_generic_field):
    def __init__(self, AP: int, BP: int, CP: int, DP: int, EP: int, FP: int, GP: int, HP: int, ICC: int):
        """
        GROMOS Cross Dihedral type for NON H Atoms

        Parameters
        ----------
        AP : int
            number of atoms forming a dihedral
        BP : int
            number of atoms forming a dihedral
        CP : int
            number of atoms forming a dihedral
        DP : int
            number of atoms forming a dihedral
        EP : int
            number of atoms forming a dihedral
        FP : int
            number of atoms forming a dihedral
        GP : int
            number of atoms forming a dihedral
        HP : int
            number of atoms forming a dihedral
        ICC : int
            dihedral type code
        """

        self.AP = AP
        self.BP = BP
        self.CP = CP
        self.DP = DP
        self.EP = EP
        self.FP = FP
        self.GP = GP
        self.HP = HP
        self.ICC = ICC

    def to_string(self):
        str_line = (
            "\t"
            + str(self.AP)
            + "\t"
            + str(self.BP)
            + "\t"
            + str(self.CP)
            + "\t"
            + str(self.DP)
            + "\t"
            + str(self.EP)
            + "\t"
            + str(self.FP)
            + "\t"
            + str(self.GP)
            + "\t"
            + str(self.HP)
            + "\t"
            + str(self.ICC)
            + "\n"
        )
        return str_line


class ljparameters_type(_generic_field):
    def __init__(self, IAC: int, JAC: int, C12: float, C6: float, CS12: float, CS6: float):
        """
               GROMOS LJ parameter pair

        Parameters
        ----------
        """
        self.IAC = IAC
        self.JAC = JAC
        self.C12 = C12
        self.C6 = C6
        self.CS12 = CS12
        self.CS6 = CS6

    def to_string(self):
        str_line = (
            "\t"
            + str(self.IAC)
            + "\t"
            + str(self.JAC)
            + "\t"
            + "{:.6e}".format(self.C12)
            + "\t"
            + "{:.6e}".format(self.C6)
            + "\t"
            + "{:.6e}".format(self.CS12)
            + "\t"
            + "{:.6e}".format(self.CS6)
            + "\n"
        )
        return str_line


class ljexception_type(_generic_field):
    def __init__(self, AT1: int, AT2: int, C12: float, C6: float):
        """
               GROMOS LJ exception pair

        Parameters
        ----------
        """
        self.AT1 = AT1
        self.AT2 = AT2
        self.C12 = C12
        self.C6 = C6

    def to_string(self):
        str_line = (
            "\t"
            + str(self.AT1)
            + "\t"
            + str(self.AT2)
            + "\t"
            + "{:.6e}".format(self.C12)
            + "\t"
            + "{:.6e}".format(self.C6)
            + "\n"
        )
        return str_line


class solventatom_type(_generic_field):
    def __init__(self, I: int, ANMS: str, IACS: int, MASS: float, CGS: float):  # noqa: E741
        """
               GROMOS solventatom line

        Parameters
        ----------
        """
        self.I = I  # noqa: E741
        self.ANMS = ANMS
        self.IACS = IACS
        self.MASS = MASS
        self.CGS = CGS

    def to_string(self):
        str_line = (
            "\t"
            + str(self.I)
            + "\t"
            + str(self.ANMS)
            + "\t"
            + str(self.IACS)
            + "\t"
            + "{:.5e}".format(self.MASS)
            + "\t"
            + "{:.5e}".format(self.CGS)
            + "\n"
        )
        return str_line


class solventconstr_type(_generic_field):
    def __init__(self, ICONS: int, JCONS: int, CONS: float):
        """
               GROMOS SOLVENTCONSTR entry

        Parameters
        ----------
        """
        self.ICONS = ICONS
        self.JCONS = JCONS
        self.CONS = CONS

    def to_string(self):
        str_line = "\t" + str(self.ICONS) + "\t" + str(self.JCONS) + "\t" + str(self.CONS) + "\n"
        return str_line


class constraint_type(_generic_field):
    def __init__(self, IC: int, JC: int, ICC: float):
        """[summary]

        Parameters
        ----------
        IC : int
            [description]
        JC : int
            [description]
        ICC : float
            [description]
        """
        self.IC = IC
        self.JC = JC
        self.ICC = ICC

    def to_string(self):
        str_line = (
            self.fieldseperator
            + str(self.IC)
            + self.fieldseperator
            + str(self.JC)
            + self.fieldseperator
            + str(self.ICC)
            + self.lineseperator
        )
        return str_line


class _topology_table_block(_iterable_topology_block):
    table_header: Iterable[str]
    table_line_type = _generic_field

    def __init__(
        self,
        content: Union[str, dict, _topology_table_block_Type],
        FORCEFIELD: FORCEFIELD = None,
        MAKETOPVERSION: MAKETOPVERSION = None,
        **kwargs
    ):
        """
        Parent class for all table like GROMOS topolgy blocks.
        Offers a standard implementation of read and write functions

        Requirements for child classes:
            table_header and table_line_type need to be defined
            additional attributes can be passed via kwargs

        Parameters
        ----------
        content : str or dict or None
            content of the table like GROMOS block
        FORCEFIELD : FORCEFIELD, optional
            [description], by default None
        MAKETOPVERSION : MAKETOPVERSION, optional
            [description], by default None
        """
        # set attributes
        for attribute in vars(__class__):
            if attribute.isupper():  # is_allcapital()):
                setattr(self, attribute, kwargs[attribute])

        # init _iterable_topology_block
        super().__init__(FORCEFIELD=FORCEFIELD, MAKETOPVERSION=MAKETOPVERSION, content=content)

        if len(kwargs.keys()) == 1:
            for key, value in kwargs.items():
                if value is None:
                    setattr(self, key, len(self.content))
                elif isinstance(value, int):
                    if value == len(self.content):  # CHECK FOR POSSIBLE ERROR
                        setattr(self, key, value)
                    else:
                        raise ValueError("In " + self.name + " is " + str(key) + " not equal to the ammount.")
                else:
                    raise IOError("I don't understand the type of " + str(key) + ": " + str(type(key)))

    def _check_import_method(self, content: str = None):
        if isinstance(content, list) and all([isinstance(x, str) for x in content]):
            self.read_content_from_str(content)
        elif isinstance(content, Iterable) and all([isinstance(x, type(self.table_line_type)) for x in content]):
            self.content = content
        elif isinstance(content, str):
            self.read_content_from_str(content.split(self.line_seperator))
        else:
            raise IOError("I don't understand the type of content: " + str(type(content)))

    def read_content_from_str(self, content: str):
        if not hasattr(self, "table_header"):
            raise Exception("Could not find table_header of " + self.name)
        elif not hasattr(self, "table_line_type"):
            raise Exception("Could not find table_line_type of " + self.name)
        else:
            if isinstance(content, str):
                lines = content.split("\n")
            else:
                lines = content
            # Table Reading:
            table_start = 0
            for line in lines:
                table_start += 1
                if all([field in line for field in self.table_header]):
                    break

            if table_start > len(lines):
                raise ValueError("Could not find the TABLE start in " + self.name)
            else:
                self._read_table(lines[table_start:])

    def _read_table(self, table: str):
        # get information of the needed structure of the sub class
        signature = inspect.signature(self.table_line_type.__init__)
        parameter_name = [name for name in signature.parameters if (name != "self")]
        parameter_type = {name: signature.parameters[name].annotation for name in parameter_name}
        # loop over the table (=content) to create the sub classes (=table_line_type)
        for table_line in table:
            if not (table_line.startswith("#") or len(table_line) == 0):
                # pre parse all non-empty non-comment lines into a list of strings
                fields = table_line.strip().split()
                if len(fields) != len(parameter_name):
                    raise Exception(
                        "Fields are not matching the ammount of needed arguments!\n "
                        "#fileds: " + str(len(fields)) + "\t#args: " + str(len(parameter_name)) + "\n\n "
                        "require: " + str(parameter_name) + "\t got: " + str(fields)
                    )
                # generate arguments dict for line parsing (= table_line_type class construction)
                kwargs = {key: parameter_type[key](field) for key, field in zip(parameter_name, fields)}
                self.content.append(self.table_line_type(**kwargs))

    def block_to_string(self) -> str:
        result = "#" + self.field_seperator + self.field_seperator.join(self.table_header) + self.line_seperator
        for x in self.content:
            result += x.to_string()
        return result


class PHYSICALCONSTANTS(_topology_block):
    def __init__(
        self,
        content: Union[str, dict, PHYSICALCONSTANTS_Type],
        FORCEFIELD: FORCEFIELD = None,
        MAKETOPVERSION: MAKETOPVERSION = None,
    ):
        # Default definition of physical constants
        self.FPEPSI = 138.9354
        self.HBAR = 0.0635078
        self.SPDL = 299792.458
        self.BOLTZ = 0.00831441
        super().__init__(FORCEFIELD=FORCEFIELD, MAKETOPVERSION=MAKETOPVERSION, content=content)

    def _check_import_method(self, content: str = None):
        # elif (type(content) == __class__):
        #    self.content = content
        if content == [[""]] or content == [""] or content is None:
            self.content = [self.FPEPSI, self.HBAR, self.SPDL, self.BOLTZ]
        elif isinstance(content, list) and all([isinstance(x, str) for x in content]):
            self.read_content_from_str(content)
        elif isinstance(content, str):
            self.read_content_from_str(content.split(self.line_seperator))
        elif isinstance(content, tuple) and len(content) == 4:
            self.FPEPSI, self.HBAR, self.SPDL, self.BOLTZ = content
            self.content = [self.FPEPSI, self.HBAR, self.SPDL, self.BOLTZ]
        else:
            raise IOError("I don't understand the type of content: " + str(type(content)))

    def read_content_from_str(self, content: str):
        if isinstance(content, str):
            lines = content.split("\n")
        else:
            lines = content

        stash = []
        for field in lines:
            if not field.strip().startswith("#") and not len(field.strip()) == 0:
                stash.append(float(field.strip()))
        if len(stash) >= 4:
            self.FPEPSI = stash[0]
            self.HBAR = stash[1]
            self.SPDL = stash[2]
            self.BOLTZ = stash[3]
            self.content = (self.FPEPSI, self.HBAR, self.SPDL, self.BOLTZ)
        elif len(stash) == 0:
            self.content = (self.FPEPSI, self.HBAR, self.SPDL, self.BOLTZ)
        else:
            raise IOError("Not enough arguments provided in PHYSICALCONSTANTS")

    def block_to_string(self) -> str:
        result = self.name + "\n"
        result += "# FPEPSI: 1.0/(4.0*PI*EPS0) (EPS0 is the permittivity of vacuum)" + self.line_seperator
        result += str(self.FPEPSI) + self.line_seperator
        result += "# HBAR: Planck's constant HBAR = H/(2* PI)" + self.line_seperator
        result += str(self.HBAR) + self.line_seperator
        result += "# SPDL: Speed of light (nm/ps)" + self.line_seperator
        result += str(self.SPDL) + self.line_seperator
        result += "# BOLTZ: Boltzmann's constant kB" + self.line_seperator
        result += str(self.BOLTZ) + self.line_seperator
        result += "END\n"
        return result


class TOPVERSION(_topology_block):
    def __init__(
        self,
        content: Union[str, Dict[str, str], TOPVERSION_Type],
        FORCEFIELD: FORCEFIELD = None,
        MAKETOPVERSION: MAKETOPVERSION = None,
    ):
        super().__init__(FORCEFIELD=FORCEFIELD, MAKETOPVERSION=MAKETOPVERSION, content=content)


class ATOMTYPENAME(_topology_block):
    def __init__(
        self,
        content: Union[str, Dict[str, str], ATOMTYPENAME_Type],
        FORCEFIELD: FORCEFIELD = None,
        MAKETOPVERSION: MAKETOPVERSION = None,
    ):
        super().__init__(FORCEFIELD=FORCEFIELD, MAKETOPVERSION=MAKETOPVERSION, content=content)


class RESNAME(_topology_block):
    def __init__(
        self,
        content: Union[str, Dict[str, str], RESNAME_Type],
        FORCEFIELD: FORCEFIELD = None,
        MAKETOPVERSION: MAKETOPVERSION = None,
    ):
        super().__init__(FORCEFIELD=FORCEFIELD, MAKETOPVERSION=MAKETOPVERSION, content=content)


class SOLUTEATOM(_iterable_topology_block):
    NRP: int
    table_header: Iterable[str] = ["IB", "JB", "ICB"]

    def __init__(
        self,
        content: Union[str, Dict[str, str], _iterable_topology_block_Type],
        FORCEFIELD: FORCEFIELD = None,
        MAKETOPVERSION: MAKETOPVERSION = None,
    ):
        super().__init__(FORCEFIELD=FORCEFIELD, MAKETOPVERSION=MAKETOPVERSION, content=content)

    def _check_import_method(self, content=None):
        if isinstance(content, list) and all([isinstance(x, str) for x in content]):
            self.read_content_from_str(content)
        elif type(content) == __class__:
            self.content = content
        elif content is __class__:
            self.content = content

    def read_content_from_str(self, content: str):
        contentLines = []
        if isinstance(content, str):
            lines = content.split("\n")
        else:
            lines = content

        for field in lines:
            if not field.strip().startswith("#") and not len(field.strip()) == 0:
                contentLines.append(field)

        # set NRP and check for sanity
        try:
            self.NRP = int(contentLines.pop(0))
        except Exception:
            self.NRP = 0
        if self.NRP < 0:
            raise IOError("NPR in SOLUTEATOM Block is " + str(self.NRP))
        elif self.NRP == 0:
            content = None
        else:
            for _ in range(self.NRP):  # main import loop
                dump1 = contentLines.pop(0).strip().split()
                if len(dump1) < 8:
                    raise IOError("Not enough arguments provided in SOLUTEATOM Block")
                else:
                    ATNM, MRES, PANM, IAC, MASS, CG, CGC, INE = dump1[0:8]
                    if 1 <= int(INE):
                        INEvalues = [int(i) for i in dump1[8:]]
                        # keep reading in lines until we have all the data needed.
                        while int(INE) > len(INEvalues):
                            try:
                                if len(contentLines) == 0:
                                    raise IOError(
                                        "Not enough lines provided for multi line INE in SOLUTEATOM Block\nATNM="
                                        + str(ATNM)
                                        + " MRES="
                                        + str(MRES)
                                    )
                                elif any(i in contentLines[0] for i in ["\t\t\t\t\t", "                   "]):
                                    INEvalues.extend([int(i) for i in contentLines.pop(0).strip().split()])
                                else:
                                    raise IOError(
                                        "no intendation detected for mult line INE in SOLUTEATOM Block or too large number for INE\nATNM="
                                        + str(ATNM)
                                        + " MRES="
                                        + str(MRES)
                                    )
                            except IOError:
                                raise IOError("Problem reading INE for ATNM=" + str(ATNM) + " MRES=" + str(MRES))
                    else:
                        INEvalues = []

                dump2 = contentLines.pop(0).strip().split()
                if len(dump2) < 1:
                    raise IOError(
                        "Not enough arguments provided in SOLUTEATOM Block\nATNM=" + str(ATNM) + " MRES=" + str(MRES)
                    )
                else:
                    INE14 = dump2[0]
                    if 1 <= int(INE14):
                        INE14values = [int(i) for i in dump2[1:]]
                        # keep reading in lines until we have all the data needed.
                        while int(INE14) > len(INE14values):
                            try:
                                if len(contentLines) == 0:
                                    raise IOError(
                                        "Not enough lines provided for multi line INE14 in SOLUTEATOM Block\nATNM="
                                        + str(ATNM)
                                        + " MRES="
                                        + str(MRES)
                                    )
                                elif any(i in contentLines[0] for i in ["\t\t\t\t\t", "                   "]):
                                    INE14values.extend([int(i) for i in contentLines.pop(0).strip().split()])
                                else:
                                    raise IOError(
                                        "no intendation detected for mult line INE14 in SOLUTEATOM Block or too large number for INE14\nATNM="
                                        + str(ATNM)
                                        + " MRES="
                                        + str(MRES)
                                    )
                            except IOError:
                                raise IOError("Problem reading INE14 for ATNM=" + str(ATNM) + " MRES=" + str(MRES))
                    else:
                        INE14values = []
                # pass everything to the subclass maker
                params = soluteatom_type(
                    int(ATNM),
                    int(MRES),
                    PANM,
                    int(IAC),
                    float(MASS),
                    float(CG),
                    int(CGC),
                    int(INE),
                    INEvalues,
                    int(INE14),
                    INE14values,
                )
                self.content.append(params)

    def block_to_string(self) -> str:
        result = self.name + self.line_seperator
        result += "#   NRP: number of solute atoms" + self.line_seperator
        result += self.field_seperator + str(self.NRP) + self.line_seperator
        result += "#  ATNM: atom number\n#  MRES: residue number\n#  PANM: atom name of solute atom\n#   IAC: integer (van der Waals) atom type code\n#  MASS: mass of solute atom\n#    CG: charge of solute atom\n#   CGC: charge group code (0 or 1)\n#   INE: number of excluded atoms\n# INE14: number of 1-4 interactions\n# ATNM MRES PANM IAC     MASS       CG  CGC INE\n#                                           INE14\n"
        for iter in self.content:
            result += iter.to_string()
        result += "END\n"
        return result


class BONDSTRETCHTYPE(_topology_table_block):
    NBTY: int
    table_header: Iterable[str] = ["CB", "CHB", "B0"]
    table_line_type = bondstretchtype_type

    def __init__(
        self,
        content: Union[Iterable[bondstretchtype_type], str],
        FORCEFIELD: FORCEFIELD = None,
        MAKETOPVERSION: MAKETOPVERSION = None,
        NBTY: int = None,
    ):
        """
                       GROMOS BONDSTRETCHTYPE block

        Parameters
        ----------
        content: Union[Iterable[bondstretchtype_type], str]
        FORCEFIELD : FORCEFIELD
        MAKETOPVERSION : MAKETOPVERSION
        NBTY : int, optional
            Number of bondstretchtypes
        """
        kwargs = {"NBTY": NBTY}
        super().__init__(content=content, FORCEFIELD=FORCEFIELD, MAKETOPVERSION=MAKETOPVERSION, **kwargs)

    def read_content_from_str(self, content: str):
        super().read_content_from_str(content)

    def block_to_string(self) -> str:
        result = self.name + "\n"
        result += "#" + self.field_seperator + "NBTY: number of covalent bond types" + self.line_seperator
        result += self.field_seperator + str(self.NBTY) + self.line_seperator
        result += super().block_to_string()
        result += "END\n"
        return result


class BOND(_topology_table_block):
    NBON: int = 1
    table_header: Iterable[str] = ["IB", "JB", "ICB"]
    table_line_type = top_bond_type

    def __init__(
        self,
        content: Union[Iterable[top_bond_type], str],
        FORCEFIELD: FORCEFIELD = None,
        MAKETOPVERSION: MAKETOPVERSION = None,
        NBON: int = None,
    ):
        """
                       GROMOS BOND block

        Parameters
        ----------
        content: Union[Iterable[top_bond_type], str]
        FORCEFIELD : FORCEFIELD
        MAKETOPVERSION : MAKETOPVERSION
        NBON : int, optional
            Number of bonds
        """
        kwargs = {"NBON": NBON}
        super().__init__(content=content, FORCEFIELD=FORCEFIELD, MAKETOPVERSION=MAKETOPVERSION, **kwargs)

    def read_content_from_str(self, content: str):
        super().read_content_from_str(content)

    def block_to_string(self) -> str:
        result = self.name + "\n"
        result += (
            "#" + self.field_seperator + "NBON: number of bonds NOT involving H atoms in solute" + self.line_seperator
        )
        result += self.field_seperator + str(self.NBON) + self.line_seperator
        result += super().block_to_string()
        result += "END\n"
        return result


class BONDH(_topology_table_block):
    NBONH: int = 1
    table_header: Iterable[str] = ["IBH", "JBH", "ICBH"]
    table_line_type = top_bond_type  # reused data type for simplicity

    def __init__(
        self,
        content: Union[Iterable[top_bond_type], str],
        FORCEFIELD: FORCEFIELD = None,
        MAKETOPVERSION: MAKETOPVERSION = None,
        NBONH: int = None,
    ):
        """
                       GROMOS BONDH block

        Parameters
        ----------
        content: Union[Iterable[top_bond_type], str]
        FORCEFIELD : FORCEFIELD
        MAKETOPVERSION : MAKETOPVERSION
        NBONH : int, optional
            Number of bonds with H
        """
        kwargs = {"NBONH": NBONH}
        super().__init__(content=content, FORCEFIELD=FORCEFIELD, MAKETOPVERSION=MAKETOPVERSION, **kwargs)

    def read_content_from_str(self, content: str):
        super().read_content_from_str(content)

    def block_to_string(self) -> str:
        result = self.name + "\n"
        result += (
            "#" + self.field_seperator + "NBONH: number of bonds involving H atoms in solute" + self.line_seperator
        )
        result += self.field_seperator + str(self.NBONH) + self.line_seperator
        result += super().block_to_string()
        result += "END\n"
        return result


class BONDANGLEBENDTYPE(_topology_table_block):
    NBTY: int
    table_header: Iterable[str] = ["CT", "CHT", "T0"]
    table_line_type = bondanglebendtype_type

    def __init__(
        self,
        content: Union[Iterable[bondanglebendtype_type], str],
        FORCEFIELD: FORCEFIELD = None,
        MAKETOPVERSION: MAKETOPVERSION = None,
        NBTY: int = None,
    ):
        """
                       GROMOS BONDSTRETCHTYPE block

        Parameters
        ----------
        content: Union[Iterable[bondanglebendtype_type], str]
        FORCEFIELD : FORCEFIELD
        MAKETOPVERSION : MAKETOPVERSION
        NBTY : int, optional
            Number of bondstretchtypes
        """
        kwargs = {"NBTY": NBTY}
        super().__init__(content=content, FORCEFIELD=FORCEFIELD, MAKETOPVERSION=MAKETOPVERSION, **kwargs)

    def read_content_from_str(self, content: str):
        super().read_content_from_str(content)

    def block_to_string(self) -> str:
        result = self.name + "\n"
        result += "#" + self.field_seperator + "NBTY: number of angle types" + self.line_seperator
        result += self.field_seperator + str(self.NBTY) + self.line_seperator
        result += super().block_to_string()
        result += "END\n"
        return result


class BONDANGLE(_topology_table_block):
    NTHE: int
    table_header: Iterable[str] = ["IT", "JT", "KT", "ICT"]
    table_line_type = bondangle_type

    def __init__(
        self,
        content: Union[Iterable[bondangle_type], str],
        FORCEFIELD: FORCEFIELD = None,
        MAKETOPVERSION: MAKETOPVERSION = None,
        NTHE: int = None,
    ):
        """
                       GROMOS BONDSTRETCHTYPE block

        Parameters
        ----------
        content: Union[Iterable[bondangle_type], str]
        FORCEFIELD : FORCEFIELD
        MAKETOPVERSION : MAKETOPVERSION
        NTHE : int, optional
            Number of bondangles
        """
        kwargs = {"NTHE": NTHE}
        super().__init__(content=content, FORCEFIELD=FORCEFIELD, MAKETOPVERSION=MAKETOPVERSION, **kwargs)

    def read_content_from_str(self, content: str):
        super().read_content_from_str(content)

    def block_to_string(self) -> str:
        result = self.name + "\n"
        result += "#" + self.field_seperator + "NTHE: number of angles" + self.line_seperator
        result += self.field_seperator + str(self.NTHE) + self.line_seperator
        result += super().block_to_string()
        result += "END\n"
        return result


class BONDANGLEH(_topology_table_block):
    NTHEH: int
    table_header: Iterable[str] = ["ITH", "JTH", "KTH", "ICTH"]
    table_line_type = bondangle_type

    def __init__(
        self,
        content: Union[Iterable[bondangle_type], str],
        FORCEFIELD: FORCEFIELD = None,
        MAKETOPVERSION: MAKETOPVERSION = None,
        NTHEH: int = None,
    ):
        """
                       GROMOS BONDANGLEH block

        Parameters
        ----------
        content: Union[Iterable[bondangle_type], str]
        FORCEFIELD : FORCEFIELD
        MAKETOPVERSION : MAKETOPVERSION
        NTHEH : int, optional
            Number of bondangles
        """
        kwargs = {"NTHEH": NTHEH}
        super().__init__(content=content, FORCEFIELD=FORCEFIELD, MAKETOPVERSION=MAKETOPVERSION, **kwargs)

    def read_content_from_str(self, content: str):
        super().read_content_from_str(content)

    def block_to_string(self) -> str:
        result = self.name + "\n"
        result += "#" + self.field_seperator + "NTHEH: number of bondangles involving a H" + self.line_seperator
        result += self.field_seperator + str(self.NTHEH) + self.line_seperator
        result += super().block_to_string()
        result += "END\n"
        return result


class IMPDIHEDRALTYPE(_topology_table_block):
    NQTY: int
    table_header: Iterable[str] = ["CQ", "Q0"]
    table_line_type = impdihedraltype_type

    def __init__(
        self,
        content: Union[Iterable[impdihedraltype_type], str],
        FORCEFIELD: FORCEFIELD = None,
        MAKETOPVERSION: MAKETOPVERSION = None,
        NQTY: int = None,
    ):
        """
                       GROMOS IMPDIHEDRALTYPE block

        Parameters
        ----------
        content: Union[Iterable[impdihedraltype_type], str]
        FORCEFIELD : FORCEFIELD
        MAKETOPVERSION : MAKETOPVERSION
        NQTY : int, optional
            Number of impdihedraltype
        """
        kwargs = {"NQTY": NQTY}
        super().__init__(content=content, FORCEFIELD=FORCEFIELD, MAKETOPVERSION=MAKETOPVERSION, **kwargs)

    def read_content_from_str(self, content: str):
        super().read_content_from_str(content)

    def block_to_string(self) -> str:
        result = self.name + "\n"
        result += "#" + self.field_seperator + "NQTY: number of improper dihedrals" + self.line_seperator
        result += self.field_seperator + str(self.NQTY) + self.line_seperator
        result += super().block_to_string()
        result += "END\n"
        return result


class IMPDIHEDRALH(_topology_table_block):
    NQHIH: int
    table_header: Iterable[str] = ["IQH", "JQH", "KQH", "LQH", "ICQH"]
    table_line_type = impdihedralh_type

    def __init__(
        self,
        content: Union[Iterable[impdihedralh_type], str],
        FORCEFIELD: FORCEFIELD = None,
        MAKETOPVERSION: MAKETOPVERSION = None,
        NQHIH: int = None,
    ):
        """
                       GROMOS IMPDIHEDRALH block

        Parameters
        ----------
        content: Union[Iterable[impdihedralh_type], str]
        FORCEFIELD : FORCEFIELD
        MAKETOPVERSION : MAKETOPVERSION
        NQHIH : int, optional
            Number of impdihedralH
        """
        kwargs = {"NQHIH": NQHIH}
        super().__init__(content=content, FORCEFIELD=FORCEFIELD, MAKETOPVERSION=MAKETOPVERSION, **kwargs)

    def read_content_from_str(self, content: str):
        super().read_content_from_str(content)

    def block_to_string(self) -> str:
        result = self.name + "\n"
        result += (
            "#" + self.field_seperator + "NQHIH: number of improper dihedrals involving H atoms" + self.line_seperator
        )
        result += self.field_seperator + str(self.NQHIH) + self.line_seperator
        result += super().block_to_string()
        result += "END\n"
        return result


class IMPDIHEDRAL(_topology_table_block):
    NQHI: int
    table_header: Iterable[str] = ["IQ", "JQ", "KQ", "LQ", "ICQ"]
    table_line_type = impdihedral_type

    def __init__(
        self,
        content: Union[Iterable[impdihedral_type], str],
        FORCEFIELD: FORCEFIELD = None,
        MAKETOPVERSION: MAKETOPVERSION = None,
        NQHI: int = None,
    ):
        """
                       GROMOS IMPDIHEDRAL block

        Parameters
        ----------
        content: Union[Iterable[impdihedral_type], str]
        FORCEFIELD : FORCEFIELD
        MAKETOPVERSION : MAKETOPVERSION
        NQHI : int, optional
            Number of impdihedral
        """
        kwargs = {"NQHI": NQHI}
        super().__init__(content=content, FORCEFIELD=FORCEFIELD, MAKETOPVERSION=MAKETOPVERSION, **kwargs)

    def read_content_from_str(self, content: str):
        super().read_content_from_str(content)

    def block_to_string(self) -> str:
        result = self.name + "\n"
        result += (
            "#"
            + self.field_seperator
            + "NQHI: number of improper dihedrals NOT involving H atoms"
            + self.line_seperator
        )
        result += self.field_seperator + str(self.NQHI) + self.line_seperator
        result += super().block_to_string()
        result += "END\n"
        return result

    # def block_to_string(self) -> str:
    #     result = self.name + "\n"
    #     result += "#" + self.field_seperator + "#  NQHI: number of improper dihedrals" + self.line_seperator
    #     result += self.field_seperator + str(self.NQHI) + self.line_seperator
    #     result += "#" + self.field_seperator + self.field_seperator.join(self.table_header) + self.line_seperator
    #     for x in self.content:
    #         result += x.to_string()
    #     result += "END\n"
    #     return result


class TORSDIHEDRALTYPE(_topology_table_block):
    NPTY: int
    table_header: Iterable[str] = ["CP", "PD", "NP"]
    table_line_type = torsdihedraltype_type

    def __init__(
        self,
        content: Union[Iterable[torsdihedraltype_type], str],
        FORCEFIELD: FORCEFIELD = None,
        MAKETOPVERSION: MAKETOPVERSION = None,
        NPTY: int = None,
    ):
        """
                       GROMOS IMPDIHEDRAL block

        Parameters
        ----------
        content: Union[Iterable[torsdihedraltype_type], str]
        FORCEFIELD : FORCEFIELD
        MAKETOPVERSION : MAKETOPVERSION
        NQTY : int, optional
            Number of torsion dihedrals
        """
        kwargs = {"NPTY": NPTY}
        super().__init__(content=content, FORCEFIELD=FORCEFIELD, MAKETOPVERSION=MAKETOPVERSION, **kwargs)

    def read_content_from_str(self, content: str):
        super().read_content_from_str(content)

    def block_to_string(self) -> str:
        result = self.name + "\n"
        result += "#" + self.field_seperator + "NPTY: number of torsion dihedrals" + self.line_seperator
        result += self.field_seperator + str(self.NPTY) + self.line_seperator
        result += super().block_to_string()
        result += "END\n"
        return result


class DIHEDRALH(_topology_table_block):
    NPHIH: int
    table_header: Iterable[str] = ["IPH", "JPH", "KPH", "LPH", "ICPH"]
    table_line_type = dihedralh_type

    def __init__(
        self,
        content: Union[Iterable[dihedralh_type], str],
        FORCEFIELD: FORCEFIELD = None,
        MAKETOPVERSION: MAKETOPVERSION = None,
        NPHIH: int = None,
    ):
        """
                       GROMOS DIHEDRAL block

        Parameters
        ----------
        content: Union[Iterable[dihedralh_type], str]
        FORCEFIELD : FORCEFIELD
        MAKETOPVERSION : MAKETOPVERSION
        NPHIH : int, optional
            Number of dihedralH
        """
        kwargs = {"NPHIH": NPHIH}
        super().__init__(content=content, FORCEFIELD=FORCEFIELD, MAKETOPVERSION=MAKETOPVERSION, **kwargs)

    def read_content_from_str(self, content: str):
        super().read_content_from_str(content)

    def block_to_string(self) -> str:
        result = self.name + "\n"
        result += (
            "#" + self.field_seperator + "NPHIH: number of torsion dihedrals involving H atoms" + self.line_seperator
        )
        result += self.field_seperator + str(self.NPHIH) + self.line_seperator
        result += super().block_to_string()
        result += "END\n"
        return result


class DIHEDRAL(_topology_table_block):
    NPHI: int
    table_header: Iterable[str] = ["IP", "JP", "KP", "LP", "ICP"]
    table_line_type = top_dihedral_type

    def __init__(
        self,
        content: Union[Iterable[top_dihedral_type], str],
        FORCEFIELD: FORCEFIELD = None,
        MAKETOPVERSION: MAKETOPVERSION = None,
        NPHI: int = None,
    ):
        """
                       GROMOS DIHEDRAL block

        Parameters
        ----------
        content: Union[Iterable[top_dihedral_type], str]
        FORCEFIELD : FORCEFIELD
        MAKETOPVERSION : MAKETOPVERSION
        NPHI : int, optional
            Number of tors dihedral
        """
        kwargs = {"NPHI": NPHI}
        super().__init__(content=content, FORCEFIELD=FORCEFIELD, MAKETOPVERSION=MAKETOPVERSION, **kwargs)

    def read_content_from_str(self, content: str):
        super().read_content_from_str(content)

    def block_to_string(self) -> str:
        result = self.name + "\n"
        result += "#" + self.field_seperator + "NPHI: number of dihedrals NOT involving H atoms" + self.line_seperator
        result += self.field_seperator + str(self.NPHI) + self.line_seperator
        result += super().block_to_string()
        result += "END\n"
        return result


class CROSSDIHEDRALH(_topology_table_block):
    NPHIH: int
    table_header: Iterable[str] = ["APH", "BPH", "CPH", "DPH", "EPH", "FPH", "GPH", "HPH", "ICCH"]
    table_line_type = crossgihedralh_type

    def __init__(
        self,
        content: Union[Iterable[crossgihedralh_type], str],
        FORCEFIELD: FORCEFIELD = None,
        MAKETOPVERSION: MAKETOPVERSION = None,
        NPHIH: int = None,
    ):
        """[summary]

        Parameters
        ----------
        content : Union[Iterable[crossgihedralh_type], str]
            [description]
        FORCEFIELD : FORCEFIELD, optional
            [description], by default None
        MAKETOPVERSION : MAKETOPVERSION, optional
            [description], by default None
        NPHIH : [type], optional
            number of cross dihedrals involving H atoms in solute, by default None
        """
        kwargs = {"NPHIH": NPHIH}
        super().__init__(content=content, FORCEFIELD=FORCEFIELD, MAKETOPVERSION=MAKETOPVERSION, **kwargs)

    def read_content_from_str(self, content: str):
        super().read_content_from_str(content)

    def block_to_string(self) -> str:
        result = self.name + "\n"
        result += "#" + self.field_seperator + "NPHIH: number of dihedrals involving H atoms" + self.line_seperator
        result += self.field_seperator + str(self.NPHIH) + self.line_seperator
        result += super().block_to_string()
        result += "END\n"
        return result


class CROSSDIHEDRAL(_topology_table_block):
    NPHI: int
    table_header: Iterable[str] = ["AP", "BP", "CP", "DP", "EP", "FP", "GP", "HP", "ICC"]
    table_line_type = crossgihedralh_type

    def __init__(
        self,
        content: Union[Iterable[crossgihedral_type], str],
        FORCEFIELD: FORCEFIELD = None,
        MAKETOPVERSION: MAKETOPVERSION = None,
        NPHI: int = None,
    ):
        """[summary]

        Parameters
        ----------
        content : Union[Iterable[crossgihedral_type], str]
            [description]
        FORCEFIELD : FORCEFIELD, optional
            [description], by default None
        MAKETOPVERSION : MAKETOPVERSION, optional
            [description], by default None
        NPHI : [type], optional
            number of cross dihedrals NOT involving H atoms in solute, by default None
        """
        kwargs = {"NPHI": NPHI}
        super().__init__(content=content, FORCEFIELD=FORCEFIELD, MAKETOPVERSION=MAKETOPVERSION, **kwargs)

    def read_content_from_str(self, content: str):
        super().read_content_from_str(content)

    def block_to_string(self) -> str:
        result = self.name + "\n"
        result += (
            "#" + self.field_seperator + "NPHI: number of cross dihedrals NOT involving H atoms" + self.line_seperator
        )
        result += self.field_seperator + str(self.NPHI) + self.line_seperator
        result += super().block_to_string()
        result += "END\n"
        return result


class LJPARAMETERS(_topology_table_block):
    NRATT2: int
    table_header: Iterable[str] = ["IAC", "JAC", "C12", "C6", "CS12", "CS6"]
    table_line_type = ljparameters_type

    def __init__(
        self,
        content: Union[Iterable[ljparameters_type], str],
        FORCEFIELD: FORCEFIELD = None,
        MAKETOPVERSION: MAKETOPVERSION = None,
        NRATT2: int = None,
    ):
        """
                       GROMOS LJPARAMETERS block

        Parameters
        ----------
        content: Union[Iterable[ljparameters_type], str]
        FORCEFIELD : FORCEFIELD
        MAKETOPVERSION : MAKETOPVERSION
        NRATT2 : int, optional
            Number of LJPARAMETERS
        """
        kwargs = {"NRATT2": NRATT2}
        super().__init__(content=content, FORCEFIELD=FORCEFIELD, MAKETOPVERSION=MAKETOPVERSION, **kwargs)

    def read_content_from_str(self, content: str):
        super().read_content_from_str(content)

    def block_to_string(self) -> str:
        result = self.name + "\n"
        result += (
            "#"
            + self.field_seperator
            + "NRATT2: number of LJ interaction types = NRATT*(NRATT+1)/2"
            + self.line_seperator
        )
        result += self.field_seperator + str(self.NRATT2) + self.line_seperator
        result += super().block_to_string()
        result += "END\n"
        return result


class _generic_topology_groups(_topology_block):
    NSM: int
    NSP: List[int]

    def __init__(
        self,
        content: Union[str, dict] = None,
        FORCEFIELD: FORCEFIELD = None,
        MAKETOPVERSION: MAKETOPVERSION = None,
        NSM: int = None,
        NSP: List[int] = None,
    ):
        if NSP is not None:
            if NSM is not None:
                if len(NSP) == NSM:
                    self.NSM = NSM
                else:
                    raise ValueError("FUn")
            else:
                self.NSM = len(NSP)
            self.NSP = list(map(int, NSP))
            super().__init__(
                FORCEFIELD=FORCEFIELD, MAKETOPVERSION=MAKETOPVERSION, content=[str(NSM)] + list(map(str, NSP))
            )
        else:
            super().__init__(FORCEFIELD=FORCEFIELD, MAKETOPVERSION=MAKETOPVERSION, content=content)

            if len(self.content) == 1 and len(self.content[0]) - 1 == int(self.content[0][0]):
                self.NSM = int(self.content[0][0])
                self.NSP = [int(x) for x in self.content[0][1:]]
            elif len(self.content) > 1:
                self.NSM = int(self.content[0][0])
                self.NSP = []
                [
                    self.NSP.extend(list(map(int, t))) if (isinstance(t, list)) else self.NSP.extend([int(t)])
                    for t in self.content[1:]
                ]
            else:
                raise ValueError("SOLUTEMOLECULES has not the correct number of fields.")

        # Clean COntent
        self.content = [[self.NSM]]
        self.content.extend([[x] for x in self.NSP])

    def block_to_string(self) -> str:
        # Clean COntent
        self.content = [[self.NSM]]
        self.content.extend([[x] for x in self.NSP])

        return super().block_to_string()


class SOLUTEMOLECULES(_generic_topology_groups):
    pass


class TEMPERATUREGROUPS(_generic_topology_groups):
    pass


class PRESSUREGROUPS(_generic_topology_groups):
    pass


class LJEXCEPTIONS(_topology_table_block):
    NEX: int
    table_header: Iterable[str] = ["AT1", "AT2", "C12", "C6"]
    table_line_type = ljexception_type

    def __init__(
        self,
        content: Union[Iterable[ljexception_type], str],
        FORCEFIELD: FORCEFIELD = None,
        MAKETOPVERSION: MAKETOPVERSION = None,
        NEX: int = None,
    ):
        """
                       GROMOS LJEXCEPTIONS block

        Parameters
        ----------
        content: Union[Iterable[ljexception_type], str]
        FORCEFIELD : FORCEFIELD
        MAKETOPVERSION : MAKETOPVERSION
        NEX : int, optional
            Number of LJEXCEPTIONS
        """
        kwargs = {"NEX": NEX}
        super().__init__(content=content, FORCEFIELD=FORCEFIELD, MAKETOPVERSION=MAKETOPVERSION, **kwargs)

    def read_content_from_str(self, content: str):
        super().read_content_from_str(content)

    def block_to_string(self) -> str:
        result = self.name + "\n"
        result += "# This block defines special LJ-interactions based on atom numbers \n# This overrules the normal LJ-parameters (including 1-4 interactions)\n"
        result += "#" + self.field_seperator + "NEX: number of exceptions" + self.line_seperator
        result += self.field_seperator + str(self.NEX) + self.line_seperator
        result += super().block_to_string()
        result += "END\n"
        return result


class SOLVENTATOM(_topology_table_block):
    NRAM: int
    table_header: Iterable[str] = ["I", "ANMS", "IACS", "MASS", "CGS"]
    table_line_type = solventatom_type

    def __init__(
        self,
        content: Union[Iterable[solventatom_type], str],
        FORCEFIELD: FORCEFIELD = None,
        MAKETOPVERSION: MAKETOPVERSION = None,
        NRAM: int = None,
    ):
        """
                       GROMOS solventatom block

        Parameters
        ----------
        content: Union[Iterable[solventatom_type], str]
        FORCEFIELD : FORCEFIELD
        MAKETOPVERSION : MAKETOPVERSION
        NRAM : int, optional
            Number of solventatom
        """
        kwargs = {"NRAM": NRAM}
        super().__init__(content=content, FORCEFIELD=FORCEFIELD, MAKETOPVERSION=MAKETOPVERSION, **kwargs)

    def read_content_from_str(self, content: str):
        super().read_content_from_str(content)

    def block_to_string(self) -> str:
        result = self.name + "\n"
        result += "#" + self.field_seperator + "NRAM: number of atoms per solvent molecule" + self.line_seperator
        result += self.field_seperator + str(self.NRAM) + self.line_seperator
        result += super().block_to_string()
        result += "END\n"
        return result


class SOLVENTCONSTR(_topology_table_block):
    NCONS: int
    table_header: Iterable[str] = ["ICONS", "JCONS", "CONS"]
    table_line_type = solventconstr_type

    def __init__(
        self,
        content: Union[Iterable[solventconstr_type], str],
        FORCEFIELD: FORCEFIELD = None,
        MAKETOPVERSION: MAKETOPVERSION = None,
        NCONS: int = None,
    ):
        """
                       GROMOS SOLVENTCONSTR block

        Parameters
        ----------
        content: Union[Iterable[solventconstr_type], str]
        FORCEFIELD : FORCEFIELD
        MAKETOPVERSION : MAKETOPVERSION
        NCONS : int, optional
            Number of SOLVENTCONSTR
        """
        kwargs = {"NCONS": NCONS}
        super().__init__(content=content, FORCEFIELD=FORCEFIELD, MAKETOPVERSION=MAKETOPVERSION, **kwargs)

    def read_content_from_str(self, content: str):
        super().read_content_from_str(content)

    def block_to_string(self) -> str:
        result = self.name + "\n"
        result += "#" + self.field_seperator + "NCONS: number of constraints" + self.line_seperator
        result += self.field_seperator + str(self.NCONS) + self.line_seperator
        result += super().block_to_string()
        result += "END\n"
        return result


class CONSTRAINT(_topology_table_block):
    NCON: int
    table_header: Iterable[str] = ["IC", "JC", "ICC"]
    table_line_type = constraint_type

    def __init__(
        self,
        content: str or dict or None,
        FORCEFIELD: FORCEFIELD = None,
        MAKETOPVERSION: MAKETOPVERSION = None,
        NCON: int = None,
    ):
        kwargs = {"NCON": NCON}
        super().__init__(content=content, FORCEFIELD=FORCEFIELD, MAKETOPVERSION=MAKETOPVERSION, **kwargs)

    def read_content_from_str(self, content: str):
        return super().read_content_from_str(content)

    def block_to_string(self) -> str:
        result = self.name + "\n"
        result += "#" + self.field_seperator + "NCON: number of constraints" + self.line_seperator
        result += self.field_seperator + str(self.NCON) + self.line_seperator
        result += super().block_to_string()
        result += "END\n"
        return result


class atom_lam_pertubation_state(_generic_field):
    state_format_pattern = " {:>5} {:>5} {:>10.5f}"

    def __init__(
        self,
        NR: int,
        RES: int,
        NAME: str,
        STATES: Dict[int, pertubation_lam_state],
        ALPHLJ: float = 1.0,
        ALPHCRF: float = 1.0,
    ):
        self.NR = int(NR)
        self.RES = int(RES)
        self.NAME = NAME
        self.STATES = STATES
        self.ALPHLJ = float(ALPHLJ)
        self.ALPHCRF = float(ALPHCRF)

    def to_string(self):
        state_str = "".join(
            [
                self.state_format_pattern.format(
                    int(self.STATES[x].IAC), float(self.STATES[x].MASS), float(self.STATES[x].CHARGE)
                )
                for x in sorted(self.STATES)
            ]
        )
        format_str = "{:>5} {:>5} {:>5}" + state_str + " {:10.5f} {:10.5f}\n"
        return format_str.format(self.NR, self.RES, self.NAME, self.ALPHLJ, self.ALPHCRF)


class PERTATOMPARAM(_generic_gromos_block):
    def __init__(
        self,
        STATEATOMS: List[atom_lam_pertubation_state] = None,
        STATEATOMHEADER: Tuple[str] = None,
        NJLA: int = None,
        STATEIDENTIFIERS: List[str] = None,
        dummy_IAC: int = 22,
        dummy_CHARGE: int = 0.0,
        content: List[str] = None,
    ):

        self.NPTB = 2
        self.dummy_IAC = dummy_IAC
        self.dummy_CHARGE = dummy_CHARGE

        if content is None:
            if STATEATOMHEADER is None:
                self.STATEATOMHEADER = [
                    "NR",
                    "RES",
                    "NAME",
                ]
                self.STATEATOMHEADER += ["ALPHLJ", "ALPHCRF"]
            else:
                self.STATEATOMHEADER = STATEATOMHEADER

            if STATEATOMS is None:
                self.STATEATOMS = []
            else:
                self.STATEATOMS = []
                self.NJLA = 0

                self.add_state_atoms(STATEATOMS)
            super().__init__(used=True, name=__class__.__name__)
        else:
            super().__init__(used=True, name=__class__.__name__, content=content)

        # You can check yourself :)
        if NJLA is not None and not len(STATEATOMS) == NJLA:
            raise ValueError(
                "NJLA must be equal to the length of STATEATOMS! NJLA="
                + str(NJLA)
                + "\t stateatoms"
                + str(len(STATEATOMS))
                + "\n\n"
                + str(self)
            )

    def read_content_from_str(self, content: List[str]):
        field = 0
        NJLA = None
        STATEIDENTIFIERS = None
        STATEATOMHEADER = None
        STATEATOMS = []
        first = True
        stdid = False
        for line in content:
            if "#" in line:
                # comment = line
                if "state_identifiers" in line:
                    stdid = True
                elif stdid:
                    STATEIDENTIFIERS = line.replace("#", "").split()
                    stdid = False
                continue
            else:
                if field > 0:
                    if first:
                        STATEATOMHEADER = [
                            "NR",
                            "RES",
                            "NAME",
                        ]
                        [
                            STATEATOMHEADER.extend(["IAC" + str(x), "MASS" + str(x), "CHARGE" + str(x)])
                            for x in range(1, 3)
                        ]
                        STATEATOMHEADER += ["ALPHLJ", "ALPHCRF"]
                        first = False

                    state_line = {key: value for key, value in zip(STATEATOMHEADER, line.split())}
                    final_state_line = {
                        key: state_line[key]
                        for key in state_line
                        if ("IAC" not in key and "CHARGE" not in key and "MASS" not in key)
                    }
                    states = {
                        x: pertubation_lam_state(
                            IAC=int(round(float(state_line["IAC" + str(x)]))),
                            MASS=float(state_line["MASS" + str(x)]),
                            CHARGE=float(state_line["CHARGE" + str(x)]),
                        )
                        for x in range(1, 3)
                    }

                    final_state_line.update({"STATES": states})
                    STATEATOMS.append(atom_lam_pertubation_state(**final_state_line))

                elif field == 0:
                    NJLA = int(line.strip())
                field += 1

        self.NJLA = NJLA
        self.STATEIDENTIFIERS = STATEIDENTIFIERS
        self.STATEATOMHEADER = STATEATOMHEADER
        self.STATEATOMS = STATEATOMS

    @property
    def nStates(self) -> int:
        return self.NPTB

    @property
    def nTotalStateAtoms(self) -> int:
        return self.NJLA

    @property
    def states(self) -> dict:
        return {
            self.STATEIDENTIFIERS[state - 1]: {
                atom.NR: atom.STATES[state] for atom in sorted(self.STATEATOMS, key=lambda x: x.NR)
            }
            for state in range(1, self.NPTB + 1)
        }

    """
    ADD FUNCTIONS
    """

    def add_state_atoms(self, state_atoms: List[atom_lam_pertubation_state]):
        """
        This function can add states and atoms, but also overwrite state values of existing atoms.
        If a new state is defined only for a subset of atoms, all other atoms are set to the default dummy.
        If a new atom misses a state definition, this state will be set to dummy.
        Parameters
        ----------
        state_atoms: List[atom_eds_pertubation_state]
        """

        # some preperations:
        pre_dummy_state = lambda atomMass: pertubation_lam_state(  # noqa: E731
            IAC=self.dummy_IAC, MASS=atomMass, CHARGE=self.dummy_CHARGE
        )
        insert_id = self.STATEATOMHEADER.index("ALPHLJ")

        # find all new states
        keys = np.array([list(natom.STATES.keys()) for natom in state_atoms], ndmin=1)
        unique_stateIDs = np.unique(np.concatenate(keys))
        # Todo: not urgent; state number adaptation ( present states 1,2,3,4 new state 8 - id should be 5 not 8)
        unique_states = list(map(str, ["state" + str(x) if isinstance(x, Number) else x for x in unique_stateIDs]))

        # insert new state IDs
        off = 0
        for unique_state in unique_stateIDs:
            self.STATEATOMHEADER.insert(insert_id + off, "IAC" + str(unique_state))
            self.STATEATOMHEADER.insert(insert_id + off + 1, "mass" + str(unique_state))
            self.STATEATOMHEADER.insert(insert_id + off + 2, "CHARGE" + str(unique_state))
            off += 3

        # add new state names
        if hasattr(self, "STATEIDENTIFIERS"):
            self.STATEIDENTIFIERS.extend(unique_states)
            self.NPTB += len(unique_states)
        else:
            self.STATEIDENTIFIERS = unique_states
            self.NPTB = len(unique_states)
        # increase the number of new states

        # 1. Update already present atoms:
        atomIDs = [atom.NR for atom in state_atoms]
        for atom in self.STATEATOMS:
            atom.STATES.update({key: val for key, val in atom.STATES.items()})
            possible_masses = [val.MASS for key, val in atom.STATES.items() if (val.MASS > 0)]
            dummy_state = pre_dummy_state(atomMass=possible_masses[0])

            if atom.NR in atomIDs:
                new_atom = state_atoms[atomIDs.index(atom.NR)]

                atom.NAME = new_atom.NAME
                atom.STATES.update({key: val for key, val in new_atom.STATES.items()})
                possible_masses = [val.MASS for key, val in new_atom.STATES.items() if (val.MASS > 0)]
                # add missing dummies
                # print(unique_stateIDs)
                atom.STATES.update({key: dummy_state for key in unique_stateIDs if key not in atom.STATES})

                # remove present atom
                del atomIDs[atomIDs.index(atom.NR)]

            else:
                # add missing dummies
                atom.STATES.update({key: dummy_state for key in unique_stateIDs if key not in atom.STATES})

        # 2. add new atoms
        new_atoms = [atom for atom in state_atoms if (atom.NR in atomIDs)]
        for atom in new_atoms:
            atom.STATES.update({key: val for key, val in atom.STATES.items()})
            possible_masses = [val.MASS for key, val in atom.STATES.items() if (val.MASS > 0)]
            dummy_state = pre_dummy_state(atomMass=possible_masses[0])

            atom.STATES.update({key: dummy_state for key in range(1, self.NPTB + 1) if (key not in atom.STATES)})
            self.STATEATOMS.append(atom)
            self.NJLA += 1

    """
    DELETING FUNCTIONS
    """

    def delete_state(self, stateIDs: Union[int, List[int]] = None, stateNames: Union[str, List[str]] = None):
        """
        This function deletes an state column.
        Parameters
        ----------
        stateIDs: int
            number of the state
        Returns
        -------
        """
        if stateIDs is not None:
            if isinstance(stateIDs, int):
                stateIDs = [stateIDs]

            for state in stateIDs:
                for atom in self.STATEATOMS:
                    if state in atom.STATES:
                        del atom.STATES[state]
                del self.STATEIDENTIFIERS[state - 1]
                self.STATEATOMHEADER = [
                    x for x in self.STATEATOMHEADER if (not x == "IAC" + str(state) and not "CHARGE" + str(state) == x)
                ]

            self.NPTB -= len(set(stateIDs))

        elif stateNames is not None:
            if isinstance(stateNames, str):
                stateNames = [stateNames]

            for stateN in stateNames:
                # print(stateN)
                stateID = self.STATEIDENTIFIERS.index(stateN) + 1

                for atom in self.STATEATOMS:
                    if stateID in atom.STATES:
                        del atom.STATES[stateID]

                del self.STATEIDENTIFIERS[stateID - 1]
                self.STATEATOMHEADER = [
                    x
                    for x in self.STATEATOMHEADER
                    if (not x == "IAC" + str(stateID) and not "CHARGE" + str(stateID) == x)
                ]
            self.NPTB -= len(set(stateNames))

        else:
            raise Exception("Please give either stateNames or stateIDs")

    def delete_atom(self, atomNR: Union[int, List[int]]):
        """
        This function removes atom lines from the ptp file.
        Parameters
        ----------
        atomNR: int
            atom to be removed.
        """
        if isinstance(atomNR, int):
            atomNR = [atomNR]

        # ind_offset = 0
        new_STATEATOMS = []
        for ind, atom in enumerate(self.STATEATOMS):
            if atom.NR in atomNR:
                continue
            else:
                new_STATEATOMS.append(atom)

        self.STATEATOMS = new_STATEATOMS
        self.NJLA -= len(atomNR)

    """
    STR FUNCTIONS
    """

    def _state_STATEATOMHEADER_str(self):
        state_format_pattern = (
            "{:>5} {:>5} {:>5}" + "".join([" {:>5}{:>5}{:>10}" for x in range(self.NPTB)]) + "    {:10} {:10}"
        )
        return state_format_pattern.format(*self.STATEATOMHEADER)

    def block_to_string(self) -> str:
        result = self.name + self.line_seperator
        result += (
            "# NJLA "
            + self.field_seperator
            + "NPTB = "
            + self.field_seperator
            + str(self.NPTB)
            + self.field_seperator
            + self.line_seperator
        )
        result += self.field_seperator + str(self.NJLA) + self.line_seperator
        result += "# state_identifiers" + self.line_seperator
        result += (
            "# "
            + self.field_seperator
            + self.field_seperator.join(map(str, self.STATEIDENTIFIERS))
            + self.line_seperator
        )
        result += "# " + self._state_STATEATOMHEADER_str() + self.line_seperator
        result += "".join(map(str, sorted(self.STATEATOMS, key=lambda x: x.NR)))
        result += "END" + self.line_seperator
        return result


class SCALEDINTERACTIONS(_generic_gromos_block):
    def __init__(self, values=None, content=None):
        """
        Not exactly sure what these parameters do
        """

        if content is None:
            super().__init__(used=True, name=__class__.__name__)
            self.values = values
        else:
            super().__init__(used=True, name="SCALEDINTERACTIONS", content=content)

    def block_to_string(self) -> str:

        result = self.name + self.line_seperator
        for i, v in enumerate(self.values):
            result += str(v) + self.field_seperator
            if not i:
                result += self.line_seperator
        result += self.line_seperator + "END\n"
        return result

    def read_content_from_str(self, content):
        # values
        values = []
        values.append(content[0])
        for v in content[1].split():
            values.append(v)

        self.values = values
