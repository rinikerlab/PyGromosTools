import re
from enum import Enum

from pygromos.files.blocks._general_blocks import _generic_gromos_block, _generic_field, _iterable_gromos_block
from pygromos.files.blocks._general_blocks import TITLE as generic_TITLE
from pygromos.files.blocks._general_blocks import TIMESTEP as generic_TIMESTEP
from pygromos.files.blocks._general_blocks import TRAJ as generic_TRAJ
from pygromos.files import blocks
from pygromos.utils.typing import List, Number

# forward declarations
TITLE: generic_TITLE = generic_TITLE
TIMESTEP: generic_TIMESTEP = generic_TIMESTEP
TRAJ: generic_TRAJ = generic_TRAJ


##########################################################################
#   ENUMS
##########################################################################


class Pbc(Enum):
    trunc_octahedron = -1
    vacuum = 0
    rectangular = 1
    triclinic = 2


##########################################################################
#   FIELDS
##########################################################################


class atomP(_generic_field):
    def __init__(self, resID: int, resName: str, atomType: str, atomID: int, xp: float, yp: float, zp: float):
        """

        Parameters
        ----------
        resID
        resName
        atomType
        atomID
        xp
        yp
        zp
        """
        self.resID = resID
        self.resName = resName
        self.atomType = atomType
        self.atomID = atomID
        self.xp = xp
        self.yp = yp
        self.zp = zp

    def to_string(self) -> str:
        return "{: >5} {: <5} {: <6} {: >5}{: 15.9f}{:15.9f}{:15.9f}\n".format(
            self.resID, self.resName, self.atomType, self.atomID, self.xp, self.yp, self.zp
        )


class atomV(_generic_field):
    def __init__(self, resID: int, resName: str, atomType: str, atomID: int, xv: float, yv: float, zv: float):
        self.resID = resID
        self.resName = resName
        self.atomType = atomType
        self.atomID = atomID
        self.xv = xv
        self.yv = yv
        self.zv = zv

    def to_string(self) -> str:
        return "{: >5} {: >3}  {: <5}{: >7}{: 15.9f}{:15.9f}{:15.9f}\n".format(
            self.resID, self.resName, self.atomType, self.atomID, self.xv, self.yv, self.zv
        )


class lattice_shift(_generic_field):
    def __init__(self, atomID: int, x: float, y: float, z: float):
        self.atomID = atomID
        self.x = x
        self.y = y
        self.z = z

    def to_string(self) -> str:
        return "{:>10}{:>10}{:>10}\n".format(self.x, self.y, self.z)


class atomSI(_generic_field):
    def __init__(self, resID: int, resName: str, atomType: str, atomID: int, sxx: float, sxy: float, sxz: float):
        self.resID = resID
        self.resName = resName
        self.atomType = atomType
        self.atomID = atomID
        self.sxx = sxx
        self.sxy = sxy
        self.sxz = sxz

    def to_string(self) -> str:
        return "{: >5} {: >3}  {: <5}{: >7}{: 15.9f}{:15.9f}{:15.9f}\n".format(
            self.resID, self.resName, self.atomType, self.atomID, self.sxx, self.sxy, self.sxz
        )


##########################################################################
#   BLOCKS
##########################################################################


class POSITION(_iterable_gromos_block):
    """
    POSITION

    Parameters
    ----------
    content: List[atomP]
        every element in this list is of atom position obj

    """

    def __init__(self, content: List[atomP]):
        super().__init__(used=True, name="POSITION", content=content)

    def _check_import_method(self, content: str):
        if isinstance(content, list) and all([isinstance(x, atomP) for x in content]):
            self.content = content
        elif isinstance(content, str):
            self.read_content_from_str(content=content.split(self.line_seperator))
        elif isinstance(content, list) and all([isinstance(x, str) for x in content]):
            self.read_content_from_str(content=content)
        elif content is None:
            self.content = []
        else:
            raise Exception("Generic Block did not understand the type of content \n content: \n" + str(content))

    def read_content_from_str(self, content: List[str]):
        self.content = [
            blocks.coords.atomP(
                resID=int(x[0]),
                resName=str(x[1]),
                atomType=str(x[2]),
                atomID=int(x[3]),
                xp=float(x[4]),
                yp=float(x[5]),
                zp=float(x[6]),
            )
            for x in list(map(lambda x: x.split(), content))
            if (not x[0] == "#")
        ]


class POSRESSPEC(_iterable_gromos_block):
    """
    POSITION

    Parameters
    ----------
    content: List[atomP]
        every element in this list is of atom position obj

    """

    def __init__(self, content: List[atomP]):
        super().__init__(used=True, name="POSRESSPEC", content=content)

    def read_content_from_str(self, content: List[str]):
        self.content = [
            blocks.coords.atomP(
                resID=int(x[0]),
                resName=str(x[1]),
                atomType=str(x[2]),
                atomID=int(x[3]),
                xp=float(x[4]),
                yp=float(x[5]),
                zp=float(x[6]),
            )
            for x in list(map(lambda x: x.split(), content))
            if (not x[0] == "#")
        ]


class VELOCITY(_iterable_gromos_block):
    """
    VELOCITY

    Parameters
    ----------
    content: List[atomV]
        every element in this list is of atom velocity obj

    """

    def __init__(self, content: List[atomV]):
        super().__init__(used=True, name="VELOCITY", content=content)

    def read_content_from_str(self, content: List[str]):
        self.content = [
            atomV(
                resID=int(x[0]),
                resName=str(x[1]),
                atomType=str(x[2]),
                atomID=int(x[3]),
                xv=float(x[4]),
                yv=float(x[5]),
                zv=float(x[6]),
            )
            for x in list(map(lambda x: x.split(), content))
            if (not x[0] == "#")
        ]


class STOCHINT(_iterable_gromos_block):
    """
    STOCHINT

    Parameters
    ----------
    content: List[atomSI]
        every element in this list is of atom stochastic interval obj

    seed: str
        contains the seed for the stochastic dynamics simulation
    """

    def __init__(self, content: List[atomSI]):
        super().__init__(used=True, name="STOCHINT", content=content)

    def block_to_string(self) -> str:
        result = self.name + self.line_seperator
        result += "#" + self.field_seperator + self.field_seperator.join(self.table_header) + "\n"
        for x in self.content:
            result += x.to_string()
        result += "# seed\n"
        result += self.seed
        if result[-1] != "\n":
            result += "\n"
        result += "END\n"
        return result

    def read_content_from_str(self, content: List[str]):
        if content is str:
            content = content.split("\n")
        self.content = [
            blocks.coords.atomSI(
                resID=int(x[0]),
                resName=str(x[1]),
                atomType=str(x[2]),
                atomID=int(x[3]),
                sxx=float(x[4]),
                sxy=float(x[5]),
                sxz=float(x[6]),
            )
            for x in list(map(lambda x: x.split(), content[:-1]))
            if (not x[0] == "#")
        ]

        # seed safety check
        x = content[-1].split()
        if len(x) == 7:
            raise ValueError(
                "The seed of STOCHINT has a length of 7 pleas use longer once of check if seed is present! \nGOT: "
                + str(x)
            )
        else:
            self.seed = content[-1]


class REFPOSITION(_iterable_gromos_block):
    """
    REFPOSITION

    Parameters
    ----------
    content: List[atomP]
        every element in this list is of atom position obj

    """

    def __init__(self, content: List[atomP]):
        """

        Parameters
        ----------
        content: List[atomP]
            every element in this list is of atom position obj

        """
        super().__init__(used=True, name="REFPOSITION", content=content)

    def _check_import_method(self, content: str):
        if isinstance(content, list) and all([isinstance(x, atomP) for x in content]):
            self.content = content
        elif isinstance(content, str):
            self.read_content_from_str(content=content.split(self.line_seperator))
        elif isinstance(content, list) and all([isinstance(x, str) for x in content]):
            self.read_content_from_str(content=content)
        elif content is None:
            self.content = []
        else:
            raise Exception("Generic Block did not understand the type of content \n content: \n" + str(content))

    def read_content_from_str(self, content: List[str]):
        self.content = [
            blocks.coords.atomP(
                resID=int(x[0]),
                resName=str(x[1]),
                atomType=str(x[2]),
                atomID=int(x[3]),
                xp=float(x[4]),
                yp=float(x[5]),
                zp=float(x[6]),
            )
            for x in list(map(lambda x: x.split(), content))
            if (not x[0] == "#")
        ]


class LATTICESHIFTS(_iterable_gromos_block):
    """
    LATTICESHIFTS

    Parameters
    ----------
    content: List[lattice_shift]
        every element in this list is a lattice shift obj

    """

    def __init__(self, content: List[lattice_shift]):
        """

        Parameters
        ----------
        content: List[lattice_shift]
            every element in this list is a lattice shift obj
        """
        super().__init__(used=True, name="LATTICESHIFTS", content=content)

    def read_content_from_str(self, content: List[str]):
        subblock1 = list(map(lambda x: re.findall(r"[\w]+", x.strip()), content))
        subblock2 = list(map(lambda x: [x.strip() for x in re.findall(r"[\W]+", x.strip())], content))

        subblock = []
        for number, sign in zip(subblock1, subblock2):
            if "#" in sign:
                continue
            elif len(sign) == 2:
                row = [number[0], sign[0] + number[1], sign[1] + number[2]]
            elif len(sign) == 3:
                row = [sign[0] + number[0], sign[1] + number[1], sign[2] + number[2]]
            else:
                raise Exception("This does not work! \n SIGN: " + str(sign) + "\n Number: " + str(number))
            subblock.append(row)

        if all(len(x) == 3 for x in subblock):
            self.content = list(
                map(
                    lambda c: blocks.coords.lattice_shift(
                        atomID=int(c[0]), x=int(c[1][0]), y=int(c[1][1]), z=int(c[1][2])
                    ),
                    enumerate(subblock),
                )
            )
        else:
            short_lines = [str(x) for x in subblock if (len(x) != 3)]
            raise IOError(
                "inconsistent Atom LatticeShifts line lenghts (have to be =3 fields!). Problem in line: "
                + "\n\t".join(short_lines)
            )


class GENBOX(_generic_gromos_block):
    """GENBOX

        This Block is representing the simulation Box in a coordinate file.

    Attributes
    ----------
    pbc: int,Pbc
        Periodic Boundary Condition
    length: List[float]
    angles: List[float]
    euler: List[float]
    origin: List[float]

    """

    def __init__(
        self,
        pbc: Pbc = Pbc(0),
        length: List[float] = [0.0, 0.0, 0.0],
        angles: List[float] = [0.0, 0.0, 0.0],
        euler: List[float] = [0.0, 0.0, 0.0],
        origin: List[float] = [0.0, 0.0, 0.0],
        content=None,
    ):
        """

        Parameters
        ----------
        pbc: int,Pbc
        length: List[float]
        angles: List[float]
        euler: List[float]
        origin: List[float]

        """

        if content is not None:
            super().__init__(used=True, name="GENBOX", content=content)
        else:
            super().__init__(used=True, name="GENBOX")
            self._pbc = Pbc(pbc)
            self._length = length
            self._angles = angles
            self._euler = euler
            self._origin = origin

    def block_to_string(self) -> str:
        result = self.name + "\n"
        result += "{:>5}\n".format(str(self.pbc.value))
        result += "{:>15.9f}{:>15.9f}{:>15.9f}\n".format(self.length[0], self.length[1], self.length[2])
        result += "{:>15.9f}{:>15.9f}{:>15.9f}\n".format(self.angles[0], self.angles[1], self.angles[2])
        result += "{:>15.9f}{:>15.9f}{:>15.9f}\n".format(self.euler[0], self.euler[1], self.euler[2])
        result += "{:>15.9f}{:>15.9f}{:>15.9f}\n".format(self.origin[0], self.origin[1], self.origin[2])
        result += "END\n"

        return result

    def read_content_from_str(self, content: List[str]):
        if len(content[0].split()) == 1:
            self._pbc = blocks.coords.Pbc(int(content[0].strip()))
        else:
            raise IOError("Could not read pbc information!")
        if len(content[1].split()) == 3:
            self._length = list(map(float, content[1].strip().split()))
        else:
            raise IOError("Could not read box length information!")
        if len(content[2].split()) == 3:
            self._angles = list(map(float, content[2].strip().split()))
        else:
            raise IOError("Could not read box angles information!")
        if len(content[3].split()) == 3:
            self._euler = list(map(float, content[3].strip().split()))
        else:
            raise IOError("Could not read box euler information!")
        if len(content[4].split()) == 3:
            self._origin = list(map(float, content[4].strip().split()))
        else:
            raise IOError("Could not read box origin information!")

    # Attributes
    @property
    def pbc(self) -> Pbc:
        return self._pbc

    @pbc.setter
    def pbc(self, pbc: Pbc):
        if isinstance(pbc, Pbc):
            self._pbc = pbc
        elif isinstance(pbc, int) or (isinstance(pbc, str) and pbc.isalnum()):
            if int(pbc) in Pbc._value2member_map_:
                self._pbc = Pbc(int(pbc))
            else:
                raise ValueError("unknown int for pbc\n Use: \n" + str(Pbc.__members__))
        else:
            raise ValueError("Periodic boundary Condition must be int or PBC-Enum")

    @property
    def length(self) -> List[float]:
        return self._length

    @length.setter
    def length(self, length: List[float]):
        if isinstance(length, List) and all([isinstance(x, Number) for x in length]):
            self._length = length
        else:
            raise ValueError("length must be List[float]")

    @property
    def angles(self) -> List[float]:
        return self._angles

    @angles.setter
    def angles(self, angles: List[float]):
        if isinstance(angles, List) and all([isinstance(x, Number) for x in angles]):
            self._angles = angles
        else:
            raise ValueError("angles must be List[float]")

    @property
    def euler(self) -> List[float]:
        return self._euler

    @euler.setter
    def euler(self, euler: List[float]):
        if isinstance(euler, List) and all([isinstance(x, Number) for x in euler]):
            self._euler = euler
        else:
            raise ValueError("euler must be List[float]")

    @property
    def origin(self) -> List[float]:
        return self._origin

    @origin.setter
    def origin(self, origin: List[float]):
        if isinstance(origin, List) and all([isinstance(x, Number) for x in origin]):
            self._angles = origin
        else:
            raise ValueError("origin must be List[float]")


class PERTDATA(_generic_gromos_block):

    content: float

    def __init__(self, content: List[str]):
        """
           This block is used for lambda-sampling and gives the lambda value of the current coordinates.

        Parameters
        ----------
        lambda_value: float
            current lambda value
        """
        super(PERTDATA, self).__init__(name=__class__.__name__, used=True, content=content)

    def read_content_from_str(self, content: List[str]):
        self.content = float("\n".join(content).strip())

    @property
    def lam(self) -> float:
        return self.content

    @lam.setter
    def lam(self, lam: float):
        self.content = float(lam)
