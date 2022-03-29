import __main__
import numpy as np
from collections import namedtuple

from pygromos.files.blocks._general_blocks import _generic_gromos_block, _generic_field
from pygromos.files.blocks._general_blocks import TITLE as generic_TITLE
from pygromos.utils.typing import Union, List, Tuple, Dict, Number

TITLE: generic_TITLE = generic_TITLE

"""
   FIELDS
"""

pertubation_eds_state = namedtuple("pertubationEdsState", ["IAC", "CHARGE"])
pertubation_lam_state_nonbonded = namedtuple("pertubationLamState", ["IAC", "MASS", "CHARGE"])

setattr(__main__, pertubation_eds_state.__name__, pertubation_eds_state)
pertubation_eds_state.__module__ = "__main__"

setattr(__main__, pertubation_lam_state_nonbonded.__name__, pertubation_lam_state_nonbonded)
pertubation_lam_state_nonbonded.__module__ = "__main__"


class atom_mass_type(_generic_field):
    def __init__(self, N: int, ATMAS: float, ATMASN: str, comment: str = ""):
        self.N = N
        self.ATMAS = ATMAS
        self.ATMASN = ATMASN
        self.comment = comment

    def to_string(self) -> str:
        return self.comment + "{:<8} {:<3.4f} {:<5}\n".format(self.N, self.ATMAS, self.ATMASN)


class atom_eds_pertubation_state(_generic_field):
    state_format_pattern = " {:>3} {:>10.5f}"

    def __init__(
        self, NR: int, NAME: str, STATES: Dict[int, pertubation_eds_state], ALPHLJ: float = 1.0, ALPHCRF: float = 1.0
    ):
        self.NR = int(NR)
        self.NAME = NAME
        self.STATES = STATES
        self.ALPHLJ = float(ALPHLJ)
        self.ALPHCRF = float(ALPHCRF)

    def to_string(self) -> str:
        state_str = "".join(
            [
                self.state_format_pattern.format(int(self.STATES[x].IAC), float(self.STATES[x].CHARGE))
                for x in sorted(self.STATES)
            ]
        )
        format_str = "{:>5} {:>5}" + state_str + " {:10.5f} {:10.5f}\n"
        return format_str.format(self.NR, self.NAME, self.ALPHLJ, self.ALPHCRF)


class atom_lam_pertubation_state(_generic_field):
    state_format_pattern = " {:>5} {:>5} {:>10.5f}"

    def __init__(
        self,
        NR: int,
        RES: int,
        NAME: str,
        STATES: Dict[int, pertubation_lam_state_nonbonded],
        ALPHLJ: float = 1.0,
        ALPHCRF: float = 1.0,
    ):
        self.NR = int(NR)
        self.RES = int(RES)
        self.NAME = NAME
        self.STATES = STATES
        self.ALPHLJ = float(ALPHLJ)
        self.ALPHCRF = float(ALPHCRF)

    def to_string(self) -> str:
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


class atom_lam_pertubation_state_bond(_generic_field):
    state_format_pattern = " {:>5}"

    def __init__(self, NR: int, atomI: int, atomJ: int, STATES: Dict[int, int]):
        self.NR = int(NR)
        self.atomI = atomI
        self.atomJ = atomJ
        self.STATES = STATES

    def to_string(self) -> str:
        state_str = "".join([self.state_format_pattern.format(int(self.STATES[x])) for x in sorted(self.STATES)])
        format_str = "{:>5} {:>5}" + state_str + "\n"
        return format_str.format(self.atomI, self.atomJ)


class atom_lam_pertubation_state_angle(_generic_field):
    state_format_pattern = " {:>5}"

    def __init__(self, NR: int, atomI: int, atomJ: int, atomK: int, STATES: Dict[int, int]):
        self.NR = int(NR)
        self.atomI = atomI
        self.atomJ = atomJ
        self.atomK = atomK
        self.STATES = STATES

    def to_string(self) -> str:
        state_str = "".join([self.state_format_pattern.format(int(self.STATES[x])) for x in sorted(self.STATES)])
        format_str = "{:>5} {:>5} {:>5}" + state_str + "\n"
        return format_str.format(self.atomI, self.atomJ, self.atomK)


class atom_lam_pertubation_state_dihedral(_generic_field):
    state_format_pattern = " {:>5}"

    def __init__(self, NR: int, atomI: int, atomJ: int, atomK: int, atomL: int, STATES: Dict[int, int]):
        self.NR = int(NR)
        self.atomI = atomI
        self.atomJ = atomJ
        self.atomK = atomK
        self.atomL = atomL
        self.STATES = STATES

    def to_string(self) -> str:
        state_str = "".join([self.state_format_pattern.format(int(self.STATES[x])) for x in sorted(self.STATES)])
        format_str = "{:>5} {:>5} {:>5} {:>5}" + state_str + "\n"
        return format_str.format(self.atomI, self.atomJ, self.atomK, self.atomL)


"""
    BLOCKS
"""
# NONBONDED


class MPERTATOM(_generic_gromos_block):
    def __init__(
        self,
        NJLA: int = None,
        NPTB: int = None,
        STATEIDENTIFIERS: List[str] = [],
        STATEATOMHEADER: Tuple[str] = ["NR", "NAME", "ALPHLJ", "ALPHCRF"],
        STATEATOMS: List[atom_eds_pertubation_state] = [],
        dummy_IAC: int = 22,
        dummy_CHARGE: float = 0.0,
        content: List[str] = None,
    ):
        """
            This block is used for lambda sampling to define the different states.

        Parameters
        ----------
        NJLA : int
            number of perturbed atoms
        NPTB : int
            number of pertubation states
        STATEIDENTIFIERS  : List[str]
            string names for states
        STATEATOMHEADER
            header for the atom description table
        STATEATOMS
            list of atoms, that shall be perturbed
        dummy_IAC
            dummy atom VdW type for perturbed atoms
        dummy_CHARGE
            dummy atom charge type for perturbed atoms
        """

        if content is None:
            super().__init__(used=True, name=__class__.__name__)
            self.NJLA = NJLA
            self.NPTB = NPTB
            self.STATEIDENTIFIERS = STATEIDENTIFIERS
            self.STATEATOMHEADER = STATEATOMHEADER
            self.STATEATOMS = STATEATOMS

        else:
            super().__init__(used=True, name=__class__.__name__, content=content)

        self.dummy_IAC = dummy_IAC
        self.dummy_CHARGE = dummy_CHARGE

    def read_content_from_str(self, content: List[str]):
        field = 0
        NJLA = None
        STATEIDENTIFIERS = []
        STATEATOMS = []

        first = True
        for line in content:
            # print(line)
            if "#" in line:
                pass
            else:
                if field > 3:
                    if first:
                        STATEATOMHEADER = [
                            "NR",
                            "NAME",
                        ]
                        [STATEATOMHEADER.extend(["IAC" + str(x), "CHARGE" + str(x)]) for x in range(1, self.NPTB + 1)]
                        STATEATOMHEADER += ["ALPHLJ", "ALPHCRF"]
                        self.STATEATOMHEADER = STATEATOMHEADER
                        first = False

                    state_line = {key: value for key, value in zip(self.STATEATOMHEADER, line.split())}
                    final_state_line = {
                        key: state_line[key] for key in state_line if ("IAC" not in key and "CHARGE" not in key)
                    }

                    states = {
                        x: pertubation_eds_state(
                            IAC=int(state_line["IAC" + str(x)]), CHARGE=float(state_line["CHARGE" + str(x)])
                        )
                        for x in range(1, 1 + self.NPTB)
                    }

                    final_state_line.update({"STATES": states})
                    STATEATOMS.append(atom_eds_pertubation_state(**final_state_line))

                elif field == 0:
                    NJLA, NPTB = tuple(map(int, line.split()))
                    self.NJLA = NJLA
                    self.NPTB = NPTB
                elif field == 1:
                    STATEIDENTIFIERS = line.split()
                    self.STATEIDENTIFIERS = STATEIDENTIFIERS
                field += 1

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

    def add_state_atoms(self, state_atoms: List[atom_eds_pertubation_state]):
        """
        This function can add states and atoms, but also overwrite state values of existing atoms.
        If a new state is defined only for a subset of atoms, all other atoms are set to the default dummy.
        If a new atom misses a state definition, this state will be set to dummy.


        Parameters
        ----------
        state_atoms: List[atom_eds_pertubation_state]

        """

        # some preperations:
        dummy_state = pertubation_eds_state(IAC=self.dummy_IAC, CHARGE=self.dummy_CHARGE)
        insert_id = self.STATEATOMHEADER.index("ALPHLJ")

        # find all new states
        unique_stateIDs = np.unique(np.concatenate([list(natom.STATES.keys()) for natom in state_atoms]))
        # TODO: not urgent; state number adaptation ( present states 1,2,3,4 new state 8 - id should be 5 not 8)
        unique_states = list(map(str, ["state" + str(x) if isinstance(x, Number) else x for x in unique_stateIDs]))

        # insert new state IDs
        off = 0
        for unique_state in unique_stateIDs:
            self.STATEATOMHEADER.insert(insert_id + off, "IAC" + str(unique_state))
            self.STATEATOMHEADER.insert(insert_id + off + 1, "CHARGE" + str(unique_state))
            off += 2

        # add new state names
        self.STATEIDENTIFIERS.extend(unique_states)

        # increase the number of new states
        self.NPTB += len(unique_states)

        # 1. Update already present atoms:
        atomIDs = [atom.NR for atom in state_atoms]
        for atom in self.STATEATOMS:
            if atom.NR in atomIDs:
                new_atom = state_atoms[atomIDs.index(atom.NR)]

                atom.NAME = new_atom.NAME
                atom.STATES.update({key: val for key, val in new_atom.STATES.items()})

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

        elif stateNames is None and stateIDs is None:
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
        state_format_pattern = "{:>3} {:>5}" + "".join([" {:>3}{:>10}" for x in range(self.NPTB)]) + "    {:10} {:10}"

        if len(self.STATEATOMHEADER) != self.NPTB * 2 + 4:
            tmp_list = " ".join(["CHARGE" + str(x) + " " + "IAC" + str(x) for x in range(self.NPTB)])
            self.STATEATOMHEADER = self.STATEATOMHEADER[:2] + tmp_list.split(" ") + self.STATEATOMHEADER[-2:]

        return state_format_pattern.format(*self.STATEATOMHEADER)

    def block_to_string(self) -> str:
        result = self.name + self.line_seperator
        result += "# NJLA " + self.field_seperator + "NPTB" + self.line_seperator
        result += self.field_seperator + str(self.NJLA) + self.field_seperator + str(self.NPTB) + self.line_seperator
        result += "# state_identifiers" + self.line_seperator
        result += (
            self.field_seperator + self.field_seperator.join(map(str, self.STATEIDENTIFIERS)) + self.line_seperator
        )
        result += "# " + self._state_STATEATOMHEADER_str() + self.line_seperator
        result += "".join(map(str, sorted(self.STATEATOMS, key=lambda x: x.NR)))
        result += "END" + self.line_seperator
        return result


class PERTATOMPARAM(_generic_gromos_block):
    def __init__(
        self,
        STATEATOMS: List[atom_lam_pertubation_state] = None,
        STATEATOMHEADER: Tuple[str] = None,
        NJLA: int = None,
        STATEIDENTIFIERS=None,
        dummy_IAC=22,
        dummy_CHARGE=0.0,
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
                for s in range(self.NPTB):
                    self.STATEATOMHEADER += [
                        "IAC",
                        "MASS",
                        "CHARGE",
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
        if NJLA is not None and len(STATEATOMS) != NJLA:
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
                        x: pertubation_lam_state_nonbonded(
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
    def states(self) -> Dict:
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
        pre_dummy_state = lambda atomMass: pertubation_lam_state_nonbonded(  # noqa: E731
            IAC=self.dummy_IAC, MASS=atomMass, CHARGE=self.dummy_CHARGE
        )
        insert_id = self.STATEATOMHEADER.index("ALPHLJ")

        # find all new states
        keys = np.array([list(natom.STATES.keys()) for natom in state_atoms], ndmin=1)
        unique_stateIDs = np.unique(np.concatenate(keys))
        # TODO: not urgent; state number adaptation ( present states 1,2,3,4 new state 8 - id should be 5 not 8)
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

        elif stateNames is None and stateIDs is None:
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


# BONDED


class PERTBONDSTRETCH(_generic_gromos_block):
    def __init__(
        self,
        STATEATOMS: List[atom_lam_pertubation_state_bond] = None,
        STATEATOMHEADER: Tuple[str] = None,
        NPB: int = None,
        dummy_BOND=22,
        content: List[str] = None,
    ):
        self.NPTB = 2
        self.dummy_BOND = dummy_BOND

        if content is None:
            if STATEATOMHEADER is None:
                self.STATEATOMHEADER = ["atomI", "atomJ", "type1", "type2"]
            else:
                self.STATEATOMHEADER = STATEATOMHEADER

            if STATEATOMS is None:
                self.STATEATOMS = []
            else:
                self.STATEATOMS = []
                self.NPB = 0

                # self.add_state_atoms(STATEATOMS)
            super().__init__(used=True, name=__class__.__name__)
        else:
            super().__init__(used=True, name=__class__.__name__, content=content)

        # You can check yourself :)
        if NPB is not None and len(STATEATOMS) != NPB:
            raise ValueError(
                "NJLA must be equal to the length of STATEATOMS! NJLA="
                + str(NPB)
                + "\t stateatoms"
                + str(len(STATEATOMS))
                + "\n\n"
                + str(self)
            )

    @property
    def nStates(self) -> int:
        return self.NPTB

    @property
    def nTotalStateAtoms(self) -> int:
        return self.NPB

    @property
    def states(self) -> dict:
        return {
            self.STATEIDENTIFIERS[state - 1]: {
                atom.NR: atom.STATES[state] for atom in sorted(self.STATEATOMS, key=lambda x: x.NR)
            }
            for state in range(1, self.NPTB + 1)
        }

    def read_content_from_str(self, content: List[str]):
        field = 0
        NPB = None
        STATEIDENTIFIERS = None
        STATEATOMHEADER = None
        STATEATOMS = []
        first = True
        stdid = False
        for line in content:
            if "#" in line:
                if "state_identifiers" in line:
                    stdid = True
                elif stdid:
                    STATEIDENTIFIERS = line.replace("#", "").split()
                    stdid = False
                continue
            else:
                if field > 0:
                    if first:
                        STATEATOMHEADER = ["atomI", "atomJ", "type1", "type2"]
                        first = False

                    state_line = {key: value for key, value in zip(STATEATOMHEADER, line.split())}
                    state_line.update({"NR": len(STATEATOMS) + 1})

                    final_state_line = {key: state_line[key] for key in state_line if ("type" not in key)}
                    states = {1: state_line["type1"], 2: state_line["type2"]}

                    final_state_line.update({"STATES": states})
                    STATEATOMS.append(atom_lam_pertubation_state_bond(**final_state_line))

                elif field == 0:
                    NPB = int(line.strip())
                field += 1

        self.NPB = NPB
        self.STATEIDENTIFIERS = STATEIDENTIFIERS
        self.STATEATOMHEADER = STATEATOMHEADER
        self.STATEATOMS = STATEATOMS

    """
    STR FUNCTIONS
    """

    def _state_STATEATOMHEADER_str(self):
        state_format_pattern = "{:>5} {:>5}" + "".join([" {:>5} " for x in range(self.NPTB)]) + ""
        return state_format_pattern.format(*self.STATEATOMHEADER)

    def block_to_string(self) -> str:
        result = self.name + self.line_seperator
        result += (
            "# NPB "
            + self.field_seperator
            + "NPTB = "
            + self.field_seperator
            + str(self.NPTB)
            + self.field_seperator
            + self.line_seperator
        )
        result += self.field_seperator + str(self.NPB) + self.line_seperator
        result += "# " + self._state_STATEATOMHEADER_str() + self.line_seperator
        result += "".join(map(str, sorted(self.STATEATOMS, key=lambda x: x.NR)))
        result += "END" + self.line_seperator
        return result


class PERTBONDSTRETCHH(_generic_gromos_block):
    def __init__(
        self,
        STATEATOMS: List[atom_lam_pertubation_state_bond] = None,
        STATEATOMHEADER: Tuple[str] = None,
        NPB: int = None,
        dummy_BOND=22,
        content: List[str] = None,
    ):
        self.NPTB = 2
        self.dummy_BOND = dummy_BOND

        if content is None:
            if STATEATOMHEADER is None:
                self.STATEATOMHEADER = ["atomI", "atomJ", "type1", "type2"]
            else:
                self.STATEATOMHEADER = STATEATOMHEADER

            if STATEATOMS is None:
                self.STATEATOMS = []
            else:
                self.STATEATOMS = []
                self.NPB = 0

                # self.add_state_atoms(STATEATOMS)
            super().__init__(used=True, name=__class__.__name__)
        else:
            super().__init__(used=True, name=__class__.__name__, content=content)

        # You can check yourself :)
        if NPB is not None and len(STATEATOMS) != NPB:
            raise ValueError(
                "NJLA must be equal to the length of STATEATOMS! NJLA="
                + str(NPB)
                + "\t stateatoms"
                + str(len(STATEATOMS))
                + "\n\n"
                + str(self)
            )

    @property
    def nStates(self) -> int:
        return self.NPTB

    @property
    def nTotalStateAtoms(self) -> int:
        return self.NPB

    @property
    def states(self) -> dict:
        return {
            self.STATEIDENTIFIERS[state - 1]: {
                atom.NR: atom.STATES[state] for atom in sorted(self.STATEATOMS, key=lambda x: x.NR)
            }
            for state in range(1, self.NPTB + 1)
        }

    def read_content_from_str(self, content: List[str]):
        field = 0
        NPB = None
        STATEIDENTIFIERS = None
        STATEATOMHEADER = None
        STATEATOMS = []
        first = True
        stdid = False
        for line in content:
            if "#" in line:
                if "state_identifiers" in line:
                    stdid = True
                elif stdid:
                    STATEIDENTIFIERS = line.replace("#", "").split()
                    stdid = False
                continue
            else:
                if field > 0:
                    if first:
                        STATEATOMHEADER = ["atomI", "atomJ", "type1", "type2"]
                        first = False

                    state_line = {key: value for key, value in zip(STATEATOMHEADER, line.split())}
                    state_line.update({"NR": len(STATEATOMS) + 1})

                    final_state_line = {key: state_line[key] for key in state_line if ("type" not in key)}
                    states = {1: state_line["type1"], 2: state_line["type2"]}

                    final_state_line.update({"STATES": states})
                    STATEATOMS.append(atom_lam_pertubation_state_bond(**final_state_line))

                elif field == 0:
                    NPB = int(line.strip())
                field += 1

        self.NPB = NPB
        self.STATEIDENTIFIERS = STATEIDENTIFIERS
        self.STATEATOMHEADER = STATEATOMHEADER
        self.STATEATOMS = STATEATOMS

    """
    STR FUNCTIONS
    """

    def _state_STATEATOMHEADER_str(self):
        state_format_pattern = "{:>5} {:>5}" + "".join([" {:>5} " for x in range(self.NPTB)]) + ""
        return state_format_pattern.format(*self.STATEATOMHEADER)

    def block_to_string(self) -> str:
        result = self.name + self.line_seperator
        result += (
            "# NPB "
            + self.field_seperator
            + "NPTB = "
            + self.field_seperator
            + str(self.NPTB)
            + self.field_seperator
            + self.line_seperator
        )
        result += self.field_seperator + str(self.NPB) + self.line_seperator
        result += "# " + self._state_STATEATOMHEADER_str() + self.line_seperator
        result += "".join(map(str, sorted(self.STATEATOMS, key=lambda x: x.NR)))
        result += "END" + self.line_seperator
        return result


# ANGLE


class PERTBONDANGLE(_generic_gromos_block):
    def __init__(
        self,
        STATEATOMS: List[atom_lam_pertubation_state_angle] = None,
        STATEATOMHEADER: Tuple[str] = None,
        NPA: int = None,
        dummy_ANGLE=22,
        content: List[str] = None,
    ):
        self.NPTB = 2
        self.dummy_ANGLE = dummy_ANGLE

        if content is None:
            if STATEATOMHEADER is None:
                self.STATEATOMHEADER = ["atomI", "atomJ", "atomK", "type1", "type2"]
            else:
                self.STATEATOMHEADER = STATEATOMHEADER

            if STATEATOMS is None:
                self.STATEATOMS = []
            else:
                self.STATEATOMS = []
                self.NPA = 0

                # self.add_state_atoms(STATEATOMS)
            super().__init__(used=True, name=__class__.__name__)
        else:
            super().__init__(used=True, name=__class__.__name__, content=content)

        # You can check yourself :)
        if NPA is not None and len(STATEATOMS) != NPA:
            raise ValueError(
                "NJLA must be equal to the length of STATEATOMS! NJLA="
                + str(NPA)
                + "\t stateatoms"
                + str(len(STATEATOMS))
                + "\n\n"
                + str(self)
            )

    @property
    def nStates(self) -> int:
        return self.NPTB

    @property
    def nTotalStateAtoms(self) -> int:
        return self.NPA

    @property
    def states(self) -> dict:
        return {
            self.STATEIDENTIFIERS[state - 1]: {
                atom.NR: atom.STATES[state] for atom in sorted(self.STATEATOMS, key=lambda x: x.NR)
            }
            for state in range(1, self.NPTB + 1)
        }

    def read_content_from_str(self, content: List[str]):
        field = 0
        NPA = None
        STATEIDENTIFIERS = None
        STATEATOMHEADER = None
        STATEATOMS = []
        first = True
        stdid = False
        for line in content:
            if "#" in line:
                if "state_identifiers" in line:
                    stdid = True
                elif stdid:
                    STATEIDENTIFIERS = line.replace("#", "").split()
                    stdid = False
                continue
            else:
                if field > 0:
                    if first:
                        STATEATOMHEADER = ["atomI", "atomJ", "atomK", "type1", "type2"]
                        first = False

                    state_line = {key: value for key, value in zip(STATEATOMHEADER, line.split())}
                    state_line.update({"NR": len(STATEATOMS) + 1})

                    final_state_line = {key: state_line[key] for key in state_line if ("type" not in key)}
                    states = {1: state_line["type1"], 2: state_line["type2"]}

                    final_state_line.update({"STATES": states})
                    STATEATOMS.append(atom_lam_pertubation_state_angle(**final_state_line))

                elif field == 0:
                    NPA = int(line.strip())
                field += 1

        self.NPA = NPA
        self.STATEIDENTIFIERS = STATEIDENTIFIERS
        self.STATEATOMHEADER = STATEATOMHEADER
        self.STATEATOMS = STATEATOMS

    """
    STR FUNCTIONS
    """

    def _state_STATEATOMHEADER_str(self):
        state_format_pattern = "{:>5} {:>5} {:>5}" + "".join([" {:>5} " for x in range(self.NPTB)]) + ""
        return state_format_pattern.format(*self.STATEATOMHEADER)

    def block_to_string(self) -> str:
        result = self.name + self.line_seperator
        result += (
            "# NPA "
            + self.field_seperator
            + "NPTB = "
            + self.field_seperator
            + str(self.NPTB)
            + self.field_seperator
            + self.line_seperator
        )
        result += self.field_seperator + str(self.NPA) + self.line_seperator
        result += "# " + self._state_STATEATOMHEADER_str() + self.line_seperator
        result += "".join(map(str, sorted(self.STATEATOMS, key=lambda x: x.NR)))
        result += "END" + self.line_seperator
        return result


class PERTBONDANGLEH(_generic_gromos_block):
    def __init__(
        self,
        STATEATOMS: List[atom_lam_pertubation_state_angle] = None,
        STATEATOMHEADER: Tuple[str] = None,
        NPA: int = None,
        dummy_ANGLE=22,
        content: List[str] = None,
    ):
        self.NPTB = 2
        self.dummy_ANGLE = dummy_ANGLE

        if content is None:
            if STATEATOMHEADER is None:
                self.STATEATOMHEADER = ["atomI", "atomJ", "atomK", "type1", "type2"]
            else:
                self.STATEATOMHEADER = STATEATOMHEADER

            if STATEATOMS is None:
                self.STATEATOMS = []
            else:
                self.STATEATOMS = []
                self.NPA = 0

                # self.add_state_atoms(STATEATOMS)
            super().__init__(used=True, name=__class__.__name__)
        else:
            super().__init__(used=True, name=__class__.__name__, content=content)

        # You can check yourself :)
        if NPA is not None and not len(STATEATOMS) == NPA:
            raise ValueError(
                "NJLA must be equal to the length of STATEATOMS! NJLA="
                + str(NPA)
                + "\t stateatoms"
                + str(len(STATEATOMS))
                + "\n\n"
                + str(self)
            )

    @property
    def nStates(self) -> int:
        return self.NPTB

    @property
    def nTotalStateAtoms(self) -> int:
        return self.NPA

    @property
    def states(self) -> dict:
        return {
            self.STATEIDENTIFIERS[state - 1]: {
                atom.NR: atom.STATES[state] for atom in sorted(self.STATEATOMS, key=lambda x: x.NR)
            }
            for state in range(1, self.NPTB + 1)
        }

    def read_content_from_str(self, content: List[str]):
        field = 0
        NPA = None
        STATEIDENTIFIERS = None
        STATEATOMHEADER = None
        STATEATOMS = []
        first = True
        stdid = False
        for line in content:
            if "#" in line:
                if "state_identifiers" in line:
                    stdid = True
                elif stdid:
                    STATEIDENTIFIERS = line.replace("#", "").split()
                    stdid = False
                continue
            else:
                if field > 0:
                    if first:
                        STATEATOMHEADER = ["atomI", "atomJ", "atomK", "type1", "type2"]
                        first = False

                    state_line = {key: value for key, value in zip(STATEATOMHEADER, line.split())}
                    state_line.update({"NR": len(STATEATOMS) + 1})

                    final_state_line = {key: state_line[key] for key in state_line if ("type" not in key)}
                    states = {1: state_line["type1"], 2: state_line["type2"]}

                    final_state_line.update({"STATES": states})
                    STATEATOMS.append(atom_lam_pertubation_state_angle(**final_state_line))

                elif field == 0:
                    NPA = int(line.strip())
                field += 1

        self.NPA = NPA
        self.STATEIDENTIFIERS = STATEIDENTIFIERS
        self.STATEATOMHEADER = STATEATOMHEADER
        self.STATEATOMS = STATEATOMS

    """
    STR FUNCTIONS
    """

    def _state_STATEATOMHEADER_str(self):
        state_format_pattern = "{:>5} {:>5} {:>5}" + "".join([" {:>5} " for x in range(self.NPTB)]) + ""
        return state_format_pattern.format(*self.STATEATOMHEADER)

    def block_to_string(self) -> str:
        result = self.name + self.line_seperator
        result += (
            "# NPA "
            + self.field_seperator
            + "NPTB = "
            + self.field_seperator
            + str(self.NPTB)
            + self.field_seperator
            + self.line_seperator
        )
        result += self.field_seperator + str(self.NPA) + self.line_seperator
        result += "# " + self._state_STATEATOMHEADER_str() + self.line_seperator
        result += "".join(map(str, sorted(self.STATEATOMS, key=lambda x: x.NR)))
        result += "END" + self.line_seperator
        return result


# DIHEDRAL


class PERTPROPERDIH(_generic_gromos_block):
    def __init__(
        self,
        STATEATOMS: List[atom_lam_pertubation_state_dihedral] = None,
        STATEATOMHEADER: Tuple[str] = None,
        NPD: int = None,
        dummy_DIH=22,
        content: List[str] = None,
    ):
        self.NPTB = 2
        self.dummy_DIH = dummy_DIH

        if content is None:
            if STATEATOMHEADER is None:
                self.STATEATOMHEADER = ["atomI", "atomJ", "atomK", "atomL", "type1", "type2"]
            else:
                self.STATEATOMHEADER = STATEATOMHEADER

            if STATEATOMS is None:
                self.STATEATOMS = []
            else:
                self.STATEATOMS = []
                self.NPD = 0

                # self.add_state_atoms(STATEATOMS)
            super().__init__(used=True, name=__class__.__name__)
        else:
            super().__init__(used=True, name=__class__.__name__, content=content)

        # You can check yourself :)
        if NPD is not None and len(STATEATOMS) != NPD:
            raise ValueError(
                "NJLA must be equal to the length of STATEATOMS! NJLA="
                + str(NPD)
                + "\t stateatoms"
                + str(len(STATEATOMS))
                + "\n\n"
                + str(self)
            )

    @property
    def nStates(self) -> int:
        return self.NPTB

    @property
    def nTotalStateAtoms(self) -> int:
        return self.NPD

    @property
    def states(self) -> dict:
        return {
            self.STATEIDENTIFIERS[state - 1]: {
                atom.NR: atom.STATES[state] for atom in sorted(self.STATEATOMS, key=lambda x: x.NR)
            }
            for state in range(1, self.NPTB + 1)
        }

    def read_content_from_str(self, content: List[str]):
        field = 0
        NPD = None
        STATEIDENTIFIERS = None
        STATEATOMHEADER = None
        STATEATOMS = []
        first = True
        stdid = False
        for line in content:
            if "#" in line:
                if "state_identifiers" in line:
                    stdid = True
                elif stdid:
                    STATEIDENTIFIERS = line.replace("#", "").split()
                    stdid = False
                continue
            else:
                if field > 0:
                    if first:
                        STATEATOMHEADER = ["atomI", "atomJ", "atomK", "atomL", "type1", "type2"]
                        first = False

                    state_line = {key: value for key, value in zip(STATEATOMHEADER, line.split())}
                    state_line.update({"NR": len(STATEATOMS) + 1})

                    final_state_line = {key: state_line[key] for key in state_line if ("type" not in key)}
                    states = {1: state_line["type1"], 2: state_line["type2"]}

                    final_state_line.update({"STATES": states})
                    STATEATOMS.append(atom_lam_pertubation_state_dihedral(**final_state_line))

                elif field == 0:
                    NPD = int(line.strip())
                field += 1

        self.NPD = NPD
        self.STATEIDENTIFIERS = STATEIDENTIFIERS
        self.STATEATOMHEADER = STATEATOMHEADER
        self.STATEATOMS = STATEATOMS

    """
    STR FUNCTIONS
    """

    def _state_STATEATOMHEADER_str(self):
        state_format_pattern = "{:>5} {:>5} {:>5} {:>5}" + "".join([" {:>5} " for x in range(self.NPTB)]) + ""
        return state_format_pattern.format(*self.STATEATOMHEADER)

    def block_to_string(self) -> str:
        result = self.name + self.line_seperator
        result += (
            "# NPD "
            + self.field_seperator
            + "NPTB = "
            + self.field_seperator
            + str(self.NPTB)
            + self.field_seperator
            + self.line_seperator
        )
        result += self.field_seperator + str(self.NPD) + self.line_seperator
        result += "# " + self._state_STATEATOMHEADER_str() + self.line_seperator
        result += "".join(map(str, sorted(self.STATEATOMS, key=lambda x: x.NR)))
        result += "END" + self.line_seperator
        return result


class PERTPROPERDIHH(_generic_gromos_block):
    def __init__(
        self,
        STATEATOMS: List[atom_lam_pertubation_state_dihedral] = None,
        STATEATOMHEADER: Tuple[str] = None,
        NPD: int = None,
        dummy_DIH=22,
        content: List[str] = None,
    ):
        self.NPTB = 2
        self.dummy_DIH = dummy_DIH

        if content is None:
            if STATEATOMHEADER is None:
                self.STATEATOMHEADER = ["atomI", "atomJ", "atomK", "atomL", "type1", "type2"]
            else:
                self.STATEATOMHEADER = STATEATOMHEADER

            if STATEATOMS is None:
                self.STATEATOMS = []
            else:
                self.STATEATOMS = []
                self.NPD = 0

                # self.add_state_atoms(STATEATOMS)
            super().__init__(used=True, name=__class__.__name__)
        else:
            super().__init__(used=True, name=__class__.__name__, content=content)

        # You can check yourself :)
        if NPD is not None and len(STATEATOMS) != NPD:
            raise ValueError(
                "NJLA must be equal to the length of STATEATOMS! NJLA="
                + str(NPD)
                + "\t stateatoms"
                + str(len(STATEATOMS))
                + "\n\n"
                + str(self)
            )

    @property
    def nStates(self) -> int:
        return self.NPTB

    @property
    def nTotalStateAtoms(self) -> int:
        return self.NPD

    @property
    def states(self) -> dict:
        return {
            self.STATEIDENTIFIERS[state - 1]: {
                atom.NR: atom.STATES[state] for atom in sorted(self.STATEATOMS, key=lambda x: x.NR)
            }
            for state in range(1, self.NPTB + 1)
        }

    def read_content_from_str(self, content: List[str]):
        field = 0
        NPD = None
        STATEIDENTIFIERS = None
        STATEATOMHEADER = None
        STATEATOMS = []
        first = True
        stdid = False
        for line in content:
            if "#" in line:
                if "state_identifiers" in line:
                    stdid = True
                elif stdid:
                    STATEIDENTIFIERS = line.replace("#", "").split()
                    stdid = False
                continue
            else:
                if field > 0:
                    if first:
                        STATEATOMHEADER = ["atomI", "atomJ", "atomK", "atomL", "type1", "type2"]
                        first = False

                    state_line = {key: value for key, value in zip(STATEATOMHEADER, line.split())}
                    state_line.update({"NR": len(STATEATOMS) + 1})

                    final_state_line = {key: state_line[key] for key in state_line if ("type" not in key)}
                    states = {1: state_line["type1"], 2: state_line["type2"]}

                    final_state_line.update({"STATES": states})
                    STATEATOMS.append(atom_lam_pertubation_state_dihedral(**final_state_line))

                elif field == 0:
                    NPD = int(line.strip())
                field += 1

        self.NPD = NPD
        self.STATEIDENTIFIERS = STATEIDENTIFIERS
        self.STATEATOMHEADER = STATEATOMHEADER
        self.STATEATOMS = STATEATOMS

    """
    STR FUNCTIONS
    """

    def _state_STATEATOMHEADER_str(self):
        state_format_pattern = "{:>5} {:>5} {:>5} {:>5}" + "".join([" {:>5} " for x in range(self.NPTB)]) + ""
        return state_format_pattern.format(*self.STATEATOMHEADER)

    def block_to_string(self) -> str:
        result = self.name + self.line_seperator
        result += (
            "# NPD "
            + self.field_seperator
            + "NPTB = "
            + self.field_seperator
            + str(self.NPTB)
            + self.field_seperator
            + self.line_seperator
        )
        result += self.field_seperator + str(self.NPD) + self.line_seperator
        result += "# " + self._state_STATEATOMHEADER_str() + self.line_seperator
        result += "".join(map(str, sorted(self.STATEATOMS, key=lambda x: x.NR)))
        result += "END" + self.line_seperator
        return result

    class PERTPROPERDIH(_generic_gromos_block):
        def __init__(
            self,
            STATEATOMS: List[atom_lam_pertubation_state_dihedral] = None,
            STATEATOMHEADER: Tuple[str] = None,
            NPD: int = None,
            dummy_DIH=22,
            content: List[str] = None,
        ):
            self.NPTB = 2
            self.dummy_DIH = dummy_DIH

            if content is None:
                if STATEATOMHEADER is None:
                    self.STATEATOMHEADER = ["atomI", "atomJ", "atomK", "atomL", "type1", "type2"]
                else:
                    self.STATEATOMHEADER = STATEATOMHEADER

                if STATEATOMS is None:
                    self.STATEATOMS = []
                else:
                    self.STATEATOMS = []
                    self.NPD = 0

                    # self.add_state_atoms(STATEATOMS)
                super().__init__(used=True, name=__class__.__name__)
            else:
                super().__init__(used=True, name=__class__.__name__, content=content)

            # You can check yourself :)
            if NPD is not None and len(STATEATOMS) != NPD:
                raise ValueError(
                    "NJLA must be equal to the length of STATEATOMS! NJLA="
                    + str(NPD)
                    + "\t stateatoms"
                    + str(len(STATEATOMS))
                    + "\n\n"
                    + str(self)
                )

        @property
        def nStates(self) -> int:
            return self.NPTB

        @property
        def nTotalStateAtoms(self) -> int:
            return self.NPD

        @property
        def states(self) -> dict:
            return {
                self.STATEIDENTIFIERS[state - 1]: {
                    atom.NR: atom.STATES[state] for atom in sorted(self.STATEATOMS, key=lambda x: x.NR)
                }
                for state in range(1, self.NPTB + 1)
            }

        def read_content_from_str(self, content: List[str]):
            field = 0
            NPD = None
            STATEIDENTIFIERS = None
            STATEATOMHEADER = None
            STATEATOMS = []
            first = True
            stdid = False
            for line in content:
                if "#" in line:
                    if "state_identifiers" in line:
                        stdid = True
                    elif stdid:
                        STATEIDENTIFIERS = line.replace("#", "").split()
                        stdid = False
                    continue
                else:
                    if field > 0:
                        if first:
                            STATEATOMHEADER = ["atomI", "atomJ", "atomK", "atomL", "type1", "type2"]
                            first = False

                        state_line = {key: value for key, value in zip(STATEATOMHEADER, line.split())}
                        state_line.update({"NR": len(STATEATOMS) + 1})

                        final_state_line = {key: state_line[key] for key in state_line if ("type" not in key)}
                        states = {1: state_line["type1"], 2: state_line["type2"]}

                        final_state_line.update({"STATES": states})
                        STATEATOMS.append(atom_lam_pertubation_state_dihedral(**final_state_line))

                    elif field == 0:
                        NPD = int(line.strip())
                    field += 1

            self.NPD = NPD
            self.STATEIDENTIFIERS = STATEIDENTIFIERS
            self.STATEATOMHEADER = STATEATOMHEADER
            self.STATEATOMS = STATEATOMS

        """
        STR FUNCTIONS
        """

        def _state_STATEATOMHEADER_str(self):
            state_format_pattern = "{:>5} {:>5} {:>5} {:>5}" + "".join([" {:>5} " for x in range(self.NPTB)]) + ""
            return state_format_pattern.format(*self.STATEATOMHEADER)

        def block_to_string(self) -> str:
            result = self.name + self.line_seperator
            result += (
                "# NPD "
                + self.field_seperator
                + "NPTB = "
                + self.field_seperator
                + str(self.NPTB)
                + self.field_seperator
                + self.line_seperator
            )
            result += self.field_seperator + str(self.NPD) + self.line_seperator
            result += "# " + self._state_STATEATOMHEADER_str() + self.line_seperator
            result += "".join(map(str, sorted(self.STATEATOMS, key=lambda x: x.NR)))
            result += "END" + self.line_seperator
            return result


# IMPROPER


class PERTIMROPERDIH(_generic_gromos_block):
    def __init__(
        self,
        STATEATOMS: List[atom_lam_pertubation_state_dihedral] = None,
        STATEATOMHEADER: Tuple[str] = None,
        NPD: int = None,
        dummy_IMP=22,
        content: List[str] = None,
    ):
        self.NPTB = 2
        self.dummy_DIH = dummy_IMP

        if content is None:
            if STATEATOMHEADER is None:
                self.STATEATOMHEADER = ["atomI", "atomJ", "atomK", "atomL", "type1", "type2"]
            else:
                self.STATEATOMHEADER = STATEATOMHEADER

            if STATEATOMS is None:
                self.STATEATOMS = []
            else:
                self.STATEATOMS = []
                self.NPD = 0

                # self.add_state_atoms(STATEATOMS)
            super().__init__(used=True, name=__class__.__name__)
        else:
            super().__init__(used=True, name=__class__.__name__, content=content)

        # You can check yourself :)
        if NPD is not None and len(STATEATOMS) != NPD:
            raise ValueError(
                "NPD must be equal to the length of STATEATOMS! NPD="
                + str(NPD)
                + "\t stateatoms"
                + str(len(STATEATOMS))
                + "\n\n"
                + str(self)
            )

    @property
    def nStates(self) -> int:
        return self.NPTB

    @property
    def nTotalStateAtoms(self) -> int:
        return self.NPD

    @property
    def states(self) -> dict:
        return {
            self.STATEIDENTIFIERS[state - 1]: {
                atom.NR: atom.STATES[state] for atom in sorted(self.STATEATOMS, key=lambda x: x.NR)
            }
            for state in range(1, self.NPTB + 1)
        }

    def read_content_from_str(self, content: List[str]):
        field = 0
        NPD = None
        STATEIDENTIFIERS = None
        STATEATOMHEADER = None
        STATEATOMS = []
        first = True
        stdid = False
        for line in content:
            if "#" in line:
                if "state_identifiers" in line:
                    stdid = True
                elif stdid:
                    STATEIDENTIFIERS = line.replace("#", "").split()
                    stdid = False
                continue
            else:
                if field > 0:
                    if first:
                        STATEATOMHEADER = ["atomI", "atomJ", "atomK", "atomL", "type1", "type2"]
                        first = False

                    state_line = {key: value for key, value in zip(STATEATOMHEADER, line.split())}
                    state_line.update({"NR": len(STATEATOMS) + 1})

                    final_state_line = {key: state_line[key] for key in state_line if ("type" not in key)}
                    states = {1: state_line["type1"], 2: state_line["type2"]}

                    final_state_line.update({"STATES": states})
                    STATEATOMS.append(atom_lam_pertubation_state_dihedral(**final_state_line))

                elif field == 0:
                    NPD = int(line.strip())
                field += 1

        self.NPD = NPD
        self.STATEIDENTIFIERS = STATEIDENTIFIERS
        self.STATEATOMHEADER = STATEATOMHEADER
        self.STATEATOMS = STATEATOMS

    """
    STR FUNCTIONS
    """

    def _state_STATEATOMHEADER_str(self):
        state_format_pattern = "{:>5} {:>5} {:>5} {:>5}" + "".join([" {:>5} " for x in range(self.NPTB)]) + ""
        return state_format_pattern.format(*self.STATEATOMHEADER)

    def block_to_string(self) -> str:
        result = self.name + self.line_seperator
        result += (
            "# NPD "
            + self.field_seperator
            + "NPTB = "
            + self.field_seperator
            + str(self.NPTB)
            + self.field_seperator
            + self.line_seperator
        )
        result += self.field_seperator + str(self.NPD) + self.line_seperator
        result += "# " + self._state_STATEATOMHEADER_str() + self.line_seperator
        result += "".join(map(str, sorted(self.STATEATOMS, key=lambda x: x.NR)))
        result += "END" + self.line_seperator
        return result


class PERTIMROPERDIHH(_generic_gromos_block):
    def __init__(
        self,
        STATEATOMS: List[atom_lam_pertubation_state_dihedral] = None,
        STATEATOMHEADER: Tuple[str] = None,
        NPD: int = None,
        dummy_IMP=22,
        content: List[str] = None,
    ):
        self.NPTB = 2
        self.dummy_DIH = dummy_IMP

        if content is None:
            if STATEATOMHEADER is None:
                self.STATEATOMHEADER = ["atomI", "atomJ", "atomK", "atomL", "type1", "type2"]
            else:
                self.STATEATOMHEADER = STATEATOMHEADER

            if STATEATOMS is None:
                self.STATEATOMS = []
            else:
                self.STATEATOMS = []
                self.NPD = 0

                # self.add_state_atoms(STATEATOMS)
            super().__init__(used=True, name=__class__.__name__)
        else:
            super().__init__(used=True, name=__class__.__name__, content=content)

        # You can check yourself :)
        if NPD is not None and len(STATEATOMS) != NPD:
            raise ValueError(
                "NPD must be equal to the length of STATEATOMS! NPD="
                + str(NPD)
                + "\t stateatoms"
                + str(len(STATEATOMS))
                + "\n\n"
                + str(self)
            )

    @property
    def nStates(self) -> int:
        return self.NPTB

    @property
    def nTotalStateAtoms(self) -> int:
        return self.NPD

    @property
    def states(self) -> dict:
        return {
            self.STATEIDENTIFIERS[state - 1]: {
                atom.NR: atom.STATES[state] for atom in sorted(self.STATEATOMS, key=lambda x: x.NR)
            }
            for state in range(1, self.NPTB + 1)
        }

    def read_content_from_str(self, content: List[str]):
        field = 0
        NPD = None
        STATEIDENTIFIERS = None
        STATEATOMHEADER = None
        STATEATOMS = []
        first = True
        stdid = False
        for line in content:
            if "#" in line:
                if "state_identifiers" in line:
                    stdid = True
                elif stdid:
                    STATEIDENTIFIERS = line.replace("#", "").split()
                    stdid = False
                continue
            else:
                if field > 0:
                    if first:
                        STATEATOMHEADER = ["atomI", "atomJ", "atomK", "atomL", "type1", "type2"]
                        first = False

                    state_line = {key: value for key, value in zip(STATEATOMHEADER, line.split())}
                    state_line.update({"NR": len(STATEATOMS) + 1})

                    final_state_line = {key: state_line[key] for key in state_line if ("type" not in key)}
                    states = {1: state_line["type1"], 2: state_line["type2"]}

                    final_state_line.update({"STATES": states})
                    STATEATOMS.append(atom_lam_pertubation_state_dihedral(**final_state_line))

                elif field == 0:
                    NPD = int(line.strip())
                field += 1

        self.NPD = NPD
        self.STATEIDENTIFIERS = STATEIDENTIFIERS
        self.STATEATOMHEADER = STATEATOMHEADER
        self.STATEATOMS = STATEATOMS

    """
    STR FUNCTIONS
    """

    def _state_STATEATOMHEADER_str(self):
        state_format_pattern = "{:>5} {:>5} {:>5} {:>5}" + "".join([" {:>5} " for x in range(self.NPTB)]) + ""
        return state_format_pattern.format(*self.STATEATOMHEADER)

    def block_to_string(self) -> str:
        result = self.name + self.line_seperator
        result += (
            "# NPD "
            + self.field_seperator
            + "NPTB = "
            + self.field_seperator
            + str(self.NPTB)
            + self.field_seperator
            + self.line_seperator
        )
        result += self.field_seperator + str(self.NPD) + self.line_seperator
        result += "# " + self._state_STATEATOMHEADER_str() + self.line_seperator
        result += "".join(map(str, sorted(self.STATEATOMS, key=lambda x: x.NR)))
        result += "END" + self.line_seperator
        return result
