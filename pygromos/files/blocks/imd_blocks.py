import warnings

from pygromos.files.blocks._general_blocks import TITLE as generic_TITLE
from pygromos.files.blocks._general_blocks import _generic_gromos_block
from pygromos.utils.typing import List, Dict, Union, Number

# forward declarations
TITLE: generic_TITLE = generic_TITLE


class _generic_imd_block(_generic_gromos_block):
    name = "genericBlock"
    _order: List[List[str]]  # contains the ordering of all fields in a block

    def __init__(self, used: bool, content: Dict = None):
        super().__init__(name=self.name, used=used, content=content)

    def block_to_string(self) -> str:
        result = ""
        result += str(self.name) + self.line_seperator
        for line in self._order:
            for field_names in line:
                result += "#" + self.field_seperator + self.field_seperator.join(field_names) + self.line_seperator
                result += self.field_seperator
                for element in field_names:
                    element = element.split("(")[0]
                    element = element.split(":")[0]
                    element = element.split(" ")[0]
                    element = element.replace(" ", "")
                    attribute = self.__getattribute__(element)

                    try:
                        if isinstance(attribute, (str, Number, bool)):  # One element field
                            if isinstance(attribute, (bool, int)):
                                attribute = int(attribute)
                            elif isinstance(attribute, float):  # supress scientific notation for floats!
                                attribute = format(attribute, "f")
                            result += str(attribute) + self.field_seperator
                        elif isinstance(attribute, List) and all(
                            [isinstance(x, (str, Number)) for x in attribute]
                        ):  # list content
                            if all([isinstance(x, str) for x in attribute]):
                                result += self.field_seperator.join(attribute) + self.field_seperator
                            elif all([isinstance(x, Number) for x in attribute]):
                                result += self.field_seperator.join(map(str, attribute)) + self.field_seperator
                            else:
                                raise ValueError(
                                    "could not Interpret list:  " + str(element) + "\t\n" + str(attribute) + "\nEOF\n"
                                )
                        elif isinstance(attribute, List) and all([isinstance(x, List) for x in attribute]):  # matrices
                            pre_parsed_rows = map(lambda x: self.field_seperator.join(map(str, x)), attribute)
                            result += (self.line_seperator + self.field_seperator).join(
                                pre_parsed_rows
                            ) + self.field_seperator
                        else:
                            result += (self.line_seperator + self.field_seperator).join(
                                map(self.field_seperator.join, attribute)
                            ) + self.field_seperator
                    except ValueError:
                        raise ValueError(
                            "Could not convert attribute "
                            + str(element)
                            + " to string!\n value of attribute was: \n"
                            + str(attribute)
                            + "\nEOF\n"
                        )

                result += self.line_seperator
        result += "END\n"
        return result

    def read_content_from_str(self, content: List[str]):
        keyLineNumb = 0
        contentlines = []
        for line in content:
            if line.startswith("# "):
                # a new key line was found add values to class
                if len(contentlines) > 0:
                    self._parse_key_content(keyLineNumb=keyLineNumb, contentlines=contentlines)
                    keyLineNumb += 1  # go next key line
                    contentlines = []  # reset content storage
            else:
                # no key found -> add data to list to be later added to key
                fields = line.split(self.field_seperator)
                while "" in fields:
                    fields.remove("")
                contentlines.append(fields)
        # add fianl key after loop
        self._parse_key_content(keyLineNumb=keyLineNumb, contentlines=contentlines)

    def _parse_key_content(self, keyLineNumb: int = 0, contentlines: List[str] = []):
        # print("----------------")
        # print(contentlines)
        # print(self._order[0])
        # print(keyLineNumb)
        if len(contentlines) == 1:  # it's a one liner
            if len(contentlines[0]) == len(self._order[0][keyLineNumb]):
                # found key-value match
                for key, field in zip(self._order[0][keyLineNumb], contentlines[0]):
                    # first bring key to attribute form
                    key = self._clean_key_from_order(key=key)
                    field = self._try_to_convert_field(field=field, key=key)
                    setattr(self, key, field)
        else:  # parse multi line
            if len(self._order[0][keyLineNumb]) == 1:
                # one key but multiple fields
                key = self._clean_key_from_order(key=self._order[0][keyLineNumb][0])
                field = self._try_to_convert_field(field=contentlines, key=key)
                setattr(self, key, field)
            else:
                # set attribute with a list over multiple lines for keys in paralles (example temp baths)
                for i, key in enumerate(self._order[0][keyLineNumb]):
                    field = [x[i] for x in contentlines]
                    key = self._clean_key_from_order(key=key)
                    field = self._try_to_convert_field(field=field, key=key)
                    setattr(self, key, field)

    def _clean_key_from_order(self, key: str) -> str:
        # key = key.replace("("," ").replace("."," ")
        # key.split(" ")
        # return key[0]
        if "(" in key:
            return key.split("(")[0]
        if "." in key:
            return key.split(".")[0]
        if " " in key:
            return key.split(" ")[0]
        return key

    def _try_to_convert_field(self, field: str, key: str) -> any:
        if isinstance(field, str):
            try:
                if self.__annotations__[key] is bool:
                    return bool(int(field))  # solves the type problem for simple types
                return self.__annotations__[key](field)  # solves the type problem for simple types
            except Exception:
                return field
        elif isinstance(field, list):
            if all([isinstance(x, str) for x in field]):
                try:
                    return [float(x) for x in field]
                except Exception:
                    return field
            elif all([isinstance(x, list) for x in field]) and all([[isinstance(x, list) for x in y] for y in field]):
                try:
                    return [[float(x) for x in y] for y in field]
                except Exception:
                    return field


class SYSTEM(_generic_imd_block):
    """System Block

        The system block defines the number of solute molecules and solvent molecules

    Attributes
    ----------
    NPM:    int
        Number of Solute Molecules
    NSM:    int
        Number of Solvent Molecules



    """

    name: str = "SYSTEM"

    # fields
    NPM: int  # number of Solute Molecules
    NSM: int  # number of Solvent Molecules

    _order = [[["NPM", "NSM"]]]

    def __init__(self, NPM: int = 0, NSM: int = 0, content: List[str] = None):
        super().__init__(used=True, content=content)
        if content is None:
            self.NPM = int(NPM)
            self.NSM = int(NSM)


class STEP(_generic_imd_block):
    """Step Block

        This Block gives the number of simulation steps,

    Attributes
    -----------
    NSTLIM: int
        number of simulations Step till terminating.
    T:  float
        Starting Time
    DT: float
        time step [fs]

    """

    name: str = "STEP"

    # fields
    NSTLIM: int
    T: float
    DT: float

    _order = [[["NSTLIM", "T", "DT"]]]

    def __init__(self, NSTLIM: int = 0, T: float = 0, DT: float = 0, content: List[str] = None):
        super().__init__(used=True, content=content)
        if content is None:
            self.NSTLIM = int(NSTLIM)
            self.T = float(T)
            self.DT = float(DT)


class REPLICA(_generic_imd_block):
    """REPLICA Block

        This block is controlling the REPLICA settings in gromos. It works both for T-REMD and H-REMD

    Attributes
    ----------
    RETL:   bool
        Do you want to rund an RE simulation or not?
        RETL=0 : turn off RE
        RETL=1 : turn on RE
    NRET:   int
        Number of temperatures
    RET:    List[float]
        List of temperature at each replica
    LRESCALE:   bool
        controls temperature scaling
        0 don't scale temperatures after exchange trial
        1 scale temperatures after exchange trial
    NRELAM: int
        number of lambda values
    RELAM:  List[float]
        lambda value for each lambda-replica
    RETS:   List[float]
        timestep of each LAMBDA replica
    NRETRIAL:   int
        How many replica exchanges trials should be executed? (NRETRIAL*STEP.NSTLIM == total simulation time)
    NREQUIL:    int
        How many equilibration runs shall be exectured? (NREQUIL*STEP.NSTLIM == total simulation time)
    CONT:   bool
        Is this a continuation run?
        0 start from one configuration file
        1 start from multiple configuration files


    Example Block
    -------------
    # RETL
    1
    # NRET
    10
    # RET(1..NRET)
    300.0  320.0  340.0 360.0 380.0
    400.0  420.0  440.0 460.0 480.0
    # LRESCALE
    1
    # NRELAM
    10
    # RELAM(1..NRELAM
    0.0    0.1    0.2   0.3   0.4
    0.5    0.6    0.7   0.8   0.9
    # RETS(1..NRELAM)
    0.002  0.001  0.001 0.001 0.002
    0.003  0.002  0.001 0.001 0.002
    # NRETRIAL
    100
    # NREQUIL
    10
    # CONT
    0
    END

    """

    name: str = "REPLICA"

    RETL: bool

    NRET: int
    RET: List[float]

    LRESCALE: bool

    NRELAM: int
    RELAM: List[float]
    RETS: List[float]

    NRETRIAL: int
    NREQUIL: int
    CONT: bool

    _order = [
        [
            ["RETL"],
            ["NRET"],
            ["RET(1 ... NRET)"],
            ["LRESCALE"],
            ["NRELAM"],
            ["RELAM(1 ... NRELAM)"],
            ["RETS(1 ... NRELAM)"],
            ["NRETRIAL", "NREQUIL", "CONT"],
        ]
    ]

    def __init__(
        self,
        RETL: bool = False,
        NRET: int = 0,
        RET: List[float] = [],
        LRESCALE: int = 0,
        NRELAM: int = 0,
        RELAM: List[float] = [],
        RETS: List[float] = [],
        NRETRIAL: int = 0,
        NREQUIL: int = 0,
        CONT: bool = False,
        content: List[str] = None,
    ):
        super().__init__(used=True, content=content)
        if content is None:
            self.RETL = bool(RETL)

            self.NRET = int(NRET)
            self.RET = RET

            self.LRESCALE = bool(LRESCALE)

            self.NRELAM = int(NRELAM)
            self.RELAM = RELAM
            self.RETS = RETS

            self.NRETRIAL = int(NRETRIAL)
            self.NREQUIL = int(NREQUIL)
            self.CONT = bool(CONT)

    def read_content_from_str(self, content: List[str]):
        try:
            setattr(self, "RETL", int(content[1].split()[0]))
            setattr(self, "NRET", int(content[3].split()[0]))
            T_values = list(map(float, content[5].split()))
            if len(T_values) == self.NRET:
                setattr(self, "RET", T_values)
            else:
                raise IOError("REPLICA: NRET was not equal to the number of Temperatures (RET) in IMD!")
            setattr(self, "LRESCALE", int(content[7].split()[0]))
            setattr(self, "NRELAM", int(content[9].split()[0]))
            lambda_val = list(map(float, content[11].split()))
            if len(lambda_val) == self.NRELAM:
                setattr(self, "RELAM", lambda_val)
            else:
                raise IOError("REPLICA: NRELAM was not equal to the number of lambda values (RELAM) in IMD!")
            lambda_timestep = list(map(float, content[13].split()))
            if len(lambda_timestep) == self.NRELAM:
                setattr(self, "RETS", lambda_timestep)
            else:
                raise IOError("REPLICA: NRELAM was not equal to the number of lambda timestep values (RETS) in IMD!")
            [setattr(self, key, int(value)) for key, value in zip(self._order[0][-1], content[-1].split())]
        except Exception as err:
            raise IOError("Could not parse block from str - " + __class__.__name__ + "\n" + str(err.args))


class NEW_REPLICA_EDS(_generic_imd_block):
    """REPLICA_EDS Block

        This block is controlling the REPLICA_EDS settings in  gromos and is basically a mixture of EDS and RE block. (Don't use them when using this block!)

    Attributes
    ----------
    REEDS:  bool
        Shall REEDS be activated?
    NRES:   int
        Number of s-Values
    NEOFF:  int
        number of replica exchange eds energy Offset vectors (only used for 2D REEDS - still needs to be present ;))
    NUMSTATES:  int
        Number of EDS-states

    RES:    List[float]
        s_values for all replicas
    EIR:    List[List[float]]
        energy offsets for all replicas and all states  List[List[float]] = REPLICA[EDS_STATE[EIR]]
    NERTRIAL: int
        How many replica exchanges trials should be executed? (NRETRIAL*STEP.NSTLIM == total simulation time)
    NREQUIL: int
        How many equilibration runs shall be exectured? (NREQUIL*STEP.NSTLIM == total simulation time)
    EDS_STAT_OUT: int
        Shall the replica exchange information be outputted? (__future__ frequency of output.)
    CONT: bool
        Is this a continuation run?
    """

    name: str = "REPLICA_EDS"

    REEDS: bool

    NRES: int
    NEOFF: int
    NUMSTATES: int

    RES: List[float]
    EIR: List[float]

    NRETRIAL: int
    NREQUIL: int
    EDS_STAT_OUT: int
    CONT: bool
    PERIODIC: int

    _order = [
        [
            ["REEDS"],
            ["NRES", "NUMSTATES", "NEOFF"],
            ["RES(1 ... NRES)"],
            ["EIR(NUMSTATES x NRES)"],
            ["NRETRIAL", "NREQUIL", "CONT", "EDS_STAT_OUT", "PERIODIC"],
        ]
    ]

    def __init__(
        self,
        REEDS: bool = False,
        NRES: int = 0,
        NUMSTATES: int = 0,
        NEOFF: int = 0,
        RES: List[float] = [],
        EIR: List[List[float]] = [[]],
        NRETRIAL: int = 0,
        NREQUIL: int = 0,
        EDS_STAT_OUT: int = 0,
        CONT: bool = False,
        PERIODIC: int = 0,
        content: List[str] = None,
    ):
        super().__init__(used=True, content=content)
        if content is None:
            self.REEDS = bool(REEDS)

            self.NRES = int(NRES)
            self.NEOFF = int(NEOFF)
            self.NUMSTATES = int(NUMSTATES)

            self.RES = RES
            self.EIR = EIR

            self.NRETRIAL = int(NRETRIAL)
            self.NREQUIL = int(NREQUIL)
            self.CONT = bool(CONT)
            self.EDS_STAT_OUT = int(EDS_STAT_OUT)
            self.PERIODIC = int(PERIODIC)

    def read_content_from_str(self, content: List[str]):
        try:
            setattr(self, "REEDS", int(content[1].split()[0]))
            setattr(self, "NRES", int(content[3].split()[0]))
            setattr(self, "NEOFF", int(content[3].split()[1]))
            setattr(self, "NUMSTATES", int(content[3].split()[2]))
            s_values = list(map(float, content[5].split()))
            if len(s_values) == self.NRES:
                setattr(self, "RES", s_values)
            else:
                raise IOError("REPLICA_EDS: NRES was not equal to the number of s-values in IMD!")
            EIR = []
            for ind in range(7, 7 + self.NUMSTATES):
                EIR_line = list(map(float, content[ind].split()))
                if len(EIR_line) != self.NRES:
                    raise IOError("REPLICA_EDS: NRES was not equal to the number of EIRs given in IMD!")
                EIR.append(EIR_line)
            setattr(self, "EIR", EIR)
            [setattr(self, key, int(value)) for key, value in zip(self._order[0][-1], content[-1].split())]

        except Exception as err:
            raise IOError("Could not parse block from str - " + __class__.__name__ + "\n" + str(err.args))


class REPLICA_EDS(_generic_imd_block):
    name: str = "REPLICA_EDS"

    REEDS: bool

    NRES: int
    NUMSTATES: int

    RES: List[float]
    EIR: List[float]

    NRETRIAL: int
    NREQUIL: int
    EDS_STAT_OUT: int
    CONT: bool

    _order = [
        [
            ["REEDS"],
            ["NRES", "NUMSTATES"],
            ["RES(1 ... NRES)"],
            ["EIR(NUMSTATES x NRES)"],
            ["NRETRIAL", "NREQUIL", "CONT", "EDS_STAT_OUT"],
        ]
    ]

    def __init__(
        self,
        REEDS: bool = True,
        NRES: int = 0,
        NUMSTATES: int = 0,
        RES: List[float] = [0],
        EIR: List[List[float]] = [[0]],
        NRETRIAL: int = 0,
        NREQUIL: int = 0,
        EDS_STAT_OUT: int = 0,
        CONT: bool = True,
        content: List[str] = None,
    ):
        """REPLICA_EDS Block

            This block is controlling the REPLICA_EDS settings in  gromos and is basically a mixture of EDS and RE block. (Don't use them when using this block!)

        Attributes
        ----------
        REEDS:  bool
            Shall REEDS be activated?
        NRES:   int
            Number of s-Values
        NUMSTATES:  int
            Number of EDS-states

        RES:    List[float]
            s_values for all replicas
        EIR:    List[List[float]]
            energy offsets for all replicas and all states  List[List[float]] = REPLICA[EDS_STATE[EIR]]
        NERTRIAL: int
            How many replica exchanges trials should be executed? (NRETRIAL*STEP.NSTLIM == total simulation time)
        NREQUIL: int
            How many equilibration runs shall be exectured? (NREQUIL*STEP.NSTLIM == total simulation time)
        EDS_STAT_OUT: int
            Shall the replica exchange information be outputted? (__future__ frequency of output.)
        CONT: bool
            Is this a continuation run?
        """
        super().__init__(used=True, content=content)
        if content is None:
            self.REEDS = REEDS

            self.NRES = NRES
            self.NUMSTATES = NUMSTATES

            self.RES = RES
            self.EIR = EIR

            self.NRETRIAL = NRETRIAL
            self.NREQUIL = NREQUIL
            self.CONT = CONT
            self.EDS_STAT_OUT = EDS_STAT_OUT

    def read_content_from_str(self, content: List[str]):
        try:
            setattr(self, "REEDS", int(content[1].split()[0]))
            setattr(self, "NRES", int(content[3].split()[0]))
            setattr(self, "NUMSTATES", int(content[3].split()[1]))
            s_values = list(map(float, content[5].split()))
            if len(s_values) == self.NRES:
                setattr(self, "RES", s_values)
            else:
                raise IOError("REPLICA_EDS: NRES was not equal to the number of s-values in IMD!")
            EIR = []
            for ind in range(7, 7 + self.NUMSTATES):
                EIR_line = list(map(float, content[ind].split()))
                if len(EIR_line) != self.NRES:
                    raise IOError("REPLICA_EDS: NRES was not equal to the number of EIRs given in IMD!")
                EIR.append(EIR_line)
            setattr(self, "EIR", EIR)
            [setattr(self, key, int(value)) for key, value in zip(self._order[0][-1], content[-1].split())]

        except Exception as err:
            raise IOError("Could not parse block from str - " + __class__.__name__ + "\n" + str(err.args))


class BOUNDCOND(_generic_imd_block):
    """Boundary Condition Block

        This block describes the boundary condition of the coordinate system.

    Attributes
    ----------
        NTB : int, optional
            Boundary conditions, by default 0
            -1 : truncated octahedral
             0 : vacuum
             1 : rectangular
             2 : triclinic
        NDFMIN : int, optional
            number of degrees of freedom subtracted for temperature, by default 0

    """

    name: str = "BOUNDCOND"

    NTB: int
    NDFMIN: int

    _order = [[["NTB", "NDFMIN"]]]

    def __init__(self, NTB: int = 0, NDFMIN: int = 0, content: List[str] = None):
        super().__init__(used=True, content=content)
        if content is None:
            self.NTB = int(NTB)
            self.NDFMIN = int(NDFMIN)


class STOCHDYN(_generic_imd_block):
    """Stochastic Dynamics block

        This block is for Stochastic Dynamics

    Attributes
    ----------
    NTSD:    int
    NTFR:    int
    NSFR:    int
    NBREF:   int
    RCUTF:   float
    CFRIC:   float
    TEMPSD:  float


    """

    name: str = "STOCHDYN"

    NTSD: int
    NTFR: int
    NSFR: int
    NBREF: int
    RCUTF: float
    CFRIC: float
    TEMPSD: float

    _order = [[["NTSD", "NTFR", "NSFR", "NBREF", "RCUTF", "CFRIC", "TEMPSD"]]]

    def __init__(
        self,
        NTSD: int = 0,
        NTFR: int = 0,
        NSFR: int = 0,
        NBREF: int = 0,
        RCUTF: float = 0,
        CFRIC: float = 0,
        TEMPSD: float = 0,
        content: List[str] = None,
    ):
        super().__init__(used=True, content=content)
        if content is None:
            self.NTSD = int(NTSD)
            self.NTFR = int(NTFR)
            self.NSFR = int(NSFR)
            self.NBREF = int(NBREF)
            self.RCUTF = float(RCUTF)
            self.CFRIC = float(CFRIC)
            self.TEMPSD = float(TEMPSD)


class PERTURBATION(_generic_imd_block):
    """Pertubation Block

        This block is for Thermodynamic integration


    Attributes
    ----------
    NTG:    int
        0..1 controls use of free-energy calculation.
              0: no free-energy calculation (default)
              1: calculate dH/dRLAM
    NRDGL: int
        0,1 controls reading of initial value for RLAM.
             0: use initial RLAM parameter from PERTURBATION block
             1: read from configuration
    RLAM: float
        0.0..1.0 initial value for lambda
    DLAMT: float
        >= 0.0 rate of lambda increase in time.
    ALPHLJ: float
        >= 0.0 Lennard-Jones soft-core parameter
    ALPHC: float
        >= 0.0 Coulomb-RF soft-core parameter
    NLAM: int
        > 0 power dependence of lambda coupling
    NSCALE: int
        0..2 controls use of interaction scaling
            0: no interaction scaling
            1: interaction scaling
            2: perturbation for all atom pairs with scaled
               interactions. No perturbation for others.

    ExampleBlock
    ____________
    #     NTG   NRDGL    RLAM   DLAMT
            1   0         0.0     0.0
    #  ALPHLJ   ALPHC    NLAM  NSCALE
        0.5     0.5         2       0

    """

    name: str = "PERTURBATION"

    NTG: int
    NRDGL: int
    RLAM: float
    DLAMT: float
    ALPHLJ: float
    ALPHC: float
    NLAM: int
    NSCALE: int

    _order = [[["NTG", "NRDGL", "RLAM", "DLAMT"], ["ALPHLJ", "ALPHC", "NLAM", "NSCALE"]]]

    def __init__(
        self,
        NTG: int = 0,
        NRDGL: int = 0,
        RLAM: float = 0,
        DLAMT: float = 0,
        ALPHLJ: float = 0,
        ALPHC: float = 0,
        NLAM: int = 0,
        NSCALE: int = 0,
        content: List[str] = None,
    ):
        super().__init__(used=True, content=content)
        if content is None:
            self.NTG = int(NTG)
            self.NRDGL = int(NRDGL)
            self.RLAM = float(RLAM)
            self.DLAMT = float(DLAMT)
            self.ALPHLJ = float(ALPHLJ)
            self.ALPHC = float(ALPHC)
            self.NLAM = int(NLAM)
            self.NSCALE = int(NSCALE)


class PRECALCLAM(_generic_imd_block):
    """
    Attributes
    ----------
    NTG:    int


    Returns
    -------

    """

    name = "PRECALCLAM"
    NRLAM: int
    MINLAM: float
    MAXLAM: float

    _order = [[["NRLAM", "MINLAM", "MAXLAM"]]]

    def __init__(self, NRLAM: int = 0, MINLAM: float = 0, MAXLAM: float = 0, content: List[str] = None):
        """
            Can be used to caluclate multiple extra lambda values
        Parameters
        ----------
        NRLAM: int
                0  : off
                >1 : precalculating energies for NRLAM extra lambda values
        MINLAM
              between 0 and 1: minimum lambda value to precalculate energies
        MAXLAM
            between MINLAM and 1: maximum lambda value to precalculate energies
        """
        super().__init__(used=True, content=content)
        if content is None:
            self.NRLAM = int(NRLAM)
            self.MINLAM = float(MINLAM)
            self.MAXLAM = float(MAXLAM)


class MULTIBATH(_generic_imd_block):
    """MULTIBATH Block

        This block controls the Thermostat of a simulation. Multiple temperature baths can be coupled.

    Attributes
    ----------
    ALGORITHM:  int
        temperature coupling algorithm
        weak-coupling(0)
        nose-hoover(1)
        nose-hoover-chains(2)	num
        (where num is the number of chains to use)
    NUM: int, optional
        Mumber of chains in Nosé Hoover chains scheme [only specify when needed]
    NBATHS: int
        Number of temperature baths
    TEMP0:  List[float]
        temperature of each bath (list len == NBATHS)
    TAU:    List[float]
            coupling of the temperaturebath to the system  (list len == NBATHS)
    DOFSET: int
        Number of set of Degrees of freedom.
    LAST:   List[int]   (list len == DOFSET)
        Last atoms of each DOFSet
    COMBATH:    List[int]   (list len == DOFSET)
        Index of the temperature baths
    IRBATH: List[int]   (list len == DOFSET)
        IRBAHT?

    See Also
    --------
        _generic_imd_block, _generic_gromos_block

    """

    name: str = "MULTIBATH"

    ALGORITHM: int
    NUM: int
    NBATHS: int
    TEMP0: List[float]
    TAU: List[float]
    DOFSET: int
    LAST: List[int]
    COMBATH: List[int]
    IRBATH: List[int]

    _order: List[List[str]] = [
        [
            ["ALGORITHM"],
            ["NUM"],
            ["NBATHS"],
            ["TEMP0(1 ... NBATHS)", "TAU(1 ... NBATHS)"],
            ["DOFSET"],
            ["LAST(1 ... DOFSET)", "COMBATH(1 ... DOFSET)", "IRBATH(1 ... DOFSET)"],
        ]
    ]
    # num is not part of the imd file!?

    def __init__(
        self,
        ALGORITHM: int = 0,
        NBATHS: int = 0,
        TEMP0: List[float] = [],
        TAU: List[float] = [],
        DOFSET: int = 0,
        LAST: List[int] = [],
        COMBATH: List[int] = [],
        IRBATH: List[int] = [],
        NUM: int = None,
        content: List[str] = None,
    ):
        if content is None:
            self.ALGORITHM = int(ALGORITHM)
            self.NUM = NUM
            self.NBATHS = int(NBATHS)
            self.TEMP0 = list(TEMP0)
            self.TAU = list(TAU)
            self.DOFSET = int(DOFSET)
            self.LAST = list(LAST)
            self.COMBATH = list(COMBATH)
            self.IRBATH = list(IRBATH)

            if not len(TEMP0) == len(TAU):
                warnings.warn("Warning in MULTIBATH block. There must be the same number of TEMP0 and TAU parameters")
            if not len(TEMP0) == int(NBATHS):
                warnings.warn("Warning in MULTIBATH block. There must be the same number of BATHS and TEMP0 parameters")
                warnings.warn(TEMP0)
                warnings.warn(NBATHS)
            if not len(LAST) == int(NBATHS):
                warnings.warn("Warning in MULTIBATH block. There must be the same number of BATHS and LAST parameters")
            if not len(LAST) == len(COMBATH):
                warnings.warn(
                    "Warning in MULTIBATH block. There must be the same number of COMBATH and LAST parameters."
                )
            if not len(LAST) == len(IRBATH):
                warnings.warn("Warning in MULTIBATH block. There must be the same number of IRBATH and LAST parameters")
            super().__init__(used=True)
        else:
            super().__init__(used=True, content=content)

    def read_content_from_str(self, content: List[str]):
        try:
            setattr(self, self._order[0][0][0], int(content[1].split()[0]))
            if "NUM" in content[0]:
                setattr(self, self._order[0][1][0], int(content[1].split()[1]))

            setattr(self, self._order[0][2][0], int(content[3].split()[0]))
            TEMP0 = []
            TAU = []
            # Read in the temperatures of each temp bath
            for ind in range(5, 5 + self.NBATHS):
                fields = list(map(float, content[ind].split()))
                TEMP0.append(fields[0])
                TAU.append(fields[1])
            setattr(self, "TEMP0", TEMP0)
            setattr(self, "TAU", TAU)
            # Read in DOFSET
            setattr(self, self._order[0][4][0], int(content[5 + self.NBATHS + 1].split()[0]))

            LAST = []
            COMBATH = []
            IRBATH = []
            # Read in the last atoms of each temp bath
            for ind in range(5 + 3 + self.NBATHS, 5 + 3 + 2 * self.NBATHS):
                fields = list(map(float, content[ind].split()))
                LAST.append(fields[0])
                COMBATH.append(fields[1])
                IRBATH.append(fields[2])
            setattr(self, "LAST", LAST)
            setattr(self, "COMBATH", COMBATH)
            setattr(self, "IRBATH", IRBATH)

        except Exception as err:
            raise IOError("Could not parse " + __class__.__name__ + " block!\n" + str(err.args))

    def adapt_multibath(
        self,
        last_atoms_bath: Dict[int, int],
        algorithm: int = None,
        num: int = None,
        T: Union[float, List[float]] = None,
        TAU: float = None,
    ) -> None:
        """adapt_multibath
                This function is adding each atom set into a single multibath.
                #TODO implementation not correct with com_bath and irbath! Works for super simple cases though

        Parameters
        ----------
        last_atoms_bath :   Dict[int,int]
            lastatom:bath
        algorithm : int
            int code for algorihtm
        T : float,List[float], optional
            temperature value
        TAU :   float, optional
            coupling var

        Returns
        -------
        None
        """

        if T is None:
            if any([self.TEMP0[0] != T_test for T_test in self.TEMP0]):
                raise ValueError("Temperatures are not the same, this is not implemented in adapt multibath!")
            T = self.TEMP0[0]
        else:
            self.TEMP0 = [float(T) for x in range(len(self.TEMP0))]

        if TAU is None:
            if any([self.TAU[0] != T_test for T_test in self.TAU]):
                raise ValueError("Temperatures are not the same, this is not implemented in adapt multibath!")

            TAU = self.TAU[0]
        else:
            self.TAU = [float(TAU) for x in range(len(self.TAU))]

        if algorithm is None:
            pass
        else:
            self.ALGORITHM = int(algorithm)

        if num is not None:
            self.NUM = int(num)

        # print("MBATH")
        # print(last_atoms_bath)
        # print("\n")
        # print(last_atoms_bath.values())
        # print("\n")
        # print(set(last_atoms_bath.values()))
        self.NBATHS = len(set(last_atoms_bath.values()))
        self.DOFSET = len(last_atoms_bath)
        self.LAST = [int(last_atom) for last_atom in last_atoms_bath]
        self.COMBATH = [int(last_atoms_bath[last_atom]) for last_atom in last_atoms_bath]
        self.IRBATH = [int(last_atoms_bath[last_atom]) for last_atom in last_atoms_bath]

        if self.NBATHS != len(self.TEMP0):
            self.TEMP0 = [float(self.TEMP0[0]) for x in range(self.NBATHS)]
            self.TAU = [float(self.TAU[0]) for x in range(self.NBATHS)]

    def block_to_string(self) -> str:
        result = ""
        result += str(self.name) + "\n"
        result += "# " + self.field_seperator.join(self._order[0][0]) + "\n"
        result += "  " + str(int(self.ALGORITHM)) + "\n"

        if self.ALGORITHM == "2":
            if self.NUM is None:
                raise Exception("You need to specify the NUM parameter for MULTIBATH if ALGORITHM is 2!")

            result += "# " + self.field_seperator.join(self._order[0][1]) + "\n"
            result += "  " + str(self.NUM) + "\n"

        result += "# " + self.field_seperator.join(self._order[0][2]) + "\n"
        result += "  " + str(int(self.NBATHS)) + "\n"
        result += "# " + self.field_seperator.join(self._order[0][3]) + "\n"
        for index in range(len(self.TEMP0)):
            result += "  " + str(float(self.TEMP0[index])) + self.field_seperator + str(float(self.TAU[index])) + "\n"
        result += "# " + self.field_seperator.join(map(str, self._order[0][4])) + "\n"
        result += "  " + str(int(self.DOFSET)) + "\n"
        result += "# " + self.field_seperator.join(map(str, self._order[0][5])) + "\n"
        for index in range(len(self.LAST)):
            result += (
                "  "
                + str(int(self.LAST[index]))
                + self.field_seperator
                + str(int(self.COMBATH[index]))
                + self.field_seperator
                + str(int(self.IRBATH[index]))
                + "\n"
            )
        result += "END\n"
        return result


class PRESSURESCALE(_generic_imd_block):
    """PRESSURESCALE Block
    This block controls the barostat of the simulation

    Attributes
    -----------
    COUPLE: int
        off(0), calc(1), scale(2)
    SCALE:  int
        off(0), iso(1), aniso(2), full(3), semianiso(4)
    COMP:   float
        compessibility
    TAUP:   float
        coupling strength
    VIRIAL: int
        none(0), atomic(1), group(2)
    SEMIANISOTROPIC:    List[int]
        (semianisotropic couplings: X, Y, Z)
        e.g. 1 1 2: x and y jointly coupled and z separately coupled
        e.g. 0 0 1: constant area (xy-plane) and z coupled to a bath
    PRES0:  List[List[float]]
        reference pressure

    """

    name: str = "PRESSURESCALE"

    COUPLE: int
    SCALE: int
    COMP: float
    TAUP: float
    VIRIAL: int
    SEMIANISOTROPIC: List[int]
    PRES0: List[List[float]]

    _order = [
        [["COUPLE", "SCALE", "COMP", "TAUP", "VIRIAL"], ["SEMIANISOTROPIC COUPLINGS(X, Y, Z)"], ["PRES0(1...3,1...3)"]]
    ]

    def __init__(
        self,
        COUPLE: int = 0,
        SCALE: int = 0,
        COMP: float = 0,
        TAUP: float = 0,
        VIRIAL: int = 0,
        SEMIANISOTROPIC: List[int] = [],
        PRES0: List[List[float]] = [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
        content: List[str] = None,
    ):
        super().__init__(used=True, content=content)
        if content is None:
            self.COUPLE = int(COUPLE)
            self.SCALE = int(SCALE)
            self.COMP = float(COMP)
            self.TAUP = float(TAUP)
            self.VIRIAL = int(VIRIAL)
            self.SEMIANISOTROPIC = list(SEMIANISOTROPIC)
            self.PRES0 = list(PRES0)

    def read_content_from_str(self, content: List[str]):
        try:
            [setattr(self, key, value) for key, value in zip(self._order[0][0], content[1].split())]
            setattr(self, "SEMIANISOTROPIC", list(map(int, content[3].split())))
            setattr(self, "PRES0", [list(map(float, x.split())) for x in content[5:]])
        except Exception as err:
            raise IOError("Could not read string for block " + __class__.__name__ + "\n" + str(err.args))


class FORCE(_generic_imd_block):
    """FORCE Block
        this Block is controlling the forcefield options. It can turn on and of terms, as well as generate force sub groups for analysis.

    Attributes
    -----------
    BONDS:  bool
    ANGLES: bool
    IMPROPER:   bool
    DIHEDRAL:   bool
    ELECTROSTATIC:  bool
    VDW:    bool
    NEGR:   int
        Number of Energy subgroups
    NRE:    List[int]
        List of last atoms for Energy subgroups. (NRE len == NEGR)

    """

    name: str = "FORCE"

    BONDS: bool
    ANGLES: bool
    IMPROPER: bool
    DIHEDRAL: bool
    ELECTROSTATIC: bool
    VDW: bool
    NEGR: int
    NRE: List[int]

    _order = [[["BONDS", "ANGLES", "IMPROPER", "DIHEDRAL", "ELECTROSTATIC", "VDW"], ["NEGR", "NRE"]]]

    def __init__(
        self,
        BONDS: bool = True,
        ANGLES: bool = True,
        IMPROPER: bool = True,
        DIHEDRAL: bool = True,
        ELECTROSTATIC: bool = True,
        VDW: bool = True,
        NEGR: int = 0,
        NRE: List[int] = [],
        content: List[str] = None,
    ):
        """
        Args:
            BONDS:
            ANGLES:
            IMPROPER:
            DIHEDRAL:
            ELECTROSTATIC:
            VDW:
            NEGR:
            NRE (list):
        """
        if content is None:
            self.BONDS = bool(int(BONDS))
            self.ANGLES = bool(int(ANGLES))
            self.IMPROPER = bool(int(IMPROPER))
            self.DIHEDRAL = bool(int(DIHEDRAL))
            self.ELECTROSTATIC = bool(int(ELECTROSTATIC))
            self.VDW = bool(int(VDW))
            try:
                self.NEGR = int(NEGR)
            except Exception:
                if isinstance(NEGR, list):
                    self.NEGR = int(NEGR[0])
                    NRE = NEGR[1:] + NRE
            self.NRE = NRE
            super().__init__(used=True)

        else:
            super().__init__(used=True, content=content)

    def _parse_key_content(self, keyLineNumb: int = 0, contentlines: List[str] = []):
        if keyLineNumb == 0:
            for key, field in zip(self._order[0][keyLineNumb], contentlines[0]):
                # first bring key to attribute form
                key = self._clean_key_from_order(key=key)
                field = self._try_to_convert_field(field=field, key=key)
                setattr(self, key, field)
        else:
            # NEGR
            key = "NEGR"
            field = self._try_to_convert_field(field=contentlines[0][0], key=key)
            setattr(self, key, field)
            # NRE
            key = "NRE"
            field = self._try_to_convert_field(field=contentlines[0][1:], key=key)
            setattr(self, key, field)

    def read_content_from_str(self, content: List[str]):
        try:
            fields = content[1].split(self.field_seperator)
            while "" in fields:
                fields.remove("")
            self._parse_key_content(keyLineNumb=0, contentlines=[fields])
            setattr(self, self._order[0][1][0], int(content[3].split()[0]))
            setattr(self, self._order[0][1][1], list(map(int, content[3].split()[1:])))
        except Exception as err:
            raise IOError("Could not parse FORCE block!\n" + str(err.args))

    def adapt_energy_groups(self, energy_groups: Dict[int, int]):
        """
            Change the Force groups.

        Parameters
        ----------
        residues : Dict[int, int]
            the dictionary contains the forceGroups with the ID of the respective last atom.

           {forceGroupID: lastAtom} = {1:12, 2:26, 3:128}
        """
        self.NEGR = len(energy_groups)
        self.NRE = [int(energy_groups[forceGroupID]) for forceGroupID in sorted(energy_groups)]

    def __adapt_energy_groups(self, residues: Dict[str, Dict[int, int]]):
        """
        Old REEDS option

        adapt_energy_groups
            This method is very "crude" and will put each residue into an own energy group.

        Parameters
        ----------
        residues :  Dict[str, Dict[int, int]]
            you can get this dict from cnf class.

        Returns
        -------
        None

        """

        # set Energy Group ammount
        self.NEGR = len(residues)  # total ammount of Engergy groups
        if "SOLV" in residues and len(residues["SOLV"]) == 0:
            self.NEGR -= 1

        # make residues sortable,
        dict_m = {}
        solvent_names = ["WAT", "SOLV"]
        for x in sorted(residues):
            if type(residues[x]) == dict and x not in dict_m and x not in solvent_names:
                dict_m.update(residues[x])
            elif x in dict_m and not x == "SOLV" and not x == "SOL":
                raise Exception("Found mutliple residues for the same residue id!")

        # build up NRE
        self.NRE = []
        count = 0
        for x in sorted(dict_m.keys()):
            val = dict_m[x]
            count += val
            self.NRE.append(count)

        # add solvents to NRE, but only if there are solvents
        present_solvents = [x for x in residues if (x in solvent_names)]
        if len(present_solvents) > 0:
            if [isinstance(residues[x], dict) and residues[x] != 0 for x in present_solvents]:
                self.NRE.append(sum(residues["SOLV"].values()) + count)
            elif [isinstance(residues[x], int) and residues[x] != 0 for x in present_solvents]:
                self.NRE.append(residues["SOLV"] + count)


class CONSTRAINT(_generic_imd_block):
    """CONSTRAINT block
        This block is controlling constraining the atoms during a simulation.

    Attributes
    ----------
    NTC:    int
    NTCP:   int
    NTCP0:  int
    NTCS:   int
    NTCS0:  int


    """

    name: str = "CONSTRAINT"

    NTC: int
    NTCP: int
    NTCP0: float
    NTCS: int
    NTCS0: float

    _order = [[["NTC"], ["NTCP", "NTCP0(1)"], ["NTCS", "NTCS0(1)"]]]

    def __init__(
        self, NTC: int = 0, NTCP: int = 0, NTCP0: float = 0, NTCS: int = 0, NTCS0: float = 0, content: List[str] = None
    ):
        """
        Args:
            NTC:
            NTCP:
            NTCP0:
            NTCS:
            NTCS0:
        """
        super().__init__(used=True, content=content)
        if content is None:
            self.NTC = int(NTC)
            self.NTCP = int(NTCP)
            self.NTCP0 = float(NTCP0)
            self.NTCS = int(NTCS)
            self.NTCS0 = float(NTCS0)


class PAIRLIST(_generic_imd_block):
    """PAIRLIST Block

            This block is controlling the pairlist control.

    Attributes
    ----------
    ALGORITHM: int
        standard(0) (gromos96 like pairlist)
        grid(1) (md++ grid pairlist)
        grid_cell(2) (creates a mask)
    NSNB:   int
        frequency (number of steps) a pairlist is constructed
    RCUTPL: float
        short-range cut-off in twin-range
    RCUTL:  float
        intermediate-range cut-off in twin-range
    SIZE:   str, float
        grid cell size (or auto = 0.5 * RCUTP)
    TYPE:   str, bool
        chargegoup(0) (chargegroup based cutoff)
        atomic(1)     (atom based cutoff)


    """

    name: str = "PAIRLIST"

    ALGORITHM: int
    NSNB: int
    RCUTPL: float
    RCUTL: float
    SIZE: Union[str, float]
    TYPE: Union[str, bool]

    _order = [[["ALGORITHM", "NSNB", "RCUTP", "RCUTL", "SIZE", "TYPE"]]]

    def __init__(
        self,
        ALGORITHM: int = 0,
        NSNB: int = 0,
        RCUTP: float = 0,
        RCUTL: float = 0,
        SIZE: Union[str, float] = 0,
        TYPE: Union[str, bool] = False,
        content: List[str] = None,
    ):
        """
        Args:
            ALGORITHM:
            NSNB:
            RCUTP:
            RCUTL:
            SIZE:
            TYPE:
        """
        super().__init__(used=True, content=content)
        if content is None:
            self.ALGORITHM = int(ALGORITHM)
            self.NSNB = int(NSNB)
            self.RCUTP = float(RCUTP)
            self.RCUTL = float(RCUTL)
            self.SIZE = SIZE
            self.TYPE = TYPE


class NONBONDED(_generic_imd_block):
    """NONBONDED block

        This block is controlling the Nonbonded term evaluation

    Attributes
    ----------
    NLRELE: int
        1-3 method to handle electrostatic interactions
        -1 : reaction-field (LSERF compatibility mode)
        0 : no electrostatic interactions
        1 : reaction-field
        2 : Ewald method
        3 : P3M method
    APPAK: float
        >= 0.0 reaction-field inverse Debye screening length
    RCRF: float
        >= 0.0 reaction-field radius
    EPSRF: float
         = 0.0 || > 1.0 reaction-field permittivity
    NSLFEXCL: bool
        contribution of excluded atoms to reaction field false=off/true=on
    NSHAPE: float
        -1..10 lattice sum charge-shaping function
        -1 : gaussian
        0..10 : polynomial
    ASHAPE: float
        > 0.0 width of the lattice sum charge-shaping function
    NA2CLC: int
        0..4 controls evaluation of lattice sum A2 term
        0 : A2 = A2~ = 0
        1 : A2~ exact, A2 = A2~
        2 : A2 numerical, A2~ = A2
        3 : A2~ exact from Ewald or from mesh and atom coords, A2 numerical
        4 : A2~ averaged from mesh only, A2 numerical
    TOLA2: float
         > 0.0 tolerance for numerical A2 evaluation
    EPSLS: float
        = 0.0 || > 1.0 lattice sum permittivity (0.0 = tinfoil)
    NKX, NKY, NKZ: float
        > 0 maximum absolute Ewald k-vector components
    KCUT: float
        > 0.0 Ewald k-space cutoff
    NGX, NGY, NGZ: float
        > 0 P3M number of grid points
    NASORD: int
        1..5 order of mesh charge assignment function
    NFDORD: int
        0..5 order of the mesh finite difference operator
        0 : ik - differentiation
        1..5 : finite differentiation
    NALIAS: float
        > 0 number of mesh alias vectors considered
    NSPORD: float
        order of SPME B-spline functions (not available)
    NQEVAL: float
        >= 0 controls accuracy reevaluation
        0 : do not reevaluate
        > 0 : evaluate every NQEVAL steps
    FACCUR: float
        > 0.0 rms force error threshold to recompute influence function
    NRDGRD: bool
        0,1 read influence function
        0 : calculate influence function at simulation start up
        1 : read influence function from file (not yet implemented)
    NWRGRD: bool
        0,1 write influence function
        0 : do not write
        1 : write at the end of the simulation (not yet implemented)
    NLRLJ: bool
        0,1 controls long-range Lennard-Jones corrections
        0 : no corrections
        1 : do corrections (not yet implemented)
    SLVDNS: float
        > 0.0 average solvent density for long-range LJ correction (ignored)
    """

    name: str = "NONBONDED"

    NLRELE: int
    APPAK: float
    RCRF: float
    EPSRF: float
    NSLFEXCL: bool
    NSHAPE: int
    ASHAPE: float
    NA2CLC: int
    TOLA2: str
    EPSLS: float
    NKX: int
    NKY: int
    NKZ: int
    KCUT: float
    NGX: int
    NGY: int
    NGZ: int
    NASORD: int
    NFDORD: int
    NALIAS: int
    NSPORD: int
    NQEVAL: int
    FACCUR: float
    NRDGRD: bool
    NWRGRD: bool
    NLRLJ: bool
    SLVDNS: float

    _order = [
        [
            ["NLRELE"],
            ["APPAK", "RCRF", "EPSRF", "NSLFEXCL"],
            ["NSHAPE", "ASHAPE", "NA2CLC", "TOLA2", "EPSLS"],
            ["NKX", "NKY", "NKZ", "KCUT"],
            ["NGX", "NGY", "NGZ", "NASORD", "NFDORD", "NALIAS", "NSPORD"],
            ["NQEVAL", "FACCUR", "NRDGRD", "NWRGRD"],
            ["NLRLJ", "SLVDNS"],
        ]
    ]

    def __init__(
        self,
        NLRELE: int = 0,
        APPAK: float = 0,
        RCRF: float = 0,
        EPSRF: float = 0,
        NSLFEXCL: bool = False,
        NSHAPE: int = 0,
        ASHAPE: float = 0,
        NA2CLC: int = 0,
        TOLA2: str = 0,
        EPSLS: float = 0,
        NKX: int = 0,
        NKY: int = 0,
        NKZ: int = 0,
        KCUT: float = 0,
        NGX: int = 0,
        NGY: int = 0,
        NGZ: int = 0,
        NASORD: int = 0,
        NFDORD: int = 0,
        NALIAS: int = 0,
        NSPORD: int = 0,
        NQEVAL: int = 0,
        FACCUR: float = 0,
        NRDGRD: bool = False,
        NWRGRD: bool = False,
        NLRLJ: bool = False,
        SLVDNS: float = 0,
        content: List[str] = None,
    ):
        super().__init__(used=True, content=content)

        if content is None:
            self.NLRELE = NLRELE
            self.APPAK = APPAK
            self.RCRF = RCRF
            self.EPSRF = EPSRF
            self.NSLFEXCL = NSLFEXCL
            self.NSHAPE = int(NSHAPE)
            self.ASHAPE = ASHAPE
            self.NA2CLC = int(NA2CLC)
            self.TOLA2 = str(TOLA2)
            self.EPSLS = EPSLS
            self.NKX = int(NKX)
            self.NKY = int(NKY)
            self.NKZ = int(NKZ)
            self.KCUT = KCUT
            self.NGX = int(NGX)
            self.NGY = int(NGY)
            self.NGZ = int(NGZ)
            self.NASORD = int(NASORD)
            self.NFDORD = int(NFDORD)
            self.NALIAS = int(NALIAS)
            self.NSPORD = int(NSPORD)
            self.NQEVAL = int(NQEVAL)
            self.FACCUR = FACCUR
            self.NRDGRD = bool(int(NRDGRD))
            self.NWRGRD = bool(int(NWRGRD))
            self.NLRLJ = bool(int(NLRLJ))
            self.SLVDNS = SLVDNS


class INITIALISE(_generic_imd_block):
    """INITIALISE block

        This block controls the Initialisation of a simulation.

    Attributes
    ----------
    NTIVEL: bool
        0,1 controls generation of initial velocities
        0: read from configuration (default)
        1: generate from Maxell distribution at temperature TEMPI
    NTISHK: int
        0..3 controls shaking of initial configuration
        0: no intial SHAKE (default)
        1: initial SHAKE on coordinates only
        2: initial SHAKE on velocities only
        3: initial SHAKE on coordinates and velocities
    NTINHT: bool
        0,1 controls generation of initial Nose-Hoover chain variables
        0: read from configuration (default)
        1: reset variables to zero.
    NTINHB: bool
        0,1 controls generation of initial Nose-Hoover (chain) barostat variables
        0: read from strartup file (if applicable) (default)
        1: reset variables to zero
    NTISHI: bool
        0,1 controls initial setting for lattice shift vectors
        0: read from configuration (default)
        1: reset shifts to zero.
    NTIRTC: bool
        0,1 controls initial setting of positions and orientations for roto-translational constraints
        0: read from configuration (default)
        1: reset based on initial configuraion of startup file
    NTICOM: int
        0,1,2 controls initial removal of COM motion
        0: no initial system COM motion removal (default)
        1: initial COM translation is removed
        2: initial COM rotation is removed
    NTISTI: bool
        0,1 controls generation of stochastic integrals
        0: read stochastic integrals and IG from configuration (default)
        1: set stochastic integrals to zero and use IG from here.
    IG: int
        random number generator seed
    TEMPI:  float
        initial temperature
    """

    name = "INITIALISE"

    NTIVEL: bool
    NTISHK: int
    NTINHT: bool
    NTINHB: bool
    NTISHI: bool
    NTIRTC: bool
    NTICOM: int
    NTISTI: bool
    IG: int
    TEMPI: float

    _order = [[["NTIVEL", "NTISHK", "NTINHT", "NTINHB"], ["NTISHI", "NTIRTC", "NTICOM"], ["NTISTI"], ["IG", "TEMPI"]]]

    def __init__(
        self,
        NTIVEL: bool = False,
        NTISHK: int = 0,
        NTINHT: bool = False,
        NTINHB: bool = False,
        NTISHI: bool = False,
        NTIRTC: bool = False,
        NTICOM: int = 0,
        NTISTI: bool = False,
        IG: int = 0,
        TEMPI: float = 0,
        content: List[str] = None,
    ):
        """
        Args:
            NTIVEL:
            NTISHK:
            NTINHT:
            NTINHB:
            NTISHI:
            NTIRTC:
            NTICOM:
            NTISTI:
            IG:
            TEMPI:
        """
        super().__init__(used=True, content=content)
        if content is None:
            self.NTIVEL = NTIVEL
            self.NTISHK = NTISHK
            self.NTINHT = NTINHT
            self.NTINHB = NTINHB
            self.NTISHI = NTISHI
            self.NTIRTC = NTIRTC
            self.NTICOM = NTICOM
            self.NTISTI = NTISTI
            self.IG = IG
            self.TEMPI = TEMPI


class COMTRANSROT(_generic_imd_block):
    """COMTRANSROT block

    This block controls the center of mass translation and rotation removal. (flying ice cube problem)

    Attributes
    ----------

    NSCM : int
        controls system centre-of-mass (com) motion removal
        0: no com motion removal (default)
        < 0: com translation and rotation are removed every abs(NSCM) steps.
        > 0: com translation is removed every NSCM steps.
    """

    name: str = "COMTRANSROT"

    NSCM: int

    _order = [[["NSCM"]]]

    def __init__(self, NSCM: int = 0, content: List[str] = None):
        super().__init__(used=True, content=content)
        if content is None:
            self.NSCM = NSCM


class QMMM(_generic_imd_block):
    """QMMM Block

        The QMMM block describes how to treat an embedded QM zone

    Attributes
    ----------
    NTQMMM:    int
        Define if QM/MM is applied (1) or not (0, default)
    NTQMSW:    int
        QM software package to use. Available options are MNDO (0, default), Turbomole (1), DFTB (2), MOPAC (3), ORCA (4), and XTB (5)
    RCUTQ:     float
        Cutoff for inclusion of MM charge groups. If the value is 0, all particles are considered.
    NTWQMMM:   int
        Write QM/MM related data to special trajectory every NTWQMMM step. Switch off with 0.
    QMLJ:      int
        Define if LJ interaction is activated in the QM zone (1) or not (0).
    MMSCAL:    float
        Scale MM charges with (2/pi)*atan(x*(r_{mm}-r_{mm})). Values > 0.0 describe the scaling factor, values < 0.0 turn off scaling (default).

    """

    name: str = "QMMM"

    # fields
    NTQMMM: int  # Define if QM/MM is applied (1) or not (0, default)
    NTQMSW: int  # QM software package to use. Available options are MNDO (0, default), Turbomole (1), DFTB (2), MOPAC (3), ORCA (4), and XTB (5)
    RCUTQ: float  # Cutoff for inclusion of MM charge groups. If the value is 0, all particles are considered.
    NTWQMMM: int  # Write QM/MM related data to special trajectory every NTWQMMM step. Switch off with 0.
    QMLJ: int  # Define if LJ interaction is activated in the QM zone (1) or not (0).
    MMSCAL: float  # Scale MM charges with (2/pi)*atan(x*(r_{mm}-r_{mm})). Values > 0.0 describe the scaling factor, values < 0.0 turn off scaling (default).

    _order = [[["NTQMMM", "NTQMSW", "RCUTQ", "NTWQMMM", "QMLJ", "MMSCAL"]]]

    def __init__(
        self,
        NTQMMM: int = 0,
        NTQMSW: int = 0,
        RCUTQ: float = 0.0,
        NTWQMMM: int = 0,
        QMLJ: int = 0,
        MMSCAL: float = -1.0,
        content: List[str] = None,
    ):
        super().__init__(used=True, content=content)
        if content is None:
            self.NTQMMM = int(NTQMMM)
            self.NTQMSW = int(NTQMSW)
            self.RCUTQ = float(RCUTQ)
            self.NTWQMMM = int(NTWQMMM)
            self.QMLJ = int(QMLJ)
            self.MMSCAL = float(MMSCAL)


class EDS(_generic_imd_block):
    """EDS block

    This block is used in an EDS simulation.

    Attributes
    -----------
    NUMSTATES:  int
        EDS-States
    S:  float
        smoothness parameter
    EIR:    List[float]
        energy offsets
    EDS:    bool, optional
        turn on EDS_simulation
    ALPHLJ: float, optional
    ALPHC:  float, optional
    FUNCTIONAL: int, optional
        1: Single s Hamiltonian (default)
        2: Hamiltonian with NUMSTATES*(NUMSTATES-1)/2 (pairwise) S parameters ==> changes type of S
        3: Hamiltonian with (NUMSTATES-1) S parameters ==> changes type of S
    """

    _order = [[["EDS"], ["ALPHLJ", "ALPHC"], ["FUNCTIONAL FORM"], ["NUMSTATES"], ["S"], ["EIR"]]]

    def __init__(
        self,
        NUMSTATES: int = 0,
        S: float = 0,
        EIR: List[float] = [],
        EDS: bool = 1,
        ALPHLJ: float = 0.0,
        ALPHC: float = 0.0,
        FUNCTIONAL: int = 1,
        content: List[str] = None,
    ):
        super().__init__(used=True, content=content)
        if content is None:
            self.name = "EDS"
            self.EDS = EDS
            self.ALPHLJ = ALPHLJ
            self.ALPHC = ALPHC
            self.FUNCTIONAL = FUNCTIONAL
            self.NUMSTATES = NUMSTATES
            self.S = S
            self.EIR = EIR

    def read_content_from_str(self, content: List[str]):
        try:
            setattr(self, "EDS", int(content[1].split()[0]))
            setattr(self, "ALPHLJ", float(content[3].split()[0]))
            setattr(self, "ALPHC", float(content[3].split()[1]))
            setattr(self, "FUNCTIONAL", int(content[5].split()[0]))
            setattr(self, "NUMSTATES", int(content[7].split()[0]))
            setattr(self, "S", list(map(float, content[9].split())))
            setattr(self, "EIR", list(map(float, content[11].split())))

        except Exception as err:
            raise IOError("Could not parse block from str - " + __class__.__name__ + "\n" + str(err.args))


class DISTANCERES(_generic_imd_block):
    """DISTANCERES Block"""

    name = "DISTANCERES"

    NTDIR: int
    NTDIRA: int
    CDIR: int
    DIR0: int
    TAUDIR: int
    FORCESCALE: int
    VDIR: int
    NTWDIR: int

    _order = [[["NTDIR", "NTDIRA", "CDIR", "DIR0", "TAUDIR", "FORCESCALE", "VDIR", "NTWDIR"]]]

    def __init__(
        self,
        NTDIR: int = 0,
        NTDIRA: int = 0,
        CDIR: int = 0,
        DIR0: int = 0,
        TAUDIR: int = 1,
        FORCESCALE: int = 0,
        VDIR: int = 0,
        NTWDIR: int = 0,
        content: List[str] = None,
    ):
        """
        Args:
            NTDIR:
            NTDIRA:
            CDIR:
            DIR0:
            TAUDIR:
            FORCESCALE:
            VDIR:
            NTWDIR:
        """
        super().__init__(used=True, content=content)
        if content is None:
            self.NTDIR = NTDIR
            self.NTDIRA = NTDIRA
            self.CDIR = CDIR
            self.DIR0 = DIR0
            self.TAUDIR = TAUDIR
            self.FORCESCALE = FORCESCALE
            self.VDIR = VDIR
            self.NTWDIR = NTWDIR


class POSITIONRES(_generic_imd_block):
    """POSITIONRES block

        This block allows the managment of the Position Restraints during the Simulation.

    Attributes
    ----------
    NTPOR:
        0..3 controls atom positions re(con)straining.
        0: no position re(con)straints (default)
        1: restraining with force constant CPOR
        2: restraining with force constant CPOR weighted by atomic B-factors
        3: constraining
    NTPORB:
        0,1 controls reading of reference positions and B-factors
        0: read reference positions from startup file.
        1: read reference positions and B-factors from special file
    NTPORS:
        0,1 controls scaling of reference positions upon pressure scaling
        0: do not scale reference positions
        1: scale reference positions
    CPOR:
        >= 0 position restraining force constant
    """

    name = "POSITIONRES"

    NTPOR: int
    NTPORB: int
    NTPORS: int
    CPOR: int

    _order = [[["NTPOR", "NTPORB", "NTPORS", "CPOR"]]]

    def __init__(self, NTPOR: int = 0, NTPORB: int = 0, NTPORS: int = 0, CPOR: int = 0, content: List[str] = None):
        super().__init__(used=True, content=content)
        if content is None:
            self.NTPOR = NTPOR
            self.NTPORB = NTPORB
            self.NTPORS = NTPORS
            self.CPOR = CPOR


class PRINTOUT(_generic_imd_block):
    """PRINTOUT block

        This Block manages the output frequency into the .omd/std-out file.

    Attributes
    ----------
    NTPR: int
        print out energies, etc. every NTPR steps
    NTPP: int
        =1 perform dihedral angle transition monitoring

    """

    name: str = "PRINTOUT"

    NTPR: int
    NTPP: int

    _order = [[["NTPR", "NTPP"]]]

    def __init__(self, NTPR: int = 0, NTPP: int = 0, content: List[str] = None):
        super().__init__(used=True, content=content)
        if content is None:
            self.name = "PRINTOUT"
            self.NTPR = NTPR
            self.NTPP = NTPP


class ENERGYMIN(_generic_imd_block):
    """ENERGYMIN block
    This block takes care of managing the Energyminimization controls.

    Attributes
    ----------
    NTEM:   int
        controls energy minimisation mode.
        0: do not do energy minimisation (default)
        1: steepest-descent minimisation
        2: Fletcher-Reeves conjugate-gradient minimisation
        3: Polak-Ribiere conjugate-gradient minimisation
    NCYC:   int
        >0 number of steps before resetting the conjugate-gradient search direction
        =0 reset only if the energy grows in the search direction
    DELE:   float
        >0.0 energy threshold for convergence
        >0.0 (conjugate-gradient) RMS force threshold for convergence
    DX0:    float
        >0.0 initial step size
    DXM:    float
        >0.0 maximum step size
    NMIN:   float
        >0 minimum number of minimisation steps
    FLIM:   float
        >=0.0 limit force to maximum value (FLIM > 0.0 is not recommended)
    """

    name = "ENERGYMIN"

    NTEM: int
    NCYC: int
    DELE: float
    DX0: float
    DXM: float
    NMIN: int
    FLIM: float

    _order = [["NTEM    NCYC    DELE    DX0     DXM   NMIN   FLIM".split()]]

    def __init__(
        self,
        NTEM: int = 0,
        NCYC: int = 0,
        DELE: float = 0,
        DX0: float = 0,
        DXM: float = 0,
        NMIN: int = 0,
        FLIM: float = 0,
        content: List[str] = None,
    ):
        super().__init__(used=True, content=content)
        if content is None:
            self.NTEM = NTEM
            self.NCYC = NCYC
            self.DELE = DELE
            self.DX0 = DX0
            self.DXM = DXM
            self.NMIN = NMIN
            self.FLIM = FLIM


class WRITETRAJ(_generic_imd_block):
    """WRITETRAJ

        The WRITETRAJ Block manages the output frequency of GROMOS. Here you can determine how often a file should be written.

    Attributes
    ----------
    NTWX:   int
        controls writing of coordinate trajectory
        0: no coordinate trajectory is written (default)
        >0: write solute and solvent coordinates every NTWX steps
        <0: write solute coordinates every |NTWX| steps
    NTWSE:   int
        >= 0 selection criteria for coordinate trajectory writing
        0: write normal coordinate trajectory
        >0: write minimum-energy coordinate and energy trajectory (based on the
        energy entry selected by NTWSE and as blocks of length NTWX)
        (see configuration/energy.cc or ene_ana library for indices)
    NTWV:   int
        controls writing of velocity trajectory
        0: no velocity trajectory is written (default)
        >0: write solute and solvent velocities every NTWV steps
        <0: write solute velocities every |NTWV| steps
    NTWF:   int
        controls writing of force trajectory
        0: no force trajectory is written (default)
        >0: write solute and solvent forces every NTWF steps
        <0: write solute forces every |NTWF| steps
    NTWE:   int
        >= 0 controls writing of energy trajectory
        0: no energy trajectory is written (default)
        >0: write energy trajectory every NTWE steps
    NTWG:   int
        >= 0 controls writing of free energy trajectory
        0: no free energy trajectory is written (default)
        >0: write free energy trajectory every NTWG steps
    NTWB:   int
        >= 0 controls writing of block-averaged energy trajectory
        0: no block averaged energy trajectory is written (default)
        >0: write block-averaged energy variables every |NTWB| steps
            (and free energies if NTWG > 0) trajectory
    """

    name = "WRITETRAJ"

    NTWX: int
    NTWSE: int
    NTWV: int
    NTWF: int
    NTWE: int
    NTWG: int
    NTWB: int

    _order = [[["NTWX", "NTWSE", "NTWV", "NTWF", "NTWE", "NTWG", "NTWB"]]]

    def __init__(
        self,
        NTWX: int = 0,
        NTWSE: int = 0,
        NTWV: int = 0,
        NTWF: int = 0,
        NTWE: int = 0,
        NTWG: int = 0,
        NTWB: int = 0,
        content: List[str] = None,
    ):
        """
        Args:
            NTWX:
            NTWSE:
            NTWV:
            NTWF:
            NTWE:
            NTWG:
            NTWB:
        """
        super().__init__(used=True, content=content)
        if content is None:
            self.NTWX = int(NTWX)
            self.NTWSE = int(NTWSE)
            self.NTWV = int(NTWV)
            self.NTWF = int(NTWF)
            self.NTWE = int(NTWE)
            self.NTWG = int(NTWG)
            self.NTWB = int(NTWB)


class AMBER(_generic_imd_block):
    """AMBER block

    Attributes
    ----------
    Amber:  bool
    AMBSCAL:    float



    """

    name: str = "AMBER"

    Amber: bool
    AMBSCAL: float

    _order = [[["AMBER", "AMBSCAL"]]]

    def __init__(self, AMBER: int = 0, AMBSCAL: int = 0, content: List[str] = None):
        super().__init__(used=True, content=content)
        if content is None:
            self.name = "AMBER"
            self.AMBER = AMBER
            self.AMBSCAL = AMBSCAL


class COVALENTFORM(_generic_imd_block):
    """COVALENTFORM Block

        The COVALENTFORM Block manages the functional form of the Bonded contributions to the force field.
        It is optional in an imd file to define the block.

    Attributes
    ----------
    NTBBH:  int
        0,1 controls bond-stretching potential
        0: quartic form (default)
        1: harmonic form
    NTBAH:  int
        0,1 controls bond-angle bending potential
        0: cosine-harmonic (default)
        1: harmonic
    NTBDN:  int
        0,1 controls torsional dihedral potential
        0: arbitrary phase shifts (default)
        1: phase shifts limited to 0 and 180 degrees.
    """

    name = "COVALENTFORM"

    NTBBH: int
    NTBAH: int
    NTBDN: int

    _order = [[["NTBBH", "NTBAH", "NTBDN"]]]

    def __init__(self, NTBBH: int = 0, NTBAH: int = 0, NTBDN: int = 0, content: List[str] = None):
        super().__init__(used=True, content=content)
        if content is None:
            self.NTBBH = NTBBH
            self.NTBAH = NTBAH
            self.NTBDN = NTBDN


class INNERLOOP(_generic_imd_block):
    """INNERLOOP block

    Attributes
    ----------
    NTILM: int
    NTILS: int
    NGPUS: int
    """

    name = "INNERLOOP"

    NTILM: int
    NTILS: int
    NGPUS: int

    _order = [[["NTILM", "NTILS", "NGPUS"]]]

    def __init__(self, NTILM: int = 0, NTILS: int = 0, NGPUS: int = 0, content: List[str] = None):
        super().__init__(used=True, content=content)
        if content is None:
            self.NTILM = NTILM
            self.NTILS = NTILS
            self.NGPUS = NGPUS


class LAMBDA(_generic_imd_block):
    """
    LAMBDA block
        This block is controlling the perturbation for free energy calculations.

    Attributes
    ----------
    NTIL:   int
    NTLI(1):   int or List[int]
    NILG1(1):  int or List[int]
    NILG2(1):  int or List[int]
    ALI(1):    float or List[float]
    BLI(1):    float or List[float]
    CLI(1):    float or List[float]
    DLI(1):    float or List[float]
    ELI(1):    float or List[float]

    """

    name: str = "LAMBDA"

    NTIL: int
    NTLI: List[int]
    NILG1: List[int]
    NILG2: List[int]
    ALI: List[float]
    BLI: List[float]
    CLI: List[float]
    DLI: List[float]
    ELI: List[float]

    _order: List[List[str]] = [
        [
            ["NTIL"],
            ["NTLI(1..)", "NILG1(1..)", "NILG2(1..)", "ALI(1..)", "BLI(1;;)", "CLI(1;;)", "DLI(1;;)", "ELI(1;;)"],
        ]
    ]

    def __init__(
        self,
        NTIL: int = 0,
        NTLI: List[int] = [],
        NILG1: List[int] = [],
        NILG2: List[int] = [],
        ALI: List[float] = [],
        BLI: List[float] = [],
        CLI: List[float] = [],
        DLI: List[float] = [],
        ELI: List[float] = [],
        content: List[str] = None,
    ):
        super().__init__(used=True, content=content)
        if content is None:
            self.NTIL = NTIL
            self.NTLI = NTLI
            self.NILG1 = NILG1
            self.NILG2 = NILG2
            self.ALI = ALI
            self.BLI = BLI
            self.CLI = CLI
            self.DLI = DLI
            self.ELI = ELI

            if not len(NTLI) == len(NILG1):
                warnings.warn("Warning in LAMBDA Block. There must be the same number of NTLI and NILG1 parameters.")
            if not len(NTLI) == len(NILG2):
                warnings.warn("Warning in LAMBDA Block. There must be the same number of NTLI and NILG2 parameters.")
            if not len(NTLI) == len(ALI):
                warnings.warn("Warning in LAMBDA Block. There must be the same number of NTLI and ALI parameters.")
            if not len(NTLI) == len(BLI):
                warnings.warn("Warning in LAMBDA Block. There must be the same number of NTLI and BLI parameters.")
            if not len(NTLI) == len(CLI):
                warnings.warn("Warning in LAMBDA Block. There must be the same number of NTLI and CLI parameters.")
            if not len(NTLI) == len(DLI):
                warnings.warn("Warning in LAMBDA Block. There must be the same number of NTLI and DLI parameters.")
            if not len(NTLI) == len(ELI):
                warnings.warn("Warning in LAMBDA Block. There must be the same number of NTLI and ELI parameters.")


class ROTTRANS(_generic_imd_block):
    """Roto-translational block

        This block is for roto-translational constraints.
        Note: use either centre of mass removal or roto-translational constraints but not both!

    Attributes
    ----------
    RTC: int
        turn roto-translational constraints on (1)
    RTCLAST: int
        last atom of subset to be roto-translationally constrained
    """

    name: str = "ROTTRANS"

    RTC: int
    RTCLAST: int

    _order = [[["RTC", "RTCLAST"]]]

    def __init__(self, RTC: int = 0, RTCLAST: int = 0, content: List[str] = None):
        super().__init__(used=True, content=content)
        if content is None:
            self.RTC = RTC
            self.RTCLAST = RTCLAST


class RANDOMNUMBERS(_generic_imd_block):
    """
    Random Numbers Block
    """

    name: str = "RANDOMNUMBERS"

    NTRNG: int
    NTGSL: int

    _order = [[["NTRNG"], ["NTGSL"]]]

    def __init__(self, NTRNG: int = 0, NTGSL: int = 0, content: List[str] = None):
        super().__init__(used=True, content=content)
        if content is None:
            self.NTRNG = int(NTRNG)
            self.NTGSL = int(NTGSL)
