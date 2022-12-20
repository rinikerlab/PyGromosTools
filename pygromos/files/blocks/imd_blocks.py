import warnings
from numbers import Number
from typing import List, Dict

from pygromos.files.blocks._general_blocks import _generic_gromos_block
from pygromos.files.blocks._general_blocks import TITLE, TIMESTEP

# forward declarations
TITLE:TITLE = TITLE

class _generic_imd_block(_generic_gromos_block):

    name = "genericBlock"
    _order: List[List[str]]  # contains the ordering of all fields in a block

    def __init__(self, used: bool):
        super().__init__(name=self.name, used=used)

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
                        if isinstance(attribute, (str, Number)):  # One element field
                            if (isinstance(attribute, float)):  # supress scientific notation for floats!
                                attribute = format(attribute, "f")
                            result += str(attribute) + self.field_seperator
                        elif isinstance(attribute, List) and all([isinstance(x, (str, Number)) for x in attribute]):  # list content
                            if (all([isinstance(x, str) for x in attribute])):
                                result += self.field_seperator.join(attribute) + self.field_seperator
                            elif (all([isinstance(x, Number) for x in attribute])):
                                if (isinstance(attribute[0], float)):
                                    result += self.field_seperator.join(map(lambda x: format(x, "f"), attribute)) + self.field_seperator
                                else:
                                    result += self.field_seperator.join(map(str, attribute)) + self.field_seperator
                            else:
                                raise ValueError("could not Interpret list:  " + str(element) + "\t\n" + str(attribute) + "\nEOF\n")
                        elif isinstance(attribute, List) and all([isinstance(x, List) for x in attribute]):  # matrices
                            pre_parsed_rows = map(lambda x: self.field_seperator.join(map(str, x)), attribute)
                            result += (self.line_seperator + self.field_seperator).join(pre_parsed_rows) + self.field_seperator
                        else:
                            result += (self.line_seperator + self.field_seperator).join(
                                map(self.field_seperator.join, attribute)) + self.field_seperator
                    except:
                        raise ValueError(
                            "Could not convert attribute " + str(element) + " to string!\n value of attribute was: \n" + str(attribute) + "\nEOF\n")

                result += self.line_seperator
        result += "END\n"
        return result


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
    NPM: int = 0  # number of Solute Molecules
    NSM: int = 0  # number of Solvent Molecules

    _order = [[["NPM", "NSM"]]]

    def __init__(self, NPM, NSM):
        super().__init__(used=True)
        self.NPM = NPM
        self.NSM = NSM


class STEP(_generic_imd_block):
    """ Step Block

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

    def __init__(self, NSTLIM, T, DT):
        super().__init__(used=True)
        self.NSTLIM = NSTLIM
        self.T = T
        self.DT = DT

class REPLICA_EDS(_generic_imd_block):
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
    name: str = "REPLICA_EDS"

    REEDS: bool

    NRES: int
    NEOFF:int
    NUMSTATES: int

    RES: List[float]
    EIR: List[float]

    NRETRIAL: int
    NREQUIL: int
    EDS_STAT_OUT: int
    CONT: bool
    PERIODIC:int

    _order = [[["REEDS"], ["NRES", "NUMSTATES", "NEOFF"], ["RES(1 ... NRES)"],
               ["EIR(NUMSTATES x NRES)"], ["NRETRIAL", "NREQUIL", "CONT", "EDS_STAT_OUT", "PERIODIC"]]]

    def __init__(self, REEDS: bool, NRES: int, NUMSTATES: int, NEOFF:int, RES: List[float], EIR: List[List[float]], NRETRIAL: int, NREQUIL: int,
                 EDS_STAT_OUT: int, CONT: bool, PERIODIC:int):
        super().__init__(used=True)
        self.REEDS = REEDS

        self.NRES = NRES
        self.NEOFF = NEOFF
        self.NUMSTATES = NUMSTATES

        self.RES = RES
        self.EIR = EIR

        self.NRETRIAL = NRETRIAL
        self.NREQUIL = NREQUIL
        self.CONT = CONT
        self.EDS_STAT_OUT = EDS_STAT_OUT
        self.PERIODIC = PERIODIC

class BOUNDCOND(_generic_imd_block):
    """Boundary Condition Block

        This block controles

    Attributes
    ----------
    NTB:    int
    NDFMIN: int

    """

    name: str = "BOUNDCOND"

    NTB: int
    NDFMIN: int

    _order = [[["NTB", "NDFMIN"]]]

    def __init__(self, NTB: int, NDFMIN: int):
        super().__init__(used=True)
        self.NTB = NTB
        self.NDFMIN = NDFMIN

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

    NTSD:    int
    NTFR:    int
    NSFR:    int
    NBREF:   int
    RCUTF:   float
    CFRIC:   float
    TEMPSD:  float

    _order = [[["NTSD", "NTFR", "NSFR", "NBREF", "RCUTF", "CFRIC", "TEMPSD"]]]

    def __init__(self, 
                NTSD:    int,
                NTFR:    int,
                NSFR:    int,
                NBREF:   int,
                RCUTF:   float,
                CFRIC:   float,
                TEMPSD:  float):
        super().__init__(used=True)
        self.NTSD = NTSD
        self.NTFR = NTFR
        self.NSFR = NSFR
        self.NBREF = NBREF
        self.RCUTF = RCUTF
        self.CFRIC = CFRIC
        self.TEMPSD = TEMPSD

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

    RTC:     int
    RTCLAST: int

    _order = [[["RTC", "RTCLAST"]]]

    def __init__(self, 
                RTC:     int,
                RTCLAST: int):
        super().__init__(used=True)
        self.RTC = RTC
        self.RTCLAST = RTCLAST

class PERTURBATION(_generic_imd_block):
    """perturbation Block

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

    _order = [[["NTG", "NRDGL", "RLAM", "DLAMT"],
               ["ALPHLJ", "ALPHC", "NLAM", "NSCALE"]]]

    def __init__(self, NTG: int, NRDGL: int, RLAM:float, DLAMT:float, ALPHLJ:float, ALPHC:float, NLAM:int, NSCALE:int):
        super().__init__(used=True)
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

    def __init__(self, NRLAM: int,MINLAM:float, MAXLAM:float):
        super().__init__(used=True)
        self.NRLAM = NRLAM
        self.MINLAM = MINLAM
        self.MAXLAM = MAXLAM

class MULTIBATH(_generic_imd_block):
    """MULTIBATH Block

        This block controls the Thermostat of a simulation. Multiple temperature baths can be coupled.

    Attributes
    ----------
    ALGORITHM:  int
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

    _order: List[List[str]] = [[["ALGORITHM"], ["NUM"], ["NBATHS"], ["TEMP0(1 ... NBATHS)", "TAU(1 ... NBATHS)"],
                                ["DOFSET"], ["LAST(1 ... DOFSET)", "COMBATH(1 ... DOFSET)", "IRBATH(1 ... DOFSET)"]]]

    def __init__(self, ALGORITHM: int, NBATHS: int, TEMP0: List[float], TAU: List[float], DOFSET: int, LAST: List[int], COMBATH: List[int],
                 IRBATH: List[int], NUM: int = None):

        super().__init__(used=True)
        self.ALGORITHM = ALGORITHM
        self.NUM = NUM
        self.NBATHS = NBATHS
        self.TEMP0 = TEMP0
        self.TAU = TAU
        self.DOFSET = DOFSET
        self.LAST = LAST
        self.COMBATH = COMBATH
        self.IRBATH = IRBATH
        
        if not len(TEMP0) == len(TAU):
            warnings.warn("Warning in MULTIBATH block. There must be the same number of TEMP0 and TAU parameters")
        if not len(TEMP0) == int(NBATHS):
            warnings.warn("Warning in MULTIBATH block. There must be the same number of BATHS and TEMP0 parameters")
            warnings.warn(TEMP0)
            warnings.warn(NBATHS)
        if not len(LAST) == int(NBATHS):
            warnings.warn("Warning in MULTIBATH block. There must be the same number of BATHS and LAST parameters")
        if not len(LAST) == len(COMBATH):
            warnings.warn("Warning in MULTIBATH block. There must be the same number of COMBATH and LAST parameters.")
        if not len(LAST) == len(IRBATH):
            warnings.warn("Warning in MULTIBATH block. There must be the same number of IRBATH and LAST parameters")

    def adapt_multibath(self, last_atoms_bath: Dict[int, int], algorithm: int = None, num: int = None, T: (float, List[float]) = None, TAU: float = None) -> None:
        """ adapt_multibath
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

        if (T == None):
            if (any([self.TEMP0[0] != T_test for T_test in self.TEMP0])):
                raise ValueError("Temperatures are not the same, this is not implemented in adapt multibath!")
            T = self.TEMP0[0]
        else:
            self.TEMP0 = [T for x in range(len(self.TEMP0))]

        if (TAU == None):
            if (any([self.TAU[0] != T_test for T_test in self.TAU])):
                raise ValueError("Temperatures are not the same, this is not implemented in adapt multibath!")

            TAU = self.TAU[0]
        else:
            self.TAU = [TAU for x in range(len(self.TAU))]

        if (algorithm == None):
            pass
        else:
            self.ALGORITHM = algorithm

        if (num != None):
            self.NUM = num

        # TODO implementation not correct with com_bath and irbath! Works for super simple cases though
        print(last_atoms_bath)
        print("\n")
        print(last_atoms_bath.values())
        print("\n")
        print(set(last_atoms_bath.values()))
        self.NBATHS = len(set(last_atoms_bath.values()))
        self.DOFSET = len(last_atoms_bath)
        self.LAST = [last_atom for last_atom in last_atoms_bath]
        self.COMBATH = [last_atoms_bath[last_atom] for last_atom in last_atoms_bath]
        self.IRBATH = [last_atoms_bath[last_atom] for last_atom in last_atoms_bath]

        if (self.NBATHS > len(self.TEMP0)):
            self.TEMP0 = [self.TEMP0[0] for x in range(self.NBATHS)]
            self.TAU = [self.TAU[0] for x in range(self.NBATHS)]

    def block_to_string(self) -> str:
        result = ""
        result += str(self.name) + "\n"
        result += "# " + "\t".join(self._order[0][0]) + "\n"
        result += "  " + str(self.ALGORITHM) + "\n"

        if(self.ALGORITHM == "2"):
            if(self.NUM is None):
                raise Exception("You need to specify the NUM parameter for MULTIBATH if ALGORITHM is 2!")
               
            result += "# " + "\t".join(self._order[0][1]) + "\n"
            result += "  " + str(self.NUM) + "\n"

        result += "# " + "\t".join(self._order[0][2]) + "\n"
        result += "  " + str(self.NBATHS) + "\n"
        result += "# " + "\t".join(self._order[0][3]) + "\n"
        for index in range(len(self.TEMP0)):
            result += "  " + str(self.TEMP0[index]) + "\t" + str(self.TAU[index]) + "\n"
        result += "# " + "\t".join(map(str, self._order[0][4])) + "\n"
        result += "  " + str(self.DOFSET) + "\n"
        result += "# " + "\t".join(map(str, self._order[0][5])) + "\n"
        for index in range(len(self.LAST)):
            result += "  " + str(self.LAST[index]) + "\t\t" + str(self.COMBATH[index]) + "\t" + str(self.IRBATH[index]) + "\n"
        result += "END\n"
        return result


class PRESSURESCALE(_generic_imd_block):
    """PRESSURESCALE Block
        This block controls the barostat of the simulation

        Attributes
        -----------
        COUPLE: int
        SCALE:  int
        COMP:   float
        TAUP:   float
        VIRIAL: int
        SEMIANISOTROPIC:    List[int]
        PRES0:  List[List[float]]

    """
    name: str = "PRESSURESCALE"

    COUPLE: int
    SCALE: int
    COMP: float
    TAUP: float
    VIRIAL: int
    SEMIANISOTROPIC: List[int]
    PRES0: List[List[float]]

    _order = [[["COUPLE", "SCALE", "COMP", "TAUP", "VIRIAL"], ["SEMIANISOTROPIC COUPLINGS(X, Y, Z)"], ["PRES0(1...3,1...3)"]]]

    def __init__(self, COUPLE:int, SCALE:int, COMP:float, TAUP:float, VIRIAL:int, SEMIANISOTROPIC:List[int], PRES0:List[List[float]]):
        super().__init__(used=True)
        self.COUPLE = COUPLE
        self.SCALE = SCALE
        self.COMP = COMP
        self.TAUP = TAUP
        self.VIRIAL = VIRIAL
        self.SEMIANISOTROPIC = SEMIANISOTROPIC
        self.PRES0 = PRES0


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

    def __init__(self, BONDS: bool, ANGLES: bool, IMPROPER: bool, DIHEDRAL: bool, ELECTROSTATIC: bool, VDW: bool, NEGR: int, NRE: List[int]):
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
        super().__init__(used=True)
        self.BONDS = BONDS
        self.ANGLES = ANGLES
        self.IMPROPER = IMPROPER
        self.DIHEDRAL = DIHEDRAL
        self.ELECTROSTATIC = ELECTROSTATIC
        self.VDW = VDW
        self.NEGR = NEGR
        self.NRE = NRE

    def adapt_energy_groups(self, residues: Dict[str, Dict[int, int]]):
        """adapt_energy_groups
            This method is very "crude" and will put each residue into an own energy group.

        Parameters
        ----------
        residues :  Dict[str, Dict[int, int]]
            you can get this dict from cnf class.

        Returns
        -------
        None

        """

        print(residues)

        # set Energy Group amount
        self.NEGR = len(residues)  # total amount of Energy groups
        if("SOLV" in residues and residues["SOLV"] == 0):
            self.NEGR -= 1
        # if("WAT" in residues):
        #     wat_key = list(residues["WAT"].keys())[0]
        #     if(len(residues["WAT"]) == 1):
        #         self.NEGR += (int)(residues["WAT"][wat_key] / 3) -1
        #     else:
        #         self.NEGR += len(residues["WAT"]) -1

        # make residues sortable,
        dict_m = {}
        for x in sorted(residues):
            if (type(residues[x]) == dict and not x in dict_m and not x =="WAT"):
                if len(residues[x]) > 1: # If multiple residues have same name in a row (typically, counter-ions)
                    dict_m.update({min(residues[x].keys()):len(residues[x])}) # Set them in a single energy group
                else:
                    dict_m.update(residues[x])
            elif (x in dict_m and not x == "SOLV" and not x == "SOL"):
                raise Exception("Found mutliple residues for the same residue id!")

        # build up NRE
        self.NRE = []
        count = 0
        for x in sorted(dict_m.keys()):
            val = dict_m[x]
            count += val
            self.NRE.append(count)

        # add solvents to NRE, but only if there are solvents
        if ("SOLV" in residues and residues["SOLV"] != 0):
            self.NRE.append(residues["SOLV"] + count)

        if ("WAT" in residues and residues["WAT"] != 0):
            self.NRE.append(sum(list(residues["WAT"].values()))+ count)

        # if("WAT" in residues and residues["WAT"][wat_key] != 0):
        #     self.NRE = self.NRE[:-1]
        #     lim = self.NRE[len(self.NRE)-1] + residues["WAT"][wat_key]+1
        #     for i in range(self.NRE[len(self.NRE)-1]+3, lim, 3):
        #         self.NRE.append(i)


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
    NTCP0: int
    NTCS: int
    NTCS0: int

    _order = [[["NTC"], ["NTCP", "NTCP0(1)"], ["NTCS", "NTCS0(1)"]]]

    def __init__(self, NTC, NTCP, NTCP0, NTCS, NTCS0):
        """
        Args:
            NTC:
            NTCP:
            NTCP0:
            NTCS:
            NTCS0:
        """
        super().__init__(used=True)
        self.NTC = NTC
        self.NTCP = NTCP
        self.NTCP0 = NTCP0
        self.NTCS = NTCS
        self.NTCS0 = NTCS0


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
    SIZE: (str or float)
    TYPE: (str, bool)

    _order = [[["ALGORITHM", "NSNB", "RCUTP", "RCUTL", "SIZE", "TYPE"]]]

    def __init__(self, ALGORITHM: int, NSNB: int, RCUTP: float, RCUTL: float, SIZE: (str or float), TYPE: (str or bool)):
        """
        Args:
            ALGORITHM:
            NSNB:
            RCUTP:
            RCUTL:
            SIZE:
            TYPE:
        """
        super().__init__(used=True)
        self.ALGORITHM = ALGORITHM
        self.NSNB = NSNB
        self.RCUTP = RCUTP
        self.RCUTL = RCUTL
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
        4 : reaction-field with cutoff
        5 : reaction-field with cutoff + extra output
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
    NSHAPE: float
    ASHAPE: float
    NA2CLC: int
    TOLA2: float
    EPSLS: float
    NKX: float
    NKY: float
    NKZ: float
    KCUT: float
    NGX: float
    NGY: float
    NGZ: float
    NASORD: int
    NFDORD: int
    NALIAS: float
    NSPORD: float
    NQEVAL: float
    FACCUR: float
    NRDGRD: bool
    NWRGRD: bool
    NLRLJ: bool
    SLVDNS: float

    _order = [[["NLRELE"], ["APPAK", "RCRF", "EPSRF", "NSLFEXCL"], ["NSHAPE", "ASHAPE", "NA2CLC", "TOLA2", "EPSLS"],
               ["NKX", "NKY", "NKZ", "KCUT"], ["NGX", "NGY", "NGZ", "NASORD", "NFDORD", "NALIAS", "NSPORD"],
               ["NQEVAL", "FACCUR", "NRDGRD", "NWRGRD"], ["NLRLJ", "SLVDNS"]]]

    def __init__(self, NLRELE: int, APPAK: float, RCRF: float, EPSRF: float, NSLFEXCL: bool, NSHAPE: float, ASHAPE: float, NA2CLC: int, TOLA2: float,
                 EPSLS: float,
                 NKX: float, NKY: float, NKZ: float, KCUT: float, NGX: float, NGY: float, NGZ: float, NASORD: int,
                 NFDORD: int, NALIAS: float, NSPORD: float, NQEVAL: float, FACCUR: float, NRDGRD: bool, NWRGRD: bool, NLRLJ: bool, SLVDNS: float):
        """
        Args:
            NLRELE:
            APPAK:
            RCRF:
            EPSRF:
            NSLFEXCL:
            NSHAPE:
            ASHAPE:
            NA2CLC:
            TOLA2:
            EPSLS:
            NKX:
            NKY:
            NKZ:
            KCUT:
            NGX:
            NGY:
            NGZ:
            NASORD:
            NFDORD:
            NALIAS:
            NSPORD:
            NQEVAL:
            FACCUR:
            NRDGRD:
            NWRGRD:
            NLRLJ:
            SLVDNS:
        """
        super().__init__(used=True)
        self.NLRELE = NLRELE
        self.APPAK = APPAK
        self.RCRF = RCRF
        self.EPSRF = EPSRF
        self.NSLFEXCL = NSLFEXCL
        self.NSHAPE = NSHAPE
        self.ASHAPE = ASHAPE
        self.NA2CLC = NA2CLC
        self.TOLA2 = TOLA2
        self.EPSLS = EPSLS
        self.NKX = NKX
        self.NKY = NKY
        self.NKZ = NKZ
        self.KCUT = KCUT
        self.NGX = NGX
        self.NGY = NGY
        self.NGZ = NGZ
        self.NASORD = NASORD
        self.NFDORD = NFDORD
        self.NALIAS = NALIAS
        self.NSPORD = NSPORD
        self.NQEVAL = NQEVAL
        self.FACCUR = FACCUR
        self.NRDGRD = NRDGRD
        self.NWRGRD = NWRGRD
        self.NLRLJ = NLRLJ
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

    def __init__(self, NTIVEL: bool, NTISHK: int, NTINHT: bool, NTINHB: bool, NTISHI: bool, NTIRTC: bool, NTICOM: int, NTISTI: bool, IG: int,
                 TEMPI: float):
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
        super().__init__(used=True)
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

    def __init__(self, NSCM):
        super().__init__(used=True)
        self.NSCM = NSCM


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
        FUNCTIONAL: int, optional
            1: Single s Hamiltonian (default)
            2: Hamiltonian with NUMSTATES*(NUMSTATES-1)/2 (pairwise) S parameters ==> changes type of S
            3: Hamiltonian with (NUMSTATES-1) S parameters ==> changes type of S
    """

    _order = [[["EDS"], ["FUNCTIONAL FORM"], ["NUMSTATES"], ["S"], ["EIR"]]]

    def __init__(self, NUMSTATES: int, S: float, EIR: List[float], EDS: bool = 1, FUNCTIONAL: int = 1):
        super().__init__(used=True)
        self.name = "EDS"
        self.EDS = EDS
        self.FUNCTIONAL = FUNCTIONAL
        self.NUMSTATES = NUMSTATES
        self.S = S
        self.EIR = EIR


class DISTANCERES(_generic_imd_block):
    """DISTANCERES Block

    """
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

    def __init__(self, NTDIR: int, NTDIRA: int, CDIR: int, DIR0: int, TAUDIR: int, FORCESCALE: int, VDIR: int, NTWDIR: int):
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
        super().__init__(used=True)
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

        Attributes
        ----------
        NTPOR:
        NTPORB:
        NTPORS:
        CPOR:
    """

    name = "POSITIONRES"

    NTPOR: int
    NTPORB: int
    NTPORS: int
    CPOR: int

    _order = [[["NTPOR", "NTPORB", "NTPORS", "CPOR"]]]

    def __init__(self, NTPOR, NTPORB, NTPORS, CPOR):
        super().__init__(used=True)
        self.NTPOR = NTPOR
        self.NTPORB = NTPORB
        self.NTPORS = NTPORS
        self.CPOR = CPOR


class PRINTOUT(_generic_imd_block):
    """PRINTOUT block

    Attributes
    ----------
    NTPR: int
    NTPP: int

    """
    name: str = "PRINTOUT"

    NTPR: int
    NTPP: int

    _order = [[["NTPR", "NTPP"]]]

    def __init__(self, NTPR, NTPP):
        super().__init__(used=True)
        self.name = "PRINTOUT"
        self.NTPR = NTPR
        self.NTPP = NTPP


class ENERGYMIN(_generic_imd_block):
    """ENERGYMIN block

        Attributes
        ----------
        NTEM:   int
        NCYC:   int
        DELE:   float
        DX0:    float
        DXM:    float
        NMIN:   float
        FLIM:   float

    """
    name = "ENERGYMIN"

    NTEM: int
    NCYC: int
    DELE: float
    DX0: float
    DXM: float
    NMIN: float
    FLIM: float

    _order = [["NTEM    NCYC    DELE    DX0     DXM   NMIN   FLIM".split()]]

    def __init__(self, NTEM, NCYC, DELE, DX0, DXM, NMIN, FLIM):
        super().__init__(used=True)
        self.NTEM = NTEM
        self.NCYC = NCYC
        self.DELE = DELE
        self.DX0 = DX0
        self.DXM = DXM
        self.NMIN = NMIN
        self.FLIM = FLIM


class WRITETRAJ(_generic_imd_block):
    """WRITETRAJ

    Attributes
    ----------
    NTWX:   int
    NTWSE:   int
    NTWV:   int
    NTWF:   int
    NTWE:   int
    NTWG:   int
    NTWB:   int
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

    def __init__(self, NTWX: int = 0, NTWSE: int = 0, NTWV: int = 0, NTWF: int = 0, NTWE: int = 0, NTWG: int = 0, NTWB: int = 0):
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
        super().__init__(used=True)
        self.NTWX = NTWX
        self.NTWSE = NTWSE
        self.NTWV = NTWV
        self.NTWF = NTWF
        self.NTWE = NTWE
        self.NTWG = NTWG
        self.NTWB = NTWB


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

    def __init__(self, AMBER, AMBSCAL):
        super().__init__(used=True)
        self.name = "AMBER"
        self.AMBER = AMBER
        self.AMBSCAL = AMBSCAL


class COVALENTFORM(_generic_imd_block):
    """COVALENTFORM Block

        Attributes
        ----------
        NTBBH:  int
        NTBAH:  int
        NTBDN:  int
    """

    name = "COVALENTFORM"

    NTBBH: int
    NTBAH: int
    NTBDN: int

    _order = [[["NTBBH", "NTBAH", "NTBDN"]]]

    def __init__(self, NTBBH, NTBAH, NTBDN):
        super().__init__(used=True)
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

    def __init__(self, NTILM, NTILS, NGPUS):
        super().__init__(used=True)
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

    NTIL:   int
    NTLI:   List[int]
    NILG1:  List[int]
    NILG2:  List[int]
    ALI:    List[float]
    BLI:    List[float]
    CLI:    List[float]
    DLI:    List[float]
    ELI:    List[float]

    _order: List[List[str]] = [[["NTIL"], ["NTLI(1..)", "NILG1(1..)", "NILG2(1..)", "ALI(1..)", "BLI(1;;)", "CLI(1;;)", "DLI(1;;)", "ELI(1;;)"]]]

    def __init__(self, NTIL: int, NTLI: List[int], NILG1: List[int], NILG2: List[int], ALI: List[float], BLI: List[float], CLI: List[float], DLI: List[float], ELI: List[float]):
        super().__init__(used=True)
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
