from pygromos.files.blocks._general_blocks import _generic_gromos_block, _generic_field
from pygromos.utils.typing import Dict


class repex_system(_generic_gromos_block):
    def __init__(self, s: list, T: (float or list), state_eir: dict):
        super().__init__(used=True, name="REPLICAEXSYSTEM")
        self.T = T
        # self.lambda_values = lambda_values
        self.s = s
        self.state_eir = state_eir

    def block_to_string(self) -> str:
        # result = self.name + "\n"
        if not (type(self.T) is float):
            result = "Number of temperatures:\t " + str(len(self.T)) + "\n"
        else:
            result = "Number of temperatures:\t 1\n"
        result += "Dimension of temperature values:\t 1\n"
        result += "Number of lambda values: \t " + str(len(self.s)) + "\n"
        result += "T \t " + str(self.T) + "\n"
        result += "lambda \t " + " ".join(map(str, self.s)) + "\n"
        result += "s (RE-EDS) \t " + " ".join(map(str, self.s)) + "\n"

        for state in self.state_eir.keys():
            result += (
                "eir(s), numstate = " + str(state) + " (RE - EDS) " + "\t".join(map(str, self.state_eir[state])) + "\n"
            )

        # result ="END\n"
        return result


class replica_stat(_generic_field):
    ID: int
    partner: int
    run: int

    li: float
    Ti: float
    Epoti: float

    lj: float
    Tj: float
    Epotj: float

    p: float
    s: bool
    si: float
    sj: float

    state_potentials: Dict[str, float]

    def __init__(
        self,
        ID: int,
        partner: int,
        run: int,
        Ti: float,
        Epoti: float,
        Tj: float,
        Epotj: float,
        p: float,
        s: bool,
        si: float,
        sj: float,
        state_potentials: dict = None,
        potentials: dict = None,
        li: float = None,
        lj: float = None,
    ):
        """

        Parameters
        ----------
        ID :
        partner :
        run :
        Ti :
        Epoti :
        Tj :
        Epotj :
        p :
        s :
        si :
        sj :
        state_potentials :
        potentials :
        li :
        lj :
        """
        if isinstance(state_potentials, type(None)) and isinstance(potentials, type(None)):
            raise Exception("Need either state_potential or potentials for replica_statistic.")
        elif isinstance(state_potentials, type(None)) and not isinstance(potentials, type(None)):
            state_potentials = potentials

        self.ID = int(ID)
        self.partner = int(partner)
        self.run = int(run)

        self.li = float(li) if (li) else float(si)
        self.Ti = float(Ti)
        self.Epoti = float(Epoti)
        self.lj = float(lj) if (lj) else float(sj)
        self.Tj = float(Tj)
        self.Epotj = float(Epotj)

        self.p = float(p)
        self.s = int(s)
        self.si = float(si)
        self.sj = float(sj)

        self.state_potentials = dict(state_potentials)

    def to_string(self) -> str:
        return (
            "\t".join(
                [
                    str(self.ID),
                    str(self.partner),
                    str(self.run),
                    str(self.li),
                    str(self.Ti),
                    str(self.Epoti),
                    str(self.lj),
                    str(self.Tj),
                    str(self.Epotj),
                    str(self.p),
                    str(self.s),
                    str(self.si),
                    str(self.sj),
                    " ".join(list(map(str, self.state_potentials.values()))),
                ]
            )
            + "\n"
        )
