import numpy as np
from scipy import constants


def get_density(mass: np.array, volume: np.array, atomu=1 / constants.Avogadro * 10e23) -> np.array:
    """
    Calculate the density

    Parameters
    ----------
    mass: np.array
    volume: np.array

    Returns
    -------
    np.array
        resulting density
    """
    return atomu * mass / volume


def get_Hvap(
    liq_totpot_energy: float, gas_totpot_energy: float, nMolecules=1, temperature=None, R=constants.R / 1000
) -> float:
    # calculate heat of vaporization
    rt_constant = R * temperature  # R in kilojoule_per_mole/kelvin * T
    heat_vap = gas_totpot_energy - liq_totpot_energy / nMolecules + rt_constant
    return heat_vap
