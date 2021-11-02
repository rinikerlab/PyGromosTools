import numpy as np



def get_density(mass:np.array, volume:np.array)->np.array:
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
    return 1.66056 * mass / volume

from scipy import constants

def get_Hvap(liq_nonbonded_energies:np.array, gas_nonbonded_energies:np.array,
             nMolecules=1, temperature=None, R=constants.R/1000) -> float:
    # calculate heat of vaporization
    rt_constant = R * temperature # R in kilojoule_per_mole/kelvin * T
    heat_vap = gas_nonbonded_energies - liq_nonbonded_energies / nMolecules + rt_constant
    return heat_vap