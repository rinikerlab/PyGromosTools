import numpy as np


def calculate_distance(atomA: np.array, atomB: np.array) -> np.array:

    return np.sqrt(np.sum((atomB - atomA) ** 2))


def rms(in_values) -> float:
    """helper function for RMSD. Calculates the root mean square for a array of (position/velocity) arrays

    Parameters
    ----------
    in_values : np.array of np.arrays

    Returns
    -------
    float
        root mean square
    """
    return np.sqrt(np.sum(np.square(in_values)) / len(in_values))


def periodic_distance(vec: np.array, grid: np.array) -> np.array:

    for i in range(3):
        if vec[i] > (grid[i] / 2):
            vec[i] = grid[i] - vec[i]
        elif vec[i] < (grid[i] / 2):
            vec[i] = grid[i] + vec[i]
    return vec


def periodic_shift(vec: np.array, grid: np.array) -> np.array:
    for i in range(3):
        if vec[i] > (grid[i] / 2):
            vec[i] = vec[i] - grid[i]
    return vec
