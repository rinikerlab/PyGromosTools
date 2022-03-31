"""
COM_TOP

Python version of the GROMOS++ program to combine two topologies.
All molecules are combined to a single topology and it can be decided where the solvent and parameters ar taken from.

This mehod can also be called directly via top1.com_top(top2), top1 + top2 or top1+=top2

Author: Marc Lehner
"""

from pygromos.utils.typing import List
from pygromos.files.topology.top import Top


def com_top(
    top1: Top,
    top2: Top,
    topo_multiplier: List[int] = [1, 1],
    solvFrom1: bool = True,
    paramFrom1: bool = True,
    verbose: bool = True,
) -> Top:
    """
        Python version of the GROMOS++ program to combine two topologies.
        All molecules are combined to a single topology and it can be decided where the solvent and parameters ar taken from.

        This mehod can also be called directly via top1.com_top(top2), top1 + top2 or top1+=top2

        Author: Marc Lehner

    Parameters
    ----------
    top1 : Top
        a topology
    top2 : Top
        another topology
    topo_multiplier : List[int]
        multiplier for topos. entries must be >=0
    solvFrom1 : bool, optional
        wether to take the solvent from topology 1, by default True
    paramFrom1 : bool, optional
        wether to use topo1 as main, by default True
    verbose : bool, optional
        additional proints for information, by default True

    Returns
    -------
    Top
        combined topology
    """
    # sanity checks
    if len(topo_multiplier) != 2:
        raise Exception("topo_multiplier length is not 2")
    if topo_multiplier[0] < 0 or topo_multiplier[1] < 0:
        raise Exception("Does not work with negative multipliers")

    # create the return top
    retTop = Top(in_value=None)

    if paramFrom1:
        retTop = top1.multiply_top(topo_multiplier[0], verbose=verbose)
        retTop = retTop._add_top(
            top=top2.multiply_top(topo_multiplier[1], verbose=verbose), solvFrom1=solvFrom1, verbose=verbose
        )
    else:
        retTop = top2.multiply_top(topo_multiplier[1], verbose=verbose)
        retTop = retTop._add_top(
            top=top1.multiply_top(topo_multiplier[0], verbose=verbose), solvFrom1=solvFrom1, verbose=verbose
        )
    return retTop
