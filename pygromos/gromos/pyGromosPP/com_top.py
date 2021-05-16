""" 
COM_TOP

Python version of the GROMOS++ program to combine two topologies. 
All molecules are combined to a single topology and it can be decided where the solvent and parameters ar taken from.

This mehod can also be called directly via top1.com_top(top2), top1 + top2 or top1+=top2

Author: Marc Lehner
"""

from copy import deepcopy
from typing import Union
from pygromos.files.topology.top import Top

def com_top(top1:Top, top2:Top, solvFrom1:bool=True, paramFrom1:bool=True, verbose:bool=True) -> Top:
    """[summary]

    Parameters
    ----------
    top1 : Top
        a topology
    top2 : Top
        another topology
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
    # create the return top
    if paramFrom1:
        return top1._add_top(top=top2, solvFrom1=solvFrom1, verbose=verbose)
    else:
        return top2._add_top(top=top1, solvFrom1=solvFrom1, verbose=verbose)

    