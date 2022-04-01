---
date: '2022-04-01T13:55:18.499Z'
docname: _source/pygromos.gromos.pyGromosPP
images: {}
path: /source-pygromos-gromos-py-gromos-pp
title: pygromos.gromos.pyGromosPP package
---

# pygromos.gromos.pyGromosPP package

## Submodules

## pygromos.gromos.pyGromosPP.com_top module

COM_TOP

Python version of the GROMOS++ program to combine two topologies.
All molecules are combined to a single topology and it can be decided where the solvent and parameters ar taken from.

This mehod can also be called directly via top1.com_top(top2), top1 + top2 or top1+=top2

Author: Marc Lehner


### pygromos.gromos.pyGromosPP.com_top.com_top(top1: [pygromos.files.topology.top.Top](#pygromos.files.topology.top.Top), top2: [pygromos.files.topology.top.Top](#pygromos.files.topology.top.Top), topo_multiplier: List[int] = [1, 1], solvFrom1: bool = True, paramFrom1: bool = True, verbose: bool = True)
> Python version of the GROMOS++ program to combine two topologies.
> All molecules are combined to a single topology and it can be decided where the solvent and parameters ar taken from.

> This mehod can also be called directly via top1.com_top(top2), top1 + top2 or top1+=top2

> Author: Marc Lehner


* **Parameters**

    
    * **top1** (*Top*) – a topology


    * **top2** (*Top*) – another topology


    * **topo_multiplier** (*List[int]*) – multiplier for topos. entries must be >=0


    * **solvFrom1** (*bool, optional*) – wether to take the solvent from topology 1, by default True


    * **paramFrom1** (*bool, optional*) – wether to use topo1 as main, by default True


    * **verbose** (*bool, optional*) – additional proints for information, by default True



* **Returns**

    combined topology



* **Return type**

    [Top](#pygromos.files.topology.top.Top)


## pygromos.gromos.pyGromosPP.ran_box module

Python implementation of the Gromos++ program ran_box which is used to generate randomized configurations for liquids (and gases)

Author: Marc Lehner


### pygromos.gromos.pyGromosPP.ran_box.ran_box(in_top_path: str, in_cnf_path: str, out_cnf_path: str = '', periodic_boundary_condition: str = 'r', nmolecule: int = 1, dens: float = 1.0, threshold: Optional[float] = None, layer: bool = False, boxsize: Optional[float] = None, fixfirst: bool = False, seed: Optional[float] = None, _binary_name: str = 'ran_box', verbose: bool = True, return_command_only: bool = False)
## Module contents
