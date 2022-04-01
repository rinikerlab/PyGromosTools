---
date: '2022-04-01T13:55:18.499Z'
docname: _source/pygromos.files.topology
images: {}
path: /source-pygromos-files-topology
title: pygromos.files.topology package
---

# pygromos.files.topology package

## Submodules

## pygromos.files.topology.disres module


### _class_ pygromos.files.topology.disres.Disres(in_value: Optional[Union[str, Dict]] = None)
Bases: `pygromos.files.topology.disres.Distance_restraints`


### _class_ pygromos.files.topology.disres.Distance_restraints(in_value: Optional[Union[str, Dict]] = None)
Bases: `pygromos.files._basics._general_gromos_file._general_gromos_file`


#### read_blocks()

#### required_blocks(_ = ['TITLE', 'DISTANCERESPEC'_ )
## pygromos.files.topology.ifp module

File:
Description:

Author:


### _class_ pygromos.files.topology.ifp.Ifp(in_value: Union[str, Dict])
Bases: `pygromos.files._basics._general_gromos_file._general_gromos_file`


#### path(_: st_ )

#### read_file()
> give back the content. WARNING DEAPRECEATED.

**WARNING**: DEAPRECEATED


* **Returns**

    key is the block name of the gromos file, any is the content of a block



* **Return type**

    Dict[str, any]


## pygromos.files.topology.mtb module

File:
Description:

Author: Marc Lehner


### _class_ pygromos.files.topology.mtb.Mtb(in_value: Union[str, Dict], _future_file: bool = False)
Bases: `pygromos.files._basics._general_gromos_file._general_gromos_file`


#### all_res_names(_: Lis_ )

#### mtb_ends(_: Dict[str, [pygromos.files.blocks.mtb_blocks.MTBUILDBLEND](#pygromos.files.blocks.mtb_blocks.MTBUILDBLEND)_ )

#### mtb_solutes(_: Dict[str, [pygromos.files.blocks.mtb_blocks.MTBUILDBLSOLUTE](#pygromos.files.blocks.mtb_blocks.MTBUILDBLSOLUTE)_ )

#### mtb_solvents(_: Dict[str, [pygromos.files.blocks.mtb_blocks.MTBUILDBLSOLVENT](#pygromos.files.blocks.mtb_blocks.MTBUILDBLSOLVENT)_ )

#### read_file()
> give back the content. WARNING DEAPRECEATED.

**WARNING**: DEAPRECEATED


* **Returns**

    key is the block name of the gromos file, any is the content of a block



* **Return type**

    Dict[str, any]



#### read_mtb_file(path: str)
## pygromos.files.topology.ptp module


### _class_ pygromos.files.topology.ptp.Pertubation_topology(in_value: Optional[Union[str, dict]] = None)
Bases: `pygromos.files._basics._general_gromos_file._general_gromos_file`


#### MPERATOM(_: [pygromos.files.blocks.pertubation_blocks.MPERTATOM](#pygromos.files.blocks.pertubation_blocks.MPERTATOM_ )

#### PERTATOMPARAM(_: [pygromos.files.blocks.pertubation_blocks.PERTATOMPARAM](#pygromos.files.blocks.pertubation_blocks.PERTATOMPARAM_ )

#### PERTBONDANGLE(_: [pygromos.files.blocks.pertubation_blocks.PERTBONDANGLE](#pygromos.files.blocks.pertubation_blocks.PERTBONDANGLE_ )

#### PERTBONDANGLEH(_: [pygromos.files.blocks.pertubation_blocks.PERTBONDANGLEH](#pygromos.files.blocks.pertubation_blocks.PERTBONDANGLEH_ )

#### PERTBONDSTRETCH(_: [pygromos.files.blocks.pertubation_blocks.PERTBONDSTRETCH](#pygromos.files.blocks.pertubation_blocks.PERTBONDSTRETCH_ )

#### PERTBONDSTRETCHH(_: [pygromos.files.blocks.pertubation_blocks.PERTBONDSTRETCHH](#pygromos.files.blocks.pertubation_blocks.PERTBONDSTRETCHH_ )

#### PERTPROPERDIH(_: [pygromos.files.blocks.pertubation_blocks.PERTPROPERDIH](#pygromos.files.blocks.pertubation_blocks.PERTPROPERDIH_ )

#### TITLE(_: [pygromos.files.blocks._general_blocks.TITLE](#pygromos.files.blocks.pertubation_blocks.TITLE_ )

#### read_blocks()

#### required_blocks(_ = ['TITLE'_ )

### _class_ pygromos.files.topology.ptp.Ptp(in_value: Optional[Union[str, dict]] = None)
Bases: `pygromos.files.topology.ptp.Pertubation_topology`


#### MPERATOM(_: [pygromos.files.blocks.pertubation_blocks.MPERTATOM](#pygromos.files.blocks.pertubation_blocks.MPERTATOM_ )

#### PERTATOMPARAM(_: [pygromos.files.blocks.pertubation_blocks.PERTATOMPARAM](#pygromos.files.blocks.pertubation_blocks.PERTATOMPARAM_ )

#### PERTBONDANGLE(_: [pygromos.files.blocks.pertubation_blocks.PERTBONDANGLE](#pygromos.files.blocks.pertubation_blocks.PERTBONDANGLE_ )

#### PERTBONDANGLEH(_: [pygromos.files.blocks.pertubation_blocks.PERTBONDANGLEH](#pygromos.files.blocks.pertubation_blocks.PERTBONDANGLEH_ )

#### PERTBONDSTRETCH(_: [pygromos.files.blocks.pertubation_blocks.PERTBONDSTRETCH](#pygromos.files.blocks.pertubation_blocks.PERTBONDSTRETCH_ )

#### PERTBONDSTRETCHH(_: [pygromos.files.blocks.pertubation_blocks.PERTBONDSTRETCHH](#pygromos.files.blocks.pertubation_blocks.PERTBONDSTRETCHH_ )

#### PERTPROPERDIH(_: [pygromos.files.blocks.pertubation_blocks.PERTPROPERDIH](#pygromos.files.blocks.pertubation_blocks.PERTPROPERDIH_ )

#### TITLE(_: [pygromos.files.blocks._general_blocks.TITLE](#pygromos.files.blocks.pertubation_blocks.TITLE_ )

#### path(_: st_ )
## pygromos.files.topology.top module

File:            gromos++ topo file functions
Warnings: this CLASS IS NOT IMPLEMENTED!
TODO:REWORK
Description:

> in this lib, gromos topo file mainpulating functions are gathered

Author: Marc Lehner, Benjamin Ries


### _class_ pygromos.files.topology.top.Top(in_value: Union[str, dict, pygromos.utils.typing.Top_Type], _future_file: bool = False)
Bases: `pygromos.files._basics._general_gromos_file._general_gromos_file`


#### _add_top(top: Optional[pygromos.utils.typing.Top_Type], solvFrom1: bool = True, verbose: bool = False)
combines two topologies. Parameters are taken from the initial topology.
But missing parameters from the second topology will be added.
Can be used like com_top from Gromos++


* **Parameters**

    
    * **top** (*TopType*) – second topology to add to the first topology


    * **solvFrom1** (*bool, optional*) – should the solvent be taken from the first topology? (else second), by default True


    * **verbose** (*bool, optional*) – extra print statements, by default True


    * **sanityCheck** (*bool, optional*) – feature and compatibility check, by default True



* **Returns**

    returns a topology made by combing two topologies



* **Return type**

    TopType



#### add_new_CONSTRAINT(IC: int, JC: int, ICC: float, verbose: bool = False)
adds a CONSTRAINT entry to the topology


* **Parameters**

    
    * **IC** (*int*) – atom index I


    * **JC** (*int*) – atom index J


    * **ICC** (*float*) – constraint length


    * **verbose** (*bool, optional*)



#### add_new_LJparameter(C6: float, C12: float, CS6: float = 0, CS12: float = 0, combination_rule: str = 'geometric', verbose=False, AddATOMTYPENAME: Optional[str] = None, lowerBound: float = 1e-100)
add a LJ entry to the LJ parameter block


* **Parameters**

    
    * **C6** (*float*)


    * **C12** (*float*)


    * **CS6** (*float, optional*)


    * **CS12** (*float, optional*)


    * **combination_rule** (*str, optional*) – no other options supported rigth now, by default “geometric”


    * **AddATOMTYPENAME** (*str, optional*) – if not None a new atomtype is made, by default None


    * **lowerBound** (*float, optional*) – saftey, by default 1e-100



#### add_new_PRESSUREGROUPS(number: int, verbose: bool = False)

#### add_new_SOLUTEMOLECULES(number: int, verbose: bool = False)

#### add_new_TEMPERATUREGROUPS(number: int, verbose: bool = False)

#### add_new_angle(k: float, kh: float, b0: float, atomI: int, atomJ: int, atomK: int, includesH: bool = False, verbose: bool = False, convertToQuartic: bool = False)
add a angle between atom I, J and K to the ANGLE block


* **Parameters**

    
    * **k** (*float*) – force konstant


    * **kh** (*float*) – force konstant harmonic


    * **b0** (*float*) – angle at which the force is 0


    * **atomI** (*int*) – atom I


    * **atomJ** (*int*) – atom J


    * **atomK** (*int*) – atom K


    * **includesH** (*bool, optional*) – ANGLE or ANGLEH, by default False


    * **convertToQuartic** (*bool, optional*) – auto convert, by default False



#### add_new_atom(ATNM: int = 0, MRES: int = 0, PANM: str = '_', IAC: int = 1, MASS: float = 1.0, CG: int = 0, CGC: int = 1, INE: list = [], INE14: list = [], verbose=False, C6: Optional[float] = None, C12: Optional[float] = None, CS6: float = 0, CS12: float = 0, IACname: Optional[str] = None)
add a atom to a system (with creating a new atomtype if needed and adding LJ parameters if needed)


* **Parameters**

    
    * **ATNM** (*int, optional*) – number of the atom in the system, by default 0


    * **MRES** (*int, optional*) – residue number, by default 0


    * **PANM** (*str, optional*) – name of the atom, by default ‘_’


    * **IAC** (*int, optional*) – atomtype number of the atom, by default 1


    * **MASS** (*float, optional*) – mass of the atom, by default 1.0


    * **CG** (*int, optional*) – charge of the atom, by default 0


    * **CGC** (*int, optional*) – charge group bool, by default 1


    * **INE** (*list, optional*) – INE list, by default []


    * **INE14** (*list, optional*) – INE14 list, by default []


    * **C6** (*float, optional*) – C6 value, by default None


    * **C12** (*float, optional*) – C12 value, by default None


    * **CS6** (*float, optional*) – CS6 value, by default 0


    * **CS12** (*float, optional*) – CS12 value, by default 0


    * **IACname** (*str, optional*) – new IACname if NONE PANM is used, by default None



#### add_new_atomtype(name: str, verbose: bool = False)
add a atomtype to ATOMTYPENAME block


* **Parameters**

    
    * **name** (*str*) – new atomtype name


    * **verbose** (*bool, optional*) – by default False



#### add_new_bond(k: float, b0: float, atomI: int, atomJ: int, includesH: bool = False, verbose: bool = False)
add a bond between atom I and J to the BOND block


* **Parameters**

    
    * **k** (*float*) – force konstant


    * **b0** (*float*) – distance at which the force is 0


    * **atomI** (*int*) – atom I


    * **atomJ** (*int*) – atom J


    * **includesH** (*bool, optional*) – wheter it should be added to BOND or BONDH, by default False



#### add_new_crossdihedral(verbose: bool = False)

#### add_new_impdihedral(CQ: float, Q0: float, atomI: int, atomJ: int, atomK: int, atomL: int, includesH: bool = False, verbose: bool = False)
add a new impdihedral


* **Parameters**

    
    * **CQ** (*float*) – force constant


    * **Q0** (*float*) – Q0


    * **atomI** (*int*) – atom I


    * **atomJ** (*int*) – atom J


    * **atomK** (*int*) – atom K


    * **atomL** (*int*) – atom L


    * **includesH** (*bool, optional*) – IMPDIHEDRALH or IMPDIHEDRAL, by default False



#### add_new_impdihedral_type(CQ: float, Q0: float, verbose: bool = False)
add a new impodihedraltype


* **Parameters**

    
    * **CQ** (*float*) – force constant


    * **Q0** (*float*) – Q0



#### add_new_resname(name: str, verbose: bool = False)
add a resname to the RESNAME block


* **Parameters**

    
    * **name** (*str*) – resname name


    * **verbose** (*bool, optional*) – by default False



#### add_new_soluteatom(ATNM: int = 0, MRES: int = 0, PANM: str = '', IAC: int = 0, MASS: float = 0, CG: float = 0, CGC: int = 0, INE: list = [], INE14: list = [], verbose: bool = False)
add a soluteatom to the SOLUTEATOM block


#### add_new_torsiondihedral(CP: float, PD: float, NP: int, atomI: int, atomJ: int, atomK: int, atomL: int, includesH: bool = False, verbose: bool = False)
add a torsiondihedral between atom I, J, K and L to the TORSIONDIHEDRAL block


* **Parameters**

    
    * **CP** (*float*) – force constant


    * **PD** (*float*) – phase


    * **NP** (*int*) – multiplicity


    * **atomI** (*int*) – atom I


    * **atomJ** (*int*) – atom J


    * **atomK** (*int*) – atom K


    * **atomL** (*int*) – atom L


    * **includesH** (*bool, optional*) – DIHEDRAL or DIHEDRALH, by default False



#### find_LJparameterNumber(C12: float, C6: float)
find the LJ parameter number


#### get_LJparameter_from_IAC(IAC: int)
get the LJ parameter from the IAC number


* **Parameters**

    **IAC** (*int*) – [description]



#### get_diff_to_top(top: pygromos.utils.typing.Top_Type)

#### get_mass()
Calculates the total mass of the solute molecule


* **Returns**

    total mass in a.u.



* **Return type**

    float



#### get_num_atomtypes()

#### harmonic2quarticAngleConversion(kh: float, b0: float)
conversion of a harmonic bondanglebending force constant to a cubic in cosine/quartic one


* **Parameters**

    
    * **kh** (*float*) – harmonic bondanglebending force constant (CHT)


    * **b0** (*float*) – bondangle 0



* **Returns**

    cubic in cosine force constant (CT)



* **Return type**

    float



#### make_ordered(orderList: Optional[list] = None)

#### multiply_top(n_muliplication: int, unifyGroups: bool = False, verbose: bool = False)

#### read_file()
> give back the content. WARNING DEAPRECEATED.

**WARNING**: DEAPRECEATED


* **Returns**

    key is the block name of the gromos file, any is the content of a block



* **Return type**

    Dict[str, any]



### pygromos.files.topology.top.check_top()

### pygromos.files.topology.top.combine_topologies()

### pygromos.files.topology.top.make_topolog(input_arg, build, param, seq, solve='H2O')
## Module contents
