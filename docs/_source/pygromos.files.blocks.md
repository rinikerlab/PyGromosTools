---
date: '2022-04-01T13:55:18.499Z'
docname: _source/pygromos.files.blocks
images: {}
path: /source-pygromos-files-blocks
title: pygromos.files.blocks package
---

# pygromos.files.blocks package

## Submodules

## pygromos.files.blocks.coord_blocks module


### _class_ pygromos.files.blocks.coord_blocks.GENBOX(pbc: pygromos.files.blocks.coord_blocks.Pbc = Pbc.vacuum, length: List[float] = [0.0, 0.0, 0.0], angles: List[float] = [0.0, 0.0, 0.0], euler: List[float] = [0.0, 0.0, 0.0], origin: List[float] = [0.0, 0.0, 0.0], content=None)
Bases: `pygromos.files.blocks._general_blocks._generic_gromos_block`

> This Block is representing the simulation Box in a coordinate file.


* **Variables**

    
    * **~GENBOX.pbc** (*int**,**Pbc*) – Periodic Boundary Condition


    * **~GENBOX.length** (*List**[**float**]*) – 


    * **~GENBOX.angles** (*List**[**float**]*) – 


    * **~GENBOX.euler** (*List**[**float**]*) – 


    * **~GENBOX.origin** (*List**[**float**]*) – 



#### \__init__(pbc: pygromos.files.blocks.coord_blocks.Pbc = Pbc.vacuum, length: List[float] = [0.0, 0.0, 0.0], angles: List[float] = [0.0, 0.0, 0.0], euler: List[float] = [0.0, 0.0, 0.0], origin: List[float] = [0.0, 0.0, 0.0], content=None)

* **Parameters**

    
    * **pbc** (*int,Pbc*)


    * **length** (*List[float]*)


    * **angles** (*List[float]*)


    * **euler** (*List[float]*)


    * **origin** (*List[float]*)



#### _property_ angles(_: List[float_ )

#### block_to_string()

#### comment(_: st_ )

#### content(_: Iterabl_ )

#### _property_ euler(_: List[float_ )

#### _property_ length(_: List[float_ )

#### _property_ origin(_: List[float_ )

#### _property_ pbc(_: pygromos.files.blocks.coord_blocks.Pb_ )

#### read_content_from_str(content: List[str])

### _class_ pygromos.files.blocks.coord_blocks.LATTICESHIFTS(content: List[pygromos.files.blocks.coord_blocks.lattice_shift])
Bases: `pygromos.files.blocks._general_blocks._iterable_gromos_block`


* **Parameters**

    **content** (*List[lattice_shift]*) – every element in this list is a lattice shift obj



#### \__init__(content: List[pygromos.files.blocks.coord_blocks.lattice_shift])

* **Parameters**

    **content** (*List[lattice_shift]*) – every element in this list is a lattice shift obj



#### comment(_: st_ )

#### read_content_from_str(content: List[str])

### _class_ pygromos.files.blocks.coord_blocks.PERTDATA(content: List[str])
Bases: `pygromos.files.blocks._general_blocks._generic_gromos_block`


#### \__init__(content: List[str])
> This block is used for lambda-sampling and gives the lambda value of the current coordinates.


* **Parameters**

    **lambda_value** (*float*) – current lambda value



#### content(_: floa_ )

#### _property_ lam(_: floa_ )

#### read_content_from_str(content: List[str])

### _class_ pygromos.files.blocks.coord_blocks.POSITION(content: List[pygromos.files.blocks.coord_blocks.atomP])
Bases: `pygromos.files.blocks._general_blocks._iterable_gromos_block`


* **Parameters**

    **content** (*List[atomP]*) – every element in this list is of atom position obj



#### comment(_: st_ )

#### read_content_from_str(content: List[str])

### _class_ pygromos.files.blocks.coord_blocks.POSRESSPEC(content: List[pygromos.files.blocks.coord_blocks.atomP])
Bases: `pygromos.files.blocks._general_blocks._iterable_gromos_block`

POSITION


* **Parameters**

    **content** (*List[atomP]*) – every element in this list is of atom position obj



#### comment(_: st_ )

#### read_content_from_str(content: List[str])

### _class_ pygromos.files.blocks.coord_blocks.Pbc(value)
Bases: `enum.Enum`

An enumeration.


#### rectangular(_ = _ )

#### triclinic(_ = _ )

#### trunc_octahedron(_ = -_ )

#### vacuum(_ = _ )

### _class_ pygromos.files.blocks.coord_blocks.REFPOSITION(content: List[pygromos.files.blocks.coord_blocks.atomP])
Bases: `pygromos.files.blocks._general_blocks._iterable_gromos_block`


* **Parameters**

    **content** (*List[atomP]*) – every element in this list is of atom position obj



#### \__init__(content: List[pygromos.files.blocks.coord_blocks.atomP])

* **Parameters**

    **content** (*List[atomP]*) – every element in this list is of atom position obj



#### comment(_: st_ )

#### read_content_from_str(content: List[str])

### _class_ pygromos.files.blocks.coord_blocks.STOCHINT(content: List[pygromos.files.blocks.coord_blocks.atomSI])
Bases: `pygromos.files.blocks._general_blocks._iterable_gromos_block`


* **Parameters**

    
    * **content** (*List[atomSI]*) – every element in this list is of atom stochastic interval obj


    * **seed** (*str*) – contains the seed for the stochastic dynamics simulation



#### block_to_string()

#### comment(_: st_ )

#### read_content_from_str(content: List[str])

### _class_ pygromos.files.blocks.coord_blocks.VELOCITY(content: List[pygromos.files.blocks.coord_blocks.atomV])
Bases: `pygromos.files.blocks._general_blocks._iterable_gromos_block`


* **Parameters**

    **content** (*List[atomV]*) – every element in this list is of atom velocity obj



#### comment(_: st_ )

#### read_content_from_str(content: List[str])

### _class_ pygromos.files.blocks.coord_blocks.atomP(resID: int, resName: str, atomType: str, atomID: int, xp: float, yp: float, zp: float)
Bases: `pygromos.files.blocks._general_blocks._generic_field`


#### \__init__(resID: int, resName: str, atomType: str, atomID: int, xp: float, yp: float, zp: float)

* **Parameters**

    
    * **resID**


    * **resName**


    * **atomType**


    * **atomID**


    * **xp**


    * **yp**


    * **zp**



#### to_string()

### _class_ pygromos.files.blocks.coord_blocks.atomSI(resID: int, resName: str, atomType: str, atomID: int, sxx: float, sxy: float, sxz: float)
Bases: `pygromos.files.blocks._general_blocks._generic_field`


#### to_string()

### _class_ pygromos.files.blocks.coord_blocks.atomV(resID: int, resName: str, atomType: str, atomID: int, xv: float, yv: float, zv: float)
Bases: `pygromos.files.blocks._general_blocks._generic_field`


#### to_string()

### _class_ pygromos.files.blocks.coord_blocks.lattice_shift(atomID: int, x: float, y: float, z: float)
Bases: `pygromos.files.blocks._general_blocks._generic_field`


#### to_string()
## pygromos.files.blocks.imd_blocks module


### _class_ pygromos.files.blocks.imd_blocks.AMBER(AMBER: int = 0, AMBSCAL: int = 0, content: Optional[List[str]] = None)
Bases: `pygromos.files.blocks.imd_blocks._generic_imd_block`

AMBER block


* **Variables**

    
    * **~AMBER.Amber** (*bool*) – 


    * **~AMBER.AMBSCAL** (*float*) – 



#### AMBSCAL(_: floa_ )

#### Amber(_: boo_ )

#### name(_: st_ _ = 'AMBER_ )

### _class_ pygromos.files.blocks.imd_blocks.BOUNDCOND(NTB: int = 0, NDFMIN: int = 0, content: Optional[List[str]] = None)
Bases: `pygromos.files.blocks.imd_blocks._generic_imd_block`

Boundary Condition Block

> This block describes the boundary condition of the coordinate system.


* **Variables**

    
    * **~BOUNDCOND.NTB** (*int**, **optional*) – Boundary conditions, by default 0
    -1 : truncated octahedral

    > 0 : vacuum
    > 1 : rectangular
    > 2 : triclinic



    * **~BOUNDCOND.NDFMIN** (*int**, **optional*) – number of degrees of freedom subtracted for temperature, by default 0



#### NDFMIN(_: in_ )

#### NTB(_: in_ )

#### name(_: st_ _ = 'BOUNDCOND_ )

### _class_ pygromos.files.blocks.imd_blocks.COMTRANSROT(NSCM: int = 0, content: Optional[List[str]] = None)
Bases: `pygromos.files.blocks.imd_blocks._generic_imd_block`

COMTRANSROT block

This block controls the center of mass translation and rotation removal. (flying ice cube problem)


* **Variables**

    **~COMTRANSROT.NSCM** (*int*) – controls system centre-of-mass (com) motion removal
    0: no com motion removal (default)
    < 0: com translation and rotation are removed every abs(NSCM) steps.
    > 0: com translation is removed every NSCM steps.



#### NSCM(_: in_ )

#### name(_: st_ _ = 'COMTRANSROT_ )

### _class_ pygromos.files.blocks.imd_blocks.CONSTRAINT(NTC: int = 0, NTCP: int = 0, NTCP0: float = 0, NTCS: int = 0, NTCS0: float = 0, content: Optional[List[str]] = None)
Bases: `pygromos.files.blocks.imd_blocks._generic_imd_block`

CONSTRAINT block

    This block is controlling constraining the atoms during a simulation.


* **Variables**

    
    * **~CONSTRAINT.NTC** (*int*) – 


    * **~CONSTRAINT.NTCP** (*int*) – 


    * **~CONSTRAINT.NTCP0** (*int*) – 


    * **~CONSTRAINT.NTCS** (*int*) – 


    * **~CONSTRAINT.NTCS0** (*int*) – 



#### NTC(_: in_ )

#### NTCP(_: in_ )

#### NTCP0(_: floa_ )

#### NTCS(_: in_ )

#### NTCS0(_: floa_ )

#### \__init__(NTC: int = 0, NTCP: int = 0, NTCP0: float = 0, NTCS: int = 0, NTCS0: float = 0, content: Optional[List[str]] = None)
Args:

    NTC:
    NTCP:
    NTCP0:
    NTCS:
    NTCS0:


#### name(_: st_ _ = 'CONSTRAINT_ )

### _class_ pygromos.files.blocks.imd_blocks.COVALENTFORM(NTBBH: int = 0, NTBAH: int = 0, NTBDN: int = 0, content: Optional[List[str]] = None)
Bases: `pygromos.files.blocks.imd_blocks._generic_imd_block`

COVALENTFORM Block

> The COVALENTFORM Block manages the functional form of the Bonded contributions to the force field.
> It is optional in an imd file to define the block.


* **Variables**

    
    * **~COVALENTFORM.NTBBH** (*int*) – 0,1 controls bond-stretching potential
    0: quartic form (default)
    1: harmonic form


    * **~COVALENTFORM.NTBAH** (*int*) – 0,1 controls bond-angle bending potential
    0: cosine-harmonic (default)
    1: harmonic


    * **~COVALENTFORM.NTBDN** (*int*) – 0,1 controls torsional dihedral potential
    0: arbitrary phase shifts (default)
    1: phase shifts limited to 0 and 180 degrees.



#### NTBAH(_: in_ )

#### NTBBH(_: in_ )

#### NTBDN(_: in_ )

#### name(_ = 'COVALENTFORM_ )

### _class_ pygromos.files.blocks.imd_blocks.DISTANCERES(NTDIR: int = 0, NTDIRA: int = 0, CDIR: int = 0, DIR0: int = 0, TAUDIR: int = 0, FORCESCALE: int = 0, VDIR: int = 0, NTWDIR: int = 0, content: Optional[List[str]] = None)
Bases: `pygromos.files.blocks.imd_blocks._generic_imd_block`

DISTANCERES Block


#### CDIR(_: in_ )

#### DIR0(_: in_ )

#### FORCESCALE(_: in_ )

#### NTDIR(_: in_ )

#### NTDIRA(_: in_ )

#### NTWDIR(_: in_ )

#### TAUDIR(_: in_ )

#### VDIR(_: in_ )

#### \__init__(NTDIR: int = 0, NTDIRA: int = 0, CDIR: int = 0, DIR0: int = 0, TAUDIR: int = 0, FORCESCALE: int = 0, VDIR: int = 0, NTWDIR: int = 0, content: Optional[List[str]] = None)
Args:

    NTDIR:
    NTDIRA:
    CDIR:
    DIR0:
    TAUDIR:
    FORCESCALE:
    VDIR:
    NTWDIR:


#### name(_ = 'DISTANCERES_ )

### _class_ pygromos.files.blocks.imd_blocks.EDS(NUMSTATES: int = 0, S: float = 0, EIR: List[float] = [], EDS: bool = 1, ALPHLJ: float = 0.0, ALPHC: float = 0.0, FUNCTIONAL: int = 1, content: Optional[List[str]] = None)
Bases: `pygromos.files.blocks.imd_blocks._generic_imd_block`

EDS block

This block is used in an EDS simulation.


* **Variables**

    
    * **~EDS.NUMSTATES** (*int*) – EDS-States


    * **~EDS.S** (*float*) – smoothness parameter


    * **~EDS.EIR** (*List**[**float**]*) – energy offsets


    * **~EDS.EDS** (*bool**, **optional*) – turn on EDS_simulation


    * **~EDS.ALPHLJ** (*float**, **optional*) – 


    * **~EDS.ALPHC** (*float**, **optional*) – 


    * **~EDS.FUNCTIONAL** (*int**, **optional*) – 1: Single s Hamiltonian (default)
    2: Hamiltonian with NUMSTATES\*(NUMSTATES-1)/2 (pairwise) S parameters ==> changes type of S
    3: Hamiltonian with (NUMSTATES-1) S parameters ==> changes type of S



#### read_content_from_str(content: List[str])

### _class_ pygromos.files.blocks.imd_blocks.ENERGYMIN(NTEM: int = 0, NCYC: int = 0, DELE: float = 0, DX0: float = 0, DXM: float = 0, NMIN: int = 0, FLIM: float = 0, content: Optional[List[str]] = None)
Bases: `pygromos.files.blocks.imd_blocks._generic_imd_block`

ENERGYMIN block
This block takes care of managing the Energyminimization controls.


* **Variables**

    
    * **~ENERGYMIN.NTEM** (*int*) – controls energy minimisation mode.
    0: do not do energy minimisation (default)
    1: steepest-descent minimisation
    2: Fletcher-Reeves conjugate-gradient minimisation
    3: Polak-Ribiere conjugate-gradient minimisation


    * **~ENERGYMIN.NCYC** (*int*) – >0 number of steps before resetting the conjugate-gradient search direction
    =0 reset only if the energy grows in the search direction


    * **~ENERGYMIN.DELE** (*float*) – >0.0 energy threshold for convergence
    >0.0 (conjugate-gradient) RMS force threshold for convergence


    * **~ENERGYMIN.DX0** (*float*) – >0.0 initial step size


    * **~ENERGYMIN.DXM** (*float*) – >0.0 maximum step size


    * **~ENERGYMIN.NMIN** (*float*) – >0 minimum number of minimisation steps


    * **~ENERGYMIN.FLIM** (*float*) – >=0.0 limit force to maximum value (FLIM > 0.0 is not recommended)



#### DELE(_: floa_ )

#### DX0(_: floa_ )

#### DXM(_: floa_ )

#### FLIM(_: floa_ )

#### NCYC(_: in_ )

#### NMIN(_: in_ )

#### NTEM(_: in_ )

#### name(_ = 'ENERGYMIN_ )

### _class_ pygromos.files.blocks.imd_blocks.FORCE(BONDS: bool = True, ANGLES: bool = True, IMPROPER: bool = True, DIHEDRAL: bool = True, ELECTROSTATIC: bool = True, VDW: bool = True, NEGR: int = 0, NRE: List[int] = [], content: Optional[List[str]] = None)
Bases: `pygromos.files.blocks.imd_blocks._generic_imd_block`

FORCE Block

    this Block is controlling the forcefield options. It can turn on and of terms, as well as generate force sub groups for analysis.


* **Variables**

    
    * **~FORCE.BONDS** (*bool*) – 


    * **~FORCE.ANGLES** (*bool*) – 


    * **~FORCE.IMPROPER** (*bool*) – 


    * **~FORCE.DIHEDRAL** (*bool*) – 


    * **~FORCE.ELECTROSTATIC** (*bool*) – 


    * **~FORCE.VDW** (*bool*) – 


    * **~FORCE.NEGR** (*int*) – Number of Energy subgroups


    * **~FORCE.NRE** (*List**[**int**]*) – List of last atoms for Energy subgroups. (NRE len == NEGR)



#### ANGLES(_: boo_ )

#### BONDS(_: boo_ )

#### DIHEDRAL(_: boo_ )

#### ELECTROSTATIC(_: boo_ )

#### IMPROPER(_: boo_ )

#### NEGR(_: in_ )

#### NRE(_: List[int_ )

#### VDW(_: boo_ )

#### \__init__(BONDS: bool = True, ANGLES: bool = True, IMPROPER: bool = True, DIHEDRAL: bool = True, ELECTROSTATIC: bool = True, VDW: bool = True, NEGR: int = 0, NRE: List[int] = [], content: Optional[List[str]] = None)
Args:

    BONDS:
    ANGLES:
    IMPROPER:
    DIHEDRAL:
    ELECTROSTATIC:
    VDW:
    NEGR:
    NRE (list):


#### adapt_energy_groups(energy_groups: Dict[int, int])
> Change the Force groups.


* **Parameters**

    **residues** (*Dict[int, int]*) –  the dictionary contains the forceGroups with the ID of the respective last atom.

    {forceGroupID: lastAtom} = {1:12, 2:26, 3:128}



#### name(_: st_ _ = 'FORCE_ )

#### read_content_from_str(content: List[str])

### _class_ pygromos.files.blocks.imd_blocks.INITIALISE(NTIVEL: bool = False, NTISHK: int = 0, NTINHT: bool = False, NTINHB: bool = False, NTISHI: bool = False, NTIRTC: bool = False, NTICOM: int = 0, NTISTI: bool = False, IG: int = 0, TEMPI: float = 0, content: Optional[List[str]] = None)
Bases: `pygromos.files.blocks.imd_blocks._generic_imd_block`

INITIALISE block

> This block controls the Initialisation of a simulation.


* **Variables**

    
    * **~INITIALISE.NTIVEL** (*bool*) – 0,1 controls generation of initial velocities
    0: read from configuration (default)
    1: generate from Maxell distribution at temperature TEMPI


    * **~INITIALISE.NTISHK** (*int*) – 0..3 controls shaking of initial configuration
    0: no intial SHAKE (default)
    1: initial SHAKE on coordinates only
    2: initial SHAKE on velocities only
    3: initial SHAKE on coordinates and velocities


    * **~INITIALISE.NTINHT** (*bool*) – 0,1 controls generation of initial Nose-Hoover chain variables
    0: read from configuration (default)
    1: reset variables to zero.


    * **~INITIALISE.NTINHB** (*bool*) – 0,1 controls generation of initial Nose-Hoover (chain) barostat variables
    0: read from strartup file (if applicable) (default)
    1: reset variables to zero


    * **~INITIALISE.NTISHI** (*bool*) – 0,1 controls initial setting for lattice shift vectors
    0: read from configuration (default)
    1: reset shifts to zero.


    * **~INITIALISE.NTIRTC** (*bool*) – 0,1 controls initial setting of positions and orientations for roto-translational constraints
    0: read from configuration (default)
    1: reset based on initial configuraion of startup file


    * **~INITIALISE.NTICOM** (*int*) – 0,1,2 controls initial removal of COM motion
    0: no initial system COM motion removal (default)
    1: initial COM translation is removed
    2: initial COM rotation is removed


    * **~INITIALISE.NTISTI** (*bool*) – 0,1 controls generation of stochastic integrals
    0: read stochastic integrals and IG from configuration (default)
    1: set stochastic integrals to zero and use IG from here.


    * **~INITIALISE.IG** (*int*) – random number generator seed


    * **~INITIALISE.TEMPI** (*float*) – initial temperature



#### IG(_: in_ )

#### NTICOM(_: in_ )

#### NTINHB(_: boo_ )

#### NTINHT(_: boo_ )

#### NTIRTC(_: boo_ )

#### NTISHI(_: boo_ )

#### NTISHK(_: in_ )

#### NTISTI(_: boo_ )

#### NTIVEL(_: boo_ )

#### TEMPI(_: floa_ )

#### \__init__(NTIVEL: bool = False, NTISHK: int = 0, NTINHT: bool = False, NTINHB: bool = False, NTISHI: bool = False, NTIRTC: bool = False, NTICOM: int = 0, NTISTI: bool = False, IG: int = 0, TEMPI: float = 0, content: Optional[List[str]] = None)
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


#### name(_ = 'INITIALISE_ )

### _class_ pygromos.files.blocks.imd_blocks.INNERLOOP(NTILM: int = 0, NTILS: int = 0, NGPUS: int = 0, content: Optional[List[str]] = None)
Bases: `pygromos.files.blocks.imd_blocks._generic_imd_block`

INNERLOOP block


* **Variables**

    
    * **~INNERLOOP.NTILM** (*int*) – 


    * **~INNERLOOP.NTILS** (*int*) – 


    * **~INNERLOOP.NGPUS** (*int*) – 



#### NGPUS(_: in_ )

#### NTILM(_: in_ )

#### NTILS(_: in_ )

#### name(_ = 'INNERLOOP_ )

### _class_ pygromos.files.blocks.imd_blocks.LAMBDA(NTIL: int = 0, NTLI: List[int] = [], NILG1: List[int] = [], NILG2: List[int] = [], ALI: List[float] = [], BLI: List[float] = [], CLI: List[float] = [], DLI: List[float] = [], ELI: List[float] = [], content: Optional[List[str]] = None)
Bases: `pygromos.files.blocks.imd_blocks._generic_imd_block`

LAMBDA block

    This block is controlling the perturbation for free energy calculations.


* **Variables**

    
    * **~LAMBDA.NTIL** (*int*) – 


    * **~LAMBDA.NTLI****(****1****)** (*int** or **List**[**int**]*) – 


    * **~LAMBDA.NILG1****(****1****)** (*int** or **List**[**int**]*) – 


    * **~LAMBDA.NILG2****(****1****)** (*int** or **List**[**int**]*) – 


    * **~LAMBDA.ALI****(****1****)** (*float** or **List**[**float**]*) – 


    * **~LAMBDA.BLI****(****1****)** (*float** or **List**[**float**]*) – 


    * **~LAMBDA.CLI****(****1****)** (*float** or **List**[**float**]*) – 


    * **~LAMBDA.DLI****(****1****)** (*float** or **List**[**float**]*) – 


    * **~LAMBDA.ELI****(****1****)** (*float** or **List**[**float**]*) – 



#### ALI(_: List[float_ )

#### BLI(_: List[float_ )

#### CLI(_: List[float_ )

#### DLI(_: List[float_ )

#### ELI(_: List[float_ )

#### NILG1(_: List[int_ )

#### NILG2(_: List[int_ )

#### NTIL(_: in_ )

#### NTLI(_: List[int_ )

#### name(_: st_ _ = 'LAMBDA_ )

### _class_ pygromos.files.blocks.imd_blocks.MULTIBATH(ALGORITHM: int = 0, NBATHS: int = 0, TEMP0: List[float] = [], TAU: List[float] = [], DOFSET: int = 0, LAST: List[int] = [], COMBATH: List[int] = [], IRBATH: List[int] = [], NUM: Optional[int] = None, content: Optional[List[str]] = None)
Bases: `pygromos.files.blocks.imd_blocks._generic_imd_block`

MULTIBATH Block

> This block controls the Thermostat of a simulation. Multiple temperature baths can be coupled.


* **Variables**

    
    * **~MULTIBATH.ALGORITHM** (*int*) – temperature coupling algorithm
    weak-coupling(0)
    nose-hoover(1)
    nose-hoover-chains(2)   num
    (where num is the number of chains to use)


    * **~MULTIBATH.NUM** (*int**, **optional*) – Mumber of chains in Nosé Hoover chains scheme [only specify when needed]


    * **~MULTIBATH.NBATHS** (*int*) – Number of temperature baths


    * **~MULTIBATH.TEMP0** (*List**[**float**]*) – temperature of each bath (list len == NBATHS)


    * **~MULTIBATH.TAU** (*List**[**float**]*) – coupling of the temperaturebath to the system  (list len == NBATHS)


    * **~MULTIBATH.DOFSET** (*int*) – Number of set of Degrees of freedom.


    * **~MULTIBATH.LAST** (*List**[**int**]   **(**list len == DOFSET**)*) – Last atoms of each DOFSet


    * **~MULTIBATH.COMBATH** (*List**[**int**]   **(**list len == DOFSET**)*) – Index of the temperature baths


    * **~MULTIBATH.IRBATH** (*List**[**int**]   **(**list len == DOFSET**)*) – IRBAHT?



#### ALGORITHM(_: in_ )

#### COMBATH(_: List[int_ )

#### DOFSET(_: in_ )

#### IRBATH(_: List[int_ )

#### LAST(_: List[int_ )

#### NBATHS(_: in_ )

#### NUM(_: in_ )

#### TAU(_: List[float_ )

#### TEMP0(_: List[float_ )

#### adapt_multibath(last_atoms_bath: Dict[int, int], algorithm: Optional[int] = None, num: Optional[int] = None, T: Optional[Union[float, List[float]]] = None, TAU: Optional[float] = None)
This function is adding each atom set into a single multibath.

    #TODO implementation not correct with com_bath and irbath! Works for super simple cases though


* **Parameters**

    
    * **last_atoms_bath** (*Dict[int,int]*) – lastatom:bath


    * **algorithm** (*int*) – int code for algorihtm


    * **T** (*float,List[float], optional*) – temperature value


    * **TAU** (*float, optional*) – coupling var



* **Return type**

    None



#### block_to_string()

#### name(_: st_ _ = 'MULTIBATH_ )

#### read_content_from_str(content: List[str])

### _class_ pygromos.files.blocks.imd_blocks.NEW_REPLICA_EDS(REEDS: bool = False, NRES: int = 0, NUMSTATES: int = 0, NEOFF: int = 0, RES: List[float] = [], EIR: List[List[float]] = [[]], NRETRIAL: int = 0, NREQUIL: int = 0, EDS_STAT_OUT: int = 0, CONT: bool = False, PERIODIC: int = 0, content: Optional[List[str]] = None)
Bases: `pygromos.files.blocks.imd_blocks._generic_imd_block`

REPLICA_EDS Block

> This block is controlling the REPLICA_EDS settings in  gromos and is basically a mixture of EDS and RE block. (Don’t use them when using this block!)


* **Variables**

    
    * **~NEW_REPLICA_EDS.REEDS** (*bool*) – Shall REEDS be activated?


    * **~NEW_REPLICA_EDS.NRES** (*int*) – Number of s-Values


    * **~NEW_REPLICA_EDS.NEOFF** (*int*) – number of replica exchange eds energy Offset vectors (only used for 2D REEDS - still needs to be present ;))


    * **~NEW_REPLICA_EDS.NUMSTATES** (*int*) – Number of EDS-states


    * **~NEW_REPLICA_EDS.RES** (*List**[**float**]*) – s_values for all replicas


    * **~NEW_REPLICA_EDS.EIR** (*List**[**List**[**float**]**]*) – energy offsets for all replicas and all states  List[List[float]] = REPLICA[EDS_STATE[EIR]]


    * **~NEW_REPLICA_EDS.NERTRIAL** (*int*) – How many replica exchanges trials should be executed? (NRETRIAL\*STEP.NSTLIM == total simulation time)


    * **~NEW_REPLICA_EDS.NREQUIL** (*int*) – How many equilibration runs shall be exectured? (NREQUIL\*STEP.NSTLIM == total simulation time)


    * **~NEW_REPLICA_EDS.EDS_STAT_OUT** (*int*) – Shall the replica exchange information be outputted? (__future__ frequency of output.)


    * **~NEW_REPLICA_EDS.CONT** (*bool*) – Is this a continuation run?



#### CONT(_: boo_ )

#### EDS_STAT_OUT(_: in_ )

#### EIR(_: List[float_ )

#### NEOFF(_: in_ )

#### NREQUIL(_: in_ )

#### NRES(_: in_ )

#### NRETRIAL(_: in_ )

#### NUMSTATES(_: in_ )

#### PERIODIC(_: in_ )

#### REEDS(_: boo_ )

#### RES(_: List[float_ )

#### name(_: st_ _ = 'REPLICA_EDS_ )

#### read_content_from_str(content: List[str])

### _class_ pygromos.files.blocks.imd_blocks.NONBONDED(NLRELE: int = 0, APPAK: float = 0, RCRF: float = 0, EPSRF: float = 0, NSLFEXCL: bool = False, NSHAPE: int = 0, ASHAPE: float = 0, NA2CLC: int = 0, TOLA2: str = 0, EPSLS: float = 0, NKX: int = 0, NKY: int = 0, NKZ: int = 0, KCUT: float = 0, NGX: int = 0, NGY: int = 0, NGZ: int = 0, NASORD: int = 0, NFDORD: int = 0, NALIAS: int = 0, NSPORD: int = 0, NQEVAL: int = 0, FACCUR: float = 0, NRDGRD: bool = False, NWRGRD: bool = False, NLRLJ: bool = False, SLVDNS: float = 0, content: Optional[List[str]] = None)
Bases: `pygromos.files.blocks.imd_blocks._generic_imd_block`

NONBONDED block

> This block is controlling the Nonbonded term evaluation


* **Variables**

    
    * **~NONBONDED.NLRELE** (*int*) – 1-3 method to handle electrostatic interactions
    -1 : reaction-field (LSERF compatibility mode)
    0 : no electrostatic interactions
    1 : reaction-field
    2 : Ewald method
    3 : P3M method


    * **~NONBONDED.APPAK** (*float*) – >= 0.0 reaction-field inverse Debye screening length


    * **~NONBONDED.RCRF** (*float*) – >= 0.0 reaction-field radius


    * **~NONBONDED.EPSRF** (*float*) – = 0.0 || > 1.0 reaction-field permittivity


    * **~NONBONDED.NSLFEXCL** (*bool*) – contribution of excluded atoms to reaction field false=off/true=on


    * **~NONBONDED.NSHAPE** (*float*) – -1..10 lattice sum charge-shaping function
    -1 : gaussian
    0..10 : polynomial


    * **~NONBONDED.ASHAPE** (*float*) – > 0.0 width of the lattice sum charge-shaping function


    * **~NONBONDED.NA2CLC** (*int*) – 0..4 controls evaluation of lattice sum A2 term
    0 : A2 = A2~ = 0
    1 : A2~ exact, A2 = A2~
    2 : A2 numerical, A2~ = A2
    3 : A2~ exact from Ewald or from mesh and atom coords, A2 numerical
    4 : A2~ averaged from mesh only, A2 numerical


    * **~NONBONDED.TOLA2** (*float*) – > 0.0 tolerance for numerical A2 evaluation


    * **~NONBONDED.EPSLS** (*float*) – = 0.0 || > 1.0 lattice sum permittivity (0.0 = tinfoil)


    * **NKZ** (*NKX**, **NKY**,*) – > 0 maximum absolute Ewald k-vector components


    * **~NONBONDED.KCUT** (*float*) – > 0.0 Ewald k-space cutoff


    * **NGZ** (*NGX**, **NGY**,*) – > 0 P3M number of grid points


    * **~NONBONDED.NASORD** (*int*) – 1..5 order of mesh charge assignment function


    * **~NONBONDED.NFDORD** (*int*) – 0..5 order of the mesh finite difference operator
    0 : ik - differentiation
    1..5 : finite differentiation


    * **~NONBONDED.NALIAS** (*float*) – > 0 number of mesh alias vectors considered


    * **~NONBONDED.NSPORD** (*float*) – order of SPME B-spline functions (not available)


    * **~NONBONDED.NQEVAL** (*float*) – >= 0 controls accuracy reevaluation
    0 : do not reevaluate
    > 0 : evaluate every NQEVAL steps


    * **~NONBONDED.FACCUR** (*float*) – > 0.0 rms force error threshold to recompute influence function


    * **~NONBONDED.NRDGRD** (*bool*) – 0,1 read influence function
    0 : calculate influence function at simulation start up
    1 : read influence function from file (not yet implemented)


    * **~NONBONDED.NWRGRD** (*bool*) – 0,1 write influence function
    0 : do not write
    1 : write at the end of the simulation (not yet implemented)


    * **~NONBONDED.NLRLJ** (*bool*) – 0,1 controls long-range Lennard-Jones corrections
    0 : no corrections
    1 : do corrections (not yet implemented)


    * **~NONBONDED.SLVDNS** (*float*) – > 0.0 average solvent density for long-range LJ correction (ignored)



#### APPAK(_: floa_ )

#### ASHAPE(_: floa_ )

#### EPSLS(_: floa_ )

#### EPSRF(_: floa_ )

#### FACCUR(_: floa_ )

#### KCUT(_: floa_ )

#### NA2CLC(_: in_ )

#### NALIAS(_: in_ )

#### NASORD(_: in_ )

#### NFDORD(_: in_ )

#### NGX(_: in_ )

#### NGY(_: in_ )

#### NGZ(_: in_ )

#### NKX(_: in_ )

#### NKY(_: in_ )

#### NKZ(_: in_ )

#### NLRELE(_: in_ )

#### NLRLJ(_: boo_ )

#### NQEVAL(_: in_ )

#### NRDGRD(_: boo_ )

#### NSHAPE(_: in_ )

#### NSLFEXCL(_: boo_ )

#### NSPORD(_: in_ )

#### NWRGRD(_: boo_ )

#### RCRF(_: floa_ )

#### SLVDNS(_: floa_ )

#### TOLA2(_: st_ )

#### name(_: st_ _ = 'NONBONDED_ )

### _class_ pygromos.files.blocks.imd_blocks.PAIRLIST(ALGORITHM: int = 0, NSNB: int = 0, RCUTP: float = 0, RCUTL: float = 0, SIZE: Union[str, float] = 0, TYPE: Union[str, bool] = False, content: Optional[List[str]] = None)
Bases: `pygromos.files.blocks.imd_blocks._generic_imd_block`

PAIRLIST Block

> This block is controlling the pairlist control.


* **Variables**

    
    * **~PAIRLIST.ALGORITHM** (*int*) – standard(0) (gromos96 like pairlist)
    grid(1) (md++ grid pairlist)
    grid_cell(2) (creates a mask)


    * **~PAIRLIST.NSNB** (*int*) – frequency (number of steps) a pairlist is constructed


    * **~PAIRLIST.RCUTPL** (*float*) – short-range cut-off in twin-range


    * **~PAIRLIST.RCUTL** (*float*) – intermediate-range cut-off in twin-range


    * **~PAIRLIST.SIZE** (*str**, **float*) – grid cell size (or auto = 0.5 \* RCUTP)


    * **~PAIRLIST.TYPE** (*str**, **bool*) – chargegoup(0) (chargegroup based cutoff)
    atomic(1)     (atom based cutoff)



#### ALGORITHM(_: in_ )

#### NSNB(_: in_ )

#### RCUTL(_: floa_ )

#### RCUTPL(_: floa_ )

#### SIZE(_: Union[str, float_ )

#### TYPE(_: Union[str, bool_ )

#### \__init__(ALGORITHM: int = 0, NSNB: int = 0, RCUTP: float = 0, RCUTL: float = 0, SIZE: Union[str, float] = 0, TYPE: Union[str, bool] = False, content: Optional[List[str]] = None)
Args:

    ALGORITHM:
    NSNB:
    RCUTP:
    RCUTL:
    SIZE:
    TYPE:


#### name(_: st_ _ = 'PAIRLIST_ )

### _class_ pygromos.files.blocks.imd_blocks.PERTURBATION(NTG: int = 0, NRDGL: int = 0, RLAM: float = 0, DLAMT: float = 0, ALPHLJ: float = 0, ALPHC: float = 0, NLAM: int = 0, NSCALE: int = 0, content: Optional[List[str]] = None)
Bases: `pygromos.files.blocks.imd_blocks._generic_imd_block`

Pertubation Block

> This block is for Thermodynamic integration


* **Variables**

    
    * **~PERTURBATION.NTG** (*int*) – 0..1 controls use of free-energy calculation.
    0: no free-energy calculation (default)
    1: calculate dH/dRLAM


    * **~PERTURBATION.NRDGL** (*int*) – 0,1 controls reading of initial value for RLAM.
    0: use initial RLAM parameter from PERTURBATION block
    1: read from configuration


    * **~PERTURBATION.RLAM** (*float*) – 0.0..1.0 initial value for lambda


    * **~PERTURBATION.DLAMT** (*float*) – >= 0.0 rate of lambda increase in time.


    * **~PERTURBATION.ALPHLJ** (*float*) – >= 0.0 Lennard-Jones soft-core parameter


    * **~PERTURBATION.ALPHC** (*float*) – >= 0.0 Coulomb-RF soft-core parameter


    * **~PERTURBATION.NLAM** (*int*) – > 0 power dependence of lambda coupling


    * **~PERTURBATION.NSCALE** (*int*) – 0..2 controls use of interaction scaling
    0: no interaction scaling
    1: interaction scaling
    2: perturbation for all atom pairs with scaled

    > interactions. No perturbation for others.



    * **~PERTURBATION.ExampleBlock** – 


    * **~PERTURBATION.____________** – 


    * **DLAMT** (*#     NTG   NRDGL    RLAM*) – 1   0         0.0     0.0


    * **NSCALE** (*#  ALPHLJ   ALPHC    NLAM*) – 0.5     0.5         2       0



#### ALPHC(_: floa_ )

#### ALPHLJ(_: floa_ )

#### DLAMT(_: floa_ )

#### NLAM(_: in_ )

#### NRDGL(_: in_ )

#### NSCALE(_: in_ )

#### NTG(_: in_ )

#### RLAM(_: floa_ )

#### name(_: st_ _ = 'PERTURBATION_ )

### _class_ pygromos.files.blocks.imd_blocks.POSITIONRES(NTPOR: int = 0, NTPORB: int = 0, NTPORS: int = 0, CPOR: int = 0, content: Optional[List[str]] = None)
Bases: `pygromos.files.blocks.imd_blocks._generic_imd_block`

POSITIONRES block

> This block allows the managment of the Position Restraints during the Simulation.


* **Variables**

    
    * **~POSITIONRES.NTPOR** (*int*) – 0..3 controls atom positions re(con)straining.
    0: no position re(con)straints (default)
    1: restraining with force constant CPOR
    2: restraining with force constant CPOR weighted by atomic B-factors
    3: constraining


    * **~POSITIONRES.NTPORB** (*int*) – 0,1 controls reading of reference positions and B-factors
    0: read reference positions from startup file.
    1: read reference positions and B-factors from special file


    * **~POSITIONRES.NTPORS** (*int*) – 0,1 controls scaling of reference positions upon pressure scaling
    0: do not scale reference positions
    1: scale reference positions


    * **~POSITIONRES.CPOR** (*int*) – >= 0 position restraining force constant



#### CPOR(_: in_ )

#### NTPOR(_: in_ )

#### NTPORB(_: in_ )

#### NTPORS(_: in_ )

#### name(_ = 'POSITIONRES_ )

### _class_ pygromos.files.blocks.imd_blocks.PRECALCLAM(NRLAM: int = 0, MINLAM: float = 0, MAXLAM: float = 0, content: Optional[List[str]] = None)
Bases: `pygromos.files.blocks.imd_blocks._generic_imd_block`


* **Variables**

    **~PRECALCLAM.NTG** (*int*) – 



#### MAXLAM(_: floa_ )

#### MINLAM(_: floa_ )

#### NRLAM(_: in_ )

#### \__init__(NRLAM: int = 0, MINLAM: float = 0, MAXLAM: float = 0, content: Optional[List[str]] = None)
> Can be used to caluclate multiple extra lambda values


* **Parameters**

    
    * **NRLAM** (*int*) – 0  : off
    >1 : precalculating energies for NRLAM extra lambda values


    * **MINLAM** – between 0 and 1: minimum lambda value to precalculate energies


    * **MAXLAM** – between MINLAM and 1: maximum lambda value to precalculate energies



#### name(_ = 'PRECALCLAM_ )

### _class_ pygromos.files.blocks.imd_blocks.PRESSURESCALE(COUPLE: int = 0, SCALE: int = 0, COMP: float = 0, TAUP: float = 0, VIRIAL: int = 0, SEMIANISOTROPIC: List[int] = [], PRES0: List[List[float]] = [[0, 0, 0], [0, 0, 0], [0, 0, 0]], content: Optional[List[str]] = None)
Bases: `pygromos.files.blocks.imd_blocks._generic_imd_block`

PRESSURESCALE Block
This block controls the barostat of the simulation


* **Variables**

    
    * **~PRESSURESCALE.COUPLE** (*int*) – off(0), calc(1), scale(2)


    * **~PRESSURESCALE.SCALE** (*int*) – off(0), iso(1), aniso(2), full(3), semianiso(4)


    * **~PRESSURESCALE.COMP** (*float*) – compessibility


    * **~PRESSURESCALE.TAUP** (*float*) – coupling strength


    * **~PRESSURESCALE.VIRIAL** (*int*) – none(0), atomic(1), group(2)


    * **~PRESSURESCALE.SEMIANISOTROPIC** (*List**[**int**]*) – (semianisotropic couplings: X, Y, Z)
    e.g. 1 1 2: x and y jointly coupled and z separately coupled
    e.g. 0 0 1: constant area (xy-plane) and z coupled to a bath


    * **~PRESSURESCALE.PRES0** (*List**[**List**[**float**]**]*) – reference pressure



#### COMP(_: floa_ )

#### COUPLE(_: in_ )

#### PRES0(_: List[List[float]_ )

#### SCALE(_: in_ )

#### SEMIANISOTROPIC(_: List[int_ )

#### TAUP(_: floa_ )

#### VIRIAL(_: in_ )

#### name(_: st_ _ = 'PRESSURESCALE_ )

#### read_content_from_str(content: List[str])

### _class_ pygromos.files.blocks.imd_blocks.PRINTOUT(NTPR: int = 0, NTPP: int = 0, content: Optional[List[str]] = None)
Bases: `pygromos.files.blocks.imd_blocks._generic_imd_block`

PRINTOUT block

> This Block manages the output frequency into the .omd/std-out file.


* **Variables**

    
    * **~PRINTOUT.NTPR** (*int*) – print out energies, etc. every NTPR steps


    * **~PRINTOUT.NTPP** (*int*) – =1 perform dihedral angle transition monitoring



#### NTPP(_: in_ )

#### NTPR(_: in_ )

#### name(_: st_ _ = 'PRINTOUT_ )

### _class_ pygromos.files.blocks.imd_blocks.QMMM(NTQMMM: int = 0, NTQMSW: int = 0, RCUTQ: float = 0.0, NTWQMMM: int = 0, QMLJ: int = 0, MMSCAL: float = - 1.0, content: Optional[List[str]] = None)
Bases: `pygromos.files.blocks.imd_blocks._generic_imd_block`

QMMM Block

> The QMMM block describes how to treat an embedded QM zone


* **Variables**

    
    * **~QMMM.NTQMMM** (*int*) – Define if QM/MM is applied (1) or not (0, default)


    * **~QMMM.NTQMSW** (*int*) – QM software package to use. Available options are MNDO (0, default), Turbomole (1), DFTB (2), MOPAC (3), ORCA (4), and XTB (5)


    * **~QMMM.RCUTQ** (*float*) – Cutoff for inclusion of MM charge groups. If the value is 0, all particles are considered.


    * **~QMMM.NTWQMMM** (*int*) – Write QM/MM related data to special trajectory every NTWQMMM step. Switch off with 0.


    * **~QMMM.QMLJ** (*int*) – Define if LJ interaction is activated in the QM zone (1) or not (0).


    * **~QMMM.MMSCAL** (*float*) – Scale MM charges with (2/pi)\*atan(x\*(r_{mm}-r_{mm})). Values > 0.0 describe the scaling factor, values < 0.0 turn off scaling (default).



#### MMSCAL(_: floa_ )

#### NTQMMM(_: in_ )

#### NTQMSW(_: in_ )

#### NTWQMMM(_: in_ )

#### QMLJ(_: in_ )

#### RCUTQ(_: floa_ )

#### name(_: st_ _ = 'QMMM_ )

### _class_ pygromos.files.blocks.imd_blocks.RANDOMNUMBERS(NTRNG: int = 0, NTGSL: int = 0, content: Optional[List[str]] = None)
Bases: `pygromos.files.blocks.imd_blocks._generic_imd_block`

Random Numbers Block


#### NTGSL(_: in_ )

#### NTRNG(_: in_ )

#### name(_: st_ _ = 'RANDOMNUMBERS_ )

### _class_ pygromos.files.blocks.imd_blocks.REPLICA(RETL: bool = False, NRET: int = 0, RET: List[float] = [], LRESCALE: int = 0, NRELAM: int = 0, RELAM: List[float] = [], RETS: List[float] = [], NRETRIAL: int = 0, NREQUIL: int = 0, CONT: bool = False, content: Optional[List[str]] = None)
Bases: `pygromos.files.blocks.imd_blocks._generic_imd_block`

REPLICA Block

> This block is controlling the REPLICA settings in gromos. It works both for T-REMD and H-REMD


* **Variables**

    
    * **~REPLICA.RETL** (*bool*) – Do you want to rund an RE simulation or not?
    RETL=0 : turn off RE
    RETL=1 : turn on RE


    * **~REPLICA.NRET** (*int*) – Number of temperatures


    * **~REPLICA.RET** (*List**[**float**]*) – List of temperature at each replica


    * **~REPLICA.LRESCALE** (*bool*) – controls temperature scaling
    0 don’t scale temperatures after exchange trial
    1 scale temperatures after exchange trial


    * **~REPLICA.NRELAM** (*int*) – number of lambda values


    * **~REPLICA.RELAM** (*List**[**float**]*) – lambda value for each lambda-replica


    * **~REPLICA.RETS** (*List**[**float**]*) – timestep of each LAMBDA replica


    * **~REPLICA.NRETRIAL** (*int*) – How many replica exchanges trials should be executed? (NRETRIAL\*STEP.NSTLIM == total simulation time)


    * **~REPLICA.NREQUIL** (*int*) – How many equilibration runs shall be exectured? (NREQUIL\*STEP.NSTLIM == total simulation time)


    * **~REPLICA.CONT** (*bool*) – Is this a continuation run?
    0 start from one configuration file
    1 start from multiple configuration files


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


#### CONT(_: boo_ )

#### LRESCALE(_: boo_ )

#### NRELAM(_: in_ )

#### NREQUIL(_: in_ )

#### NRET(_: in_ )

#### NRETRIAL(_: in_ )

#### RELAM(_: List[float_ )

#### RET(_: List[float_ )

#### RETL(_: boo_ )

#### RETS(_: List[float_ )

#### name(_: st_ _ = 'REPLICA_ )

#### read_content_from_str(content: List[str])

### _class_ pygromos.files.blocks.imd_blocks.REPLICA_EDS(REEDS: bool = True, NRES: int = 0, NUMSTATES: int = 0, RES: List[float] = [0], EIR: List[List[float]] = [[0]], NRETRIAL: int = 0, NREQUIL: int = 0, EDS_STAT_OUT: int = 0, CONT: bool = True, content: Optional[List[str]] = None)
Bases: `pygromos.files.blocks.imd_blocks._generic_imd_block`


#### CONT(_: boo_ )

#### EDS_STAT_OUT(_: in_ )

#### EIR(_: List[float_ )

#### NREQUIL(_: in_ )

#### NRES(_: in_ )

#### NRETRIAL(_: in_ )

#### NUMSTATES(_: in_ )

#### REEDS(_: boo_ )

#### RES(_: List[float_ )

#### \__init__(REEDS: bool = True, NRES: int = 0, NUMSTATES: int = 0, RES: List[float] = [0], EIR: List[List[float]] = [[0]], NRETRIAL: int = 0, NREQUIL: int = 0, EDS_STAT_OUT: int = 0, CONT: bool = True, content: Optional[List[str]] = None)
REPLICA_EDS Block

> This block is controlling the REPLICA_EDS settings in  gromos and is basically a mixture of EDS and RE block. (Don’t use them when using this block!)


* **Variables**

    
    * **~REPLICA_EDS.__init__.REEDS** (*bool*) – Shall REEDS be activated?


    * **~REPLICA_EDS.__init__.NRES** (*int*) – Number of s-Values


    * **~REPLICA_EDS.__init__.NUMSTATES** (*int*) – Number of EDS-states


    * **~REPLICA_EDS.__init__.RES** (*List**[**float**]*) – s_values for all replicas


    * **~REPLICA_EDS.__init__.EIR** (*List**[**List**[**float**]**]*) – energy offsets for all replicas and all states  List[List[float]] = REPLICA[EDS_STATE[EIR]]


    * **~REPLICA_EDS.__init__.NERTRIAL** (*int*) – How many replica exchanges trials should be executed? (NRETRIAL\*STEP.NSTLIM == total simulation time)


    * **~REPLICA_EDS.__init__.NREQUIL** (*int*) – How many equilibration runs shall be exectured? (NREQUIL\*STEP.NSTLIM == total simulation time)


    * **~REPLICA_EDS.__init__.EDS_STAT_OUT** (*int*) – Shall the replica exchange information be outputted? (__future__ frequency of output.)


    * **~REPLICA_EDS.__init__.CONT** (*bool*) – Is this a continuation run?



#### name(_: st_ _ = 'REPLICA_EDS_ )

#### read_content_from_str(content: List[str])

### _class_ pygromos.files.blocks.imd_blocks.ROTTRANS(RTC: int = 0, RTCLAST: int = 0, content: Optional[List[str]] = None)
Bases: `pygromos.files.blocks.imd_blocks._generic_imd_block`

Roto-translational block

> This block is for roto-translational constraints.
> Note: use either centre of mass removal or roto-translational constraints but not both!


* **Variables**

    
    * **~ROTTRANS.RTC** (*int*) – turn roto-translational constraints on (1)


    * **~ROTTRANS.RTCLAST** (*int*) – last atom of subset to be roto-translationally constrained



#### RTC(_: in_ )

#### RTCLAST(_: in_ )

#### name(_: st_ _ = 'ROTTRANS_ )

### _class_ pygromos.files.blocks.imd_blocks.STEP(NSTLIM: int = 0, T: float = 0, DT: float = 0, content: Optional[List[str]] = None)
Bases: `pygromos.files.blocks.imd_blocks._generic_imd_block`

Step Block

> This Block gives the number of simulation steps,


* **Variables**

    
    * **~STEP.NSTLIM** (*int*) – number of simulations Step till terminating.


    * **~STEP.T** (*float*) – Starting Time


    * **~STEP.DT** (*float*) – time step [fs]



#### DT(_: floa_ )

#### NSTLIM(_: in_ )

#### T(_: floa_ )

#### name(_: st_ _ = 'STEP_ )

### _class_ pygromos.files.blocks.imd_blocks.STOCHDYN(NTSD: int = 0, NTFR: int = 0, NSFR: int = 0, NBREF: int = 0, RCUTF: float = 0, CFRIC: float = 0, TEMPSD: float = 0, content: Optional[List[str]] = None)
Bases: `pygromos.files.blocks.imd_blocks._generic_imd_block`

Stochastic Dynamics block

> This block is for Stochastic Dynamics


* **Variables**

    
    * **~STOCHDYN.NTSD** (*int*) – 


    * **~STOCHDYN.NTFR** (*int*) – 


    * **~STOCHDYN.NSFR** (*int*) – 


    * **~STOCHDYN.NBREF** (*int*) – 


    * **~STOCHDYN.RCUTF** (*float*) – 


    * **~STOCHDYN.CFRIC** (*float*) – 


    * **~STOCHDYN.TEMPSD** (*float*) – 



#### CFRIC(_: floa_ )

#### NBREF(_: in_ )

#### NSFR(_: in_ )

#### NTFR(_: in_ )

#### NTSD(_: in_ )

#### RCUTF(_: floa_ )

#### TEMPSD(_: floa_ )

#### name(_: st_ _ = 'STOCHDYN_ )

### _class_ pygromos.files.blocks.imd_blocks.SYSTEM(NPM: int = 0, NSM: int = 0, content: Optional[List[str]] = None)
Bases: `pygromos.files.blocks.imd_blocks._generic_imd_block`

System Block

> The system block defines the number of solute molecules and solvent molecules


* **Variables**

    
    * **~SYSTEM.NPM** (*int*) – Number of Solute Molecules


    * **~SYSTEM.NSM** (*int*) – Number of Solvent Molecules



#### NPM(_: in_ )

#### NSM(_: in_ )

#### name(_: st_ _ = 'SYSTEM_ )

### _class_ pygromos.files.blocks.imd_blocks.WRITETRAJ(NTWX: int = 0, NTWSE: int = 0, NTWV: int = 0, NTWF: int = 0, NTWE: int = 0, NTWG: int = 0, NTWB: int = 0, content: Optional[List[str]] = None)
Bases: `pygromos.files.blocks.imd_blocks._generic_imd_block`

> The WRITETRAJ Block manages the output frequency of GROMOS. Here you can determine how often a file should be written.


* **Variables**

    
    * **~WRITETRAJ.NTWX** (*int*) – controls writing of coordinate trajectory
    0: no coordinate trajectory is written (default)
    >0: write solute and solvent coordinates every NTWX steps
    <0: write solute coordinates every 

    ```
    |NTWX|
    ```

     steps


    * **~WRITETRAJ.NTWSE** (*int*) – >= 0 selection criteria for coordinate trajectory writing
    0: write normal coordinate trajectory
    >0: write minimum-energy coordinate and energy trajectory (based on the
    energy entry selected by NTWSE and as blocks of length NTWX)
    (see configuration/energy.cc or ene_ana library for indices)


    * **~WRITETRAJ.NTWV** (*int*) – controls writing of velocity trajectory
    0: no velocity trajectory is written (default)
    >0: write solute and solvent velocities every NTWV steps
    <0: write solute velocities every 

    ```
    |NTWV|
    ```

     steps


    * **~WRITETRAJ.NTWF** (*int*) – controls writing of force trajectory
    0: no force trajectory is written (default)
    >0: write solute and solvent forces every NTWF steps
    <0: write solute forces every 

    ```
    |NTWF|
    ```

     steps


    * **~WRITETRAJ.NTWE** (*int*) – >= 0 controls writing of energy trajectory
    0: no energy trajectory is written (default)
    >0: write energy trajectory every NTWE steps


    * **~WRITETRAJ.NTWG** (*int*) – >= 0 controls writing of free energy trajectory
    0: no free energy trajectory is written (default)
    >0: write free energy trajectory every NTWG steps


    * **~WRITETRAJ.NTWB** (*int*) – >= 0 controls writing of block-averaged energy trajectory
    0: no block averaged energy trajectory is written (default)
    >0: write block-averaged energy variables every 

    ```
    |NTWB|
    ```

     steps

    > (and free energies if NTWG > 0) trajectory




#### NTWB(_: in_ )

#### NTWE(_: in_ )

#### NTWF(_: in_ )

#### NTWG(_: in_ )

#### NTWSE(_: in_ )

#### NTWV(_: in_ )

#### NTWX(_: in_ )

#### \__init__(NTWX: int = 0, NTWSE: int = 0, NTWV: int = 0, NTWF: int = 0, NTWE: int = 0, NTWG: int = 0, NTWB: int = 0, content: Optional[List[str]] = None)
Args:

    NTWX:
    NTWSE:
    NTWV:
    NTWF:
    NTWE:
    NTWG:
    NTWB:


#### name(_ = 'WRITETRAJ_ )
## pygromos.files.blocks.miscBlocks module


### _class_ pygromos.files.blocks.miscBlocks.ATOMNAMELIB(content: Union[List[str], Dict[str, Dict[str, str]]])
Bases: `pygromos.files.blocks._general_blocks._generic_gromos_block`


#### \__init__(content: Union[List[str], Dict[str, Dict[str, str]]])
> atom name translation lib part of the resn lib


* **Parameters**

    **content** (*Union[List[str], Dict[str, Dict[str, str]]]*)



#### block_to_string()

#### fields(_ = ['top_resn', 'atom_pdb', 'atom_top'_ )

#### pdb_top(_: defaultdict(<class 'dict'>, {}_ )

#### read_content_from_str(content: Union[List[str], Dict[str, Dict[str, str]]])

### _class_ pygromos.files.blocks.miscBlocks.RESIDUENAMELIB(content: Union[str, List[str], Dict[str, Dict[str, str]]])
Bases: `pygromos.files.blocks._general_blocks._generic_gromos_block`


#### \__init__(content: Union[str, List[str], Dict[str, Dict[str, str]]])
content is the content for the  residuenamelib


* **Parameters**

    **content** (*Union[List[str], Dict[str, Dict[str, str]]]*)



#### block_to_string()

#### fields(_ = ['pdb_name', 'g96top_name'_ )

#### pdb_top(_: defaultdict(<class 'list'>, {}_ )

#### read_content_from_str(content: List[str])
## pygromos.files.blocks.mtb_blocks module


### _class_ pygromos.files.blocks.mtb_blocks.LINKEXCLUSIONS(FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None, content=None)
Bases: `pygromos.files.blocks.mtb_blocks.mtb_blocks`


#### FORCEFIELD(_: pygromos.files.blocks.topology_blocks.FORCEFIEL_ )

#### MAKETOPVERSION(_: pygromos.files.blocks.topology_blocks.MAKETOPVERSIO_ )

#### NRNE(_: in_ )

#### block_to_string()

#### read_content_from_str(content: str)

### _class_ pygromos.files.blocks.mtb_blocks.MTBUILDBLEND(FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None, content=None)
Bases: `pygromos.files.blocks.mtb_blocks.mtb_blocks`


#### FORCEFIELD(_: pygromos.files.blocks.topology_blocks.FORCEFIEL_ )

#### MAKETOPVERSION(_: pygromos.files.blocks.topology_blocks.MAKETOPVERSIO_ )

#### angles(_: List[pygromos.files.blocks.mtb_blocks.mtb_angles_field_ )

#### atoms(_: List[pygromos.files.blocks.mtb_blocks.mtb_atoms_field_ )

#### block_to_string()

#### bonds(_: List[pygromos.files.blocks.mtb_blocks.mtb_bonds_field_ )

#### dihedrals(_: List[pygromos.files.blocks.mtb_blocks.mtb_dihedral_field_ )

#### improper_dihedrals(_: List[pygromos.files.blocks.mtb_blocks.mtb_dihedral_field_ )

#### lj_exceptions(_: List[pygromos.files.blocks.mtb_blocks.mtb_lj_exceptions_field_ )

#### read_content_from_str(content: str)

#### replacing_atoms(_: List[pygromos.files.blocks.mtb_blocks.mtb_trailing_atoms_field_ )

#### skip_lines(_: in_ _ = _ )

### _class_ pygromos.files.blocks.mtb_blocks.MTBUILDBLSOLUTE(FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None, content=None)
Bases: `pygromos.files.blocks.mtb_blocks.mtb_blocks`


#### FORCEFIELD(_: pygromos.files.blocks.topology_blocks.FORCEFIEL_ )

#### MAKETOPVERSION(_: pygromos.files.blocks.topology_blocks.MAKETOPVERSIO_ )

#### angles(_: List[pygromos.files.blocks.mtb_blocks.mtb_angles_field_ )

#### atoms(_: List[pygromos.files.blocks.mtb_blocks.mtb_atoms_field_ )

#### block_to_string()

#### bonds(_: List[pygromos.files.blocks.mtb_blocks.mtb_bonds_field_ )

#### dihedrals(_: List[pygromos.files.blocks.mtb_blocks.mtb_dihedral_field_ )

#### improper_dihedrals(_: List[pygromos.files.blocks.mtb_blocks.mtb_dihedral_field_ )

#### lj_exceptions(_: List[pygromos.files.blocks.mtb_blocks.mtb_lj_exceptions_field_ )

#### preceding_exclusions(_: List[pygromos.files.blocks.mtb_blocks.mtb_preceding_exclusions_field_ )

#### read_content_from_str(content: str)

#### skip_lines(_: in_ _ = _ )

#### trailing_atoms(_: List[pygromos.files.blocks.mtb_blocks.mtb_trailing_atoms_field_ )

### _class_ pygromos.files.blocks.mtb_blocks.MTBUILDBLSOLVENT(FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None, content=None)
Bases: `pygromos.files.blocks.mtb_blocks.mtb_blocks`


#### FORCEFIELD(_: pygromos.files.blocks.topology_blocks.FORCEFIEL_ )

#### MAKETOPVERSION(_: pygromos.files.blocks.topology_blocks.MAKETOPVERSIO_ )

#### atoms(_: List[pygromos.files.blocks.mtb_blocks.mtb_atoms_field_ )

#### block_to_string()

#### constraints(_: List[pygromos.files.blocks.mtb_blocks.mtb_constraints_field_ )

#### read_content_from_str(content: str)

### _class_ pygromos.files.blocks.mtb_blocks.mtb_angles_field(IB: int, JB: int, KB: int, MCB: int)
Bases: `pygromos.files.blocks.mtb_blocks.mtb_fields`


#### to_string()

### _class_ pygromos.files.blocks.mtb_blocks.mtb_atoms_field(ATOM: int, ANM: str, IACM: int, MASS: int, CGMI: float, CGM: int, MAE: int, MSAE: List[int])
Bases: `pygromos.files.blocks.mtb_blocks.mtb_fields`


#### to_string()

### _class_ pygromos.files.blocks.mtb_blocks.mtb_atoms_solvent_field(ATOM: int, ANM: str, IACM: int, MASS: int, CG: float)
Bases: `pygromos.files.blocks.mtb_blocks.mtb_fields`


#### to_string()

### _class_ pygromos.files.blocks.mtb_blocks.mtb_blocks(name: Optional[str] = None, used: Optional[bool] = None, content: Optional[str] = None)
Bases: `pygromos.files.blocks._general_blocks._generic_gromos_block`


#### comment(_: st_ )

#### content(_: Iterabl_ )

### _class_ pygromos.files.blocks.mtb_blocks.mtb_bonds_field(IB: int, JB: int, MCB: int)
Bases: `pygromos.files.blocks.mtb_blocks.mtb_fields`


#### to_string()

### _class_ pygromos.files.blocks.mtb_blocks.mtb_constraints_field(IB: int, JB: int, LENGTH: float)
Bases: `pygromos.files.blocks.mtb_blocks.mtb_fields`


#### to_string()

### _class_ pygromos.files.blocks.mtb_blocks.mtb_dihedral_field(IB: int, JB: int, KB: int, LB: int, MCB: int)
Bases: `pygromos.files.blocks.mtb_blocks.mtb_fields`


#### to_string()

### _class_ pygromos.files.blocks.mtb_blocks.mtb_fields()
Bases: `pygromos.files.blocks._general_blocks._generic_field`


#### field_continue_next_line(_ = '\\n\\t\\t\\t\\t\\t\\t\\t\\t\\t\\t_ )

#### fieldseperator(_ = '  _ )

#### lineseperator(_ = ' \\n_ )

### _class_ pygromos.files.blocks.mtb_blocks.mtb_lj_exceptions_field(iac: int, jac: int, mcb: int)
Bases: `pygromos.files.blocks.mtb_blocks.mtb_fields`


#### to_string()

### _class_ pygromos.files.blocks.mtb_blocks.mtb_preceding_exclusions_field(ATOM: int, MAE: int, MSAE: List[int])
Bases: `pygromos.files.blocks.mtb_blocks.mtb_fields`


#### to_string()

### _class_ pygromos.files.blocks.mtb_blocks.mtb_trailing_atoms_field(ATOM: int, ANM: str, IACM: int, MASS: int, CGMI: float, CGM: int)
Bases: `pygromos.files.blocks.mtb_blocks.mtb_fields`


#### to_string()
## pygromos.files.blocks.pertubation_blocks module


### _class_ pygromos.files.blocks.pertubation_blocks.MPERTATOM(NJLA: Optional[int] = None, NPTB: Optional[int] = None, STATEIDENTIFIERS: List[str] = [], STATEATOMHEADER: Tuple[str] = ['NR', 'NAME', 'ALPHLJ', 'ALPHCRF'], STATEATOMS: List[pygromos.files.blocks.pertubation_blocks.atom_eds_pertubation_state] = [], dummy_IAC: int = 22, dummy_CHARGE: float = 0.0, content: Optional[List[str]] = None)
Bases: `pygromos.files.blocks._general_blocks._generic_gromos_block`


#### \__init__(NJLA: Optional[int] = None, NPTB: Optional[int] = None, STATEIDENTIFIERS: List[str] = [], STATEATOMHEADER: Tuple[str] = ['NR', 'NAME', 'ALPHLJ', 'ALPHCRF'], STATEATOMS: List[pygromos.files.blocks.pertubation_blocks.atom_eds_pertubation_state] = [], dummy_IAC: int = 22, dummy_CHARGE: float = 0.0, content: Optional[List[str]] = None)
> This block is used for lambda sampling to define the different states.


* **Parameters**

    
    * **NJLA** (*int*) – number of perturbed atoms


    * **NPTB** (*int*) – number of pertubation states


    * **STATEIDENTIFIERS** (*List[str]*) – string names for states


    * **STATEATOMHEADER** – header for the atom description table


    * **STATEATOMS** – list of atoms, that shall be perturbed


    * **dummy_IAC** – dummy atom VdW type for perturbed atoms


    * **dummy_CHARGE** – dummy atom charge type for perturbed atoms



#### add_state_atoms(state_atoms: List[pygromos.files.blocks.pertubation_blocks.atom_eds_pertubation_state])
This function can add states and atoms, but also overwrite state values of existing atoms.
If a new state is defined only for a subset of atoms, all other atoms are set to the default dummy.
If a new atom misses a state definition, this state will be set to dummy.


* **Parameters**

    **state_atoms** (*List[atom_eds_pertubation_state]*)



#### block_to_string()

#### comment(_: st_ )

#### content(_: Iterabl_ )

#### delete_atom(atomNR: Union[int, List[int]])
This function removes atom lines from the ptp file.


* **Parameters**

    **atomNR** (*int*) – atom to be removed.



#### delete_state(stateIDs: Optional[Union[int, List[int]]] = None, stateNames: Optional[Union[str, List[str]]] = None)
This function deletes an state column.


* **Parameters**

    **stateIDs** (*int*) – number of the state



#### _property_ nStates(_: in_ )

#### _property_ nTotalStateAtoms(_: in_ )

#### read_content_from_str(content: List[str])

#### _property_ states(_: dic_ )

### _class_ pygromos.files.blocks.pertubation_blocks.PERTATOMPARAM(STATEATOMS: Optional[List[pygromos.files.blocks.pertubation_blocks.atom_lam_pertubation_state]] = None, STATEATOMHEADER: Optional[Tuple[str]] = None, NJLA: Optional[int] = None, STATEIDENTIFIERS=None, dummy_IAC=22, dummy_CHARGE=0.0, content: Optional[List[str]] = None)
Bases: `pygromos.files.blocks._general_blocks._generic_gromos_block`


#### add_state_atoms(state_atoms: List[pygromos.files.blocks.pertubation_blocks.atom_lam_pertubation_state])
This function can add states and atoms, but also overwrite state values of existing atoms.
If a new state is defined only for a subset of atoms, all other atoms are set to the default dummy.
If a new atom misses a state definition, this state will be set to dummy.


* **Parameters**

    **state_atoms** (*List[atom_eds_pertubation_state]*)



#### block_to_string()

#### comment(_: st_ )

#### content(_: Iterabl_ )

#### delete_atom(atomNR: Union[int, List[int]])
This function removes atom lines from the ptp file.


* **Parameters**

    **atomNR** (*int*) – atom to be removed.



#### delete_state(stateIDs: Optional[Union[int, List[int]]] = None, stateNames: Optional[Union[str, List[str]]] = None)
This function deletes an state column.


* **Parameters**

    **stateIDs** (*int*) – number of the state



#### _property_ nStates(_: in_ )

#### _property_ nTotalStateAtoms(_: in_ )

#### read_content_from_str(content: List[str])

#### _property_ states(_: Dic_ )

### _class_ pygromos.files.blocks.pertubation_blocks.PERTBONDANGLE(STATEATOMS: Optional[List[pygromos.files.blocks.pertubation_blocks.atom_lam_pertubation_state_angle]] = None, STATEATOMHEADER: Optional[Tuple[str]] = None, NPA: Optional[int] = None, dummy_ANGLE=22, content: Optional[List[str]] = None)
Bases: `pygromos.files.blocks._general_blocks._generic_gromos_block`


#### block_to_string()

#### comment(_: st_ )

#### content(_: Iterabl_ )

#### _property_ nStates(_: in_ )

#### _property_ nTotalStateAtoms(_: in_ )

#### read_content_from_str(content: List[str])

#### _property_ states(_: dic_ )

### _class_ pygromos.files.blocks.pertubation_blocks.PERTBONDANGLEH(STATEATOMS: Optional[List[pygromos.files.blocks.pertubation_blocks.atom_lam_pertubation_state_angle]] = None, STATEATOMHEADER: Optional[Tuple[str]] = None, NPA: Optional[int] = None, dummy_ANGLE=22, content: Optional[List[str]] = None)
Bases: `pygromos.files.blocks._general_blocks._generic_gromos_block`


#### block_to_string()

#### comment(_: st_ )

#### content(_: Iterabl_ )

#### _property_ nStates(_: in_ )

#### _property_ nTotalStateAtoms(_: in_ )

#### read_content_from_str(content: List[str])

#### _property_ states(_: dic_ )

### _class_ pygromos.files.blocks.pertubation_blocks.PERTBONDSTRETCH(STATEATOMS: Optional[List[pygromos.files.blocks.pertubation_blocks.atom_lam_pertubation_state_bond]] = None, STATEATOMHEADER: Optional[Tuple[str]] = None, NPB: Optional[int] = None, dummy_BOND=22, content: Optional[List[str]] = None)
Bases: `pygromos.files.blocks._general_blocks._generic_gromos_block`


#### block_to_string()

#### comment(_: st_ )

#### content(_: Iterabl_ )

#### _property_ nStates(_: in_ )

#### _property_ nTotalStateAtoms(_: in_ )

#### read_content_from_str(content: List[str])

#### _property_ states(_: dic_ )

### _class_ pygromos.files.blocks.pertubation_blocks.PERTBONDSTRETCHH(STATEATOMS: Optional[List[pygromos.files.blocks.pertubation_blocks.atom_lam_pertubation_state_bond]] = None, STATEATOMHEADER: Optional[Tuple[str]] = None, NPB: Optional[int] = None, dummy_BOND=22, content: Optional[List[str]] = None)
Bases: `pygromos.files.blocks._general_blocks._generic_gromos_block`


#### block_to_string()

#### comment(_: st_ )

#### content(_: Iterabl_ )

#### _property_ nStates(_: in_ )

#### _property_ nTotalStateAtoms(_: in_ )

#### read_content_from_str(content: List[str])

#### _property_ states(_: dic_ )

### _class_ pygromos.files.blocks.pertubation_blocks.PERTIMROPERDIH(STATEATOMS: Optional[List[pygromos.files.blocks.pertubation_blocks.atom_lam_pertubation_state_dihedral]] = None, STATEATOMHEADER: Optional[Tuple[str]] = None, NPD: Optional[int] = None, dummy_IMP=22, content: Optional[List[str]] = None)
Bases: `pygromos.files.blocks._general_blocks._generic_gromos_block`


#### block_to_string()

#### comment(_: st_ )

#### content(_: Iterabl_ )

#### _property_ nStates(_: in_ )

#### _property_ nTotalStateAtoms(_: in_ )

#### read_content_from_str(content: List[str])

#### _property_ states(_: dic_ )

### _class_ pygromos.files.blocks.pertubation_blocks.PERTIMROPERDIHH(STATEATOMS: Optional[List[pygromos.files.blocks.pertubation_blocks.atom_lam_pertubation_state_dihedral]] = None, STATEATOMHEADER: Optional[Tuple[str]] = None, NPD: Optional[int] = None, dummy_IMP=22, content: Optional[List[str]] = None)
Bases: `pygromos.files.blocks._general_blocks._generic_gromos_block`


#### block_to_string()

#### comment(_: st_ )

#### content(_: Iterabl_ )

#### _property_ nStates(_: in_ )

#### _property_ nTotalStateAtoms(_: in_ )

#### read_content_from_str(content: List[str])

#### _property_ states(_: dic_ )

### _class_ pygromos.files.blocks.pertubation_blocks.PERTPROPERDIH(STATEATOMS: Optional[List[pygromos.files.blocks.pertubation_blocks.atom_lam_pertubation_state_dihedral]] = None, STATEATOMHEADER: Optional[Tuple[str]] = None, NPD: Optional[int] = None, dummy_DIH=22, content: Optional[List[str]] = None)
Bases: `pygromos.files.blocks._general_blocks._generic_gromos_block`


#### block_to_string()

#### comment(_: st_ )

#### content(_: Iterabl_ )

#### _property_ nStates(_: in_ )

#### _property_ nTotalStateAtoms(_: in_ )

#### read_content_from_str(content: List[str])

#### _property_ states(_: dic_ )

### _class_ pygromos.files.blocks.pertubation_blocks.PERTPROPERDIHH(STATEATOMS: Optional[List[pygromos.files.blocks.pertubation_blocks.atom_lam_pertubation_state_dihedral]] = None, STATEATOMHEADER: Optional[Tuple[str]] = None, NPD: Optional[int] = None, dummy_DIH=22, content: Optional[List[str]] = None)
Bases: `pygromos.files.blocks._general_blocks._generic_gromos_block`


#### _class_ PERTPROPERDIH(STATEATOMS: Optional[List[pygromos.files.blocks.pertubation_blocks.atom_lam_pertubation_state_dihedral]] = None, STATEATOMHEADER: Optional[Tuple[str]] = None, NPD: Optional[int] = None, dummy_DIH=22, content: Optional[List[str]] = None)
Bases: `pygromos.files.blocks._general_blocks._generic_gromos_block`


#### block_to_string()

#### comment(_: st_ )

#### content(_: Iterabl_ )

#### _property_ nStates(_: in_ )

#### _property_ nTotalStateAtoms(_: in_ )

#### read_content_from_str(content: List[str])

#### _property_ states(_: dic_ )

#### block_to_string()

#### comment(_: st_ )

#### content(_: Iterabl_ )

#### _property_ nStates(_: in_ )

#### _property_ nTotalStateAtoms(_: in_ )

#### read_content_from_str(content: List[str])

#### _property_ states(_: dic_ )

### _class_ pygromos.files.blocks.pertubation_blocks.TITLE(content: str, field_seperator: str = '\\t', line_seperator: str = '\\n', name: str = 'TITLE', used: bool = True)
Bases: `pygromos.files.blocks._general_blocks._generic_gromos_block`


#### block_to_string()

#### content(_: st_ )

#### field_seperator(_: st_ _ = '\\t_ )

#### line_seperator(_: st_ _ = '\\n_ )

#### order(_ = [[['content']]_ )

#### pyGromosWatermark(_: st_ _ = '>>> Generated with PyGromosTools (riniker group) <<<_ )

#### read_content_from_str(content: List[str])

### _class_ pygromos.files.blocks.pertubation_blocks.atom_eds_pertubation_state(NR: int, NAME: str, STATES: Dict[int, __main__.pertubationEdsState], ALPHLJ: float = 1.0, ALPHCRF: float = 1.0)
Bases: `pygromos.files.blocks._general_blocks._generic_field`


#### state_format_pattern(_ = ' {:>3} {:>10.5f}_ )

#### to_string()

### _class_ pygromos.files.blocks.pertubation_blocks.atom_lam_pertubation_state(NR: int, RES: int, NAME: str, STATES: Dict[int, __main__.pertubationLamState], ALPHLJ: float = 1.0, ALPHCRF: float = 1.0)
Bases: `pygromos.files.blocks._general_blocks._generic_field`


#### state_format_pattern(_ = ' {:>5} {:>5} {:>10.5f}_ )

#### to_string()

### _class_ pygromos.files.blocks.pertubation_blocks.atom_lam_pertubation_state_angle(NR: int, atomI: int, atomJ: int, atomK: int, STATES: Dict[int, int])
Bases: `pygromos.files.blocks._general_blocks._generic_field`


#### state_format_pattern(_ = ' {:>5}_ )

#### to_string()

### _class_ pygromos.files.blocks.pertubation_blocks.atom_lam_pertubation_state_bond(NR: int, atomI: int, atomJ: int, STATES: Dict[int, int])
Bases: `pygromos.files.blocks._general_blocks._generic_field`


#### state_format_pattern(_ = ' {:>5}_ )

#### to_string()

### _class_ pygromos.files.blocks.pertubation_blocks.atom_lam_pertubation_state_dihedral(NR: int, atomI: int, atomJ: int, atomK: int, atomL: int, STATES: Dict[int, int])
Bases: `pygromos.files.blocks._general_blocks._generic_field`


#### state_format_pattern(_ = ' {:>5}_ )

#### to_string()

### _class_ pygromos.files.blocks.pertubation_blocks.atom_mass_type(N: int, ATMAS: float, ATMASN: str, comment: str = '')
Bases: `pygromos.files.blocks._general_blocks._generic_field`


#### to_string()
## pygromos.files.blocks.qmmm_blocks module


### _class_ pygromos.files.blocks.qmmm_blocks.DFTBELEMENTS(content: Union[str, dict, pygromos.utils.typing.DFTBELEMENTS_Type])
Bases: `pygromos.files.blocks.topology_blocks._topology_block`


#### FORCEFIELD(_: pygromos.files.blocks.topology_blocks.FORCEFIEL_ )

#### MAKETOPVERSION(_: pygromos.files.blocks.topology_blocks.MAKETOPVERSIO_ )

### _class_ pygromos.files.blocks.qmmm_blocks.MNDOELEMENTS(content: Union[str, dict, pygromos.utils.typing.MNDOELEMENTS_Type])
Bases: `pygromos.files.blocks.topology_blocks._topology_block`


#### FORCEFIELD(_: pygromos.files.blocks.topology_blocks.FORCEFIEL_ )

#### MAKETOPVERSION(_: pygromos.files.blocks.topology_blocks.MAKETOPVERSIO_ )

### _class_ pygromos.files.blocks.qmmm_blocks.MOPACELEMENTS(content: Union[str, dict, pygromos.utils.typing.MOPACELEMENTS_Type])
Bases: `pygromos.files.blocks.topology_blocks._topology_block`


#### FORCEFIELD(_: pygromos.files.blocks.topology_blocks.FORCEFIEL_ )

#### MAKETOPVERSION(_: pygromos.files.blocks.topology_blocks.MAKETOPVERSIO_ )

### _class_ pygromos.files.blocks.qmmm_blocks.ORCAELEMENTS(content: Union[str, dict, pygromos.utils.typing.ORCAELEMENTS_Type])
Bases: `pygromos.files.blocks.topology_blocks._topology_block`


#### FORCEFIELD(_: pygromos.files.blocks.topology_blocks.FORCEFIEL_ )

#### MAKETOPVERSION(_: pygromos.files.blocks.topology_blocks.MAKETOPVERSIO_ )

### _class_ pygromos.files.blocks.qmmm_blocks.QMUNIT(content: Union[str, dict, pygromos.utils.typing.QMUNIT_Type], QLGL: float = 0.052918, QEGE: float = 2625.5, QCGC: float = 1.0, QIGI: float = 0.1)
Bases: `pygromos.files.blocks.topology_blocks._topology_block`


#### block_to_string()

#### read_content_from_str(content: str)

#### table_header(_: Iterable[str_ _ = ['QLGL', 'QEGE', 'QCGC', 'QIGI'_ )
The QMUNIT block


* **Parameters**

    
    * **QLGL** (*float*) – QM length to Gromos length (e.g. Bohr to nm)


    * **QEGE** (*float*) – QM energy to Gromos energy (e.g. Hartree to kJ / mol)


    * **QCGC** (*float*) – Gromos charge to QM charge


    * **QIGI** (*float*) – QM input units to Gromos input units (e.g. Angstrom to nm)



### _class_ pygromos.files.blocks.qmmm_blocks.QMZONE(content: Union[Iterable[pygromos.files.blocks.qmmm_blocks.qmzone_field], str])
Bases: `pygromos.files.blocks.topology_blocks._topology_table_block`


#### NBON(_: in_ _ = _ )

#### block_to_string()

#### table_header(_: Iterable[str_ _ = ['QMEN', 'QMEI', 'QMEZ', 'QMEB'_ )

#### table_line_type()
alias of `pygromos.files.blocks.qmmm_blocks.qmzone_field`


### _class_ pygromos.files.blocks.qmmm_blocks.TURBOMOLEELEMENTS(content: Union[str, dict, pygromos.utils.typing.TURBOMOLEELEMENTS_Type])
Bases: `pygromos.files.blocks.topology_blocks._topology_block`


#### FORCEFIELD(_: pygromos.files.blocks.topology_blocks.FORCEFIEL_ )

#### MAKETOPVERSION(_: pygromos.files.blocks.topology_blocks.MAKETOPVERSIO_ )

### _class_ pygromos.files.blocks.qmmm_blocks.XTBELEMENTS(content: Union[str, dict, pygromos.utils.typing.XTBELEMENTS_Type])
Bases: `pygromos.files.blocks.topology_blocks._topology_block`


#### FORCEFIELD(_: pygromos.files.blocks.topology_blocks.FORCEFIEL_ )

#### MAKETOPVERSION(_: pygromos.files.blocks.topology_blocks.MAKETOPVERSIO_ )

### _class_ pygromos.files.blocks.qmmm_blocks.qmzone_field(QMEN: str, QMEI: int, QMEZ: int, QMEB: int)
Bases: `pygromos.files.blocks._general_blocks._generic_field`


#### \__init__(QMEN: str, QMEI: int, QMEZ: int, QMEB: int)
One line of the QMZONE block


* **Parameters**

    
    * **QMEN** (*str*) – (QM element name), indicates the atom identifier for this position


    * **QMEI** (*int*) – (QM element iterator), specifies an iterator over the atom positions


    * **QMEZ** (*int*) – (QM element Z), specifies the nuclear charge Z of the atom position


    * **QMEB** (*int*) – (QM element bond), specifies whether a bond can be broken or not (== 0)



#### to_string()
## pygromos.files.blocks.replica_exchange_blocks module


### _class_ pygromos.files.blocks.replica_exchange_blocks.repex_system(s: list, T: float, state_eir: dict)
Bases: `pygromos.files.blocks._general_blocks._generic_gromos_block`


#### block_to_string()

#### comment(_: st_ )

#### content(_: Iterabl_ )

### _class_ pygromos.files.blocks.replica_exchange_blocks.replica_stat(ID: int, partner: int, run: int, Ti: float, Epoti: float, Tj: float, Epotj: float, p: float, s: bool, si: float, sj: float, state_potentials: Optional[dict] = None, potentials: Optional[dict] = None, li: Optional[float] = None, lj: Optional[float] = None)
Bases: `pygromos.files.blocks._general_blocks._generic_field`


#### Epoti(_: floa_ )

#### Epotj(_: floa_ )

#### ID(_: in_ )

#### Ti(_: floa_ )

#### Tj(_: floa_ )

#### \__init__(ID: int, partner: int, run: int, Ti: float, Epoti: float, Tj: float, Epotj: float, p: float, s: bool, si: float, sj: float, state_potentials: Optional[dict] = None, potentials: Optional[dict] = None, li: Optional[float] = None, lj: Optional[float] = None)

* **Parameters**

    
    * **ID**


    * **partner**


    * **run**


    * **Ti**


    * **Epoti**


    * **Tj**


    * **Epotj**


    * **p**


    * **s**


    * **si**


    * **sj**


    * **state_potentials**


    * **potentials**


    * **li**


    * **lj**



#### li(_: floa_ )

#### lj(_: floa_ )

#### p(_: floa_ )

#### partner(_: in_ )

#### run(_: in_ )

#### s(_: boo_ )

#### si(_: floa_ )

#### sj(_: floa_ )

#### state_potentials(_: Dict[str, float_ )

#### to_string()
## pygromos.files.blocks.topology_blocks module


### _class_ pygromos.files.blocks.topology_blocks.ATOMTYPENAME(content: Union[str, Dict[str, str], pygromos.utils.typing.ATOMTYPENAME_Type], FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None)
Bases: `pygromos.files.blocks.topology_blocks._topology_block`


#### FORCEFIELD(_: pygromos.files.blocks.topology_blocks.FORCEFIEL_ )

#### MAKETOPVERSION(_: pygromos.files.blocks.topology_blocks.MAKETOPVERSIO_ )

### _class_ pygromos.files.blocks.topology_blocks.BOND(content: Union[Iterable[pygromos.files.blocks.topology_blocks.top_bond_type], str], FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None, NBON: Optional[int] = None)
Bases: `pygromos.files.blocks.topology_blocks._topology_table_block`


#### NBON(_: in_ _ = _ )

#### \__init__(content: Union[Iterable[pygromos.files.blocks.topology_blocks.top_bond_type], str], FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None, NBON: Optional[int] = None)
> GROMOS BOND block


* **Parameters**

    
    * **content** (*Union[Iterable[top_bond_type], str]*)


    * **FORCEFIELD** (*FORCEFIELD*)


    * **MAKETOPVERSION** (*MAKETOPVERSION*)


    * **NBON** (*int, optional*) – Number of bonds



#### block_to_string()

#### read_content_from_str(content: str)

#### table_header(_: Iterable[str_ _ = ['IB', 'JB', 'ICB'_ )

#### table_line_type()
alias of `pygromos.files.blocks.topology_blocks.top_bond_type`


### _class_ pygromos.files.blocks.topology_blocks.BONDANGLE(content: Union[Iterable[pygromos.files.blocks.topology_blocks.bondangle_type], str], FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None, NTHE: Optional[int] = None)
Bases: `pygromos.files.blocks.topology_blocks._topology_table_block`


#### NTHE(_: in_ )

#### \__init__(content: Union[Iterable[pygromos.files.blocks.topology_blocks.bondangle_type], str], FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None, NTHE: Optional[int] = None)
> GROMOS BONDSTRETCHTYPE block


* **Parameters**

    
    * **content** (*Union[Iterable[bondangle_type], str]*)


    * **FORCEFIELD** (*FORCEFIELD*)


    * **MAKETOPVERSION** (*MAKETOPVERSION*)


    * **NTHE** (*int, optional*) – Number of bondangles



#### block_to_string()

#### read_content_from_str(content: str)

#### table_header(_: Iterable[str_ _ = ['IT', 'JT', 'KT', 'ICT'_ )

#### table_line_type()
alias of `pygromos.files.blocks.topology_blocks.bondangle_type`


### _class_ pygromos.files.blocks.topology_blocks.BONDANGLEBENDTYPE(content: Union[Iterable[pygromos.files.blocks.topology_blocks.bondanglebendtype_type], str], FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None, NBTY: Optional[int] = None)
Bases: `pygromos.files.blocks.topology_blocks._topology_table_block`


#### NBTY(_: in_ )

#### \__init__(content: Union[Iterable[pygromos.files.blocks.topology_blocks.bondanglebendtype_type], str], FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None, NBTY: Optional[int] = None)
> GROMOS BONDSTRETCHTYPE block


* **Parameters**

    
    * **content** (*Union[Iterable[bondanglebendtype_type], str]*)


    * **FORCEFIELD** (*FORCEFIELD*)


    * **MAKETOPVERSION** (*MAKETOPVERSION*)


    * **NBTY** (*int, optional*) – Number of bondstretchtypes



#### block_to_string()

#### read_content_from_str(content: str)

#### table_header(_: Iterable[str_ _ = ['CT', 'CHT', 'T0'_ )

#### table_line_type()
alias of `pygromos.files.blocks.topology_blocks.bondanglebendtype_type`


### _class_ pygromos.files.blocks.topology_blocks.BONDANGLEBENDTYPECODE(content: Union[Iterable[pygromos.files.blocks.topology_blocks.angle_type], str], FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None, NRTTY: Optional[int] = None, NTTY: Optional[int] = None)
Bases: `pygromos.files.blocks.topology_blocks._iterable_topology_block`


#### NRTTY(_: in_ )

#### NTTY(_: in_ )

#### \__init__(content: Union[Iterable[pygromos.files.blocks.topology_blocks.angle_type], str], FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None, NRTTY: Optional[int] = None, NTTY: Optional[int] = None)
> GROMOS bond-stretching parameters


* **Parameters**

    
    * **content** (*Union[Iterable[atom_mass_type], str]*)


    * **FORCEFIELD** (*FORCEFIELD*)


    * **MAKETOPVERSION** (*MAKETOPVERSION*)


    * **NRBTY** (*int, optional*) – Number of bond types


    * **NBTY** (*int, optional*) – Number of maximal bond index.rst



#### block_to_string()

#### read_content_from_str(content: str)

#### table_header(_: Iterable[str_ _ = ['ICT(H)[N]', 'CT[N]', 'CHT[N]', '(T0[N])'_ )

### _class_ pygromos.files.blocks.topology_blocks.BONDANGLEH(content: Union[Iterable[pygromos.files.blocks.topology_blocks.bondangle_type], str], FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None, NTHEH: Optional[int] = None)
Bases: `pygromos.files.blocks.topology_blocks._topology_table_block`


#### NTHEH(_: in_ )

#### \__init__(content: Union[Iterable[pygromos.files.blocks.topology_blocks.bondangle_type], str], FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None, NTHEH: Optional[int] = None)
> GROMOS BONDANGLEH block


* **Parameters**

    
    * **content** (*Union[Iterable[bondangle_type], str]*)


    * **FORCEFIELD** (*FORCEFIELD*)


    * **MAKETOPVERSION** (*MAKETOPVERSION*)


    * **NTHEH** (*int, optional*) – Number of bondangles



#### block_to_string()

#### read_content_from_str(content: str)

#### table_header(_: Iterable[str_ _ = ['ITH', 'JTH', 'KTH', 'ICTH'_ )

#### table_line_type()
alias of `pygromos.files.blocks.topology_blocks.bondangle_type`


### _class_ pygromos.files.blocks.topology_blocks.BONDH(content: Union[Iterable[pygromos.files.blocks.topology_blocks.top_bond_type], str], FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None, NBONH: Optional[int] = None)
Bases: `pygromos.files.blocks.topology_blocks._topology_table_block`


#### NBONH(_: in_ _ = _ )

#### \__init__(content: Union[Iterable[pygromos.files.blocks.topology_blocks.top_bond_type], str], FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None, NBONH: Optional[int] = None)
> GROMOS BONDH block


* **Parameters**

    
    * **content** (*Union[Iterable[top_bond_type], str]*)


    * **FORCEFIELD** (*FORCEFIELD*)


    * **MAKETOPVERSION** (*MAKETOPVERSION*)


    * **NBONH** (*int, optional*) – Number of bonds with H



#### block_to_string()

#### read_content_from_str(content: str)

#### table_header(_: Iterable[str_ _ = ['IBH', 'JBH', 'ICBH'_ )

#### table_line_type()
alias of `pygromos.files.blocks.topology_blocks.top_bond_type`


### _class_ pygromos.files.blocks.topology_blocks.BONDSTRETCHTYPE(content: Union[Iterable[pygromos.files.blocks.topology_blocks.bondstretchtype_type], str], FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None, NBTY: Optional[int] = None)
Bases: `pygromos.files.blocks.topology_blocks._topology_table_block`


#### NBTY(_: in_ )

#### \__init__(content: Union[Iterable[pygromos.files.blocks.topology_blocks.bondstretchtype_type], str], FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None, NBTY: Optional[int] = None)
> GROMOS BONDSTRETCHTYPE block


* **Parameters**

    
    * **content** (*Union[Iterable[bondstretchtype_type], str]*)


    * **FORCEFIELD** (*FORCEFIELD*)


    * **MAKETOPVERSION** (*MAKETOPVERSION*)


    * **NBTY** (*int, optional*) – Number of bondstretchtypes



#### block_to_string()

#### read_content_from_str(content: str)

#### table_header(_: Iterable[str_ _ = ['CB', 'CHB', 'B0'_ )

#### table_line_type()
alias of `pygromos.files.blocks.topology_blocks.bondstretchtype_type`


### _class_ pygromos.files.blocks.topology_blocks.BONDSTRETCHTYPECODE(content: Union[Iterable[pygromos.files.blocks.topology_blocks.bond_type], str], FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None, NRBTY: Optional[int] = None, NBTY: Optional[int] = None)
Bases: `pygromos.files.blocks.topology_blocks._iterable_topology_block`


#### NBTY(_: in_ )

#### NRBTY(_: in_ )

#### \__init__(content: Union[Iterable[pygromos.files.blocks.topology_blocks.bond_type], str], FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None, NRBTY: Optional[int] = None, NBTY: Optional[int] = None)
> GROMOS bond-stretching parameters


* **Parameters**

    
    * **content** (*Union[Iterable[atom_mass_type], str]*)


    * **FORCEFIELD** (*FORCEFIELD*)


    * **MAKETOPVERSION** (*MAKETOPVERSION*)


    * **NRBTY** (*int, optional*) – Number of bond types


    * **NBTY** (*int, optional*) – Number of maximal bond index.rst



#### block_to_string()

#### read_content_from_str(content: str)

#### table_header(_: Iterable[str_ _ = ['ICB(H)[N]', 'CB[N]', 'HB[N]', 'B0[N]'_ )

### _class_ pygromos.files.blocks.topology_blocks.CONSTRAINT(content: str, FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None, NCON: Optional[int] = None)
Bases: `pygromos.files.blocks.topology_blocks._topology_table_block`


#### NCON(_: in_ )

#### block_to_string()

#### read_content_from_str(content: str)

#### table_header(_: Iterable[str_ _ = ['IC', 'JC', 'ICC'_ )

#### table_line_type()
alias of `pygromos.files.blocks.topology_blocks.constraint_type`


### _class_ pygromos.files.blocks.topology_blocks.CROSSDIHEDRAL(content: Union[Iterable[pygromos.files.blocks.topology_blocks.crossgihedral_type], str], FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None, NPHI: Optional[int] = None)
Bases: `pygromos.files.blocks.topology_blocks._topology_table_block`


#### NPHI(_: in_ )

#### \__init__(content: Union[Iterable[pygromos.files.blocks.topology_blocks.crossgihedral_type], str], FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None, NPHI: Optional[int] = None)
[summary]


* **Parameters**

    
    * **content** (*Union[Iterable[crossgihedral_type], str]*) – [description]


    * **FORCEFIELD** (*FORCEFIELD, optional*) – [description], by default None


    * **MAKETOPVERSION** (*MAKETOPVERSION, optional*) – [description], by default None


    * **NPHI** (*[type], optional*) – number of cross dihedrals NOT involving H atoms in solute, by default None



#### block_to_string()

#### read_content_from_str(content: str)

#### table_header(_: Iterable[str_ _ = ['AP', 'BP', 'CP', 'DP', 'EP', 'FP', 'GP', 'HP', 'ICC'_ )

#### table_line_type()
alias of `pygromos.files.blocks.topology_blocks.crossgihedralh_type`


### _class_ pygromos.files.blocks.topology_blocks.CROSSDIHEDRALH(content: Union[Iterable[pygromos.files.blocks.topology_blocks.crossgihedralh_type], str], FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None, NPHIH: Optional[int] = None)
Bases: `pygromos.files.blocks.topology_blocks._topology_table_block`


#### NPHIH(_: in_ )

#### \__init__(content: Union[Iterable[pygromos.files.blocks.topology_blocks.crossgihedralh_type], str], FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None, NPHIH: Optional[int] = None)
[summary]


* **Parameters**

    
    * **content** (*Union[Iterable[crossgihedralh_type], str]*) – [description]


    * **FORCEFIELD** (*FORCEFIELD, optional*) – [description], by default None


    * **MAKETOPVERSION** (*MAKETOPVERSION, optional*) – [description], by default None


    * **NPHIH** (*[type], optional*) – number of cross dihedrals involving H atoms in solute, by default None



#### block_to_string()

#### read_content_from_str(content: str)

#### table_header(_: Iterable[str_ _ = ['APH', 'BPH', 'CPH', 'DPH', 'EPH', 'FPH', 'GPH', 'HPH', 'ICCH'_ )

#### table_line_type()
alias of `pygromos.files.blocks.topology_blocks.crossgihedralh_type`


### _class_ pygromos.files.blocks.topology_blocks.DIHEDRAL(content: Union[Iterable[pygromos.files.blocks.topology_blocks.top_dihedral_type], str], FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None, NPHI: Optional[int] = None)
Bases: `pygromos.files.blocks.topology_blocks._topology_table_block`


#### NPHI(_: in_ )

#### \__init__(content: Union[Iterable[pygromos.files.blocks.topology_blocks.top_dihedral_type], str], FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None, NPHI: Optional[int] = None)
> GROMOS DIHEDRAL block


* **Parameters**

    
    * **content** (*Union[Iterable[top_dihedral_type], str]*)


    * **FORCEFIELD** (*FORCEFIELD*)


    * **MAKETOPVERSION** (*MAKETOPVERSION*)


    * **NPHI** (*int, optional*) – Number of tors dihedral



#### block_to_string()

#### read_content_from_str(content: str)

#### table_header(_: Iterable[str_ _ = ['IP', 'JP', 'KP', 'LP', 'ICP'_ )

#### table_line_type()
alias of `pygromos.files.blocks.topology_blocks.top_dihedral_type`


### _class_ pygromos.files.blocks.topology_blocks.DIHEDRALH(content: Union[Iterable[pygromos.files.blocks.topology_blocks.dihedralh_type], str], FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None, NPHIH: Optional[int] = None)
Bases: `pygromos.files.blocks.topology_blocks._topology_table_block`


#### NPHIH(_: in_ )

#### \__init__(content: Union[Iterable[pygromos.files.blocks.topology_blocks.dihedralh_type], str], FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None, NPHIH: Optional[int] = None)
> GROMOS DIHEDRAL block


* **Parameters**

    
    * **content** (*Union[Iterable[dihedralh_type], str]*)


    * **FORCEFIELD** (*FORCEFIELD*)


    * **MAKETOPVERSION** (*MAKETOPVERSION*)


    * **NPHIH** (*int, optional*) – Number of dihedralH



#### block_to_string()

#### read_content_from_str(content: str)

#### table_header(_: Iterable[str_ _ = ['IPH', 'JPH', 'KPH', 'LPH', 'ICPH'_ )

#### table_line_type()
alias of `pygromos.files.blocks.topology_blocks.dihedralh_type`


### _class_ pygromos.files.blocks.topology_blocks.DISTANCERESSPEC(KDISH: Optional[int] = None, KDISC: Optional[int] = None, RESTRAINTHEADER: Optional[list] = None, RESTRAINTS: Optional[list] = None, content: Optional[List[str]] = None)
Bases: `pygromos.files.blocks._general_blocks._generic_gromos_block`


#### \__init__(KDISH: Optional[int] = None, KDISC: Optional[int] = None, RESTRAINTHEADER: Optional[list] = None, RESTRAINTS: Optional[list] = None, content: Optional[List[str]] = None)

* **Parameters**

    
    * **KDISH**


    * **KDISC**


    * **RESTRAINTHEADER**


    * **RESTRAINTS**



#### block_to_string()

#### comment(_: st_ )

#### content(_: Iterabl_ )

#### read_content_from_str(content: List[str])

### _class_ pygromos.files.blocks.topology_blocks.FORCEFIELD(NAME: Optional[str] = None, content: Optional[List[str]] = None)
Bases: `pygromos.files.blocks._general_blocks._generic_gromos_block`


#### NAME(_: st_ )

#### block_to_string()

### _class_ pygromos.files.blocks.topology_blocks.IMPDIHEDRAL(content: Union[Iterable[pygromos.files.blocks.topology_blocks.impdihedral_type], str], FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None, NQHI: Optional[int] = None)
Bases: `pygromos.files.blocks.topology_blocks._topology_table_block`


#### NQHI(_: in_ )

#### \__init__(content: Union[Iterable[pygromos.files.blocks.topology_blocks.impdihedral_type], str], FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None, NQHI: Optional[int] = None)
> GROMOS IMPDIHEDRAL block


* **Parameters**

    
    * **content** (*Union[Iterable[impdihedral_type], str]*)


    * **FORCEFIELD** (*FORCEFIELD*)


    * **MAKETOPVERSION** (*MAKETOPVERSION*)


    * **NQHI** (*int, optional*) – Number of impdihedral



#### block_to_string()

#### read_content_from_str(content: str)

#### table_header(_: Iterable[str_ _ = ['IQ', 'JQ', 'KQ', 'LQ', 'ICQ'_ )

#### table_line_type()
alias of `pygromos.files.blocks.topology_blocks.impdihedral_type`


### _class_ pygromos.files.blocks.topology_blocks.IMPDIHEDRALH(content: Union[Iterable[pygromos.files.blocks.topology_blocks.impdihedralh_type], str], FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None, NQHIH: Optional[int] = None)
Bases: `pygromos.files.blocks.topology_blocks._topology_table_block`


#### NQHIH(_: in_ )

#### \__init__(content: Union[Iterable[pygromos.files.blocks.topology_blocks.impdihedralh_type], str], FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None, NQHIH: Optional[int] = None)
> GROMOS IMPDIHEDRALH block


* **Parameters**

    
    * **content** (*Union[Iterable[impdihedralh_type], str]*)


    * **FORCEFIELD** (*FORCEFIELD*)


    * **MAKETOPVERSION** (*MAKETOPVERSION*)


    * **NQHIH** (*int, optional*) – Number of impdihedralH



#### block_to_string()

#### read_content_from_str(content: str)

#### table_header(_: Iterable[str_ _ = ['IQH', 'JQH', 'KQH', 'LQH', 'ICQH'_ )

#### table_line_type()
alias of `pygromos.files.blocks.topology_blocks.impdihedralh_type`


### _class_ pygromos.files.blocks.topology_blocks.IMPDIHEDRALTYPE(content: Union[Iterable[pygromos.files.blocks.topology_blocks.impdihedraltype_type], str], FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None, NQTY: Optional[int] = None)
Bases: `pygromos.files.blocks.topology_blocks._topology_table_block`


#### NQTY(_: in_ )

#### \__init__(content: Union[Iterable[pygromos.files.blocks.topology_blocks.impdihedraltype_type], str], FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None, NQTY: Optional[int] = None)
> GROMOS IMPDIHEDRALTYPE block


* **Parameters**

    
    * **content** (*Union[Iterable[impdihedraltype_type], str]*)


    * **FORCEFIELD** (*FORCEFIELD*)


    * **MAKETOPVERSION** (*MAKETOPVERSION*)


    * **NQTY** (*int, optional*) – Number of impdihedraltype



#### block_to_string()

#### read_content_from_str(content: str)

#### table_header(_: Iterable[str_ _ = ['CQ', 'Q0'_ )

#### table_line_type()
alias of `pygromos.files.blocks.topology_blocks.impdihedraltype_type`


### _class_ pygromos.files.blocks.topology_blocks.IMPDIHEDRALTYPECODE(content: Union[Iterable[pygromos.files.blocks.topology_blocks.improper_dihedral_type], str], FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None, NRQTY: Optional[int] = None, NQTY: Optional[int] = None)
Bases: `pygromos.files.blocks.topology_blocks._iterable_topology_block`


#### NQTY(_: in_ )

#### NRQTY(_: in_ )

#### \__init__(content: Union[Iterable[pygromos.files.blocks.topology_blocks.improper_dihedral_type], str], FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None, NRQTY: Optional[int] = None, NQTY: Optional[int] = None)
> GROMOS improper dihedral type parameters


* **Parameters**

    
    * **content** (*Union[Iterable[atom_mass_type], str]*)


    * **FORCEFIELD** (*FORCEFIELD*)


    * **MAKETOPVERSION** (*MAKETOPVERSION*)


    * **NRQTY** (*int, optional*) – Number of improperDihedrals types


    * **NQTY** (*int, optional*) – Number of maximal improperDihedrals index.rst



#### block_to_string()

#### read_content_from_str(content: str)

#### table_header(_: Iterable[str_ _ = ['ICQ', 'CQ', 'Q0'_ )

### _class_ pygromos.files.blocks.topology_blocks.LJEXCEPTIONS(content: Union[Iterable[pygromos.files.blocks.topology_blocks.ljexception_type], str], FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None, NEX: Optional[int] = None)
Bases: `pygromos.files.blocks.topology_blocks._topology_table_block`


#### NEX(_: in_ )

#### \__init__(content: Union[Iterable[pygromos.files.blocks.topology_blocks.ljexception_type], str], FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None, NEX: Optional[int] = None)
> GROMOS LJEXCEPTIONS block


* **Parameters**

    
    * **content** (*Union[Iterable[ljexception_type], str]*)


    * **FORCEFIELD** (*FORCEFIELD*)


    * **MAKETOPVERSION** (*MAKETOPVERSION*)


    * **NEX** (*int, optional*) – Number of LJEXCEPTIONS



#### block_to_string()

#### read_content_from_str(content: str)

#### table_header(_: Iterable[str_ _ = ['AT1', 'AT2', 'C12', 'C6'_ )

#### table_line_type()
alias of `pygromos.files.blocks.topology_blocks.ljexception_type`


### _class_ pygromos.files.blocks.topology_blocks.LJPARAMETERS(content: Union[Iterable[pygromos.files.blocks.topology_blocks.ljparameters_type], str], FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None, NRATT2: Optional[int] = None)
Bases: `pygromos.files.blocks.topology_blocks._topology_table_block`


#### NRATT2(_: in_ )

#### \__init__(content: Union[Iterable[pygromos.files.blocks.topology_blocks.ljparameters_type], str], FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None, NRATT2: Optional[int] = None)
> GROMOS LJPARAMETERS block


* **Parameters**

    
    * **content** (*Union[Iterable[ljparameters_type], str]*)


    * **FORCEFIELD** (*FORCEFIELD*)


    * **MAKETOPVERSION** (*MAKETOPVERSION*)


    * **NRATT2** (*int, optional*) – Number of LJPARAMETERS



#### block_to_string()

#### read_content_from_str(content: str)

#### table_header(_: Iterable[str_ _ = ['IAC', 'JAC', 'C12', 'C6', 'CS12', 'CS6'_ )

#### table_line_type()
alias of `pygromos.files.blocks.topology_blocks.ljparameters_type`


### _class_ pygromos.files.blocks.topology_blocks.MAKETOPVERSION(VERSION: Optional[str] = None, content: Optional[List[str]] = None)
Bases: `pygromos.files.blocks._general_blocks._generic_gromos_block`


#### VERSION(_: st_ )

#### block_to_string()

### _class_ pygromos.files.blocks.topology_blocks.MASSATOMTYPECODE(content: Union[Iterable[pygromos.files.blocks.topology_blocks.atom_mass_type], Iterable[str]], FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None, NRMATY: Optional[int] = None, NMATY: Optional[int] = None)
Bases: `pygromos.files.blocks.topology_blocks._iterable_topology_block`


#### NMATY(_: in_ )

#### NRMATY(_: in_ )

#### block_to_string()

#### read_content_from_str(content: Union[str, List[str]])

#### table_header(_: Iterable[str_ _ = ['N', 'ATMAS', 'ATMASN'_ )

### _class_ pygromos.files.blocks.topology_blocks.MIXEDATOMLJPAIR(content: List[str], NRMTT: Optional[int] = None, FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None)
Bases: `pygromos.files.blocks.topology_blocks._iterable_topology_block`


#### NRMTT(_: in_ )

#### \__init__(content: List[str], NRMTT: Optional[int] = None, FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None)

* **Parameters**

    
    * **content**


    * **NRATT**


    * **FORCEFIELD**


    * **MAKETOPVERSION**



#### block_to_string()

#### read_content_from_str(content: str)

#### table_header(_: Iterable[str_ _ = ['IACI', 'IACJ', 'C6', 'C12_1', 'C12_2', 'C12_3'_ )
GROMOS 43A1 normal van der Waals parameters for mixed atom type pairs (I,J)


### _class_ pygromos.files.blocks.topology_blocks.PERTATOMPARAM(STATEATOMS: Optional[List[pygromos.files.blocks.topology_blocks.atom_lam_pertubation_state]] = None, STATEATOMHEADER: Optional[Tuple[str]] = None, NJLA: Optional[int] = None, STATEIDENTIFIERS: Optional[List[str]] = None, dummy_IAC: int = 22, dummy_CHARGE: int = 0.0, content: Optional[List[str]] = None)
Bases: `pygromos.files.blocks._general_blocks._generic_gromos_block`


#### add_state_atoms(state_atoms: List[pygromos.files.blocks.topology_blocks.atom_lam_pertubation_state])
This function can add states and atoms, but also overwrite state values of existing atoms.
If a new state is defined only for a subset of atoms, all other atoms are set to the default dummy.
If a new atom misses a state definition, this state will be set to dummy.
:Parameters: **state_atoms** (*List[atom_eds_pertubation_state]*)


#### block_to_string()

#### comment(_: st_ )

#### content(_: Iterabl_ )

#### delete_atom(atomNR: Union[int, List[int]])
This function removes atom lines from the ptp file.
:Parameters: **atomNR** (*int*) – atom to be removed.


#### delete_state(stateIDs: Optional[Union[int, List[int]]] = None, stateNames: Optional[Union[str, List[str]]] = None)
This function deletes an state column.
:Parameters: **stateIDs** (*int*) – number of the state


#### _property_ nStates(_: in_ )

#### _property_ nTotalStateAtoms(_: in_ )

#### read_content_from_str(content: List[str])

#### _property_ states(_: dic_ )

### _class_ pygromos.files.blocks.topology_blocks.PHYSICALCONSTANTS(content: Union[str, dict, pygromos.utils.typing.PHYSICALCONSTANTS_Type], FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None)
Bases: `pygromos.files.blocks.topology_blocks._topology_block`


#### FORCEFIELD(_: pygromos.files.blocks.topology_blocks.FORCEFIEL_ )

#### MAKETOPVERSION(_: pygromos.files.blocks.topology_blocks.MAKETOPVERSIO_ )

#### block_to_string()

#### read_content_from_str(content: str)

### _class_ pygromos.files.blocks.topology_blocks.PRESSUREGROUPS(content: Optional[Union[str, dict]] = None, FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None, NSM: Optional[int] = None, NSP: Optional[List[int]] = None)
Bases: `pygromos.files.blocks.topology_blocks._generic_topology_groups`


#### NSM(_: in_ )

#### NSP(_: List[int_ )

### _class_ pygromos.files.blocks.topology_blocks.RESNAME(content: Union[str, Dict[str, str], pygromos.utils.typing.RESNAME_Type], FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None)
Bases: `pygromos.files.blocks.topology_blocks._topology_block`


#### FORCEFIELD(_: pygromos.files.blocks.topology_blocks.FORCEFIEL_ )

#### MAKETOPVERSION(_: pygromos.files.blocks.topology_blocks.MAKETOPVERSIO_ )

### _class_ pygromos.files.blocks.topology_blocks.SCALEDINTERACTIONS(values=None, content=None)
Bases: `pygromos.files.blocks._general_blocks._generic_gromos_block`


#### \__init__(values=None, content=None)
Not exactly sure what these parameters do


#### block_to_string()

#### comment(_: st_ )

#### content(_: Iterabl_ )

#### read_content_from_str(content)

### _class_ pygromos.files.blocks.topology_blocks.SINGLEATOMLJPAIR(content: List[str], NRATT: Optional[int] = None, FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None)
Bases: `pygromos.files.blocks.topology_blocks._iterable_topology_block`


#### NRATT(_: in_ )

#### \__init__(content: List[str], NRATT: Optional[int] = None, FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None)

* **Parameters**

    
    * **content**


    * **NRATT**


    * **FORCEFIELD**


    * **MAKETOPVERSION**



#### block_to_string()

#### read_content_from_str(content: str)

#### table_header(_ = ['IAC', 'TYPE', 'C6', 'C12(1)', 'C12(2)', 'C12(3)'_ )

### _class_ pygromos.files.blocks.topology_blocks.SOLUTEATOM(content: Union[str, Dict[str, str], pygromos.utils.typing._iterable_topology_block_Type], FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None)
Bases: `pygromos.files.blocks.topology_blocks._iterable_topology_block`


#### NRP(_: in_ )

#### block_to_string()

#### read_content_from_str(content: str)

#### table_header(_: Iterable[str_ _ = ['IB', 'JB', 'ICB'_ )

### _class_ pygromos.files.blocks.topology_blocks.SOLUTEMOLECULES(content: Optional[Union[str, dict]] = None, FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None, NSM: Optional[int] = None, NSP: Optional[List[int]] = None)
Bases: `pygromos.files.blocks.topology_blocks._generic_topology_groups`


#### NSM(_: in_ )

#### NSP(_: List[int_ )

### _class_ pygromos.files.blocks.topology_blocks.SOLVENTATOM(content: Union[Iterable[pygromos.files.blocks.topology_blocks.solventatom_type], str], FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None, NRAM: Optional[int] = None)
Bases: `pygromos.files.blocks.topology_blocks._topology_table_block`


#### NRAM(_: in_ )

#### \__init__(content: Union[Iterable[pygromos.files.blocks.topology_blocks.solventatom_type], str], FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None, NRAM: Optional[int] = None)
> GROMOS solventatom block


* **Parameters**

    
    * **content** (*Union[Iterable[solventatom_type], str]*)


    * **FORCEFIELD** (*FORCEFIELD*)


    * **MAKETOPVERSION** (*MAKETOPVERSION*)


    * **NRAM** (*int, optional*) – Number of solventatom



#### block_to_string()

#### read_content_from_str(content: str)

#### table_header(_: Iterable[str_ _ = ['I', 'ANMS', 'IACS', 'MASS', 'CGS'_ )

#### table_line_type()
alias of `pygromos.files.blocks.topology_blocks.solventatom_type`


### _class_ pygromos.files.blocks.topology_blocks.SOLVENTCONSTR(content: Union[Iterable[pygromos.files.blocks.topology_blocks.solventconstr_type], str], FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None, NCONS: Optional[int] = None)
Bases: `pygromos.files.blocks.topology_blocks._topology_table_block`


#### NCONS(_: in_ )

#### \__init__(content: Union[Iterable[pygromos.files.blocks.topology_blocks.solventconstr_type], str], FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None, NCONS: Optional[int] = None)
> GROMOS SOLVENTCONSTR block


* **Parameters**

    
    * **content** (*Union[Iterable[solventconstr_type], str]*)


    * **FORCEFIELD** (*FORCEFIELD*)


    * **MAKETOPVERSION** (*MAKETOPVERSION*)


    * **NCONS** (*int, optional*) – Number of SOLVENTCONSTR



#### block_to_string()

#### read_content_from_str(content: str)

#### table_header(_: Iterable[str_ _ = ['ICONS', 'JCONS', 'CONS'_ )

#### table_line_type()
alias of `pygromos.files.blocks.topology_blocks.solventconstr_type`


### _class_ pygromos.files.blocks.topology_blocks.SPECATOMLJPAIR(content: List[str], NRST: Optional[int] = None, FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None)
Bases: `pygromos.files.blocks.topology_blocks._iterable_topology_block`


#### NRST(_: in_ )

#### \__init__(content: List[str], NRST: Optional[int] = None, FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None)

* **Parameters**

    
    * **content**


    * **NRATT**


    * **FORCEFIELD**


    * **MAKETOPVERSION**



#### read_content_from_str(content: str)

#### table_header(_: Iterable[str_ _ = ['???'_ )

### _class_ pygromos.files.blocks.topology_blocks.TEMPERATUREGROUPS(content: Optional[Union[str, dict]] = None, FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None, NSM: Optional[int] = None, NSP: Optional[List[int]] = None)
Bases: `pygromos.files.blocks.topology_blocks._generic_topology_groups`


#### NSM(_: in_ )

#### NSP(_: List[int_ )

### _class_ pygromos.files.blocks.topology_blocks.TOPVERSION(content: Union[str, Dict[str, str], pygromos.utils.typing.TOPVERSION_Type], FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None)
Bases: `pygromos.files.blocks.topology_blocks._topology_block`


#### FORCEFIELD(_: pygromos.files.blocks.topology_blocks.FORCEFIEL_ )

#### MAKETOPVERSION(_: pygromos.files.blocks.topology_blocks.MAKETOPVERSIO_ )

### _class_ pygromos.files.blocks.topology_blocks.TORSDIHEDRALTYPE(content: Union[Iterable[pygromos.files.blocks.topology_blocks.torsdihedraltype_type], str], FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None, NPTY: Optional[int] = None)
Bases: `pygromos.files.blocks.topology_blocks._topology_table_block`


#### NPTY(_: in_ )

#### \__init__(content: Union[Iterable[pygromos.files.blocks.topology_blocks.torsdihedraltype_type], str], FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None, NPTY: Optional[int] = None)
> GROMOS IMPDIHEDRAL block


* **Parameters**

    
    * **content** (*Union[Iterable[torsdihedraltype_type], str]*)


    * **FORCEFIELD** (*FORCEFIELD*)


    * **MAKETOPVERSION** (*MAKETOPVERSION*)


    * **NQTY** (*int, optional*) – Number of torsion dihedrals



#### block_to_string()

#### read_content_from_str(content: str)

#### table_header(_: Iterable[str_ _ = ['CP', 'PD', 'NP'_ )

#### table_line_type()
alias of `pygromos.files.blocks.topology_blocks.torsdihedraltype_type`


### _class_ pygromos.files.blocks.topology_blocks.TORSDIHEDRALTYPECODE(content: Union[Iterable[pygromos.files.blocks.topology_blocks.dihedral_type], str], FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None, NRPTY: Optional[int] = None, NPTY: Optional[int] = None)
Bases: `pygromos.files.blocks.topology_blocks._iterable_topology_block`


#### NPTY(_: in_ )

#### NRPTY(_: in_ )

#### \__init__(content: Union[Iterable[pygromos.files.blocks.topology_blocks.dihedral_type], str], FORCEFIELD: Optional[pygromos.files.blocks.topology_blocks.FORCEFIELD] = None, MAKETOPVERSION: Optional[pygromos.files.blocks.topology_blocks.MAKETOPVERSION] = None, NRPTY: Optional[int] = None, NPTY: Optional[int] = None)
> GROMOS (trigonometric) dihedral torsional angle parameters


* **Parameters**

    
    * **content** (*Union[Iterable[atom_mass_type], str]*)


    * **FORCEFIELD** (*FORCEFIELD*)


    * **MAKETOPVERSION** (*MAKETOPVERSION*)


    * **NRPTY** (*int, optional*) – Number of dihedral-angle types


    * **NPTY** (*int, optional*) – Number of maximal dihedral-angle index.rst



#### block_to_string()

#### read_content_from_str(content: str)

#### table_header(_: Iterable[str_ _ = ['ICP(H)[N]', 'CP[N]', 'PD', 'NP'_ )

### _class_ pygromos.files.blocks.topology_blocks.angle_type(ICT: int, CT: float, CHT: float, T0: float, atomI: Union[str, Iterable[str]], atomJ: Union[str, Iterable[str]], atomK: Union[str, Iterable[str]], specialNumber: int)
Bases: `pygromos.files.blocks._general_blocks._generic_field`


#### \__init__(ICT: int, CT: float, CHT: float, T0: float, atomI: Union[str, Iterable[str]], atomJ: Union[str, Iterable[str]], atomK: Union[str, Iterable[str]], specialNumber: int)
> GROMOS bond-stretching parameters for one possible bond


* **Parameters**

    
    * **ICT** (*int*) – Bond-angle type code


    * **CT** (*float*) – Non-harmonic force constant


    * **CHT** (*float*) – Harmonic force constant


    * **T0** (*float*) – Ideal bond angle


    * **atomI** (*str, Iterable[str]*) – possible atom I useages


    * **atomJ** (*str, Iterable[str]*) – possible atom J useages


    * **atomK** (*str, Iterable[str]*) – possible atom K useages


    * **specialNumber** (*int*) – No Idea @Todo: refactor later!



#### to_string()

### _class_ pygromos.files.blocks.topology_blocks.atom_lam_pertubation_state(NR: int, RES: int, NAME: str, STATES: Dict[int, pygromos.files.blocks.topology_blocks.pertubationLamState], ALPHLJ: float = 1.0, ALPHCRF: float = 1.0)
Bases: `pygromos.files.blocks._general_blocks._generic_field`


#### state_format_pattern(_ = ' {:>5} {:>5} {:>10.5f}_ )

#### to_string()

### _class_ pygromos.files.blocks.topology_blocks.atom_mass_type(N: int, ATMAS: float, ATMASN: str, comment: str = '')
Bases: `pygromos.files.blocks._general_blocks._generic_field`


#### to_string()

### _class_ pygromos.files.blocks.topology_blocks.atom_pair_distanceRes(i1: int, j1: int, k1: int, l1: int, type1: Union[pygromos.files.blocks.topology_blocks.geometric_code, int], i2: int, j2: int, k2: int, l2: int, type2: Union[pygromos.files.blocks.topology_blocks.geometric_code, int], r0: float, w0: float, rah: Union[pygromos.files.blocks.topology_blocks.distant_Restraint_Type, int], comment: str = '')
Bases: `pygromos.files.blocks._general_blocks._generic_field`


#### \__init__(i1: int, j1: int, k1: int, l1: int, type1: Union[pygromos.files.blocks.topology_blocks.geometric_code, int], i2: int, j2: int, k2: int, l2: int, type2: Union[pygromos.files.blocks.topology_blocks.geometric_code, int], r0: float, w0: float, rah: Union[pygromos.files.blocks.topology_blocks.distant_Restraint_Type, int], comment: str = '')

* **Parameters**

    
    * **i1** (*int*) – id of atom i of first molecule


    * **j1** (*int*) – id of atom j of first molecule


    * **k1** (*int*) – id of atom k of first molecule


    * **l1** (*int*) – id of atom l of first molecule


    * **type1** (*int*) – geometric restraintype of first molecule


    * **i2** (*int*) – id of atom i of second molecule


    * **j2** (*int*) – id of atom j of second molecule


    * **k2** (*int*) – id of atom k of second molecule


    * **l2** (*int*) – id of atom l of second molecule


    * **type2** (*int*) – geometric restraintype of second molecule


    * **r0** (*float*) – radius_0 of restraint


    * **w0** (*float*) – weighting of restraint


    * **rah** (*int*) – restraint_type


    * **comment** (*str, optional*) – comment for this restraint



#### to_string()

### _class_ pygromos.files.blocks.topology_blocks.bond_type(ICB: int, CB: float, HB: float, B0: float, atomI: Union[str, Iterable[str]], atomJ: Union[str, Iterable[str]], specialNumber: int)
Bases: `pygromos.files.blocks._general_blocks._generic_field`


#### \__init__(ICB: int, CB: float, HB: float, B0: float, atomI: Union[str, Iterable[str]], atomJ: Union[str, Iterable[str]], specialNumber: int)
> GROMOS bond-stretching parameters for one possible bond


* **Parameters**

    
    * **ICB** (*int*) – Bond type code


    * **CB** (*float*) – Quartic bond-stretch force constant


    * **HB** (*float*) – Harmonic bond-stretch force constant


    * **B0** (*float*) – Ideal bond length


    * **atomI** (*str, Iterable[str]*) – possible atom I useages


    * **atomJ** (*str, Iterable[str]*) – possible atom J useages


    * **specialNumber** (*int*) – No Idea @Todo: refactor later!



#### to_string()

### _class_ pygromos.files.blocks.topology_blocks.bondangle_type(IT: int, JT: int, KT: int, ICT: int)
Bases: `pygromos.files.blocks._general_blocks._generic_field`


#### \__init__(IT: int, JT: int, KT: int, ICT: int)
> GROMOS bondangleype for a single pair


* **Parameters**

    
    * **IT** (*int*) – atom number in angle


    * **JT** (*int*) – atom number in angle


    * **KT** (*int*) – atom number in angle


    * **ICT** (*int*) – bond angle type code



#### to_string()

### _class_ pygromos.files.blocks.topology_blocks.bondanglebendtype_type(CT: float, CHT: float, T0: float)
Bases: `pygromos.files.blocks._general_blocks._generic_field`


#### \__init__(CT: float, CHT: float, T0: float)
> GROMOS bondanglebendtype for a single angle


* **Parameters**

    
    * **CT** (*float*) – quartic force constant


    * **CHT** (*float*) – harmonic force constant


    * **T0** (*float*) – angle at minimum energy



#### to_string()

### _class_ pygromos.files.blocks.topology_blocks.bondstretchtype_type(CB: float, CHB: float, B0: float)
Bases: `pygromos.files.blocks._general_blocks._generic_field`


#### \__init__(CB: float, CHB: float, B0: float)
> GROMOS bondstretchtype for a single pair


* **Parameters**

    
    * **CB** (*float*) – quartic force constant


    * **CHB** (*float*) – harmonic force constant


    * **B0** (*float*) – bond length at minimum energy



#### to_string()

### _class_ pygromos.files.blocks.topology_blocks.constraint_type(IC: int, JC: int, ICC: float)
Bases: `pygromos.files.blocks._general_blocks._generic_field`


#### \__init__(IC: int, JC: int, ICC: float)
[summary]


* **Parameters**

    
    * **IC** (*int*) – [description]


    * **JC** (*int*) – [description]


    * **ICC** (*float*) – [description]



#### to_string()

### _class_ pygromos.files.blocks.topology_blocks.crossgihedral_type(AP: int, BP: int, CP: int, DP: int, EP: int, FP: int, GP: int, HP: int, ICC: int)
Bases: `pygromos.files.blocks._general_blocks._generic_field`


#### \__init__(AP: int, BP: int, CP: int, DP: int, EP: int, FP: int, GP: int, HP: int, ICC: int)
GROMOS Cross Dihedral type for NON H Atoms


* **Parameters**

    
    * **AP** (*int*) – number of atoms forming a dihedral


    * **BP** (*int*) – number of atoms forming a dihedral


    * **CP** (*int*) – number of atoms forming a dihedral


    * **DP** (*int*) – number of atoms forming a dihedral


    * **EP** (*int*) – number of atoms forming a dihedral


    * **FP** (*int*) – number of atoms forming a dihedral


    * **GP** (*int*) – number of atoms forming a dihedral


    * **HP** (*int*) – number of atoms forming a dihedral


    * **ICC** (*int*) – dihedral type code



#### to_string()

### _class_ pygromos.files.blocks.topology_blocks.crossgihedralh_type(APH: int, BPH: int, CPH: int, DPH: int, EPH: int, FPH: int, GPH: int, HPH: int, ICCH: int)
Bases: `pygromos.files.blocks._general_blocks._generic_field`


#### \__init__(APH: int, BPH: int, CPH: int, DPH: int, EPH: int, FPH: int, GPH: int, HPH: int, ICCH: int)
GROMOS Cross Dihedral type for H


* **Parameters**

    
    * **APH** (*int*) – number of atoms forming a dihedral


    * **BPH** (*int*) – number of atoms forming a dihedral


    * **CPH** (*int*) – number of atoms forming a dihedral


    * **DPH** (*int*) – number of atoms forming a dihedral


    * **EPH** (*int*) – number of atoms forming a dihedral


    * **FPH** (*int*) – number of atoms forming a dihedral


    * **GPH** (*int*) – number of atoms forming a dihedral


    * **HPH** (*int*) – number of atoms forming a dihedral


    * **ICCH** (*int*) – dihedral type code



#### to_string()

### _class_ pygromos.files.blocks.topology_blocks.dihedral_type(ICP: int, CP: float, PD: float, NP: int, atomI: str, atomJ: str, atomK: str, atomL: str, special_number: float, concrete_example: str = '')
Bases: `pygromos.files.blocks._general_blocks._generic_field`


#### \__init__(ICP: int, CP: float, PD: float, NP: int, atomI: str, atomJ: str, atomK: str, atomL: str, special_number: float, concrete_example: str = '')
> GROMOS improper (harmonic) dihedral angle parameters


* **Parameters**

    
    * **ICQ** (*int*) – Dihedral-angle type code


    * **CQ** (*float*) – Force constant


    * **Q0** (*float*) – Phase shift


    * **NP** (*float*) – Multiplicity


    * **atomI** (*str*) – Examples for atomI


    * **atomJ** (*str*) – Examples for atomJ


    * **atomK** (*str*) – Examples for atomK


    * **atomL** (*str*) – Examples for atomL


    * **special_number** (*int*) – No Idea @Todo: refactor later!


    * **concrete_example** (*str, optional*) – giving a concrete example



#### to_string()

### _class_ pygromos.files.blocks.topology_blocks.dihedralh_type(IPH: int, JPH: int, KPH: int, LPH: int, ICPH: int)
Bases: `pygromos.files.blocks._general_blocks._generic_field`


#### \__init__(IPH: int, JPH: int, KPH: int, LPH: int, ICPH: int)
> GROMOS dihedral for a single pair


#### to_string()

### _class_ pygromos.files.blocks.topology_blocks.distant_Restraint_Type(value)
Bases: `enum.Enum`

An enumeration.


#### full_harmonic_distance_restraint(_ = _ )

#### half_harmonic_attractive(_ = _ )

#### half_harmonic_repulsive(_ = -_ )

### _class_ pygromos.files.blocks.topology_blocks.geometric_code(value)
Bases: `enum.Enum`

An enumeration.


#### pseudo_H_atom_goc_H_atoms_CH3(_ = _ )

#### pseudo_H_atoms_goc_of_three_CH3(_ = _ )

#### pseudo_H_atoms_goc_of_two_CH3(_ = _ )

#### real_atom(_ = _ )

#### virtual_H_atom_aliphaticC(_ = _ )

#### virtual_H_atom_aliphaticC_strange(_ = _ )

#### virtual_H_atom_aromaticC(_ = _ )

#### virtual_H_atoms_goc_aliph(_ = _ )

#### virtual_atoms_cog(_ = -_ )

#### virtual_atoms_com(_ = -_ )

### _class_ pygromos.files.blocks.topology_blocks.impdihedral_type(IQ: int, JQ: int, KQ: int, LQ: int, ICQ: int)
Bases: `pygromos.files.blocks._general_blocks._generic_field`


#### \__init__(IQ: int, JQ: int, KQ: int, LQ: int, ICQ: int)
> GROMOS impdihedral for a single pair


* **Parameters**

    
    * **IQ** (*int*)


    * **JQ** (*int*)


    * **KQ** (*int*)


    * **LQ** (*int*) – IQ,JQ,KQ,LQ: atom sequence numbers of atoms forming an improper dihedral


    * **ICQ** (*int*) – improper dihedral type code



#### to_string()

### _class_ pygromos.files.blocks.topology_blocks.impdihedralh_type(IQH: int, JQH: int, KQH: int, LQH: int, ICQH: int)
Bases: `pygromos.files.blocks._general_blocks._generic_field`


#### \__init__(IQH: int, JQH: int, KQH: int, LQH: int, ICQH: int)
> GROMOS impdihedralH for a single pair


* **Parameters**

    
    * **IQH** (*int*)


    * **JQH** (*int*)


    * **KQH** (*int*)


    * **LQH** (*int*) – IQH,JQH,KQH,LQH: atom sequence numbers of atoms forming an improper dihedral


    * **ICQH** (*int*) – improper dihedral type code



#### to_string()

### _class_ pygromos.files.blocks.topology_blocks.impdihedraltype_type(CQ: float, Q0: float)
Bases: `pygromos.files.blocks._general_blocks._generic_field`


#### \__init__(CQ: float, Q0: float)
> GROMOS impdihedraltype for a single pair


* **Parameters**

    
    * **CQ** (*float*) – force constant of improper dihedral per degrees square


    * **Q0** (*float*) – improper dihedral angle at minimum energy in degrees



#### to_string()

### _class_ pygromos.files.blocks.topology_blocks.improper_dihedral_type(ICQ: int, CQ: float, Q0: float, group_type: str, special_number: int)
Bases: `pygromos.files.blocks._general_blocks._generic_field`


#### \__init__(ICQ: int, CQ: float, Q0: float, group_type: str, special_number: int)
> GROMOS improper (harmonic) dihedral angle parameters


* **Parameters**

    
    * **ICQ** (*int*) – Improper dihedral-angle type code


    * **CQ** (*float*) – Force constant


    * **Q0** (*float*) – Ideal improper dihedral angle


    * **group_type** (*str*) – Example Useage


    * **special_number** (*int*) – No Idea @Todo: refactor later!



#### to_string()

### _class_ pygromos.files.blocks.topology_blocks.ljexception_type(AT1: int, AT2: int, C12: float, C6: float)
Bases: `pygromos.files.blocks._general_blocks._generic_field`


#### \__init__(AT1: int, AT2: int, C12: float, C6: float)
> GROMOS LJ exception pair


#### to_string()

### _class_ pygromos.files.blocks.topology_blocks.ljparameters_type(IAC: int, JAC: int, C12: float, C6: float, CS12: float, CS6: float)
Bases: `pygromos.files.blocks._general_blocks._generic_field`


#### \__init__(IAC: int, JAC: int, C12: float, C6: float, CS12: float, CS6: float)
> GROMOS LJ parameter pair


#### to_string()

### _class_ pygromos.files.blocks.topology_blocks.mixed_atom_lj_pair_type(IACI: int, IACJ: int, C6: float, C12_1: float, C12_2: float, C12_3: float)
Bases: `pygromos.files.blocks._general_blocks._generic_field`


#### to_string()

### pygromos.files.blocks.topology_blocks.pertubation_lam_state()
alias of `pygromos.files.blocks.topology_blocks.pertubationLamState`


### _class_ pygromos.files.blocks.topology_blocks.single_atom_lj_pair_type(IAC: int, TYPE: str, C6: float, C12_1: float, C12_2: float, C12_3: float, CS6: float, CS12: float, LJ14PAIR: Iterable[float])
Bases: `pygromos.files.blocks._general_blocks._generic_field`


#### \__init__(IAC: int, TYPE: str, C6: float, C12_1: float, C12_2: float, C12_3: float, CS6: float, CS12: float, LJ14PAIR: Iterable[float])

* **Parameters**

    
    * **IAC** (*int*) – vander wals type index.rst


    * **TYPE** (*str*) – atom type


    * **C6** (*float*) – square-root of C6


    * **C12_2** (*float*) – square-root of C6


    * **C12_3** (*float*) – square-root of C6


    * **CS6** (*float*) – No Idea @Todo: refactor later!


    * **CS12** (*float*) – No Idea @Todo: refactor later!


    * **LJ14PAIR** (*Iterable[float]*) – No Idea @Todo: refactor later!



#### to_string()

### _class_ pygromos.files.blocks.topology_blocks.soluteatom_type(ATNM: int, MRES: int, PANM: str, IAC: int, MASS: float, CG: float, CGC: int, INE: int, INEvalues: List[int], INE14: int, INE14values: List[int])
Bases: `pygromos.files.blocks._general_blocks._generic_field`


#### \__init__(ATNM: int, MRES: int, PANM: str, IAC: int, MASS: float, CG: float, CGC: int, INE: int, INEvalues: List[int], INE14: int, INE14values: List[int])
soluteatom_type


* **Parameters**

    
    * **ATNM** (*int*) – [description]


    * **MRES** (*int*) – [description]


    * **PANM** (*str*) – [description]


    * **IAC** (*int*) – [description]


    * **MASS** (*float*) – [description]


    * **CG** (*float*) – [description]


    * **CGC** (*int*) – [description]


    * **INE** (*int*) – [description]


    * **INEvalues** (*[type]*) – [description]


    * **INE14** (*int*) – [description]


    * **INE14values** (*[type]*) – [description]



#### to_string()

### _class_ pygromos.files.blocks.topology_blocks.solventatom_type(I: int, ANMS: str, IACS: int, MASS: float, CGS: float)
Bases: `pygromos.files.blocks._general_blocks._generic_field`


#### \__init__(I: int, ANMS: str, IACS: int, MASS: float, CGS: float)
> GROMOS solventatom line


#### to_string()

### _class_ pygromos.files.blocks.topology_blocks.solventconstr_type(ICONS: int, JCONS: int, CONS: float)
Bases: `pygromos.files.blocks._general_blocks._generic_field`


#### \__init__(ICONS: int, JCONS: int, CONS: float)
> GROMOS SOLVENTCONSTR entry


#### to_string()

### _class_ pygromos.files.blocks.topology_blocks.special_atom_lj_pair_type(c)
Bases: `pygromos.files.blocks._general_blocks._generic_field`


#### \__init__(c)

* **Parameters**

    
    * **IAC** (*int*) – vander wals type index.rst


    * **TYPE** (*str*) – atom type


    * **C6** (*float*) – square-root of C6


    * **C12_2** (*float*) – square-root of C6


    * **C12_3** (*float*) – square-root of C6


    * **CS6** (*float*) – No Idea @Todo: refactor later!


    * **CS12** (*float*) – No Idea @Todo: refactor later!


    * **LJ14PAIR** (*Iterable[float]*) – No Idea @Todo: refactor later!



#### to_string()

### _class_ pygromos.files.blocks.topology_blocks.top_bond_type(IB: int, JB: int, ICB: int)
Bases: `pygromos.files.blocks._general_blocks._generic_field`


#### \__init__(IB: int, JB: int, ICB: int)
> GROMOS bond definition for a single pair


* **Parameters**

    
    * **IB** (*int*) – Atom number i


    * **JB** (*int*) – Atom number j


    * **ICB** (*int*) – Bond type code



#### to_string()

### _class_ pygromos.files.blocks.topology_blocks.top_dihedral_type(IP: int, JP: int, KP: int, LP: int, ICP: int)
Bases: `pygromos.files.blocks._general_blocks._generic_field`


#### \__init__(IP: int, JP: int, KP: int, LP: int, ICP: int)
> GROMOS dihedral for a single pair


#### to_string()

### _class_ pygromos.files.blocks.topology_blocks.torsdihedraltype_type(CP: float, PD: float, NP: int)
Bases: `pygromos.files.blocks._general_blocks._generic_field`


#### \__init__(CP: float, PD: float, NP: int)
> GROMOS dihedraltype for a single pair


* **Parameters**

    
    * **CP** (*float*) – force constant


    * **PD** (*float*) – phase-shift angle


    * **NP** (*int*) – multiplicity



#### to_string()
## Module contents


### _class_ pygromos.files.blocks.all_blocks_class()
Bases: `object`

This class is a very simplistic block library containing all present gromos blocks.


#### get_all_blocks()
