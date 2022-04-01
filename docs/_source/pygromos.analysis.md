---
date: '2022-04-01T13:55:18.499Z'
docname: _source/pygromos.analysis
images: {}
path: /source-pygromos-analysis
title: pygromos.analysis package
---

# pygromos.analysis package

## Submodules

## pygromos.analysis.coordinate_analysis module


### pygromos.analysis.coordinate_analysis.calculate_distance(atomA: numpy.array, atomB: numpy.array)

### pygromos.analysis.coordinate_analysis.periodic_distance(vec: numpy.array, grid: numpy.array)

### pygromos.analysis.coordinate_analysis.periodic_shift(vec: numpy.array, grid: numpy.array)

### pygromos.analysis.coordinate_analysis.rms(in_values)
helper function for RMSD. Calculates the root mean square for a array of (position/velocity) arrays


* **Parameters**

    **in_values** (*np.array of np.arrays*)



* **Returns**

    root mean square



* **Return type**

    float


## pygromos.analysis.energy_analysis module


### pygromos.analysis.energy_analysis.get_Hvap(liq_totpot_energy: float, gas_totpot_energy: float, nMolecules=1, temperature=None, R=0.008314462618)

### pygromos.analysis.energy_analysis.get_density(mass: numpy.array, volume: numpy.array, atomu=1.6605390671738465)
Calculate the density


* **Parameters**

    
    * **mass** (*np.array*)


    * **volume** (*np.array*)



* **Returns**

    resulting density



* **Return type**

    np.array


## pygromos.analysis.error_estimate module

File: calculation of error estimates
Warnings: this class is WIP!

Description:

    Implementation of functions estimating the Error for calculated values using Gromos

Author: Paul Katzberger


### _class_ pygromos.analysis.error_estimate.error_estimator(values: List[float])
Bases: `object`

Class to calculate the Error Estimate as implemented in ene_ana


#### \__init__(values: List[float])
Initialize calculation
:Parameters: **values** (*List[float]*) – list of ordered values for which the error should be estimated


#### calculate_error_estimate()
Calculation of the Error Estimate as in ene_ana
:returns: **d_ee** – Error Estimate for provided list
:rtype: float


#### calculate_rmsd()
Calculate rmsd
:returns: **rmsd** – rmsd of values
:rtype: float

## pygromos.analysis.free_energy_calculation module

Free Energy Calculations:

    This module contains functions for free energy calculations
    author: Gerhard König, Benjamin Schroeder

    THIS file was pirated from here: [https://github.com/rinikerlab/Ensembler/blob/master/ensembler/analysis/freeEnergyCalculation.py](https://github.com/rinikerlab/Ensembler/blob/master/ensembler/analysis/freeEnergyCalculation.py)


### _class_ pygromos.analysis.free_energy_calculation.bar(C: float = 0.0, T: float = 298, k: float = 8.31446261815324, kT: bool = False, kJ: bool = False, kCal: bool = False, convergence_radius: float = 1e-05, max_iterations: int = 500, min_iterations: int = 1)
Bases: `pygromos.analysis.free_energy_calculation.bennetAcceptanceRatio`


#### convergence_radius(_: floa_ )

#### max_iterations(_: in_ )

#### min_iterations(_: in_ )

### _class_ pygromos.analysis.free_energy_calculation.bennetAcceptanceRatio(C: float = 0.0, T: float = 298, k: float = 8.31446261815324, kT: bool = False, kJ: bool = False, kCal: bool = False, convergence_radius: float = 1e-05, max_iterations: int = 500, min_iterations: int = 1)
Bases: `pygromos.analysis.free_energy_calculation._FreeEnergyCalculator`

> This class implements the BAR method.
> $dF = -

rac{1}{eta} \* ln(
rac{flangle(V_j-V_i+C)
angle_i}{flangle(V_i-V_j-C)
angle_j})+C$

> with :
> $ f(x) =

rac{1}{1+e^(eta x)}$ - fermi function


#### C(_ = _ )

#### T(_ = _ )

#### Vi_i(_ = Vi__ )

#### Vi_j(_ = Vi__ )

#### Vj_i(_ = Vj__ )

#### Vj_j(_ = Vj__ )

#### \__init__(C: float = 0.0, T: float = 298, k: float = 8.31446261815324, kT: bool = False, kJ: bool = False, kCal: bool = False, convergence_radius: float = 1e-05, max_iterations: int = 500, min_iterations: int = 1)
Here you can set Class wide the parameters T and k for the bennet acceptance ration (BAR) Equation
:Parameters: \* **C** (*float, optional*) – is the initial guess of the free Energy.

> 
> * **T** (*float, optional*) – Temperature in Kelvin, defaults to 398


> * **k** (*float, optional*) – boltzmann Constant, defaults to const.k\*const.Avogadro


> * **kT** (*bool, optional*) – overwrites T and k to set all results in units of $k_bT$


> * **kJ** (*bool, optional*) – overwrites k to get the Boltzman constant with units kJ/(mol\*K)


> * **kCal** (*bool, optional*) – overwrites k to get the Boltzman constant with units kcal/(mol\*K)


> * **convergence_radius** (*float, optional*) – when is the result converged? if the deviation of one to another iteration is below the convergence radius.


> * **max_iterations** (*int, optional*) – maximal number of iterations.


> * **min_iterations** (*int, optional*) – minimal number of iterations


#### _calc_bar(C: numbers.Number, Vj_i: numpy.array, Vi_i: numpy.array, Vi_j: numpy.array, Vj_j: numpy.array)
> _calc_bar

>     this function is calculating the free energy difference of two states for one iteration of the BAR method.

>         It is implemented straight forwad, but therefore not very numerical stable.


* **Parameters**

    
    * **Vi_i** (*np.array*) – potential energies of stateI while sampling stateI


    * **Vj_i** (*np.array*) – potential energies of stateJ while sampling stateI


    * **Vi_j** (*np.array*) – potential energies of stateI while sampling stateJ


    * **Vj_j** (*np.array*) – potential energies of stateJ while sampling stateJ



* **Returns**

    free energy difference



* **Return type**

    float



#### _calc_bar_mpmath(C: numbers.Number, Vj_i: numpy.array, Vi_i: numpy.array, Vi_j: numpy.array, Vj_j: numpy.array)
> _calc_bar

>     this function is calculating the free energy difference of two states for one iteration of the BAR method.

>         It is implemented straight forwad, but therefore not very numerical stable.


* **Parameters**

    
    * **Vi_i** (*np.array*) – potential energies of stateI while sampling stateI


    * **Vj_i** (*np.array*) – potential energies of stateJ while sampling stateI


    * **Vi_j** (*np.array*) – potential energies of stateI while sampling stateJ


    * **Vj_j** (*np.array*) – potential energies of stateJ while sampling stateJ



* **Returns**

    free energy difference



* **Return type**

    float



#### _calculate_optimize(Vi_i: (typing.Iterable[numbers.Number], <class 'numbers.Number'>), Vj_i: (typing.Iterable[numbers.Number], <class 'numbers.Number'>), Vi_j: (typing.Iterable[numbers.Number], <class 'numbers.Number'>), Vj_j: (typing.Iterable[numbers.Number], <class 'numbers.Number'>), C0: float = 0, verbose: bool = False)
this function is calculating the free energy difference of two states with the BAR method.

    it iterates over the _calc_bar method and determines the convergence and the result.


* **Parameters**

    
    * **Vi_i** (*np.array*) – potential energies of stateI while sampling stateI


    * **Vj_i** (*np.array*) – potential energies of stateJ while sampling stateI


    * **Vi_j** (*np.array*) – potential energies of stateI while sampling stateJ


    * **Vj_j** (*np.array*) – potential energies of stateJ while sampling stateJ



* **Returns**

    free energy difference



* **Return type**

    float



#### beta(_ = bet_ )

#### calculate(Vi_i: Iterable[numbers.Number], Vj_i: Iterable[numbers.Number], Vi_j: Iterable[numbers.Number], Vj_j: Iterable[numbers.Number], verbose: bool = False)
> calculate

>     this function is calculating the free energy difference of two states with the BAR method.


* **Parameters**

    
    * **Vi_i** (*np.array*) – potential energies of stateI while sampling stateI


    * **Vj_i** (*np.array*) – potential energies of stateJ while sampling stateI


    * **Vi_j** (*np.array*) – potential energies of stateI while sampling stateJ


    * **Vj_j** (*np.array*) – potential energies of stateJ while sampling stateJ



* **Returns**

    free energy difference



* **Return type**

    float



#### constants(_: dic_ _ = {T: 298, k: 8.31446261815324, C: <class 'numbers.Number'>_ )

#### convergence_radius(_: floa_ )

#### equation(_: sympy.core.function.Functio_ _ = (-log(exp((C - Vi_i + Vj_i)/(T\*k))) + log(exp((C + Vi_j - Vj_j)/(T\*k))))/(T\*k_ )

#### k(_ = _ )

#### max_iterations(_: in_ )

#### min_iterations(_: in_ )

#### set_parameters(C: Optional[float] = None, T: Optional[float] = None, k: Optional[float] = None)
> set_parameters setter for the parameters T and k


* **Parameters**

    
    * **T** (*float, optional*) – Temperature in Kelvin, defaults to 398


    * **k** (*float, optional*) – boltzmann Constant, defaults to const.k\*const.Avogadro


    * **C** (*float, optional*) – C is the initial guess of the free energy difference.



### _class_ pygromos.analysis.free_energy_calculation.dfEDS(kCal: bool = False, T: float = 298, k: float = 8.31446261815324, kT: bool = False, kJ: bool = False)
Bases: `pygromos.analysis.free_energy_calculation.threeStateZwanzig`


### _class_ pygromos.analysis.free_energy_calculation.threeStateZwanzig(kCal: bool = False, T: float = 298, k: float = 8.31446261815324, kT: bool = False, kJ: bool = False)
Bases: `pygromos.analysis.free_energy_calculation.zwanzigEquation`

> this class provides the implementation for the Free energy calculation with EDS.
> It calculates the free energy via the reference state.
> $dF = dF_{BR}-dF_{AR} =

rac{1}{eta} \* ( ln(langle e^{-eta \* (V_j-V_R)}
angle) - ln(langle e^{-eta \* (V_i-V_R)}
angle))$


#### T(_ = _ )

#### Vi(_ = V_ )

#### Vj(_ = V_ )

#### Vr(_ = V_ )

#### \__init__(kCal: bool = False, T: float = 298, k: float = 8.31446261815324, kT: bool = False, kJ: bool = False)
this class provides the implementation for the Free energy calculation with EDS.

    It calculates the free energy via the reference state.


* **Parameters**

    
    * **T** (*float, optional*) – Temperature in Kelvin, defaults to 398


    * **k** (*float, optional*) – boltzmann Constant, defaults to const.k\*const.Avogadro


    * **kT** (*bool, optional*) – overwrites T and k to set all results in units of $k_bT$


    * **kJ** (*bool, optional*) – overwrites k to get the Boltzman constant with units kJ/(mol\*K)


    * **kCal** (*bool, optional*) – overwrites k to get the Boltzman constant with units kcal/(mol\*K)



#### _calculate_implementation_useZwanzig(Vi: (typing.Iterable, <class 'numbers.Number'>), Vj: (typing.Iterable, <class 'numbers.Number'>), Vr: (typing.Iterable[numbers.Number], <class 'numbers.Number'>))
> calculate

>     this method calculates the zwanzig equation via the intermediate reference state R using the Zwanzig equation.
>     it directly accesses the zwanzig implmentation.


* **Parameters**

    
    * **Vi** (*np.array*) – the potential energy of stateI while sampling stateR


    * **Vj** (*np.array*) – the potential energy of stateJ while sampling stateR


    * **Vr** (*np.array*) – the potential energy of stateR while sampling stateR



* **Returns**

    free energy difference



* **Return type**

    float



#### calculate(Vi: (typing.Iterable[numbers.Number], <class 'numbers.Number'>), Vj: (typing.Iterable[numbers.Number], <class 'numbers.Number'>), Vr: (typing.Iterable[numbers.Number], <class 'numbers.Number'>))
> calculate

>     this method calculates the zwanzig equation via the intermediate reference state R using the Zwanzig equation.


* **Parameters**

    
    * **Vi** (*np.array*) – the potential energy of stateI while sampling stateR


    * **Vj** (*np.array*) – the potential energy of stateJ while sampling stateR


    * **Vr** (*np.array*) – the potential energy of stateR while sampling stateR



* **Returns**

    free energy difference



* **Return type**

    float



#### equation(_: sympy.core.function.Functio_ _ = -(log(exp(-(Vi - Vr)/(T\*k))) - log(exp(-(Vj - Vr)/(T\*k))))/(T\*k_ )

#### k(_ = _ )

### _class_ pygromos.analysis.free_energy_calculation.zwanzig(T: float = 298, k: float = 8.31446261815324, kT: bool = False, kJ: bool = False, kCal: bool = False)
Bases: `pygromos.analysis.free_energy_calculation.zwanzigEquation`


#### constants(_: dic_ )

### _class_ pygromos.analysis.free_energy_calculation.zwanzigEquation(T: float = 298, k: float = 8.31446261815324, kT: bool = False, kJ: bool = False, kCal: bool = False)
Bases: `pygromos.analysis.free_energy_calculation._FreeEnergyCalculator`

Zwanzig Equation

    This class is a nice wrapper for the zwanzig Equation.
    dF = - eta ln(langle e^(-beta(V_j-V_i))

angle)


#### T(_ = _ )

#### Vi(_ = V_ )

#### Vj(_ = V_ )

#### \__init__(T: float = 298, k: float = 8.31446261815324, kT: bool = False, kJ: bool = False, kCal: bool = False)
Here you can set Class wide the parameters T and k for the Zwanzig Equation
:Parameters: \* **T** (*float, optional*) – Temperature in Kelvin, defaults to 398

> 
> * **k** (*float, optional*) – boltzmann Constant, defaults to const.k\*const.Avogadro


> * **kT** (*bool, optional*) – overwrites T and k to set all results in units of $k_bT$


> * **kJ** (*bool, optional*) – overwrites k to get the Boltzman constant with units kJ/(mol\*K)


> * **kCal** (*bool, optional*) – overwrites k to get the Boltzman constant with units kcal/(mol\*K)


#### _calculate_efficient(Vi: (typing.Iterable, <class 'numbers.Number'>), Vj: (typing.Iterable, <class 'numbers.Number'>))
> _calculate_efficient
> Calculate a free energy difference with the Zwanzig equation (aka exponential formula or thermodynamic perturbation).
> The initial state of the free energy difference is denoted as 0, the final state is called 1.
> The function expects two arrays of size n with potential energies. The first array, u00, contains the potential energies of a set
> of Boltzmann-weighted conformations of an MD or MC trajectory of the initial state, analyzed with the Hamiltonian of the
> initial state. The second array, u01 , contains the potential energies of a trajectory of the initial state that was
> analyzed with the potential energy function of the final state. The variable kT expects the product of the Boltzmann
> constant with the temperature that was used to generate the trajectory in the respective units of the potential energies.
> This is an efficient more overflow robust implementation of the Zwanzig Equation.
> @Author: Gerhard König
> See Zwanzig, R. W. J. Chem. Phys. 1954, 22, 1420-1426. doi:10.1063/1.1740409
> $dF =

rac{1}{eta} \* ln(langlee^{-eta \* (V_j-V_i)}
angle)$

> Vi

>     Potential energies of state I

> Vj

>     Potential energies of state J

> float

>     free energy difference


#### _calculate_implementation_bruteForce(Vi: (typing.Iterable, <class 'numbers.Number'>), Vj: (typing.Iterable, <class 'numbers.Number'>))
> _calculate_implementation_bruteForce

>     This is a plain implementation of the zwanzig equation. It is not very numerical robust


* **Parameters**

    
    * **Vi**


    * **Vj**



* **Parameters**

    
    * **Vi** (*np.array*) – Potential energies of state I


    * **Vj** (*np.array*) – Potential energies of state J



* **Returns**

    free energy difference



* **Return type**

    float



#### _calculate_meanEfficient(Vi: (typing.Iterable, <class 'numbers.Number'>), Vj: (typing.Iterable, <class 'numbers.Number'>))
Calculate a free energy difference with the Zwanzig equation (aka exponential formula or thermodynamic perturbation).
The initial state of the free energy difference is denoted as 0, the final state is called 1.
The function expects two arrays of size n with potential energies. The first array, u00, contains the potential energies of a set
of Boltzmann-weighted conformations of an MD or MC trajectory of the initial state, analyzed with the Hamiltonian of the
initial state. The second array, u01 , contains the potential energies of a trajectory of the initial state that was
analyzed with the potential energy function of the final state. The variable kT expects the product of the Boltzmann
constant with the temperature that was used to generate the trajectory in the respective units of the potential energies.
This is an efficient more overflow robust implementation of the Zwanzig Equation.
@Author: Gerhard König
See Zwanzig, R. W. J. Chem. Phys. 1954, 22, 1420-1426. doi:10.1063/1.1740409

> _calculate_implementation_bruteForce

>     This is a plain implementation of the zwanzig equation. It is not very numerical robust


* **Parameters**

    
    * **Vi** (*np.array*) – Potential energies of state I


    * **Vj** (*np.array*) – Potential energies of state J



* **Returns**

    free energy difference



* **Return type**

    float



#### _calculate_mpmath(Vi: (typing.Iterable, <class 'numbers.Number'>), Vj: (typing.Iterable, <class 'numbers.Number'>))
> implementation of zwanzig with mpmath package, another way of having a robust variant,
> but this one is very close to the initial equation thanks to the mpmath package.
> $dF =

rac{1}{eta} \* ln(langlee^{-eta \* (V_j-V_i)}
angle)$

> Vi

>     Potential energies of state I

> Vj

>     Potential energies of state J

> float

>     free energy difference


#### calculate(Vi: (typing.Iterable[numbers.Number], <class 'numbers.Number'>), Vj: (typing.Iterable[numbers.Number], <class 'numbers.Number'>))
zwanzig
Calculate a free energy difference with the Zwanzig equation (aka exponential formula or thermodynamic perturbation).
The initial state of the free energy difference is denoted as 0, the final state is called 1.
The function expects two arrays of size n with potential energies. The first array, u00, contains the potential energies of a set
of Boltzmann-weighted conformations of an MD or MC trajectory of the initial state, analyzed with the Hamiltonian of the
initial state. The second array, u01 , contains the potential energies of a trajectory of the initial state that was
analyzed with the potential energy function of the final state. The variable kT expects the product of the Boltzmann
constant with the temperature that was used to generate the trajectory in the respective units of the potential energies.
See Zwanzig, R. W. J. Chem. Phys. 1954, 22, 1420-1426. doi:10.1063/1.1740409
:Parameters: \* **Vi** (*np.array*) – Potential energies of state I

> 
> * **Vj** (*np.array*) – Potential energies of state J


* **Returns**

    free energy difference



* **Return type**

    float



#### constants(_: dic_ )

#### equation(_: sympy.core.function.Functio_ _ = -log(exp(-(-Vi + Vj)/(T\*k)))/(T\*k_ )

#### k(_ = _ )

#### set_parameters(T: Optional[float] = None, k: Optional[float] = None)
> set_parameters setter for the parameters T and k


* **Parameters**

    
    * **T** (*float, optional*) – Temperature in Kelvin, defaults to 398


    * **k** (*float, optional*) – boltzmann Constant, defaults to const.k\*const.Avogadro


## Module contents

This is a future folder
