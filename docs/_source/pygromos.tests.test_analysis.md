---
date: '2022-04-01T13:55:18.499Z'
docname: _source/pygromos.tests.test_analysis
images: {}
path: /source-pygromos-tests-test-analysis
title: pygromos.tests.test_analysis package
---

# pygromos.tests.test_analysis package

## Submodules

## pygromos.tests.test_analysis.test_error_estimate module


### _class_ pygromos.tests.test_analysis.test_error_estimate.test_ee(methodName='runTest')
Bases: `unittest.case.TestCase`


#### error_estimate_class()
alias of [`pygromos.analysis.error_estimate.error_estimator`](#pygromos.analysis.error_estimate.error_estimator)


#### test_array(_ = array([   0,    1,    2, ..., 9997, 9998, 9999]_ )

#### test_constructor()

#### test_error_estimate()
## pygromos.tests.test_analysis.test_freeEnergy module


### _class_ pygromos.tests.test_analysis.test_freeEnergy.test_BAR(methodName='runTest')
Bases: `pygromos.tests.test_analysis.test_freeEnergy.test_ZwanzigEquation`


#### feCalculation()
alias of [`pygromos.analysis.free_energy_calculation.bennetAcceptanceRatio`](#pygromos.analysis.free_energy_calculation.bennetAcceptanceRatio)


#### test_free_Energy1()

### _class_ pygromos.tests.test_analysis.test_freeEnergy.test_ZwanzigEquation(methodName='runTest')
Bases: `unittest.case.TestCase`


#### feCalculation()
alias of [`pygromos.analysis.free_energy_calculation.zwanzigEquation`](#pygromos.analysis.free_energy_calculation.zwanzigEquation)


#### test_constructor()

#### test_free_Energy1()

### _class_ pygromos.tests.test_analysis.test_freeEnergy.test_threeStateZwanzigReweighting(methodName='runTest')
Bases: `pygromos.tests.test_analysis.test_freeEnergy.test_ZwanzigEquation`


#### feCalculation()
alias of [`pygromos.analysis.free_energy_calculation.threeStateZwanzig`](#pygromos.analysis.free_energy_calculation.threeStateZwanzig)


#### test_free_Energy1()
## Module contents
