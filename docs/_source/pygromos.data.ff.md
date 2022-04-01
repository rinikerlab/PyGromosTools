---
date: '2022-04-01T13:55:18.499Z'
docname: _source/pygromos.data.ff
images: {}
path: /source-pygromos-data-ff
title: pygromos.data.ff package
---

# pygromos.data.ff package

## Subpackages


* [pygromos.data.ff.Gromos2016H66 package]()


    * [Module contents](#module-pygromos.data.ff.Gromos2016H66)


* [pygromos.data.ff.Gromos54A7 package]()


    * [Module contents](#module-pygromos.data.ff.Gromos54A7)


* [pygromos.data.ff.MixHexGromos54A7 package]()


    * [Module contents](#module-pygromos.data.ff.MixHexGromos54A7)


## Submodules

## pygromos.data.ff.addAtomType module


### pygromos.data.ff.addAtomType.addAtomType(atomTypes, N)

* **Parameters**

    
    * **atomTypes**


    * **N**



### pygromos.data.ff.addAtomType.check(atomTypes)

* **Parameters**

    **atomTypes**



### pygromos.data.ff.addAtomType.getAtomType(block, row)

* **Parameters**

    
    * **block**


    * **row**



### pygromos.data.ff.addAtomType.getBlock(lines, row)

* **Parameters**

    
    * **lines**


    * **row**



### pygromos.data.ff.addAtomType.main()

### pygromos.data.ff.addAtomType.readBlock(block)

* **Parameters**

    **block**



### pygromos.data.ff.addAtomType.readFile(fileName)

* **Parameters**

    **fileName**



### pygromos.data.ff.addAtomType.writeAtomTypes(atomTypes, out)

* **Parameters**

    
    * **atomTypes**


    * **out**



### pygromos.data.ff.addAtomType.writeLine(line, file)

* **Parameters**

    
    * **line**


    * **file**



### pygromos.data.ff.addAtomType.writeOut(row1, row2, atomTypes, fileName, outName)

* **Parameters**

    
    * **row1**


    * **row2**


    * **atomTypes**


    * **fileName**


    * **outName**


## pygromos.data.ff.atomType module


### _class_ pygromos.data.ff.atomType.AtomType(atomNum, atomName, sqrtC06, sqrtC12_1, sqrtC12_2, sqrtC12_3, sqrtC06NB, sqrtC12NB, matrix)
Bases: `object`


#### \__init__(atomNum, atomName, sqrtC06, sqrtC12_1, sqrtC12_2, sqrtC12_3, sqrtC06NB, sqrtC12NB, matrix)

* **Parameters**

    
    * **atomNum**


    * **atomName**


    * **sqrtC06**


    * **sqrtC12_1**


    * **sqrtC12_2**


    * **sqrtC12_3**


    * **sqrtC06NB**


    * **sqrtC12NB**


    * **matrix**



#### addToMatrix(n)

* **Parameters**

    **n**


## Module contents

Here we have additional files that might be usefull.
