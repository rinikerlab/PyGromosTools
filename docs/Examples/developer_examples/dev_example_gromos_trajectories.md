---
date: '2022-04-01T13:55:18.499Z'
docname: Examples/developer_examples/dev_example_gromos_trajectories
images: {}
path: /examples-developer-examples-dev-example-gromos-trajectories
title: Gromos Trajectory evaluation with Pygromos and Pandas
---

# Gromos Trajectory evaluation with Pygromos and Pandas

## Example file for the evaluation of GROMOS trajectory files in pygromos

## 1) TRC

### 1.1) TRC import

The TRC class bridges the GROMOS TRC file based strucutres and mdtraj.Trajectories.

One can read in TRCs and use them just like mdtraj.Trajectories. One can than write them out both as .h5 as well as GROMOS .trc/trc.gz files.

If you have a function that’s generally useful, please contact the developers to possibly add it to the pygromos code to help other people :)

### 1.2) TRC file handling

#### Save as .h5

##### Save as trc

##### export last frame as cnf

##### Get selection of trajectories

##### Save out selection

### 1.3) TRC Calculate

### 1.4) TRC Visualization

#### Do some changes and see how the visualization changes

To jump to a certain frame, use the slider in the widget or select frames with the following cell:

#### let’s find hydrogenbonds

#### add surface

#### remove surface

#### add transperent cartoon

#### show distance between atom pair

#### once you are happy with the result render and download image

#### Go crazy and make a movie to show off

## 2) TRE

### 2.1) Tre import and structure

Tre files contain all energy related data (like split up energy terms, temperature, pressure, …..). In PyGromos they generally share the same block structure as other files, but all the data inside the specific timesteps is stored efficiently inside a pandas DataFrame, here called tre.database . This database offers manipulation with all pandas functions. Alternatively many common functions are provided inside the Tre class.

This class should in principle replace further usage of the gromos++ ene_ana function, since all these operation can be done efficiently on the pandas DataFrame.

We are currently working on adding more common functions to the Tre class. If you find a useful function please contact the developers so the function can be added for general usage :)

### 2.2) Common Tre functions

Tables and lists inside the database are stored in numpy arrays. For example the two temperatures from the previous example are stored in a numpy array of size 2 since it has two temperature baths

Specific values inside a tre file can also be directly accessed with numpy and pandas syntax

### $\lambda$-Sampling & TREs

### EDS in TREs

## Concatenate and Copy multiple Trajectories

Trajectories offer a wide range of additional file manipulations. Trajectory classes can be copied (deep) and added to each other to concatenate multiple small simulation pieces into one large trajectory.

In the new combined trajectory we have one long trajectory made from the two smaller ones. The length is one element shorter, since normally the last element of the first trajectory and the first element of the second trajectory is the same element. This can be controlled via the option “skip_new_0=True” in the add_traj() function which is the core of the “+” operator for trajectories. In the following line the default behavior can be seen as a smooth numbering in the TIMESTEPs.
