---
date: '2022-04-01T13:55:18.499Z'
docname: Examples/index
images: {}
path: /examples-index
title: Examples
---

# Examples


* [TI Calculation]()


    * [Imports](#Imports)


    * [Input files](#Input-files)


    * [Vacuum Simulation](#Vacuum-Simulation)


        * [Direction A->B](#Direction-A->B)


            * [Setup:](#Setup:)


        * [RUN Emin](#RUN-Emin)


    * [RUN Test SD EQ](#RUN-Test-SD-EQ)


    * [Further Analysis:](#Further-Analysis:)


    * [Lambda Sampling](#Lambda-Sampling)


        * [Setup again](#Setup-again)


    * [Submission](#Submission)


* [Amber2Gromos Example]()


    * [create top + cnf for aniline](#create-top-+-cnf-for-aniline)


    * [create top + cnf for aniline, solvated in CHCl3](#create-top-+-cnf-for-aniline,-solvated-in-CHCl3)


    * [create top + cnf for aniline, solvated in TIP3P water](#create-top-+-cnf-for-aniline,-solvated-in-TIP3P-water)


* [Usage of the PyGromosTools on HPC-cluster]()


    * [Initialization](#Initialization)


    * [Energy minimization](#Energy-minimization)


        * [Print out the new file paths](#Print-out-the-new-file-paths)


        * [Visualize the minimized system](#Visualize-the-minimized-system)


    * [Equilibration followed by subsequent production runs](#Equilibration-followed-by-subsequent-production-runs)


        * [Visualize the last configuration](#Visualize-the-last-configuration)


* [QM/MM in Gromos]()


    * [Support for QMMM functionality in `GROMOS` input files](#Support-for-QMMM-functionality-in-GROMOS-input-files)


        * [Load an `imd` file containing a QMMM block](#Load-an-imd-file-containing-a-QMMM-block)


        * [Print out different sections of the QMMM block](#Print-out-different-sections-of-the-QMMM-block)


        * [Change a block value and print again](#Change-a-block-value-and-print-again)


    * [Directly manipulate a QMMM specification file](#Directly-manipulate-a-QMMM-specification-file)


    * [Print out and change some blocks in the QMMM specification file](#Print-out-and-change-some-blocks-in-the-QMMM-specification-file)


        * [Title block](#Title-block)


        * [QMZONE block](#QMZONE-block)


        * [QMUNIT block](#QMUNIT-block)


        * [XTBELEMENTS block](#XTBELEMENTS-block)


        * [A helper function that returns all QM engines specified in the QM/MM specification file](#A-helper-function-that-returns-all-QM-engines-specified-in-the-QM/MM-specification-file)


        * [Store your QMMM specification file with all your other simulation files in a `Gromos_System` object](#Store-your-QMMM-specification-file-with-all-your-other-simulation-files-in-a-Gromos_System-object)


    * [Run QM/MM Simulations](#Run-QM/MM-Simulations)


* [Calculation of free energy of evaporization]()


* [OpenFF Automatic Parametrization]()


* [Calculation of Self Solvation Free Energy]()


    * [Do Imports](#Do-Imports)


    * [Choose Molecule to run calculation for](#Choose-Molecule-to-run-calculation-for)


        * [create the gromos_system from a smile and get the number of atoms](#create-the-gromos_system-from-a-smile-and-get-the-number-of-atoms)


    * [Create The Solvation free energy calculation system](#Create-The-Solvation-free-energy-calculation-system)


        * [Create Liquid](#Create-Liquid)


        * [Minimize Liquid](#Minimize-Liquid)


        * [Change the number of cores for longer runs](#Change-the-number-of-cores-for-longer-runs)


    * [Equilibrate System](#Equilibrate-System)


    * [Do TI calculation](#Do-TI-calculation)


    * [Read out results](#Read-out-results)


* [Developer Examples]()


    * [Implemented pyGromosPP programs]()


        * [Ran_Box](#Ran_Box)


    * [PyGromos File Examples]()


        * [IMD - Simulation Paramter File](#IMD---Simulation-Paramter-File)


        * [CNF - Coordinate File](#CNF---Coordinate-File)


            * [Deleting Residues](#Deleting-Residues)


            * [Generate position Restraint Files](#Generate-position-Restraint-Files)


            * [Compact all](#Compact-all)


        * [TOP - Topology File](#TOP---Topology-File)


        * [TRC - Coordinate Trajectory File](#TRC---Coordinate-Trajectory-File)


        * [Other Files](#Other-Files)


            * [MTB - topology building block file](#MTB---topology-building-block-file)


            * [IFP - topology parameter file](#IFP---topology-parameter-file)


            * [disres - distance restraint file](#disres---distance-restraint-file)


        * [PTP - Pertubations for free energy calculations](#PTP---Pertubations-for-free-energy-calculations)


            * [defining some state types for later use :)](#defining-some-state-types-for-later-use-:))


            * [Read in ptp file](#Read-in-ptp-file)


            * [delete full state](#delete-full-state)


            * [delete specific atoms](#delete-specific-atoms)


            * [Add atom or state or overwrite atominformation (except atom.NR)](#Add-atom-or-state-or-overwrite-atominformation-(except-atom.NR))


    * [SD Simulations with Submission System]()


        * [Imports](#Imports)


        * [Input files](#Input-files)


        * [RUN Emin](#RUN-Emin)


        * [RUN SD Simulation](#RUN-SD-Simulation)


        * [Further Analysis:](#Further-Analysis:)


            * [Coordinate Analysis](#Coordinate-Analysis)


            * [Energy Analysis](#Energy-Analysis)


    * [Gromos Trajectory evaluation with Pygromos and Pandas]()


        * [Example file for the evaluation of GROMOS trajectory files in pygromos](#Example-file-for-the-evaluation-of-GROMOS-trajectory-files-in-pygromos)


        * [1) TRC](#1)-TRC)


            * [1.1) TRC import](#1.1)-TRC-import)


            * [1.2) TRC file handling](#1.2)-TRC-file-handling)


            * [1.3) TRC Calculate](#1.3)-TRC-Calculate)


            * [1.4) TRC Visualization](#1.4)-TRC-Visualization)


        * [2) TRE](#2)-TRE)


            * [2.1) Tre import and structure](#2.1)-Tre-import-and-structure)


            * [2.2) Common Tre functions](#2.2)-Common-Tre-functions)


            * [$\lambda$-Sampling & TREs](#\lambda-Sampling-&-TREs)


            * [EDS in TREs](#EDS-in-TREs)


        * [Concatenate and Copy multiple Trajectories](#Concatenate-and-Copy-multiple-Trajectories)


    * [Gromos Tutorial Pipeline - OLD]()


        * [Build initial files](#Build-initial-files)


            * [generate Topology](#generate-Topology)


            * [generate coordinates](#generate-coordinates)


            * [energy minimization - Vacuum](#energy-minimization---Vacuum)


        * [Solvatistation and Solvent Energy Minimization](#Solvatistation-and-Solvent-Energy-Minimization)


            * [build box system](#build-box-system)


            * [Add Ions](#Add-Ions)


            * [Energy Minimization BOX](#Energy-Minimization-BOX)


        * [Simulation](#Simulation)


            * [Equilibration NVP](#Equilibration-NVP)


            * [MD NVP](#MD-NVP)


        * [Analysis](#Analysis)


    * [Submission Systems]()


        * [Local Submission](#Local-Submission)


        * [LSF Submission](#LSF-Submission)


            * [Queue Checking:](#Queue-Checking:)


            * [Submission:](#Submission:)


            * [Submitting multiple jobs](#Submitting-multiple-jobs)


            * [Killing a jobs](#Killing-a-jobs)


    * [Topology creation and modifications with PyGromosTools]()


        * [Automatic System generation and combining of topologies](#Automatic-System-generation-and-combining-of-topologies)


        * [Manual topology modifications](#Manual-topology-modifications)
