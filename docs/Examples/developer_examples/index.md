---
date: '2022-04-01T13:55:18.499Z'
docname: Examples/developer_examples/index
images: {}
path: /examples-developer-examples-index
title: Developer Examples
---

# Developer Examples


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


            * [Save as .h5](#Save-as-.h5)


        * [1.3) TRC Calculate](#1.3)-TRC-Calculate)


        * [1.4) TRC Visualization](#1.4)-TRC-Visualization)


            * [Do some changes and see how the visualization changes](#Do-some-changes-and-see-how-the-visualization-changes)


            * [let’s find hydrogenbonds](#let’s-find-hydrogenbonds)


            * [add surface](#add-surface)


            * [remove surface](#remove-surface)


            * [add transperent cartoon](#add-transperent-cartoon)


            * [show distance between atom pair](#show-distance-between-atom-pair)


            * [once you are happy with the result render and download image](#once-you-are-happy-with-the-result-render-and-download-image)


            * [Go crazy and make a movie to show off](#Go-crazy-and-make-a-movie-to-show-off)


    * [2) TRE](#2)-TRE)


        * [2.1) Tre import and structure](#2.1)-Tre-import-and-structure)


        * [2.2) Common Tre functions](#2.2)-Common-Tre-functions)


        * [$\lambda$-Sampling & TREs](#\lambda-Sampling-&-TREs)


        * [EDS in TREs](#EDS-in-TREs)


    * [Concatenate and Copy multiple Trajectories](#Concatenate-and-Copy-multiple-Trajectories)


* [Gromos Tutorial Pipeline - OLD]()


    * [Build initial files](#Build-initial-files)


        * [generate Topology](#generate-Topology)


            * [build single topologies](#build-single-topologies)


            * [combine topology](#combine-topology)


        * [generate coordinates](#generate-coordinates)


            * [add hydrogens](#add-hydrogens)


            * [cnf to pdb](#cnf-to-pdb)


        * [energy minimization - Vacuum](#energy-minimization---Vacuum)


    * [Solvatistation and Solvent Energy Minimization](#Solvatistation-and-Solvent-Energy-Minimization)


        * [build box system](#build-box-system)


            * [to pdb](#to-pdb)


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
