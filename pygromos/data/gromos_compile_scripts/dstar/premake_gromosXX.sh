#!/usr/bin/env bash

#PATHS
SRCPATH=your/path/to/gromosXX
buildPATH=you/want/gromos/here/build_gromosXX

gslLIBRARY=/opt/dstar/progs/gcc-4.8.5/gsl-2.1/
fftwLIBRARY=/opt/dstar/progs/gcc-4.8.5/fftw-3.3.4/

CXXFLAG=/opt/dstar/progs/gcc-4.8.5/openmpi-1.10.2/bin/mpiCC
CCFLAG=/opt/dstar/progs/gcc-4.8.5/openmpi-1.10.2/bin/mpicc

#DO
cd $SRCPATH
echo "${PWD}"
mkdir -p $buildPATH
./Config.sh
mkdir LINUX -p
cd LINUX
CXX=${CXXFLAG}
CC=${CCFLAG}
#OPTIONS
../configure --disable-shared --enable-debug --with-gsl=/opt/dstar/progs/gcc-4.8.5/gsl-2.1/ --with-fftw=/opt/dstar/progs/gcc-4.8.5/fftw-3.3.4/ --disable-openmp --enable-mpi CXX=/opt/dstar/progs/gcc-4.8.5/openmpi-1.10.2/bin/mpiCC CC=/opt/dstar/progs/gcc-4.8.5/openmpi-1.10.2/bin/mpicc --prefix=${buildPATH}
