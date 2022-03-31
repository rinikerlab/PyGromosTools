#PATH
GXXPATH=your/path/to/gromosXX
GXXBUILDPATH=you/want/gromos/here/build_gromosXX

#MODULES:
module load open_mpi/1.6.5
module load gsl
module load gcc
module load fftw
module load openblas/0.2.8_seq

#DO
cd ${GXXPATH}
mkdir LINUX -p
cd LINUX

../configure --prefix=${BUILDXX} --disable-debug  --enable-mpi --disable-openmp CXX=/cluster/apps/openmpi/1.6.5/x86_64/gcc_4.8.2/bin/mpiCC CC=/cluster/apps/openmpi/1.6.5/x86_64/gcc_4.8.2/bin/mpiCC
