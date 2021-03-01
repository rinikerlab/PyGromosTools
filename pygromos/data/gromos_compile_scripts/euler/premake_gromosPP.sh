#PATH
PLSPATH=your/path/to/gromosXX
PLSBUILDPATH=you/want/gromos/here/build_gromosXX

#MODULE
module load open_mpi/1.6.5
module load gsl
module load gcc
module load fftw
module load openblas/0.2.8_seq

#DO
cd $PLSPATH
mkdir LINUX -p
cd LINUX
../configure --disable-shared --disable-debug --enable-openmp \
             --prefix=${GPPPATHBUILD} || exit 1
