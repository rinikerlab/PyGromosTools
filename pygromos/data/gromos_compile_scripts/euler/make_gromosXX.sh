#PATH
GXXPATH=your/path/to/gromosXX

#MODULES
module load open_mpi/1.6.5
module load gsl
module load gcc
module load fftw
module load openblas/0.2.8_seq

#DO
cd ${GXXPATH}
mkdir LINUX
cd LINUX
make -j8
make install
