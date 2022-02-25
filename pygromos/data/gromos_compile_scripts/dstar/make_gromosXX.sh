#!/usr/bin/env bash
SRCPATH=your/path/to/gromosXX
buildPATH=you/want/gromos/here/build_gromosXX

cd $SRCPATH
cd LINUX
make -j8
make install # -d
