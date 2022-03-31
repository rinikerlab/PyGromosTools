#!/usr/bin/env bash

#PATHS
SRCPATH=your/path/to/gromosXX
buildPATH=you/want/gromos/here/build_gromosXX

#DO
rm -r ${SRCPATH}/LINUX/*
rmdir ${SRCPATH}/LINUX
rm -r ${SRCPATH}/autom4te.cache

rm -r ${buildPATH}
