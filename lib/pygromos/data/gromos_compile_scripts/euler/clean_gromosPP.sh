SRCPATH=your/path/to/gromosXX
buildPATH=you/want/gromos/here/build_gromosXX


cd ${SRCPATH}
rm LINUX/* -r
rmdir LINUX
rm ${buildPATH}/* -r
rmdir ${buildPATH}
