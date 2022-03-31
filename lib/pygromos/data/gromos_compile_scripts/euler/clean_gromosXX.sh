GXXPATH=your/path/to/gromosXX
GXXBUILDPATH=you/want/gromos/here/build_gromosXX

cd ${GXXPATH}
rm LINUX/* -r
rmdir LINUX
rm ${GXXBUILDPATH}/* -r
rmdir ${GXXBUILDPATH}
