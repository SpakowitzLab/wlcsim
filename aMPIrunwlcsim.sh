#/bin/bash

make dataclean
make
MPIEXE=
if [[ -f /usr/local/bin/mpirun ]] ; then
    MPIEXE=/usr/local/bin/mpirun
else
    MPIEXE=/usr/bin/mpirun
fi
$MPIEXE -np 1 wlcsim.exe -i input/params
