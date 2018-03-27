#/bin/bash

make dataclean
make
/usr/bin/mpirun -np 10 wlcsim.exe -i input/params
