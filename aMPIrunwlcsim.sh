#/bin/bash

make dataclean
make
mpirun -np 10 wlcsim.exe -i input/params
