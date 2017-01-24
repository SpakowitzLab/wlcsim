#/bin/bash

echo "Compile"
#rm MCparrll_out
cd code
# compile with mpi's fortran compiler
mpifort -o  MCparrll_out SIMcode/*f90  SIMcode/*f95 MCcode/*f95 /usr/lib/liblapack.so -O5 
cd ..
mv code/MCparrll_out .
#mkdir -p data
#mv data/* trash/

#touch data/error
echo "Now run"
# now run the output
# --prefix used to avoid changing path
mpirun -np 3 xterm -e gdb MCparrll_out

