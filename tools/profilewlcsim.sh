#/bin/bash
#sleep 30

rm gmon.out
echo "Compile"

cd code
# compile with mpi's fortran compiler
mpifort -c mersenne_twister.f90
mpifort -fbounds-check -O5 -g -pg  mersenne_twister.o SIMcode/* DATAcode/* MCcode/* -o MCparrll_out 
cd ..
mv code/MCparrll_out .
mv data/* trash/


echo "Now run"
# now run the output
# --prefix used to avoid changing path
mpirun --prefix ~/openmpi/ -np 1 MCparrll_out

gprof -l MCparrll_out

# cd code
#gfortran -O3 -fbounds-check -o wlcsim SIMcode/* BDcode/* DATAcode/* MCcode/*
# cd ..
# mv code/wlcsim .
# mv data/* trash/.
# ./wlcsim
