#/bin/bash

echo "Compile ..."
gfortran -c -fbounds-check -Wall -W -fmax-errors=5 -O1 ../util/binning.f95 test_binning.f90
gfortran *.o -o testBin

echo "now Run ..."
./testBin

rm *.o


