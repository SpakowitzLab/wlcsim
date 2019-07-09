#/bin/bash

# Generate LEF and LEF dependent precalculations
cd input
    python3 makeBindPairs.py
    python3 generateSpiders.py
cd ..

# Clear old stuff
rm -f data/*
rm -f wlcsim.exe
make clean

# Compile fortran
make

# Run fortran
./wlcsim.exe -i input/params
