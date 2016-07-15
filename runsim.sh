#!/bin/bash
set -o pipefail
set -eu

cd code
gfortran -cpp -c SIMcode/mt19937.f90 -J. -I.
gfortran -o ../wlcsim SIMcode/* BDcode/* DATAcode/* MCcode/* -I.
cd ..
mkdir -p data trash
mv data/* trash || true
./wlcsim
