#!/bin/bash
set -o pipefail
set -eu

cd code
gfortran -c SIMcode/mt19937.f90
gfortran -o ../wlcsim SIMcode/* BDcode/* DATAcode/* MCcode/*
cd ..
mkdir -p data trash
mv data/* trash || true
./wlcsim
