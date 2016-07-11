#!/bin/bash
set -o pipefail
set -eu

gfortran -o wlcsim code/SIMcode/* code/BDcode/* code/DATAcode/* code/MCcode/* -Icode
mkdir -p data trash
mv data/* trash || true
./wlcsim
