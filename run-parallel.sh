#!/bin/bash

set -eu
set -o pipefail

codedir=`pwd`

COUNTER=0
until [ ! -d par-run-dir.$COUNTER ]; do
    let COUNTER+=1
done
pardir="par-run-dir.$COUNTER"
echo "Running simulations in ${COUNTER}!"
nprocs=`grep -c ^processor /proc/cpuinfo`
let nprocs-=1
echo "Using $nprocs processors!"
mkdir -p "$pardir"
cd "$pardir"
for core in `seq 1 $nprocs`; do
    rundir=run.$core
    mkdir -p "$rundir"
    cp -r "${codedir}/input" "${codedir}/run_parameter.pl" "$rundir"
    cd "$rundir"
    mkdir -p data savedata
    echo "#!/bin/bash"  >> runsim.sh
    echo "rm data/*"    >> runsim.sh
    echo "../../wlcsim" >> runsim.sh
    chmod a+x runsim.sh
    ./run_parameter.pl &
    cd ..
done

