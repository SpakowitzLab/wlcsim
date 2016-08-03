#!/bin/bash
# this script runs parallel simulations using run_parameter.pl
# to vary parameters.
# to start it on several computers, simply run it on them individually
# each computer should create its own subfolder to save its results
#
# you can pass a name for the run to save it in a unique folder, otherwise
# the script will just use par-run-dir/run.$runnumber as the folder to save in

set -eu
set -o pipefail

codedir=`pwd`
compname=`hostname`
#nprocs=`grep -c ^processor /proc/cpuinfo`
nprocs=8
commit=`git rev-parse HEAD`
if [ $# -eq 0 ]; then
    run_name="par-run-dir/run"
else
    run_name="par-run-dir/$1"
fi

# first, make sure we've build the latest version of the executable possible
make clean && make

# now, get the next available directory name to save output in given
# run name
COUNTER=0
until [ ! -d ${run_name}.$COUNTER/$compname ]; do
    let COUNTER+=1
done
pardir="${run_name}.$COUNTER/$compname"
echo "Running simulations in ${pardir}!"
mkdir -p "$pardir"
echo "Logging current commit hash: $commit"
echo "$commit" > "$pardir/commit_hash"
let nprocs-=1
echo "Using $nprocs processors!"
cd "$pardir"
for core in `seq 1 $nprocs`; do
    rundir=core.$core
    mkdir -p "$rundir"
    cp -r "${codedir}/input" "${codedir}/run_parameter.pl" "${codedir}/wlcsim.exe" "$rundir"
    cd "$rundir"
    mkdir -p data savedata
    echo "#!/bin/bash"  >> runsim.sh
    echo "rm data/* || true"    >> runsim.sh
    echo "./wlcsim.exe >/dev/null 2>&1" >> runsim.sh
    chmod a+x runsim.sh
    ./run_parameter.pl &
    cd ..
done

