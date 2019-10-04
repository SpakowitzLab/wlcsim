#! /bin/bash
## Number of nodes
# SBATCH -N 1
## Number of cores
# SBATCH -n 2

echo This is job $SLURM_JOB_ID

# conda env
source activate wlcsim
echo `which python`

rm data/*
rm wlcsim.exe

make

mkdir allData

for iter in {1..50..1}
  do
      mpirun --prefix /usr/mpi/gcc/openmpi-1.10.3a1/ wlcsim.exe 
      newfile="data$iter"
      mv data $newfile
      mv $newfile allData
      echo "done with iteration $iter out of 50"
      mkdir data
  done

