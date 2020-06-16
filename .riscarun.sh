#! /bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p risc

echo This is job $SLURM_JOB_ID
echo has $SLURM_NTASKS cores

# conda env
source activate wlcsim
echo `which python`

rm data/*
rm wlcsim.exe

make

#mkdir allData

#for iter in {1..50..1}
#  do
mpirun --prefix /usr/mpi/gcc/openmpi-1.10.3a1/ wlcsim.exe 
#      newfile="data$iter"
#      mv data $newfile
#      mv $newfile allData
#      echo "done with iteration $iter out of 50"
#      mkdir data
#  done

