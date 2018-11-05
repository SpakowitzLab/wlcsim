#/bin/bash
python unstructure.py src/wlcsim/wlcsim.f03
sed -i -e "s/wlc_d%/wlc__/g" src/wlcsim/wlcsim.f03
sed -i -e "s/, wlc_d//g" src/wlcsim/wlcsim.f03
sed -i -e "s/,wlc_d//g" src/wlcsim/wlcsim.f03
sed -i -e "s/wlc_d,//g" src/wlcsim/wlcsim.f03
sed -i -e "s/wlc_d ,//g" src/wlcsim/wlcsim.f03
sed -i -e "s/wlc_d//g" src/wlcsim/wlcsim.f03
sed -i -e "s/wlc__/wlc_/g" src/wlcsim/wlcsim.f03
sed -i -e "s/wlcsim_data,//g" src/wlcsim/wlcsim.f03
sed -i -e "s/wlcsim_data ,//g" src/wlcsim/wlcsim.f03
sed -i -e "s/,wlcsim_data//g" src/wlcsim/wlcsim.f03
sed -i -e "s/, wlcsim_data//g" src/wlcsim/wlcsim.f03
sed -i -e "s/use params, only: wlcsim_data//g" src/wlcsim/wlcsim.f03
python unstructure.py src/wlcsim/params.f03
sed -i -e "s/wlc_d%/wlc__/g" src/wlcsim/params.f03
sed -i -e "s/, wlc_d//g" src/wlcsim/params.f03
sed -i -e "s/,wlc_d//g" src/wlcsim/params.f03
sed -i -e "s/wlc_d,//g" src/wlcsim/params.f03
sed -i -e "s/wlc_d ,//g" src/wlcsim/params.f03
sed -i -e "s/wlc_d//g" src/wlcsim/params.f03
sed -i -e "s/wlc__/wlc_/g" src/wlcsim/params.f03
sed -i -e "s/wlcsim_data,//g" src/wlcsim/params.f03
sed -i -e "s/wlcsim_data ,//g" src/wlcsim/params.f03
sed -i -e "s/,wlcsim_data//g" src/wlcsim/params.f03
sed -i -e "s/, wlcsim_data//g" src/wlcsim/params.f03
sed -i -e "s/use params, only: wlcsim_data//g" src/wlcsim/params.f03
python unstructure.py src/wlcsim/initcond.f03
sed -i -e "s/wlc_d%/wlc__/g" src/wlcsim/initcond.f03
sed -i -e "s/, wlc_d//g" src/wlcsim/initcond.f03
sed -i -e "s/,wlc_d//g" src/wlcsim/initcond.f03
sed -i -e "s/wlc_d,//g" src/wlcsim/initcond.f03
sed -i -e "s/wlc_d ,//g" src/wlcsim/initcond.f03
sed -i -e "s/wlc_d//g" src/wlcsim/initcond.f03
sed -i -e "s/wlc__/wlc_/g" src/wlcsim/initcond.f03
sed -i -e "s/wlcsim_data,//g" src/wlcsim/initcond.f03
sed -i -e "s/wlcsim_data ,//g" src/wlcsim/initcond.f03
sed -i -e "s/,wlcsim_data//g" src/wlcsim/initcond.f03
sed -i -e "s/, wlcsim_data//g" src/wlcsim/initcond.f03
sed -i -e "s/use params, only: wlcsim_data//g" src/wlcsim/initcond.f03
python unstructure.py src/wlcsim/restart_mpi.f90
sed -i -e "s/wlc_d%/wlc__/g" src/wlcsim/restart_mpi.f90
sed -i -e "s/, wlc_d//g" src/wlcsim/restart_mpi.f90
sed -i -e "s/,wlc_d//g" src/wlcsim/restart_mpi.f90
sed -i -e "s/wlc_d,//g" src/wlcsim/restart_mpi.f90
sed -i -e "s/wlc_d ,//g" src/wlcsim/restart_mpi.f90
sed -i -e "s/wlc_d//g" src/wlcsim/restart_mpi.f90
sed -i -e "s/wlc__/wlc_/g" src/wlcsim/restart_mpi.f90
sed -i -e "s/wlcsim_data,//g" src/wlcsim/restart_mpi.f90
sed -i -e "s/wlcsim_data ,//g" src/wlcsim/restart_mpi.f90
sed -i -e "s/,wlcsim_data//g" src/wlcsim/restart_mpi.f90
sed -i -e "s/, wlcsim_data//g" src/wlcsim/restart_mpi.f90
sed -i -e "s/use params, only: wlcsim_data//g" src/wlcsim/restart_mpi.f90
python unstructure.py src/wlcsim/wlcsim_bruno_mc.f03
sed -i -e "s/wlc_d%/wlc__/g" src/wlcsim/wlcsim_bruno_mc.f03
sed -i -e "s/, wlc_d//g" src/wlcsim/wlcsim_bruno_mc.f03
sed -i -e "s/,wlc_d//g" src/wlcsim/wlcsim_bruno_mc.f03
sed -i -e "s/wlc_d,//g" src/wlcsim/wlcsim_bruno_mc.f03
sed -i -e "s/wlc_d ,//g" src/wlcsim/wlcsim_bruno_mc.f03
sed -i -e "s/wlc_d//g" src/wlcsim/wlcsim_bruno_mc.f03
sed -i -e "s/wlc__/wlc_/g" src/wlcsim/wlcsim_bruno_mc.f03
sed -i -e "s/wlcsim_data,//g" src/wlcsim/wlcsim_bruno_mc.f03
sed -i -e "s/wlcsim_data ,//g" src/wlcsim/wlcsim_bruno_mc.f03
sed -i -e "s/,wlcsim_data//g" src/wlcsim/wlcsim_bruno_mc.f03
sed -i -e "s/, wlcsim_data//g" src/wlcsim/wlcsim_bruno_mc.f03
sed -i -e "s/use params, only: wlcsim_data//g" src/wlcsim/wlcsim_bruno_mc.f03
python unstructure.py src/wlcsim/wlcsim_bruno_looping_events.f03
sed -i -e "s/wlc_d%/wlc__/g" src/wlcsim/wlcsim_bruno_looping_events.f03
sed -i -e "s/, wlc_d//g" src/wlcsim/wlcsim_bruno_looping_events.f03
sed -i -e "s/,wlc_d//g" src/wlcsim/wlcsim_bruno_looping_events.f03
sed -i -e "s/wlc_d,//g" src/wlcsim/wlcsim_bruno_looping_events.f03
sed -i -e "s/wlc_d ,//g" src/wlcsim/wlcsim_bruno_looping_events.f03
sed -i -e "s/wlc_d//g" src/wlcsim/wlcsim_bruno_looping_events.f03
sed -i -e "s/wlc__/wlc_/g" src/wlcsim/wlcsim_bruno_looping_events.f03
sed -i -e "s/wlcsim_data,//g" src/wlcsim/wlcsim_bruno_looping_events.f03
sed -i -e "s/wlcsim_data ,//g" src/wlcsim/wlcsim_bruno_looping_events.f03
sed -i -e "s/,wlcsim_data//g" src/wlcsim/wlcsim_bruno_looping_events.f03
sed -i -e "s/, wlcsim_data//g" src/wlcsim/wlcsim_bruno_looping_events.f03
sed -i -e "s/use params, only: wlcsim_data//g" src/wlcsim/wlcsim_bruno_looping_events.f03
python unstructure.py src/wlcsim/init_MPI.f03
sed -i -e "s/wlc_d%/wlc__/g" src/wlcsim/init_MPI.f03
sed -i -e "s/, wlc_d//g" src/wlcsim/init_MPI.f03
sed -i -e "s/,wlc_d//g" src/wlcsim/init_MPI.f03
sed -i -e "s/wlc_d,//g" src/wlcsim/init_MPI.f03
sed -i -e "s/wlc_d ,//g" src/wlcsim/init_MPI.f03
sed -i -e "s/wlc_d//g" src/wlcsim/init_MPI.f03
sed -i -e "s/wlc__/wlc_/g" src/wlcsim/init_MPI.f03
sed -i -e "s/wlcsim_data,//g" src/wlcsim/init_MPI.f03
sed -i -e "s/wlcsim_data ,//g" src/wlcsim/init_MPI.f03
sed -i -e "s/,wlcsim_data//g" src/wlcsim/init_MPI.f03
sed -i -e "s/, wlcsim_data//g" src/wlcsim/init_MPI.f03
sed -i -e "s/use params, only: wlcsim_data//g" src/wlcsim/init_MPI.f03
python unstructure.py src/wlcsim/wlcsim_quinn.f03
sed -i -e "s/wlc_d%/wlc__/g" src/wlcsim/wlcsim_quinn.f03
sed -i -e "s/, wlc_d//g" src/wlcsim/wlcsim_quinn.f03
sed -i -e "s/,wlc_d//g" src/wlcsim/wlcsim_quinn.f03
sed -i -e "s/wlc_d,//g" src/wlcsim/wlcsim_quinn.f03
sed -i -e "s/wlc_d ,//g" src/wlcsim/wlcsim_quinn.f03
sed -i -e "s/wlc_d//g" src/wlcsim/wlcsim_quinn.f03
sed -i -e "s/wlc__/wlc_/g" src/wlcsim/wlcsim_quinn.f03
sed -i -e "s/wlcsim_data,//g" src/wlcsim/wlcsim_quinn.f03
sed -i -e "s/wlcsim_data ,//g" src/wlcsim/wlcsim_quinn.f03
sed -i -e "s/,wlcsim_data//g" src/wlcsim/wlcsim_quinn.f03
sed -i -e "s/, wlcsim_data//g" src/wlcsim/wlcsim_quinn.f03
sed -i -e "s/use params, only: wlcsim_data//g" src/wlcsim/wlcsim_quinn.f03
python unstructure.py src/wlcsim/wlcsim_brad.f03
sed -i -e "s/wlc_d%/wlc__/g" src/wlcsim/wlcsim_brad.f03
sed -i -e "s/, wlc_d//g" src/wlcsim/wlcsim_brad.f03
sed -i -e "s/,wlc_d//g" src/wlcsim/wlcsim_brad.f03
sed -i -e "s/wlc_d,//g" src/wlcsim/wlcsim_brad.f03
sed -i -e "s/wlc_d ,//g" src/wlcsim/wlcsim_brad.f03
sed -i -e "s/wlc_d//g" src/wlcsim/wlcsim_brad.f03
sed -i -e "s/wlc__/wlc_/g" src/wlcsim/wlcsim_brad.f03
sed -i -e "s/wlcsim_data,//g" src/wlcsim/wlcsim_brad.f03
sed -i -e "s/wlcsim_data ,//g" src/wlcsim/wlcsim_brad.f03
sed -i -e "s/,wlcsim_data//g" src/wlcsim/wlcsim_brad.f03
sed -i -e "s/, wlcsim_data//g" src/wlcsim/wlcsim_brad.f03
sed -i -e "s/use params, only: wlcsim_data//g" src/wlcsim/wlcsim_brad.f03
python unstructure.py src/wlcsim/wlcsim_bruno.f03
sed -i -e "s/wlc_d%/wlc__/g" src/wlcsim/wlcsim_bruno.f03
sed -i -e "s/, wlc_d//g" src/wlcsim/wlcsim_bruno.f03
sed -i -e "s/,wlc_d//g" src/wlcsim/wlcsim_bruno.f03
sed -i -e "s/wlc_d,//g" src/wlcsim/wlcsim_bruno.f03
sed -i -e "s/wlc_d ,//g" src/wlcsim/wlcsim_bruno.f03
sed -i -e "s/wlc_d//g" src/wlcsim/wlcsim_bruno.f03
sed -i -e "s/wlc__/wlc_/g" src/wlcsim/wlcsim_bruno.f03
sed -i -e "s/wlcsim_data,//g" src/wlcsim/wlcsim_bruno.f03
sed -i -e "s/wlcsim_data ,//g" src/wlcsim/wlcsim_bruno.f03
sed -i -e "s/,wlcsim_data//g" src/wlcsim/wlcsim_bruno.f03
sed -i -e "s/, wlcsim_data//g" src/wlcsim/wlcsim_bruno.f03
sed -i -e "s/use params, only: wlcsim_data//g" src/wlcsim/wlcsim_bruno.f03
python unstructure.py src/wlcsim/MCparrll_mpi.f90
sed -i -e "s/wlc_d%/wlc__/g" src/wlcsim/MCparrll_mpi.f90
sed -i -e "s/, wlc_d//g" src/wlcsim/MCparrll_mpi.f90
sed -i -e "s/,wlc_d//g" src/wlcsim/MCparrll_mpi.f90
sed -i -e "s/wlc_d,//g" src/wlcsim/MCparrll_mpi.f90
sed -i -e "s/wlc_d ,//g" src/wlcsim/MCparrll_mpi.f90
sed -i -e "s/wlc_d//g" src/wlcsim/MCparrll_mpi.f90
sed -i -e "s/wlc__/wlc_/g" src/wlcsim/MCparrll_mpi.f90
sed -i -e "s/wlcsim_data,//g" src/wlcsim/MCparrll_mpi.f90
sed -i -e "s/wlcsim_data ,//g" src/wlcsim/MCparrll_mpi.f90
sed -i -e "s/,wlcsim_data//g" src/wlcsim/MCparrll_mpi.f90
sed -i -e "s/, wlcsim_data//g" src/wlcsim/MCparrll_mpi.f90
sed -i -e "s/use params, only: wlcsim_data//g" src/wlcsim/MCparrll_mpi.f90
