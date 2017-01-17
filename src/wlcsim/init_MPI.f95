!subroutine for initializing MPI
!Brad Krajina
!Last edited 2017/1/7
#if MPI_VERSION
subroutine init_mpi(wlc_d)
  
  use mpi
  use params, only: wlcsim_params, wlcsim_data
  implicit none

  ! Simulation state
  type(wlcsim_data)   :: wlc_d

  !Initialize MPI
  call MPI_Init ( wlc_d%error )
  if (wlc_d%error.ne.0) then
     print*, "MPI_Init", wlc_d%error
  endif

  ! Get number of processes
  call MPI_Comm_size ( MPI_COMM_WORLD, wlc_d%numProcesses, wlc_d%error )
  if (wlc_d%error.ne.0) then
     print*, "MPI_Comm_size", wlc_d%error
  endif

  ! Get individual process id
  call MPI_Comm_rank ( MPI_COMM_WORLD, wlc_d%id, wlc_d%error )
  if (wlc_d%error.ne.0) then
     print*, "MPI_Comm_rank", wlc_d%error
  endif
end subroutine 
#endif
