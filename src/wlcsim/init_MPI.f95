!subroutine for initializing MPI
!Brad Krajina
!Last edited 2017/1/7
#if MPI_VERSION
subroutine init_mpi(wlc_d)

  use mpi
  use params, only: wlcsim_params, wlcsim_data
  implicit none

  ! Simulation state
  type(wlcsim_data), intent(inout) :: wlc_d

  !Initialize MPI
  call MPI_Init(wlc_d%error)
  call stop_if_err(wlc_d%error, "Failed to MPI_Init.")

  ! Get number of processes
  call MPI_Comm_size ( MPI_COMM_WORLD, wlc_d%numProcesses, wlc_d%error )
  call stop_if_err(wlc_d%error, "Failed to get num_processes.")

  ! Get individual process id
  call MPI_Comm_rank ( MPI_COMM_WORLD, wlc_d%id, wlc_d%error )
  call stop_if_err(wlc_d%error, "Failed to get num_processes.")
end subroutine
#endif
