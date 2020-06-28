!subroutine for initializing MPI
!Brad Krajina
!Last edited 2017/1/7
#if MPI_VERSION
subroutine init_mpi()
! values from wlcsim_data
   use params, only: wlc_repSuffix, wlc_error, wlc_numProcesses, wlc_id

   use mpi
   use params, only: wlcsim_params, MAXFILENAMELEN

   implicit none

   ! Simulation state

   ! string to work with
   character(MAXFILENAMELEN) iostrg    ! for file naming

   !Initialize MPI
   call MPI_Init(wlc_error)
   call stop_if_err(wlc_error, "Failed to MPI_Init.")

   ! Get number of processes
   call MPI_Comm_size(MPI_COMM_WORLD, wlc_numProcesses, wlc_error)
   call stop_if_err(wlc_error, "Failed to get num_processes.")

   ! Get individual process id
   call MPI_Comm_rank(MPI_COMM_WORLD, wlc_id, wlc_error)
   call stop_if_err(wlc_error, "Failed to get num_processes.")

   ! give each replica a different suffix to append to files
   write (iostrg, "(I4)") wlc_id
   iostrg = adjustL(iostrg)
   iostrg = trim(iostrg)
   iostrg = "v"//trim(iostrg)
   iostrg = trim(iostrg)
   wlc_repSuffix = iostrg
end subroutine
#endif
