! use mpi

subroutine wlcsim_quinn(save_ind, mc, md)
    use params
    implicit none
    integer, intent(in) :: save_ind ! 1, 2, ...
    type(wlcsim_params), intent(inout) :: mc
    type(wlcsim_data), intent(inout) :: md
    integer, save :: i, id, num_processes
    integer (kind=4) error

    if (save_ind == 1) then
        call MPI_Init(error)
        call stop_if_err(error, "Failed to MPI_Init.")
        call MPI_Comm_size(MPI_COMM_WORLD, num_processes, error)
        call stop_if_err(error, "Failed to get num_processes.")
        call MPI_Comm_rank(MPI_COMM_WORLD, id, error)
        call stop_if_err(error, "Failed to get num_processes.")
        if (num_processes > 1) then
            call setup_mpi_stuff()
        endif
    endif

     ! for changing constants durring run
     if (mc%useSchedule) then
         call strength_schedule(mc,mc%field_int_on)
     endif


    do i=1,mc%numReplicaExchangesBetweenSaves
!   * Perform a MC simulation *
        call VerifyEnegiesFromScratch(mc, md)
        call MCsim(mc,md,mc%numStepsBetweenExchanges,mc%field_int_on,md%rand_stat)
        print*, '________________________________________'
        print*, 'Time point ',save_ind, ' out of', save_indMAX
        call wlcsim_params_printEnergies(mc)
        call wlcsim_params_printWindowStats(mc)
        !call wlcsim_params_printPhi(mc,md)
    enddo
end subroutine wlcsim_quinn

subroutine setup_mpi()
  ! use mpi
  use params, only: dp

  implicit none
  integer ( kind = 4 ) error
  integer ( kind = 4 ) id
  integer ( kind = 4 ) p
!
!  Initialize MPI.
!
  call MPI_Init ( error )
  if (error.ne.0) then
      print*, "MPI_Init", error
  endif
!
!  Get the number of processes.
!
  call MPI_Comm_size ( MPI_COMM_WORLD, p, error )
  if (error.ne.0) then
      print*, "MPI_Comm_size", error
  endif
!
!  Get the individual process ID.
!
  call MPI_Comm_rank ( MPI_COMM_WORLD, id, error )
  if (error.ne.0) then
      print*, "MPI_Comm_rank", error
  endif
!
!  print a message.
!
  if ( id == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  WLC MC sim with tempering using MPI'
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  The number of threads being used is ', p
    write ( *, '(a,i8)' ) '  The number of replicas is ', p-1


  end if
  if (p.gt.1) then
      call paraTemp ( p, id )
  else
      call singleCall()
  endif
end subroutine


subroutine wlcsim_quinn_orig(rand_stat)

!
!     This simulation tracks the dynamics of a single polymer
!     chain modeled as a discrete wormlike chain with bending
!     and stretching energy.
!
!     Andrew Spakowitz
!     Written 9-2-13
!
!     Edited by Shifan  prior to 2016
!     Edited heavily by Quinn in spring of 2016
!
!     Variables within the simulation

  use mersenne_twister  ! so that we know the size of rand_stat
  use params

  IMPLICIT NONE
  type(random_stat), intent(inout) :: rand_stat ! state of random number generator

  ! miscellaneous
  integer I
  character(4) fileind       ! index of output
  character(MAXFILENAMELEN) iostr       ! file for output

!     Simulation input variables

  integer INTON             ! Include polymer interactions

! simulation data strucutres
  type(wlcsim_params) mc
  type(wlcsim_data) md

!-------------------------------------------------------
!
!    Set simulation parameers
!
!--------------------------------------------------------
  iostr='input/params'
  print*, "setting parameters from: ", iostr
  call wlcsim_params_setparams(mc,iostr)
  call wlcsim_params_allocate(mc,md)

!-----------------------------------------------------
!
!   Initial condition / Restart
!
!----------------------------------------------------
  !  Calculate volume of bins
  if (mc%confinetype.eq.3) then
      call MC_caclVolume(mc%confinetype,mc%NBINX,mc%dbin, mc%lbox(1), &
                         md%Vol,rand_stat)  ! calculate partial volumes
  else
      do I=1,mc%NBIN
           md%Vol(I)=mc%dbin**3
      enddo
  endif

  if (mc%restart) then
      call pt_restart(mc,md)
      !print*, '-----load simulation-----'
      !iostr='BinaryfileName'
      !stop 1
      !call wlcsim_params_readBindary(mc,md,iostr)
  else


    !   Setup the initial condition
      call initcond(md%R,md%U,md%AB,mc%NT,mc%NB,mc%NP,mc%FRMfile,pack_as_para(mc),mc%lbox, &
                    mc%initCondType,rand_stat)

    !   Load in AB sequence
      if (mc%FRMCHEM) then
          iostr='input/ab'
          call wlcsim_params_loadAB(mc,md,iostr)
      else
          call initchem(md%AB,mc%NT,mc%NB,mc%nBpM,mc%NP,mc%FA,mc%LAM,rand_stat)
      endif


    !   Load methalation sequence
      if (mc%FRMchem) then
          ! more to come here ...
          print*, "wlcsim: FRMMETH not fininshed"
          stop 1
      else
          call initchem(md%METH,mc%NT,mc%NB,mc%nBpM,mc%NP,mc%F_METH,mc%LAM_METH,md%rand_stat)
      endif

    ! Load External field
      if (mc%FRMField) then
          iostr='input/h_A'
          Call wlcsim_params_LoadField(mc,md,iostr)
      else
          Call wlcsim_params_MakeField(mc,md)
      endif

      ! Get assignement from head node
      call PT_override(mc,md)

      iostr='data/r0'
      I=0;
      call wlcsim_params_saveR(mc,md,iostr,.false.)

      iostr='data/params'
      call wlcsim_params_saveparameters(mc,iostr)

      iostr='data/u0'
      call wlcsim_params_saveU(mc,md,iostr)

  endif

 ! call wlcsim_params_printDescription(mc)  ! output simulation parameters

!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!              Begin simulation
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  print*, 'Beginning simulation: rep', md%rep, " id", mc%id

  do WHILE ((save_ind).LE.save_indMAX)

     ! for changing constants durring run
     if (mc%useSchedule) then
         call strength_schedule(mc,INTON)
     endif


!   * Perform a MC simulation *
    call MCsim(mc,md,mc%NSTEP,INTON,rand_stat)

!    Save the conformation and the metrics
    write (fileind,'(I4)'), save_ind

    !Save various energy contiributions to file
    iostr='data/out1'
    call wlcsim_params_appendEnergyData(mc,iostr)

    !part 2.5 - adaptations
    iostr='data/out3'
    call wlcsim_params_appendAdaptData(mc,iostr)

    if (mc%savePhi) then
        write(iostr,"(I6)"), save_ind
        iostr='data/phi' // trim(adjustL(iostr))
        call wlcsim_params_savePHI(mc,md,iostr)
    endif

    write(iostr,"(I6)"), save_ind
    iostr='data/r' // trim(adjustL(iostr))
    call wlcsim_params_saveR(mc,md,iostr,0)

    if (mc%saveU) then
        write(iostr,"(I6)"), save_ind
        iostr='data/u' // trim(adjustL(iostr))
        call wlcsim_params_saveU(mc,md,iostr)
    endif


    print*, '________________________________________'
    print*, 'Time point ',save_ind, ' out of', save_indMAX
    call wlcsim_params_printEnergies(mc)
    call wlcsim_params_printWindowStats(mc)
    !call wlcsim_params_printPhi(mc,md)
    save_ind=save_ind+1
  enddo

end
