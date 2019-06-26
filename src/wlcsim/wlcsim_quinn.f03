#include "../defines.inc"

subroutine wlcsim_quinn(save_ind, wlc_p)
! values from wlcsim_data
use params, only: wlc_mc_ind, wlc_numProcesses, wlc_id
#if MPI_VERSION
    use mpi
#endif
    use params
    implicit none
    integer, intent(in) :: save_ind ! 1, 2, ...
    type(wlcsim_params), intent(inout) :: wlc_p

    ! to minimize code rewriting, we use our old name for save_ind internally
    wlc_mc_ind = save_ind



#if MPI_VERSION
    if (wlc_numProcesses == 1) then
        call onlyNode(wlc_p)
    elseif (wlc_id == 0) then
        call head_node(wlc_numProcesses)
    else
        call worker_node(wlc_p)
    endif
#else
    call onlyNode(wlc_p)
#endif



    print*, '________________________________________'
    print*, 'Time point ',save_ind, ' out of', WLC_P__NUMSAVEPOINTS, 'Thread id', wlc_id
    call printEnergies()
    call printWindowStats(wlc_p)
    !call wlcsim_params_printPhi(wlc_p)

end subroutine wlcsim_quinn

#if MPI_VERSION
subroutine head_node(process)
! values from wlcsim_data
use params, only: wlc_mc_ind, wlc_rand_stat
    use mersenne_twister
    use mpi
    use params
    use energies
!   MPI variables
    implicit none
    integer, intent(in) :: process ! number of therads
    integer ( kind = 4 ) dest   !destination id for messages
    integer ( kind = 4 ) source  !source id for messages
    integer ( kind = 4 ) error  ! error id for MIP functions
    integer ( kind = 4 ) status(MPI_status_SIZE) ! MPI stuff

    !   variable for random number generator seeding
    real(dp) urand(1)

!   for head node use only variables
    integer rep ! physical replica number, for loops
    integer temp ! for castling
    logical keepGoing   ! set to false when NaN encountered
    real(dp) x(NUMBER_OF_ENERGY_TYPES) ! slice of xMtrx
    real(dp) cof(NUMBER_OF_ENERGY_TYPES) ! slice of cofMtrx
    integer N_average      ! number of attempts since last average
    integer upSuccess(process-1)  ! number of successes since last average
    integer downSuccess(process-1) ! number of successes since last average
    integer nExchange ! total number of exchanges attemted
    real(dp) energy ! for deciding to accept exchange
    integer term ! for loopin over terms
    real(dp) cof_path_by_energy_type ! function
    integer nPTReplicas

    !   Quinn's parallel tempering head node variables
    integer, allocatable :: nodeNumber(:)  ! list of which nodes are which
    real(dp), allocatable :: xMtrx(:,:)  ! sum of bound states
    real(dp), allocatable :: cofMtrx(:,:) ! mu or chi or whatever
    real(dp), allocatable :: s_vals(:) ! path parameter
    integer ii

    nPTReplicas=process-1
    nExchange = 0

    ! Allocate the head node variables for keeping track of which node is which
    allocate( xMtrx(nPTReplicas,NUMBER_OF_ENERGY_TYPES))
    allocate( cofMtrx(nPTReplicas,NUMBER_OF_ENERGY_TYPES))
    allocate( nodeNumber(nPTReplicas))
    allocate( s_vals(nPTReplicas))

    ! Keep track of up and down successes (head node only needs this)
    do rep = 1,nPTReplicas
        upSuccess(rep) = 0
        downSuccess(rep) = 0
        s_vals(rep) = WLC_P__INITIAL_MAX_S*(dble(rep)-1.0_dp)/(dble(nPTReplicas)-1.0_dp)
    enddo

    ! Set initial values for parrallel tempering values
    do rep = 1,nPTReplicas
        do ii = 1,NUMBER_OF_ENERGY_TYPES
            cofMtrx(rep,ii) = cof_path_by_energy_type(ii,s_vals(rep))
        enddo
    enddo

    N_average = 0

    ! Initially replica numbers are same as nodes
    do rep = 1,nPTReplicas
        nodeNumber(rep) = rep
    enddo
    !-----------------------------------------------------
    !
    !    Begin Main loop
    !
    !------------------------------------------------------
    keepGoing = .True.
    do while(keepGoing)
        ! give workers thier jobs
        do rep = 1,nPTReplicas
            dest = nodeNumber(rep)
            call MPI_Send (rep,1, MPI_integer, dest,   0, &
                            MPI_COMM_WORLD,error )
            if (WLC_P__RESTART.and.wlc_mc_ind.eq.1) then
                source = dest
                call MPI_Recv (cof, NUMBER_OF_ENERGY_TYPES, MPI_doUBLE_PRECISION, source, 0, &
                               MPI_COMM_WORLD, status, error )
                cofMtrx(rep,:) = cof
            else
                cof = cofMtrx(rep,:)
                call MPI_Send (cof,NUMBER_OF_ENERGY_TYPES, MPI_doUBLE_PRECISION, dest,   0, &
                                MPI_COMM_WORLD,error )
            endif
        enddo
        ! get results from workers

        do rep = 1,nPTReplicas
            source = nodeNumber(rep)
            call MPI_Recv ( x, NUMBER_OF_ENERGY_TYPES, MPI_doUBLE_PRECISION, source, 0, &
                           MPI_COMM_WORLD, status, error )
            xMtrx(rep,:) = x
            if(isnan(x(1))) then ! endo of program
                keepGoing = .false.
                return
            endif
        enddo

        source = 1
        call MPI_Recv (wlc_mc_ind, 1, MPI_integer, source, 0, &
                       MPI_COMM_WORLD, status, error )

        ! do replica exchange
        do rep = 1,(nPTReplicas-1)
            energy = 0.0_dp
            do term = 1,NUMBER_OF_ENERGY_TYPES
                energy = energy-(xMtrx(rep + 1,term)-xMtrx(rep,term))*&
                              (cofMtrx(rep + 1,term)-cofMtrx(rep,term))
            enddo
            call random_number(urand,wlc_rand_stat)
            if (exp(-1.0_dp*energy).gt.urand(1)) then
                if (WLC_P__PTON) then
                    temp = nodeNumber(rep)
                    nodeNumber(rep) = nodeNumber(rep + 1)
                    nodeNumber(rep + 1) = temp
                endif
                upSuccess(rep) = upSuccess(rep) + 1
                downSuccess(rep + 1) = downSuccess(rep + 1) + 1
                x = xMtrx(rep,:)
                xMtrx(rep,:) = xMtrx(rep + 1,:)
                xMtrx(rep + 1,:) = x
                nExchange = nExchange + 1
            endif
        enddo


        ! track/adapt acceptance rates
        N_average = N_average + 1
        if (N_average.ge.WLC_P__NREPADAPT) then
            call save_repHistory(upSuccess,downSuccess,nPTReplicas, &
                                 cofMtrx,xMtrx,nodeNumber,N_average,&
                                 nExchange,wlc_mc_ind,NUMBER_OF_ENERGY_TYPES,s_vals)

            if ((wlc_mc_ind.ge.WLC_P__INDSTARTREPADAPT).and. &
                (wlc_mc_ind.lt.WLC_P__INDENDREPADAPT)) then ! insert input defined location here
                call adaptCof(downSuccess,nPTReplicas,s_vals,N_average,&
                               WLC_P__LOWERREPEXE,WLC_P__UPPERREPEXE,&
                               WLC_P__LOWERCOFRAIL,WLC_P__UPPERCOFRAIL,&
                               WLC_P__REPANNEALSPEED,WLC_P__REPLICABOUNDS)
                do rep = 1,nPTReplicas
                    do ii = 1,NUMBER_OF_ENERGY_TYPES
                        cofMtrx(rep,ii) = cof_path_by_energy_type(ii,s_vals(rep))
                    enddo
                enddo
            endif
            N_average = 0
            do rep = 1,nPTReplicas
                upSuccess(rep) = 0
                downSuccess(rep) = 0
            enddo
        endif
    enddo

    deallocate(xMtrx)
    deallocate(cofMtrx)
    deallocate(nodeNumber)
    deallocate(s_vals)

end subroutine head_node
#endif

function cof_path_by_energy_type(energy_type, s) result(cof)
    use energies, only: energyOf, mu_, umbrella_, umbrellaQuadratic_
    use umbrella, only: setUmbrellaCof, setUmbrellaQuadraticCof
    use params, only: dp, nan
    implicit none
    integer, intent(in) :: energy_type
    real(dp), intent(in) :: s
    real(dp) cof

    cof = energyOf(energy_type)%cof ! Use Default value

    if (energyOf(energy_type)%parallel_temper) then
        ! special instructions for specified types
        if (energy_type == umbrella_ .and. WLC_P__UMBRELLA) then
            cof = setUmbrellaCof(s)
        elseif (energy_type == umbrellaQuadratic_ .and. WLC_P__UMBRELLA) then
            cof = setUmbrellaQuadraticCof(s)
        elseif (energy_type == mu_ ) then
            cof = s-2.5_dp ! set for Quinn's chromatin problem
        else
            ! Parallel temper from 0 to default value
            cof = energyOf(energy_type)%cof*s/WLC_P__INITIAL_MAX_S
        endif
    endif
    return
end function cof_path_by_energy_type


#if MPI_VERSION
subroutine worker_node(wlc_p)
! values from wlcsim_data
use params, only: wlc_ind_exchange, wlc_mc_ind
    use energies, only: energyOf, NUMBER_OF_ENERGY_TYPES
    use mpi
    use params
    use mersenne_twister
    implicit none
    integer ( kind = 4 ), save :: id = -1     ! which processor I am
    integer ( kind = 4 ) error  ! error id for MIP functions
    type(wlcsim_params), intent(inout) :: wlc_p
    integer i, ii
    logical system_has_been_changed
    real :: start, finish

    if (id == -1) then
        call MPI_Comm_rank(MPI_COMM_WORLD, id, error)
        call stop_if_err(error, "Failed to get num_processes.")
    endif
    if (wlc_mc_ind == 1) then
        if (WLC_P__PT_CHI .or. WLC_P__PT_H .or. WLC_P__PT_KAP .or. WLC_P__PT_MU .or. WLC_P__PT_COUPLE .or. WLC_P__ENSEMBLE_BIND) then
            call startWorker()
        endif
    endif

    call schedule(wlc_p,system_has_been_changed)

    if (system_has_been_changed) then
        call CalculateEnergiesFromScratch(wlc_p) ! Calculate dE and dx from scratch
        do ii = 1, NUMBER_OF_ENERGY_TYPES
            energyOf(ii)%E = energyOf(ii)%dE
            energyOf(ii)%x = energyOf(ii)%dx
        enddo
    else
        call VerifyEnergiesFromScratch(wlc_p)
    endif

    ! ------------------------------
    !
    ! call main simulation code
    !
    !  --------------------------------

    call cpu_time(start)
    do i = 1,WLC_P__NREPLICAEXCHANGEPERSAVEPOINT
        wlc_ind_exchange=i
        !   * Perform a MC simulation *
        call MCsim(wlc_p,WLC_P__STEPSPEREXCHANGE)

        !   * Replica Exchange *
        call replicaExchange()
    enddo
    call cpu_time(finish)
    print*, "Save Point time", finish-start, " seconds"
end subroutine worker_node
#endif

subroutine onlyNode(wlc_p)
! values from wlcsim_data
    use params, only: wlcsim_params
    use energies, only: energyOf, NUMBER_OF_ENERGY_TYPES
    implicit none
    type(wlcsim_params), intent(inout) :: wlc_p
    logical system_has_been_changed
    real :: start, finish
    integer ii
    !   * Perform a MC simulation *
    call schedule(wlc_p,system_has_been_changed)
    if (system_has_been_changed) then
        call CalculateEnergiesFromScratch(wlc_p)
        do ii = 1, NUMBER_OF_ENERGY_TYPES
            energyOf(ii)%E = energyOf(ii)%dE
            energyOf(ii)%x = energyOf(ii)%dx
        enddo
    else
        call VerifyEnergiesFromScratch(wlc_p)
    endif
    call cpu_time(start)
    call MCsim(wlc_p,WLC_P__NREPLICAEXCHANGEPERSAVEPOINT*WLC_P__STEPSPEREXCHANGE)
    call cpu_time(finish)
    print*, "Save Point time", finish-start, " seconds"
end subroutine onlyNode
subroutine schedule(wlc_p,system_has_been_changed)
! values from wlcsim_data
use params, only: wlc_mc_ind, eps
    use energies, only: energyOf, NUMBER_OF_ENERGY_TYPES
    use params
    implicit none
    type(wlcsim_params), intent(inout) :: wlc_p
    logical, intent(out) :: system_has_been_changed
    integer ii

    system_has_been_changed = .False.
    ! ------------------------------
    !
    ! Different instructions for each save point
    !
    !  --------------------------------
    if (wlc_mc_ind <= 1) system_has_been_changed = .TRUE.
    if (wlc_mc_ind < WLC_P__NNOINT) then
        wlc_p%field_int_on_currently = .false.
    elseif (WLC_P__FIELD_INT_ON) then
        if (.not.wlc_p%field_int_on_currently)  system_has_been_changed = .TRUE.
        wlc_p%field_int_on_currently = .true.
    endif

    do ii = 1,NUMBER_OF_ENERGY_TYPES
        if(wlc_mc_ind >= energyOf(ii)%ind_on) then
            if (energyOf(ii)%isOn) cycle
            system_has_been_changed = .TRUE.
            energyOf(ii)%isOn =  .TRUE.
        endif
    enddo

    if ((wlc_mc_ind.gt.WLC_P__INDSTARTREPADAPT).and. &
        (wlc_mc_ind.le.WLC_P__INDENDREPADAPT)) then ! addapt Cof was run
        system_has_been_changed = .TRUE.
    endif
end subroutine
