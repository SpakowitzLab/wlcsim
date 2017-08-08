
subroutine wlcsim_quinn(save_ind, wlc_d, wlc_p)
#if MPI_VERSION
    use mpi
#endif
    use params
    implicit none
    integer, intent(in) :: save_ind ! 1, 2, ...
    type(wlcsim_params), intent(inout) :: wlc_p
    type(wlcsim_data), intent(inout) :: wlc_d
    integer (kind = 4) error

    ! to minimize code rewriting, we use our old name for save_ind internally
    wlc_d%mc_ind = save_ind



#if MPI_VERSION
    if (wlc_d%numProcesses == 1) then
        call onlyNode(wlc_p, wlc_d)
    elseif (wlc_d%id == 0) then
        call head_node(wlc_p, wlc_d,wlc_d%numProcesses)
    else
        call worker_node(wlc_p, wlc_d)
    endif
#else
    call onlyNode(wlc_p, wlc_d)
#endif



    print*, '________________________________________'
    print*, 'Time point ',save_ind, ' out of', wlc_p%numSavePoints, 'Thread id', wlc_d%id
    call printEnergies(wlc_d)
    call printWindowStats(wlc_p, wlc_d)
    !call wlcsim_params_printPhi(wlc_p, wlc_d)

end subroutine wlcsim_quinn

#if MPI_VERSION
subroutine head_node(wlc_p, wlc_d,process)
    use mersenne_twister
    use mpi
    use params
!   MPI variables
    implicit none
    integer, intent(in) :: process ! number of therads
    integer ( kind = 4 ) dest   !destination id for messages
    integer ( kind = 4 ) source  !source id for messages
    integer ( kind = 4 ) id     ! which processor I am
    integer ( kind = 4 ) error  ! error id for MIP functions
    integer ( kind = 4 ) status(MPI_status_SIZE) ! MPI stuff
    type(wlcsim_params), intent(inout) :: wlc_p
    type(wlcsim_data), intent(inout) :: wlc_d

    !   variable for random number generator seeding
    type(random_stat) rand_stat  ! state of random number chain
    integer Irand     ! Seed
    character(8) datedum  ! trash
    character(10) timedum ! trash
    character(5) zonedum  ! trash
    integer seedvalues(8) ! clock readings
    real urand(1)

!   for head node use only variables
    integer rep ! physical replica number, for loops
    integer temp ! for castling
    logical keepGoing   ! set to false when NaN encountered
    integer, parameter :: nTerms = 9  ! number of energy terms
    real(dp) x(nTerms) ! slice of xMtrx
    real(dp) cof(nTerms) ! slice of cofMtrx
    integer N_average      ! number of attempts since last average
    integer upSuccess(process-1)  ! number of successes since last average
    integer downSuccess(process-1) ! number of successes since last average
    integer nExchange ! total number of exchanges attemted
    real(dp) energy ! for deciding to accept exchange
    integer term ! for loopin over terms
    real(dp) h_path,chi_path,mu_path,kap_path,HP1_Bind_path,maierSaupe_path ! functions
    integer nPTReplicas

    !   Quinn's parallel tempering head node variables
    integer, allocatable :: nodeNumber(:)  ! list of which nodes are which
    real(dp), allocatable :: xMtrx(:,:)  ! sum of bound states
    real(dp), allocatable :: cofMtrx(:,:) ! mu or chi or whatever
    real(dp), allocatable :: s_vals(:) ! path parameter

    nPTReplicas=process-1

    ! Allocate the head node variables for keeping track of which node is which
    allocate( xMtrx(nPTReplicas,nTerms))
    allocate( cofMtrx(nPTReplicas,nTerms))
    allocate( nodeNumber(nPTReplicas))
    allocate( s_vals(nPTReplicas))

    ! Keep track of up and down successes (head node only needs this)
    do rep = 1,nPTReplicas
        upSuccess(rep) = 0
        downSuccess(rep) = 0
        s_vals(rep) = wlc_p%inITIAL_MAX_S*(dble(rep)-1.0_dp)/(dble(nPTReplicas)-1.0_dp)
    enddo

    ! Set initial values for parrallel tempering values
    do rep = 1,nPTReplicas
        if (wlc_p%PT_chi) then
            cofMtrx(rep,1) = chi_path(s_vals(rep))
        else
            cofMtrx(rep,1) = wlc_p%chi
        endif
        if (wlc_p%PT_mu) then
            cofMtrx(rep,2) = mu_path(s_vals(rep))
        else
            cofMtrx(rep,2) = wlc_p%mu
        endif
        if (wlc_p%PT_h) then
            cofMtrx(rep,3) = h_path(s_vals(rep))
        else
            cofMtrx(rep,3) = wlc_p%hA
        endif
        if (wlc_p%PT_couple) then
            cofMtrx(rep,4) = HP1_Bind_path(s_vals(rep))
        else
            cofMtrx(rep,4) = wlc_p%HP1_Bind
        endif
        if (wlc_p%PT_Kap) then
            cofMtrx(rep,5) = kap_path(s_vals(rep))
        else
            cofMtrx(rep,5) = wlc_p%KAP
        endif
        cofMtrx(rep,6) = 0
        cofMtrx(rep,7) = 0
        cofMtrx(rep,8) = 0
        if (wlc_p%PT_MaierSaupe) then
            cofMtrx(rep,9) = maierSaupe_path(s_vals(rep))
        else
            cofMtrx(rep,9) = wlc_p%chi_l2
        endif
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
            if (wlc_p%restart.and.wlc_d%mc_ind.eq.1) then
                source = dest
                call MPI_Recv (cof, nTerms, MPI_doUBLE_PRECISION, source, 0, &
                               MPI_COMM_WORLD, status, error )
                cofMtrx(rep,:) = cof
            else
                cof = cofMtrx(rep,:)
                call MPI_Send (cof,nTerms, MPI_doUBLE_PRECISION, dest,   0, &
                                MPI_COMM_WORLD,error )
            endif
        enddo
        ! get results from workers

        do rep = 1,nPTReplicas
            source = nodeNumber(rep)
            call MPI_Recv ( x, nTerms, MPI_doUBLE_PRECISION, source, 0, &
                           MPI_COMM_WORLD, status, error )
            xMtrx(rep,:) = x
            if(isnan(x(1))) then ! endo of program
                keepGoing = .false.
                return
            endif
        enddo

        source = 1
        call MPI_Recv (wlc_d%mc_ind, 1, MPI_integer, source, 0, &
                       MPI_COMM_WORLD, status, error )

        ! do replica exchange
        do rep = 1,(nPTReplicas-1)
            energy = 0.0_dp
            do term = 1,nTerms
                energy = energy-(xMtrx(rep + 1,term)-xMtrx(rep,term))*&
                              (cofMtrx(rep + 1,term)-cofMtrx(rep,term))
            enddo
            call random_number(urand,rand_stat)
            if (exp(-1.0_dp*energy).gt.urand(1)) then
                if (wlc_p%PTON) then
                    temp = nodeNumber(rep)
                    nodeNumber(rep) = nodeNumber(rep + 1)
                    nodeNumber(rep + 1) = temp
                endif
                upSuccess(rep) = upSuccess(rep) + 1
                downSuccess(rep + 1) = downSuccess(rep + 1) + 1
                x = xMtrx(rep,:)
                xMtrx(rep,:) = xMtrx(rep + 1,:)
                xMtrx(rep + 1,:) = x
            endif
        enddo


        ! track/adapt acceptance rates
        N_average = N_average + 1
        if (N_average.ge.wlc_p%NRepAdapt) then
            call save_repHistory(upSuccess,downSuccess,nPTReplicas, &
                                 cofMtrx,xMtrx,nodeNumber,N_average,&
                                 nExchange,wlc_d%mc_ind,nTerms,s_vals)

            if ((wlc_d%mc_ind.ge.wlc_p%indStartRepAdapt).and. &
                (wlc_d%mc_ind.lt.wlc_p%indendRepAdapt)) then ! insert input defined location here
                call adaptCof(downSuccess,nPTReplicas,s_vals,N_average,&
                               wlc_p%lowerRepExe,wlc_p%upperRepExe,&
                               wlc_p%lowerCofRail,wlc_p%upperCofRail,&
                               wlc_p%RepAnnealSpeed,wlc_p%replicaBounds)
                do rep = 1,nPTReplicas
                    if (wlc_p%PT_chi) then
                        cofMtrx(rep,1) = chi_path(s_vals(rep))
                    endif
                    if (wlc_p%PT_mu) then
                        cofMtrx(rep,2) = mu_path(s_vals(rep))
                    endif
                    if (wlc_p%PT_h) then
                        cofMtrx(rep,3) = h_path(s_vals(rep))
                    endif
                    if (wlc_p%PT_couple) then
                        cofMtrx(rep,4) = HP1_Bind_path(s_vals(rep))
                    endif
                    if (wlc_p%PT_Kap) then
                        cofMtrx(rep,5) = kap_path(s_vals(rep))
                    endif
                enddo
            endif
            N_average = 0
            do rep = 1,nPTReplicas
                upSuccess(rep) = 0
                downSuccess(rep) = 0
            enddo
        endif
        nExchange = nExchange + 1
    enddo

    deallocate(xMtrx)
    deallocate(cofMtrx)
    deallocate(nodeNumber)
    deallocate(s_vals)

end subroutine head_node
#endif

function chi_path(s) result(chi)
    use params, only: dp
    implicit none
    real(dp), intent(in) :: s
    real(dp) chi
    real(dp) chi_max
    if (.false.) then
        chi_max = 2.00_dp
        if (s.lt.0.5_dp) then
            chi = 0.0_dp
        else
            chi = chi_max*2.0_dp*(s-0.5_dp)
        endif
    else
        chi = s
    endif
end function chi_path
function h_path(s) result(h)
    use params, only: dp
    implicit none
    real(dp), intent(in) :: s
    real(dp) h
    real(dp) h_max
    if (.false.) then
        h_max = 10.0_dp
        if (s.lt.0.5_dp) then
             h = h_max*s*2.0_dp
        else
             h = h_max*(1.0_dp-s)*2.0_dp
        endif
    else
        h = 0.0_dp
    endif
end function h_path
function mu_path(s) result(mu)
    use params, only: dp
    implicit none
    real(dp), intent(in) :: s
    real(dp) mu
    mu = s
end function mu_path
function kap_path(s) result(kap)
    use params, only: dp
    implicit none
    real(dp), intent(in) :: s
    real(dp) kap
    kap = s*10.0_dp
end function kap_path
function maierSaupe_path(s) result(output)
    use params, only: dp
    implicit none
    real(dp), intent(in) :: s
    real(dp) output
    output = s*1.0_dp
end function maierSaupe_path
function hp1_bind_path(s) result(hp1_bind)
    use params, only: dp
    implicit none
    real(dp), intent(in) :: s
    real(dp) hp1_bind
    hp1_bind = s
end function hp1_bind_path

#if MPI_VERSION
subroutine worker_node(wlc_p, wlc_d)
    use mpi
    use params
    use mersenne_twister
    implicit none
    integer ( kind = 4 ), save :: id = -1     ! which processor I am
    integer ( kind = 4 ) error  ! error id for MIP functions
    integer ( kind = 4 ) status(MPI_status_SIZE) ! MPI stuff
    type(wlcsim_params), intent(inout) :: wlc_p
    type(wlcsim_data), intent(inout) :: wlc_d
    type(random_stat) rand_stat  ! state of random number chain
    integer i
    logical system_has_been_changed
    system_has_been_changed = .False.
    if (id == -1) then
        call MPI_Comm_rank(MPI_COMM_WORLD, id, error)
        call stop_if_err(error, "Failed to get num_processes.")
    endif
    if (wlc_d%mc_ind == 1) then
        if (wlc_p%PT_chi .or. wlc_p%PT_h .or. wlc_p%PT_kap .or. wlc_p%PT_mu .or. wlc_p%PT_couple) then
            call startWorker(wlc_p, wlc_d)
        endif
    endif

    ! ------------------------------
    !
    ! Different instructions for each save point
    !
    !  --------------------------------
    if (wlc_d%mc_ind <= 1) system_has_been_changed = .TRUE.
    if (wlc_d%mc_ind <= wlc_p%NNOinT) then
        wlc_p%field_int_on = .false.
    else
        if (.not.wlc_p%field_int_on)  system_has_been_changed = .TRUE.
        wlc_p%field_int_on = .true.
    endif
    if(wlc_d%mc_ind.lt.wlc_p%N_KAP_ON) then
        wlc_p%KAP_ON = 0.0_dp
    else
        if (wlc_p%KAP_ON.eq.0.0_dp) system_has_been_changed = .TRUE.
        wlc_p%KAP_ON = 1.0_dp
    endif

    if(wlc_d%mc_ind.lt.wlc_p%N_CHI_ON) then
        wlc_p%CHI_ON = 0.0_dp
    else
        if (wlc_p%CHI_ON.eq.0.0_dp) system_has_been_changed = .TRUE.
        wlc_p%CHI_ON = 1.0_dp
    endif

    if ((wlc_d%mc_ind.gt.wlc_p%indStartRepAdapt).and. &
        (wlc_d%mc_ind.le.wlc_p%indendRepAdapt)) then ! addapt Cof was run
        system_has_been_changed = .TRUE.
    endif

    if (system_has_been_changed) then
        call CalculateEnergiesFromScratch(wlc_p, wlc_d)
        if (wlc_p%field_int_on) then
            wlc_d%ECouple =wlc_d%dECouple
            wlc_d%EKap    =wlc_d%dEKap
            wlc_d%ECHI    =wlc_d%dECHI
            wlc_d%EField  =wlc_d%dEField
            wlc_d%x_Field =wlc_d%dx_Field
            wlc_d%x_couple = wlc_d%dx_couple
            wlc_d%x_Kap   =wlc_d%dx_Kap
            wlc_d%x_Chi   =wlc_d%dx_Chi
        else
            wlc_d%ECouple =0.0_dp
            wlc_d%EKap    =0.0_dp
            wlc_d%ECHI    =0.0_dp
            wlc_d%EField  =0.0_dp
            wlc_d%x_Field =0.0_dp
            wlc_d%x_couple = 0.0_dp
            wlc_d%x_Kap   =0.0_dp
            wlc_d%x_Chi   =0.0_dp
        endif
        if (wlc_p%bind_On) then
            wlc_d%ebind   =wlc_d%debind
            wlc_d%x_mu    =wlc_d%dx_mu
        else
            wlc_d%ebind   =0.0_dp
            wlc_d%x_mu    =0.0_dp
        endif
    else
        call VerifyEnergiesFromScratch(wlc_p, wlc_d)
    endif

    ! ------------------------------
    !
    ! call main simulation code
    !
    !  --------------------------------

    do i = 1,wlc_p%nReplicaExchangePerSavePoint
        wlc_d%ind_exchange=i
        !   * Perform a MC simulation *
        call MCsim(wlc_p, wlc_d,wlc_p%stepsPerExchange)

        !   * Replica Exchange *
        call replicaExchange(wlc_p,wlc_d)

    enddo
end subroutine worker_node
#endif

subroutine onlyNode(wlc_p, wlc_d)
    use params
    implicit none
    type(wlcsim_params), intent(inout) :: wlc_p
    type(wlcsim_data), intent(inout) :: wlc_d
    !   * Perform a MC simulation *
    call VerifyEnergiesFromScratch(wlc_p, wlc_d)
    call MCsim(wlc_p, wlc_d,wlc_p%nReplicaExchangePerSavePoint*wlc_p%stepsPerExchange)
end subroutine onlyNode
