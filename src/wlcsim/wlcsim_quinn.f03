
subroutine wlcsim_quinn(save_ind, mc, md)
    use mpi
    use params
    implicit none
    integer, intent(in) :: save_ind ! 1, 2, ...
    type(wlcsim_params), intent(inout) :: mc
    type(wlcsim_data), intent(inout) :: md
    integer, save :: id, num_processes
    integer (kind=4) error

    if (save_ind == 1) then
        call MPI_Init(error)
        call stop_if_err(error, "Failed to MPI_Init.")
        call MPI_Comm_size(MPI_COMM_WORLD, num_processes, error)
        call stop_if_err(error, "Failed to get num_processes.")
        call MPI_Comm_rank(MPI_COMM_WORLD, id, error)
        call stop_if_err(error, "Failed to get num_processes.")
    endif

    md%mc_ind=save_ind

    if (num_processes == 1) then
        call onlyNode(mc,md)
    elseif (id == 0) then
        call head_node(mc,md,num_processes)
    else
        call worker_node(mc,md)
    endif



    print*, '________________________________________'
    print*, 'Time point ',save_ind, ' out of', mc%numSavePoints, 'Thread id', id
    call printEnergies(md)
    call printWindowStats(mc,md)
    !call wlcsim_params_printPhi(mc,md)


end subroutine wlcsim_quinn

subroutine head_node(mc,md,p)
    use mersenne_twister
    use mpi
    use params
!   MPI variables
    implicit none
    integer, intent(in) :: p ! number of therads
    integer ( kind = 4 ) dest   !destination id for messages
    integer ( kind = 4 ) source  !source id for messages
    integer ( kind = 4 ) id     ! which processor I am
    integer ( kind = 4 ) error  ! error id for MIP functions
    integer ( kind = 4 ) status(MPI_status_SIZE) ! MPI stuff
    type(wlcsim_params), intent(inout) :: mc
    type(wlcsim_data), intent(inout) :: md
    integer, allocatable :: nodeNumber(:)  ! list of which nodes are which
    real(dp), allocatable :: xMtrx(:,:)  ! sum of bound states
    real(dp), allocatable :: cofMtrx(:,:) ! mu or chi or whatever
    real(dp), allocatable :: s_vals(:) ! path parameter

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
    integer, parameter :: nTerms=8  ! number of energy terms
    double precision x(nTerms) ! slice of xMtrx
    double precision cof(nTerms) ! slice of cofMtrx
    integer N_average      ! number of attempts since last average
    integer upSuccess(p-1)  ! number of successes since last average
    integer downSuccess(p-1) ! number of successes since last average
    integer nExchange ! total number of exchanges attemted
    double precision energy ! for deciding to accept exchange
    integer term ! for loopin over terms
    double precision h_path,chi_path,mu_path,kap_path,HP1_Bind_path ! functions
    integer nPTReplicas

    id=0 ! head node
    nPTReplicas=p-1
    ! -----------------------------------------------
    !
    !   Generate thread safe random number seeds
    !
    !--------------------------------------------
    if (.false.) then ! set spedific seed
        Irand=7171
    else ! seed from clock
        call date_and_time(datedum,timedum,zonedum,seedvalues)
        Irand=int(-seedvalues(5)*1E7-seedvalues(6)*1E5 &
                  -seedvalues(7)*1E3-seedvalues(8))
        Irand=mod(Irand,10000)
        print*, "Random Intiger seed:",Irand
    endif
    call random_setseed(Irand*(id+1),rand_stat) ! random seed for head node
    do dest=1,nPTReplicas ! send out the others
        call MPI_Send (Irand,1, MPI_integer, dest,   0, &
                        MPI_COMM_WORLD,error )
    enddo
    ! -------------------------
    !
    !   innitialize
    !
    ! --------------------------
    nPTReplicas = p-1;

    allocate( xMtrx(nPTReplicas,nTerms))
    allocate( cofMtrx(nPTReplicas,nTerms))
    allocate( nodeNumber(nPTReplicas))
    allocate( s_vals(nPTReplicas))

    do rep=1,nPTReplicas
        upSuccess(rep)=0
        downSuccess(rep)=0
        s_vals(rep)=mc%INITIAL_MAX_S*dble(rep)/dble(nPTReplicas)
    enddo

    do rep=1,nPTReplicas
        if (mc%PT_chi) then
            cofMtrx(rep,1)=chi_path(s_vals(rep))
        else
            cofMtrx(rep,1)=mc%chi
        endif
        if (mc%PT_mu) then
            cofMtrx(rep,2)=mu_path(s_vals(rep))
        else
            cofMtrx(rep,2)=mc%mu
        endif
        if (mc%PT_h) then
            cofMtrx(rep,3)=h_path(s_vals(rep))
        else
            cofMtrx(rep,3)=mc%hA
        endif
        if (mc%PT_couple) then
            cofMtrx(rep,4)=HP1_Bind_path(s_vals(rep))
        else
            cofMtrx(rep,4)=mc%HP1_Bind
        endif
        if (mc%PT_Kap) then
            cofMtrx(rep,5)=kap_path(s_vals(rep))
        else
            cofMtrx(rep,5)=mc%KAP
        endif
        cofMtrx(rep,6)=0
        cofMtrx(rep,7)=0
        cofMtrx(rep,8)=0
    enddo

    N_average=0

    ! Initially replica numbers are same as nodes
    do rep=1,nPTReplicas
        nodeNumber(rep)=rep
    enddo
    !-----------------------------------------------------
    !
    !    Begin Main loop
    !
    !------------------------------------------------------
    keepGoing=.True.
    do while(keepGoing)
        ! give workers thier jobs
        do rep=1,nPTReplicas
            dest=nodeNumber(rep)
            call MPI_Send (rep,1, MPI_integer, dest,   0, &
                            MPI_COMM_WORLD,error )
            if (mc%restart.and.md%mc_ind.eq.1) then
                source=dest
                call MPI_Recv (cof, nTerms, MPI_doUBLE_PRECISION, source, 0, &
                               MPI_COMM_WORLD, status, error )
                cofMtrx(rep,:)=cof
            else
                cof=cofMtrx(rep,:)
                call MPI_Send (cof,nTerms, MPI_doUBLE_PRECISION, dest,   0, &
                                MPI_COMM_WORLD,error )
            endif
        enddo
        ! get results from workers

        do rep=1,nPTReplicas
            source=nodeNumber(rep)
            call MPI_Recv ( x, nTerms, MPI_doUBLE_PRECISION, source, 0, &
                           MPI_COMM_WORLD, status, error )
            xMtrx(rep,:)=x
            if(isnan(x(1))) then ! endo of program
                keepGoing=.false.
                return
            endif
        enddo

        source=1
        call MPI_Recv (md%mc_ind, 1, MPI_integer, source, 0, &
                       MPI_COMM_WORLD, status, error )

        ! do replica exchange
        do rep=1,(nPTReplicas-1)
            energy=0.0_dp
            do term=1,nTerms
                energy=energy-(xMtrx(rep+1,term)-xMtrx(rep,term))*&
                              (cofMtrx(rep+1,term)-cofMtrx(rep,term))
            enddo
            call random_number(urand,rand_stat)
            if (exp(-1.0_dp*energy).gt.urand(1)) then
                if (mc%PTON) then
                    temp=nodeNumber(rep)
                    nodeNumber(rep)=nodeNumber(rep+1)
                    nodeNumber(rep+1)=temp
                endif
                upSuccess(rep)=upSuccess(rep)+1
                downSuccess(rep+1)=downSuccess(rep+1)+1
                x=xMtrx(rep,:)
                xMtrx(rep,:)=xMtrx(rep+1,:)
                xMtrx(rep+1,:)=x
            endif
        enddo


        ! track/adapt acceptance rates
        N_average=N_average+1
        if (N_average.ge.mc%NRepAdapt) then
            call save_repHistory(upSuccess,downSuccess,nPTReplicas, &
                                 cofMtrx,xMtrx,nodeNumber,N_average,&
                                 nExchange,md%mc_ind,nTerms,s_vals)

            if ((md%mc_ind.ge.mc%indStartRepAdapt).and. &
                (md%mc_ind.lt.mc%indendRepAdapt)) then ! insert input defined location here
                call adaptCof(downSuccess,nPTReplicas,s_vals,N_average,&
                               mc%lowerRepExe,mc%upperRepExe,&
                               mc%lowerCofRail,mc%upperCofRail,&
                               mc%RepAnnealSpeed,mc%replicaBounds)
                do rep=1,nPTReplicas
                    if (mc%PT_chi) then
                        cofMtrx(rep,1)=chi_path(s_vals(rep))
                    endif
                    if (mc%PT_mu) then
                        cofMtrx(rep,2)=mu_path(s_vals(rep))
                    endif
                    if (mc%PT_h) then
                        cofMtrx(rep,3)=h_path(s_vals(rep))
                    endif
                    if (mc%PT_couple) then
                        cofMtrx(rep,4)=HP1_Bind_path(s_vals(rep))
                    endif
                    if (mc%PT_Kap) then
                        cofMtrx(rep,5)=kap_path(s_vals(rep))
                    endif
                enddo
            endif
            N_average=0
            do rep=1,nPTReplicas
                upSuccess(rep)=0
                downSuccess(rep)=0
            enddo
        endif
        nExchange=nExchange+1
    enddo

    deallocate(xMtrx)
    deallocate(cofMtrx)
    deallocate(nodeNumber)
    deallocate(s_vals)

end subroutine head_node
function chi_path(s) result(chi)
    use params, only: dp
    implicit none
    real(dp), intent(in) :: s
    real(dp) chi
    real(dp) chi_max
    if (.false.) then
        chi_max=2.00_dp
        if (s.lt.0.5_dp) then
            chi=0.0_dp
        else
            chi=chi_max*2.0_dp*(s-0.5_dp)
        endif
    else
        chi=s
    endif
end function chi_path
function h_path(s) result(h)
    use params, only: dp
    implicit none
    real(dp), intent(in) :: s
    real(dp) h
    real(dp) h_max
    if (.false.) then
        h_max=10.0_dp
        if (s.lt.0.5_dp) then
             h=h_max*s*2.0_dp
        else
             h=h_max*(1.0_dp-s)*2.0_dp
        endif
    else
        h=0.0_dp
    endif
end function h_path
function mu_path(s) result(mu)
    use params, only: dp
    implicit none
    real(dp), intent(in) :: s
    real(dp) mu
    mu=s
end function mu_path
function kap_path(s) result(kap)
    use params, only: dp
    implicit none
    real(dp), intent(in) :: s
    real(dp) kap
    kap=s*10.0_dp
end function kap_path
function hp1_bind_path(s) result(hp1_bind)
    use params, only: dp
    implicit none
    real(dp), intent(in) :: s
    real(dp) hp1_bind
    hp1_bind=s
end function hp1_bind_path

subroutine worker_node(mc,md)
    use mpi
    use params
    use mersenne_twister
    integer ( kind = 4 ) dest   !destination id for messages
    integer ( kind = 4 ) source  !source id for messages
    integer ( kind = 4 ), save :: id = -1     ! which processor I am
    integer ( kind = 4 ) ierror  ! error id for MIP functions
    integer ( kind = 4 ) status(MPI_status_SIZE) ! MPI stuff
    type(wlcsim_params), intent(inout) :: mc
    type(wlcsim_data), intent(inout) :: md
    type(random_stat) rand_stat  ! state of random number chain

    logical system_has_been_changed
    system_has_been_changed=.False.
    if (id == -1) then
        call MPI_Comm_rank(MPI_COMM_WORLD, id, ierror)
        call stop_if_err(error, "Failed to get num_processes.")
    endif

    source = 0
    dest = 0
    ! -----------------------------------------------
    !
    !   Generate thread safe random number chain: rand_stat
    !
    !--------------------------------------------
    if (md%mc_ind==1) then
        call MPI_Recv ( Irand, 1, MPI_integer, source, 0, &
                        MPI_COMM_WORLD, status, ierror )
        call random_setseed(Irand*(id+1),rand_stat) ! random seed for head node
        if (mc%restart) then
            call pt_restart(mc,md)
        endif
        system_has_been_changed=.TRUE.
    endif
    ! ------------------------------
    !
    ! Different instructions for each save point
    !
    !  --------------------------------

    if (md%mc_ind.LE.mc%NNOINT) then
        mc%field_int_on=.false.
    else
        if (.not.mc%field_int_on)  system_has_been_changed=.TRUE.
        mc%field_int_on=.true.
    endif
    if(md%mc_ind.lt.mc%N_KAP_ON) then
        mc%KAP_ON=0.0_dp
    else
        if (mc%KAP_ON.eq.0.0_dp) system_has_been_changed=.TRUE.
        mc%KAP_ON=1.0_dp
    endif

    if(md%mc_ind.lt.mc%N_CHI_ON) then
        mc%CHI_ON=0.0_dp
    else
        if (mc%CHI_ON.eq.0.0_dp) system_has_been_changed=.TRUE.
        mc%CHI_ON=1.0_dp
    endif

    if ((md%mc_ind.gt.mc%indStartRepAdapt).and. &
        (md%mc_ind.le.mc%indendRepAdapt)) then ! addapt Cof was run
        system_has_been_changed=.TRUE.
    endif

    if (system_has_been_changed) then
        call CalculateEnergiesFromScratch(mc,md)
        if (mc%field_int_on) then
            md%ECouple =md%dECouple
            md%EKap    =md%dEKap
            md%ECHI    =md%dECHI
            md%EField  =md%dEField
            md%x_Field =md%dx_Field
            md%x_couple=md%dx_couple
            md%x_Kap   =md%dx_Kap
            md%x_Chi   =md%dx_Chi
        else
            md%ECouple =0.0_dp
            md%EKap    =0.0_dp
            md%ECHI    =0.0_dp
            md%EField  =0.0_dp
            md%x_Field =0.0_dp
            md%x_couple=0.0_dp
            md%x_Kap   =0.0_dp
            md%x_Chi   =0.0_dp
        endif
        if (mc%bind_On) then
            md%ebind   =md%debind
            md%x_mu    =md%dx_mu
        else
            md%ebind   =0.0_dp
            md%x_mu    =0.0_dp
        endif
    else
        call VerifyEnegiesFromScratch(mc, md)
    endif

    ! ------------------------------
    !
    ! call main simulation code
    !
    !  --------------------------------

    do i=1,mc%nReplicaExchangePerSavePoint

        !   * Perform a MC simulation *
        call MCsim(mc,md,mc%stepsPerExchange)

        !   * Replica Exchange *
        call replicaExchange(mc)

    enddo
end subroutine worker_node

subroutine onlyNode(mc,md)
    use params
    type(wlcsim_params), intent(inout) :: mc
    type(wlcsim_data), intent(inout) :: md
    !   * Perform a MC simulation *
    call VerifyEnegiesFromScratch(mc, md)
    call MCsim(mc,md,mc%nReplicaExchangePerSavePoint*mc%stepsPerExchange)
end subroutine onlyNode
