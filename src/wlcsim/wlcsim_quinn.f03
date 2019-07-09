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
        call head_node(wlc_p,wlc_numProcesses)
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
subroutine head_node(wlc_p,process)
! values from wlcsim_data
use params, only: wlc_mc_ind, wlc_rand_stat
    use mersenne_twister
    use mpi
    use params
!   MPI variables
    implicit none
    integer, intent(in) :: process ! number of therads
    integer ( kind = 4 ) dest   !destination id for messages
    integer ( kind = 4 ) source  !source id for messages
    integer ( kind = 4 ) error  ! error id for MIP functions
    integer ( kind = 4 ) status(MPI_status_SIZE) ! MPI stuff
    type(wlcsim_params), intent(inout) :: wlc_p

    !   variable for random number generator seeding
    real(dp) urand(1)

!   for head node use only variables
    integer rep ! physical replica number, for loops
    integer temp ! for castling
    logical keepGoing   ! set to false when NaN encountered
    integer, parameter :: nTerms = 10  ! number of energy terms
    real(dp) x(nTerms) ! slice of xMtrx
    real(dp) cof(nTerms) ! slice of cofMtrx
    integer N_average      ! number of attempts since last average
    integer upSuccess(process-1)  ! number of successes since last average
    integer downSuccess(process-1) ! number of successes since last average
    integer nExchange ! total number of exchanges attemted
    real(dp) energy ! for deciding to accept exchange
    integer term ! for loopin over terms
    real(dp) h_path,chi_path,mu_path,kap_path,HP1_Bind_path,maierSaupe_path,AEF_path ! functions
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
        s_vals(rep) = WLC_P__INITIAL_MAX_S*(dble(rep)-1.0_dp)/(dble(nPTReplicas)-1.0_dp)
    enddo

    ! Set initial values for parrallel tempering values
    do rep = 1,nPTReplicas
        if (WLC_P__PT_CHI) then
            cofMtrx(rep,1) = chi_path(s_vals(rep))
        else
            cofMtrx(rep,1) = wlc_p%CHI
        endif
        if (WLC_P__PT_MU) then
            cofMtrx(rep,2) = mu_path(s_vals(rep))
        else
            cofMtrx(rep,2) = wlc_p%MU
        endif
        if (WLC_P__PT_H) then
            cofMtrx(rep,3) = h_path(s_vals(rep))
        else
            cofMtrx(rep,3) = wlc_p%HA
        endif
        if (WLC_P__PT_COUPLE) then
            cofMtrx(rep,4) = HP1_Bind_path(s_vals(rep))
        else
            cofMtrx(rep,4) = wlc_p%HP1_BIND
        endif
        if (WLC_P__PT_KAP) then
            cofMtrx(rep,5) = kap_path(s_vals(rep))
        else
            cofMtrx(rep,5) = wlc_p%KAP
        endif
        cofMtrx(rep,6) = 0.0_dp
        cofMtrx(rep,7) = 0.0_dp
        cofMtrx(rep,8) = 0.0_dp
        if (WLC_P__PT_MAIERSAUPE) then
            cofMtrx(rep,9) = maierSaupe_path(s_vals(rep))
        else
            cofMtrx(rep,9) = wlc_p%CHI_L2
        endif
        if (WLC_P__PT_AEF) then
            cofMtrx(rep,10) = AEF_path(s_vals(rep))
        else
            cofMtrx(rep,10) = wlc_p%AEF
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
            if (WLC_P__RESTART.and.wlc_mc_ind.eq.1) then
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
        call MPI_Recv (wlc_mc_ind, 1, MPI_integer, source, 0, &
                       MPI_COMM_WORLD, status, error )

        ! do replica exchange
        do rep = 1,(nPTReplicas-1)
            energy = 0.0_dp
            do term = 1,nTerms
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
            endif
        enddo


        ! track/adapt acceptance rates
        N_average = N_average + 1
        if (N_average.ge.WLC_P__NREPADAPT) then
            call save_repHistory(upSuccess,downSuccess,nPTReplicas, &
                                 cofMtrx,xMtrx,nodeNumber,N_average,&
                                 nExchange,wlc_mc_ind,nTerms,s_vals)

            if ((wlc_mc_ind.ge.WLC_P__INDSTARTREPADAPT).and. &
                (wlc_mc_ind.lt.WLC_P__INDENDREPADAPT)) then ! insert input defined location here
                call adaptCof(downSuccess,nPTReplicas,s_vals,N_average,&
                               WLC_P__LOWERREPEXE,WLC_P__UPPERREPEXE,&
                               WLC_P__LOWERCOFRAIL,WLC_P__UPPERCOFRAIL,&
                               WLC_P__REPANNEALSPEED,WLC_P__REPLICABOUNDS)
                do rep = 1,nPTReplicas
                    if (WLC_P__PT_CHI) then
                        cofMtrx(rep,1) = chi_path(s_vals(rep))
                    endif
                    if (WLC_P__PT_MU) then
                        cofMtrx(rep,2) = mu_path(s_vals(rep))
                    endif
                    if (WLC_P__PT_H) then
                        cofMtrx(rep,3) = h_path(s_vals(rep))
                    endif
                    if (WLC_P__PT_COUPLE) then
                        cofMtrx(rep,4) = HP1_Bind_path(s_vals(rep))
                    endif
                    if (WLC_P__PT_KAP) then
                        cofMtrx(rep,5) = kap_path(s_vals(rep))
                    endif
                    if (WLC_P__PT_MAIERSAUPE) then
                        cofMtrx(rep,9) = maiersaupe_path(s_vals(rep))
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
        chi = -s
    endif
end function chi_path
function h_path(s) result(h)
    use params, only: dp
    implicit none
    real(dp), intent(in) :: s
    real(dp) h
    h = 1.0_dp*s
end function h_path
function mu_path(s) result(mu)
    use params, only: dp
    implicit none
    real(dp), intent(in) :: s
    real(dp) mu
    mu = s-2.0_dp
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
    output = s*(-1.0_dp)
end function maierSaupe_path
function hp1_bind_path(s) result(hp1_bind)
    use params, only: dp
    implicit none
    real(dp), intent(in) :: s
    real(dp) hp1_bind
    hp1_bind = s
end function hp1_bind_path
function AEF_path(s) result(AEF)
    use params, only: dp
    implicit none
    real(dp), intent(in) :: s
    real(dp) AEF
    AEF = -1.0_dp*s
end function AEF_path

#if MPI_VERSION
subroutine worker_node(wlc_p)
! values from wlcsim_data
use params, only: wlc_x_ExternalField, wlc_EmaierSaupe, wlc_deelas, wlc_dx_couple, wlc_x_Chi&
    , wlc_dECouple, wlc_debind, wlc_x_externalField, wlc_dEExternalField, wlc_deMaierSaupe, wlc_x_couple&
    , wlc_EField, wlc_ebind, wlc_dx_Kap, wlc_dx_externalField, wlc_ind_exchange, wlc_x_Field&
    , wlc_dx_Chi, wlc_EKap, wlc_x_maierSaupe, wlc_dECHI, wlc_dx_maierSaupe, wlc_dEKap&
    , wlc_dEField, wlc_eExternalField, wlc_ECouple, wlc_ECHI, wlc_x_mu, wlc_eMu&
    , wlc_mc_ind, wlc_deMu, wlc_dx_Field, wlc_eelas, wlc_dx_mu, wlc_x_Kap&
    , wlc_EMaierSaupe
    use mpi
    use params
    use mersenne_twister
    implicit none
    integer ( kind = 4 ), save :: id = -1     ! which processor I am
    integer ( kind = 4 ) error  ! error id for MIP functions
    type(wlcsim_params), intent(inout) :: wlc_p
    integer i
    logical system_has_been_changed
    real :: start, finish

    if (id == -1) then
        call MPI_Comm_rank(MPI_COMM_WORLD, id, error)
        call stop_if_err(error, "Failed to get num_processes.")
    endif
    if (wlc_mc_ind == 1) then
        if (WLC_P__PT_CHI .or. WLC_P__PT_H .or. WLC_P__PT_KAP .or. WLC_P__PT_MU .or. WLC_P__PT_COUPLE .or. WLC_P__ENSEMBLE_BIND) then
            call startWorker(wlc_p)
        endif
    endif

    call schedule(wlc_p,system_has_been_changed)

    if (system_has_been_changed) then
        call CalculateEnergiesFromScratch(wlc_p)
        wlc_eelas = wlc_deelas
        if (wlc_p%field_int_on_currently) then
            wlc_ECouple =wlc_dECouple
            wlc_EKap    =wlc_dEKap
            wlc_ECHI    =wlc_dECHI
            wlc_EField  =wlc_dEField
            wlc_EMaierSaupe = wlc_deMaierSaupe
            wlc_x_Field =wlc_dx_Field
            wlc_x_couple = wlc_dx_couple
            wlc_x_Kap   =wlc_dx_Kap
            wlc_x_Chi   =wlc_dx_Chi
            wlc_x_maierSaupe = wlc_dx_maierSaupe
        else
            wlc_ECouple =0.0_dp
            wlc_EKap    =0.0_dp
            wlc_ECHI    =0.0_dp
            wlc_EField  =0.0_dp
            wlc_EmaierSaupe = 0.0_dp
            wlc_x_Field =0.0_dp
            wlc_x_couple = 0.0_dp
            wlc_x_Kap   =0.0_dp
            wlc_x_Chi   =0.0_dp
            wlc_x_maierSaupe = 0.0_dp
        endif
        if (WLC_P__VARIABLE_CHEM_STATE) then
            wlc_ebind   =wlc_debind
            wlc_eMu     =wlc_deMu
            wlc_x_mu    =wlc_dx_mu
        else
            wlc_ebind   =0.0_dp
            wlc_eMu     =0.0_dp
            wlc_x_mu    =0.0_dp
        endif
        if(WLC_P__APPLY_EXTERNAL_FIELD) then
            wlc_eExternalField = wlc_dEExternalField
            wlc_x_externalField = wlc_dx_externalField
            if (abs(wlc_eExternalField-wlc_p%AEF*wlc_x_ExternalField).gt.0.00001) then
                print*, "error in wlcsim_quinn"
                stop
            endif
        else
            wlc_eExternalField = 0.0_dp
            wlc_x_externalField = 0.0_dp
        endif
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
        call replicaExchange(wlc_p)
    enddo
    call cpu_time(finish)
    print*, "Save Point time", finish-start, " seconds"
end subroutine worker_node
#endif

subroutine onlyNode(wlc_p)
! values from wlcsim_data
use params, only: wlc_x_maiersaupe, wlc_dx_couple, wlc_x_Chi, wlc_dECouple, wlc_debind&
    , wlc_dEExternalField, wlc_x_couple, wlc_ebind, wlc_EField, wlc_dx_Kap, wlc_dEmaiersaupe&
    , wlc_x_Field, wlc_dx_Chi, wlc_EKap, wlc_dx_maiersaupe, wlc_dECHI, wlc_dEKap&
    , wlc_dEField, wlc_eExternalField, wlc_ECouple, wlc_ECHI, wlc_x_mu, wlc_eMu&
    , wlc_Emaiersaupe, wlc_deMu, wlc_dx_Field, wlc_dx_mu, wlc_x_Kap
    use params
    implicit none
    type(wlcsim_params), intent(inout) :: wlc_p
    logical system_has_been_changed
    real :: start, finish
    !   * Perform a MC simulation *
    call schedule(wlc_p,system_has_been_changed)
    if (system_has_been_changed) then
        call CalculateEnergiesFromScratch(wlc_p)
        if (wlc_p%field_int_on_currently) then
            wlc_ECouple =wlc_dECouple
            wlc_EKap    =wlc_dEKap
            wlc_ECHI    =wlc_dECHI
            wlc_EField  =wlc_dEField
            wlc_Emaiersaupe = wlc_dEmaiersaupe
            wlc_x_Field =wlc_dx_Field
            wlc_x_maiersaupe = wlc_dx_maiersaupe
            wlc_x_couple = wlc_dx_couple
            wlc_x_Kap   =wlc_dx_Kap
            wlc_x_Chi   =wlc_dx_Chi
        else
            wlc_ECouple =0.0_dp
            wlc_EKap    =0.0_dp
            wlc_ECHI    =0.0_dp
            wlc_EField  =0.0_dp
            wlc_Emaiersaupe = 0.0_dp
            wlc_x_Field =0.0_dp
            wlc_x_couple = 0.0_dp
            wlc_x_Kap   =0.0_dp
            wlc_x_Chi   =0.0_dp
            wlc_x_maiersaupe = 0.0_dp
        endif
        if (WLC_P__VARIABLE_CHEM_STATE) then
            wlc_ebind   =wlc_debind
            wlc_eMu     =wlc_deMu
            wlc_x_mu    =wlc_dx_mu
        else
            wlc_ebind   =0.0_dp
            wlc_eMu     =0.0_dp
            wlc_x_mu    =0.0_dp
        endif
        if(WLC_P__APPLY_EXTERNAL_FIELD) then
            wlc_eExternalField = wlc_dEExternalField
        else
            wlc_eExternalField = 0.0_dp
        endif
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
    use params
    implicit none
    type(wlcsim_params), intent(inout) :: wlc_p
    logical, intent(out) :: system_has_been_changed

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
    if(wlc_mc_ind.lt.WLC_P__N_KAP_ON) then
        wlc_p%KAP_ON = 0.0_dp
    else
        if (abs(wlc_p%KAP_ON) < eps) system_has_been_changed = .TRUE.
        wlc_p%KAP_ON = 1.0_dp
    endif

    if(wlc_mc_ind.lt.WLC_P__N_CHI_ON) then
        wlc_p%CHI_ON = 0.0_dp
    else
        if (abs(wlc_p%CHI_ON) < eps) system_has_been_changed = .TRUE.
        wlc_p%CHI_ON = 1.0_dp
    endif

    if(wlc_mc_ind.lt.WLC_P__N_EXTERNAL_ON) then
        wlc_p%AEF = 0.0_dp
    else
        if (wlc_mc_ind == WLC_P__N_EXTERNAL_ON) system_has_been_changed = .TRUE.
        wlc_p%AEF = WLC_P__AmplitudeExternalField
    endif

    if(wlc_mc_ind.lt.WLC_P__N_CHI_L2_ON) then
        wlc_p%CHI_L2_ON = .False.
    else
        if (.not. wlc_p%CHI_L2_ON) system_has_been_changed = .TRUE.
        wlc_p%CHI_L2_ON = .True.
    endif

    if ((wlc_mc_ind.gt.WLC_P__INDSTARTREPADAPT).and. &
        (wlc_mc_ind.le.WLC_P__INDENDREPADAPT)) then ! addapt Cof was run
        system_has_been_changed = .TRUE.
    endif
end subroutine
