#include "../defines.inc"
! file defines routines useful for parallel tempering continuous variables over
! MPI. Thus, only compiled if MPI is available.
#if MPI_VERSION

subroutine startWorker(wlc_p, wlc_d)
use params
use mpi
! Override initialization with parallel setup parameters
!  In particualar it changes: wlc_p%AB, wlc_p%rep, wlc_p%MU, wlc_d%repSuffix
    Implicit none
    type(wlcsim_params), intent(inout) :: wlc_p
    type(wlcsim_data), intent(inout) :: wlc_d
    integer (kind = 4) dest ! message destination
    integer (kind = 4) source ! message source
    integer (kind = 4) id, nThreads
    integer (kind = 4) error  ! error id for MIP functions
    character(MAXFILENAMELEN) iostrg    ! for file naming
    integer ( kind = 4 ) status(MPI_status_SIZE) ! MPI stuff
    integer, parameter :: nTerms = 9  ! number of energy terms
    real(dp) cof(nTerms)

    call MPI_COMM_SIZE(MPI_COMM_WORLD,nThreads,error)
    call MPI_COMM_RANK(MPI_COMM_WORLD,id,error)
    if (wlc_d%id .ne. id) then
        print*, "ID mismatch in PT_override!"
        stop
    endif
    if (nThreads.lt.3) then
        wlc_d%repSuffix = "v1"
        wlc_d%rep = 1
        !wlc_d%id = int(id)
        print*, "No PT_override. Input values used."
        return
    endif

    !---------------------------------------------
    !
    !     Quenched Disorder must be same!
    !     Copy from replica 1 to others.
    !
    !----------------------------------------------
    if (WLC_P__FIELD_INT_ON) then
        if (wlc_d%id.eq.1) then
            do dest = 2,nThreads-1
                if(WLC_P__VARIABLE_CHEM_STATE) then
                    if (.not. WLC_P__CHEM_SEQ_FROM_FILE) then
                        call MPI_Send (wlc_d%METH,wlc_p%NT, MPI_integer, dest,   0, &
                                       MPI_COMM_WORLD,error )
                    endif
                else
                    call MPI_Send (wlc_d%AB,wlc_p%NT, MPI_integer, dest,   0, &
                                   MPI_COMM_WORLD,error )
                endif
            enddo
        elseif (wlc_d%id.gt.1) then
            source = 1
            if(WLC_P__VARIABLE_CHEM_STATE) then
                if (.not. WLC_P__CHEM_SEQ_FROM_FILE) then
                    call MPI_Recv (wlc_d%METH, wlc_p%NT, MPI_integer, source, 0, &
                                   MPI_COMM_WORLD, status, error )
                endif
            else
                call MPI_Recv (wlc_d%AB, wlc_p%NT, MPI_integer, source, 0, &
                               MPI_COMM_WORLD, status, error )
            endif
        endif
    endif

    !---------------------------------------------------
    !
    !   Receive instructions from head node
    !
    !----------------------------------------------------
    source = 0
    dest = 0
    call MPI_Recv ( wlc_d%rep, 1, MPI_integer, source, 0, &
      MPI_COMM_WORLD, status, error )

    call MPI_Recv ( cof, nTerms, MPI_doUBLE_PRECISION, source, 0, &
      MPI_COMM_WORLD, status, error )

    wlc_p%CHI      =cof(1)
    wlc_p%MU       =cof(2)
    wlc_p%HA      =cof(3)
    wlc_p%HP1_BIND =cof(4)
    wlc_p%KAP      =cof(5)
    !wlc_p%para(1)  =cof(6) ! eb, eperp, epar
    !wlc_p%para(2)  =cof(7)
    !wlc_p%para(3)  =cof(8)
    wlc_p%CHI_L2   =cof(9)

    write(iostrg,"(I4)") wlc_d%rep
    iostrg = adjustL(iostrg)
    iostrg = trim(iostrg)
    iostrg = "v"//trim(iostrg)
    iostrg = trim(iostrg)
    wlc_d%repSuffix = iostrg

end subroutine
subroutine replicaExchange(wlc_p, wlc_d)
! This checks in with the mpi head node to
! For parallel tempering of the form:  E = cof*x
! 1: Tell head node the x value
! 2: Recive replica assignment from head node
! 3: Recive assigned cof value
use params, only : wlcsim_params, wlcsim_data, dp, MAXFILENAMELEN, epsApprox
use mpi
    implicit none
    integer, parameter :: nTerms = 9  ! number of energy terms
    integer (kind = 4) id, error
    type(wlcsim_params), intent(inout) :: wlc_p
    type(wlcsim_data), intent(inout) :: wlc_d
    integer (kind = 4) dest ! message destination
    integer (kind = 4) source ! message source
    integer (kind = 4) nThreads
    character(MAXFILENAMELEN) iostr ! for handling sufix string
    integer status(MPI_status_SIZE)  ! MPI status
    real(dp) cof(nTerms)
    real(dp) chi_Old
    real(dp) mu_old
    real(dp) hA_Old
    real(dp) HP1_Bind_Old
    real(dp) Kap_Old
    real(dp) chi_l2_old
    real(dp) x(nTerms)
    real(dp) test(5)

    call MPI_COMM_SIZE(MPI_COMM_WORLD,nThreads,error)
    if (nThreads.lt.3) return

    x(1) = wlc_d%x_Chi
    x(2) = wlc_d%x_mu
    x(3) = wlc_d%x_Field
    x(4) = wlc_d%x_couple
    x(5) = wlc_d%x_kap
    x(6) = 0.0_dp !x(6) = wlc_p%EElas(1)/wlc_p%para(1)
    x(7) = 0.0_dp !x(7) = wlc_p%EElas(2)/wlc_p%para(2)
    x(8) = 0.0_dp !x(8) = wlc_p%EElas(3)/wlc_p%para(3)
    x(9) = wlc_d%x_maierSaupe

    test(1) = wlc_d%EChi/wlc_p%CHI
    test(3) = wlc_d%EField/wlc_p%HA
    test(4) = wlc_d%ECouple/wlc_p%HP1_BIND
    test(5) = wlc_d%EKap/wlc_p%KAP

    chi_Old = wlc_p%CHI
    mu_old = wlc_p%MU
    hA_Old = wlc_p%HA
    HP1_Bind_Old = wlc_p%HP1_BIND
    Kap_Old = wlc_p%KAP
    chi_l2_old = wlc_p%CHI_L2

    ! send number bound to head node
    dest = 0
    call MPI_Send(x,nTerms,MPI_doUBLE_PRECISION,dest,0,MPI_COMM_WORLD,error)
    call MPI_COMM_RANK(MPI_COMM_WORLD,id,error)
    ! send ind to head node
    if (id.eq.1) then
        call MPI_Send(wlc_d%mc_ind,1,MPI_integer,dest,0,MPI_COMM_WORLD,error)
    endif
    ! hear back on which replica and it's mu value
    source = 0
    ! get new replica number
    call MPI_Recv(wlc_d%rep,1,MPI_integer,source,0, &
                  MPI_COMM_WORLD,status,error)
    ! get new mu value
    call MPI_Recv(cof,nTerms,MPI_doUBLE_PRECISION,source,0,&
                  MPI_COMM_WORLD,status,error)

    wlc_p%CHI      =cof(1)
    wlc_p%MU       =cof(2)
    wlc_p%HA       =cof(3)
    wlc_p%HP1_BIND =cof(4)
    wlc_p%KAP      =cof(5)
    !wlc_p%para(1)  =cof(6)
    !wlc_p%para(2)  =cof(7)
    !wlc_p%para(3)  =cof(8)
    wlc_p%CHI_L2 = cof(9)

    if (abs(wlc_d%EChi-wlc_d%x_chi*chi_old).gt.epsApprox) then
        print*, "Error in replicaExchange"
        print*, "wlc_p%EChi",wlc_d%EChi,"x(1)*CofOld(1)",wlc_d%x_chi*chi_old
        stop 1
    endif

    wlc_d%EChi    =wlc_d%EChi    +wlc_d%x_chi      *(wlc_p%CHI      -chi_old)
    wlc_d%EMu     =wlc_d%EMu     +wlc_d%x_mu       *(wlc_p%MU       -mu_old)
    wlc_d%EField  =wlc_d%EField  +wlc_d%x_field    *(wlc_p%HA       -hA_old)
    wlc_d%ECouple =wlc_d%ECouple +wlc_d%x_couple   *(wlc_p%HP1_BIND -HP1_Bind_Old)
    wlc_d%EKap    =wlc_d%EKap    +wlc_d%x_Kap      *(wlc_p%KAP      -Kap_Old)
   ! wlc_p%EElas(1) = wlc_p%EElas(1) + x(6)*(Cof(6)-CofOld(6))
   ! wlc_p%EElas(2) = wlc_p%EElas(2) + x(7)*(Cof(7)-CofOld(7))
   ! wlc_p%EElas(3) = wlc_p%EElas(3) + x(8)*(Cof(8)-CofOld(8))
    wlc_d%eMaierSaupe    =wlc_d%eMaierSaupe    +wlc_d%x_maierSaupe*(wlc_p%CHI_L2-chi_l2_old)

    if (abs(wlc_d%EChi-wlc_p%CHI*wlc_d%x_chi).gt.epsApprox) then
        print*, "Error in replicaExchange"
        print*, "wlc_p%EChi",wlc_d%EChi,"x(1)*Cof(1)",wlc_p%CHI*wlc_d%x_chi
        stop 1
    endif

    ! change output file sufix
    write(iostr,"(I4)") wlc_d%rep
    iostr = adjustL(iostr)
    iostr = trim(iostr)
    iostr = "v"//trim(iostr)
    iostr = trim(iostr)
    wlc_d%repSuffix = iostr

end subroutine

#endif
