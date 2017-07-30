! file defines routines useful for parallel tempering continuous variables over
! MPI. Thus, only compiled if MPI is available.
#if MPI_VERSION

subroutine startWorker(wlc_p, wlc_d)
use params
use mpi
! Override initialization with parallel setup parameters
!  In particualar it changes: wlc_p%AB, wlc_p%rep, wlc_p%mu, wlc_d%repSuffix
    Implicit none
    type(wlcsim_params), intent(inout) :: wlc_p
    type(wlcsim_data), intent(inout) :: wlc_d
    integer (kind = 4) dest ! message destination
    integer (kind = 4) source ! message source
    integer (kind = 4) id, nThreads,ierror
    integer (kind = 4) error  ! error id for MIP functions
    character(MAXFILENAMELEN) iostrg    ! for file naming
    integer ( kind = 4 ) status(MPI_status_SIZE) ! MPI stuff
    integer, parameter :: nTerms = 8  ! number of energy terms
    real(dp) cof(nTerms)

    call MPI_COMM_SIZE(MPI_COMM_WORLD,nThreads,ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD,id,ierror)
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
    if (wlc_d%id.eq.1) then
        do dest = 2,nThreads-1
            if(wlc_p%bind_On) then
                call MPI_Send (wlc_d%METH,wlc_p%NT, MPI_integer, dest,   0, &
                               MPI_COMM_WORLD,error )
            else
                call MPI_Send (wlc_d%AB,wlc_p%NT, MPI_integer, dest,   0, &
                               MPI_COMM_WORLD,error )
            endif
        enddo
    elseif (wlc_d%id.gt.1) then
        source = 1
        if(wlc_p%bind_On) then
            call MPI_Recv (wlc_d%METH, wlc_p%NT, MPI_integer, source, 0, &
                           MPI_COMM_WORLD, status, error )
        else
            call MPI_Recv (wlc_d%AB, wlc_p%NT, MPI_integer, source, 0, &
                           MPI_COMM_WORLD, status, error )
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

    wlc_p%chi      =cof(1)
    wlc_p%mu       =cof(2)
    wlc_p%hA      =cof(3)
    wlc_p%HP1_Bind =cof(4)
    wlc_p%Kap      =cof(5)
    !wlc_p%para(1)  =cof(6)
    !wlc_p%para(2)  =cof(7)
    !wlc_p%para(3)  =cof(8)

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
use params
use mpi
    implicit none
    integer, parameter :: nTerms = 8  ! number of energy terms
    integer (kind = 4) id, ierror
    type(wlcsim_params), intent(inout) :: wlc_p
    type(wlcsim_data), intent(inout) :: wlc_d
    integer i  ! working intiger
    integer (kind = 4) dest ! message destination
    integer (kind = 4) source ! message source
    integer (kind = 4) nThreads
    character(MAXFILENAMELEN) iostr ! for handling sufix string
    integer status(MPI_status_SIZE)  ! MPI status
    real(dp) cof(nTerms)
    real(dp) cofOld(nTerms)
    real(dp) x(nTerms)
    real(dp) test(5)
    logical isfile

    call MPI_COMM_SIZE(MPI_COMM_WORLD,nThreads,ierror)
    if (nThreads.lt.3) return

    x(1) = wlc_d%x_Chi
    x(2) = wlc_d%x_mu
    x(3) = wlc_d%x_Field
    x(4) = wlc_d%x_couple
    x(5) = wlc_d%x_kap
    x(6) = 0.0_dp !x(6) = wlc_p%EElas(1)/wlc_p%para(1)
    x(7) = 0.0_dp !x(7) = wlc_p%EElas(2)/wlc_p%para(2)
    x(8) = 0.0_dp !x(8) = wlc_p%EElas(3)/wlc_p%para(3)

    test(1) = wlc_d%EChi/wlc_p%Chi
    test(3) = wlc_d%EField/wlc_p%hA
    test(4) = wlc_d%ECouple/wlc_p%HP1_Bind
    test(5) = wlc_d%EKap/wlc_p%Kap

    cofOld(1) = wlc_p%chi
    cofOld(2) = wlc_p%mu
    cofOld(3) = wlc_p%hA
    cofOld(4) = wlc_p%HP1_Bind
    cofOld(5) = wlc_p%Kap
    cofOld(6) = 0!wlc_p%para(1)
    cofOld(7) = 0!wlc_p%para(2)
    cofOld(8) = 0!wlc_p%para(3)

    do i = 1,5
        if (i.eq.2) cycle ! doesn't work for mu
        if (abs(cofOld(I)*(test(I)-x(I))).lt.0.0001) cycle
        if (abs(cofOld(I)).lt.0.00000001) cycle
        inquire(file = "data/error", exist = isfile)
        if (isfile) then
            open (unit = 1, file = "data/error", status ='OLD', POSITION = "append")
        else
            open (unit = 1, file = "data/error", status = 'new')
        endif
        write(1,*) "Error in replicaExchange"
        write(1,*) "I",I," test",test(I)," x",x(I)," cof",cofOld(I)
        print*, "Error in replicaExchange"
        print*, "I",I," test",test(I)," x",x(I)," cof",cofOld(I)
        close (1)
    enddo

    ! send number bound to head node
    dest = 0
    call MPI_Send(x,nTerms,MPI_doUBLE_PRECISION,dest,0,MPI_COMM_WORLD,ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD,id,ierror)
    ! send ind to head node
    if (id.eq.1) then
        call MPI_Send(wlc_d%mc_ind,1,MPI_integer,dest,0,MPI_COMM_WORLD,ierror)
    endif
    ! hear back on which replica and it's mu value
    source = 0
    ! get new replica number
    call MPI_Recv(wlc_d%rep,1,MPI_integer,source,0, &
                  MPI_COMM_WORLD,status,ierror)
    ! get new mu value
    call MPI_Recv(cof,nTerms,MPI_doUBLE_PRECISION,source,0,&
                  MPI_COMM_WORLD,status,ierror)

    wlc_p%chi      =cof(1)
    wlc_p%mu       =cof(2)
    wlc_p%hA       =cof(3)
    wlc_p%HP1_Bind =cof(4)
    wlc_p%Kap      =cof(5)
    !wlc_p%para(1)  =cof(6)
    !wlc_p%para(2)  =cof(7)
    !wlc_p%para(3)  =cof(8)

    if (abs(wlc_d%EChi-x(1)*CofOld(1)).gt.0.0000001_dp) then
        print*, "Error in replicaExchange"
        print*, "wlc_p%EChi",wlc_d%EChi,"x(1)*CofOld(1)",x(1)*CofOld(1)
        stop 1
    endif

    wlc_d%EChi    =wlc_d%EChi    +x(1)*(Cof(1)-CofOld(1))
    wlc_d%ebind   =wlc_d%ebind   +x(2)*(Cof(2)-CofOld(2))
    wlc_d%EField  =wlc_d%EField  +x(3)*(Cof(3)-CofOld(3))
    wlc_d%ECouple =wlc_d%ECouple +x(4)*(Cof(4)-CofOld(4))
    wlc_d%EKap    =wlc_d%EKap    +x(5)*(Cof(5)-CofOld(5))
   ! wlc_p%EElas(1) = wlc_p%EElas(1) + x(6)*(Cof(6)-CofOld(6))
   ! wlc_p%EElas(2) = wlc_p%EElas(2) + x(7)*(Cof(7)-CofOld(7))
   ! wlc_p%EElas(3) = wlc_p%EElas(3) + x(8)*(Cof(8)-CofOld(8))

    if (abs(wlc_d%EChi-x(1)*Cof(1)).gt.0.000001_dp) then
        print*, "Error in replicaExchange"
        print*, "wlc_p%EChi",wlc_d%EChi,"x(1)*Cof(1)",x(1)*Cof(1)
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
