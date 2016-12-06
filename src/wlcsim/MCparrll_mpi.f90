subroutine PT_override(mc,md)
use params
use mpi
! Override initialization with parallel setup parameters
!  In particualar it changes: mc%AB, mc%rep, mc%mu, mc%repSuffix
    Implicit none
    type(wlcsim_params), intent(inout) :: mc
    type(wlcsim_data), intent(inout) :: md
    integer (kind=4) dest ! message destination
    integer (kind=4) source ! message source
    integer (kind=4) id, nThreads,ierror
    integer (kind=4) error  ! error id for MIP functions
    character(MAXFILENAMELEN) iostrg    ! for file naming
    integer ( kind = 4 ) status(MPI_status_SIZE) ! MPI stuff
    integer, parameter :: nTerms=8  ! number of energy terms
    real(dp) cof(nTerms)

    call MPI_COMM_SIZE(MPI_COMM_WORLD,nThreads,ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD,id,ierror)
    if (nThreads.lt.3) then
        mc%repSuffix="v1"
        md%rep=1
        md%id=int(id)
        print*, "No PT_override. Input values used."
        return
    endif

    !---------------------------------------------
    !
    !     Quenched Disorder must be same!
    !     Copy from replica 1 to others.
    !
    !----------------------------------------------
    if (id.eq.1) then
        do dest=2,nThreads-1
            if(mc%simtype.eq.1) then
                call MPI_Send (md%METH,mc%NT, MPI_integer, dest,   0, &
                               MPI_COMM_WORLD,error )
            elseif(mc%simtype.eq.0) then
                call MPI_Send (md%AB,mc%NT, MPI_integer, dest,   0, &
                               MPI_COMM_WORLD,error )
            else
                print*, "Error in PT_override. simtype doesn't exist."
                stop 1
            endif
        enddo
    else
        source=1
        if(mc%simtype.eq.1) then
            call MPI_Recv (md%METH, mc%NT, MPI_integer, source, 0, &
                           MPI_COMM_WORLD, status, error )
        elseif(mc%simtype.eq.0) then
            call MPI_Recv (md%AB, mc%NT, MPI_integer, source, 0, &
                           MPI_COMM_WORLD, status, error )
        else
            print*, "Error in PT_override. simtype doesn't exist."
            stop 1
        endif
    endif

    !---------------------------------------------------
    !
    !   Receive instructions from head node
    !
    !----------------------------------------------------
    source=0
    dest=0
    call MPI_Recv ( md%rep, 1, MPI_integer, source, 0, &
      MPI_COMM_WORLD, status, error )

    call MPI_Recv ( cof, nTerms, MPI_doUBLE_PRECISION, source, 0, &
      MPI_COMM_WORLD, status, error )

    mc%chi      =cof(1)
    mc%mu       =cof(2)
    mc%hA      =cof(3)
    mc%HP1_Bind =cof(4)
    mc%Kap      =cof(5)
    !mc%para(1)  =cof(6)
    !mc%para(2)  =cof(7)
    !mc%para(3)  =cof(8)

    write(iostrg,"(I4)") md%rep
    iostrg=adjustL(iostrg)
    iostrg=trim(iostrg)
    iostrg="v"//trim(iostrg)
    iostrg=trim(iostrg)
    mc%repSuffix=iostrg

    ! keep track of which thread you are
    md%id=int(id)
end subroutine
subroutine replicaExchange(mc,md)
! This checks in with the mpi head node to
! For parallel tempering of the form:  E=cof*x
! 1: Tell head node the x value
! 2: Recive replica assignment from head node
! 3: Recive assigned cof value
use params
use mpi
    IMPLICIT NONE
    integer, parameter :: nTerms=8  ! number of energy terms
    integer (kind=4) id, ierror
    type(wlcsim_params), intent(inout) :: mc
    type(wlcsim_data), intent(inout) :: md
    integer i  ! working intiger
    integer (kind=4) dest ! message destination
    integer (kind=4) source ! message source
    integer (kind=4) nThreads
    character(MAXFILENAMELEN) iostr ! for handling sufix string
    integer status(MPI_status_SIZE)  ! MPI status
    real(dp) cof(nTerms)
    real(dp) cofOld(nTerms)
    real(dp) x(nTerms)
    real(dp) test(5)
    logical isfile

    call MPI_COMM_SIZE(MPI_COMM_WORLD,nThreads,ierror)
    if (nThreads.lt.3) return

    x(1)=md%x_Chi
    x(2)=md%x_mu
    x(3)=md%x_Field
    x(4)=md%x_couple
    x(5)=md%x_kap
    x(6)=0.0_dp !x(6)=mc%EElas(1)/mc%para(1)
    x(7)=0.0_dp !x(7)=mc%EElas(2)/mc%para(2)
    x(8)=0.0_dp !x(8)=mc%EElas(3)/mc%para(3)

    test(1)=md%EChi/mc%Chi
    test(3)=md%EField/mc%hA
    test(4)=md%ECouple/mc%HP1_Bind
    test(5)=md%EKap/mc%Kap

    cofOld(1)=mc%chi
    cofOld(2)=mc%mu
    cofOld(3)=mc%hA
    cofOld(4)=mc%HP1_Bind
    cofOld(5)=mc%Kap
    cofOld(6)=0!mc%para(1)
    cofOld(7)=0!mc%para(2)
    cofOld(8)=0!mc%para(3)

    do i=1,5
        if (i.eq.2) cycle ! doesn't work for mu
        if (abs(cofOld(I)*(test(I)-x(I))).lt.0.0001) cycle
        if (abs(cofOld(I)).lt.0.00000001) cycle
        inquire(file = "data/error", exist=isfile)
        if (isfile) then
            open (unit = 1, file = "data/error", status ='OLD', POSITION="append")
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
    dest=0
    call MPI_Send(x,nTerms,MPI_doUBLE_PRECISION,dest,0,MPI_COMM_WORLD,ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD,id,ierror)
    ! send ind to head node
    if (id.eq.1) then
        call MPI_Send(md%mc_ind,1,MPI_integer,dest,0,MPI_COMM_WORLD,ierror)
    endif
    ! hear back on which replica and it's mu value
    source=0
    ! get new replica number
    call MPI_Recv(md%rep,1,MPI_integer,source,0, &
                  MPI_COMM_WORLD,status,ierror)
    ! get new mu value
    call MPI_Recv(cof,nTerms,MPI_doUBLE_PRECISION,source,0,&
                  MPI_COMM_WORLD,status,ierror)

    mc%chi      =cof(1)
    mc%mu       =cof(2)
    mc%hA       =cof(3)
    mc%HP1_Bind =cof(4)
    mc%Kap      =cof(5)
    !mc%para(1)  =cof(6)
    !mc%para(2)  =cof(7)
    !mc%para(3)  =cof(8)

    if (abs(md%EChi-x(1)*CofOld(1)).gt.0.0000001_dp) then
        print*, "Error in replicaExchange"
        print*, "mc%EChi",md%EChi,"x(1)*CofOld(1)",x(1)*CofOld(1)
        stop 1
    endif

    md%EChi    =md%EChi    +x(1)*(Cof(1)-CofOld(1))
    md%ebind   =md%ebind   +x(2)*(Cof(2)-CofOld(2))
    md%EField  =md%EField  +x(3)*(Cof(3)-CofOld(3))
    md%ECouple =md%ECouple +x(4)*(Cof(4)-CofOld(4))
    md%EKap    =md%EKap    +x(5)*(Cof(5)-CofOld(5))
   ! mc%EElas(1)=mc%EElas(1)+x(6)*(Cof(6)-CofOld(6))
   ! mc%EElas(2)=mc%EElas(2)+x(7)*(Cof(7)-CofOld(7))
   ! mc%EElas(3)=mc%EElas(3)+x(8)*(Cof(8)-CofOld(8))

    if (abs(md%EChi-x(1)*Cof(1)).gt.0.000001_dp) then
        print*, "Error in replicaExchange"
        print*, "mc%EChi",md%EChi,"x(1)*Cof(1)",x(1)*Cof(1)
        stop 1
    endif

    ! change output file sufix
    write(iostr,"(I4)") md%rep
    iostr=adjustL(iostr)
    iostr=trim(iostr)
    iostr="v"//trim(iostr)
    iostr=trim(iostr)
    mc%repSuffix=iostr


    ! keep track of which thread you are
    md%id=int(id)
end subroutine
