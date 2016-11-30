program main

!*****************************************************************************
!
!  This program is a MTMC with parallel tempering.
!  The MC code is in wlcsim
!
!
!  Modified:   5/6/2016
!
!  Author:  Quinn MacPherson
!
!*****************************************************************************
  use mpi
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
!
!  Shut down MPI.
!
  call MPI_Finalize ( error )
!
!  Terminate.
!
  if ( id == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Normal end of execution.'
    write ( *, '(a)' ) ' '
  end if

  stop
end
subroutine singleCall()
    use mersenne_twister
    implicit none
!   variable for random number generator seeding
    type(random_stat) rand_stat  ! state of random number chain
    integer Irand     ! Seed
    character*8 datedum  ! trash
    character*10 timedum ! trash
    character*5 zonedum  ! trash
    integer seedvalues(8) ! clock readings
    !real urand(1) ! use this to generate a random number
    if (.false.) then ! set spedific seed
        Irand=7171
    else ! seed from clock
        call date_and_time(datedum,timedum,zonedum,seedvalues)
        Irand=int(-seedvalues(5)*1E7-seedvalues(6)*1E5&
                  -seedvalues(7)*1E3-seedvalues(8))
        Irand=mod(Irand,10000)
        print*, "Random Intiger seed:",Irand
    endif
    call random_setseed(Irand,rand_stat) ! random seed for head node
    print*, "calling single wlcsim"
    call wlcsim(rand_stat)
   !call simpleSim(rand_stat)
end subroutine
subroutine paraTemp ( p, id)

!*****************************************************************************
!
!
  use mersenne_twister
    implicit none
!   MPI variables
    integer ( kind = 4 ) dest   !destination id for messages
    integer ( kind = 4 ) source  !source id for messages
    integer ( kind = 4 ) error  ! error id for MIP functions
    integer ( kind = 4 ), intent(in) :: id     ! which processor I am
    integer ( kind = 4 ), intent(in) :: p ! number of threads
    integer ( kind = 4 ) status(MPI_status_SIZE) ! MPI stuff
    integer nPTReplicas     ! number of replicas.

!   variable for random number generator seeding
    type(random_stat) rand_stat  ! state of random number chain
    integer Irand     ! Seed
    character*8 datedum  ! trash
    character*10 timedum ! trash
    character*5 zonedum  ! trash
    integer seedvalues(8) ! clock readings
    real urand(1)

!   worker node only variables
    real(dp) nan_dp !chemical potential

!   for head node use only variables
    integer rep ! physical replica number, for loops
    integer temp ! for castling
    logical keepGoing   ! set to false when NaN encountered
    integer, allocatable :: nodeNumber(:)  ! list of which nodes are which
    real(dp), allocatable :: xMtrx(:,:)  ! sum of bound states
    real(dp), allocatable :: cofMtrx(:,:) ! mu or chi or whatever
    real(dp), allocatable :: s_vals(:) ! path parameter
    integer, parameter :: nTerms=8  ! number of energy terms
    real(dp) x(nTerms) ! slice of xMtrx
    real(dp) cof(nTerms) ! slice of cofMtrx
    integer N_average      ! number of attempts since last average
    integer upSuccess(p-1)  ! number of successes since last average
    integer downSuccess(p-1) ! number of successes since last average
    integer nExchange ! total number of exchanges attemted
    type(wlcsim_params) mc ! genaral symulation parameters
    character*16 fileName ! ouput filename
    real(dp) energy ! for deciding to accept exchange
    integer term ! for loopin over terms
    real(dp) h_path,chi_path,mu_path,kap_path,HP1_Bind_path ! functions

    nExchange=0
    nPTReplicas=p-1
    fileName='input/params'
!  -----------------------------------------------------------------
!
!      Head Node
!
!  -----------------------------------------------------------------
    if ( id == 0 ) then
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

        call wlcsim_params_setparams(mc,fileName) ! so that the head thread knows the  parameters

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
                cofMtrx(rep,3)=mc%h_A
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
            cofMtrx(rep,6)=mc%para(1)
            cofMtrx(rep,7)=mc%para(2)
            cofMtrx(rep,8)=mc%para(3)
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
        mc%ind=0
        keepGoing=.True.
        do while(keepGoing)
            ! give workers thier jobs
            do rep=1,nPTReplicas
                dest=nodeNumber(rep)
                call MPI_Send (rep,1, MPI_integer, dest,   0, &
                                MPI_COMM_WORLD,error )
                if (mc%restart.and.mc%ind.eq.0) then
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
            call MPI_Recv (mc%ind, 1, MPI_integer, source, 0, &
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
                                     nExchange,mc%ind,nTerms,s_vals)

                if ((mc%ind.ge.mc%indStartRepAdapt).and. &
                    (mc%ind.lt.mc%indendRepAdapt)) then ! insert input defined location here
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
    else
!  -----------------------------------------------------------------
!
!      Worker Node
!
!  -----------------------------------------------------------------

        source = 0
        dest = 0
        ! -----------------------------------------------
        !
        !   Generate thread safe random number chain: rand_stat
        !
        !--------------------------------------------

        call MPI_Recv ( Irand, 1, MPI_integer, source, 0, &
                        MPI_COMM_WORLD, status, error )
        call random_setseed(Irand*(id+1),rand_stat) ! random seed for head node
        ! ------------------------------
        !
        ! call main simulation code
        !
        !  --------------------------------
        call wlcsim(rand_stat)
        !call simpleSim(rand_stat)
        nan_dp=0; nan_dp=nan_dp/nan_dp !NaN
        x(1)=nan_dp
        print*, "Node ",id," sending normal exit code."
        call MPI_Send(x,nTerms,MPI_doUBLE_PRECISION,0,0,MPI_COMM_WORLD,error)
    end if


    print*, "Node ",id," exiting normally."
    return
end
function chi_path(s) result(chi)
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
    implicit none
    real(dp), intent(in) :: s
    real(dp) mu
    mu=s
end function mu_path
function kap_path(s) result(kap)
    implicit none
    real(dp), intent(in) :: s
    real(dp) kap
    kap=s*10.0_dp
end function kap_path
function hp1_bind_path(s) result(hp1_bind)
    implicit none
    real(dp), intent(in) :: s
    real(dp) hp1_bind
    hp1_bind=s
end function hp1_bind_path
subroutine PT_override(mc,md)
! Override initialization with parallel setup parameters
!  In particualar it changes: mc%AB, mc%rep, mc%mu, mc%repSufix
    Implicit none
    type(wlcsim_params), intent(inout) :: mc
    type(wlcsim_data), intent(inout) :: md
    integer (kind=4) dest ! message destination
    integer (kind=4) source ! message source
    integer (kind=4) id, nThreads,ierror
    integer (kind=4) error  ! error id for MIP functions
    character*16 iostrg    ! for file naming
    integer ( kind = 4 ) status(MPI_status_SIZE) ! MPI stuff
    integer, parameter :: nTerms=8  ! number of energy terms
    real(dp) cof(nTerms)

    call MPI_COMM_SIZE(MPI_COMM_WORLD,nThreads,ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD,id,ierror)
    if (nThreads.lt.3) then
        mc%repSufix="v1"
        mc%rep=1
        mc%id=int(id)
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
    call MPI_Recv ( mc%rep, 1, MPI_integer, source, 0, &
      MPI_COMM_WORLD, status, error )

    call MPI_Recv ( cof, nTerms, MPI_doUBLE_PRECISION, source, 0, &
      MPI_COMM_WORLD, status, error )

    mc%chi      =cof(1)
    mc%mu       =cof(2)
    mc%h_A      =cof(3)
    mc%HP1_Bind =cof(4)
    mc%Kap      =cof(5)
    !mc%para(1)  =cof(6)
    !mc%para(2)  =cof(7)
    !mc%para(3)  =cof(8)

    write(iostrg,"(I4)"), mc%rep
    iostrg=adjustL(iostrg)
    iostrg=trim(iostrg)
    iostrg="v"//trim(iostrg)
    iostrg=trim(iostrg)
    mc%repSufix=iostrg

    ! keep track of which thread you are
    mc%id=int(id)
end subroutine
subroutine replicaExchange(mc)
! This checks in with the mpi head node to
! For parallel tempering of the form:  E=cof*x
! 1: Tell head node the x value
! 2: Recive replica assignment from head node
! 3: Recive assigned cof value
    IMPLICIT NONE
    integer, parameter :: nTerms=8  ! number of energy terms
    integer (kind=4) id, ierror
    type(wlcsim_params), intent(inout) :: mc
    integer i  ! working intiger
    integer (kind=4) dest ! message destination
    integer (kind=4) source ! message source
    integer (kind=4) nThreads
    character*16 iostr ! for handling sufix string
    integer status(MPI_status_SIZE)  ! MPI status
    real(dp) cof(nTerms)
    real(dp) cofOld(nTerms)
    real(dp) x(nTerms)
    real(dp) test(5)
    logical isfile

    call MPI_COMM_SIZE(MPI_COMM_WORLD,nThreads,ierror)
    if (nThreads.lt.3) return

    x(1)=mc%x_Chi
    x(2)=mc%x_mu
    x(3)=mc%x_Field
    x(4)=mc%x_couple
    x(5)=mc%x_kap
    x(6)=0.0_dp !x(6)=mc%EElas(1)/mc%para(1)
    x(7)=0.0_dp !x(7)=mc%EElas(2)/mc%para(2)
    x(8)=0.0_dp !x(8)=mc%EElas(3)/mc%para(3)

    test(1)=mc%EChi/mc%Chi
    test(3)=mc%EField/mc%h_A
    test(4)=mc%ECouple/mc%HP1_Bind
    test(5)=mc%EKap/mc%Kap

    cofOld(1)=mc%chi
    cofOld(2)=mc%mu
    cofOld(3)=mc%h_A
    cofOld(4)=mc%HP1_Bind
    cofOld(5)=mc%Kap
    cofOld(6)=mc%para(1)
    cofOld(7)=mc%para(2)
    cofOld(8)=mc%para(3)

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
        write(1,*), "Error in replicaExchange"
        write(1,*), "I",I," test",test(I)," x",x(I)," cof",cofOld(I)
        print*, "Error in replicaExchange"
        print*, "I",I," test",test(I)," x",x(I)," cof",cofOld(I)
        close (1)
    enddo

    ! send number bound to head node
    dest=0
    call MPI_Send(x,nTerms,MPI_doUBLE_PRECISION,dest,0,MPI_COMM_WORLD,mc%error)
    call MPI_COMM_RANK(MPI_COMM_WORLD,id,ierror)
    ! send ind to head node
    if (id.eq.1) then
        call MPI_Send(mc%ind,1,MPI_integer,dest,0,MPI_COMM_WORLD,mc%error)
    endif
    ! hear back on which replica and it's mu value
    source=0
    ! get new replica number
    call MPI_Recv(mc%rep,1,MPI_integer,source,0, &
                  MPI_COMM_WORLD,status,mc%error)
    ! get new mu value
    call MPI_Recv(cof,nTerms,MPI_doUBLE_PRECISION,source,0,&
                  MPI_COMM_WORLD,status,mc%error)

    mc%chi      =cof(1)
    mc%mu       =cof(2)
    mc%h_A      =cof(3)
    mc%HP1_Bind =cof(4)
    mc%Kap      =cof(5)
    !mc%para(1)  =cof(6)
    !mc%para(2)  =cof(7)
    !mc%para(3)  =cof(8)

    if (abs(mc%EChi-x(1)*CofOld(1)).gt.0.0000001_dp) then
        print*, "Error in replicaExchange"
        print*, "mc%EChi",mc%EChi,"x(1)*CofOld(1)",x(1)*CofOld(1)
        stop 1
    endif

    mc%EChi    =mc%EChi    +x(1)*(Cof(1)-CofOld(1))
    mc%ebind   =mc%ebind   +x(2)*(Cof(2)-CofOld(2))
    mc%EField  =mc%EField  +x(3)*(Cof(3)-CofOld(3))
    mc%ECouple =mc%ECouple +x(4)*(Cof(4)-CofOld(4))
    mc%EKap    =mc%EKap    +x(5)*(Cof(5)-CofOld(5))
   ! mc%EElas(1)=mc%EElas(1)+x(6)*(Cof(6)-CofOld(6))
   ! mc%EElas(2)=mc%EElas(2)+x(7)*(Cof(7)-CofOld(7))
   ! mc%EElas(3)=mc%EElas(3)+x(8)*(Cof(8)-CofOld(8))

    if (abs(mc%EChi-x(1)*Cof(1)).gt.0.000001_dp) then
        print*, "Error in replicaExchange"
        print*, "mc%EChi",mc%EChi,"x(1)*Cof(1)",x(1)*Cof(1)
        stop 1
    endif

    ! change output file sufix
    write(iostr,"(I4)"), mc%rep
    iostr=adjustL(iostr)
    iostr=trim(iostr)
    iostr="v"//trim(iostr)
    iostr=trim(iostr)
    mc%repSufix=iostr


    ! keep track of which thread you are
    mc%id=int(id)
end subroutine
