#include "../defines.inc"
! file defines routines useful for parallel tempering continuous variables over
! MPI. Thus, only compiled if MPI is available.
#if MPI_VERSION

subroutine startWorker()
! values from wlcsim_data
use params, only: wlc_id, wlc_AB, wlc_rep, wlc_repSuffix, wlc_METH
use params
use mpi
use energies, only: NUMBER_OF_ENERGY_TYPES, energyOf
! Override initialization with parallel setup parameters energyOf(*_)%cof
    Implicit none
    integer (kind = 4) dest ! message destination
    integer (kind = 4) source ! message source
    integer (kind = 4) id, nThreads
    integer (kind = 4) error  ! error id for MIP functions
    character(MAXFILENAMELEN) iostrg    ! for file naming
    integer ( kind = 4 ) status(MPI_status_SIZE) ! MPI stuff
    real(dp) cof(NUMBER_OF_ENERGY_TYPES)
    integer ii

    call MPI_COMM_SIZE(MPI_COMM_WORLD,nThreads,error)
    call MPI_COMM_RANK(MPI_COMM_WORLD,id,error)
    if (wlc_id .ne. id) then
        print*, "ID mismatch in PT_override!"
        stop
    endif
    if (nThreads.lt.3) then
        wlc_repSuffix = "v1"
        wlc_rep = 1
        !wlc_id = int(id)
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
        if (wlc_id.eq.1) then
            do dest = 2,nThreads-1
                if(WLC_P__VARIABLE_CHEM_STATE) then
                    if (.not. WLC_P__CHEM_SEQ_FROM_FILE) then
                        call MPI_Send (wlc_METH,WLC_P__NT, MPI_integer, dest,   0, &
                                       MPI_COMM_WORLD,error )
                    endif
                else
                    call MPI_Send (wlc_AB,WLC_P__NT, MPI_integer, dest,   0, &
                                   MPI_COMM_WORLD,error )
                endif
            enddo
        elseif (wlc_id.gt.1) then
            source = 1
            if(WLC_P__VARIABLE_CHEM_STATE) then
                if (.not. WLC_P__CHEM_SEQ_FROM_FILE) then
                    call MPI_Recv (wlc_METH, WLC_P__NT, MPI_integer, source, 0, &
                                   MPI_COMM_WORLD, status, error )
                endif
            else
                call MPI_Recv (wlc_AB, WLC_P__NT, MPI_integer, source, 0, &
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
    call MPI_Recv ( wlc_rep, 1, MPI_integer, source, 0, &
      MPI_COMM_WORLD, status, error )

    call MPI_Recv ( cof, NUMBER_OF_ENERGY_TYPES, MPI_doUBLE_PRECISION, source, 0, &
      MPI_COMM_WORLD, status, error )

    ! set cof values
    do ii = 1, NUMBER_OF_ENERGY_TYPES
        energyOf(ii)%cof = cof(ii)
    enddo

    write(iostrg,"(I4)") wlc_rep
    iostrg = adjustL(iostrg)
    iostrg = trim(iostrg)
    iostrg = "v"//trim(iostrg)
    iostrg = trim(iostrg)
    wlc_repSuffix = iostrg

end subroutine
subroutine replicaExchange()
! values from wlcsim_data
use params, only: wlc_mc_ind, wlc_rep,  wlc_repSuffix
! This checks in with the mpi head node to
! For parallel tempering of the form:  E = cof*x
! 1: Tell head node the x value
! 2: Recive replica assignment from head node
! 3: Recive assigned cof value
use params, only : dp, MAXFILENAMELEN, epsApprox
use mpi
use energies, only: NUMBER_OF_ENERGY_TYPES, energyOf
    implicit none
    integer (kind = 4) id, error
    integer (kind = 4) dest ! message destination
    integer (kind = 4) source ! message source
    integer (kind = 4) nThreads
    character(MAXFILENAMELEN) iostr ! for handling sufix string
    integer status(MPI_status_SIZE)  ! MPI status
    real(dp) cof(NUMBER_OF_ENERGY_TYPES)
    real(dp) cof_old(NUMBER_OF_ENERGY_TYPES)
    real(dp) x(NUMBER_OF_ENERGY_TYPES)
    integer ii

    call MPI_COMM_SIZE(MPI_COMM_WORLD,nThreads,error)
    if (nThreads.lt.3) return

    do ii = 1,NUMBER_OF_ENERGY_TYPES
        x(ii) = energyOf(ii)%x  ! package x values for sending
        cof_old(ii) = energyOf(ii)%cof ! save old cof values for below
    enddo

    ! send x to head node
    dest = 0
    call MPI_Send(x,NUMBER_OF_ENERGY_TYPES,MPI_doUBLE_PRECISION,dest,0,MPI_COMM_WORLD,error)
    call MPI_COMM_RANK(MPI_COMM_WORLD,id,error)
    ! send ind to head node
    if (id.eq.1) then
        call MPI_Send(wlc_mc_ind,1,MPI_integer,dest,0,MPI_COMM_WORLD,error)
    endif
    ! hear back on which replica and it's cof value
    source = 0
    ! get new replica number
    call MPI_Recv(wlc_rep,1,MPI_integer,source,0, &
                  MPI_COMM_WORLD,status,error)
    ! get new cof value
    call MPI_Recv(cof,NUMBER_OF_ENERGY_TYPES,MPI_doUBLE_PRECISION,source,0,&
                  MPI_COMM_WORLD,status,error)

    do ii = 1,NUMBER_OF_ENERGY_TYPES
        ! update cof
        energyOf(ii)%cof = cof(ii)
        ! Verify that the old E=cof*x
        if (abs(energyOf(ii)%E-energyOf(ii)%x*cof_old(ii)).gt.epsApprox) then
            print*, "Error in replicaExchange of ", energyOf(ii)%name_str, " before"
            print*, "E=",energyOf(ii)%E," != x*CofOld",energyOf(ii)%x*cof_old(ii)
            stop 1
        endif
        ! Update energy useing x and change in cof
        energyOf(ii)%E=energyOf(ii)%E+energyOf(ii)%x*(energyOf(ii)%cof-cof_old(ii))
        ! Verify that the new E=cof*x
        if (abs(energyOf(ii)%E-energyOf(ii)%cof*energyOf(ii)%x).gt.epsApprox) then
            print*, "Error in replicaExchange of ", energyOf(ii)%name_str, " after"
            print*, "E=",energyOf(ii)%E," but  cof*x)=",energyOf(ii)%cof*energyOf(ii)%x
            print*, "cof",energyOf(ii)%cof,"- cof_old",cof_old(ii),"=",energyOf(ii)%cof-cof_old(ii)
            print*, "x",energyOf(ii)%x
            stop 1
        endif
    enddo
    ! change output file sufix
    write(iostr,"(I4)") wlc_rep
    iostr = adjustL(iostr)
    iostr = trim(iostr)
    iostr = "v"//trim(iostr)
    iostr = trim(iostr)
    wlc_repSuffix = iostr

end subroutine

#endif
