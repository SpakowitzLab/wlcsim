#include "../defines.inc"
#if MPI_VERSION
subroutine pt_restart(wlc_p)
! values from wlcsim_data
use params, only: wlc_WindoW, wlc_PHIT, wlc_repSuffix,  wlc_rep&
    , wlc_MCAMP&
    , wlc_U, wlc_R, wlc_id&
    , wlc_mc_ind, wlc_AB 
use energies
! Takes wlcsim_params and wlcsim_data and restarts the MPI workers for running
! parallel-tempered MC simulations.
!
! This function takes the place of PT_override in the case of restart
! This will read from a output directory and restart multiple replicas
! Override initialization with parallel setup parameters
!  In particualar it changes: wlc_p%AB, wlc_p%rep, energyOf(mu_)%cof, wlc_p%repSuffix
    use mpi
    use params
    Implicit none
    type(wlcsim_params), intent(inout) :: wlc_p
    integer (kind = 4) dest ! message destination
    integer (kind = 4) source ! message source
    integer (kind = 4) id, nThreads,ierror
    integer (kind = 4) error  ! error id for MIP functions
    character(MAXFILENAMELEN) iostrg    ! for file naming
    character(16) vNum    ! for file naming
    character(MAXFILENAMELEN) dir
    integer ( kind = 4 ) status(MPI_status_SIZE) ! MPI stuff
    integer, parameter :: nTerms = 8  ! number of energy terms
    real(dp) mag ! magnitude for renormalizing U
    real(dp) cof(nTerms)
    integer I ! bead index

    ! file parsing
    integer ios ! read status (detect end of file)
    real(dp) temp(28)  ! values

    ! Which replica am I?
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nThreads,ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD,id,ierror)
    source = 0;
    call MPI_Recv ( wlc_rep, 1, MPI_integer, source, 0, &
                   MPI_COMM_WORLD, status, error )
    if (wlc_rep.ne.id) then
        print*, "That's not what I expected! see restart"
    endif
    if (nThreads.lt.3) then
        print*, "don't use pt_restart for fewer than 3 treads"
        stop 1
    endif
    write(vNum,'(I4)') wlc_rep
    vNum = adJustL(vNum)
    vNum = "v"//trim(vNum)

    ! Where to read from
    dir = "data/"

    ! read Some operation variables
    ! Many of these aren't necessary but a few are
    iostrg = trim(dir)//"out1"
    iostrg = trim(iostrg)//trim(vNum)
    print*, "reading", iostrg
    open(unit = 1, file = iostrg, status ='OLD')
    read(1,*)
    do while (.TRUE.)
        read(1,*,IOSTAT = ios) temp(1), temp(2), temp(3), temp(4), temp(5), &
                             temp(6), temp(7), temp(8), temp(9), temp(10), &
                             temp(11), temp(12), temp(13), temp(14), &
                             temp(15), temp(16)
        if (ios.eq.0) then
            wlc_mc_ind = nint(temp(1))
            energyOf(bend_)%E = temp(3)
            energyOf(stretch_)%E = temp(4)
            energyOf(shear_)%E = temp(5)
            energyOf(couple_)%E = temp(6)
            energyOf(kap_)%E = temp(7)
            energyOf(chi_)%E = temp(8)
            energyOf(field_)%E = temp(9)
            energyOf(bind_)%E = temp(10)
            energyOf(mu_)%x = temp(11)
            energyOf(couple_)%cof = temp(12)
            energyOf(chi_)%cof = temp(13)
            energyOf(mu_)%cof = temp(14)
            energyOf(kap_)%cof = temp(15)
            energyOf(field_)%cof = temp(16)
            ! x_ms
            ! chi_l2
            ! E_mu
        else
            Exit
        endif
    enddo
    close(1)
    print*, "first set from file", iostrg
    print*, temp
    ! not sure if the following if statments are necessary
    if (energyOf(chi_)%cof.ne.0.0) then
        energyOf(chi_)%x = energyOf(chi_)%E/energyOf(chi_)%cof
    endif
    if (energyOf(chi_)%cof.ne.0.0) then
        energyOf(couple_)%x = energyOf(couple_)%E/energyOf(couple_)%cof
    endif
    if (energyOf(kap_)%cof.ne.0) then
        energyOf(kap_)%x = energyOf(kap_)%E/energyOf(kap_)%cof
    endif
    if (energyOf(field_)%x.ne.0.0) then
        energyOf(field_)%x = energyOf(field_)%E/energyOf(field_)%cof
    endif
    if (energyOf(mu_)%cof.ne.0.0) then
        energyOf(mu_)%x = energyOf(mu_)%E/energyOf(mu_)%cof
    endif

    ! read back in addaptation stuff, May make slight difference
    iostrg = trim(dir)//"out3"
    iostrg = trim(iostrg)//trim(vNum)
    print*, iostrg
    open(unit = 1, file = iostrg, status ='OLD')
    read(1,*)
    do while (.TRUE.)
        read(1,*,IOSTAT = ios)  temp(1), temp(2), temp(3), temp(4), temp(5), &
                              temp(6), temp(7), temp(8), temp(9), temp(10), &
                              temp(11), temp(12), temp(13), temp(14), &
                              temp(15), temp(16), temp(17), temp(18), &
                              temp(19), temp(20), temp(21), temp(22), &
                              temp(23), temp(24), temp(25), temp(26), &
                              temp(27), temp(28)
        if (ios.eq.0) then
            wlc_WindoW(1) = temp(3); wlc_MCAMP(1) = temp(4); wlc_PHIT(1) = temp(5);
            wlc_WindoW(2) = temp(6); wlc_MCAMP(2) = temp(7); wlc_PHIT(2) = temp(8);
            wlc_WindoW(3) = temp(9); wlc_MCAMP(3) = temp(10); wlc_PHIT(3) = temp(11);
            wlc_p%MOVEON(4) = nint(temp(12)); wlc_MCAMP(4) = temp(13); wlc_PHIT(4) = temp(14);
            wlc_p%MOVEON(5) = nint(temp(15)); wlc_MCAMP(5) = temp(16); wlc_PHIT(5) = temp(17);
            wlc_p%MOVEON(6) = nint(temp(18)); wlc_MCAMP(6) = temp(19); wlc_PHIT(6) = temp(20);
            wlc_p%MOVEON(7) = nint(temp(21)); wlc_PHIT(7) = temp(22);
            wlc_p%MOVEON(8) = nint(temp(23)); wlc_PHIT(8) = temp(24);
            wlc_p%MOVEON(9) = nint(temp(25)); wlc_PHIT(9) = temp(26);
            wlc_p%MOVEON(10) = nint(temp(27)); wlc_PHIT(10) = temp(28)
        else
            Exit
        endif
    enddo
    close(1)
    print*, "second set from file", iostrg
    print*, temp



    ! read R and AB from file
    write(iostrg,"(I8)") wlc_mc_ind
    iostrg = adjustL(iostrg)
    iostrg = "r"//trim(iostrg)
    iostrg = trim(dir)//trim(iostrg)
    iostrg = trim(iostrg)//trim(vNum)
    print*, "reading", iostrg
    open (unit = 5, file = iostrg, status = 'OLD')
    print*, "NT = ",WLC_P__NT
    ios = 0;
    do I = 1,WLC_P__NT
       if (ios.ne.0) then
           print*, "Problem while reading R, Possible incomplete file"
           stop 1
       endif
       read(5,*) wlc_R(1,I),wlc_R(2,I),wlc_R(3,I),wlc_AB(I)
    enddo
    close(5)

    ! read U
    write(iostrg,"(I8)") wlc_mc_ind
    iostrg = adjustL(iostrg)
    iostrg = "u"//trim(iostrg)
    iostrg = trim(dir)//trim(iostrg)
    iostrg = trim(iostrg)//trim(vNum)
    ! read U from file
    open (unit = 5, file = iostrg, status = 'OLD')
    do I = 1,WLC_P__NT
       read(5,*) wlc_U(1,I),wlc_U(2,I),wlc_U(3,I)
       mag = sqrt(wlc_U(1,I)**2 + wlc_U(2,I)**2 + wlc_U(3,I)**2)
       wlc_U(1,I) = wlc_U(1,I)/mag
       wlc_U(2,I) = wlc_U(2,I)/mag
       wlc_U(3,I) = wlc_U(3,I)/mag
    enddo
    close(5)

    ! Let head node know what cof values you read
    cof(1) = energyOf(chi_)%cof
    cof(2) = energyOf(mu_)%cof
    cof(3) = energyOf(field_)%cof
    cof(4) = energyOf(couple_)%cof
    cof(5) = energyOf(kap_)%cof
    cof(6) = 0
    cof(7) = 0
    cof(8) = 0
    dest = 0
    call MPI_Send (cof,nTerms, MPI_doUBLE_PRECISION, dest,   0, &
                    MPI_COMM_WORLD,error )

    ! Make repsuffix
    write(iostrg,"(I4)") wlc_rep
    iostrg = adjustL(iostrg)
    iostrg = trim(iostrg)
    iostrg = "v"//trim(iostrg)
    iostrg = trim(iostrg)
    wlc_repSuffix = trim(iostrg)

    ! keep track of which thread you are
    wlc_id = int(id)
end subroutine
#endif
