#if MPI_VERSION
subroutine pt_restart(wlc_p,wlc_d)
! Takes wlcsim_params and wlcsim_data and restarts the MPI workers for running
! parallel-tempered MC simulations.
!
! This function takes the place of PT_override in the case of restart
! This will read from a output directory and restart multiple replicas
! Override initialization with parallel setup parameters
!  In particualar it changes: wlc_p%AB, wlc_p%rep, wlc_p%mu, wlc_p%repSuffix
    use mpi
    use params
    Implicit none
    type(wlcsim_params), intent(inout) :: wlc_p
    type(wlcsim_data), intent(inout) :: wlc_d
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
    call MPI_Recv ( wlc_d%rep, 1, MPI_integer, source, 0, &
                   MPI_COMM_WORLD, status, error )
    if (wlc_d%rep.ne.id) then
        print*, "That's not what I expected! see restart"
    endif
    if (nThreads.lt.3) then
        print*, "don't use pt_restart for fewer than 3 treads"
        stop 1
    endif
    write(vNum,'(I4)') wlc_d%rep
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
            wlc_d%mc_ind = nint(temp(1))
            wlc_d%EElas(1) = temp(3)
            wlc_d%EElas(2) = temp(4)
            wlc_d%EElas(3) = temp(5)
            wlc_d%ECouple = temp(6)
            wlc_d%EKap = temp(7)
            wlc_d%EChi = temp(8)
            wlc_d%EField = temp(9)
            wlc_d%ebind = temp(10)
            wlc_d%x_Mu = temp(11)
            wlc_p%HP1_Bind = temp(12)
            wlc_p%chi = temp(13)
            wlc_p%mu = temp(14)
            wlc_p%Kap = temp(15)
            wlc_p%hA = temp(16)
        else
            Exit
        endif
    enddo
    close(1)
    print*, "first set from file", iostrg
    print*, temp
    ! not sure if the following if statments are necessary
    if (wlc_p%Chi.ne.0.0) then
        wlc_d%x_Chi = wlc_d%EChi/wlc_p%Chi
    endif
    if (wlc_p%Chi.ne.0.0) then
        wlc_d%x_Couple = wlc_d%ECouple/wlc_p%HP1_Bind
    endif
    if (wlc_p%Kap.ne.0) then
        wlc_d%x_Kap = wlc_d%EKap/wlc_p%Kap
    endif
    if (wlc_d%x_Field.ne.0.0) then
        wlc_d%x_Field = wlc_d%EField/wlc_p%hA
    endif
    if (wlc_p%Mu.ne.0.0) then
        wlc_d%x_Mu = wlc_d%ebind/wlc_p%Mu
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
            wlc_d%WindoW(1) = temp(3); wlc_d%MCAMP(1) = temp(4); wlc_d%PHIT(1) = temp(5);
            wlc_d%WindoW(2) = temp(6); wlc_d%MCAMP(2) = temp(7); wlc_d%PHIT(2) = temp(8);
            wlc_d%WindoW(3) = temp(9); wlc_d%MCAMP(3) = temp(10); wlc_d%PHIT(3) = temp(11);
            wlc_p%MOVEON(4) = nint(temp(12)); wlc_d%MCAMP(4) = temp(13); wlc_d%PHIT(4) = temp(14);
            wlc_p%MOVEON(5) = nint(temp(15)); wlc_d%MCAMP(5) = temp(16); wlc_d%PHIT(5) = temp(17);
            wlc_p%MOVEON(6) = nint(temp(18)); wlc_d%MCAMP(6) = temp(19); wlc_d%PHIT(6) = temp(20);
            wlc_p%MOVEON(7) = nint(temp(21)); wlc_d%PHIT(7) = temp(22);
            wlc_p%MOVEON(8) = nint(temp(23)); wlc_d%PHIT(8) = temp(24);
            wlc_p%MOVEON(9) = nint(temp(25)); wlc_d%PHIT(9) = temp(26);
            wlc_p%MOVEON(10) = nint(temp(27)); wlc_d%PHIT(10) = temp(28)
        else
            Exit
        endif
    enddo
    close(1)
    print*, "second set from file", iostrg
    print*, temp



    ! read R and AB from file
    write(iostrg,"(I8)") wlc_d%mc_ind
    iostrg = adjustL(iostrg)
    iostrg = "r"//trim(iostrg)
    iostrg = trim(dir)//trim(iostrg)
    iostrg = trim(iostrg)//trim(vNum)
    print*, "reading", iostrg
    open (unit = 5, file = iostrg, status = 'OLD')
    print*, "NT = ",wlc_p%NT
    ios = 0;
    do I = 1,wlc_p%NT
       if (ios.ne.0) then
           print*, "Problem while reading R, Possible incomplete file"
           stop 1
       endif
       read(5,*) wlc_d%R(I,1),wlc_d%R(I,2),wlc_d%R(I,3),wlc_d%AB(I)
    enddo
    close(5)

    ! read U
    write(iostrg,"(I8)") wlc_d%mc_ind
    iostrg = adjustL(iostrg)
    iostrg = "u"//trim(iostrg)
    iostrg = trim(dir)//trim(iostrg)
    iostrg = trim(iostrg)//trim(vNum)
    ! read U from file
    open (unit = 5, file = iostrg, status = 'OLD')
    do I = 1,wlc_p%NT
       read(5,*) wlc_d%U(I,1),wlc_d%U(I,2),wlc_d%U(I,3)
       mag = sqrt(wlc_d%U(I,1)**2 + wlc_d%U(I,2)**2 + wlc_d%U(I,3)**2)
       wlc_d%U(I,1) = wlc_d%U(I,1)/mag
       wlc_d%U(I,2) = wlc_d%U(I,2)/mag
       wlc_d%U(I,3) = wlc_d%U(I,3)/mag
    enddo
    close(5)

    ! Let head node know what cof values you read
    cof(1) = wlc_p%chi
    cof(2) = wlc_p%mu
    cof(3) = wlc_p%hA
    cof(4) = wlc_p%HP1_Bind
    cof(5) = wlc_p%KAP
    cof(6) = 0
    cof(7) = 0
    cof(8) = 0
    dest = 0
    call MPI_Send (cof,nTerms, MPI_doUBLE_PRECISION, dest,   0, &
                    MPI_COMM_WORLD,error )

    ! Make repsuffix
    write(iostrg,"(I4)") wlc_d%rep
    iostrg = adjustL(iostrg)
    iostrg = trim(iostrg)
    iostrg = "v"//trim(iostrg)
    iostrg = trim(iostrg)
    wlc_d%repSuffix = trim(iostrg)

    ! keep track of which thread you are
    wlc_d%id = int(id)
end subroutine
#endif
