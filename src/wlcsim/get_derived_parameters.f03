! *---------------------------------------------------------------*
!
!     subroutine getpara.f95
!     Setup the parameters for the simulation
!
!     1. Determine the simulation type
!     2. Evaluate the polymer elastic parameters
!     3. Determine the parameters for Brownian dynamics simulation
!
!     Andrew Spakowitz
!     8/17/15
!

subroutine get_derived_parameters(wlc_p)

    use params

    implicit none

    integer i,ind
    real(dp) m

    type(wlcsim_params), intent(inout) :: wlc_p
    REAL(dp) :: pvec(679, 8) ! array holding dssWLC params calculated by Elena

    !Calculate total number of beads
    wlc_p%nT = wlc_p%nB*wlc_p%nP

    if (wlc_p%NB == 1.0d0) then
        ! since we use "DEL" as an intermediate, we need at least two beads
        PRinT*, 'Some intermediate calculations used require at least two beads, 1 requested.'
        STOP 1
    endif

    ! calculate metrics that don't change between WLC, ssWLC, GC
    if (wlc_p%RinG) then
        wlc_p%DEL = wlc_p%L/wlc_p%LP/(wlc_p%NB)
    else
        wlc_p%DEL = wlc_p%L/wlc_p%LP/(wlc_p%NB-1.0_dp)
    ENDif
    ! std dev of interbead distribution of GC, used to initialize
    wlc_p%SIGMA = sqrt(2.0_dp*wlc_p%LP*wlc_p%L/3.0_dp/wlc_p%NB)

!     Load the tabulated parameters

    open (UNIT = 5,FILE = 'input/dssWLCparams',STATUS = 'OLD')
    do I = 1,679
        READ(5,*) PVEC(I,1),PVEC(I,2),PVEC(I,3),PVEC(I,4),PVEC(I,5),PVEC(I,6),PVEC(I,7),PVEC(I,8)
    ENDdo
    CLOSE(5)


!     Setup the parameters for WLC simulation

    ! if del < 0.01
    if (wlc_p%DEL < PVEC(1,1)) then
        PRinT*, 'It has never been known if the WLC code actually works.'
        PRinT*, 'An entire summer student (Luis Nieves) was thrown at this'
        PRinT*, 'problem and it is still not solved.'
        stop 1
        wlc_p%EB = wlc_p%LP/wlc_p%DEL
        wlc_p%GAM = wlc_p%DEL
        wlc_p%XIR = wlc_p%L/wlc_p%LP/wlc_p%NB
        wlc_p%SIMTYPE = 1

!    Setup the parameters for GC simulation

    ! if del > 10
    elseif (wlc_p%DEL > PVEC(679,1)) then
        wlc_p%EPAR = 1.5/wlc_p%DEL
        wlc_p%GAM = 0.0_dp
        wlc_p%SIMTYPE = 3
        wlc_p%XIR = wlc_p%L/wlc_p%NB/wlc_p%LP

!    Setup the parameters for ssWLC simulation
    ! if 0.01 <= del <= 10
    else !  if (DEL >= PVEC(1,1).AND.DEL <= PVEC(679,1)) then
        wlc_p%SIMTYPE = 2

    ! find(del < pvec, 1, 'first')
    inD = 1
    do while (wlc_p%DEL > PVEC(inD,1))
        inD = inD + 1
    enddo

    !     Perform linear interpolations
    I = 2
    M = (PVEC(inD,I)-PVEC(inD-1,I))/(PVEC(inD,1)-PVEC(inD-1,1))
    wlc_p%EB = M*(wlc_p%DEL-PVEC(inD,1)) + PVEC(inD,I)

    I = 3
    M = (PVEC(inD,I)-PVEC(inD-1,I))/(PVEC(inD,1)-PVEC(inD-1,1))
    wlc_p%GAM = M*(wlc_p%DEL-PVEC(inD,1)) + PVEC(inD,I)

    I = 4
    M = (PVEC(inD,I)-PVEC(inD-1,I))/(PVEC(inD,1)-PVEC(inD-1,1))
    wlc_p%EPAR = M*(wlc_p%DEL-PVEC(inD,1)) + PVEC(inD,I)

    I = 5
    M = (PVEC(inD,I)-PVEC(inD-1,I))/(PVEC(inD,1)-PVEC(inD-1,1))
    wlc_p%EPERP = M*(wlc_p%DEL-PVEC(inD,1)) + PVEC(inD,I)

    I = 6
    M = (PVEC(inD,I)-PVEC(inD-1,I))/(PVEC(inD,1)-PVEC(inD-1,1))
    wlc_p%ETA = M*(wlc_p%DEL-PVEC(inD,1)) + PVEC(inD,I)

    I = 7
    M = (PVEC(inD,I)-PVEC(inD-1,I))/(PVEC(inD,1)-PVEC(inD-1,1))
    wlc_p%XIU = M*(wlc_p%DEL-PVEC(inD,1)) + PVEC(inD,I)

    ! The values read in from file are all non-dimentionalized by the
    ! persistance length.  We now re-dimentionalize them.
    ! We also divied by DEL which is also re-dimentionalized.

    wlc_p%EB = wlc_p%LP*wlc_p%EB/(wlc_p%DEL*wlc_p%LP)
    wlc_p%EPAR = wlc_p%EPAR/(wlc_p%DEL*wlc_p%LP*wlc_p%LP)
    wlc_p%EPERP = wlc_p%EPERP/(wlc_p%DEL*wlc_p%LP*wlc_p%LP)
    wlc_p%GAM = wlc_p%DEL*wlc_p%LP*wlc_p%GAM
    wlc_p%ETA = wlc_p%ETA/wlc_p%LP
    wlc_p%XIU = wlc_p%XIU*wlc_p%L/wlc_p%NB/wlc_p%LP
    wlc_p%XIR = wlc_p%L/wlc_p%LP/wlc_p%NB
    wlc_p%DT = 0.5*wlc_p%XIU/(wlc_p%EPERP*wlc_p%GAM**2.)

    ! wlc_p%L0 = wlc_p%GAM  ! not sure why this was included
    endif

    return
end subroutine get_derived_parameters

!---------------------------------------------------------------*
