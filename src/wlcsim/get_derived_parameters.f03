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

    if (wlc_p%NB.EQ.1.0d0) then
        PRINT*, 'Non-dimensionalization used requires at least two beads, 1 requested.'
        STOP 1
    endif
    wlc_p%REND=wlc_p%REND*wlc_p%L/wlc_p%LP

    !TODO code smell, ask brad what he needs here
    IF (wlc_p%RING) THEN
        wlc_p%DEL=wlc_p%L/wlc_p%LP/(wlc_p%NB)
    ELSE
        wlc_p%DEL=wlc_p%L/wlc_p%LP/(wlc_p%NB-1.0_dp)
    ENDIF

!     Load the tabulated parameters

    OPEN (UNIT=5,FILE='input/dssWLCparams',STATUS='OLD')
    DO I=1,679
        READ(5,*) PVEC(I,1),PVEC(I,2),PVEC(I,3),PVEC(I,4),PVEC(I,5),PVEC(I,6),PVEC(I,7),PVEC(I,8)
    ENDDO
    CLOSE(5)


!     Setup the parameters for WLC simulation

    if (wlc_p%DEL.LT.PVEC(1,1)) then
        PRINT*, 'It has never been known if the WLC code actually works.'
        PRINT*, 'An entire summer student (Luis Nieves) was thrown at this'
        PRINT*, 'problem and it is still not solved.'
        stop 1
        wlc_p%EB=wlc_p%LP/wlc_p%DEL
        wlc_p%GAM=wlc_p%DEL
        wlc_p%XIR=wlc_p%L/wlc_p%LP/wlc_p%NB
        wlc_p%SIMTYPE=1

!    Setup the parameters for GC simulation

    elseif (wlc_p%DEL.GT.PVEC(679,1)) then
        wlc_p%EPAR=1.5/wlc_p%DEL
        wlc_p%GAM=0.0_dp
        wlc_p%SIMTYPE=3
        wlc_p%XIR=wlc_p%L/wlc_p%NB/wlc_p%LP

!    Setup the parameters for ssWLC simulation

    else !  if (DEL.GE.PVEC(1,1).AND.DEL.LE.PVEC(679,1)) then
        wlc_p%SIMTYPE=2

    ! find(del < pvec, 1, 'first')
    IND=1
    do while (wlc_p%DEL.GT.PVEC(IND,1))
        IND=IND+1
    enddo

    !     Perform linear interpolations
    I=2
    M=(PVEC(IND,I)-PVEC(IND-1,I))/(PVEC(IND,1)-PVEC(IND-1,1))
    wlc_p%EB=M*(wlc_p%DEL-PVEC(IND,1))+PVEC(IND,I)

    I=3
    M=(PVEC(IND,I)-PVEC(IND-1,I))/(PVEC(IND,1)-PVEC(IND-1,1))
    wlc_p%GAM=M*(wlc_p%DEL-PVEC(IND,1))+PVEC(IND,I)

    I=4
    M=(PVEC(IND,I)-PVEC(IND-1,I))/(PVEC(IND,1)-PVEC(IND-1,1))
    wlc_p%EPAR=M*(wlc_p%DEL-PVEC(IND,1))+PVEC(IND,I)

    I=5
    M=(PVEC(IND,I)-PVEC(IND-1,I))/(PVEC(IND,1)-PVEC(IND-1,1))
    wlc_p%EPERP=M*(wlc_p%DEL-PVEC(IND,1))+PVEC(IND,I)

    I=6
    M=(PVEC(IND,I)-PVEC(IND-1,I))/(PVEC(IND,1)-PVEC(IND-1,1))
    wlc_p%ETA=M*(wlc_p%DEL-PVEC(IND,1))+PVEC(IND,I)

    I=7
    M=(PVEC(IND,I)-PVEC(IND-1,I))/(PVEC(IND,1)-PVEC(IND-1,1))
    wlc_p%XIU=M*(wlc_p%DEL-PVEC(IND,1))+PVEC(IND,I)

    wlc_p%EB=wlc_p%LP*wlc_p%EB/(wlc_p%DEL*wlc_p%LP)
    wlc_p%EPAR=wlc_p%EPAR/(wlc_p%DEL*wlc_p%LP*wlc_p%LP)
    wlc_p%EPERP=wlc_p%EPERP/(wlc_p%DEL*wlc_p%LP*wlc_p%LP)
    wlc_p%GAM=wlc_p%DEL*wlc_p%LP*wlc_p%GAM
    wlc_p%ETA=wlc_p%ETA/wlc_p%LP
    wlc_p%XIU=wlc_p%XIU*wlc_p%L/wlc_p%NB/wlc_p%LP
    wlc_p%XIR=wlc_p%L/wlc_p%LP/wlc_p%NB
    wlc_p%DT=0.5*wlc_p%XIU/(wlc_p%EPERP*wlc_p%GAM**2.)

    wlc_p%L0 = wlc_p%GAM
    endif

    return
end subroutine get_derived_parameters

!---------------------------------------------------------------*
