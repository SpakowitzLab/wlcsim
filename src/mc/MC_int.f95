!---------------------------------------------------------------!

!
!     This subroutine calculates the change in the self energy for
!     a small Monte Carlo move in the position.
!
!     Corrections to force magnitude made 6-3-04.
!
!     Andrew Spakowitz
!     Written 6-29-04
!
!     Edited by Shifan
!
!     Edited by Quinn in 2016

SUBROUTINE MC_int(mc,md,I1,I2,initialize)
use params, only: dp, wlcsim_params, wlcsim_data
IMPLICIT NONE

TYPE(wlcsim_params), intent(in) :: mc   ! <---- Contains output
TYPE(wlcsim_data), intent(inout) :: md
LOGICAL, intent(in) :: initialize   ! if true, calculate absolute energy
INTEGER, intent(in) :: I1           ! Test bead position 1
INTEGER, intent(in) :: I2           ! Test bead position 2

!   Internal variables
INTEGER I                 ! For looping over bins
INTEGER IB                ! Bead index
INTEGER rrdr ! -1 if r, 1 if r+dr
INTEGER IX(2),IY(2),IZ(2)
DOUBLE PRECISION WX(2),WY(2),WZ(2)
DOUBLE PRECISION WTOT       ! total weight ascribed to bin
DOUBLE PRECISION RBIN(3)    ! bead position
INTEGER INDBIN              ! index of bin
INTEGER ISX,ISY,ISZ
LOGICAL isA   ! The bead is of type A

! Copy so I don't have to type mc% everywhere
INTEGER NBINX(3)
double precision temp
NBINX=mc%NBINX

if ((mc%solType.ne.0).or.&
    ((mc%confineType.ne.0).and.(mc%confineType.ne.4))) then
    print*, "change back to not using temp"
    stop 1
endif


! -------------------------------------------------------------
!
!  Calculate change (or value if initialize) of phi for A and B
!
!--------------------------------------------------------------
if (initialize) then
    do I=1,mc%NBIN
       md%PHIA(I)=0.0_dp
       md%DPHIA(I)=0.0_dp
       md%PHIB(I)=0.0_dp
       md%DPHIB(I)=0.0_dp
       md%INDPHI(I)=0
    enddo
    md%DEKap=0
    md%DECouple=0
    md%DEChi=0
endif

md%NPHI=0
do IB=I1,I2
  do rrdr=-1,1,2
   ! on initialize only add current position
   ! otherwise subract current and add new
   if (initialize.and.(rrdr.eq.-1)) CYCLE
   if ((rrdr.eq.-1).or.initialize) then
       RBIN(1)=md%R(IB,1)
       RBIN(2)=md%R(IB,2)
       RBIN(3)=md%R(IB,3)
   else
       RBIN(1)=md%RP(IB,1)
       RBIN(2)=md%RP(IB,2)
       RBIN(3)=md%RP(IB,3)
   endif
   isA=md%AB(IB).eq.1
   temp=rrdr*mc%beadVolume/(mc%dbin**3)
   ! --------------------------------------------------
   !
   !  Interpolate beads into bins
   !
   ! --------------------------------------------------
   call interp(mc%confineType,RBIN,mc%LBOX,mc%NBINX,mc%dbin,IX,IY,IZ,WX,WY,WZ)

   ! -------------------------------------------------------
   !
   ! Count beads in bins
   !
   ! ------------------------------------------------------
   !   Add or Subtract volume fraction with weighting from each bin
   !   I know that it looks bad to have this section of code twice but it
   !   makes it faster.
   if (isA) then
       do ISX=1,2
          if ((IX(ISX).le.0).OR.(IX(ISX).ge.(NBINX(1)+1))) CYCLE
          do ISY=1,2
             if ((IY(ISY).le.0).OR.(IY(ISY).ge.(NBINX(2)+1))) CYCLE
             do ISZ=1,2
                if ((IZ(ISZ).le.0).OR.(IZ(ISZ).ge.(NBINX(3)+1))) cycle
                WTOT=WX(ISX)*WY(ISY)*WZ(ISZ)
                INDBIN=IX(ISX)+(IY(ISY)-1)*NBINX(1)+(IZ(ISZ)-1)*NBINX(1)*NBINX(2)
                if (initialize) then
                    ! Set all phi values on initialize
                    md%PHIA(INDBIN)=md%PHIA(INDBIN)+WTOT*mc%beadVolume/md%Vol(INDBIN)
                else
                    ! Generate list of which phi's change and by how much
                    I=md%NPHI
                    do
                       if (I.eq.0) then
                          md%NPHI=md%NPHI+1
                          md%INDPHI(md%NPHI)=INDBIN
                          !md%DPHIA(mc%NPHI)=rrdr*WTOT*mc%V/md%Vol(INDBIN)
                          md%DPHIA(md%NPHI)=temp*WTOT
                          md%DPHIB(md%NPHI)=0.0_dp
                          exit
                       elseif (INDBIN.EQ.md%INDPHI(I)) then
                          !md%DPHIA(I)=md%DPHIA(I)+rrdr*WTOT*mc%V/md%Vol(INDBIN)
                          md%DPHIA(I)=md%DPHIA(I)+temp*WTOT
                          exit
                       else
                          I=I-1
                       endif
                    enddo
                endif
             enddo
          enddo
       enddo
   else
       do ISX=1,2
          if ((IX(ISX).le.0).OR.(IX(ISX).ge.(NBINX(1)+1))) CYCLE
          do ISY=1,2
             if ((IY(ISY).le.0).OR.(IY(ISY).ge.(NBINX(2)+1))) CYCLE
             do ISZ=1,2
                if ((IZ(ISZ).le.0).OR.(IZ(ISZ).ge.(NBINX(3)+1))) cycle
                WTOT=WX(ISX)*WY(ISY)*WZ(ISZ)
                INDBIN=IX(ISX)+(IY(ISY)-1)*NBINX(1)+(IZ(ISZ)-1)*NBINX(1)*NBINX(2)
                if (initialize) then
                    ! Set all phi values on initialize
                    md%PHIB(INDBIN)=md%PHIB(INDBIN)+WTOT*mc%beadVolume/md%Vol(INDBIN)
                else
                    ! Generate list of which phi's change and by how much
                    I=md%NPHI
                    do
                       if (I.eq.0) then
                          md%NPHI=md%NPHI+1
                          md%INDPHI(md%NPHI)=INDBIN
                          md%DPHIA(md%NPHI)=0.0_dp
                          !md%DPHIB(mc%NPHI)=rrdr*WTOT*mc%V/md%Vol(INDBIN)
                          md%DPHIB(md%NPHI)=temp*WTOT
                          exit
                       elseif (INDBIN.EQ.md%INDPHI(I)) then
                          !md%DPHIB(I)=md%DPHIB(I)+rrdr*WTOT*mc%V/md%Vol(INDBIN)
                          md%DPHIB(I)=md%DPHIB(I)+temp*WTOT
                          exit
                       else
                          I=I-1
                       endif
                    enddo
                endif
             enddo !ISZ
          enddo !ISY
       enddo !ISX
   endif
 enddo ! loop over rrdr.  A.k.a new and old
enddo ! loop over IB  A.k.a. beads
! ---------------------------------------------------------------------
!
! Calcualte change in (or abosolute) energy
!
!---------------------------------------------------------------------
call hamiltonian(mc,md,initialize)

RETURN
END

!---------------------------------------------------------------!
