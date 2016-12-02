!---------------------------------------------------------------!
!     
!     This subroutine calculates the field energy from scratch.
!     Also sets PhiA and PhiB to 0.
!     Set x_chi, x
!
!     Andrew Spakowitz
!     Written 6-29-04
!
!     Edited by Shifan
!     Edited by Quinn in 2016

SUBROUTINE MC_int_initialize(mc,md)
use params, only: dp, wlcsim_params, wlcsim_data
IMPLICIT NONE

TYPE(wlcsim_params), intent(in) :: mc   ! <---- Contains output
TYPE(wlcsim_data), intent(inout) :: md
LOGICAL initialize   ! if true, calculate absolute energy

!   Internal variables
INTEGER I                 ! For looping over bins
INTEGER IB                ! Bead index
INTEGER IX(2),IY(2),IZ(2)
real(dp) WX(2),WY(2),WZ(2)
real(dp) WTOT       ! total weight ascribed to bin
real(dp) RBIN(3)    ! bead position
INTEGER INDBIN              ! index of bin
INTEGER ISX,ISY,ISZ
LOGICAL isA   ! The bead is of type A

! Copy so I don't have to type mc% everywhere
INTEGER NBINX(3)
NBINX=mc%NBINX

! -------------------------------------------------------------
!
!  Calculate change (or value if initialize) of phi for A and B
!
!--------------------------------------------------------------
do I=1,mc%NBIN
   md%PHIA(I)=0.0_dp
   md%DPHIA(I)=0.0_dp
   md%PHIB(I)=0.0_dp
   md%DPHIB(I)=0.0_dp
   md%INDPHI(I)=0
enddo

md%NPHI=0
do IB=1,mc%NT
   RBIN(1)=md%R(IB,1)
   RBIN(2)=md%R(IB,2)
   RBIN(3)=md%R(IB,3)

   isA=md%AB(IB).eq.1
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

                ! Set all phi values on initialize
                md%PHIA(INDBIN)=md%PHIA(INDBIN)+WTOT*mc%beadVolume/md%Vol(INDBIN)

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

                ! Set all phi values on initialize
                md%PHIB(INDBIN)=md%PHIB(INDBIN)+WTOT*mc%beadVolume/md%Vol(INDBIN)

             enddo !ISZ
          enddo !ISY
       enddo !ISX
   endif
enddo ! loop over IB  A.k.a. beads
! ---------------------------------------------------------------------
!
! Calcualte abosolute energy
!
!---------------------------------------------------------------------
initialize=.True.
call hamiltonian(mc,md,initialize)

RETURN
END

!---------------------------------------------------------------!
!---------------------------------------------------------------!

!
!     This subroutine calculates the change in the self energy for
!     a small Monte Carlo move in the position.
!
!     Andrew Spakowitz
!     Written 6-29-04
!
!     Edited by Shifan
!     Edited by Quinn in 2016
!--------------------------------------------------------------!
SUBROUTINE MC_int_update(mc,md,I1,I2)
use params, only: dp, wlcsim_params, wlcsim_data
IMPLICIT NONE

TYPE(wlcsim_params), intent(in) :: mc   ! <---- Contains output
TYPE(wlcsim_data), intent(inout) :: md
LOGICAL initialize   ! if true, calculate absolute energy
INTEGER, intent(in) :: I1           ! Test bead position 1
INTEGER, intent(in) :: I2           ! Test bead position 2

!   Internal variables
INTEGER I                 ! For looping over bins
INTEGER IB                ! Bead index
INTEGER rrdr ! -1 if r, 1 if r+dr
INTEGER IX(2),IY(2),IZ(2)
real(dp) WX(2),WY(2),WZ(2)
real(dp) WTOT       ! total weight ascribed to bin
real(dp) RBIN(3)    ! bead position
INTEGER INDBIN              ! index of bin
INTEGER ISX,ISY,ISZ
LOGICAL isA   ! The bead is of type A

! Copy so I don't have to type mc% everywhere
INTEGER NBINX(3)
double precision temp

NBINX=mc%NBINX

md%NPHI=0
do IB=I1,I2
  do rrdr=-1,1,2
   ! on initialize only add current position
   ! otherwise subract current and add new
   if (rrdr.eq.-1) then
       RBIN(1)=md%R(IB,1)
       RBIN(2)=md%R(IB,2)
       RBIN(3)=md%R(IB,3)
   else
       RBIN(1)=md%RP(IB,1)
       RBIN(2)=md%RP(IB,2)
       RBIN(3)=md%RP(IB,3)
   endif
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
   isA=md%AB(IB).eq.1
   if (mc%confineType==0 .or. mc%confineType==4) then
       ! If periodic than you can assume that all bins are included and have a volume
       ! of dbin**3
       temp=rrdr*mc%beadVolume/(mc%dbin**3)
       if (isA) then
           do ISX=1,2
              do ISY=1,2
                 do ISZ=1,2
                    WTOT=WX(ISX)*WY(ISY)*WZ(ISZ)
                    INDBIN=IX(ISX)+(IY(ISY)-1)*NBINX(1)+(IZ(ISZ)-1)*NBINX(1)*NBINX(2)

                    ! Generate list of which phi's change and by how much
                    I=md%NPHI
                    do
                       if (I.eq.0) then
                          md%NPHI=md%NPHI+1
                          md%INDPHI(md%NPHI)=INDBIN
                          md%DPHIA(md%NPHI)=temp*WTOT
                          md%DPHIB(md%NPHI)=0.0_dp
                          exit
                       elseif (INDBIN.EQ.md%INDPHI(I)) then
                          md%DPHIA(I)=md%DPHIA(I)+temp*WTOT
                          exit
                       else
                          I=I-1
                       endif
                    enddo
                 enddo
              enddo
           enddo
       else
           do ISX=1,2
              do ISY=1,2
                 do ISZ=1,2
                    WTOT=WX(ISX)*WY(ISY)*WZ(ISZ)
                    INDBIN=IX(ISX)+(IY(ISY)-1)*NBINX(1)+(IZ(ISZ)-1)*NBINX(1)*NBINX(2)

                    ! Generate list of which phi's change and by how much
                    I=md%NPHI
                    do
                       if (I.eq.0) then
                          md%NPHI=md%NPHI+1
                          md%INDPHI(md%NPHI)=INDBIN
                          md%DPHIA(md%NPHI)=0.0_dp
                          md%DPHIB(md%NPHI)=temp*WTOT
                          exit
                       elseif (INDBIN.EQ.md%INDPHI(I)) then
                          md%DPHIB(I)=md%DPHIB(I)+temp*WTOT
                          exit
                       else
                          I=I-1
                       endif
                    enddo
                 enddo !ISZ
              enddo !ISY
           enddo !ISX
       endif
    else ! not periodic
       ! if not periodic you need to check wheter bin is outside and
       ! also need to divide by volume of bin 
       if (isA) then
           do ISX=1,2
              if ((IX(ISX).le.0).OR.(IX(ISX).ge.(NBINX(1)+1))) CYCLE
              do ISY=1,2
                 if ((IY(ISY).le.0).OR.(IY(ISY).ge.(NBINX(2)+1))) CYCLE
                 do ISZ=1,2
                    if ((IZ(ISZ).le.0).OR.(IZ(ISZ).ge.(NBINX(3)+1))) cycle
                    WTOT=WX(ISX)*WY(ISY)*WZ(ISZ)
                    INDBIN=IX(ISX)+(IY(ISY)-1)*NBINX(1)+(IZ(ISZ)-1)*NBINX(1)*NBINX(2)

                    ! Generate list of which phi's change and by how much
                    I=md%NPHI
                    do
                       if (I.eq.0) then
                          md%NPHI=md%NPHI+1
                          md%INDPHI(md%NPHI)=INDBIN
                          md%DPHIA(md%NPHI)=rrdr*WTOT*mc%beadVolume/md%Vol(INDBIN)
                          !md%DPHIA(md%NPHI)=temp*WTOT
                          md%DPHIB(md%NPHI)=0.0_dp
                          exit
                       elseif (INDBIN.EQ.md%INDPHI(I)) then
                          md%DPHIA(I)=md%DPHIA(I)+rrdr*WTOT*mc%beadVolume/md%Vol(INDBIN)
                          !md%DPHIA(I)=md%DPHIA(I)+temp*WTOT
                          exit
                       else
                          I=I-1
                       endif
                    enddo
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

                    ! Generate list of which phi's change and by how much
                    I=md%NPHI
                    do
                       if (I.eq.0) then
                          md%NPHI=md%NPHI+1
                          md%INDPHI(md%NPHI)=INDBIN
                          md%DPHIA(md%NPHI)=0.0_dp
                          md%DPHIB(md%NPHI)=rrdr*WTOT*mc%beadVolume/md%Vol(INDBIN)
                          !md%DPHIB(md%NPHI)=temp*WTOT
                          exit
                       elseif (INDBIN.EQ.md%INDPHI(I)) then
                          md%DPHIB(I)=md%DPHIB(I)+rrdr*WTOT*mc%beadVolume/md%Vol(INDBIN)
                          !md%DPHIB(I)=md%DPHIB(I)+temp*WTOT
                          exit
                       else
                          I=I-1
                       endif
                    enddo
                 enddo !ISZ
              enddo !ISY
           enddo !ISX
       endif
    endif

 enddo ! loop over rrdr.  A.k.a new and old
enddo ! loop over IB  A.k.a. beads
! ---------------------------------------------------------------------
!
! Calcualte change in energy
!
!---------------------------------------------------------------------
initialize=.False.
call hamiltonian(mc,md,initialize)

RETURN
END

!---------------------------------------------------------------!
