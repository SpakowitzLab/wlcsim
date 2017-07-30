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

subroutine MC_int_initialize(wlc_p,wlc_d)
use params, only: dp, wlcsim_params, wlcsim_data
implicit none

TYPE(wlcsim_params), intent(in) :: wlc_p   ! <---- Contains output
TYPE(wlcsim_data), intent(inout) :: wlc_d
LOGICAL initialize   ! if true, calculate absolute energy

!   Internal variables
integer I                 ! For looping over bins
integer IB                ! Bead index
integer IX(2),IY(2),IZ(2)
real(dp) WX(2),WY(2),WZ(2)
real(dp) WTOT       ! total weight ascribed to bin
real(dp) RBin(3)    ! bead position
integer inDBin              ! index of bin
integer ISX,ISY,ISZ
LOGICAL isA   ! The bead is of type A

! Copy so I don't have to type wlc_p% everywhere
integer NBinX(3)
NBinX = wlc_p%NBinX

! -------------------------------------------------------------
!
!  Calculate change (or value if initialize) of phi for A and B
!
!--------------------------------------------------------------
do I = 1,wlc_p%NBin
   wlc_d%PHIA(I) = 0.0_dp
   wlc_d%DPHIA(I) = 0.0_dp
   wlc_d%PHIB(I) = 0.0_dp
   wlc_d%DPHIB(I) = 0.0_dp
   wlc_d%inDPHI(I) = 0
enddo

wlc_d%NPHI = 0
do IB = 1,wlc_p%NT
   RBin(1) = wlc_d%R(IB,1)
   RBin(2) = wlc_d%R(IB,2)
   RBin(3) = wlc_d%R(IB,3)

   isA = wlc_d%AB(IB).eq.1
   ! --------------------------------------------------
   !
   !  Interpolate beads into bins
   !
   ! --------------------------------------------------
   call interp(wlc_p%confineType,RBin,wlc_p%LBOX,wlc_p%NBinX,wlc_p%dbin,IX,IY,IZ,WX,WY,WZ)

   ! -------------------------------------------------------
   !
   ! Count beads in bins
   !
   ! ------------------------------------------------------
   !   Add or Subtract volume fraction with weighting from each bin
   !   I know that it looks bad to have this section of code twice but it
   !   makes it faster.
   if (isA) then
       do ISX = 1,2
          if ((IX(ISX).le.0).OR.(IX(ISX).ge.(NBinX(1) + 1))) CYCLE
          do ISY = 1,2
             if ((IY(ISY).le.0).OR.(IY(ISY).ge.(NBinX(2) + 1))) CYCLE
             do ISZ = 1,2
                if ((IZ(ISZ).le.0).OR.(IZ(ISZ).ge.(NBinX(3) + 1))) cycle
                WTOT = WX(ISX)*WY(ISY)*WZ(ISZ)
                inDBin = IX(ISX) + (IY(ISY)-1)*NBinX(1) + (IZ(ISZ)-1)*NBinX(1)*NBinX(2)

                ! Set all phi values on initialize
                wlc_d%PHIA(inDBin) = wlc_d%PHIA(inDBin) + WTOT*wlc_p%beadVolume/wlc_d%Vol(inDBin)

             enddo
          enddo
       enddo
   else
       do ISX = 1,2
          if ((IX(ISX).le.0).OR.(IX(ISX).ge.(NBinX(1) + 1))) CYCLE
          do ISY = 1,2
             if ((IY(ISY).le.0).OR.(IY(ISY).ge.(NBinX(2) + 1))) CYCLE
             do ISZ = 1,2
                if ((IZ(ISZ).le.0).OR.(IZ(ISZ).ge.(NBinX(3) + 1))) cycle
                WTOT = WX(ISX)*WY(ISY)*WZ(ISZ)
                inDBin = IX(ISX) + (IY(ISY)-1)*NBinX(1) + (IZ(ISZ)-1)*NBinX(1)*NBinX(2)

                ! Set all phi values on initialize
                wlc_d%PHIB(inDBin) = wlc_d%PHIB(inDBin) + WTOT*wlc_p%beadVolume/wlc_d%Vol(inDBin)

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
initialize = .True.
call hamiltonian(wlc_p,wlc_d,initialize)

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
subroutine MC_int_update(wlc_p,wlc_d,I1,I2)
use params, only: dp, wlcsim_params, wlcsim_data
implicit none

TYPE(wlcsim_params), intent(in) :: wlc_p   ! <---- Contains output
TYPE(wlcsim_data), intent(inout) :: wlc_d
LOGICAL initialize   ! if true, calculate absolute energy
integer, intent(in) :: I1           ! Test bead position 1
integer, intent(in) :: I2           ! Test bead position 2

!   Internal variables
integer I                 ! For looping over bins
integer IB                ! Bead index
integer rrdr ! -1 if r, 1 if r + dr
integer IX(2),IY(2),IZ(2)
real(dp) WX(2),WY(2),WZ(2)
real(dp) WTOT       ! total weight ascribed to bin
real(dp) RBin(3)    ! bead position
integer inDBin              ! index of bin
integer ISX,ISY,ISZ
LOGICAL isA   ! The bead is of type A

! Copy so I don't have to type wlc_p% everywhere
integer NBinX(3)
real(dp) temp

NBinX = wlc_p%NBinX

wlc_d%NPHI = 0
do IB = I1,I2
  do rrdr = -1,1,2
   ! on initialize only add current position
   ! otherwise subract current and add new
   if (rrdr.eq.-1) then
       RBin(1) = wlc_d%R(IB,1)
       RBin(2) = wlc_d%R(IB,2)
       RBin(3) = wlc_d%R(IB,3)
   else
       RBin(1) = wlc_d%RP(IB,1)
       RBin(2) = wlc_d%RP(IB,2)
       RBin(3) = wlc_d%RP(IB,3)
   endif
   ! --------------------------------------------------
   !
   !  Interpolate beads into bins
   !
   ! --------------------------------------------------
   call interp(wlc_p%confineType,RBin,wlc_p%LBOX,wlc_p%NBinX,wlc_p%dbin,IX,IY,IZ,WX,WY,WZ)

   ! -------------------------------------------------------
   !
   ! Count beads in bins
   !
   ! ------------------------------------------------------
   !   Add or Subtract volume fraction with weighting from each bin
   !   I know that it looks bad to have this section of code twice but it
   !   makes it faster.
   isA = wlc_d%AB(IB).eq.1
   if (wlc_p%confineType == 0 .or. wlc_p%confineType == 4) then
       ! If periodic than you can assume that all bins are included and have a volume
       ! of dbin**3
       temp = rrdr*wlc_p%beadVolume/(wlc_p%dbin**3)
       if (isA) then
           do ISX = 1,2
              do ISY = 1,2
                 do ISZ = 1,2
                    WTOT = WX(ISX)*WY(ISY)*WZ(ISZ)
                    inDBin = IX(ISX) + (IY(ISY)-1)*NBinX(1) + (IZ(ISZ)-1)*NBinX(1)*NBinX(2)

                    ! Generate list of which phi's change and by how much
                    I = wlc_d%NPHI
                    do
                       if (I.eq.0) then
                          wlc_d%NPHI = wlc_d%NPHI + 1
                          wlc_d%inDPHI(wlc_d%NPHI) = inDBin
                          wlc_d%DPHIA(wlc_d%NPHI) = temp*WTOT
                          wlc_d%DPHIB(wlc_d%NPHI) = 0.0_dp
                          exit
                       elseif (inDBin == wlc_d%inDPHI(I)) then
                          wlc_d%DPHIA(I) = wlc_d%DPHIA(I) + temp*WTOT
                          exit
                       else
                          I = I-1
                       endif
                    enddo
                 enddo
              enddo
           enddo
       else
           do ISX = 1,2
              do ISY = 1,2
                 do ISZ = 1,2
                    WTOT = WX(ISX)*WY(ISY)*WZ(ISZ)
                    inDBin = IX(ISX) + (IY(ISY)-1)*NBinX(1) + (IZ(ISZ)-1)*NBinX(1)*NBinX(2)

                    ! Generate list of which phi's change and by how much
                    I = wlc_d%NPHI
                    do
                       if (I.eq.0) then
                          wlc_d%NPHI = wlc_d%NPHI + 1
                          wlc_d%inDPHI(wlc_d%NPHI) = inDBin
                          wlc_d%DPHIA(wlc_d%NPHI) = 0.0_dp
                          wlc_d%DPHIB(wlc_d%NPHI) = temp*WTOT
                          exit
                       elseif (inDBin == wlc_d%inDPHI(I)) then
                          wlc_d%DPHIB(I) = wlc_d%DPHIB(I) + temp*WTOT
                          exit
                       else
                          I = I-1
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
           do ISX = 1,2
              if ((IX(ISX).le.0).OR.(IX(ISX).ge.(NBinX(1) + 1))) CYCLE
              do ISY = 1,2
                 if ((IY(ISY).le.0).OR.(IY(ISY).ge.(NBinX(2) + 1))) CYCLE
                 do ISZ = 1,2
                    if ((IZ(ISZ).le.0).OR.(IZ(ISZ).ge.(NBinX(3) + 1))) cycle
                    WTOT = WX(ISX)*WY(ISY)*WZ(ISZ)
                    inDBin = IX(ISX) + (IY(ISY)-1)*NBinX(1) + (IZ(ISZ)-1)*NBinX(1)*NBinX(2)

                    ! Generate list of which phi's change and by how much
                    I = wlc_d%NPHI
                    do
                       if (I.eq.0) then
                          wlc_d%NPHI = wlc_d%NPHI + 1
                          wlc_d%inDPHI(wlc_d%NPHI) = inDBin
                          wlc_d%DPHIA(wlc_d%NPHI) = rrdr*WTOT*wlc_p%beadVolume/wlc_d%Vol(inDBin)
                          !wlc_d%DPHIA(wlc_d%NPHI) = temp*WTOT
                          wlc_d%DPHIB(wlc_d%NPHI) = 0.0_dp
                          exit
                       elseif (inDBin == wlc_d%inDPHI(I)) then
                          wlc_d%DPHIA(I) = wlc_d%DPHIA(I) + rrdr*WTOT*wlc_p%beadVolume/wlc_d%Vol(inDBin)
                          !wlc_d%DPHIA(I) = wlc_d%DPHIA(I) + temp*WTOT
                          exit
                       else
                          I = I-1
                       endif
                    enddo
                 enddo
              enddo
           enddo
       else
           do ISX = 1,2
              if ((IX(ISX).le.0).OR.(IX(ISX).ge.(NBinX(1) + 1))) CYCLE
              do ISY = 1,2
                 if ((IY(ISY).le.0).OR.(IY(ISY).ge.(NBinX(2) + 1))) CYCLE
                 do ISZ = 1,2
                    if ((IZ(ISZ).le.0).OR.(IZ(ISZ).ge.(NBinX(3) + 1))) cycle
                    WTOT = WX(ISX)*WY(ISY)*WZ(ISZ)
                    inDBin = IX(ISX) + (IY(ISY)-1)*NBinX(1) + (IZ(ISZ)-1)*NBinX(1)*NBinX(2)

                    ! Generate list of which phi's change and by how much
                    I = wlc_d%NPHI
                    do
                       if (I.eq.0) then
                          wlc_d%NPHI = wlc_d%NPHI + 1
                          wlc_d%inDPHI(wlc_d%NPHI) = inDBin
                          wlc_d%DPHIA(wlc_d%NPHI) = 0.0_dp
                          wlc_d%DPHIB(wlc_d%NPHI) = rrdr*WTOT*wlc_p%beadVolume/wlc_d%Vol(inDBin)
                          !wlc_d%DPHIB(wlc_d%NPHI) = temp*WTOT
                          exit
                       elseif (inDBin == wlc_d%inDPHI(I)) then
                          wlc_d%DPHIB(I) = wlc_d%DPHIB(I) + rrdr*WTOT*wlc_p%beadVolume/wlc_d%Vol(inDBin)
                          !wlc_d%DPHIB(I) = wlc_d%DPHIB(I) + temp*WTOT
                          exit
                       else
                          I = I-1
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
initialize = .False.
call hamiltonian(wlc_p,wlc_d,initialize)

RETURN
END

!---------------------------------------------------------------!
