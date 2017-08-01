!---------------------------------------------------------------!

!
!     This subroutine calculates the change in the self energy for
!     a small Monte Carlo move in the position.
!
!     Andrew Spakowitz
!     Written 6-29-04
!
!     Edited by Shifan
!
!     Edited by Quinn in 2016

subroutine MC_int_rep(wlc_p,wlc_d,I1,I2,forward)
use params
implicit none

!   iputs
TYPE(wlcsim_params), intent(in) :: wlc_p   ! <---- Contains output
TYPE(wlcsim_data), intent(inout) :: wlc_d
integer, intent(in) :: I1  ! Test bead position 1
integer, intent(in) :: I2  ! Test bead position 2

!   Internal variables
integer I                 ! For looping over bins
integer II                ! For looping over IB
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
real(dp) temp    !for speeding up code
LOGICAL, intent(in) :: forward ! move forward

NBinX = wlc_p%NBinX

wlc_d%NPHI = 0
! -------------------------------------------------------------
!
!  Calculate end beads
!
!--------------------------------------------------------------

do II = 1,2
  if (II.eq.1) then
      IB = I1
      if (forward) then
          rrdr = -1
      else
          rrdr = 1
      endif
  elseif (II.eq.2) then
      IB = I2
      if (forward) then
          rrdr = 1
      else
          rrdr = -1
      endif
  else
      print*, "Error in MC_int_rep, II = {1,2}"
      stop 1
  endif
   ! subract current and add new
   if (rrdr.eq.-1) then
       RBin(1) = wlc_d%R(1,IB)
       RBin(2) = wlc_d%R(2,IB)
       RBin(3) = wlc_d%R(3,IB)
   else
       RBin(1) = wlc_d%RP(IB,1)
       RBin(2) = wlc_d%RP(IB,2)
       RBin(3) = wlc_d%RP(IB,3)
   endif
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
                ! Generate list of which phi's change and by how much
                I = wlc_d%NPHI
                do
                   if (I.eq.0) then
                      wlc_d%NPHI = wlc_d%NPHI + 1
                      wlc_d%inDPHI(wlc_d%NPHI) = inDBin
                      wlc_d%DPHIA(wlc_d%NPHI) = rrdr*WTOT*wlc_p%beadVolume/wlc_d%Vol(inDBin)
                      wlc_d%DPHIB(wlc_d%NPHI) = 0.0_dp
                      exit
                   elseif (inDBin == wlc_d%inDPHI(I)) then
                      wlc_d%DPHIA(I) = wlc_d%DPHIA(I) + rrdr*WTOT*wlc_p%beadVolume/wlc_d%Vol(inDBin)
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
                      exit
                   elseif (inDBin == wlc_d%inDPHI(I)) then
                      wlc_d%DPHIB(I) = wlc_d%DPHIB(I) + rrdr*WTOT*wlc_p%beadVolume/wlc_d%Vol(inDBin)
                      exit
                   else
                      I = I-1
                   endif
                enddo
             enddo !ISZ
          enddo !ISY
       enddo !ISX
   endif
enddo ! loop over IB  A.k.a. beads
! ---------------------------------------------------------------------
!



!----------------------------------------------------------
!
!  Now do intermediate Beads
!
!-----------------------------------------------------------
do IB = I1,I2-1
   if (wlc_d%AB(IB).eq.wlc_d%AB(IB + 1)) CYCLE
   if (forward) then
       RBin(1) = wlc_d%RP(IB,1)
       RBin(2) = wlc_d%RP(IB,2)
       RBin(3) = wlc_d%RP(IB,3)
       isA = wlc_d%AB(IB).eq.1
   else
       RBin(1) = wlc_d%R(1,IB)
       RBin(2) = wlc_d%R(2,IB)
       RBin(3) = wlc_d%R(3,IB)
       isA = wlc_d%AB(IB + 1).eq.1
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
                      temp = WTOT*wlc_p%beadVolume/wlc_d%Vol(inDBin)
                      wlc_d%DPHIA(wlc_d%NPHI) = temp
                      wlc_d%DPHIB(wlc_d%NPHI) = -temp
                      exit
                   elseif (inDBin == wlc_d%inDPHI(I)) then
                      temp = WTOT*wlc_p%beadVolume/wlc_d%Vol(inDBin)
                      wlc_d%DPHIA(I) = wlc_d%DPHIA(I) + temp
                      wlc_d%DPHIB(I) = wlc_d%DPHIB(I)-temp
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
                      temp = WTOT*wlc_p%beadVolume/wlc_d%Vol(inDBin)
                      wlc_d%DPHIA(wlc_d%NPHI) = -temp
                      wlc_d%DPHIB(wlc_d%NPHI) = temp
                      exit
                   elseif (inDBin == wlc_d%inDPHI(I)) then
                      temp = WTOT*wlc_p%beadVolume/wlc_d%Vol(inDBin)
                      wlc_d%DPHIA(I) = wlc_d%DPHIA(I)-temp
                      wlc_d%DPHIB(I) = wlc_d%DPHIB(I) + temp
                      exit
                   else
                      I = I-1
                   endif
                enddo
             enddo !ISZ
          enddo !ISY
       enddo !ISX
   endif
enddo ! loop over IB  A.k.a. beads
! ---------------------------------------------------------------------
!
! Calcualte change in energy
!
!---------------------------------------------------------------------
call hamiltonian(wlc_p,wlc_d,.false.)

RETURN
END

!---------------------------------------------------------------!
