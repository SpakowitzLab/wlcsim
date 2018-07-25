#include "../defines.inc"
!---------------------------------------------------------------!
!
!
!     This subroutine calculates the change in field interaction
!     energy for a MC move which swaps two chains.
!
!
!     Written by Quinn in 2016 based on code from Andrew Spakowitz and Shifan
!
!---------------------------------------------------------------!
subroutine MC_int_swap(wlc_p,wlc_d,I1,I2,I3,I4)
use params,only:dp,wlcsim_params,wlcsim_data
implicit none

TYPE(wlcsim_params), intent(inout) :: wlc_p   ! <---- Contains output
TYPE(wlcsim_data), intent(inout) :: wlc_d
integer, intent(in) :: I1      ! Test bead position 1
integer, intent(in) :: I2      ! Test bead position 2
integer, intent(in) :: I3      ! Test bead, first bead of second section
integer, intent(in) :: I4      ! Test bead, second bead of second section

!   Internal variables
integer I                 ! For looping over bins
integer IB                ! Bead index
integer IB2               ! Index you are swapping with
integer rrdr ! -1 if r, 1 if r + dr
integer IX(2),IY(2),IZ(2)
real(dp) WX(2),WY(2),WZ(2)
real(dp) WTOT       ! total weight ascribed to bin
real(dp) RBin(3)    ! bead position
integer inDBin              ! index of bin
integer ISX,ISY,ISZ
!LOGICAL isA   ! The bead is of type A
integer AminusB
integer m_index ! m quantum number for sphrical harmonics
integer NBinX(3)
real(dp) temp    !for speeding up code
real(dp), dimension(-2:2) :: phi2
NBinX = wlc_p%NBINX


if (WLC_P__TWO_TAIL) then
    print*, "The swap move is not currently set up for two Tail"
    stop 1
endif
if (I2-I1 + 1.ne.WLC_P__NB) then
    print*, "Error in MC_int_swap. I2-I1 + 1.ne.NB"
    stop 1
endif
if (I4-I3 + 1.ne.WLC_P__NB) then
    print*, "Error in MC_int_swap. I2-I1 + 1.ne.NB"
    stop 1
endif
if (.not.(min(I1,I2)>max(I3,I4) .or. min(I3,I4)>max(I1,I2))) then
    print*, "Error in MC_int_swap. Overlappling regions"
    print*, I1,I2,"|",I3,I4
endif


! -------------------------------------------------------------
!
!  Calculate change of phi for A and B
!
!--------------------------------------------------------------

wlc_d%NPHI = 0
do IB = I1,I2
  IB2 = IB + I3-I1
  !No need to do calculation if identities are the same
  if (wlc_d%AB(IB).eq.wlc_d%AB(IB2)) cycle
  do rrdr = -1,1,2
   ! on initialize only add current position
   ! otherwise subract current and add new
   if (rrdr.eq.-1) then
       RBin(1) = wlc_d%R(1,IB)
       RBin(2) = wlc_d%R(2,IB)
       RBin(3) = wlc_d%R(3,IB)
   else
       RBin(1) = wlc_d%RP(1,IB)
       RBin(2) = wlc_d%RP(2,IB)
       RBin(3) = wlc_d%RP(3,IB)
   endif
   AminusB = -1+2*wlc_d%AB(IB) ! -1 if B and +1 if A
   if (wlc_p%CHI_L2_ON) then
       if (rrdr == -1) then
           call Y2calc(wlc_d%U(:,IB),phi2)
       else
           call Y2calc(wlc_d%UP(:,IB),phi2)
       endif
   else
       ! You could give some MS parameter to B as well if you wanted
       phi2=0.0_dp
   endif

   ! --------------------------------------------------
   !
   !  Interpolate beads into bins
   !
   ! --------------------------------------------------
   call interp(wlc_p,RBin,IX,IY,IZ,WX,WY,WZ)

   ! -------------------------------------------------------
   !
   ! Count beads in bins
   !
   ! ------------------------------------------------------
   !   Add or Subtract volume fraction with weighting from each bin
   !   I know that it looks bad to have this section of code twice but it
   !   makes it faster.
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
                  temp = AminusB*rrdr*WTOT*WLC_P__BEADVOLUME/(WLC_P__DBIN**3)
                  wlc_d%DPHIA(wlc_d%NPHI) = temp
                  wlc_d%DPHIB(wlc_d%NPHI) = -temp
                  if(wlc_p%CHI_L2_ON) then
                      do m_index = -2,2
                          wlc_d%dPHI_l2(m_index,wlc_d%NPHI) = phi2(m_index)*temp
                      enddo
                  endif
                  exit
               elseif (inDBin == wlc_d%inDPHI(I)) then
                  temp = AminusB*rrdr*WTOT*WLC_P__BEADVOLUME/(WLC_P__DBIN**3)
                  wlc_d%DPHIA(I) = wlc_d%DPHIA(I) + temp
                  wlc_d%DPHIB(I) = wlc_d%DPHIB(I)-temp
                  if(wlc_p%CHI_L2_ON) then
                      do m_index = -2,2
                          wlc_d%dPHI_l2(m_index,I) = wlc_d%dPHI_l2(m_index,I) + &
                                      phi2(m_index)*temp
                      enddo
                  endif
                  exit
               else
                  I = I-1
               endif
            enddo
         enddo
      enddo
   enddo
 enddo ! loop over rrdr.  A.k.a new and old
enddo ! loop over IB  A.k.a. beads
call hamiltonian(wlc_p,wlc_d,.false.)

if (abs(wlc_d%DEKap).gt.0.0001) then
    print*, "Error in MC_int_swap.  Kappa energy shouldn't change on move 9"
    print*, "DEKap", wlc_d%DEKap
    stop 1
endif
RETURN
END

!---------------------------------------------------------------!
