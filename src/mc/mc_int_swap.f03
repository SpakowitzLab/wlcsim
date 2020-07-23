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
subroutine mc_int_swap(wlc_p,I1,I2,I3,I4)
! values from wlcsim_data
use params, only: wlc_UP, wlc_DPHIA, wlc_RP, wlc_NPHI&
    , wlc_AB, wlc_R, wlc_dPHI_l2, wlc_U, wlc_DPHIB, wlc_inDPHI&
    , wlc_ind_in_list
use params,only:dp,wlcsim_params
use polydispersity, only: length_of_chain_containing
use energies, only: energyOf, kap_, maierSaupe_
implicit none

TYPE(wlcsim_params), intent(inout) :: wlc_p   ! <---- Contains output
integer, intent(in) :: I1      ! Test bead position 1
integer, intent(in) :: I2      ! Test bead position 2
integer, intent(in) :: I3      ! Test bead, first bead of second section
integer, intent(in) :: I4      ! Test bead, second bead of second section

!   Internal variables
integer I, J                 ! For looping over bins
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
real(dp) temp    !for speeding up code
real(dp), dimension(-2:2) :: phi2
real(dp) change


if (WLC_P__TWO_TAIL) then
    print*, "The swap move is not currently set up for two Tail"
    stop 1
endif
if (I2-I1 + 1.ne.length_of_chain_containing(I1)) then
    print*, "Error in mc_int_swap. I2-I1 + 1.ne.NB"
    stop 1
endif
if (I4-I3 + 1.ne.length_of_chain_containing(I3)) then
    print*, "Error in mc_int_swap. I2-I1 + 1.ne.NB"
    stop 1
endif
if (length_of_chain_containing(I1).ne.length_of_chain_containing(I3)) then
    print*, "Error, polymers must be of same length"
    stop 1
endif
if (.not.(min(I1,I2)>max(I3,I4) .or. min(I3,I4)>max(I1,I2))) then
    print*, "Error in mc_int_swap. Overlappling regions"
    print*, I1,I2,"|",I3,I4
endif
if (WLC_P__FIELDINTERACTIONTYPE == 'chromatin2') then
    print*, "swap move not set up for chromatin2."
    stop
endif


! -------------------------------------------------------------
!
!  Calculate change of phi for A and B
!
!--------------------------------------------------------------

wlc_NPHI = 0
do IB = I1,I2
  IB2 = IB + I3-I1
  !No need to do calculation if identities are the same
  if (wlc_AB(IB).eq.wlc_AB(IB2)) cycle
  do rrdr = -1,1,2
   ! on initialize only add current position
   ! otherwise subract current and add new
   if (rrdr.eq.-1) then
       RBin(1) = wlc_R(1,IB)
       RBin(2) = wlc_R(2,IB)
       RBin(3) = wlc_R(3,IB)
   else
       RBin(1) = wlc_RP(1,IB)
       RBin(2) = wlc_RP(2,IB)
       RBin(3) = wlc_RP(3,IB)
   endif
   AminusB = -1+2*wlc_AB(IB) ! -1 if B and +1 if A
   if (WLC_P__CHI_L2_ABLE .and. energyOf(maierSaupe_)%isOn) then
       if (rrdr == -1) then
           call y2_calc(wlc_U(:,IB),phi2)
       else
           call y2_calc(wlc_UP(:,IB),phi2)
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
   change = real(AminusB*rrdr,dp)*WLC_P__BEADVOLUME/(WLC_P__DBIN**3)
   do ISZ = 1,2
      do ISY = 1,2
         do ISX = 1,2
            WTOT = WX(ISX)*WY(ISY)*WZ(ISZ)
            inDBin = IX(ISX) + (IY(ISY)-1)*WLC_P__NBIN_X + (IZ(ISZ)-1)*WLC_P__NBIN_X*WLC_P__NBIN_Y
            temp = WTOT*change
            ! Generate list of which phi's change and by how much
            I = wlc_ind_in_list(indBin)
            if (I == -1) then
                wlc_NPHI = wlc_NPHI + 1
                wlc_ind_in_list(indBin) = wlc_NPHI
                wlc_inDPHI(wlc_NPHI) = inDBin
                wlc_DPHIA(wlc_NPHI) = temp
                wlc_DPHIB(wlc_NPHI) = -temp
                if(WLC_P__CHI_L2_ABLE .and. energyOf(maierSaupe_)%isOn) then
                    do m_index = -2,2
                        wlc_dPHI_l2(m_index,wlc_NPHI) = phi2(m_index)*temp
                    enddo
                endif
            else
                wlc_DPHIA(I) = wlc_DPHIA(I) + temp
                wlc_DPHIB(I) = wlc_DPHIB(I)-temp
                if(WLC_P__CHI_L2_ABLE .and. energyOf(maierSaupe_)%isOn) then
                    do m_index = -2,2
                        wlc_dPHI_l2(m_index,I) = wlc_dPHI_l2(m_index,I) + &
                                    phi2(m_index)*temp
                    enddo
                endif
            endif
         enddo
      enddo
   enddo
 enddo ! loop over rrdr.  A.k.a new and old
enddo ! loop over IB  A.k.a. beads
do I = 1,wlc_NPHI
   J = wlc_inDPHI(I)
   wlc_ind_in_list(J) = -1
enddo
call hamiltonian(wlc_p,.false.)

if (abs(energyOf(kap_)%dx).gt.0.0001_dp) then
    print*, "Error in mc_int_swap.  Kappa energy shouldn't change on move 9"
    print*, "DEKap", energyOf(kap_)%dx
    stop 1
endif
RETURN
END

!---------------------------------------------------------------!
