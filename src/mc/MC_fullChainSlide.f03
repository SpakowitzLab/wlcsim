#include "../defines.inc"
!--------------------------------------------------------------*
!
!           Makes Monti Carlo Moves
!
!    Quinn separated out this file on 8/9/17
!
!---------------------------------------------------------------

! variables that need to be allocated only on certain branches moved into MD to prevent segfaults
! please move other variables in as you see fit
subroutine MC_fullChainSlide(wlc_d,IB1,IB2,IT1,IT2,MCAMP,rand_stat)

use mersenne_twister
use params, only: dp, wlcsim_data

implicit none
type(wlcsim_data), intent(inout) :: wlc_d
integer, intent(out) :: IB1   ! Test bead position 1
integer, intent(out) :: IT1   ! Index of test bead 1
integer, intent(out) :: IB2   ! Test bead position 2
integer, intent(out) :: IT2   ! Index of test bead 2

integer IP    ! Test polymer
integer I ! Test indices
! Things for random number generator
type(random_stat), intent(inout) :: rand_stat  ! status of random number generator
real urand(3)  ! random vector
integer irnd(1)
! Variables for the crank-shaft move

real(dp) P1(3)    ! Point on rotation line
real(dp), intent(in) :: MCAMP ! Amplitude of random change
real(dp) DR(3)    ! Displacement for slide move


!TOdo saving RP is not actually needed, even in these cases, but Brad's code assumes that we have RP.
if (WLC_P__RING .OR. WLC_P__INTERP_BEAD_LENNARD_JONES) then
    wlc_d%RP = wlc_d%R
    wlc_d%UP = wlc_d%U
    P1 = 0.0_dp
endif

!     Perform full chain slide move (MCTYPE 6)

call random_index(WLC_P__NP,irnd,rand_stat)
IP=irnd(1)
IB1 = 1
IB2 = WLC_P__NB
IT1 = WLC_P__NB*(IP-1) + IB1
IT2 = WLC_P__NB*(IP-1) + IB2

call random_number(urand,rand_stat)
DR(1) = MCAMP*(urand(1)-0.5_dp)
DR(2) = MCAMP*(urand(2)-0.5_dp)
DR(3) = MCAMP*(urand(3)-0.5_dp)

do I = IT1,IT2
    wlc_d%RP(:,I) = wlc_d%R(:,I) + DR
enddo
wlc_d%UP(:,IT1:IT2) = wlc_d%U(:,IT1:IT2)
if (WLC_P__LOCAL_TWIST)  wlc_d%VP(:,IT1:IT2) = wlc_d%V(:,IT1:IT2)
end subroutine
