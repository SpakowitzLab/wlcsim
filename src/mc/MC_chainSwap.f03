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
subroutine MC_chainSwap(IB1,IB2,IT1,IT2,rand_stat,IT3,IT4)
! values from wlcsim_data
use params, only: wlc_U, wlc_RP, wlc_VP, wlc_UP, wlc_R&
    , wlc_V

use mersenne_twister


implicit none
integer IP    ! Test polymer
integer IP2   ! Second Test polymer
integer, intent(out) :: IB1   ! Test bead position 1
integer, intent(out) :: IT1   ! Index of test bead 1
integer, intent(out) :: IB2   ! Test bead position 2
integer, intent(out) :: IT2   ! Index of test bead 2
integer, intent(out) :: IT3   ! Test bead position 3 if applicable
integer, intent(out) :: IT4   ! Test bead position 4 if applicable

! Things for random number generator
type(random_stat), intent(inout) :: rand_stat  ! status of random number generator
integer irnd(1)
integer I


!TOdo saving RP is not actually needed, even in these cases, but Brad's code assumes that we have RP.
if (WLC_P__RING .OR. WLC_P__INTERP_BEAD_LENNARD_JONES) then
    wlc_RP = wlc_R
    wlc_UP = wlc_U
endif

! switch two chains
call random_index(WLC_P__NP,irnd,rand_stat)
IP=irnd(1)
call random_index(WLC_P__NP,irnd,rand_stat)
IP2=irnd(1)
! Don't switch a chain with itself
if (IP.eq.IP2) then
    IP2 = IP-1
    if (IP2.eq.0) then
        IP2 = 2
    endif
endif
IT1 = WLC_P__NB*(IP-1) + 1
IT2 = WLC_P__NB*(IP-1) + WLC_P__NB
IT3 = WLC_P__NB*(IP2-1) + 1
IT4 = WLC_P__NB*(IP2-1) + WLC_P__NB
do I = 0,WLC_P__NB-1
    wlc_RP(:,IT1 + I) = wlc_R(:,IT3 + I)
    wlc_UP(:,IT1 + I) = wlc_U(:,IT3 + I)
    wlc_RP(:,IT3 + I) = wlc_R(:,IT1 + I)
    wlc_UP(:,IT3 + I) = wlc_U(:,IT1 + I)
    if (WLC_P__LOCAL_TWIST) then
        wlc_VP(:,IT1 + I) = wlc_V(:,IT3 + I)
        wlc_VP(:,IT3 + I) = wlc_V(:,IT1 + I)
    endif
ENDdo
IB1 = -2000000
IB2 = -2000000
end subroutine
