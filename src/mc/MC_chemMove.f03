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
subroutine MC_chemMove(IB1,IB2,IT1,IT2,WindoW,rand_stat,success)
! values from wlcsim_data
use params, only: wlc_V, wlc_R, wlc_RP, wlc_AB, wlc_U&
    , wlc_UP, wlc_ABP, wlc_VP

use mersenne_twister
use params, only: dp
use windowTools, only: drawWindow

implicit none
integer, intent(out) :: IB1   ! Test bead position 1
integer, intent(out) :: IT1   ! Index of test bead 1
integer, intent(out) :: IB2   ! Test bead position 2
integer, intent(out) :: IT2   ! Index of test bead 2
logical, intent(out) ::success
integer dib   ! number of beads moved by move

integer IP    ! Test polymer
integer J  ! Test indices
! Things for random number generator
type(random_stat), intent(inout) :: rand_stat  ! status of random number generator
real(dp) urnd(1) ! single random number
real(dp), intent(in) :: WindoW ! Size of window for bead selection

integer, parameter, dimension(0:3) :: changeBoth = [3, 2, 1, 0]
integer, parameter, dimension(0:3) :: changeFirst = [2, 3, 0, 1]
integer, parameter, dimension(0:3) :: changeSecond = [1,0, 3, 2]

!TOdo saving RP is not actually needed, even in these cases, but Brad's code assumes that we have RP.
if (WLC_P__RING .OR. WLC_P__INTERP_BEAD_LENNARD_JONES) then
    wlc_RP = wlc_R
    wlc_UP = wlc_U
endif

! Change wlc_AB (a.k.a HP1 binding type fore section of polymer)

call drawWindow(window,WLC_P__MAXWINDOW_SLIDE_MOVE,.false.,rand_stat,&
                IT1,IT2,IB1,IB2,IP,DIB,success)
if (success .eqv.  .false.) return

!keep binding constant within monomers
IT1 = IT1-MOD(IT1-1,WLC_P__NBPM)
IT2 = IT2-MOD(IT2-1,WLC_P__NBPM) + WLC_P__NBPM-1
IB1 = MOD(IT1,WLC_P__NB)
IB2 = MOD(IT2,WLC_P__NB)


if (WLC_P__TWO_TAIL) then
    call random_number(urnd,rand_stat)
    if (urnd(1)>WLC_P__PROBSINGLESWAP) then
        do J = IT1, IT2
            wlc_ABP(J) = changeBoth(wlc_AB(J))
        enddo
    elseif (urnd(1)<WLC_P__PROBSINGLESWAP/2.0_dp) then
        do J = IT1, IT2
            wlc_ABP(J) = changeFirst(wlc_AB(J))
        enddo
    else
        do J = IT1, IT2
            wlc_ABP(J) = changeSecond(wlc_AB(J))
        enddo
    endif
else
    do J = IT1,IT2
        wlc_ABP(J) = 1-wlc_AB(J)
    ENDdo
endif

! Copy over R,U,V to RP, UP, and VP
wlc_RP(:,IT1:IT2) = wlc_R(:,IT1:IT2)
wlc_UP(:,IT1:IT2) = wlc_U(:,IT1:IT2)
if (WLC_P__LOCAL_TWIST) wlc_VP(:,IT1:IT2) = wlc_V(:,IT1:IT2)
end subroutine
