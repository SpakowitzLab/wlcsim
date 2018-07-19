#include "../defines.inc"
!--------------------------------------------------------------*
!
!           Makes Monti Carlo Moves
!
!    Quinn split out this file on 8/9/17
!
!---------------------------------------------------------------

! variables that need to be allocated only on certain branches moved into MD to prevent segfaults
! please move other variables in as you see fit
subroutine MC_slide(wlc_p,wlc_d,IB1,IB2,IT1,IT2 &
                  ,MCAMP,WindoW,rand_stat &
                  ,dib,success)

use mersenne_twister
use params, only: dp,wlcsim_params, wlcsim_data
use windowTools, only: drawWindow

implicit none
type(wlcsim_params), intent(in) :: wlc_p
type(wlcsim_data), intent(inout) :: wlc_d
integer, intent(out) :: IB1   ! Test bead position 1
integer, intent(out) :: IT1   ! Index of test bead 1
integer, intent(out) :: IB2   ! Test bead position 2
integer, intent(out) :: IT2   ! Index of test bead 2
integer, intent(out) :: dib   ! number of beads moved by move (plus or minus a few)
logical, intent(out) :: success

integer IP    ! Test polymer
integer I,J  ! Test indices
! Things for random number generator
type(random_stat), intent(inout) :: rand_stat  ! status of random number generator
real urand(3)  ! random vector
real(dp), intent(in) :: MCAMP ! Amplitude of random change
!integer, intent(in) :: WLC_P__WINTYPE
real(dp), intent(in) :: WindoW ! Size of window for bead selection
real(dp) DR(3)    ! Displacement for slide move


!TOdo saving RP is not actually needed, even in these cases, but Brad's code assumes that we have RP.
if (WLC_P__RING .OR. WLC_P__INTERP_BEAD_LENNARD_JONES) then
    wlc_d%RP = wlc_d%R
    wlc_d%UP = wlc_d%U
endif

!     Perform slide move (MCTYPE 2)
call drawWindow(wlc_d,window,WLC_P__MAXWINDOW_SLIDE_MOVE,.true.,rand_stat,&
                IT1,IT2,IB1,IB2,IP,DIB,success)
if (success .eqv. .false.) return

call random_number(urand,rand_stat)
DR(1) = MCAMP*(urand(1)-0.5)
DR(2) = MCAMP*(urand(2)-0.5)
DR(3) = MCAMP*(urand(3)-0.5)

I = IT1
do  J = 0,DIB

    if (I == (WLC_P__NB*IP + 1).AND.WLC_P__RING) then
       I = WLC_P__NB*(IP-1) + 1
    endif

    wlc_d%RP(:,I) = wlc_d%R(:,I) + DR
    wlc_d%UP(:,I) = wlc_d%U(:,I)
    if (WLC_P__LOCAL_TWIST) wlc_d%VP(:,I) = wlc_d%V(:,I)
    I = I + 1

ENDdo
end subroutine
