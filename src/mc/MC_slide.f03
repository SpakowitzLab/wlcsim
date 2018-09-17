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
subroutine MC_slide(IB1,IB2,IT1,IT2,MCAMP,WindoW,rand_stat,dib,success)
! values from wlcsim_data
use params, only: wlc_RP, wlc_R, wlc_U, wlc_VP, wlc_V, wlc_UP, &
    wlc_nBend, wlc_bendPoints, wlc_pointsMoved, wlc_nPointsMoved, &
    wlc_ExplicitBindingPair

use mersenne_twister
use params, only: dp
use windowTools, only: drawWindow

implicit none
integer, intent(out) :: IB1   ! Test bead position 1
integer, intent(out) :: IT1   ! Index of test bead 1
integer, intent(out) :: IB2   ! Test bead position 2
integer, intent(out) :: IT2   ! Index of test bead 2
integer, intent(out) :: dib   ! number of beads moved by move (plus or minus a few)
logical, intent(out) :: success

integer IP    ! Test polymer
integer I,J  ! Test indices
integer ii
integer otherEnd
! Things for random number generator
type(random_stat), intent(inout) :: rand_stat  ! status of random number generator
real(dp) urand(3)  ! random vector
real(dp), intent(in) :: MCAMP ! Amplitude of random change
!integer, intent(in) :: WLC_P__WINTYPE
real(dp), intent(in) :: WindoW ! Size of window for bead selection
real(dp) DR(3)    ! Displacement for slide move


!TOdo saving RP is not actually needed, even in these cases, but Brad's code assumes that we have RP.
if (WLC_P__RING .OR. WLC_P__INTERP_BEAD_LENNARD_JONES) then
    wlc_RP = wlc_R
    wlc_UP = wlc_U
endif

!     Perform slide move (MCTYPE 2)
call drawWindow(window,WLC_P__MAXWINDOW_SLIDE_MOVE,.false.,rand_stat,&
                IT1,IT2,IB1,IB2,IP,DIB,success)

call random_number(urand,rand_stat)
DR(1) = MCAMP*(urand(1)-0.5_dp)
DR(2) = MCAMP*(urand(2)-0.5_dp)
DR(3) = MCAMP*(urand(3)-0.5_dp)

!  Which elastic segments change
wlc_nBend = 0
if (IB1>1) then
    wlc_nBend = wlc_nBend + 1
    wlc_bendPoints(wlc_nBend)=IT1-1
    I=IT1-1
    wlc_RP(:,I)=wlc_R(:,I)
    wlc_UP(:,I)=wlc_U(:,I)
    if (WLC_P__LOCAL_TWIST) wlc_VP(:,I) = wlc_V(:,I)
endif
if (IB2<WLC_P__NB) then
    wlc_nBend = wlc_nBend + 1
    wlc_bendPoints(wlc_nBend)=IT2
    I=IT2+1
    wlc_RP(:,I)=wlc_R(:,I)
    wlc_UP(:,I)=wlc_U(:,I)
    if (WLC_P__LOCAL_TWIST) wlc_VP(:,I) = wlc_V(:,I)
endif
! Move Explicity Bound points along with they are bound to
if (WLC_P__EXPLICIT_BINDING) then
    do ii=IT1,IT2
        otherEnd = wlc_ExplicitBindingPair(ii)
        if (otherEnd == -1) cycle

        if (otherEnd >= IT1 .and. otherEnd <= IT2) cycle

        ! Move Explicitly bound point along with section
        I=otherEnd
        wlc_RP(:,I)=wlc_R(:,I) + DR
        wlc_UP(:,I)=wlc_U(:,I)
        if (WLC_P__LOCAL_TWIST) wlc_VP(:,I) = wlc_V(:,I)
        wlc_nPointsMoved=wlc_nPointsMoved+1
        wlc_pointsMoved(wlc_nPointsMoved)=otherEnd
        ! Add ajacent points to RP and bendPoints
        if (otherEnd .ne. IT1-1 .and. mod(otherEnd,WLC_P__NB) .ne. 0 ) then
            wlc_nBend=wlc_nBend+1
            wlc_bendPoints(wlc_nBend)=otherEnd
            I=otherEnd+1
            wlc_RP(:,I)=wlc_R(:,I)
            wlc_UP(:,I)=wlc_U(:,I)
            if (WLC_P__LOCAL_TWIST) wlc_VP(:,I) = wlc_V(:,I)
        endif
        if (otherEnd .ne. IT2+1 .and. mod(otherEnd,WLC_P__NB) .ne. 1) then
            wlc_nBend=wlc_nBend+1
            wlc_bendPoints(wlc_nBend)=otherEnd-1
            I=otherEnd-1
            wlc_RP(:,I)=wlc_R(:,I)
            wlc_UP(:,I)=wlc_U(:,I)
            if (WLC_P__LOCAL_TWIST) wlc_VP(:,I) = wlc_V(:,I)
        endif
    enddo
endif

I = IT1
do  J = 0,DIB

    if (I == (WLC_P__NB*IP + 1).AND.WLC_P__RING) then
       I = WLC_P__NB*(IP-1) + 1
    endif

    wlc_RP(:,I) = wlc_R(:,I) + DR
    wlc_UP(:,I) = wlc_U(:,I)
    if (WLC_P__LOCAL_TWIST) wlc_VP(:,I) = wlc_V(:,I)
    wlc_nPointsMoved=wlc_nPointsMoved+1
    wlc_pointsMoved(wlc_nPointsMoved)=I
    I = I + 1

ENDdo
end subroutine
