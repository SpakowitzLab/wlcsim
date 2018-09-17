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
subroutine MC_pivot(IB1,IB2,IT1,IT2,MCAMP,WindoW,rand_stat,success)
! values from wlcsim_data
use params, only: wlc_RP, wlc_R, wlc_UP, wlc_U, wlc_V&
    , wlc_VP, wlc_bendPoints, wlc_nBend, wlc_nPointsMoved, wlc_pointsMoved

use mersenne_twister
use params, only: dp
use vector_utils, only: randomUnitVec, rotateR, rotateU, axisAngle
use windowTools, only: exponential_random_int, enforceBinding

implicit none
integer, intent(out) :: IB1   ! Test bead position 1
integer, intent(out) :: IT1   ! Index of test bead 1
integer, intent(out) :: IB2   ! Test bead position 2
integer, intent(out) :: IT2   ! Index of test bead 2
logical, intent(out) :: success

integer IP    ! Test polymer
integer I  ! Test indices
! Things for random number generator
type(random_stat), intent(inout) :: rand_stat  ! status of random number generator
real(dp) urnd(1) ! single random number
integer irnd(1)

real(dp) TA(3)    ! Axis of rotation
real(dp) P1(3)    ! Point on rotation line
real(dp) ROT(3,4) ! Rotation matrix
real(dp) ALPHA    ! Angle of move

!     MC adaptation variables

real(dp), intent(in) :: MCAMP ! Amplitude of random change
real(dp), intent(in) :: WindoW ! Size of window for bead selection

!TOdo saving RP is not actually needed, even in these cases, but Brad's code assumes that we have RP.
if (WLC_P__RING .OR. WLC_P__INTERP_BEAD_LENNARD_JONES) then
    wlc_RP = wlc_R
    wlc_UP = wlc_U
    P1 = 0.0_dp
endif

! We don't have to protect moves 4-10 with if WLC_P__RING because the code is identical in both cases
!     Perform pivot move (MCTYPE 3)
    call random_index(WLC_P__NP,irnd,rand_stat)
    IP=irnd(1)
    call random_index(WLC_P__NB,irnd,rand_stat)
    IB1=irnd(1)
    call random_number(urnd,rand_stat)
    if (urnd(1).gt.0.5_dp) then
        IB2 = exponential_random_int(window,rand_stat) + 1
        if (IB2 > WLC_P__NB) then
            IB2 = WLC_P__NB
        endif
        IB1 = 1
        IT1 = WLC_P__NB*(IP-1) + IB1
        IT2 = WLC_P__NB*(IP-1) + IB2
        P1 = wlc_R(:,IT2)
        if (IB2<WLC_P__NB) then
            wlc_nBend = wlc_nBend + 1
            wlc_bendPoints(wlc_nBend)=IT2
            I=IT2+1
            wlc_RP(:,I)=wlc_R(:,I)
            wlc_UP(:,I)=wlc_U(:,I)
            if (WLC_P__LOCAL_TWIST) wlc_VP(:,I) = wlc_V(:,I)
        endif
    else
        IB1 = WLC_P__NB-exponential_random_int(window,rand_stat)
        if (IB1 < 1) then
            IB1 = 1
        endif
        IB2 = WLC_P__NB
        IT1 = WLC_P__NB*(IP-1) + IB1
        IT2 = WLC_P__NB*(IP-1) + IB2
        P1 = wlc_R(:,IT1)
        if (IB1>1) then
            wlc_nBend = wlc_nBend + 1
            wlc_bendPoints(wlc_nBend)=IT1-1
            I=IT1-1
            wlc_RP(:,I)=wlc_R(:,I)
            wlc_UP(:,I)=wlc_U(:,I)
            if (WLC_P__LOCAL_TWIST) wlc_VP(:,I) = wlc_V(:,I)
        endif
    endif

   call randomUnitVec(TA,rand_stat)
   call random_number(urnd,rand_stat)
   ALPHA = MCAMP*(urnd(1)-0.5_dp)
   call axisAngle(ROT,alpha,TA,P1)

    if (WLC_P__EXPLICIT_BINDING) then
        call enforceBinding(rand_stat,IB1,IB2,IT1,IT2,WLC_P__MAXWINDOW_PIVOT_MOVE,success)
        if (success .eqv. .False.) return
    else
        success = .TRUE.
    endif

    do I=IT1,IT2
        wlc_RP(:,I) = rotateR(ROT,wlc_R(:,I))
        wlc_UP(:,I) = rotateU(ROT,wlc_U(:,I))
        if (WLC_P__LOCAL_TWIST) then
            wlc_VP(:,I) = rotateU(ROT,wlc_V(:,I))
        endif
        wlc_nPointsMoved=wlc_nPointsMoved+1
        wlc_pointsMoved(wlc_nPointsMoved)=I
    enddo
end subroutine
