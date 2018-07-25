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
subroutine MC_pivot(wlc_d,IB1,IB2,IT1,IT2,MCAMP,WindoW,rand_stat,success)

use mersenne_twister
use params, only: dp, wlcsim_data
use vector_utils, only: randomUnitVec, rotateR, rotateU, axisAngle
use windowTools, only: exponential_random_int, enforceBinding

implicit none
type(wlcsim_data), intent(inout) :: wlc_d
integer, intent(out) :: IB1   ! Test bead position 1
integer, intent(out) :: IT1   ! Index of test bead 1
integer, intent(out) :: IB2   ! Test bead position 2
integer, intent(out) :: IT2   ! Index of test bead 2
logical, intent(out) :: success

integer IP    ! Test polymer
integer I  ! Test indices
! Things for random number generator
type(random_stat), intent(inout) :: rand_stat  ! status of random number generator
real urnd(1) ! single random number
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
    wlc_d%RP = wlc_d%R
    wlc_d%UP = wlc_d%U
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
        P1 = wlc_d%R(:,IT2)
    else
        IB1 = WLC_P__NB-exponential_random_int(window,rand_stat)
        if (IB1 < 1) then
            IB1 = 1
        endif
        IB2 = WLC_P__NB
        IT1 = WLC_P__NB*(IP-1) + IB1
        IT2 = WLC_P__NB*(IP-1) + IB2
        P1 = wlc_d%R(:,IT1)
    endif

   call randomUnitVec(TA,rand_stat)
   call random_number(urnd,rand_stat)
   ALPHA = MCAMP*(urnd(1)-0.5_dp)
   call axisAngle(ROT,alpha,TA,P1)

    if (WLC_P__EXPLICIT_BINDING) then
        call enforceBinding(rand_stat,IB1,IB2,IT1,IT2,wlc_d,WLC_P__MAXWINDOW_PIVOT_MOVE,success)
    else
        success = .TRUE.
    endif

    do I=IT1,IT2
        wlc_d%RP(:,I) = rotateR(ROT,wlc_d%R(:,I))
        wlc_d%UP(:,I) = rotateU(ROT,wlc_d%U(:,I))
        if (WLC_P__LOCAL_TWIST) then
            wlc_d%VP(:,I) = rotateU(ROT,wlc_d%V(:,I))
        endif
    enddo
end subroutine
