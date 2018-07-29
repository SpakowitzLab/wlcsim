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
subroutine MC_fullChainRotation(IB1,IB2,IT1,IT2,MCAMP,rand_stat)
! values from wlcsim_data
use params, only: wlc_R, wlc_V, wlc_RP, wlc_VP, wlc_UP, wlc_U
use mersenne_twister
use params, only: dp
use vector_utils, only: rotateR, rotateU, axisAngle, randomUnitVec
implicit none
integer, intent(out) :: IB1   ! Test bead position 1
integer, intent(out) :: IT1   ! Index of test bead 1
integer, intent(out) :: IB2   ! Test bead position 2
integer, intent(out) :: IT2   ! Index of test bead 2

integer IP    ! Test polymer
integer I  ! Test indices
! Things for random number generator
type(random_stat), intent(inout) :: rand_stat  ! status of random number generator
real(dp) urnd(1) ! single random number
integer irnd(1)
! Variables for the crank-shaft move

real(dp) TA(3)    ! Axis of rotation
real(dp) P1(3)    ! Point on rotation line
real(dp) ROT(3,4) ! Rotation matrix

real(dp) ALPHA    ! Angle of move
real(dp), intent(in) :: MCAMP ! Amplitude of random change


!TOdo saving RP is not actually needed, even in these cases, but Brad's code assumes that we have RP.
if (WLC_P__RING .OR. WLC_P__INTERP_BEAD_LENNARD_JONES) then
    wlc_RP = wlc_R
    wlc_UP = wlc_U
    P1 = 0.0_dp
endif

!     Perform a full chain rotation

call random_index(WLC_P__NP,irnd,rand_stat)
IP=irnd(1)
IB1 = 1
IB2 = WLC_P__NB
IT1 = WLC_P__NB*(IP-1) + IB1
IT2 = WLC_P__NB*(IP-1) + IB2

call randomUnitVec(TA,rand_stat)

! use ~central bead to put axes through
! you could also use center of mass if you wanted
P1 = wlc_R(:,(IT1 + IT2)/2)

call random_number(urnd,rand_stat)
ALPHA = MCAMP*(urnd(1)-0.5)

call axisAngle(ROT,alpha,TA,P1)
do I = IT1,IT2
    wlc_RP(:,I) = rotateR(ROT,wlc_R(:,I))
    wlc_UP(:,I) = rotateU(ROT,wlc_U(:,I))
    if (WLC_P__LOCAL_TWIST) wlc_VP(:,I) = rotateU(ROT,wlc_V(:,I))
enddo
end subroutine
