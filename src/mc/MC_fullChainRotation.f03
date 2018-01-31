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
subroutine MC_fullChainRotation(wlc_p,R,U,RP,UP,IP,IB1,IB2,IT1,IT2 &
                  ,MCAMP,rand_stat)

use mersenne_twister
use params, only: dp, pi,wlcsim_params

!TODO: replace R,U,RP,UP .... with wlc_d

implicit none
type(wlcsim_params), intent(in) :: wlc_p
real(dp), intent(in) :: R(3,wlc_p%NT)  ! Bead positions
real(dp), intent(in) :: U(3,wlc_p%NT)  ! Tangent vectors
real(dp), intent(out) :: RP(3,wlc_p%NT)  ! Bead positions
real(dp), intent(out) :: UP(3,wlc_p%NT)  ! Tangent vectors
integer, intent(out) :: IP    ! Test polymer
integer, intent(out) :: IB1   ! Test bead position 1
integer, intent(out) :: IT1   ! Index of test bead 1
integer, intent(out) :: IB2   ! Test bead position 2
integer, intent(out) :: IT2   ! Index of test bead 2

integer I  ! Test indices
! Things for random number generator
type(random_stat), intent(inout) :: rand_stat  ! status of random number generator
real urand(3)  ! random vector
real urnd(1) ! single random number
integer irnd(1)
! Variables for the crank-shaft move

real(dp) TA(3)    ! Axis of rotation
real(dp) P1(3)    ! Point on rotation line
real(dp) ROT(4,4) ! Rotation matrix

real(dp) ALPHA    ! Angle of move
real(dp) BETA     ! Angle of move

real(dp), intent(in) :: MCAMP ! Amplitude of random change


!TOdo saving RP is not actually needed, even in these cases, but Brad's code assumes that we have RP.
if (WLC_P__RING .OR. WLC_P__INTERP_BEAD_LENNARD_JONES) then
    RP = R
    UP = U
    P1 = 0.0_dp
endif

!     Perform a full chain rotation

call random_index(WLC_P__NP,irnd,rand_stat)
IP=irnd(1)
IB1 = 1
IB2 = WLC_P__NB
IT1 = WLC_P__NB*(IP-1) + IB1
IT2 = WLC_P__NB*(IP-1) + IB2

call random_number(urand,rand_stat)
ALPHA = 2.0_dp*PI*urand(2)
BETA = acos(2.0_dp*urand(3)-1.0_dp)
TA(1) = sin(BETA)*cos(ALPHA)
TA(2) = sin(BETA)*sin(ALPHA)
TA(3) = cos(BETA)

! use ~central bead to put axes through
! you could also use center of mass if you wanted
P1(1) = R(1,(IT1 + IT2)/2)
P1(2) = R(2,(IT1 + IT2)/2)
P1(3) = R(3,(IT1 + IT2)/2)

call random_number(urnd,rand_stat)
ALPHA = MCAMP*(urnd(1)-0.5)

ROT(1,1) = TA(1)**2. + (TA(2)**2. + TA(3)**2.)*cos(ALPHA)
ROT(1,2) = TA(1)*TA(2)*(1.-cos(ALPHA))-TA(3)*sin(ALPHA)
ROT(1,3) = TA(1)*TA(3)*(1.-cos(ALPHA)) + TA(2)*sin(ALPHA)
ROT(1,4) = (P1(1)*(1.-TA(1)**2.) &
-TA(1)*(P1(2)*TA(2) + P1(3)*TA(3)))*(1.-cos(ALPHA)) + (P1(2)*TA(3)-P1(3)*TA(2))*sin(ALPHA)

ROT(2,1) = TA(1)*TA(2)*(1.-cos(ALPHA)) + TA(3)*sin(ALPHA)
ROT(2,2) = TA(2)**2. + (TA(1)**2. + TA(3)**2.)*cos(ALPHA)
ROT(2,3) = TA(2)*TA(3)*(1.-cos(ALPHA))-TA(1)*sin(ALPHA)
ROT(2,4) = (P1(2)*(1.-TA(2)**2.) &
-TA(2)*(P1(1)*TA(1) + P1(3)*TA(3)))*(1.-cos(ALPHA)) + (P1(3)*TA(1)-P1(1)*TA(3))*sin(ALPHA)

ROT(3,1) = TA(1)*TA(3)*(1.-cos(ALPHA))-TA(2)*sin(ALPHA)
ROT(3,2) = TA(2)*TA(3)*(1.-cos(ALPHA)) + TA(1)*sin(ALPHA)
ROT(3,3) = TA(3)**2. + (TA(1)**2. + TA(2)**2.)*cos(ALPHA)
ROT(3,4) = (P1(3)*(1.-TA(3)**2.) &
-TA(3)*(P1(1)*TA(1) + P1(2)*TA(2)))*(1.-cos(ALPHA)) + (P1(1)*TA(2)-P1(2)*TA(1))*sin(ALPHA)

do I = IT1,IT2
   RP(1,I) = ROT(1,4) + ROT(1,1)*R(1,I) + ROT(1,2)*R(2,I) + ROT(1,3)*R(3,I)
   RP(2,I) = ROT(2,4) + ROT(2,1)*R(1,I) + ROT(2,2)*R(2,I) + ROT(2,3)*R(3,I)
   RP(3,I) = ROT(3,4) + ROT(3,1)*R(1,I) + ROT(3,2)*R(2,I) + ROT(3,3)*R(3,I)
   UP(1,I) = ROT(1,1)*U(1,I) + ROT(1,2)*U(2,I) + ROT(1,3)*U(3,I)
   UP(2,I) = ROT(2,1)*U(1,I) + ROT(2,2)*U(2,I) + ROT(2,3)*U(3,I)
   UP(3,I) = ROT(3,1)*U(1,I) + ROT(3,2)*U(2,I) + ROT(3,3)*U(3,I)
enddo
end subroutine
