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
subroutine MC_rotate(wlc_p,R,U,RP,UP,IB1,IB2,IT1,IT2 &
                  ,MCAMP,rand_stat)

use mersenne_twister
use params, only: dp, pi,wlcsim_params
use vector_utils, only: randomUnitVec

!TODO: replace R,U,RP,UP .... with wlc_d

implicit none
type(wlcsim_params), intent(in) :: wlc_p
real(dp), intent(in) :: R(3,WLC_P__NT)  ! Bead positions
real(dp), intent(in) :: U(3,WLC_P__NT)  ! Tangent vectors
real(dp), intent(out) :: RP(3,WLC_P__NT)  ! Bead positions
real(dp), intent(out) :: UP(3,WLC_P__NT)  ! Tangent vectors
integer, intent(out) :: IB1   ! Test bead position 1
integer, intent(out) :: IT1   ! Index of test bead 1
integer, intent(out) :: IB2   ! Test bead position 2
integer, intent(out) :: IT2   ! Index of test bead 2

integer IP    ! Test polymer
integer I  ! Test indices
! Things for random number generator
type(random_stat), intent(inout) :: rand_stat  ! status of random number generator
integer irnd(1)  ! random vector
! Variables for the crank-shaft move

real(dp) TA(3)    ! Axis of rotation
real(dp) ROT(3,4) ! Rotation matrix
real(dp) ALPHA    ! Angle of move
real urnd(1) ! single random number
real(dp), parameter, dimension(3) ::  P1 = [0.0_dp, 0.0_dp, 0.0_dp]

!     MC adaptation variables

real(dp), intent(in) :: MCAMP ! Amplitude of random change


!TOdo saving RP is not actually needed, even in these cases, but Brad's code assumes that we have RP.
if (WLC_P__RING .OR. WLC_P__INTERP_BEAD_LENNARD_JONES) then
    RP = R
    UP = U
endif

!     Perform rotate move (MCTYPE 4)
!     a.k.a. rotate a single bead
call random_index(WLC_P__NP,irnd,rand_stat)
IP=irnd(1)
call random_index(WLC_P__NB,irnd,rand_stat)
IB1=irnd(1)
IB2 = IB1
IT1 = WLC_P__NB*(IP-1) + IB1
IT2 = WLC_P__NB*(IP-1) + IB2

call randomUnitVec(TA,rand_stat)
call random_number(urnd,rand_stat)
ALPHA = MCAMP*(urnd(1)-0.5_dp)
call axisAngle(ROT,alpha,TA,P1)

I = IT1
UP(1,I) = ROT(1,1)*U(1,I) + ROT(1,2)*U(2,I) + ROT(1,3)*U(3,I)
UP(2,I) = ROT(2,1)*U(1,I) + ROT(2,2)*U(2,I) + ROT(2,3)*U(3,I)
UP(3,I) = ROT(3,1)*U(1,I) + ROT(3,2)*U(2,I) + ROT(3,3)*U(3,I)
RP(1,I) = R(1,I)
RP(2,I) = R(2,I)
RP(3,I) = R(3,I)
end subroutine
