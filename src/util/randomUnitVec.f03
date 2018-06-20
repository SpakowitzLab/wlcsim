!! ---------------------------------------------------------------
!
!    Generates random unit vector
!
!------------------------------------------------------------
subroutine randomUnitVec(U,rand_stat)

use mersenne_twister
use params, only: dp, pi

implicit none
type(random_stat), intent(inout) :: rand_stat  ! status of random number generator
real(dp), intent(out) :: U(3)
real(dp) z
real(dp) theta
real(dp) temp
real urand(2)

call random_number(urand,rand_stat)
theta = urand(1)*2.0_dp*PI
z = urand(2)*2.0_dp-1.0_dp
temp = sqrt(1.0_dp-z**2)
U(1) = temp*cos(theta)
U(2) = temp*sin(theta)
U(3) = z

! Alturnative algorithem
!call random_number(urand,rand_stat)
!ALPHA = 2.0_dp*PI*urand(1)
!BETA = acos(2.0_dp*urand(2)-1.0_dp)
!U(1) = sin(BETA)*cos(ALPHA)
!U(2) = sin(BETA)*sin(ALPHA)
!U(3) = cos(BETA)

! Alturnative algorithem
!real urand(3)
!call random_gauss(urand, rand_stat)
!U = urand
!U = U/norm2(U)
end subroutine
