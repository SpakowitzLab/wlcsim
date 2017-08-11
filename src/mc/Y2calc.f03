!---------------------------------------------------------------!
!
!     Written by Quinn MacPherson 8/4/17
!---------------------------------------------------------------!
subroutine Y2calc(u,Y2real)
use params, only: dp, pi
implicit none

!real(dp), parameter :: pi=3.14159265359_dp

real(dp), intent(in) :: u(3)
real(dp), intent(out) :: Y2real(5)
real(dp), parameter, dimension(5) :: coefficient = (/0.5_dp*sqrt(15.0_dp/pi), &
                                     0.5_dp*sqrt(15.0_dp/pi), &
                                     0.25_dp*sqrt(5.0_dp/pi), &
                                     0.5_dp*sqrt(15.0_dp/pi), &
                                     0.25_dp*sqrt(15.0_dp/pi)/)
Y2real(1)=coefficient(1)*u(1)*u(2)
Y2real(2)=coefficient(2)*u(2)*u(3)
Y2real(3)=coefficient(3)*(2*u(3)**2 - u(1)**2 - u(2)**2)
Y2real(4)=coefficient(4)*u(3)*u(1)
Y2real(5)=coefficient(5)*(u(1)**2 - u(2)**2)
RETURN
END

!---------------------------------------------------------------!
