  !Calculate the Gauss integral between a pair of line segments (half the contribution to the writhe)
  !Based on methods by Klenin (2000)

subroutine GAUSSPAIR(R1,R2,e1,e2,s1,s2,dWr)
  use params, only : dp
  implicit none

  real(dp), PARAMETER :: PI = 3.141592654d0 ! Value of pi


  real(dp), intent(in) :: R1(3),R2(3)  ! PositionS
  real(dp), intent(in) :: e1(3),e2(3)  ! segment orientations
  real(dp), intent(in) :: s1,s2        ! segment lenghts
  !Geometric variables
  real(dp)  r12(3)                 ! Relative position vector
  real(dp)  beta                   ! Angle between tangents
  real(dp)  cosB
  real(dp)  sin2B
  real(dp)  e3(3)
  real(dp)  a0
  real(dp)  a1
  real(dp)  a2
  real(dp)  a3
  real(dp)  b
  !Variables for writhe integral
  real(dp)  t1
  real(dp)  t2
  real(dp)  F1
  real(dp)  F2
  real(dp)  F3
  real(dp)  F4

  !Counter variables
  integer  I,J
  real(dp), intent(out) :: dWr      ! Writhe

  r12 = R2-R1
  dWr = 0.

  cosB = doT_PRODUCT(e1,e2)
  B = ACOS(cosB)

  sin2B = sin(B)**2.
  e3(1) = e1(2)*e2(3)-e1(3)*e2(2)
  e3(2) = e1(3)*e2(1)-e1(1)*e2(3)
  e3(3) = e1(1)*e2(2)-e1(2)*e2(1)

  if (abs(sin2B).lt.(10.**(-15.)))then
     sin2b = 0.
  endif

  a1 = doT_PRODUCT(r12,e2*cosB-e1)/(sin2B)
  a2 = doT_PRODUCT(r12,e2-e1*cosB)/sin2B
  a0 = doT_PRODUCT(r12,e3)/sin2B

  t1 = a1 + s1
  t2 = a2 + s2
  F1 = -ATAN((t1*t2 + (a0**2.)*cosB)/(a0*SQRT((t1**2. + t2**2.-2.*t1*t2*cosB + (a0**2.)*sin2B))))/(4.*PI)
  t1 = a1 + s1
  t2 = a2
  F2 = -ATAN((t1*t2 + (a0**2.)*cosB)/(a0*SQRT((t1**2. + t2**2.-2.*t1*t2*cosB + (a0**2.)*sin2B))))/(4.*PI)
  t1 = a1
  t2 = a2 + s2
  F3 = -ATAN((t1*t2 + (a0**2.)*cosB)/(a0*SQRT((t1**2. + t2**2.-2.*t1*t2*cosB + (a0**2.)*sin2B))))/(4.*PI)
  t1 = a1
  t2 = a2
  F4 = -ATAN((t1*t2 + (a0**2.)*cosB)/(a0*SQRT((t1**2. + t2**2.-2.*t1*t2*cosB + (a0**2.)*sin2B))))/(4.*PI)
  dWr = F1-F2-F3 + F4

  !Gauss integral = 0 for segments in plane (dWr = NaN)
  if (dWr /= dWr) then
     dWr = 0.
  endif
  !Set dwr = 0 for segments nearly in-plane to avoid errors due to floating point
  if (abs(a0) < 10.**(-12.)) then
     dWr = 0.
  endif


  RETURN
END subroutine GAUSSPAIR


