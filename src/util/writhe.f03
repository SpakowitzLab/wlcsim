  !Calculate the Writhe of a curve using the method of Klenin (2000)

subroutine WRITHE(R,N,Wr)
  use params, only : dp, pi

  implicit none


  integer, intent(in) :: N                 ! Number of beads
  real(dp), intent(in) :: R(3,N)  ! PositionS
    !Geometric variables
  real(dp)  r1(3)                  ! Position bead 1
  real(dp)  r2(3)                  ! Position bead 2
  real(dp)  r12(3)                 ! Relative position vector
  real(dp)  s1                     ! Length of segment 1
  real(dp)  s2                     ! Length of segment 2
  real(dp)  beta                   ! Angle between tangents
  real(dp)  e1(3)                  ! Tangent of first segment
  real(dp)  e2(3)                  ! Tangent of second segment
  real(dp)  cosB
  real(dp)  sin2B
  real(dp)  e3(3)
  real(dp)  OMEGA(N,N)             ! Matrix of solid angles
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
  integer  I,IP1
  integer  J,JP1

  real(dp)  dWr
  real(dp), intent(out) :: Wr      ! Writhe


  OMEGA = 0
  Wr = 0.

  do  I = 2,N
      if (I == N) then
           IP1 = 1
        else
           IP1 = I + 1
        ENDif
     do J = 1,I-1



        if (J == N) then
           JP1 = 1
        else
           JP1 = J + 1
        ENDif



        r1 = R(:,I)
        r2 = R(:,J)
        r12 = r2-r1

        s2 = SQRT(SUM((R(:,JP1)-R(:,J))**2))
        s1 = SQRT(SUM((R(:,IP1)-R(:,I))**2))
        e2 = (R(:,JP1)-R(:,J))/s2
        e1 = (R(:,IP1)-R(:,I))/s1
        cosB = doT_PRODUCT(e1,e2)
        B = ACOS(cosB)

       ! sin2B = 1.-(cosB**2.)
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

           if (dWr /= dWr) then
              dWr = 0.
           endif
           if (JP1 == I.OR.IP1 == J) then
              dWr = 0.
           endif

           if (abs(a0) < 10.**(-10.)) then
              dWr = 0.
           endif


           Wr = Wr + dWr


     ENDdo
  ENDdo

  Wr = 2.*Wr

  RETURN
END subroutine WRITHE










