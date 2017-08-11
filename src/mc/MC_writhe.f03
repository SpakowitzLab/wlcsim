  !Calculate the writhe between segment that was moved and segment that was not moved using the method of Klenin (2000)

subroutine WRITHE_MOVE(R,RP,IB1,IB2,IT1,IT2,IP,N,Wr)
  use params, only : dp, pi
  implicit none

  integer, intent(in) :: N                 ! Number of beads per polymer
  real(dp), intent(in) :: R(3,N),RP(3,N)  ! Position
  integer, intent(in) :: IB1,IB2,IT1,IT2   ! Bead indices for terminal points of segment moved
  integer, intent(in) :: IP                ! Test polymer index
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
  real(dp)  a0
  real(dp)  a1
  real(dp)  a2
  real(dp)  a3

  !Variables for writhe integral
  real(dp)  t1
  real(dp)  t2
  real(dp)  F1
  real(dp)  F2
  real(dp)  F3
  real(dp)  F4

  !Counter variables
  integer  I
  integer  J
  integer NBI,NBO,II,IO,IIP1,IOP1
  real(dp)  dWr
  real(dp), intent(out) :: Wr      ! Writhe


  if (IB2 >= IB1) then
     NBI = IB2-IB1
     NBO = N-NBI
  else
     NBO = IB1-IB2
     NBI = N-NBO
  ENDif


  Wr = 0.
  II = IB1

  do  I = 1,NBI-1
     IO = IB1
     if (II.gt.N) then
        II = 1
        IIP1 = II + 1
     elseif (II == N) then
        IIP1 = 1
     else
        IIP1 = II + 1
     ENDif

     do J = 1,NBI-1


        if (IO > N) then
           IO = 1
           IOP1 = IO + 1
        elseif (IO == N) then
           IOP1 = 1
        else
           IOP1 = IO + 1
        ENDif

        r1 = R(:,II)
        r2 = R(:,IO)
        r12 = r2-r1

        if (RP(1,Io) /= R(1,Io).OR.RP(2,Io) /= R(2,Io).OR.RP(3,Io) /= R(3,Io)) then
           ! PRinT *, "POSITIONS NOT EQUAL", IO,"of", NBO
        ENDif

        if (RP(1,IOP1) /= R(1,IOP1).OR.RP(2,IOP1) /= R(2,IOP1).OR.RP(3,IOP1) /= RP(3,IOP1)) then
           !  PRinT *, "POSITIONS NOT EQUAL AT IOP1", IOP1
        endif


        s2 = SQRT(SUM((R(:,IOP1)-R(:,IO))**2))
        s1 = SQRT(SUM((R(:,IIP1)-R(:,II))**2))
        e2 = (R(:,IOP1)-R(:,IO))/s2
        e1 = (R(:,IIP1)-R(:,II))/s1
        cosB = doT_PRODUCT(e1,e2)

        sin2B = 1-(cosB**2)

        e3(1) = e1(2)*e2(3)-e1(3)*e2(2)
        e3(2) = e1(3)*e2(1)-e1(1)*e2(3)
        e3(3) = e1(1)*e2(2)-e1(2)*e2(1)

        a1 = doT_PRODUCT(r12,e2*cosB-e1)/(sin2B)
        a2 = doT_PRODUCT(r12,e2-e1*cosB)/sin2B
        a0 = doT_PRODUCT(r12,e3)/sin2B



        t1 = a1 + s1
        t2 = a2 + s2
        F1 = -ATAN((t1*t2 + (a0**2)*cosB)/(a0*SQRT((t1**2 + t2**2-2.*t1*t2*cosB + (a0**2)*sin2B))))/(4.*PI)
        t1 = a1 + s1
        t2 = a2
        F2 = -ATAN((t1*t2 + (a0**2)*cosB)/(a0*SQRT((t1**2 + t2**2-2*t1*t2*cosB + (a0**2)*sin2B))))/(4.*PI)
        t1 = a1
        t2 = a2 + s2
        F3 = -ATAN((t1*t2 + (a0**2)*cosB)/(a0*SQRT((t1**2 + t2**2-2*t1*t2*cosB + (a0**2)*sin2B))))/(4.*PI)
        t1 = a1
        t2 = a2
        F4 = -ATAN((t1*t2 + (a0**2)*cosB)/(a0*SQRT((t1**2 + t2**2-2*t1*t2*cosB + (a0**2)*sin2B))))/(4.*PI)
        dWr = F1-F2-F3 + F4

        if (dWr /= dWr) then
           dWr = 0.
        endif

        Wr = Wr + dWr

        IO = IO + 1
     ENDdo
     II = II + 1
  ENDdo

  Wr = 2.*Wr

  RETURN
END subroutine WRITHE_MOVE










