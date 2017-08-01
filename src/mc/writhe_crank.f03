  !Calculate the portion of the writhe that changes during a crankshaft move

subroutine WRITHECRANK(R,IT1,IT2,N,Wr)
  use params, only : dp
  implicit none

  real(dp), PARAMETER :: PI = 3.141592654d0 ! Value of pi


  integer, intent(in) :: N                 ! Number of beads
  real(dp), intent(in) :: R(3,N)  ! Positions
  !Geometric variables
  real(dp)  r1(3)                  ! Position bead 1
  real(dp)  r2(3)                  ! Position bead 2
  real(dp)  r12(3)                 ! Relative position vector
  real(dp)  s1                     ! Length of segment 1
  real(dp)  s2                     ! Length of segment 2
  real(dp)  beta                   ! Angle between tangents
  real(dp)  e1(3)                  ! Tangent of first segment
  real(dp)  e2(3)                  ! Tangent of second segment
  real(dp)  e3(3)
  integer  DIB,DIO                         !Number of segments inside section rotated,outside


  !Counter variables
  integer IT1,IT2                         !Indices of endpoints in segment rotated
  integer II,IO                           !Index of bead inside segment rotated,outside
  integer  I,IP1
  integer  J,JP1

  real(dp)  dWr
  real(dp), intent(out) :: Wr      ! Writhe



  Wr = 0.
  if (IT2 >= IT1) then
     DIB = IT2-IT1
  else
     DIB = (N-IT1) + IT2
  ENDif

  DIO = N-DIB
  II = IT1



  do  I = 1,DIB

     if (II == N) then
        IP1 = 1
     elseif (II == N + 1) then
        II = 1
        IP1 = II + 1
     else
        IP1 = II + 1
     ENDif
     IO = IT2
     do J = 1,DIO



        if (IO == N) then
           JP1 = 1
        elseif (IO == N + 1) then
           IO = 1
           JP1 = IO + 1
        else
           JP1 = IO + 1
        ENDif



        r1 = R(:,II)
        r2 = R(:,IO)
        r12 = r2-r1

        s2 = SQRT(SUM((R(:,JP1)-R(:,IO))**2))
        s1 = SQRT(SUM((R(:,IP1)-R(:,II))**2))
        e2 = (R(:,JP1)-R(:,IO))/s2
        e1 = (R(:,IP1)-R(:,II))/s1

        CALL GAUSSPAIR(R1,R2,e1,e2,s1,s2,dWr)
        !Wr = 0 for segments in plane (Wr = NaN)
        if (dWr /= dWr) then
           dWr = 0.
        endif
        !Wr = 0 for adjacent segments
        if (JP1 == II.OR.IP1 == IO) then
           dWr = 0.
        endif


        Wr = Wr + dWr

        IO = IO + 1

     ENDdo
     II = II + 1
  ENDdo


  Wr = 2.*Wr

  RETURN
END subroutine WRITHECRANK










