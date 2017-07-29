  !Calculate the portion of the writhe that changes during a segment slide move

subroutine WRITHESLIDE(R,IT1,IT2,N,Wr)
  use params, only : dp
  implicit none

  real(dp), PARAMETER :: PI = 3.141592654d0 ! Value of pi


  integer, intent(in) :: N                 ! Number of beads
  real(dp), intent(in) :: R(N,3)  ! PositionS
  !Geometric variables
  real(dp)  r1(3)                  ! Position bead 1
  real(dp)  r2(3)                  ! Position bead 2
  real(dp)  r12(3)                 ! Relative position vector
  real(dp)  s1                     ! Length of segment 1
  real(dp)  s2                     ! Length of segment 2
  real(dp)  e1(3)                  ! Tangent of first segment
  real(dp)  e2(3)                  ! Tangent of second segment
  real(dp)  e3(3)
  integer  DIB                             !Number of segments inside segment moved
  integer DIO                              !Number of segments unchanged by move

  !Counter variables
  integer IT1,IT2,II,IO,IS1,IS2
  integer  I,IP1,IS1P1,IS2P1
  integer  J,JP1

  real(dp)  dWr
  real(dp), intent(out) :: Wr      ! Writhe



  Wr = 0.
  if (IT2 >= IT1) then
     DIB = IT2-IT1
  else
     DIB = (N-IT1) + IT2
  ENDif
  DIO = N-DIB-2
  IS1P1 = IT1                                ! Index of first stretched (or compressed) segment (not inside beads moved)
  IS2 = IT2                                  ! Index of second stretched (or compressed) segment (not inside beads moved)

  if (IS1P1 == 1) then
     IS1 = N
  else
     IS1 = IS1P1-1
  ENDif

  if (IS2 == N) then
     IS2P1 = 1
  else
     IS2P1 = IS2 + 1
  ENDif

  II = IT1


  if (IT1 /= IT2) then
     !Sum over segments inside segment slid (inner segments)
     do  I = 1,DIB

        if (II == N) then
           IP1 = 1
        elseif (II == N + 1) then
           II = 1
           IP1 = II + 1
        else
           IP1 = II + 1
        ENDif
        IO = IT2 + 1
        !Sum over segments unchanged by slide (outer segments)
        do J = 1,DIO

           if (IO == N) then
              JP1 = 1
           elseif (IO == N + 1) then
              IO = 1
              JP1 = IO + 1
           else
              JP1 = IO + 1
           ENDif


           !Calculate Gauss integeral between outer segment and inner segment

           r2 = R(IO,:)
           r1 = R(II,:)
           r12 = r2-r1

           s2 = SQRT(SUM((R(JP1,:)-R(IO,:))**2))
           s1 = SQRT(SUM((R(IP1,:)-R(II,:))**2))
           e2 = (R(JP1,:)-R(IO,:))/s2
           e1 = (R(IP1,:)-R(II,:))/s1

           CALL GAUSSPAIR(R1,R2,e1,e2,s1,s2,dWr)
           !Writhe = 0 between segments in plane (dWr = NaN)
           if (dWr /= dWr) then
              dWr = 0.
           endif
           !Writhe = 0 between adjacent segments
           if (JP1 == II.OR.IP1 == IO) then
              dWr = 0.
           endif

           Wr = Wr + dWr

           if (I == 1) then

              !Calculate Gauss integeral between outer segment and first stretched segment

              r2 = R(IO,:)
              r1 = R(IS1,:)

              r12 = r2-r1

              s2 = SQRT(SUM((R(JP1,:)-R(IO,:))**2))
              s1 = SQRT(SUM((R(IS1P1,:)-R(IS1,:))**2))
              e2 = (R(JP1,:)-R(IO,:))/s2
              e1 = (R(IS1P1,:)-R(IS1,:))/s1

              CALL GAUSSPAIR(R1,R2,e1,e2,s1,s2,dWr)
              !Writhe = 0 between segments in plane (dWr = NaN)
              if (dWr /= dWr) then
                 dWr = 0.
              endif
              !Writhe = 0 between adjacent segments
              if (JP1 == IS1.OR.IS1P1 == IO) then
                 dWr = 0.
              endif

              Wr = Wr + dWr

              !Calculate Gauss integeral between outer segment and second stretched segment
              r2 = R(IO,:)
              r1 = R(IS2,:)

              r12 = r2-r1

              s2 = SQRT(SUM((R(JP1,:)-R(IO,:))**2))
              s1 = SQRT(SUM((R(IS2P1,:)-R(IS2,:))**2))
              e2 = (R(JP1,:)-R(IO,:))/s2
              e1 = (R(IS2P1,:)-R(IS2,:))/s1

              CALL GAUSSPAIR(R1,R2,e1,e2,s1,s2,dWr)
              !Writhe = 0 between segments in plane (dWr = NaN)
              if (dWr /= dWr) then
                 dWr = 0.
              endif
              !Writhe = 0 between adjacent segments
              if (JP1 == IS2.OR.IS2P1 == IO) then
                 dWr = 0.
              endif
              Wr = Wr + dWr
           ENDif



           IO = IO + 1

        ENDdo

        !Calculate Gauss integral between inner segment and first stretched segment

        r2 = R(II,:)
        r1 = R(IS1,:)

        r12 = r2-r1

        s2 = SQRT(SUM((R(IP1,:)-R(II,:))**2))
        s1 = SQRT(SUM((R(IS1P1,:)-R(IS1,:))**2))
        e2 = (R(IP1,:)-R(II,:))/s2
        e1 = (R(IS1P1,:)-R(IS1,:))/s1

        CALL GAUSSPAIR(R1,R2,e1,e2,s1,s2,dWr)
        !Writhe = 0 between segments in plane (dWr = NaN)
        if (dWr /= dWr) then
           dWr = 0.
        endif
        !Writhe = 0 between adjacent segments
        if (IP1 == IS1.OR.IS1P1 == II) then
           dWr = 0.
        endif

        Wr = Wr + dWr

        !Calculate Gauss integral between inner segment and second stretched segment
        r2 = R(II,:)
        r1 = R(IS2,:)

        r12 = r2-r1

        s2 = SQRT(SUM((R(IP1,:)-R(II,:))**2))
        s1 = SQRT(SUM((R(IS2P1,:)-R(IS2,:))**2))
        e2 = (R(IP1,:)-R(II,:))/s2
        e1 = (R(IS2P1,:)-R(IS2,:))/s1

        CALL GAUSSPAIR(R1,R2,e1,e2,s1,s2,dWr)
        !Writhe = 0 between segments in plane (dWr = NaN)
        if (dWr /= dWr) then
           dWr = 0.
        endif
        !Writhe = 0 between adjacent segments
        if (IP1 == IS2.OR.IS2P1 == II) then
           dWr = 0.
        endif
        Wr = Wr + dWr

        II = II + 1
     ENDdo
  ENDif

  ! calculate Gauss integral  between two stretched segments
  r2 = R(IS2,:)
  r1 = R(IS1,:)

  r12 = r2-r1

  s2 = SQRT(SUM((R(IS2P1,:)-R(IS2,:))**2))
  s1 = SQRT(SUM((R(IS1P1,:)-R(IS1,:))**2))
  e2 = (R(IS2P1,:)-R(IS2,:))/s2
  e1 = (R(IS1P1,:)-R(IS1,:))/s1

  CALL  GAUSSPAIR(R1,R2,e1,e2,s1,s2,dWr)
  !Writhe = 0 between segments in plane (dWr = NaN)
  if (dWr /= dWr) then
     dWr = 0.
  endif
  !Writhe = 0 between adjacent segments
  if (IS2P1 == IS1.OR.IS1P1 == IS2) then
     dWr = 0.
  endif
  WR = WR + DWR




  if (IT1 == IT2) then
     !If IT1 == IT2, Get writhe between segment unaffected by move and two segments affected

     IO = IT1 + 1
     !Sum over segments unchanged by slide (outer segments)
     do J = 1,DIO

        if (IO == N) then
           JP1 = 1
        elseif (IO == N + 1) then
           IO = 1
           JP1 = IO + 1
        else
           JP1 = IO + 1
        ENDif


        !Calculate Gauss integeral between outer segment and first stretched segment

        r2 = R(IO,:)
        r1 = R(IS1,:)

        r12 = r2-r1

        s2 = SQRT(SUM((R(JP1,:)-R(IO,:))**2))
        s1 = SQRT(SUM((R(IS1P1,:)-R(IS1,:))**2))
        e2 = (R(JP1,:)-R(IO,:))/s2
        e1 = (R(IS1P1,:)-R(IS1,:))/s1

        CALL GAUSSPAIR(R1,R2,e1,e2,s1,s2,dWr)
        !Writhe = 0 between segments in plane (dWr = NaN)
        if (dWr /= dWr) then
           dWr = 0.
        endif
        !Writhe = 0 between adjacent segments
        if (JP1 == IS1.OR.IS1P1 == IO) then
           dWr = 0.
        endif


        Wr = Wr + dWr

        !Calculate Gauss integeral between outer segment and second stretched segment
        r2 = R(IO,:)
        r1 = R(IS2,:)

        r12 = r2-r1

        s2 = SQRT(SUM((R(JP1,:)-R(IO,:))**2))
        s1 = SQRT(SUM((R(IS2P1,:)-R(IS2,:))**2))
        e2 = (R(JP1,:)-R(IO,:))/s2
        e1 = (R(IS2P1,:)-R(IS2,:))/s1

        CALL GAUSSPAIR(R1,R2,e1,e2,s1,s2,dWr)
        !Writhe = 0 between segments in plane (dWr = NaN)
        if (dWr /= dWr) then
           dWr = 0.
        endif
        !Writhe = 0 between adjacent segments
        if (JP1 == IS2.OR.IS2P1 == IO) then
           dWr = 0.
        endif

        Wr = Wr + dWr
        IO = IO + 1
     ENDdo
  ENDif



  Wr = 2.*Wr

  RETURN
END subroutine WRITHESLIDE










