  !Calculate the portion of the writhe that changes during a segment slide move

SUBROUTINE WRITHESLIDE(R,IT1,IT2,N,Wr)
  IMPLICIT NONE

  DOUBLE PRECISION, PARAMETER :: PI=3.141592654d0 ! Value of pi


  INTEGER, INTENT(IN) :: N                 ! Number of beads
  DOUBLE PRECISION, INTENT(IN) :: R(N,3)  ! PositionS
  !Geometric variables
  DOUBLE PRECISION  r1(3)                  ! Position bead 1
  DOUBLE PRECISION  r2(3)                  ! Position bead 2
  DOUBLE PRECISION  r12(3)                 ! Relative position vector
  DOUBLE PRECISION  s1                     ! Length of segment 1
  DOUBLE PRECISION  s2                     ! Length of segment 2
  DOUBLE PRECISION  e1(3)                  ! Tangent of first segment
  DOUBLE PRECISION  e2(3)                  ! Tangent of second segment
  DOUBLE PRECISION  e3(3)
  INTEGER  DIB                             !Number of segments inside segment moved
  INTEGER DIO                              !Number of segments unchanged by move

  !Counter variables
  INTEGER IT1,IT2,II,IO,IS1,IS2
  INTEGER  I,IP1,IS1P1,IS2P1
  INTEGER  J,JP1

  DOUBLE PRECISION  dWr
  DOUBLE PRECISION, INTENT(OUT) :: Wr      ! Writhe



  Wr=0.
  IF (IT2.GE.IT1) THEN
     DIB=IT2-IT1
  ELSE
     DIB=(N-IT1)+IT2
  ENDIF
  DIO=N-DIB-2
  IS1P1=IT1                                ! Index of first stretched (or compressed) segment (not inside beads moved)
  IS2=IT2                                  ! Index of second stretched (or compressed) segment (not inside beads moved)

  IF (IS1P1.EQ.1) THEN
     IS1=N
  ELSE
     IS1=IS1P1-1
  ENDIF

  IF (IS2.EQ.N) THEN
     IS2P1=1
  ELSE
     IS2P1=IS2+1
  ENDIF

  II=IT1


  IF (IT1.NE.IT2) then
     !Sum over segments inside segment slid (inner segments)
     DO  I=1,DIB

        IF (II.EQ.N) THEN
           IP1=1
        ELSEIF (II.EQ.N+1) THEN
           II=1
           IP1=II+1
        ELSE
           IP1=II+1
        ENDIF
        IO=IT2+1
        !Sum over segments unchanged by slide (outer segments)
        DO J=1,DIO

           IF (IO.EQ.N) THEN
              JP1=1
           ELSEIF (IO.EQ.N+1) THEN
              IO=1
              JP1=IO+1
           ELSE
              JP1=IO+1
           ENDIF


           !Calculate Gauss integeral between outer segment and inner segment

           r2=R(IO,:)
           r1=R(II,:)
           r12=r2-r1

           s2=SQRT(SUM((R(JP1,:)-R(IO,:))**2))
           s1=SQRT(SUM((R(IP1,:)-R(II,:))**2))
           e2=(R(JP1,:)-R(IO,:))/s2
           e1=(R(IP1,:)-R(II,:))/s1

           CALL GAUSSPAIR(R1,R2,e1,e2,s1,s2,dWr)
           !Writhe=0 between segments in plane (dWr=NaN)
           if (dWr.NE.dWr) then
              dWr=0.
           endif
           !Writhe=0 between adjacent segments
           if (JP1.EQ.II.OR.IP1.EQ.IO) THEN
              dWr=0.
           endif

           Wr=Wr+dWr

           IF (I.EQ.1) THEN

              !Calculate Gauss integeral between outer segment and first stretched segment

              r2=R(IO,:)
              r1=R(IS1,:)

              r12=r2-r1

              s2=SQRT(SUM((R(JP1,:)-R(IO,:))**2))
              s1=SQRT(SUM((R(IS1P1,:)-R(IS1,:))**2))
              e2=(R(JP1,:)-R(IO,:))/s2
              e1=(R(IS1P1,:)-R(IS1,:))/s1

              CALL GAUSSPAIR(R1,R2,e1,e2,s1,s2,dWr)
              !Writhe=0 between segments in plane (dWr=NaN)
              if (dWr.NE.dWr) then
                 dWr=0.
              endif
              !Writhe=0 between adjacent segments
              if (JP1.EQ.IS1.OR.IS1P1.EQ.IO) THEN
                 dWr=0.
              endif

              Wr=Wr+dWr

              !Calculate Gauss integeral between outer segment and second stretched segment
              r2=R(IO,:)
              r1=R(IS2,:)

              r12=r2-r1

              s2=SQRT(SUM((R(JP1,:)-R(IO,:))**2))
              s1=SQRT(SUM((R(IS2P1,:)-R(IS2,:))**2))
              e2=(R(JP1,:)-R(IO,:))/s2
              e1=(R(IS2P1,:)-R(IS2,:))/s1

              CALL GAUSSPAIR(R1,R2,e1,e2,s1,s2,dWr)
              !Writhe=0 between segments in plane (dWr=NaN)
              if (dWr.NE.dWr) then
                 dWr=0.
              endif
              !Writhe=0 between adjacent segments
              if (JP1.EQ.IS2.OR.IS2P1.EQ.IO) THEN
                 dWr=0.
              endif
              Wr=Wr+dWr
           ENDIF



           IO=IO+1

        ENDDO

        !Calculate Gauss integral between inner segment and first stretched segment

        r2=R(II,:)
        r1=R(IS1,:)

        r12=r2-r1

        s2=SQRT(SUM((R(IP1,:)-R(II,:))**2))
        s1=SQRT(SUM((R(IS1P1,:)-R(IS1,:))**2))
        e2=(R(IP1,:)-R(II,:))/s2
        e1=(R(IS1P1,:)-R(IS1,:))/s1

        CALL GAUSSPAIR(R1,R2,e1,e2,s1,s2,dWr)
        !Writhe=0 between segments in plane (dWr=NaN)
        if (dWr.NE.dWr) then
           dWr=0.
        endif
        !Writhe=0 between adjacent segments
        if (IP1.EQ.IS1.OR.IS1P1.EQ.II) THEN
           dWr=0.
        endif

        Wr=Wr+dWr

        !Calculate Gauss integral between inner segment and second stretched segment
        r2=R(II,:)
        r1=R(IS2,:)

        r12=r2-r1

        s2=SQRT(SUM((R(IP1,:)-R(II,:))**2))
        s1=SQRT(SUM((R(IS2P1,:)-R(IS2,:))**2))
        e2=(R(IP1,:)-R(II,:))/s2
        e1=(R(IS2P1,:)-R(IS2,:))/s1

        CALL GAUSSPAIR(R1,R2,e1,e2,s1,s2,dWr)
        !Writhe=0 between segments in plane (dWr=NaN)
        if (dWr.NE.dWr) then
           dWr=0.
        endif
        !Writhe=0 between adjacent segments
        if (IP1.EQ.IS2.OR.IS2P1.EQ.II) THEN
           dWr=0.
        endif
        Wr=Wr+dWr

        II=II+1
     ENDDO
  ENDIF

  ! calculate Gauss integral  between two stretched segments
  r2=R(IS2,:)
  r1=R(IS1,:)

  r12=r2-r1

  s2=SQRT(SUM((R(IS2P1,:)-R(IS2,:))**2))
  s1=SQRT(SUM((R(IS1P1,:)-R(IS1,:))**2))
  e2=(R(IS2P1,:)-R(IS2,:))/s2
  e1=(R(IS1P1,:)-R(IS1,:))/s1

  CALL  GAUSSPAIR(R1,R2,e1,e2,s1,s2,dWr)
  !Writhe=0 between segments in plane (dWr=NaN)
  if (dWr.NE.dWr) then
     dWr=0.
  endif
  !Writhe=0 between adjacent segments
  if (IS2P1.EQ.IS1.OR.IS1P1.EQ.IS2) THEN
     dWr=0.
  endif
  WR=WR+DWR




  IF (IT1.EQ.IT2) THEN
     !If IT1.EQ.IT2, Get writhe between segment unaffected by move and two segments affected

     IO=IT1+1
     !Sum over segments unchanged by slide (outer segments)
     DO J=1,DIO

        IF (IO.EQ.N) THEN
           JP1=1
        ELSEIF (IO.EQ.N+1) THEN
           IO=1
           JP1=IO+1
        ELSE
           JP1=IO+1
        ENDIF


        !Calculate Gauss integeral between outer segment and first stretched segment

        r2=R(IO,:)
        r1=R(IS1,:)

        r12=r2-r1

        s2=SQRT(SUM((R(JP1,:)-R(IO,:))**2))
        s1=SQRT(SUM((R(IS1P1,:)-R(IS1,:))**2))
        e2=(R(JP1,:)-R(IO,:))/s2
        e1=(R(IS1P1,:)-R(IS1,:))/s1

        CALL GAUSSPAIR(R1,R2,e1,e2,s1,s2,dWr)
        !Writhe=0 between segments in plane (dWr=NaN)
        if (dWr.NE.dWr) then
           dWr=0.
        endif
        !Writhe=0 between adjacent segments
        if (JP1.EQ.IS1.OR.IS1P1.EQ.IO) THEN
           dWr=0.
        endif


        Wr=Wr+dWr

        !Calculate Gauss integeral between outer segment and second stretched segment
        r2=R(IO,:)
        r1=R(IS2,:)

        r12=r2-r1

        s2=SQRT(SUM((R(JP1,:)-R(IO,:))**2))
        s1=SQRT(SUM((R(IS2P1,:)-R(IS2,:))**2))
        e2=(R(JP1,:)-R(IO,:))/s2
        e1=(R(IS2P1,:)-R(IS2,:))/s1

        CALL GAUSSPAIR(R1,R2,e1,e2,s1,s2,dWr)
        !Writhe=0 between segments in plane (dWr=NaN)
        if (dWr.NE.dWr) then
           dWr=0.
        endif
        !Writhe=0 between adjacent segments
        if (JP1.EQ.IS2.OR.IS2P1.EQ.IO) THEN
           dWr=0.
        endif

        Wr=Wr+dWr
        IO=IO+1
     ENDDO
  ENDIF



  Wr=2.*Wr

  RETURN
END SUBROUTINE WRITHESLIDE










