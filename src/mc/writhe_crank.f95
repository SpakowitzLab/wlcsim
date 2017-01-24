  !Calculate the portion of the writhe that changes during a crankshaft move

SUBROUTINE WRITHECRANK(R,IT1,IT2,N,Wr)
  IMPLICIT NONE

  DOUBLE PRECISION, PARAMETER :: PI=3.141592654d0 ! Value of pi


  INTEGER, INTENT(IN) :: N                 ! Number of beads
  DOUBLE PRECISION, INTENT(IN) :: R(N,3)  ! Positions
  !Geometric variables
  DOUBLE PRECISION  r1(3)                  ! Position bead 1
  DOUBLE PRECISION  r2(3)                  ! Position bead 2
  DOUBLE PRECISION  r12(3)                 ! Relative position vector
  DOUBLE PRECISION  s1                     ! Length of segment 1
  DOUBLE PRECISION  s2                     ! Length of segment 2
  DOUBLE PRECISION  beta                   ! Angle between tangents
  DOUBLE PRECISION  e1(3)                  ! Tangent of first segment
  DOUBLE PRECISION  e2(3)                  ! Tangent of second segment
  DOUBLE PRECISION  e3(3)
  INTEGER  DIB,DIO                         !Number of segments inside section rotated,outside


  !Counter variables
  INTEGER IT1,IT2                         !Indices of endpoints in segment rotated
  INTEGER II,IO                           !Index of bead inside segment rotated,outside
  INTEGER  I,IP1
  INTEGER  J,JP1

  DOUBLE PRECISION  dWr
  DOUBLE PRECISION, INTENT(OUT) :: Wr      ! Writhe



  Wr=0.
  IF (IT2.GE.IT1) THEN
     DIB=IT2-IT1
  ELSE
     DIB=(N-IT1)+IT2
  ENDIF

  DIO=N-DIB
  II=IT1



  DO  I=1,DIB

     IF (II.EQ.N) THEN
        IP1=1
     ELSEIF (II.EQ.N+1) THEN
        II=1
        IP1=II+1
     ELSE
        IP1=II+1
     ENDIF
     IO=IT2
     DO J=1,DIO



        IF (IO.EQ.N) THEN
           JP1=1
        ELSEIF (IO.EQ.N+1) THEN
           IO=1
           JP1=IO+1
        ELSE
           JP1=IO+1
        ENDIF



        r1=R(II,:)
        r2=R(IO,:)
        r12=r2-r1

        s2=SQRT(SUM((R(JP1,:)-R(IO,:))**2))
        s1=SQRT(SUM((R(IP1,:)-R(II,:))**2))
        e2=(R(JP1,:)-R(IO,:))/s2
        e1=(R(IP1,:)-R(II,:))/s1

        CALL GAUSSPAIR(R1,R2,e1,e2,s1,s2,dWr)
        !Wr=0 for segments in plane (Wr=NaN)
        if (dWr.NE.dWr) then
           dWr=0.
        endif
        !Wr=0 for adjacent segments
        if (JP1.EQ.II.OR.IP1.EQ.IO) THEN
           dWr=0.
        endif


        Wr=Wr+dWr

        IO=IO+1

     ENDDO
     II=II+1
  ENDDO


  Wr=2.*Wr

  RETURN
END SUBROUTINE WRITHECRANK










