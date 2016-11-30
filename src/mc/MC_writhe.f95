  !Calculate the writhe between segment that was moved and segment that was not moved using the method of Klenin (2000)

SUBROUTINE WRITHE_MOVE(R,RP,IB1,IB2,IT1,IT2,IP,N,Wr)
  IMPLICIT NONE

  DOUBLE PRECISION, PARAMETER :: PI=3.141592654d0 ! Value of pi


  INTEGER, INTENT(IN) :: N                 ! Number of beads per polymer
  DOUBLE PRECISION, INTENT(IN) :: R(N,3),RP(N,3)  ! Position
  INTEGER, INTENT(IN) :: IB1,IB2,IT1,IT2   ! Bead indices for terminal points of segment moved
  INTEGER, INTENT(IN) :: IP                ! Test polymer index
  !Geometric variables
  DOUBLE PRECISION  r1(3)                  ! Position bead 1
  DOUBLE PRECISION  r2(3)                  ! Position bead 2
  DOUBLE PRECISION  r12(3)                 ! Relative position vector
  DOUBLE PRECISION  s1                     ! Length of segment 1
  DOUBLE PRECISION  s2                     ! Length of segment 2
  DOUBLE PRECISION  beta                   ! Angle between tangents
  DOUBLE PRECISION  e1(3)                  ! Tangent of first segment
  DOUBLE PRECISION  e2(3)                  ! Tangent of second segment
  DOUBLE PRECISION  cosB
  DOUBLE PRECISION  sin2B
  DOUBLE PRECISION  e3(3)
  DOUBLE PRECISION  a0
  DOUBLE PRECISION  a1
  DOUBLE PRECISION  a2
  DOUBLE PRECISION  a3

  !Variables for writhe integral
  DOUBLE PRECISION  t1
  DOUBLE PRECISION  t2
  DOUBLE PRECISION  F1
  DOUBLE PRECISION  F2
  DOUBLE PRECISION  F3
  DOUBLE PRECISION  F4

  !Counter variables
  INTEGER  I
  INTEGER  J
  INTEGER NBI,NBO,II,IO,IIP1,IOP1
  DOUBLE PRECISION  dWr
  DOUBLE PRECISION, INTENT(OUT) :: Wr      ! Writhe


  IF (IB2.GE.IB1) THEN
     NBI=IB2-IB1
     NBO=N-NBI
  ELSE
     NBO=IB1-IB2
     NBI=N-NBO
  ENDIF


  Wr=0.
  II=IB1

  DO  I=1,NBI-1
     IO=IB1
     IF (II.gt.N) THEN
        II=1
        IIP1=II+1
     ELSEIF (II.EQ.N) THEN
        IIP1=1
     ELSE
        IIP1=II+1
     ENDIF

     DO J=1,NBI-1


        IF (IO.GT.N) THEN
           IO=1
           IOP1=IO+1
        ELSEIF (IO.EQ.N) THEN
           IOP1=1
        ELSE
           IOP1=IO+1
        ENDIF

        r1=R(II,:)
        r2=R(IO,:)
        r12=r2-r1

        if (RP(Io,1).NE.R(Io,1).OR.RP(Io,2).NE.R(Io,2).OR.RP(Io,3).NE.R(Io,3)) THEN
           ! PRINT *, "POSITIONS NOT EQUAL", IO,"of", NBO
        ENDIF

        if (RP(IOP1,1).NE.R(IOP1,1).OR.RP(IOP1,2).NE.R(IOP1,2).OR.RP(IOP1,3).NE.RP(IOP1,3)) THEN
           !  PRINT *, "POSITIONS NOT EQUAL AT IOP1", IOP1
        endif


        s2=SQRT(SUM((R(IOP1,:)-R(IO,:))**2))
        s1=SQRT(SUM((R(IIP1,:)-R(II,:))**2))
        e2=(R(IOP1,:)-R(IO,:))/s2
        e1=(R(IIP1,:)-R(II,:))/s1
        cosB=DOT_PRODUCT(e1,e2)

        sin2B=1-(cosB**2)

        e3(1)=e1(2)*e2(3)-e1(3)*e2(2)
        e3(2)=e1(3)*e2(1)-e1(1)*e2(3)
        e3(3)=e1(1)*e2(2)-e1(2)*e2(1)

        a1=DOT_PRODUCT(r12,e2*cosB-e1)/(sin2B)
        a2=DOT_PRODUCT(r12,e2-e1*cosB)/sin2B
        a0=DOT_PRODUCT(r12,e3)/sin2B



        t1=a1+s1
        t2=a2+s2
        F1=-ATAN((t1*t2+(a0**2)*cosB)/(a0*SQRT((t1**2+t2**2-2.*t1*t2*cosB+(a0**2)*sin2B))))/(4.*PI)
        t1=a1+s1
        t2=a2
        F2=-ATAN((t1*t2+(a0**2)*cosB)/(a0*SQRT((t1**2+t2**2-2*t1*t2*cosB+(a0**2)*sin2B))))/(4.*PI)
        t1=a1
        t2=a2+s2
        F3=-ATAN((t1*t2+(a0**2)*cosB)/(a0*SQRT((t1**2+t2**2-2*t1*t2*cosB+(a0**2)*sin2B))))/(4.*PI)
        t1=a1
        t2=a2
        F4=-ATAN((t1*t2+(a0**2)*cosB)/(a0*SQRT((t1**2+t2**2-2*t1*t2*cosB+(a0**2)*sin2B))))/(4.*PI)
        dWr=F1-F2-F3+F4

        if (dWr.NE.dWr) then
           dWr=0.
        endif

        Wr=Wr+dWr

        IO=IO+1
     ENDDO
     II=II+1
  ENDDO

  Wr=2.*Wr

  RETURN
END SUBROUTINE WRITHE_MOVE










