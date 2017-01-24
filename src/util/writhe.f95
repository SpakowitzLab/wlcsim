  !Calculate the Writhe of a curve using the method of Klenin (2000)

SUBROUTINE WRITHE(R,N,Wr)
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
  DOUBLE PRECISION  beta                   ! Angle between tangents
  DOUBLE PRECISION  e1(3)                  ! Tangent of first segment
  DOUBLE PRECISION  e2(3)                  ! Tangent of second segment
  DOUBLE PRECISION  cosB
  DOUBLE PRECISION  sin2B
  DOUBLE PRECISION  e3(3)
  DOUBLE PRECISION  OMEGA(N,N)             ! Matrix of solid angles
  DOUBLE PRECISION  a0
  DOUBLE PRECISION  a1
  DOUBLE PRECISION  a2
  DOUBLE PRECISION  a3
  DOUBLE PRECISION  b
   !Variables for writhe integral
  DOUBLE PRECISION  t1
  DOUBLE PRECISION  t2
  DOUBLE PRECISION  F1
  DOUBLE PRECISION  F2
  DOUBLE PRECISION  F3
  DOUBLE PRECISION  F4

  !Counter variables
  INTEGER  I,IP1
  INTEGER  J,JP1

  DOUBLE PRECISION  dWr
  DOUBLE PRECISION, INTENT(OUT) :: Wr      ! Writhe


  OMEGA=0
  Wr=0.

  DO  I=2,N
      IF (I.EQ.N) THEN
           IP1=1
        ELSE
           IP1=I+1
        ENDIF
     DO J=1,I-1



        IF (J.EQ.N) THEN
           JP1=1
        ELSE
           JP1=J+1
        ENDIF



        r1=R(I,:)
        r2=R(J,:)
        r12=r2-r1

        s2=SQRT(SUM((R(JP1,:)-R(J,:))**2))
        s1=SQRT(SUM((R(IP1,:)-R(I,:))**2))
        e2=(R(JP1,:)-R(J,:))/s2
        e1=(R(IP1,:)-R(I,:))/s1
        cosB=DOT_PRODUCT(e1,e2)
        B=ACOS(cosB)

       ! sin2B=1.-(cosB**2.)
        sin2B=sin(B)**2.
        e3(1)=e1(2)*e2(3)-e1(3)*e2(2)
        e3(2)=e1(3)*e2(1)-e1(1)*e2(3)
        e3(3)=e1(1)*e2(2)-e1(2)*e2(1)

        if (abs(sin2B).lt.(10.**(-15.)))then
           sin2b=0.
        endif

        a1=DOT_PRODUCT(r12,e2*cosB-e1)/(sin2B)
        a2=DOT_PRODUCT(r12,e2-e1*cosB)/sin2B
        a0=DOT_PRODUCT(r12,e3)/sin2B



           t1=a1+s1
           t2=a2+s2
           F1=-ATAN((t1*t2+(a0**2.)*cosB)/(a0*SQRT((t1**2.+t2**2.-2.*t1*t2*cosB+(a0**2.)*sin2B))))/(4.*PI)
           t1=a1+s1
           t2=a2
           F2=-ATAN((t1*t2+(a0**2.)*cosB)/(a0*SQRT((t1**2.+t2**2.-2.*t1*t2*cosB+(a0**2.)*sin2B))))/(4.*PI)
           t1=a1
           t2=a2+s2
           F3=-ATAN((t1*t2+(a0**2.)*cosB)/(a0*SQRT((t1**2.+t2**2.-2.*t1*t2*cosB+(a0**2.)*sin2B))))/(4.*PI)
           t1=a1
           t2=a2
           F4=-ATAN((t1*t2+(a0**2.)*cosB)/(a0*SQRT((t1**2.+t2**2.-2.*t1*t2*cosB+(a0**2.)*sin2B))))/(4.*PI)
           dWr=F1-F2-F3+F4

           if (dWr.NE.dWr) then
              dWr=0.
           endif
           if (JP1.EQ.I.OR.IP1.EQ.J) THEN
              dWr=0.
           endif

           if (abs(a0).LT.10.**(-10.)) then
              dWr=0.
           endif


           Wr=Wr+dWr


     ENDDO
  ENDDO

  Wr=2.*Wr

  RETURN
END SUBROUTINE WRITHE










