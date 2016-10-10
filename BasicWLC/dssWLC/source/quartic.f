C     ***QUARTIC************************************************25.03.98
C     Solution of a quartic equation
C     ref.: J. E. Hacke, Amer. Math. Monthly, Vol. 48, 327-328, (1941)
C     NO WARRANTY, ALWAYS TEST THIS SUBROUTINE AFTER DOWNLOADING
C Downloaded from: http://van-der-waals.pc.uni-koeln.de/quartic/quartic.f
C     ******************************************************************
C     dd(0:4)     (i)  vector containing the polynomial coefficients
C     sol(1:4)    (o)  results, real part
C     soli(1:4)   (o)  results, imaginary part
C     Nsol        (o)  number of real solutions 
C     ==================================================================
      subroutine quartic(dd,sol,soli,Nsol)
      implicit double precision (a-h,o-z)
      dimension dd(0:4),sol(4),soli(4)
      dimension AA(0:3),z(3)
C
      Nsol = 0
      a = dd(4)
      b = dd(3)
      c = dd(2)
      d = dd(1)
      e = dd(0)
C
      if (dd(4).eq.0.d+0) then
	write(6,*)'ERROR: NOT A QUARTIC EQUATION'
	return
      endif
C
      p = (-3.d+0*b**2 + 8.d+0*a*c)/(8.d+0*a**2)
      q = (b**3 - 4.d+0*a*b*c + 8.d+0*d*a**2)/(8.d+0*a**3)
      r = (-3.d+0*b**4 + 16.d+0*a*b**2*c - 64.d+0*a**2*b*d + 
     &      256.d+0*a**3*e)/(256.d+0*a**4)
C
C     solve cubic resolvent
      AA(3) =  8.d+0
      AA(2) = -4.d+0*p 
      AA(1) = -8.d+0*r
      AA(0) =  4.d+0*p*r - q**2
      call cubic(AA,z,ncube)
C      
      zsol = -1.d+99
      do 5 i=1,ncube
 5      zsol = max(zsol,z(i))
      z(1) = zsol
      xK2 = 2.d+0 * z(1) - p
      xK  = sqrt(xK2)
C-----------------------------------------------
      if (xK.eq.0.d+0) then
        xL2 = z(1)**2 - r
	if (xL2.lt.0.d+0) then
	  write(6,*)'Sorry, no solution'
	  return
        endif
	xL  = sqrt(xL2)
      else
        xL = q/(2.d+0 * xK)
      endif
C-----------------------------------------------
      sqp = xK2 - 4.d+0*(z(1) + xL)
      sqm = xK2 - 4.d+0*(z(1) - xL)
C
      do 10 i=1,4
 10     soli(i) = 0.d+0
      if       (sqp.ge.0.d+0 .and. sqm.ge.0.d+0) then
	sol(1) = 0.5d+0*( xK + sqrt(sqp))
	sol(2) = 0.5d+0*( xK - sqrt(sqp))
	sol(3) = 0.5d+0*(-xK + sqrt(sqm))
	sol(4) = 0.5d+0*(-xK - sqrt(sqm))
	Nsol = 4
      else if  (sqp.ge.0.d+0 .and. sqm.lt.0.d+0) then
	sol(1) =  0.5d+0*(xK + sqrt(sqp))
	sol(2) =  0.5d+0*(xK - sqrt(sqp))
	sol(3) = -0.5d+0*xK 
	sol(4) = -0.5d+0*xK 
	soli(3) =  sqrt(-0.25d+0 * sqm)
	soli(4) = -sqrt(-0.25d+0 * sqm)
	Nsol = 2
      else if  (sqp.lt.0.d+0 .and. sqm.ge.0.d+0) then
	sol(1) = 0.5d+0*(-xK + sqrt(sqm))
	sol(2) = 0.5d+0*(-xK - sqrt(sqm))
	sol(3) =  0.5d+0*xK 
	sol(4) =  0.5d+0*xK 
	soli(3) =  sqrt(-0.25d+0 * sqp)
	soli(4) = -sqrt(-0.25d+0 * sqp)
	Nsol = 2
      else if  (sqp.lt.0.d+0 .and. sqm.lt.0.d+0) then
	sol(1) = -0.5d+0*xK 
	sol(2) = -0.5d+0*xK 
	soli(1) =  sqrt(-0.25d+0 * sqm)
	soli(2) = -sqrt(-0.25d+0 * sqm)
	sol(3) =  0.5d+0*xK 
	sol(4) =  0.5d+0*xK 
	soli(3) =  sqrt(-0.25d+0 * sqp)
	soli(4) = -sqrt(-0.25d+0 * sqp)
	Nsol = 0
      endif
      do 20 i=1,4
 20     sol(i) = sol(i) - b/(4.d+0*a)
C
      return
      END

C     ***CUBIC************************************************08.11.1986
C     Solution of a cubic equation
C     Equations of lesser degree are solved by the appropriate formulas.
C     The solutions are arranged in ascending order.
C     NO WARRANTY, ALWAYS TEST THIS SUBROUTINE AFTER DOWNLOADING
C     ******************************************************************
C     A(0:3)      (i)  vector containing the polynomial coefficients
C     X(1:L)      (o)  results
C     L           (o)  number of valid solutions (beginning with X(1))
C     ==================================================================
      SUBROUTINE CUBIC(A,X,L)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(0:3),X(3),U(3)
      PARAMETER(PI=3.1415926535897932D+0,THIRD=1.D+0/3.D+0)
      INTRINSIC MIN,MAX,ACOS
C
C     define cubic root as statement function
      CBRT(Z)=SIGN(ABS(Z)**THIRD,Z)
C
C     ==== determine the degree of the polynomial ====
C
      IF (A(3).NE.0.D+0) THEN
C
C       cubic problem
	W=A(2)/A(3)*THIRD
	P=(A(1)/A(3)*THIRD-W**2)**3
	Q=-.5D+0*(2.D+0*W**3-(A(1)*W-A(0))/A(3))
	DIS=Q**2+P
	IF (DIS.LT.0.D+0) THEN
C         three real solutions!
C         Confine the argument of ACOS to the interval [-1;1]!
	  PHI=ACOS(MIN(1.D+0,MAX(-1.D+0,Q/SQRT(-P))))
	  P=2.D+0*(-P)**(5.D-1*THIRD)
	  DO 100 I=1,3
  100       U(I)=P*COS((PHI+DBLE(2*I)*PI)*THIRD)-W
	  X(1)=MIN(U(1),U(2),U(3))
	  X(2)=MAX(MIN(U(1),U(2)),MIN(U(1),U(3)),MIN(U(2),U(3)))
	  X(3)=MAX(U(1),U(2),U(3))
	  L=3
	ELSE
C         only one real solution!
	  DIS=SQRT(DIS)
	  X(1)=CBRT(Q+DIS)+CBRT(Q-DIS)-W
	  L=1
	END IF
C
      ELSE IF (A(2).NE.0.D+0) THEN
C
C       quadratic problem
	P=5.D-1*A(1)/A(2)
	DIS=P**2-A(0)/A(2)
	IF (DIS.GE.0.D+0) THEN
C         two real solutions!
	  X(1)=-P-SQRT(DIS)
	  X(2)=-P+SQRT(DIS)
	  L=2
	ELSE
C         no real solution!
	  L=0
	END IF
C
      ELSE IF (A(1).NE.0.D+0) THEN
C
C       linear equation
	X(1)=-A(0)/A(1)
	L=1
C
      ELSE
C       no equation
	L=0
      END IF
C
C     ==== perform one step of a newton iteration in order to minimize
C          round-off errors ====
      DO 110 I=1,L
	X(I)=X(I)-(A(0)+X(I)*(A(1)+X(I)*(A(2)+X(I)*A(3))))
     *  /(A(1)+X(I)*(2.D+0*A(2)+X(I)*3.D+0*A(3)))
  110 CONTINUE
      RETURN
      END
