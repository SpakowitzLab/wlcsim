SUBROUTINE A0BOUNDFUNC(R,PARAM,FU,DU)
! function used to find enveloping Lorentzian for cylindrical gaussian
! PARAM is array of parameters (A,B,R0,M0)
! FU returns function values, DU returns derivative
! sign flipped upside down to do minimization rather than maximization
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: R, PARAM(4)
  DOUBLE PRECISION, INTENT(OUT) :: FU, DU
  DOUBLE PRECISION :: A, B, R0, M0

  PRINT*, 'TESTX0:', PARAM
  A = PARAM(1); B = PARAM(2); R0 = PARAM(3); M0 = PARAM(4)

  IF (ABS(R-R0).LT.1D-15) THEN
     ! use asymptotic form near r0
     DU =  (-B+SQRT(2/A+B**2))/(6+3*A*B**2)
     FU = 1 + B*SQRT(A/(2+A*B**2))/(2*A) + DU*(R-R0)
  ELSE
     FU = R*(R-R0)**2/(M0*EXP(A*(R-B)**2) - R)
     DU = (R-R0)*(-2*R**2-EXP(A*(B-R)**2)*M0*(R*(-3-2*A*(B-R)*(R-R0))+R0)) / &
          & (R-M0*EXP(A*(B-R)**2))**2
  END IF

  FU = -FU; DU = -DU;

  PRINT*, 'TESTX1:', R, PARAM, FU, DU
END SUBROUTINE A0BOUNDFUNC
