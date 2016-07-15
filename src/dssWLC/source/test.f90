PROGRAM MAIN

  call testbessel
  !CALL TESTMULTINORM
  !CALL TESTCHOLESKY
  !CALL TESTMULTINORM
CONTAINS
  
  SUBROUTINE TESTBESSEL
    ! test bessel function
    implicit none
    DOUBLE PRECISION :: X, ANS

    X = 5D0
    CALL BESSF_I0(X,ANS)
    PRINT*, ANS
  END SUBROUTINE TESTBESSEL

SUBROUTINE TESTMULTINORM
  ! test multivariate normal random number generator
  USE MT19937, ONLY : RNORM, MVNORM
  IMPLICIT NONE
  INTEGER, PARAMETER :: N = 2, NDEV = 10000
  DOUBLE PRECISION :: COV(N,N), CHOLCOV(N,N)
  DOUBLE PRECISION :: Y(N), XDEV(NDEV,N), MU(N)
  INTEGER :: INFO, I, J

  MU = (/1D0,-2D0/)

  ! covariance matrix
  COV(1,:) = (/4D0,0D0/)
  COV(2,:) = (/2D0,3D0/)

  CHOLCOV = COV
  CALL DPOTRF('L',N,CHOLCOV,N,INFO)
  
 ! PRINT*, 'INFO:', INFO
  
  DO I = 1,NDEV
     XDEV(I,:) = MVNORM(N,MU,CHOLCOV)
      PRINT*, XDEV(I,:)
  END DO

  ! ! Get standard normal random deviates
  ! DO I = 1,NDEV
  !    DO J = 1,N
  !       Y(J) = RNORM()
  !    ENDDO
     
  !    ! convert to multivariate normal deviates
      CALL DTRMV('L','N','N',N,CHOLCOV,N,Y,1)
  !    XDEV(I,:) = Y + MU
  !    PRINT*, XDEV(I,:)
  ! END DO
  
  

  END SUBROUTINE TESTMULTINORM

SUBROUTINE TESTCHOLESKY
  ! test cholesky decomposition with lapack
  IMPLICIT NONE
  INTEGER, PARAMETER :: N = 3, LDA = N
  DOUBLE PRECISION :: A(N,N), XVEC(N)
  INTEGER :: INFO, I

  A(1,:) = (/6D0,0D0,0D0/)
  A(2,:) = (/2D0,4D0,0D0/)
  A(3,:) = (/3D0, 0D0, 5D0/)

  CALL DPOTRF('L',N,A,LDA,INFO)

  PRINT*, 'INFO:', INFO

  DO I = 1,N
     PRINT*, A(I,:)
  END DO

  ! test matrix-vector multiple
  XVEC = (/1D0,2D0,-3D0/)
  CALL DTRMV('L','N','N',N,A,N,XVEC,1)

  PRINT*, 'XVEC:', XVEC
END SUBROUTINE TESTCHOLESKY

! SUBROUTINE TESTQUARTIC
!   ! test the quartic equation solver

!   IMPLICIT NONE
!   DOUBLE PRECISION :: COEFF(5), SOL(4),SOLI(4)
!   INTEGER :: NSOL

!   !COEFF = (/2.,5.1,-1.,0.,-1./)
!   COEFF = (/-1d0,0d0,-1d0,5.1d0,2d0/)

!   CALL QUARTIC(COEFF,SOL,SOLI,NSOL)

!   PRINT*, NSOL
!   PRINT*, 'SOL', SOL
!   PRINT*, 'SOLI', SOLI
! END SUBROUTINE TESTQUARTIC
END PROGRAM MAIN
