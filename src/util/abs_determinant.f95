! This subroutine calculates the absolute value of the determinant of a matrix
! using a LAPACK Q-R factorization. The absolute value of the determinant is
! the sum of the absolute values of the elements on the diagonal of the R
! matrix
! Calculate the absolute value of the determinant of a square NXN
! matrix using an LU factorization by LAPACK


SUBROUTINE abs_determinant(A,N,det)
  INTEGER N       !Size of the matrix A
  DOUBLE PRECISION A(N,N)
  DOUBLE PRECISION det
  INTEGER LDA
  INTEGER INFO
  INTEGER I
  DOUBLE PRECISION LU(N,N)
  INTEGER IPIV(N,N)

  LDA=N
  LWORK=N
  LU=A

  CALL DGETF2(N,N,A,LDA,IPIV,INFO)
  !Returned LU  matrix contains the U matrix from the LU factorization
  !on the upper diagonal. Take the absolute value of the product of the diagonals
  !to get the absolute value of the determinant

  det=1.

  DO I=1,N
     det=det*A(I,I)
  ENDDO
  det=ABS(det)


END SUBROUTINE abs_determinant

