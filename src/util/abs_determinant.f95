! This subroutine calculates the absolute value of the determinant of a matrix
! using a LAPACK Q-R factorization. The absolute value of the determinant is
! the sum of the absolute values of the elements on the diagonal of the R
! matrix
! Calculate the absolute value of the determinant of a square NXN
! matrix using an LU factorization by LAPACK


subroutine abs_determinant(A,N,det)

  use params, only : dp

  integer N       !Size of the matrix A
  real(dp) A(N,N)
  real(dp) det
  integer LDA
  integer inFO
  integer I
  real(dp) LU(N,N)
  integer IPIV(N,N)

  LDA = N
  LWORK = N
  LU = A

  CALL DGETF2(N,N,A,LDA,IPIV,inFO)
  !Returned LU  matrix contains the U matrix from the LU factorization
  !on the upper diagonal. Take the absolute value of the product of the diagonals
  !to get the absolute value of the determinant

  det = 1.

  do I = 1,N
     det = det*A(I,I)
  ENDdo
  det = ABS(det)


END subroutine abs_determinant

