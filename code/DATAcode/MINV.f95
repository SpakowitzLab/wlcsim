!---------------------------------------------------------------*

!
! This subroutine finds the inverse of a 3x3 matrix and finds
! its determinant
!
! Andrew Spakowitz
! Written 9-1-04

      SUBROUTINE MINV(A,B,DET)

      DOUBLE PRECISION A(3,3)        ! Original matrix
      DOUBLE PRECISION B(3,3)        ! Inverse matrix
      DOUBLE PRECISION DET        ! Determinant
      DOUBLE PRECISION COA(3,3)        ! Cofactor of A


! First find the determinant

      DET = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))

! Find the Cofactor of A

      COA(1,1) = A(2,2)*A(3,3)-A(2,3)*A(3,2)
      COA(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
      COA(1,3) = A(2,1)*A(3,2)-A(2,2)*A(3,1)
      COA(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
      COA(2,2) = A(1,1)*A(3,3)-A(1,3)*A(3,1)
      COA(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
      COA(3,1) = A(1,2)*A(2,3)-A(1,3)*A(2,2)
      COA(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
      COA(3,3) = A(1,1)*A(2,2)-A(1,2)*A(2,1)

! Find the inverse of A

      B(1,1) = COA(1,1)/DET
      B(1,2) = COA(2,1)/DET
      B(1,3) = COA(3,1)/DET
      B(2,1) = COA(1,2)/DET
      B(2,2) = COA(2,2)/DET
      B(2,3) = COA(3,2)/DET
      B(3,1) = COA(1,3)/DET
      B(3,2) = COA(2,3)/DET
      B(3,3) = COA(3,3)/DET

      RETURN
      END

!---------------------------------------------------------------*