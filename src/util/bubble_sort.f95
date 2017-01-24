!Sort an NXM matrix by a given column using a
!Bubble sort algorithm (slow)

SUBROUTINE bubble_sort(A,N,M,isort)

  INTEGER N !Number of rows
  INTEGER M !Number of columns
  DOUBLE PRECISION A(N,M)
  INTEGER isort !column on which to sort rows of matrix
  INTEGER I,J
  DOUBLE PRECISION Temp(M)
  DO I=1,N
     DO J=1,N-1
        IF (A(J,isort).GT.A(J+1,isort)) THEN
           Temp=A(J,:)
           A(J,:)=A(J+1,:)
           A(J+1,:)=Temp
        ENDIF
     ENDDO
  ENDDO

  RETURN
ENDSUBROUTINE bubble_sort

