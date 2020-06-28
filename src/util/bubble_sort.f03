!Sort an NXM matrix by a given column using a
!Bubble sort algorithm (slow)

subroutine bubble_sort(A, N, M, isort)
   use params, only: dp
   implicit none
   integer N !Number of rows
   integer M !Number of columns
   real(dp) A(N, M)
   integer isort !column on which to sort rows of matrix
   integer I, J
   real(dp) Temp(M)
   do I = 1, N
      do J = 1, N - 1
         if (A(J, isort) > A(J + 1, isort)) then
            Temp = A(J, :)
            A(J, :) = A(J + 1, :)
            A(J + 1, :) = Temp
         ENDif
      ENDdo
   ENDdo

   RETURN
ENDsubroutine bubble_sort
