!---------------------------------------------------------------
! Sort an index-array and an index-index array by a double array target
!---------------------------------------------------------------

! use quicksort
subroutine qcolsort(array_size, indexi, index, value)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer, parameter :: QSORT_THRESHOLD = 32
    integer, intent(in) :: array_size
    integer, intent(inout) :: index(array_size)
    integer, intent(inout) :: indexi(array_size)
    real(dp), intent(in) :: value(array_size)

    include "qsort_inline.inc"
contains
    include "colsort.inc"
end subroutine qcolsort

! use insertion sort
subroutine icolsort(array_size, indexi, index, value)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer, intent(in) :: array_size
    integer, intent(inout) :: index(array_size)
    integer, intent(inout) :: indexi(array_size)
    real(dp), intent(in) :: value(array_size)

    integer :: left, right
    do right = 2, array_size
        left = right - 1
        if (less_than(right, left)) then
            do  ! need two separate if's since fortran has no short circuit "&&"
                if (left < 2) exit
                if (.NOT.less_than(right, left - 1)) exit
                left = left - 1
            enddo
            call rshift(left, right)
        endif
    enddo
contains
    include "colsort.inc"
end subroutine icolsort
