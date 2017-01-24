logical function isanynan(arr)
    use params, only : dp

    implicit none

    real(dp), dimension(:,:) :: arr
    integer, dimension(2) :: shapes
    integer :: i, j

    shapes = shape(arr)
    isanynan = .false.
    do i = 1, shapes(1)
        do j = 1, shapes(2)
            if (isnan(arr(i,j))) then
                isanynan = .true.
                return
            endif
        enddo
    enddo
endfunction

logical function b_any(arr)
    implicit none
    logical :: arr(:,:)
    b_any = any(arr)
end function
    ! integer :: r

    ! r = rank(arr)
    ! isanynan = isanynan_helper(arr, r)
! end function

! logical isanynan_helper(arr, rank)
    ! real(REAL64), intent(in), dimension(*) :: arr
    ! integer, intent(in) :: rank
    ! integer, dimension(rank) :: sizes

    ! sizes = shape(arr)
    ! do i=1,
