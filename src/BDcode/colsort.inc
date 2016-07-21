! helper routines for sorting an array ind based on an array r, while keeping a
! third array indi such that:
! r(ind(:,k),k) is in order for k in 1,2,3    i.e. [~,ind(:,1)] = sort(r(:,1)
! ind(indi(i,k),k) == i for k in 1,2,3
! (indi(:,k) == indexi, ind(:,k) == index, r(:,k) == value
!
! in reasoning about the following, it is useful to note that ind and indi as
! functions from 1,n to 1,n are inverses: indi(ind(i)) == i && ind(indi(i)) == i

    ! do-nothing initialization routine required for qsort_inline.inc
    subroutine init()
    end subroutine init

    ! swap indices a,b
    subroutine swap(a,b)
    ! less_than(index(a),index(b)) will lead to a swap
        integer, intent(in) :: a,b
        integer :: hold
        hold=index(a)
        index(a)=index(b)
        index(b)=hold
        ! easy to see this is what you want if you know you want to make
        ! indexi the inverse of index (as a function from N to N)
        indexi(index(a)) = a
        indexi(index(b)) = b
    end subroutine swap

    ! circular shift-right by one:
    subroutine rshift(left,right)
        implicit none
        integer, intent(in) :: left, right
        integer :: hold, i
        hold=index(right)
        ! This syntax is valid, but has poor optimization in GFortran:
        ! index(left+1:right)=index(left:right-1)
        do i=right,left+1,-1
            index(i)=index(i-1)
        end do
        index(left)=hold
        do i=left,right
            indexi(index(i)) = indexi(index(i)) + 1
        end do
        indexi(hold) = left
    end subroutine rshift

    ! fuzzy comparator prevents unecessary switching in super-fine time steps
    logical &
    function less_than(a,b)
        integer, intent(in) :: a,b
!       real(dp), parameter :: small=1.0e-15
!       if ( abs(value(index(a))-value(index(b))) < small ) then
!           less_than = index(a) < index(b)
!       else
            less_than = value(index(a)) < value(index(b))
!       end if
   end function less_than
