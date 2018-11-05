module holds
    real, allocatable :: R(:,:)
    real, allocatable :: U(:,:)
contains
    subroutine alloc()
        allocate(R(3,3))
        allocate(U(3,3))
    end subroutine
end module


program main

    use holds, only : R, U, alloc
    use arrayf, only : f

    integer I

    call alloc()

    I = 1

    R(:,:) = 1.0
    U(:,:) = 1.0

    PRINT *, R

    U(:,I) = f(R(:,I))

    PRINT *, U

end program main

