module holds
    type h1
        real, allocatable :: R(:,:)
        real, allocatable :: U(:,:)
    end type
contains
    subroutine alloc(h)
        type(h1), intent(out) :: h
        allocate(h%R(3,3))
        allocate(h%U(3,3))
    end subroutine
end module


program main

    use holds, only : h1, alloc
    use arrayf, only : rotateR

    integer I

    type(h1) h
    real ROT(3,3)

    ROT(:,:) = 1.0

    call alloc(h)

    I = 1

    h%R(:,:) = 1.0
    h%U(:,:) = 1.0

    PRINT *, h%R

    h%U(:,I) = rotateR(ROT, h%R(:,I))

    PRINT *, h%U

end program main

