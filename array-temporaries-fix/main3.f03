program main

    use arrayf, only : f

    real, allocatable :: R(:,:)
    real, allocatable :: U(:,:)

    integer I

    allocate(R(3,3))
    allocate(U(3,3))
    I = 1

    R(:,:) = 1
    U(:,:) = 1

    PRINT *, R

    U(:,I) = f(R(:,I))

    PRINT *, U

end program main
