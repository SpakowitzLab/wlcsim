module holds
    real :: R(3,3)
    real :: U(3,3)
end module


program main

    use holds, only : R, U
    use arrayf, only : rotateR

    integer I

    real ROT(3,3)

    ROT(:,:) = 1.0

    I = 1

    R(:,:) = 1.0
    U(:,:) = 1.0

    PRINT *, R

    U(:,I) = rotateR(ROT, R(:,I))

    PRINT *, U

end program main

