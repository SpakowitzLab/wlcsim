module holds
    type h1
        real :: R(3,3)
        real :: U(3,3)
    end type
end module


program main

    use holds, only : h1
    use arrayf, only : rotateR

    integer I

    type(h1) h
    real ROT(3,3)

    ROT(:,:) = 1.0

    I = 1

    h%R(:,:) = 1.0
    h%U(:,:) = 1.0

    PRINT *, h%R

    h%U(:,I) = rotateR(ROT, h%R(:,I))

    PRINT *, h%U

end program main

