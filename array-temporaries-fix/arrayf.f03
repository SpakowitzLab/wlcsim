module arrayf
contains
    function f(R)
        real, intent(in) :: R(3)
        real f(3)
        f = R + 1
    end function

    function rotateR(ROT,R)
        implicit none
        real, intent(in) :: ROT(3,3) ! Rotation matrix
        real, intent(in) :: R(3)
        real :: rotateR(3)
        rotateR(1) = ROT(1,1)*R(1) + ROT(1,2)*R(2) + ROT(1,3)*R(3)
        rotateR(2) = ROT(2,1)*R(1) + ROT(2,2)*R(2) + ROT(2,3)*R(3)
        rotateR(3) = ROT(3,1)*R(1) + ROT(3,2)*R(2) + ROT(3,3)*R(3)
    end function rotateR
    function E_SSWLC(R)
        implicit none
        real, intent(in), dimension(3) :: R ! R of bead i
        real, dimension(4) :: E_SSWLC
        !real DR(3)
        !real DRPAR
        !real DRPERP(3)
        !real GI(3)

        !DR = R-RM1
        !DRPAR = dot_product(DR,UM1)
        !DRPERP = DR - DRPAR*UM1
        !GI = U - UM1 - ETA*DRPERP
        !E_SSWLC(1)=0.5*EB*dot_product(GI,GI)
        !E_SSWLC(2)=0.5*EPAR*(DRPAR-GAM)**2
        !E_SSWLC(3)=0.5*EPERP*dot_product(DRPERP,DRPERP)
        !E_SSWLC(4)=0.0

        !DR = R-RM1
        E_SSWLC(1)=0.0
        E_SSWLC(2)=0.5
        E_SSWLC(3)=0.5
        E_SSWLC(4)=0.0
    end function E_SSWLC
    function simple()
        implicit none
        real, dimension(4) :: simple
        simple(1)=0.0
        simple(2)=0.5
        simple(3)=0.5
        simple(4)=0.0
    end function
end module

