!-----------------------------------------------------------------
!
!     This program linearly interpolates a bead at RBin into
!     8 bins indexed by IX, IY, IZ with weights WX, WY, WZ
!
!        Addapted from MC_int.f95
!        By Quinn MacPherson in summer 2016
!
!---------------------------------------------------------------
subroutine interp(confineType,RBin,LBOX,NBinX,dbin,IX,IY,IZ,WX,WY,WZ)
use params, only: dp
implicit none
integer, intent (in) :: confineType
real(dp), intent(inout) :: RBin(3) ! position or posiion within bin
real(dp), intent(in) :: LBOX(3) ! Side length of box
integer, intent(in) :: NBinX(3)      ! number of discritations in each direction
real(dp), intent(in) :: dbin  ! size of discritation
integer, intent(out) :: IX(2)  ! Output
integer, intent(out) :: IY(2) ! Output
integer, intent(out) :: IZ(2)  ! Output
real(dp), intent(out) :: WX(2) ! Output
real(dp), intent(out) :: WY(2) ! Output
real(dp), intent(out) :: WZ(2) ! Output
SELECT CASE (confineType)
CASE (0) ! Box from 0-LBOX, Bins split by boundaries
    ! Periodic BC
    RBin(1) = RBin(1)-floor(RBin(1)/LBOX(1))*LBOX(1)
    RBin(2) = RBin(2)-floor(RBin(2)/LBOX(2))*LBOX(2)
    RBin(3) = RBin(3)-floor(RBin(3)/LBOX(3))*LBOX(3)

    ! Binning
    IX(1) = ceiling(RBin(1)/dbin)
    IY(1) = ceiling(RBin(2)/dbin)
    IZ(1) = ceiling(RBin(3)/dbin)

    IX(2) = IX(1)-1
    IY(2) = IY(1)-1
    IZ(2) = IZ(1)-1

    ! Calculate the bin weighting
    WX(2) = (dbin*IX(1)-RBin(1))/dbin
    WX(1) = 1.0_dp-WX(2)
    WY(2) = (dbin*IY(1)-RBin(2))/dbin
    WY(1) = 1.0_dp-WY(2)
    WZ(2) = (dbin*IZ(1)-RBin(3))/dbin
    WZ(1) = 1.0_dp-WZ(2)

    ! Periodic BC on Bins:
    IX(1) = IX(1)-floor(REAL((IX(1)-1))/REAL(NBinX(1))) * NBinX(1)
    IX(2) = IX(2)-floor(REAL((IX(2)-1))/REAL(NBinX(1))) * NBinX(1)
    IY(1) = IY(1)-floor(REAL((IY(1)-1))/REAL(NBinX(2))) * NBinX(2)
    IY(2) = IY(2)-floor(REAL((IY(2)-1))/REAL(NBinX(2))) * NBinX(2)
    IZ(1) = IZ(1)-floor(REAL((IZ(1)-1))/REAL(NBinX(3))) * NBinX(3)
    IZ(2) = IZ(2)-floor(REAL((IZ(2)-1))/REAL(NBinX(3))) * NBinX(3)
CASE (1)
    ! Periodic BC
    RBin(1) = RBin(1)-floor(RBin(1)/LBOX(2))*LBOX(1)
    RBin(2) = RBin(2)-floor(RBin(2)/LBOX(1))*LBOX(2)

    ! Binning
    IX(1) = ceiling(RBin(1)/dbin)
    IY(1) = ceiling(RBin(2)/dbin)
    IZ(1) = nint(RBin(3)/dbin) + 1 ! Note 1.0 so that box centers are on half intigers

    IX(2) = IX(1)-1
    IY(2) = IY(1)-1
    IZ(2) = IZ(1)-1

    ! Calculate the bin weighting
    WX(2) = (dbin*IX(1)-RBin(1))/dbin   ! WX(2) = (RBin(1)-IX(1)*dbin)/(IX(2)*dbin-IX(1)*dbin)
    WX(1) = 1.0_dp-WX(2)              ! WX(1) = (IX(2)*dbin-RBin(1))/(IX(2)*dbin-IX(1)*dbin)
    WY(2) = (dbin*IY(1)-RBin(2))/dbin   ! WY(2) = (RBin(2)-IY(1)*dbin)/(IY(2)*dbin-IY(1)*dbin)
    WY(1) = 1.0_dp-WY(2)              ! WY(1) = (IY(2)*dbin-RBin(2))/(IY(2)*dbin-IY(1)*dbin)
    WZ(2) = (dbin*IZ(1)-0.5_dp*dbin-RBin(3))/dbin   ! WZ(2) = (RBin(3)-IZ(1)*dbin)/(IZ(2)*dbin-IZ(1)*dbin)
    WZ(1) = 1.0_dp-WZ(2)                   ! WZ(1) = (IZ(2)*dbin-RBin(3))/(IZ(2)*dbin-IZ(1)*dbin)

    if ((WZ(1).lt.0).OR.(WZ(2).lt.0)) then
        print*, "negitive W"
        stop 1
    endif

    ! Periodic BC on Bins:
    IX(1) = IX(1)-floor(REAL((IX(1)-1))/REAL(NBinX(1))) * NBinX(1)
    IX(2) = IX(2)-floor(REAL((IX(2)-1))/REAL(NBinX(1))) * NBinX(1)
    IY(1) = IY(1)-floor(REAL((IY(1)-1))/REAL(NBinX(2))) * NBinX(2)
    IY(2) = IY(2)-floor(REAL((IY(2)-1))/REAL(NBinX(2))) * NBinX(2)
CASE (2) ! Box confinement
    ! Binning
    IX(1) = nint(RBin(1)/dbin) + 1
    IY(1) = nint(RBin(2)/dbin) + 1
    IZ(1) = nint(RBin(3)/dbin) + 1 ! Note +1 because fortran starts a 1

    IX(2) = IX(1)-1
    IY(2) = IY(1)-1
    IZ(2) = IZ(1)-1

    ! Calculate the bin weighting
    WX(2) = (dbin*IX(1)-0.5_dp*dbin-RBin(1))/dbin
    WX(1) = 1.0_dp-WX(2)
    WY(2) = (dbin*IY(1)-0.5_dp*dbin-RBin(2))/dbin
    WY(1) = 1.0_dp-WY(2)
    WZ(2) = (dbin*IZ(1)-0.5_dp*dbin-RBin(3))/dbin
    WZ(1) = 1.0_dp-WZ(2)
CASE (3)
    ! Binning
    IX(1) = nint(RBin(1)/dbin) + 1
    IY(1) = nint(RBin(2)/dbin) + 1
    IZ(1) = nint(RBin(3)/dbin) + 1 ! Note +1 because fortran starts a 1

    IX(2) = IX(1)-1
    IY(2) = IY(1)-1
    IZ(2) = IZ(1)-1

    ! Calculate the bin weighting
    WX(2) = (dbin*IX(1)-0.5_dp*dbin-RBin(1))/dbin
    WX(1) = 1.0_dp-WX(2)
    WY(2) = (dbin*IY(1)-0.5_dp*dbin-RBin(2))/dbin
    WY(1) = 1.0_dp-WY(2)
    WZ(2) = (dbin*IZ(1)-0.5_dp*dbin-RBin(3))/dbin
    WZ(1) = 1.0_dp-WZ(2)
CASE (4) ! Box from 0-LBOX, Bins split by boundaries
    ! Periodic BC
    RBin(1) = RBin(1)-floor(RBin(1)/LBOX(1))*LBOX(1)
    RBin(2) = RBin(2)-floor(RBin(2)/LBOX(2))*LBOX(2)
    RBin(3) = RBin(3)-floor(RBin(3)/LBOX(3))*LBOX(3)

    ! Binning
    IX(1) = ceiling(RBin(1)/dbin)
    IY(1) = ceiling(RBin(2)/dbin)
    IZ(1) = ceiling(RBin(3)/dbin)

    IX(2) = IX(1)-1
    IY(2) = IY(1)-1
    IZ(2) = IZ(1)-1

    ! Calculate the bin weighting
    WX(2) = (dbin*IX(1)-RBin(1))/dbin
    WX(1) = 1.0_dp-WX(2)
    WY(2) = (dbin*IY(1)-RBin(2))/dbin
    WY(1) = 1.0_dp-WY(2)
    WZ(2) = (dbin*IZ(1)-RBin(3))/dbin
    WZ(1) = 1.0_dp-WZ(2)

    ! Periodic BC on Bins:
    IX(1) = IX(1)-floor(REAL((IX(1)-1))/REAL(NBinX(1))) * NBinX(1)
    IX(2) = IX(2)-floor(REAL((IX(2)-1))/REAL(NBinX(1))) * NBinX(1)
    IY(1) = IY(1)-floor(REAL((IY(1)-1))/REAL(NBinX(2))) * NBinX(2)
    IY(2) = IY(2)-floor(REAL((IY(2)-1))/REAL(NBinX(2))) * NBinX(2)
    IZ(1) = IZ(1)-floor(REAL((IZ(1)-1))/REAL(NBinX(3))) * NBinX(3)
    IZ(2) = IZ(2)-floor(REAL((IZ(2)-1))/REAL(NBinX(3))) * NBinX(3)
END SELECT
return
end subroutine
