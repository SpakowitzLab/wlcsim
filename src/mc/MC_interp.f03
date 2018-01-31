#include "../defines.inc"
!-----------------------------------------------------------------
!
!     This program linearly interpolates a bead at RBin into
!     8 bins indexed by IX, IY, IZ with weights WX, WY, WZ
!
!        Addapted from MC_int.f95
!        By Quinn MacPherson in summer 2016
!
!---------------------------------------------------------------
subroutine interp(wlc_p,RBin,IX,IY,IZ,WX,WY,WZ)
use params, only: dp, wlcsim_params
use inputparams, only: MAXPARAMLEN
implicit none
Type(wlcsim_params), intent(in) :: wlc_p      ! system varibles

real(dp), intent(inout) :: RBin(3) ! position or posiion within bin
integer, intent(out) :: IX(2)  ! Output
integer, intent(out) :: IY(2) ! Output
integer, intent(out) :: IZ(2)  ! Output
real(dp), intent(out) :: WX(2) ! Output
real(dp), intent(out) :: WY(2) ! Output
real(dp), intent(out) :: WZ(2) ! Output
SELECT CASE (WLC_P__CONFINETYPE)
CASE ('none') ! Box from 0-wlc_p%LBOX, Bins split by boundaries
    ! Periodic BC
    RBin(1) = RBin(1)-floor(RBin(1)/WLC_P__LBOX_X)*WLC_P__LBOX_X
    RBin(2) = RBin(2)-floor(RBin(2)/WLC_P__LBOX_Y)*WLC_P__LBOX_Y
    RBin(3) = RBin(3)-floor(RBin(3)/WLC_P__LBOX_Z)*WLC_P__LBOX_Z

    ! Binning
    IX(1) = ceiling(RBin(1)/WLC_P__DBIN)
    IY(1) = ceiling(RBin(2)/WLC_P__DBIN)
    IZ(1) = ceiling(RBin(3)/WLC_P__DBIN)

    IX(2) = IX(1)-1
    IY(2) = IY(1)-1
    IZ(2) = IZ(1)-1

    ! Calculate the bin weighting
    WX(2) = (WLC_P__DBIN*IX(1)-RBin(1))/WLC_P__DBIN
    WX(1) = 1.0_dp-WX(2)
    WY(2) = (WLC_P__DBIN*IY(1)-RBin(2))/WLC_P__DBIN
    WY(1) = 1.0_dp-WY(2)
    WZ(2) = (WLC_P__DBIN*IZ(1)-RBin(3))/WLC_P__DBIN
    WZ(1) = 1.0_dp-WZ(2)

    ! Periodic BC on Bins:
    IX(1) = IX(1)-floor(REAL((IX(1)-1))/REAL(wlc_p%NBINX(1))) * wlc_p%NBINX(1)
    IX(2) = IX(2)-floor(REAL((IX(2)-1))/REAL(wlc_p%NBINX(1))) * wlc_p%NBINX(1)
    IY(1) = IY(1)-floor(REAL((IY(1)-1))/REAL(wlc_p%NBINX(2))) * wlc_p%NBINX(2)
    IY(2) = IY(2)-floor(REAL((IY(2)-1))/REAL(wlc_p%NBINX(2))) * wlc_p%NBINX(2)
    IZ(1) = IZ(1)-floor(REAL((IZ(1)-1))/REAL(wlc_p%NBINX(3))) * wlc_p%NBINX(3)
    IZ(2) = IZ(2)-floor(REAL((IZ(2)-1))/REAL(wlc_p%NBINX(3))) * wlc_p%NBINX(3)
CASE ('platesInZperiodicXY')
    ! Periodic BC
    RBin(1) = RBin(1)-floor(RBin(1)/WLC_P__LBOX_Y)*WLC_P__LBOX_X
    RBin(2) = RBin(2)-floor(RBin(2)/WLC_P__LBOX_X)*WLC_P__LBOX_Y

    ! Binning
    IX(1) = ceiling(RBin(1)/WLC_P__DBIN)
    IY(1) = ceiling(RBin(2)/WLC_P__DBIN)
    IZ(1) = nint(RBin(3)/WLC_P__DBIN) + 1 ! Note 1.0 so that box centers are on half intigers

    IX(2) = IX(1)-1
    IY(2) = IY(1)-1
    IZ(2) = IZ(1)-1

    ! Calculate the bin weighting
    WX(2) = (WLC_P__DBIN*IX(1)-RBin(1))/WLC_P__DBIN   ! WX(2) = (RBin(1)-IX(1)*WLC_P__DBIN)/(IX(2)*WLC_P__DBIN-IX(1)*WLC_P__DBIN)
    WX(1) = 1.0_dp-WX(2)              ! WX(1) = (IX(2)*WLC_P__DBIN-RBin(1))/(IX(2)*WLC_P__DBIN-IX(1)*WLC_P__DBIN)
    WY(2) = (WLC_P__DBIN*IY(1)-RBin(2))/WLC_P__DBIN   ! WY(2) = (RBin(2)-IY(1)*WLC_P__DBIN)/(IY(2)*WLC_P__DBIN-IY(1)*WLC_P__DBIN)
    WY(1) = 1.0_dp-WY(2)              ! WY(1) = (IY(2)*WLC_P__DBIN-RBin(2))/(IY(2)*WLC_P__DBIN-IY(1)*WLC_P__DBIN)
    WZ(2) = (WLC_P__DBIN*IZ(1)-0.5_dp*WLC_P__DBIN-RBin(3))/WLC_P__DBIN   ! WZ(2) = (RBin(3)-IZ(1)*WLC_P__DBIN)/(IZ(2)*WLC_P__DBIN-IZ(1)*WLC_P__DBIN)
    WZ(1) = 1.0_dp-WZ(2)                   ! WZ(1) = (IZ(2)*WLC_P__DBIN-RBin(3))/(IZ(2)*WLC_P__DBIN-IZ(1)*WLC_P__DBIN)

    if ((WZ(1).lt.0).OR.(WZ(2).lt.0)) then
        print*, "negitive W"
        stop 1
    endif

    ! Periodic BC on Bins:
    IX(1) = IX(1)-floor(REAL((IX(1)-1))/REAL(wlc_p%NBINX(1))) * wlc_p%NBINX(1)
    IX(2) = IX(2)-floor(REAL((IX(2)-1))/REAL(wlc_p%NBINX(1))) * wlc_p%NBINX(1)
    IY(1) = IY(1)-floor(REAL((IY(1)-1))/REAL(wlc_p%NBINX(2))) * wlc_p%NBINX(2)
    IY(2) = IY(2)-floor(REAL((IY(2)-1))/REAL(wlc_p%NBINX(2))) * wlc_p%NBINX(2)
CASE ('cube') ! Box confinement
    ! Binning
    IX(1) = nint(RBin(1)/WLC_P__DBIN) + 1
    IY(1) = nint(RBin(2)/WLC_P__DBIN) + 1
    IZ(1) = nint(RBin(3)/WLC_P__DBIN) + 1 ! Note +1 because fortran starts a 1

    IX(2) = IX(1)-1
    IY(2) = IY(1)-1
    IZ(2) = IZ(1)-1

    ! Calculate the bin weighting
    WX(2) = (WLC_P__DBIN*IX(1)-0.5_dp*WLC_P__DBIN-RBin(1))/WLC_P__DBIN
    WX(1) = 1.0_dp-WX(2)
    WY(2) = (WLC_P__DBIN*IY(1)-0.5_dp*WLC_P__DBIN-RBin(2))/WLC_P__DBIN
    WY(1) = 1.0_dp-WY(2)
    WZ(2) = (WLC_P__DBIN*IZ(1)-0.5_dp*WLC_P__DBIN-RBin(3))/WLC_P__DBIN
    WZ(1) = 1.0_dp-WZ(2)
CASE ('sphere')
    ! Binning
    IX(1) = nint(RBin(1)/WLC_P__DBIN) + 1
    IY(1) = nint(RBin(2)/WLC_P__DBIN) + 1
    IZ(1) = nint(RBin(3)/WLC_P__DBIN) + 1 ! Note +1 because fortran starts a 1

    IX(2) = IX(1)-1
    IY(2) = IY(1)-1
    IZ(2) = IZ(1)-1

    ! Calculate the bin weighting
    WX(2) = (WLC_P__DBIN*IX(1)-0.5_dp*WLC_P__DBIN-RBin(1))/WLC_P__DBIN
    WX(1) = 1.0_dp-WX(2)
    WY(2) = (WLC_P__DBIN*IY(1)-0.5_dp*WLC_P__DBIN-RBin(2))/WLC_P__DBIN
    WY(1) = 1.0_dp-WY(2)
    WZ(2) = (WLC_P__DBIN*IZ(1)-0.5_dp*WLC_P__DBIN-RBin(3))/WLC_P__DBIN
    WZ(1) = 1.0_dp-WZ(2)
CASE ('periodicUnequal') ! Box from 0-wlc_p%LBOX, Bins split by boundaries
    ! Periodic BC
    RBin(1) = RBin(1)-floor(RBin(1)/WLC_P__LBOX_X)*WLC_P__LBOX_X
    RBin(2) = RBin(2)-floor(RBin(2)/WLC_P__LBOX_Y)*WLC_P__LBOX_Y
    RBin(3) = RBin(3)-floor(RBin(3)/WLC_P__LBOX_Z)*WLC_P__LBOX_Z

    ! Binning
    IX(1) = ceiling(RBin(1)/WLC_P__DBIN)
    IY(1) = ceiling(RBin(2)/WLC_P__DBIN)
    IZ(1) = ceiling(RBin(3)/WLC_P__DBIN)

    IX(2) = IX(1)-1
    IY(2) = IY(1)-1
    IZ(2) = IZ(1)-1

    ! Calculate the bin weighting
    WX(2) = (WLC_P__DBIN*IX(1)-RBin(1))/WLC_P__DBIN
    WX(1) = 1.0_dp-WX(2)
    WY(2) = (WLC_P__DBIN*IY(1)-RBin(2))/WLC_P__DBIN
    WY(1) = 1.0_dp-WY(2)
    WZ(2) = (WLC_P__DBIN*IZ(1)-RBin(3))/WLC_P__DBIN
    WZ(1) = 1.0_dp-WZ(2)

    ! Periodic BC on Bins:
    IX(1) = IX(1)-floor(REAL((IX(1)-1))/REAL(wlc_p%NBINX(1))) * wlc_p%NBINX(1)
    IX(2) = IX(2)-floor(REAL((IX(2)-1))/REAL(wlc_p%NBINX(1))) * wlc_p%NBINX(1)
    IY(1) = IY(1)-floor(REAL((IY(1)-1))/REAL(wlc_p%NBINX(2))) * wlc_p%NBINX(2)
    IY(2) = IY(2)-floor(REAL((IY(2)-1))/REAL(wlc_p%NBINX(2))) * wlc_p%NBINX(2)
    IZ(1) = IZ(1)-floor(REAL((IZ(1)-1))/REAL(wlc_p%NBINX(3))) * wlc_p%NBINX(3)
    IZ(2) = IZ(2)-floor(REAL((IZ(2)-1))/REAL(wlc_p%NBINX(3))) * wlc_p%NBINX(3)
END SELECT
return
end subroutine
