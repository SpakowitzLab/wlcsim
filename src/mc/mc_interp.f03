#include "../defines.inc"
!-----------------------------------------------------------------
!
!        Addapted from MC_int.f95
!        By Quinn MacPherson in summer 2016
!
!---------------------------------------------------------------
subroutine interp(wlc_p,RBin,IX,IY,IZ,WX,WY,WZ)
!This program linearly interpolates a bead at RBin into
!8 bins indexed by IX, IY, IZ with weights WX, WY, WZ.
!Interpolation is based on Pike (2009), "Theoretically informed
!coarse grain simulation of polymeric systems".
!
!For exampe a if a bead is closer to the center of IX(1) then IX(2)
!then WX(1) will be proportionally higher then WX(2).  The total
!weight in any of the eight bins is given by WX*WY*WZ.

use params, only: dp, wlcsim_params
implicit none
Type(wlcsim_params), intent(in) :: wlc_p      ! system varibles

real(dp), intent(inout) :: RBin(3) ! position of bead
integer, intent(out) :: IX(2)  ! Range of bins in x direction
integer, intent(out) :: IY(2) ! Range of bins in y direction
integer, intent(out) :: IZ(2)  ! Range of bins in z direction
real(dp), intent(out) :: WX(2) ! Fractional contributions in x
real(dp), intent(out) :: WY(2) ! Fractional contributions in y
real(dp), intent(out) :: WZ(2) ! Fractional contributions in z
SELECT CASE (WLC_P__CONFINETYPE)
CASE ('none') ! Box from 0-wlc_p%LBOX, Bins split by boundaries
    ! Periodic BC
    RBin(1) = MODULO(RBin(1),WLC_P__LBOX_X)
    RBin(2) = MODULO(RBin(2),WLC_P__LBOX_Y)
    RBin(3) = MODULO(RBin(3),WLC_P__LBOX_Z)

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
    IX(1) = MODULO(IX(1)-1,WLC_P__NBIN_X) + 1
    IX(2) = MODULO(IX(2)-1,WLC_P__NBIN_X) + 1
    IY(1) = MODULO(IY(1)-1,WLC_P__NBIN_Y) + 1
    IY(2) = MODULO(IY(2)-1,WLC_P__NBIN_Y) + 1
    IZ(1) = MODULO(IZ(1)-1,WLC_P__NBIN_Z) + 1
    IZ(2) = MODULO(IZ(2)-1,WLC_P__NBIN_Z) + 1
CASE ('platesInZperiodicXY')
    ! Periodic BC
    RBin(1) = MODULO(RBin(1),WLC_P__LBOX_X)
    RBin(2) = MODULO(RBin(2),WLC_P__LBOX_Y)

    ! Binning
    IX(1) = ceiling(RBin(1)/WLC_P__DBIN)
    IY(1) = ceiling(RBin(2)/WLC_P__DBIN)

    IX(2) = IX(1)-1
    IY(2) = IY(1)-1

    ! Calculate the bin weighting
    WX(2) = (WLC_P__DBIN*IX(1)-RBin(1))/WLC_P__DBIN   ! WX(2) = (RBin(1)-IX(1)*WLC_P__DBIN)/(IX(2)*WLC_P__DBIN-IX(1)*WLC_P__DBIN)
    WX(1) = 1.0_dp-WX(2)              ! WX(1) = (IX(2)*WLC_P__DBIN-RBin(1))/(IX(2)*WLC_P__DBIN-IX(1)*WLC_P__DBIN)
    WY(2) = (WLC_P__DBIN*IY(1)-RBin(2))/WLC_P__DBIN   ! WY(2) = (RBin(2)-IY(1)*WLC_P__DBIN)/(IY(2)*WLC_P__DBIN-IY(1)*WLC_P__DBIN)
    WY(1) = 1.0_dp-WY(2)              ! WY(1) = (IY(2)*WLC_P__DBIN-RBin(2))/(IY(2)*WLC_P__DBIN-IY(1)*WLC_P__DBIN)

    ! Periodic BC on Bins:
    IX(1) = MODULO(IX(1)-1,WLC_P__NBIN_X) + 1
    IX(2) = MODULO(IX(2)-1,WLC_P__NBIN_X) + 1
    IY(1) = MODULO(IY(1)-1,WLC_P__NBIN_Y) + 1
    IY(2) = MODULO(IY(2)-1,WLC_P__NBIN_Y) + 1

    if (WLC_P__BOUNDARY_TYPE == "ExtendBinsPast") then
        IZ(1) = nint(RBin(3)/WLC_P__DBIN) + 1 ! Note 1.0 so that box centers are on half intigers
        IZ(2) = IZ(1)-1
        WZ(2) = (WLC_P__DBIN*IZ(1)-0.5_dp*WLC_P__DBIN-RBin(3))/WLC_P__DBIN
        WZ(1) = 1.0_dp-WZ(2)
    elseif (WLC_P__BOUNDARY_TYPE == "SolidEdgeBin") then
        if (RBin(3) <= WLC_P__DBIN*0.5_dp) then
            IZ(1) = 2
            IZ(2) = 1
            WZ(1) = 0.0_dp
            WZ(2) = 1.0_dp
        elseif (RBin(3) >= WLC_P__LBOX_Z - WLC_P__DBIN*0.5_dp) then
            IZ(2) = WLC_P__NBIN_Z - 1
            IZ(1) = WLC_P__NBIN_Z
            WZ(1) = 1.0_dp
            WZ(2) = 0.0_dp
        else
            IZ(1) = nint(RBin(3)/WLC_P__DBIN) + 1
            IZ(2) = IZ(1)-1
            WZ(2) = (WLC_P__DBIN*IZ(1)-0.5_dp*WLC_P__DBIN-RBin(3))/WLC_P__DBIN
            WZ(1) = 1.0_dp-WZ(2)
        endif
    endif
CASE ('cube') ! Box confinement
    if (WLC_P__BOUNDARY_TYPE == "ExtendBinsPast") then
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
    elseif (WLC_P__BOUNDARY_TYPE == "SolidEdgeBin") then
        if (RBin(1) <= WLC_P__DBIN*0.5_dp) then
            IX(1) = 2
            IX(2) = 1
            WX(1) = 0.0_dp
            WX(2) = 1.0_dp
        elseif (RBin(1) >= WLC_P__LBOX_X - WLC_P__DBIN*0.5_dp) then
            IX(2) = WLC_P__NBIN_X - 1
            IX(1) = WLC_P__NBIN_X
            WX(1) = 1.0_dp
            WX(2) = 0.0_dp
        else
            IX(1) = nint(RBin(1)/WLC_P__DBIN) + 1
            IX(2) = IX(1)-1
            WX(2) = (WLC_P__DBIN*IX(1)-0.5_dp*WLC_P__DBIN-RBin(1))/WLC_P__DBIN
            WX(1) = 1.0_dp-WX(2)
        endif
        !  -----
        if (RBin(2) <= WLC_P__DBIN*0.5_dp) then
            IY(1) = 2
            IY(2) = 1
            WY(1) = 0.0_dp
            WY(2) = 1.0_dp
        elseif (RBin(2) >= WLC_P__LBOX_Y - WLC_P__DBIN*0.5_dp) then
            IY(2) = WLC_P__NBIN_Y - 1
            IY(1) = WLC_P__NBIN_Y
            WY(1) = 1.0_dp
            WY(2) = 0.0_dp
        else
            IY(1) = nint(RBin(2)/WLC_P__DBIN) + 1
            IY(2) = IY(1)-1
            WY(2) = (WLC_P__DBIN*IY(1)-0.5_dp*WLC_P__DBIN-RBin(2))/WLC_P__DBIN
            WY(1) = 1.0_dp-WY(2)
        endif
        !  -----
        if (RBin(3) <= WLC_P__DBIN*0.5_dp) then
            IZ(1) = 2
            IZ(2) = 1
            WZ(1) = 0.0_dp
            WZ(2) = 1.0_dp
        elseif (RBin(3) >= WLC_P__LBOX_Z - WLC_P__DBIN*0.5_dp) then
            IZ(2) = WLC_P__NBIN_Z - 1
            IZ(1) = WLC_P__NBIN_Z
            WZ(1) = 1.0_dp
            WZ(2) = 0.0_dp
        else
            IZ(1) = nint(RBin(3)/WLC_P__DBIN) + 1
            IZ(2) = IZ(1)-1
            WZ(2) = (WLC_P__DBIN*IZ(1)-0.5_dp*WLC_P__DBIN-RBin(3))/WLC_P__DBIN
            WZ(1) = 1.0_dp-WZ(2)
        endif
    endif
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
    print*, "CONFINETYPE = periodicUnequal has be depreciated"
    print*, "Please use 'none'"
    stop 1
END SELECT
return
end subroutine
