#include "../defines.inc"
!---------------------------------------------------------------*
!
!          For calculating total energy due to twist by subtracting
!          from writh.  From Brad's Code
!
!--------------------------------------------------------------
subroutine MC_global_twist(IT1,IT2,MCTYPE)
! values from wlcsim_data
use params, only: wlc_R, wlc_RP
use energies, only: energyOf, global_twistLiner_, global_twistquadratic_

use params, only: dp, pi, wlcsim_params, nan
implicit none
integer, intent(in) :: IT1
integer, intent(in) :: IT2

!     Polymer properties

real(dp) tw       ! Twist
real(dp) twP      ! Twist of test structure
real(dp) DWR      ! Change in Writhe
real(dp) WRM,WRMP ! Component of writhe affected by move
integer MCTYPE            ! MC move type

! Setup parameters

    if (MCTYPE == 1) then
        CALL WRITHECRANK(wlc_R,IT1,IT2,WLC_P__NB,WRM)
        CALL WRITHECRANK(wlc_RP,IT1,IT2,WLC_P__NB,WRMP)
        DWR = WRMP-WRM
    elseif (MCTYPE == 2) then
        CALL WRITHESLIDE(wlc_R,IT1,IT2,WLC_P__NB,WRM)
        CALL WRITHESLIDE(wlc_RP,IT1,IT2,WLC_P__NB,WRMP)
        DWR = WRMP-WRM
    else
        DWR = 0.0_dp
    ENDif
    energyOf(global_twistLiner_)%dx = DWR
    energyOf(global_twistQuadratic_)%dx = DWR**2 + 2*DWR*energyOf(global_twistLiner_)%x

    RETURN
END

!---------------------------------------------------------------*
