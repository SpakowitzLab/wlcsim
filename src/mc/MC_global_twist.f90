#include "../defines.inc"
!---------------------------------------------------------------*
!
!          For calculating total energy due to twist by subtracting
!          from writh.  From Brad's Code
!
!--------------------------------------------------------------
subroutine MC_global_twist(wlc_p,IT1,IT2,MCTYPE,WRP,DETWIST)
! values from wlcsim_data
use params, only: wlc_R, wlc_RP, wlc_WR

use params, only: dp, pi,wlcsim_params
implicit none
type(wlcsim_params), intent(in) :: wlc_p
integer, intent(in) :: IT1
integer, intent(in) :: IT2
real(dp), intent(out) :: DETWIST

!     Polymer properties

real(dp) tw       ! Twist
real(dp) twP      ! Twist of test structure
real(dp), intent(out) :: WRP      ! Writhe of test structure
real(dp) DWR      ! Change in Writhe
real(dp) WRM,WRMP ! Component of writhe affected by move
integer MCTYPE            ! MC move type

! Setup parameters

    DETWIST = 0.0_dp ! WLC_P__TWIST energy

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
    WRP = wlc_WR + DWR
    tw = REAL(wlc_p%LK,dp)-wlc_WR
    twP = REAL(wlc_p%LK,dp)-WRP
    DETWIST = DETWIST + &
                      (((2.0_dp*pi*twP)**2)*WLC_P__LT/(2.0_dp*WLC_P__L))&
                      -(((2.0_dp*pi*TW)**2)*WLC_P__LT/(2.0_dp*WLC_P__L))

    RETURN
END

!---------------------------------------------------------------*
