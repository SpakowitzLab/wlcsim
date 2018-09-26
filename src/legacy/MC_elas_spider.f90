#include "../defines.inc"
!---------------------------------------------------------------*
!
!     Calculate the change in the polymer elastic energy
!     due to the displacement from a MC move
!
subroutine MC_eelas_spider(wlc_p,DEELAS,spider_id,EB,EPAR,EPERP,GAM,ETA)
! values from wlcsim_data
use params, only: wlc_spiders, wlc_U, wlc_R, wlc_UP, wlc_RP

use params, only: dp, wlcsim_params
use MC_wlc, only: E_wlc, E_SSWLC, E_GAUSS
implicit none
type(wlcsim_params), intent(in) :: wlc_p
integer, intent(in) :: spider_id
real(dp), intent(out) :: DEELAS(4)   ! Change in ECOM

!     Polymer properties

real(dp), intent(in) :: EB
real(dp), intent(in) :: EPAR
real(dp), intent(in) :: EPERP
real(dp), intent(in) :: GAM
real(dp), intent(in) :: ETA

integer hip, knee, toe, leg_n, ii, I
integer left_ends(3)  ! left end bead of the three kninks
real(dp) energy_change(4)

! Setup parameters

DEELAS(1) = 0.0_dp ! bending energy
DEELAS(2) = 0.0_dp ! parallel stretch energy
DEELAS(3) = 0.0_dp ! perpendicular stretch energy
DEELAS(4) = 0.0_dp ! WLC_P__TWIST energy

if (WLC_P__RING) then
    print*, "Ring not set up for spider move"
    stop 1
endif

!     Calculate the change in the energy
do leg_n = 1,wlc_spiders(spider_id)%nLegs
    hip = wlc_spiders(spider_id)%legs(1,leg_n)
    knee= wlc_spiders(spider_id)%legs(2,leg_n)
    toe = wlc_spiders(spider_id)%legs(3,leg_n)

    ! note that hip and toe U vectors don't rotate
    ! note also that the knee U vector rotates with the hip
    if (toe>hip) then
        left_ends(1) = hip 
        left_ends(2) = knee
        left_ends(3) = toe-1
    else
        left_ends(1) = toe
        left_ends(2) = knee-1
        left_ends(3) = hip-1
    endif
    do ii = 1,3
        I = left_ends(ii)
        ! MC move only affects energy if it's interior to the polymer, since there are only crankshaft moves, and end
        ! crankshafts don't affect polymer
        if (wlc_p%SIMTYPE == 1) then
            print*, "You will need to update this section before use."
            print*, "Finish implementing IT1 and IT2"
            stop 1
            !DEELAS(1) = E_wlc(RP(:,IT1M1), RP(:,IT1), RP(:,IT1P1), EB ) 
            !DEELAS(1) = DEELAS(1) - E_wlc(R(:,IT1M1), R(:,IT1), R(:,IT1P1), EB) 
        
        elseif (wlc_p%SIMTYPE == 2) then
        
            !function E_SSWLC(R,RM1,U,UM1,EB,EPAR,EPERP,ETA,GAM)
            energy_change = E_SSWLC(wlc_RP(:,I+1),wlc_RP(:,I),&
                                      wlc_UP(:,I+1),wlc_UP(:,I),&
                                      EB,EPAR,EPERP,ETA,GAM)
            DEELAS = DEELAS + energy_change
            energy_change = E_SSWLC( wlc_R(:,I+1), wlc_R(:,I),&
                                       wlc_U(:,I+1), wlc_U(:,I),&
                                       EB,EPAR,EPERP,ETA,GAM)
            DEELAS = DEELAS - energy_change
        
        elseif (wlc_p%SIMTYPE == 3) then
            DEELAS(2) = DEELAS(2) + E_GAUSS(wlc_RP(:,I+1),wlc_RP(:,I),EPAR)
            DEELAS(2) = DEELAS(2) - E_GAUSS( wlc_R(:,I+1), wlc_R(:,I),EPAR)
        endif
    enddo
enddo
RETURN
END

!---------------------------------------------------------------*
