#include "../defines.inc"
!---------------------------------------------------------------*

!     subroutine MC_eelas
!
!     Calculate the change in the polymer elastic energy
!     due to the displacement from a MC move
!
subroutine MC_eelas_spider(wlc_p,wlc_d,DEELAS,spider_id,EB,EPAR,EPERP,GAM,ETA)

use params, only: dp, pi,wlcsim_params, wlcsim_data
use MC_wlc, only: E_wlc, E_SSWLC, E_GAUSS
implicit none
type(wlcsim_params), intent(in) :: wlc_p
type(wlcsim_data), intent(in) :: wlc_d
integer, intent(in) :: spider_id
real(dp), intent(out) :: DEELAS(4)   ! Change in ECOM

!     Polymer properties

real(dp), intent(in) :: EB
real(dp), intent(in) :: EPAR
real(dp), intent(in) :: EPERP
real(dp), intent(in) :: GAM
real(dp), intent(in) :: ETA

integer hip, knee, toe, leg_n, ii, I, last_bead_on_polymer
integer left_ends(3)  ! left end bead of the three kninks

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
do leg_n = 1,wlc_d%spiders(spider_id)%nLegs
    hip = wlc_d%spiders(spider_id)%legs(1,leg_n)
    knee= wlc_d%spiders(spider_id)%legs(2,leg_n)
    toe = wlc_d%spiders(spider_id)%legs(3,leg_n)
    last_bead_on_polymer = (knee/WLC_P__NB + 1)*WLC_P__NB

    ! note that hip and toe U vectors don't rotate
    ! note also that the knee U vector rotates with the hip
    if (toe>hip) then
        left_ends = [hip,knee,toe-1]
    else
        left_ends = [toe,knee-1,hip-1]
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
            DEELAS = DEELAS + E_SSWLC(wlc_d%RP(:,I+1),wlc_d%RP(:,I),&
                                      wlc_d%UP(:,I+1),wlc_d%UP(:,I),&
                                      EB,EPAR,EPERP,ETA,GAM)
            DEELAS = DEELAS - E_SSWLC( wlc_d%R(:,I+1), wlc_d%R(:,I),&
                                       wlc_d%U(:,I+1), wlc_d%U(:,I),&
                                       EB,EPAR,EPERP,ETA,GAM)
        
        elseif (wlc_p%SIMTYPE == 3) then
            DEELAS(2) = DEELAS(2) + E_GAUSS(wlc_d%RP(:,I+1),wlc_d%RP(:,I),EPAR)
            DEELAS(2) = DEELAS(2) - E_GAUSS( wlc_d%R(:,I+1), wlc_d%R(:,I),EPAR)
        endif
    enddo
enddo
RETURN
END

!---------------------------------------------------------------*
