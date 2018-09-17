#include "../defines.inc"
!---------------------------------------------------------------*

!     subroutine MC_eelas
!
!     Calculate the change in the polymer elastic energy
!     due to the displacement from a MC move
!
subroutine MC_eelas(wlc_p,EB,EPAR,EPERP,GAM,ETA)
! values from wlcsim_data
use params, only: wlc_U, wlc_nucleosomeWrap, wlc_VP, wlc_V, wlc_DEELAS&
    , wlc_R, wlc_UP, wlc_basepairs, wlc_RP, wlc_bendPoints, wlc_nBend

use params, only: dp, wlcsim_params
use MC_wlc, only: E_wlc, E_SSWLC, E_GAUSS
use nucleosome, only: nucleosome_energy
implicit none
type(wlcsim_params), intent(in) :: wlc_p

!     Polymer properties

real(dp), intent(in) :: EB
real(dp), intent(in) :: EPAR
real(dp), intent(in) :: EPERP
real(dp), intent(in) :: GAM
real(dp), intent(in) :: ETA
integer IT2
integer IT2P1
integer IT2M1
real(dp) energy_change(4)
integer ii

! Setup parameters

wlc_DEELAS(1) = 0.0_dp ! bending energy
wlc_DEELAS(2) = 0.0_dp ! parallel stretch energy
wlc_DEELAS(3) = 0.0_dp ! perpendicular stretch energy
wlc_DEELAS(4) = 0.0_dp ! WLC_P__TWIST energy

!     Calculate the change in the energy
do ii=1,wlc_nBend
    IT2=wlc_bendPoints(ii)

    if (WLC_P__RING .and. MOD(IT2,WLC_P__NB) == 0) then
        ! loop to beginning of polymer
        IT2P1 = (IT2/WLC_P__NB -1)*WLC_P__NB + 1
    else
        IT2P1 = IT2 + 1
    endif

    ! if we're talking about a WLC, if we crankshaft a single bead, that's a no-op, since the u's are directly
    ! determined by the r's. Thus we're not worried about double counting the energy change here since the energy change
    ! should be zero by definition if IB1 = =IB2.

    if (WLC_P__ELASTICITY_TYPE == "constant") then
        if (wlc_p%SIMTYPE == 1.AND.((MOD(IT2,WLC_P__NB) /= 1).OR.(WLC_P__RING))) then
            if (MOD(IT2,WLC_P__NB) == 1) then
                IT2M1 = WLC_P__NB
            else
                IT2M1 = IT2 - 1
            endif
            Print*, "This section is out of date"
            print*, "The variable IT2M1 is never used!"
            stop
            wlc_DEELAS(1) = wlc_DEELAS(1) - E_wlc(wlc_RP(:,IT2M1),wlc_RP(:,IT2),wlc_R(:,IT2P1),EB)
            wlc_DEELAS(1) = wlc_DEELAS(1) - E_wlc(wlc_R(:,IT2M1),wlc_R(:,IT2),wlc_R(:,IT2P1),EB)

        elseif (wlc_p%SIMTYPE == 2) then
            !function E_SSWLC(R,RM1,U,UM1,EB,EPAR,EPERP,ETA,GAM)
            energy_change = E_SSWLC(wlc_RP(:,IT2P1),wlc_RP(:,IT2),&
                                                  wlc_UP(:,IT2P1),wlc_UP(:,IT2),&
                                                  EB,EPAR,EPERP,ETA,GAM)
            wlc_DEELAS = wlc_DEELAS + energy_change
            energy_change = E_SSWLC(wlc_R(:,IT2P1), wlc_R(:,IT2),&
                                                  wlc_U(:,IT2P1), wlc_U(:,IT2),&
                                                  EB,EPAR,EPERP,ETA,GAM)
            wlc_DEELAS = wlc_DEELAS - energy_change

        elseif (wlc_p%SIMTYPE == 3) then
            wlc_DEELAS(2) = wlc_DEELAS(2) + E_GAUSS(wlc_R(:,IT2P1),wlc_RP(:,IT2),EPAR)
            wlc_DEELAS(2) = wlc_DEELAS(2) - E_GAUSS(wlc_R(:,IT2P1), wlc_R(:,IT2),EPAR)
        endif
    elseif (WLC_P__ELASTICITY_TYPE == "nucleosomes") then
            energy_change = nucleosome_energy(wlc_RP(:,IT2P1),wlc_RP(:,IT2)&
                                            ,wlc_UP(:,IT2P1),wlc_UP(:,IT2)&
                                            ,wlc_VP(:,IT2P1),wlc_VP(:,IT2)&
                                            ,wlc_basepairs(IT2)&
                                            ,wlc_nucleosomeWrap(IT2))
            wlc_DEELAS = wlc_DEELAS + energy_change
            energy_change = nucleosome_energy(wlc_R(:,IT2P1),wlc_R(:,IT2)&
                                            ,wlc_U(:,IT2P1),wlc_U(:,IT2)&
                                            ,wlc_V(:,IT2P1),wlc_V(:,IT2)&
                                            ,wlc_basepairs(IT2)&
                                            ,wlc_nucleosomeWrap(IT2))
            wlc_DEELAS = wlc_DEELAS - energy_change
    endif
enddo


RETURN
END

!---------------------------------------------------------------*
