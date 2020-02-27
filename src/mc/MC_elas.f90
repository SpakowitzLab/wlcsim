#include "../defines.inc"
!---------------------------------------------------------------*

!     subroutine MC_eelas
!
!     Calculate the change in the polymer elastic energy
!     due to the displacement from a MC move
!
subroutine MC_eelas(wlc_p)
! values from wlcsim_data
use params, only: wlc_U, wlc_nucleosomeWrap, wlc_VP, wlc_V&
    , wlc_R, wlc_UP, wlc_basepairs, wlc_RP, wlc_bendPoints, wlc_nBend

use params, only: dp, wlcsim_params, NAN
use MC_wlc, only: E_wlc, E_SSWLC, E_SSWLCWT, E_GAUSS
use nucleosome, only: nucleosome_energy, internucleosome_energy
use polydispersity, only: is_right_end, leftmost_from, is_left_end, rightmost_from
use energies, only: energyOf, bend_, stretch_, shear_, twist_, internucleosome_

implicit none
type(wlcsim_params), intent(in) :: wlc_p

!     Polymer properties

integer IT2
integer IT2P1
integer IT2M1
real(dp) energy_change(5)
real(dp) interNucEnergy_change(1)
integer ii, j
real(dp), dimension(3) :: tempR, tempU, tempV ! vectors of j nucleosomes
integer nucPlus ! how many nucs j is from i
real(dp), parameter :: tau = 5.0_DP ! 0E distance between nucs (ideally this will be data to read in)

! Setup parameters
energy_change = 0.0_dp

!     Calculate the change in the energy
do ii=1,wlc_nBend
    IT2=wlc_bendPoints(ii)

    if (WLC_P__RING .and. is_right_end(IT2)) then
        ! loop to beginning of polymer
        IT2P1 = leftmost_from(IT2)
    else
        IT2P1 = IT2 + 1
    endif

    ! if we're talking about a WLC, if we crankshaft a single bead, that's a no-op, since the u's are directly
    ! determined by the r's. Thus we're not worried about double counting the energy change here since the energy change
    ! should be zero by definition if IB1 = =IB2.

    if (WLC_P__ELASTICITY_TYPE == "constant") then
        if (wlc_p%SIMTYPE == 1.AND.(.not. is_left_end(IT2))) then
            if (is_left_end(IT2)) then
                IT2M1 = rightmost_from(IT2)
            else
                IT2M1 = IT2 - 1
            endif
            Print*, "This section is out of date"
            print*, "The variable IT2M1 is never used!"
            stop
            energyOf(bend_)%dx = energyOf(bend_)%dx - E_wlc(wlc_RP(:,IT2M1),wlc_RP(:,IT2),wlc_R(:,IT2P1),wlc_p%EB)
            energyOf(bend_)%dx = energyOf(bend_)%dx - E_wlc(wlc_R(:,IT2M1),wlc_R(:,IT2),wlc_R(:,IT2P1),wlc_p%EB)

        elseif (wlc_p%SIMTYPE == 2) then
            !function E_SSWLC(R,RM1,U,UM1,wlc_p%EB,wlc_p%EPAR,wlc_p%EPERP,wlc_p%ETA,wlc_p%GAM)
            if (WLC_P__LOCAL_TWIST) then
                energy_change(1:4) = energy_change(1:4) + E_SSWLCWT(wlc_RP(:,IT2P1),wlc_RP(:,IT2),&
                                                  wlc_UP(:,IT2P1),wlc_UP(:,IT2),&
                                                  wlc_VP(:,IT2P1),wlc_VP(:,IT2),&
                                                  wlc_p%EB,wlc_p%EPAR,wlc_p%EPERP,wlc_p%ETA,wlc_p%GAM, wlc_p%ETWIST)
                energy_change(1:4) = energy_change(1:4) - E_SSWLCWT(wlc_R(:,IT2P1), wlc_R(:,IT2),&
                                                  wlc_U(:,IT2P1), wlc_U(:,IT2),&
                                                  wlc_V(:,IT2P1), wlc_V(:,IT2),&
                                                  wlc_p%EB,wlc_p%EPAR,wlc_p%EPERP,wlc_p%ETA,wlc_p%GAM, wlc_p%ETWIST)
            else
                energy_change(1:4) = energy_change(1:4) + E_SSWLC(wlc_RP(:,IT2P1),wlc_RP(:,IT2),&
                                                  wlc_UP(:,IT2P1),wlc_UP(:,IT2),&
                                                  wlc_p%EB,wlc_p%EPAR,wlc_p%EPERP,wlc_p%ETA,wlc_p%GAM)
                energy_change(1:4) = energy_change(1:4) - E_SSWLC(wlc_R(:,IT2P1), wlc_R(:,IT2),&
                                                  wlc_U(:,IT2P1), wlc_U(:,IT2),&
                                                  wlc_p%EB,wlc_p%EPAR,wlc_p%EPERP,wlc_p%ETA,wlc_p%GAM)
            endif

        elseif (wlc_p%SIMTYPE == 3) then
            energyOf(stretch_)%dx = energyOf(stretch_)%dx + E_GAUSS(wlc_R(:,IT2P1),wlc_RP(:,IT2),wlc_p%EPAR)
            energyOf(stretch_)%dx = energyOf(stretch_)%dx - E_GAUSS(wlc_R(:,IT2P1), wlc_R(:,IT2),wlc_p%EPAR)
        endif
    elseif (WLC_P__ELASTICITY_TYPE == "nucleosomes") then
            energy_change(1:4) = energy_change(1:4) + nucleosome_energy(wlc_RP(:,IT2P1),wlc_RP(:,IT2)&
                                            ,wlc_UP(:,IT2P1),wlc_UP(:,IT2)&
                                            ,wlc_VP(:,IT2P1),wlc_VP(:,IT2)&
                                            ,wlc_basepairs(IT2)&
                                            ,wlc_nucleosomeWrap(IT2))
            energy_change(1:4) = energy_change(1:4) - nucleosome_energy(wlc_R(:,IT2P1),wlc_R(:,IT2)&
                                            ,wlc_U(:,IT2P1),wlc_U(:,IT2)&
                                            ,wlc_V(:,IT2P1),wlc_V(:,IT2)&
                                            ,wlc_basepairs(IT2)&
                                            ,wlc_nucleosomeWrap(IT2))
            if (WLC_P__INTERNUCLEOSOME /= 0) then 
                if (wlc_nucleosomeWrap(IT2) /= 1 ) then 
                    nucPlus = 0
                    do j = IT2, WLC_P__NT
                        if (wlc_nucleosomeWrap(j) /= 1) then
                            nucPlus = nucPlus + 1
                            ! will eventually want to use specific taus for i and j pairs but for now just look at next nearest
                            if (nucPlus == 2) then 
                                if (isnan(wlc_RP(1,j)) .eqv. .FALSE.) then 
                                    tempR = wlc_RP(:,j)
                                    tempU = wlc_UP(:,j)
                                    tempV = wlc_VP(:,j)
                                else
                                    tempR = wlc_R(:,j)
                                    tempU = wlc_U(:,j)
                                    tempV = wlc_V(:,j)
                                endif
                                energy_change(5) = energy_change(5) + &
                        internucleosome_energy(wlc_RP(:,IT2),tempR,wlc_UP(:,IT2),tempU,wlc_VP(:,IT2),tempV,tau)
                                energy_change(5) = energy_change(5) - &
                        internucleosome_energy(wlc_R(:,IT2),wlc_R(:,j),wlc_U(:,IT2),wlc_U(:,j),wlc_V(:,IT2),wlc_V(:,j),tau)
                                exit
                            endif
                        endif
                    enddo
                endif
            endif
    endif
enddo
if (wlc_p%SIMTYPE == 2 .or. WLC_P__ELASTICITY_TYPE == "nucleosomes") then
    energyOf(bend_)%dx = energy_change(1)
    energyOf(stretch_)%dx = energy_change(2)
    energyOf(shear_)%dx = energy_change(3)
    energyOf(twist_)%dx = energy_change(4)
    if (WLC_P__INTERNUCLEOSOME /= 0) then 
        energyOf(internucleosome_)%dx = energy_change(5)
    endif
endif
RETURN
END

!---------------------------------------------------------------*
