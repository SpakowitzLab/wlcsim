#include "../defines.inc"
!---------------------------------------------------------------*

!     subroutine MC_eelas
!
!     Calculate the change in the polymer elastic energy
!     due to the displacement from a MC move
!
subroutine MC_eelas(wlc_p,IB1,IB2,IT1,IT2,EB,EPAR,EPERP,GAM,ETA,MCTYPE,WRP)
! values from wlcsim_data
use params, only: wlc_U, wlc_nucleosomeWrap, wlc_VP, wlc_V, wlc_DEELAS&
    , wlc_R, wlc_UP, wlc_basepairs, wlc_RP, wlc_WR

use params, only: dp, pi,wlcsim_params
use MC_wlc, only: E_wlc, E_SSWLC, E_GAUSS
use nucleosome, only: nucleosome_energy
implicit none
type(wlcsim_params), intent(in) :: wlc_p
integer, intent(in) :: IB1               ! Test bead position 1
integer, intent(in) :: IT1               ! Index of test bead 1
integer, intent(in) :: IB2               ! Test bead position 2
integer, intent(in) :: IT2               ! Index of test bead 2

!     Polymer properties

real(dp), intent(in) :: EB
real(dp), intent(in) :: EPAR
real(dp), intent(in) :: EPERP
real(dp), intent(in) :: GAM
real(dp), intent(in) :: ETA
real(dp) tw       ! Twist
real(dp) twP      ! Twist of test structure
real(dp), intent(out) :: WRP      ! Writhe of test structure
real(dp) DWR      ! Change in Writhe
real(dp) WRM,WRMP ! Component of writhe affected by move
integer IT1M1
integer IT1P1
integer IT2M1
integer IT2P1
integer MCTYPE            ! MC move type
real(dp) energy_change(4)

! Setup parameters

wlc_DEELAS(1) = 0.0_dp ! bending energy
wlc_DEELAS(2) = 0.0_dp ! parallel stretch energy
wlc_DEELAS(3) = 0.0_dp ! perpendicular stretch energy
wlc_DEELAS(4) = 0.0_dp ! WLC_P__TWIST energy

!     Calculate the change in the energy

if ((IB1 /= 1).or.(WLC_P__RING)) then
    ! so if WLC_P__RING == 1, i.e.
    if (IB1 == 1) then
        ! then the bead to the left of IB1 is actually the end bead due to the WLC_P__RING inducing periodic boundaries on the array
        IT1M1 = WLC_P__NB
    else
        IT1M1 = IT1 - 1
    endif

    ! MC move only affects energy if it's interior to the polymer, since there are only crankshaft moves, and end
    ! crankshafts don't affect polymer
    if (WLC_P__ELASTICITY_TYPE == "constant") then
        if (wlc_p%SIMTYPE == 1.AND.(IB1 /= WLC_P__NB.OR.WLC_P__RING)) then
            if (IB1 == WLC_P__NB) then
                IT1P1 = 1
            else
                IT1P1 = IT1 + 1
            endif
            print*, "You will need to update this section before use."
            print*, "Finish implementing IT1 and IT2"
            stop 1
            wlc_DEELAS(1) = wlc_DEELAS(1) - E_wlc(wlc_R(:,IT1M1),wlc_RP(:,IT1),&
                                          wlc_RP(:,IT1P1), EB )
            wlc_DEELAS(1) = wlc_DEELAS(1) - E_wlc(wlc_R(:,IT1M1),wlc_R(:,IT1),&
                                          wlc_R(:,IT1P1), EB)

        elseif (wlc_p%SIMTYPE == 2) then
            energy_change = E_SSWLC(wlc_RP(:,IT1),wlc_R(:,IT1M1),&
                                      wlc_UP(:,IT1),wlc_U(:,IT1M1),&
                                      EB,EPAR,EPERP,ETA,GAM)
            wlc_DEELAS = wlc_DEELAS + energy_change
            energy_change =  E_SSWLC( wlc_R(:,IT1),wlc_R(:,IT1M1),&
                                       wlc_U(:,IT1),wlc_U(:,IT1M1)&
                                       ,EB,EPAR,EPERP,ETA,GAM)
            wlc_DEELAS = wlc_DEELAS - energy_change

        elseif (wlc_p%SIMTYPE == 3) then
            wlc_DEELAS(2) = wlc_DEELAS(2) + E_GAUSS(wlc_RP(:,IT1),wlc_R(:,IT1M1),EPAR)
            wlc_DEELAS(2) = wlc_DEELAS(2) - E_GAUSS( wlc_R(:,IT1),wlc_R(:,IT1M1),EPAR)
        endif
    elseif (WLC_P__ELASTICITY_TYPE == "nucleosomes") then
            energy_change = nucleosome_energy(wlc_RP(:,IT1),wlc_R(:,IT1M1)&
                                            ,wlc_UP(:,IT1),wlc_U(:,IT1M1)&
                                            ,wlc_VP(:,IT1),wlc_V(:,IT1M1)&
                                            ,wlc_basepairs(IT1M1)&
                                            ,wlc_nucleosomeWrap(IT1M1))
            wlc_DEELAS = wlc_DEELAS + energy_change
            energy_change = nucleosome_energy(wlc_R(:,IT1),wlc_R(:,IT1M1)&
                                            ,wlc_U(:,IT1),wlc_U(:,IT1M1)&
                                            ,wlc_V(:,IT1),wlc_V(:,IT1M1)&
                                            ,wlc_basepairs(IT1M1)&
                                            ,wlc_nucleosomeWrap(IT1M1))
            wlc_DEELAS = wlc_DEELAS - energy_change
    endif
endif

if ((IB2 /= WLC_P__NB).or.(WLC_P__RING)) then
    if (IB2 == WLC_P__NB) then
        IT2P1 = 1
    else
        IT2P1 = IT2 + 1
    endif

    ! if we're talking about a WLC, if we crankshaft a single bead, that's a no-op, since the u's are directly
    ! determined by the r's. Thus we're not worried about double counting the energy change here since the energy change
    ! should be zero by definition if IB1 = =IB2.
    if (WLC_P__ELASTICITY_TYPE == "constant") then
        if (wlc_p%SIMTYPE == 1.AND.((IB2 /= 1).OR.(WLC_P__RING))) then
            if (IB2 == 1) then
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
            energy_change = E_SSWLC(wlc_R(:,IT2P1),wlc_RP(:,IT2),&
                                                  wlc_U(:,IT2P1),wlc_UP(:,IT2),&
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
            energy_change = nucleosome_energy(wlc_R(:,IT2P1),wlc_RP(:,IT2)&
                                            ,wlc_U(:,IT2P1),wlc_UP(:,IT2)&
                                            ,wlc_V(:,IT2P1),wlc_VP(:,IT2)&
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
endif

if (WLC_P__RING.AND.WLC_P__TWIST) then
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
    wlc_DEELAS(4) = wlc_DEELAS(4) + &
                      (((2.0_dp*pi*twP)**2)*WLC_P__LT/(2.0_dp*WLC_P__L))&
                      -(((2.0_dp*pi*TW)**2)*WLC_P__LT/(2.0_dp*WLC_P__L))
ENDif

RETURN
END

!---------------------------------------------------------------*
