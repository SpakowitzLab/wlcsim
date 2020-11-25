#include "../defines.inc"
subroutine mc_eelas(wlc_p)
!Calculate the change in the polymer elastic energy
!due to the displacement from a MC move.  Only calculate for
!bonds that have been changed which are listed in wlc_bendPints.

! values from wlcsim_data
   use params, only: wlc_U, wlc_nucleosomeWrap, wlc_VP, wlc_V &
                     , wlc_R, wlc_UP, wlc_basepairs, wlc_RP, wlc_bendPoints, wlc_nBend, &
                     wlc_basepairs_prop, wlc_nucleosomeWrap_prop

   use params, only: dp, wlcsim_params
   use MC_wlc, only: E_wlc, E_SSWLC, E_SSWLCWT, E_GAUSS
   use nucleosome, only: nucleosome_energy
   use polydispersity, only: is_right_end, leftmost_from, is_left_end, rightmost_from
   use energies, only: energyOf, bend_, stretch_, shear_, twist_

   implicit none
   type(wlcsim_params), intent(in) :: wlc_p

!     Polymer properties

   integer IT2
   integer IT2P1
   integer IT2M1
   real(dp) energy_change(4)
   real(dp) basepairs(WLC_P__NT)
   real(dp) wrapping(WLC_P__NT)
   integer ii

! Setup parameters
   energy_change = 0.0_dp

!     Calculate the change in the energy
   do ii = 1, wlc_nBend
      IT2 = wlc_bendPoints(ii)

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
         if (wlc_p%SIMTYPE == 1 .AND. (.not. is_left_end(IT2))) then
            if (is_left_end(IT2)) then
               IT2M1 = rightmost_from(IT2)
            else
               IT2M1 = IT2 - 1
            endif
            Print *, "This section is out of date"
            print *, "The variable IT2M1 is never used!"
            stop
            energyOf(bend_)%dx = energyOf(bend_)%dx - E_wlc(wlc_RP(:, IT2M1), wlc_RP(:, IT2), wlc_R(:, IT2P1), wlc_p%EB)
            energyOf(bend_)%dx = energyOf(bend_)%dx - E_wlc(wlc_R(:, IT2M1), wlc_R(:, IT2), wlc_R(:, IT2P1), wlc_p%EB)

         elseif (wlc_p%SIMTYPE == 2) then
            !function E_SSWLC(R,RM1,U,UM1,wlc_p%EB,wlc_p%EPAR,wlc_p%EPERP,wlc_p%ETA,wlc_p%GAM)
            if (WLC_P__LOCAL_TWIST) then
               energy_change = energy_change + E_SSWLCWT(wlc_RP(:, IT2P1), wlc_RP(:, IT2), &
                                                         wlc_UP(:, IT2P1), wlc_UP(:, IT2), &
                                                         wlc_VP(:, IT2P1), wlc_VP(:, IT2), &
                                                         wlc_p%EB, wlc_p%EPAR, wlc_p%EPERP, wlc_p%ETA, wlc_p%GAM, wlc_p%ETWIST)
               energy_change = energy_change - E_SSWLCWT(wlc_R(:, IT2P1), wlc_R(:, IT2), &
                                                         wlc_U(:, IT2P1), wlc_U(:, IT2), &
                                                         wlc_V(:, IT2P1), wlc_V(:, IT2), &
                                                         wlc_p%EB, wlc_p%EPAR, wlc_p%EPERP, wlc_p%ETA, wlc_p%GAM, wlc_p%ETWIST)
            else
               energy_change = energy_change + E_SSWLC(wlc_RP(:, IT2P1), wlc_RP(:, IT2), &
                                                       wlc_UP(:, IT2P1), wlc_UP(:, IT2), &
                                                       wlc_p%EB, wlc_p%EPAR, wlc_p%EPERP, wlc_p%ETA, wlc_p%GAM)
               energy_change = energy_change - E_SSWLC(wlc_R(:, IT2P1), wlc_R(:, IT2), &
                                                       wlc_U(:, IT2P1), wlc_U(:, IT2), &
                                                       wlc_p%EB, wlc_p%EPAR, wlc_p%EPERP, wlc_p%ETA, wlc_p%GAM)
            endif

         elseif (wlc_p%SIMTYPE == 3) then
            energyOf(stretch_)%dx = energyOf(stretch_)%dx + E_GAUSS(wlc_R(:, IT2P1), wlc_RP(:, IT2), wlc_p%EPAR)
            energyOf(stretch_)%dx = energyOf(stretch_)%dx - E_GAUSS(wlc_R(:, IT2P1), wlc_R(:, IT2), wlc_p%EPAR)
         endif
      elseif (WLC_P__ELASTICITY_TYPE == "nucleosomes") then
         ! set basepairs and wrapping vector
         basepairs = wlc_basepairs
         wrapping = wlc_nucleosomeWrap
         if (WLC_P__MOVEON_NUCLEOSOME_SLIDE == 1 .AND. WLC_P__MOVEON_NUCLEOSOME_BREATHE == 0) then 
            basepairs = wlc_basepairs_prop
         endif 
         if (WLC_P__MOVEON_NUCLEOSOME_BREATHE == 1) then 
            basepairs = wlc_basepairs_prop
            wrapping = wlc_nucleosomeWrap_prop
         endif 
         energy_change = energy_change + nucleosome_energy(wlc_RP(:, IT2P1), wlc_RP(:, IT2) &
                                                           , wlc_UP(:, IT2P1), wlc_UP(:, IT2) &
                                                           , wlc_VP(:, IT2P1), wlc_VP(:, IT2) &
                                                           , basepairs(IT2) &
                                                           , wrapping(IT2))
         energy_change = energy_change - nucleosome_energy(wlc_R(:, IT2P1), wlc_R(:, IT2) &
                                                           , wlc_U(:, IT2P1), wlc_U(:, IT2) &
                                                           , wlc_V(:, IT2P1), wlc_V(:, IT2) &
                                                           , wlc_basepairs(IT2) &
                                                           , wlc_nucleosomeWrap(IT2))
      endif
   enddo
   if (wlc_p%SIMTYPE == 2 .or. WLC_P__ELASTICITY_TYPE == "nucleosomes") then
      energyOf(bend_)%dx = energy_change(1)
      energyOf(stretch_)%dx = energy_change(2)
      energyOf(shear_)%dx = energy_change(3)
      energyOf(twist_)%dx = energy_change(4)
   endif
   RETURN
END

!---------------------------------------------------------------*
