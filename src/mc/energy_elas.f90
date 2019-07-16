#include "../defines.inc"
!---------------------------------------------------------------*

!
!     This subroutine calculates the elastic energies for a wormlike
!     chain with a stretching potential. The stretch and bend
!     moduli are fed along with the bead positions.


      subroutine energy_elas(EELAS,wlc_p)
      ! values from wlcsim_data
      use params, only: wlc_nucleosomeWrap, wlc_basepairs, wlc_V, wlc_R, wlc_U
      use params, only: dp, pi,  wlcsim_params, nan
      use MC_wlc, only: E_SSWLC, E_SSWLCWT
      use nucleosome, only: nucleosome_energy
      use polydispersity, only: first_bead_of_chain, length_of_chain
      implicit none
      type(wlcsim_params),intent(in) :: wlc_p
      real(dp), intent(out):: EELAS(4) ! Elastic force
      integer WR,TW ! writhe, twist
      integer I,J,IB,ibp1            ! Index holders
      real(dp) energy_change(4)

      EELAS = 0.0_dp
      IB = 1
      do I = 1,WLC_P__NP
         do J = 1,length_of_chain(I)
            if (WLC_P__RING) then
                if (J == length_of_chain(I)) then
                    IBP1 = first_bead_of_chain(I)
                else
                    IBP1 = IB + 1
                ENDif
            elseif (J == length_of_chain(I)) then
                CYCLE
            else
                IBP1 = IB + 1
            ENDif
            if (WLC_P__ELASTICITY_TYPE == "constant") then
                if (WLC_P__LOCAL_TWIST) then
                    energy_change =  E_SSWLCWT( wlc_R(:,IBP1), wlc_R(:,IB),&
                                           wlc_U(:,IBP1), wlc_U(:,IB),&
                                           wlc_V(:,IBP1), wlc_V(:,IB),&
                                           wlc_p%EB, wlc_p%EPAR, wlc_p%EPERP,wlc_p%ETA, wlc_p%GAM, wlc_p%ETWIST)
                    EELAS = EELAS + energy_change
                else
                    energy_change =  E_SSWLC( wlc_R(:,IBP1), wlc_R(:,IB),&
                                           wlc_U(:,IBP1), wlc_U(:,IB),&
                                           wlc_p%EB, wlc_p%EPAR, wlc_p%EPERP,wlc_p%ETA, wlc_p%GAM)
                    EELAS = EELAS + energy_change
                endif
            elseif (WLC_P__ELASTICITY_TYPE == "nucleosomes") then
                 energy_change = nucleosome_energy(wlc_R(:,IBP1),wlc_R(:,IB)&
                                                  ,wlc_U(:,IBP1),wlc_U(:,IB)&
                                                  ,wlc_V(:,IBP1),wlc_V(:,IB)&
                                                  ,wlc_basepairs(IB)&
                                                  ,wlc_nucleosomeWrap(IB))
                 EELAS = EELAS + energy_change
            endif
            IB = IB + 1
         ENDdo
         IB = IB + 1
      ENDdo

      ! Note Global twist is calculated elesweare, this only calculates local twist!

      RETURN
      END

!---------------------------------------------------------------*
