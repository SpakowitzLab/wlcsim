#include "../defines.inc"
!---------------------------------------------------------------*

!
!     This subroutine calculates the elastic energies for a wormlike
!     chain with a stretching potential. The stretch and bend
!     moduli are fed along with the bead positions.


      subroutine energy_elas(EELAS,wlc_d,wlc_p)
      use params, only: dp, pi, wlcsim_data, wlcsim_params
      use MC_wlc, only: E_SSWLC
      use nucleosome, only: nucleosome_energy
      implicit none
      type(wlcsim_data), intent(in) :: wlc_d
      type(wlcsim_params),intent(in) :: wlc_p
      real(dp), intent(out):: EELAS(4) ! Elastic force
      integer WR,TW ! writhe, twist
      integer I,J,IB,ibp1            ! Index holders

      EELAS = 0.0_dp
      IB = 1
      do I = 1,WLC_P__NP
         do J = 1,WLC_P__NB
            if (WLC_P__RING) then
                if (J == WLC_P__NB) then
                    IBP1 = 1 + (I-1)*WLC_P__NB
                else
                    IBP1 = IB + 1
                ENDif
            elseif (J == WLC_P__NB) then
                CYCLE
            else
                IBP1 = IB + 1
            ENDif
            if (WLC_P__ELASTICITY_TYPE == "constant") then
                 EELAS = EELAS +  E_SSWLC( wlc_d%R(:,IBP1), wlc_d%R(:,IB),&
                                           wlc_d%U(:,IBP1), wlc_d%U(:,IB),&
                                           wlc_p%EB, wlc_p%EPAR, wlc_p%EPERP,wlc_p%ETA, wlc_p%GAM)
            elseif (WLC_P__ELASTICITY_TYPE == "nucleosomes") then
                 EELAS = EELAS + nucleosome_energy(wlc_d%R(:,IBP1),wlc_d%R(:,IB)&
                                                  ,wlc_d%U(:,IBP1),wlc_d%U(:,IB)&
                                                  ,wlc_d%V(:,IBP1),wlc_d%V(:,IB)&
                                                  ,wlc_d%basepairs(IB)&
                                                  ,wlc_d%nucleosomeWrap(IB))
            endif
            IB = IB + 1
         ENDdo
         IB = IB + 1
      ENDdo

      ! Get Twist Energy
      if (WLC_P__TWIST) then
          call WRITHE(wlc_d%R,WLC_P__NB,Wr)
          Tw = wlc_p%Lk-Wr
          EELAS(4) = ((2*PI*Tw)**2)*WLC_P__LT/(2*WLC_P__L)
      ENDif

      RETURN
      END

!---------------------------------------------------------------*
