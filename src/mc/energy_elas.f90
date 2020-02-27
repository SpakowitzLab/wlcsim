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
      use nucleosome, only: nucleosome_energy, internucleosome_energy
      use polydispersity, only: first_bead_of_chain, length_of_chain
      implicit none
      type(wlcsim_params),intent(in) :: wlc_p
      real(dp), intent(out):: EELAS(5) ! Elastic force
      integer WR,TW ! writhe, twist
      integer I,J,IB,ibp1,jj            ! Index holders
      real(dp) energy_change(5)
      integer nucPlus ! how many nucs j is from i
      real(dp), parameter :: tau = 5.0_DP ! 0E distance between nucs (ideally this will be data to read in)


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
                IB = IB + 1
                CYCLE
            else
                IBP1 = IB + 1
            ENDif
            if (WLC_P__ELASTICITY_TYPE == "constant") then
                if (WLC_P__LOCAL_TWIST) then
                    energy_change(1:4) =  E_SSWLCWT( wlc_R(:,IBP1), wlc_R(:,IB),&
                                           wlc_U(:,IBP1), wlc_U(:,IB),&
                                           wlc_V(:,IBP1), wlc_V(:,IB),&
                                           wlc_p%EB, wlc_p%EPAR, wlc_p%EPERP,wlc_p%ETA, wlc_p%GAM, wlc_p%ETWIST)
                    EELAS(1:4) = EELAS(1:4) + energy_change(1:4)
                else
                    energy_change(1:4) =  E_SSWLC( wlc_R(:,IBP1), wlc_R(:,IB),&
                                           wlc_U(:,IBP1), wlc_U(:,IB),&
                                           wlc_p%EB, wlc_p%EPAR, wlc_p%EPERP,wlc_p%ETA, wlc_p%GAM)
                    EELAS(1:4) = EELAS(1:4) + energy_change(1:4)
                endif
            elseif (WLC_P__ELASTICITY_TYPE == "nucleosomes") then
                 energy_change(1:4) = nucleosome_energy(wlc_R(:,IBP1),wlc_R(:,IB)&
                                                  ,wlc_U(:,IBP1),wlc_U(:,IB)&
                                                  ,wlc_V(:,IBP1),wlc_V(:,IB)&
                                                  ,wlc_basepairs(IB)&
                                                  ,wlc_nucleosomeWrap(IB))
                 if (WLC_P__INTERNUCLEOSOME /= 0) then 
                    if (wlc_nucleosomeWrap(IB) /= 1 ) then 
                        nucPlus = 0
                        do jj = IB, WLC_P__NT
                            if (wlc_nucleosomeWrap(jj) /= 1) then
                                nucPlus = nucPlus + 1
                                ! will eventually want to use specific taus for i and j pairs but for now just look at next nearest
                                if (nucPlus == 2) then 
                                    energy_change(5) = energy_change(5) + &
                            internucleosome_energy(wlc_R(:,IB),wlc_R(:,jj),wlc_U(:,IB),wlc_U(:,jj),wlc_V(:,IB),wlc_V(:,jj),tau)
                                    exit
                                endif
                            endif
                        enddo
                    endif
                    EELAS(5) = EELAS(5) + energy_change(5)
                 endif
                 EELAS(1:4) = EELAS(1:4) + energy_change(1:4)
            endif
            IB = IB + 1
         ENDdo
      ENDdo

      ! Note Global twist is calculated elesweare, this only calculates local twist!

      RETURN
      END

!---------------------------------------------------------------*
