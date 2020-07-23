#include "../defines.inc"

subroutine hamiltonian(wlc_p, initialize)
!--------------------------------------------------------------------
! This subroutine calculates the field energy (Hamiltonian) based on
! volume fractions phi_A and phi_B.  The expresion used to calculate the
! energy is different for different systems.  The setting
! WLC_P__FIELDINTERACTIONTYPE determines which expression to use.
!
! If `initialize` is true then the energy for all bins is calculated.
! Otherwise calculate only for bins specified by wlc_inDPHI.
!-------------------------------------------------------------------

! values from wlcsim_data
   use params, only: wlc_NPHI, wlc_inDPHI &
                     , wlc_PHIB, wlc_PHIH, wlc_PHI_l2 &
                     , wlc_DPHIB, wlc_DPHI_l2 &
                     , wlc_Vol, wlc_DPHIA &
                     , wlc_PHIA, wlc_phiH_l2
   use energies, only: energyOf, chi_, couple_, kap_, field_, maierSaupe_
   use params, only: dp, wlcsim_params
   implicit none
   TYPE(wlcsim_params), intent(inout) :: wlc_p ! data
   logical, intent(in) :: initialize ! Need to do all beads?
   real(dp) PHIPoly ! fraction polymer
   real(dp) phi_A ! demsotu of A
   real(dp) phi_B ! density of B
   real(dp) phi_h ! strength of field
   real(dp) phi_l2 ! strength of field
   real(dp) VV ! volume of bin
   integer I, J, m_index ! for looping

   VV = WLC_P__DBIN**3
   if (initialize) then  ! calculate absolute energy

      select case (WLC_P__FIELDINTERACTIONTYPE) ! pick which keyword, case matchign string must be all uppercase

         !-------------------------------------------------------
         !
         !   Applied aligning field for Paul
         !
         !-------------------------------------------------------
      case ('AppliedAligningField')
         if (energyOf(maierSaupe_)%isOn) then
            do I = 1, wlc_p%NBIN
               if (WLC_P__FRACTIONAL_BIN) VV = wlc_Vol(I)
               energyOf(field_)%dx = energyOf(field_)%dx + VV*wlc_PHI_l2(0, I)
            enddo
         endif
         !-------------------------------------------------------
         !
         !   Applied aligning field for Luke
         !
         !-------------------------------------------------------
      case ('AppliedAligningFieldMelt')
         if (energyOf(maierSaupe_)%isOn) then
            do I = 1, wlc_p%NBIN
               if (WLC_P__FRACTIONAL_BIN) VV = wlc_Vol(I)
               do m_index = -2, 2
                  energyOf(field_)%dx = energyOf(field_)%dx + &
                                        VV*wlc_PHI_l2(m_index, I)*wlc_PHIH_l2(m_index, I)
               enddo
            enddo
         endif
         do I = 1, wlc_p%NBIN
            if (WLC_P__FRACTIONAL_BIN) then
               VV = wlc_Vol(I)
               if (VV .le. 0.1_dp) CYCLE
            endif
            energyOf(kap_)%dx = energyOf(kap_)%dx + VV*((wlc_PHIA(I) + wlc_PHIB(I) - 1.0_dp)**2)
         enddo
         !-------------------------------------------------------
         !
         !   Maier Saupe Melt Interaction
         !
         !-------------------------------------------------------
         ! In this problem Kap and chi are in units of kT/(simulation units cubed)
         ! If VV=1.0 than this is just kT/(bin volume)
      case ('MaierSaupe')
         if (energyOf(maierSaupe_)%isOn) then
            do I = 1, wlc_p%NBIN
               do m_index = -2, 2
                  if (WLC_P__FRACTIONAL_BIN) VV = wlc_Vol(I)
                  energyOf(maierSaupe_)%dx = energyOf(maierSaupe_)%dx + VV*wlc_PHI_l2(m_index, I)**2
               enddo
            enddo
         endif
         do I = 1, wlc_p%NBIN
            if (WLC_P__FRACTIONAL_BIN) then
               VV = wlc_Vol(I)
               if (VV .le. 0.1_dp) CYCLE
            endif
            energyOf(kap_)%dx = energyOf(kap_)%dx + VV*((wlc_PHIA(I) + wlc_PHIB(I) - 1.0_dp)**2)
         enddo
         !------------------------------------------------------------
         !
         !    A-B melt hamiltonian
         !
         !--------------------------------------------------------------
         ! Here Chi and Kap are in units of KT/beadVolume
      case ('ABmelt') ! Melt Hamiltonian
         do I = 1, wlc_p%NBIN
            if (WLC_P__FRACTIONAL_BIN) then
               VV = wlc_Vol(I)
               if (VV .le. 0.1_dp) CYCLE
            endif
            energyOf(chi_)%dx = energyOf(chi_)%dx + (VV/WLC_P__BEADVOLUME)*(wlc_PHIA(I)*wlc_PHIB(I))
            energyOf(kap_)%dx = energyOf(kap_)%dx + (VV/WLC_P__BEADVOLUME)*((wlc_PHIA(I) + wlc_PHIB(I) - 1.0_dp)**2)
            energyOf(field_)%dx = energyOf(field_)%dx - wlc_PHIH(I)*wlc_PHIA(I)
         enddo
         !------------------------------------------------------------
         !
         !    A-B solution hamiltonian
         !
         !--------------------------------------------------------------
         ! Here Chi and Kap are in units of KT/beadVolume
      case ('ABsolution') ! A,B,Solvent Hamiltonian
         do I = 1, wlc_p%NBIN
            if (WLC_P__FRACTIONAL_BIN) then
               VV = wlc_Vol(I)
               if (VV .le. 0.1_dp) CYCLE
            endif
            energyOf(chi_)%dx = energyOf(chi_)%dx + (VV/WLC_P__BEADVOLUME)*(wlc_PHIA(I)*wlc_PHIB(I))
            PHIPoly = wlc_PHIA(I) + wlc_PHIB(I)
            if (PHIPoly > 1.0_dp) then
               energyOf(kap_)%dx = energyOf(kap_)%dx + (VV/WLC_P__BEADVOLUME)*(PHIPoly - 1.0_dp)**2
            endif
         enddo
         !---------------------------------------------------------
         !
         !    Chromatin Hamiltonian
         !
         ! ---------------------------------------------------------
         ! Here Chi and Kap are in units of KT/beadVolume
      case ('chromatin')
         do I = 1, wlc_p%NBIN
            if (WLC_P__FRACTIONAL_BIN) then
               VV = wlc_Vol(I)
               if (VV .le. 0.1_dp) CYCLE
            endif
            PHIPoly = wlc_PHIA(I) + wlc_PHIB(I)
            energyOf(chi_)%dx = energyOf(chi_)%dx + (VV/WLC_P__BEADVOLUME)*PHIPoly*(1.0_dp - PHIPoly)
            energyOf(couple_)%dx = energyOf(couple_)%dx + VV*(wlc_PHIA(I))**2
            if (PHIPoly > 0.5_dp) then
               energyOf(kap_)%dx = energyOf(kap_)%dx + (VV/WLC_P__BEADVOLUME)*(PHIPoly - 0.5_dp)**2
            endif
         enddo
      case ('chromatin2')
         do I = 1, wlc_p%NBIN
            if (WLC_P__FRACTIONAL_BIN) then
               VV = wlc_Vol(I)
               if (VV .le. 0.1_dp) CYCLE
            endif
            energyOf(chi_)%dx = energyOf(chi_)%dx + VV*wlc_PHIB(I)*(1.0_dp - wlc_PHIB(I))
            energyOf(couple_)%dx = energyOf(couple_)%dx + VV*(wlc_PHIA(I))**2
            if (wlc_PHIB(I) > 0.5_dp) then
               energyOf(kap_)%dx = energyOf(kap_)%dx + VV*(wlc_PHIB(I) - 0.5_dp)**2
            endif
         enddo
      end select

   else ! Calculate change in energy

      select case (WLC_P__FIELDINTERACTIONTYPE) ! pick which keyword, case matchign string must be all uppercase
         !-------------------------------------------------------
         !
         !   Applied aligning field for Paul
         !
         !-------------------------------------------------------
      case ('AppliedAligningField')
         if (WLC_P__FRACTIONAL_BIN) VV = wlc_Vol(I)
         if (energyOf(maierSaupe_)%isOn) then
            do I = 1, wlc_NPHI
               energyOf(field_)%dx = energyOf(field_)%dx + VV*wlc_DPHI_l2(0, I)
            enddo
         endif
         !-------------------------------------------------------
         !
         !   Applied aligning field for Luke
         !
         !-------------------------------------------------------
      case ('AppliedAligningFieldMelt')
         do I = 1, wlc_NPHI
            J = wlc_inDPHI(I)
            if (energyOf(maierSaupe_)%isOn) then
               do m_index = -2, 2
                  energyOf(field_)%dx = energyOf(field_)%dx + &
                                        VV*wlc_DPHI_l2(m_index, I)*wlc_PHIH_l2(m_index, J)
               enddo
            endif
            phi_A = wlc_PHIA(J)
            phi_B = wlc_PHIB(J)
            energyOf(kap_)%dx = energyOf(kap_)%dx - VV*((phi_A + phi_B - 1.0_dp)**2)
            phi_A = phi_A + wlc_DPHIA(I)
            phi_B = phi_B + wlc_DPHIB(I)
            energyOf(kap_)%dx = energyOf(kap_)%dx + VV*((phi_A + phi_B - 1.0_dp)**2)
         enddo
         !-------------------------------------------------------
         !
         !   Maier Saupe Melt Interaction
         !
         !-------------------------------------------------------
         ! In this problem Kap and chi are in units of kT/binVolume
      case ('MaierSaupe')
         if (energyOf(maierSaupe_)%isOn) then
            do I = 1, wlc_NPHI
               if (WLC_P__FRACTIONAL_BIN) then
                  VV = wlc_Vol(I)
                  if (VV .le. 0.1_dp) CYCLE
               endif
               J = wlc_inDPHI(I)
               ! minus old
               phi_A = wlc_PHIA(J)
               phi_B = wlc_PHIB(J)
               energyOf(kap_)%dx = energyOf(kap_)%dx - VV*((phi_A + phi_B - 1.0_dp)**2)

               ! plus new
               phi_A = phi_A + wlc_DPHIA(I)
               phi_B = phi_B + wlc_DPHIB(I)
               energyOf(kap_)%dx = energyOf(kap_)%dx + VV*((phi_A + phi_B - 1.0_dp)**2)

               do m_index = -2, 2
                  ! minus old
                  phi_l2 = wlc_PHI_l2(m_index, J)
                  energyOf(maierSaupe_)%dx = energyOf(maierSaupe_)%dx - VV*phi_l2**2
                  ! plus new
                  phi_l2 = phi_l2 + wlc_DPHI_l2(m_index, I)
                  energyOf(maierSaupe_)%dx = energyOf(maierSaupe_)%dx + VV*phi_l2**2
               enddo
            enddo
         else
            do I = 1, wlc_NPHI
               if (WLC_P__FRACTIONAL_BIN) then
                  VV = wlc_Vol(I)
                  if (VV .le. 0.1_dp) CYCLE
               endif
               J = wlc_inDPHI(I)
               ! minus old
               phi_A = wlc_PHIA(J)
               phi_B = wlc_PHIB(J)
               energyOf(kap_)%dx = energyOf(kap_)%dx - VV*((phi_A + phi_B - 1.0_dp)**2)

               ! plus new
               phi_A = phi_A + wlc_DPHIA(I)
               phi_B = phi_B + wlc_DPHIB(I)
               energyOf(kap_)%dx = energyOf(kap_)%dx + VV*((phi_A + phi_B - 1.0_dp)**2)
            enddo
         endif

         !------------------------------------------------------------
         !
         !    A-B melt hamiltonian
         !
         !--------------------------------------------------------------
         ! Here Chi and Kap are in units of KT/beadVolume
      case ('ABmelt') ! Melt Hamiltonian
         do I = 1, wlc_NPHI
            J = wlc_inDPHI(I)
            if (WLC_P__FRACTIONAL_BIN) then
               VV = wlc_Vol(I)
               if (VV .le. 0.1_dp) CYCLE
            endif
            ! new
            phi_A = wlc_PHIA(J) + wlc_DPHIA(I)
            phi_B = wlc_PHIB(J) + wlc_DPHIB(I)
            phi_h = wlc_PHIH(J)
            energyOf(chi_)%dx = energyOf(chi_)%dx + (VV/WLC_P__BEADVOLUME)*phi_A*phi_B
            energyOf(kap_)%dx = energyOf(kap_)%dx + (VV/WLC_P__BEADVOLUME)*((phi_A + phi_B - 1.0_dp)**2)
            energyOf(field_)%dx = energyOf(field_)%dx - phi_h*phi_A
            ! minus old
            energyOf(chi_)%dx = energyOf(chi_)%dx - (VV/WLC_P__BEADVOLUME)*(wlc_PHIA(J)*wlc_PHIB(J))
            energyOf(kap_)%dx = energyOf(kap_)%dx - (VV/WLC_P__BEADVOLUME)*((wlc_PHIA(J) + wlc_PHIB(J) - 1.0_dp)**2)
            energyOf(field_)%dx = energyOf(field_)%dx + phi_h*wlc_PHIA(J)
         enddo
         !------------------------------------------------------------
         !
         !    A-B solution hamiltonian
         !
         !--------------------------------------------------------------
         ! Here Chi and Kap are in units of KT/beadVolume
      case ('ABsolution') ! Melt Hamiltonian
         do I = 1, wlc_NPHI
            J = wlc_inDPHI(I)
            if (WLC_P__FRACTIONAL_BIN) then
               VV = wlc_Vol(I)
               if (VV .le. 0.1_dp) CYCLE
            endif
            ! new
            phi_A = wlc_PHIA(J) + wlc_DPHIA(I)
            phi_B = wlc_PHIB(J) + wlc_DPHIB(I)
            phi_h = wlc_PHIH(J)
            energyOf(chi_)%dx = energyOf(chi_)%dx + (VV/WLC_P__BEADVOLUME)*phi_A*phi_B
            PHIPoly = phi_A + phi_B
            if (PHIPoly > 1.0_dp) then
               energyOf(kap_)%dx = energyOf(kap_)%dx + (VV/WLC_P__BEADVOLUME)*(PHIPoly - 1.0_dp)**2
            endif
            ! minus old
            phi_A = wlc_PHIA(J)
            phi_B = wlc_PHIB(J)
            PHIPoly = phi_A + phi_B
            if (PHIPoly > 1.0_dp) then
               energyOf(kap_)%dx = energyOf(kap_)%dx - (VV/WLC_P__BEADVOLUME)*(PHIPoly - 1.0_dp)**2
            endif
            energyOf(chi_)%dx = energyOf(chi_)%dx - (VV/WLC_P__BEADVOLUME)*(phi_A + phi_B)
         enddo
         !---------------------------------------------------------
         !
         !    Chromatin Hamiltonian
         !
         ! ---------------------------------------------------------
         ! Here Chi and Kap are in units of KT/beadVolume
      case ('chromatin')
         do I = 1, wlc_NPHI
            J = wlc_inDPHI(I)
            if (WLC_P__FRACTIONAL_BIN) then
               VV = wlc_Vol(I)
               if (VV .le. 0.1_dp) CYCLE
            endif
            ! new ...
            PHIPoly = wlc_PHIA(J) + wlc_DPHIA(I) + wlc_PHIB(J) + wlc_DPHIB(I)
            energyOf(chi_)%dx = energyOf(chi_)%dx + (VV/WLC_P__BEADVOLUME)*PHIPoly*(1.0_dp - PHIPoly)
            energyOf(couple_)%dx = energyOf(couple_)%dx + VV*(wlc_PHIA(J) + wlc_DPHIA(I))**2
            if (PHIPoly > 0.5_dp) then
               energyOf(kap_)%dx = energyOf(kap_)%dx + (VV/WLC_P__BEADVOLUME)*(PHIPoly - 0.5_dp)**2
            endif
            ! minus old
            PHIPoly = wlc_PHIA(J) + wlc_PHIB(J)
            energyOf(chi_)%dx = energyOf(chi_)%dx - (VV/WLC_P__BEADVOLUME)*PHIPoly*(1.0_dp - PHIPoly)
            energyOf(couple_)%dx = energyOf(couple_)%dx - VV*(wlc_PHIA(J))**2
            if (PHIPoly > 0.5_dp) then
               energyOf(kap_)%dx = energyOf(kap_)%dx - (VV/WLC_P__BEADVOLUME)*(PHIPoly - 0.5_dp)**2
            endif
         enddo
      case ('chromatin2')
         do I = 1, wlc_NPHI
            J = wlc_inDPHI(I)
            if (WLC_P__FRACTIONAL_BIN) then
               VV = wlc_Vol(I)
               if (VV .le. 0.1_dp) CYCLE
            endif
            ! new ...
            PHI_A = wlc_PHIA(J) + wlc_DPHIA(I)
            PHI_B = wlc_PHIB(J) + wlc_DPHIB(I)
            energyOf(chi_)%dx = energyOf(chi_)%dx + VV*PHI_B*(1.0_dp - PHI_B)
            energyOf(couple_)%dx = energyOf(couple_)%dx + VV*(PHI_A**2)
            if (PHI_B > 0.5_dp) then
               energyOf(kap_)%dx = energyOf(kap_)%dx + VV*(PHI_B - 0.5_dp)**2
            endif
            ! minus old
            PHI_A = wlc_PHIA(J)
            PHI_B = wlc_PHIB(J)
            energyOf(chi_)%dx = energyOf(chi_)%dx - VV*PHI_B*(1.0_dp - PHI_B)
            energyOf(couple_)%dx = energyOf(couple_)%dx - VV*(PHI_A**2)
            if (PHI_B > 0.5_dp) then
               energyOf(kap_)%dx = energyOf(kap_)%dx - VV*(PHI_B - 0.5_dp)**2
            endif
         enddo
      end select
   endif
   energyOf(chi_)%dx = energyOf(chi_)%dx
   energyOf(couple_)%dx = energyOf(couple_)%dx
   energyOf(kap_)%dx = energyOf(kap_)%dx
   RETURN
END subroutine

!---------------------------------------------------------------!
