#include "../defines.inc"
!--------------------------------------------------------------------
!
!
! This subroutine calculates the field Hamiltonian from the phi values.
! It puts the output in: wlc_dx_chi, wlc_dx_chi, wlc_dx_couple, wlc_dx_couple, wlc_dx_Kap, wlc_dx_Kap, wlc_DEChi,wlc_dx_chi, wlc_DECouple,
! wlc_dx_couple, wlc_DEKap, wlc_dx_Kap, wlc_DEField, wlc_dx_Field
!      by Quinn MacPherson based on code from Shifan Mao
!       Made a separate function on 7/8/16
!
!   If initialize then calculate all bins.
!   Otherwise calcualte only specified bins.
!-------------------------------------------------------------------

subroutine hamiltonian(wlc_p,initialize)
! values from wlcsim_data
use params, only: wlc_deMaierSaupe, wlc_dx_Field, wlc_NPHI, wlc_inDPHI, wlc_DEKap&
    , wlc_PHIB, wlc_PHIH, wlc_Dx_Field, wlc_dx_couple, wlc_PHI_l2, wlc_Dx_Couple&
    , wlc_DPHIB, wlc_DPHI_l2, wlc_dx_Kap, wlc_DEField, wlc_dx_Chi, wlc_DEChi&
    , wlc_Vol, wlc_Dx_Chi, wlc_DECouple, wlc_DPHIA, wlc_dx_maierSaupe, wlc_dx_chi&
    , wlc_Dx_Kap, wlc_PHIA, wlc_phiH_l2
use params,only: dp,wlcsim_params
implicit none
TYPE(wlcsim_params), intent(inout) :: wlc_p
logical, intent(in) :: initialize ! Need to do all beads
real(dp) PHIPoly ! fraction polymer
real(dp) phi_A ! demsotu of A
real(dp) phi_B ! density of B
real(dp) phi_h ! strength of field
real(dp) phi_l2 ! strength of field
real(dp) VV ! volume of bin
integer I,J,m_index ! for looping

wlc_dx_Chi = 0.0_dp
wlc_Dx_Couple = 0.0_dp
wlc_Dx_Kap = 0.0_dp
wlc_Dx_Field = 0.0_dp
wlc_dx_maierSaupe = 0.0_dp
VV = WLC_P__DBIN**3
if (initialize) then  ! calculate absolute energy

    select case(WLC_P__FIELDINTERACTIONTYPE) ! pick which keyword, case matchign string must be all uppercase

    !-------------------------------------------------------
    !
    !   Applied aligning field for Paul
    !
    !-------------------------------------------------------
    case('AppliedAligningField')
        if (wlc_p%CHI_L2_ON) then
            do I = 1,wlc_p%NBIN
                if (WLC_P__FRACTIONAL_BIN) VV = wlc_Vol(I)
                wlc_dx_Field =  wlc_dx_Field + VV*wlc_PHI_l2(0,I)
            enddo
        endif
    !-------------------------------------------------------
    !
    !   Applied aligning field for Luke
    !
    !-------------------------------------------------------
    case('AppliedAligningFieldMelt')
        if (wlc_p%CHI_L2_ON) then
            do I = 1,wlc_p%NBIN
                if (WLC_P__FRACTIONAL_BIN) VV = wlc_Vol(I)
                do m_index = -2,2
                    wlc_dx_Field = wlc_dx_Field + &
                                        VV*wlc_PHI_l2(m_index,I)*wlc_PHIH_l2(m_index,I)
                enddo
            enddo
        endif
        do I = 1,wlc_p%NBIN
            if (WLC_P__FRACTIONAL_BIN) then
                VV = wlc_Vol(I)
                if (VV.le.0.1_dp) CYCLE
            endif
            wlc_Dx_Kap = wlc_dx_Kap + VV*((wlc_PHIA(I) + wlc_PHIB(I)-1.0_dp)**2)
        enddo
    !-------------------------------------------------------
    !
    !   Maier Saupe Melt Interaction
    !
    !-------------------------------------------------------
    ! In this problem Kap and chi are in units of kT/(simulation units cubed)
    ! If VV=1.0 than this is just kT/(bin volume)
    case('MaierSaupe')
        if (wlc_p%CHI_L2_ON) then
            do I = 1,wlc_p%NBIN
                do m_index = -2,2
                    if (WLC_P__FRACTIONAL_BIN) VV = wlc_Vol(I)
                    wlc_dx_maierSaupe =  wlc_dx_maierSaupe + VV*wlc_PHI_l2(m_index,I)**2
                enddo
            enddo
        endif
        do I = 1,wlc_p%NBIN
            if (WLC_P__FRACTIONAL_BIN) then
                VV = wlc_Vol(I)
                if (VV.le.0.1_dp) CYCLE
            endif
            wlc_Dx_Kap = wlc_dx_Kap + VV*((wlc_PHIA(I) + wlc_PHIB(I)-1.0_dp)**2)
        enddo
    !------------------------------------------------------------
    !
    !    A-B melt hamiltonian
    !
    !--------------------------------------------------------------
    ! Here Chi and Kap are in units of KT/beadVolume
    case('ABmelt') ! Melt Hamiltonian
        do I = 1,wlc_p%NBIN
            if (WLC_P__FRACTIONAL_BIN) then
                VV = wlc_Vol(I)
                if (VV.le.0.1_dp) CYCLE
            endif
            wlc_Dx_Chi = wlc_Dx_Chi + (VV/WLC_P__BEADVOLUME)*(wlc_PHIA(I)*wlc_PHIB(I))
            wlc_Dx_Kap = wlc_dx_Kap + (VV/WLC_P__BEADVOLUME)*((wlc_PHIA(I) + wlc_PHIB(I)-1.0_dp)**2)
            wlc_Dx_Field = wlc_dx_Field-wlc_PHIH(I)*wlc_PHIA(I)
        enddo
    !------------------------------------------------------------
    !
    !    A-B solution hamiltonian
    !
    !--------------------------------------------------------------
    ! Here Chi and Kap are in units of KT/beadVolume
    case('ABsolution') ! A,B,Solvent Hamiltonian
        do I = 1,wlc_p%NBIN
            if (WLC_P__FRACTIONAL_BIN) then
                VV = wlc_Vol(I)
                if (VV.le.0.1_dp) CYCLE
            endif
            wlc_Dx_Chi = wlc_Dx_Chi + (VV/WLC_P__BEADVOLUME)*(wlc_PHIA(I)*wlc_PHIB(I))
            PHIPoly = wlc_PHIA(I) + wlc_PHIB(I)
            if(PHIPoly > 1.0_dp) then
                wlc_Dx_Kap = wlc_Dx_Kap + (VV/WLC_P__BEADVOLUME)*(PHIPoly-1.0_dp)**2
            endif
        enddo
    !---------------------------------------------------------
    !
    !    Chromatin Hamiltonian
    !
    ! ---------------------------------------------------------
    ! Here Chi and Kap are in units of KT/beadVolume
    case('chromatin')
        do I = 1,wlc_p%NBIN
            if (WLC_P__FRACTIONAL_BIN) then
                VV = wlc_Vol(I)
                if (VV.le.0.1_dp) CYCLE
            endif
            PHIPoly = wlc_PHIA(I) + wlc_PHIB(I)
            wlc_Dx_Chi = wlc_Dx_Chi + (VV/WLC_P__BEADVOLUME)*PHIPoly*(1.0_dp-PHIPoly)
            wlc_Dx_Couple = wlc_Dx_Couple + VV*(wlc_PHIA(I))**2
            if(PHIPoly > 0.5_dp) then
                wlc_Dx_Kap = wlc_Dx_Kap + (VV/WLC_P__BEADVOLUME)*(PHIPoly-0.5_dp)**2
            endif
        enddo
    case('chromatin2')
        do I = 1,wlc_p%NBIN
            if (WLC_P__FRACTIONAL_BIN) then
                VV = wlc_Vol(I)
                if (VV.le.0.1_dp) CYCLE
            endif
            wlc_Dx_Chi = wlc_Dx_Chi + VV*wlc_PHIB(I)*(1.0_dp-wlc_PHIB(I))
            wlc_Dx_Couple = wlc_Dx_Couple + VV*(wlc_PHIA(I))**2
            if(wlc_PHIB(I) > 0.5_dp) then
                wlc_Dx_Kap = wlc_Dx_Kap + VV*(wlc_PHIB(I) - 0.5_dp)**2
            endif
        enddo
    end select


else ! Calculate change in energy

    select case(WLC_P__FIELDINTERACTIONTYPE) ! pick which keyword, case matchign string must be all uppercase
    !-------------------------------------------------------
    !
    !   Applied aligning field for Paul
    !
    !-------------------------------------------------------
    case('AppliedAligningField')
        if (WLC_P__FRACTIONAL_BIN) VV = wlc_Vol(I)
        if (wlc_p%CHI_L2_ON) then
            do I = 1,wlc_NPHI
                wlc_dx_Field =  wlc_dx_Field + VV*wlc_DPHI_l2(0,I)
            enddo
        endif
    !-------------------------------------------------------
    !
    !   Applied aligning field for Luke
    !
    !-------------------------------------------------------
    case('AppliedAligningFieldMelt')
        do I = 1,wlc_NPHI
            if (wlc_p%CHI_L2_ON) then
                do m_index = -2,2
                    wlc_dx_Field = wlc_dx_Field + &
                                        VV*wlc_DPHI_l2(m_index,I)*wlc_PHIH_l2(m_index,I)
                enddo
            endif
            phi_A = wlc_PHIA(J)
            phi_B = wlc_PHIB(J)
            wlc_Dx_Kap = wlc_Dx_Kap - VV*((phi_A + phi_B-1.0_dp)**2)
            phi_A = phi_A + wlc_DPHIA(I)
            phi_B = phi_B + wlc_DPHIB(I)
            wlc_Dx_Kap = wlc_Dx_Kap + VV*((phi_A + phi_B-1.0_dp)**2)
        enddo
    !-------------------------------------------------------
    !
    !   Maier Saupe Melt Interaction
    !
    !-------------------------------------------------------
    ! In this problem Kap and chi are in units of kT/binVolume
    case('MaierSaupe')
        if (wlc_p%CHI_L2_ON) then
            do I = 1,wlc_NPHI
                if (WLC_P__FRACTIONAL_BIN) then
                    VV = wlc_Vol(I)
                    if (VV.le.0.1_dp) CYCLE
                endif
                J = wlc_inDPHI(I)
                ! minus old
                phi_A = wlc_PHIA(J)
                phi_B = wlc_PHIB(J)
                wlc_Dx_Kap = wlc_Dx_Kap - VV*((phi_A + phi_B-1.0_dp)**2)


                ! plus new
                phi_A = phi_A + wlc_DPHIA(I)
                phi_B = phi_B + wlc_DPHIB(I)
                wlc_Dx_Kap = wlc_Dx_Kap + VV*((phi_A + phi_B-1.0_dp)**2)

                do m_index = -2,2
                        ! minus old
                        phi_l2 = wlc_PHI_l2(m_index,J)
                        wlc_dx_maierSaupe =  wlc_dx_maierSaupe - VV*phi_l2**2
                        ! plus new
                        phi_l2 = phi_l2 + wlc_DPHI_l2(m_index,I)
                        wlc_dx_maierSaupe =  wlc_dx_maierSaupe + VV*phi_l2**2
                enddo
            enddo
        else
            do I = 1,wlc_NPHI
                if (WLC_P__FRACTIONAL_BIN) then
                    VV = wlc_Vol(I)
                    if (VV.le.0.1_dp) CYCLE
                endif
                J = wlc_inDPHI(I)
                ! minus old
                phi_A = wlc_PHIA(J)
                phi_B = wlc_PHIB(J)
                wlc_Dx_Kap = wlc_Dx_Kap - VV*((phi_A + phi_B-1.0_dp)**2)

                ! plus new
                phi_A = phi_A + wlc_DPHIA(I)
                phi_B = phi_B + wlc_DPHIB(I)
                wlc_Dx_Kap = wlc_Dx_Kap + VV*((phi_A + phi_B-1.0_dp)**2)
            enddo
        endif

    !------------------------------------------------------------
    !
    !    A-B melt hamiltonian
    !
    !--------------------------------------------------------------
    ! Here Chi and Kap are in units of KT/beadVolume
    case('ABmelt') ! Melt Hamiltonian
        do I = 1,wlc_NPHI
            J = wlc_inDPHI(I)
            if (WLC_P__FRACTIONAL_BIN) then
                VV = wlc_Vol(I)
                if (VV.le.0.1_dp) CYCLE
            endif
            ! new
            phi_A = wlc_PHIA(J) + wlc_DPHIA(I)
            phi_B = wlc_PHIB(J) + wlc_DPHIB(I)
            phi_h = wlc_PHIH(J)
            wlc_Dx_Chi = wlc_Dx_Chi + (VV/WLC_P__BEADVOLUME)*phi_A*phi_B
            wlc_Dx_Kap = wlc_Dx_Kap + (VV/WLC_P__BEADVOLUME)*((phi_A + phi_B-1.0_dp)**2)
            wlc_Dx_Field = wlc_Dx_Field-phi_h*phi_A
            ! minus old
            wlc_Dx_Chi = wlc_Dx_Chi-(VV/WLC_P__BEADVOLUME)*(wlc_PHIA(J)*wlc_PHIB(J))
            wlc_Dx_Kap = wlc_Dx_Kap-(VV/WLC_P__BEADVOLUME)*((wlc_PHIA(J) + wlc_PHIB(J)-1.0_dp)**2)
            wlc_Dx_Field = wlc_Dx_Field + phi_h*wlc_PHIA(J)
        enddo
    !------------------------------------------------------------
    !
    !    A-B solution hamiltonian
    !
    !--------------------------------------------------------------
    ! Here Chi and Kap are in units of KT/beadVolume
    case('ABsolution') ! Melt Hamiltonian
        do I = 1,wlc_NPHI
            J = wlc_inDPHI(I)
            if (WLC_P__FRACTIONAL_BIN) then
                VV = wlc_Vol(I)
                if (VV.le.0.1_dp) CYCLE
            endif
            ! new
            phi_A = wlc_PHIA(J) + wlc_DPHIA(I)
            phi_B = wlc_PHIB(J) + wlc_DPHIB(I)
            phi_h = wlc_PHIH(J)
            wlc_Dx_Chi = wlc_Dx_Chi + (VV/WLC_P__BEADVOLUME)*phi_A*phi_B
            PHIPoly = phi_A + phi_B
            if(PHIPoly > 1.0_dp) then
                wlc_Dx_Kap = wlc_Dx_Kap + (VV/WLC_P__BEADVOLUME)*(PHIPoly-1.0_dp)**2
            endif
            ! minus old
            phi_A = wlc_PHIA(J)
            phi_B = wlc_PHIB(J)
            PHIPoly = phi_A + phi_B
            if(PHIPoly > 1.0_dp) then
                wlc_Dx_Kap = wlc_Dx_Kap - (VV/WLC_P__BEADVOLUME)*(PHIPoly-1.0_dp)**2
            endif
            wlc_Dx_Chi = wlc_Dx_Chi-(VV/WLC_P__BEADVOLUME)*(phi_A+phi_B)
        enddo
    !---------------------------------------------------------
    !
    !    Chromatin Hamiltonian
    !
    ! ---------------------------------------------------------
    ! Here Chi and Kap are in units of KT/beadVolume
    case('chromatin')
        do I = 1,wlc_NPHI
            J = wlc_inDPHI(I)
            if (WLC_P__FRACTIONAL_BIN) then
                VV = wlc_Vol(I)
                if (VV.le.0.1_dp) CYCLE
            endif
            ! new ...
            PHIPoly = wlc_PHIA(J) + wlc_DPHIA(I) + wlc_PHIB(J) + wlc_DPHIB(I)
            wlc_Dx_Chi = wlc_Dx_Chi + (VV/WLC_P__BEADVOLUME)*PHIPoly*(1.0_dp-PHIPoly)
            wlc_Dx_Couple = wlc_Dx_Couple + VV*(wlc_PHIA(J) + wlc_DPHIA(I))**2
            if(PHIPoly > 0.5_dp) then
               wlc_Dx_Kap = wlc_Dx_Kap + (VV/WLC_P__BEADVOLUME)*(PHIPoly-0.5_dp)**2
            endif
            ! minus old
            PHIPoly = wlc_PHIA(J) + wlc_PHIB(J)
            wlc_Dx_Chi = wlc_Dx_Chi-(VV/WLC_P__BEADVOLUME)*PHIPoly*(1.0_dp-PHIPoly)
            wlc_Dx_Couple = wlc_Dx_Couple-VV*(wlc_PHIA(J))**2
            if(PHIPoly > 0.5_dp) then
               wlc_Dx_Kap = wlc_Dx_Kap-(VV/WLC_P__BEADVOLUME)*(PHIPoly-0.5_dp)**2
            endif
        enddo
    case('chromatin2')
        do I = 1,wlc_NPHI
            J = wlc_inDPHI(I)
            if (WLC_P__FRACTIONAL_BIN) then
                VV = wlc_Vol(I)
                if (VV.le.0.1_dp) CYCLE
            endif
            ! new ...
            PHI_A = wlc_PHIA(J) + wlc_DPHIA(I)
            PHI_B = wlc_PHIB(J) + wlc_DPHIB(I)
            wlc_Dx_Chi = wlc_Dx_Chi + VV*PHI_B*(1.0_dp-PHI_B)
            wlc_Dx_Couple = wlc_Dx_Couple + VV*(PHI_A**2)
            if(PHI_B > 0.5_dp) then
               wlc_Dx_Kap = wlc_Dx_Kap + VV*(PHI_B-0.5_dp)**2
            endif
            ! minus old
            PHI_A = wlc_PHIA(J)
            PHI_B = wlc_PHIB(J)
            wlc_Dx_Chi = wlc_Dx_Chi - VV*PHI_B*(1.0_dp-PHI_B)
            wlc_Dx_Couple = wlc_Dx_Couple - VV*(PHI_A**2)
            if(PHI_B > 0.5_dp) then
               wlc_Dx_Kap = wlc_Dx_Kap - VV*(PHI_B-0.5_dp)**2
            endif
        enddo
    end select
endif
wlc_dx_chi = wlc_dx_chi*wlc_p%CHI_ON
wlc_dx_couple = wlc_dx_couple*wlc_p%COUPLE_ON
wlc_dx_Kap = wlc_dx_Kap*wlc_p%KAP_ON

wlc_DEChi = wlc_p%CHI*        wlc_dx_chi
wlc_DECouple = wlc_p%HP1_BIND*wlc_dx_couple
wlc_DEKap = wlc_p%KAP*        wlc_dx_Kap
wlc_DEField = wlc_p%HA*       wlc_dx_Field
wlc_deMaierSaupe = wlc_p%CHI_L2*wlc_dx_maierSaupe
RETURN
END subroutine

!---------------------------------------------------------------!
