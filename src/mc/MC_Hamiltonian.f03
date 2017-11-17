#include "../defines.inc"
!--------------------------------------------------------------------
!
!
! This subroutine calculates the field Hamiltonian from the phi values.
! It puts the output in: wlc_d%dx_chi, wlc_d%dx_chi, wlc_d%dx_couple, wlc_d%dx_couple, wlc_d%dx_Kap, wlc_d%dx_Kap, wlc_d%DEChi,wlc_d%dx_chi, wlc_d%DECouple,
! wlc_d%dx_couple, wlc_d%DEKap, wlc_d%dx_Kap, wlc_d%DEField, wlc_d%dx_Field
!      by Quinn MacPherson based on code from Shifan Mao
!       Made a separate function on 7/8/16
!
!   If initialize then calculate all bins.
!   Otherwise calcualte only specified bins.
!-------------------------------------------------------------------

subroutine hamiltonian(wlc_p,wlc_d,initialize)
use params,only: dp,wlcsim_params,wlcsim_data
implicit none
TYPE(wlcsim_params), intent(inout) :: wlc_p   ! <---- Contains output
TYPE(wlcsim_data), intent(inout) :: wlc_d
logical, intent(in) :: initialize ! Need to do all beads
real(dp) PHIPoly ! fraction polymer
real(dp) phi_A ! demsotu of A
real(dp) phi_B ! density of B
real(dp) phi_h ! strength of field
real(dp) phi_l2 ! strength of field
real(dp) VV ! volume of bin
integer I,J,m_index ! for looping

wlc_d%dx_Chi = 0.0_dp
wlc_d%Dx_Couple = 0.0_dp
wlc_d%Dx_Kap = 0.0_dp
wlc_d%Dx_Field = 0.0_dp
wlc_d%dx_maierSaupe = 0.0_dp
if (initialize) then  ! calculate absolute energy

    select case(WLC_P__FIELDINTERACTIONTYPE) ! pick which keyword, case matchign string must be all uppercase
    !-------------------------------------------------------
    !
    !   Maier Saupe Melt Interaction
    !
    !-------------------------------------------------------
    ! In this problem Kap and chi are in units of kT/(simulation units cubed)
    ! If VV=1.0 than this is just kT/(bin volume)
    case('MaierSaupe') 
        VV = WLC_P__DBIN**3

        if (wlc_p%CHI_L2_ON) then
            do I = 1,wlc_p%NBIN
                do m_index = -2,2
                    wlc_d%dx_maierSaupe =  wlc_d%dx_maierSaupe + VV*wlc_d%PHI_l2(m_index,I)**2
                enddo
            enddo
        endif
        do I = 1,wlc_p%NBIN
            VV = wlc_d%Vol(I)
            if (VV.le.0.1_dp) CYCLE
            wlc_d%Dx_Kap = wlc_d%dx_Kap + VV*((wlc_d%PHIA(I) + wlc_d%PHIB(I)-1.0_dp)**2)
        enddo
    !------------------------------------------------------------
    !
    !    A-B melt hamiltonian
    !
    !--------------------------------------------------------------
    ! Here Chi and Kap are in units of KT/beadVolume
    case('ABmelt') ! Melt Hamiltonian
        do I = 1,wlc_p%NBIN
            VV = wlc_d%Vol(I)
            if (VV.le.0.1_dp) CYCLE
            wlc_d%Dx_Chi = wlc_d%Dx_Chi + (VV/WLC_P__BEADVOLUME)*(wlc_d%PHIA(I)*wlc_d%PHIB(I))
            wlc_d%Dx_Kap = wlc_d%dx_Kap + (VV/WLC_P__BEADVOLUME)*((wlc_d%PHIA(I) + wlc_d%PHIB(I)-1.0_dp)**2)
            wlc_d%Dx_Field = wlc_d%dx_Field-wlc_d%PHIH(I)*wlc_d%PHIA(I)
        enddo
    !------------------------------------------------------------
    !
    !    A-B solution hamiltonian
    !
    !--------------------------------------------------------------
    ! Here Chi and Kap are in units of KT/beadVolume
    case('ABsolution') ! A,B,Solvent Hamiltonian
        do I = 1,wlc_p%NBIN
            VV = wlc_d%Vol(I)
            if (VV.le.0.1_dp) CYCLE
            wlc_d%Dx_Chi = wlc_d%Dx_Chi + (VV/WLC_P__BEADVOLUME)*(wlc_d%PHIA(I)*wlc_d%PHIB(I))
            PHIPoly = wlc_d%PHIA(I) + wlc_d%PHIB(I)
            if(PHIPoly > 1.0_dp) then
                wlc_d%Dx_Kap = wlc_d%Dx_Kap + (VV/WLC_P__BEADVOLUME)*(PHIPoly-1.0_dp)**2
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
            !VV = wlc_d%Vol(I)
            VV = WLC_P__DBIN**3
            if (VV.le.0.1_dp) CYCLE
            PHIPoly = wlc_d%PHIA(I) + wlc_d%PHIB(I)
            wlc_d%Dx_Chi = wlc_d%Dx_Chi + (VV/WLC_P__BEADVOLUME)*PHIPoly*(1.0_dp-PHIPoly)
            wlc_d%Dx_Couple = wlc_d%Dx_Couple + VV*(wlc_d%PHIA(I))**2
            if(PHIPoly > 1.0_dp) then
                wlc_d%Dx_Kap = wlc_d%Dx_Kap + (VV/WLC_P__BEADVOLUME)*(PHIPoly-1.0_dp)**2
            endif
        enddo
    end select


else ! Calculate change in energy

    select case(WLC_P__FIELDINTERACTIONTYPE) ! pick which keyword, case matchign string must be all uppercase
    !-------------------------------------------------------
    !
    !   Maier Saupe Melt Interaction
    !
    !-------------------------------------------------------
    ! In this problem Kap and chi are in units of kT/binVolume
    case('MaierSaupe') 
        VV = WLC_P__DBIN**3
        if (wlc_p%CHI_L2_ON) then
            do I = 1,wlc_d%NPHI
                J = wlc_d%inDPHI(I)
                ! minus old
                phi_A = wlc_d%PHIA(J) 
                phi_B = wlc_d%PHIB(J) 
                wlc_d%Dx_Kap = wlc_d%Dx_Kap - VV*((phi_A + phi_B-1.0_dp)**2)


                ! plus new
                phi_A = phi_A + wlc_d%DPHIA(I)
                phi_B = phi_B + wlc_d%DPHIB(I)
                wlc_d%Dx_Kap = wlc_d%Dx_Kap + VV*((phi_A + phi_B-1.0_dp)**2)

                do m_index = -2,2
                        ! minus old
                        phi_l2 = wlc_d%PHI_l2(m_index,J)
                        wlc_d%dx_maierSaupe =  wlc_d%dx_maierSaupe - VV*phi_l2**2
                        ! plus new
                        phi_l2 = phi_l2 + wlc_d%DPHI_l2(m_index,I)
                        wlc_d%dx_maierSaupe =  wlc_d%dx_maierSaupe + VV*phi_l2**2
                enddo
            enddo
        else
            do I = 1,wlc_d%NPHI
                J = wlc_d%inDPHI(I)
                ! minus old
                phi_A = wlc_d%PHIA(J) 
                phi_B = wlc_d%PHIB(J) 
                wlc_d%Dx_Kap = wlc_d%Dx_Kap - VV*((phi_A + phi_B-1.0_dp)**2)

                ! plus new
                phi_A = phi_A + wlc_d%DPHIA(I)
                phi_B = phi_B + wlc_d%DPHIB(I)
                wlc_d%Dx_Kap = wlc_d%Dx_Kap + VV*((phi_A + phi_B-1.0_dp)**2)
            enddo
        endif

    !------------------------------------------------------------
    !
    !    A-B melt hamiltonian
    !
    !--------------------------------------------------------------
    ! Here Chi and Kap are in units of KT/beadVolume
    case('ABmelt') ! Melt Hamiltonian
        do I = 1,wlc_d%NPHI
            J = wlc_d%inDPHI(I)
            VV = wlc_d%Vol(J)
            if (VV.le.0.1_dp) CYCLE
            ! new
            phi_A = wlc_d%PHIA(J) + wlc_d%DPHIA(I)
            phi_B = wlc_d%PHIB(J) + wlc_d%DPHIB(I)
            phi_h = wlc_d%PHIH(J)
            wlc_d%Dx_Chi = wlc_d%Dx_Chi + (VV/WLC_P__BEADVOLUME)*phi_A*phi_B
            wlc_d%Dx_Kap = wlc_d%Dx_Kap + (VV/WLC_P__BEADVOLUME)*((phi_A + phi_B-1.0_dp)**2)
            wlc_d%Dx_Field = wlc_d%Dx_Field-phi_h*phi_A
            ! minus old
            wlc_d%Dx_Chi = wlc_d%Dx_Chi-(VV/WLC_P__BEADVOLUME)*(wlc_d%PHIA(J)*wlc_d%PHIB(J))
            wlc_d%Dx_Kap = wlc_d%Dx_Kap-(VV/WLC_P__BEADVOLUME)*((wlc_d%PHIA(J) + wlc_d%PHIB(J)-1.0_dp)**2)
            wlc_d%Dx_Field = wlc_d%Dx_Field + phi_h*wlc_d%PHIA(J)
        enddo
    !------------------------------------------------------------
    !
    !    A-B solution hamiltonian
    !
    !--------------------------------------------------------------
    ! Here Chi and Kap are in units of KT/beadVolume
    case('ABsolution') ! Melt Hamiltonian
        do I = 1,wlc_d%NPHI
            J = wlc_d%inDPHI(I)
            VV = wlc_d%Vol(J)
            if (VV.le.0.1_dp) CYCLE
            ! new
            phi_A = wlc_d%PHIA(J) + wlc_d%DPHIA(I)
            phi_B = wlc_d%PHIB(J) + wlc_d%DPHIB(I)
            phi_h = wlc_d%PHIH(J)
            wlc_d%Dx_Chi = wlc_d%Dx_Chi + (VV/WLC_P__BEADVOLUME)*phi_A*phi_B
            PHIPoly = phi_A + phi_B
            if(PHIPoly > 1.0_dp) then
                wlc_d%Dx_Kap = wlc_d%Dx_Kap + (VV/WLC_P__BEADVOLUME)*(PHIPoly-1.0_dp)**2
            endif
            ! minus old
            phi_A = wlc_d%PHIA(J)
            phi_B = wlc_d%PHIB(J)
            PHIPoly = phi_A + phi_B
            if(PHIPoly > 1.0_dp) then
                wlc_d%Dx_Kap = wlc_d%Dx_Kap - (VV/WLC_P__BEADVOLUME)*(PHIPoly-1.0_dp)**2
            endif
            wlc_d%Dx_Chi = wlc_d%Dx_Chi-(VV/WLC_P__BEADVOLUME)*(phi_A+phi_B)
        enddo
    !---------------------------------------------------------
    !
    !    Chromatin Hamiltonian
    !
    ! ---------------------------------------------------------
    ! Here Chi and Kap are in units of KT/beadVolume
    case('chromatin')
        do I = 1,wlc_d%NPHI
            J = wlc_d%inDPHI(I)
            !VV = wlc_d%Vol(J)
            VV = WLC_P__DBIN**3
            if (VV.le.0.1_dp) CYCLE
            ! new ...
            PHIPoly = wlc_d%PHIA(J) + wlc_d%DPHIA(I) + wlc_d%PHIB(J) + wlc_d%DPHIB(I)
            wlc_d%Dx_Chi = wlc_d%Dx_Chi + (VV/WLC_P__BEADVOLUME)*PHIPoly*(1.0_dp-PHIPoly)
            wlc_d%Dx_Couple = wlc_d%Dx_Couple + VV*(wlc_d%PHIA(J) + wlc_d%DPHIA(I))**2
            if(PHIPoly > 1.0_dp) then
               wlc_d%Dx_Kap = wlc_d%Dx_Kap + (VV/WLC_P__BEADVOLUME)*(PHIPoly-1.0_dp)**2
            endif
            ! minus old
            PHIPoly = wlc_d%PHIA(J) + wlc_d%PHIB(J)
            wlc_d%Dx_Chi = wlc_d%Dx_Chi-(VV/WLC_P__BEADVOLUME)*PHIPoly*(1.0_dp-PHIPoly)
            wlc_d%Dx_Couple = wlc_d%Dx_Couple-VV*(wlc_d%PHIA(J))**2
            if(PHIPoly > 1.0_dp) then
               wlc_d%Dx_Kap = wlc_d%Dx_Kap-(VV/WLC_P__BEADVOLUME)*(PHIPoly-1.0_dp)**2
            endif
        enddo
    end select
endif
wlc_d%dx_chi = wlc_d%dx_chi*wlc_p%CHI_ON
wlc_d%dx_couple = wlc_d%dx_couple*wlc_p%COUPLE_ON
wlc_d%dx_Kap = wlc_d%dx_Kap*wlc_p%KAP_ON

wlc_d%DEChi = wlc_p%CHI*        wlc_d%dx_chi
wlc_d%DECouple = wlc_p%HP1_BIND*wlc_d%dx_couple
wlc_d%DEKap = wlc_p%KAP*        wlc_d%dx_Kap
wlc_d%DEField = wlc_p%HA*       wlc_d%dx_Field
wlc_d%deMaierSaupe = wlc_p%CHI_L2*wlc_d%dx_maierSaupe
RETURN
END subroutine

!---------------------------------------------------------------!
