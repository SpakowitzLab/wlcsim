!--------------------------------------------------------------------
!
!
! This subroutine calculates the field Hamiltonian from the phi values.
! It puts the output in: md%dx_chi, md%dx_chi, md%dx_couple, md%dx_couple, md%dx_Kap, md%dx_Kap, md%DEChi,md%dx_chi, md%DECouple,
! md%dx_couple, md%DEKap, md%dx_Kap, md%DEField, md%dx_Field
!      by Quinn MacPherson based on code from Shifan Mao
!       Made a separate function on 7/8/16
!
!   If initialize then calculate all bins.
!   Otherwise calcualte only specified bins.
!-------------------------------------------------------------------

subroutine hamiltonian(mc,md,initialize)
use params
implicit none
TYPE(wlcsim_params), intent(inout) :: mc   ! <---- Contains output
TYPE(wlcsim_data), intent(inout) :: md
logical, intent(in) :: initialize ! Need to do all beads
double precision PHIPoly ! fraction polymer
double precision phi_A ! demsotu of A
double precision phi_B ! density of B
double precision phi_h ! strength of field
double precision VV ! volume of bin
integer I,J ! for looping


md%dx_Chi=0.0_dp
md%Dx_Couple=0.0_dp
md%Dx_Kap=0.0_dp
md%Dx_Field=0.0_dp
if (initialize) then  ! calculate absolute energy
    if (mc%solType.eq.0) then ! Melt Hamiltonian
        do I=1,mc%NBIN
            VV=md%Vol(I)
            if (VV.le.0.1_dp) CYCLE
            md%Dx_Chi=md%Dx_Chi+(VV/mc%beadVolume)*(md%PHIA(I)*md%PHIB(I))
            md%Dx_Kap=md%dx_Kap+(VV/mc%beadVolume)*((md%PHIA(I)+md%PHIB(I)-1.0_dp)**2)
            md%Dx_Field=md%dx_Field-md%PHIH(I)*md%PHIA(I)
        enddo
    elseif(mc%solType.eq.1) then ! Chromatin Hamiltonian
        do I=1,mc%NBIN
            VV=md%Vol(I)
            if (VV.le.0.1_dp) CYCLE
            PHIPoly=md%PHIA(I)+md%PHIB(I)
            md%Dx_Chi=md%Dx_Chi+(VV/mc%beadVolume)*PHIPoly*(1.0_dp-PHIPoly)
            md%Dx_Couple=md%Dx_Couple+VV*(md%PHIA(I))**2
            if(PHIPoly.GT.1.0_dp) then
                md%Dx_Kap=md%Dx_Kap+(VV/mc%beadVolume)*(PHIPoly-1.0_dp)**2
            endif
        enddo
    else
        print*, "Error in MC_int, solType",mc%solType, &
                " notdefined"
    endif
else ! Calculate change in energy
    if (mc%solType.eq.0) then ! Melt Hamiltonian
        do I=1,md%NPHI
            J=md%INDPHI(I)
            VV=md%Vol(J)
            if (VV.le.0.1_dp) CYCLE
            ! new
            phi_A=md%PHIA(J)+md%DPHIA(I)
            phi_B=md%PHIB(J)+md%DPHIB(I)
            phi_h=md%PHIH(J)
            md%Dx_Chi=md%Dx_Chi+(VV/mc%beadVolume)*phi_A*phi_B
            md%Dx_Kap=md%Dx_Kap+(VV/mc%beadVolume)*((phi_A+phi_B-1.0_dp)**2)
            md%Dx_Field=md%Dx_Field-phi_h*phi_A
            ! minus old
            md%Dx_Chi=md%Dx_Chi-(VV/mc%beadVolume)*(md%PHIA(J)*md%PHIB(J))
            md%Dx_Kap=md%Dx_Kap-(VV/mc%beadVolume)*((md%PHIA(J)+md%PHIB(J)-1.0_dp)**2)
            md%Dx_Field=md%Dx_Field+phi_h*md%PHIA(J)
        enddo
    elseif(mc%solType.eq.1) then ! Chromatin Hamiltonian
        do I=1,md%NPHI
            J=md%INDPHI(I)
            VV=md%Vol(J)
            if (VV.le.0.1_dp) CYCLE
            ! new ...
            PHIPoly=md%PHIA(J)+md%DPHIA(I)+md%PHIB(J)+md%DPHIB(I)
            md%Dx_Chi=md%Dx_Chi+(VV/mc%beadVolume)*PHIPoly*(1.0_dp-PHIPoly)
            md%Dx_Couple=md%Dx_Couple+VV*(md%PHIA(J)+md%DPHIA(I))**2
            if(PHIPoly.GT.1.0_dp) then
               md%Dx_Kap=md%Dx_Kap+(VV/mc%beadVolume)*(PHIPoly-1.0_dp)**2
            endif
            ! minus old
            PHIPoly=md%PHIA(J)+md%PHIB(J)
            md%Dx_Chi=md%Dx_Chi-(VV/mc%beadVolume)*PHIPoly*(1.0_dp-PHIPoly)
            md%Dx_Couple=md%Dx_Couple-VV*(md%PHIA(J))**2
            if(PHIPoly.GT.1.0_dp) then
               md%Dx_Kap=md%Dx_Kap-(VV/mc%beadVolume)*(PHIPoly-1.0_dp)**2
            endif
        enddo
    endif
endif
md%dx_chi=md%dx_chi*mc%CHI_ON
md%dx_couple=md%dx_couple*mc%Couple_ON
md%dx_Kap=md%dx_Kap*mc%KAP_ON

md%DEChi=mc%Chi*        md%dx_chi
md%DECouple=mc%HP1_Bind*md%dx_couple
md%DEKap=mc%Kap*        md%dx_Kap
md%DEField=mc%hA*       md%dx_Field
RETURN
END subroutine

!---------------------------------------------------------------!
