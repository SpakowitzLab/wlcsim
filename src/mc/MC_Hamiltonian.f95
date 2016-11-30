!--------------------------------------------------------------------
!
!
! This subroutine calculates the field Hamiltonian from the phi values.
!
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
TYPE(wlcsim_data), intent(in) :: md
logical, intent(in) :: initialize ! Need to do all beads
double precision PHIPoly ! fraction polymer
double precision phi_A ! demsotu of A
double precision phi_B ! density of B
double precision phi_h ! strength of field
double precision VV ! volume of bin
integer I,J ! for looping


mc%dx_Chi=0.0_dp
mc%Dx_Couple=0.0_dp
mc%Dx_Kap=0.0_dp
mc%Dx_Field=0.0_dp
if (initialize) then  ! calculate absolute energy
    if (mc%solType.eq.0) then ! Melt Hamiltonian
        do I=1,mc%NBIN
            VV=md%Vol(I)
            if (VV.le.0.1_dp) CYCLE
            mc%Dx_Chi=mc%Dx_Chi+(VV/mc%V)*(md%PHIA(I)*md%PHIB(I))
            mc%Dx_Kap=mc%dx_Kap+(VV/mc%V)*((md%PHIA(I)+md%PHIB(I)-1.0_dp)**2)
            mc%Dx_Field=mc%dx_Field-md%PHIH(I)*md%PHIA(I)
        enddo
    elseif(mc%solType.eq.1) then ! Chromatin Hamiltonian
        do I=1,mc%NBIN
            VV=md%Vol(I)
            if (VV.le.0.1_dp) CYCLE
            PHIPoly=md%PHIA(I)+md%PHIB(I)
            mc%Dx_Chi=mc%Dx_Chi+(VV/mc%V)*PHIPoly*(1.0_dp-PHIPoly)
            mc%Dx_Couple=mc%Dx_Couple+VV*(md%PHIA(I))**2
            if(PHIPoly.GT.1.0_dp) then
               mc%Dx_Kap=mc%Dx_Kap+(VV/mc%V)*(PHIPoly-1.0_dp)**2
            endif
        enddo
    else
        print*, "Error in MC_int, solType",mc%solType, &
                " notdefined"
    endif
else ! Calculate change in energy
    if (mc%solType.eq.0) then ! Melt Hamiltonian
        do I=1,mc%NPHI
            J=md%INDPHI(I)
            VV=md%Vol(J)
            if (VV.le.0.1_dp) CYCLE
            ! new
            phi_A=md%PHIA(J)+md%DPHIA(I)
            phi_B=md%PHIB(J)+md%DPHIB(I)
            phi_h=md%PHIH(J)
            mc%Dx_Chi=mc%Dx_Chi+(VV/mc%V)*phi_A*phi_B
            mc%Dx_Kap=mc%Dx_Kap+(VV/mc%V)*((phi_A+phi_B-1.0_dp)**2)
            mc%Dx_Field=mc%Dx_Field-phi_h*phi_A
            ! minus old
            mc%Dx_Chi=mc%Dx_Chi-(VV/mc%V)*(md%PHIA(J)*md%PHIB(J))
            mc%Dx_Kap=mc%Dx_Kap-(VV/mc%V)*((md%PHIA(J)+md%PHIB(J)-1.0_dp)**2)
            mc%Dx_Field=mc%Dx_Field+phi_h*md%PHIA(J)
        enddo
    elseif(mc%solType.eq.1) then ! Chromatin Hamiltonian
        do I=1,mc%NPHI
            J=md%INDPHI(I)
            VV=md%Vol(J)
            if (VV.le.0.1_dp) CYCLE
            ! new ...
            PHIPoly=md%PHIA(J)+md%DPHIA(I)+md%PHIB(J)+md%DPHIB(I)
            mc%Dx_Chi=mc%Dx_Chi+(VV/mc%V)*PHIPoly*(1.0_dp-PHIPoly)
            mc%Dx_Couple=mc%Dx_Couple+VV*(md%PHIA(J)+md%DPHIA(I))**2
            if(PHIPoly.GT.1.0_dp) then
               mc%Dx_Kap=mc%Dx_Kap+(VV/mc%V)*(PHIPoly-1.0_dp)**2
            endif
            ! minus old
            PHIPoly=md%PHIA(J)+md%PHIB(J)
            mc%Dx_Chi=mc%Dx_Chi-(VV/mc%V)*PHIPoly*(1.0_dp-PHIPoly)
            mc%Dx_Couple=mc%Dx_Couple-VV*(md%PHIA(J))**2
            if(PHIPoly.GT.1.0_dp) then
               mc%Dx_Kap=mc%Dx_Kap-(VV/mc%V)*(PHIPoly-1.0_dp)**2
            endif
        enddo
    endif
endif
mc%dx_chi=mc%dx_chi*mc%CHI_ON
mc%dx_couple=mc%dx_couple*mc%Couple_ON
mc%dx_Kap=mc%dx_Kap*mc%KAP_ON

mc%DEChi=mc%Chi*        mc%dx_chi
mc%DECouple=mc%HP1_Bind*mc%dx_couple
mc%DEKap=mc%Kap*        mc%dx_Kap
mc%DEField=mc%h_A*      mc%dx_Field
RETURN
END subroutine

!---------------------------------------------------------------!
