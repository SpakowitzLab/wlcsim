#include "../defines.inc"
!
! This file is where quinn hardcodes certain parameters of the MC that he wants
! to change depending on what step he's on.
!---------------------------------------------------------------*
subroutine strength_schedule(wlc_p,wlc_d)
    use params
    implicit none
    type(wlcsim_params), intent(inout) :: wlc_p
    type(wlcsim_data), intent(inout) :: wlc_d

    if (wlc_d%mc_ind <= WLC_P__NNOINT) then
        wlc_p%FIELD_INT_ON = .true.
    else
        wlc_p%FIELD_INT_ON = .true.
    endif
    if(wlc_d%mc_ind.lt.WLC_P__N_KAP_ON) then
        wlc_p%KAP_ON = 0.0_dp
    else
        wlc_p%KAP_ON = 1.0_dp
    endif

    if(wlc_d%mc_ind.lt.WLC_P__N_CHI_ON) then
!        PTON = .False.
        wlc_p%CHI_ON = 0.0_dp
    else
!        PTON = .True.
        wlc_p%CHI_ON = 1.0_dp
    endif

    return
end subroutine
!---------------------------------------------------------------*
