!
! This file is where quinn hardcodes certain parameters of the MC that he wants
! to change depending on what step he's on.
!---------------------------------------------------------------*
subroutine strength_schedule(mc,md)
    use params
    implicit none
    type(wlcsim_params), intent(inout) :: mc
    type(wlcsim_data), intent(inout) :: md

    if (wlc_d%mc_ind <= wlc_p%NNOinT) then
        wlc_p%field_int_on = .true.
    else
        wlc_p%field_int_on = .true.
    endif
    if(wlc_d%mc_ind.lt.wlc_p%N_KAP_ON) then
        wlc_p%KAP_ON = 0.0_dp
    else
        wlc_p%KAP_ON = 1.0_dp
    endif

    if(wlc_d%mc_ind.lt.wlc_p%N_CHI_ON) then
!        PTON = .False.
        wlc_p%CHI_ON = 0.0_dp
    else
!        PTON = .True.
        wlc_p%CHI_ON = 1.0_dp
    endif

    return
end subroutine
!---------------------------------------------------------------*
