!
! This file is where quinn hardcodes certain parameters of the MC that he wants
! to change depending on what step he's on.
!---------------------------------------------------------------*
subroutine strength_schedule(mc,md)
    use params
    implicit none
    type(wlcsim_params), intent(inout) :: mc
    type(wlcsim_data), intent(inout) :: md

    if (md%mc_ind.LE.mc%NNOINT) then
        mc%field_int_on=.true.
    else
        mc%field_int_on=.true.
    endif
    if(md%mc_ind.lt.mc%N_KAP_ON) then
        mc%KAP_ON=0.0_dp
    else
        mc%KAP_ON=1.0_dp
    endif

    if(md%mc_ind.lt.mc%N_CHI_ON) then
!        PTON=.False.
        mc%CHI_ON=0.0_dp
    else
!        PTON=.True.
        mc%CHI_ON=1.0_dp
    endif

    return
end subroutine
!---------------------------------------------------------------*
