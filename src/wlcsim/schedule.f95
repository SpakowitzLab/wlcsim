!
! This file is where quinn hardcodes certain parameters of the MC that he wants
! to change depending on what step he's on.
!---------------------------------------------------------------*
subroutine strength_schedule(mc,inton)
    use setPrecision
    use simMod
    implicit none
    type(MCvar), intent(inout) :: mc
    integer, intent(out) :: inton

    if (mc%ind.LE.mc%NNOINT) then
        INTON=0
    else
        INTON=1
    endif
    if(mc%ind.lt.mc%N_KAP_ON) then
        mc%KAP_ON=0.0_dp
    else
        mc%KAP_ON=1.0_dp
    endif

    if(mc%ind.lt.mc%N_CHI_ON) then
!        PTON=.False.
        mc%CHI_ON=0.0_dp
    else
!        PTON=.True.
        mc%CHI_ON=1.0_dp
    endif


! -----------------------------------------------
!
!      Shut off moves for testing purposes
!
! -----------------------------------------------

    if (.true.) then
        return
    endif

    if (mc%ind.lt.30) then
        mc%moveon(1) = 1;
        mc%moveon(2) = 1;
        mc%moveon(3) = 1;
        mc%moveon(4) = 1;
        mc%moveon(5) = 1;
        mc%moveon(6) = 1;
        mc%moveon(9) = 1;
        mc%moveon(10) = 1;
    elseif (mc%ind.lt.40) then
        mc%moveon(1) = 1;
        mc%moveon(2) = 1;
        mc%moveon(3) = 1;
        mc%moveon(4) = 1;
        mc%moveon(5) = 1;
        mc%moveon(6) = 1;
        mc%moveon(9) = 0;
        mc%moveon(10) = 0;
    elseif (mc%ind.lt.50) then
        mc%moveon(1) = 1;
        mc%moveon(2) = 1;
        mc%moveon(3) = 1;
        mc%moveon(4) = 1;
        mc%moveon(5) = 1;
        mc%moveon(6) = 1;
        mc%moveon(9) = 1;
        mc%moveon(10) = 0;
    elseif (mc%ind.lt.60) then
        mc%moveon(1) = 1;
        mc%moveon(2) = 1;
        mc%moveon(3) = 1;
        mc%moveon(4) = 1;
        mc%moveon(5) = 1;
        mc%moveon(6) = 1;
        mc%moveon(9) = 0;
        mc%moveon(10) = 1;
    elseif (mc%ind.lt.70) then
        mc%moveon(1) = 1;
        mc%moveon(2) = 1;
        mc%moveon(3) = 1;
        mc%moveon(4) = 1;
        mc%moveon(5) = 1;
        mc%moveon(6) = 0;
        mc%moveon(9) = 1;
        mc%moveon(10) = 1;
    elseif (mc%ind.lt.80) then
        mc%moveon(1) = 1;
        mc%moveon(2) = 1;
        mc%moveon(3) = 1;
        mc%moveon(4) = 1;
        mc%moveon(5) = 0;
        mc%moveon(6) = 1;
        mc%moveon(9) = 1;
        mc%moveon(10) = 1;
    elseif (mc%ind.lt.90) then
        mc%moveon(1) = 1;
        mc%moveon(2) = 1;
        mc%moveon(3) = 1;
        mc%moveon(4) = 0;
        mc%moveon(5) = 1;
        mc%moveon(6) = 1;
        mc%moveon(9) = 1;
        mc%moveon(10) = 1;
    elseif (mc%ind.lt.100) then
        mc%moveon(1) = 1;
        mc%moveon(2) = 1;
        mc%moveon(3) = 0;
        mc%moveon(4) = 1;
        mc%moveon(5) = 1;
        mc%moveon(6) = 1;
        mc%moveon(9) = 1;
        mc%moveon(10) = 1;
    elseif (mc%ind.lt.110) then
        mc%moveon(1) = 1;
        mc%moveon(2) = 0;
        mc%moveon(3) = 1;
        mc%moveon(4) = 1;
        mc%moveon(5) = 1;
        mc%moveon(6) = 1;
        mc%moveon(9) = 1;
        mc%moveon(10) = 1;
    elseif (mc%ind.lt.120) then
        mc%moveon(1) = 0;
        mc%moveon(2) = 1;
        mc%moveon(3) = 1;
        mc%moveon(4) = 1;
        mc%moveon(5) = 1;
        mc%moveon(6) = 1;
        mc%moveon(9) = 1;
        mc%moveon(10) = 1;
    else
        mc%moveon(1) = 1;
        mc%moveon(2) = 1;
        mc%moveon(3) = 1;
        mc%moveon(4) = 1;
        mc%moveon(5) = 1;
        mc%moveon(6) = 1;
        mc%moveon(9) = 1;
        mc%moveon(10) = 1;
    endif



    return
end subroutine
!---------------------------------------------------------------*
