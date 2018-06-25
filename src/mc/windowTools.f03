#include "../defines.inc"
module windowToos
implicit none
contains

function exponential_random_int(window,rand_stat) result(output)
    ! this function gives a random exponentially distributed intiger
    ! the most likely outcome is 0
    use params, only: dp
    use mersenne_twister
    implicit none
    type(random_stat), intent(inout) :: rand_stat  ! status of random number generator
    real urnd(1) ! single random number
    real(dp), intent(in) :: window
    integer output
    call random_number(urnd,rand_stat)
    output  = nint(-1.0_dp*log(urnd(1)+0.000001_dp)*window+0.0001_dp)
    output = abs(output)
end function exponential_random_int

subroutine enforceBinding(rand_stat,IB1,IB2,IT1,IT2,wlc_d,max_window,success)
use params, only: dp, wlcsim_data
use mersenne_twister
implicit none
type(random_stat), intent(inout) :: rand_stat  ! status of random number generator
type(wlcsim_data), intent(in) :: wlc_d
integer, intent(inout) :: IT1,IT2,IB1,IB2
real(dp), intent(in) :: max_window
logical, intent(inout) :: success
integer otherEnd
real urnd(1) ! single random number
integer IT1_temp, IT2_temp, IB1_temp, IB2_temp, I
call random_number(urnd,rand_stat)
if (WLC_P__PROB_BIND_RESPECTING_MOVE > urnd(1)) then
    IB1_temp=IB1
    IT1_temp=IT1
    IB2_temp=IB2
    IT2_temp=IT2
    do I =IT1,IT2
        otherEnd=wlc_d%ExplicitBindingPair(I)
        if (WLC_P__NP>1) then
            ! make sure the other end is on the same polymer
            if ((IT1-1)/WLC_P__NB .ne. (otherEnd-1)/WLC_P__NB) cycle
        endif
        if (otherEnd < 1) cycle
        if (otherEnd < IT1) then  ! Loop to point before IT1
            IB1_temp=IB1-IT1+otherEnd
            IT1_temp=otherEnd
        elseif (otherEnd > IT2) then ! Loop to point after IT2
            IB2_temp=IB2-IT2+otherEnd
            IT2_temp=otherEnd
            exit
        endif
    enddo
    if (IB2_temp-IB1_temp<max_window) then
        ! prevent extremely long crank shaft moves
        IB1=IB1_temp
        IT1=IT1_temp
        IB2=IB2_temp
        IT2=IT2_temp
    else
        if (WLC_P__PROB_BIND_RESPECTING_MOVE > 0.9999_dp)  success=.FALSE.
    endif
endif
end subroutine

end module windowToos
