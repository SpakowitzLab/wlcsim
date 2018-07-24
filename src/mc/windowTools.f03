#include "../defines.inc"
module windowTools
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
logical, intent(out) :: success
integer otherEnd
real urnd(1) ! single random number
integer IT1_temp, IT2_temp, IB1_temp, IB2_temp, I
success = .True.
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
        return
    endif
endif
end subroutine

subroutine drawWindow(wlc_d,window,maxWindow,enforceBind,rand_stat,IT1,IT2,IB1,IB2,IP,DIB,success)
use params, only: dp, wlcsim_data
use mersenne_twister
implicit none
type(wlcsim_data), intent(in) :: wlc_d
real(dp), intent(in) :: WindoW ! Size of window for bead selection
real(dp), intent(in) :: maxWindow
logical, intent(in) :: enforceBind  ! don't do move with internal bind
type(random_stat), intent(inout) :: rand_stat  ! status of random number generator
integer, intent(out) :: IP    ! Test polymer
integer, intent(out) :: IB1   ! Test bead position 1
integer, intent(out) :: IT1   ! Index of test bead 1
integer, intent(out) :: IB2   ! Test bead position 2
integer, intent(out) :: IT2   ! Index of test bead 2
integer, intent(out) :: dib   ! number of beads moved by move (plus or minus a few)
logical, intent(out) :: success
integer irnd(1)
real urnd(1) ! single random number
integer TEMP

success = .TRUE.  ! True unless set to false
call random_index(WLC_P__NP,irnd,rand_stat)
IP=irnd(1)
call random_index(WLC_P__NB,irnd,rand_stat)
IB1=irnd(1)
if (WLC_P__WINTYPE.eq.0) then
    IB2 = IB1 +exponential_random_int(window,rand_stat)
elseif (WLC_P__WINTYPE.eq.1.and..not.WLC_P__RING) then
    call random_number(urnd,rand_stat)
    IB2 = IB1 + (2*nint(urnd(1))-1)* &
           exponential_random_int(window,rand_stat)
elseif (WLC_P__WINTYPE.eq.1.and.WLC_P__RING) then
    IB2 = IB1 + exponential_random_int(window,rand_stat)
else
    call stop_if_err(1, "Warning: WLC_P__WINTYPE not recognized")
endif

DIB = IB2-IB1

if (WLC_P__RING) then
    if (IB2 > WLC_P__NB) then
        IB2 = DIB-(WLC_P__NB-IB1)
    endif
    if (WLC_P__EXPLICIT_BINDING) then
        print*, "Ring polymer not set up to use explicit binding"
        print*, "Need to write special loop skiping code"
        stop
    endif
else
    if (IB2 > WLC_P__NB) then
        IB2 = WLC_P__NB
    endif
    if (IB2 < 1) then
       IB2 = 1
    endif
    if (IB2 < IB1) then
        TEMP = IB1
        IB1 = IB2
        IB2 = TEMP
    endif
    IT2 = WLC_P__NB*(IP-1) + IB2
    IT1 = WLC_P__NB*(IP-1) + IB1
    if (WLC_P__EXPLICIT_BINDING .and. enforceBind) then
        call enforceBinding(rand_stat,IB1,IB2,IT1,IT2,wlc_d,maxWindow,success)
        if (success .eqv. .False.) return
    endif
endif

IT1 = WLC_P__NB*(IP-1) + IB1
IT2 = WLC_P__NB*(IP-1) + IB2

DIB = IB2-IB1
end subroutine drawWindow


end module windowTools
