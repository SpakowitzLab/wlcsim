#include "../defines.inc"
!---------------------------------------------------------------!

!
!     This subroutine calculates the change in the self energy for
!     a small Monte Carlo move in the position.
!
!     Andrew Spakowitz
!     Written 6-29-04
!
!     Edited by Shifan
!
!     Edited by Quinn in 2016

subroutine MC_explicit_binding(wlc_p,wlc_d,IT1,IT2,IT3,IT4,MCTYPE)
use params
implicit none

!   iputs
TYPE(wlcsim_params), intent(in) :: wlc_p   ! <---- Contains output
TYPE(wlcsim_data), intent(inout) :: wlc_d
integer, intent(in) :: IT1  ! Test bead position 1
integer, intent(in) :: IT2  ! Test bead position 2
integer, intent(in) :: IT3  ! Test bead position 1
integer, intent(in) :: IT4  ! Test bead position 2
integer, intent(in) :: MCTYPE

!   Internal variables

integer ii
integer otherEnd
real(dp) r(3)
real(dp) before, after, delta

wlc_d%DEExplicitBinding = 0.0_dp

if(MCTYPE==4 .or. MCTYPE==7 .or. MCTYPE==12) RETURN
if(MCTYPE==9) then
    print*, "swap move not finished for explicit binding"
    stop 1
endif

do ii = IT1, IT2
    otherEnd = wlc_d%ExplicitBindingPair(ii) 
    if(otherEnd .le. 0) cycle

    ! plus new
    if (otherEnd < IT1 .or. otherEnd > IT2) then
        r = wlc_d%RP(:,ii) - wlc_d%R(:,otherEnd)
    else
        if (otherEnd .le. ii) cycle ! don't doule count
        r = wlc_d%RP(:,ii) - wlc_d%RP(:,otherEnd)
    endif
    wlc_d%DEExplicitBinding = wlc_d%DEExplicitBinding + &
        WLC_P__EXPLICIT_BIND_ENERGY*(r(1)**2 + r(2)**2 + r(3)**2 )

    ! minus old
    r = wlc_d%R(:,ii) - wlc_d%R(:,otherEnd)
    wlc_d%DEExplicitBinding = wlc_d%DEExplicitBinding - &
        WLC_P__EXPLICIT_BIND_ENERGY*(r(1)**2 + r(2)**2 + r(3)**2 )

enddo


RETURN
END

!---------------------------------------------------------------!
subroutine MC_explicit_binding_from_scratch(wlc_p,wlc_d)
use params
implicit none

!   iputs
TYPE(wlcsim_params), intent(in) :: wlc_p   ! <---- Contains output
TYPE(wlcsim_data), intent(inout) :: wlc_d

!   Internal variables

integer ii
integer otherEnd
real(dp) r(3)

wlc_d%DEExplicitBinding = 0.0_dp

do ii = 1,WLC_P__NT
    otherEnd = wlc_d%ExplicitBindingPair(ii)
    if (otherEnd .le. 0) cycle
    if (wlc_d%ExplicitBindingPair(otherEnd).ne.ii) then
        print*, "Reflexive not found",ii,otherEnd
    endif
    if(otherEnd .le. ii) cycle ! Exclude <0 and don't double count

    ! plus new
    r = wlc_d%R(:,ii) - wlc_d%R(:,otherEnd)
    wlc_d%DEExplicitBinding = wlc_d%DEExplicitBinding + &
        WLC_P__EXPLICIT_BIND_ENERGY*(r(1)**2 + r(2)**2 + r(3)**2 )

enddo


RETURN
END

!---------------------------------------------------------------!
