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

subroutine MC_explicit_binding()
! values from wlcsim_data
use params, only: wlc_DEExplicitBinding, wlc_R, wlc_RP, &
    wlc_ExplicitBindingPair,  wlc_nPointsMoved, wlc_pointsMoved
use params, only:  dp
implicit none

!   Internal variables

integer ii,jj
integer otherEnd
real(dp) r(3)

wlc_DEExplicitBinding = 0.0_dp


do jj = 1,wlc_nPointsMoved
    ii = wlc_pointsMoved(jj)
    otherEnd = wlc_ExplicitBindingPair(ii) 
    if(otherEnd .le. 0) cycle

    ! plus new
    if (isnan(wlc_RP(1,otherEnd))) then  ! assumes that RP is up to date or nan
        r = wlc_RP(:,ii) - wlc_R(:,otherEnd)
    else
        if (otherEnd .le. ii) cycle ! don't doule count
        r = wlc_RP(:,ii) - wlc_RP(:,otherEnd)
    endif
    wlc_DEExplicitBinding = wlc_DEExplicitBinding + &
        WLC_P__EXPLICIT_BIND_ENERGY*(r(1)**2 + r(2)**2 + r(3)**2 )

    ! minus old
    r = wlc_R(:,ii) - wlc_R(:,otherEnd)
    wlc_DEExplicitBinding = wlc_DEExplicitBinding - &
        WLC_P__EXPLICIT_BIND_ENERGY*(r(1)**2 + r(2)**2 + r(3)**2 )

enddo


RETURN
END

!---------------------------------------------------------------!
subroutine MC_explicit_binding_from_scratch()
! values from wlcsim_data
use params, only: wlc_DEExplicitBinding, wlc_R, wlc_ExplicitBindingPair
use params, only:  dp
implicit none

!   iputs

!   Internal variables

integer ii
integer otherEnd
real(dp) r(3)

wlc_DEExplicitBinding = 0.0_dp

do ii = 1,WLC_P__NT
    otherEnd = wlc_ExplicitBindingPair(ii)
    if (otherEnd .le. 0) cycle
    if (wlc_ExplicitBindingPair(otherEnd).ne.ii) then
        print*, "Reflexive not found",ii,otherEnd
    endif
    if(otherEnd .le. ii) cycle ! Exclude <0 and don't double count

    ! plus new
    r = wlc_R(:,ii) - wlc_R(:,otherEnd)
    wlc_DEExplicitBinding = wlc_DEExplicitBinding + &
        WLC_P__EXPLICIT_BIND_ENERGY*(r(1)**2 + r(2)**2 + r(3)**2 )

enddo


RETURN
END

!---------------------------------------------------------------!
