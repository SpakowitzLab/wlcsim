#include "../defines.inc"
module updateRU
! use updateRU, only: updateR, checkR
contains

subroutine updateR(I)
! values from wlcsim_data
use params, only: wlc_bin, wlc_R_period, wlc_R, wlc_UP, wlc_VP&
    , wlc_U, wlc_V, wlc_RP
use params, only:  dp
use binning, only: removeBead, addBead
implicit none
integer, intent(in) :: I

if (WLC_P__NEIGHBOR_BINS) then
    if (WLC_P__CONFINETYPE == 'excludedShpereInPeriodic') then
        call removeBead(wlc_bin,wlc_R_period(:,I),I)
    elseif (WLC_P__CONFINETYPE == 'none') then
        call removeBead(wlc_bin,wlc_R_period(:,I),I)
    elseif (WLC_P__CONFINETYPE == 'sphere') then
        call removeBead(wlc_bin,wlc_R(:,I),I)
    else
        print*, "Not an option yet.  See MCsim."
        stop 1
    endif
endif
wlc_R(:,I) = wlc_RP(:,I)
wlc_U(:,I) = wlc_UP(:,I)
if (WLC_P__LOCAL_TWIST) then
    wlc_V(:,I) = wlc_VP(:,I)/norm2(wlc_VP(:,I))
endif
if (WLC_P__NEIGHBOR_BINS) then
    if (WLC_P__CONFINETYPE == 'excludedShpereInPeriodic' .or. WLC_P__CONFINETYPE=='none') then
        wlc_R_period(1,I)=modulo(wlc_R(1,I),WLC_P__LBOX_X)
        wlc_R_period(2,I)=modulo(wlc_R(2,I),WLC_P__LBOX_Y)
        wlc_R_period(3,I)=modulo(wlc_R(3,I),WLC_P__LBOX_Z)
        call addBead(wlc_bin,wlc_R_period,WLC_P__NT,I)
    elseif (WLC_P__CONFINETYPE == 'sphere') then
        call addBead(wlc_bin,wlc_R,WLC_P__NT,I)
    else
        print*, "Not an option yet.  See MCsim."
        stop 1
    endif
endif

end subroutine updateR


function checkR(I,IT1,IT2)
! values from wlcsim_data
use params, only: wlc_ExplicitBindingPair, wlc_R, wlc_RP
use params, only:  eps
implicit none
integer, intent(in) :: IT1,IT2
integer, intent(in) :: I
logical checkR

integer otherEnd
otherEnd= wlc_ExplicitBindingPair(I)
if (WLC_P__NETWORK) then
    print*, "checkR not implemented for Network"
    stop 1
endif
if (otherEnd == -1) then
    checkR = .False.
    return
endif

if (IT1 <= otherend .and. otherEnd<=IT2) then

    if (norm2(wlc_RP(:,I) - wlc_RP(:,otherEnd)) > eps) then
        print*, "moved apart"
        print*, "I",I,"otherEnd",otherEnd
        print*, " R(I)",wlc_R(:,I)
        print*, " R(o)",wlc_R(:,otherEnd)
        print*, "RP(I)",wlc_RP(:,I)
        print*, "RP(o)",wlc_RP(:,otherEnd)
        print*, "error", norm2(wlc_RP(:,I) - wlc_RP(:,otherEnd))
        checkR = .True.
    endif
else
    if (norm2(wlc_RP(:,I) - wlc_R(:,otherEnd)) > eps) then
        print*, "moved apart"
        print*, "I",I,"otherEnd",otherEnd
        print*, " R(I)",wlc_R(:,I)
        print*, " R(o)",wlc_R(:,otherEnd)
        print*, "RP(I)",wlc_RP(:,I)
        print*, "R(o)",wlc_R(:,otherEnd)
        print*, "error", norm2(wlc_RP(:,I) - wlc_RP(:,otherEnd))
        checkR = .True.
    endif
endif
end function checkR

end module
