#include "../defines.inc"
subroutine updateR(wlc_d,I)
use params, only: wlcsim_data, dp
implicit none
Type(wlcsim_data), intent(inout) :: wlc_d     ! system allocated data
integer, intent(in) :: I

if (WLC_P__NEIGHBOR_BINS) then
    if (WLC_P__CONFINETYPE == 'excludedShpereInPeriodic') then
        call removeBead(wlc_d%bin,wlc_d%R_period(:,I),I)
    elseif (WLC_P__CONFINETYPE == 'none') then
        ! call removeBead(wlc_d%bin,wlc_d%R(:,I),I)
    else
        print*, "Not an option yet.  See MCsim."
    endif
endif
wlc_d%R(:,I) = wlc_d%RP(:,I)
wlc_d%U(:,I) = wlc_d%UP(:,I)
if (WLC_P__NEIGHBOR_BINS) then
    if (WLC_P__CONFINETYPE == 'excludedShpereInPeriodic') then
        wlc_d%R_period(1,I)=modulo(wlc_d%R(1,I),WLC_P__LBOX_X)
        wlc_d%R_period(2,I)=modulo(wlc_d%R(2,I),WLC_P__LBOX_Y)
        wlc_d%R_period(3,I)=modulo(wlc_d%R(3,I),WLC_P__LBOX_Z)
        call addBead(wlc_d%bin,wlc_d%R_period,WLC_P__NT,I)
    elseif (WLC_P__CONFINETYPE == 'none') then
        call addBead(wlc_d%bin,wlc_d%R,WLC_P__NT,I)
    else
        print*, "Not an option yet.  See MCsim."
    endif
endif

end subroutine updateR
