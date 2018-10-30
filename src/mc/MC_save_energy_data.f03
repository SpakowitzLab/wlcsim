#include "../defines.inc"

subroutine MC_save_energy_data(MCTYPE)
! values from wlcsim_data
use params, only: wlc_DEElas, wlc_DEKap, wlc_DECouple, wlc_DEChi, wlc_DEBind,wlc_DEExplicitBinding, wlc_DEExternalField
implicit none
integer MCTYPE

if (MCTYPE == 12) then
    print*, "---------------------------------"
    print*, "DEElas",wlc_DEElas
    print*, "DEKap",wlc_DEKap
    print*, "DECouple", wlc_DECouple
    print*, "DEChi", wlc_DEChi
    print*, "DEBind", wlc_DEBind
    print*, "DEExplicitBinding", wlc_DEExplicitBinding
    print*, "DEExternalField", wlc_DEExternalField


endif

end subroutine

subroutine MC_discribe_spider(spider_id,dr,PRINT_UNIT)
use params, only: wlc_R, wlc_spiders,dp
implicit none
integer spider_id
real(dp) dr(3)
integer hip, knee, toe, leg_n,section
real(dp) hipR(3),toeR(3),kneeR(3),thigh,shin
integer PRINT_UNIT

write(PRINT_UNIT,*) "Number of sections",wlc_spiders(spider_id)%nSections
do section = 1,wlc_spiders(spider_id)%nSections
    print*, wlc_spiders(spider_id)%sections(:,section)
enddo
write(PRINT_UNIT,*) "dr", dr
write(PRINT_UNIT,*) wlc_spiders(spider_id)%nLegs, "legs"
do leg_n = 1,wlc_spiders(spider_id)%nLegs
    hip = wlc_spiders(spider_id)%legs(1,leg_n)
    knee= wlc_spiders(spider_id)%legs(2,leg_n)
    toe = wlc_spiders(spider_id)%legs(3,leg_n)
    hipR=wlc_R(:,hip)
    kneeR=wlc_R(:,knee)
    toeR=wlc_R(:,toe)
    thigh = norm2(kneeR-hipR)
    shin = norm2(kneeR-toeR)
    write(PRINT_UNIT,*) "hip",hip,"knee",knee,"toe",toe
    write(PRINT_UNIT,*) "hipR",hipR,"kneeR",kneeR,"toeR",toeR
    write(PRINT_UNIT,*) "thigh",thigh,"shin",shin,"extent",norm2(hipR-toeR)
enddo

    


end subroutine
