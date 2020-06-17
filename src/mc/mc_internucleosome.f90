#include "../defines.inc"
subroutine mc_internucleosome(left,right)
!Calculate the change in the polymer internucleosome attraction energy
! the actual energy function is in the nucleosome module file and is 
! based off of work from the de pablo group

! values from wlcsim_data
use params, only: dp, wlc_U, wlc_nucleosomeWrap, wlc_VP, wlc_V&
    , wlc_R, wlc_UP, wlc_RP
use nucleosome, only: internucleosome_energy
use energies, only: energyOf, internucleosome_

implicit none
integer, intent(in) :: left, right ! leftmost and rightmost moved beads
real(dp) delInt     ! change in internucleosome attraction energy
integer i, j

delInt = 0
! will only do internucleosome if sterics is on (probably)
do i = left,right
    if (wlc_nucleosomeWrap(i)==1) cycle
    do j = 1,WLC_P__NT
        if (wlc_nucleosomeWrap(j)==1 .or. (j>=left .and. j<=i)) cycle
            ! old config
            if (j > i) then
                delInt = delInt - internucleosome_energy(wlc_R(:,i),wlc_R(:,j),&
                                                        wlc_U(:,i),wlc_U(:,j),&
                                                        wlc_V(:,i),wlc_V(:,j))
            else
                delInt = delInt - internucleosome_energy(wlc_R(:,j),wlc_R(:,i),&
                                                        wlc_U(:,j),wlc_U(:,i),&
                                                        wlc_V(:,j),wlc_V(:,i))
            endif
        ! new config
        if (i >= left .AND. i <= right) then 
            if (j >= left .AND. j <= right) then ! i in moved, j in moved
                if (j > i) then 
                    delInt = delInt + internucleosome_energy(wlc_RP(:,i),wlc_RP(:,j),&
                                                            wlc_UP(:,i),wlc_UP(:,j),&
                                                            wlc_VP(:,i),wlc_VP(:,j))
                else
                    delInt = delInt + internucleosome_energy(wlc_RP(:,j),wlc_RP(:,i),&
                                                            wlc_UP(:,j),wlc_UP(:,i),&
                                                            wlc_VP(:,j),wlc_VP(:,i))
                endif
            else ! i in moved, j not in moved
                if (j > i) then 
                    delInt = delInt + internucleosome_energy(wlc_RP(:,i),wlc_R(:,j),&
                                                            wlc_UP(:,i),wlc_U(:,j),&
                                                            wlc_VP(:,i),wlc_V(:,j))
                else
                    delInt = delInt + internucleosome_energy(wlc_R(:,j),wlc_RP(:,i),&
                                                            wlc_U(:,j),wlc_UP(:,i),&
                                                            wlc_V(:,j),wlc_VP(:,i))
                endif
            endif
        else 
            if (j >= left .AND. j <= right) then ! i not in moved, j in moved
                if (j > i) then 
                    delInt = delInt + internucleosome_energy(wlc_R(:,i),wlc_RP(:,j),&
                                                            wlc_U(:,i),wlc_UP(:,j),&
                                                            wlc_V(:,i),wlc_VP(:,j))
                else
                    delInt = delInt + internucleosome_energy(wlc_RP(:,j),wlc_R(:,i),&
                                                            wlc_UP(:,j),wlc_U(:,i),&
                                                            wlc_VP(:,j),wlc_V(:,i))
                endif
            else ! i not in moved, j not in moved
                if (j > i ) then 
                    delInt = delInt + internucleosome_energy(wlc_R(:,i),wlc_R(:,j),&
                                                            wlc_U(:,i),wlc_U(:,j),&
                                                            wlc_V(:,i),wlc_V(:,j))
                else
                    delInt = delInt + internucleosome_energy(wlc_R(:,j),wlc_R(:,i),&
                                                            wlc_U(:,j),wlc_U(:,i),&
                                                            wlc_V(:,j),wlc_V(:,i))
                endif
            endif
        endif
    enddo
enddo
energyOf(internucleosome_)%dx = delInt

END subroutine 