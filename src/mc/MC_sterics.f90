#include "../defines.inc"
!--------------------------------------------------------------
!
!     subroutine MC_sterics
!
!     Determine if moved beads intersect old beads.
!     Written by Nicole Pagane, Dec 2019
!     Derived from Quinn's cylinder collision code
!     GJK algorithm used for collision detection
!--------------------------------------------------------------
subroutine MC_sterics(collisions,left,right,ii,nn,neighbors,checkType)
! values from wlcsim_data
use binning
use GJKAlgorithm, only: GJK, constructPolygonPrism
use params, only: dp, wlc_RP, wlc_R, wlc_UP, wlc_U, wlc_VP, wlc_V, &
    wlc_basepairs, wlc_nucleosomeWrap
use nucleosome, only: nucleosomeProp

implicit none
integer, intent(inout) :: collisions
integer, intent(in) :: left, right 
integer, intent(in) :: ii
integer, intent(in) :: nn ! number of neighbors
integer, intent(in) :: neighbors(nn) ! ID of neighboring beads
integer, intent(in) :: checkType ! 0 for not moved beads, 1 for proposed move
logical :: isNucleosome ! whether or not the moved bead is a nucleosome
logical :: isM1DNA ! whether or not the moved bead -1 is DNA
integer, parameter :: s = 12 ! num sides of desired polygon
real(dp), dimension(3) :: tempR, tempU, tempV
real(dp), dimension(s,3) :: poly1Plus, poly2Plus, poly1Minus
integer jj, i

!collisions = 0

! determine identity of moving bead
if (wlc_nucleosomeWrap(ii) /= 1) then ! is nucleosome
    isNucleosome = .TRUE.
else
    isNucleosome = .FALSE.
endif

! determine idenity of neighboring beads (-1)
if (ii-1 >= 1) then ! on chain
    if (wlc_nucleosomeWrap(ii-1) == 1) then ! M1 is DNA
        isM1DNA = .TRUE.
        if (checkType == 1) then 
            ! construct polygon for i-1 to i bead
            if (ii > left ) then 
                poly1Minus = constructPolygonPrism(wlc_RP(:,ii-1), wlc_RP(:,ii), wlc_nucleosomeWrap(ii-1), &
                            wlc_UP(:,ii-1), wlc_VP(:,ii-1), s)
            else
                poly1Minus = constructPolygonPrism(wlc_R(:,ii-1), wlc_RP(:,ii), wlc_nucleosomeWrap(ii-1), &
                            wlc_U(:,ii-1), wlc_V(:,ii-1), s)
            endif
        else
            poly1Minus = constructPolygonPrism(wlc_R(:,ii-1), wlc_R(:,ii), wlc_nucleosomeWrap(ii-1), &
                            wlc_U(:,ii-1), wlc_V(:,ii-1), s)
        endif
    else ! M1 is nucloeosme
        isM1DNA = .FALSE.
        ! construct polygon for i-1 (end of nuc) to i bead
        if ( ii > left .AND. checkType == 1) then 
            call nucleosomeProp(wlc_UP(:,ii-1), wlc_VP(:,ii-1), wlc_RP(:,ii-1), &
                    wlc_basepairs(ii-1),wlc_nucleosomeWrap(ii-1), &
                    tempU, tempV, tempR)
        else
            call nucleosomeProp(wlc_U(:,ii-1), wlc_V(:,ii-1), wlc_R(:,ii-1), &
                    wlc_basepairs(ii-1),wlc_nucleosomeWrap(ii-1), &
                    tempU, tempV, tempR)
        endif
        ! construct polygon from end of nuc to i bead
        if (checkType == 1) then 
            poly1Minus = constructPolygonPrism(tempR, wlc_RP(:,ii), wlc_nucleosomeWrap(ii), &
                        tempU, tempV, s)
        else
            poly1Minus = constructPolygonPrism(tempR, wlc_R(:,ii), wlc_nucleosomeWrap(ii), &
                        tempU, tempV, s)
        endif
    endif
else
    isM1DNA = .FALSE.
endif

! construct polygon for i to i+1 bead
if (checkType == 1) then 
    if (ii < right .AND. ii < WLC_P__NT ) then 
        poly1Plus = constructPolygonPrism(wlc_RP(:,ii), wlc_RP(:,ii+1), wlc_nucleosomeWrap(ii), &
                    wlc_UP(:,ii), wlc_VP(:,ii), s)
    else if (ii < WLC_P__NT) then 
        poly1Plus = constructPolygonPrism(wlc_RP(:,ii), wlc_R(:,ii+1), wlc_nucleosomeWrap(ii), &
                    wlc_UP(:,ii), wlc_VP(:,ii), s)
    endif
else 
    if (ii < WLC_P__NT) then 
        poly1Plus = constructPolygonPrism(wlc_R(:,ii), wlc_R(:,ii+1), wlc_nucleosomeWrap(ii), &
                    wlc_U(:,ii), wlc_V(:,ii), s)
    endif
endif


! iterate through all possible interactions 
do jj = 1, nn
    if (neighbors(jj) < WLC_P__NT) then
        if (checkType == 1) then 
            if (jj+1 < left .OR. jj > right) then 
                poly2Plus = constructPolygonPrism(wlc_R(:,neighbors(jj)), wlc_R(:,neighbors(jj)+1), &
                    wlc_nucleosomeWrap(neighbors(jj)), wlc_U(:,neighbors(jj)), wlc_V(:,neighbors(jj)), s)
            else if (jj+1 == left) then 
                poly2Plus = constructPolygonPrism(wlc_R(:,neighbors(jj)), wlc_RP(:,neighbors(jj)+1), &
                    wlc_nucleosomeWrap(neighbors(jj)), wlc_U(:,neighbors(jj)), wlc_V(:,neighbors(jj)), s)
            else if (jj >= left .AND. jj < right .AND. jj+1 > right) then 
                poly2Plus = constructPolygonPrism(wlc_RP(:,neighbors(jj)), wlc_R(:,neighbors(jj)+1), &
                    wlc_nucleosomeWrap(neighbors(jj)), wlc_UP(:,neighbors(jj)), wlc_VP(:,neighbors(jj)), s)
            else ! both in RP
                poly2Plus = constructPolygonPrism(wlc_RP(:,neighbors(jj)), wlc_RP(:,neighbors(jj)+1), &
                    wlc_nucleosomeWrap(neighbors(jj)), wlc_UP(:,neighbors(jj)), wlc_VP(:,neighbors(jj)), s)
            endif
        else
            poly2Plus = constructPolygonPrism(wlc_R(:,neighbors(jj)), wlc_R(:,neighbors(jj)+1), &
                    wlc_nucleosomeWrap(neighbors(jj)), wlc_U(:,neighbors(jj)), wlc_V(:,neighbors(jj)), s)
        endif
        ! check identity of all other beads in chain 
        if (isNucleosome .AND. wlc_nucleosomeWrap(neighbors(jj)) /= 1 .AND. ii < WLC_P__NT ) then ! sphere-sphere collision
            ! check for collision
            collisions = collisions + 30*GJK(poly1Plus, poly2Plus, s)
        ! moved bead nuc + DNA
        else if (isNucleosome .AND. wlc_nucleosomeWrap(neighbors(jj)) == 1 .AND. ii < WLC_P__NT ) then ! sphere-line collision
            ! ignore 10bp nearest nuc
            if ( (neighbors(jj) < ii .AND. sum(wlc_basepairs(neighbors(jj):ii-1)) >= 10) .OR. &
                (neighbors(jj) >= ii .AND. sum(wlc_basepairs(ii:neighbors(jj))) >= 10) ) then 
                if (wlc_nucleosomeWrap(neighbors(jj)+1) == 1) then
                    ! check for collision
                    collisions = collisions + GJK(poly1Plus, poly2Plus, s)
                endif
            endif 
        ! moved bead DNA + nuc
        else if ( (isNucleosome .EQV. .FALSE.) .AND. (wlc_nucleosomeWrap(neighbors(jj)) /= 1) ) then ! sphere-line collision
            ! ignore 10bp nearest nuc
            if ( (ii < neighbors(jj) .AND. sum(wlc_basepairs(ii:neighbors(jj)-1)) >= 10) .OR. &
                    (ii >= neighbors(jj) .AND. sum(wlc_basepairs(neighbors(jj):ii)) >= 10) ) then 
                ! P1 segment is start of either DNA or nucleosome regardless, can complete line segment
                if (ii < WLC_P__NT ) then 
                    ! check for collision
                    collisions = collisions + GJK(poly1Plus, poly2Plus, s)
                endif
                if (ii > 1 .AND. ii == left) then ! only the first moved bead is allowed to check backwards
                    ! check for collision
                    collisions = collisions + GJK(poly1Minus, poly2Plus, s)
                endif
            endif
        else ! line-line collision
            ! P1 segment is start of either DNA or nucleosome regardless, can complete line segment
            if ((wlc_nucleosomeWrap(neighbors(jj)+1) == 1) .AND. ((ii+1 < neighbors(jj)) &
                .OR. (neighbors(jj)+1 < ii)) .AND. (ii < WLC_P__NT)) then
                ! check for collision
                collisions = collisions + GJK(poly1Plus, poly2Plus, s)
            endif
            ! only the first moved bead is allowed to check backwards in addition to making sure the logic of line segments
            if ((wlc_nucleosomeWrap(neighbors(jj)+1) == 1) .AND. ((neighbors(jj)+1 < ii-1) &
                .OR. (ii < neighbors(jj))) .AND. (ii > 1) .AND. ii == left) then 
                ! check for collision
                collisions = collisions + GJK(poly1Minus, poly2Plus, s)
            endif
        endif
    endif
enddo
END
! ---------------------------------------------------------------*