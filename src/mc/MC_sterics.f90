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
subroutine MC_sterics(collisions,IB1,IB2,IT1,IT2,MCTYPE,forward)
! values from wlcsim_data
use binning, only: addBead, removeBead, findNeighbors
use params, only: dp, wlc_RP, wlc_R, wlc_UP, wlc_U, wlc_VP, wlc_V, &
    wlc_bin, wlc_bendPoints, wlc_nBend, wlc_R_GJK, wlc_nucleosomeWrap
use polydispersity, only: length_of_chain, chain_ID, leftmost_from, is_right_end, rightmost_from
use GJKAlgorithm, only: constructPolygonPrism
implicit none
integer, intent(out) :: collisions
integer, intent(in) :: IB1               ! Test bead position 1
integer, intent(in) :: IT1               ! Index of test bead 1
integer, intent(in) :: IB2               ! Test bead position 2
integer, intent(in) :: IT2               ! Index of test bead 2
integer, intent(in) :: MCTYPE            ! MC move type
logical, intent(in) :: forward           ! direction of reptation move
real(dp) R(3,WLC_P__NT-1) ! all bead locations
real(dp) distances(1000) ! Returned distances
real(dp) :: radius = 2.0*WLC_P__NUCLEOSOME_RADIUS ! nm
integer neighbors(1000) ! ID of neighboring beads
integer nn ! number of neighbors
integer left, right, leftExclude, rightExclude, i
real(dp) poly1M(2,3)
real(dp) poly1P(2,3)
integer :: s = 2 ! this is just to get the center of the shape

if ( MCTYPE == 7 .or. MCTYPE == 8 .or. MCTYPE == 9) then
    left = -1
elseif (MCTYPE == 4) then ! rotate  (matters since not spherically symetric)
    left = IT1
    right = left
elseif (MCTYPE == 1 .or. MCTYPE == 3) then ! Crank-shaft or pivot
    left = IT1
    right = IT2-1
    if (IB1 == 1) then
        leftExclude = left
    else
        leftExclude = left -1
    endif
    if (is_right_end(IT2)) then
        rightExclude = right
    else
        rightExclude = right + 1
    endif
elseif (MCTYPE == 2) then ! Slide
    if (IB1 == 1) then
        left = IT1
    else
        left = IT1 - 1
        !wlc_RP(:,IT1-1) = wlc_R(:,IT1-1) ! need to extend RP
    endif
    if (IB1 <= 2) then
        leftExclude = left
    else
        leftExclude = left -1
    endif
    if (IB2 == WLC_P__NB) then
        right = IT2 - 1
    else
        right = IT2
        !wlc_RP(:,IT2+1) = wlc_R(:,IT2+1) ! need to extend RP
    endif
    if (IB2 >= rightmost_from(IT2)-1) then
        rightExclude = right
    else
        rightExclude = right + 1
    endif
    !leftExclude = -1
    !rightExclude = -1
elseif (MCTYPE == 5 .or. MCTYPE == 6) then ! Full chain move
    left = IT1
    right = IT2-1
    leftExclude = left
    rightExclude = right
elseif (MCTYPE == 10 .or. MCTYPE == 11) then ! reptation or super rep.
    if (forward) then
        left = IT2-1
    else
        left = IT1
    endif
    leftExclude = left ! Only exclude self
    rightExclude = left
    right = left
else
    print*, "collision not set up for this movetype ", MCTYPE
    stop
endif

if (left /= -1 ) then 
    !left = wlc_bendPoints(1)
    !right = wlc_bendPoints(wlc_nBend)
    R = wlc_R_GJK
    collisions = 0
    ! check for neighbors on old beads
    ! do i = left, right
    !     nn = 0
    !     call findNeighbors(wlc_bin,wlc_R_GJK(:,i),radius,wlc_R_GJK,WLC_P__NT-1,1000,neighbors,distances,nn)
    !     ! check for collisions
    !     call sterics_check(collisions,left,right,i,nn,neighbors(1:nn),0)
    ! enddo
    ! replace old beads with new moved beads
    do i = left, right
        !call removeBead(wlc_bin,wlc_R(:,i),i)
        ! if bead i moves, then remove virtual beads i-1 and i
        if (i > 1 .AND. i == left) then 
            call removeBead(wlc_bin,wlc_R_GJK(:,i-1),i-1)
        endif
        if (i < WLC_P__NT) then 
            call removeBead(wlc_bin,wlc_R_GJK(:,i),i)
        endif
        !R(:,i) = wlc_RP(:,i)
        !call addBead(wlc_bin,R,WLC_P__NT,i)
        ! add back in virtual beads i-1 and i for moved bead i
        if (i > 1 .AND. i == left) then 
            poly1M = constructPolygonPrism(wlc_R(:,i-1), wlc_RP(:,i), &
                wlc_nucleosomeWrap(i-1), wlc_U(:,i-1), wlc_V(:,i-1), s)
            R(1,i-1) = sum(poly1M(:,1))/s
            R(2,i-1) = sum(poly1M(:,2))/s
            R(3,i-1) = sum(poly1M(:,3))/s
            call addBead(wlc_bin,R,WLC_P__NT-1,i-1)
        else if (i > 1) then 
            poly1M = constructPolygonPrism(wlc_RP(:,i-1), wlc_RP(:,i), &
                wlc_nucleosomeWrap(i-1), wlc_UP(:,i-1), wlc_VP(:,i-1), s)
        endif
        if (i < WLC_P__NT) then
            if (i == right) then 
                poly1P = constructPolygonPrism(wlc_RP(:,i), wlc_R(:,i+1), &
                    wlc_nucleosomeWrap(i),wlc_UP(:,i), wlc_VP(:,i), s)
            else
                poly1P = constructPolygonPrism(wlc_RP(:,i), wlc_RP(:,i+1), &
                    wlc_nucleosomeWrap(i),wlc_UP(:,i), wlc_VP(:,i), s)
            endif
            R(1,i) = sum(poly1P(:,1))/s
            R(2,i) = sum(poly1P(:,2))/s
            R(3,i) = sum(poly1P(:,3))/s
            call addBead(wlc_bin,R,WLC_P__NT-1,i)
        endif
    enddo
    !collisions = -collisions
    ! check for neighbors on new beads
    do i = left, right
        nn = 0
        call findNeighbors(wlc_bin,R(:,i),radius,R,WLC_P__NT-1,1000,neighbors,distances,nn)
        ! check for collisions
        call sterics_check(collisions,left,right,i,nn,neighbors(1:nn),1)
    enddo
    ! add back beads here in case move is rejected
    do i = left, right
        if (i > 1 .AND. i == left) then 
            call removeBead(wlc_bin,R(:,i-1),i-1)
            call addBead(wlc_bin,wlc_R_GJK,WLC_P__NT-1,i-1)
        endif
        if (i < WLC_P__NT) then 
            call removeBead(wlc_bin,R(:,i),i)
            call addBead(wlc_bin,wlc_R_GJK,WLC_P__NT-1,i)
        endif
    enddo
endif

END subroutine MC_sterics

! sterics check subroutine to check for different types of collisions
subroutine sterics_check(collisions,left,right,ii,nn,neighbors,checkType)
! values from wlcsim_data
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
integer :: isM1DNA ! 0 is DNA, 1 is nucleosome, 2 is not on chain
integer, parameter :: s = WLC_P__GJK_POLYGON ! num sides of desired polygon
real(dp), dimension(3) :: tempR, tempU, tempV
real(dp), dimension(s,3) :: poly1Plus, poly2Plus, poly1Minus
integer jj, i

! determine identity of moving bead
if (wlc_nucleosomeWrap(ii) /= 1) then ! is nucleosome
    isNucleosome = .TRUE.
else
    isNucleosome = .FALSE.
endif

! determine idenity of neighboring beads (-1)
if (ii-1 >= 1) then ! -1 on chain
    if (wlc_nucleosomeWrap(ii-1) == 1) then ! M1 is DNA
        isM1DNA = 0
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
        isM1DNA = 1
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
    isM1DNA = 2
endif

! construct polygon for i to i+1 bead
if ( ii < WLC_P__NT ) then 
    if (checkType == 1) then 
        if (ii+1 <= right) then 
           poly1Plus = constructPolygonPrism(wlc_RP(:,ii), wlc_RP(:,ii+1), wlc_nucleosomeWrap(ii), &
                       wlc_UP(:,ii), wlc_VP(:,ii), s)
        else 
           poly1Plus = constructPolygonPrism(wlc_RP(:,ii), wlc_R(:,ii+1), wlc_nucleosomeWrap(ii), &
                       wlc_UP(:,ii), wlc_VP(:,ii), s)
        endif
    else 
        poly1Plus = constructPolygonPrism(wlc_R(:,ii), wlc_R(:,ii+1), wlc_nucleosomeWrap(ii), &
                    wlc_U(:,ii), wlc_V(:,ii), s)
    endif
else
    return ! on last bead which SHOULD NOT extend forward
endif

! iterate through all possible interactions 
do jj = 1, nn
    if (neighbors(jj) >= left .and. neighbors(jj) <= ii) cycle
    if (neighbors(jj) < WLC_P__NT) then ! last bead does not extend forward
        if (checkType == 1) then 
            if (jj+1 < left .OR. jj > right) then 
                poly2Plus = constructPolygonPrism(wlc_R(:,neighbors(jj)), wlc_R(:,neighbors(jj)+1), &
                    wlc_nucleosomeWrap(neighbors(jj)), wlc_U(:,neighbors(jj)), wlc_V(:,neighbors(jj)), s)
            else if (jj+1 == left) then 
                poly2Plus = constructPolygonPrism(wlc_R(:,neighbors(jj)), wlc_RP(:,neighbors(jj)+1), &
                    wlc_nucleosomeWrap(neighbors(jj)), wlc_U(:,neighbors(jj)), wlc_V(:,neighbors(jj)), s)
            else if (jj == right) then 
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
        if (isNucleosome .AND. wlc_nucleosomeWrap(neighbors(jj)) /= 1 ) then ! sphere-sphere collision
            ! check for collision
            collisions = collisions + 10*GJK(poly1Plus, poly2Plus, s)
            ! if (GJK(poly1Plus, poly2Plus, s) > 0 ) then 
            !     print*, GJK(poly1Plus, poly2Plus, s), 'nuc-nuc', ii, neighbors(jj)
            ! endif
        ! moved bead nuc + DNA
        else if (isNucleosome .AND. wlc_nucleosomeWrap(neighbors(jj)) == 1) then ! nuc-DNA 
            ! ignore 10bp nearest nuc
            if ( (neighbors(jj) < ii .AND. sum(wlc_basepairs(neighbors(jj):ii-1)) > 10) .OR. &
                (neighbors(jj) >= ii .AND. sum(wlc_basepairs(ii:neighbors(jj))) > 10) ) then 
                ! check for collision
                collisions = collisions + GJK(poly1Plus, poly2Plus, s)
                ! if (GJK(poly1Plus, poly2Plus, s) > 0) then 
                !     print*, GJK(poly1Plus, poly2Plus, s), 'nuc-dna', ii, neighbors(jj)
                ! endif
            endif 
        ! moved bead DNA + nuc
        else if ( (isNucleosome .EQV. .FALSE.) .AND. (wlc_nucleosomeWrap(neighbors(jj)) /= 1) ) then ! nuc-DNA 
            ! ignore 10bp nearest nuc
            if ( (ii < neighbors(jj) .AND. sum(wlc_basepairs(ii:neighbors(jj)-1)) > 10) .OR. &
                    (ii >= neighbors(jj) .AND. sum(wlc_basepairs(neighbors(jj):ii)) > 10) ) then 
                ! P1 segment is start of either DNA or nucleosome regardless, can complete line segment
                ! check for collision
                collisions = collisions + GJK(poly1Plus, poly2Plus, s)
                !if (GJK(poly1Plus, poly2Plus, s) > 0) then 
                !    print*, GJK(poly1Plus, poly2Plus, s), 'dna-nuc forwards', ii, neighbors(jj)
                !endif
                if (ii > 1 .AND. (ii == left .OR. isM1DNA == 1)) then ! only the first moved bead or -1 nuc is allowed to check backwards
                    ! check for collision
                    collisions = collisions + GJK(poly1Minus, poly2Plus, s)
                    !if (GJK(poly1Minus, poly2Plus, s) > 0) then 
                    !   print*, GJK(poly1Minus, poly2Plus, s), 'dna-nuc backwards', ii, neighbors(jj)
                    !endif
                endif
            endif
        else ! DNA-DNA collision
            ! P1 segment is start of either DNA or nucleosome regardless, can complete line segment
            if (((ii+1 < neighbors(jj)) .OR. (neighbors(jj)+1 < ii)) ) then
                ! check for collision
                collisions = collisions + GJK(poly1Plus, poly2Plus, s)
                ! if (GJK(poly1Plus, poly2Plus, s) > 0) then 
                !     print*, GJK(poly1Plus, poly2Plus, s), 'dna-dna forwards', ii, neighbors(jj)
                !     print*, poly1Plus(:,1)
                !     print*, poly1Plus(:,2)
                !     print*, poly1Plus(:,3)
                !     print*, poly2Plus(:,1)
                !     print*, poly2Plus(:,2)
                !     print*, poly2Plus(:,3)
                ! endif
            endif
            ! only the first moved bead is allowed to check backwards in addition to making sure the logic of line segments
            if (((neighbors(jj)+1 < ii-1) .OR. (ii < neighbors(jj))) .AND. (ii > 1) .AND. ii == left) then 
                ! check for collision
                collisions = collisions + GJK(poly1Minus, poly2Plus, s)
                ! if (GJK(poly1Minus, poly2Plus, s) > 0) then 
                !     print*, GJK(poly1Minus, poly2Plus, s), 'dna-dna backwards', ii,neighbors(jj)
                ! endif
            endif
        endif
    endif
enddo
END subroutine sterics_check 
! ---------------------------------------------------------------*