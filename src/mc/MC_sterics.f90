#include "../defines.inc"
!--------------------------------------------------------------
!
!     subroutine MC_sterics
!
!     Determine if moved beads intersect old beads.
!     Written by Nicole Pagane, Mar 2020
!     Inspired from Quinn's cylinder collision code
!     GJK algorithm used for collision detection
!--------------------------------------------------------------


subroutine MC_sterics(collisions,IB1,IB2,IT1,IT2,MCTYPE,forward)
! values from wlcsim_data
use binning, only: addBead, removeBead, findNeighbors
use params, only: dp, wlc_RP, wlc_R, wlc_UP, wlc_U, wlc_VP, wlc_V, &
    wlc_bin, wlc_R_GJK, wlc_nucleosomeWrap, wlc_nBend, wlc_bendPoints, wlc_pointsMoved, wlc_nPointsMoved
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
real(dp) :: radius = 2*WLC_P__NUCLEOSOME_RADIUS ! nm
integer neighbors(1000) ! ID of neighboring beads
integer nn ! number of neighbors
integer left, right, i, offset
real(dp) poly(WLC_P__GJK_POLYGON,3)
integer :: s = WLC_P__GJK_POLYGON ! this is just to get the center of the shape

! only if the MC move moved a bead
if (wlc_nPointsMoved>0) then
    left = minval(wlc_pointsMoved(1:wlc_nPointsMoved))
    right = maxval(wlc_pointsMoved(1:wlc_nPointsMoved))
    ! offset -1 left for rotate
    if (MCTYPE == 4) then 
        if (left > 1 ) then 
            left = left - 1
        endif
        if (right < WLC_P__NT) then 
            right = right + 1
        endif
    endif
    ! adjust extension to left
    offset = 0
    if (left > 1) then 
        offset = -1
        wlc_RP(:,left-1) = wlc_R(:,left-1)
        wlc_UP(:,left-1) = wlc_U(:,left-1)
        wlc_VP(:,left-1) = wlc_V(:,left-1)
    endif
    ! set up for virtual bead binning
    R = wlc_R_GJK
    collisions = 0
    ! check for neighbors on old beads
    do i = left+offset, right
        nn = 0
        call findNeighbors(wlc_bin,wlc_R_GJK(:,i),radius,wlc_R_GJK,WLC_P__NT-1,1000,neighbors,distances,nn)
        ! check for collisions
        call sterics_check(collisions,left+offset,right,i,nn,neighbors(1:nn),distances(1:nn),0,.FALSE.)
    enddo
    ! replace old beads with new moved beads
    do i = left, right
        ! if bead i moves, then remove virtual beads i-1 and i
        ! add back in virtual beads i-1 and i for moved bead i
        if (i > 1 .AND. i == left) then 
            call removeBead(wlc_bin,wlc_R_GJK(:,i-1),i-1)
            poly = constructPolygonPrism(wlc_R(:,i-1), wlc_RP(:,i), &
                wlc_nucleosomeWrap(i-1), wlc_U(:,i-1), wlc_V(:,i-1), s)
            R(1,i-1) = sum(poly(:,1))/s
            R(2,i-1) = sum(poly(:,2))/s
            R(3,i-1) = sum(poly(:,3))/s
            call addBead(wlc_bin,R,WLC_P__NT-1,i-1)
        endif
        if (i < WLC_P__NT) then 
            call removeBead(wlc_bin,wlc_R_GJK(:,i),i)
            if (i == right) then 
                poly = constructPolygonPrism(wlc_RP(:,i), wlc_R(:,i+1), &
                    wlc_nucleosomeWrap(i),wlc_UP(:,i), wlc_VP(:,i), s)
            else
                poly = constructPolygonPrism(wlc_RP(:,i), wlc_RP(:,i+1), &
                    wlc_nucleosomeWrap(i),wlc_UP(:,i), wlc_VP(:,i), s)
            endif
            R(1,i) = sum(poly(:,1))/s
            R(2,i) = sum(poly(:,2))/s
            R(3,i) = sum(poly(:,3))/s
            call addBead(wlc_bin,R,WLC_P__NT-1,i)
        endif
    enddo
    collisions = -collisions
    ! check for neighbors on new beads
    do i = left+offset, right
        nn = 0
        call findNeighbors(wlc_bin,R(:,i),radius,R,WLC_P__NT-1,1000,neighbors,distances,nn)
        ! check for collisions
        call sterics_check(collisions,left+offset,right,i,nn,neighbors(1:nn),distances(1:nn),1,.FALSE.)
    enddo
    ! add back beads here if move is rejected
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
subroutine sterics_check(collisions,left,right,ii,nn,neighbors,distances,checkType,debug)
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
real(dp), intent(in) :: distances(nn) ! ID of neighboring beads
integer, intent(in) :: checkType ! 0 for not moved beads, 1 for proposed move
logical :: isNucleosome ! whether or not the moved bead is a nucleosome
integer :: isM1ii ! 0 is DNA, 1 is nucleosome, 2 is not on chain
logical :: iiNotLast ! whether or not the moved bead is the last bead or not
integer, parameter :: s = WLC_P__GJK_POLYGON ! num sides of desired polygon
real(dp), dimension(3) :: tempR, tempU, tempV
real(dp), dimension(s,3) :: poly1Plus, poly2Plus, poly1Minus
integer jj, i
logical, intent(in) :: debug ! printing debugging statements

! determine identity of moving bead
if (wlc_nucleosomeWrap(ii) /= 1) then ! is nucleosome
    isNucleosome = .TRUE.
else
    isNucleosome = .FALSE.
endif

! determine idenity of neighboring beads (-1)
if (ii > 1) then ! -1 on chain
    if (wlc_nucleosomeWrap(ii-1) == 1) then ! M1 is DNA
        isM1ii = 0
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
        isM1ii = 1
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
    isM1ii = 2
endif

! construct polygon for i to i+1 bead
if ( ii < WLC_P__NT ) then 
    iiNotLast = .TRUE.
    if (checkType == 1) then 
        if (ii < right) then 
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
    iiNotLast = .FALSE.
endif

! iterate through all possible interactions 
do jj = 1, nn
    if (neighbors(jj) >= left .and. neighbors(jj) <= ii) cycle
    ! neighbors(jj) will always be one less than NT since it is the index of the nearby virtual bead
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
    if (iiNotLast .AND. isNucleosome .AND. wlc_nucleosomeWrap(neighbors(jj)) /= 1 ) then ! nuc-nuc collision
        ! check for collision
        collisions = collisions + 30*GJK(poly1Plus, poly2Plus, s)
        if (GJK(poly1Plus, poly2Plus, s) > 0 .AND. debug) then 
            print*, GJK(poly1Plus, poly2Plus, s), 'nuc-nuc', ii, neighbors(jj), left, right
        endif
    ! moved bead nuc + DNA
    else if (iiNotLast .AND. isNucleosome .AND. wlc_nucleosomeWrap(neighbors(jj)) == 1 .AND. & 
        distances(jj) < 2*wlc_nucleosomeWrap(jj)*WLC_P__LENGTH_PER_BP+WLC_P__NUCLEOSOME_RADIUS) then ! nuc-DNA 
        ! ignore 10bp nearest nuc
        if ( (neighbors(jj) < ii .AND. sum(wlc_basepairs(neighbors(jj):ii-1)) > 10) .OR. &
            (neighbors(jj) > ii .AND. sum(wlc_basepairs(ii:neighbors(jj)-1)) > 10) ) then 
            ! check for collision
            collisions = collisions + GJK(poly1Plus, poly2Plus, s)
            if (GJK(poly1Plus, poly2Plus, s) > 0 .AND. debug) then 
                print*, GJK(poly1Plus, poly2Plus, s), 'nuc-dna', ii, neighbors(jj), left, right
            endif
        endif 
    ! moved bead DNA + nuc
    else if ( (isNucleosome .EQV. .FALSE.) .AND. (wlc_nucleosomeWrap(neighbors(jj)) /= 1) .AND. &
        distances(jj) < 2*wlc_nucleosomeWrap(jj)*WLC_P__LENGTH_PER_BP+WLC_P__NUCLEOSOME_RADIUS) then ! DNA-nuc
        ! ignore 10bp nearest nuc
        if ( (ii < neighbors(jj) .AND. sum(wlc_basepairs(ii:neighbors(jj)-1)) > 10) .OR. &
                (ii > neighbors(jj) .AND. sum(wlc_basepairs(neighbors(jj):ii-1)) > 10) ) then 
            ! P1 segment is start of either DNA or nucleosome regardless, can complete line segment
            ! check for collision
            if (iiNotLast) then 
                collisions = collisions + GJK(poly1Plus, poly2Plus, s)
                if (GJK(poly1Plus, poly2Plus, s) > 0 .AND. debug) then 
                   print*, GJK(poly1Plus, poly2Plus, s), 'dna-nuc forwards', ii, neighbors(jj), left, right
                endif
            endif
            ! only the -1 nuc or last bead can check back
            if (isM1ii == 1 .OR. (iiNotLast .EQV. .FALSE. .AND. left == i)) then 
                ! check for collision
                collisions = collisions + GJK(poly1Minus, poly2Plus, s)
                if (GJK(poly1Minus, poly2Plus, s) > 0 .AND. debug) then 
                  print*, GJK(poly1Minus, poly2Plus, s), 'dna-nuc backwards', ii, neighbors(jj), left, right
                endif
            endif
        endif
    else if (distances(jj) < 3*wlc_nucleosomeWrap(jj)*WLC_P__LENGTH_PER_BP) then ! DNA-DNA collision
        ! P1 segment is start of either DNA or nucleosome regardless, can complete line segment
        if (iiNotLast .AND. ((ii+1 < neighbors(jj)) .OR. (neighbors(jj)+1 < ii)) ) then
            ! check for collision
            collisions = collisions + GJK(poly1Plus, poly2Plus, s)
            if (GJK(poly1Plus, poly2Plus, s) > 0 .AND. debug) then 
                print*, GJK(poly1Plus, poly2Plus, s), 'dna-dna forwards', ii, neighbors(jj), left, right
                print*, poly1Plus(:,1)
                print*, poly1Plus(:,2)
                print*, poly1Plus(:,3)
                print*, poly2Plus(:,1)
                print*, poly2Plus(:,2)
                print*, poly2Plus(:,3)
            endif
        endif
        ! only the -1 nuc or last bead can check back in addition to making sure the logic of line segments
        if (((neighbors(jj)+1 < ii-1) .OR. (ii < neighbors(jj))) .AND. &
                (isM1ii == 1 .OR. (iiNotLast .EQV. .FALSE. .AND. left == i)) ) then 
            ! check for collision
            collisions = collisions + GJK(poly1Minus, poly2Plus, s)
            if (GJK(poly1Minus, poly2Plus, s) > 0 .AND. debug) then 
                print*, GJK(poly1Minus, poly2Plus, s), 'dna-dna backwards', ii,neighbors(jj), left, right
                print*, poly1Minus(:,1)
                print*, poly1Minus(:,2)
                print*, poly1Minus(:,3)
                print*, poly2Plus(:,1)
                print*, poly2Plus(:,2)
                print*, poly2Plus(:,3)
            endif
        endif
    endif
enddo
END subroutine sterics_check 
! ---------------------------------------------------------------*