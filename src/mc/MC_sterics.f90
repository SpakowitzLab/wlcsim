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


subroutine MC_sterics(collisions,left,right)
! values from wlcsim_data
use params, only: dp, wlc_RP, wlc_UP, wlc_VP, wlc_R, wlc_U, wlc_V, &
    wlc_R_GJK, wlc_nucleosomeWrap, wlc_nPointsMoved, wlc_bin
use GJKAlgorithm, only: constructPolygonPrism
! if using binning, uncomment the next line
!use binnning, only: addBead, removeBead, findNeighbors
implicit none

integer, intent(out) :: collisions
integer, intent(in) :: left, right
real(dp) RGJK(3,WLC_P__NT-1) ! all virtual bead locations
real(dp) RALL(3,WLC_P__NT) ! all bead R
real(dp) UALL(3,WLC_P__NT) ! all bead U
real(dp) VALL(3,WLC_P__NT) ! all bead V
real(dp) distances(WLC_P__NT) ! Returned distances
integer neighbors(WLC_P__NT) ! ID of neighboring beads
integer nn ! number of neighbors
integer i, offset1, offset2, oldCollisions
real(dp) poly(WLC_P__GJK_POLYGON,3)
integer :: s = WLC_P__GJK_POLYGON 

! only if the MC move moved a bead
if (wlc_nPointsMoved>0) then
    ! adjust extension to left
    offset1 = -1
    if (left == 1) then
        offset1 = 0
    endif
    ! adjust extension to right
    offset2 = 0
    if (right == WLC_P__NT) then 
        offset2 = -1
    endif
    ! set up for collision searching
    RGJK = wlc_R_GJK
    RALL = wlc_R
    UALL = wlc_U
    VALL = wlc_V
    collisions = 0
    ! check for neighbors on old beads
    print*, 'OLD START'
    do i = left+offset1, right+offset2
        call findNeighbors(RGJK(:,i),2*WLC_P__GJK_RADIUS,RGJK,WLC_P__NT-1,WLC_P__NT,neighbors,distances,nn)
        ! check for collisions
        call sterics_check(collisions,RALL,UALL,VALL,left+offset1,i,nn,neighbors(1:nn),distances(1:nn),.TRUE.)
    enddo
    print*, 'OLD END'
    if (collisions > 0) then 
        print*, 'HERE', left, right, collisions
        stop
    endif
    ! replace old beads with new moved beads
    do i = left, right
        ! update real bead locations
        RALL(:,i) = wlc_RP(:,i)
        UALL(:,i) = wlc_UP(:,i)
        VALL(:,i) = wlc_VP(:,i)/norm2(wlc_VP(:,i))
        ! if bead i moves, then remove virtual beads i-1 and i
        ! add back in virtual beads i-1 and i for moved bead i
        if (i > 1 .AND. i == left) then 
            poly = constructPolygonPrism(wlc_R(:,i-1), wlc_RP(:,i), &
                wlc_nucleosomeWrap(i-1), wlc_U(:,i-1), wlc_V(:,i-1), s)
            RGJK(1,i-1) = sum(poly(:,1))/s
            RGJK(2,i-1) = sum(poly(:,2))/s
            RGJK(3,i-1) = sum(poly(:,3))/s
        endif
        if (i < WLC_P__NT) then 
            if (i == right) then 
                poly = constructPolygonPrism(wlc_RP(:,i), wlc_R(:,i+1), &
                    wlc_nucleosomeWrap(i),wlc_UP(:,i), wlc_VP(:,i)/norm2(wlc_VP(:,i)), s)
            else
                poly = constructPolygonPrism(wlc_RP(:,i), wlc_RP(:,i+1), &
                    wlc_nucleosomeWrap(i),wlc_UP(:,i), wlc_VP(:,i)/norm2(wlc_VP(:,i)), s)
            endif
            RGJK(1,i) = sum(poly(:,1))/s
            RGJK(2,i) = sum(poly(:,2))/s
            RGJK(3,i) = sum(poly(:,3))/s
        endif
    enddo
    collisions = -collisions
    ! check for neighbors on new beads
    print*, 'NEW START'
    do i = left+offset1, right+offset2
        call findNeighbors(RGJK(:,i),2*WLC_P__GJK_RADIUS,RGJK,WLC_P__NT-1,WLC_P__NT,neighbors,distances,nn)
        ! check for collisions
        call sterics_check(collisions,RALL,UALL,VALL,left+offset1,i,nn,neighbors(1:nn),distances(1:nn),.TRUE.)
    enddo
    print*, 'NEW END', left, right
endif

END subroutine MC_sterics

! sterics check subroutine to check for different types of collisions
subroutine sterics_check(collisions,RALL,UALL,VALL,left,ii,nn,neighbors,distances,debug)
! values from wlcsim_data
use GJKAlgorithm, only: GJK, constructPolygonPrism
use params, only: dp, wlc_basepairs, wlc_nucleosomeWrap, wlc_pointsMoved, wlc_nPointsMoved
use nucleosome, only: nucleosomeProp
implicit none

integer, intent(inout) :: collisions
real(dp), intent(in) :: RALL(3,WLC_P__NT) ! all bead R
real(dp), intent(in) :: UALL(3,WLC_P__NT) ! all bead U
real(dp), intent(in) :: VALL(3,WLC_P__NT) ! all bead V
integer, intent(in) :: left ! leftmost check bead
integer, intent(in) :: ii ! index of moved bead
integer, intent(in) :: nn ! number of neighbors
integer, intent(in) :: neighbors(nn) ! ID of neighboring beads
real(dp), intent(in) :: distances(nn) ! ID of neighboring beads
logical :: isNucleosome ! whether or not the moved bead is a nucleosome
logical :: isM1ii ! whether or not the -1 ii moved bead is a nucleosome
logical :: isM1jj ! whether or not the -1 ii moved bead is a nucleosome
integer, parameter :: s = WLC_P__GJK_POLYGON ! num sides of desired polygon
real(dp), dimension(3) :: tempR, tempU, tempV
real(dp), dimension(s,3) :: poly1Plus, poly2Plus, poly1Minus, poly2Minus
integer jj
logical, intent(in) :: debug ! printing debugging statements

! determine identity of moving bead
if (wlc_nucleosomeWrap(ii) /= 1) then ! is nucleosome
    isNucleosome = .TRUE.
else
    isNucleosome = .FALSE.
endif

! determine idenity of neighboring beads (-1)
isM1ii = .false.
if (ii > 1) then ! -1 on chain
    if (wlc_nucleosomeWrap(ii-1) /= 1) then ! M1 is nucleosome
        isM1ii = .true.
        ! construct polygon for i-1 (end of nuc) to i bead
        call nucleosomeProp(UALL(:,ii-1), VALL(:,ii-1), RALL(:,ii-1), &
                wlc_basepairs(ii-1),wlc_nucleosomeWrap(ii-1), &
                tempU, tempV, tempR)
        ! construct polygon from end of nuc to i bead
        poly1Minus = constructPolygonPrism(tempR, RALL(:,ii), wlc_nucleosomeWrap(ii), &
                    tempU, tempV, s)
    endif
endif

! construct polygon for i to i+1 bead (should be on virtual bead, so i+1 is safe near end of chain)
poly1Plus = constructPolygonPrism(RALL(:,ii), RALL(:,ii+1), wlc_nucleosomeWrap(ii), &
            UALL(:,ii), VALL(:,ii), s)

! iterate through all possible interactions 
do jj = 1, nn
    if (neighbors(jj) >= left .and. neighbors(jj) <= ii ) cycle
    ! neighbors(jj) will always be one less than NT since it is the index of the nearby virtual bead
    ! determine idenity of neighboring beads (-1)
    isM1jj = .false.
    if (neighbors(jj) > 1) then ! -1 on chain
        if (wlc_nucleosomeWrap(neighbors(jj)-1) /= 1) then ! M1 is nucleosome
            isM1jj = .true.
            ! construct polygon for j-1 (end of nuc) to j bead
            call nucleosomeProp(UALL(:,neighbors(jj)-1), VALL(:,neighbors(jj)-1), RALL(:,neighbors(jj)-1), &
                    wlc_basepairs(neighbors(jj)-1),wlc_nucleosomeWrap(neighbors(jj)-1), &
                    tempU, tempV, tempR)
            ! construct polygon from end of nuc to j bead
            poly2Minus = constructPolygonPrism(tempR, RALL(:,neighbors(jj)), wlc_nucleosomeWrap(neighbors(jj)), &
                        tempU, tempV, s)
        endif
    endif
    poly2Plus = constructPolygonPrism(RALL(:,neighbors(jj)), RALL(:,neighbors(jj)+1), &
            wlc_nucleosomeWrap(neighbors(jj)), UALL(:,neighbors(jj)), VALL(:,neighbors(jj)), s)
    ! check identity of all other beads in chain 
    if (isNucleosome .AND. wlc_nucleosomeWrap(neighbors(jj)) /= 1 ) then ! nuc-nuc collision
        ! check for collision
        collisions = collisions + 10*GJK(poly1Plus, poly2Plus, s)
        print*, 'tried nuc-nuc',  ii, neighbors(jj), left, distances(jj)
        !print*, poly1Plus
        !print*, poly2Plus
        if (GJK(poly1Plus, poly2Plus, s) > 0 .AND. debug) then 
            print*, 'nuc-nuc', ii, neighbors(jj), left, distances(jj)
        endif
    ! moved bead nuc + DNA
    else if (isNucleosome .AND. wlc_nucleosomeWrap(neighbors(jj)) == 1 .AND. & 
        distances(jj) < (2*wlc_basepairs(jj)*WLC_P__LENGTH_PER_BP)+WLC_P__GJK_RADIUS) then ! nuc-DNA 
        ! ignore 10bp nearest nuc
        if ( (neighbors(jj) < ii .AND. sum(wlc_basepairs(neighbors(jj):ii-1)) > 10) .OR. &
            (neighbors(jj) > ii .AND. sum(wlc_basepairs(ii:neighbors(jj)-1)) > 10) ) then 
            ! check for collision
            print*, 'tried nuc-dna forwards',  ii, neighbors(jj), left, distances(jj)
            !print*, poly1Plus
            !print*, poly2Plus
            collisions = collisions + GJK(poly1Plus, poly2Plus, s)
            if (GJK(poly1Plus, poly2Plus, s) > 0 .AND. debug) then 
                print*, 'nuc-dna forwards', ii, neighbors(jj), left, distances(jj)
            endif
            ! only the -1 j nuc can check back
            if (isM1jj) then 
                ! check for collision
                print*, 'tried nuc-dna backwards',  ii, neighbors(jj), left, distances(jj)
                !print*, poly1Plus
                !print*, poly2Minus
                collisions = collisions + GJK(poly1Plus, poly2Minus, s)
                if (GJK(poly1Plus, poly2Minus, s) > 0 .AND. debug) then 
                  print*, 'nuc-dna backwards', ii, neighbors(jj), left, distances(jj)
                endif
            endif
        endif 
    ! moved bead DNA + nuc
    else if ( (isNucleosome .EQV. .FALSE.) .AND. (wlc_nucleosomeWrap(neighbors(jj)) /= 1) .AND. &
        distances(jj) < (2*wlc_basepairs(jj)*WLC_P__LENGTH_PER_BP)+WLC_P__GJK_RADIUS) then ! DNA-nuc
        ! ignore 10bp nearest nuc
        if ( (ii < neighbors(jj) .AND. sum(wlc_basepairs(ii:neighbors(jj)-1)) > 10) .OR. &
                (ii > neighbors(jj) .AND. sum(wlc_basepairs(neighbors(jj):ii-1)) > 10) ) then 
            ! P1 segment is start of either DNA or nucleosome regardless, can complete line segment
            ! check for collision
            print*, 'tried dna-nuc forwards',  ii, neighbors(jj), left, distances(jj)
            !print*, poly1Plus
            !print*, poly2Plus
            collisions = collisions + GJK(poly1Plus, poly2Plus, s)
            if (GJK(poly1Plus, poly2Plus, s) > 0 .AND. debug) then 
                print*, 'dna-nuc forwards', ii, neighbors(jj), left, distances(jj)
            endif
            ! only the -1 i nuc can check back
            if (isM1ii) then 
                ! check for collision
                print*, 'tried dna-nuc backwards',  ii, neighbors(jj), left, distances(jj)
                !print*, poly1Minus
                !print*, poly2Plus
                collisions = collisions + GJK(poly1Minus, poly2Plus, s)
                if (GJK(poly1Minus, poly2Plus, s) > 0 .AND. debug) then 
                  print*, 'dna-nuc backwards', ii, neighbors(jj), left, distances(jj)
                endif
            endif
        endif
    else if (distances(jj) < 3*wlc_basepairs(jj)*WLC_P__LENGTH_PER_BP) then ! DNA-DNA collision
        ! P1 segment is start of either DNA or nucleosome regardless, can complete line segment
        if ((ii+1 < neighbors(jj)) .OR. (neighbors(jj)+1 < ii) ) then
            ! check for collision
            print*, 'tried dna-dna i,j forwards',  ii, neighbors(jj), left, distances(jj)
            !print*, poly1Plus
            !print*, poly2Plus
            collisions = collisions + GJK(poly1Plus, poly2Plus, s)
            if (GJK(poly1Plus, poly2Plus, s) > 0 .AND. debug) then 
                print*, 'dna-dna i,j forwards', ii, neighbors(jj), left, distances(jj)
            endif
        endif
        ! only the -1 nuc i can check back in addition to making sure the logic of line segments
        if (isM1ii .AND. ((neighbors(jj)+1 < ii-1) .OR. (ii < neighbors(jj)))) then 
            ! check for collision
            print*, 'tried dna-dna i back, j for',  ii, neighbors(jj), left, distances(jj)
            !print*, poly1Minus
            !print*, poly2Plus
            collisions = collisions + GJK(poly1Minus, poly2Plus, s)
            if (GJK(poly1Minus, poly2Plus, s) > 0 .AND. debug) then 
                print*, 'dna-dna i back, j for', ii,neighbors(jj), left, distances(jj)
            endif
        endif
        ! only the -1 nuc j can check back in addition to making sure the logic of line segments
        if (isM1jj .AND. ((neighbors(jj) < ii) .OR. (ii+1 < neighbors(jj)-1))) then 
            ! check for collision
            print*, 'tried dna-dna i for, j back',  ii, neighbors(jj), left, distances(jj)
            !print*, poly1Plus
            !print*, poly2Minus
            collisions = collisions + GJK(poly1Plus, poly2Minus, s)
            if (GJK(poly1Plus, poly2Minus, s) > 0 .AND. debug) then 
                print*, 'dna-dna i for, j back', ii,neighbors(jj), left, distances(jj)
            endif
        endif
        ! both the -1 nuc j,i can check back in addition to making sure the logic of line segments
        if (isM1jj .AND. isM1ii .AND. ((neighbors(jj) < ii-1) .OR. (ii < neighbors(jj)-1))) then 
            ! check for collision
            print*, 'tried dna-dna i,j backwards',  ii, neighbors(jj), left, distances(jj)
            !print*, poly1Minus
            !print*, poly2Minus
            collisions = collisions + GJK(poly1Minus, poly2Minus, s)
            if (GJK(poly1Minus, poly2Minus, s) > 0 .AND. debug) then 
                print*, 'dna-dna i,j backwards', ii,neighbors(jj), left, distances(jj)
            endif
        endif
    endif
enddo
END subroutine sterics_check 

! using this instead of quinn's binning code
subroutine findNeighbors(pos,radius,beads,nBeads,neighboringMax,neighbors,distances,nn)
use params, only: dp
use vector_utils, only: distance
implicit none

real(dp), intent(in) :: pos(3)
real(dp), intent(in) :: radius
integer, intent(in) :: nBeads
real(dp), intent(in) :: beads(3,nBeads)
integer, intent(in) :: neighboringMax
integer, intent(out) :: neighbors(neighboringMax)
real(dp), intent(out) :: distances(neighboringMax)
integer, intent(out) :: nn
integer i
real(dp) dist

nn = 0
do i = 1,nBeads
    dist = distance(pos,beads(:,i))
    if (dist <= radius) then 
        nn = nn + 1
        neighbors(nn) = i
        distances(nn) = dist
    endif
enddo

end subroutine findNeighbors
! ---------------------------------------------------------------*