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
integer i, offset1, offset2
real(dp) poly(WLC_P__GJK_POLYGON,3)
integer :: s = WLC_P__GJK_POLYGON 

! only if the MC move moved a bead
if (wlc_nPointsMoved>0) then
    ! adjust extension to left
    if (left == 1) then
        offset1 = 0
    else
        offset1 = -1
    endif
    ! adjust extension to right
    if (right == WLC_P__NT) then 
        offset2 = -1
    else
        offset2 = 0
    endif
    ! set up for collision searching
    RGJK = wlc_R_GJK
    RALL = wlc_R
    UALL = wlc_U
    VALL = wlc_V
    collisions = 0
    ! check for neighbors on old beads
    do i = left+offset1, right+offset2
        call findNeighbors(RGJK(:,i),2*WLC_P__GJK_RADIUS,RGJK,WLC_P__NT-1,WLC_P__NT,neighbors,distances,nn)
        ! check for collisions
        call sterics_check(collisions,RALL,UALL,VALL,left+offset1,i,nn,neighbors(1:nn),distances(1:nn))
    enddo
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
    do i = left+offset1, right+offset2
        call findNeighbors(RGJK(:,i),2*WLC_P__GJK_RADIUS,RGJK,WLC_P__NT-1,WLC_P__NT,neighbors,distances,nn)
        ! check for collisions
        call sterics_check(collisions,RALL,UALL,VALL,left+offset1,i,nn,neighbors(1:nn),distances(1:nn))
    enddo
endif

END subroutine MC_sterics

! sterics check subroutine to check for different types of collisions
subroutine sterics_check(collisions,RALL,UALL,VALL,left,ii,nn,neighbors,distances)
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
logical :: jjGreaterThanii ! as it says
integer, parameter :: s = WLC_P__GJK_POLYGON ! num sides of desired polygon
real(dp), dimension(3) :: tempR, tempU, tempV
real(dp), dimension(s,3) :: poly1Plus, poly2Plus, poly1Minus, poly2Minus
integer jj

! determine identity of moving bead
if (wlc_nucleosomeWrap(ii) /= 1) then ! is nucleosome
    isNucleosome = .TRUE.
else
    isNucleosome = .FALSE.
endif

! construct polygon for i to i+1 bead (should be on virtual bead, so i+1 is safe near end of chain)
poly1Plus = constructPolygonPrism(RALL(:,ii), RALL(:,ii+1), wlc_nucleosomeWrap(ii), &
            UALL(:,ii), VALL(:,ii), s)

! iterate through all possible interactions 
do jj = 1, nn
    if (neighbors(jj) >= left .and. neighbors(jj) <= ii ) cycle
    poly2Plus = constructPolygonPrism(RALL(:,neighbors(jj)), RALL(:,neighbors(jj)+1), &
            wlc_nucleosomeWrap(neighbors(jj)), UALL(:,neighbors(jj)), VALL(:,neighbors(jj)), s)
    if (neighbors(jj) > ii) then 
        jjGreaterThanii = .true.
    else
        jjGreaterThanii = .false.
    endif
    ! check identity of all other beads in chain 
    if (isNucleosome .AND. wlc_nucleosomeWrap(neighbors(jj)) /= 1 ) then ! nuc-nuc collision
        ! check for collision
        if (jjGreaterThanii) then ! ORDER MATTERS
            collisions = collisions + 10*GJK(poly1Plus, poly2Plus, s)
        else
            collisions = collisions + 10*GJK(poly2Plus, poly1Plus, s)
        endif
    ! moved bead nuc + DNA
    else if (isNucleosome .AND. wlc_nucleosomeWrap(neighbors(jj)) == 1 .AND. & 
        distances(jj) < (2*wlc_basepairs(jj)*WLC_P__LENGTH_PER_BP)+WLC_P__GJK_RADIUS) then ! nuc-DNA 
        ! ignore 10bp nearest nuc
        if ( (neighbors(jj) < ii .AND. sum(wlc_basepairs(neighbors(jj):ii-1)) > 10) .OR. &
            (neighbors(jj) > ii .AND. sum(wlc_basepairs(ii:neighbors(jj)-1)) > 10) ) then 
            ! check for collision
            if (jjGreaterThanii) then ! ORDER MATTERS
                collisions = collisions + GJK(poly1Plus, poly2Plus, s)
            else
                collisions = collisions + GJK(poly2Plus, poly1Plus, s)
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
            if (jjGreaterThanii) then ! ORDER MATTERS
                collisions = collisions + GJK(poly1Plus, poly2Plus, s)
            else
                collisions = collisions + GJK(poly2Plus, poly1Plus, s)
            endif
        endif
    else if (distances(jj) < 3*wlc_basepairs(jj)*WLC_P__LENGTH_PER_BP) then ! DNA-DNA collision
        ! P1 segment is start of either DNA or nucleosome regardless, can complete line segment
        if ((ii+1 < neighbors(jj)) .OR. (neighbors(jj)+1 < ii) ) then
            ! check for collision
            if (jjGreaterThanii) then ! ORDER MATTERS
                collisions = collisions + GJK(poly1Plus, poly2Plus, s)
            else
                collisions = collisions + GJK(poly2Plus, poly1Plus, s)
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