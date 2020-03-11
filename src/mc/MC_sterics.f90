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


subroutine MC_sterics(collisions)
! values from wlcsim_data
use binning, only: addBead, removeBead, findNeighbors
use params, only: dp, wlc_RP, wlc_R, wlc_UP, wlc_U, wlc_VP, wlc_V, &
    wlc_bin, wlc_R_GJK, wlc_nucleosomeWrap, wlc_pointsMoved, wlc_nPointsMoved
use GJKAlgorithm, only: constructPolygonPrism
implicit none

integer, intent(out) :: collisions
real(dp) RGJK(3,WLC_P__NT-1) ! all virtual bead locations
real(dp) RALL(3,WLC_P__NT) ! all bead R
real(dp) UALL(3,WLC_P__NT) ! all bead U
real(dp) VALL(3,WLC_P__NT) ! all bead V
real(dp) distances(1000) ! Returned distances
integer neighbors(1000) ! ID of neighboring beads
integer nn ! number of neighbors
integer left, right, i, offset1, offset2
real(dp) poly(WLC_P__GJK_POLYGON,3)
integer :: s = WLC_P__GJK_POLYGON 

! only if the MC move moved a bead
if (wlc_nPointsMoved>0) then
    left = minval(wlc_pointsMoved(1:wlc_nPointsMoved))
    right = maxval(wlc_pointsMoved(1:wlc_nPointsMoved))
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
    ! set up for virtual bead binning
    RGJK = wlc_R_GJK
    RALL = wlc_R
    collisions = 0
    ! check for neighbors on old beads
    ! do i = left+offset, right
    !     nn = 0
    !     call findNeighbors(wlc_bin,wlc_R_GJK(:,i),2*WLC_P__GJK_RADIUS,wlc_R_GJK,WLC_P__NT-1,1000,neighbors,distances,nn)
    !     ! check for collisions
    !     call sterics_check(collisions,wlc_R,wlc_U,wlc_V,left+offset,right,i,nn,neighbors(1:nn),distances(1:nn),.FALSE.)
    ! enddo
    ! replace old beads with new moved beads
    do i = left, right
        ! if bead i moves, then remove virtual beads i-1 and i
        ! add back in virtual beads i-1 and i for moved bead i
        if (i > 1 .AND. i == left) then 
            call removeBead(wlc_bin,wlc_R_GJK(:,i-1),i-1)
            poly = constructPolygonPrism(wlc_R(:,i-1), wlc_RP(:,i), &
                wlc_nucleosomeWrap(i-1), wlc_U(:,i-1), wlc_V(:,i-1), s)
            RGJK(1,i-1) = sum(poly(:,1))/s
            RGJK(2,i-1) = sum(poly(:,2))/s
            RGJK(3,i-1) = sum(poly(:,3))/s
            call addBead(wlc_bin,RGJK,WLC_P__NT-1,i-1)
            ! update bead locations
            RALL(:,i-1) = wlc_R(:,i-1)
            UALL(:,i-1) = wlc_U(:,i-1)
            VALL(:,i-1) = wlc_V(:,i-1)
        endif
        if (i < WLC_P__NT) then 
            call removeBead(wlc_bin,wlc_R_GJK(:,i),i)
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
            call addBead(wlc_bin,RGJK,WLC_P__NT-1,i)
            ! update bead locations
            RALL(:,i) = wlc_RP(:,i)
            UALL(:,i) = wlc_UP(:,i)
            VALL(:,i) = wlc_VP(:,i)/norm2(wlc_VP(:,i))
        endif
    enddo
    collisions = -collisions
    ! check for neighbors on new beads
    do i = left+offset1, right+offset2
        nn = 0
        call findNeighbors(wlc_bin,RGJK(:,i),2*WLC_P__GJK_RADIUS,RGJK,WLC_P__NT-1,1000,neighbors,distances,nn)
        ! check for collisions
        call sterics_check(collisions,RALL,UALL,VALL,left+offset1,i,nn,neighbors(1:nn),distances(1:nn),.FALSE.)
    enddo
    ! add back beads here if move is rejected
    do i = left, right
        if (i > 1 .AND. i == left) then 
            call removeBead(wlc_bin,RGJK(:,i-1),i-1)
            call addBead(wlc_bin,wlc_R_GJK,WLC_P__NT-1,i-1)
        endif
        if (i < WLC_P__NT) then 
            call removeBead(wlc_bin,RGJK(:,i),i)
            call addBead(wlc_bin,wlc_R_GJK,WLC_P__NT-1,i)
        endif
    enddo
endif

END subroutine MC_sterics

! sterics check subroutine to check for different types of collisions
subroutine sterics_check(collisions,rall,uall,vall,left,ii,nn,neighbors,distances,debug)
! values from wlcsim_data
use GJKAlgorithm, only: GJK, constructPolygonPrism
use params, only: dp, wlc_basepairs, wlc_nucleosomeWrap, wlc_pointsMoved, wlc_nPointsMoved
use nucleosome, only: nucleosomeProp
implicit none

integer, intent(inout) :: collisions
real(dp), intent(in) :: RALL(3,WLC_P__NT) ! all bead R
real(dp), intent(in) :: UALL(3,WLC_P__NT) ! all bead U
real(dp), intent(in) :: VALL(3,WLC_P__NT) ! all bead V
integer, intent(in) :: left ! leftmost moved bead
integer, intent(in) :: ii ! index of moved bead
integer, intent(in) :: nn ! number of neighbors
integer, intent(in) :: neighbors(nn) ! ID of neighboring beads
real(dp), intent(in) :: distances(nn) ! ID of neighboring beads
logical :: isNucleosome ! whether or not the moved bead is a nucleosome
logical :: isM1ii ! whether or not the -1 moved bead is a nucleosome
logical :: iiNotLast ! whether or not the moved bead is the last bead or not
integer, parameter :: s = WLC_P__GJK_POLYGON ! num sides of desired polygon
real(dp), dimension(3) :: tempR, tempU, tempV
real(dp), dimension(s,3) :: poly1Plus, poly2Plus, poly1Minus
integer jj, i
logical, intent(in) :: debug ! printing debugging statements

if (nn <= 1) then 
    print*, nn, left, ii, collisions
    stop 
endif

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

! construct polygon for i to i+1 bead
if ( ii < WLC_P__NT ) then 
    iiNotLast = .TRUE.
    poly1Plus = constructPolygonPrism(RALL(:,ii), RALL(:,ii+1), wlc_nucleosomeWrap(ii), &
                UALL(:,ii), VALL(:,ii), s)
else
    iiNotLast = .FALSE.
endif

! iterate through all possible interactions 
do jj = 1, nn
    if (neighbors(jj) >= left .and. neighbors(jj) <= ii ) cycle
    ! neighbors(jj) will always be one less than NT since it is the index of the nearby virtual bead
    poly2Plus = constructPolygonPrism(RALL(:,neighbors(jj)), RALL(:,neighbors(jj)+1), &
            wlc_nucleosomeWrap(neighbors(jj)), UALL(:,neighbors(jj)), VALL(:,neighbors(jj)), s)
    ! check identity of all other beads in chain 
    if (isNucleosome .AND. wlc_nucleosomeWrap(neighbors(jj)) /= 1 ) then ! nuc-nuc collision
        ! check for collision
        collisions = collisions + 10*GJK(poly1Plus, poly2Plus, s)
        if (GJK(poly1Plus, poly2Plus, s) > 0 .AND. debug) then 
            print*, 'nuc-nuc', ii, neighbors(jj), left
        endif
    ! moved bead nuc + DNA
    else if (isNucleosome .AND. wlc_nucleosomeWrap(neighbors(jj)) == 1 .AND. & 
        distances(jj) < 1.5*wlc_nucleosomeWrap(jj)*WLC_P__LENGTH_PER_BP+WLC_P__NUCLEOSOME_RADIUS) then ! nuc-DNA 
        ! ignore 10bp nearest nuc
        if ( (neighbors(jj) < ii .AND. sum(wlc_basepairs(neighbors(jj):ii-1)) > 10) .OR. &
            (neighbors(jj) > ii .AND. sum(wlc_basepairs(ii:neighbors(jj)-1)) > 10) ) then 
            ! check for collision
            collisions = collisions + GJK(poly1Plus, poly2Plus, s)
            if (GJK(poly1Plus, poly2Plus, s) > 0 .AND. debug) then 
                print*, 'nuc-dna', ii, neighbors(jj), left
            endif
        endif 
    ! moved bead DNA + nuc
    else if ( (isNucleosome .EQV. .FALSE.) .AND. (wlc_nucleosomeWrap(neighbors(jj)) /= 1) .AND. &
        distances(jj) < 2*wlc_nucleosomeWrap(jj)*WLC_P__LENGTH_PER_BP+WLC_P__GJK_RADIUS) then ! DNA-nuc
        ! ignore 10bp nearest nuc
        if ( (ii < neighbors(jj) .AND. sum(wlc_basepairs(ii:neighbors(jj)-1)) > 10) .OR. &
                (ii > neighbors(jj) .AND. sum(wlc_basepairs(neighbors(jj):ii-1)) > 10) ) then 
            ! P1 segment is start of either DNA or nucleosome regardless, can complete line segment
            ! check for collision
            collisions = collisions + GJK(poly1Plus, poly2Plus, s)
            if (GJK(poly1Plus, poly2Plus, s) > 0 .AND. debug) then 
                print*, 'dna-nuc forwards', ii, neighbors(jj), left
            endif
            ! only the -1 nuc can check back
            if (isM1ii) then 
                ! check for collision
                collisions = collisions + GJK(poly1Minus, poly2Plus, s)
                if (GJK(poly1Minus, poly2Plus, s) > 0 .AND. debug) then 
                  print*, 'dna-nuc backwards', ii, neighbors(jj), left
                endif
            endif
        endif
    else if (distances(jj) < 2.5*wlc_nucleosomeWrap(jj)*WLC_P__LENGTH_PER_BP) then ! DNA-DNA collision
        ! P1 segment is start of either DNA or nucleosome regardless, can complete line segment
        if ((ii+1 < neighbors(jj)) .OR. (neighbors(jj)+1 < ii) ) then
            ! check for collision
            collisions = collisions + GJK(poly1Plus, poly2Plus, s)
            if (GJK(poly1Plus, poly2Plus, s) > 0 .AND. debug) then 
                print*, 'dna-dna forwards', ii, neighbors(jj), left
            endif
        endif
        ! only the -1 nuc can check back in addition to making sure the logic of line segments
        if (((neighbors(jj)+1 < ii-1) .OR. (ii < neighbors(jj))) .AND. isM1ii ) then 
            ! check for collision
            collisions = collisions + GJK(poly1Minus, poly2Plus, s)
            if (GJK(poly1Minus, poly2Plus, s) > 0 .AND. debug) then 
                print*, 'dna-dna backwards', ii,neighbors(jj), left
            endif
        endif
    endif
enddo
END subroutine sterics_check 
! ---------------------------------------------------------------*