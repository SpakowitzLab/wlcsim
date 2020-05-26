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


subroutine MC_sterics(collisions,left,right,netSterics)
! values from wlcsim_data
use params, only: dp, wlc_RP, wlc_UP, wlc_VP, wlc_R, wlc_U, wlc_V, wlc_GJK, &
    wlc_R_GJK, wlc_nucleosomeWrap, wlc_nPointsMoved, wlc_bin, wlc_basepairs, &
    wlc_basepairs_prop
use GJKAlgorithm, only: constructPolygonPrism
use polydispersity, only: get_IP, first_bead_of_chain, last_bead_of_chain
! if using binning, uncomment the next line
!use binning, only: addBead, removeBead, findNeighbors
implicit none

integer, intent(out) :: collisions
integer, intent(in) :: left, right
logical, intent(in) :: netSterics
real(dp) RALL(3,WLC_P__NT) ! all bead R
real(dp) UALL(3,WLC_P__NT) ! all bead U
real(dp) VALL(3,WLC_P__NT) ! all bead V
real(dp) SGJK(WLC_P__GJK_POLYGON,3,WLC_P__NT) ! all vertices for GJK
real(dp) RGJK(3,WLC_P__NT) ! all centers for GJK
real(dp) distances(WLC_P__NT) ! Returned distances
integer neighbors(WLC_P__NT) ! ID of neighboring beads
integer nn ! number of neighbors
integer i, offset1, offset2
integer IP ! chain index 
real(dp) poly(WLC_P__GJK_POLYGON,3)

! only if the MC move moved a bead
! i have commented out the quinn binning implementation and also the probabilistic collision
! acceptance stuff in case we want to switch back. if switch back to quinn code, then you will 
! need to replace all of the instances of findNeighbors to his function and not mine (see EOF)
if (wlc_nPointsMoved>0) then
    IP = get_IP(left) ! can assume left and right are on the same chain
    ! adjust extension to left
    if (left == first_bead_of_chain(IP)) then
        offset1 = 0
    else
        offset1 = -1
    endif
    ! adjust extension to right
    if (right == last_bead_of_chain(IP)) then 
        offset2 = -1
    else
        offset2 = 0
    endif
    ! set up for collision searching
    RGJK = wlc_R_GJK
    SGJK = wlc_GJK
    RALL = wlc_R
    UALL = wlc_U
    VALL = wlc_V
    collisions = 0
    if (netSterics) then 
        ! check for neighbors on old beads
        do i = left+offset1, right+offset2
        !nn = 0
        call findNeighbors(RGJK(:,i),2*WLC_P__GJK_RADIUS,RGJK,WLC_P__NT,WLC_P__NT,neighbors,distances,nn)
        ! check for collisions
        call sterics_check(collisions,RALL,UALL,VALL,SGJK,wlc_basepairs,left+offset1,i,&
                nn,neighbors(1:nn),distances(1:nn),.true.)
        enddo
    endif
    ! replace old beads with new moved beads
    do i = left, right
        ! update real bead locations
        RALL(:,i) = wlc_RP(:,i)
        UALL(:,i) = wlc_UP(:,i)
        VALL(:,i) = wlc_VP(:,i)/norm2(wlc_VP(:,i))
        ! if bead i moves, then remove virtual beads i-1 and i
        ! add back in virtual beads i-1 and i for moved bead i
        if (i > first_bead_of_chain(IP) .AND. i == left) then 
            poly = constructPolygonPrism(wlc_R(:,i-1), wlc_RP(:,i), &
                wlc_nucleosomeWrap(i-1), wlc_U(:,i-1), wlc_V(:,i-1),WLC_P__GJK_POLYGON)
            SGJK(:,:,i-1) = poly
            RGJK(1,i-1) = sum(poly(:,1)/WLC_P__GJK_POLYGON)
            RGJK(2,i-1) = sum(poly(:,2)/WLC_P__GJK_POLYGON)
            RGJK(3,i-1) = sum(poly(:,3)/WLC_P__GJK_POLYGON)
            !if (WLC_P__NEIGHBOR_BINS) then
            !   call removeBead(wlc_bin,wlc_R_GJK(:,i-1),i-1)
            !   call addBead(wlc_bin,RGJK,WLC_P__NT-1,i-1)
            !endif
        endif
        if (i < last_bead_of_chain(IP)) then 
            if (i == right) then 
                poly = constructPolygonPrism(wlc_RP(:,i), wlc_R(:,i+1), &
                    wlc_nucleosomeWrap(i),wlc_UP(:,i), wlc_VP(:,i)/norm2(wlc_VP(:,i)),WLC_P__GJK_POLYGON)
            else
                poly = constructPolygonPrism(wlc_RP(:,i), wlc_RP(:,i+1), &
                    wlc_nucleosomeWrap(i),wlc_UP(:,i), wlc_VP(:,i)/norm2(wlc_VP(:,i)),WLC_P__GJK_POLYGON)
            endif
            SGJK(:,:,i) = poly
            RGJK(1,i) = sum(poly(:,1)/WLC_P__GJK_POLYGON)
            RGJK(2,i) = sum(poly(:,2)/WLC_P__GJK_POLYGON)
            RGJK(3,i) = sum(poly(:,3)/WLC_P__GJK_POLYGON)
            !if (WLC_P__NEIGHBOR_BINS) then
            !   call removeBead(wlc_bin,wlc_R_GJK(:,i),i)
            !   call addBead(wlc_bin,RGJK,WLC_P__NT-1,i)
            !endif
        endif
    enddo
    collisions = -collisions
    ! check for neighbors on new beads
    do i = left+offset1, right+offset2
        !nn = 0
        call findNeighbors(RGJK(:,i),2*WLC_P__GJK_RADIUS,RGJK,WLC_P__NT,WLC_P__NT,neighbors,distances,nn)
        ! check for collisions
        call sterics_check(collisions,RALL,UALL,VALL,SGJK,wlc_basepairs_prop,left+offset1,i,&
                nn,neighbors(1:nn),distances(1:nn),netSterics)
    enddo
    ! if (WLC_P__NEIGHBOR_BINS) then
    !     ! add back in beads if move is rejected
    !     do i = left,right
    !         if (i > 1 .AND. i == left) then 
    !             call removeBead(wlc_bin,RGJK(:,i-1),i-1)
    !             call addBead(wlc_bin,wlc_R_GJK,WLC_P__NT-1,i-1)
    !         endif
    !         if (i < WLC_P__NT) then 
    !             call removeBead(wlc_bin,RGJK(:,i),i)
    !             call addBead(wlc_bin,wlc_R_GJK,WLC_P__NT-1,i)
    !         endif
    !     enddo
    ! endif
endif
END subroutine MC_sterics

! sterics check subroutine to check for different types of collisions
! this can be optimized by more often checking for sterics to quit early
! leaving as is to allow for probabilistic acceptances if desired 
subroutine sterics_check(collisions,RALL,UALL,VALL,SGJK,basepairs,left,ii,nn,neighbors,distances,checkAll)
! values from wlcsim_data
use GJKAlgorithm, only: GJK, constructPolygonPrism
use params, only: dp, wlc_basepairs, wlc_nucleosomeWrap, wlc_pointsMoved, wlc_nPointsMoved
use nucleosome, only: nucleosomeProp
use polydispersity, only: get_IP, first_bead_of_chain, last_bead_of_chain
implicit none

integer, intent(inout) :: collisions
real(dp), intent(in) :: RALL(3,WLC_P__NT) ! all bead R
real(dp), intent(in) :: UALL(3,WLC_P__NT) ! all bead U
real(dp), intent(in) :: VALL(3,WLC_P__NT) ! all bead V
real(dp), intent(in) :: SGJK(WLC_P__GJK_POLYGON,3,WLC_P__NT) ! all vertices for GJK
real(dp), intent(in) :: basepairs(WLC_P__NT) ! basepair discretization
integer, intent(in) :: left ! leftmost check bead
integer, intent(in) :: ii ! index of moved bead
integer, intent(in) :: nn ! number of neighbors
integer, intent(in) :: neighbors(nn) ! ID of neighboring beads
real(dp), intent(in) :: distances(nn) ! ID of neighboring beads
logical, intent(in) :: checkAll ! whether or not to check all collisions
logical :: iiIsNucleosome, jjIsNucleosome ! whether or not the moved bead is a nucleosome
logical :: jjGreaterThanii ! as it says
integer, parameter :: s = WLC_P__GJK_POLYGON ! num sides of desired polygon
real(dp), dimension(3) :: tempR, tempU, tempV
real(dp), dimension(s,3) :: poly1Plus, poly2Plus, poly1ExitDNA, poly2ExitDNA
integer jj

! determine identity of moving bead
if (wlc_nucleosomeWrap(ii) /= 1) then ! is nucleosome
    iiIsNucleosome = .TRUE.
    call nucleosomeProp(UALL(:,ii),VALL(:,ii),RALL(:,ii),basepairs(ii),wlc_nucleosomeWrap(ii),&
            tempU,tempV,tempR)
    poly1ExitDNA = constructPolygonPrism(tempR, RALL(:,ii+1), wlc_nucleosomeWrap(ii), &
            tempU, tempV, s)
else
    iiIsNucleosome = .FALSE.
endif

! construct polygon for i to i+1 bead (should be on virtual bead, so i+1 is safe near end of chain)
poly1Plus = SGJK(:,:,ii)!constructPolygonPrism(RALL(:,ii), RALL(:,ii+1), wlc_nucleosomeWrap(ii), &
            !UALL(:,ii), VALL(:,ii), s)

! iterate through all possible collisions 
do jj = 1, nn
    if (neighbors(jj) >= left .and. neighbors(jj) <= ii &
        .or. neighbors(jj) == last_bead_of_chain(get_IP(neighbors(jj))) ) cycle
    ! determine identity of potentially collided bead
    if (wlc_nucleosomeWrap(neighbors(jj)) /= 1) then ! is nucleosome
        jjIsNucleosome = .TRUE.
        call nucleosomeProp(UALL(:,neighbors(jj)),VALL(:,neighbors(jj)),RALL(:,neighbors(jj)),&
                basepairs(neighbors(jj)),wlc_nucleosomeWrap(neighbors(jj)),tempU,tempV,tempR)
        poly2ExitDNA = constructPolygonPrism(tempR, RALL(:,neighbors(jj)+1), &
                wlc_nucleosomeWrap(neighbors(jj)), tempU, tempV, s)
    else
        jjIsNucleosome = .FALSE.
    endif
    ! construct forward polygon
    poly2Plus = SGJK(:,:,neighbors(jj))!constructPolygonPrism(RALL(:,neighbors(jj)), RALL(:,neighbors(jj)+1), &
            !wlc_nucleosomeWrap(neighbors(jj)), UALL(:,neighbors(jj)), VALL(:,neighbors(jj)), s)
    ! figure out ordering, THIS IS IMPORTANT FOR CONSISTENT GJK CONVERSION
    if (neighbors(jj) > ii) then 
        jjGreaterThanii = .true.
    else
        jjGreaterThanii = .false.
    endif
    ! check identity of all other beads in chain 
    if (iiIsNucleosome .AND. jjIsNucleosome) then ! nuc i + nuc j 
        ! check for collision of nucleosomes
        if (jjGreaterThanii) then ! ORDER MATTERS
            collisions = collisions + 10*GJK(poly1Plus, poly2Plus, s)
        else
            collisions = collisions + 10*GJK(poly2Plus, poly1Plus, s)
        endif
        ! check for exit DNA collisions
        if (distances(jj) < (2*basepairs(jj)*WLC_P__LENGTH_PER_BP)+WLC_P__GJK_RADIUS) then 
            ! check for collision of nucleosome i with exit DNA j ,
            ! collision of exit DNA i with nuclesome j, and
            ! colllision of exit DNA i and j
            if (jjGreaterThanii) then ! ORDER MATTERS
                collisions = collisions + GJK(poly1Plus, poly2ExitDNA, s)
                collisions = collisions + GJK(poly1ExitDNA, poly2Plus, s)
                collisions = collisions + GJK(poly1ExitDNA, poly2ExitDNA, s)
            else
                collisions = collisions + GJK(poly2ExitDNA, poly1Plus, s)
                collisions = collisions + GJK(poly2Plus, poly1ExitDNA, s)
                collisions = collisions + GJK(poly2ExitDNA, poly1ExitDNA, s)
            endif
        endif
    ! moved bead nuc + DNA
    else if (iiIsNucleosome .AND. (jjIsNucleosome .EQV. .FALSE.) .AND. & 
        distances(jj) < (2*basepairs(jj)*WLC_P__LENGTH_PER_BP)+WLC_P__GJK_RADIUS) then ! nuc i + DNA j 
        ! ignore 10bp nearest nuc
        if ( (neighbors(jj) < ii-1 .AND. sum(basepairs(neighbors(jj)+1:ii-1)) > 10) .OR. &
            (neighbors(jj)-1 > ii .AND. sum(basepairs(ii+1:neighbors(jj)-1)) > 10) ) then 
            ! check for collision of nucleosome i with DNA j,
            ! and for collision of exit DNA i with DNA j
            if (jjGreaterThanii) then ! ORDER MATTERS
                collisions = collisions + GJK(poly1Plus, poly2Plus, s)
                collisions = collisions + GJK(poly1ExitDNA, poly2Plus, s)
            else
                collisions = collisions + GJK(poly2Plus, poly1Plus, s)
                collisions = collisions + GJK(poly2Plus, poly1ExitDNA, s)
            endif
        endif 
    ! moved bead DNA + nuc
    else if ( (iiIsNucleosome .EQV. .FALSE.) .AND. jjIsNucleosome .AND. &
        distances(jj) < (2*basepairs(jj)*WLC_P__LENGTH_PER_BP)+WLC_P__GJK_RADIUS) then ! DNA i  + nuc j
        ! ignore 10bp nearest nuc
        if ( (ii < neighbors(jj)-1 .AND. sum(basepairs(ii+1:neighbors(jj)-1)) > 10) .OR. &
                (ii-1 > neighbors(jj) .AND. sum(basepairs(neighbors(jj)+1:ii-1)) > 10) ) then 
            ! check for collision of DNA i with nucleosome j,
            ! and for collision DNA i with exit DNA j
            if (jjGreaterThanii) then ! ORDER MATTERS
                collisions = collisions + GJK(poly1Plus, poly2Plus, s)
                collisions = collisions + GJK(poly1Plus, poly2ExitDNA, s)
            else
                collisions = collisions + GJK(poly2Plus, poly1Plus, s)
                collisions = collisions + GJK(poly2ExitDNA, poly1Plus, s)
            endif
        endif
    else if (distances(jj) < 3*basepairs(jj)*WLC_P__LENGTH_PER_BP) then ! DNA i  + DNA j
        if ((ii+1 < neighbors(jj)) .OR. (neighbors(jj)+1 < ii) ) then
            ! check for collision
            if (jjGreaterThanii) then ! ORDER MATTERS
                collisions = collisions + GJK(poly1Plus, poly2Plus, s)
            else
                collisions = collisions + GJK(poly2Plus, poly1Plus, s)
            endif
        endif
    endif
    if (collisions > 0 .and. (checkAll .eqv. .false.)) return 
enddo
END subroutine sterics_check 

! ! using this instead of quinn's binning code
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