#include "../defines.inc"
!--------------------------------------------------------------
!
!     subroutine MC_sterics
!
!     Determine if moved beads intersect old beads.
!     Written by Nicole Pagane, Dec 2019
!     Sphere (nucleosomes) and line (DNA) geometries used
!--------------------------------------------------------------
subroutine MC_sterics(collide,IB1,IB2,IT1,IT2,MCTYPE,forward)
! values from wlcsim_data
use LineLineIntersection, only: LineLineIntersectionCalculation
use SphereLineIntersection, only: SphereLineIntersectionCalculation
use SphereSphereIntersection, only: SphereSphereIntersectionCalculation
use params, only: dp, wlc_RP, wlc_R, wlc_UP, wlc_U, wlc_VP, wlc_V, &
    wlc_basepairs, wlc_nucleosomeWrap
    use nucleosome, only: nucleosomeProp

implicit none
integer, intent(in) :: IB1               ! Test bead position 1
integer, intent(in) :: IT1               ! Index of test bead 1
integer, intent(in) :: IB2               ! Test bead position 2
integer, intent(in) :: IT2               ! Index of test bead 2
integer, intent(in) :: MCTYPE            ! MC move type
logical, intent(in) :: forward           ! direction of reptation move
logical, intent(out) :: collide ! collision variable
logical :: collisionVal
real(dp), parameter :: radius = 5.2 ! nm 
logical :: isNucleosome ! whether or not the moved bead is a nucleosome
logical :: isP1DNA ! whether or not the moved bead +1 is DNA
logical :: isM1DNA ! whether or not the moved bead -1 is DNA
real(dp), dimension(3) :: tempU, tempV, tempR, sphere1, sphere2
integer left, right, ii, jj

collide = .FALSE.

! slide move moves all beads from IT1-IT2
! crank shaft donesn't translate IT2 and IT2
! pivot translates one bead but not the other
! all segments will be kept track of by the bead at their left end
! left = bead to left that doesn't move
! right = right most bead that does move
if (MCTYPE == 4 .or. MCTYPE == 7 .or. MCTYPE == 8 .or. MCTYPE == 9) then
    return ! no collision
elseif (MCTYPE == 1 .or. MCTYPE == 3 ) then ! Crank-shaft or pivot
    if (IT1 == IT2) then
        return ! Bead only rotated; no collision
    endif
    left = IT1
    right = IT2-1
elseif (MCTYPE == 2) then ! Slide
    if (IB1 == 1) then
        left = IT1
    else
        left = IT1 - 1
        !wlc_RP(:,IT1-1) = wlc_R(:,IT1-1) ! need to extend RP
    endif
    if (IB2 == WLC_P__NB) then
        right = IT2 - 1
    else
        right = IT2
        !wlc_RP(:,IT2+1) = wlc_R(:,IT2+1) ! need to extend RP
    endif
elseif (MCTYPE == 5 .or. MCTYPE == 6) then ! Full chain move
    left = IT1
    right = IT2-1
elseif (MCTYPE == 10 .or. MCTYPE == 11) then ! reptation or super rep.
    if (forward) then
        left = IT2-1
    else
        left = IT1
    endif
    right = left
else
    print*, "collision not set up for this movetype ", MCTYPE
    stop
endif

! check collisions for ALL moved beads
do ii = left,right
    ! determine identity of moving bead
    if (wlc_nucleosomeWrap(ii) /= 1) then ! is nucleosome
        isNucleosome = .TRUE.
        ! find end of nucleosome to then find actual midpoint
        !call nucleosomeProp(wlc_UP(:,ii), wlc_VP(:,ii), wlc_RP(:,ii), &
        !            wlc_basepairs(ii),wlc_nucleosomeWrap(ii), &
        !            tempU, tempV, tempR)
        !sphere1 = (wlc_RP(:,ii) + tempR+tempU*WLC_P__LENGTH_PER_BP*wlc_basepairs(ii)) / 2.0
        if (ii /= WLC_P__NT) then 
            sphere1 = (wlc_RP(:,ii) + wlc_RP(:,ii+1)) / 2.0
        else
            sphere1 = wlc_RP(:,ii)
        endif
    else
        isNucleosome = .FALSE.
    endif

    ! determine idenity of neighboring beads (+1)
    if (ii+1 <= WLC_P__NT) then ! on chain
        if (wlc_nucleosomeWrap(ii+1) == 1) then ! P1 is DNA
            isP1DNA = .TRUE.
        else ! P1 is nucloeosme
            isP1DNA = .FALSE.
        endif
    else
        isP1DNA = .FALSE.
    endif

    ! determine idenity of neighboring beads (-1)
    if (ii-1 >= 1) then ! on chain
        if (wlc_nucleosomeWrap(ii-1) == 1) then ! M1 is DNA
            isM1DNA = .TRUE.
        else ! M1 is nucloeosme
            isM1DNA = .FALSE.
        endif
    else
        isM1DNA = .FALSE.
    endif

    ! iterate through all possible interactions 
    ! THIS IS BAD IF CHAIN IS BIG!!!! Then would need to transition to Quinn's code to 
    ! findNeighbors and only check distance of close beads (good for our size chain rn)
    do jj = 1, WLC_P__NT 
        if ( (jj /= ii) ) then !.AND. (jj /= ii-1) .AND. (jj /= ii+1) .AND. (jj /= ii-2) .AND. (jj /= ii+2)) then ! ignore self and neigboring beads
            ! check identity of all other beads in chain 
            if (isNucleosome .AND. wlc_nucleosomeWrap(jj) /= 1) then ! sphere-sphere collision
                ! find end of nucleosome to then find actual midpoint
                !call nucleosomeProp(wlc_U(:,jj), wlc_V(:,jj), wlc_R(:,jj), &
                !            wlc_basepairs(jj),wlc_nucleosomeWrap(jj), &
                !            tempU, tempV, tempR)
                !sphere2 = (wlc_R(:,jj) + tempR+tempU*WLC_P__LENGTH_PER_BP*wlc_basepairs(jj)) / 2.0
                if (jj /= WLC_P__NT ) then 
                    sphere2 = (wlc_RP(:,jj) + wlc_RP(:,jj+1)) / 2.0
                else
                    sphere2 = wlc_RP(:,jj)
                endif
                collisionVal = SphereSphereIntersectionCalculation(sphere1, radius, sphere2, radius)
            ! moved bead nuc + DNA
            else if (isNucleosome .AND. wlc_nucleosomeWrap(jj) == 1) then ! sphere-line collision
                if (jj < WLC_P__NT) then 
                    if (wlc_nucleosomeWrap(jj+1) == 1 .AND. (jj /= ii-1) .AND. (jj /= ii-2)) then 
                        collisionVal = SphereLineIntersectionCalculation(wlc_R(:,jj), wlc_R(:,jj+1), sphere1, radius)
                    endif
                endif
                if (jj > 2) then 
                    if (wlc_nucleosomeWrap(jj-1) == 1 .AND. (jj /= ii+1)) then !.AND. (jj /= i+2)) then
                        collisionVal = SphereLineIntersectionCalculation(wlc_R(:,jj-1), wlc_R(:,jj), sphere1, radius)
                    endif
                endif
            ! moved bead DNA + nuc
            else if ( (isNucleosome .EQV. .FALSE.) .AND. (wlc_nucleosomeWrap(jj) /= 1) ) then ! sphere-line collision
                ! find end of nucleosome to then find actual midpoint
                !call nucleosomeProp(wlc_U(:,jj), wlc_V(:,jj), wlc_R(:,jj), &
                !            wlc_basepairs(jj),wlc_nucleosomeWrap(jj), &
                !            tempU, tempV, tempR)
                !sphere2 = (wlc_R(:,jj) + tempR+tempU*WLC_P__LENGTH_PER_BP*wlc_basepairs(jj)) / 2.0
                if (jj /= WLC_P__NT ) then 
                    sphere2 = (wlc_RP(:,jj) + wlc_RP(:,jj+1)) / 2.0
                else
                    sphere2 = wlc_RP(:,jj)
                endif
                if (isP1DNA .AND. (ii /= jj-1) .AND. (ii /= jj-2)) then 
                    collisionVal = SphereLineIntersectionCalculation(wlc_RP(:,ii), wlc_RP(:,ii+1), sphere2, radius)
                endif
                if (isM1DNA .AND. (ii /= jj+1)) then 
                    collisionVal = SphereLineIntersectionCalculation(wlc_RP(:,ii-1), wlc_RP(:,ii), sphere2, radius)
                endif
            else ! line-line collision
                if (isP1DNA .AND. (jj /= WLC_P__NT)) then 
                    if ((wlc_nucleosomeWrap(jj+1) == 1) .AND. (ii+1 < jj) .AND. (jj+1 < ii)) then
                        collisionVal = LineLineIntersectionCalculation(wlc_RP(:,ii), wlc_RP(:,ii+1), wlc_R(:,jj), wlc_R(:,jj+1))
                    endif
                endif
                if (isP1DNA .AND. (jj /= 1)) then
                    if  ((wlc_nucleosomeWrap(jj-1) == 1) .AND. (ii+1 < jj-1)) then
                        collisionVal = LineLineIntersectionCalculation(wlc_RP(:,ii), wlc_RP(:,ii+1), wlc_R(:,jj-1), wlc_R(:,jj))
                    endif
                endif
                if (isM1DNA .AND. (jj /= WLC_P__NT)) then
                    if ((wlc_nucleosomeWrap(jj+1) == 1) .AND. (jj+1 < ii-1)) then
                        collisionVal = LineLineIntersectionCalculation(wlc_RP(:,ii-1), wlc_RP(:,ii), wlc_R(:,jj), wlc_R(:,jj+1))
                    endif
                endif
                if (isM1DNA .AND. (jj /= 1)) then
                    if ((wlc_nucleosomeWrap(jj-1) == 1) .AND. (ii < jj-1) .AND. (jj < ii-1)) then
                        collisionVal = LineLineIntersectionCalculation(wlc_RP(:,ii-1), wlc_RP(:,ii), wlc_R(:,jj-1), wlc_R(:,jj))
                    endif
                endif
            endif
            if (collisionVal) then 
                collide = .TRUE.
                return ! exit and report collision
            endif
        endif
    enddo
enddo
END
!---------------------------------------------------------------*