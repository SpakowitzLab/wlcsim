#include "../defines.inc"
!--------------------------------------------------------------
!
!     subroutine MC_cylinder
!
!     Determine if moved beads intersect old beads.
!     Written by Quinn MacPherson, Dec 2017
!--------------------------------------------------------------
subroutine MC_cylinder(wlc_p,wlc_d,collide,IB1,IB2,&
                    IT1,IT2, &
                    MCTYPE,forward)

use params, only: dp, pi,wlcsim_params, wlcsim_data
use binning, only: binType, findNeighbors 
implicit none
type(wlcsim_params), intent(in) :: wlc_p
Type(wlcsim_data), intent(inout) :: wlc_d
integer, intent(in) :: IB1               ! Test bead position 1
integer, intent(in) :: IT1               ! Index of test bead 1
integer, intent(in) :: IB2               ! Test bead position 2
integer, intent(in) :: IT2               ! Index of test bead 2
integer, intent(in) :: MCTYPE            ! MC move type
logical, intent(in) :: forward           ! direction of reptation move
logical, intent(out) :: collide   ! Change in ECOM
logical cylinders ! function for collision checking

integer, parameter  :: maxNeighbors = 256 ! equal length of lists
integer nNeighbors !number of neighbors found (so far)
integer neighbors(maxNeighbors) ! list of bead ID's
real(dp) distances(maxNeighbors) ! list of |r-r| values
integer II,jj
real(dp) radius 
integer left, right
integer leftExclude, rightExclude
radius = 2.0_dp*WLC_P__CHAIN_D + 2.0_dp*wlc_p%l0

! slide move moves all beads from IT1-IT2
! crank shaft donesn't translate IT2 and IT2
! pivot translates one bead but not the other
! all segments will be kept track of by the bead at their left end
! left = bead to left that doesn't move
! right = right most bead that does move
! left/right-exclude don't check collision with this range of beads

if (MCTYPE == 1 .or. MCTYPE == 3 ) then ! Crank-shaft or pivot 
    if (IT1 == IT2) then
        collide = .FALSE. ! Bead only rotated
        return
    endif
    left = IT1
    right = IT2-1
    if (IB1 == 1) then
        leftExclude = left
    else
        leftExclude = left -1
    endif
    if (IB2 == WLC_P__NBPM) then
        rightExclude = right
    else
        rightExclude = right + 1
    endif
elseif (MCTYPE == 2) then ! Slide
    if (IB1 == 1) then
        left = IT1
    else
        left = IT1 - 1
        wlc_d%RP(:,IT1-1) = wlc_d%R(:,IT1-1) ! need to extend RP
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
        wlc_d%RP(:,IT1+1) = wlc_d%R(:,IT1+1) ! need to extend RP
    endif
    if (IB2 >= WLC_P__NBPM-1) then
        rightExclude = right
    else
        rightExclude = right + 1
    endif
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
    print*, "cylinder collision not set up for this movetype ", MCTYPE
    stop
endif


do ii = left,right
    nNeighbors = 0 ! Clear list of neighbors
    call findNeighbors(wlc_d%bin,wlc_d%RP(:,II),radius,wlc_d%R,wlc_p%NT, &
                   maxNeighbors,neighbors,distances,nNeighbors)
    do jj = 1,nNeighbors
        ! Exclude self interactions
        if (neighbors(jj) >= leftExclude .and. neighbors(jj) <= rightExclude) cycle

        if (WLC_P__RING) then
            print*, "Warning, cyclinders not working for ring"
        endif

        collide = cylinders(wlc_d%RP(:,ii),wlc_d%RP(:,ii+1),&
                            wlc_d%R(:,neighbors(jj)),&
                            wlc_d%R(:,neighbors(jj+1)))
        if (collide) return
    enddo
enddo
collide = .FALSE.
RETURN
END

!---------------------------------------------------------------*
