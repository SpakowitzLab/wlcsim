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
use binning, only: binType, findNeighbors, countBeads
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

integer, parameter  :: maxNeighbors = WLC_P__NT ! equal length of lists
integer nNeighbors !number of neighbors found (so far)
integer neighbors(maxNeighbors) ! list of bead ID's
real(dp) distances(maxNeighbors) ! list of |r-r| values
integer II,jj
real(dp) radius
real(dp) R_11(3), R_12(3), R_21(3), R_22(3), new_origin(3), R_test(3), R_test2(3)
integer ix, iy, iz
integer left, right
integer leftExclude, rightExclude
integer numberInBin
radius = 1.0_dp*WLC_P__CHAIN_D + 2.2_dp*wlc_p%l0

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
        wlc_d%RP(:,IT2+1) = wlc_d%R(:,IT2+1) ! need to extend RP
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

    if (WLC_P__CONFINETYPE == 'excludedShpereInPeriodic') then
        R_test(1) = modulo(wlc_d%RP(1,II),WLC_P__LBOX_X)
        R_test(2) = modulo(wlc_d%RP(2,II),WLC_P__LBOX_Y)
        R_test(3) = modulo(wlc_d%RP(3,II),WLC_P__LBOX_Z)
        ! Loop over current an nearest periods to check for neighbors
        do ix=-1,1
            if (ix==-1 .and. R_test(1) > radius) cycle
            if (ix==1 .and. R_test(1) < WLC_P__LBOX_X-radius) cycle
            do iy=-1,1
                if (iy==-1 .and. R_test(2) > radius) cycle
                if (iy==1 .and. R_test(2) < WLC_P__LBOX_Y-radius) cycle
                do iz=-1,1
                    if (iz==-1 .and. R_test(3) > radius) cycle
                    if (iz==1 .and. R_test(3) < WLC_P__LBOX_Y-radius) cycle
                    R_test2(1)=R_test(1)+ix*WLC_P__LBOX_X
                    R_test2(2)=R_test(2)+iy*WLC_P__LBOX_Y
                    R_test2(3)=R_test(3)+iz*WLC_P__LBOX_Z
                    call findNeighbors(wlc_d%bin,R_test2,radius,wlc_d%R_period,&
                        WLC_P__NT,maxNeighbors,neighbors,distances,nNeighbors)

                enddo
            enddo
        enddo
    elseif (WLC_P__CONFINETYPE == 'none') then
        call findNeighbors(wlc_d%bin,wlc_d%RP(:,II),radius,wlc_d%R, &
            WLC_P__NT,maxNeighbors,neighbors,distances,nNeighbors)
    else
        print*, "you need to decide whether to use R or R_period"
        stop
    endif

    if (WLC_P__CONFINETYPE == 'excludedShpereInPeriodic') then
        !Put one end of one segment at the center of the repeating box
        !The other end translated from the first
        new_origin(1)=wlc_d%RP(1,ii)-WLC_P__LBOX_X/2.0_DP
        new_origin(2)=wlc_d%RP(2,ii)-WLC_P__LBOX_Y/2.0_DP
        new_origin(3)=wlc_d%RP(3,ii)-WLC_P__LBOX_Z/2.0_DP
        R_11(1)=modulo(wlc_d%RP(1,ii),new_origin(1))
        R_11(2)=modulo(wlc_d%RP(2,ii),new_origin(2))
        R_11(3)=modulo(wlc_d%RP(3,ii),new_origin(3))
        R_12=R_11+wlc_d%RP(:,ii+1)-wlc_d%RP(:,ii)
    endif
    !if (nNeighbors>256) then
    !    print*, "nNeighbors",nNeighbors,": ", distances(1), distances(2), distances(3), "..."
    !    print*, "R_test",R_test
    !    print*, "Radius ",sqrt((R_test(1)-WLC_P__LBOX_X/2.0_DP)**2+&
    !                           (R_test(2)-WLC_P__LBOX_Y/2.0_DP)**2+&
    !                           (R_test(3)-WLC_P__LBOX_Z/2.0_DP)**2)
    !    call countBeads(wlc_d%bin,numberInBin)
    !    print*, "numberInBin",numberInBin
    !endif
    do jj = 1,nNeighbors
        ! Exclude self interactions
        if (neighbors(jj) >= leftExclude .and. neighbors(jj) <= rightExclude) cycle

        if (WLC_P__RING) then
            print*, "Warning, cyclinders not working for ring"
        endif

        if (WLC_P__CONFINETYPE == 'excludedShpereInPeriodic') then
            if (neighbors(jj)==WLC_P__NT) cycle ! no cylinder comming out of last bead

            !Choose one end of the other segment that is in the same repeat
            !Translate to get the other end (which may be in a different repeat
            !This methode only works if segments are <1/4 of box size
            R_21(1)=modulo(wlc_d%R(1,neighbors(jj)),new_origin(1))
            R_21(2)=modulo(wlc_d%R(2,neighbors(jj)),new_origin(2))
            R_21(3)=modulo(wlc_d%R(3,neighbors(jj)),new_origin(3))
            R_22=R_21+wlc_d%R(:,neighbors(jj)+1)-wlc_d%R(:,neighbors(jj))

            collide = cylinders(R_11,R_12,R_21,R_22)
        elseif (WLC_P__CONFINETYPE == 'none') then
            collide = cylinders(wlc_d%RP(:,ii),wlc_d%RP(:,ii+1),&
                                wlc_d%R(:,neighbors(jj)),&
                                wlc_d%R(:,neighbors(jj)+1))
        else
            print*, "you need to decide whether to use R or R_period"
            stop
        endif
        if (collide) return
    enddo
enddo
collide = .FALSE.
RETURN
END

!---------------------------------------------------------------*
