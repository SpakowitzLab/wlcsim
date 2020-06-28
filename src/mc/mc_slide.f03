#include "../defines.inc"
!--------------------------------------------------------------*
!
!           Makes Monti Carlo Moves
!
!    Quinn split out this file on 8/9/17
!
!---------------------------------------------------------------

! variables that need to be allocated only on certain branches moved into MD to prevent segfaults
! please move other variables in as you see fit
subroutine mc_slide(IB1,IB2,IT1,IT2,MCAMP,WindoW,rand_stat,dib,success)
! values from wlcsim_data
use params, only: wlc_RP, wlc_R, wlc_U, wlc_VP, wlc_V, wlc_UP, &
    wlc_nBend, wlc_bendPoints, wlc_pointsMoved, wlc_nPointsMoved, &
    wlc_ExplicitBindingPair, wlc_network_start_index, wlc_other_beads

use mersenne_twister
use params, only: dp
use windowTools, only: draw_window
use polydispersity, only: is_right_end, is_left_end, first_bead_of_chain, last_bead_of_chain, length_of_chain
use ringHelper, only: bend_points_left_ring, bend_points_right_ring, &
    points_moved_left_ring, points_moved_right_ring


implicit none
integer, intent(out) :: IB1   ! Test bead position 1
integer, intent(out) :: IT1   ! Index of test bead 1
integer, intent(out) :: IB2   ! Test bead position 2
integer, intent(out) :: IT2   ! Index of test bead 2
integer, intent(out) :: dib   ! number of beads moved by move (plus or minus a few)
logical, intent(out) :: success

integer IP    ! Test polymer
integer I,J  ! Test indices
integer ii
integer otherEnd
! Things for random number generator
type(random_stat), intent(inout) :: rand_stat  ! status of random number generator
real(dp) urand(3)  ! random vector
real(dp), intent(in) :: MCAMP ! Amplitude of random change
!integer, intent(in) :: WLC_P__WINTYPE
real(dp), intent(in) :: WindoW ! Size of window for bead selection
real(dp) DR(3)    ! Displacement for slide move
integer indx
logical duplicates

!     Perform slide move (MCTYPE 2)
call draw_window(window,WLC_P__MAXWINDOW_SLIDE_MOVE,.false.,rand_stat,&
                IT1,IT2,IB1,IB2,IP,DIB,success)

call random_number(urand,rand_stat)
DR(1) = MCAMP*(urand(1)-0.5_dp)
DR(2) = MCAMP*(urand(2)-0.5_dp)
DR(3) = MCAMP*(urand(3)-0.5_dp)

!  Which elastic segments change
wlc_nBend = 0
if (.not. WLC_P__RING) then                 ! polymer is not a ring
    if (IB1>1) then
        wlc_nBend = wlc_nBend + 1
        wlc_bendPoints(wlc_nBend)=IT1-1
        I=IT1-1
        wlc_RP(:,I)=wlc_R(:,I)
        wlc_UP(:,I)=wlc_U(:,I)
        if (WLC_P__LOCAL_TWIST) wlc_VP(:,I) = wlc_V(:,I)
        wlc_nPointsMoved=wlc_nPointsMoved+1
        wlc_pointsMoved(wlc_nPointsMoved)=I
    endif
    if (IB2<length_of_chain(IP)) then
        wlc_nBend = wlc_nBend + 1
        wlc_bendPoints(wlc_nBend)=IT2
        I=IT2+1
        wlc_RP(:,I)=wlc_R(:,I)
        wlc_UP(:,I)=wlc_U(:,I)
        if (WLC_P__LOCAL_TWIST) wlc_VP(:,I) = wlc_V(:,I)
        wlc_nPointsMoved=wlc_nPointsMoved+1
        wlc_pointsMoved(wlc_nPointsMoved)=I
    endif
else                                        ! polymer is a ring
    if (length_of_chain(IP) - DIB > 2) then
        call bend_points_left_ring(IT1)
        call points_moved_left_ring(IT1)
        call bend_points_right_ring(IT2)
        call points_moved_right_ring(IT2)
    elseif (length_of_chain(IP) - DIB == 2) then
        call bend_points_left_ring(IT1)
        call points_moved_left_ring(IT1)
        call bend_points_right_ring(IT2)
    elseif (length_of_chain(IP) - DIB == 1) then
        call bend_points_left_ring(IT1)
    else
        print*, 'DIB should not take this value', DIB
        stop 1
    endif
endif

! Move Explicity Bound points along with they are bound to
if (WLC_P__EXPLICIT_BINDING) then
    do ii=IT1,IT2
        if (WLC_P__NETWORK) then
            do indx = wlc_network_start_index(ii), &
                          wlc_network_start_index(ii+1)-1
                otherEnd = wlc_other_beads(indx)
                call slide_another_bead(IT1, IT2, otherEnd, DR)
            enddo
        else
            otherEnd = wlc_ExplicitBindingPair(ii)
            if (otherEnd == -1) cycle ! not a loop
            if (otherEnd >= IT1 .and. otherEnd <= IT2) cycle ! internal loop
            call slide_another_bead(IT1, IT2, otherEnd, DR)
        endif
    enddo
endif

I = IT1
do  J = 0,DIB

    if (I == (last_bead_of_chain(IP) + 1).AND.WLC_P__RING) then
       I = first_bead_of_chain(IP)
    endif

    if (WLC_P__WARNING_LEVEL >= 2) then
        if (.not. isnan(wlc_RP(1,I))) then
            print*, "Should be nan still!!"
            stop
        endif
    endif
    wlc_RP(:,I) = wlc_R(:,I) + DR
    wlc_UP(:,I) = wlc_U(:,I)
    if (WLC_P__LOCAL_TWIST) wlc_VP(:,I) = wlc_V(:,I)
    wlc_nPointsMoved=wlc_nPointsMoved+1
    wlc_pointsMoved(wlc_nPointsMoved)=I
    I = I + 1
    if (WLC_P__WARNING_LEVEL >= 2) then
        if (duplicates(wlc_pointsMoved(1:wlc_nPointsMoved), wlc_nPointsMoved)) then
            print*, "Duplicates"
            stop
        endif
        if (duplicates(wlc_bendPoints(1:wlc_nBend), wlc_nBend)) then
            print*, "Duplicates bend"
            stop
        endif
    endif

ENDdo
end subroutine

function duplicates(vec, length) result(output)
    implicit none
    integer length
    integer vec(length)
    integer ii
    logical output
    output = .FALSE.
    if (length <2) then
        return
    endif
    do ii = 1,length-1
        if (ANY( vec(ii+1:length) == vec(ii)) ) then
            output = .TRUE.
            return
        endif
    enddo
end function

subroutine slide_another_bead(IT1,IT2, otherEnd, DR)
use polydispersity, only: is_right_end, is_left_end
use params, only: wlc_RP, wlc_R, wlc_U, wlc_VP, wlc_V, wlc_UP, &
    wlc_nBend, wlc_bendPoints, wlc_pointsMoved, wlc_nPointsMoved, dp
implicit none
integer, intent(in) :: IT1
integer, intent(in) :: IT2
integer, intent(in) :: otherEnd
real(dp), intent(in) :: DR(3)    ! Displacement for slide move
integer I

if (otherEnd >= IT1 .and. otherEnd <= IT2) return

! Move Explicitly bound point along with section
I=otherEnd
if (isnan(wlc_RP(1,I))) then
    wlc_nPointsMoved=wlc_nPointsMoved+1
    wlc_pointsMoved(wlc_nPointsMoved)=I
endif
wlc_RP(:,I)=wlc_R(:,I) + DR
wlc_UP(:,I)=wlc_U(:,I)
if (WLC_P__LOCAL_TWIST) wlc_VP(:,I) = wlc_V(:,I)
! Add ajacent points to RP and bendPoints
if (otherEnd .ne. IT1-1 .and. (.not. is_right_end(otherEnd))) then
    if (.not. ANY(wlc_bendPoints(1:wlc_nBend)==otherEnd)) then
        wlc_nBend=wlc_nBend+1
        wlc_bendPoints(wlc_nBend)=otherEnd
    endif
    I=otherEnd+1
    if (isnan(wlc_RP(1,I))) then
        wlc_RP(:,I)=wlc_R(:,I)
        wlc_UP(:,I)=wlc_U(:,I)
        if (WLC_P__LOCAL_TWIST) wlc_VP(:,I) = wlc_V(:,I)
        wlc_nPointsMoved=wlc_nPointsMoved+1
        wlc_pointsMoved(wlc_nPointsMoved)=I
    endif
endif
if (otherEnd .ne. IT2+1 .and. (.not. is_left_end(otherEnd))) then
    if (.not. ANY(wlc_bendPoints(1:wlc_nBend)==otherEnd-1)) then
        wlc_nBend=wlc_nBend+1
        wlc_bendPoints(wlc_nBend)=otherEnd-1
    endif
    I=otherEnd-1
    if (isnan(wlc_RP(1,I))) then
        wlc_RP(:,I)=wlc_R(:,I)
        wlc_UP(:,I)=wlc_U(:,I)
        if (WLC_P__LOCAL_TWIST) wlc_VP(:,I) = wlc_V(:,I)
        wlc_nPointsMoved=wlc_nPointsMoved+1
        wlc_pointsMoved(wlc_nPointsMoved)=I
    endif
endif

end subroutine


