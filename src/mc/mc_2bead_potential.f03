#include "../defines.inc"
!---------------------------------------------------------------!
!
!     This subroutine calculates the change in the energy resulting
!     from interbead potentials.
!
!     Written by Quinn 2/27/19
!------------------------------------------------

function two_bead_potential(x,bead1,bead2) result(potential)
    use params, only: dp
    use polydispersity, only: chain_ID
    implicit none
    real(dp) x
    real(dp) potential
    integer, intent(in) :: bead1
    integer, intent(in) :: bead2
    if (abs(bead1 - bead2) <= WLC_P__EXCLUDE_NEIGHBORS_IN_POTENTIAL &
        .and. chain_ID(bead1)==chain_ID(bead2)) then
        potential = 0.0_dp
        return
    endif
    !if (chain_ID(bead1).ne.chain_ID(bead2)) then
    !    potential = 0.0_dp
    !    return
    !endif
    if (x>1.9999999_dp*WLC_P__CHAIN_D) then
        potential = 0.0_dp
        return
    endif


    potential = exp(-x**2/(WLC_P__CHAIN_D**2))
    return
end function

subroutine MC_2bead_potential(MCTYPE)
! values from wlcsim_data
use params, only: wlc_R, wlc_RP, dp, &
     wlc_nPointsMoved, wlc_pointsMoved, nMoveTypes, wlc_bin
use binning, only: binType, findNeighbors, countBeads
use energies, only: energyOf, twoBody_
implicit none
real(dp) two_bead_potential
integer, intent(in) :: MCTYPE
integer ii,jj,kk,ll
integer, parameter  :: maxNeighbors = WLC_P__NT ! equal length of lists
integer nNeighbors !number of neighbors found (so far)
integer neighbors(maxNeighbors) ! list of bead ID's
real(dp) distances(maxNeighbors) ! list of |r-r| values
real(dp) radius, distance
integer otherBead
integer temp

radius = WLC_P__CHAIN_D*2.0_dp  ! 2 because I want to go to exp(-4)

do jj = 1,wlc_nPointsMoved
    ii = wlc_pointsMoved(jj)
    !  ---   Energy between moved beads and unmoved beads ---
    ! plus new
    nNeighbors=0
    call findNeighbors(wlc_bin, wlc_RP(:,ii), radius, wlc_R, &
        WLC_P__NT, maxNeighbors, neighbors, distances, nNeighbors)
    do kk = 1,nNeighbors
        otherBead=neighbors(kk)

        ! vvvvvvvvvvvvvvvvvvvvvvv Testing vvvvvvvvvvvvvvvvv
        if (.not. isnan(wlc_RP(1,otherBead)) .and. .FALSE.) then
            temp=0
            do ll = 1,wlc_nPointsMoved
                if (otherBead == wlc_pointsMoved(ll)) then
                    temp=1
                    exit
                endif
            enddo
            if (temp==0) then
                print*, "Points Moved", wlc_pointsMoved(1:wlc_nPointsMoved)
                print*, "otherBead",otherBead," RP", wlc_RP(1,otherBead), " R", wlc_R(1,otherBead)
                print*, "found error"
                stop
            endif
        endif
        !^^^^^^^^^^^^^^^^^^^^^ end testing ^^^^^^^^^^^^^^^^^^^

        if (.not. isnan(wlc_RP(1,otherBead))) cycle ! exclude moved beads
        distance = distances(kk)
        energyOf(twoBody_)%dx = energyOf(twoBody_)%dx + &
            two_bead_potential(distance, ii, otherBead)
    enddo
    ! minus old
    nNeighbors=0
    call findNeighbors(wlc_bin, wlc_R(:,ii), radius, wlc_R, &
        WLC_P__NT, maxNeighbors, neighbors, distances, nNeighbors)
    do kk = 1,nNeighbors
        otherBead=neighbors(kk)
        if (.not. isnan(wlc_RP(1,otherBead))) cycle ! exclude moved beads
        distance = distances(kk)
        energyOf(twoBody_)%dx = energyOf(twoBody_)%dx - &
            two_bead_potential(distance, ii, otherBead)
    enddo

    !  ---  Between moved beads ---
    do kk = jj, wlc_nPointsMoved
        ll = wlc_pointsMoved(kk)
        distance = norm2(wlc_RP(:,ll)-wlc_RP(:,ii))
        if (distance<=radius) then
            energyOf(twoBody_)%dx = energyOf(twoBody_)%dx +&
                two_bead_potential(distance,ii,ll)
        endif
        distance = norm2(wlc_R(:,ll)-wlc_R(:,ii))
        if (distance<=radius) then
            energyOf(twoBody_)%dx = energyOf(twoBody_)%dx -&
                two_bead_potential(distance,ii,ll)
        endif
    enddo
enddo

RETURN
END

!---------------------------------------------------------------!
subroutine MC_2bead_potential_from_scratch()
! values from wlcsim_data
use params, only: wlc_R
use energies, only: energyOf, twoBody_
use params, only:  dp
use polydispersity, only: chain_ID
implicit none
real(dp) two_bead_potential
integer ii, jj
real(dp) distance, radius

radius = WLC_P__CHAIN_D*2.0_dp  ! 2 because I want to go to exp(-4)

do ii = 1,WLC_P__NT
    do jj = ii, WLC_P__NT
        distance = norm2(wlc_R(:,ii)-wlc_R(:,jj))
        if (distance > radius) cycle
        energyOf(twoBody_)%dx = energyOf(twoBody_)%dx +&
            two_bead_potential(distance, ii, jj)
    enddo
enddo

RETURN
END

!---------------------------------------------------------------!
