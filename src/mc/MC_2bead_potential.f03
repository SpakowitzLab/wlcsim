#include "../defines.inc"
!---------------------------------------------------------------!
!
!     This subroutine calculates the change in the energy resulting
!     from interbead potentials.
!
!     Written by Quinn 2/27/19
!------------------------------------------------

function two_bead_potential(x) result(potential)
    use params, only: dp
    implicit none
    real(dp) x
    real(dp) potential
    potential = exp(-x**2)
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
logical, parameter, dimension(nMoveTypes) :: &
    moved_beads_change_wrt_eachother(nMoveTypes) = &
    (/.False., .False.,&
      .False., .False.,&
      .False., .False.,&
      .True., .True.,&
      .True., .True.,&
      .True., .True./)


radius = WLC_P__CHAIN_D*2.0_dp  ! 2 because I want to go to exp(-4)

do jj = 1,wlc_nPointsMoved
    ii = wlc_pointsMoved(jj)

    !  ---   Energy between moved beads and unmoved beads ---
    ! plus new
    nNeighbors=0
    call findNeighbors(wlc_bin, wlc_RP(:,ii), radius, wlc_R, &
        WLC_P__NT, maxNeighbors, neighbors, distances, nNeighbors)
    do kk = 1,nNeighbors
        if (.not. isnan(wlc_RP(1,neighbors(kk)))) cycle ! exclude moved beads
        if (abs(ii-neighbors(kk)) <= WLC_P__EXCLUDE_NEIGHBORS_IN_POTENTIAL) cycle
        distance = distances(kk)
        energyOf(twoBody_)%dx = energyOf(twoBody_)%dx + two_bead_potential(distance)
    enddo
    ! minus old
    nNeighbors=0
    call findNeighbors(wlc_bin, wlc_R(:,ii), radius, wlc_R, &
        WLC_P__NT, maxNeighbors, neighbors, distances, nNeighbors)
    do kk = 1,nNeighbors
        if (.not. isnan(wlc_RP(1,neighbors(kk)))) cycle ! exclude moved beads
        distance = distances(kk)
        energyOf(twoBody_)%dx = energyOf(twoBody_)%dx - two_bead_potential(distance)
    enddo

    !  ---  Between moved beads ---
    if (moved_beads_change_wrt_eachother(MCTYPE)) then
        do kk = jj, wlc_nPointsMoved
            ll = wlc_pointsMoved(kk)
            if (abs(ii-ll) <= WLC_P__EXCLUDE_NEIGHBORS_IN_POTENTIAL) cycle
            distance = norm2(wlc_RP(:,ll)-wlc_RP(:,ii))
            if (distance<=radius) then
                energyOf(twoBody_)%dx = energyOf(twoBody_)%dx + two_bead_potential(distance)
            endif
            distance = norm2(wlc_R(:,ll)-wlc_R(:,ii))
            if (distance<=radius) then
                energyOf(twoBody_)%dx = energyOf(twoBody_)%dx - two_bead_potential(distance)
            endif
        enddo
    endif
enddo

RETURN
END

!---------------------------------------------------------------!
subroutine MC_2bead_potential_from_scratch()
! values from wlcsim_data
use params, only: wlc_R
use energies, only: energyOf, twoBody_
use params, only:  dp
implicit none
real(dp) two_bead_potential
integer ii, jj
real(dp) distance

do ii = 1,WLC_P__NT-WLC_P__EXCLUDE_NEIGHBORS_IN_POTENTIAL-1
    do jj = ii+WLC_P__EXCLUDE_NEIGHBORS_IN_POTENTIAL+1, WLC_P__NT
        distance = norm2(wlc_R(:,ii)-wlc_R(:,jj))
        energyOf(twoBody_)%dx = energyOf(twoBody_)%dx + two_bead_potential(distance)
    enddo
enddo

RETURN
END

!---------------------------------------------------------------!
