#include "../defines.inc"
module umbrella
use precision, only: dp

implicit none
contains


! ------------------
!  Umbrella sample is of the form:
!
!   E = spring_k * ( rxn_coordinate - x_0 )**2
!   E = spring_k * rxn_coordinate**2  -  spring_k * x_0 * rxn_coordinate
!   E = UmbrellaQuardaticCof * rxn_coordinate**2 + UmbrellaCof * rxn_coordinate
!
! --------------------
function setUmbrellaCof(s) result(cof)
    use params, only: wlc_numProcesses
    implicit none
    real(dp) cof, s
    real(dp) Delta_rxn_coord, spring_k, x_0
    Delta_rxn_coord = (WLC_P__MAX_RXN_COORD-WLC_P__MIN_RXN_COORD)/(wlc_numProcesses-1.0_dp)
    spring_k = 2.0/(Delta_rxn_coord**2) ! for good mixing
    x_0 = WLC_P__MAX_RXN_COORD*s + WLC_P__MIN_RXN_COORD*(1.0_dp - s)
    cof = -spring_k*x_0
end function
function setUmbrellaQuadraticCof(s) result(cof)
    use params, only: wlc_numProcesses
    implicit none
    real(dp) cof, s
    real(dp) Delta_rxn_coord, spring_k
    Delta_rxn_coord = (WLC_P__MAX_RXN_COORD-WLC_P__MIN_RXN_COORD)/(wlc_numProcesses-1.0_dp)
    spring_k = 2.0/(Delta_rxn_coord**2) ! for good mixing
    cof = spring_k
end function


function boundary_weight(x) result(output)
    implicit none
    real(dp) x, output
    if (x<1.0_dp) then
        output = 1.0_dp
    else if (x<17.0_dp) then
        output = (17.0_dp-x)/16.0_dp
    else
        output = 0.0_dp
    endif
end function
subroutine umbrella_energy_from_scratch()
    use params, only: wlc_R
    use energies, only: energyOf, umbrella_, umbrellaQuadratic_
    implicit none
    real(dp) rxn_coordinate
    integer ii
    ! -------- Calculate rxn coordinate from scratch --------
    rxn_coordinate=0.0_dp
    do ii = 1, WLC_P__NT
        rxn_coordinate = rxn_coordinate + boundary_weight(wlc_R(1,ii))
    enddo
    rxn_coordinate = rxn_coordinate/WLC_P__NT
    ! -------- End caluclate rxn coordinate from scratch ---
    energyOf(umbrella_)%dx = rxn_coordinate
    energyOf(umbrellaQuadratic_)%dx = rxn_coordinate**2
end subroutine
subroutine umbrella_energy()
    use energies, only: energyOf, umbrella_, umbrellaQuadratic_
    use params, only: wlc_R, wlc_RP, wlc_nPointsMoved, wlc_pointsMoved
    implicit none
    real(dp) drxn_coordinate
    integer ii
    ! --------- Calculate change in rxn coordinate ---------
    drxn_coordinate=0.0_dp
    do ii = 1, wlc_nPointsMoved
        drxn_coordinate = drxn_coordinate - boundary_weight(wlc_R(1,wlc_pointsMoved(ii)))
        drxn_coordinate = drxn_coordinate + boundary_weight(wlc_RP(1,wlc_pointsMoved(ii)))
    enddo
    drxn_coordinate = drxn_coordinate/WLC_P__NT
    ! --------- End calculate change in rxn coordinate ----
    energyOf(umbrella_)%dx = drxn_coordinate
    ! d(x^2) = dx^2 + 2*x*dx
    energyOf(umbrellaQuadratic_)%dx = drxn_coordinate**2 + 2.0_dp*energyOf(umbrella_)%x*drxn_coordinate
end subroutine


end module umbrella
