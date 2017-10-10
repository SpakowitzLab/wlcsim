!--------------------------------------------------------------*
!
!           Makes Monti Carlo Moves
!
!    Quinn separated out this file on 8/9/17
!
!---------------------------------------------------------------

! variables that need to be allocated only on certain branches moved into MD to prevent segfaults
! please move other variables in as you see fit
subroutine MC_chainSwap(wlc_p,R,U,RP,UP,IP,IB1,IB2,IT1,IT2 &
                  ,rand_stat &
                  ,IT3,IT4)

use mersenne_twister
use params, only: dp,wlcsim_params

implicit none
type(wlcsim_params),intent(in) :: wlc_p
!integer, intent(in) :: wlc_p%NB     ! Number of beads on a polymer
!integer, intent(in) :: wlc_p%NP     ! Number of polymers
!integer, intent(in) :: wlc_p%NT     ! Total beads in simulation
real(dp), intent(in) :: R(3,wlc_p%NT)  ! Bead positions
real(dp), intent(in) :: U(3,wlc_p%NT)  ! Tangent vectors
real(dp), intent(out) :: RP(3,wlc_p%NT)  ! Bead positions
real(dp), intent(out) :: UP(3,wlc_p%NT)  ! Tangent vectors
integer, intent(out) :: IP    ! Test polymer
integer IP2   ! Second Test polymer if applicable
integer, intent(out) :: IB1   ! Test bead position 1
integer, intent(out) :: IT1   ! Index of test bead 1
integer, intent(out) :: IB2   ! Test bead position 2
integer, intent(out) :: IT2   ! Index of test bead 2
integer, intent(out) :: IT3   ! Test bead position 3 if applicable
integer, intent(out) :: IT4   ! Test bead position 4 if applicable
!logical, intent(in) :: wlc_p%ring
!logical, intent(in) :: wlc_p%interp_Bead_lennard_jones

! Things for random number generator
type(random_stat), intent(inout) :: rand_stat  ! status of random number generator
integer irnd(1)
integer I


!TOdo saving RP is not actually needed, even in these cases, but Brad's code assumes that we have RP.
if (wlc_p%ring .OR. wlc_p%interp_Bead_lennard_jones) then
    RP = R
    UP = U
endif

! switch two chains
call random_index(wlc_p%NP,irnd,rand_stat)
IP=irnd(1)
call random_index(wlc_p%NP,irnd,rand_stat)
IP2=irnd(1)
! Don't switch a chain with itself
if (IP.eq.IP2) then
    IP2 = IP-1
    if (IP2.eq.0) then
        IP2 = 2
    endif
endif
IT1 = wlc_p%NB*(IP-1) + 1
IT2 = wlc_p%NB*(IP-1) + wlc_p%NB
IT3 = wlc_p%NB*(IP2-1) + 1
IT4 = wlc_p%NB*(IP2-1) + wlc_p%NB
do I = 0,wlc_p%NB-1
   RP(1,IT1 + I) = R(1,IT3 + I)
   RP(2,IT1 + I) = R(2,IT3 + I)
   RP(3,IT1 + I) = R(3,IT3 + I)
   UP(1,IT1 + I) = U(1,IT3 + I)
   UP(2,IT1 + I) = U(2,IT3 + I)
   UP(3,IT1 + I) = U(3,IT3 + I)
   RP(1,IT3 + I) = R(1,IT1 + I)
   RP(2,IT3 + I) = R(2,IT1 + I)
   RP(3,IT3 + I) = R(3,IT1 + I)
   UP(1,IT3 + I) = U(1,IT1 + I)
   UP(2,IT3 + I) = U(2,IT1 + I)
   UP(3,IT3 + I) = U(3,IT1 + I)
ENDdo
IB1 = -2000000
IB2 = -2000000
end subroutine
