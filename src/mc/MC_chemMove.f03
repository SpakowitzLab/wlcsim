!--------------------------------------------------------------*
!
!           Makes Monti Carlo Moves
!
!    Quinn separated out this file on 8/9/17
!
!---------------------------------------------------------------

! variables that need to be allocated only on certain branches moved into MD to prevent segfaults
! please move other variables in as you see fit
subroutine MC_chemMove(wlc_p,R,U,RP,UP,AB,ABP,IP,IB1,IB2,IT1,IT2 &
                  ,WindoW,rand_stat)

use mersenne_twister
use params, only: dp,wlcsim_params

implicit none
type(wlcsim_params), intent(in) :: wlc_p
!integer, intent(in) :: wlc_p%NB     ! Number of beads on a polymer
!integer, intent(in) :: wlc_p%NP     ! Number of polymers
!integer, intent(in) :: wlc_p%NT     ! Total beads in simulation
real(dp), intent(in) :: R(3,wlc_p%NT)  ! Bead positions
real(dp), intent(in) :: U(3,wlc_p%NT)  ! Tangent vectors
integer, intent(in) :: AB(wlc_p%NT)  ! Tangent vectors
real(dp), intent(out) :: RP(3,wlc_p%NT)  ! Bead positions
real(dp), intent(out) :: UP(3,wlc_p%NT)  ! Tangent vectors
integer, intent(out) :: ABP(wlc_p%NT)  ! Tangent vectors
!integer, intent(in) :: wlc_p%nBPM    ! Beads per monomer, aka G
integer, intent(out) :: IP    ! Test polymer
integer, intent(out) :: IB1   ! Test bead position 1
integer, intent(out) :: IT1   ! Index of test bead 1
integer, intent(out) :: IB2   ! Test bead position 2
integer, intent(out) :: IT2   ! Index of test bead 2
!logical, intent(in) :: wlc_p%ring
!logical, intent(in) :: wlc_p%interp_bead_lennard_jones

integer I,J  ! Test indices
! Things for random number generator
type(random_stat), intent(inout) :: rand_stat  ! status of random number generator
real urand(3)  ! random vector
real urnd(1) ! single random number
integer irnd(1)
real(dp), intent(in) :: WindoW ! Size of window for bead selection
integer TEMP
integer exponential_random_int

!TOdo saving RP is not actually needed, even in these cases, but Brad's code assumes that we have RP.
if (wlc_p%ring .OR. wlc_p%interp_bead_lennard_jones) then
    RP = R
    UP = U
endif

! Change wlc_d%AB (a.k.a HP1 binding type fore section of polymer)
! Move amplitude is ignored for this move type
call random_index(wlc_p%NP,irnd,rand_stat)
IP=irnd(1)
call random_index(wlc_p%NB,irnd,rand_stat)
IB1=irnd(1)
call random_number(urand,rand_stat)
IB2 = IB1 + (2*nint(urand(3))-1)* &
        exponential_random_int(window,rand_stat)

if (IB2 < 1) then
   IB2 = 1
endif
if (IB2 > wlc_p%NB) then
   IB2 = wlc_p%NB
endif

if (IB2 < IB1) then
   TEMP = IB1
   IB1 = IB2
   IB2 = TEMP
endif
IT1 = wlc_p%NB*(IP-1) + IB1
IT2 = wlc_p%NB*(IP-1) + IB2

!keep binding constant within monomers
IT1 = IT1-MOD(IT1-1,wlc_p%nBPM)
IT2 = IT2-MOD(IT2-1,wlc_p%nBPM) + wlc_p%nBPM-1

do J = IT1,IT2
    ABP(J) = 1-AB(J)
ENDdo

!This loop may not be necessary
do I = IT1,IT2
   RP(1,I) = R(1,I)
   RP(2,I) = R(2,I)
   RP(3,I) = R(3,I)
   UP(1,I) = U(1,I)
   UP(2,I) = U(2,I)
   UP(3,I) = U(3,I)
ENDdo
end subroutine
