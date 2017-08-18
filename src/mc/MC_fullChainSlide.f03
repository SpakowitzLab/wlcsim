!--------------------------------------------------------------*
!
!           Makes Monti Carlo Moves
!
!    Quinn separated out this file on 8/9/17
!
!---------------------------------------------------------------

! variables that need to be allocated only on certain branches moved into MD to prevent segfaults
! please move other variables in as you see fit
subroutine MC_fullChainSlide(R,U,RP,UP,NT,NB,NP,IP,IB1,IB2,IT1,IT2 &
                  ,MCAMP,rand_stat &
                  ,ring,inTERP_BEAD_LENNARD_JONES)

use mersenne_twister
use params, only: dp

!TODO: replace R,U,RP,UP .... with wlc_d

implicit none

integer, intent(in) :: NB     ! Number of beads on a polymer
integer, intent(in) :: NP     ! Number of polymers
integer, intent(in) :: NT     ! Total beads in simulation
real(dp), intent(in) :: R(3,NT)  ! Bead positions
real(dp), intent(in) :: U(3,NT)  ! Tangent vectors
real(dp), intent(out) :: RP(3,NT)  ! Bead positions
real(dp), intent(out) :: UP(3,NT)  ! Tangent vectors
integer, intent(out) :: IP    ! Test polymer
integer, intent(out) :: IB1   ! Test bead position 1
integer, intent(out) :: IT1   ! Index of test bead 1
integer, intent(out) :: IB2   ! Test bead position 2
integer, intent(out) :: IT2   ! Index of test bead 2
logical, intent(in) :: ring
logical, intent(in) :: inTERP_BEAD_LENNARD_JONES

integer I ! Test indices
! Things for random number generator
type(random_stat), intent(inout) :: rand_stat  ! status of random number generator
real urand(3)  ! random vector
real urnd(1) ! single random number
integer irnd(1)
! Variables for the crank-shaft move

real(dp) P1(3)    ! Point on rotation line
real(dp), intent(in) :: MCAMP ! Amplitude of random change
real(dp) DR(3)    ! Displacement for slide move


!TOdo saving RP is not actually needed, even in these cases, but Brad's code assumes that we have RP.
if (RinG .OR. inTERP_BEAD_LENNARD_JONES) then
    RP = R
    UP = U
    P1 = 0.0_dp
endif

!     Perform full chain slide move (MCTYPE 6)

call random_index(NP,irnd,rand_stat)
IP=irnd(1)
IB1 = 1
IB2 = NB
IT1 = NB*(IP-1) + IB1
IT2 = NB*(IP-1) + IB2

call random_number(urand,rand_stat)
DR(1) = MCAMP*(urand(1)-0.5_dp)
DR(2) = MCAMP*(urand(2)-0.5_dp)
DR(3) = MCAMP*(urand(3)-0.5_dp)

do I = IT1,IT2
   RP(1,I) = R(1,I) + DR(1)
   RP(2,I) = R(2,I) + DR(2)
   RP(3,I) = R(3,I) + DR(3)
   UP(1,I) = U(1,I)
   UP(2,I) = U(2,I)
   UP(3,I) = U(3,I)
enddo
end subroutine
