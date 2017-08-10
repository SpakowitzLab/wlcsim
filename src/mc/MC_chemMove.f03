!--------------------------------------------------------------*
!
!           Makes Monti Carlo Moves
!
!    Quinn separated out this file on 8/9/17
!
!---------------------------------------------------------------

! variables that need to be allocated only on certain branches moved into MD to prevent segfaults
! please move other variables in as you see fit
subroutine MC_chemMove(R,U,RP,UP,AB,ABP,NT,NB,NP,IP,IB1,IB2,IT1,IT2,MCTYPE &
                  ,WindoW,BPM,rand_stat &
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
integer, intent(in) :: AB(NT)  ! Tangent vectors
real(dp), intent(out) :: RP(3,NT)  ! Bead positions
real(dp), intent(out) :: UP(3,NT)  ! Tangent vectors
integer, intent(out) :: ABP(NT)  ! Tangent vectors
integer, intent(in) :: BPM    ! Beads per monomer, aka G
integer, intent(out) :: IP    ! Test polymer
integer, intent(out) :: IB1   ! Test bead position 1
integer, intent(out) :: IT1   ! Index of test bead 1
integer, intent(out) :: IB2   ! Test bead position 2
integer, intent(out) :: IT2   ! Index of test bead 2
logical, intent(in) :: ring
logical, intent(in) :: inTERP_BEAD_LENNARD_JONES

integer I,J  ! Test indices
! Things for random number generator
type(random_stat), intent(inout) :: rand_stat  ! status of random number generator
real urand(3)  ! random vector
real urnd(1) ! single random number
! Variables for the crank-shaft move

real(dp) P1(3)    ! Point on rotation line


!     MC adaptation variables

integer, PARAMETER :: moveTypes = 10 ! Number of different move types
integer, intent(in) :: MCTYPE            ! Type of MC move
real(dp), intent(in) :: WindoW(moveTypes) ! Size of window for bead selection
integer TEMP



!TOdo saving RP is not actually needed, even in these cases, but Brad's code assumes that we have RP.
if (RinG .OR. inTERP_BEAD_LENNARD_JONES) then
    RP = R
    UP = U
    P1 = 0.0_dp
endif

! Change wlc_d%AB (a.k.a HP1 binding type fore section of polymer)
! Move amplitude is ignored for this move type
call random_number(urand,rand_stat)
IP = ceiling(urand(1)*NP)
IB1 = ceiling(urand(2)*NB)
call random_number(urnd,rand_stat)
IB2 = IB1 + (2*nint(urand(3))-1)* &
        nint(-1.0*log(urnd(1))*WindoW(MCTYPE))

if (IB2 < 1) then
   IB2 = 1
endif
if (IB2 > NB) then
   IB2 = NB
endif

if (IB2 < IB1) then
   TEMP = IB1
   IB1 = IB2
   IB2 = TEMP
endif
IT1 = NB*(IP-1) + IB1
IT2 = NB*(IP-1) + IB2

!keep binding constant within monomers
IT1 = IT1-MOD(IT1-1,BPM)
IT2 = IT2-MOD(IT2-1,BPM) + BPM-1

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
