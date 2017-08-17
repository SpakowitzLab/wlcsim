!--------------------------------------------------------------*
!
!           Makes Monti Carlo Moves
!
!    Quinn split out this file on 8/9/17
!
!---------------------------------------------------------------

! variables that need to be allocated only on certain branches moved into MD to prevent segfaults
! please move other variables in as you see fit
subroutine MC_slide(R,U,RP,UP,NT,NB,NP,IP,IB1,IB2,IT1,IT2 &
                  ,MCAMP,WindoW,rand_stat,winType &
                  ,dib,ring,inTERP_BEAD_LENNARD_JONES)

use mersenne_twister
use params, only: dp

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
integer, intent(out) :: dib   ! number of beads moved by move (plus or minus a few)
logical, intent(in) :: ring
logical, intent(in) :: inTERP_BEAD_LENNARD_JONES

integer I,J  ! Test indices
! Things for random number generator
type(random_stat), intent(inout) :: rand_stat  ! status of random number generator
real urand(3)  ! random vector
real urnd(1) ! single random number
real(dp), intent(in) :: MCAMP ! Amplitude of random change
integer, intent(in) :: winType
real(dp), intent(in) :: WindoW ! Size of window for bead selection
real(dp) DR(3)    ! Displacement for slide move
integer TEMP


!TOdo saving RP is not actually needed, even in these cases, but Brad's code assumes that we have RP.
if (RinG .OR. inTERP_BEAD_LENNARD_JONES) then
    RP = R
    UP = U
endif

!     Perform slide move (MCTYPE 2)

call random_number(urand,rand_stat)
IP = ceiling(urand(1)*NP)
IB1 = ceiling(urand(2)*NB)
! again, we use a window
if (winType.eq.0) then
    IB2 = IB1 + nint((urand(3)-0.5_dp)*(2.0_dp*WindoW + 1.0))
elseif (winType.eq.1.and..not.RinG) then
    call random_number(urnd,rand_stat)
    IB2 = IB1 + (2*nint(urand(3))-1)* &
            nint(-1.0*log(urnd(1))*WindoW)
elseif (winType.eq.1.and.RinG) then
    call random_number(urnd,rand_stat)
    IB2 = IB1 + nint(-1.0*log(urnd(1))*WindoW)

endif

DIB = IB2-IB1

if (RinG) then
 if (IB2 > NB) then
     IB2 = DIB-(NB-IB1)
 endif
else
 if (IB2 > NB) then
     IB2 = NB
 endif
 if (IB2 < 1) then
    IB2 = 1
 endif
 if (IB2 < IB1) then
     TEMP = IB1
     IB1 = IB2
     IB2 = TEMP
 endif
 IT2 = NB*(IP-1) + IB2
 DIB = IB2-IB1
endif

IT1 = NB*(IP-1) + IB1
IT2 = NB*(IP-1) + IB2

call random_number(urand,rand_stat)
DR(1) = MCAMP*(urand(1)-0.5)
DR(2) = MCAMP*(urand(2)-0.5)
DR(3) = MCAMP*(urand(3)-0.5)

I = IT1
do  J = 0,DIB

   if (I == (NB*IP + 1).AND.RinG) then
      I = NB*(IP-1) + 1
   endif

   RP(1,I) = R(1,I) + DR(1)
   RP(2,I) = R(2,I) + DR(2)
   RP(3,I) = R(3,I) + DR(3)
   UP(1,I) = U(1,I)
   UP(2,I) = U(2,I)
   UP(3,I) = U(3,I)
   I = I + 1

ENDdo
end subroutine
