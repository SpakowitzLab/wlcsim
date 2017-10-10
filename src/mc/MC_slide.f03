!--------------------------------------------------------------*
!
!           Makes Monti Carlo Moves
!
!    Quinn split out this file on 8/9/17
!
!---------------------------------------------------------------

! variables that need to be allocated only on certain branches moved into MD to prevent segfaults
! please move other variables in as you see fit
subroutine MC_slide(wlc_p,R,U,RP,UP,IP,IB1,IB2,IT1,IT2 &
                  ,MCAMP,WindoW,rand_stat &
                  ,dib)

use mersenne_twister
use params, only: dp,wlcsim_params

implicit none
type(wlcsim_params), intent(in) :: wlc_p
!integer, intent(in) :: wlc_p%NB     ! Number of beads on a polymer
!integer, intent(in) :: wlc_p%NP     ! Number of polymers
!integer, intent(in) :: wlc_p%NT     ! Total beads in simulation
real(dp), intent(in) :: R(3,wlc_p%NT)  ! Bead positions
real(dp), intent(in) :: U(3,wlc_p%NT)  ! Tangent vectors
real(dp), intent(out) :: RP(3,wlc_p%NT)  ! Bead positions
real(dp), intent(out) :: UP(3,wlc_p%NT)  ! Tangent vectors
integer, intent(out) :: IP    ! Test polymer
integer, intent(out) :: IB1   ! Test bead position 1
integer, intent(out) :: IT1   ! Index of test bead 1
integer, intent(out) :: IB2   ! Test bead position 2
integer, intent(out) :: IT2   ! Index of test bead 2
integer, intent(out) :: dib   ! number of beads moved by move (plus or minus a few)
!logical, intent(in) :: wlc_p%ring
!logical, intent(in) :: wlc_p%interp_bead_lennard_jones

integer I,J  ! Test indices
! Things for random number generator
type(random_stat), intent(inout) :: rand_stat  ! status of random number generator
real urand(3)  ! random vector
real urnd(1) ! single random number
integer irnd(1)
real(dp), intent(in) :: MCAMP ! Amplitude of random change
!integer, intent(in) :: wlc_p%wintype
real(dp), intent(in) :: WindoW ! Size of window for bead selection
real(dp) DR(3)    ! Displacement for slide move
integer TEMP
integer exponential_random_int


!TOdo saving RP is not actually needed, even in these cases, but Brad's code assumes that we have RP.
if (wlc_p%ring .OR. wlc_p%interp_bead_lennard_jones) then
    RP = R
    UP = U
endif

!     Perform slide move (MCTYPE 2)

call random_index(wlc_p%NP,irnd,rand_stat)
IP=irnd(1)
call random_index(wlc_p%NB,irnd,rand_stat)
IB1=irnd(1)
if (wlc_p%wintype.eq.0) then
    call random_number(urnd,rand_stat)
    IB2 = IB1 +exponential_random_int(window,rand_stat)
elseif (wlc_p%wintype.eq.1.and..not.wlc_p%ring) then
    call random_number(urand,rand_stat)
    IB2 = IB1 + (2*nint(urand(3))-1)* &
           exponential_random_int(window,rand_stat)
elseif (wlc_p%wintype.eq.1.and.wlc_p%ring) then
    call random_number(urnd,rand_stat)
    IB2 = IB1 + exponential_random_int(window,rand_stat)
else
    call stop_if_err(1, "Warning: wlc_p%wintype not recognized")
endif

DIB = IB2-IB1

if (wlc_p%ring) then
 if (IB2 > wlc_p%NB) then
     IB2 = DIB-(wlc_p%NB-IB1)
 endif
else
 if (IB2 > wlc_p%NB) then
     IB2 = wlc_p%NB
 endif
 if (IB2 < 1) then
    IB2 = 1
 endif
 if (IB2 < IB1) then
     TEMP = IB1
     IB1 = IB2
     IB2 = TEMP
 endif
 IT2 = wlc_p%NB*(IP-1) + IB2
 DIB = IB2-IB1
endif

IT1 = wlc_p%NB*(IP-1) + IB1
IT2 = wlc_p%NB*(IP-1) + IB2

call random_number(urand,rand_stat)
DR(1) = MCAMP*(urand(1)-0.5)
DR(2) = MCAMP*(urand(2)-0.5)
DR(3) = MCAMP*(urand(3)-0.5)

I = IT1
do  J = 0,DIB

   if (I == (wlc_p%NB*IP + 1).AND.wlc_p%ring) then
      I = wlc_p%NB*(IP-1) + 1
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
