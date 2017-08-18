!--------------------------------------------------------------*
!
!           Makes Monti Carlo Moves
!
!    Quinn split out this file on 8/9/17
!
!---------------------------------------------------------------

! variables that need to be allocated only on certain branches moved into MD to prevent segfaults
! please move other variables in as you see fit
subroutine MC_pivot(R,U,RP,UP,NT,NB,NP,IP,IB1,IB2,IT1,IT2 &
                  ,MCAMP,WindoW,rand_stat,winType &
                  ,ring,inTERP_BEAD_LENNARD_JONES)

use mersenne_twister
use params, only: dp, pi

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

integer I  ! Test indices
! Things for random number generator
type(random_stat), intent(inout) :: rand_stat  ! status of random number generator
real urand(3)  ! random vector
real urnd(1) ! single random number
integer irnd(1)

real(dp) TA(3)    ! Axis of rotation
real(dp) P1(3)    ! Point on rotation line
real(dp) ROT(4,4) ! Rotation matrix
real(dp) ALPHA    ! Angle of move
real(dp) BETA     ! Angle of move

!     MC adaptation variables

real(dp), intent(in) :: MCAMP ! Amplitude of random change
integer, intent(in) :: winType
real(dp), intent(in) :: WindoW ! Size of window for bead selection
integer exponential_random_int

!TOdo saving RP is not actually needed, even in these cases, but Brad's code assumes that we have RP.
if (RinG .OR. inTERP_BEAD_LENNARD_JONES) then
    RP = R
    UP = U
    P1 = 0.0_dp
endif

! We don't have to protect moves 4-10 with if ring because the code is identical in both cases
!     Perform pivot move (MCTYPE 3)
    call random_index(NP,irnd,rand_stat)
    IP=irnd(1)
    call random_index(NB,irnd,rand_stat)
    IB1=irnd(1)
    if (urnd(1).gt.0.5_dp) then
        IB2 = exponential_random_int(window,rand_stat) + 1
        if (IB2 > NB) then
            IB2 = NB
        endif
        IB1 = 1
        IT1 = NB*(IP-1) + IB1
        IT2 = NB*(IP-1) + IB2
        P1(1) = R(1,IT2)
        P1(2) = R(2,IT2)
        P1(3) = R(3,IT2)
    else
        IB1 = NB-exponential_random_int(window,rand_stat)
        if (IB1 < 1) then
            IB1 = 1
        endif
        IB2 = NB
        IT1 = NB*(IP-1) + IB1
        IT2 = NB*(IP-1) + IB2
        P1(1) = R(1,IT1)
        P1(2) = R(2,IT1)
        P1(3) = R(3,IT1)
    endif

   call random_number(urand,rand_stat)
   ALPHA = 2.*PI*urand(1)
   BETA = acos(2.*urand(2)-1.)
   TA(1) = sin(BETA)*cos(ALPHA)
   TA(2) = sin(BETA)*sin(ALPHA)
   TA(3) = cos(BETA)

   ALPHA = MCAMP*(urand(3)-0.5)

   ROT(1,1) = TA(1)**2. + (TA(2)**2. + TA(3)**2.)*cos(ALPHA)
   ROT(1,2) = TA(1)*TA(2)*(1.-cos(ALPHA))-TA(3)*sin(ALPHA)
   ROT(1,3) = TA(1)*TA(3)*(1.-cos(ALPHA)) + TA(2)*sin(ALPHA)
   ROT(1,4) = (P1(1)*(1.-TA(1)**2.) &
   -TA(1)*(P1(2)*TA(2) + P1(3)*TA(3)))*(1.-cos(ALPHA)) + (P1(2)*TA(3)-P1(3)*TA(2))*sin(ALPHA)

   ROT(2,1) = TA(1)*TA(2)*(1.-cos(ALPHA)) + TA(3)*sin(ALPHA)
   ROT(2,2) = TA(2)**2. + (TA(1)**2. + TA(3)**2.)*cos(ALPHA)
   ROT(2,3) = TA(2)*TA(3)*(1.-cos(ALPHA))-TA(1)*sin(ALPHA)
   ROT(2,4) = (P1(2)*(1.-TA(2)**2.) &
   -TA(2)*(P1(1)*TA(1) + P1(3)*TA(3)))*(1.-cos(ALPHA)) + (P1(3)*TA(1)-P1(1)*TA(3))*sin(ALPHA)

   ROT(3,1) = TA(1)*TA(3)*(1.-cos(ALPHA))-TA(2)*sin(ALPHA)
   ROT(3,2) = TA(2)*TA(3)*(1.-cos(ALPHA)) + TA(1)*sin(ALPHA)
   ROT(3,3) = TA(3)**2. + (TA(1)**2. + TA(2)**2.)*cos(ALPHA)
   ROT(3,4) = (P1(3)*(1.-TA(3)**2.) &
   -TA(3)*(P1(1)*TA(1) + P1(2)*TA(2)))*(1.-cos(ALPHA)) + (P1(1)*TA(2)-P1(2)*TA(1))*sin(ALPHA)

   do I = IT1,IT2
      RP(1,I) = ROT(1,4) + ROT(1,1)*R(1,I) + ROT(1,2)*R(2,I) + ROT(1,3)*R(3,I)
      RP(2,I) = ROT(2,4) + ROT(2,1)*R(1,I) + ROT(2,2)*R(2,I) + ROT(2,3)*R(3,I)
      RP(3,I) = ROT(3,4) + ROT(3,1)*R(1,I) + ROT(3,2)*R(2,I) + ROT(3,3)*R(3,I)
      UP(1,I) = ROT(1,1)*U(1,I) + ROT(1,2)*U(2,I) + ROT(1,3)*U(3,I)
      UP(2,I) = ROT(2,1)*U(1,I) + ROT(2,2)*U(2,I) + ROT(2,3)*U(3,I)
      UP(3,I) = ROT(3,1)*U(1,I) + ROT(3,2)*U(2,I) + ROT(3,3)*U(3,I)
   enddo
end subroutine
