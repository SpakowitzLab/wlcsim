#include "../defines.inc"
!--------------------------------------------------------------*
!
!           Makes Monti Carlo Moves
!
!    Quinn split out this file on 8/9/17
!
!---------------------------------------------------------------

! variables that need to be allocated only on certain branches moved into MD to prevent segfaults
! please move other variables in as you see fit
subroutine MC_crank(wlc_p,wlc_d,R,U,RP,UP,IB1,IB2,IT1,IT2 &
                  ,MCAMP,WindoW,rand_stat  &
                  ,dib,success)

use mersenne_twister
use params, only: dp,wlcsim_params, wlcsim_data
use vector_utils, only: rotateR, rotateU, axisAngle
use windowTools, only: drawWindow

implicit none
type(wlcsim_params),intent(in) :: wlc_p
type(wlcsim_data), intent(inout) :: wlc_d
!integer, intent(in) :: ExplicitBindingPair(WLC_P__NT)
real(dp), intent(in) :: R(3,WLC_P__NT)  ! Bead positions
real(dp), intent(in) :: U(3,WLC_P__NT)  ! Tangent vectors
real(dp), intent(out) :: RP(3,WLC_P__NT)  ! Bead positions
real(dp), intent(out) :: UP(3,WLC_P__NT)  ! Tangent vectors
integer, intent(out) :: IB1   ! Test bead position 1
integer, intent(out) :: IT1   ! Index of test bead 1
integer, intent(out) :: IB2   ! Test bead position 2
integer, intent(out) :: IT2   ! Index of test bead 2
integer, intent(out) :: dib   ! number of beads moved by move
logical, intent(inout) :: success

integer IP    ! Test polymer
integer I,J  ! Test indices
! Things for random number generator
type(random_stat), intent(inout) :: rand_stat  ! status of random number generator
real urnd(1) ! single random number
! Variables for the crank-shaft move

real(dp) TA(3)    ! Axis of rotation
real(dp) P1(3)    ! Point on rotation line
real(dp) ROT(3,4) ! Rotation matrix

real(dp) ALPHA    ! Angle of move

!     MC adaptation variables

real(dp), intent(in) :: MCAMP ! Amplitude of random change
!integer, intent(in) :: winType
real(dp), intent(in) :: WindoW ! Size of window for bead selection

! Variables for change of binding state move
real(dp) d1,d2  !for testing

!TOdo saving RP is not actually needed, even in these cases, but Brad's code assumes that we have RP.
if (WLC_P__RING .OR.WLC_P__INTERP_BEAD_LENNARD_JONES) then
    RP = R
    UP = U
    P1 = 0.0_dp
endif

!     Perform crank-shaft move (MCTYPE 1)
call drawWindow(wlc_d,window,WLC_P__MAXWINDOW_CRANK_SHAFT,.true.,rand_stat,&
                IT1,IT2,IB1,IB2,IP,DIB,success)
if (success .eqv. .false.) return


if (WLC_P__RING) then                    !Polymer is a ring
   if (IB1 == IB2.AND.IB1 == 1) then
      TA(1) = R(1,IT1 + 1)-R(1,WLC_P__NB*IP)
      TA(2) = R(2,IT1 + 1)-R(2,WLC_P__NB*IP)
      TA(3) = R(3,IT1 + 1)-R(3,WLC_P__NB*IP)
   elseif (IB1 == IB2.AND.IB1 == WLC_P__NB) then
      TA(1) = R(1,WLC_P__NB*(IP-1) + 1)-R(1,IT1-1)
      TA(2) = R(2,WLC_P__NB*(IP-1) + 1)-R(2,IT1-1)
      TA(3) = R(3,WLC_P__NB*(IP-1) + 1)-R(3,IT1-1)
   elseif (IB1 == IB2.AND.IB1 /= 1.AND.IB2 /=WLC_P__NB) then
      TA(1) = R(1,IT1 + 1)-R(1,IT1-1)
      TA(2) = R(2,IT1 + 1)-R(2,IT1-1)
      TA(3) = R(3,IT1 + 1)-R(3,IT1-1)
   else
      TA(1) = R(1,IT2)-R(1,IT1)
      TA(2) = R(2,IT2)-R(2,IT1)
      TA(3) = R(3,IT2)-R(3,IT1)
   endif
else                                 !Polymer is not a ring
  if (IB1 == IB2.AND.IB1 == 1) then
      TA(1) = R(1,IT1 + 1)-R(1,IT1)
      TA(2) = R(2,IT1 + 1)-R(2,IT1)
      TA(3) = R(3,IT1 + 1)-R(3,IT1)
   elseif (IB1 == IB2.AND.IB1 == WLC_P__NB) then
      TA(1) = R(1,WLC_P__NB*IP)-R(1,WLC_P__NB*IP-1)
      TA(2) = R(2,WLC_P__NB*IP)-R(2,WLC_P__NB*IP-1)
      TA(3) = R(3,WLC_P__NB*IP)-R(3,WLC_P__NB*IP-1)
   elseif (IB1 == IB2.AND.IB1 /= 1.AND.IB2 /= WLC_P__NB) then
      TA(1) = R(1,IT1 + 1)-R(1,IT1-1)
      TA(2) = R(2,IT1 + 1)-R(2,IT1-1)
      TA(3) = R(3,IT1 + 1)-R(3,IT1-1)
   else
      TA(1) = R(1,IT2)-R(1,IT1)
      TA(2) = R(2,IT2)-R(2,IT1)
      TA(3) = R(3,IT2)-R(3,IT1)
   endif
endif

  P1(1) = R(1,IT1)
  P1(2) = R(2,IT1)
  P1(3) = R(3,IT1)
  call random_number(urnd,rand_stat)
  ALPHA = MCAMP*(urnd(1)-0.5_dp)

  call axisAngle(ROT,alpha,TA,P1)

  I = IT1
   do J = 0,DIB
     if (I == (WLC_P__NB*IP+1).and.WLC_P__RING) then
          I = WLC_P__NB*(IP-1)+1
     endif
     RP(:,I) = rotateR(ROT,R(:,I))
     UP(:,I) = rotateU(ROT,U(:,I))
     I = I + 1
  ENDdo

!  ------begining testing---------
if(.false.) then
    ! This is a code block for testing
    if (abs(RP(1,IT1)-R(1,IT1)).gt.0.000001) then
        print*, "error in crank-shaft move"
        print*, RP(1,IT1), R(1,IT1)
        stop 1
    endif
    if (abs(RP(1,IT2)-R(1,IT2)).gt.0.000001) then
        print*, "error in crank-shaft move"
        print*, RP(1,IT1), R(1,IT1)
        stop 1
    endif
    if(IT1.ne.IT2) then
        d1 = (R(1,IT1 + 1)-R(1,IT1))**2 + &
           (R(2,IT1 + 1)-R(2,IT1))**2 + &
           (R(3,IT1 + 1)-R(3,IT1))**2
        d2 = (RP(1,IT1 + 1)-RP(1,IT1))**2 + &
           (RP(2,IT1 + 1)-RP(2,IT1))**2 + &
           (RP(3,IT1 + 1)-RP(3,IT1))**2
        if (abs(d1-d2).gt.0.000001) then
            print*, "error in crank-shaft move"
            print*, "distance change in 1"
            print*, "IT1",IT1," IT2",IT2
            print*, d1,d2
            stop 1
        endif
        d1 = (R(1,IT2-1)-R(1,IT2))**2 + &
           (R(2,IT2-1)-R(2,IT2))**2 + &
           (R(3,IT2-1)-R(3,IT2))**2
        d2 = (RP(1,IT2-1)-RP(1,IT2))**2 + &
           (RP(2,IT2-1)-RP(2,IT2))**2 + &
           (RP(3,IT2-1)-RP(3,IT2))**2
        if (abs(d1-d2).gt.0.000001) then
            print*, "error in crank-shaft move"
            print*, "distance change in 2"
            print*, d1,d2
            stop 1
        endif
    endif
endif
! --------end testing--------
end subroutine
