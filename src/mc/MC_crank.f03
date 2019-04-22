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
subroutine MC_crank(IB1,IB2,IT1,IT2,MCAMP,WindoW,rand_stat,dib,success)
! values from wlcsim_data
use params, only: wlc_VP, wlc_UP, wlc_V, wlc_RP, wlc_U, wlc_R, wlc_nBend, wlc_bendPoints, wlc_nPointsMoved, wlc_pointsMoved

use mersenne_twister
use params, only: dp,  eps
use vector_utils, only: rotateR, rotateU, axisAngle, randomUnitVec
use windowTools, only: drawWindow
use polydispersity, only: length_of_chain, first_bead_of_chain, last_bead_of_chain

implicit none
integer, intent(out) :: IB1   ! Test bead position 1
integer, intent(out) :: IT1   ! Index of test bead 1
integer, intent(out) :: IB2   ! Test bead position 2
integer, intent(out) :: IT2   ! Index of test bead 2
integer, intent(out) :: dib   ! number of beads moved by move
logical, intent(out) :: success

integer IP    ! Test polymer
integer I,J  ! Test indices
! Things for random number generator
type(random_stat), intent(inout) :: rand_stat  ! status of random number generator
real(dp) urnd(1) ! single random number
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
    wlc_RP = wlc_R
    wlc_UP = wlc_U
    P1 = 0.0_dp
endif

!     Perform crank-shaft move (MCTYPE 1)
call drawWindow(window,WLC_P__MAXWINDOW_CRANK_SHAFT,.true.,rand_stat,&
                IT1,IT2,IB1,IB2,IP,DIB,success)
if (success .eqv. .false.) return

!  Which elastic segments change
wlc_nBend = 0
if (IB1>1) then
    wlc_nBend = wlc_nBend + 1
    wlc_bendPoints(wlc_nBend)=IT1-1
    I=IT1-1
    wlc_RP(:,I)=wlc_R(:,I)
    wlc_UP(:,I)=wlc_U(:,I)
    if (WLC_P__LOCAL_TWIST) wlc_VP(:,I) = wlc_V(:,I)
    wlc_nPointsMoved=wlc_nPointsMoved+1
    wlc_pointsMoved(wlc_nPointsMoved)=I
endif
if (IB2<length_of_chain(IP)) then
    wlc_nBend = wlc_nBend + 1
    wlc_bendPoints(wlc_nBend)=IT2
    I=IT2+1
    wlc_RP(:,I)=wlc_R(:,I)
    wlc_UP(:,I)=wlc_U(:,I)
    if (WLC_P__LOCAL_TWIST) wlc_VP(:,I) = wlc_V(:,I)
    wlc_nPointsMoved=wlc_nPointsMoved+1
    wlc_pointsMoved(wlc_nPointsMoved)=I
endif

if (WLC_P__RING) then                    !Polymer is a ring
   if (IB1 == IB2.AND.IB1 == 1) then
      TA = wlc_R(:,IT1 + 1)-wlc_R(:,last_bead_of_chain(IP))
   elseif (IB1 == IB2.AND.IT1 == last_bead_of_chain(IP)) then
      TA = wlc_R(:,first_bead_of_chain(IP))-wlc_R(:,IT1-1)
   elseif (IB1 == IB2.AND.IB1 /= 1.AND.IB2 /= length_of_chain(IP)) then
      TA = wlc_R(:,IT1 + 1)-wlc_R(:,IT1-1)
   else
      TA = wlc_R(:,IT2)-wlc_R(:,IT1)
   endif
else                                 !Polymer is not a ring
  if (IB1 == IB2.AND.IB1 == 1) then
      TA = wlc_R(:,IT1 + 1)-wlc_R(:,IT1)
   elseif (IB1 == IB2.AND.IT1 == last_bead_of_chain(IP)) then
      TA = wlc_R(:,IT1)-wlc_R(:,IT1-1)
   elseif (IB1 == IB2.AND.IB1 /= 1.AND.IB2 /= length_of_chain(IP)) then
      TA = wlc_R(:,IT1 + 1)-wlc_R(:,IT1-1)
   else
      TA = wlc_R(:,IT2)-wlc_R(:,IT1)
   endif
endif
  if (norm2(TA)<eps) then
      call randomUnitVec(TA,rand_stat)
  endif

  P1 = wlc_R(:,IT1)
  call random_number(urnd,rand_stat)
  ALPHA = MCAMP*(urnd(1)-0.5_dp)

  call axisAngle(ROT,alpha,TA,P1)

  I = IT1
   do J = 0,DIB
     if (I > last_bead_of_chain(IP).and.WLC_P__RING) then
          I = first_bead_of_chain(IP)
     endif
     wlc_RP(:,I) = rotateR(ROT,wlc_R(:,I))
     wlc_UP(:,I) = rotateU(ROT,wlc_U(:,I))
     if (WLC_P__LOCAL_TWIST) wlc_VP(:,I) = rotateU(ROT,wlc_V(:,I))
     wlc_nPointsMoved=wlc_nPointsMoved+1
     wlc_pointsMoved(wlc_nPointsMoved)=I
     I = I + 1
  ENDdo

!  ------begining testing---------
if(.false.) then
    ! This is a code block for testing
    if (abs(wlc_RP(1,IT1)-wlc_R(1,IT1)).gt.0.000001_dp) then
        print*, "error in crank-shaft move"
        print*, wlc_RP(1,IT1), wlc_R(1,IT1)
        stop 1
    endif
    if (abs(wlc_RP(1,IT2)-wlc_R(1,IT2)).gt.0.000001_dp) then
        print*, "error in crank-shaft move"
        print*, wlc_RP(1,IT1), wlc_R(1,IT1)
        stop 1
    endif
    if(IT1.ne.IT2) then
        d1 = (wlc_R(1,IT1 + 1)-wlc_R(1,IT1))**2 + &
           (wlc_R(2,IT1 + 1)-wlc_R(2,IT1))**2 + &
           (wlc_R(3,IT1 + 1)-wlc_R(3,IT1))**2
        d2 = (wlc_RP(1,IT1 + 1)-wlc_RP(1,IT1))**2 + &
           (wlc_RP(2,IT1 + 1)-wlc_RP(2,IT1))**2 + &
           (wlc_RP(3,IT1 + 1)-wlc_RP(3,IT1))**2
        if (abs(d1-d2).gt.0.000001_dp) then
            print*, "error in crank-shaft move"
            print*, "distance change in 1"
            print*, "IT1",IT1," IT2",IT2
            print*, d1,d2
            stop 1
        endif
        d1 = (wlc_R(1,IT2-1)-wlc_R(1,IT2))**2 + &
           (wlc_R(2,IT2-1)-wlc_R(2,IT2))**2 + &
           (wlc_R(3,IT2-1)-wlc_R(3,IT2))**2
        d2 = (wlc_RP(1,IT2-1)-wlc_RP(1,IT2))**2 + &
           (wlc_RP(2,IT2-1)-wlc_RP(2,IT2))**2 + &
           (wlc_RP(3,IT2-1)-wlc_RP(3,IT2))**2
        if (abs(d1-d2).gt.0.000001_dp) then
            print*, "error in crank-shaft move"
            print*, "distance change in 2"
            print*, d1,d2
            stop 1
        endif
    endif
endif
! --------end testing--------
end subroutine
