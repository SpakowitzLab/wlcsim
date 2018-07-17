#include "../defines.inc"
!--------------------------------------------------------------*
!
!           Makes Monti Carlo Moves
!
!    Quinn separated out this file on 8/9/17
!
!---------------------------------------------------------------

! variables that need to be allocated only on certain branches moved into MD to prevent segfaults
! please move other variables in as you see fit
subroutine MC_chemMove(wlc_p,wlc_d,R,U,RP,UP,AB,ABP,IB1,IB2,IT1,IT2 &
                  ,WindoW,rand_stat,success)

use mersenne_twister
use params, only: dp,wlcsim_params,wlcsim_data
use windowTools, only: drawWindow

implicit none
type(wlcsim_params), intent(in) :: wlc_p
type(wlcsim_data), intent(in) :: wlc_d
real(dp), intent(in) :: R(3,WLC_P__NT)  ! Bead positions
real(dp), intent(in) :: U(3,WLC_P__NT)  ! Tangent vectors
integer, intent(in) :: AB(WLC_P__NT)  ! Tangent vectors
real(dp), intent(out) :: RP(3,WLC_P__NT)  ! Bead positions
real(dp), intent(out) :: UP(3,WLC_P__NT)  ! Tangent vectors
integer, intent(out) :: ABP(WLC_P__NT)  ! Tangent vectors
!integer, intent(in) :: WLC_P__NBPM    ! Beads per monomer, aka G
integer, intent(out) :: IB1   ! Test bead position 1
integer, intent(out) :: IT1   ! Index of test bead 1
integer, intent(out) :: IB2   ! Test bead position 2
integer, intent(out) :: IT2   ! Index of test bead 2
logical, intent(out) ::success
integer dib   ! number of beads moved by move

integer IP    ! Test polymer
integer I,J  ! Test indices
! Things for random number generator
type(random_stat), intent(inout) :: rand_stat  ! status of random number generator
real urnd(1) ! single random number
real(dp), intent(in) :: WindoW ! Size of window for bead selection
integer TEMP

integer, parameter, dimension(0:3) :: changeBoth = [3, 2, 1, 0]
integer, parameter, dimension(0:3) :: changeFirst = [2, 3, 0, 1]
integer, parameter, dimension(0:3) :: changeSecond = [1,0, 3, 2]

!TOdo saving RP is not actually needed, even in these cases, but Brad's code assumes that we have RP.
if (WLC_P__RING .OR. WLC_P__INTERP_BEAD_LENNARD_JONES) then
    RP = R
    UP = U
endif

! Change wlc_d%AB (a.k.a HP1 binding type fore section of polymer)

call drawWindow(wlc_d,window,WLC_P__MAXWINDOW_SLIDE_MOVE,.false.,rand_stat,&
                IT1,IT2,IB1,IB2,IP,DIB,success)
if (success .eqv.  .false.) return

!keep binding constant within monomers
IT1 = IT1-MOD(IT1-1,WLC_P__NBPM)
IT2 = IT2-MOD(IT2-1,WLC_P__NBPM) + WLC_P__NBPM-1
IB1 = MOD(IT1,WLC_P__NB)
IB2 = MOD(IT2,WLC_P__NB)


if (WLC_P__TWO_TAIL) then
    call random_number(urnd,rand_stat)
    if (urnd(1)>WLC_P__PROBSINGLESWAP) then
        do J = IT1, IT2
            ABP(J) = changeBoth(AB(J))
        enddo
    elseif (urnd(1)<WLC_P__PROBSINGLESWAP/2.0_dp) then
        do J = IT1, IT2
            ABP(J) = changeFirst(AB(J))
        enddo
    else
        do J = IT1, IT2
            ABP(J) = changeSecond(AB(J))
        enddo
    endif
else
    do J = IT1,IT2
        ABP(J) = 1-AB(J)
    ENDdo
endif



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
