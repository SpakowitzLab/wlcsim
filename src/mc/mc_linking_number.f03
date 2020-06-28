#include "../defines.inc"
!---------------------------------------------------!
!
!   This module calculates linking number (Lk), twist (Tw)
!   , and writhe (Wr) after each move and from scratch by
!   Lk = Tw + Wr
!
!   This module does not support multi-chain simulation
!
!           Kao 08/22/2019
!
!---------------------------------------------------!

module linkingNumber
    use precision, only: dp, pi, eps
    implicit none
contains

! Calculate linking number of the whole ring from scratch
subroutine link_twist_writhe_from_scratch()
use params, only: wlc_R, wlc_U, wlc_V, wlc_LkScratch, wlc_TwScratch,&
    wlc_WrScratch
implicit none
call calc_tw(wlc_U, wlc_V, wlc_TwScratch)
call WRITHE(wlc_R, WLC_P__NB, wlc_WrScratch)
wlc_LkScratch = wlc_TwScratch + wlc_WrScratch
end subroutine

! Calculate change in twist, writhe, and linking number
subroutine get_del_tw_wr_lk(IB1, IB2, MCTYPE, delTw, delWr, delLk)
implicit none
integer, intent(in) :: IB1
integer, intent(in) :: IB2
integer, intent(in) :: MCTYPE
real(dp), intent(out) :: delTw
real(dp), intent(out) :: delWr
real(dp), intent(out) :: delLk

if (MCTYPE == 2) then
    ! mc_slide doesn't change twist
    delTw = 0.0_dp
else
    delTw = getDelTw()
endif
if (MCTYPE == 4) then
    ! mc_rotate doesn't change writhe
    delWr = 0.0_dp
else
    delWr = getDelWr(IB1, IB2, WLC_P__NB)
endif
delLk = delTw + delWr
end subroutine

! ----------------- Helper functions/subroutine -----------------!
! Calculate twist of the whole ring using triad twist approach
! This is a helper subroutine of link_twist_writhe_from_scratch
subroutine calc_tw(U, V, Tw)
implicit none
real(dp), intent(in), dimension(:, :) :: U  ! tangent vectors
real(dp), intent(in), dimension(:, :) :: V  ! binormal vectors
real(dp), intent(out) :: Tw     ! result twist (turns)
integer IT          ! bead index that goes around a ring
integer ITP1        ! index of bead next to IT

Tw = 0.0_dp
do IT = 1, WLC_P__NB
    if (IT < WLC_P__NB) then
        ITP1 = IT + 1
    else
        ITP1 = 1
    endif
    Tw = Tw + oneLinkTriadTwist(U(:, IT), V(:, IT), U(:, ITP1), V(:, ITP1))
enddo
end subroutine

! Calculate change in twist after a move using triad twist approach
! This is a helper function of delTw_Wr_Lk
function getDelTw() result(delTw)
use params, only: wlc_nBend, wlc_bendPoints, wlc_U, wlc_V, wlc_UP, wlc_VP
use ringHelper, only: nextBead
implicit none
real(dp) delTw      ! result change in twist
real(dp) oldTw      ! twist of bendPoints before the move
real(dp) newTw      ! twist of bendPoints after the move
integer IT          ! bead index
integer ITP1        ! index of the bead next to IT
integer I

oldTw = 0.0_dp
newTw = 0.0_dp
do I = 1, wlc_nBend
    IT = wlc_bendPoints(I)
    ITP1 = nextBead(IT)
    oldTw = oldTw + oneLinkTriadTwist(wlc_U(:, IT), wlc_V(:, IT), wlc_U(:, ITP1), wlc_V(:, ITP1))
    newTw = newTw + oneLinkTriadTwist(wlc_UP(:, IT), wlc_VP(:, IT), wlc_UP(:, ITP1), wlc_VP(:, ITP1))
enddo
delTw = newTw - oldTw
end function

! Calculate twist between two triads
! Helper function to calculate twist with triad twist approach
function oneLinkTriadTwist(U, V, UP1, VP1) result(Tw)
use vector_utils, only: cross
implicit none
real(dp) Tw     ! result twist (turns)
real(dp), intent(in), dimension(3) :: U     ! tangent of i-th bead
real(dp), intent(in), dimension(3) :: V     ! binormal of i-th bead
real(dp), intent(in), dimension(3) :: UP1   ! tangent of i+1-th bead
real(dp), intent(in), dimension(3) :: VP1   ! binormal of i+1-th bead

real(dp) F(3)       ! normal of i-th bead
real(dp) FP1(3)   ! normal of i+1-th bead

F = cross(V, U)
FP1 = cross(VP1, UP1)

Tw = atan2(dot_product(V, FP1) - dot_product(F, VP1), &
            dot_product(F, FP1) + dot_producT(V, VP1))
Tw = Tw / (2.0_dp*pi)
end function

! Calculate change in writhe after each move by evaluating the
! changed integrand of the double integral described in Klenin (2000)
! Method 1 before and after a move.
function getDelWr(IB1, IB2, N) result(delWr)
use ringHelper, only: nextBead, prevBead
implicit none
real(dp) delWr      ! the result total change in writhe
integer, intent(in) :: IB1      ! local index of first bead that moved
integer, intent(in) :: IB2      ! local index of last bead that moved
integer, intent(in) :: N        ! number of beads in the ring
integer NBI     ! Number of inner segments
integer NBO     ! Number of outer segments
integer I      ! bead index for iteration
integer II
integer J      ! bead index for iteration
integer JJ

! Denote a chain segment by the starting bead
! e.g. segment IB1 means R(:, IB1 + 1) - R(:, IB1)
! Call the segments that moved as inner segments
! Call the segments that don't move as outer segments
! To iterate through every changed integrand (Omega),
! we nead to pair the inner segments up among themselves, and pair
! each inner segment with all outer segments since writhe calculation
! is the double integral of Omega.
! Note that all move segments are prevBead(IB1) through IB2
delWr = 0.0_dp
! Find the length of inner and outer segments
if (IB2 >= IB1) then
    NBI = IB2 - IB1 + 2
    NBO = N - NBI
else
    NBO = IB1 - IB2 - 2
    NBI = N - NBO
endif

! Iterate through inner segments
I = prevBead(IB1)
do II = 1, NBI
    J = nextBead(I)
    do JJ = 1, NBI - II
        delWr = delWr + delWrithePair(I, J)
        J = nextBead(J)
    enddo
    I = nextBead(I)
enddo

! Iterate through outer segments
I = nextBead(IB2)
do II = 1, NBO
    J = prevBead(IB1)
    do JJ = 1, NBI
        delWr = delWr + delWrithePair(I, J)
        J = nextBead(J)
    enddo
    I = nextBead(I)
enddo

delWr = 2.0_dp * delWr
end function

! This is a helper function of getDelWr.
! Evaluate the change in writhe integrand (Omega) of a pair
! of chain segments after a move.
function delWrithePair(I1, I2) result(delOmega)
use params, only: wlc_R, wlc_RP
use ringHelper, only: nextBead
implicit none
real(dp) delOmega      ! the result change in integrand
integer, intent(in) :: I1        ! index of first segment
integer, intent(in) :: I2        ! index of second segment
real(dp) Omega      ! Integrand of a pair of segment before move
real(dp) OmegaP     ! Integrand of a pair of segment after move
integer I1P1        ! index of bead next to I1
integer I2P1        ! index of bead next to I2
real(dp) RPI1(3)    ! position of proposed position at I1
real(dp) RPI1P1(3)  ! position of proposed position at I1P1
real(dp) RPI2(3)    ! position of proposed position at I2
real(dp) RPI2P1(3)  ! position of proposed position at I2P1

I1P1 = nextBead(I1)
I2P1 = nextBead(I2)
if (I1P1 == I2 .or. I2P1 == I1) then
    ! Since Omega(i, i) = Omega(i, i+1) = 0
    delOmega = 0.0_dp
else
    Omega = dWrithe_Klenin1b(wlc_R(:, I1), wlc_R(:, I1P1), &
        wlc_R(:, I2), wlc_R(:, I2P1))
    ! If RP(:, I) is nan, that bead doesn't move; we use R(:, I) instead.
    if (isnan(wlc_RP(1, I1))) then
        RPI1 = wlc_R(:, I1)
    else
        RPI1 = wlc_RP(:, I1)
    endif
    if (isnan(wlc_RP(1, I1P1))) then
        RPI1P1 = wlc_R(:, I1P1)
    else
        RPI1P1 = wlc_RP(:, I1P1)
    endif
    if (isnan(wlc_RP(1, I2))) then
        RPI2 = wlc_R(:, I2)
    else
        RPI2 = wlc_RP(:, I2)
    endif
    if (isnan(wlc_RP(1, I2P1))) then
        RPI2P1 = wlc_R(:, I2P1)
    else
        RPI2P1 = wlc_RP(:, I2P1)
    endif
    OmegaP = dWrithe_Klenin1b(RPI1, RPI1P1, RPI2, RPI2P1)
    delOmega = OmegaP - Omega
endif
end function

! This function find the value of the integrand for each pair of
! chain segment using Method 1b described  in Klenin (2000)
! Look at Brad's util/writhe.f03 for reference
function dWrithe_Klenin1b(r1, r1P1, r2, r2P1) result(Omega)
use vector_utils, only: cross
implicit none
real(dp) Omega      ! the result integrand (solid angle)
real(dp), intent(in), dimension(3) ::  r1(3)        ! Position bead 1
real(dp), intent(in), dimension(3) ::  r1P1(3)      ! Position bead next to bead 1
real(dp), intent(in), dimension(3) ::  r2(3)        ! Position bead 2
real(dp), intent(in), dimension(3) ::  r2P1(3)      ! Position bead next to bead 2
real(dp)  r12(3)            ! Relative position vector
real(dp)  s1                ! Length of segment 1
real(dp)  s2                ! Length of segment 2
real(dp)  beta              ! Angle between tangents
real(dp)  e1(3)             ! Tangent of first segment
real(dp)  e2(3)             ! Tangent of second segment
real(dp)  cosB
real(dp)  sin2B
real(dp)  unit_cross(3)
real(dp)  a0
real(dp)  a1
real(dp)  a2
real(dp)  a3
!Variables for writhe integral
real(dp)  t1
real(dp)  t2
real(dp)  F1
real(dp)  F2
real(dp)  F3
real(dp)  F4

r12 = r2-r1
s1 = SQRT(SUM((r1P1 - r1)**2))
s2 = SQRT(SUM((r2P1 - r2)**2))
e1 = (r1P1 - r1) / s1
e2 = (r2P1 - r2) / s2

cosB = doT_PRODUCT(e1,e2)
sin2B = 1.0_dp - (cosB**2.0_dp)
!B = ACOS(cosB)
!sin2B = sin(B)**2.

unit_cross = cross(e1, e2)

if (abs(sin2B).lt.(eps))then
    sin2b = 0.0_dp
endif

a1 = doT_PRODUCT(r12,e2*cosB-e1)/sin2B
a2 = doT_PRODUCT(r12,e2-e1*cosB)/sin2B
a0 = doT_PRODUCT(r12, unit_cross)/sin2B

t1 = a1 + s1
t2 = a2 + s2
F1 = -ATAN((t1*t2 + (a0**2.)*cosB)/(a0*SQRT((t1**2. + t2**2.-2.*t1*t2*cosB + (a0**2.)*sin2B))))/(4.*PI)
t1 = a1 + s1
t2 = a2
F2 = -ATAN((t1*t2 + (a0**2.)*cosB)/(a0*SQRT((t1**2. + t2**2.-2.*t1*t2*cosB + (a0**2.)*sin2B))))/(4.*PI)
t1 = a1
t2 = a2 + s2
F3 = -ATAN((t1*t2 + (a0**2.)*cosB)/(a0*SQRT((t1**2. + t2**2.-2.*t1*t2*cosB + (a0**2.)*sin2B))))/(4.*PI)
t1 = a1
t2 = a2
F4 = -ATAN((t1*t2 + (a0**2.)*cosB)/(a0*SQRT((t1**2. + t2**2.-2.*t1*t2*cosB + (a0**2.)*sin2B))))/(4.*PI)
Omega = F1 - F2 - F3 + F4

if (Omega /= Omega) then
    Omega = 0.0_dp
endif

if (abs(a0) < 10.0_dp**(-10.0_dp)) then
    Omega = 0.0_dp
endif
end function

! --------------- Test functions/subroutine ---------------!
! Calculate delLk with total chain calculation
function getDelLkTest() result(delLk)
implicit none
real(dp) delLk
real(dp) delTw
real(dp) delWr

delTw = getDelTwTest()
delWr = getDelWrTest()
delLk = delTw + delWr
end function

! Calculate delTw with calculating total twist
function getDelTwTest() result(delTw)
use params, only: wlc_U, wlc_V, wlc_UP, wlc_VP
implicit none
real(dp) oldTw
real(dp) newTw
real(dp) delTw

call calc_tw(wlc_U, wlc_V, oldTw)
call calc_tw(wlc_UP, wlc_VP, newTw)
delTw = newTw - oldTw
end function

! Calculate change in writhe after each move by evaluating Wr before
! and after the move. Note that to use this function, all wlc_R
! must be copied to wlc_RP so that there is no nan in the unmoved RP
function getDelWrTest() result(delWr)
use params, only: wlc_R, wlc_RP
implicit none
real(dp) oldWr
real(dp) newWr
real(dp) delWr

call WRITHE(wlc_R, WLC_P__NB, oldWr)
call WRITHE(wlc_RP, WLC_P__NB, newWr)
delWr = newWr - oldWr
end  function

end module linkingNumber
