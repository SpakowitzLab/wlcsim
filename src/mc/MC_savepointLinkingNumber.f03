#include "../defines.inc"
!----------------------------------------------------!
!
!   This module calculate linking number
!   at each save point, including at the beginning
!   of the simulation
!
!   This module exists to avoid cyclic requirement
!   of use statement in MC_linkingNumber
!
!   This module doesn't support multi-chain simulation
!   Kao 08/20/2019
!
!-----------------------------------------------------!

module savepointLinkingNumber
     use precision, only: dp, pi
     implicit none
contains

! Calculate linking number of the whole ring by Tw + Wr
subroutine calcTwWrLk(R, U, V, Tw, Wr, Lk)
implicit none
real(dp), intent(in), dimension(:, :) :: R  ! positions
real(dp), intent(in), dimension(:, :) :: U  ! tangent vectors
real(dp), intent(in), dimension(:, :) :: V  ! binormal vectors
real(dp), intent(out) :: Tw     ! result twist
real(dp), intent(out) :: Wr     ! result writhe
real(dp), intent(out) :: Lk     ! result linking number

call calcTw(U, V, Tw) 
call WRITHE(R, WLC_P__NB, Wr)
Lk = Tw + Wr

end subroutine

! Calculate twist of the whole ring using triad twist approach
subroutine calcTw(U, V, Tw)
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

! ------------------- Helper functions -------------------------!

! Twist between two triads
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

end module savepointLinkingNumber
