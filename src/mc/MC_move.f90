!--------------------------------------------------------------*
!
!           Makes Monti Carlo Moves
!
!    Quinn Made Changes to this file starting on 12/15/15
!
!---------------------------------------------------------------

! variables that need to be allocated only on certain branches moved into MD to prevent segfaults
! please move other variables in as you see fit
subroutine MC_move(wlc_p,wlc_d,IB1,IB2,IT1,IT2,IT3,IT4,IP,MCTYPE,forward,rand_stat,dib)
use mersenne_twister
use params, only: dp, pi, wlcsim_data, wlcsim_params
implicit none

integer, intent(out) :: IT1, IT2, IT3, IT4, IP, IB1, IB2, dib
integer, intent(in) :: MCTYPE
logical, intent(out) :: forward
type(random_stat), intent(inout) :: rand_stat  ! status of random number generator


type(wlcsim_params), intent(in) :: wlc_p
type(wlcsim_data), intent(inout) :: wlc_d
!     Perform crank-shaft move (MCTYPE 1)

select case(MCTYPE) ! pick which keyword, case matchign string must be all uppercase
case(1) 
call MC_crank(wlc_d%R,wlc_d%U,wlc_d%RP,wlc_d%UP,wlc_p%NT,wlc_p%NB,wlc_p%NP,&
        IP,IB1,IB2,IT1,IT2,MCTYPE &
       ,wlc_d%MCAMP,wlc_d%Window,rand_stat,wlc_p%winType &
       ,dib,wlc_p%ring,wlc_p%inTERP_BEAD_LENNARD_JONES)
case(2)
call MC_slide(wlc_d%R,wlc_d%U,wlc_d%RP,wlc_d%UP,wlc_p%NT,wlc_p%NB,wlc_p%NP,&
        IP,IB1,IB2,IT1,IT2,MCTYPE &
       ,wlc_d%MCAMP,wlc_d%Window,rand_stat,wlc_p%winType &
       ,dib,wlc_p%ring,wlc_p%inTERP_BEAD_LENNARD_JONES)
case(3)
call MC_pivot(wlc_d%R,wlc_d%U,wlc_d%RP,wlc_d%UP,wlc_p%NT,wlc_p%NB,wlc_p%NP,&
        IP,IB1,IB2,IT1,IT2,MCTYPE &
       ,wlc_d%MCAMP,wlc_d%Window,rand_stat,wlc_p%winType &
       ,wlc_p%ring,wlc_p%inTERP_BEAD_LENNARD_JONES)
case(4)
call MC_rotate(wlc_d%R,wlc_d%U,wlc_d%RP,wlc_d%UP,wlc_p%NT,wlc_p%NB,wlc_p%NP,&
        IP,IB1,IB2,IT1,IT2 &
       ,wlc_d%MCAMP,rand_stat &
       ,wlc_p%ring,wlc_p%inTERP_BEAD_LENNARD_JONES)
case(5)
call MC_fullChainRotation(wlc_d%R,wlc_d%U,wlc_d%RP,wlc_d%UP,wlc_p%NT,wlc_p%NB,&
        wlc_p%NP,IP,IB1,IB2,IT1,IT2 &
       ,wlc_d%MCAMP,rand_stat &
       ,wlc_p%ring,wlc_p%inTERP_BEAD_LENNARD_JONES)
case(6)
call MC_fullChainSlide(wlc_d%R,wlc_d%U,wlc_d%RP,wlc_d%UP,wlc_p%NT,wlc_p%NB,&
        wlc_p%NP,IP,IB1,IB2,IT1,IT2 &
       ,wlc_d%MCAMP,rand_stat &
       ,wlc_p%ring,wlc_p%inTERP_BEAD_LENNARD_JONES)
case(7)
call MC_chemMove(wlc_d%R,wlc_d%U,wlc_d%RP,wlc_d%UP,wlc_d%AB,wlc_d%ABP,wlc_p%NT,&
         wlc_p%NB,wlc_p%NP,IP,IB1,IB2,IT1,IT2,MCTYPE &
       ,wlc_d%Window,wlc_p%nBPM,rand_stat &
       ,wlc_p%ring,wlc_p%inTERP_BEAD_LENNARD_JONES)
case(8)
case(9)
call MC_chainSwap(wlc_d%R,wlc_d%U,wlc_d%RP,wlc_d%UP,wlc_p%NT,wlc_p%NB,&
        wlc_p%NP,IP,IB1,IB2,IT1,IT2 &
       ,rand_stat &
       ,IT3,IT4,wlc_p%ring,wlc_p%inTERP_BEAD_LENNARD_JONES)
case(10)
call MC_reptation(wlc_d%R,wlc_d%U,wlc_d%RP,wlc_d%UP,wlc_p%NT,wlc_p%NB,&
        wlc_p%NP,IP,IT1,IT2&
       ,rand_stat &
       ,forward,wlc_p%ring,wlc_p%inTERP_BEAD_LENNARD_JONES)
end select 
RETURN
END
subroutine test_equiv_forward(U,R,UP,RP,NT,IT1,IT2,RparaMag,RperpMag)
use params, only: dp
implicit none
! inputs
integer NT,IT1,IT2
real(dp) R(3,NT)  ! Bead positions
real(dp) U(3,NT)  ! Tangent vectors
real(dp) RP(3,NT)  ! Bead positions
real(dp) UP(3,NT)  ! Tangent vectors
real(dp) RparaMag, RperpMag

!defined
real(dp) drOld(3)
real(dp) drNew(3)
real(dp) drParOld, drParNew
real(dp) drPerpOld(3)
real(dp) drPerpNew(3)
real(dp) Eta
real(dp) GIOld(3)
real(dp) Ginew(3)
Eta = 1.89756278_dp

drOld(1) = R(1,IT1 + 1)-R(1,IT1)
drOld(2) = R(2,IT1 + 1)-R(2,IT1)
drOld(3) = R(3,IT1 + 1)-R(3,IT1)
DRPAROld = DROld(1)*U(1,IT1) + DROld(2)*U(2,IT1) + DROld(3)*U(3,IT1)
drNew(1) = RP(1,IT2)-RP(1,IT2-1)
drNew(2) = RP(2,IT2)-RP(2,IT2-1)
drNew(3) = RP(3,IT2)-RP(3,IT2-1)
DRPARNew = DRNew(1)*UP(1,IT2-1) + &
         DRNew(2)*UP(2,IT2-1) + &
         DRNew(3)*UP(3,IT2-1)
if (abs(drOld(1)**2 + drOld(2)**2 + drOld(3)**2&
      -(drNew(1)**2 + drNew(2)**2 + drNew(3)**2)).gt.0.000001) then
      print*, "drOld",drOld, " mag^2 = ",drOld(1)**2 + drOld(2)**2 + drOld(3)**2
      print*, "drNew",drNew, " mag^2 = ",drNew(1)**2 + drNew(2)**2 + drNew(3)**2
      print*, "Difference detected in test_equiv, 0"
      stop 1
endif

if (abs(drParOld-drParNew).gt.0.0000001_dp) then
    print*, "DRParOld",DRParOld,"DRParNew",DRParNew
    print*, "Difference detected in test_equiv, 1"
    stop 1
endif

drPerpOld(1) = drOld(1)-drParOld*U(1,IT1)
drPerpOld(2) = drOld(2)-drParOld*U(2,IT1)
drPerpOld(3) = drOld(3)-drParOld*U(3,IT1)
drPerpNew(1) = drNew(1)-drParNew*UP(1,IT2-1)
drPerpNew(2) = drNew(2)-drParNew*UP(2,IT2-1)
drPerpNew(3) = drNew(3)-drParNew*UP(3,IT2-1)

if (abs(drPerpOld(1)**2 + drPerpOld(2)**2 + drPerpOld(3)**2 &
      -(drPerpNew(1)**2 + drPerpNew(2)**2 + drPerpNew(3)**2)).gt.0.000001_dp) then
  print*, "drOld",sqrt(drOld(1)**2 + drOld(2)**2 + drOld(3)**2)
  print*, "drNew",sqrt(drNew(1)**2 + drNew(2)**2 + drNew(3)**2)
  print*, "dRparOld",dRparOld,"dRparNew",drParNew
  print*, "perp Old:", drPerpOld(1)**2 + drPerpOld(2)**2 + drPerpOld(3)**2
  print*, "perp New:", drPerpNew(1)**2 + drPerpNew(2)**2 + drPerpNew(3)**2
  print*, "RparaMag",RparaMag,"RperpMag",RperpMag
  print*, "Difference detected in test_equiv, 2"
  stop 1
endif

GIOld(1) = U(1,IT1 + 1)-U(1,IT1)-Eta*dRperpOld(1)
GIOld(2) = U(2,IT1 + 1)-U(2,IT1)-Eta*dRperpOld(2)
GIOld(3) = U(3,IT1 + 1)-U(3,IT1)-Eta*dRperpOld(3)
Ginew(1) = UP(1,IT2)-UP(1,IT2-1)-Eta*dRperpNew(1)
Ginew(2) = UP(2,IT2)-UP(2,IT2-1)-Eta*dRperpNew(2)
Ginew(3) = UP(3,IT2)-UP(3,IT2-1)-Eta*dRperpNew(3)

if (abs(GIOld(1)**2 + GIOld(2)**2 + GIOld(3)**2&
      -(Ginew(1)**2 + Ginew(2)**2 + Ginew(3)**2)).gt.0.000001_dp) then
  print*, "Difference detected in test_equiv, 3"
  print*, "GIOld(1)**2 + GIOld(2)**2 + GIOld(3)**2", &
           GIOld(1)**2 + GIOld(2)**2 + GIOld(3)**2
  print*, "Ginew(1)**2 + Ginew(2)**2 + Ginew(3)**2", &
          Ginew(1)**2 + Ginew(2)**2 + Ginew(3)**2
  print*, "RparaMag",RparaMag,"RperpMag",RperpMag
  stop 1
endif

return
end subroutine
subroutine random_perp(u,p,t,rand_stat)
! The subroutine generates the second two vectors in a unit triad
! The output vectors, p and t, are perpendicular to eachother and u
! The triad is randomly left or right handed
use mersenne_twister
use params, only: dp
implicit none
real(dp), PARAMETER :: PI = 3.141592654 ! Value of pi
type(random_stat) rand_stat  ! status of random number generator
real urnd(1) ! single random number

real(dp) v(2) ! random 2-vec
real(dp), intent(in) :: u(3) ! input
real(dp), intent(out) :: p(3) ! output: random perpendicular to u
real(dp), intent(out) :: t(3) ! orthogonal to p and u
real(dp) f

if (abs(u(1)**2 + u(2)**2 + u(3)**2-1.0_dp) .gt. 0.0000001_dp) then
    print*, u
    print*, "Error in random_perp, please give me a unit vector"
    stop 1
endif

call random_number(urnd,rand_stat)
v(1) = cos(2*PI*urnd(1))
v(2) = sin(2*PI*urnd(1))

if (u(3).gt.0.0) then
    f = 1.0_dp/(1 + u(3))
    p(1) = (u(3) + f*u(2)**2)*v(1) - u(2)*u(1)*v(2)*f
    p(2) = (u(3) + f*u(1)**2)*v(2) - u(2)*u(1)*v(1)*f
    p(3) = -1.0_dp*(u(2)*v(2) + u(1)*v(1))
else
    f = 1.0_dp/(1-u(3))
    p(1) = (-u(3) + f*u(2)**2)*v(1) - u(2)*u(1)*v(2)*f
    p(2) = (-u(3) + f*u(1)**2)*v(2) - u(2)*u(1)*v(1)*f
    p(3) = (u(2)*v(2) + u(1)*v(1))

endif

t(1) = u(2)*p(3)-u(3)*p(2)
t(2) = u(3)*p(1)-u(1)*p(3)
t(3) = u(1)*p(2)-u(2)*p(1)

! random sign
call random_number(urnd,rand_stat)
if (urnd(1).lt.0.5_dp) then
    t(1) = -1.0_dp*t(1)
    t(2) = -1.0_dp*t(2)
    t(3) = -1.0_dp*t(3)
endif

! Testing
!if (abs(p(1)*u(1) + p(2)*u(2) + p(3)*u(3)).gt.0.000001_dp) then
!    print*, "Error in random_perp, 1"
!    stop 1
!endif
!if (abs(p(1)**2 + p(2)**2 + p(3)**2-1) .gt. 0.0000001_dp) then
!    print*, "Error in random_perp, 2"
!    stop 1
!endif
!if (abs(t(1)**2 + t(2)**2 + t(3)**2 -1).gt.0.000001_dp) then
!    print*, "Error in random_perp, 3"
!    stop 1
!endif
!if (abs(t(1)*p(1) + t(2)*p(2) + t(3)*p(3)).gt.0.0000001_dp) then
!    print*, "Error in random_perp, 4"
!    stop 1
!endif
!if (abs(t(1)*u(1) + t(2)*u(2) + t(3)*u(3)).gt.0.0000001_dp) then
!    print*, "Error in random_perp, 5"
!    stop 1
!endif
! END Testing

return
end subroutine
!---------------------------------------------------------------!
