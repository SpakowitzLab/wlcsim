!--------------------------------------------------------------*
!
!           Makes Monti Carlo Moves
!
!    Quinn Made Changes to this file starting on 12/15/15
!
!---------------------------------------------------------------

! variables that need to be allocated only on certain branches moved into MD to prevent segfaults
! please move other variables in as you see fit
subroutine MC_move(wlc_p,wlc_d,IB1,IB2,IT1,IT2,IT3,IT4,MCTYPE,forward,rand_stat,dib,spider_id,success)
use mersenne_twister
use params, only: dp, pi, wlcsim_data, wlcsim_params
implicit none

integer, intent(out) :: IT1, IT2, IT3, IT4, IB1, IB2, dib,spider_id
logical, intent(out) :: success
integer, intent(in) :: MCTYPE
logical, intent(out) :: forward
type(random_stat), intent(inout) :: rand_stat  ! status of random number generator


type(wlcsim_params), intent(in) :: wlc_p
type(wlcsim_data), intent(inout) :: wlc_d

success = .TRUE.
spider_id = 0
select case(MCTYPE) ! pick which keyword, case matchign string must be all uppercase
case(1)
call MC_crank(wlc_p,wlc_d,wlc_d%R,wlc_d%U,wlc_d%RP,wlc_d%UP&
       ,IB1,IB2,IT1,IT2 &
       ,wlc_d%MCAMP(MCTYPE),wlc_d%Window(MCTYPE),rand_stat &
       ,dib,success)
case(2)
call MC_slide(wlc_p,wlc_d,wlc_d%R,wlc_d%U,wlc_d%RP,wlc_d%UP&
       ,IB1,IB2,IT1,IT2 &
       ,wlc_d%MCAMP(MCTYPE),wlc_d%Window(MCTYPE),rand_stat &
       ,dib,success)
case(3)
call MC_pivot(wlc_p,wlc_d,wlc_d%R,wlc_d%U,wlc_d%RP,wlc_d%UP&
       ,IB1,IB2,IT1,IT2 &
       ,wlc_d%MCAMP(MCTYPE),wlc_d%Window(MCTYPE),rand_stat,success)
case(4)
call MC_rotate(wlc_p,wlc_d%R,wlc_d%U,wlc_d%RP,wlc_d%UP&
       ,IB1,IB2,IT1,IT2 &
       ,wlc_d%MCAMP(MCTYPE),rand_stat)
case(5)
call MC_fullChainRotation(wlc_p,wlc_d%R,wlc_d%U,wlc_d%RP,wlc_d%UP&
       ,IB1,IB2,IT1,IT2 &
       ,wlc_d%MCAMP(MCTYPE),rand_stat)
case(6)
call MC_fullChainSlide(wlc_p,wlc_d%R,wlc_d%U,wlc_d%RP,wlc_d%UP&
       ,IB1,IB2,IT1,IT2 &
       ,wlc_d%MCAMP(MCTYPE),rand_stat)
case(7)
call MC_chemMove(wlc_p,wlc_d%R,wlc_d%U,wlc_d%RP,wlc_d%UP,wlc_d%AB,wlc_d%ABP,IB1,IB2,IT1,IT2 &
       ,wlc_d%Window(MCTYPE),rand_stat)
case(8)
case(9)
call MC_chainSwap(wlc_p,wlc_d%R,wlc_d%U,wlc_d%RP,wlc_d%UP,IB1,IB2,IT1,IT2 &
       ,rand_stat &
       ,IT3,IT4)
case(10)
call MC_reptation(wlc_p,wlc_d%R,wlc_d%U,wlc_d%RP,wlc_d%UP,IT1,IT2,IB1,IB2&
       ,rand_stat &
       ,forward)
case(11)
call MC_superReptation(wlc_p,wlc_d%R,wlc_d%U,wlc_d%RP,wlc_d%UP,wlc_d%AB,wlc_d%ABP&
        ,IT1,IT2,IB1,IB2,rand_stat &
       ,forward)
case(12)
call MC_spider(wlc_d,wlc_d%MCAMP,rand_stat,success,spider_id)
end select
RETURN
END
subroutine test_equiv_forward(U,R,UP,RP,NT,IT1,IT2,RparaMag,RperpMag)
use params, only: dp, eps
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
      -(drNew(1)**2 + drNew(2)**2 + drNew(3)**2)).gt.eps) then
      print*, "drOld",drOld, " mag^2 = ",drOld(1)**2 + drOld(2)**2 + drOld(3)**2
      print*, "drNew",drNew, " mag^2 = ",drNew(1)**2 + drNew(2)**2 + drNew(3)**2
      print*, "Difference detected in test_equiv, 0"
      stop 1
endif

if (abs(drParOld-drParNew).gt.eps) then
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
      -(drPerpNew(1)**2 + drPerpNew(2)**2 + drPerpNew(3)**2)).gt.eps) then
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
      -(Ginew(1)**2 + Ginew(2)**2 + Ginew(3)**2)).gt.eps) then
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
!---------------------------------------------------------------!
