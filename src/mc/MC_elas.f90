#include "../defines.inc"
!---------------------------------------------------------------*

!     subroutine MC_eelas
!
!     Calculate the change in the polymer elastic energy
!     due to the displacement from a MC move
!
subroutine MC_eelas(wlc_p,DEELAS,R,U,RP,UP,IB1,IB2,&
                    IT1,IT2,EB,EPAR,EPERP,GAM,ETA, &
                    MCTYPE,WR,WRP)

use params, only: dp, pi,wlcsim_params
use MC_wlc, only: E_wlc, E_SSWLC, E_GAUSS
implicit none
type(wlcsim_params), intent(in) :: wlc_p
real(dp), intent(in) :: R(3,WLC_P__NT)  ! Bead positions
real(dp), intent(in) :: U(3,WLC_P__NT)  ! Tangent vectors
real(dp), intent(in) :: RP(3,WLC_P__NT)  ! Bead positions
real(dp), intent(in) :: UP(3,WLC_P__NT)  ! Tangent vectors
integer, intent(in) :: IB1               ! Test bead position 1
integer, intent(in) :: IT1               ! Index of test bead 1
integer, intent(in) :: IB2               ! Test bead position 2
integer, intent(in) :: IT2               ! Index of test bead 2

real(dp), intent(out) :: DEELAS(4)   ! Change in ECOM

!     Polymer properties

real(dp), intent(in) :: EB
real(dp), intent(in) :: EPAR
real(dp), intent(in) :: EPERP
real(dp), intent(in) :: GAM
real(dp), intent(in) :: ETA
real(dp) tw       ! Twist
real(dp) twP      ! Twist of test structure
real(dp) WR       ! Writhe
real(dp), intent(out) :: WRP      ! Writhe of test structure
real(dp) DWR      ! Change in Writhe
real(dp) WRM,WRMP ! Component of writhe affected by move
integer IT1M1
integer IT1P1
integer IT2M1
integer IT2P1
integer MCTYPE            ! MC move type

! Setup parameters

      DEELAS(1) = 0.0_dp ! bending energy
      DEELAS(2) = 0.0_dp ! parallel stretch energy
      DEELAS(3) = 0.0_dp ! perpendicular stretch energy
      DEELAS(4) = 0.0_dp ! WLC_P__TWIST energy

!     Calculate the change in the energy

      if ((IB1 /= 1).or.(WLC_P__RING)) then
          ! so if WLC_P__RING == 1, i.e.
          if (IB1 == 1) then
              ! then the bead to the left of IB1 is actually the end bead due to the WLC_P__RING inducing periodic boundaries on the array
              IT1M1 = WLC_P__NB
          else
              IT1M1 = IT1 - 1
          endif

         ! MC move only affects energy if it's interior to the polymer, since there are only crankshaft moves, and end
         ! crankshafts don't affect polymer
         if (wlc_p%SIMTYPE == 1.AND.(IB1 /= WLC_P__NB.OR.WLC_P__RING)) then
             if (IB1 == WLC_P__NB) then
                 IT1P1 = 1
             else
                 IT1P1 = IT1 + 1
             endif
             print*, "You will need to update this section before use."
             print*, "Finish implementing IT1 and IT2"
             stop 1
             DEELAS(1) = DEELAS(1) - E_wlc(R(:,IT1M1), RP(:,IT1), RP(:,IT1P1), EB )
             DEELAS(1) = DEELAS(1) - E_wlc(R(:,IT1M1), R(:,IT1), R(:,IT1P1), EB)

         elseif (wlc_p%SIMTYPE == 2) then
             !function E_SSWLC(R,RM1,U,UM1,EB,EPAR,EPERP,ETA,GAM)
             DEELAS = DEELAS + E_SSWLC(RP(:,IT1),R(:,IT1M1),UP(:,IT1),U(:,IT1M1),EB,EPAR,EPERP,ETA,GAM)
             DEELAS = DEELAS - E_SSWLC( R(:,IT1),R(:,IT1M1), U(:,IT1),U(:,IT1M1),EB,EPAR,EPERP,ETA,GAM)

         elseif (wlc_p%SIMTYPE == 3) then
             DEELAS(2) = DEELAS(2) + E_GAUSS(RP(:,IT1),R(:,IT1M1),EPAR)
             DEELAS(2) = DEELAS(2) - E_GAUSS( R(:,IT1),R(:,IT1M1),EPAR)
         endif
      endif

      if ((IB2 /= WLC_P__NB).or.(WLC_P__RING)) then
         if (IB2 == WLC_P__NB) then
             IT2P1 = 1
         else
             IT2P1 = IT2 + 1
         endif

         ! if we're talking about a WLC, if we crankshaft a single bead, that's a no-op, since the u's are directly
         ! determined by the r's. Thus we're not worried about double counting the energy change here since the energy change
         ! should be zero by definition if IB1 = =IB2.
         if (wlc_p%SIMTYPE == 1.AND.((IB2 /= 1).OR.(WLC_P__RING))) then
             if (IB2 == 1) then
                 IT2M1 = WLC_P__NB
             else
                 IT2M1 = IT2 - 1
             endif
             Print*, "This section is out of date"
             print*, "The variable IT2M1 is never used!"
             stop
             DEELAS(1) = DEELAS(1) - E_wlc(RP(:,IT2M1),RP(:,IT2),R(:,IT2P1),EB)
             DEELAS(1) = DEELAS(1) - E_wlc(R(:,IT2M1),R(:,IT2),R(:,IT2P1),EB)

         elseif (wlc_p%SIMTYPE == 2) then
             !function E_SSWLC(R,RM1,U,UM1,EB,EPAR,EPERP,ETA,GAM)
             DEELAS = DEELAS + E_SSWLC(R(:,IT2P1),RP(:,IT2),U(:,IT2P1),UP(:,IT2),EB,EPAR,EPERP,ETA,GAM)
             DEELAS = DEELAS - E_SSWLC(R(:,IT2P1), R(:,IT2),U(:,IT2P1), U(:,IT2),EB,EPAR,EPERP,ETA,GAM)

         elseif (wlc_p%SIMTYPE == 3) then
             DEELAS(2) = DEELAS(2) + E_GAUSS(R(:,IT2P1),RP(:,IT2),EPAR)
             DEELAS(2) = DEELAS(2) - E_GAUSS(R(:,IT2P1), R(:,IT2),EPAR)
         endif

      endif

      if (WLC_P__RING.AND.WLC_P__TWIST) then
          if (MCTYPE == 1) then
              CALL WRITHECRANK(R,IT1,IT2,WLC_P__NB,WRM)
              CALL WRITHECRANK(RP,IT1,IT2,WLC_P__NB,WRMP)
              DWR = WRMP-WRM
          elseif (MCTYPE == 2) then
              CALL WRITHESLIDE(R,IT1,IT2,WLC_P__NB,WRM)
              CALL WRITHESLIDE(RP,IT1,IT2,WLC_P__NB,WRMP)
              DWR = WRMP-WRM
          else
              DWR = 0.
          ENDif
          WRP = WR + DWR
          tw = REAL(wlc_p%LK)-WR
          twP = REAL(wlc_p%LK)-WRP
          DEELAS(4) = DEELAS(4) + (((2.*pi*twP)**2.)*WLC_P__LT/(2.*WLC_P__L))-(((2.*pi*TW)**2.)*WLC_P__LT/(2.*WLC_P__L))
      ENDif

      RETURN
      END

!---------------------------------------------------------------*
