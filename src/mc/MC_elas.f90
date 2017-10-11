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
implicit none
type(wlcsim_params), intent(in) :: wlc_p
real(dp), intent(in) :: R(3,wlc_p%NT)  ! Bead positions
real(dp), intent(in) :: U(3,wlc_p%NT)  ! Tangent vectors
real(dp), intent(in) :: RP(3,wlc_p%NT)  ! Bead positions
real(dp), intent(in) :: UP(3,wlc_p%NT)  ! Tangent vectors
integer, intent(in) :: IB1               ! Test bead position 1
integer, intent(in) :: IT1               ! Index of test bead 1
integer, intent(in) :: IB2               ! Test bead position 2
integer, intent(in) :: IT2               ! Index of test bead 2

real(dp) :: U0(3,wlc_p%NT)  ! dummy variable for calculating u in wlc case
real(dp) :: UP0(3,wlc_p%NT)  ! dummy variable for calculating u in wlc case
real(dp) :: UNORM,U1U2

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

real(dp) DR(3),DRPAR,DRPERP(3)
real(dp) GI(3)

! Setup parameters

      DEELAS(1) = 0.0_dp ! bending energy
      DEELAS(2) = 0.0_dp ! parallel stretch energy
      DEELAS(3) = 0.0_dp ! perpendicular stretch energy
      DEELAS(4) = 0.0_dp ! wlc_p%twist energy

!     Calculate the change in the energy

      if ((IB1 /= 1).or.(wlc_p%ring)) then
          ! so if wlc_p%ring == 1, i.e.
          if (IB1 == 1) then
              ! then the bead to the left of IB1 is actually the end bead due to the wlc_p%ring inducing periodic boundaries on the array
              IT1M1 = wlc_p%NB
          else
              IT1M1 = IT1 - 1
          endif

         ! MC move only affects energy if it's interior to the polymer, since there are only crankshaft moves, and end
         ! crankshafts don't affect polymer
         if (wlc_p%SIMTYPE == 1.AND.(IB1 /= wlc_p%NB.OR.wlc_p%ring)) then
             if (IB1 == wlc_p%NB) then
                 IT1P1 = 1
             else
                 IT1P1 = IT1 + 1
             endif
             print*, "You will need to update this section before use."
             print*, "Finish implementing IT1 and IT2"
             stop 1

            U0(1,IT1M1) = R(1,IT1)-R(1,IT1M1)
            U0(2,IT1M1) = R(2,IT1)-R(2,IT1M1)
            U0(3,IT1M1) = R(3,IT1)-R(3,IT1M1)
            UNORM = sqrt(U0(1,IT1M1)**2. + U0(2,IT1M1)**2. + U0(3,IT1M1)**2.)
            U0(1,IT1M1) = U0(1,IT1M1)/UNORM
            U0(2,IT1M1) = U0(2,IT1M1)/UNORM
            U0(3,IT1M1) = U0(3,IT1M1)/UNORM

            U0(1,IT1) = R(1,IT1 + 1)-R(1,IT1)
            U0(2,IT1) = R(2,IT1 + 1)-R(2,IT1)
            U0(3,IT1) = R(3,IT1 + 1)-R(3,IT1)
            UNORM = sqrt(U0(1,IT1)**2. + U0(2,IT1)**2. + U0(3,IT1)**2.)
            U0(1,IT1) = U0(1,IT1)/UNORM
            U0(2,IT1) = U0(2,IT1)/UNORM
            U0(3,IT1) = U0(3,IT1)/UNORM

            UP0(1,IT1) = RP(1,IT1 + 1)-RP(1,IT1)
            UP0(2,IT1) = RP(2,IT1 + 1)-RP(2,IT1)
            UP0(3,IT1) = RP(3,IT1 + 1)-RP(3,IT1)
            UNORM = sqrt(UP0(1,IT1)**2. + UP0(2,IT1)**2. + UP0(3,IT1)**2.)
            UP0(1,IT1) = UP0(1,IT1)/UNORM
            UP0(2,IT1) = UP0(2,IT1)/UNORM
            UP0(3,IT1) = UP0(3,IT1)/UNORM

            U1U2 = U0(1,IT1M1)*U0(1,IT1) + U0(2,IT1M1)*U0(2,IT1) + U0(3,IT1M1)*U0(3,IT1)
            ! only bending energy in WLC, others are init to zero, so sum in total energy calculation later will be no-op
            DEELAS(1) = DEELAS(1) + EB*U1U2
            U1U2 = U(1,IT1M1)*UP0(1,IT1) + U(2,IT1M1)*UP(2,IT1) + U(3,IT1M1)*UP(3,IT1)
            DEELAS(1) = DEELAS(1)-EB*U1U2

         elseif (wlc_p%SIMTYPE == 2) then

            DR(1) = R(1,IT1)-R(1,IT1M1)
            DR(2) = R(2,IT1)-R(2,IT1M1)
            DR(3) = R(3,IT1)-R(3,IT1M1)
            DRPAR = DR(1)*U(1,IT1M1) + DR(2)*U(2,IT1M1) + DR(3)*U(3,IT1M1)

            DRPERP(1) = DR(1)-DRPAR*U(1,IT1M1)
            DRPERP(2) = DR(2)-DRPAR*U(2,IT1M1)
            DRPERP(3) = DR(3)-DRPAR*U(3,IT1M1)
            !U1U2 = U(1,IT1M1)*U(1,IT1) + U(2,IT1M1)*U(2,IT1) + U(3,IT1M1)*U(3,IT1)

            GI(1) = (U(1,IT1)-U(1,IT1M1)-ETA*DRPERP(1))
            GI(2) = (U(2,IT1)-U(2,IT1M1)-ETA*DRPERP(2))
            GI(3) = (U(3,IT1)-U(3,IT1M1)-ETA*DRPERP(3))

            DEELAS(1) = DEELAS(1)-0.5_dp*EB*(GI(1)**2 + GI(2)**2 + GI(3)**2)
            DEELAS(2) = DEELAS(2)-0.5_dp*EPAR*(DRPAR-GAM)**2
            DEELAS(3) = DEELAS(3)-0.5_dp*EPERP*(DRPERP(1)**2 + DRPERP(2)**2. + DRPERP(3)**2)

            DR(1) = RP(1,IT1)-R(1,IT1M1)
            DR(2) = RP(2,IT1)-R(2,IT1M1)
            DR(3) = RP(3,IT1)-R(3,IT1M1)
            DRPAR = DR(1)*U(1,IT1M1) + DR(2)*U(2,IT1M1) + DR(3)*U(3,IT1M1)

            DRPERP(1) = DR(1)-DRPAR*U(1,IT1M1)
            DRPERP(2) = DR(2)-DRPAR*U(2,IT1M1)
            DRPERP(3) = DR(3)-DRPAR*U(3,IT1M1)
            !U1U2 = U(1,IT1M1)*UP(1,IT1) + U(2,IT1M1)*UP(2,IT1) + U(3,IT1M1)*UP(3,IT1)

            GI(1) = (UP(1,IT1)-U(1,IT1M1)-ETA*DRPERP(1))
            GI(2) = (UP(2,IT1)-U(2,IT1M1)-ETA*DRPERP(2))
            GI(3) = (UP(3,IT1)-U(3,IT1M1)-ETA*DRPERP(3))

            ! in ssWLC, we have, bending, parallel extension, and perpendicular extension energy, which we must sum later to get
            ! total energy
            DEELAS(1) = DEELAS(1) + 0.5_dp*EB*(GI(1)**2 + GI(2)**2 + GI(3)**2)
            DEELAS(2) = DEELAS(2) + 0.5_dp*EPAR*(DRPAR-GAM)**2
            DEELAS(3) = DEELaS(3) + 0.5_dp*EPERP*(DRPERP(1)**2 + DRPERP(2)**2 + DRPERP(3)**2)

         elseif (wlc_p%SIMTYPE == 3) then

            DR(1) = R(1,IT1)-R(1,IT1M1)
            DR(2) = R(2,IT1)-R(2,IT1M1)
            DR(3) = R(3,IT1)-R(3,IT1M1)
            ! in gaussian chain, there's only parallel stretching energy. DEELAS init'd to zeros, so sum(DEELAS) == DEELAS(2) later
            DEELAS(2) = DEELAS(2)-0.5*EPAR*(DR(1)**2. + DR(2)**2. + DR(3)**2.)
            DR(1) = RP(1,IT1)-R(1,IT1M1)
            DR(2) = RP(2,IT1)-R(2,IT1M1)
            DR(3) = RP(3,IT1)-R(3,IT1M1)
            DEELAS(2) = DEELAS(2) + 0.5*EPAR*(DR(1)**2. + DR(2)**2. + DR(3)**2.)

         endif
      endif

      if ((IB2 /= wlc_p%NB).or.(wlc_p%ring)) then
         if (IB2 == wlc_p%NB) then
             IT2P1 = 1
         else
             IT2P1 = IT2 + 1
         endif

         ! if we're talking about a WLC, if we crankshaft a single bead, that's a no-op, since the u's are directly
         ! determined by the r's. Thus we're not worried about double counting the energy change here since the energy change
         ! should be zero by definition if IB1 = =IB2.
         if (wlc_p%SIMTYPE == 1.AND.((IB2 /= 1).OR.(wlc_p%ring))) then
            if (IB2 == 1) then
                IT2M1 = wlc_p%NB
            else
                IT2M1 = IT2 - 1
            endif
            Print*, "This section is out of date"
            print*, "The variable IT2M1 is never used!"
            stop
            U0(1,IT2-1) = R(1,IT2)-R(1,IT2-1)
            U0(2,IT2-1) = R(2,IT2)-R(2,IT2-1)
            U0(3,IT2-1) = R(3,IT2)-R(3,IT2-1)
            UNORM = sqrt(U(1,IT2-1)**2. + U(2,IT2-1)**2. + U(3,IT2-1)**2.)
            U0(1,IT2-1) = U0(1,IT2-1)/UNORM
            U0(2,IT2-1) = U0(2,IT2-1)/UNORM
            U0(3,IT2-1) = U0(3,IT2-1)/UNORM

            U0(1,IT2) = R(1,IT2 + 1)-R(1,IT2)
            U0(2,IT2) = R(2,IT2 + 1)-R(2,IT2)
            U0(3,IT2) = R(3,IT2 + 1)-R(3,IT2)
            UNORM = sqrt(U0(1,IT2)**2. + U0(2,IT2)**2. + U0(3,IT2)**2.)
            U0(1,IT2) = U0(1,IT2)/UNORM
            U0(2,IT2) = U0(2,IT2)/UNORM
            U0(3,IT2) = U0(3,IT2)/UNORM

            UP0(1,IT2-1) = RP(1,IT2)-RP(1,IT2-1)
            UP0(2,IT2-1) = RP(2,IT2)-RP(2,IT2-1)
            UP0(3,IT2-1) = RP(3,IT2)-RP(3,IT2-1)
            UNORM = sqrt(UP0(1,IT2-1)**2. + UP0(2,IT2-1)**2. + UP0(3,IT2-1)**2.)
            UP0(1,IT2-1) = UP0(1,IT2-1)/UNORM
            UP0(2,IT2-1) = UP0(2,IT2-1)/UNORM
            UP0(3,IT2-1) = UP0(3,IT2-1)/UNORM

            U1U2 = U0(1,IT2-1)*U0(1,IT2) + U0(2,IT2-1)*U0(2,IT2) + U0(3,IT2-1)*U0(3,IT2)
            DEELAS(1) = DEELAS(1) + EB*U1U2
            U1U2 = UP0(1,IT2-1)*U0(1,IT2) + UP0(2,IT2-1)*U0(2,IT2) + UP0(3,IT2-1)*U0(3,IT2)
            DEELAS(1) = DEELAS(1)-EB*U1U2

         elseif (wlc_p%SIMTYPE == 2) then

            DR(1) = R(1,IT2P1)-R(1,IT2)
            DR(2) = R(2,IT2P1)-R(2,IT2)
            DR(3) = R(3,IT2P1)-R(3,IT2)
            DRPAR = DR(1)*U(1,IT2) + DR(2)*U(2,IT2) + DR(3)*U(3,IT2)

            DRPERP(1) = DR(1)-DRPAR*U(1,IT2)
            DRPERP(2) = DR(2)-DRPAR*U(2,IT2)
            DRPERP(3) = DR(3)-DRPAR*U(3,IT2)
            !U1U2 = U(1,IT2)*U(1,IT2 + 1) + U(2,IT2)*U(2,IT2 + 1) + U(3,IT2)*U(3,IT2 + 1)

            GI(1) = (U(1,IT2P1)-U(1,IT2)-ETA*DRPERP(1))
            GI(2) = (U(2,IT2P1)-U(2,IT2)-ETA*DRPERP(2))
            GI(3) = (U(3,IT2P1)-U(3,IT2)-ETA*DRPERP(3))

            DEELAS(1) = DEELAS(1)-0.5_dp*EB*(GI(1)**2. + GI(2)**2. + GI(3)**2.)
            DEELAS(2) = DEELAS(2)-0.5_dp*EPAR*(DRPAR-GAM)**2.
            DEELAS(3) = DEELAS(3)-0.5_dp*EPERP*(DRPERP(1)**2. + DRPERP(2)**2. + DRPERP(3)**2.)

            DR(1) = R(1,IT2P1)-RP(1,IT2)
            DR(2) = R(2,IT2P1)-RP(2,IT2)
            DR(3) = R(3,IT2P1)-RP(3,IT2)
            DRPAR = DR(1)*UP(1,IT2) + DR(2)*UP(2,IT2) + DR(3)*UP(3,IT2)

            DRPERP(1) = DR(1)-DRPAR*UP(1,IT2)
            DRPERP(2) = DR(2)-DRPAR*UP(2,IT2)
            DRPERP(3) = DR(3)-DRPAR*UP(3,IT2)
            !U1U2 = UP(1,IT2)*U(1,IT2 + 1) + UP(2,IT2)*U(2,IT2 + 1) + UP(3,IT2)*U(3,IT2 + 1)

            GI(1) = (U(1,IT2P1)-UP(1,IT2)-ETA*DRPERP(1))
            GI(2) = (U(2,IT2P1)-UP(2,IT2)-ETA*DRPERP(2))
            GI(3) = (U(3,IT2P1)-UP(3,IT2)-ETA*DRPERP(3))

            DEELAS(1) = DEELAS(1) + 0.5_dp*EB*(GI(1)**2. + GI(2)**2 + GI(3)**2)
            DEELAS(2) = DEELAS(2) + 0.5_dp*EPAR*(DRPAR-GAM)**2.
            DEELAS(3) = DEELAS(3) + 0.5_dp*EPERP*(DRPERP(1)**2. + DRPERP(2)**2. + DRPERP(3)**2.)

         elseif (wlc_p%SIMTYPE == 3) then

            DR(1) = R(1,IT2P1)-R(1,IT2)
            DR(2) = R(2,IT2P1)-R(2,IT2)
            DR(3) = R(3,IT2P1)-R(3,IT2)
            DEELAS(2) = DEELAS(2)-0.5*EPAR*(DR(1)**2. + DR(2)**2. + DR(3)**2.)
            DR(1) = R(1,IT2P1)-RP(1,IT2)
            DR(2) = R(2,IT2P1)-RP(2,IT2)
            DR(3) = R(3,IT2P1)-RP(3,IT2)
            DEELAS(2) = DEELAS(2) + 0.5*EPAR*(DR(1)**2. + DR(2)**2. + DR(3)**2.)

         endif

      endif

      if (wlc_p%ring.AND.wlc_p%twIST) then
          if (MCTYPE == 1) then
              CALL WRITHECRANK(R,IT1,IT2,wlc_p%NB,WRM)
              CALL WRITHECRANK(RP,IT1,IT2,wlc_p%NB,WRMP)
              DWR = WRMP-WRM
          elseif (MCTYPE == 2) then
              CALL WRITHESLIDE(R,IT1,IT2,wlc_p%NB,WRM)
              CALL WRITHESLIDE(RP,IT1,IT2,wlc_p%NB,WRMP)
              DWR = WRMP-WRM
          else
              DWR = 0.
          ENDif
          WRP = WR + DWR
          tw = REAL(wlc_p%LK)-WR
          twP = REAL(wlc_p%LK)-WRP
          DEELAS(4) = DEELAS(4) + (((2.*pi*twP)**2.)*wlc_p%lt/(2.*wlc_p%L))-(((2.*pi*TW)**2.)*wlc_p%LT/(2.*wlc_p%L))
      ENDif

      RETURN
      END

!---------------------------------------------------------------*
