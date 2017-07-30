!---------------------------------------------------------------*

!
!     This subroutine calculates the elastic forces for a wormlike
!     chain with a stretching potential.  The stretch and bend
!     moduli are fed along with the bead positions.
!
!     Andrew Spakowitz
!     Written 9-1-04

      subroutine force_elas(FELAS,TELAS,R,U,NT,N,NP, EB,EPAR,EPERP,GAM,ETA,SIMTYPE)
      use params, only : dp
      implicit none
      real(dp), intent(out) :: FELAS(NT,3) ! Elastic force
      real(dp), intent(out) :: TELAS(NT,3) ! Elastic force
      real(dp), intent(in) ::  R(NT,3)  ! Bead positions
      ! u will be modified iff simtype = 1, if want this to be intent(in), just
      ! modify the below to use a separate variable for this
      real(dp), intent(inout) :: U(NT,3)  ! Tangent vectors
      real(dp) B(NT)    ! Bond lengths
      integer I,J,IB            ! Index holders
      integer, intent(in) :: N,NT,NP  ! Number of bead
      integer SIMTYPE           ! Simulation method (WLC = 1,SSWLC = 2,GC = 3)

!     Polymer properties

      real(dp), intent(in) :: EB,EPAR,EPERP
      real(dp), intent(in) :: GAM,ETA

!     Variables for force and torque calculations

      real(dp) DR(3),DRPAR,DRPERP(3)
      real(dp) FI(3),TI1(3),TI2(3)
      real(dp) U1U2,GI(3)

      IB = 1
      do 10 I = 1,NP
         do 20 J = 1,N
            FELAS(IB,1) = 0.
            FELAS(IB,2) = 0.
            FELAS(IB,3) = 0.
            TELAS(IB,1) = 0.
            TELAS(IB,2) = 0.
            TELAS(IB,3) = 0.
            if (SIMTYPE == 1.AND.J <= (N-1)) then
               U(IB,1) = R(IB + 1,1)-R(IB,1)
               U(IB,2) = R(IB + 1,2)-R(IB,2)
               U(IB,3) = R(IB + 1,3)-R(IB,3)
               B(IB) = sqrt(U(IB,1)**2. + U(IB,2)**2. + U(IB,3)**2.)
               U(IB,1) = U(IB,1)/B(IB)
               U(IB,2) = U(IB,2)/B(IB)
               U(IB,3) = U(IB,3)/B(IB)
            endif
            IB = IB + 1

 20      continue
 10   continue

!     Calculate the forces and torques

      do 30 I = 1,NP
         do 40 J = 1,(N-1)
            IB = J + N*(I-1)
            if (SIMTYPE == 1) then

               if (J <= (N-2)) then

                  U1U2 = U(IB,1)*U(IB + 1,1) + U(IB,2)*U(IB + 1,2) + U(IB,3)*U(IB + 1,3)

                  GI(1) = EB*(U(IB + 1,1)-U1U2*U(IB,1))/B(IB)
                  GI(2) = EB*(U(IB + 1,2)-U1U2*U(IB,2))/B(IB)
                  GI(3) = EB*(U(IB + 1,3)-U1U2*U(IB,3))/B(IB)

                  FELAS(IB,1) = FELAS(IB,1)-GI(1)
                  FELAS(IB,2) = FELAS(IB,2)-GI(2)
                  FELAS(IB,3) = FELAS(IB,3)-GI(3)
                  FELAS(IB + 1,1) = FELAS(IB + 1,1) + GI(1)
                  FELAS(IB + 1,2) = FELAS(IB + 1,2) + GI(2)
                  FELAS(IB + 1,3) = FELAS(IB + 1,3) + GI(3)

                  GI(1) = EB*(U(IB,1)-U1U2*U(IB + 1,1))/B(IB + 1)
                  GI(2) = EB*(U(IB,2)-U1U2*U(IB + 1,2))/B(IB + 1)
                  GI(3) = EB*(U(IB,3)-U1U2*U(IB + 1,3))/B(IB + 1)

                  FELAS(IB + 1,1) = FELAS(IB + 1,1)-GI(1)
                  FELAS(IB + 1,2) = FELAS(IB + 1,2)-GI(2)
                  FELAS(IB + 1,3) = FELAS(IB + 1,3)-GI(3)
                  FELAS(IB + 2,1) = FELAS(IB + 2,1) + GI(1)
                  FELAS(IB + 2,2) = FELAS(IB + 2,2) + GI(2)
                  FELAS(IB + 2,3) = FELAS(IB + 2,3) + GI(3)

               endif
               TELAS(IB,1) = 0.
               TELAS(IB,2) = 0.
               TELAS(IB,3) = 0.
               TELAS(IB + 1,1) = 0.
               TELAS(IB + 1,2) = 0.
               TELAS(IB + 1,3) = 0.

            elseif (SIMTYPE == 2) then
               DR(1) = R(IB + 1,1)-R(IB,1)
               DR(2) = R(IB + 1,2)-R(IB,2)
               DR(3) = R(IB + 1,3)-R(IB,3)
               DRPAR = DR(1)*U(IB,1) + DR(2)*U(IB,2) + DR(3)*U(IB,3)

               DRPERP(1) = DR(1)-DRPAR*U(IB,1)
               DRPERP(2) = DR(2)-DRPAR*U(IB,2)
               DRPERP(3) = DR(3)-DRPAR*U(IB,3)
               U1U2 = U(IB,1)*U(IB + 1,1) + U(IB,2)*U(IB + 1,2) + U(IB,3)*U(IB + 1,3)

               GI(1) = U(IB + 1,1)-U1U2*U(IB,1)-ETA*DRPERP(1)
               GI(2) = U(IB + 1,2)-U1U2*U(IB,2)-ETA*DRPERP(2)
               GI(3) = U(IB + 1,3)-U1U2*U(IB,3)-ETA*DRPERP(3)

               FI(1) = -ETA*EB*GI(1) + EPAR*(DRPAR-GAM)*U(IB,1) + EPERP*DRPERP(1)
               FI(2) = -ETA*EB*GI(2) + EPAR*(DRPAR-GAM)*U(IB,2) + EPERP*DRPERP(2)
               FI(3) = -ETA*EB*GI(3) + EPAR*(DRPAR-GAM)*U(IB,3) + EPERP*DRPERP(3)

               FELAS(IB,1) = FELAS(IB,1) + FI(1)
               FELAS(IB,2) = FELAS(IB,2) + FI(2)
               FELAS(IB,3) = FELAS(IB,3) + FI(3)
               FELAS(IB + 1,1) = FELAS(IB + 1,1)-FI(1)
               FELAS(IB + 1,2) = FELAS(IB + 1,2)-FI(2)
               FELAS(IB + 1,3) = FELAS(IB + 1,3)-FI(3)

               GI(1) = U(IB + 1,1)-U(IB,1)-ETA*DRPERP(1)
               GI(2) = U(IB + 1,2)-U(IB,2)-ETA*DRPERP(2)
               GI(3) = U(IB + 1,3)-U(IB,3)-ETA*DRPERP(3)

               TI1(1) = EB*GI(1)
               TI1(2) = EB*GI(2)
               TI1(3) = EB*GI(3)

               TI2(1) = -ETA*EB*DRPAR*GI(1) + ETA*EB*(1-U1U2)*DR(1)-EPAR*(DRPAR-GAM)*DR(1) + EPERP*DRPAR*DRPERP(1)
               TI2(2) = -ETA*EB*DRPAR*GI(2) + ETA*EB*(1-U1U2)*DR(2)-EPAR*(DRPAR-GAM)*DR(2) + EPERP*DRPAR*DRPERP(2)
               TI2(3) = -ETA*EB*DRPAR*GI(3) + ETA*EB*(1-U1U2)*DR(3)-EPAR*(DRPAR-GAM)*DR(3) + EPERP*DRPAR*DRPERP(3)

               TELAS(IB,1) = TELAS(IB,1) + TI1(1) + TI2(1)
               TELAS(IB,2) = TELAS(IB,2) + TI1(2) + TI2(2)
               TELAS(IB,3) = TELAS(IB,3) + TI1(3) + TI2(3)
               TELAS(IB + 1,1) = TELAS(IB + 1,1)-TI1(1)
               TELAS(IB + 1,2) = TELAS(IB + 1,2)-TI1(2)
               TELAS(IB + 1,3) = TELAS(IB + 1,3)-TI1(3)

            elseif (SIMTYPE == 3) then

               DR(1) = R(IB + 1,1)-R(IB,1)
               DR(2) = R(IB + 1,2)-R(IB,2)
               DR(3) = R(IB + 1,3)-R(IB,3)

               FELAS(IB,1) = FELAS(IB,1) + EPAR*DR(1)
               FELAS(IB,2) = FELAS(IB,2) + EPAR*DR(2)
               FELAS(IB,3) = FELAS(IB,3) + EPAR*DR(3)
               FELAS(IB + 1,1) = FELAS(IB + 1,1)-EPAR*DR(1)
               FELAS(IB + 1,2) = FELAS(IB + 1,2)-EPAR*DR(2)
               FELAS(IB + 1,3) = FELAS(IB + 1,3)-EPAR*DR(3)

               TELAS(IB,1) = 0.
               TELAS(IB,2) = 0.
               TELAS(IB,3) = 0.
               TELAS(IB + 1,1) = 0.
               TELAS(IB + 1,2) = 0.
               TELAS(IB + 1,3) = 0.

            endif

 40      continue
 30   continue

      RETURN
      END

!---------------------------------------------------------------*
