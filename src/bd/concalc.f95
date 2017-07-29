!---------------------------------------------------------------*

!
!   This subroutine performs the constraint forces
!
!   Andrew Spakowitz
!   Written 9-8-04

      subroutine concalc(R,DRDT,NT,N,NP,XI,L0,DT,RK,BROWN)

      use params, only : dp

! Variables from the simulation

      real(dp) DRDT(NT,3,4)  ! Rate of change
      real(dp) R(NT,3)       ! Bead positions
      real(dp) U(N-1,3)      ! Unit tangent vector
      real(dp) B(N-1)         ! Bond length
      integer N,NT,NP                 ! Number of beads
      real(dp) DT         ! Time step size
      real(dp) XI            ! Drag coefficient
      real(dp) L0            ! Bond distances
      integer I,J,IB                        ! Index holder
      integer RK                    ! RK integer
      integer BROWN             ! Logic for BD forces

! Variables for use in the constraint calculation

      real(dp) BLAM(N-1),ADIAG(N-1)
      real(dp) ASUPER(N-1),ASUB(N-1)
      real(dp) C
      integer inFO

! Variables of use in the pseudopotential calculation

      real(dp) FPS(N,3),APS(N-1),BPS(N-2)
      real(dp) DETLT(N-1),DETGT(N-1),DETER

! Relaxation parameter

      C = 0.1/DT

! Calculate the bond length and tangent vectors

      do 5 J = 1,NP
         IB = N*(J-1)

         do 10 I = 1,(N-1)
            U(I,1) = R(I + IB + 1,1)-R(I + IB,1)
            U(I,2) = R(I + IB + 1,2)-R(I + IB,2)
            U(I,3) = R(I + IB + 1,3)-R(I + IB,3)
            B(I) = sqrt(U(I,1)**2. + U(I,2)**2. + U(I,3)**2.)
            U(I,1) = U(I,1)/B(I)
            U(I,2) = U(I,2)/B(I)
            U(I,3) = U(I,3)/B(I)
 10      continue

! Setup the A component

         do 15 I = 1,(N-1)
            ADIAG(I) = 2.*B(I)*B(I)
            if (BROWN == 1) then
               APS(I) = ADIAG(I)
            endif
            if (I >= 2) then
               ASUB(I) = -B(I)*B(I-1)*(U(I,1)*U(I-1,1)+ &
                    U(I,2)*U(I-1,2) + U(I,3)*U(I-1,3))
               ASUPER(I-1) = ASUB(I)
               if (BROWN == 1) then
                  BPS(I-1) = ASUPER(I-1)
               endif
            endif
15       continue
         ASUB(1) = 0.
         ASUPER(N-1) = 0.

! Calculate psuedo-potential (if Brownian forces on)

         if (BROWN == 1) then
            DETLT(1) = 1.0
            DETLT(2) = APS(1)
            DETGT(N-1) = 1.0
            DETGT(N-2) = APS(N-1)
            do 20 I = 3,(N-1)
               DETLT(I) = APS(I-1)*DETLT(I-1)- &
                    BPS(I-2)*BPS(I-2)*DETLT(I-2)
               DETGT(N-I) = APS(N-I + 1)*DETGT(N-I + 1)- &
                    BPS(N-I + 1)*BPS(N-I + 1)*DETGT(N-I + 2)
 20         continue
            DETER = APS(N-1)*DETLT(N-1)-BPS(N-2)*BPS(N-2)*DETLT(N-2)

            FPS(1,1) = DETGT(2)*BPS(1)*B(2)*U(2,1)/DETER
            FPS(1,2) = DETGT(2)*BPS(1)*B(2)*U(2,2)/DETER
            FPS(1,3) = DETGT(2)*BPS(1)*B(2)*U(2,3)/DETER
            FPS(2,1) = DETGT(2)*BPS(1)*(B(1)*U(1,1)- &
                 B(2)*U(2,1))/DETER+ &
                 DETGT(3)*APS(1)*BPS(2)*B(3)*U(3,1)/DETER
            FPS(2,2) = DETGT(2)*BPS(1)*(B(1)*U(1,2)- &
                 B(2)*U(2,2))/DETER+ &
                 DETGT(3)*APS(1)*BPS(2)*B(3)*U(3,2)/DETER
            FPS(2,3) = DETGT(2)*BPS(1)*(B(1)*U(1,3)- &
                 B(2)*U(2,3))/DETER+ &
                 DETGT(3)*APS(1)*BPS(2)*B(3)*U(3,3)/DETER
            do 30 I = 3,(N-2)
               FPS(I,1) = (DETLT(I-1)*DETGT(I)*BPS(I-1)*(-B(I)*U(I,1) + B(I-1)*U(I-1,1))+ &
                    DETLT(I-1)*DETGT(I + 1)*APS(I-1)*BPS(I)*B(I + 1)*U(I + 1,1)- &
                    DETLT(I-2)*DETGT(I)*APS(I)*BPS(I-2)*B(I-2)*U(I-2,1)- &
                    DETLT(I-2)*DETGT(I + 1)*( &
                    (BPS(I-2)**2.)*BPS(I)*B(I + 1)*U(I + 1,1)- &
                    (BPS(I)**2.)*BPS(I-2)*B(I-2)*U(I-2,1)))/DETER
               FPS(I,2) = (DETLT(I-1)*DETGT(I)*BPS(I-1)*(-B(I)*U(I,2) + B(I-1)*U(I-1,2))+ &
                    DETLT(I-1)*DETGT(I + 1)*APS(I-1)*BPS(I)*B(I + 1)*U(I + 1,2)- &
                    DETLT(I-2)*DETGT(I)*APS(I)*BPS(I-2)*B(I-2)*U(I-2,2)- &
                    DETLT(I-2)*DETGT(I + 1)*( &
                    (BPS(I-2)**2.)*BPS(I)*B(I + 1)*U(I + 1,2)- &
                    (BPS(I)**2.)*BPS(I-2)*B(I-2)*U(I-2,2)))/DETER
               FPS(I,3) = (DETLT(I-1)*DETGT(I)*BPS(I-1)*(-B(I)*U(I,3) + B(I-1)*U(I-1,3))+ &
                    DETLT(I-1)*DETGT(I + 1)*APS(I-1)*BPS(I)*B(I + 1)*U(I + 1,3)- &
                    DETLT(I-2)*DETGT(I)*APS(I)*BPS(I-2)*B(I-2)*U(I-2,3)- &
                    DETLT(I-2)*DETGT(I + 1)*( &
                    (BPS(I-2)**2.)*BPS(I)*B(I + 1)*U(I + 1,3)- &
                    (BPS(I)**2.)*BPS(I-2)*B(I-2)*U(I-2,3)))/DETER
 30         continue
            FPS(N-1,1) = DETLT(N-2)*BPS(N-2)*(-B(N-1)*U(N-1,1)+ &
                 B(N-2)*U(N-2,1))/DETER- &
                 DETLT(N-3)*APS(N-1)*BPS(N-3)*B(N-3)*U(N-3,1)/DETER
            FPS(N-1,2) = DETLT(N-2)*BPS(N-2)*(-B(N-1)*U(N-1,2)+ &
                 B(N-2)*U(N-2,2))/DETER- &
                 DETLT(N-3)*APS(N-1)*BPS(N-3)*B(N-3)*U(N-3,2)/DETER
            FPS(N-1,3) = DETLT(N-2)*BPS(N-2)*(-B(N-1)*U(N-1,3)+ &
                 B(N-2)*U(N-2,3))/DETER- &
                 DETLT(N-3)*APS(N-1)*BPS(N-3)*B(N-3)*U(N-3,3)/DETER
            FPS(N,1) = -DETLT(N-2)*BPS(N-2)*B(N-2)*U(N-2,1)/DETER
            FPS(N,2) = -DETLT(N-2)*BPS(N-2)*B(N-2)*U(N-2,2)/DETER
            FPS(N,3) = -DETLT(N-2)*BPS(N-2)*B(N-2)*U(N-2,3)/DETER

            do 35 I = 1,N
               DRDT(I + IB,1,RK) = DRDT(I + IB,1,RK) + FPS(I,1)/XI
               DRDT(I + IB,2,RK) = DRDT(I + IB,2,RK) + FPS(I,2)/XI
               DRDT(I + IB,3,RK) = DRDT(I + IB,3,RK) + FPS(I,3)/XI
 35         continue
         endif


! Setup the B component

         do 40 I = 1,(N-1)
            BLAM(I) = B(I)*XI*(U(I,1)*(DRDT(I + 1 + IB,1,RK)-DRDT(I + IB,1,RK))+ &
                 U(I,2)*(DRDT(I + 1 + IB,2,RK)-DRDT(I + IB,2,RK))+ &
                 U(I,3)*(DRDT(I + 1 + IB,3,RK)-DRDT(I + IB,3,RK)))+ &
                 XI*C*(B(I)**.2-L0**2.)
 40      continue

! Calculate the new rates of change

         call DGTSV((N-1),1,ASUB,ADIAG,ASUPER,BLAM,(N-1),inFO)

         DRDT(1 + IB,1,RK) = DRDT(1 + IB,1,RK) + BLAM(1)*B(1)*U(1,1)/XI
         DRDT(1 + IB,2,RK) = DRDT(1 + IB,2,RK) + BLAM(1)*B(1)*U(1,2)/XI
         DRDT(1 + IB,3,RK) = DRDT(1 + IB,3,RK) + BLAM(1)*B(1)*U(1,3)/XI
         do 50 I = 2,(N-1)
            DRDT(I + IB,1,RK) = DRDT(I + IB,1,RK) + (BLAM(I)*B(I)*U(I,1)- &
                 BLAM(I-1)*B(I-1)*U(I-1,1))/XI
            DRDT(I + IB,2,RK) = DRDT(I + IB,2,RK) + (BLAM(I)*B(I)*U(I,2)- &
                 BLAM(I-1)*B(I-1)*U(I-1,2))/XI
            DRDT(I + IB,3,RK) = DRDT(I + IB,3,RK) + (BLAM(I)*B(I)*U(I,3)- &
                 BLAM(I-1)*B(I-1)*U(I-1,3))/XI
 50      continue
         DRDT(N + IB,1,RK) = DRDT(N + IB,1,RK)-BLAM(N-1)*B(N-1)*U(N-1,1)/XI
         DRDT(N + IB,2,RK) = DRDT(N + IB,2,RK)-BLAM(N-1)*B(N-1)*U(N-1,2)/XI
         DRDT(N + IB,3,RK) = DRDT(N + IB,3,RK)-BLAM(N-1)*B(N-1)*U(N-1,3)/XI

 5    continue

      RETURN
      END

!---------------------------------------------------------------*
