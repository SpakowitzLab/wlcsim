!---------------------------------------------------------------*

!
!     This subroutine performs the RK step for user inputted timestep
!     size
!
!     Andrew Spakowitz
!     Written 6-6-04

      subroutine RKstep(RS,R,US,U,DRDT,DUDT,NT,N,NP,RK,DT)

      use params, only : dp, pi

      real(dp) RS(NT,3)  ! Saved bead positions
      real(dp) R(NT,3)  ! Temp bead positions
      real(dp) US(NT,3) ! Unit tangent
      real(dp) U(NT,3) ! Unit tangent
      real(dp) DRDT(NT,3,4) ! Change rate of beads
      real(dp) DUDT(NT,3,4) ! Change rate of beads
      real(dp) DT       ! Time step size
      integer N,NT              ! Bead numbers
      integer RK                ! RK number
      integer I,J,IB         ! Index number
      real(dp) MAGU

      IB = 1
      do 10 I = 1,NP
         do 20 J = 1,N
            if(RK == 1) then
               R(IB,1) = RS(IB,1) + DT*DRDT(IB,1,RK)/2.
               R(IB,2) = RS(IB,2) + DT*DRDT(IB,2,RK)/2.
               R(IB,3) = RS(IB,3) + DT*DRDT(IB,3,RK)/2.
               U(IB,1) = US(IB,1) + DT*DUDT(IB,1,RK)/2.
               U(IB,2) = US(IB,2) + DT*DUDT(IB,2,RK)/2.
               U(IB,3) = US(IB,3) + DT*DUDT(IB,3,RK)/2.
            elseif(RK == 2) then
               R(IB,1) = RS(IB,1) + DT*DRDT(IB,1,RK)/2.
               R(IB,2) = RS(IB,2) + DT*DRDT(IB,2,RK)/2.
               R(IB,3) = RS(IB,3) + DT*DRDT(IB,3,RK)/2.
               U(IB,1) = US(IB,1) + DT*DUDT(IB,1,RK)/2.
               U(IB,2) = US(IB,2) + DT*DUDT(IB,2,RK)/2.
               U(IB,3) = US(IB,3) + DT*DUDT(IB,3,RK)/2.
            elseif(RK == 3) then
               R(IB,1) = RS(IB,1) + DT*DRDT(IB,1,RK)
               R(IB,2) = RS(IB,2) + DT*DRDT(IB,2,RK)
               R(IB,3) = RS(IB,3) + DT*DRDT(IB,3,RK)
               U(IB,1) = US(IB,1) + DT*DUDT(IB,1,RK)
               U(IB,2) = US(IB,2) + DT*DUDT(IB,2,RK)
               U(IB,3) = US(IB,3) + DT*DUDT(IB,3,RK)
            elseif(RK == 4) then
               R(IB,1) = RS(IB,1) + DT*(DRDT(IB,1,1)/6. + DRDT(IB,1,2)/3. + DRDT(IB,1,3)/3. + DRDT(IB,1,4)/6.)
               R(IB,2) = RS(IB,2) + DT*(DRDT(IB,2,1)/6. + DRDT(IB,2,2)/3. + DRDT(IB,2,3)/3. + DRDT(IB,2,4)/6.)
               R(IB,3) = RS(IB,3) + DT*(DRDT(IB,3,1)/6. + DRDT(IB,3,2)/3. + DRDT(IB,3,3)/3. + DRDT(IB,3,4)/6.)
               U(IB,1) = US(IB,1) + DT*(DUDT(IB,1,1)/6. + DUDT(IB,1,2)/3. + DUDT(IB,1,3)/3. + DUDT(IB,1,4)/6.)
               U(IB,2) = US(IB,2) + DT*(DUDT(IB,2,1)/6. + DUDT(IB,2,2)/3. + DUDT(IB,2,3)/3. + DUDT(IB,2,4)/6.)
               U(IB,3) = US(IB,3) + DT*(DUDT(IB,3,1)/6. + DUDT(IB,3,2)/3. + DUDT(IB,3,3)/3. + DUDT(IB,3,4)/6.)
            endif

            MAGU = sqrt(U(IB,1)**2. + U(IB,2)**2. + U(IB,3)**2.)
            U(IB,1) = U(IB,1)/MAGU
            U(IB,2) = U(IB,2)/MAGU
            U(IB,3) = U(IB,3)/MAGU

            IB = IB + 1
 20      continue
 10   continue

      RETURN
      END

!---------------------------------------------------------------*

