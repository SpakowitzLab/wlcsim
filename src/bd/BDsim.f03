!--------------------------------------------------------------*

      subroutine BDsim(R,U,NT,N,NP,TIME,TTOT,DT,BROWN, &
           inTON,IDUM,PARA,SIMTYPE,COLLISION_TIME,COL_DIST, &
           COL_TYPE)
!
!     External subroutine to perform a Brownian dynamics simulation.
!
!     Andrew Spakowitz
!     Written 11-11-13

      use mt19937, only : rnorm
      use params, only : dp
      implicit none
      real(dp) R(3,NT)  ! Bead positions
      real(dp) U(3,NT)  ! Tangent vectors
      real(dp) TIME     ! Time of BD simulation
      real(dp) TTOT     ! Final time of BD simulation
      integer N,NP,NT           ! Number of beads

!     Variables in the simulation

      real(dp) RS(3,NT) ! R during the step
      real(dp) US(3,NT) ! R during the step
      real(dp) DT,DT0       ! Time step size
      integer RK                ! Runge-Kutta index
      real(dp) DRDT(3,NT,4) ! Position rate of change
      real(dp) DUDT(3,NT,4) ! Position rate of change
      integer I,J,IB            ! Index Holders
      real(dp) doTU

!     Variables for use in the force calculations

      real(dp) FELAS(3,NT) ! Elastic force
      real(dp) FPONP(3,NT) ! self-int force
      real(dp) TELAS(3,NT) ! Elastic force
      real(dp) TPONP(3,NT) ! self-int force
      real(dp) FORCE    ! External force

!     Variables in the simulation

      real(dp) EB,EPAR,EPERP
      real(dp) GAM,ETA
      real(dp) XIR,XIU
      real(dp) LBOX     ! Box edge length
      real(dp) LHC      ! Length of HC int
      real(dp) VHC      ! HC strength
      real(dp) PARA(10)
      integer SIMTYPE           ! Simulation method (WLC = 1,SSWLC = 2,GC = 3)

!     Variables used for the Brownian forces

      real(dp) FRAND(3,NT) ! Random force
      real(dp) TRAND(3,NT) ! Random force
      real(dp) MAGR,MAGU ! Mag of Brownian forces
      integer BROWN             ! Logic for BD forces
      integer inTON             ! Include polymer interactions
      integer IDUM              ! Seed for the generator

!     Variables for the timestep switch

      integer SWDT

!     Variable to hold time of first collisions between each bead
      real(dp) COLLISION_TIME(NT,NT)
      real(dp) COL_DIST ! l1 dist to trigger collision
      integer COL_TYPE ! algorithm to use for collision detection

!     Load the input parameters

      EB = PARA(1)
      EPAR = PARA(2)
      EPERP = PARA(3)
      GAM = PARA(4)
      ETA = PARA(5)
      XIR = PARA(6)
      XIU = PARA(7)
      LBOX = PARA(8)
      LHC = PARA(9)
      VHC = PARA(10)

      MAGR = sqrt(XIR*2.0/DT)
      MAGU = sqrt(XIU*2.0/DT)
      DT0 = DT
      SWDT = 0

!     Setup the geometric parameters and initialize random forces

      IB = 1
      do 10 I = 1,NP
         do 20 J = 1,N
            RS(1,IB) = R(1,IB)
            RS(2,IB) = R(2,IB)
            RS(3,IB) = R(3,IB)
            US(1,IB) = U(1,IB)
            US(2,IB) = U(2,IB)
            US(3,IB) = U(3,IB)
            FELAS(1,IB) = 0.
            FELAS(2,IB) = 0.
            FELAS(3,IB) = 0.
            FRAND(1,IB) = 0.
            FRAND(2,IB) = 0.
            FRAND(3,IB) = 0.
            FPONP(1,IB) = 0.
            FPONP(2,IB) = 0.
            FPONP(3,IB) = 0.
            TELAS(1,IB) = 0.
            TELAS(2,IB) = 0.
            TELAS(3,IB) = 0.
            TRAND(1,IB) = 0.
            TRAND(2,IB) = 0.
            TRAND(3,IB) = 0.
            TPONP(1,IB) = 0.
            TPONP(2,IB) = 0.
            TPONP(3,IB) = 0.
            IB = IB + 1
 20      continue
 10   continue

!     Begin the time integration

      do while (TIME < TTOT)

         call CHECK_COLLISIONS(R, NT, COLLISION_TIME, COL_DIST, TIME, COL_TYPE)

!     Calculate the random forces and torques for use in this
!     timestep calculation if BROWN = 1

         RK = 1
         do while (RK <= 4)

 130        continue

            if (BROWN == 1.AND.RK == 1) then
               IB = 1
               do 30 I = 1,NP
                  do 40 J = 1,N
                     FRAND(1,IB) = MAGR*rnorm()
                     FRAND(2,IB) = MAGR*rnorm()
                     FRAND(3,IB) = MAGR*rnorm()
                     TRAND(1,IB) = MAGU*rnorm()
                     TRAND(2,IB) = MAGU*rnorm()
                     TRAND(3,IB) = MAGU*rnorm()
                     IB = IB + 1
 40               continue
 30            continue
            endif

!     Calculate the four Runge-Kutta derivatives


!     Calculate the elastic forces (same as free chain)

            call force_elas(FELAS,TELAS,R,U,NT,N,NP,EB,EPAR,EPERP,GAM,ETA,SIMTYPE)

!     Calculate the self forces

            if (inTON == 1) then
               call force_ponp(FPONP,R,NT,N,NP,LHC,VHC,LBOX,GAM,DT,XIR,SWDT)

!     If timestep is switch, reset coords and redo step

               if (SWDT == 1) then
                  print*, "Time-step switch", DT,RK,TIME
                  SWDT = 0
                  MAGR = sqrt(XIR*2.0/DT)
                  MAGU = sqrt(XIU*2.0/DT)
                  IB = 1
                  do 60 I = 1,NP
                     do 65 J = 1,N
                        R(1,IB) = RS(1,IB)
                        R(2,IB) = RS(2,IB)
                        R(3,IB) = RS(3,IB)
                        U(1,IB) = US(1,IB)
                        U(2,IB) = US(2,IB)
                        U(3,IB) = US(3,IB)
                        IB = IB + 1
 65                  continue
 60               continue
                  RK = 1
                  goto 130
               endif
            endif


!     Calculate the change in the position vector

            IB = 1
            do 70 I = 1,NP
               do 80 J = 1,N
                  DRDT(1,IB,RK) = (FELAS(1,IB) + FPONP(1,IB))/XIR
                  DRDT(2,IB,RK) = (FELAS(2,IB) + FPONP(2,IB))/XIR
                  DRDT(3,IB,RK) = (FELAS(3,IB) + FPONP(3,IB))/XIR
                  DUDT(1,IB,RK) = (TELAS(1,IB) + TPONP(1,IB))/XIU
                  DUDT(2,IB,RK) = (TELAS(2,IB) + TPONP(2,IB))/XIU
                  DUDT(3,IB,RK) = (TELAS(3,IB) + TPONP(3,IB))/XIU

                  if (BROWN == 0) then
                     doTU = DUDT(1,IB,RK)*U(1,IB) + DUDT(2,IB,RK)*U(2,IB) + DUDT(3,IB,RK)*U(3,IB)
                     DUDT(1,IB,RK) = DUDT(1,IB,RK)-doTU*U(1,IB)
                     DUDT(2,IB,RK) = DUDT(2,IB,RK)-doTU*U(2,IB)
                     DUDT(3,IB,RK) = DUDT(3,IB,RK)-doTU*U(3,IB)
                  endif
                  IB = IB + 1
 80            continue
 70         continue

            if (BROWN == 1) then
               IB = 1
               do 90 I = 1,NP
                  do 100 J = 1,N
                     DRDT(1,IB,RK) = DRDT(1,IB,RK) + FRAND(1,IB)/XIR
                     DRDT(2,IB,RK) = DRDT(2,IB,RK) + FRAND(2,IB)/XIR
                     DRDT(3,IB,RK) = DRDT(3,IB,RK) + FRAND(3,IB)/XIR
                     DUDT(1,IB,RK) = DUDT(1,IB,RK) + TRAND(1,IB)/XIU
                     DUDT(2,IB,RK) = DUDT(2,IB,RK) + TRAND(2,IB)/XIU
                     DUDT(3,IB,RK) = DUDT(3,IB,RK) + TRAND(3,IB)/XIU

                     doTU = DUDT(1,IB,RK)*U(1,IB) + DUDT(2,IB,RK)*U(2,IB) + DUDT(3,IB,RK)*U(3,IB)
                     DUDT(1,IB,RK) = DUDT(1,IB,RK)-doTU*U(1,IB)
                     DUDT(2,IB,RK) = DUDT(2,IB,RK)-doTU*U(2,IB)
                     DUDT(3,IB,RK) = DUDT(3,IB,RK)-doTU*U(3,IB)

                     IB = IB + 1
 100               continue
 90            continue
            endif

!     If SIMTYPE = 1 (WLC), calculate the constraint forces

            if (SIMTYPE == 1) then
               call concalc(R,DRDT,NT,N,NP,XIR,GAM,DT,RK,BROWN)
            endif

!     Step forward using the RK algorithm

            call RKstep(RS,R,US,U,DRDT,DUDT,NT,N,NP,RK,DT)

            RK = RK + 1

         ENDdo

         TIME = TIME + DT

!     Swap old variables for new ones

         DT = DT0
         MAGR = sqrt(XIR*2.0/DT)
         MAGU = sqrt(XIU*2.0/DT)

         IB = 1
         do 110 I = 1,NP
            do 120 J = 1,N
               RS(1,IB) = R(1,IB)
               RS(2,IB) = R(2,IB)
               RS(3,IB) = R(3,IB)
               US(1,IB) = U(1,IB)
               US(2,IB) = U(2,IB)
               US(3,IB) = U(3,IB)
               IB = IB + 1
 120        continue
 110     continue

      ENDdo

      RETURN
      END

!---------------------------------------------------------------*
