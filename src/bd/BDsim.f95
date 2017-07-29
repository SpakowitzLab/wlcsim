!---------------------------------------------------------------*

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

      real(dp) R(NT,3)  ! Bead positions
      real(dp) U(NT,3)  ! Tangent vectors
      real(dp) TIME     ! Time of BD simulation
      real(dp) TTOT     ! Final time of BD simulation
      integer N,NP,NT           ! Number of beads

!     Variables in the simulation

      real(dp) RS(NT,3) ! R during the step
      real(dp) US(NT,3) ! R during the step
      real(dp) DT,DT0       ! Time step size
      integer RK                ! Runge-Kutta index
      real(dp) DRDT(NT,3,4) ! Position rate of change
      real(dp) DUDT(NT,3,4) ! Position rate of change
      integer I,J,IB            ! Index Holders
      real(dp) doTU

!     Variables for use in the force calculations

      real(dp) FELAS(NT,3) ! Elastic force
      real(dp) FPONP(NT,3) ! self-int force
      real(dp) TELAS(NT,3) ! Elastic force
      real(dp) TPONP(NT,3) ! self-int force
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

      real(dp) FRAND(NT,3) ! Random force
      real(dp) TRAND(NT,3) ! Random force
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
            RS(IB,1) = R(IB,1)
            RS(IB,2) = R(IB,2)
            RS(IB,3) = R(IB,3)
            US(IB,1) = U(IB,1)
            US(IB,2) = U(IB,2)
            US(IB,3) = U(IB,3)
            FELAS(IB,1) = 0.
            FELAS(IB,2) = 0.
            FELAS(IB,3) = 0.
            FRAND(IB,1) = 0.
            FRAND(IB,2) = 0.
            FRAND(IB,3) = 0.
            FPONP(IB,1) = 0.
            FPONP(IB,2) = 0.
            FPONP(IB,3) = 0.
            TELAS(IB,1) = 0.
            TELAS(IB,2) = 0.
            TELAS(IB,3) = 0.
            TRAND(IB,1) = 0.
            TRAND(IB,2) = 0.
            TRAND(IB,3) = 0.
            TPONP(IB,1) = 0.
            TPONP(IB,2) = 0.
            TPONP(IB,3) = 0.
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
                     FRAND(IB,1) = MAGR*rnorm()
                     FRAND(IB,2) = MAGR*rnorm()
                     FRAND(IB,3) = MAGR*rnorm()
                     TRAND(IB,1) = MAGU*rnorm()
                     TRAND(IB,2) = MAGU*rnorm()
                     TRAND(IB,3) = MAGU*rnorm()
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
                        R(IB,1) = RS(IB,1)
                        R(IB,2) = RS(IB,2)
                        R(IB,3) = RS(IB,3)
                        U(IB,1) = US(IB,1)
                        U(IB,2) = US(IB,2)
                        U(IB,3) = US(IB,3)
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
                  DRDT(IB,1,RK) = (FELAS(IB,1) + FPONP(IB,1))/XIR
                  DRDT(IB,2,RK) = (FELAS(IB,2) + FPONP(IB,2))/XIR
                  DRDT(IB,3,RK) = (FELAS(IB,3) + FPONP(IB,3))/XIR
                  DUDT(IB,1,RK) = (TELAS(IB,1) + TPONP(IB,1))/XIU
                  DUDT(IB,2,RK) = (TELAS(IB,2) + TPONP(IB,2))/XIU
                  DUDT(IB,3,RK) = (TELAS(IB,3) + TPONP(IB,3))/XIU

                  if (BROWN == 0) then
                     doTU = DUDT(IB,1,RK)*U(IB,1) + DUDT(IB,2,RK)*U(IB,2) + DUDT(IB,3,RK)*U(IB,3)
                     DUDT(IB,1,RK) = DUDT(IB,1,RK)-doTU*U(IB,1)
                     DUDT(IB,2,RK) = DUDT(IB,2,RK)-doTU*U(IB,2)
                     DUDT(IB,3,RK) = DUDT(IB,3,RK)-doTU*U(IB,3)
                  endif
                  IB = IB + 1
 80            continue
 70         continue

            if (BROWN == 1) then
               IB = 1
               do 90 I = 1,NP
                  do 100 J = 1,N
                     DRDT(IB,1,RK) = DRDT(IB,1,RK) + FRAND(IB,1)/XIR
                     DRDT(IB,2,RK) = DRDT(IB,2,RK) + FRAND(IB,2)/XIR
                     DRDT(IB,3,RK) = DRDT(IB,3,RK) + FRAND(IB,3)/XIR
                     DUDT(IB,1,RK) = DUDT(IB,1,RK) + TRAND(IB,1)/XIU
                     DUDT(IB,2,RK) = DUDT(IB,2,RK) + TRAND(IB,2)/XIU
                     DUDT(IB,3,RK) = DUDT(IB,3,RK) + TRAND(IB,3)/XIU

                     doTU = DUDT(IB,1,RK)*U(IB,1) + DUDT(IB,2,RK)*U(IB,2) + DUDT(IB,3,RK)*U(IB,3)
                     DUDT(IB,1,RK) = DUDT(IB,1,RK)-doTU*U(IB,1)
                     DUDT(IB,2,RK) = DUDT(IB,2,RK)-doTU*U(IB,2)
                     DUDT(IB,3,RK) = DUDT(IB,3,RK)-doTU*U(IB,3)

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
               RS(IB,1) = R(IB,1)
               RS(IB,2) = R(IB,2)
               RS(IB,3) = R(IB,3)
               US(IB,1) = U(IB,1)
               US(IB,2) = U(IB,2)
               US(IB,3) = U(IB,3)
               IB = IB + 1
 120        continue
 110     continue

      ENDdo

      RETURN
      END

!---------------------------------------------------------------*
