!---------------------------------------------------------------*

      SUBROUTINE BDsim(R,U,NT,N,NP,TIME,TTOT,DT,BROWN, &
           INTON,IDUM,PARA,SIMTYPE,HAS_COLLIDED,FPT_DIST, &
           COL_TYPE)

!
!     External subroutine to perform a Brownian dynamics simulation.
!
!     Andrew Spakowitz
!     Written 11-11-13

      use mt19937, only : rnorm

      DOUBLE PRECISION R(NT,3)  ! Bead positions
      DOUBLE PRECISION U(NT,3)  ! Tangent vectors
      DOUBLE PRECISION TIME     ! Time of BD simulation
      DOUBLE PRECISION TTOT     ! Final time of BD simulation
      INTEGER N,NP,NT           ! Number of beads

!     Variables in the simulation

      DOUBLE PRECISION B(NT-1)  ! Bond length
      DOUBLE PRECISION RS(NT,3) ! R during the step
      DOUBLE PRECISION US(NT,3) ! R during the step
      DOUBLE PRECISION L0       ! Bond distances
      DOUBLE PRECISION DT,DT0       ! Time step size
      INTEGER RK                ! Runge-Kutta index
      DOUBLE PRECISION DRDT(NT,3,4) ! Position rate of change
      DOUBLE PRECISION DUDT(NT,3,4) ! Position rate of change
      INTEGER I,J,IB            ! Index Holders
      DOUBLE PRECISION DOTU
      DOUBLE PRECISION R0(3)

!     Variables for use in the force calculations

      DOUBLE PRECISION FELAS(NT,3) ! Elastic force
      DOUBLE PRECISION FPONP(NT,3) ! self-int force
      DOUBLE PRECISION TELAS(NT,3) ! Elastic force
      DOUBLE PRECISION TPONP(NT,3) ! self-int force
      DOUBLE PRECISION FORCE    ! External force
      INTEGER FON               ! Is force on?

!     Variables in the simulation

      DOUBLE PRECISION EB,EPAR,EPERP
      DOUBLE PRECISION GAM,ETA
      DOUBLE PRECISION XIR,XIU
      DOUBLE PRECISION LBOX     ! Box edge length
      DOUBLE PRECISION LHC      ! Length of HC int
      DOUBLE PRECISION VHC      ! HC strength
      DOUBLE PRECISION PARA(10)
      INTEGER SIMTYPE           ! Simulation method (WLC=1,SSWLC=2,GC=3)

!     Variables used for the Brownian forces

      DOUBLE PRECISION FRAND(NT,3) ! Random force
      DOUBLE PRECISION TRAND(NT,3) ! Random force
      DOUBLE PRECISION MAGR,MAGU ! Mag of Brownian forces
      INTEGER BROWN             ! Logic for BD forces
      INTEGER INTON             ! Include polymer interactions
      REAL ran1                 ! Random number generator
      INTEGER IDUM              ! Seed for the generator
      INTEGER NOW(3)            ! Time now (hr,min,sec)

!     Variables for the timestep switch

      INTEGER SWDT

!     Variable to hold time of first collisions between each bead
      DOUBLE PRECISION HAS_COLLIDED(NT,NT)
      DOUBLE PRECISION FPT_DIST ! l1 dist to trigger collision
      INTEGER COL_TYPE ! algorithm to use for collision detection

!     Load the input parameters

      EB=PARA(1)
      EPAR=PARA(2)
      EPERP=PARA(3)
      GAM=PARA(4)
      ETA=PARA(5)
      XIR=PARA(6)
      XIU=PARA(7)
      LBOX=PARA(8)
      LHC=PARA(9)
      VHC=PARA(10)

      MAGR=sqrt(XIR*2.0/DT)
      MAGU=sqrt(XIU*2.0/DT)
      DT0=DT
      SWDT=0

!     Setup the geometric parameters and initialize random forces

      IB=1
      DO 10 I=1,NP
         DO 20 J=1,N
            RS(IB,1)=R(IB,1)
            RS(IB,2)=R(IB,2)
            RS(IB,3)=R(IB,3)
            US(IB,1)=U(IB,1)
            US(IB,2)=U(IB,2)
            US(IB,3)=U(IB,3)
            FELAS(IB,1)=0.
            FELAS(IB,2)=0.
            FELAS(IB,3)=0.
            FRAND(IB,1)=0.
            FRAND(IB,2)=0.
            FRAND(IB,3)=0.
            FPONP(IB,1)=0.
            FPONP(IB,2)=0.
            FPONP(IB,3)=0.
            TELAS(IB,1)=0.
            TELAS(IB,2)=0.
            TELAS(IB,3)=0.
            TRAND(IB,1)=0.
            TRAND(IB,2)=0.
            TRAND(IB,3)=0.
            TPONP(IB,1)=0.
            TPONP(IB,2)=0.
            TPONP(IB,3)=0.
            IB=IB+1
 20      CONTINUE
 10   CONTINUE

!     Begin the time integration

      DO WHILE (TIME.LT.TTOT)

         call CHECK_COLLISIONS(R, NT, HAS_COLLIDED, FPT_DIST, TIME, COL_TYPE)


!     Calculate the random forces and torques for use in this
!     timestep calculation if BROWN=1

         RK=1
         DO WHILE (RK.LE.4)

 130        CONTINUE

            if (BROWN.EQ.1.AND.RK.EQ.1) then
               IB=1
               DO 30 I=1,NP
                  DO 40 J=1,N
                     FRAND(IB,1)=MAGR*rnorm()
                     FRAND(IB,2)=MAGR*rnorm()
                     FRAND(IB,3)=MAGR*rnorm()
                     TRAND(IB,1)=MAGU*rnorm()
                     TRAND(IB,2)=MAGU*rnorm()
                     TRAND(IB,3)=MAGU*rnorm()
                     IB=IB+1
 40               CONTINUE
 30            CONTINUE
            endif

!     Calculate the four Runge-Kutta derivatives


!     Calculate the elastic forces (same as free chain)

            call force_elas(FELAS,TELAS,R,U,NT,N,NP,EB,EPAR,EPERP,GAM,ETA,SIMTYPE)

!     Calculate the self forces

            if (INTON.EQ.1) then
               call force_ponp(FPONP,R,NT,N,NP,LHC,VHC,LBOX,GAM,DT,XIR,SWDT)

!     If timestep is switch, reset coords and redo step

               if (SWDT.EQ.1) then
                  print*, "Time-step switch", DT,RK,TIME
                  SWDT=0
                  MAGR=sqrt(XIR*2.0/DT)
                  MAGU=sqrt(XIU*2.0/DT)
                  IB=1
                  DO 60 I=1,NP
                     DO 65 J=1,N
                        R(IB,1)=RS(IB,1)
                        R(IB,2)=RS(IB,2)
                        R(IB,3)=RS(IB,3)
                        U(IB,1)=US(IB,1)
                        U(IB,2)=US(IB,2)
                        U(IB,3)=US(IB,3)
                        IB=IB+1
 65                  CONTINUE
 60               CONTINUE
                  RK=1
                  goto 130
               endif
            endif


!     Calculate the change in the position vector

            IB=1
            DO 70 I=1,NP
               DO 80 J=1,N
                  DRDT(IB,1,RK)=(FELAS(IB,1)+FPONP(IB,1))/XIR
                  DRDT(IB,2,RK)=(FELAS(IB,2)+FPONP(IB,2))/XIR
                  DRDT(IB,3,RK)=(FELAS(IB,3)+FPONP(IB,3))/XIR
                  DUDT(IB,1,RK)=(TELAS(IB,1)+TPONP(IB,1))/XIU
                  DUDT(IB,2,RK)=(TELAS(IB,2)+TPONP(IB,2))/XIU
                  DUDT(IB,3,RK)=(TELAS(IB,3)+TPONP(IB,3))/XIU

                  if (BROWN.EQ.0) then
                     DOTU=DUDT(IB,1,RK)*U(IB,1)+DUDT(IB,2,RK)*U(IB,2)+DUDT(IB,3,RK)*U(IB,3)
                     DUDT(IB,1,RK)=DUDT(IB,1,RK)-DOTU*U(IB,1)
                     DUDT(IB,2,RK)=DUDT(IB,2,RK)-DOTU*U(IB,2)
                     DUDT(IB,3,RK)=DUDT(IB,3,RK)-DOTU*U(IB,3)
                  endif
                  IB=IB+1
 80            CONTINUE
 70         CONTINUE

            if (BROWN.EQ.1) then
               IB=1
               DO 90 I=1,NP
                  DO 100 J=1,N
                     DRDT(IB,1,RK)=DRDT(IB,1,RK)+FRAND(IB,1)/XIR
                     DRDT(IB,2,RK)=DRDT(IB,2,RK)+FRAND(IB,2)/XIR
                     DRDT(IB,3,RK)=DRDT(IB,3,RK)+FRAND(IB,3)/XIR
                     DUDT(IB,1,RK)=DUDT(IB,1,RK)+TRAND(IB,1)/XIU
                     DUDT(IB,2,RK)=DUDT(IB,2,RK)+TRAND(IB,2)/XIU
                     DUDT(IB,3,RK)=DUDT(IB,3,RK)+TRAND(IB,3)/XIU

                     DOTU=DUDT(IB,1,RK)*U(IB,1)+DUDT(IB,2,RK)*U(IB,2)+DUDT(IB,3,RK)*U(IB,3)
                     DUDT(IB,1,RK)=DUDT(IB,1,RK)-DOTU*U(IB,1)
                     DUDT(IB,2,RK)=DUDT(IB,2,RK)-DOTU*U(IB,2)
                     DUDT(IB,3,RK)=DUDT(IB,3,RK)-DOTU*U(IB,3)

                     IB=IB+1
 100               CONTINUE
 90            CONTINUE
            endif

!     If SIMTYPE=1 (WLC), calculate the constraint forces

            if (SIMTYPE.EQ.1) then
               call concalc(R,DRDT,NT,N,NP,XIR,GAM,DT,RK,BROWN)
            endif

!     Step forward using the RK algorithm

            call RKstep(RS,R,US,U,DRDT,DUDT,NT,N,NP,RK,DT)

            RK=RK+1

         ENDDO

         TIME=TIME+DT

!     Swap old variables for new ones

         DT=DT0
         MAGR=sqrt(XIR*2.0/DT)
         MAGU=sqrt(XIU*2.0/DT)

         IB=1
         DO 110 I=1,NP
            DO 120 J=1,N
               RS(IB,1)=R(IB,1)
               RS(IB,2)=R(IB,2)
               RS(IB,3)=R(IB,3)
               US(IB,1)=U(IB,1)
               US(IB,2)=U(IB,2)
               US(IB,3)=U(IB,3)
               IB=IB+1
 120        CONTINUE
 110     CONTINUE

      ENDDO

      RETURN
      END

!---------------------------------------------------------------*
