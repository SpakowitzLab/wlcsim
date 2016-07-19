!---------------------------------------------------------------*

      PROGRAM wlcsim

!
!     WLC Simulation Package:
!     Simulation Package for Brownian dynamics and
!     Monte Carlo Simulation
!
!     Andrew Spakowitz
!     Version 1.0
!     8/17/2015
!

!     Variables within the simulation

      use mt19937, only : grnd, sgrnd, rnorm, mt, mti

      implicit none

      external WRITE_COLTIMES ! must declare sighandlers as external

      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:):: R     ! Conformation of polymer chains
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:):: U     ! Conformation of polymer chains
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:):: R0 ! Conformation of polymer chains
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:):: U0 ! Conformation of polymer chains

      INTEGER NT                 ! Number of beads in simulation
      INTEGER N                 ! Number of beads in simulation
      INTEGER NP                ! Number of polymers in simulation
      DOUBLE PRECISION L0       ! Equilibrium segment length
      DOUBLE PRECISION ENERGY   ! Total energy
      DOUBLE PRECISION TIME     ! Current time
      DOUBLE PRECISION TSAVE     ! Time of save point
      DOUBLE PRECISION T0,TF    ! Initial/final times
      DOUBLE PRECISION DT       ! Time step size
      INTEGER I,J,IB            ! Index
      INTEGER INDMAX            ! Maximum index in series
      INTEGER IND               ! Ind in series
      INTEGER TENS              ! Decimal of index
      character*4 fileind       ! Index of output
      character*16 snapnm       ! File for output

!     Simulation input variables

      INTEGER FRMFILE           ! Initial condition
      INTEGER BROWN             ! Include Brownian forces
      INTEGER INTON             ! Include polymer interactions
      INTEGER LOGTIME           ! Is data recorded in log time?
      DOUBLE PRECISION DT0      ! Initial time step size
      INTEGER NSTEP,NINIT

!     Monte Carlo variables

      DOUBLE PRECISION MCAMP(6) ! Amplitude of random change
      INTEGER MOVEON(6)            ! Is the move active
      INTEGER WINDOW(6)            ! Size of window for bead selection
      INTEGER SUCCESS(6)        ! Number of successes

!     Energy variables

      DOUBLE PRECISION EELAS(3) ! Elastic energy
      DOUBLE PRECISION EPONP    ! Poly-poly energy

!     Structure analysis

      DOUBLE PRECISION RCOM(3)  ! Center of mass
      DOUBLE PRECISION DELR(3)  ! Mag of gyration tensor
      DOUBLE PRECISION RCOM0(3) ! Init val RCOM
      DOUBLE PRECISION DELR0(3) ! Init val DELR
      DOUBLE PRECISION DRCOM    ! Change in RCOM
      DOUBLE PRECISION SIG(3,3)
      DOUBLE PRECISION COR

!     Variables in the simulation

      DOUBLE PRECISION PARA(10)
      INTEGER SIMTYPE           ! Simulation method (WLC=1,SSWLC=2,GC=3)

!     Variables for the random number generators

      INTEGER IDUM              ! Seed for the generator
      DOUBLE PRECISION MOM(6)

!     Variable to hold time of first collisions between each bead
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:):: HAS_COLLIDED
      DOUBLE PRECISION FPT_DIST ! l1 dist to trigger collision
      INTEGER COL_TYPE ! what kind of collision checking to use


!     Load in the parameters for the simulation

      open (unit=5, file='input/input')
      read (unit=5, fmt='(24(/))')
      read (unit=5, fmt=*) N
      read (unit=5, fmt='(2(/))')
      read (unit=5, fmt=*) NP
      read (unit=5, fmt='(2(/))')
      read (unit=5, fmt=*) TF
      read (unit=5, fmt='(2(/))')
      read (unit=5, fmt=*) INDMAX
      read (unit=5, fmt='(2(/))')
      read (unit=5, fmt=*) DT
      read (unit=5, fmt='(2(/))')
      read (unit=5, fmt=*) FRMFILE
      read (unit=5, fmt='(2(/))')
      read (unit=5, fmt=*) BROWN
      read (unit=5, fmt='(2(/))')
      read (unit=5, fmt=*) INTON
      read (unit=5, fmt='(2(/))')
      read (unit=5, fmt=*) LOGTIME
      read (unit=5, fmt='(2(/))')
      read (unit=5, fmt=*) NINIT
      read (unit=5, fmt='(2(/))')
      read (unit=5, fmt=*) NSTEP
      read (unit=5, fmt='(2(/))')
      read (unit=5, fmt=*) FPT_DIST
      read (unit=5, fmt='(2(/))')
      read (unit=5, fmt=*) COL_TYPE
      close(5)
      call getpara(PARA,DT,SIMTYPE)
      DT0=DT

      NT=N*NP
      ALLOCATE(R(NT,3))
      ALLOCATE(U(NT,3))
      ALLOCATE(R0(NT,3))
      ALLOCATE(U0(NT,3))
      if (COL_TYPE.NE.0) then
         ALLOCATE(HAS_COLLIDED(NT,NT))
         HAS_COLLIDED = -1.0d+0
      endif

!     Setup the initial condition

      call initcond(R,U,NT,N,NP,IDUM,FRMFILE,PARA)

!     Turn on moves for each simulation type

      if (SIMTYPE.EQ.1) then
         MCAMP(1)=1.
         MCAMP(2)=1.
         MCAMP(3)=1.
         MCAMP(4)=1.
         MCAMP(5)=1.
         MCAMP(6)=1.
         MOVEON(1)=1
         MOVEON(2)=0
         MOVEON(3)=1
         MOVEON(4)=0
      elseif (SIMTYPE.EQ.2) then
         MCAMP(1)=1.
         MCAMP(2)=1.
         MCAMP(3)=1.
         MCAMP(4)=1.
         MCAMP(5)=1.
         MCAMP(6)=1.
         MOVEON(1)=1
         MOVEON(2)=1
         MOVEON(3)=1
         MOVEON(4)=1
      elseif (SIMTYPE.EQ.3) then
         MCAMP(1)=1.
         MCAMP(2)=1.
         MCAMP(3)=1.
         MCAMP(4)=1.
         MCAMP(5)=1.
         MCAMP(6)=1.
         MOVEON(1)=1
         MOVEON(2)=1
         MOVEON(3)=1
         MOVEON(4)=0
      endif

!     Turn off whole chain rotation and translation if interactions are off

      if (INTON.EQ.1) then
         MOVEON(5)=1
         MOVEON(6)=1
      else
         MOVEON(5)=0
         MOVEON(6)=0
      endif

!     Perform an initialization MC simulation

      call MCsim(R,U,NT,N,NP,NINIT,BROWN,INTON,IDUM,PARA,MCAMP, &
           SUCCESS,MOVEON,WINDOW,SIMTYPE)

!     Save the conformation and PSI angles

      OPEN (UNIT = 1, FILE = 'data/r0', STATUS = 'NEW')
      IB=1
      DO 10 I=1,NP
         DO 20 J=1,N
            R0(IB,1)=R(IB,1)
            R0(IB,2)=R(IB,2)
            R0(IB,3)=R(IB,3)
            U0(IB,1)=U(IB,1)
            U0(IB,2)=U(IB,2)
            U0(IB,3)=U(IB,3)
            WRITE(1,*) R(IB,1),R(IB,2),R(IB,3)
            IB=IB+1
 20      CONTINUE
 10   CONTINUE
      CLOSE(1)

      OPEN (UNIT = 1, FILE = 'data/u0', STATUS = 'NEW')
      IB=1
      DO 30 I=1,NP
         DO 40 J=1,N
            WRITE(1,*) U(IB,1),U(IB,2),U(IB,3)
            IB=IB+1
 40      CONTINUE
 30   CONTINUE
      CLOSE(1)

!     Begin simulation

      IND=1
      TIME=0.

!     Open the output files

      OPEN (UNIT = 2, FILE = 'data/out1', STATUS = 'NEW')
      OPEN (UNIT = 3, FILE = 'data/out2', STATUS = 'NEW')
      OPEN (UNIT = 4, FILE = 'data/out3', STATUS = 'NEW')

      call stress(SIG,R,U,NT,N,NP,PARA,INTON)

      WRITE(3,*) real(SIG(1,1)),real(SIG(1,2)),real(SIG(1,3)),real(SIG(2,1)),real(SIG(2,2))
      WRITE(4,*) real(SIG(2,3)),real(SIG(3,1)),real(SIG(3,2)),real(SIG(3,3))

      DO WHILE (IND.LE.INDMAX)

!     Perform a MC simulation, only if NSTEP.NE.0

         call MCsim(R,U,NT,N,NP,NSTEP,BROWN,INTON,IDUM,PARA,MCAMP, &
              SUCCESS,MOVEON,WINDOW,SIMTYPE)

!     Perform a Brownian dynamics simulation over time step

         if (LOGTIME.EQ.0) then
            TSAVE = TF*IND/INDMAX
         else
            TSAVE = DT0*exp((IND-1.)/(INDMAX-1.)*log(TF/DT0))
         endif
         if (NSTEP.EQ.0) then
            call BDsim(R,U,NT,N,NP,TIME,TSAVE,DT,BROWN,INTON,IDUM, &
                       PARA,SIMTYPE,HAS_COLLIDED,FPT_DIST,COL_TYPE)
         endif

!     Save the conformation and the metrics

         TENS=nint(log10(1.*IND)-0.4999)+1
         write (fileind,'(I4)'), IND
         snapnm= 'data/r'//fileind((4-TENS+1):4)
         OPEN (UNIT = 1, FILE = snapnm, STATUS = 'NEW')
         IB=1
         DO 50 I=1,NP
            DO 60 J=1,N
               WRITE(1,*) R(IB,1),R(IB,2),R(IB,3)
               IB=IB+1
 60         CONTINUE
 50      CONTINUE
         CLOSE(1)

         snapnm= 'data/u'//fileind((4-TENS+1):4)
         OPEN (UNIT = 1, FILE = snapnm, STATUS = 'NEW')
         IB=1
         DO 70 I=1,NP
            DO 80 J=1,N
               WRITE(1,*) U(IB,1),U(IB,2),U(IB,3)
               IB=IB+1
 80         CONTINUE
 70      CONTINUE
         CLOSE(1)

         snapnm='data/coltimes'
         OPEN (UNIT=1, FILE=snapnm, STATUS='REPLACE')
         DO, I=1,NT
             WRITE(1,*) ( HAS_COLLIDED(i,j), j=1,NT )
         ENDDO
         CLOSE(1)

         call stress(SIG,R,U,NT,N,NP,PARA,INTON,SIMTYPE)
         call stressp(COR,R,U,R0,U0,NT,N,NP,PARA,INTON,SIMTYPE)

         call energy_elas(EELAS,R,U,NT,N,NP,PARA)
         EPONP=0.
         if (INTON.EQ.1) then
            call energy_ponp(EPONP,R,NT,N,NP,PARA)
         endif
         WRITE(2,*) real(TIME),real(EELAS(1)),real(EELAS(2)),real(EELAS(3)),real(EPONP),real(COR)


         WRITE(3,*) real(SIG(1,1)),real(SIG(1,2)),real(SIG(1,3)),real(SIG(2,1)),real(SIG(2,2))
         WRITE(4,*) real(SIG(2,3)),real(SIG(3,1)),real(SIG(3,2)),real(SIG(3,3))

         PRINT*, '________________________________________'
         PRINT*, 'Time point ',IND, ' out of', INDMAX
         PRINT*, 'Current time ',TIME
         PRINT*, 'Bending energy ', EELAS(1)
         PRINT*, 'Par compression energy ', EELAS(2)
         PRINT*, 'Perp compression energy ', EELAS(3)
         PRINT*, 'Polymer-polymer energy ', EPONP
         PRINT*, 'Current number of beads ', N
         PRINT*, 'Time step ', DT
         print*, 'End-to-end distance poly 1 ', &
              sqrt((R(N,1)-R(1,1))**2.+(R(N,2)-R(1,2))**2.+(R(N,3)-R(1,3))**2.)
         PRINT*, 'Simulation type ', SIMTYPE

         IND=IND+1

      ENDDO
      END


!---------------------------------------------------------------*

