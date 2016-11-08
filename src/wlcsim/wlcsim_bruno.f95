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

      implicit none

      external WRITE_COLTIMES ! must declare sighandlers as external

      real(dp), allocatable, dimension(:,:):: R     ! Conformation of polymer chains
      real(dp), allocatable, dimension(:,:):: U     ! Conformation of polymer chains
      real(dp), allocatable, dimension(:,:):: R0 ! Conformation of polymer chains
      real(dp), allocatable, dimension(:,:):: U0 ! Conformation of polymer chains

      integer NT                 ! Number of beads in simulation
      integer N                 ! Number of beads in simulation
      integer NP                ! Number of polymers in simulation
      real(dp) L0       ! Equilibrium segment length
      real(dp) ENERGY   ! Total energy
      real(dp) TIME     ! Current time
      real(dp) TSAVE     ! Time of save point
      real(dp) T0,TF    ! Initial/final times
      real(dp) DT       ! Time step size
      integer I,J,IB            ! index
      integer indMAX            ! Maximum index in series
      integer ind               ! ind in series
      integer TENS              ! Decimal of index
      character*4 fileind       ! index of output
      character*16 snapnm       ! file for output

!     Simulation input variables

      integer FRMfile           ! Initial condition
      integer BROWN             ! Include Brownian forces
      integer INTON             ! Include polymer interactions
      integer LOGTIME           ! Is data recorded in log time?
      real(dp) DT0      ! Initial time step size
      integer NSTEP,NINIT

!     Monte Carlo variables

      real(dp) MCAMP(6) ! Amplitude of random change
      integer MOVEON(6)            ! Is the move active
      integer WindoW(6)            ! Size of window for bead selection
      integer SUCCESS(6)        ! Number of successes

!     Energy variables

      real(dp) EELAS(3) ! Elastic energy
      real(dp) EPONP    ! Poly-poly energy

!     Structure analysis

      real(dp) RCOM(3)  ! Center of mass
      real(dp) delR(3)  ! Mag of gyration tensor
      real(dp) RCOM0(3) ! Init val RCOM
      real(dp) delR0(3) ! Init val delR
      real(dp) DRCOM    ! Change in RCOM
      real(dp) SIG(3,3)
      real(dp) COR

!     Variables in the simulation

      real(dp) para(10)
      integer SIMtype           ! Simulation method (WLC=1,SSWLC=2,GC=3)

!     Variables for the random number generators

      integer IDUM              ! Seed for the generator
      real(dp) MOM(6)

!     Variable to hold time of first collisions between each bead
      real(dp), allocatable, dimension(:,:):: HAS_COLLIDED
      integer NUM_POSSIBLE_COLLISIONS
      real(dp) FPT_DIST ! l1 dist to trigger collision
      integer COL_type ! what kind of collision checking to use

!     Variables to control simulation output
      integer SAVE_RU
      integer EXIT_WHEN_COLLIDED

!     Load in the parameters for the simulation

      open (unit=5, file='input/input')
      read (unit=5, fmt='(24(/))')
      read (unit=5, fmt=*) N
      read (unit=5, fmt='(2(/))')
      read (unit=5, fmt=*) NP
      read (unit=5, fmt='(2(/))')
      read (unit=5, fmt=*) TF
      read (unit=5, fmt='(2(/))')
      read (unit=5, fmt=*) indMAX
      read (unit=5, fmt='(2(/))')
      read (unit=5, fmt=*) DT
      read (unit=5, fmt='(2(/))')
      read (unit=5, fmt=*) FRMfile
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
      read (unit=5, fmt=*) COL_type
      read (unit=5, fmt='(2(/))')
      read (unit=5, fmt=*) SAVE_RU
      read (unit=5, fmt='(2(/))')
      read (unit=5, fmt=*) EXIT_WHEN_COLLIDED
      close(5)
      NUM_POSSIBLE_COLLISIONS = N*N - N
      call getpara(para,DT,SIMtype)
      DT0=DT

      NT=N*NP
      ALLOCATE(R(NT,3))
      ALLOCATE(U(NT,3))
      ALLOCATE(R0(NT,3))
      ALLOCATE(U0(NT,3))
      if (COL_type.NE.0) then
         ALLOCATE(HAS_COLLIDED(NT,NT))
         HAS_COLLIDED = -1.0d+0
      else if (COL_type.EQ.2) then
          WRITE(*,*) "Not yet implemented: KD-tree based collision detection."
          stop 1
      else
         ALLOCATE(HAS_COLLIDED(1,1))
         HAS_COLLIDED = -1.0d+0
      endif

!     Setup the initial condition

      call initcond(R,U,NT,N,NP,IDUM,FRMfile,para)

!     Turn on moves for each simulation type

      if (SIMtype.EQ.1) then
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
      elseif (SIMtype.EQ.2) then
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
      elseif (SIMtype.EQ.3) then
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

      call MCsim(R,U,NT,N,NP,NINIT,BROWN,INTON,IDUM,para,MCAMP, &
           SUCCESS,MOVEON,WindoW,SIMtype)

!     Save the conformation and PSI angles

      open (unit = 1, file = 'data/r0', status = 'NEW')
      IB=1
      do 10 I=1,NP
         do 20 J=1,N
            R0(IB,1)=R(IB,1)
            R0(IB,2)=R(IB,2)
            R0(IB,3)=R(IB,3)
            U0(IB,1)=U(IB,1)
            U0(IB,2)=U(IB,2)
            U0(IB,3)=U(IB,3)
            WRITE(1,*) R(IB,1),R(IB,2),R(IB,3)
            IB=IB+1
 20      continue
 10   continue
      close(1)

      open (unit = 1, file = 'data/u0', status = 'NEW')
      IB=1
      do 30 I=1,NP
         do 40 J=1,N
            WRITE(1,*) U(IB,1),U(IB,2),U(IB,3)
            IB=IB+1
 40      continue
 30   continue
      close(1)

!     Begin simulation

      ind=1
      TIME=0.

!     open the output files

      open (unit = 2, file = 'data/out1', status = 'NEW')
      open (unit = 3, file = 'data/out2', status = 'NEW')
      open (unit = 4, file = 'data/out3', status = 'NEW')

      call stress(SIG,R,U,NT,N,NP,para,INTON,SIMtype)

      WRITE(3,*) real(SIG(1,1)),real(SIG(1,2)),real(SIG(1,3)),real(SIG(2,1)),real(SIG(2,2))
      WRITE(4,*) real(SIG(2,3)),real(SIG(3,1)),real(SIG(3,2)),real(SIG(3,3))

      do WHILE (ind.LE.indMAX)

!     Perform a MC simulation, only if NSTEP.NE.0

         call MCsim(R,U,NT,N,NP,NSTEP,BROWN,INTON,IDUM,para,MCAMP, &
              SUCCESS,MOVEON,WindoW,SIMtype)

!     Perform a Brownian dynamics simulation over time step

         if (LOGTIME.EQ.0) then
            TSAVE = TF*ind/indMAX
         else
            TSAVE = DT0*exp((ind-1.)/(indMAX-1.)*log(TF/DT0))
         endif
         if (NSTEP.EQ.0) then
            call BDsim(R,U,NT,N,NP,TIME,TSAVE,DT,BROWN,INTON,IDUM, &
                       para,SIMtype,HAS_COLLIDED,FPT_DIST,COL_type)
         endif

!     Save the conformation and the metrics

         TENS=nint(log10(1.*ind)-0.4999)+1
         write (fileind,'(I4)'), ind

         if (SAVE_RU.NE.0) then

            if (SAVE_RU.EQ.1) then
                snapnm= 'data/r'//fileind((4-TENS+1):4)
            else if (SAVE_RU.EQ.2) then
                snapnm= 'data/r'
                open (unit = 1, file = 'data/t', status = 'REPLACE')
                    WRITE(1,*) fileind((4-TENS+1):4)
                close(1)
            endif
            open (unit = 1, file = snapnm, status = 'REPLACE')
                IB=1
                do 50 I=1,NP
                    do 60 J=1,N
                    WRITE(1,*) R(IB,1),R(IB,2),R(IB,3)
                    IB=IB+1
60                  continue
50              continue
            close(1)

            if (SAVE_RU.EQ.1) then
                snapnm= 'data/u'//fileind((4-TENS+1):4)
            else if (SAVE_RU.EQ.2) then
                snapnm= 'data/u'
            endif
            open (unit = 1, file = snapnm, status = 'REPLACE')
                IB=1
                do 70 I=1,NP
                    do 80 J=1,N
                    WRITE(1,*) U(IB,1),U(IB,2),U(IB,3)
                    IB=IB+1
80                  continue
70              continue
            close(1)

         endif

         snapnm= 'data/coltimes'
         if (COL_type.NE.0) then
            open (unit=1, file=snapnm, status='REPLACE')
                do, I=1,NT
                    WRITE(1,*) ( HAS_COLLIDED(i,j), j=1,NT )
                enddo
            close(1)
         endif

         call stress(SIG,R,U,NT,N,NP,para,INTON,SIMtype)
         call stressp(COR,R,U,R0,U0,NT,N,NP,para,INTON,SIMtype)

         call energy_elas(EELAS,R,U,NT,N,NP,para)
         EPONP=0.
         if (INTON.EQ.1) then
            call energy_ponp(EPONP,R,NT,N,NP,para)
         endif
         WRITE(2,*) real(TIME),real(EELAS(1)),real(EELAS(2)),real(EELAS(3)),real(EPONP),real(COR)


         WRITE(3,*) real(SIG(1,1)),real(SIG(1,2)),real(SIG(1,3)),real(SIG(2,1)),real(SIG(2,2))
         WRITE(4,*) real(SIG(2,3)),real(SIG(3,1)),real(SIG(3,2)),real(SIG(3,3))

         print*, '________________________________________'
         print*, 'Time point ',ind, ' out of', indMAX
         print*, 'Current time ',TIME
         print*, 'Bending energy ', EELAS(1)
         print*, 'Par compression energy ', EELAS(2)
         print*, 'Perp compression energy ', EELAS(3)
         print*, 'Polymer-polymer energy ', EPONP
         print*, 'Current number of beads ', N
         print*, 'Time step ', DT
         print*, 'end-to-end distance poly 1 ', &
              sqrt((R(N,1)-R(1,1))**2.+(R(N,2)-R(1,2))**2.+(R(N,3)-R(1,3))**2.)
         print*, 'Simulation type ', SIMtype

         if (COUNT(HAS_COLLIDED.NE.-1.0d+0) == NUM_POSSIBLE_COLLISIONS) then
             if (EXIT_WHEN_COLLIDED.NE.0) then
                EXIT
             endif
         endif

         ind=ind+1

      enddo
      end


!---------------------------------------------------------------*

