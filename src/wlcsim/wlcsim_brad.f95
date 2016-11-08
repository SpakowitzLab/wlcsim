!---------------------------------------------------------------*
      
PROGRAM wlcsim
      
!     This simulation models the equilibrium conformation statistics of DNA molecules modeled
!     as a discrete,stretchable, shearable worm-like chain (see Koslover, 2013) using a Metropolis
!     Monte Carlo procedure. Brownian dynamics simulations are also included.

!     The program involves performing a series of monte carlo simulations with a specified number of monte
!     carlo steps. After each monte carlo simulation, the configuration, energy, and structural parameters
!     of the resulting polymer configuration are saved for future analysis.

!     This version of the simulation allows for a parallel tempering scheme in which multiple "replicas"
!     with different linking numbers are simulated in parallel, periodically testing for exchange between
!     replicas. Note that this parallel tempering scheme applies only to closed-circular DNA.

!     DNA molecules may be linear or circular. Topology of circular DNA molecules are specified
!     by their linking number. This allows the properties of supercoiled DNA molecules to be studied.

!     Brad Krajina, Andrew J. Spakowitz.
!     Last modified: 10/2/2014     



  use mpi
  use mersenne_twister 
  IMPLICIT NONE

  double precision, PARAMETER :: PI=3.141592654d0 ! Value of pi

  !Variables used in simulation 

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:):: R	 ! Conformation of polymer chains
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:):: U	 ! Conformation of polymer chains
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:):: R0	 ! Conformation of polymer chains
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:):: U0	 ! Conformation of polymer chains


  DOUBLE PRECISION L0       ! Equilibrium segment length
  DOUBLE PRECISION DEL      ! Segment contour length
  DOUBLE PRECISION ENERGY   ! Total energy
  DOUBLE PRECISION TIME     ! Current time
  DOUBLE PRECISION TSAVE     ! Time of save point
  DOUBLE PRECISION T0,TF    ! Initial/final times
  DOUBLE PRECISION DT       ! Time step size
  INTEGER I,J,IB            ! Index
  INTEGER INDMAX            ! Maximum index in series
  INTEGER IND               ! Ind in series
  INTEGER TENS              ! Decimal of index
  character*5 fileind       ! Index of output
  character*16 snapnm       ! File for output


  !     Simulation input variables

  INTEGER NT                 ! Number of beads in simulation
  INTEGER N                 ! Number of beads in simulation
  INTEGER NP                ! Number of polymers in simulation
  DOUBLE PRECISION L        ! Total contour length
  DOUBLE PRECISION LP       ! Bending persistencce length
  DOUBLE PRECISION LT       ! Twist persistence length
  INTEGER LK                ! Linking number
  INTEGER FRMFILE           ! Initial condition
  INTEGER BROWN             ! Include Brownian forces
  INTEGER INTON             ! Include polymer interactions
  INTEGER LOGTIME           ! Is data recorded in log time?
  INTEGER RING              ! Is polymer a ring?
  INTEGER TWIST             ! Include twist?
  DOUBLE PRECISION DT0      ! Initial time step size
  INTEGER NSTEP,NINIT

  !     Monte Carlo variables

  DOUBLE PRECISION MCAMP(6) ! Amplitude of random change
  INTEGER MOVEON(6)			! Is the move active
  INTEGER WINDOW(6)			! Size of window for bead selection
  INTEGER SUCCESS(6),SUCCESS_TOTAL(6)        ! Number of successes

  !MPI variables
  integer IOStatus
  integer (kind = 4) error
  integer (kind = 4) p
  integer (kind = 4) id
  integer (kind = 4) status(MPI_STATUS_SIZE)
  integer rep
  integer dest
  integer source
  integer, allocatable :: nodeNumber(:)
  !   variable for random number generator 
  type(random_stat) rand_stat  ! state of random number chain
  integer Irand     ! Seed
  character*8 datedum  ! trash
  character*10 timedum ! trash
  character*5 zonedum  ! trash
  integer seedvalues(8) ! clock readings
  real urand(1)
  type(random_stat) rstat  ! state of random number chain
  type(random_stat) rstati  ! state of random number chain
  type(random_stat), allocatable ::  stat(:)

!x  type(random_stat), allocatable:: stat(:) !for random numer generator

  !     Parallel tempering variables

  INTEGER ParTempOn         !Use parallel tempering?
  INTEGER, allocatable :: Lks(:)               !Vector of Lks for difference replicas
  DOUBLE PRECISION, allocatable :: Es(:)       !Vector of total energies for different replicas
  DOUBLE PRECISION, allocatable :: EElasRep(:,:) !Matrix of elastic energies for replicas
  DOUBLE PRECISION, allocatable :: Wrs(:)        !Vector of writhes for different replicas
  DOUBLE PRECISION, allocatable :: ETwists(:)    !Vector of twist energies for different replicas
  DOUBLE PRECISION ETwist
  DOUBLE PRECISION dE_exchange                 !Change in energy for replica exchange
  INTEGER TempLk                               !Temporary Lk value for performing exchange of Lk
  INTEGER NLks              !Number of LKs for replica parallel tempering
  INTEGER LkMin             !Minimum LK for parallel tempering
  INTEGER LkMax             !Maximum LK for parallel tempering
  INTEGER LkStep            !Step size between LK values
  CHARACTER*4 LkStr         !String variable for Lk
  CHARACTER*4 LkPlusStr     !String variable for Lk replica one Lk step above current replica
  CHARACTER*4 LkMinusStr      !String variable for Lk replica one LK step below current replica
  CHARACTER*4 LkSwapStr     !String variable for Lk replica to test for exchange
  DOUBLE PRECISION Prob     !Probability of swapping current replica with test replica 
  DOUBLE PRECISION Test     !Test number for determining whether to swap (random number between 0. and 1.)
  INTEGER, allocatable :: NSwapPlus(:)         !Number of times replica with Lk one step above is swapped
  INTEGER, allocatable :: NSwapMinus(:)        !Number of times replica with Lk one step below is swapped
  INTEGER, allocatable :: NTrialPlus(:)        !Number of trials for exchanging with the replica on Lk step above 
  INTEGER, allocatable :: NTrialMinus(:)       !Number of trials for exchanging with the replica on Lk step below
  INTEGER Restart              !Is simulation restarting from a previously interrupted run (1 or 0)?
  DOUBLE PRECISION, ALLOCATABLE ::  Swap(:,:)     !Matrix of swaps with other replicas (+1,-1,or 0). Each column corresponds to a replica
  DOUBLE PRECISION, ALLOCATABLE ::  PSwapPlus(:)  !Vector of swap probabilities of replicas with Lk above
  DOUBLE PRECISION, ALLOCATABLE :: PSwapMinus(:)  !Vector of swap probabilities of replicas with Lk below

  !     Energy variables

  DOUBLE PRECISION EELAS(4) ! Elastic energy
  DOUBLE PRECISION EPONP    ! Poly-poly energy
  DOUBLE PRECISION ETOT     ! Total chain energy
  DOUBLE PRECISION, ALLOCATABLE :: EELASALL(:,:,:)
  DOUBLE PRECISION, ALLOCATABLE :: EPONPALL(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: ETotAll(:,:)
  DOUBLE PRECISION  ETotAvg
  DOUBLE PRECISION  ETotStdev
  DOUBLE PRECISION  EtotStderr

  !     Structure analysis

  DOUBLE PRECISION, ALLOCATABLE :: RCOM(:,:) ! Center of mass
  DOUBLE PRECISION, ALLOCATABLE :: RGYR(:)
  DOUBLE PRECISION, ALLOCATABLE :: R2(:)  !Second moment of end-to-end displacement 
  DOUBLE PRECISION, ALLOCATABLE :: R4(:)  !Fourth moment of end-to-end displacement 
  DOUBLE PRECISION, ALLOCATABLE :: R6(:)  !Sixth moment of end-to-end displacement 
  DOUBLE PRECISION, ALLOCATABLE :: DR(:)  !End-to-end displacement 

  DOUBLE PRECISION, ALLOCATABLE :: RGYRALL(:,:) !Trajectory of the radius of gyrationover the simulation for all replicas
  DOUBLE PRECISION, ALLOCATABLE :: WrAll(:,:) !Trajectory of writhe over the simulation for all replicas
  DOUBLE PRECISION, ALLOCATABLE ::  RGYRSQALL(:,:) ! Trajectory of the Radius of gyration squared over the simulation for all replicas
  DOUBLE PRECISION, ALLOCATABLE ::  RCOMSQ(:)  ! Center of mass squared
  DOUBLE PRECISION, ALLOCATABLE ::  RGYRSQ(:)  ! Radius of gyration squared


  DOUBLE PRECISION RGYRSQ_AVG
  DOUBLE PRECISION RGYRSQ_STDEV
  DOUBLE PRECISION RGYRSQ_STDER
  DOUBLE PRECISION RGYR_AVG
  DOUBLE PRECISION RGYR_STDEV
  DOUBLE PRECISION RGYR_STDER
  DOUBLE PRECISION R2_AVG
  DOUBLE PRECISION R2_STDEV
  DOUBLE PRECISION R2_STDER
  DOUBLE PRECISION R4_AVG
  DOUBLE PRECISION R4_STDEV
  DOUBLE PRECISION R4_STDER
  DOUBLE PRECISION R6_AVG
  DOUBLE PRECISION R6_STDEV
  DOUBLE PRECISION R6_STDER
  DOUBLE PRECISION DR_AVG
  DOUBLE PRECISION DR_STDEV
  DOUBLE PRECISION DR_STDER
  DOUBLE PRECISION WR_AVG
  DOUBLE PRECISION WR_STDEV
  DOUBLE PRECISION WR_STDER
  DOUBLE PRECISION Etot_avg
  DOUBLE PRECISION EPonp_avg
  DOUBLE PRECISION Etot_stdev
  DOUBLE PRECISION EPonp_stdev
  DOUBLE PRECISION Etot_stder
  DOUBLE PRECISION EPonp_Stder
  DOUBLE PRECISION EElas_avg(4)
  DOUBLE PRECISION EElas_stdev(4)
  DOUBLE PRECISION EElas_stder(4)

  DOUBLE PRECISION Wr       ! Writhe
  DOUBLE PRECISION Tw       ! Twist
  DOUBLE PRECISION DELR(3)  ! Mag of gyration tensor
  DOUBLE PRECISION RCOM0(3) ! Init val RCOM
  DOUBLE PRECISION DELR0(3) ! Init val DELR
  DOUBLE PRECISION DRCOM    ! Change in RCOM
  DOUBLE PRECISION SIG(3,3)
  DOUBLE PRECISION COR

  !     Algorithm analysis variables
  DOUBLE PRECISION, ALLOCATABLE ::  MCAMP1ALL(:)
  DOUBLE PRECISION, ALLOCATABLE :: MCAMP2ALL(:)
  DOUBLE PRECISION, ALLOCATABLE :: WINDOW1ALL(:)
  DOUBLE PRECISION, ALLOCATABLE :: WINDOW2ALL(:)
  DOUBLE PRECISION, ALLOCATABLE :: MCSTEPCUM(:)
  DOUBLE PRECISION, ALLOCATABLE :: SUCCESSALL(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: PHITALL(:,:)
  DOUBLE PRECISION, ALLOCATABLE ::  RGYRSQ_AUTO(:)  !auto-correlation of radius of gyration
  DOUBLE PRECISION, ALLOCATABLE :: RSQ_AUTO(:)     !auto-correlation of RSQ
  DOUBLE PRECISION, ALLOCATABLE :: Wr_AUTO(:)                      !auto-correlation of Wr
  DOUBLE PRECISION, ALLOCATABLE :: Energy_AUTO(:)                      !auto-correlation of Wr

  INTEGER N_auto !maximum difference in indices between elements on which auto-correlation is computed

  !     Variables in the simulation

  DOUBLE PRECISION PARA(10)

  !     Variables for the random number generators

  INTEGER IDUM              ! Seed for the generator
  DOUBLE PRECISION MOM(6)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Initialize MPI
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call MPI_Init ( error )
  if (error.ne.0) then
     print*, "MPI_Init", error
  endif

  !     Get number of processes
  call MPI_Comm_size ( MPI_COMM_WORLD, p, error )
  if (error.ne.0) then
     print*, "MPI_Comm_size", error
  endif

  !     Get individual process id
  call MPI_Comm_rank ( MPI_COMM_WORLD, id, error )
  if (error.ne.0) then
     print*, "MPI_Comm_rank", error
  endif

  !  Print a message.
  if ( id == 0 ) then
     write ( *, '(a)' ) ' '
     write ( *, '(a)' ) ' basicPT_mpi:'
     write ( *, '(a)' ) '  FORTRAN90/MPI version'
     write ( *, '(a)' ) '  This is a basic parallel tempering MC sim'
     write ( *, '(a)' ) ' '
     write ( *, '(a,i8)' ) '  The number of threads being used is ', p
     write ( *, '(a,i8)' ) '  The number of replicas is ', p-1


  endif



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Load in the parameters for the simulation
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  open (unit=5, file='input/input')
  read (unit=5, fmt='(32(/))')
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
  read (unit=5, fmt=*) RING
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) TWIST
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) LOGTIME
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) NINIT
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) NSTEP
  close(5)

  !Read from the parallel tempering input file to 
  !get the Lks for the simulation

  NLKs = 0

   open (unit = 1, file = 'input/Lks')
   do
      read(unit = 1, fmt = *,iostat=IOstatus) TempLk
      if (IOstatus /= 0) exit
      NLks = NLKs + 1
   end do
   close(unit = 1)

   ALLOCATE(Lks(NLks))

   open(unit = 1, file = 'input/Lks')
   do i = 1, NLks
      read(unit = 1,fmt = *) Lks(i)
   enddo
   close(unit = 1)

  !Exit if the number of Lks is not the same as the number of replicas

  if (NLks.NE.(p-1)) then
     if (id == 0) then
        print *, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        print *, 'Number of Lks must be equal to number of worker nodes!'
        print *, 'Exiting...'
        print *, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     endif
     call exit()
  endif
  !     Get elastic parameters and allocate vectors

  call getpara(PARA,DT,DEL,L,LP,LT,Lk,RING)
  DT0=DT
  NT=N*NP
  ALLOCATE(R(NT,3))
  ALLOCATE(U(NT,3))
  ALLOCATE(R0(NT,3))
  ALLOCATE(U0(NT,3))
  ALLOCATE(RCOM(NP,3))
  ALLOCATE(RCOMSQ(NP))
  ALLOCATE(RGYRSQ(NP))
  ALLOCATE(RGYR(NP))
  ALLOCATE(RGYRALL(INDMAX,NP))
  ALLOCATE(R2(NP))
  ALLOCATE(R4(NP))
  ALLOCATE(R6(NP))
  ALLOCATE(DR(NP))

  ALLOCATE(Swap(NLks,INDMAX))

  N_auto=INDMAX/10

  ALLOCATE(RGYRSQ_AUTO(N_auto))
  ALLOCATE(RSQ_AUTO(N_auto))
  ALLOCATE(Wr_AUTO(N_auto))
  ALLOCATE(Energy_AUTO(N_auto))
  ALLOCATE(stat(NLks))
  ALLOCATE(PSwapPlus(NLks))
  ALLOCATE(PSwapMinus(NLks))

  !Set which MC moves to use

  MOVEON(1)=1
  MOVEON(2)=1
  MOVEON(3)=1
  MOVEON(4)=1

  if (RING.EQ.1) then
     MOVEON(3)=0
  endif

  if (INTON.EQ.1.AND.NP.GT.1) then
     MOVEON(5)=1
     MOVEON(6)=1
  else
     MOVEON(5)=0
     MOVEON(6)=0
  endif

  !Allocate vectors for replicas
  ALLOCATE(Wrs(NLks))
  ALLOCATE(Es(NLks))
  ALLOCATE(EElasRep(4,NLks))

  ALLOCATE(ETwists(NLks))
  ALLOCATE(nodeNumber(NLks))
  ALLOCATE(NSwapPlus(NLks))
  ALLOCATE(NSwapMinus(NLks))
  ALLOCATE(NTrialPlus(NLks))
  ALLOCATE(NTrialMinus(NLks))

  ALLOCATE(WrAll(INDMAX,NLks))
  ALLOCATE(MCAMP1ALL(INDMAX))
  ALLOCATE(MCAMP2ALL(INDMAX))
  ALLOCATE(WINDOW1ALL(INDMAX))
  ALLOCATE(WINDOW2ALL(INDMAX))
  ALLOCATE(MCSTEPCUM(INDMAX))
  ALLOCATE(SUCCESSALL(INDMAX,6))
  ALLOCATE(PHITALL(INDMAX,6))
  ALLOCATE(EELASALL(INDMAX,4,NLks))
  ALLOCATE(EPONPALL(INDMAX,NLks))
  ALLOCATE(ETotAll(INDMAX,NLks))
  ALLOCATE(RGYRSQALL(INDMAX,NLks))

  ! Initially replica numbers are same as nodes
  do rep=1,NLks
     nodeNumber(rep)=rep
  enddo

  ! Initialize values 
  WINDOW = N
  NSwapPlus = 0
  NSwapMinus = 0
  NTrialPlus = 0
  NTrialMinus = 0


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Begin Simulation
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Check process id number to see if you are a worker or 
  !a head node

  if (id == 0) then
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !Generate the random seeds
     !and perform initialization Monte Carlo
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     call date_and_time(datedum,timedum,zonedum,seedvalues)
     Irand=int(-seedvalues(5)*1E7-seedvalues(6)*1E5 &
          -seedvalues(7)*1E3-seedvalues(8))
     Irand=mod(Irand,10000)
     call random_setseed(Irand,rand_stat)
     call random_number(urand,rand_stat)

     do rep=1,NLks
        dest=nodeNumber(rep)

        !Send out the random seed to the worker nodes
        call MPI_Send (Irand,1, MPI_INTEGER, dest,   0, &
             MPI_COMM_WORLD,error )

        !Send out the Lk to the worker nodes
        Lk = Lks(rep)
        call MPI_Send (Lk,1, MPI_INTEGER, dest,   0, &
             MPI_COMM_WORLD,error )
    
     enddo
     !Wait for worker nodes to send back Lk to indicate seed set
     !and initialization complete
     do rep = 1,NLks
        source=nodeNumber(rep)
        call MPI_Recv(Lk,1, MPI_INTEGER, source,   0, &
             MPI_COMM_WORLD,status,error )
     enddo

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !Begin Head Node Monte Carlo Simulation Loop
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     do ind = 1,indmax

        print *, '*******************************************'
        print *, 'MC Step ', ind
        !Send out the Lk to the worker nodes to begin their
        !simulations
        do rep = 1,NLks
           dest = nodeNumber(rep)
           Lk = Lks(rep)
           call MPI_Send (Lk,1, MPI_INTEGER, dest,   0, &
                MPI_COMM_WORLD,error )
        enddo

        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !Hear back from workers
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !Loop over Lk replicas

        do rep = 1,NLks
           Lk = Lks(rep)
           !Determine which node this Lk is running on
           source = nodeNumber(rep)
           !Get the energy and writhe from this Lk
           call MPI_RECV(Wr,1, MPI_DOUBLE_PRECISION, source, 0, &
                MPI_COMM_WORLD,status,error )
           call MPI_RECV(RGYRSQ(1),1, MPI_DOUBLE_PRECISION, source, 0, &
                MPI_COMM_WORLD,status,error )
           call MPI_RECV(Etot,1, MPI_DOUBLE_PRECISION, source, 0, &
                MPI_COMM_WORLD,status,error )
           call MPI_RECV(EPonp,1, MPI_DOUBLE_PRECISION, source, 0, &
                MPI_COMM_WORLD,status,error )
           call MPI_RECV(EELAS,4, MPI_DOUBLE_PRECISION, source, 0, &
                MPI_COMM_WORLD,status,error )

           Wrs(rep) = Wr
           
           !Update the trajectory matrices
           WrAll(ind,rep) = Wr
           RGyrSqAll(ind,rep) = RGYRSQ(1)
           EELASALL(ind,:,rep) = EELAS
           EPONPALL(ind,rep) = EPONP
           ETotAll(ind,rep) = Etot

        enddo

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !Perform Replica Exchange
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !Loop over Lk replicas. Each replica first looks to the replica below, then the one above
        
        do rep = 1,NLks-1
           !We are asking each polymer configuration to swap its Lk. The only contribution to the change in energy
           !is the change in twist energy. The new twist energy is the twist associated with the current writhe
           !evaluated under the trial Lk

           NTrialPlus(rep) = NTrialPlus(rep) + 1

           NTrialMinus(rep+1) = NTrialMinus(rep+1) + 1

           dE_exchange = 0.0d0
           
           dE_exchange = -ETwists(rep) - Etwists(rep+1)

           !Change in twist energy due to rep taking Lk above it
           dE_exchange = dE_exchange +((2*PI*(Lks(rep + 1)-Wrs(rep))**2.d0)*LT/(2.d0*L))
           
           !Change in twist energy due to rep + 1 taking Lk below it
           dE_exchange = dE_exchange +(2*PI*(Lks(rep)-Wrs(rep + 1))**2.d0)*LT/(2.d0*L)

           !Generate a random number for the exchange
           call random_number(urand,rand_stat)

           !Determine whether to exchange
           Prob = exp(-dE_exchange)
           if (urand(1).le.Prob) then
              TempLk = Lks(rep)
              Lks(rep) = Lks(rep + 1)
              Lks(rep + 1) = TempLk
              TempLk = nodeNumber(rep)
              nodeNumber(rep) = nodeNumber(rep+1)
              nodenumber(rep+1) = TempLk
              NSwapPlus(rep) = NSwapPlus(rep) +1
              NSwapMinus(rep+1) = NSwapMinus(rep +1) + 1
           endif
        enddo

     enddo

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !Head node finished with MC loop
     !Calculate and save simulation statistics
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     !Calculate swap prbabilities Note that lowest Lk does not swap below
     !and largest Lk does not swap above
     do rep = 1,NLks-1
        PSwapPlus(rep) = real(NSwapPlus(rep))/real(NTrialPlus(rep))
        PSwapMinus(rep+1) = real(NSwapMinus(rep+1))/real(NTrialMinus(rep+1))
     enddo

     !Calculate averages, standard deviations, and autocorrelations of strucutral quantities
     !and energies

     do rep = 1,NLks

        write(LkStr,'(I4)') Lks(rep)

        !Calculate average, standard deviation, and standard error
        Wr_Avg = sum(WrAll(:,rep))/real(indmax)
        RGyrSq_Avg = sum(RGyrSqAll(:,rep))/real(indmax)
        Etot_Avg = sum(EtotAll(:,rep))/real(indmax)
        EElas_Avg = sum(EElasAll(:,:,rep),dim = 1)/real(indmax)


        !Write to file
        open(1,file = 'data/LK_' //TRIM(ADJUSTL(LkStr)) //'/Wr_Avg', status = 'replace')
        open(2,file = 'data/LK_' //TRIM(ADJUSTL(LkStr)) //'/RGYRSQ_Avg', status = 'replace')
        open(3,file = 'data/LK_' //TRIM(ADJUSTL(LkStr)) //'/ETOT_Avg', status = 'replace')
        open(4,file = 'data/LK_' //TRIM(ADJUSTL(LkStr)) //'/EELAS_Avg', status = 'replace')
        
        write(1,*) Wr_Avg
        write(2,*) RGYRSQ_Avg
        write(3,*) ETOT_Avg
        write(4,*) EELAS_Avg

        close(1)
        close(2)
        close(3)
        close(4)
        
        Wr_stdev = sqrt(sum((WrAll(:,rep) - Wr_Avg)**2.d0)/real(indmax))
        RGyrSq_stdev = sqrt(sum((RGyrSqAll(:,rep) - RGyrSq_Avg)**2.d0)/real(indmax))
        Etot_stdev = sqrt(sum((EtotAll(:,rep) - Etot_Avg)**2.d0)/real(indmax))
        EElas_stdev = sqrt(sum((EElasAll(:,:,rep) - Etot_Avg)**2.d0,dim = 1)/real(indmax))

        !Write to file
        open(1,file = 'data/LK_' //TRIM(ADJUSTL(LkStr)) //'/Wr_stdev', status = 'replace')
        open(2,file = 'data/LK_' //TRIM(ADJUSTL(LkStr)) //'/RGYRSQ_stdev', status = 'replace')
        open(3,file = 'data/LK_' //TRIM(ADJUSTL(LkStr)) //'/ETOT_stdev', status = 'replace')
        open(4,file = 'data/LK_' //TRIM(ADJUSTL(LkStr)) //'/EELAS_stdev', status = 'replace')
        
        write(1,*) Wr_stdev
        write(2,*) RGYRSQ_stdev
        write(3,*) ETOT_stdev
        write(4,*) EELAS_stdev

        close(1)
        close(2)
        close(3)
        close(4)
      
        Wr_stder = Wr_stdev/(sqrt(real(indmax)))
        RGyrSq_stder = RGyrSq_stdev/(sqrt(real(indmax)))
        Etot_stder = Etot_stdev/(sqrt(real(indmax)))
        EElas_stder = EElas_stdev/sqrt(real(indmax))

        !Write to file
        open(1,file = 'data/LK_' //TRIM(ADJUSTL(LkStr)) //'/Wr_stder', status = 'replace')
        open(2,file = 'data/LK_' //TRIM(ADJUSTL(LkStr)) //'/RGYRSQ_stder', status = 'replace')
        open(3,file = 'data/LK_' //TRIM(ADJUSTL(LkStr)) //'/ETOT_stder', status = 'replace')
        open(4,file = 'data/LK_' //TRIM(ADJUSTL(LkStr)) //'/EELAS_stder', status = 'replace')
        
        write(1,*) Wr_Avg
        write(2,*) RGYRSQ_Avg
        write(3,*) ETOT_Avg
        write(4,*) EELAS_Avg

        close(1)
        close(2)
        close(3)
        close(4)

        !Calculate autocorrelations
        CALL auto_correlation_vector(RGyrSqAll(:,rep),indmax,1,N_auto,RGYRSQ_AUTO)
        CALL auto_correlation_vector(WrAll(:,rep),indmax,1,N_auto,Wr_AUTO)

        open(1,file = 'data/LK_' //TRIM(ADJUSTL(LkStr)) //'/Wr_Auto', status = 'replace')
        open(2,file = 'data/LK_' //TRIM(ADJUSTL(LkStr)) //'/RGYRSQ_Auto', status = 'replace')
        
        do i =1,N_auto
           write(1,*) Wr_auto(i)
           write(2,*) RGyrSq_auto(i)
        enddo

        close(1)
        close(2)

        !Write trajectories to file

        open(1,file = 'data/LK_'//TRIM(ADJUSTL(LkStr))//'/Wr_all', status = 'replace')
        open(2,file = 'data/LK_'//TRIM(ADJUSTL(LkStr))//'/RGYRSQ_all', status = 'replace')
        open(3,file = 'data/LK_'//TRIM(ADJUSTL(LkStr))//'/RGYR_all', status = 'replace')
        open(4,file = 'data/LK_'//TRIM(ADJUSTL(LkStr))//'/ETOT_all', status = 'replace')
        open(5,file = 'data/LK_'//TRIM(ADJUSTL(LkStr))//'/EPONP_all', status = 'replace')
        open(6,file = 'data/LK_'//TRIM(ADJUSTL(LkStr))//'/EELAS_all', status = 'replace')
        

        do i = 1,indmax
           write(1,*) WrAll(i,rep)
           write(2,*) RGyrSqAll(i,rep)
           write(3,*) sqrt(RGyrSqAll(i,rep))
           write(4,*) ETotAll(i,rep)
           write(5,*) EPonpAll(i,rep)
           write(6,*) EElasAll(i,:,rep)
        enddo

        close(1)
        close(2)
        close(3)
        close(4)
        close(5)
        close(6)

     enddo


     !Save swap probabilities to file
     open(unit = 1, file = 'data/PSwapPlus', status = 'replace')
     open(unit = 2, file = 'data/PSwapMinus', status = 'replace')

     do rep=1,NLks -1
        write(1,*) PSwapPlus(rep)
        write(2,*) PSwapMinus(rep+1)
     enddo
     close(1)
     close(2)
   else
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Worker nodes
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      source = 0
      !Genderate a random seed for this thread
      call MPI_Recv(Irand,1, MPI_INTEGER, source,   0, &
           MPI_COMM_WORLD,status,error )
      call MPI_Recv(Lk,1, MPI_INTEGER, source,   0, &
           MPI_COMM_WORLD,status,error )

      call random_setseed(Irand*(id+1),rand_stat)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Perform initialization Monte Carlo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Initialize configuration
      call initcond(R,U,NT,N,NP,FRMFILE,PARA,RING,rand_stat)
      !Initialize window size
      WINDOW = N
      !Perform initalization Monte Carlo simulation
      call MCsim(R,U,NT,N,NP,NINIT,BROWN,INTON,IDUM,PARA, &
           MCAMP,SUCCESS,SUCCESS_TOTAL,MOVEON,WINDOW,RING,TWIST,Lk,LT,LP,L,rand_stat)

      !Save configuration to HD
      write(LkStr, '(I4)') Lk

      open(unit = 1, file = 'data/LK_'//TRIM(ADJUSTL(LKStr))//'/r0',status = 'replace')
      do i = 1,N
         write(1,*) R(i,:)
      enddo
      close(unit = 1)

      open(unit = 1, file = 'data/LK_'//TRIM(ADJUSTL(LKStr))//'/u0', status = 'replace')
      do i = 1,N
         write(1,*) U(i,:)
      enddo
      close(unit = 1)
      
      !Communicate with head node after seed set and initialization complete
      call MPI_SEND(Lk,1, MPI_INTEGER, source, 0, &
             MPI_COMM_WORLD,error )

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Worker Monte Carlo Simulation loop
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do ind = 1,indmax

         !Wait for head node to send out Lk value
         source = 0
         call MPI_Recv(Lk,1, MPI_INTEGER, source,   0, &
              MPI_COMM_WORLD,status,error )
         !Perform a Monte Carlo simulation
         call MCsim(R,U,NT,N,NP,NSTEP,BROWN,INTON,IDUM,PARA, &
           MCAMP,SUCCESS,SUCCESS_TOTAL,MOVEON,WINDOW,RING,TWIST,Lk,LT,LP,L,rand_stat)
         
         !Calculate dimensions of polymer, energy, and writhe
         CALL getdim(N,NP,NT,R,RCOM,RCOMSQ,RGYRSQ,R2,R4,R6,DR)
        
         call WRITHE(R,N, Wr)

         !Update trajectory vectors of structural quantities
    !     WrAll(ind) = Wr
    !     RGyrSqAll(ind,:) = RGYRSQ
    !     R2All(ind,:) = R2
                 
         CALL energy_elas(EELAS,R,U,NT,N,NP,PARA,RING,TWIST,Lk,lt,LP,L)
                        
         EPONP=0.
         if (INTON.EQ.1) then
            call  ENERGY_SELF_CHAIN(EPONP,R,NT,N,NP,PARA,RING)
         endif

         !Update energy trajectory vectors

         ETOT=EPONP +SUM(EELAS)

         !Save the polymer configuration
         write(fileInd,'(I5)') ind
         open(unit = 1, file = 'data/LK_'//TRIM(ADJUSTL(LKStr))//'/r'&
              //TRIM(ADJUSTL(fileInd)),status = 'replace')
         do i = 1,N
            write(1,*) R(i,:)
         enddo
         close(unit = 1)

         open(unit = 1, file = 'data/LK_'//TRIM(ADJUSTL(LKStr))//'/u'&
              //TRIM(ADJUSTL(fileInd)),status = 'replace')
         do i = 1,N
            write(1,*) U(i,:)
         enddo
         close(unit = 1)
         

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !Communicate with the head node
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         !Send back Writhe to the head node
         call MPI_SEND(Wr,1, MPI_DOUBLE_PRECISION, source, 0, &
             MPI_COMM_WORLD,error )
         call MPI_SEND(RGYRSQ(1),1, MPI_DOUBLE_PRECISION, source, 0, &
             MPI_COMM_WORLD,error )
         call MPI_SEND(Etot,1, MPI_DOUBLE_PRECISION, source, 0, &
             MPI_COMM_WORLD,error )
         call MPI_SEND(EPonp,1, MPI_DOUBLE_PRECISION, source, 0, &
             MPI_COMM_WORLD,error )
         call MPI_SEND(EELAS,4, MPI_DOUBLE_PRECISION, source, 0, &
             MPI_COMM_WORLD,error )

      enddo
      
   endif

  call MPI_finalize(error)

END PROGRAM wlcsim
      
    

     

      
!---------------------------------------------------------------*
      
