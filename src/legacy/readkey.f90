subroutine READKEY
  ! this subroutine reads in keywords from a parameter file
  ! it sets the various global variables defined in KEYS module
  ! name of the parameter file is param.* where * is a keyword argument
  ! if no keyword argument is supplied, the default is just a file called param
  ! The EXTRAPARAMFILES keyword will allow extra parameter files to be
  ! read in as well

  use KEYS
  use inPUTPARAMS, only : READLinE, READA, READF, READI, REAdo
  use GENUTIL

  implicit none

  ! ---- stuff for inputing the parameter file in free format --------
  CHARACTER*100 :: ARG ! command line argument
  integer :: NUMARG ! number of command line arguments
  integer :: NITEMS ! number of items on the line in the parameter file
  integer :: PF ! input file unit
  LOGICAL :: FILEEND = .FALSE. ! done reading file?
  CHARACTER*100 :: WORD ! keyword
  ! -------------- for reading multiple parameter files --------
  integer, PARAMETER :: MAXNFILES = 10
  CHARACTER*100 :: PARAMFILES(MAXNFILES)
  integer :: NPARAMFILES, NPARAMREAD
  ! ------ for initializing stuff
  integer :: TIMEVAL(8), SEED
  !real(dp) :: ROTMAT(3,3)
  ! ---------------- temporary variables ---------------
  integer :: DUMI, I, TMPI, DUMI1, DUMI2, DUMI3
  CHARACTER*100 :: DUMSTR
  LOGICAL :: LDUM, TRACKDISTSET

  ! ------------------------
  ! set variable defaults
  ! ------------------------
  ACTION = 'none'
  RNGSEED = 0
  VERBOSE = .FALSE.

  ! geometry and energy parameters
  SHEARABLE = .TRUE.
  STRETCHABLE = .TRUE.
  COUPLED = .TRUE.
  LS = 0.1D0;
  LP = 1;
  EC = 0;
  EPERP = 1D3;
  EPAR = 1D3;
  GAM = 1D0;
  STARTNPT = 100;
  MAXNPT = 100;
  forCE = 0D0
  FinITEXT = .FALSE.
  FinITSHEAR = 1D-3
  NEDGESEG = 0
  EDGELS = 0.1D0
  EDGELP = 1;
  EDGEGAM = 1;
  EDGEEPAR = 1D3;
  EDGEEPERP = 1D3;
  EDGEEC = 0;

  ! input/output
  outFILE = '*.out'
  DUMPSNAPSHOTS = .FALSE.
  SNAPSHOTEVERY = 1
  SNAPSHOTFILE = '*.snap.out'
  RESTART = .FALSE.
  RESTARTFILE = '*.snap.out'
  APPENDSNAPSHOTS = .FALSE.
  SKIPREAD = 1
  STARTEQUIL = .FALSE.
  EQUILSAMPLETYPE = 1

  ! monte carlo
  MCPRinTFREQ = 100
  MCTOTSTEPS = 1000
  MCinITSTEPS = 100
  MCSTATSTEPS = 100
  MCoutPUTFREQ = 100

  ADJUSTEVERY = 1000
  FACCTARGET = 0.5D0
  FACCTOL = 0.1D0
  ADJUSTSCL = 2D0
  doREDISC = .FALSE.
  doLOCALMOVES = .FALSE.
  outPUTBEADWEIGHT = .FALSE.
  inTuWEIGHTNPT = 500
  inTRWEIGHTNPT = 50

  ! brownian dynamics
  DELTSCL = 1D-4
  FRICTR = 1D0
  FRICTU = 1D0
  FRICTPERLEN = .FALSE.
  FRICTOB = 10D0
  RAdoB = 1D0
  MOdoB = 1D3
  BDSTEPS = 1000
  BDPRinTEVERY = 1
  BDPRinTLOG = .FALSE.
  LOGRTERM = .FALSE.
  FIXBEAD1 = .FALSE.
  FIXBEADMID = .FALSE.
  RUNGEKUTTA = 4
  STRESSFILE = '*.stress.out'
  GAUSSIANCHAin = .FALSE.
  doBROWN = .TRUE.
  ! coefficient for the relaxation force in the bead-rod brownian dynamics
  ! that keeps the segment length more or less constant
  BRCRELAX = 0.1;
  usePSEUdoforCE = .TRUE.
  CONSTMOD = 1D4
  MU = 0D0
  ! tracking looping first passage times
  TRACKLOOPinG = .FALSE.
  LOOPRAD = 0.1
  LOOPFILE= "*.loop.out"

  inITRANGE = 1D0

  useSTERICS = .FALSE.
  STERRAD = 1D0
  STERSKIP = 1
  STERMOD = 1D3

  MinSEGLEN = 0.1D0
  MAXSEGLEN = 5D0

  ! groups of chains
  PARAMFROMSNAPSHOT = .FALSE.
  NCONNECT = 0
  NCHAin = 1
  SQUARELATTICE = .FALSE.
  NforCE = 0
  forCE = 0D0
  CONNECTIONS = 0
  CONNECTPOS = .TRUE.
  CONNECTUVEC = .FALSE.
  CONPOSMOD = 1D3
  CONUVECMOD = 1D3
  TRACKDISTSET = .FALSE.
  TRACKDIST = 0
  FIXCONNECT = .FALSE.
  NFIXBEAD = 0
  FIXBEAD = 0
  FIXBOUNDARY = .FALSE.
  SETSHEAR = .FALSE.
  SHEARGAMMA = 0D0
  DIAMONDLATTICE = .FALSE.
  NDIAMOND = (/1,1/)
  WIDTHDIAMOND = -1D0
  LENDIAMOND = 1
  STARTCOLLAPSE = .FALSE.
  useBDENERGY = .FALSE.

  RESTART = .FALSE.
  RESTARTFILE = 'start.out'
  SKIPREAD = 0

  EQUILBEADROD = .FALSE.
  STARTEQUILLP = 1D0

  ! -------------------------
  ! Read in all parameter files, starting with the ones specified on command line
  ! --------------------------

  PF = 55 ! i/o unit number to be used for parameter files

  ! get input parameter files from command line
  NPARAMFILES = 0
  NUMARG = COMMAND_ARGUMENT_COUNT()
  if (NUMARG = =0) then
     NPARAMFILES = 1
     PARAMFILES(1) = 'param'
     ARG = ''
  else
     do I = 1,NUMARG
        CALL GETARG(I, ARG)
        NPARAMFILES = NPARAMFILES + 1
        WRITE(DUMSTR,'(A)') 'param.' //TRIM(ADJUSTL(ARG))
        PARAMFILES(NPARAMFILES) = DUMSTR
     ENDdo
     ! reset arg to its original value
     if (NUMARG > 1) CALL GETARG(1,ARG)
  ENDif

  NPARAMREAD = 0 ! keep track of how many files have been read
  do while (NPARAMREAD < NPARAMFILES)
     NPARAMREAD = NPARAMREAD + 1

     PRinT*, 'Reading parameter file: ', PARAMFILES(NPARAMREAD)
     inQUIRE(FILE = PARAMFILES(NPARAMREAD),EXIST = LDUM)
     if (.NOT.LDUM) then
        PRinT*, 'ERROR in READKEY: Parameter file ', TRIM(ADJUSTL(PARAMFILES(NPARAMREAD))), ' does not exist.'
        STOP 1
     ENDif
     open(UNIT = PF, FILE = PARAMFILES(NPARAMREAD), STATUS = 'OLD')

     ! read in the keywords one line at a time
     do
        CALL READLinE(PF,FILEEND,NITEMS)
        if (FILEEND.and.nitems.eq.0) EXIT

        ! skip empty lines
        if (NITEMS == 0) CYCLE

        ! Read in the keyword for this line
        CALL READA(WORD,CASESET = 1)

        ! Skip any empty lines or any comment lines
        if (WORD(1:1) == '#') CYCLE

        SELECT CASE(WORD) ! pick which keyword
        CASE('ACTION')
           CALL READA(ACTION, CASESET = 1)
        CASE('ADJUSTRANGE')
           CALL READI(ADJUSTEVERY)
           if (NITEMS > 2) CALL READF(FACCTARGET)
           if (NITEMS > 3) CALL READF(FACCTOL)
           if (NITEMS > 4) CALL READF(ADJUSTSCL)
        CASE('BDSTEPS')
           CALL READI(BDSTEPS)
           if (NITEMS > 2) CALL READF(BDPRinTEVERY)
           if (NITEMS > 3) CALL REAdo(BDPRinTLOG)
        CASE('BRCRELAX')
           CALL READF(BRCRELAX)
        CASE('CONNECT')
           NCONNECT = NCONNECT + 1
           if (NCONNECT > MAXNCONNECT) then
              PRinT*, 'TOO MANY EXPLICIT CONNECTIONS. RAISE MAXNCONNECT'
              STOP 1
           ENDif
           do DUMI = 1,4
              CALL READI(CONNECTIONS(NCONNECT,DUMI))
           ENDdo
        CASE('CONNECTMOD')
           CALL READF(CONPOSMOD)
           CALL READF(CONUVECMOD)
        CASE('CONNECTTYPE')
           CALL REAdo(CONNECTPOS)
           CALL REAdo(CONNECTUVEC)
        CASE('CONSTMOD')
           CALL READF(CONSTMOD)
        CASE('COUPLED')
           CALL REAdo(COUPLED)
        CASE('DELTSCL')
           CALL READF(DELTSCL)
        CASE('DIAMONDLATTICE')
           DIAMONDLATTICE = .TRUE.
           do DUMI = 1,2
              CALL READI(NDIAMOND(DUMI))
           ENDdo
           CALL READi(LENDIAMOND)
           if (NITEMS > 4) CALL READF(WIDTHDIAMOND)
        CASE('doLOCALMOVES')
           doLOCALMOVES = .TRUE. ! do single bead moves for 1-chain MC
        CASE('EC')
           CALL READF(EC)
        CASE('EDGESEGS')
           CALL READI(NEDGESEG)
           CALL READF(EDGELS)
           CALL READF(EDGELP)
           CALL READF(EDGEGAM)
           CALL READF(EDGEEPAR)
           CALL READF(EDGEEPERP)
           CALL READF(EDGEEC)
        CASE('EPAR')
           CALL READF(EPAR)
        CASE('EPERP')
           CALL READF(EPERP)
        CASE('FinITEXT')
           if (NITEMS > 1) then
              CALL READF(FinITSHEAR)
           ENDif
           FinITEXT = .TRUE.
        CASE('FIXBEAD')
           NFIXBEAD = NFIXBEAD + 1
           if (NFIXBEAD > MAXFIXBEAD) then
              PRinT*, 'ERROR: too many fixed bead lines'
              STOP 1
           ENDif
           CALL READI(FIXBEAD(NFIXBEAD,1))
           if (NITEMS > 2) then
              CALL READI(FIXBEAD(NFIXBEAD,2))
              CALL REAdo(LDUM)
              if (LDUM) FIXBEAD(NFIXBEAD,3) = 1
              CALL REAdo(LDUM)
              if (LDUM) FIXBEAD(NFIXBEAD,4) = 1
           else
              FIXBEAD(NFIXBEAD,2) = 1
           ENDif
        CASE('FIXBEAD1')
           FIXBEAD1 = .TRUE.
        CASE('FIXBEADMID')
           FIXBEADMID = .TRUE.
        CASE('FIXBOUNDARY')
           if (NITEMS > 1) then
              CALL REAdo(FIXBOUNDARY(1))
              CALL REAdo(FIXBOUNDARY(2))
           ENDif
           if (NITEMS > 3) then
              CALL REAdo(FIXBOUNDARY(3))
              CALL REAdo(FIXBOUNDARY(4))
           ENDif
        CASE('FIXCONNECT')
           FIXCONNECT = .TRUE.
        CASE('forCE')
           NforCE = NforCE + 1
           if (NforCE > MAXNforCE) then
              PRinT*, 'TOO MANY forCE! RAISE MAXNforCE.'
              stop 1
           ENDif
           CALL READI(forCEBEAD(NforCE,1))
           CALL READI(forCEBEAD(NforCE,2))
           do DUMI = 1,3
              CALL READF(forCE(NforCE,DUMI))
           ENDdo
        CASE('FRICT')
           CALL READF(FRICTR)
           CALL READF(FRICTU)
           if (NITEMS > 3) then
              CALL REAdo(FRICTPERLEN)
           ENDif
        CASE('GAM')
           CALL READF(GAM)
        CASE('GAUSSIANCHAin')
           GAUSSIANCHAin = .TRUE.
        CASE('inITRANGE')
           do DUMI = 1,4
              CALL READF(inITRANGE(DUMI))
           ENDdo
        CASE('LOGRTERM')
           LOGRTERM = .TRUE.
        CASE('LOOPinG')
           TRACKLOOPinG = .TRUE.
           if (NITEMS > 1) CALL READF(LOOPRAD)
           if (NITEMS > 2) CALL READA(LOOPFILE)
        CASE('LP')
           CALL READF(LP)
        CASE('LS')
           CALL READF(LS)
        CASE('MCPRinTFREQ')
           CALL READI(MCPRinTFREQ)
           if (NITEMS > 2) then
              CALL READI(MCoutPUTFREQ)
           else
              MCoutPUTFREQ = MCPRinTFREQ
           ENDif
        CASE('MCSTEPS')
           CALL READI(MCTOTSTEPS)
           if (NITEMS > 2) then
              CALL READI(MCSTATSTEPS)
           endif
           if (NITEMS > 3) then
              CALL READI(MCinITSTEPS)
           ENDif
        CASE('MU')
           CALL READF(MU)
        CASE('NCHAin')
           CALL READI(NCHAin)
        CASE('NOBROWN')
           doBROWN = .FALSE.
        CASE('NPT')
           ! starting number of points; maximal number
           ! if not specified, assuming maximal number is the starting number
           CALL READI(STARTNPT)
           if (NITEMS > 2) then
              CALL READI(MAXNPT)
           else
              MAXNPT = STARTNPT
           ENDif
        CASE('OBSTACLE')
           CALL READF(RAdoB)
           CALL READF(MOdoB)
           CALL READF(FRICTOB)
        CASE('outFILE')
           CALL READA(outFILE)
        CASE('outPUTBEADWEIGHT')
           ! output the partition function associated with each mobile bead
           ! integrating over R and U vecs separately
           outPUTBEADWEIGHT = .TRUE.
           if (NITEMS > 1) then
              ! number of integration points in each dim when integrating over u vector
              CALL READI(inTUWEIGHTNPT)
           ENDif
           if (NITEMS > 2) then
              CALL READI(inTRWEIGHTNPT)
           ENDif
        CASE('PARAMFROMSNAPSHOT')
           if (NITEMS > 1) then
              CALL REAdo(PARAMFROMSNAPSHOT)
           else
              PARAMFROMSNAPSHOT = .TRUE.
           ENDif
        CASE('REDISCRETIZE')
           doREDISC = .TRUE.
           if (NITEMS > 1) CALL READF(MinSEGLEN)
           if (NITEMS > 2) CALL READF(MAXSEGLEN)
        CASE('RESTART')
           RESTART = .TRUE.
           if (NITEMS > 1) CALL READA(RESTARTFILE)
           if (NITEMS > 2) CALL READI(SKIPREAD)
        CASE('RNGSEED')
           CALL READI(RNGSEED)
        CASE('RUNGEKUTTA')
           CALL READI(RUNGEKUTTA)
        CASE('SETSHEAR')
           SETSHEAR = .TRUE.
           CALL READF(SHEARGAMMA)
        CASE('SHEARABLE')
           CALL REAdo(SHEARABLE)
        CASE('SNAPSHOTS')
           DUMPSNAPSHOTS = .TRUE.
           if (NITEMS > 1) CALL READI(SNAPSHOTEVERY)
           if (NITEMS > 2) CALL READA(SNAPSHOTFILE)
           if (NITEMS > 3) CALL REAdo(APPENDSNAPSHOTS)
        CASE('STARTEQUIL')
           ! start with properly sampled equilibrium conformations
           STARTEQUIL = .TRUE.
           if (NITEMS > 1) CALL READI(EQUILSAMPLETYPE)
           if (NITEMS > 2) then
              EQUILBEADROD = .TRUE.
              CALL READF(STARTEQUILLP)
           ENDif
        CASE('SQUARELATTICE')
           SQUARELATTICE = .TRUE.
        CASE('STARTCOLLAPSE')
           STARTCOLLAPSE = .TRUE.
        CASE('STERICS')
           useSTERICS = .TRUE.
           CALL READF(STERRAD)
           if (NITEMS > 2) CALL READI(STERSKIP)
           if (NITEMS > 3) CALL READF(STERMOD)
        CASE('STRESSFILE')
           CALL READA(STRESSFILE)
        CASE('STRETCHABLE')
           CALL REAdo(STRETCHABLE)
        CASE('TRACKDIST')
           TRACKDISTSET = .TRUE.
           do DUMI = 1,4
              CALL READI(TRACKDIST(DUMI))
           ENDdo
        CASE('useBDENERGY')
           useBDENERGY = .TRUE. ! use BD energy for MC calculations
        CASE('usePSEUdoforCE')
           ! use pseudo-potential force for bead-rod BD simulations?
           CALL REAdo(usePSEUdoforCE)
        CASE('VERBOSE')
           CALL REAdo(VERBOSE)
        CASE DEFAULT
           print*, 'ERROR: unidentified keyword ', TRIM(WORD), " Will ignore."
        END SELECT
     ENDdo
     CLOSE(PF)
  ENDdo

  ! ----- set some more defaults -----
  if (.NOT.TRACKDISTSET) then
     TRACKDIST = (/1,1,STARTNPT,1/)
  ENDif

  ! -----------------
  ! check validity of some values, raise errors or adjust as necessary
  ! -----------------

  if (STARTNPT <= 0.OR.MAXNPT < STARTNPT) then
     PRinT*, 'ERROR in NPT VALUES',STARTNPT,MAXNPT
     STOP 1
  ENDif
  if (EPERP < 0) then
     PRinT*, 'ERROR in EPERP VALUE', EPERP
     STOP 1
  ENDif
  if (EPAR < 0) then
     PRinT*, 'ERROR in EPAR VALUE', EPAR
     STOP 1
  ENDif
  if (LS < 0) then
     PRinT*, 'ERROR in LS VALUE', LS
     STOP 1
  ENDif
  if (LP < 0) then
     PRinT*, 'ERROR in LP VALUE', LP
     STOP 1
  ENDif

  if (DIAMONDLATTICE) then
     ! reset number of chains and length of chains based on diamond lattice
     NCHAin = 2*(NDIAMOND(1) + NDIAMOND(2)-1)
     MAXNPT = 2*MinVAL(NDIAMOND)*LENDIAMOND + 1
     if (WIDTHDIAMOND < 0) then
        WIDTHDIAMOND = GAM*LS*LENDIAMOND/SQRT(2D0)*2
     ENDif
     PRinT*, 'Recalculating nchain and maxnpt for diamond lattice:', NCHAin, MAXNPT, WIDTHDIAMOND
  ENDif

  if (TRACKDIST(1) <= 0.OR.TRACKDIST(1) > MAXNPT&
       & .OR.TRACKDIST(3) <= 0.OR.TRACKDIST(3) > MAXNPT&
       & .OR.TRACKDIST(2) <= 0.OR.TRACKDIST(2) > NCHAin &
       & .OR.TRACKDIST(4) <= 0.OR.TRACKDIST(4) > NCHAin) then
     PRinT*, 'ERROR: BAD TRACKDIST', TRACKDIST
     STOP 1
  ENDif

  do DUMI = 1,NforCE
     if (forCEBEAD(DUMI,1) <= 0.OR.forCEBEAD(DUMI,1) > MAXNPT &
          & .OR.forCEBEAD(DUMI,2) <= 0.OR.forCEBEAD(DUMI,2) > NCHAin) then
        PRinT*, 'ERROR: BAD forCE', forCEBEAD(DUMI,:)
        STOP 1
     ENDif
  ENDdo

  ! ----------- fix file names -----------
  CALL REPLACESUBSTR(outFILE,'*',TRIM(ADJUSTL(ARG)))
  CALL REPLACESUBSTR(SNAPSHOTFILE,'*',TRIM(ADJUSTL(ARG)))
  CALL REPLACESUBSTR(RESTARTFILE,'*',TRIM(ADJUSTL(ARG)))
  CALL REPLACESUBSTR(STRESSFILE,'*',TRIM(ADJUSTL(ARG)))
  CALL REPLACESUBSTR(LOOPFILE,'*',TRIM(ADJUSTL(ARG)))
  ! ---------------------------

  ! Initiate random number generator
  if (RNGSEED == 0) then
     ! use the current time of day in milliseconds
     CALL DATE_AND_TIME(VALUES = TIMEVAL)
     SEED = TIMEVAL(5)*3600*1000 + TIMEVAL(6)*60*1000 + TIMEVAL(7)*1000 + TIMEVAL(8)
  elseif (RNGSEED == -1) then
     ! use the last 5 characters in the command-line argument
     SEED = STRinG2NUM(TRIM(ADJUSTL(ARG)))
  elseif (RNGSEED == -2) then
     ! use the last 4 characters in the command-line argument
     ! and additionally the millisecond time
     CALL DATE_AND_TIME(VALUES = TIMEVAL)
     SEED = STRinG2NUM(TRIM(ADJUSTL(ARG)),TIMEVAL(8))
  else
     ! use this seed directly
     SEED = RNGSEED
  ENDif

  print*, 'Initiating Mersenne twister random number generator with seed:', SEED
  CALL SGRND(SEED)

  print*, '------------Parameter values : -------------------'
  print*, 'ACTION: ', TRIM(ADJUSTL(ACTION))
  print*, 'Output file: ', TRIM(outFILE)
  if (DUMPSNAPSHOTS) then
     PRinT*, 'Dumping snapshot every', SNAPSHOTEVERY,'steps. In file:', TRIM(ADJUSTL(SNAPSHOTFILE))
  ENDif
  if (RESTART) then
     PRinT*, 'Restarting from file:', trim(adjustl(RESTARTFILE))
  ENDif
  print*, 'Number of chains:', NCHAin
  print*, 'STARTNPT, MAXNPT, LS:', STARTNPT, MAXNPT,LS
  PRinT*, 'LP, GAM, EPAR, EPERP, EC:', LP, GAM, EPAR, EPERP, EC
  PRinT*, 'FinITE EXTENSION?:', FinITEXT, FinITSHEAR
  PRinT*, 'FRICTION COEFFICIENTS:', FRICTR, FRICTU
  PRinT*, 'OBSTACLE:', RAdoB,MOdoB,FRICTOB
  PRinT*, 'CONSTRAinT, STERIC MODULUS, mu:', CONSTMOD, STERMOD, MU
  if (useSTERICS) then
     PRinT*, 'Using sterics, with radius:', STERRAD
  ENDif
  PRinT*, 'NUMBER OF CONNECTIONS:', NCONNECT, SQUARELATTICE
  PRinT*, 'FIXED CONNECTIONS?:', FIXCONNECT
  print*, 'Tracking distance btwn points:', TRACKDIST
  if (NFIXBEAD > 0) then
     PRinT*, 'FIXED BEADS:'
     do DUMI = 1,NFIXBEAD
        PRinT*, FIXBEAD(DUMI,:)
     ENDdo
  ENDif
  if (ANY(FIXBOUNDARY)) PRinT*, 'FIXinG BOUNDARIES.', FIXBOUNDARY
  if (GAUSSIANCHAin) PRinT*, 'Treating chain as a plain gaussian with stretch modulus EPAR'
  if (STARTEQUIL) PRinT*, 'Starting from equilibrated configurations.'
  if (NEDGESEG > 0) PRinT*, 'For ', NEDGESEG, ' edge segments parameters are:', &
       & EDGELS, EDGELP, EDGEGAM, EDGEEPAR, EDGEEPERP, EDGEEC
  print*, '----------------------------------------------------'


END subroutine READKEY
