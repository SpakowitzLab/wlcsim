SUBROUTINE READKEY
  ! this subroutine reads in keywords from a parameter file
  ! it sets the various global variables defined in KEYS module
  ! name of the parameter file is param.* where * is a keyword argument
  ! if no keyword argument is supplied, the default is just a file called param
  ! The EXTRAPARAMFILES keyword will allow extra parameter files to be
  ! read in as well

  USE KEYS
  USE INPUTPARAMS, ONLY : READLINE, READA, READF, READI, READO
  USE GENUTIL

  IMPLICIT NONE

  ! ---- stuff for inputing the parameter file in free format --------
  CHARACTER*100 :: ARG ! command line argument
  INTEGER :: NUMARG ! number of command line arguments
  INTEGER :: NITEMS ! number of items on the line in the parameter file
  INTEGER :: PF ! input file unit
  LOGICAL :: FILEEND=.FALSE. ! done reading file?
  CHARACTER*100 :: WORD ! keyword
  ! -------------- for reading multiple parameter files --------
  INTEGER, PARAMETER :: MAXNFILES = 10
  CHARACTER*100 :: PARAMFILES(MAXNFILES)
  INTEGER :: NPARAMFILES, NPARAMREAD
  ! ------ for initializing stuff
  INTEGER :: TIMEVAL(8), SEED
  !DOUBLE PRECISION :: ROTMAT(3,3)
  ! ---------------- temporary variables ---------------
  INTEGER :: DUMI, I, TMPI, DUMI1, DUMI2, DUMI3
  CHARACTER*100 :: DUMSTR
  LOGICAL :: LDUM, TRACKDISTSET

  ! ------------------------
  ! set variable defaults
  ! ------------------------
  ACTION = 'NONE'
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
  FORCE = 0D0
  FINITEXT = .FALSE.
  FINITSHEAR = 1D-3
  NEDGESEG = 0
  EDGELS = 0.1D0
  EDGELP = 1;
  EDGEGAM = 1;
  EDGEEPAR = 1D3;
  EDGEEPERP = 1D3;
  EDGEEC = 0;

  ! input/output
  OUTFILE = '*.out'
  DUMPSNAPSHOTS = .FALSE.
  SNAPSHOTEVERY = 1
  SNAPSHOTFILE = '*.snap.out'
  RESTART = .FALSE.
  RESTARTFILE = '*.snap.out'
  APPENDSNAPSHOTS = .FALSE.
  SKIPREAD=1
  STARTEQUIL = .FALSE.
  EQUILSAMPLETYPE = 1

  ! monte carlo
  MCPRINTFREQ = 100
  MCTOTSTEPS = 1000
  MCINITSTEPS = 100
  MCSTATSTEPS = 100
  MCOUTPUTFREQ = 100

  ADJUSTEVERY = 1000
  FACCTARGET = 0.5D0
  FACCTOL = 0.1D0
  ADJUSTSCL = 2D0
  DOREDISC = .FALSE.
  DOLOCALMOVES = .FALSE.
  OUTPUTBEADWEIGHT = .FALSE.
  INTuWEIGHTNPT = 500
  INTRWEIGHTNPT = 50

  ! brownian dynamics
  DELTSCL = 1D-4
  FRICTR = 1D0
  FRICTU = 1D0
  FRICTPERLEN = .FALSE.
  FRICTOB = 10D0
  RADOB = 1D0
  MODOB = 1D3
  BDSTEPS = 1000
  BDPRINTEVERY = 1
  BDPRINTLOG = .FALSE.
  LOGRTERM = .FALSE.
  FIXBEAD1 = .FALSE.
  FIXBEADMID = .FALSE.
  RUNGEKUTTA = 4
  STRESSFILE = '*.stress.out'
  GAUSSIANCHAIN = .FALSE.
  DOBROWN = .TRUE.
  ! coefficient for the relaxation force in the bead-rod brownian dynamics
  ! that keeps the segment length more or less constant
  BRCRELAX = 0.1;
  USEPSEUDOFORCE = .TRUE.
  CONSTMOD = 1D4
  MU = 0D0
  ! tracking looping first passage times
  TRACKLOOPING=.FALSE.
  LOOPRAD = 0.1
  LOOPFILE= "*.loop.out"

  INITRANGE = 1D0

  USESTERICS = .FALSE.
  STERRAD = 1D0
  STERSKIP=1
  STERMOD = 1D3

  MINSEGLEN = 0.1D0
  MAXSEGLEN = 5D0

  ! groups of chains
  PARAMFROMSNAPSHOT = .FALSE.
  NCONNECT = 0
  NCHAIN = 1
  SQUARELATTICE = .FALSE.
  NFORCE = 0
  FORCE = 0D0
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
  USEBDENERGY=.FALSE.

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
  IF (NUMARG==0) THEN
     NPARAMFILES = 1
     PARAMFILES(1) = 'param'
     ARG = ''
  ELSE
     DO I = 1,NUMARG
        CALL GETARG(I, ARG)
        NPARAMFILES = NPARAMFILES + 1
        WRITE(DUMSTR,'(A)') 'param.' //TRIM(ADJUSTL(ARG))
        PARAMFILES(NPARAMFILES) = DUMSTR
     ENDDO
     ! reset arg to its original value
     IF (NUMARG.GT.1) CALL GETARG(1,ARG)
  ENDIF

  NPARAMREAD = 0 ! keep track of how many files have been read
  DO WHILE (NPARAMREAD.LT.NPARAMFILES)
     NPARAMREAD = NPARAMREAD + 1

     PRINT*, 'Reading parameter file: ', PARAMFILES(NPARAMREAD)
     INQUIRE(FILE=PARAMFILES(NPARAMREAD),EXIST=LDUM)
     IF (.NOT.LDUM) THEN
        PRINT*, 'ERROR in READKEY: Parameter file ', TRIM(ADJUSTL(PARAMFILES(NPARAMREAD))), ' does not exist.'
        STOP 1
     ENDIF
     OPEN(UNIT=PF, FILE=PARAMFILES(NPARAMREAD), STATUS='OLD')

     ! read in the keywords one line at a time
     DO
        CALL READLINE(PF,FILEEND,NITEMS)
        IF (FILEEND.and.nitems.eq.0) EXIT

        ! skip empty lines
        IF (NITEMS.EQ.0) CYCLE

        ! Read in the keyword for this line
        CALL READA(WORD,CASESET=1)

        ! Skip any empty lines or any comment lines
        IF (WORD(1:1).EQ.'#') CYCLE

        SELECT CASE(WORD) ! pick which keyword
        CASE('ACTION')
           CALL READA(ACTION, CASESET=1)
        CASE('ADJUSTRANGE')
           CALL READI(ADJUSTEVERY)
           IF (NITEMS.GT.2) CALL READF(FACCTARGET)
           IF (NITEMS.GT.3) CALL READF(FACCTOL)
           IF (NITEMS.GT.4) CALL READF(ADJUSTSCL)
        CASE('BDSTEPS')
           CALL READI(BDSTEPS)
           IF (NITEMS.GT.2) CALL READF(BDPRINTEVERY)
           IF (NITEMS.GT.3) CALL READO(BDPRINTLOG)
        CASE('BRCRELAX')
           CALL READF(BRCRELAX)
        CASE('CONNECT')
           NCONNECT = NCONNECT + 1
           IF (NCONNECT.GT.MAXNCONNECT) THEN
              PRINT*, 'TOO MANY EXPLICIT CONNECTIONS. RAISE MAXNCONNECT'
              STOP 1
           ENDIF
           DO DUMI = 1,4
              CALL READI(CONNECTIONS(NCONNECT,DUMI))
           ENDDO
        CASE('CONNECTMOD')
           CALL READF(CONPOSMOD)
           CALL READF(CONUVECMOD)
        CASE('CONNECTTYPE')
           CALL READO(CONNECTPOS)
           CALL READO(CONNECTUVEC)
        CASE('CONSTMOD')
           CALL READF(CONSTMOD)
        CASE('COUPLED')
           CALL READO(COUPLED)
        CASE('DELTSCL')
           CALL READF(DELTSCL)
        CASE('DIAMONDLATTICE')
           DIAMONDLATTICE = .TRUE.
           DO DUMI = 1,2
              CALL READI(NDIAMOND(DUMI))
           ENDDO
           CALL READi(LENDIAMOND)
           IF (NITEMS.GT.4) CALL READF(WIDTHDIAMOND)
        CASE('DOLOCALMOVES')
           DOLOCALMOVES = .TRUE. ! do single bead moves for 1-chain MC
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
        CASE('FINITEXT')
           IF (NITEMS.GT.1) THEN
              CALL READF(FINITSHEAR)
           ENDIF
           FINITEXT = .TRUE.
        CASE('FIXBEAD')
           NFIXBEAD = NFIXBEAD + 1
           IF (NFIXBEAD.GT.MAXFIXBEAD) THEN
              PRINT*, 'ERROR: too many fixed bead lines'
              STOP 1
           ENDIF
           CALL READI(FIXBEAD(NFIXBEAD,1))
           IF (NITEMS.GT.2) THEN
              CALL READI(FIXBEAD(NFIXBEAD,2))
              CALL READO(LDUM)
              IF (LDUM) FIXBEAD(NFIXBEAD,3) = 1
              CALL READO(LDUM)
              IF (LDUM) FIXBEAD(NFIXBEAD,4) = 1
           ELSE
              FIXBEAD(NFIXBEAD,2) = 1
           ENDIF
        CASE('FIXBEAD1')
           FIXBEAD1 = .TRUE.
        CASE('FIXBEADMID')
           FIXBEADMID = .TRUE.
        CASE('FIXBOUNDARY')
           IF (NITEMS.GT.1) THEN
              CALL READO(FIXBOUNDARY(1))
              CALL READO(FIXBOUNDARY(2))
           ENDIF
           IF (NITEMS.GT.3) THEN
              CALL READO(FIXBOUNDARY(3))
              CALL READO(FIXBOUNDARY(4))
           ENDIF
        CASE('FIXCONNECT')
           FIXCONNECT=.TRUE.
        CASE('FORCE')
           NFORCE = NFORCE + 1
           IF (NFORCE.GT.MAXNFORCE) THEN
              PRINT*, 'TOO MANY FORCE! RAISE MAXNFORCE.'
              stop 1
           ENDIF
           CALL READI(FORCEBEAD(NFORCE,1))
           CALL READI(FORCEBEAD(NFORCE,2))
           DO DUMI = 1,3
              CALL READF(FORCE(NFORCE,DUMI))
           ENDDO
        CASE('FRICT')
           CALL READF(FRICTR)
           CALL READF(FRICTU)
           IF (NITEMS.GT.3) THEN
              CALL READO(FRICTPERLEN)
           ENDIF
        CASE('GAM')
           CALL READF(GAM)
        CASE('GAUSSIANCHAIN')
           GAUSSIANCHAIN = .TRUE.
        CASE('INITRANGE')
           DO DUMI = 1,4
              CALL READF(INITRANGE(DUMI))
           ENDDO
        CASE('LOGRTERM')
           LOGRTERM = .TRUE.
        CASE('LOOPING')
           TRACKLOOPING = .TRUE.
           IF (NITEMS.GT.1) CALL READF(LOOPRAD)
           IF (NITEMS.GT.2) CALL READA(LOOPFILE)
        CASE('LP')
           CALL READF(LP)
        CASE('LS')
           CALL READF(LS)
        CASE('MCPRINTFREQ')
           CALL READI(MCPRINTFREQ)
           IF (NITEMS.GT.2) THEN
              CALL READI(MCOUTPUTFREQ)
           ELSE
              MCOUTPUTFREQ = MCPRINTFREQ
           ENDIF
        CASE('MCSTEPS')
           CALL READI(MCTOTSTEPS)
           IF (NITEMS.GT.2) THEN
              CALL READI(MCSTATSTEPS)
           endif
           IF (NITEMS.GT.3) THEN
              CALL READI(MCINITSTEPS)
           ENDIF
        CASE('MU')
           CALL READF(MU)
        CASE('NCHAIN')
           CALL READI(NCHAIN)
        CASE('NOBROWN')
           DOBROWN = .FALSE.
        CASE('NPT')
           ! starting number of points; maximal number
           ! if not specified, assuming maximal number is the starting number
           CALL READI(STARTNPT)
           IF (NITEMS.GT.2) THEN
              CALL READI(MAXNPT)
           ELSE
              MAXNPT = STARTNPT
           ENDIF
        CASE('OBSTACLE')
           CALL READF(RADOB)
           CALL READF(MODOB)
           CALL READF(FRICTOB)
        CASE('OUTFILE')
           CALL READA(OUTFILE)
        CASE('OUTPUTBEADWEIGHT')
           ! output the partition function associated with each mobile bead
           ! integrating over R and U vecs separately
           OUTPUTBEADWEIGHT = .TRUE.
           IF (NITEMS.GT.1) THEN
              ! number of integration points in each dim when integrating over u vector
              CALL READI(INTUWEIGHTNPT)
           ENDIF
           IF (NITEMS.GT.2) THEN
              CALL READI(INTRWEIGHTNPT)
           ENDIF
        CASE('PARAMFROMSNAPSHOT')
           IF (NITEMS.GT.1) THEN
              CALL READO(PARAMFROMSNAPSHOT)
           ELSE
              PARAMFROMSNAPSHOT = .TRUE.
           ENDIF
        CASE('REDISCRETIZE')
           DOREDISC = .TRUE.
           IF (NITEMS.GT.1) CALL READF(MINSEGLEN)
           IF (NITEMS.GT.2) CALL READF(MAXSEGLEN)
        CASE('RESTART')
           RESTART = .TRUE.
           IF (NITEMS.GT.1) CALL READA(RESTARTFILE)
           IF (NITEMS.GT.2) CALL READI(SKIPREAD)
        CASE('RNGSEED')
           CALL READI(RNGSEED)
        CASE('RUNGEKUTTA')
           CALL READI(RUNGEKUTTA)
        CASE('SETSHEAR')
           SETSHEAR = .TRUE.
           CALL READF(SHEARGAMMA)
        CASE('SHEARABLE')
           CALL READO(SHEARABLE)
        CASE('SNAPSHOTS')
           DUMPSNAPSHOTS = .TRUE.
           IF (NITEMS.GT.1) CALL READI(SNAPSHOTEVERY)
           IF (NITEMS.GT.2) CALL READA(SNAPSHOTFILE)
           IF (NITEMS.GT.3) CALL READO(APPENDSNAPSHOTS)
        CASE('STARTEQUIL')
           ! start with properly sampled equilibrium conformations
           STARTEQUIL = .TRUE.
           IF (NITEMS.GT.1) CALL READI(EQUILSAMPLETYPE)
           IF (NITEMS.GT.2) THEN
              EQUILBEADROD = .TRUE.
              CALL READF(STARTEQUILLP)
           ENDIF
        CASE('SQUARELATTICE')
           SQUARELATTICE = .TRUE.
        CASE('STARTCOLLAPSE')
           STARTCOLLAPSE = .TRUE.
        CASE('STERICS')
           USESTERICS = .TRUE.
           CALL READF(STERRAD)
           IF (NITEMS.GT.2) CALL READI(STERSKIP)
           IF (NITEMS.GT.3) CALL READF(STERMOD)
        CASE('STRESSFILE')
           CALL READA(STRESSFILE)
        CASE('STRETCHABLE')
           CALL READO(STRETCHABLE)
        CASE('TRACKDIST')
           TRACKDISTSET = .TRUE.
           DO DUMI = 1,4
              CALL READI(TRACKDIST(DUMI))
           ENDDO
        CASE('USEBDENERGY')
           USEBDENERGY = .TRUE. ! use BD energy for MC calculations
        CASE('USEPSEUDOFORCE')
           ! use pseudo-potential force for bead-rod BD simulations?
           CALL READO(USEPSEUDOFORCE)
        CASE('VERBOSE')
           CALL READO(VERBOSE)
        CASE DEFAULT
           print*, 'ERROR: unidentified keyword ', TRIM(WORD), " Will ignore."
        END SELECT
     ENDDO
     CLOSE(PF)
  ENDDO

  ! ----- set some more defaults -----
  IF (.NOT.TRACKDISTSET) THEN
     TRACKDIST = (/1,1,STARTNPT,1/)
  ENDIF

  ! -----------------
  ! check validity of some values, raise errors or adjust as necessary
  ! -----------------

  IF (STARTNPT.LE.0.OR.MAXNPT.LT.STARTNPT) THEN
     PRINT*, 'ERROR IN NPT VALUES',STARTNPT,MAXNPT
     STOP 1
  ENDIF
  IF (EPERP.LT.0) THEN
     PRINT*, 'ERROR IN EPERP VALUE', EPERP
     STOP 1
  ENDIF
  IF (EPAR.LT.0) THEN
     PRINT*, 'ERROR IN EPAR VALUE', EPAR
     STOP 1
  ENDIF
  IF (LS.LT.0) THEN
     PRINT*, 'ERROR IN LS VALUE', LS
     STOP 1
  ENDIF
  IF (LP.LT.0) THEN
     PRINT*, 'ERROR IN LP VALUE', LP
     STOP 1
  ENDIF

  IF (DIAMONDLATTICE) THEN
     ! reset number of chains and length of chains based on diamond lattice
     NCHAIN = 2*(NDIAMOND(1)+NDIAMOND(2)-1)
     MAXNPT = 2*MINVAL(NDIAMOND)*LENDIAMOND + 1
     IF (WIDTHDIAMOND.LT.0) THEN
        WIDTHDIAMOND = GAM*LS*LENDIAMOND/SQRT(2D0)*2
     ENDIF
     PRINT*, 'Recalculating nchain and maxnpt for diamond lattice:', NCHAIN, MAXNPT, WIDTHDIAMOND
  ENDIF

  IF (TRACKDIST(1).LE.0.OR.TRACKDIST(1).GT.MAXNPT&
       & .OR.TRACKDIST(3).LE.0.OR.TRACKDIST(3).GT.MAXNPT&
       & .OR.TRACKDIST(2).LE.0.OR.TRACKDIST(2).GT.NCHAIN &
       & .OR.TRACKDIST(4).LE.0.OR.TRACKDIST(4).GT.NCHAIN) THEN
     PRINT*, 'ERROR: BAD TRACKDIST', TRACKDIST
     STOP 1
  ENDIF

  DO DUMI = 1,NFORCE
     IF (FORCEBEAD(DUMI,1).LE.0.OR.FORCEBEAD(DUMI,1).GT.MAXNPT &
          & .OR.FORCEBEAD(DUMI,2).LE.0.OR.FORCEBEAD(DUMI,2).GT.NCHAIN) THEN
        PRINT*, 'ERROR: BAD FORCE', FORCEBEAD(DUMI,:)
        STOP 1
     ENDIF
  ENDDO

  ! ----------- fix file names -----------
  CALL REPLACESUBSTR(OUTFILE,'*',TRIM(ADJUSTL(ARG)))
  CALL REPLACESUBSTR(SNAPSHOTFILE,'*',TRIM(ADJUSTL(ARG)))
  CALL REPLACESUBSTR(RESTARTFILE,'*',TRIM(ADJUSTL(ARG)))
  CALL REPLACESUBSTR(STRESSFILE,'*',TRIM(ADJUSTL(ARG)))
  CALL REPLACESUBSTR(LOOPFILE,'*',TRIM(ADJUSTL(ARG)))
  ! ---------------------------

  ! Initiate random number generator
  IF (RNGSEED.EQ.0) THEN
     ! use the current time of day in milliseconds
     CALL DATE_AND_TIME(VALUES=TIMEVAL)
     SEED = TIMEVAL(5)*3600*1000 + TIMEVAL(6)*60*1000 + TIMEVAL(7)*1000 + TIMEVAL(8)
  ELSEIF (RNGSEED.EQ.-1) THEN
     ! use the last 5 characters in the command-line argument
     SEED = STRING2NUM(TRIM(ADJUSTL(ARG)))
  ELSEIF (RNGSEED.EQ.-2) THEN
     ! use the last 4 characters in the command-line argument
     ! and additionally the millisecond time
     CALL DATE_AND_TIME(VALUES=TIMEVAL)
     SEED = STRING2NUM(TRIM(ADJUSTL(ARG)),TIMEVAL(8))
  ELSE
     ! use this seed directly
     SEED = RNGSEED
  ENDIF

  print*, 'Initiating Mersenne twister random number generator with seed:', SEED
  CALL SGRND(SEED)

  print*, '------------Parameter values : -------------------'
  print*, 'ACTION: ', TRIM(ADJUSTL(ACTION))
  print*, 'Output file: ', TRIM(OUTFILE)
  IF (DUMPSNAPSHOTS) THEN
     PRINT*, 'Dumping snapshot every', SNAPSHOTEVERY,'steps. In file:', TRIM(ADJUSTL(SNAPSHOTFILE))
  ENDIF
  IF (RESTART) THEN
     PRINT*, 'Restarting from file:', trim(adjustl(RESTARTFILE))
  ENDIF
  print*, 'Number of chains:', NCHAIN
  print*, 'STARTNPT, MAXNPT, LS:', STARTNPT, MAXNPT,LS
  PRINT*, 'LP, GAM, EPAR, EPERP, EC:', LP, GAM, EPAR, EPERP, EC
  PRINT*, 'FINITE EXTENSION?:', FINITEXT, FINITSHEAR
  PRINT*, 'FRICTION COEFFICIENTS:', FRICTR, FRICTU
  PRINT*, 'OBSTACLE:', RADOB,MODOB,FRICTOB
  PRINT*, 'CONSTRAINT, STERIC MODULUS, mu:', CONSTMOD, STERMOD, MU
  IF (USESTERICS) THEN
     PRINT*, 'Using sterics, with radius:', STERRAD
  ENDIF
  PRINT*, 'NUMBER OF CONNECTIONS:', NCONNECT, SQUARELATTICE
  PRINT*, 'FIXED CONNECTIONS?:', FIXCONNECT
  print*, 'Tracking distance btwn points:', TRACKDIST
  IF (NFIXBEAD.GT.0) THEN
     PRINT*, 'FIXED BEADS:'
     DO DUMI = 1,NFIXBEAD
        PRINT*, FIXBEAD(DUMI,:)
     ENDDO
  ENDIF
  IF (ANY(FIXBOUNDARY)) PRINT*, 'FIXING BOUNDARIES.', FIXBOUNDARY
  IF (GAUSSIANCHAIN) PRINT*, 'Treating chain as a plain gaussian with stretch modulus EPAR'
  IF (STARTEQUIL) PRINT*, 'Starting from equilibrated configurations.'
  IF (NEDGESEG.GT.0) PRINT*, 'For ', NEDGESEG, ' edge segments parameters are:', &
       & EDGELS, EDGELP, EDGEGAM, EDGEEPAR, EDGEEPERP, EDGEEC
  print*, '----------------------------------------------------'


END SUBROUTINE READKEY
