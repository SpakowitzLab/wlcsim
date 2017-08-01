!   --------------------------------------------------------------
!
!   A module containing all global constants, and which defines two types, one
!   for simulation constants (parameters) and another for all the variables that
!   can change during the simulation, to prevent each function from requiring
!   dozens of parameters.
!
!   The various functions for loading/saving/modifying the simulation parameters
!   as a whole should also go here.
!
!   --------------------------------------------------------------
module params
    use, intrinsic :: iso_fortran_env
    use, intrinsic :: IEEE_ARITHMETIC
    use mersenne_twister
    use precision, only: dp
    use inputparams, only: MAXPARAMLEN

    implicit none

    public

    !!!     hardcoded params. will need to change if certain parts of code change
    ! number of wlc_p move types
    integer, parameter :: nMoveTypes = 10

    !!!     arbitrary technical choices
    ! used for all character buffers holding filenames
    integer, parameter :: MAXFILENAMELEN = 500
    ! unique file "units" to use for each file
    integer, parameter :: inFileUnit = 51
    integer, parameter :: outFileUnit = 52
    ! number of digits in max allowed number of save points, also change below
    integer, parameter :: MAX_LOG10_SAVES = 10
    ! format string to use for num2str(save_ind)
    character(5), parameter :: num2strFormatString = '(I10)'

    !!!     universal constants
    ! fully accurate, adaptive precision
    real(dp) :: pi = 4 * atan(1.0_dp)
    real(dp) :: nan
    ! ! won't get optimized away by compiler, see e.g.
    ! ! https://software.intel.com/en-us/forums/intel-visual-fortran-compiler-for-windows/topic/294680
    ! real(dp) :: nan = transfer((/ Z'00000000', Z'7FF80000' /),1.0_dp)
    ! the following would be preferred, but generated compilation errors...
    !real(dp) :: nan = IEEE_VALUE(nan, IEEE_QUIET_NAN)
    real(dp) :: inf
    ! ! doesn't work, inf needs to happen at runtime in fortran
    ! real(dp) :: one = 1.0_dp
    ! real(dp) :: inf = -log(one - one)
    ! the following would be preferred, but generated compilation errors...
    !real(dp) :: inf = ieee_value(inf, ieee_positive_inf)

    ! for all parameters that cannot change during individual simulations
    ! these are documented more thoroughly where they are read in (see the
    ! subroutine get_input_from_file), in the docs (TOdo), and often the default
    ! values will help with understanding how the variable is used.
    !
    ! many of these variables are used only in certain kinds of simulations
    type wlcsim_params
        character(MAXPARAMLEN) codeName ! which simulation code to run
    !   Simulation parameters
        integer simType           ! whether to use WLC, ssWLC, or Gaussian Chain
        integer nT                ! Total number of beads  NT = NP*N*G
        integer nB                ! Number of beads in a polymer NB = N*G
        integer nMpP              ! Number of monomers (NOT BEADS!) in a polymer
        integer nBpM              ! number beads per monomer
        integer nP                ! Number of polymers
        real(dp) dt ! sets time scale of simulation
        real(dp) l  ! length of each polymer in simulation
        real(dp) lp       ! persistence length
        real(dp) lt       ! twist persistence length
        real(dp) l0       ! Equilibrium segment length (same as gam)
        real(dp) beadVolume        ! Bead volume
        real(dp) fPoly    ! volume fraction of Polymer in simulation volume
        real(dp) fA       ! Fraction of A beads
        real(dp) lam      ! Chemical correlation parameter (eigenvalue of transition matrix that generates A/B's)
        real(dp) eb     ! Energy of bending
        real(dp) epar   ! energy "parallel" i.e. psring energy
        real(dp) eperp  ! energy "perpendicular" i.e. shear energy
        real(dp) eself  ! energy "polymer on polymer" (self-interaction)
        real(dp) gam    ! average equilibrium interbead spacing
        real(dp) eta    ! bend-shear coupling parameter
        real(dp) xir    ! drag per unit persistence length
        real(dp) sigma  ! variance of interbead position distribution of appropriately renormalized gaussian chain
        real(dp) xiu    ! rotational drag
        real(dp) eps      ! number of kuhn lengths between beads
        real(dp) del      ! number of persistence lengths between beads
        real(dp) chi      ! Chi parameter value (solvent-polymer) (Flory-Huggins separation constant (how much A/B's hate each))
        real(dp) kap      ! Incompressibility parameter of the melt
        real(dp) collisionRadius ! radius triggering collisions to be recorded in "coltimes"
        real(dp) lhc    !TOdo something to do with intrapolymer interaction strength
        real(dp) vhc    !TOdo something to do with intrapolymer interaction strength, fill in defaults, etc

    !   for passing 1st order phase transition in (quinn/shifan's) random copolymer wlc_p sims
        real(dp) k_field   ! wave vector of applied sinusoidal field (used in PT to step around 1st order phase transition)
        real(dp) hA       ! strength of applied sinusoidal field (used in PT to step around 1st order phase transition)
        real(dp) rend   ! initial end-to-end distance (if applicable in initialization)

    !   for simulating chromatin methylation
        real(dp) EU       ! Energy of HP1 binding for unmethalated sites
        real(dp) EM       ! Energy of HP1 binding for methalated sites
        real(dp) HP1_Bind ! Energy of binding of HP1 to eachother
        real(dp) F_METH   ! Fraction methalated is using option 2 ( for initialization )
        real(dp) LAM_METH ! eigenvalue of transition matrix generating initial methalation
        real(dp) mu       ! chemical potential of HP1

    !   box size things
        integer NBin     ! Number of bins
        integer NBinX(3) ! Number of MC bins on an edge
        integer nColBin  ! Number of collision-detection bins on each edge
        real(dp) lbox(3)  ! monte carlo field bin total box length (approximate)
        real(dp) dbin      ! monte carlo field bin discretization size (approximate)

    !   Monte Carlo Variables (for adaptation)
        integer movetypes
        real(dp) PDesire(nMoveTypes) ! desired hit rate
        real(dp) MAXWindoW(nMoveTypes)         ! Max Size of window for bead selection
        real(dp) MinWindoW(nMoveTypes)         ! Min Size of window for bead selection
        real(dp) MinAMP(nMoveTypes) ! minium amplitude
        real(dp) MaxAMP(nMoveTypes) ! maximum amplitude
        integer MOVEON(nMoveTypes)         ! Is the move active
        real(dp) winTarget(nMoveTypes) ! target for ratio of window to anmplitude
        integer NADAPT(nMoveTypes) ! Nunber of steps between adapt
        real(dp) min_accept  ! threshold for deciding to usually not use a move
        integer reduce_move  ! whether or not to stop usuing a move when it goes below min_accept success
        integer winType      ! distributionof segment size in crankshaft move (unif = 0, exp = 1)

    !   Timing variables
        integer stepsPerExchange   ! number of steps between parallel tempering
        integer nReplicaExchangePerSavePoint ! read teh variable name
        integer numSavePoints      ! total number of save points
        integer stepsPerSave       ! steps per save point
        integer NNoInt             ! save points before turning on NNoInt
        integer N_KAP_ON           ! when to turn KAP energy on
        integer N_CHI_ON           ! when to turn CHI energy on
        integer nInitMCSteps       ! number of mc steps before starting BD

    !   Switches
        logical ring              ! whether the polymer is a ring
        logical twist             ! whether to include twist (wlc_p only for now)
        integer LK                ! Linking number
        integer confinetype       ! type of Boundary Conditions
        integer initCondType           ! initial condition type
        logical field_interactions ! field-based self interactions on
        logical intrapolymer_stick_crossing_enforced ! field-based self interactions on
        logical FRwlc_pHEM           ! read initial chemical sequence from file
        logical FRMchem           ! read initial chemical/methylation state from file
        logical FRMfile           ! read initial condition R from file
        logical FRMField          ! read initial field from file
        integer collisionDetectionType      ! save first passage time vectors to file
        logical exitWhenCollided  ! stop sim with coltimes is full
        logical saveR             ! save R vectors to file
        logical saveU             ! save U vectors to file
        logical saveAB            ! save AB (chemical identity) to file
        logical savePhi           ! save Phi vectors to file
        integer solType           ! Melt vs. Solution, Choose hamiltonian
        logical recenter_on       ! recenter in "quasi"-periodic boundary, should be off in BD
        logical useSchedule       ! use scheduled change in interactions strength(s)
        real(dp) KAP_ON     ! fraction of KAP energy contributing to "calculated" energy
        real(dp) CHI_ON     ! fraction of CHI energy contributing to "calculated" energy
        real(dp) Couple_ON  ! fraction of Coupling energy contributing to "calculated" energy
        logical restart     ! whether we are restarting from a previous sim or not
        logical inTERP_BEAD_LENNARD_JONES ! whether to have inter bead lennard jones energies
        logical field_int_on ! include field interactions (e.g. A/B interactions)
                             ! uses many of the same functions as the chemical
                             ! identity/"meth"ylation code, but energies are calcualted via a field-based approach
        logical bind_On ! chemical identities of each bead are tracked in the "meth" variable

    !   parallel Tempering parameters
        logical PTON    ! whether or not to parallel temper
        logical PT_twist
        logical PT_chi
        logical PT_h
        logical PT_kap
        logical PT_mu
        logical PT_couple

    !   Replica Dynamic Cof choice
        integer NRepAdapt ! number of exchange attemts between adapt
        real(dp) lowerRepExe ! when to decrese cof spacing
        real(dp) upperRepExe ! when to increase cof spacing
        real(dp) lowerCofRail ! minumum acceptable Cof
        real(dp) upperCofRail ! maximum acceptable Cof
        integer indStartRepAdapt
        integer indendRepAdapt
        real(dp) repAnnealSpeed  ! for annealing
        logical replicaBounds
        real(dp) inITIAL_MAX_S

    end type

    ! for variables that can change during the simulation
    type wlcsim_data
        ! a
        ! really
        ! long
        ! comment
        ! for wlcsim

        ! for R
        ! i coudl
        ! do this
        real(dp), allocatable, dimension(:,:):: R   ! Conformation of polymer chains to asldkfjalsdkfj askldf aklsjf aklsf alskf alskdfj alskdfj asldkf asldkf alsdkfj
        real(dp), allocatable, dimension(:,:):: U   ! Conformation of polymer chains
        real(dp), allocatable, dimension(:,:):: RP !Test Bead positions - only valid from IT1 to IT2
        real(dp), allocatable, dimension(:,:):: UP !Test target vectors - only valid from IT1 to IT2
        real(dp), allocatable, dimension(:):: PHIA ! Volume fraction of A
        real(dp), allocatable, dimension(:):: PHIB ! Volume fraction of B
        real(dp), allocatable, dimension(:):: PHIH ! Quinn's sinusoidal field for passing 1st order phase transitions
        real(dp), allocatable, dimension(:):: Vol  ! Volume fraction of A
        integer, allocatable, dimension(:):: AB    ! Chemical identity of beads
        integer, allocatable, dimension(:):: ABP   ! Test Chemical identity of beads
        integer, allocatable, dimension(:):: METH  ! Methalation state of beads
        real(dp), allocatable, dimension(:):: DPHIA    ! Change in phi A
        real(dp), allocatable, dimension(:):: DPHIB    ! Change in phi A
        integer, allocatable, dimension(:) :: indPHI   ! indices of the phi
        ! simulation times at which (i,j)th bead pair first collided
        real(dp), allocatable, dimension(:,:) :: coltimes
        real(dp) :: wr

    !   Twist variables
        real(dp), ALLOCATABLE :: CROSS(:,:)   !Matrix of information for crossings in a 2-D projection of the polymer
        real(dp), ALLOCATABLE :: CROSSP(:,:)  !Matrix of crossings for the trial configuration
        integer NCROSS
        integer NCROSSP
        integer CrossSize


    !   Monte Carlo Variables (for adaptation)
        real(dp) MCAMP(nMoveTypes) ! Amplitude of random change
        real(dp) WindoW(nMoveTypes)         ! Size of window for bead selection
        integer SUCCESS(nMoveTypes)        ! Number of successes
        integer successTOTAL(nMoveTypes)               !Total number of successes
        real(dp) PHit(nMoveTypes) ! hit rate

    !   Energys
        !real(dp) Eint     ! running Eint
        real(dp) eElas(4) ! Elastic force
        real(dp) eChi     ! CHI energy
        real(dp) eKap     ! KAP energy
        real(dp) eCouple  ! Coupling
        real(dp) eBind    ! binding energy
        real(dp) eField   ! Field energy
        real(dp) eSelf    ! repulsive lennard jones on closest approach self-interaction energy (polymer on polymer)
        real(dp) eKnot    ! 0-inf potential for if chain crossed itself

    !   Congigate Energy variables (needed to avoid NaN when cof-> 0 in rep exchange)
        real(dp) x_Chi,   dx_Chi
        real(dp) x_Couple,dx_Couple
        real(dp) x_Kap,   dx_Kap
        real(dp) x_Field, dx_Field
        real(dp) x_Mu,    dx_Mu

    !   Move Variables
        real(dp) DEELAS(4) ! Change in bending energy
    !    real(dp) DEinT    ! Change in self energy
        real(dp) DECouple ! Coupling energy
        real(dp) DEChi    ! chi interaction energy
        real(dp) DEKap    ! compression energy
        real(dp) Debind   ! Change in binding energy
        real(dp) DEField  ! Change in field energy
        real(dp) DESelf   ! change in self interaction energy
        real(dp) ECon     ! Confinement Energy
        integer NPHI  ! NUMBER o phi values that change, i.e. number of bins that were affected

    !   Parallel tempering variables
        integer numProcesses !number of MPI processes running
        integer rep  ! which replica am I
        integer id   ! which thread am I
        integer error  ! MPI error
        integer, allocatable, dimension(:) ::  LKs    !Vector of linking numbers for replicas
        integer nLKs      !Number of linking number replicas to parallel temper over
        real(dp), allocatable, dimension(:) :: Wrs !Vector of writhe for each replica
        real(dp), allocatable, dimension(:,:) :: eelasREPLICAS !elastic energies of replicas
        integer replicaSTART !index for replica to start with for exchange loop
        integer replicaEND   !index for replica to end at for exchange loop
        integer, allocatable, dimension(:) :: nTRIALup !number of times this replica has attempted to swap with replica above
        integer, allocatable, dimension(:) :: nTRIALdown !number of times this replica has attempted to swap with replica below
        integer, allocatable, dimension(:) :: nSWAPup !number of times this replica has swapped with replica above
        integer, allocatable, dimension(:) :: nSWAPdown !number of times this replica has swapped with replica below
        integer, allocatable, dimension(:) :: nodeNUMBER !vector of replicas indices for nodes
        character(MAXFILENAMELEN) repSuffix    ! prefix for writing files


    !   random number generator state
        type(random_stat) rand_stat
        integer rand_seed

    !   indices
        integer mc_ind                  ! current save point index for mc
        integer time_ind                ! current time point
        real(dp) time
    end type


contains

    subroutine set_param_defaults(wlc_p,wlc_d)
        implicit none
        ! WARNinG: changing this to intent(out) means that unassigned values
        ! here will become undefined upon return, due to Fortran's weird
        ! intent(out) semantics for records, this would require that a default
        ! value always be given to new parameters in wlc_p, else we would get a
        ! compile time catchable runtime error that is not caught by gcc as of v5.0
        !
        ! this is almost definitely undesireable, since "undefined" means the
        ! behavior will depend on which compiler is used
        type(wlcsim_params), intent(inout) :: wlc_p
        type(wlcsim_data), intent(inout) :: wlc_d
        integer mctype
        ! file IO
        wlc_p%FRMfile = .FALSE.      ! don't load initial bead positions from file
        wlc_p%FRMCHEM = .FALSE.      ! don't load initial "chem" status from file
        wlc_p%FRMFIELD = .FALSE.     ! don't load initial field values from file
        wlc_p%saveR = .TRUE.         ! do save orientation vectors (makes restart of ssWLC possible)
        wlc_p%saveU = .TRUE.         ! do save orientation vectors (makes restart of ssWLC possible)
        wlc_p%saveAB = .False.     ! dont' save AB by default, almost nobody uses this
        wlc_p%collisionDetectionType = 0         ! don't track first passage time collisions between beads
        wlc_p%collisionRadius = 0    ! never collide except on floating-point coincidence
        wlc_p%savePhi = .FALSE.      ! don't save A/B density per bin (not needed for restart)
        wlc_p%FRwlc_pHEM = .FALSE.      ! don't load initial a/b states from file
        wlc_p%restart = .FALSE.      ! don't restart from previously saved simulation

        ! geometry options
        wlc_p%NP  =1               ! one polymer
        wlc_p%nB  =200             ! 200 beads per polymer
        wlc_p%nBpM = 10
        wlc_p%lp = 1                ! units of lp by default
        wlc_p%lt = 1                ! twist persistence length equals persistence length by default
        wlc_p%nMpP = wlc_p%nB / wlc_p%nBpM
        wlc_p%nT = wlc_p%nP * wlc_p%nB
        wlc_p%lbox(1) = nan     ! box size, *MUST* be set by user
        wlc_p%lbox(2) = nan
        wlc_p%lbox(3) = nan
        wlc_p%nColBin = 1   ! equivalent to collisionDetectionType = 1
        wlc_p%dbin = NaN ! set in tweak_param_defaults
        ! wlc_p%l0  =1.25_dp         ! TOdo: not input
        wlc_p%beadVolume  = 0.1_dp ! much smaller than space between beads
        wlc_p%fA  =0.5_dp  ! half A, half B by default
        wlc_p%LAM =0.0_dp  ! perfectly random sequence  (see generating_sequences.rst for details)
        wlc_p%F_METH = 0.5_dp ! half beads methylated by default
        wlc_p%LAM_METH = 0.9_dp ! highly alternating sequence by default
        wlc_p%fPoly = 0.025_dp   ! volume fraction of plymer corresponding to HELA DNA in cytoplasm
        wlc_p%k_field = 0.0_dp ! some previous values: !1.5708_dp !0.3145_dp

        ! energy parameters
        wlc_p%EPS =0.3_dp ! TOdo: not input
        wlc_p%CHI =0.0_dp ! don't use chi by default
        wlc_p%hA =0.0_dp  ! don't use weird artificial field by default
        wlc_p%KAP =10.0_dp ! "fairly" incompressible --Quinn
        wlc_p%EU  =0.0_dp ! a function of coarse graining. This should be set by hand if needed.
        wlc_p%EM  =0.0_dp ! by default, no hp1 binding energy included
        wlc_p%mu  =0.0_dp ! by default, no hp1 binding included
        wlc_p%HP1_Bind = 0.0_dp ! by default, no binding of HP1 to each other

        ! options
        wlc_p%codeName= "brad" ! not bruno, brad, or quinn, so will error unless specified elsewehre
        wlc_p%movetypes = nMoveTypes
        wlc_p%initCondType = 0 ! 0 for initializing polymer in non-random straight line
        wlc_p%confineType = 0 ! 0 for no confinement
        wlc_p%solType = 1       ! solution, not melt, by default
        wlc_p%ring = .false.    ! not a ring by default
        wlc_p%twist = .false.    ! don't include twist by default
        wlc_p%lk = 0    ! no linking number (lays flat) by default
        wlc_p%min_accept = 0.05 ! if a move succeeds < 5% of the time, start using it only every reduce_move cycles
        wlc_p%exitWhenCollided = .FALSE. ! stop sim when coltimes is full
        wlc_p%field_int_on = .FALSE. ! no field interactions by default
        wlc_p%bind_On = .FALSE. ! no binding energy by default
        wlc_p%inTERP_BEAD_LENNARD_JONES = .FALSE. ! no intrapolymer interactions by default

        ! timing options
        wlc_p%dt  = 1              ! set time scale to unit
        wlc_p%nInitMCSteps = 4000  ! number of initilizing mc steps. 1000s x num polymers is good
        wlc_p%stepsPerSave = 2000  ! number of simulation steps to take
        wlc_p%numSavePoints = 200    ! 200 total save points, i.e. 2000 steps per save point
        wlc_p%NNoInt = 100    ! number of simulation steps before turning on interactions in Quinn's wlc_p scheduler
        wlc_p%reduce_move = 10 ! use moves that fall below the min_accept threshold only once every 10 times they would otherwise be used
        wlc_p%winType = 1   ! exponential fragment sizes mix better
        wlc_p%useSchedule = .False. ! use Quinn's scheduler to modify wlc_p params halfway through the simulation
        wlc_p%KAP_ON = 1.0_dp ! use full value of compression energy
        wlc_p%CHI_ON = 1.0_dp ! use full value of chi energy
        wlc_p%Couple_ON = 1.0_dp ! use full value for coupling energy
        wlc_p%N_KAP_ON = 1 ! turn on compression energy immediately
        wlc_p%N_CHI_ON = 1 ! turn on chi energy immediately
        wlc_p%recenter_on = .TRUE. ! recenter the polymer in the box if it exists the boundary
        wlc_p%inITIAL_MAX_S = 0.0_dp !TOdo: for now must be set explicitly, was 0.1, Quinn, what is this value?

        ! replica options
        wlc_p%PTON = .FALSE.  ! use parallel if applicable
        wlc_p%stepsPerExchange = 100      ! 100 steps between parallel tempering is pretty frequent
        wlc_p%nReplicaExchangePerSavePoint = 1000      ! make this large
        wlc_p%NRepAdapt = 1000  ! 1000 exchange attempts between adaptations
        wlc_p%lowerRepExe = 0.09 ! TOdo: enter justification for these defaults, if any.
        wlc_p%upperRepExe = 0.18 ! TOdo: fine if the only justification is "these just work"
        wlc_p%lowerCofRail = 0.005
        wlc_p%upperCofRail = 0.1
        wlc_p%indStartRepAdapt = 10
        wlc_p%indendRepAdapt = 20
        wlc_p%repAnnealSpeed = 0.01
        wlc_p%replicaBounds = .TRUE.
        wlc_p%PT_twist =.False. ! don't parallel temper linking number (twist) by default
        wlc_p%PT_chi =.False. ! don't parallel temper chi by default
        wlc_p%PT_h =.False. ! don't parallel temper h by default
        wlc_p%PT_kap =.False. ! don't parallel temper kap by default
        wlc_p%PT_mu =.False. ! don't parallel temper mu by default
        wlc_p%PT_couple =.False. ! don't parallel temper HP1 binding by default


        !switches to turn on various types of moves
        wlc_p%MOVEON(1) = 1  ! crank-shaft move
        wlc_p%MOVEON(2) = 1  ! slide move
        wlc_p%MOVEON(3) = 1  ! pivot move
        wlc_p%MOVEON(4) = 1  ! rotate move
        wlc_p%MOVEON(5) = 1  ! full chain rotation
        wlc_p%MOVEON(6) = 1  ! full chain slide
        wlc_p%MOVEON(7) = 1  ! Change in Binding state
        wlc_p%MOVEON(8) = 0  ! Chain flip ! TOdo not working
        wlc_p%MOVEON(9) = 1  ! Chain exchange
        wlc_p%MOVEON(10) = 1 ! Reptation

        ! Balance move amplitude and number of beads
        do mctype = 1,wlc_p%movetypes
            wlc_p%winTarget(mctype) = 8.0_dp
            wlc_p%MinWindoW(mctype) = nan
            wlc_d%MCAMP(mctype) = nan
        enddo
        do mctype = 1,wlc_p%movetypes
            wlc_p%NADAPT(mctype) = 1000 ! adapt after at most 1000 steps
            wlc_p%PDESIRE(mctype) = 0.5_dp ! Target
            wlc_d%SUCCESS(mctype) = 0
            wlc_d%SUCCESStotal(mctype) = 0
            wlc_d%PHIT(mctype) = 0.0_dp
        enddo

    end subroutine

    subroutine read_input_file(infile, wlc_p)
        use inPUTparaMS, only : readLinE, readA, readF, readI, reado
        implicit none
        type(wlcsim_params), intent(inout) :: wlc_p
        character(MAXFILENAMELEN), intent(in) :: infile
        character(MAXPARAMLEN) :: WORD ! parameter name currently being read in
        logical fileend ! have we reached end of file?
        integer nitems  ! number of items read from line
        integer pf      ! file unit for input file

        pf = inFileUnit
        open(unit = PF,file = infile,status = 'OLD')

        ! read in the keywords one line at a time
        do
        call READLinE(PF, fileend, NITEMS)
        if (fileend .and. nitems == 0) exit

        ! skip empty lines
        if (NITEMS == 0) cycle

        ! read in the keyword for this line, convert to upper case for matching
        call readA(WORD, CASESET = 1)

        ! Skip any empty lines or any comment lines
        if (WORD(1:1) == '#') cycle

        select case(WORD) ! pick which keyword, case matchign string must be all uppercase
        case('CODENAME') ! select version of wlcsim to run
            call readA(wlc_p%codeName)
        case('INITCONDTYPE')
            call readI(wlc_p%initCondType)
            ! initCondType      |  Discription
            ! _____________|_________________________________
            !    1         |   straight line in y direction with random starting
            !    2         |   rerandomize when reaching boundary, slit in z dir
            !    3         |   rerandomize when reaching boundary, cube boundary
            !    4         |   rerandomize when reaching boundary, shpere
            !    7         |   initialize as gaussian chain, redraw if outside bdry
        case('CONFINETYPE')
            call readI(wlc_p%confinetype)
            ! confinetype  |  Discription
            ! _____________|_________________________________
            !    0         |  No confinement, periodic cube
            !    1         |  Between two plates in Z direction at 0 and lbox
            !    2         |  Cube of size lbox**3,  range: 0-lbox
            !    3         |  Circle of radius lbox, centered at lbox/2
            !    4         |  Periodic, unequal dimensions
        case('RECENTERON')
            call reado(wlc_p%recenter_on) ! recenter in periodic boundary
        case('SOLTYPE')
            call readI(wlc_p%solType)
            ! solType      | Discription
            !______________|_________________________________
            !    0         | Melt density fluctuates around fixed mean
            !    1         | Solution (For DNA)
        case('FRMCHEM')
            call reado(wlc_p%FRMCHEM) ! Initial chemical/methylation sequence from file
        case('FRMFILE')
            call reado(wlc_p%FRMfile) ! read configuration from file
        case('TWIST')
            call reado(wlc_p%twist) ! whether to include twist energies in wlc_p
        case('RING')
            call reado(wlc_p%ring) ! whether polymer is a ring or not
        case('INTERPBEADLENNARDJONES')
            call reado(wlc_p%inTERP_BEAD_LENNARD_JONES) ! whether polymer is a ring or not
        case('FIELDINTON')
            call reado(wlc_p%field_int_on) !include field interactions
        case('BINDON')
            call reado(wlc_p%bind_on) ! Whether to include a binding state model
        case('LK')
            call readi(wlc_p%lk) ! linking number
        case('PTON')
            call reado(wlc_p%PTON) ! parallel Tempering on
        case('SAVER')
            call reado(wlc_p%saveR)  ! save R vectors to file (every savepoint)
        case('SAVEU')
            call reado(wlc_p%saveU)  ! save R vectors to file (every savepoint)
        case('SAVEAB')
            call reado(wlc_p%saveAB)  ! save AB vectors to file (every savepoint)
        case('SAVEPHI')
            call reado(wlc_p%savePhi) ! save Phi vectors to file (every savepoint)
        case('COLLISIONRADIUS')
            call readf(wlc_p%collisionRadius)
        case('COLLISIONDETECTIONTYPE')
            call readi(wlc_p%collisionDetectionType)
            ! collisionDetectionType   |  Description
            ! _____________|_________________________________
            !    0         |  No tracking fpt
            !    1         |  Use naive, O(n^2) algo for collision checking
            !    2         |  KDtree-based col checking (not implemented)
            !    3         |  custom, fast col checker written by bruno
            !    4         |  bin-based collision pruning, recentering via CoM
        case('EXITWHENCOLLIDED')
            call reado(wlc_p%exitWhenCollided)  ! save u vectors to file (every savepoint)
        case('NB')
            call readi(wlc_p%nb)  ! number of beads in the polymer
        case('DT')
            call readF(wlc_p%dt)  ! time step of simulation. scaled non-dimensionalized time
        case('L')
            call readF(wlc_p%l)  ! actual length in AU of polymer we want to simulate
        case('LT')
            call readF(wlc_p%lt)  ! persistence length
        case('LP')
            call readF(wlc_p%lp)  ! twist persistence length
        case('DBIN')
            call readF(wlc_p%dbin) ! spaitial descretation length, not tested
        case('LBOX')
            call readF(wlc_p%lbox(1)) ! side length of box
            wlc_p%lbox(2) = wlc_p%lbox(1)
            wlc_p%lbox(3) = wlc_p%lbox(1)
        case('LBOXX')
            call readF(wlc_p%lbox(1)) ! side length of box in x direction
        case('LBOXY')
            call readF(wlc_p%lbox(2)) ! side length of box in y direction
        case('LBOXZ')
            call readF(wlc_p%lbox(3)) ! side length of box in z direction
        case('NCOLBIN')
            call readI(wlc_p%nColBin) ! number of bins in each dimension
        case('NP')
            call readI(wlc_p%NP)  ! Number of polymers
        case('NBPM')
            call readI(wlc_p%nBpM) ! Beads per monomer
        case('NMPP')
            call readI(wlc_p%nMpP) ! Number of monomers in a polymer
        case('NNOINT')
            call readI(wlc_p%NNoInt) ! save points before turning on interaction
        case('NKAPON')
            call readI(wlc_p%N_KAP_ON) ! when to turn compression energy on
        case('NCHION')
            call readI(wlc_p%N_CHI_ON) ! when to turn CHI energy on
        case('NUMSAVEPOINTS')
            call readI(wlc_p%numSavePoints) ! total number of save points
        case('NINITMCSTEPS')
            call readI(wlc_p%nInitMCSteps) ! num initial mc steps
        case('STEPSPERSAVE')
            call readI(wlc_p%stepsPerSave) ! steps per save point
        case('STEPSPEREXCHANGE')
            call readI(wlc_p%stepsPerExchange) ! number of steps between parallel tempering
        case('NREPLICAEXCHANGEPERSAVEPOINT')
            call readI(wlc_p%nReplicaExchangePerSavePoint) ! read the variable
        case('FPOLY')
            call readF(wlc_p%Fpoly) ! Fraction Polymer
        case('BEADVOLUME')
            call readF(wlc_p%beadVolume) ! Bead volume
        case('FA')
            call readF(wlc_p%FA) ! Fraction of A beads (fraction bound)
        case('LAM')
            call readF(wlc_p%LAM) ! Chemical correlation parameter
        case('EPS')
            call readF(wlc_p%EPS) ! Elasticity l0/(2lp)
        case('VHC')
            call readF(wlc_p%VHC) ! hard-core lennard jones potential strength
        case('LHC')
            call readF(wlc_p%LHC) ! hard-core lennard jones diameter
        case('CHI')
            call readF(wlc_p%CHI) ! CHI parameter (definition depends on  hamiltoniaon
        case('HA')
            call readF(wlc_p%hA) ! strength of externally applied field
        case('KAP')
            call readF(wlc_p%KAP) !  Incompressibility parameter
        case('EU')
            call readF(wlc_p%EU) ! Energy of binding for unmethalated
        case('EM')
            call readF(wlc_p%EM) ! Energy of binding for methalated
        case('MU')
            call readF(wlc_p%MU) ! chemical potential of HP1
        case('ENERGYCOUPLE')
            call readF(wlc_p%HP1_Bind) ! Energy of binding of HP1 to eachother
        case('FMETH')
            call readF(wlc_p%F_METH) ! Fraction methalated
        case('LAMMETH')
            call readF(wlc_p%LAM_METH) ! eigenvalue of methalation setup
        case('CRANKSHAFTON')
            call readI(wlc_p%MOVEON(1)) ! is Crank shaft move on 1/0
        case('SLIDEON')
            call readI(wlc_p%MOVEON(2)) ! is Slide move on 1/0
        case('PIVOTON')
            call readI(wlc_p%MOVEON(3)) ! is Pivot move on 1/0
        case('ROTATEON')
            call readI(wlc_p%MOVEON(4)) ! is single bead rotate on 1/0
        case('FULLCHAINROTATIONON')
            call readI(wlc_p%MOVEON(5)) ! is full chain rotate on 1/0
        case('FULLCHAINSLIDEON')
            call readI(wlc_p%MOVEON(6)) ! is full chain slide on 1/0
        case('BINDMOVEON')
            call readI(wlc_p%MOVEON(7)) ! is bind/unbind move on 1/0
        case('CHAINFLIPMOVEON')
            call readI(wlc_p%MOVEON(8)) ! is flip move move on 1/0
        case('CHAINSWAPMOVEON')
            call readI(wlc_p%MOVEON(9)) ! is chain swap move on 1/0
        case('REPTATIONMOVEON')
            call readI(wlc_p%MOVEON(10)) ! is reptation move on 1/0
        case('MINCRANKSHAFTWIN')
            call readF(wlc_p%MinWindoW(1)) ! min mean window size
        case('MINSLIDEWIN')
            call readF(wlc_p%MinWindoW(2))
        case('MINPIVOTWIN')
            call readF(wlc_p%MinWindoW(3))
        case('MINBINDWIN')
            call readF(wlc_p%MinWindoW(7))
        case('REDUCEMOVE')
            call readI(wlc_p%reduce_move) !  only exicute unlikely movetypes every ____ cycles
        case('WINTYPE')
            call readI(wlc_p%winType)   ! fragment size distribution for crankshaft move
        case('MINACCEPT')
            call readF(wlc_p%Min_ACCEPT) ! below which moves are turned off
        case('CRANKSHAFTTARGET')
            call readF(wlc_p%winTarget(1)) ! target window size for crank shaft move
        case('SLIDETARGET')
            call readF(wlc_p%winTarget(2)) ! target window size for slide move
        case('PIVOTTARGET')
            call readF(wlc_p%winTarget(3)) ! target window size for Pivot move
        case('STRENGTHSCHEDULE')
            call reado(wlc_p%useSchedule) ! use scheduled ramp in interaction strength(s)
        case('NREPADAPT')
            call readI(wlc_p%NRepAdapt)  ! number of exchange attemts between adapt
        case('LOWERREPEXE')
            call readF(wlc_p%lowerRepExe) ! when to decrease cof spacing
        case('UPPERREPEXE')
            call readF(wlc_p%upperRepExe) ! when to increase cof spacing
        case('LOWERCOFRAIL')
            call readF(wlc_p%lowerCofRail) ! minumum acceptable Cof
        case('UPPERCOFRAIL')
            call readF(wlc_p%upperCofRail) ! maximum acceptable Cof
        case('INDSTARTREPADAPT')
            call readI(wlc_p%indStartRepAdapt) ! ind to start rep. cof. adaptiation on
        case('INDENDREPADAPT')
            call readI(wlc_p%indendRepAdapt) ! turn off rep adapt
        case('REPANNEALSPEED')
            call readF(wlc_p%repAnnealSpeed)  ! max change in cof. every adjust
        case('FRMFIELD')
            call reado(wlc_p%FRMFIELD)  ! read field from file
        case('KFIELD')
            call readF(wlc_p%k_field)  ! wave mode for default field
        case('REPLICABOUNDS')
            call reado(wlc_p%replicaBounds) ! insure that 0 < s < 1
        case('INITIALMAXS')
            call readF(wlc_p%inITIAL_MAX_S) ! inital s of rep with highest s
        case('PARALLELTEMPCHI')
            call reado(wlc_p%PT_chi) ! parallel temper chi
        case('PARALLELTEMPH')
            call reado(wlc_p%PT_h) ! parallel temper h
        case('PARALLELTEMPKAP')
            call reado(wlc_p%PT_kap) ! parallel temper kap
        case('PARALLELTEMPMU')
            call reado(wlc_p%PT_mu) ! parallel temper mu
        case('PARALLELTEMPCOUPLE')
            call reado(wlc_p%PT_couple) ! parallel temper HP1_bind
        case('PARALLELTEMPTWIST')
            call reado(wlc_p%pt_twist)  ! parallel temper over linking numbers
        case('RESTART')
            call reado(wlc_p%restart) ! Restart from parallel tempering
        case default
            print *, "Warning, the input key ", trim(Word), " was not recognized."
            print *, "    ...Checking deprecated variable names"
            select case(WORD) ! pick which keyword
            case('RECENTER_ON')
                call reado(wlc_p%recenter_on) ! recenter in periodic boundary
            case('INTERP_BEAD_LENNARD_JONES')
                call reado(wlc_p%inTERP_BEAD_LENNARD_JONES) ! whether polymer is a ring or not
            case('FIELD_INT_ON')
                call reado(wlc_p%field_int_on) !include field interactions
            case('BIND_ON')
                call reado(wlc_p%bind_on) ! Whether to include a binding state model
            case('LK')
                call readi(wlc_p%lk) ! linking number
            case('PTON')
                call reado(wlc_p%PTON) ! parallel Tempering on
            case('SAVE_R')
                call reado(wlc_p%saveR)  ! save u vectors to file (every savepoint)
            case('FPT_COL_TYPE')
                call readi(wlc_p%collisionDetectionType)  ! save u vectors to file (every savepoint)
                ! collisionDetectionType   |  Description
                ! _____________|_________________________________
                !    0         |  No tracking fpt
                !    1         |  Use naive, O(n^2) algo for collision checking
                !    2         |  KDtree-based col checking (not implemented)
                !    3         |  custom, fast col checker written by bruno
            case('SAVE_U')
                call reado(wlc_p%saveU)  ! save u vectors to file (every savepoint)
            case('SAVE_PHI')
                call reado(wlc_p%savePhi) ! save Phi vectors to file (every savepoint)
            case('EXIT_WHEN_COLLIDED')
                call reado(wlc_p%exitWhenCollided)  ! save u vectors to file (every savepoint)
            case('N_KAP_ON')
                call readI(wlc_p%N_KAP_ON) ! when to turn compression energy on
            case('N_CHI_ON')
                call readI(wlc_p%N_CHI_ON) ! when to turn CHI energy on
            case('H_A')
                call readF(wlc_p%hA) ! strength of externally applied field
            case('HP1_BinD')
                call readF(wlc_p%HP1_Bind) ! Energy of binding of HP1 to eachother
            case('F_METH')
                call readF(wlc_p%F_METH) ! Fraction methalated
            case('LAM_METH')
                call readF(wlc_p%LAM_METH) ! eigenvalue of methalation setup
            case('CRANK_SHAFT_ON')
                call readI(wlc_p%MOVEON(1)) ! is Crank shaft move on 1/0
            case('SLIDE_ON')
                call readI(wlc_p%MOVEON(2)) ! is Slide move on 1/0
            case('PIVOT_ON')
                call readI(wlc_p%MOVEON(3)) ! is Pivot move on 1/0
            case('ROTATE_ON')
                call readI(wlc_p%MOVEON(4)) ! is single bead rotate on 1/0
            case('FULL_CHAIN_ROTATION_ON')
                call readI(wlc_p%MOVEON(5)) ! is full chain rotate on 1/0
            case('FULL_CHAIN_SLIDE_ON')
                call readI(wlc_p%MOVEON(6)) ! is full chain slide on 1/0
            case('BIND_MOVE_ON')
                call readI(wlc_p%MOVEON(7)) ! is bind/unbind move on 1/0
            case('CHAIN_FLIP_MOVE_ON')
                call readI(wlc_p%MOVEON(8)) ! is flip move move on 1/0
            case('CHAIN_SWAP_MOVE_ON')
                call readI(wlc_p%MOVEON(9)) ! is chain swap move on 1/0
            case('REPTATION_MOVE_ON')
                call readI(wlc_p%MOVEON(10)) ! is reptation move on 1/0
            case('MIN_CRANK_SHAFT_WIN')
                call readF(wlc_p%MinWindoW(1)) ! min mean window size
            case('MIN_SLIDE_WIN')
                call readF(wlc_p%MinWindoW(2))
            case('MIN_PIVOT_WIN')
                call readF(wlc_p%MinWindoW(3))
            case('MIN_BIND_WIN')
                call readF(wlc_p%MinWindoW(7))
            case('REDUCE_MOVE')
                call readI(wlc_p%reduce_move) !  only exicute unlikely movetypes every ____ cycles
            case('WIN_TYPE')
                call readI(wlc_p%winType)   ! fragment size distribution for crankshaft move
            case('MIN_ACCEPT')
                call readF(wlc_p%Min_ACCEPT) ! below which moves are turned off
            case('CRANK_SHAFT_TARGET')
                call readF(wlc_p%winTarget(1)) ! target window size for crank shaft move
            case('SLIDE_TARGET')
                call readF(wlc_p%winTarget(2)) ! target window size for slide move
            case('PIVOT_TARGET')
                call readF(wlc_p%winTarget(3)) ! target window size for Pivot move
            case('STRENGTH_SCHEDULE')
                call reado(wlc_p%useSchedule) ! use scheduled ramp in interaction strength(s)
            case('N_REP_ADAPT')
                call readI(wlc_p%NRepAdapt)  ! number of exchange attemts between adapt
            case('LOWER_REP_EXE')
                call readF(wlc_p%lowerRepExe) ! when to decrease cof spacing
            case('UPPER_REP_EXE')
                call readF(wlc_p%upperRepExe) ! when to increase cof spacing
            case('LOWER_COF_RAIL')
                call readF(wlc_p%lowerCofRail) ! minumum acceptable Cof
            case('UPPER_COF_RAIL')
                call readF(wlc_p%upperCofRail) ! maximum acceptable Cof
            case('IND_START_REP_ADAPT')
                call readI(wlc_p%indStartRepAdapt) ! ind to start rep. cof. adaptiation on
            case('IND_END_REP_ADAPT')
                call readI(wlc_p%indendRepAdapt) ! turn off rep adapt
            case('REP_ANNEAL_SPEED')
                call readF(wlc_p%repAnnealSpeed)  ! max change in cof. every adjust
            case('FRMFIELD')
                call reado(wlc_p%FRMFIELD)  ! read field from file
            case('K_FIELD')
                call readF(wlc_p%k_field)  ! wave mode for default field
            case('REPLICA_BOUNDS')
                call reado(wlc_p%replicaBounds) ! insure that 0 < s < 1
            case('INITIAL_MAX_S')
                call readF(wlc_p%inITIAL_MAX_S) ! inital s of rep with highest s
            case('PT_CHI')
                call reado(wlc_p%PT_chi) ! parallel temper chi
            case('PT_H')
                call reado(wlc_p%PT_h) ! parallel temper h
            case('PT_KAP')
                call reado(wlc_p%PT_kap) ! parallel temper kap
            case('PT_MU')
                call reado(wlc_p%PT_mu) ! parallel temper mu
            case('PT_COUPLE')
                call reado(wlc_p%PT_couple) ! parallel temper HP1_bind
            case('PT_TWIST')
                call reado(wlc_p%pt_twist)  ! parallel temper over linking numbers
            case default
                print*, "params%read_input_file: ERROR: Unidentified keyword:", &
                        TRIM(WORD)
                stop 1
            end select
        end select
        end do
        close(PF)
    end subroutine read_input_file


    subroutine idiot_checks(wlc_p, wlc_d)
#if MPI_VERSION
        use mpi
#endif
        implicit none
        type(wlcsim_params), intent(inout) :: wlc_p
        type(wlcsim_data), intent(inout) :: wlc_d
        logical err
        integer (kind = 4) mpi_err

        if (wlc_p%ring) then
            if (wlc_p%NP .gt. 1) then
                print*, "As of the writing of this error message"
                print*, "MC_eelals and possible energy_elas are"
                print*, "not capable of more than one rings"
                stop
            endif
            if (wlc_p%initCondType == 7) then
                print*, "initCondType = 7 doesn't know how to make a ring."
                stop
            endif
        endif
        if ((wlc_p%NT > 200000).OR.(wlc_p%NT.lt.1)) then
            print*, "ERROR: Requested ", wlc_p%NT," beads."
            stop 1
        endif

        if (wlc_p%lBox(1) .ne. wlc_p%lBox(1)) then
            print*, "No box size set.  If you need a box please specify it."
            call stop_if_err(wlc_p%initCondType /= 0, &
                'Only one initial polymer config supported if you''re not '//&
                'using LBOX to define a MC simulation box.')
        else
            if ((wlc_p%NBin > 20000).or.(wlc_p%NBin.lt.1)) then
                print*, "ERROR: Requested ", wlc_p%NBin," bins."
                print*, "You probably don't want this."
                print*, "Comment me out if you do."
                stop 1
            endif
            ! we no longer specify fPoly, it is set in tweak_param_defaults
            ! if (wlc_p%confineType.eq.3) then
            !     if (abs((wlc_p%fPoly*(1.0/6.0_dp)*PI*wlc_p%LBOX(1)**3 / &
            !             (wlc_p%beadVolume*wlc_p%NT)) - 1)>0.02) then
            !          print*, "Error: volume fraction incorrect"
            !          stop
            !      endif
            ! else
            !     if (abs(wlc_p%fPoly*wlc_p%LBOX(1)**3/(wlc_p%beadVolume*wlc_p%NT) - 1)>0.02) then
            !          print*, "Error: volume fraction incorrect"
            !          stop
            !     endif
            ! endif
        endif

        call stop_if_err(wlc_p%collisionDetectionType == 2, &
            'KD-tree based collision detection not yet implemented.')

        call stop_if_err(wlc_p%REND > wlc_p%L, &
            "Requesting initial end-to-end distance larger than polymer length.")

        if (wlc_p%codeName == 'quinn') then
           if ((wlc_p%NBinX(1)-wlc_p%NBinX(2).ne.0).or. &
                (wlc_p%NBinX(1)-wlc_p%NBinX(3).ne.0)) then
              err = wlc_p%solType.eq.1
              call stop_if_err(err, "Solution not tested with non-cube box, more coding needed")
              err = wlc_p%confinetype.ne.4
              call stop_if_err(err, "Unequal boundaries require confinetype = 4")
              err = wlc_p%initCondType.eq.4
              call stop_if_err(err, "You shouldn't put a sphere in and unequal box!")
           endif

           err = wlc_p%NBinX(1)*wlc_p%NBinX(2)*wlc_p%NBinX(3).ne.wlc_p%NBin
           call stop_if_err(err, "error in mcsim. Wrong number of bins")

           !TOdo: replace with semantic descriptions of error encountered, instead
           ! of simply outputting the input that the user put in
           if (wlc_p%NT.ne.wlc_p%nMpP*wlc_p%NP*wlc_p%nBpM) then
              print*, "error in mcsim.  NT = ",wlc_p%NT," nMpP = ",wlc_p%nMpP," NP = ",wlc_p%NP," nBpM = ",wlc_p%nBpM
              stop 1
           endif

           if (wlc_p%NB.ne.wlc_p%nMpP*wlc_p%nBpM) then
              print*, "error in mcsim.  NB = ",wlc_p%NB," nMpP = ",wlc_p%nMpP," nBpM = ",wlc_p%nBpM
              stop 1
           endif

           err = wlc_p%NNoInt.gt.wlc_p%indStartRepAdapt
           call stop_if_err(err, "error in mcsim. don't run adapt without int on")

           err = wlc_p%NNoInt.gt.wlc_p%N_CHI_ON
           call stop_if_err(err, "error in mcsim. Can't have chi without int on")

           err = wlc_p%NNoInt.gt.wlc_p%N_KAP_ON
           call stop_if_err(err, "error in mcsim. Can't have kap without int on")

        endif


#if MPI_VERSION
    if (wlc_p%pt_twist) then
        if (.NOT.wlc_p%twist) then
            print *, 'parallel tempering on twist, but twist off'
            stop
        endif
        if (wlc_d%nLKs + 1.ne.wlc_d%numProcesses) then
            print *, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            print *, 'number of threads not equal to number of replicas!'
            print *, 'exiting...'
            print *, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            stop
        endif
  endif
#endif
    end subroutine


    subroutine get_input_from_file(infile, wlc_d, wlc_p)
        ! Based on Elena's readkeys subroutine
        implicit none
        type(wlcsim_params), intent(out) :: wlc_p
        type(wlcsim_data), intent(out) :: wlc_d
        character(MAXFILENAMELEN), intent(in) :: infile  ! file with parameters

        ! baseline defaults
        call set_param_defaults(wlc_p, wlc_d)

        call read_input_file(infile, wlc_p)

        ! advanced defaults that require some inputs to specify
        call tweak_param_defaults(wlc_p, wlc_d)

        ! get derived parameters that aren't directly input from file
        call get_derived_parameters(wlc_p)

        !If parallel tempering is on, read the Lks
        if (wlc_p%pt_twist) then
            call get_LKs_from_file(wlc_d)
        endif

        call printDescription(wlc_p)
        call idiot_checks(wlc_p, wlc_d)

    end subroutine



    subroutine initialize_wlcsim_data(wlc_d, wlc_p)
#if MPI_VERSION
        use mpi
#endif
        implicit none
        type(wlcsim_data), intent(inout)   :: wlc_d
        type(wlcsim_params), intent(in)    :: wlc_p
        character(8) datedum  ! trash
        character(10) timedum ! trash
        character(5) zonedum  ! trash
        integer seedvalues(8) ! clock readings
        integer NT  ! total number of beads
        integer NBin ! total number of bins
        integer i
        integer irand
        integer ( kind = 4 ) dest   !destination id for messages
        integer ( kind = 4 ) source  !source id for messages
        integer ( kind = 4 ) status(MPI_status_SIZE) ! MPI stuff
        integer ( kind = 4 ) error  ! error id for MIP functions
        nt = wlc_p%nt
        nbin = wlc_p%nbin

#if MPI_VERSION
        call init_MPI(wlc_d)
#endif
        allocate(wlc_d%R(3,NT))
        allocate(wlc_d%U(3,NT))
        if (wlc_p%codeName /= 'bruno' .OR. wlc_p%nInitMCSteps /= 0) then
            allocate(wlc_d%RP(NT,3))
            allocate(wlc_d%UP(NT,3))
        endif
        !TOdo these should in principle be inside the following if statement,
        !but it's not clear if that's possible without adding a bunch of dirty
        !if statements deep inside mc_move. which is fine, but I would want to
        !check with quinn *exactly* in which cases they're needed if i do that
        if (wlc_p%field_int_on) then
            allocate(wlc_d%AB(NT))   !Chemical identity aka binding state
            allocate(wlc_d%ABP(NT))   !Chemical identity aka binding state
            allocate(wlc_d%PHIA(NBin))
            allocate(wlc_d%PHIB(NBin))
            allocate(wlc_d%DPHIA(NBin))
            allocate(wlc_d%DPHIB(NBin))
            allocate(wlc_d%indPHI(NBin))
            allocate(wlc_d%PhiH(NBin))
            allocate(wlc_d%Vol(NBin))
            do I = 1,NBin
                wlc_d%PHIA(I) = 0.0_dp
                wlc_d%PHIB(I) = 0.0_dp
            enddo
        endif
        if (wlc_p%bind_on) then
            allocate(wlc_d%METH(NT)) !Underlying methalation profile
        endif
        !Allocate vector of writhe and elastic energies for replicas
        if (wlc_p%pt_twist) then
            allocate(wlc_d%Wrs(wlc_d%nLKs))
            allocate(wlc_d%eelasREPLICAS(wlc_d%nLKs,4))

        endif
        if (wlc_p%ring) then !TOdo this should be if ("knot")
            wlc_d%NCross = 0
            wlc_d%CrossSize = wlc_p%NB**2
            allocate(wlc_d%Cross(wlc_d%CrossSize,6))
            allocate(wlc_d%CrossP(wlc_d%CrossSize,6))
        endif
        !If parallel tempering is on, initialize the nodeNumbers

        if (wlc_p%pt_twist) then

            !Allocate node numbers
            allocate(wlc_d%nodeNUMBER(wlc_d%nLKs))
            do i = 1,wlc_d%nLKs
                wlc_d%nodeNUMBER(i) = i
            enddo

            !Initially, replica start and replica end are the first and second to last replicas for even
            !nLKs and the first and second to last for odd nLKs
            if (mod(wlc_d%nLKs,2).eq.0) then
                wlc_d%replicaSTART = 1
                wlc_d%replicaEND = wlc_d%nLKs - 1
            else
                wlc_d%replicaSTART = 1
                wlc_d%replicaEND = wlc_d%nLKs - 2
            endif

            !Allocate the number of replica exchange trials and successes and initialize to zero
            allocate(wlc_d%nSWAPup(wlc_d%nLKs))
            allocate(wlc_d%nSWAPdown(wlc_d%nLKs))
            allocate(wlc_d%nTRIALup(wlc_d%nLKs))
            allocate(wlc_d%nTRIALdown(wlc_d%nLKs))

            wlc_d%nSWAPup = 0
            wlc_d%nSWAPdown = 0
            wlc_d%nTRIALup = 0
            wlc_d%nTRIALdown = 0

        endif

        if (wlc_p%collisionDetectionType /= 0) then
            allocate(wlc_d%coltimes(NT,NT))
            wlc_d%coltimes = -1.0_dp
        endif

#if MPI_VERSION
        ! -----------------------------------------------
        !
        !   Generate thread safe random number seeds
        !
        !--------------------------------------------
        if (wlc_d%id .eq. 0) then ! head node
            if (.false.) then ! set spedific seed
                Irand = 7171
            else ! seed from clock
                call date_and_time(datedum,timedum,zonedum,seedvalues)
                Irand = int(-seedvalues(5)*1E7-seedvalues(6)*1E5 &
                          -seedvalues(7)*1E3-seedvalues(8))
                Irand = mod(Irand,10000)
                print*, "Random Intiger seed:",Irand
            endif
            call random_setseed(Irand*(wlc_d%id + 1),wlc_d%rand_stat) ! random seed for head node
            do dest = 1,wlc_d%numProcesses-1 ! send out the others
                call MPI_Send (Irand,1, MPI_integer, dest,   0, &
                                MPI_COMM_WORLD,error )
            enddo
        else ! worker node
            source = 0
            call MPI_Recv ( Irand, 1, MPI_integer, source, 0, &
                            MPI_COMM_WORLD, status, error )
            call random_setseed(Irand*(wlc_d%id + 1),wlc_d%rand_stat)
            !if (wlc_d%restart) then
            !    call pt_restart(wlc_p,wlc_d)
            !endif
        endif
#else
        if (.false.) then ! if you wanted to set specific seed
            wlc_d%rand_seed = 7171
        else ! seed from clock
            call date_and_time(datedum,timedum,zonedum,seedvalues)
            ! funny business
            wlc_d%rand_seed = int(-seedvalues(5)*1E7-seedvalues(6)*1E5 &
                      -seedvalues(7)*1E3-seedvalues(8))
            wlc_d%rand_seed = mod(wlc_d%rand_seed,10000)
            ! print*, "Random Intiger seed:",wlc_d%rand_seed
        endif

        call random_setseed(wlc_d%rand_seed, wlc_d%rand_stat)
#endif
        call initcond(wlc_d%R, wlc_d%U, wlc_p%NT, wlc_p%NB, &
            wlc_p%NP, wlc_p%frmfile, pack_as_para(wlc_p), wlc_p%lbox, &
            wlc_p%initCondType, wlc_d%rand_stat, wlc_p%ring, wlc_p)

        if (wlc_p%field_int_on) then
            call initchem(wlc_d%AB, wlc_p%nT, wlc_p%nMpP, wlc_p%nBpM, wlc_p%nP, wlc_p%fA, wlc_p%lam, wlc_d%rand_stat)
            ! calculate volumes of bins
            if (wlc_p%confineType.eq.3) then
                call MC_calcVolume(wlc_p%confinetype, wlc_p%NBinX, wlc_p%dBin, &
                                wlc_p%LBox(1), wlc_d%Vol, wlc_d%rand_stat)
            else
                do I = 1,NBin
                    wlc_d%Vol(I) = wlc_p%dBin**3
                enddo
            endif
        endif

        if (wlc_p%bind_on) then
            call initchem(wlc_d%meth, wlc_p%nt, wlc_p%nMpP, wlc_p%nBpM, wlc_p%nP, wlc_p%fA, wlc_p%lam, wlc_d%rand_stat)
        endif

        ! initialize energies
!        call CalculateEnergiesFromScratch(wlc_p,wlc_d)
        wlc_d%EElas   = wlc_d%dEElas
        if (wlc_p%field_int_on) then
            wlc_d%ECouple =wlc_d%dECouple
            wlc_d%EKap    =wlc_d%dEKap
            wlc_d%ECHI    =wlc_d%dECHI
            wlc_d%EField  =wlc_d%dEField
            wlc_d%x_Field =wlc_d%dx_Field
            wlc_d%x_couple = wlc_d%dx_couple
            wlc_d%x_Kap   =wlc_d%dx_Kap
            wlc_d%x_Chi   =wlc_d%dx_Chi
        else
            wlc_d%ECouple =0.0_dp
            wlc_d%EKap    =0.0_dp
            wlc_d%ECHI    =0.0_dp
            wlc_d%EField  =0.0_dp
            wlc_d%x_Field =0.0_dp
            wlc_d%x_couple = 0.0_dp
            wlc_d%x_Kap   =0.0_dp
            wlc_d%x_Chi   =0.0_dp
        endif
        if (wlc_p%bind_On) then
            wlc_d%ebind   =wlc_d%debind
            wlc_d%x_mu    =wlc_d%dx_mu
        else
            wlc_d%ebind   =0.0_dp
            wlc_d%x_mu    =0.0_dp
        endif
        if(wlc_p%Ring) then
            wlc_d%eKnot   =1.0
        endif

        wlc_d%time = 0
        wlc_d%time_ind = 0
        wlc_d%mc_ind = 0



    end subroutine initialize_wlcsim_data

    function pack_as_para(wlc_p) result(para)
        implicit none
        type(wlcsim_params), intent(in) :: wlc_p
        real(dp) para(10)
        para(1) = wlc_p%EB
        para(2) = wlc_p%EPAR
        para(3) = wlc_p%EPERP
        para(4) = wlc_p%GAM
        para(5) = wlc_p%ETA
        para(6) = wlc_p%XIR
        para(7) = wlc_p%XIU
        para(8) = wlc_p%LBOX(1)
        para(9) = wlc_p%lhc
        para(10) = wlc_p%VHC
    end function pack_as_para


    subroutine printDescription(wlc_p)
        implicit none
        type(wlcsim_params), intent(in) :: wlc_p
        print*, "---------------System Description---------------"
        print*, " type of simulation, codeName", wlc_p%codeName
        print*, " WLC, DSSWLC, GC, simType", wlc_p%simType
        print*, "Bead variables:"
        print*, " Total number of beads, NT = ", wlc_p%NT
        print*, " Number of beads in a polymer, NB = ", wlc_p%NB
        print*, " Number of monomers in a polymer, nMpP = ", wlc_p%nMpP
        print*, " Number of polymers, NP = ",wlc_p%NP
        print*, " Number of beads in a monomer, nBpM = ", wlc_p%nBpM
        print*, " fraction Methalated", wlc_p%F_METH
        print*, " LAM_METH", wlc_p%LAM_METH
        print*, " "
        print*, "Length and volume Variables:"
        print*, " persistance length =",(wlc_p%L0/(2.0_dp*wlc_p%EPS))
        print*, " length of each polymer in simulation, l = ",wlc_p%l
        print*, " twist persistence length, lt", wlc_p%lt
        print*, " lbox = ", wlc_p%lbox(1), wlc_p%lbox(2), wlc_p%lbox(3)
        print*, " Number of bins in x direction", &
                   wlc_p%NBinX(1), wlc_p%NBinX(2),wlc_p%NBinX(3)
        print*, " Number of bins", wlc_p%NBin
        print*, " spatial descritation dbin = ",wlc_p%dbin
        print*, " L0 = ", wlc_p%L0
        print*, " volume fraction polymer =", wlc_p%Fpoly
        print*, " bead volume V = ", wlc_p%beadVolume
        print*, " number of kuhn lengths between beads, eps ", wlc_p%eps
        print*, " "
        print*, "Energy Variables"
        print*, " elasticity EPS =", wlc_p%EPS
        print*, " solvent-polymer CHI =",wlc_p%CHI
        print*, " compression cof, KAP =", wlc_p%KAP
        print*, " field strength, hA =", wlc_p%hA
        print*, " -energy of binding unmethalated ", wlc_p%EU," more positive for favorable binding"
        print*, " -energy of binding methalated",wlc_p%EM
        print*, " HP1_Binding energy parameter", wlc_p%HP1_Bind
        print*, " chemical potential of HP1, mu", wlc_p%mu
        print*, " bend-shear coupling parameter, eta ", wlc_p%eta
        print*, " "
        print*, "Time Variables"
        print*, " stepsPerExchange", wlc_p%stepsPerExchange
        print*, " nReplicaExchangePerSavePoint", wlc_p%nReplicaExchangePerSavePoint
        print*, " numSavePoints", wlc_p%numSavePoints
        print*, " stepsPerSave", wlc_p%stepsPerSave
        print*, " "
        print*, "Switches:"
        print*, " confinetype:",wlc_p%confinetype
        print*, " initCondType:",wlc_p%initCondType
        print*, " ring:", wlc_p%ring
        print*, " twist:", wlc_p%twist
        print*, " "
        print*, "---------------------------------------------"

    end subroutine

    subroutine tweak_param_defaults(wlc_p, wlc_d)
        implicit none
        type(wlcsim_params), intent(inout) :: wlc_p
        type(wlcsim_data), intent(inout) :: wlc_d
        real(dp) :: default_window

        if (wlc_p%dbin /= wlc_p%dbin) then
            ! discretizing at 1 persistence length seems to be a reasonable default
            wlc_p%dbin = wlc_p%lp
        endif

        if (wlc_p%confineType.eq.3) then
            wlc_p%fPoly = 6.0_dp*wlc_p%beadVolume*wlc_p%NT &
                /PI/wlc_p%LBOX(1)/wlc_p%LBOX(1)/wlc_p%LBOX(1)
        else
            wlc_p%fPoly = wlc_p%beadVolume*wlc_p%NT &
                /wlc_p%LBOX(1)/wlc_p%LBOX(2)/wlc_p%LBOX(3)
        endif

        wlc_p%L0 = wlc_p%l/real(wlc_p%nb)
        !  Edit the following to optimize wlc_p performance
        !  Monte-Carlo simulation parameters
        wlc_d%MCAMP(1) = 0.5_dp*PI
        wlc_d%MCAMP(2) = 0.3_dp*wlc_p%L0
        wlc_d%MCAMP(3) = 0.5_dp*PI
        wlc_d%MCAMP(4) = 0.5_dp*PI
        wlc_d%MCAMP(5) = 0.5_dp*PI
        wlc_d%MCAMP(6) = 5.0_dp*wlc_p%L0
        wlc_d%MCAMP(7) = nan
        wlc_d%MCAMP(8) = nan
        wlc_d%MCAMP(9) = nan
        wlc_d%MCAMP(10) = nan

        !switches to turn on various types of moves
        wlc_p%MOVEON(1) = 1  ! crank-shaft move
        wlc_p%MOVEON(2) = 1  ! slide move
        wlc_p%MOVEON(3) = 1  ! pivot move
        wlc_p%MOVEON(4) = 1  ! rotate move
        wlc_p%MOVEON(5) = 1  ! full chain rotation
        wlc_p%MOVEON(6) = 1  ! full chain slide
        ! if we're not using field interactions
        ! energies, then this should never be on
        if (wlc_p%field_int_on) then
            wlc_p%MOVEON(7) = 1  ! Change in Binding state
        else
            wlc_p%MOVEON(7) = 0
        endif
        wlc_p%MOVEON(8) = 0  ! Chain flip ! TOdo not working
        ! if number of polymers is 1, this should never be on
        if (wlc_p%np < 2) then
            wlc_p%moveon(9) = 0
        else
            wlc_p%MOVEON(9) = 1  ! Chain exchange
        endif
        wlc_p%MOVEON(10) = 1 ! Reptation
        if (wlc_p%codeName == 'quinn') then
            !switches to turn on various types of moves
            wlc_p%MOVEON(1) = 1  ! crank-shaft move
            wlc_p%MOVEON(2) = 1  ! slide move
            wlc_p%MOVEON(3) = 1  ! pivot move
            wlc_p%MOVEON(4) = 1  ! rotate move
            wlc_p%MOVEON(5) = 0  ! full chain rotation
            wlc_p%MOVEON(6) = 0  ! full chain slide
            wlc_p%MOVEON(7) = 1  ! Change in Binding state
            wlc_p%MOVEON(8) = 0  ! Chain flip
            wlc_p%MOVEON(9) = 0  ! Chain exchange
            wlc_p%MOVEON(10) = 0 ! Reptation
        endif

        !     Initial segment window for wlc_p moves
        default_window = wlc_p%nB
        default_window = max(default_window, 1.0_dp*wlc_p%nMpP*wlc_p%nBpM)
        default_window = default_window/5.0_dp
        wlc_d%Window(1) = default_window ! 15.0_dp ! used to be N*G
        wlc_d%Window(2) = default_window ! 15.0_dp ! used to be N*G
        wlc_d%Window(3) = default_window ! 15.0_dp ! used to be N*G
        wlc_d%Window(4) = default_window ! 1.0_dp
        wlc_d%Window(5) = default_window ! dble(wlc_p%nMpP*wlc_p%nBpM)
        wlc_d%Window(6) = default_window ! dble(wlc_p%nMpP*wlc_p%nBpM)
        wlc_d%Window(7) = default_window ! 15.0_dp ! used to be N*G
        wlc_d%Window(8) = default_window ! dble(wlc_p%nMpP*wlc_p%nBpM)
        wlc_d%Window(9) = default_window ! dble(wlc_p%nMpP*wlc_p%nBpM)
        wlc_d%Window(9) = 1.0_dp

        !    Maximum window size (large windows are expensive)
        wlc_p%MAXWindoW(1) = dble(min(150,wlc_p%NB))
        wlc_p%MAXWindoW(2) = dble(min(150,wlc_p%NB))
        wlc_p%MAXWindoW(3) = dble(min(150,wlc_p%NB))
        wlc_p%MAXWindoW(4) = nan
        wlc_p%MAXWindoW(5) = nan
        wlc_p%MAXWindoW(6) = nan
        wlc_p%MAXWindoW(7) = dble(min(4,wlc_p%NB))
        wlc_p%MAXWindoW(8) = nan
        wlc_p%MAXWindoW(9) = nan
        wlc_p%MAXWindoW(9) = nan ! need to chaige code to allow >1


        if (wlc_p%MinWindow(1).ne.wlc_p%MinWindow(1)) wlc_p%MinWindoW(1) = dble(min(10,wlc_p%NB))
        if (wlc_p%MinWindow(2).ne.wlc_p%MinWindow(2)) wlc_p%MinWindoW(2) = dble(min(10,wlc_p%NB))
        if (wlc_p%MinWindow(3).ne.wlc_p%MinWindow(3)) wlc_p%MinWindoW(3) = dble(min(10,wlc_p%NB))
        if (wlc_p%MinWindow(7).ne.wlc_p%MinWindow(7)) wlc_p%MinWindoW(7) = dble(min(10,wlc_p%NB))

        wlc_p%MinAMP(1) = 0.1_dp*PI
        wlc_p%MinAMP(2) = 0.2_dp*wlc_p%L0
        wlc_p%MinAMP(3) = 0.2_dp*PI
        wlc_p%MinAMP(4) = 0.2_dp*PI
        wlc_p%MinAMP(5) = 0.05_dp*PI
        wlc_p%MinAMP(6) = 0.2_dp*wlc_p%L0
        wlc_p%MinAMP(7) = nan
        wlc_p%MinAMP(8) = nan
        wlc_p%MinAMP(9) = nan
        wlc_p%MinAMP(10) = nan

        wlc_p%MAXAMP(1) = 1.0_dp*PI
        wlc_p%MAXAMP(2) = 1.0_dp*wlc_p%L0
        wlc_p%MAXAMP(3) = 1.0_dp*PI
        wlc_p%MAXAMP(4) = 1.0_dp*PI
        wlc_p%MAXAMP(5) = 1.0_dp*PI
        wlc_p%MAXAMP(6) = 0.1*wlc_p%lbox(1)
        wlc_p%MAXAMP(7) = nan
        wlc_p%MAXAMP(8) = nan
        wlc_p%MAXAMP(9) = nan
        wlc_p%MAXAMP(10) = nan

        ! Solution
        wlc_p%NBinX(1) = nint(wlc_p%LBOX(1)/wlc_p%dBin)
        wlc_p%NBinX(2) = nint(wlc_p%LBOX(2)/wlc_p%dBin)
        wlc_p%NBinX(3) = nint(wlc_p%LBOX(3)/wlc_p%dBin)
        wlc_p%NBin = wlc_p%NBinX(1)*wlc_p%NBinX(2)*wlc_p%NBinX(3)
        if (abs(wlc_p%dBin*sum(wlc_p%NBinX) /sum(wlc_p%LBOX) -1)>0.001) then
            print*, "Warning: I have changed the descritation lenth from:", wlc_p%dBin
            wlc_p%dBin = wlc_p%lBox(1)*dble(wlc_p%NBinX(1))
            print*, "to ", wlc_p%dBin
        endif

        if (wlc_p%codeName == 'brad') then
            ! initialize windows to number of beads
            wlc_p%MAXWindoW = wlc_p%nB         ! Max Size of window for bead selection
            wlc_p% MinWindoW  = 1         ! Min Size of window for bead selection

            ! Window amplitudes
            wlc_p%MinAMP = 0.0_dp ! minium amplitude
            wlc_p%MinAMP(1) = 0.07_dp*pi
            wlc_p%MinAMP(2) = 0.01_dp*wlc_p%l/wlc_p%nB
            wlc_p%MaxAMP = 2.0_dp*pi
            wlc_p%MaxAMP(2) = wlc_p%lbox(1)
            wlc_p%MaxAMP(6) = wlc_p%lbox(1)

            ! Turn off saving AB
            wlc_p%saveAB = .FALSE.

            !Set which moves are used
            wlc_p%MOVEON(1) = 1  ! crank-shaft move
            wlc_p%MOVEON(2) = 1  ! slide move
            wlc_p%MOVEON(3) = 1  ! pivot move
            wlc_p%MOVEON(4) = 1  ! rotate move
            wlc_p%MOVEON(5) = 0  ! full chain rotation
            wlc_p%MOVEON(6) = 0  ! full chain slide
            wlc_p%MOVEON(7) = 0  ! Change in Binding state
            wlc_p%MOVEON(8) = 0  ! Chain flip
            wlc_p%MOVEON(9) = 0  ! Chain exchange
            wlc_p%MOVEON(10) = 0 ! Reptation

        endif

        ! If ring is on, turn off the pivot move
        if (wlc_p%ring) then
            wlc_p%moveON(3) = 0
        endif

    end subroutine

    subroutine wlcsim_params_recenter(wlc_p,wlc_d)
    !  Prevents drift in periodic BC
        implicit none
        type(wlcsim_params), intent(in) :: wlc_p
        type(wlcsim_data), intent(inout) :: wlc_d
        integer IB, I, J   ! Couners
        real(dp) R0(3)  ! Offset to move by
        IB = 1
        do I = 1,wlc_p%NP
        R0(1) = nint(wlc_d%R(1,IB)/wlc_p%lbox(1)-0.5_dp)*wlc_p%lbox(1)
        R0(2) = nint(wlc_d%R(2,IB)/wlc_p%lbox(2)-0.5_dp)*wlc_p%lbox(2)
        R0(3) = nint(wlc_d%R(3,IB)/wlc_p%lbox(3)-0.5_dp)*wlc_p%lbox(3)
        if (abs(R0(1)*R0(2)*R0(3)) .gt. 0.0001_dp) then
            do J = 1,wlc_p%NB
                wlc_d%R(1,IB) = wlc_d%R(1,IB)-R0(1)
                wlc_d%R(2,IB) = wlc_d%R(2,IB)-R0(2)
                wlc_d%R(3,IB) = wlc_d%R(3,IB)-R0(3)
                IB = IB + 1
            enddo
        endif
        enddo
    end subroutine

    subroutine printSimInfo(i, wlc_p, wlc_d)
    ! print out current simulation metainformation
        implicit none
        type(wlcsim_params), intent(in) :: wlc_p
        type(wlcsim_data), intent(in) :: wlc_d
        integer, intent(in) :: i
        print*, 'Current time ', wlc_d%time
        print*, 'Time point ', wlc_d%time_ind, ' out of ', wlc_p%stepsPerSave*wlc_p%numSavePoints
        print*, 'Save point ', i, ' out of ', wlc_p%numSavePoints
    end subroutine

    subroutine printEnergies(wlc_d)
    ! For realtime feedback on wlc_p simulation
        implicit none
        type(wlcsim_data), intent(in) :: wlc_d
        print*, "ECouple:", wlc_d%ECouple
        print*, "Bending energy", wlc_d%EELAS(1)
        print*, "Par compression energy", wlc_d%EELAS(2)
        print*, "Shear energy", wlc_d%EELAS(3)
        print*, "ECHI", wlc_d%ECHI
        print*, "EField", wlc_d%EField
        print*, "EKAP", wlc_d%EKAP
        print*, "ebind", wlc_d%ebind
    end subroutine

    subroutine wlcsim_params_printPhi(wlc_p,wlc_d)
    ! prints densities for trouble shooting
        implicit none
        type(wlcsim_params), intent(in) :: wlc_p
        type(wlcsim_data), intent(in) :: wlc_d
        integer I
        real(dp) EKap, ECouple, EChi,VV, PHIPOly
        print*,"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        print*, " PHIA  | PHIB  | PPoly |  Vol  | EKap  | EChi  |ECouple|"
        do I = 1,wlc_p%NBin
            VV = wlc_d%Vol(I)
            if (VV.le.0.1_dp) cycle
            PHIPOLY = wlc_d%PHIA(I) + wlc_d%PHIB(I)
            EChi = VV*(wlc_p%CHI/wlc_p%beadVolume)*PHIPoly*(1.0_dp-PHIPoly)
            ECouple = VV*wlc_p%HP1_Bind*(wlc_d%PHIA(I))**2
            if(PHIPoly > 1.0_dp) then
            EKap = VV*(wlc_p%KAP/wlc_p%beadVolume)*(PHIPoly-1.0_dp)**2
            else
            cycle
            EKap = 0.0_dp
            endif
            write(*,"(4f8.4,3f8.1)") wlc_d%PHIA(I), wlc_d%PHIB(I), &
                                wlc_d%PHIA(I) + wlc_d%PHIB(I),wlc_d%Vol(I),&
                                EKap,EChi,ECouple
        enddo
        print*,"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    end subroutine

    subroutine printWindowStats(wlc_p, wlc_d)
    ! For realtime feedback on adaptation
        implicit none
        type(wlcsim_params), intent(in) :: wlc_p
        type(wlcsim_data), intent(in) :: wlc_d
        integer I ! counter
        I = 0
        print*, "Succes | MCAMP | WindoW| type "
        do I = 1,wlc_p%movetypes
            if (wlc_p%MOVEON(i).eq.1) then
                write(*,"(f8.5,2f8.2,1I8)") wlc_d%phit(i), wlc_d%MCAMP(i),  wlc_d%WindoW(i), i
            endif
        enddo
        return
    end subroutine

    subroutine wlcsim_params_LoadField(wlc_p,wlc_d,fileName)
        implicit none
        type(wlcsim_params), intent(in) :: wlc_p
        type(wlcsim_data), intent(inout) :: wlc_d
        integer I
        character(MAXFILENAMELEN) fileName ! file name to load from
        open (unit = inFileUnit, file = fileName, status = 'OLD')
        do I = 1,wlc_p%NBin
            read(inFileUnit,*) wlc_d%PHIH(I)
        enddo
        return
    end subroutine

    subroutine wlcsim_params_MakeField(wlc_p,wlc_d)
        implicit none
        type(wlcsim_params), intent(in) :: wlc_p
        type(wlcsim_data), intent(inout) :: wlc_d
        integer indBin  ! index of bin
        integer IX,IY,IZ ! bin corrdinates

        do IX = 1,wlc_p%NBinX(1)
            do IY = 1,wlc_p%NBinX(2)
                do IZ = 1,wlc_p%NBinX(3)
                    indBin = IX + &
                        (IY-1)*wlc_p%NBinX(1) + &
                        (IZ-1)*wlc_p%NBinX(1)*wlc_p%NBinX(2)
                    wlc_d%PHIH(indBin) = dsin(wlc_p%k_field*wlc_p%dbin*dble(IX))
                enddo
            enddo
        enddo
        return
    end subroutine

    subroutine wlcsim_params_loadAB(wlc_p,wlc_d,fileName)
    ! Loads AB for file...has not been tested
        implicit none
        type(wlcsim_params), intent(in) :: wlc_p
        type(wlcsim_data), intent(inout) :: wlc_d
        character(MAXFILENAMELEN), intent(in) :: fileName ! file name to load from
        integer IB, I, J ! counters
        open (unit = inFileUnit, file = fileName, status = 'OLD')
        IB = 1
        do I = 1,wlc_p%NP
        do J = 1,wlc_p%NB
            read(inFileUnit,"(I2)") wlc_d%AB(IB)
            IB = IB + 1
            enddo
        enddo
        close(inFileUnit)
    end subroutine

    subroutine wlcsim_params_saveR(wlc_p,wlc_d,fileName,repeatingBC,stat)
    ! Writes R and AB to file for analysis
    ! Rx  Ry  Rz AB
        implicit none
        logical, intent(in) :: repeatingBC  ! 1 for reapeating boundary conditions
        integer I,J,IB  ! counters
        type(wlcsim_params), intent(in) :: wlc_p
        type(wlcsim_data), intent(in) :: wlc_d
        character(MAXFILENAMELEN), intent(in) :: fileName
        character(MAXFILENAMELEN) fullName
        character(len = *), intent(in) :: stat
        fullName = trim(fileName) // trim(wlc_d%repSuffix)
        fullName = trim(fullName)
        open (unit = outFileUnit, file = fullName, status = stat)
        IB = 1
        if (repeatingBC) then
           do I = 1,wlc_p%NP
              do J = 1,wlc_p%NB
                 if (wlc_p%saveAB) then
                    write(outFileUnit,"(3f10.3,I2)") &
                         wlc_d%R(1,IB)-0.*nint(wlc_d%R(1,IB)/wlc_p%lbox(1)-0.5_dp)*wlc_p%lbox(1), &
                         wlc_d%R(2,IB)-0.*nint(wlc_d%R(2,IB)/wlc_p%lbox(2)-0.5_dp)*wlc_p%lbox(2), &
                         wlc_d%R(3,IB)-0.*nint(wlc_d%R(3,IB)/wlc_p%lbox(3)-0.5_dp)*wlc_p%lbox(3), &
                         wlc_d%AB(IB)
                 else
                    write(outFileUnit,"(3f10.3)") &
                         wlc_d%R(1,IB)-0.*nint(wlc_d%R(1,IB)/wlc_p%lbox(1)-0.5_dp)*wlc_p%lbox(1), &
                         wlc_d%R(2,IB)-0.*nint(wlc_d%R(2,IB)/wlc_p%lbox(2)-0.5_dp)*wlc_p%lbox(2), &
                         wlc_d%R(3,IB)-0.*nint(wlc_d%R(3,IB)/wlc_p%lbox(3)-0.5_dp)*wlc_p%lbox(3)
                 endif
                 IB = IB + 1
              enddo
           enddo
           print*, "Error in wlcsim_params_saveR"
           print*, "Are you sure you want repeating BC?"
           print*, "Quinn put this in ages ago but never implemented it...."
           stop 1
        else
           do I = 1,wlc_p%NP
              do J = 1,wlc_p%NB
                  if (wlc_p%saveAB) then
                     write(outFileUnit,"(3f10.3,I2)") &
                            wlc_d%R(1,IB),wlc_d%R(2,IB),wlc_d%R(3,IB),wlc_d%AB(IB)
                  else
                     write(outFileUnit,"(3f10.3)") &
                           wlc_d%R(1,IB),wlc_d%R(2,IB),wlc_d%R(3,IB)
                  endif
                  IB = IB + 1
              enddo
           enddo
        endif
        close(outFileUnit)

      end subroutine wlcsim_params_saveR

    subroutine wlcsim_params_savePHI(wlc_p,wlc_d,fileName)
    ! Saves PHIA and PHIB to file for analysis
        implicit none
        integer I  ! counters
        type(wlcsim_params), intent(in) :: wlc_p
        type(wlcsim_data), intent(in) :: wlc_d
        character(MAXFILENAMELEN), intent(in) :: fileName
        character(MAXFILENAMELEN) fullName
        fullName=  trim(fileName) // trim(wlc_d%repSuffix)
        open (unit = outFileUnit, file = fullName, status = 'NEW')
        do I = 1,wlc_p%NBin
            write(outFileUnit,"(2f7.2)") wlc_d%PHIA(I),wlc_d%PHIB(I)
        enddo
        close(outFileUnit)
    end subroutine

    subroutine wlcsim_params_saveU(wlc_p,wlc_d,fileName,stat)
    ! Saves U to ASCII file for analisys
        implicit none
        integer I,J,IB  ! counters
        type(wlcsim_params), intent(in) :: wlc_p
        type(wlcsim_data), intent(in) :: wlc_d
        character(MAXFILENAMELEN), intent(in) :: fileName
        character(MAXFILENAMELEN) fullName
        character(len = *), intent(in) :: stat
        fullName=  trim(fileName) // trim(wlc_d%repSuffix)
        open (unit = outFileUnit, file = fullName, status = stat)
        IB = 1
        do I = 1,wlc_p%NP
            do J = 1,wlc_p%NB
                write(outFileUnit,"(3f8.3,2I2)") wlc_d%U(1,IB),wlc_d%U(2,IB),wlc_d%U(3,IB)
                IB = IB + 1
            enddo
        enddo
        close(outFileUnit)
    end subroutine

    subroutine save_parameters(wlc_p,fileName)
        ! Write a number of parameters ASCII variables to file for reccords
        implicit none
        type(wlcsim_params), intent(in) :: wlc_p
        character(len=*), intent(in) :: fileName
        open (unit =outFileUnit, file = fileName, status = 'NEW')
            write(outFileUnit,"(I8)") wlc_p%NT ! 1 Number of beads in simulation
            write(outFileUnit,"(I8)") wlc_p%nMpP  ! 2 Number of monomers in a polymer
            write(outFileUnit,"(I8)") wlc_p%NB ! 3 Number of beads in a polymer
            write(outFileUnit,"(I8)") wlc_p%NP ! 4 Number of polymers in simulation
            write(outFileUnit,"(I8)") wlc_p%NT ! 5 Number of beads in simulation
            write(outFileUnit,"(I8)") wlc_p%nBpM  ! 6 Number of beads per monomer

            write(outFileUnit,"(f10.5)") wlc_p%L0    ! Equilibrium segment length
            write(outFileUnit,"(f10.5)") wlc_p%CHI  ! 8  initail CHI parameter value
            write(outFileUnit,"(f10.5)") wlc_p%Fpoly ! Fraction polymer
            write(outFileUnit,"(f10.5)") wlc_p%lbox(1)  ! 10 Lenth of box
            write(outFileUnit,"(f10.5)") wlc_p%EU    ! Energy unmethalated
            write(outFileUnit,"(f10.5)") wlc_p%EM    ! 12 Energy methalated
            write(outFileUnit,"(f10.5)") wlc_p%HP1_Bind ! Energy of HP1 binding
            write(outFileUnit,"(f10.5)") (wlc_p%L0/wlc_p%EPS) ! 14 Khun lenth
            write(outFileUnit,"(A)") "-999"  ! for historic reasons
            write(outFileUnit,"(f10.5)") wlc_p%F_METH  ! methalation fraction
            write(outFileUnit,"(f10.5)") wlc_p%LAM_METH  ! methalation lambda
        close(outFileUnit)
    end subroutine

    subroutine wlcsim_params_appendEnergyData(save_ind, wlc_p, wlc_d, fileName)
    ! print Energy data
        implicit none
        type(wlcsim_params), intent(in) :: wlc_p
        type(wlcsim_data), intent(in) :: wlc_d
        integer, intent(in) :: save_ind
        character(MAXFILENAMELEN), intent(in) :: fileName
        LOGICAL isfile
        character(MAXFILENAMELEN) fullName
        fullName=  trim(fileName) // trim(wlc_d%repSuffix)
        inquire(file = fullName, exist = isfile)
        if (isfile) then
            open (unit = outFileUnit, file = fullName, status ='OLD', POSITION = "append")
        else
            open (unit = outFileUnit, file = fullName, status = 'new')
            write(outFileUnit,*) "ind | id |",&
                       " ebend  | eparll | EShear | ECoupl | E Kap  | E Chi  |",&
                       " EField | ebind  |  x_Mu  | Couple |  Chi   |  mu    |",&
                       "  Kap   | Field  |"
        endif
        write(outFileUnit,"(2I5, 9f9.1,5f9.4)") save_ind, wlc_d%id, &
            wlc_d%EELAS(1), wlc_d%EELAS(2), wlc_d%EELAS(3), wlc_d%ECouple, &
            wlc_d%EKap, wlc_d%ECHI, wlc_d%EField, wlc_d%ebind, wlc_d%x_Mu, &
            wlc_p%HP1_Bind*wlc_p%Couple_on, wlc_p%CHI*wlc_p%CHI_ON, wlc_p%mu, wlc_p%KAP*wlc_p%KAP_ON,&
            wlc_p%hA
        close(outFileUnit)
    end subroutine

    subroutine wlcsim_params_appendAdaptData(save_ind, wlc_p, wlc_d, fileName)
    ! Appends wlc_p move adaptation data to the file
        implicit none
        type(wlcsim_params), intent(in) :: wlc_p
        type(wlcsim_data), intent(in) :: wlc_d
        integer, intent(in) :: save_ind
        LOGICAL isfile
        character(MAXFILENAMELEN), intent(in) :: fileName
        character(MAXFILENAMELEN) fullName
        fullName=  trim(fileName) // trim(wlc_d%repSuffix)
        inquire(file = fullName, exist = isfile)
        if (isfile) then
            open (unit = outFileUnit, file = fullName, status ='OLD', POSITION = "append")
        else
            open (unit = outFileUnit, file = fullName, status = 'new')
            write(outFileUnit,*) "ind| id|",&
                       " Win 1 | AMP 1 | SUC 1 | Win 2 | AMP 2 | SUC 2 |",&
                       " Win 3 | AMP 3 | SUC 3 | ON  4 | AMP 4 | SUC 4 |",&
                       " ON  5 | AMP 5 | SUC 5 | ON  6 | AMP 6 | SUC 6 |",&
                       " ON  7 | SUC 7 | ON  8 | SUC 8 |", &
                       " ON  9 | SUC 9 | ON 10 | SUC 10|"
        endif
        write(outFileUnit,"(2I4,26f8.3)") save_ind,wlc_d%id,&
            real(wlc_d%WindoW(1)),wlc_d%MCAMP(1),wlc_d%PHIT(1), &
            real(wlc_d%WindoW(2)),wlc_d%MCAMP(2),wlc_d%PHIT(2), &
            real(wlc_d%WindoW(3)),wlc_d%MCAMP(3),wlc_d%PHIT(3), &
            real(wlc_p%MOVEON(4)),wlc_d%MCAMP(4),wlc_d%PHIT(4), &
            real(wlc_p%MOVEON(5)),wlc_d%MCAMP(5),wlc_d%PHIT(5), &
            real(wlc_p%MOVEON(6)),wlc_d%MCAMP(6),wlc_d%PHIT(6), &
            real(wlc_p%MOVEON(7)),wlc_d%PHIT(7), &
            real(wlc_p%MOVEON(8)),wlc_d%PHIT(8), &
            real(wlc_p%MOVEON(9)),wlc_d%PHIT(9), &
            real(wlc_p%MOVEON(10)),wlc_d%PHIT(10)
        close(outFileUnit)
    end subroutine
    subroutine wlcsim_params_writebinary(wlc_p,wlc_d,baseName)
    !    This function writes the contence of the structures wlc_p and wlc_d
    !  to a binary file.  if you add more variables to wlc_d you need to
    !  a seperate write command for them as it is not possible to write
    !  a structure with allocatables to a binar file.
    !    The contence are stored in
    !     baseName//'R'
    !     baseName//'U'
    !     etc.
        implicit none
        integer sizeOftype         ! for binary saving
        type(wlcsim_params), intent(in) :: wlc_p             ! to be save or filled
        type(wlcsim_data), intent(in) :: wlc_d             ! to be save or filled
        CHARACTER(LEN = 16), intent(in) :: baseName ! for example 'record/'
        CHARACTER(LEN = 16) fileName ! fileName
        CHARACTER(LEN = 16) suffix    ! end of file name
        LOGICAL exists    ! does file already exist?

        !  ------parameters -----

        sizeOftype = int(SIZEOF(wlc_p))
        suffix = 'parameters'
        fileName = trim(baseName) // trim(suffix)
        inquire(file = fileName,exist = exists)
        if(exists) then
            open(unit = outFileUnit,file = fileName, status = 'old', &
                form = 'unformatted',access = 'direct',recl = sizeOftype)
        else
            open(unit = outFileUnit,file = fileName, status = 'new', &
                form = 'unformatted',access = 'direct',recl = sizeOftype)
        endif
        write(outFileUnit,rec = 1) wlc_p
        close(outFileUnit)

        ! -------- R --------

        sizeOftype = int(SIZEOF(wlc_d%R))
        suffix = 'R'
        fileName = trim(baseName) // trim(suffix)
        inquire(file = fileName,exist = exists)
        if(exists) then
            open(unit = outFileUnit,file = fileName, status = 'old', &
                form = 'unformatted',access = 'direct',recl = sizeOftype)
        else
            open(unit = outFileUnit,file = fileName, status = 'new', &
                form = 'unformatted',access = 'direct',recl = sizeOftype)
        endif
        write(outFileUnit,rec = 1) wlc_d%R
        close(outFileUnit)

        ! -------- U --------

        sizeOftype = int(SIZEOF(wlc_d%U))
        suffix = 'U'
        fileName = trim(baseName) // trim(suffix)
        inquire(file = fileName,exist = exists)
        if(exists) then
            open(unit = outFileUnit,file = fileName, status = 'old', &
                form = 'unformatted',access = 'direct',recl = sizeOftype)
        else
            open(unit = outFileUnit,file = fileName, status = 'new', &
                form = 'unformatted',access = 'direct',recl = sizeOftype)
        endif
        write(outFileUnit,rec = 1) wlc_d%U
        close(outFileUnit)

        ! -------- AB --------

        sizeOftype = int(SIZEOF(wlc_d%AB))
        suffix = 'AB'
        fileName = trim(baseName) // trim(suffix)
        inquire(file = fileName,exist = exists)
        if(exists) then
            open(unit = outFileUnit,file = fileName, status = 'old', &
                form = 'unformatted',access = 'direct',recl = sizeOftype)
        else
            open(unit = outFileUnit,file = fileName, status = 'new', &
                form = 'unformatted',access = 'direct',recl = sizeOftype)
        endif
        write(outFileUnit,rec = 1) wlc_d%AB
        close(outFileUnit)

        ! -------- Vol --------

        sizeOftype = int(SIZEOF(wlc_d%Vol))
        suffix = 'Vol'
        fileName = trim(baseName) // trim(suffix)
        inquire(file = fileName,exist = exists)
        if(exists) then
            open(unit = outFileUnit,file = fileName, status = 'old', &
                form = 'unformatted',access = 'direct',recl = sizeOftype)
        else
            open(unit = outFileUnit,file = fileName, status = 'new', &
                form = 'unformatted',access = 'direct',recl = sizeOftype)
        endif
        write(outFileUnit,rec = 1) wlc_d%Vol
        close(outFileUnit)
    end subroutine

    subroutine wlcsim_params_readBinary(wlc_p, wlc_d, baseName)
    ! This function reads what wlcsim_params_writebinary writes and
    ! stores it to wlc_p and wlc_d.  Be sure to allocate wlc_d before
    ! calling this command.
        implicit none
        integer sizeOftype         ! for binary saving
        type(wlcsim_params) wlc_p             ! to be save or filled
        type(wlcsim_data) wlc_d             ! to be save or filled
        CHARACTER(LEN = 16) baseName ! for example 'record/'
        CHARACTER(LEN = 16) fileName ! fileName
        CHARACTER(LEN = 16) suffix    ! end of file name
        LOGICAL exists    ! does file already exist?

        !  ------parameters -----

        sizeOftype = int(SIZEOF(wlc_p))
        suffix = 'parameters'
        fileName = trim(baseName) // trim(suffix)
        inquire(file = fileName,exist = exists)
        if(exists) then
            open(unit = inFileUnit,file = fileName, status = 'old', &
                form = 'unformatted',access = 'direct',recl = sizeOftype)
        else
            print*, 'Error in wlcsim_params_readBinary. file ',fileName,'does not exist'
            stop 1
        endif
        read(inFileUnit,rec = 1) wlc_p
        close(inFileUnit)

        ! -------- R --------

        sizeOftype = int(SIZEOF(wlc_d%R))
        suffix = 'R'
        fileName = trim(baseName) // trim(suffix)
        inquire(file = fileName,exist = exists)
        if(exists) then
            open(unit = inFileUnit,file = fileName, status = 'old', &
                form = 'unformatted',access = 'direct',recl = sizeOftype)
        else
            print*, 'Error in wlcsim_params_readBinary. file ',fileName,'does not exist'
            stop 1
        endif
        read(inFileUnit,rec = 1) wlc_d%R
        close(inFileUnit)

        ! -------- U --------

        sizeOftype = int(SIZEOF(wlc_d%U))
        suffix = 'U'
        fileName = trim(baseName) // trim(suffix)
        inquire(file = fileName,exist = exists)
        if(exists) then
            open(unit = inFileUnit,file = fileName, status = 'old', &
                form = 'unformatted',access = 'direct',recl = sizeOftype)
        else
            print*, 'Error in wlcsim_params_readBinary. file ',fileName,'does not exist'
            stop 1
        endif
        read(inFileUnit,rec = 1) wlc_d%U
        close(inFileUnit)

        ! -------- AB --------

        sizeOftype = int(SIZEOF(wlc_d%AB))
        suffix = 'AB'
        fileName = trim(baseName) // trim(suffix)
        inquire(file = fileName,exist = exists)
        if(exists) then
            open(unit = inFileUnit,file = fileName, status = 'old', &
                form = 'unformatted',access = 'direct',recl = sizeOftype)
        else
            print*, 'Error in wlcsim_params_readBinary. file ',fileName,'does not exist'
            stop 1
        endif
        read(inFileUnit,rec = 1) wlc_d%AB
        close(inFileUnit)

        ! -------- Vol --------

        sizeOftype = int(SIZEOF(wlc_d%Vol))
        suffix = 'Vol'
        fileName = trim(baseName) // trim(suffix)
        inquire(file = fileName,exist = exists)
        if(exists) then
            open(unit = inFileUnit,file = fileName, status = 'old', &
                form = 'unformatted',access = 'direct',recl = sizeOftype)
        else
            print*, 'Error in wlcsim_params_readBinary. file ',fileName,'does not exist'
            stop 1
        endif
        read(inFileUnit,rec = 1) wlc_d%Vol
        close(inFileUnit)
    end subroutine

    subroutine save_simulation_state(save_ind, wlc_d, wlc_p, outfile_base, stat)
        implicit none
        type(wlcsim_data), intent(in) :: wlc_d
        type(wlcsim_params), intent(in) :: wlc_p
        integer, intent(in) :: save_ind
        character(MAX_LOG10_SAVES) :: fileind ! num2str(i)
        character(MAXFILENAMELEN) :: filename
        character(MAXFILENAMELEN) :: outfile_base
        character(len = *) :: stat
        integer :: ind, j

        write (fileind,num2strFormatString) save_ind

        !Save various energy contiributions to file
        filename = trim(adjustL(outfile_base)) // 'energies'
        call wlcsim_params_appendEnergyData(save_ind, wlc_p, wlc_d, filename)

        !part 2.5 - adaptations
        filename = trim(adjustL(outfile_base)) // 'adaptations'
        call wlcsim_params_appendAdaptData(save_ind, wlc_p, wlc_d, filename)

        if (wlc_p%savePhi) then
            write(filename,num2strFormatString) save_ind
            filename = trim(adjustL(outfile_base)) // 'phi' // trim(adjustL(filename))
            call wlcsim_params_savePHI(wlc_p,wlc_d,filename)
        endif

        if (wlc_p%saveR) then
            write(filename,num2strFormatString) save_ind
            filename = trim(adjustL(outfile_base)) // 'r' // trim(adjustL(filename))
            call wlcsim_params_saveR(wlc_p,wlc_d,filename,.false.,stat)
        endif

        if (wlc_p%saveU) then
            write(filename,num2strFormatString) save_ind
            filename = trim(adjustL(outfile_base)) // 'u' // trim(adjustL(filename))
            call wlcsim_params_saveU(wlc_p,wlc_d,filename,stat)
        endif

        if (wlc_p%collisionDetectionType /= 0) then
            filename = trim(adjustL(outfile_base)) // 'coltimes'
            open(unit = outFileUnit, file = filename, status = 'REPLACE')
            do ind = 1,wlc_p%nt
                write(outFileUnit,*) (wlc_d%coltimes(ind,j), j = 1,wlc_p%nt)
            enddo
            close(outFileUnit)
        endif
    end subroutine save_simulation_state

    subroutine setup_runtime_floats()
        inf = ieee_value(inf, ieee_positive_inf)
        nan = ieee_value(nan, ieee_quiet_nan)
    end subroutine

    !Get Lks for parallel tempering from file
    subroutine get_LKs_from_file(wlc_d)
    implicit none
    type(wlcsim_data), intent(inout) :: wlc_d
    integer nLKs !number of linking numbers
    integer IOstatus
    integer TempLk
    integer i
    nLKs = 0
    open (unit = 1, file = 'input/LKs')
    do
        read(unit = 1, fmt = *,iostat = IOstatus) TempLk
        if (IOstatus /= 0) exit
        nLKs = nLKs + 1
    end do
    close(unit = 1)

    wlc_d%nLKs = nLKs
    allocate(wlc_d%LKs(nLks))

    open(unit = 1, file = 'input/LKs')
    do i = 1, nLks
        read(unit = 1,fmt = *) wlc_d%Lks(i)
    enddo
    close(unit = 1)
    end subroutine get_LKs_from_file
end module params
