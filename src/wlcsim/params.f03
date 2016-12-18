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
    use inputparams, only: MAXPARAMLEN

    IMPLICIT NONE

    public

    !!!     hardcoded params. will need to change if certain parts of code change
    ! number of wlc_p move types
    integer, parameter :: nMoveTypes = 10

    !!!     arbitrary technical choices
    ! precision of simulations
    integer, parameter :: dp = real64 ! preferred over SELECTED_real_Kind(15,307)
                                      ! only available as of fortran 2008
    ! used for all character buffers holding filenames
    integer, parameter :: MAXFILENAMELEN = 500
    ! unique file "units" to use for each file
    integer, parameter :: inFileUnit = 51
    integer, parameter :: outFileUnit = 52
    ! number of digits in max allowed number of save points, also change below
    integer, parameter :: MAX_LOG10_SAVES = 4
    ! format string to use for num2str(save_ind)
    character(4), parameter :: num2strFormatString = '(I4)'

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
    ! subroutine get_input_from_file), in the docs (TODO), and often the default
    ! values will help with understanding how the variable is used.
    !
    ! many of these variables are used only in certain kinds of simulations
    type wlcsim_params
        character(MAXPARAMLEN) codeName ! which simulation code to run
    !   Simulation parameters
        integer simType           ! whether to use WLC, ssWLC, or Gaussian Chain
        integer nT                ! Total number of beads  NT=NP*N*G
        integer nB                ! Number of beads in a polymer NB=N*G
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
        real(dp) xiu    ! rotational drag
        real(dp) eps      ! number of kuhn lengths between beads
        real(dp) del      ! number of persistence lengths between beads
        real(dp) chi      ! Chi parameter value (solvent-polymer) (Flory-Huggins separation constant (how much A/B's hate each))
        real(dp) kap      ! Incompressibility parameter of the melt
        real(dp) fpt_dist ! radius triggering collisions to be recorded in "coltimes"
        real(dp) lhc    !TODO something to do with intrapolymer interaction strength
        real(dp) vhc    !TODO something to do with intrapolymer interaction strength, fill in defaults, etc

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
        integer NBIN     ! Number of bins
        integer NBINX(3) ! Number of bin on an edge
        real(dp) lbox(3)  ! Box length (approximate)
        real(dp) dbin      ! Discretization size (approximate)

    !   Monte Carlo Variables (for adaptation)
        integer movetypes
        real(dp) PDesire(nMoveTypes) ! desired hit rate
        real(dp) MAXWindoW(nMoveTypes)         ! Max Size of window for bead selection
        real(dp) MINWindoW(nMoveTypes)         ! Min Size of window for bead selection
        real(dp) MinAMP(nMoveTypes) ! minium amplitude
        real(dp) MaxAMP(nMoveTypes) ! maximum amplitude
        integer MOVEON(nMoveTypes)         ! Is the move active
        real(dp) winTarget(nMoveTypes) ! target for ratio of window to anmplitude
        integer NADAPT(nMoveTypes) ! Nunber of steps between adapt
        real(dp) min_accept  ! threshold for deciding to usually not use a move
        integer reduce_move  ! whether or not to stop usuing a move when it goes below min_accept success
        integer winType      ! distributionof segment size in crankshaft move (unif=0, exp=1)

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
        integer fptColType      ! save first passage time vectors to file
        logical exitWhenCollided  ! stop sim with coltimes is full
        logical saveR             ! save R vectors to file
        logical saveU             ! save U vectors to file
        logical savePhi           ! save Phi vectors to file
        integer solType           ! Melt vs. Solution, Choose hamiltonian
        logical recenter_on       ! recenter in "quasi"-periodic boundary, should be off in BD
        logical useSchedule       ! use scheduled change in interactions strength(s)
        real(dp) KAP_ON     ! fraction of KAP energy contributing to "calculated" energy
        real(dp) CHI_ON     ! fraction of CHI energy contributing to "calculated" energy
        real(dp) Couple_ON  ! fraction of Coupling energy contributing to "calculated" energy
        logical restart     ! whether we are restarting from a previous sim or not
        logical INTERP_BEAD_LENNARD_JONES ! whether to have inter bead lennard jones energies
        logical field_int_on ! include field interactions (e.g. A/B interactions)
        logical bind_On

    !   parallel Tempering parameters
        character(MAXFILENAMELEN) repSuffix    ! prefix for writing files
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
        real(dp) INITIAL_MAX_S

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
        DOUBLE PRECISION, ALLOCATABLE :: CROSS(:,:)   !Matrix of information for crossings in a 2-D projection of the polymer
        DOUBLE PRECISION, ALLOCATABLE :: CROSSP(:,:)  !Matrix of crossings for the trial configuration
        INTEGER NCROSS
        INTEGER NCROSSP
        INTEGER CrossSize


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
    !    real(dp) DEINT    ! Change in self energy
        real(dp) DECouple ! Coupling energy
        real(dp) DEChi    ! chi interaction energy
        real(dp) DEKap    ! compression energy
        real(dp) Debind   ! Change in binding energy
        real(dp) DEField  ! Change in field energy
        real(dp) DESelf   ! change in self interaction energy
        real(dp) ECon     ! Confinement Energy
        integer NPHI  ! NUMBER o phi values that change, i.e. number of bins that were affected

    !   Parallel tempering variables
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
        ! WARNING: changing this to intent(out) means that unassigned values
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
        wlc_p%FRMfile=.FALSE.      ! don't load initial bead positions from file
        wlc_p%FRMCHEM=.FALSE.      ! don't load initial "chem" status from file
        wlc_p%FRMFIELD=.FALSE.     ! don't load initial field values from file
        wlc_p%saveR=.TRUE.         ! do save orientation vectors (makes restart of ssWLC possible)
        wlc_p%saveU=.TRUE.         ! do save orientation vectors (makes restart of ssWLC possible)
        wlc_p%fptColType=0         ! don't track first passage time collisions between beads
        wlc_p%savePhi=.FALSE.      ! don't save A/B density per bin (not needed for restart)
        wlc_p%FRwlc_pHEM=.FALSE.      ! don't load initial a/b states from file
        wlc_p%restart=.FALSE.      ! don't restart from previously saved simulation

        ! geometry options
        wlc_p%NP  =1               ! one polymer
        wlc_p%nB  =200             ! 200 beads per polymer
        wlc_p%nBpM = 10
        wlc_p%lp = 1                ! units of lp by default
        wlc_p%lt = 1                ! twist persistence length equals persistence length by default
        wlc_p%nMpP = wlc_p%nB / wlc_p%nBpM
        wlc_p%nT = wlc_p%nP * wlc_p%nB
        wlc_p%lbox(1)=nan     ! box size
        wlc_p%lbox(2)=nan
        wlc_p%lbox(3)=nan
        wlc_p%dbin =1.0_dp         ! unit bin size
        ! wlc_p%l0  =1.25_dp         ! TODO: not input
        wlc_p%beadVolume  = 0.1_dp ! much smaller than space between beads
        wlc_p%fA  =0.5_dp  ! half A, half B by default
        wlc_p%LAM =0.0_dp  ! perfectly random sequence  (see generating_sequences.rst for details)
        wlc_p%F_METH=0.5_dp ! half beads methylated by default
        wlc_p%LAM_METH=0.9_dp ! highly alternating sequence by default
        wlc_p%fPoly=0.025_dp   ! volume fraction of plymer corresponding to HELA DNA in cytoplasm
        wlc_p%k_field=0.0_dp ! some previous values: !1.5708_dp !0.3145_dp

        ! energy parameters
        wlc_p%EPS =0.3_dp ! TODO: not input
        wlc_p%CHI =0.0_dp ! don't use chi by default
        wlc_p%hA =0.0_dp  ! don't use weird artificial field by default
        wlc_p%KAP =10.0_dp ! "fairly" incompressible --Quinn
        wlc_p%EU  =0.0_dp ! a function of coarse graining. This should be set by hand if needed.
        wlc_p%EM  =0.0_dp ! by default, no hp1 binding energy included
        wlc_p%mu  =0.0_dp ! by default, no hp1 binding included
        wlc_p%HP1_Bind=0.0_dp ! by default, no binding of HP1 to each other

        ! options
        wlc_p%codeName= "" ! not bruno, brad, or quinn, so will error unless specified elsewehre
        wlc_p%movetypes=nMoveTypes
        wlc_p%initCondType = 4 ! 4 for sphereical
        wlc_p%confineType = 3 ! 3 for spherical
        wlc_p%solType=0    ! you better at least know whether you want a melt or solution
        wlc_p%ring=.false.    ! not a ring by default
        wlc_p%twist=.false.    ! don't include twist by default
        wlc_p%lk=0    ! no linking number (lays flat) by default
        wlc_p%min_accept=0.05 ! if a move succeeds < 5% of the time, start using it only every reduce_move cycles
        wlc_p%exitWhenCollided = .FALSE. ! stop sim when coltimes is full
        wlc_p%field_int_on = .FALSE. ! no field interactions by default
        wlc_p%bind_On = .FALSE. ! no binding energy by default
        wlc_p%INTERP_BEAD_LENNARD_JONES = .FALSE. ! no intrapolymer interactions by default

        ! timing options
        wlc_p%dt  = 1              ! set time scale to unit
        wlc_p%nInitMCSteps=4000  ! number of initilizing mc steps. 1000s x num polymers is good
        wlc_p%stepsPerSave=2000  ! number of simulation steps to take
        wlc_p%numSavePoints=200    ! 200 total save points, i.e. 2000 steps per save point
        wlc_p%NNoInt=100    ! number of simulation steps before turning on interactions in Quinn's wlc_p scheduler
        wlc_p%reduce_move=10 ! use moves that fall below the min_accept threshold only once every 10 times they would otherwise be used
        wlc_p%winType = 1   ! exponential fragment sizes mix better
        wlc_p%useSchedule=.False. ! use Quinn's scheduler to modify wlc_p params halfway through the simulation
        wlc_p%KAP_ON=1.0_dp ! use full value of compression energy
        wlc_p%CHI_ON=1.0_dp ! use full value of chi energy
        wlc_p%Couple_ON=1.0_dp ! use full value for coupling energy
        wlc_p%N_KAP_ON=1 ! turn on compression energy immediately
        wlc_p%N_CHI_ON=1 ! turn on chi energy immediately
        wlc_p%recenter_on=.TRUE. ! recenter the polymer in the box if it exists the boundary
        wlc_p%INITIAL_MAX_S=0.0_dp !TODO: for now must be set explicitly, was 0.1, Quinn, what is this value?

        ! replica options
        wlc_p%PTON=.FALSE.  ! use parallel if applicable
        wlc_p%stepsPerExchange=100      ! 100 steps between parallel tempering is pretty frequent
        wlc_p%nReplicaExchangePerSavePoint=1000      ! make this large
        wlc_p%NRepAdapt=1000  ! 1000 exchange attempts between adaptations
        wlc_p%lowerRepExe=0.09 ! TODO: enter justification for these defaults, if any.
        wlc_p%upperRepExe=0.18 ! TODO: fine if the only justification is "these just work"
        wlc_p%lowerCofRail=0.005
        wlc_p%upperCofRail=0.1
        wlc_p%indStartRepAdapt=10
        wlc_p%indendRepAdapt=20
        wlc_p%repAnnealSpeed=0.01
        wlc_p%replicaBounds=.TRUE.
        wlc_p%PT_twist =.False. ! don't parallel temper chi by default
        wlc_p%PT_chi =.False. ! don't parallel temper chi by default
        wlc_p%PT_h =.False. ! don't parallel temper h by default
        wlc_p%PT_kap =.False. ! don't parallel temper kap by default
        wlc_p%PT_mu =.False. ! don't parallel temper mu by default
        wlc_p%PT_couple =.False. ! don't parallel temper HP1 binding by default
        wlc_p%repSuffix = '' !default no suffix


        !switches to turn on various types of moves
        wlc_p%MOVEON(1)=1  ! crank-shaft move
        wlc_p%MOVEON(2)=1  ! slide move
        wlc_p%MOVEON(3)=1  ! pivot move
        wlc_p%MOVEON(4)=1  ! rotate move
        wlc_p%MOVEON(5)=1  ! full chain rotation
        wlc_p%MOVEON(6)=1  ! full chain slide
        wlc_p%MOVEON(7)=1  ! Change in Binding state
        wlc_p%MOVEON(8)=0  ! Chain flip ! TODO not working
        wlc_p%MOVEON(9)=1  ! Chain exchange
        wlc_p%MOVEON(10)=1 ! Reptation

        ! Balance move amplitude and number of beads
        do mctype=1,wlc_p%movetypes
            wlc_p%winTarget(mctype)=8.0_dp
            wlc_p%MINWindoW(mctype)=nan
            wlc_d%MCAMP(mctype)=nan 
        enddo
        do mctype=1,wlc_p%movetypes
            wlc_p%NADAPT(mctype)=1000 ! adapt after at most 1000 steps
            wlc_p%PDESIRE(mctype)=0.5_dp ! Target
            wlc_d%SUCCESS(mctype)=0
            wlc_d%SUCCESStotal(mctype)=0
            wlc_d%PHIT(mctype)=0.0_dp
        enddo

    end subroutine

    subroutine read_input_file(infile, wlc_p)
        use INPUTparaMS, only : readLINE, readA, readF, readI, reado
        implicit none
        type(wlcsim_params), intent(inout) :: wlc_p
        character(MAXFILENAMELEN), intent(in) :: infile
        character(MAXPARAMLEN) :: WORD ! parameter name currently being read in
        logical fileend ! have we reached end of file?
        integer nitems  ! number of items read from line
        integer pf      ! file unit for input file
        integer codeNameNumber

        pf = inFileUnit
        open(unit=PF,file=infile,status='OLD')

        ! read in the keywords one line at a time
        do
        CALL READLINE(PF, fileend, NITEMS)
        if (fileend .and. nitems == 0) exit

        ! skip empty lines
        if (NITEMS.EQ.0) cycle

        ! read in the keyword for this line
        CALL readA(WORD) ! you can call readA(WORD,CASESET=1)

        ! Skip any empty lines or any comment lines
        if (WORD(1:1).EQ.'#') cycle

        SELECT CASE(WORD) ! pick which keyword
        CASE('codeName') ! select version of wlcsim to run
            call readi(codeNameNumber)
            SELECT CASE(codeNameNumber)
            case(1)
                wlc_p%codeName='brad'
            case(2)
                wlc_p%codeName='bruno'
            case(3)
                wlc_p%codeName='quinn'
            case default
                print*, "1 for brad, 2 for bruno, 3 for quinn"
                print*, "unrecognized codeName"
                stop
            endselect
        CASE('initCondType')
            Call readI(wlc_p%initCondType)
            ! initCondType      |  Discription
            ! _____________|_________________________________
            !    1         |   straight line in y direction with random starting
            !    2         |   rerandomize when reaching boundary, slit in z dir
            !    3         |   rerandomize when reaching boundary, cube boundary
            !    4         |   rerandomize when reaching boundary, shpere
        CASE('confineType')
            Call readI(wlc_p%confinetype)
            ! confinetype  |  Discription
            ! _____________|_________________________________
            !    0         |  No confinement, periodic cube
            !    1         |  Between two plates in Z direction at 0 and lbox
            !    2         |  Cube of size lbox**3,  range: 0-lbox
            !    3         |  Circle of radius lbox, centered at lbox/2
            !    4         |  Periodic, unequal dimensions
        CASE('recenterOn')
            Call reado(wlc_p%recenter_on) ! recenter in periodic boundary
        CASE('solType')
            Call readI(wlc_p%solType)
            ! solType      | Discription
            !______________|_________________________________
            !    0         | Melt density fluctuates around fixed mean
            !    1         | Solution (For DNA)
        CASE('frmChem')
            Call reado(wlc_p%FRMCHEM) ! Initial chemical/methylation sequence from file
        CASE('frmFile')
            call reado(wlc_p%FRMfile) ! read configuration from file
        CASE('twist')
            CALL reado(wlc_p%twist) ! whether to include twist energies in wlc_p
        CASE('ring')
            CALL reado(wlc_p%ring) ! whether polymer is a ring or not
        CASE('interpBeadLennardJones')
            CALL reado(wlc_p%INTERP_BEAD_LENNARD_JONES) ! whether polymer is a ring or not
        CASE('fieldIntOn')
            CALL reado(wlc_p%field_int_on) !include field interactions
        CASE('bindOn')
            CALL reado(wlc_p%bind_on) ! Whether to include a binding state model
        CASE('lk')
            CALL readi(wlc_p%lk) ! linking number
        CASE('ptOn')
            CALL reado(wlc_p%PTON) ! parallel Tempering on
        CASE('saveR')
            Call reado(wlc_p%saveR)  ! save r vectors to file (every savepoint)
        CASE('fptColType')
            Call readi(wlc_p%fptColType)  
            ! fptColType   |  Description
            ! _____________|_________________________________
            !    0         |  No tracking fpt
            !    1         |  Use naive, O(n^2) algo for collision checking
            !    2         |  KDtree-based col checking (not implemented)
            !    3         |  custom, fast col checker written by bruno
        CASE('saveU')
            Call reado(wlc_p%saveU)  ! save u vectors to file (every savepoint)
        CASE('savePhi')
            Call reado(wlc_p%savePhi) ! save Phi vectors to file (every savepoint)
        CASE('exitWhenCollided')
            Call reado(wlc_p%exitWhenCollided) ! for finding looping time
        CASE('nb')
            Call readi(wlc_p%nb)  ! number of beads in the polymer
        CASE('dt')
            Call readF(wlc_p%dt)  ! time step of simulation. scaled non-dimensionalized time
        CASE('l')
            Call readF(wlc_p%l)  ! actual length in AU of polymer we want to simulate
        CASE('lt')
            Call readF(wlc_p%lt)  ! twist persistence length
        CASE('lp')
            Call readF(wlc_p%lp)  ! persistence length
        CASE('dBin')
            Call readF(wlc_p%dbin) ! spaitial descretation length, not tested
        CASE('lBox')
            Call readF(wlc_p%lbox(1)) ! side length of box
            wlc_p%lbox(2)=wlc_p%lbox(1)
            wlc_p%lbox(3)=wlc_p%lbox(1)
        CASE('lBoxX')
            Call readF(wlc_p%lbox(1)) ! side length of box in x direction
        CASE('lBoxY')
            Call readF(wlc_p%lbox(2)) ! side length of box in y direction
        CASE('lBoxZ')
            Call readF(wlc_p%lbox(3)) ! side length of box in z direction
        CASE('np')
            CALL readI(wlc_p%NP)  ! Number of polymers
        CASE('nBpM')
            Call readI(wlc_p%nBpM) ! Beads per monomer
        CASE('nMpP')
            CALL readI(wlc_p%nMpP) ! Number of monomers in a polymer
        CASE('nNoInt')
            Call readI(wlc_p%NNoInt) ! save points before turning on interaction
        CASE('nKapOn')
            call readI(wlc_p%N_KAP_ON) ! when to turn compression energy on
        CASE('nChiOn')
            call readI(wlc_p%N_CHI_ON) ! when to turn CHI energy on
        CASE('numSavePoints')
            Call readI(wlc_p%numSavePoints) ! total number of save points
        CASE('nInitMCSteps')
            Call readI(wlc_p%nInitMCSteps) ! num initial mc steps
        CASE('stepsPerSave')
            Call readI(wlc_p%stepsPerSave) ! steps per save point
        CASE('stepsPerExchange')
            Call readI(wlc_p%stepsPerExchange) ! number of steps between parallel tempering
        CASE('nReplicaExchangePerSavePoint')
            call readI(wlc_p%nReplicaExchangePerSavePoint) ! read the variable
        CASE('fPoly')
            Call readF(wlc_p%Fpoly) ! Fraction Polymer
        CASE('beadVolume')
            Call readF(wlc_p%beadVolume) ! Bead volume
        CASE('fA')
            Call readF(wlc_p%FA) ! Fraction of A beads (fraction bound)
        CASE('lam')
            Call readF(wlc_p%LAM) ! Chemical correlation parameter
        CASE('EPS')
            Call readF(wlc_p%EPS) ! Elasticity l0/(2lp)
            print*, 'EPS depricated'
            stop 
        CASE('chi')
            Call readF(wlc_p%CHI) ! CHI parameter (definition depends on  hamiltoniaon
        CASE('hA')
            Call readF(wlc_p%hA) ! strength of externally applied field
        CASE('kap')
            Call readF(wlc_p%KAP) !  Incompressibility parameter
        CASE('EU')
            Call readF(wlc_p%EU) ! Energy of binding for unmethalated
        CASE('EM')
            Call readF(wlc_p%EM) ! Energy of binding for methalated
        CASE('mu')
            Call readF(wlc_p%MU) ! chemical potential of HP1
        CASE('energyCouple')
            Call readF(wlc_p%HP1_Bind) ! Energy of binding of HP1 to eachother
        CASE('fMeth')
            Call readF(wlc_p%F_METH) ! Fraction methalated
        CASE('lamMeth')
            Call readF(wlc_p%LAM_METH) ! eigenvalue of methalation setup
        CASE('crankShaftOn')
            Call readI(wlc_p%MOVEON(1)) ! is Crank shaft move on 1/0
        CASE('slideOn')
            Call readI(wlc_p%MOVEON(2)) ! is Slide move on 1/0
        CASE('pivotOn')
            Call readI(wlc_p%MOVEON(3)) ! is Pivot move on 1/0
        CASE('rotateOn')
            Call readI(wlc_p%MOVEON(4)) ! is single bead rotate on 1/0
        CASE('fullChainRotationOn')
            Call readI(wlc_p%MOVEON(5)) ! is full chain rotate on 1/0
        CASE('fullChainSlideOn')
            Call readI(wlc_p%MOVEON(6)) ! is full chain slide on 1/0
        CASE('bindMoveOn')
            Call readI(wlc_p%MOVEON(7)) ! is bind/unbind move on 1/0
        CASE('chainFlipMoveOn')
            Call readI(wlc_p%MOVEON(8)) ! is flip move move on 1/0
        CASE('chainSwapMoveOn')
            Call readI(wlc_p%MOVEON(9)) ! is chain swap move on 1/0
        CASE('reptationMoveOn')
            Call readI(wlc_p%MOVEON(10)) ! is reptation move on 1/0
        CASE('minCrankShaftWin')
            Call readF(wlc_p%MINWindoW(1)) ! min mean window size
        CASE('minSlideWin')
            Call readF(wlc_p%MINWindoW(2))
        CASE('minPivotWin')
            Call readF(wlc_p%MINWindoW(3))
        CASE('minBindWin')
            Call readF(wlc_p%MINWindoW(7))
        CASE('reduceMove')
            Call readI(wlc_p%reduce_move) !  only exicute unlikely movetypes every ____ cycles
        CASE('winType')
            call readI(wlc_p%winType)   ! fragment size distribution for crankshaft move
        CASE('minAccept')
            Call readF(wlc_p%MIN_ACCEPT) ! below which moves are turned off
        CASE('crankShaftTarget')
            Call readF(wlc_p%winTarget(1)) ! target window size for crank shaft move
        CASE('slideTarget')
            Call readF(wlc_p%winTarget(2)) ! target window size for slide move
        CASE('pivotTarget')
            Call readF(wlc_p%winTarget(3)) ! target window size for Pivot move
        CASE('strengthSchedule')
            Call reado(wlc_p%useSchedule) ! use scheduled ramp in interaction strength(s)
        CASE('nRepAdapt')
            Call readI(wlc_p%NRepAdapt)  ! number of exchange attemts between adapt
        CASE('lowerRepExe')
            Call readF(wlc_p%lowerRepExe) ! when to decrease cof spacing
        CASE('upperRepExe')
            Call readF(wlc_p%upperRepExe) ! when to increase cof spacing
        CASE('lowerCofRail')
            Call readF(wlc_p%lowerCofRail) ! minumum acceptable Cof
        CASE('upperCofRail')
            Call readF(wlc_p%upperCofRail) ! maximum acceptable Cof
        CASE('indStartRepAdapt')
            Call readI(wlc_p%indStartRepAdapt) ! ind to start rep. cof. adaptiation on
        CASE('indEndRepAdapt')
            Call readI(wlc_p%indendRepAdapt) ! turn off rep adapt
        CASE('repAnnealSpeed')
            Call readF(wlc_p%repAnnealSpeed)  ! max change in cof. every adjust
        CASE('frmField')
            Call reado(wlc_p%FRMFIELD)  ! read field from file
        CASE('kField')
            Call readF(wlc_p%k_field)  ! wave mode for default field
        CASE('replicaBounds')
            Call reado(wlc_p%replicaBounds) ! insure that 0 < s < 1
        CASE('initialMaxS')
            call readF(wlc_p%INITIAL_MAX_S) ! inital s of rep with highest s
        CASE('parallelTempChi')
            call reado(wlc_p%PT_chi) ! parallel temper chi
        CASE('parallelTempH')
            call reado(wlc_p%PT_h) ! parallel temper h
        CASE('parallelTempKap')
            call reado(wlc_p%PT_kap) ! parallel temper kap
        CASE('parallelTempMu')
            call reado(wlc_p%PT_mu) ! parallel temper mu
        CASE('parallelTempCouple')
            call reado(wlc_p%PT_couple) ! parallel temper HP1_bind
        CASE('parallelTempTwist')
            call reado(wlc_p%pt_twist)  ! parallel temper over linking numbers
        CASE('restart')
            call reado(wlc_p%restart) ! Restart from parallel tempering
        CASE DEFAULT
            print*, "Warning, the input key ", trim(Word), " hase been changed to camel case"
           SELECT CASE(WORD) ! pick which keyword
            CASE('initCondType')
                Call readI(wlc_p%initCondType)
                ! initCondType      |  Discription
                ! _____________|_________________________________
                !    1         |   straight line in y direction with random starting
                !    2         |   rerandomize when reaching boundary, slit in z dir
                !    3         |   rerandomize when reaching boundary, cube boundary
                !    4         |   rerandomize when reaching boundary, shpere
            CASE('CONFINEtype')
                Call readI(wlc_p%confinetype)
                ! confinetype  |  Discription
                ! _____________|_________________________________
                !    0         |  No confinement, periodic cube
                !    1         |  Between two plates in Z direction at 0 and lbox
                !    2         |  Cube of size lbox**3,  range: 0-lbox
                !    3         |  Circle of radius lbox, centered at lbox/2
                !    4         |  Periodic, unequal dimensions
            CASE('RECENTER_ON')
                Call reado(wlc_p%recenter_on) ! recenter in periodic boundary
            CASE('solType')
                Call readI(wlc_p%solType)
                ! solType      | Discription
                !______________|_________________________________
                !    0         | Melt density fluctuates around fixed mean
                !    1         | Solution (For DNA)
            CASE('FRMCHEM')
                Call reado(wlc_p%FRMCHEM) ! Initial chemical/methylation sequence from file
            CASE('FRMfile')
                call reado(wlc_p%FRMfile) ! read configuration from file
            CASE('TWIST')
                CALL reado(wlc_p%twist) ! whether to include twist energies in wlc_p
            CASE('RING')
                CALL reado(wlc_p%ring) ! whether polymer is a ring or not
            CASE('INTERP_BEAD_LENNARD_JONES')
                CALL reado(wlc_p%INTERP_BEAD_LENNARD_JONES) ! whether polymer is a ring or not
            CASE('field_int_on')
                CALL reado(wlc_p%field_int_on) !include field interactions
            CASE('BIND_ON')
                CALL reado(wlc_p%bind_on) ! Whether to include a binding state model
            CASE('LK')
                CALL readi(wlc_p%lk) ! linking number
            CASE('PTON')
                CALL reado(wlc_p%PTON) ! parallel Tempering on
            CASE('SAVE_R')
                Call reado(wlc_p%saveR)  ! save u vectors to file (every savepoint)
            CASE('FPT_COL_TYPE')
                Call readi(wlc_p%fptColType)  ! save u vectors to file (every savepoint)
                ! fptColType   |  Description
                ! _____________|_________________________________
                !    0         |  No tracking fpt
                !    1         |  Use naive, O(n^2) algo for collision checking
                !    2         |  KDtree-based col checking (not implemented)
                !    3         |  custom, fast col checker written by bruno
            CASE('SAVE_U')
                Call reado(wlc_p%saveU)  ! save u vectors to file (every savepoint)
            CASE('SAVE_PHI')
                Call reado(wlc_p%savePhi) ! save Phi vectors to file (every savepoint)
            CASE('EXIT_WHEN_COLLIDED')
                Call reado(wlc_p%exitWhenCollided)  ! save u vectors to file (every savepoint)
            CASE('nb')
                Call readi(wlc_p%nb)  ! number of beads in the polymer
            CASE('dt')
                Call readF(wlc_p%dt)  ! time step of simulation. scaled non-dimensionalized time
            CASE('l')
                Call readF(wlc_p%l)  ! actual length in AU of polymer we want to simulate
            CASE('lt')
                Call readF(wlc_p%lt)  ! persistence length
            CASE('lp')
                Call readF(wlc_p%lp)  ! twist persistence length
            CASE('dbin')
                Call readF(wlc_p%dbin) ! spaitial descretation length, not tested
            CASE('lbox')
                Call readF(wlc_p%lbox(1)) ! side length of box
                wlc_p%lbox(2)=wlc_p%lbox(1)
                wlc_p%lbox(3)=wlc_p%lbox(1)
            CASE('lboxX')
                Call readF(wlc_p%lbox(1)) ! side length of box in x direction
            CASE('lboxY')
                Call readF(wlc_p%lbox(2)) ! side length of box in y direction
            CASE('lboxZ')
                Call readF(wlc_p%lbox(3)) ! side length of box in z direction
            CASE('NP')
                CALL readI(wlc_p%NP)  ! Number of polymers
            CASE('nBpM')
                Call readI(wlc_p%nBpM) ! Beads per monomer
            CASE('nMpP')
                CALL readI(wlc_p%nMpP) ! Number of monomers in a polymer
            CASE('NNOINT')
                Call readI(wlc_p%NNoInt) ! save points before turning on interaction
            CASE('N_KAP_ON')
                call readI(wlc_p%N_KAP_ON) ! when to turn compression energy on
            CASE('N_CHI_ON')
                call readI(wlc_p%N_CHI_ON) ! when to turn CHI energy on
            CASE('numSavePoints')
                Call readI(wlc_p%numSavePoints) ! total number of save points
            CASE('nInitMCSteps')
                Call readI(wlc_p%nInitMCSteps) ! num initial mc steps
            CASE('stepsPerSave')
                Call readI(wlc_p%stepsPerSave) ! steps per save point
            CASE('stepsPerExchange')
                Call readI(wlc_p%stepsPerExchange) ! number of steps between parallel tempering
            CASE('nReplicaExchangePerSavePoint')
                call readI(wlc_p%nReplicaExchangePerSavePoint) ! read the variable
            CASE('FPOLY')
                Call readF(wlc_p%Fpoly) ! Fraction Polymer
            CASE('V')
                Call readF(wlc_p%beadVolume) ! Bead volume
            CASE('FA')
                Call readF(wlc_p%FA) ! Fraction of A beads (fraction bound)
            CASE('LAM')
                Call readF(wlc_p%LAM) ! Chemical correlation parameter
            CASE('EPS')
                Call readF(wlc_p%EPS) ! Elasticity l0/(2lp)
            CASE('CHI')
                Call readF(wlc_p%CHI) ! CHI parameter (definition depends on  hamiltoniaon
            CASE('H_A')
                Call readF(wlc_p%hA) ! strength of externally applied field
            CASE('KAP')
                Call readF(wlc_p%KAP) !  Incompressibility parameter
            CASE('EU')
                Call readF(wlc_p%EU) ! Energy of binding for unmethalated
            CASE('EM')
                Call readF(wlc_p%EM) ! Energy of binding for methalated
            CASE('MU')
                Call readF(wlc_p%MU) ! chemical potential of HP1
            CASE('HP1_Bind')
                Call readF(wlc_p%HP1_Bind) ! Energy of binding of HP1 to eachother
            CASE('F_METH')
                Call readF(wlc_p%F_METH) ! Fraction methalated
            CASE('LAM_METH')
                Call readF(wlc_p%LAM_METH) ! eigenvalue of methalation setup
            CASE('CRANK_SHAFT_ON')
                Call readI(wlc_p%MOVEON(1)) ! is Crank shaft move on 1/0
            CASE('SLIDE_ON')
                Call readI(wlc_p%MOVEON(2)) ! is Slide move on 1/0
            CASE('PIVOT_ON')
                Call readI(wlc_p%MOVEON(3)) ! is Pivot move on 1/0
            CASE('ROTATE_ON')
                Call readI(wlc_p%MOVEON(4)) ! is single bead rotate on 1/0
            CASE('FULL_CHAIN_ROTATION_ON')
                Call readI(wlc_p%MOVEON(5)) ! is full chain rotate on 1/0
            CASE('FULL_CHAIN_SLIDE_ON')
                Call readI(wlc_p%MOVEON(6)) ! is full chain slide on 1/0
            CASE('Bind_MOVE_ON')
                Call readI(wlc_p%MOVEON(7)) ! is bind/unbind move on 1/0
            CASE('CHAIN_FLIP_MOVE_ON')
                Call readI(wlc_p%MOVEON(8)) ! is flip move move on 1/0
            CASE('CHAIN_SWAP_MOVE_ON')
                Call readI(wlc_p%MOVEON(9)) ! is chain swap move on 1/0
            CASE('REPTATION_MOVE_ON')
                Call readI(wlc_p%MOVEON(10)) ! is reptation move on 1/0
            CASE('MIN_CRANK_SHAFT_WIN')
                Call readF(wlc_p%MINWindoW(1)) ! min mean window size
            CASE('MIN_SLIDE_WIN')
                Call readF(wlc_p%MINWindoW(2))
            CASE('MIN_PIVOT_WIN')
                Call readF(wlc_p%MINWindoW(3))
            CASE('MIN_Bind_WIN')
                Call readF(wlc_p%MINWindoW(7))
            CASE('REDUCE_MOVE')
                Call readI(wlc_p%reduce_move) !  only exicute unlikely movetypes every ____ cycles
            CASE('WIN_TYPE')
                call readI(wlc_p%winType)   ! fragment size distribution for crankshaft move
            CASE('MIN_ACCEPT')
                Call readF(wlc_p%MIN_ACCEPT) ! below which moves are turned off
            CASE('CRANK_SHAFT_TARGET')
                Call readF(wlc_p%winTarget(1)) ! target window size for crank shaft move
            CASE('SLIDE_TARGET')
                Call readF(wlc_p%winTarget(2)) ! target window size for slide move
            CASE('PIVOT_TARGET')
                Call readF(wlc_p%winTarget(3)) ! target window size for Pivot move
            CASE('STRENGTH_SCHEDULE')
                Call reado(wlc_p%useSchedule) ! use scheduled ramp in interaction strength(s)
            CASE('N_REP_ADAPT')
                Call readI(wlc_p%NRepAdapt)  ! number of exchange attemts between adapt
            CASE('LOWER_REP_EXE')
                Call readF(wlc_p%lowerRepExe) ! when to decrease cof spacing
            CASE('UPPER_REP_EXE')
                Call readF(wlc_p%upperRepExe) ! when to increase cof spacing
            CASE('LOWER_COF_RAIL')
                Call readF(wlc_p%lowerCofRail) ! minumum acceptable Cof
            CASE('UPPER_COF_RAIL')
                Call readF(wlc_p%upperCofRail) ! maximum acceptable Cof
            CASE('ind_START_REP_ADAPT')
                Call readI(wlc_p%indStartRepAdapt) ! ind to start rep. cof. adaptiation on
            CASE('ind_end_REP_ADAPT')
                Call readI(wlc_p%indendRepAdapt) ! turn off rep adapt
            CASE('REP_ANNEAL_SPEED')
                Call readF(wlc_p%repAnnealSpeed)  ! max change in cof. every adjust
            CASE('FRMFIELD')
                Call reado(wlc_p%FRMFIELD)  ! read field from file
            CASE('K_FIELD')
                Call readF(wlc_p%k_field)  ! wave mode for default field
            CASE('REPLICA_BOUNDS')
                Call reado(wlc_p%replicaBounds) ! insure that 0 < s < 1
            CASE('INITIAL_MAX_S')
                call readF(wlc_p%INITIAL_MAX_S) ! inital s of rep with highest s
            CASE('PT_CHI')
                call reado(wlc_p%PT_chi) ! parallel temper chi
            CASE('PT_H')
                call reado(wlc_p%PT_h) ! parallel temper h
            CASE('PT_KAP')
                call reado(wlc_p%PT_kap) ! parallel temper kap
            CASE('PT_MU')
                call reado(wlc_p%PT_mu) ! parallel temper mu
            CASE('PT_COUPLE')
                call reado(wlc_p%PT_couple) ! parallel temper HP1_bind
            CASE('PT_TWIST')
                call reado(wlc_p%pt_twist)  ! parallel temper over linking numbers
            CASE('RESTART')
                call reado(wlc_p%restart) ! Restart from parallel tempering
            CASE DEFAULT
                print*, "Error in params%read_input_file.  Unidentified keyword:", &
                        TRIM(WORD)
                stop 1
            endSELECT

            !print*, "Error in params%read_input_file.  Unidentified keyword:", &
            !        TRIM(WORD)
            !stop 1
        endSELECT
        enddo
        close(PF)
    end subroutine read_input_file


    subroutine idiot_checks(wlc_p, wlc_d)
#if MPI_VERSION
        use mpi
#endif
        IMPLICIT NONE
        type(wlcsim_params), intent(out) :: wlc_p
        type(wlcsim_data), intent(out) :: wlc_d
        logical err
        integer numProcesses ! number of threads running
        integer (kind=4) mpi_err

        if ((wlc_p%NT.GT.200000).OR.(wlc_p%NT.lt.1)) then
            print*, "Requested ", wlc_p%NT," beads."
            stop 1
        endif
        if ((wlc_p%NBIN.GT.20000).or.(wlc_p%NBIN.lt.1)) then
            print*, "Requested ", wlc_p%NBIN," bins."
            stop 1
        endif

        if (wlc_p%lBox(1) .ne. wlc_p%lBox(1)) then
            print*, "No box size set.  If you need a box please specify it."
        else
            if (wlc_p%confineType.eq.3) then
                if (abs((wlc_p%fPoly*(1.0/6.0_dp)*PI*wlc_p%LBOX(1)**3 / &
                        (wlc_p%beadVolume*wlc_p%NT)) - 1)>0.02) then
                     print*, "Error: volume fraction incorrect"
                     stop
                 endif
            else
                if (abs(wlc_p%fPoly*wlc_p%LBOX(1)**3/(wlc_p%beadVolume*wlc_p%NT)-1)>0.02) then
                     print*, "Error: volume fraction incorrect"
                     stop
                endif
            endif
        endif

        call stop_if_err(wlc_p%fptColType == 2, &
            'KD-tree based collision detection not yet implemented.')

        call stop_if_err(wlc_p%REND > wlc_p%L, &
            "Requesting initial end-to-end distance larger than polymer length.")

        if ((wlc_p%NBINX(1)-wlc_p%NBINX(2).ne.0).or. &
            (wlc_p%NBINX(1)-wlc_p%NBINX(3).ne.0)) then
            err = wlc_p%solType.eq.1
            call stop_if_err(err, "Solution not tested with non-cube box, more coding needed")
            err = wlc_p%confinetype.ne.4
            call stop_if_err(err, "Unequal boundaries require confinetype=4")
            err = wlc_p%initCondType.eq.4
            call stop_if_err(err, "You shouldn't put a shpere in and unequal box!")
        endif

        err = wlc_p%NBINX(1)*wlc_p%NBINX(2)*wlc_p%NBINX(3).ne.wlc_p%NBIN
        call stop_if_err(err, "error in mcsim. Wrong number of bins")

        if (wlc_p%NT.ne.wlc_p%nMpP*wlc_p%NP*wlc_p%nBpM) then
            print*, "Number of beads doesn't add up.  You should have NT=nMpP*NP*nBpM"
            print*, "Insteead you have:  NT=",wlc_p%NT," nMpP=",wlc_p%nMpP," NP=",wlc_p%NP," nBpM=",wlc_p%nBpM
            stop 1
        endif

        if (wlc_p%NB.ne.wlc_p%nMpP*wlc_p%nBpM) then
            print*, "error in mcsim.  NB=",wlc_p%NB," nMpP=",wlc_p%nMpP," nBpM=",wlc_p%nBpM
            stop 1
        endif

        err = wlc_p%NNoInt.gt.wlc_p%indStartRepAdapt
        call stop_if_err(err, "error in mcsim. don't run adapt without int on")

        err = wlc_p%NNoInt.gt.wlc_p%N_CHI_ON
        call stop_if_err(err, "error in mcsim. Can't have chi without int on")

        err = wlc_p%NNoInt.gt.wlc_p%N_KAP_ON
        call stop_if_err(err, "error in mcsim. Can't have kap without int on")
#if MPI_VERSION
    if (wlc_p%pt_twist) then
        if (.NOT.wlc_p%twist) then
            print *, 'parallel tempering on twist, but twist off'
            stop
        endif
        call MPI_COMM_SIZE(MPI_COMM_WORLD,numProcesses,mpi_err)
        if (wlc_d%nLKs+1.ne.numProcesses) then
            print *, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            print *, 'number of threads not equal to number of replicas!'
            print *, 'exiting...'
            print *, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            stop
        endif
  endif
#endif
    end subroutine


    subroutine get_input_from_file(infile, wlc_p, wlc_d)
        ! Based on Elena's readkeys subroutine
        IMPLICIT NONE
        type(wlcsim_params), intent(out) :: wlc_p
        type(wlcsim_data), intent(out) :: wlc_d
        character(MAXFILENAMELEN), intent(in) :: infile  ! file with parameters

        call set_param_defaults(wlc_p,wlc_d)

        call read_input_file(infile, wlc_p)

        ! for derived parameters other than wlc ones
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
        implicit none
        type(wlcsim_data), intent(inout)   :: wlc_d
        type(wlcsim_params), intent(in)    :: wlc_p
        character(8) datedum  ! trash
        character(10) timedum ! trash
        character(5) zonedum  ! trash
        integer seedvalues(8) ! clock readings
        integer NT  ! total number of beads
        integer NBIN ! total number of bins
        integer i
        nt = wlc_p%nt
        nbin = wlc_p%nbin

        ALLOCATE(wlc_d%R(NT,3))
        ALLOCATE(wlc_d%U(NT,3))
        if (wlc_p%codeName /= 'bruno' .OR. wlc_p%nInitMCSteps /= 0) then
            Allocate(wlc_d%RP(NT,3))
            Allocate(wlc_d%UP(NT,3))
        endif
        if (wlc_p%codeName == 'quinn') then
            ALLOCATE(wlc_d%AB(NT))   !Chemical identity aka binding state
            ALLOCATE(wlc_d%ABP(NT))   !Chemical identity aka binding state
            ALLOCATE(wlc_d%PHIA(NBIN))
            ALLOCATE(wlc_d%PHIB(NBIN))
            ALLOCATE(wlc_d%DPHIA(NBIN))
            ALLOCATE(wlc_d%DPHIB(NBIN))
            Allocate(wlc_d%indPHI(NBIN))
            Allocate(wlc_d%PhiH(NBIN))
            ALLOCATE(wlc_d%Vol(NBIN))
            ALLOCATE(wlc_d%METH(NT)) !Underlying methalation profile
            do I=1,wlc_p%NBIN
                wlc_d%PHIA(I)=0.0_dp
                wlc_d%PHIB(I)=0.0_dp
            enddo
        endif
        !Allocate vector of writhe and elastic energies for replicas
        if (wlc_p%pt_twist) then
            allocate(wlc_d%Wrs(wlc_d%nLKs))
            allocate(wlc_d%eelasREPLICAS(wlc_d%nLKs,4))

        endif
        if (wlc_p%ring) then !TODO this should be if ("knot")
            wlc_d%NCross=0
            wlc_d%CrossSize=wlc_p%NB**2
            ALLOCATE(wlc_d%Cross(wlc_d%CrossSize,6))
            ALLOCATE(wlc_d%CrossP(wlc_d%CrossSize,6))
        endif
        !If parallel tempering is on, initialize the nodeNumbers
        !and initialize MPI
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

        if (wlc_p%fptColType /= 0) then
            allocate(wlc_d%coltimes(NT,NT))
            wlc_d%coltimes = -1.0_dp
        endif

        if (.false.) then ! if you wanted to set specific seed
            wlc_d%rand_seed=7171
        else ! seed from clock
            call date_and_time(datedum,timedum,zonedum,seedvalues)
            ! funny business
            wlc_d%rand_seed=int(-seedvalues(5)*1E7-seedvalues(6)*1E5 &
                      -seedvalues(7)*1E3-seedvalues(8))
            wlc_d%rand_seed=mod(wlc_d%rand_seed,10000)
            ! print*, "Random Intiger seed:",wlc_d%rand_seed
        endif

        call random_setseed(wlc_d%rand_seed, wlc_d%rand_stat)

        call initcond(wlc_d%R, wlc_d%U, wlc_d%AB, wlc_p%NT, wlc_p%NB, &
            wlc_p%NP, wlc_p%frmfile, pack_as_para(wlc_p), wlc_p%lbox, &
            wlc_p%initCondType, wlc_d%rand_stat)


        ! initialize energies
        call CalculateEnergiesFromScratch(wlc_p,wlc_d)
        wlc_d%EElas   =wlc_d%dEElas
        if (wlc_p%field_int_on) then
            wlc_d%ECouple =wlc_d%dECouple
            wlc_d%EKap    =wlc_d%dEKap
            wlc_d%ECHI    =wlc_d%dECHI
            wlc_d%EField  =wlc_d%dEField
            wlc_d%x_Field =wlc_d%dx_Field
            wlc_d%x_couple=wlc_d%dx_couple
            wlc_d%x_Kap   =wlc_d%dx_Kap
            wlc_d%x_Chi   =wlc_d%dx_Chi
        else
            wlc_d%ECouple =0.0_dp
            wlc_d%EKap    =0.0_dp
            wlc_d%ECHI    =0.0_dp
            wlc_d%EField  =0.0_dp
            wlc_d%x_Field =0.0_dp
            wlc_d%x_couple=0.0_dp
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
        PARA(1) = wlc_p%EB
        PARA(2) = wlc_p%EPAR
        PARA(3) = wlc_p%EPERP
        PARA(4) = wlc_p%GAM
        PARA(5) = wlc_p%ETA
        PARA(6) = wlc_p%XIR
        PARA(7) = wlc_p%XIU
        PARA(8) = wlc_p%LBOX(1)
        PARA(9) = wlc_p%lhc
        PARA(10) = wlc_p%VHC
    end function pack_as_para


    subroutine printDescription(wlc_p)
        IMPLICIT NONE
        type(wlcsim_params), intent(in) :: wlc_p
        print*, "---------------System Description---------------"
        print*, "Bead variables:"
        print*, " Total number of beads, NT=", wlc_p%NT
        print*, " Number of beads in a polymer, NB=", wlc_p%NB
        print*, " Number of monomers in a polymer, N=", wlc_p%nMpP
        print*, " Number of polymers, NP=",wlc_p%NP
        print*, " Number of beads in a monomer, nBpM=", wlc_p%nBpM
        print*, " fraction Methalated", wlc_p%F_METH
        print*, " LAM_METH", wlc_p%LAM_METH
        print*, "Length and volume Variables:"
        print*, " persistance length =",(wlc_p%L0/(2.0_dp*wlc_p%EPS))
        print*, " lbox=", wlc_p%lbox(1), wlc_p%lbox(2), wlc_p%lbox(3)
        print*, " Number of bins in x direction", &
                wlc_p%NBINX(1), wlc_p%NBINX(2),wlc_p%NBINX(3)
        print*, " Number of bins", wlc_p%NBIN
        print*, " spatial descritation dbin=",wlc_p%dbin
        print*, " L0=", wlc_p%L0
        print*, " volume fraction polymer =", wlc_p%Fpoly
        print*, " bead volume V=", wlc_p%beadVolume
        print*, "Energy Variables"
        print*, " elasticity EPS =", wlc_p%EPS
        print*, " solvent-polymer CHI =",wlc_p%CHI
        print*, " compression cof, KAP =", wlc_p%KAP
        print*, " field strength, hA =", wlc_p%hA

        print*, " -energy of binding unmethalated ", wlc_p%EU," more positive for favorable binding"
        print*, " -energy of binding methalated",wlc_p%EM
        print*, " HP1_Binding energy parameter", wlc_p%HP1_Bind
        print*, " chemical potential of HP1", wlc_p%mu
        print*, "Other:"
        print*, " confinetype:",wlc_p%confinetype
        print*, " initCondType:",wlc_p%initCondType
        print*, "---------------------------------------------"

    end subroutine

    subroutine tweak_param_defaults(wlc_p, wlc_d)
        IMPLICIT NONE
        type(wlcsim_params), intent(inout) :: wlc_p
        type(wlcsim_data), intent(inout) :: wlc_d


        wlc_p%L0 = wlc_p%l/wlc_p%nb
        !  Edit the following to optimize wlc_p performance
        !  Monte-Carlo simulation parameters
        if (wlc_d%MCAMP(1) .ne. wlc_d%MCAMP(1)) wlc_d%MCAMP(1)=0.5_dp*PI
        if (wlc_d%MCAMP(2) .ne. wlc_d%MCAMP(2)) wlc_d%MCAMP(2)=0.3_dp*wlc_p%L0
        if (wlc_d%MCAMP(3) .ne. wlc_d%MCAMP(3)) wlc_d%MCAMP(3)=0.5_dp*PI
        if (wlc_d%MCAMP(4) .ne. wlc_d%MCAMP(4)) wlc_d%MCAMP(4)=0.5_dp*PI
        if (wlc_d%MCAMP(5) .ne. wlc_d%MCAMP(5)) wlc_d%MCAMP(5)=0.5_dp*PI
        if (wlc_d%MCAMP(6) .ne. wlc_d%MCAMP(6)) wlc_d%MCAMP(6)=5.0_dp*wlc_p%L0

        !     Initial segment window for wlc_p moves
        wlc_d%Window(1)=15.0_dp ! used to be N*G
        wlc_d%Window(2)=15.0_dp ! used to be N*G
        wlc_d%Window(3)=15.0_dp ! used to be N*G
        wlc_d%Window(4)=1.0_dp
        wlc_d%Window(5)=dble(wlc_p%nMpP*wlc_p%nBpM)
        wlc_d%Window(6)=dble(wlc_p%nMpP*wlc_p%nBpM)
        wlc_d%Window(7)=15.0_dp ! used to be N*G
        wlc_d%Window(8)=dble(wlc_p%nMpP*wlc_p%nBpM)
        wlc_d%Window(9)=dble(wlc_p%nMpP*wlc_p%nBpM)
        wlc_d%Window(9)=1.0_dp

        !    Maximum window size (large windows are expensive)
        wlc_p%MAXWindoW(1)=dble(min(150,wlc_p%NB))
        wlc_p%MAXWindoW(2)=dble(min(150,wlc_p%NB))
        wlc_p%MAXWindoW(3)=dble(min(150,wlc_p%NB))
        wlc_p%MAXWindoW(4)=nan
        wlc_p%MAXWindoW(5)=nan
        wlc_p%MAXWindoW(6)=nan
        wlc_p%MAXWindoW(7)=dble(min(4,wlc_p%NB))
        wlc_p%MAXWindoW(8)=nan
        wlc_p%MAXWindoW(9)=nan
        wlc_p%MAXWindoW(9)=nan ! need to chaige code to allow >1

        
        if (wlc_p%MINWindow(1).ne.wlc_p%MINWindow(1)) wlc_p%MINWindoW(1)=dble(min(10,wlc_p%NB))
        if (wlc_p%MINWindow(2).ne.wlc_p%MINWindow(2)) wlc_p%MINWindoW(2)=dble(min(10,wlc_p%NB))
        if (wlc_p%MINWindow(3).ne.wlc_p%MINWindow(3)) wlc_p%MINWindoW(3)=dble(min(10,wlc_p%NB))
        if (wlc_p%MINWindow(7).ne.wlc_p%MINWindow(7)) wlc_p%MINWindoW(7)=dble(min(10,wlc_p%NB))

        wlc_p%MINAMP(1)=0.1_dp*PI
        wlc_p%MINAMP(2)=0.2_dp*wlc_p%L0
        wlc_p%MINAMP(3)=0.2_dp*PI
        wlc_p%MINAMP(4)=0.2_dp*PI
        wlc_p%MINAMP(5)=0.05_dp*PI
        wlc_p%MINAMP(6)=0.2_dp*wlc_p%L0
        wlc_p%MINAMP(7)=nan
        wlc_p%MINAMP(8)=nan
        wlc_p%MINAMP(9)=nan
        wlc_p%MINAMP(10)=nan

        wlc_p%MAXAMP(1)=1.0_dp*PI
        wlc_p%MAXAMP(2)=1.0_dp*wlc_p%L0
        wlc_p%MAXAMP(3)=1.0_dp*PI
        wlc_p%MAXAMP(4)=1.0_dp*PI
        wlc_p%MAXAMP(5)=1.0_dp*PI
        wlc_p%MAXAMP(6)=0.1*wlc_p%lbox(1)
        wlc_p%MAXAMP(7)=nan
        wlc_p%MAXAMP(8)=nan
        wlc_p%MAXAMP(9)=nan
        wlc_p%MAXAMP(10)=nan

        ! Solution 
        wlc_p%NBINX(1)=nint(wlc_p%LBOX(1)/wlc_p%dBin)
        wlc_p%NBINX(2)=nint(wlc_p%LBOX(2)/wlc_p%dBin)
        wlc_p%NBINX(3)=nint(wlc_p%LBOX(3)/wlc_p%dBin)
        wlc_p%NBIN=wlc_p%NBINX(1)*wlc_p%NBINX(2)*wlc_p%NBINX(3)
        if (abs(wlc_p%dBin*sum(wlc_p%NBINX) /sum(wlc_p%LBOX) -1)>0.001) then
            print*, "Warning: I have changed the descritation lenth from:", wlc_p%dBin
            wlc_p%dBin=wlc_p%lBox(1)*dble(wlc_p%NBINX(1))
            print*, "to ", wlc_p%dBin
        endif

        if (wlc_p%codeName == 'brad') then
            ! initialize windows to number of beads
            wlc_p%MAXWindoW = wlc_p%nB         ! Max Size of window for bead selection
            wlc_p% MINWindoW  = 1         ! Min Size of window for bead selection

            ! Window amplitudes
            wlc_p%MinAMP = 0.0_dp ! minium amplitude
            wlc_p%MinAMP(1) = 0.07_dp*pi
            wlc_p%MinAMP(2) = 0.01_dp*wlc_p%l/wlc_p%nB
            wlc_p%MaxAMP = 2.0_dp*pi
            wlc_p%MaxAMP(2) = wlc_p%lbox(1)
            wlc_p%MaxAMP(6) = wlc_p%lbox(1)
        endif

        ! If ring is on, turn off the pivot move
        if (wlc_p%ring) then
            wlc_p%moveON(3) = 0
        endif

    end subroutine

    subroutine wlcsim_params_recenter(wlc_p,wlc_d)
    !  Prevents drift in periodic BC
        IMPLICIT NONE
        type(wlcsim_params), intent(in) :: wlc_p
        type(wlcsim_data), intent(inout) :: wlc_d
        integer IB, I, J   ! Couners
        real(dp) R0(3)  ! Offset to move by
        IB=1
        do I=1,wlc_p%NP
        R0(1)=nint(wlc_d%R(IB,1)/wlc_p%lbox(1)-0.5_dp)*wlc_p%lbox(1)
        R0(2)=nint(wlc_d%R(IB,2)/wlc_p%lbox(2)-0.5_dp)*wlc_p%lbox(2)
        R0(3)=nint(wlc_d%R(IB,3)/wlc_p%lbox(3)-0.5_dp)*wlc_p%lbox(3)
        if (abs(R0(1)*R0(2)*R0(3)) .gt. 0.0001_dp) then
            do J=1,wlc_p%NB
                wlc_d%R(IB,1)=wlc_d%R(IB,1)-R0(1)
                wlc_d%R(IB,2)=wlc_d%R(IB,2)-R0(2)
                wlc_d%R(IB,3)=wlc_d%R(IB,3)-R0(3)
                IB=IB+1
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
        IMPLICIT NONE
        type(wlcsim_params), intent(in) :: wlc_p
        type(wlcsim_data), intent(in) :: wlc_d
        integer I
        real(dp) EKap, ECouple, EChi,VV, PHIPOly
        print*,"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        print*, " PHIA  | PHIB  | PPoly |  Vol  | EKap  | EChi  |ECouple|"
        do I=1,wlc_p%NBIN
            VV=wlc_d%Vol(I)
            if (VV.le.0.1_dp) cycle
            PHIPOLY=wlc_d%PHIA(I)+wlc_d%PHIB(I)
            EChi=VV*(wlc_p%CHI/wlc_p%beadVolume)*PHIPoly*(1.0_dp-PHIPoly)
            ECouple=VV*wlc_p%HP1_Bind*(wlc_d%PHIA(I))**2
            if(PHIPoly.GT.1.0_dp) then
            EKap=VV*(wlc_p%KAP/wlc_p%beadVolume)*(PHIPoly-1.0_dp)**2
            else
            cycle
            EKap=0.0_dp
            endif
            write(*,"(4f8.4,3f8.1)") wlc_d%PHIA(I), wlc_d%PHIB(I), &
                                wlc_d%PHIA(I)+wlc_d%PHIB(I),wlc_d%Vol(I),&
                                EKap,EChi,ECouple
        enddo
        print*,"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    end subroutine

    subroutine printWindowStats(wlc_p, wlc_d)
    ! For realtime feedback on adaptation
        IMPLICIT NONE
        type(wlcsim_params), intent(in) :: wlc_p
        type(wlcsim_data), intent(in) :: wlc_d
        integer I ! counter
        I=0
        print*, "Succes | MCAMP | WindoW| type "
        do I=1,wlc_p%movetypes
            if (wlc_p%MOVEON(i).eq.1) then
                write(*,"(f8.5,2f8.2,1I8)") wlc_d%phit(i), wlc_d%MCAMP(i),  wlc_d%WindoW(i), i
            endif
        enddo
        return
    end subroutine

    subroutine wlcsim_params_LoadField(wlc_p,wlc_d,fileName)
        IMPLICIT NONE
        type(wlcsim_params), intent(in) :: wlc_p
        type(wlcsim_data), intent(inout) :: wlc_d
        integer I
        character(MAXFILENAMELEN) fileName ! file name to load from
        open (unit = inFileUnit, file = fileName, status = 'OLD')
        do I=1,wlc_p%NBIN
            read(inFileUnit,*) wlc_d%PHIH(I)
        enddo
        return
    end subroutine

    subroutine wlcsim_params_MakeField(wlc_p,wlc_d)
        IMPLICIT NONE
        type(wlcsim_params), intent(in) :: wlc_p
        type(wlcsim_data), intent(inout) :: wlc_d
        integer indBIN  ! index of bin
        integer IX,IY,IZ ! bin corrdinates

        do IX=1,wlc_p%NBINX(1)
            do IY=1,wlc_p%NBINX(2)
                do IZ=1,wlc_p%NBINX(3)
                    indBIN=IX+&
                        (IY-1)*wlc_p%NBINX(1)+&
                        (IZ-1)*wlc_p%NBINX(1)*wlc_p%NBINX(2)
                    wlc_d%PHIH(indBIN)=dsin(wlc_p%k_field*wlc_p%dbin*dble(IX))
                enddo
            enddo
        enddo
        return
    end subroutine

    subroutine wlcsim_params_loadAB(wlc_p,wlc_d,fileName)
    ! Loads AB for file...has not been tested
        IMPLICIT NONE
        type(wlcsim_params), intent(in) :: wlc_p
        type(wlcsim_data), intent(inout) :: wlc_d
        character(MAXFILENAMELEN), intent(in) :: fileName ! file name to load from
        integer IB, I, J ! counters
        open (unit = inFileUnit, file = fileName, status = 'OLD')
        IB=1
        do I=1,wlc_p%NP
        do J=1,wlc_p%NB
            read(inFileUnit,"(I2)") wlc_d%AB(IB)
            IB=IB+1
            enddo
        enddo
        close(inFileUnit)
    end subroutine

    subroutine wlcsim_params_saveR(wlc_p,wlc_d,fileName,repeatingBC)
    ! Writes R and AB to file for analysis
    ! Rx  Ry  Rz AB
        IMPLICIT NONE
        logical, intent(in) :: repeatingBC  ! 1 for reapeating boundary conditions
        integer I,J,IB  ! counters
        type(wlcsim_params), intent(in) :: wlc_p
        type(wlcsim_data), intent(in) :: wlc_d
        character(MAXFILENAMELEN), intent(in) :: fileName
        character(MAXFILENAMELEN) fullName
        fullName=  trim(fileName) // trim(wlc_p%repSuffix)
        fullName = trim(fullName)
        open (unit = outFileUnit, file = fullName, status = 'NEW')
        IB=1
        if (repeatingBC) then
            do I=1,wlc_p%NP
                do J=1,wlc_p%NB
                        WRITE(outFileUnit,"(3f10.3,I2)") &
                            wlc_d%R(IB,1)-0.*nint(wlc_d%R(IB,1)/wlc_p%lbox(1)-0.5_dp)*wlc_p%lbox(1), &
                            wlc_d%R(IB,2)-0.*nint(wlc_d%R(IB,2)/wlc_p%lbox(2)-0.5_dp)*wlc_p%lbox(2), &
                            wlc_d%R(IB,3)-0.*nint(wlc_d%R(IB,3)/wlc_p%lbox(3)-0.5_dp)*wlc_p%lbox(3), &
                            wlc_d%AB(IB)
                    IB=IB+1
                enddo
            enddo
            print*, "Error in wlcsim_params_saveR"
            print*, "Are you sure you want reapeating BC"
            stop 1
        else
            do I=1,wlc_p%NP
                do J=1,wlc_p%NB
                    if (wlc_p%solType.eq.0) then
                        WRITE(outFileUnit,"(3f10.3,I2)") &
                            wlc_d%R(IB,1),wlc_d%R(IB,2),wlc_d%R(IB,3),wlc_d%AB(IB)
                    else
                        WRITE(outFileUnit,"(3f10.3,I2)") &
                            wlc_d%R(IB,1),wlc_d%R(IB,2),wlc_d%R(IB,3),wlc_d%AB(IB), wlc_d%METH(IB)
                    endif
                    IB=IB+1
                enddo
            enddo
        endif
        close(outFileUnit)
    end subroutine

    subroutine wlcsim_params_savePHI(wlc_p,wlc_d,fileName)
    ! Saves PHIA and PHIB to file for analysis
        IMPLICIT NONE
        integer I  ! counters
        type(wlcsim_params), intent(in) :: wlc_p
        type(wlcsim_data), intent(in) :: wlc_d
        character(MAXFILENAMELEN), intent(in) :: fileName
        character(MAXFILENAMELEN) fullName
        fullName=  trim(fileName) // trim(wlc_p%repSuffix)
        open (unit = outFileUnit, file = fullName, status = 'NEW')
        do I=1,wlc_p%NBIN
            WRITE(outFileUnit,"(2f7.2)") wlc_d%PHIA(I),wlc_d%PHIB(I)
        enddo
        close(outFileUnit)
    end subroutine

    subroutine wlcsim_params_saveU(wlc_p,wlc_d,fileName)
    ! Saves U to ASCII file for analisys
        IMPLICIT NONE
        integer I,J,IB  ! counters
        type(wlcsim_params), intent(in) :: wlc_p
        type(wlcsim_data), intent(in) :: wlc_d
        character(MAXFILENAMELEN), intent(in) :: fileName
        character(MAXFILENAMELEN) fullName
        fullName=  trim(fileName) // trim(wlc_p%repSuffix)
        open (unit = outFileUnit, file = fullName, status = 'NEW')
        IB=1
        do I=1,wlc_p%NP
            do J=1,wlc_p%NB
                WRITE(outFileUnit,"(3f8.3,2I2)") wlc_d%U(IB,1),wlc_d%U(IB,2),wlc_d%U(IB,3)
                IB=IB+1
            enddo
        enddo
        close(outFileUnit)
    end subroutine

    subroutine save_parameters(wlc_p,fileName)
    ! Write a number of parameters ASCII variables to file for reccords
        IMPLICIT NONE
        type(wlcsim_params), intent(in) :: wlc_p
        character(MAXFILENAMELEN), intent(in) :: fileName
        character(MAXFILENAMELEN) fullName
        fullName=  trim(fileName) // trim(wlc_p%repSuffix)
        open (unit =outFileUnit, file = fullName, status = 'NEW')
            WRITE(outFileUnit,"(I8)") wlc_p%NT ! 1 Number of beads in simulation
            WRITE(outFileUnit,"(I8)") wlc_p%nMpP  ! 2 Number of monomers in a polymer
            WRITE(outFileUnit,"(I8)") wlc_p%NB ! 3 Number of beads in a polymer
            WRITE(outFileUnit,"(I8)") wlc_p%NP ! 4 Number of polymers in simulation
            WRITE(outFileUnit,"(I8)") wlc_p%NT ! 5 Number of beads in simulation
            WRITE(outFileUnit,"(I8)") wlc_p%nBpM  ! 6 Number of beads per monomer

            WRITE(outFileUnit,"(f10.5)") wlc_p%L0    ! Equilibrium segment length
            WRITE(outFileUnit,"(f10.5)") wlc_p%CHI  ! 8  initail CHI parameter value
            WRITE(outFileUnit,"(f10.5)") wlc_p%Fpoly ! Fraction polymer
            WRITE(outFileUnit,"(f10.5)") wlc_p%lbox(1)  ! 10 Lenth of box
            WRITE(outFileUnit,"(f10.5)") wlc_p%EU    ! Energy unmethalated
            WRITE(outFileUnit,"(f10.5)") wlc_p%EM    ! 12 Energy methalated
            WRITE(outFileUnit,"(f10.5)") wlc_p%HP1_Bind ! Energy of HP1 binding
            WRITE(outFileUnit,"(f10.5)") (wlc_p%L0/wlc_p%EPS) ! 14 Khun lenth
            WRITE(outFileUnit,"(A)") "-999"  ! for historic reasons
            WRITE(outFileUnit,"(f10.5)") wlc_p%F_METH  ! methalation fraction
            WRITE(outFileUnit,"(f10.5)") wlc_p%LAM_METH  ! methalation lambda
        close(outFileUnit)
    end subroutine

    subroutine wlcsim_params_appendEnergyData(save_ind, wlc_p, wlc_d, fileName)
    ! print Energy data
        IMPLICIT NONE
        type(wlcsim_params), intent(in) :: wlc_p
        type(wlcsim_data), intent(in) :: wlc_d
        integer, intent(in) :: save_ind
        character(MAXFILENAMELEN), intent(in) :: fileName
        LOGICAL isfile
        character(MAXFILENAMELEN) fullName
        fullName=  trim(fileName) // trim(wlc_p%repSuffix)
        inquire(file = fullName, exist=isfile)
        if (isfile) then
            open (unit = outFileUnit, file = fullName, status ='OLD', POSITION="append")
        else
            open (unit = outFileUnit, file = fullName, status = 'new')
            WRITE(outFileUnit,*) "ind | id |",&
                       " ebend  | eparll | EShear | ECoupl | E Kap  | E Chi  |",&
                       " EField | ebind  |  x_Mu  | Couple |  Chi   |  mu    |",&
                       "  Kap   | Field  |"
        endif
        WRITE(outFileUnit,"(2I5, 9f9.1,5f9.4)") save_ind, wlc_d%id, &
            wlc_d%EELAS(1), wlc_d%EELAS(2), wlc_d%EELAS(3), wlc_d%ECouple, &
            wlc_d%EKap, wlc_d%ECHI, wlc_d%EField, wlc_d%ebind, wlc_d%x_Mu, &
            wlc_p%HP1_Bind*wlc_p%Couple_on, wlc_p%CHI*wlc_p%CHI_ON, wlc_p%mu, wlc_p%KAP*wlc_p%KAP_ON,&
            wlc_p%hA
        close(outFileUnit)
    end subroutine

    subroutine wlcsim_params_appendAdaptData(save_ind, wlc_p, wlc_d, fileName)
    ! Appends wlc_p move adaptation data to the file
        IMPLICIT NONE
        type(wlcsim_params), intent(in) :: wlc_p
        type(wlcsim_data), intent(in) :: wlc_d
        integer, intent(in) :: save_ind
        LOGICAL isfile
        character(MAXFILENAMELEN), intent(in) :: fileName
        character(MAXFILENAMELEN) fullName
        fullName=  trim(fileName) // trim(wlc_p%repSuffix)
        inquire(file = fullName, exist=isfile)
        if (isfile) then
            open (unit = outFileUnit, file = fullName, status ='OLD', POSITION="append")
        else
            open (unit = outFileUnit, file = fullName, status = 'new')
            WRITE(outFileUnit,*) "ind| id|",&
                       " WIN 1 | AMP 1 | SUC 1 | WIN 2 | AMP 2 | SUC 2 |",&
                       " WIN 3 | AMP 3 | SUC 3 | ON  4 | AMP 4 | SUC 4 |",&
                       " ON  5 | AMP 5 | SUC 5 | ON  6 | AMP 6 | SUC 6 |",&
                       " ON  7 | SUC 7 | ON  8 | SUC 8 |", &
                       " ON  9 | SUC 9 | ON 10 | SUC 10|"
        endif
        WRITE(outFileUnit,"(2I4,26f8.3)") save_ind,wlc_d%id,&
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
        IMPLICIT NONE
        integer sizeOftype         ! for binary saving
        type(wlcsim_params), intent(in) :: wlc_p             ! to be save or filled
        type(wlcsim_data), intent(in) :: wlc_d             ! to be save or filled
        CHARACTER(LEN=16), intent(in) :: baseName ! for example 'record/'
        CHARACTER(LEN=16) fileName ! fileName
        CHARACTER(LEN=16) suffix    ! end of file name
        LOGICAL exists    ! does file already exist?

        !  ------parameters -----

        sizeOftype=int(SIZEOF(wlc_p))
        suffix='parameters'
        fileName=trim(baseName) // trim(suffix)
        inquire(file=fileName,exist=exists)
        if(exists) then
            open(unit=outFileUnit,file=fileName, status='old', &
                form='unformatted',access='direct',recl=sizeOftype)
        else
            open(unit=outFileUnit,file=fileName, status='new', &
                form='unformatted',access='direct',recl=sizeOftype)
        endif
        write(outFileUnit,rec=1) wlc_p
        close(outFileUnit)

        ! -------- R --------

        sizeOftype=int(SIZEOF(wlc_d%R))
        suffix='R'
        fileName=trim(baseName) // trim(suffix)
        inquire(file=fileName,exist=exists)
        if(exists) then
            open(unit=outFileUnit,file=fileName, status='old', &
                form='unformatted',access='direct',recl=sizeOftype)
        else
            open(unit=outFileUnit,file=fileName, status='new', &
                form='unformatted',access='direct',recl=sizeOftype)
        endif
        write(outFileUnit,rec=1) wlc_d%R
        close(outFileUnit)

        ! -------- U --------

        sizeOftype=int(SIZEOF(wlc_d%U))
        suffix='U'
        fileName=trim(baseName) // trim(suffix)
        inquire(file=fileName,exist=exists)
        if(exists) then
            open(unit=outFileUnit,file=fileName, status='old', &
                form='unformatted',access='direct',recl=sizeOftype)
        else
            open(unit=outFileUnit,file=fileName, status='new', &
                form='unformatted',access='direct',recl=sizeOftype)
        endif
        write(outFileUnit,rec=1) wlc_d%U
        close(outFileUnit)

        ! -------- AB --------

        sizeOftype=int(SIZEOF(wlc_d%AB))
        suffix='AB'
        fileName=trim(baseName) // trim(suffix)
        inquire(file=fileName,exist=exists)
        if(exists) then
            open(unit=outFileUnit,file=fileName, status='old', &
                form='unformatted',access='direct',recl=sizeOftype)
        else
            open(unit=outFileUnit,file=fileName, status='new', &
                form='unformatted',access='direct',recl=sizeOftype)
        endif
        write(outFileUnit,rec=1) wlc_d%AB
        close(outFileUnit)

        ! -------- Vol --------

        sizeOftype=int(SIZEOF(wlc_d%Vol))
        suffix='Vol'
        fileName=trim(baseName) // trim(suffix)
        inquire(file=fileName,exist=exists)
        if(exists) then
            open(unit=outFileUnit,file=fileName, status='old', &
                form='unformatted',access='direct',recl=sizeOftype)
        else
            open(unit=outFileUnit,file=fileName, status='new', &
                form='unformatted',access='direct',recl=sizeOftype)
        endif
        write(outFileUnit,rec=1) wlc_d%Vol
        close(outFileUnit)
    end subroutine

    subroutine wlcsim_params_readBinary(wlc_p, wlc_d, baseName)
    ! This function reads what wlcsim_params_writebinary writes and
    ! stores it to wlc_p and wlc_d.  Be sure to allocate wlc_d before
    ! calling this command.
        IMPLICIT NONE
        integer sizeOftype         ! for binary saving
        type(wlcsim_params) wlc_p             ! to be save or filled
        type(wlcsim_data) wlc_d             ! to be save or filled
        CHARACTER(LEN=16) baseName ! for example 'record/'
        CHARACTER(LEN=16) fileName ! fileName
        CHARACTER(LEN=16) suffix    ! end of file name
        LOGICAL exists    ! does file already exist?

        !  ------parameters -----

        sizeOftype=int(SIZEOF(wlc_p))
        suffix='parameters'
        fileName=trim(baseName) // trim(suffix)
        inquire(file=fileName,exist=exists)
        if(exists) then
            open(unit=inFileUnit,file=fileName, status='old', &
                form='unformatted',access='direct',recl=sizeOftype)
        else
            print*, 'Error in wlcsim_params_readBinary. file ',fileName,'does not exist'
            stop 1
        endif
        read(inFileUnit,rec=1) wlc_p
        close(inFileUnit)

        ! -------- R --------

        sizeOftype=int(SIZEOF(wlc_d%R))
        suffix='R'
        fileName=trim(baseName) // trim(suffix)
        inquire(file=fileName,exist=exists)
        if(exists) then
            open(unit=inFileUnit,file=fileName, status='old', &
                form='unformatted',access='direct',recl=sizeOftype)
        else
            print*, 'Error in wlcsim_params_readBinary. file ',fileName,'does not exist'
            stop 1
        endif
        read(inFileUnit,rec=1) wlc_d%R
        close(inFileUnit)

        ! -------- U --------

        sizeOftype=int(SIZEOF(wlc_d%U))
        suffix='U'
        fileName=trim(baseName) // trim(suffix)
        inquire(file=fileName,exist=exists)
        if(exists) then
            open(unit=inFileUnit,file=fileName, status='old', &
                form='unformatted',access='direct',recl=sizeOftype)
        else
            print*, 'Error in wlcsim_params_readBinary. file ',fileName,'does not exist'
            stop 1
        endif
        read(inFileUnit,rec=1) wlc_d%U
        close(inFileUnit)

        ! -------- AB --------

        sizeOftype=int(SIZEOF(wlc_d%AB))
        suffix='AB'
        fileName=trim(baseName) // trim(suffix)
        inquire(file=fileName,exist=exists)
        if(exists) then
            open(unit=inFileUnit,file=fileName, status='old', &
                form='unformatted',access='direct',recl=sizeOftype)
        else
            print*, 'Error in wlcsim_params_readBinary. file ',fileName,'does not exist'
            stop 1
        endif
        read(inFileUnit,rec=1) wlc_d%AB
        close(inFileUnit)

        ! -------- Vol --------

        sizeOftype=int(SIZEOF(wlc_d%Vol))
        suffix='Vol'
        fileName=trim(baseName) // trim(suffix)
        inquire(file=fileName,exist=exists)
        if(exists) then
            open(unit=inFileUnit,file=fileName, status='old', &
                form='unformatted',access='direct',recl=sizeOftype)
        else
            print*, 'Error in wlcsim_params_readBinary. file ',fileName,'does not exist'
            stop 1
        endif
        read(inFileUnit,rec=1) wlc_d%Vol
        close(inFileUnit)
    end subroutine

    subroutine save_simulation_state(save_ind, wlc_d, wlc_p, outfile_base)
        implicit none
        type(wlcsim_data), intent(in) :: wlc_d
        type(wlcsim_params), intent(in) :: wlc_p
        integer, intent(in) :: save_ind
        character(MAX_LOG10_SAVES) :: fileind ! num2str(i)
        character(MAXFILENAMELEN) :: filename
        character(MAXFILENAMELEN) :: outfile_base
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
            call wlcsim_params_saveR(wlc_p,wlc_d,filename,.false.)
        endif

        if (wlc_p%saveU) then
            write(filename,num2strFormatString) save_ind
            filename = trim(adjustL(outfile_base)) // 'u' // trim(adjustL(filename))
            call wlcsim_params_saveU(wlc_p,wlc_d,filename)
        endif

        if (wlc_p%fptColType /= 0) then
            filename = trim(adjustL(outfile_base)) // 'coltimes'
            open(unit=outFileUnit, file=filename, status='REPLACE')
            do ind=1,wlc_p%nt
                write(outFileUnit,*) (wlc_d%coltimes(ind,j), j=1,wlc_p%nt)
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
    type(wlcsim_data), intent(inout) :: wlc_d
    integer nLKs !number of linking numbers
    integer IOstatus
    integer TempLk
    integer i
    nLKs = 0
    open (unit = 1, file = 'input/LKs')
    do
        read(unit = 1, fmt = *,iostat=IOstatus) TempLk
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
