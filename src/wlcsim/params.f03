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
    use inputparams, only: MAXPARAMLEN

    IMPLICIT NONE

    private

    !!!     hardcoded params. will need to change if certain parts of code change
    ! number of wlcsim_p move types
    integer, parameter :: nMoveTypes = 10

    !!!     arbitrary technical choices
    ! precision of simulations
    integer, parameter :: dp = real64 ! preferred over SELECTED_real_Kind(15,307)
                                      ! only available as of fortran 2008
    ! used for all character buffers holding filenames
    integer, parameter :: MAXFILENAMELEN = 500
    ! unique file "units" to use for each file
    integer, parameter :: inFileUnit = 51

    !!!     universal constants
    ! fully accurate, adaptive precision
    real(dp) :: pi = 4 * atan(1.0_dp)
    ! won't get optimized away by compiler, see e.g.
    ! https://software.intel.com/en-us/forums/intel-visual-fortran-compiler-for-windows/topic/294680
    real(dp) :: nan = transfer((/ Z'00000000', Z'7FF80000' /),1.0_dp)
    ! the following would be preferred, but generated compilation errors...
    !real(dp) :: nan = IEEE_VALUE(IEEE_QUIET_NAN)

    ! for all parameters that cannot change during individual simulations
    ! these are documented more thoroughly where they are read in (see the
    ! subroutine get_input_from_file), in the docs (TODO), and often the default
    ! values will help with understanding how the variable is used.
    !
    ! many of these variables are used only in certain kinds of simulations
    type wlcsim_params
    !   Simulation parameters
        integer nT                ! Total number of beads  NT=NP*N*G
        integer nB                ! Number of beads in a polymer NB=N*G
        integer nMpP              ! Number of monomers (NOT BEADS!) in a polymer
        integer nBpM              ! number beads per monomer
        integer nP                ! Number of polymers
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
        real(dp) gam    ! average equilibrium interbead spacing
        real(dp) eta    ! bend-shear coupling parameter
        real(dp) xir    ! drag per unit persistence length
        real(dp) xiu    ! rotational drag
        real(dp) eps      ! number of kuhn lengths between beads
        real(dp) del      ! number of persistence lengths between beads
        real(dp) chi      ! Chi parameter value (solvent-polymer) (Flory-Huggins separation constant (how much A/B's hate each))
        real(dp) kap      ! Incompressibility parameter of the melt

    !   for passing 1st order phase transition in (quinn/shifan's) random copolymer wlcsim_p sims
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

    !   Timing variables
        integer NPT                ! number of steps between parallel tempering
        integer indMAX             ! total number of save points
        integer NSTEP              ! steps per save point
        integer NNoInt             ! save points before turning on NNoInt
        integer N_KAP_ON           ! when to turn KAP energy on
        integer N_CHI_ON           ! when to turn CHI energy on

    !   Switches
        logical ring              ! whether the polymer is a ring
        logical twist             ! whether to include twist (wlcsim_p only for now)
        integer LK                ! Linking number
        integer confinetype       ! type of Boundary Conditions
        integer initCondType           ! initial condition type
        logical field_interactions ! field-based self interactions on
        logical intrapolymer_stick_crossing_enforced ! field-based self interactions on
        logical FRwlcsim_pHEM           ! read initial chemical sequence from file
        logical FRMchem           ! read initial chemical/methylation state from file
        logical FRMfile           ! read initial condition R from file
        logical FRMField          ! read initial field from file
        logical saveU             ! save U vectors to file
        logical savePhi           ! save Phi vectors to file
        integer simtype           ! Melt vs. Solution, Choose hamiltonian
        logical recenter_on       ! recenter in "quasi"-periodic boundary, should be off in BD
        logical useSchedule       ! use scheduled change in interactions strength(s)
        real(dp) KAP_ON     ! fraction of KAP energy contributing to "calculated" energy
        real(dp) CHI_ON     ! fraction of CHI energy contributing to "calculated" energy
        real(dp) Couple_ON  ! fraction of Coupling energy contributing to "calculated" energy
        logical restart     ! whether we are restarting from a previous sim or not

    !   parallel Tempering parameters
        character(16) repSuffix    ! prefix for writing files
        logical PTON    ! whether or not to parallel temper
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
        real(dp), allocatable, dimension(:,:):: R   ! Conformation of polymer chains
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

    !   Monte Carlo Variables (for adaptation)
        real(dp) MCAMP(nMoveTypes) ! Amplitude of random change
        real(dp) WindoW(nMoveTypes)         ! Size of window for bead selection
        integer SUCCESS(nMoveTypes)        ! Number of successes
        real(dp) PHit(nMoveTypes) ! hit rate

    !   Energys
        !real(dp) Eint     ! running Eint
        real(dp) EELAS(3) ! Elastic force
        real(dp) ECHI     ! CHI energy
        real(dp) EKAP     ! KAP energy
        real(dp) ECouple  ! Coupling
        real(dp) ebind    ! binding energy
        real(dp) EField   ! Field energy

    !   Congigate Energy variables (needed to avoid NaN when cof-> 0 in rep exchange)
        real(dp) x_Chi,   dx_Chi
        real(dp) x_Couple,dx_Couple
        real(dp) x_Kap,   dx_Kap
        real(dp) x_Field, dx_Field
        real(dp) x_Mu,    dx_Mu

    !   Move Variables
        real(dp) DEELAS(3)   ! Change in bending energy
    !    real(dp) DEINT    ! Change in self energy
        real(dp) DECouple ! Coupling energy
        real(dp) DEChi    ! chi interaction energy
        real(dp) DEKap    ! compression energy
        real(dp) Debind   ! Change in binding energy
        real(dp) DEField  ! Change in field energy
        real(dp) ECon     ! Confinement Energy
        integer NPHI  ! NUMBER o phi values that change, i.e. number of bins that were affected

    !   Parallel tempering variables
        integer rep  ! which replica am I
        integer id   ! which thread am I
        integer error  ! MPI error

    !   indices
        integer ind                ! current save point
    end type

    public :: nMoveTypes, dp, pi, wlcsim_params, wlcsim_data

contains

    subroutine set_param_defaults(wlc_p)
        implicit none
        type(wlcsim_params), intent(out) :: wlc_p
        ! file IO
        wlc_p%FRMfile=.FALSE.      ! don't load initial bead positions from file
        wlc_p%FRMCHEM=.FALSE.      ! don't load initial "chem" status from file
        wlc_p%FRMFIELD=.false.     ! don't load initial field values from file
        wlc_p%saveU=.true.         ! do save orientation vectors (makes restart of ssWLC possible)
        wlc_p%savePhi=.FALSE.      ! don't save A/B density per bin (not needed for restart)
        wlc_p%FRwlcsim_pHEM=.FALSE.      ! don't load initial a/b states from file
        wlc_p%restart=.FALSE.      ! don't restart from previously saved simulation

        ! geometry options
        wlc_p%NP  =1               ! one polymer
        wlc_p%nB  =200             ! 200 beads per polymer
        wlc_p%nBpM = 10
        wlc_p%lp = 1                ! units of lp by default
        wlc_p%lt = 1                ! twist persistence length equals persistence length by default
        wlc_p%nMpP = wlc_p%nB / wlc_p%nBpM
        wlc_p%nT = wlc_p%nP * wlc_p%nB
        wlc_p%lbox(1)=25.0_dp      ! arbitrary box size
        wlc_p%lbox(2)=25.0_dp
        wlc_p%lbox(3)=25.0_dp
        wlc_p%dbin =1.0_dp         ! unit bin size
        wlc_p%l0  =1.25_dp         ! TODO: not input
        wlc_p%beadVolume   =0.1_dp ! much smaller than space between beads
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
        wlc_p%movetypes=nMoveTypes
        wlc_p%initCondType = 4 ! 4 for sphereical
        wlc_p%confineType = 3 ! 3 for spherical
        wlc_p%simtype=0    ! you better at least choose what kind of simulation you want
        wlc_p%ring=.false.    ! not a ring by default
        wlc_p%twist=.false.    ! don't include twist by default
        wlc_p%lk=0    ! no linking number (lays flat) by default
        wlc_p%min_accept=0.05 ! if a move succeeds < 5% of the time, start using it only every reduce_move cycles

        ! timing options
        wlc_p%NStep=400000  ! number of simulation steps to take
        wlc_p%NNoInt=100    ! number of simulation steps before turning on interactions in Quinn's wlcsim_p scheduler
        wlc_p%indMAX=200    ! 2000 steps per save point
        wlc_p%reduce_move=10 ! use moves that fall below the min_accept threshold only once every 10 times they would otherwise be used
        wlc_p%useSchedule=.False. ! use Quinn's scheduler to modify wlcsim_p params halfway through the simulation
        wlc_p%KAP_ON=1.0_dp ! use full value of compression energy
        wlc_p%CHI_ON=1.0_dp ! use full value of chi energy
        wlc_p%Couple_ON=1.0_dp ! use full value for coupling energy
        wlc_p%N_KAP_ON=1 ! turn on compression energy immediately
        wlc_p%N_CHI_ON=1 ! turn on chi energy immediately
        wlc_p%recenter_on=.TRUE. ! recenter the polymer in the box if it exists the boundary
        wlc_p%INITIAL_MAX_S=0.0_dp !TODO: for now must be set explicitly, was 0.1, Quinn, what is this value?

        ! replica options
        wlc_p%PTON=.TRUE.  ! use parallel if applicable
        wlc_p%NPT=100      ! 100 steps between parallel tempering is pretty frequent
        wlc_p%NRepAdapt=1000  ! 1000 exchange attempts between adaptations
        wlc_p%lowerRepExe=0.09 ! TODO: enter justification for these defaults, if any.
        wlc_p%upperRepExe=0.18 ! TODO: fine if the only justification is "these just work"
        wlc_p%lowerCofRail=0.005
        wlc_p%upperCofRail=0.1
        wlc_p%indStartRepAdapt=10
        wlc_p%indendRepAdapt=20
        wlc_p%repAnnealSpeed=0.01
        wlc_p%replicaBounds=.TRUE.
        wlc_p%PT_chi =.False. ! don't parallel temper chi by default
        wlc_p%PT_h =.False. ! don't parallel temper h by default
        wlc_p%PT_kap =.False. ! don't parallel temper kap by default
        wlc_p%PT_mu =.False. ! don't parallel temper mu by default
        wlc_p%PT_couple =.False. ! don't parallel temper HP1 binding by default
    end subroutine

    subroutine read_from_file(infile, wlc_p)
        use INPUTparaMS, only : readLINE, readA, readF, readI, reado
        implicit none
        type(wlcsim_params), intent(inout) :: wlc_p
        character(MAXFILENAMELEN), intent(in) :: infile
        character(MAXPARAMLEN) :: WORD ! parameter name currently being read in
        logical fileend ! have we reached end of file?
        integer nitems  ! number of items read from line
        integer pf      ! file unit for input file

        pf = inFileUnit
        open(unit=PF,file=infile,status='OLD')

        ! read in the keywords one line at a time
        do
        CALL READLINE(PF, fileend, NITEMS)
        if (fileend .and. nitems == 0) exit

        ! skip empty lines
        if (NITEMS.EQ.0) cycle

        ! read in the keyword for this line
        CALL readA(WORD,CASESET=1)

        ! Skip any empty lines or any comment lines
        if (WORD(1:1).EQ.'#') cycle

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
        CASE('SIMtype')
            Call readI(wlc_p%simtype)
            ! simtype      | Discription
            !______________|_________________________________
            !    0         | Melt density fluctuates around fixed mean
            !    1         | Solution (For DNA)
        CASE('FRMCHEM')
            Call reado(wlc_p%FRMCHEM) ! Initial chemical/methylation sequence from file
        CASE('FRMfile')
            call reado(wlc_p%FRMfile) ! read configuration from file
        CASE('TWIST')
            CALL reado(wlc_p%twist) ! whether to include twist energies in wlcsim_p
        CASE('RING')
            CALL reado(wlc_p%ring) ! whether polymer is a ring or not
        CASE('LK')
            CALL readi(wlc_p%lk) ! linking number
        CASE('PTON')
            CALL reado(wlc_p%PTON) ! parallel Tempering on
        CASE('SAVE_U')
            Call reado(wlc_p%saveU)  ! save u vectors to file (every savepoint)
        CASE('SAVE_PHI')
            Call reado(wlc_p%savePhi) ! save Phi vectors to file (every savepoint)
        CASE('nb')
            Call readi(wlc_p%nb)  ! actual length in AU of polymer we want to simulate
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
        CASE('indMAX')
            Call readI(wlc_p%indMAX) ! total number of save points
        CASE('NSTEP')
            Call readI(wlc_p%NStep) ! steps per save point
        CASE('NPT')
            Call readI(wlc_p%NPT) ! number of steps between parallel tempering
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
        CASE('RESTART')
            call reado(wlc_p%restart) ! Restart from parallel tempering
        CASE DEFAULT
            print*, "Error in params%read_from_file.  Unidentified keyword:", &
                    TRIM(WORD)
            stop 1
        endSELECT
        enddo
        close(PF)
    end subroutine read_from_file

    subroutine idiot_checks(wlc_p)
        IMPLICIT NONE
        type(wlcsim_params), intent(out) :: wlc_p

        if ((wlc_p%NBINX(1)-wlc_p%NBINX(2).ne.0).or. &
            (wlc_p%NBINX(1)-wlc_p%NBINX(3).ne.0)) then
            if (wlc_p%simtype.eq.1) then
                print*, "Solution not tested with non-cube box, more coding needed"
                stop 1
            endif
            if (wlc_p%confinetype.ne.4) then
                print*, "Unequal boundaries require confinetype=4"
                stop 1
            endif
            if (wlc_p%initCondType.eq.4) then
                print*, "You shouldn't put a shpere in and unequal box!"
                stop 1
            endif
        endif
        if (wlc_p%NBINX(1)*wlc_p%NBINX(2)*wlc_p%NBINX(3).ne.wlc_p%NBIN) then
            print*, "error in mcsim. Wrong number of bins"
            stop 1
        endif
        if (wlc_p%NT.ne.wlc_p%nMpP*wlc_p%NP*wlc_p%nBpM) then
            print*, "error in mcsim.  NT=",wlc_p%NT," nMpP=",wlc_p%nMpP," NP=",wlc_p%NP," nBpM=",wlc_p%nBpM
            stop 1
        endif
        if (wlc_p%NB.ne.wlc_p%nMpP*wlc_p%nBpM) then
            print*, "error in mcsim.  NB=",wlc_p%NB," nMpP=",wlc_p%nMpP," nBpM=",wlc_p%nBpM
            stop 1
        endif
        if (wlc_p%NNoInt.gt.wlc_p%indStartRepAdapt) then
            print*, "error in mcsim. don't run adapt without int on"
            stop 1
        endif
        if (wlc_p%NNoInt.gt.wlc_p%N_CHI_ON) then
            print*, "error in mcsim. Can't have chi without int on"
            stop 1
        endif
        if (wlc_p%NNoInt.gt.wlc_p%N_KAP_ON) then
            print*, "error in mcsim. Can't have kap without int on"
            stop 1
        endif
    end subroutine


    subroutine get_input_from_file(infile, wlc_p, wlc_d)
    ! Based on Elena's readkeys subroutine
        IMPLICIT NONE
        type(wlcsim_params), intent(out) :: wlc_p
        type(wlcsim_data), intent(out) :: wlc_d
        character(1024), intent(in) :: infile  ! file with parameters

        call set_param_defaults(wlc_p)

        call tweak_param_defaults(wlc_p, wlc_d)

        call read_from_file(infile, wlc_p)

        ! get derived parameters that aren't directly input from file
        call getpara(wlc_p)

        call idiot_checks(wlc_p)

    end subroutine


    subroutine initialize_wlcsim_data(wlc_d)
        implicit none
        type(wlcsim_data), intent(out)   :: wlc_d

        ! initialize energies to zero
        wlc_d%EElas=0.0_dp
        wlc_d%ECouple=0.0_dp
        wlc_d%ebind=0.0_dp
        wlc_d%EKap=0.0_dp
        wlc_d%ECHI=0.0_dp
        wlc_d%EField=0.0_dp
        wlc_d%x_mu=0.0_dp
        wlc_d%x_Field=0.0_dp
        wlc_d%x_couple=0.0_dp
        wlc_d%x_Kap=0.0_dp
        wlc_d%x_Chi=0.0_dp


    end subroutine initialize_wlcsim_data


    subroutine printDescription(wlcsim_p)
        IMPLICIT NONE
        type(wlcsim_params), intent(in) :: wlcsim_p
        print*, "---------------System Description---------------"
        print*, "Bead variables:"
        print*, " Total number of beads, NT=", wlcsim_p%NT
        print*, " Number of beads in a polymer, NB=", wlcsim_p%NB
        print*, " Number of monomers in a polymer, N=", wlcsim_p%nMpP
        print*, " Number of polymers, NP=",wlcsim_p%NP
        print*, " Number of beads in a monomer, nBpM=", wlcsim_p%nBpM
        print*, " fraction Methalated", wlcsim_p%F_METH
        print*, " LAM_METH", wlcsim_p%LAM_METH
        print*, "Length and volume Variables:"
        print*, " persistance length =",(wlcsim_p%L0/(2.0_dp*wlcsim_p%EPS))
        print*, " lbox=", wlcsim_p%lbox(1), wlcsim_p%lbox(2), wlcsim_p%lbox(3)
        print*, " Number of bins in x direction", &
                wlcsim_p%NBINX(1), wlcsim_p%NBINX(2),wlcsim_p%NBINX(3)
        print*, " Number of bins", wlcsim_p%NBIN
        print*, " spatial descritation dbin=",wlcsim_p%dbin
        print*, " L0=", wlcsim_p%L0
        print*, " volume fraction polymer =", wlcsim_p%Fpoly
        print*, " bead volume V=", wlcsim_p%beadVolume
        print*, "Energy Variables"
        print*, " elasticity EPS =", wlcsim_p%EPS
        print*, " solvent-polymer CHI =",wlcsim_p%CHI
        print*, " compression cof, KAP =", wlcsim_p%KAP
        print*, " field strength, hA =", wlcsim_p%hA

        print*, " -energy of binding unmethalated ", wlcsim_p%EU," more positive for favorable binding"
        print*, " -energy of binding methalated",wlcsim_p%EM
        print*, " HP1_Binding energy parameter", wlcsim_p%HP1_Bind
        print*, " chemical potential of HP1", wlcsim_p%mu
        print*, "Other:"
        print*, " confinetype:",wlcsim_p%confinetype
        print*, " initCondType:",wlcsim_p%initCondType
        print*, "---------------------------------------------"

    end subroutine
    subroutine wlcsim_params_allocate(wlcsim_p,wlcsim_d)
        IMPLICIT NONE
        type(wlcsim_params), intent(in) :: wlcsim_p
        type(wlcsim_data), intent(out) :: wlcsim_d
        integer NT  ! total number of beads
        integer NBIN ! total number of bins
        NT=wlcsim_p%NT
        NBIN=wlcsim_p%NBIN

        ! Quinn's sanity checks are removed. We should allow user to allocate as much as they want...
!         if ((NT.GT.200000).OR.(NT.lt.1)) then
!             print*, "Tried to allocate ", NT," beads in wlcsim_params_allocate"
!             stop 1
!         endif
!         if ((NBIN.GT.20000).or.(NBIN.lt.1)) then
!             print*, "Tried to allocate ",NBIN," bins in wlcsim_params_allocate"
!             stop 1
!         endif
        ALLOCATE(wlcsim_d%R(NT,3))
        ALLOCATE(wlcsim_d%U(NT,3))
        Allocate(wlcsim_d%RP(NT,3))
        Allocate(wlcsim_d%UP(NT,3))
        ALLOCATE(wlcsim_d%AB(NT))   !Chemical identity aka binding state
        ALLOCATE(wlcsim_d%ABP(NT))   !Chemical identity aka binding state
        ALLOCATE(wlcsim_d%METH(NT)) !Underlying methalation profile
        ALLOCATE(wlcsim_d%PHIA(NBIN))
        ALLOCATE(wlcsim_d%PHIB(NBIN))
        ALLOCATE(wlcsim_d%DPHIA(NBIN))
        ALLOCATE(wlcsim_d%DPHIB(NBIN))
        ALLOCATE(wlcsim_d%Vol(NBIN))
        Allocate(wlcsim_d%indPHI(NBIN))
        Allocate(wlcsim_d%PhiH(NBIN))

    end subroutine

    subroutine tweak_param_defaults(wlcsim_p, wlcsim_d)
        IMPLICIT NONE
        type(wlcsim_params), intent(inout) :: wlcsim_p
        type(wlcsim_data), intent(inout) :: wlcsim_d
        integer mctype ! type of move
    !   Edit the following to optimize wlcsim_p performance
        !  Monte-Carlo simulation parameters
        wlcsim_d%MCAMP(1)=0.5_dp*PI
        wlcsim_d%MCAMP(2)=0.3_dp*wlcsim_p%L0
        wlcsim_d%MCAMP(3)=0.5_dp*PI
        wlcsim_d%MCAMP(4)=0.5_dp*PI
        wlcsim_d%MCAMP(5)=0.5_dp*PI
        wlcsim_d%MCAMP(6)=5.0_dp*wlcsim_p%L0
        wlcsim_d%MCAMP(7)=nan
        wlcsim_d%MCAMP(8)=nan
        wlcsim_d%MCAMP(9)=nan
        wlcsim_d%MCAMP(10)=nan
        !switches to turn on various types of moves
        wlcsim_p%MOVEON(1)=1  ! crank-shaft move
        wlcsim_p%MOVEON(2)=1  ! slide move
        wlcsim_p%MOVEON(3)=1  ! pivot move
        wlcsim_p%MOVEON(4)=1  ! rotate move
        wlcsim_p%MOVEON(5)=0  ! full chain rotation
        wlcsim_p%MOVEON(6)=0  ! full chain slide
        wlcsim_p%MOVEON(7)=1  ! Change in Binding state
        wlcsim_p%MOVEON(8)=0  ! Chain flip
        wlcsim_p%MOVEON(9)=0  ! Chain exchange
        wlcsim_p%MOVEON(10)=0 ! Reptation

        !     Initial segment window for wlcsim_p moves
        wlcsim_d%WindoW(1)=15.0_dp ! used to be N*G
        wlcsim_d%WindoW(2)=15.0_dp ! used to be N*G
        wlcsim_d%WindoW(3)=15.0_dp ! used to be N*G
        wlcsim_d%WindoW(4)=1.0_dp
        wlcsim_d%WindoW(5)=dble(wlcsim_p%nMpP*wlcsim_p%nBpM)
        wlcsim_d%WindoW(6)=dble(wlcsim_p%nMpP*wlcsim_p%nBpM)
        wlcsim_d%WindoW(7)=15.0_dp ! used to be N*G
        wlcsim_d%WindoW(8)=dble(wlcsim_p%nMpP*wlcsim_p%nBpM)
        wlcsim_d%WindoW(9)=dble(wlcsim_p%nMpP*wlcsim_p%nBpM)
        wlcsim_d%WindoW(9)=1.0_dp

        !    Maximum window size (large windows are expensive)
        wlcsim_p%MAXWindoW(1)=dble(min(150,wlcsim_p%NB))
        wlcsim_p%MAXWindoW(2)=dble(min(150,wlcsim_p%NB))
        wlcsim_p%MAXWindoW(3)=dble(min(150,wlcsim_p%NB))
        wlcsim_p%MAXWindoW(4)=nan
        wlcsim_p%MAXWindoW(5)=nan
        wlcsim_p%MAXWindoW(6)=nan
        wlcsim_p%MAXWindoW(7)=dble(min(4,wlcsim_p%NB))
        wlcsim_p%MAXWindoW(8)=nan
        wlcsim_p%MAXWindoW(9)=nan
        wlcsim_p%MAXWindoW(9)=nan ! need to chaige code to allow >1

        wlcsim_p%MINWindoW(1)=dble(min(4,wlcsim_p%NB))
        wlcsim_p%MINWindoW(2)=dble(min(4,wlcsim_p%NB))
        wlcsim_p%MINWindoW(3)=dble(min(4,wlcsim_p%NB))
        wlcsim_p%MINWindoW(4)=nan
        wlcsim_p%MINWindoW(5)=nan
        wlcsim_p%MINWindoW(6)=nan
        wlcsim_p%MINWindoW(7)=dble(min(4,wlcsim_p%NB))
        wlcsim_p%MINWindoW(8)=nan
        wlcsim_p%MINWindoW(9)=nan
        wlcsim_p%MINWindoW(10)=nan
        do mctype=1,wlcsim_p%movetypes
            wlcsim_p%winTarget(mctype)=8.0_dp
        enddo

        wlcsim_p%MINAMP(1)=0.1_dp*PI
        wlcsim_p%MINAMP(2)=0.2_dp*wlcsim_p%L0
        wlcsim_p%MINAMP(3)=0.2_dp*PI
        wlcsim_p%MINAMP(4)=0.2_dp*PI
        wlcsim_p%MINAMP(5)=0.05_dp*PI
        wlcsim_p%MINAMP(6)=0.2_dp*wlcsim_p%L0
        wlcsim_p%MINAMP(7)=nan
        wlcsim_p%MINAMP(8)=nan
        wlcsim_p%MINAMP(9)=nan
        wlcsim_p%MINAMP(10)=nan

        wlcsim_p%MAXAMP(1)=1.0_dp*PI
        wlcsim_p%MAXAMP(2)=1.0_dp*wlcsim_p%L0
        wlcsim_p%MAXAMP(3)=1.0_dp*PI
        wlcsim_p%MAXAMP(4)=1.0_dp*PI
        wlcsim_p%MAXAMP(5)=1.0_dp*PI
        wlcsim_p%MAXAMP(6)=0.1*wlcsim_p%lbox(1)
        wlcsim_p%MAXAMP(7)=nan
        wlcsim_p%MAXAMP(8)=nan
        wlcsim_p%MAXAMP(9)=nan
        wlcsim_p%MAXAMP(9)=nan

        do mctype=1,wlcsim_p%movetypes
            wlcsim_p%NADAPT(mctype)=1000 ! adapt after at most 1000 steps
            wlcsim_p%PDESIRE(mctype)=0.5_dp ! Target
            wlcsim_d%SUCCESS(mctype)=0
            wlcsim_d%PHIT(mctype)=0.0_dp
        enddo
    end subroutine

    subroutine wlcsim_params_recenter(wlcsim_p,wlcsim_d)
    !  Prevents drift in periodic BC
        IMPLICIT NONE
        type(wlcsim_params), intent(in) :: wlcsim_p
        type(wlcsim_data), intent(inout) :: wlcsim_d
        integer IB, I, J   ! Couners
        real(dp) R0(3)  ! Offset to move by
        IB=1
        do I=1,wlcsim_p%NP
        R0(1)=nint(wlcsim_d%R(IB,1)/wlcsim_p%lbox(1)-0.5_dp)*wlcsim_p%lbox(1)
        R0(2)=nint(wlcsim_d%R(IB,2)/wlcsim_p%lbox(2)-0.5_dp)*wlcsim_p%lbox(2)
        R0(3)=nint(wlcsim_d%R(IB,3)/wlcsim_p%lbox(3)-0.5_dp)*wlcsim_p%lbox(3)
        if (abs(R0(1)*R0(2)*R0(3)) .gt. 0.0001_dp) then
            do J=1,wlcsim_p%NB
                wlcsim_d%R(IB,1)=wlcsim_d%R(IB,1)-R0(1)
                wlcsim_d%R(IB,2)=wlcsim_d%R(IB,2)-R0(2)
                wlcsim_d%R(IB,3)=wlcsim_d%R(IB,3)-R0(3)
                IB=IB+1
            enddo
        endif
        enddo
    end subroutine
    subroutine wlcsim_params_printEnergies(wlcsim_d)
    ! For realtime feedback on wlcsim_p simulation
        IMPLICIT NONE
        type(wlcsim_data), intent(in) :: wlcsim_d
        print*, "ECouple:", wlcsim_d%ECouple
        print*, "Bending energy", wlcsim_d%EELAS(1)
        print*, "Par compression energy", wlcsim_d%EELAS(2)
        print*, "Shear energy", wlcsim_d%EELAS(3)
        print*, "ECHI", wlcsim_d%ECHI
        print*, "EField", wlcsim_d%EField
        print*, "EKAP", wlcsim_d%EKAP
        print*, "ebind", wlcsim_d%ebind
    end subroutine
    subroutine wlcsim_params_printPhi(wlcsim_p,wlcsim_d)
    ! prints densities for trouble shooting
        IMPLICIT NONE
        type(wlcsim_params), intent(in) :: wlcsim_p
        type(wlcsim_data), intent(in) :: wlcsim_d
        integer I
        real(dp) EKap, ECouple, EChi,VV, PHIPOly
        print*,"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        print*, " PHIA  | PHIB  | PPoly |  Vol  | EKap  | EChi  |ECouple|"
        do I=1,wlcsim_p%NBIN
            VV=wlcsim_d%Vol(I)
            if (VV.le.0.1_dp) cycle
            PHIPOLY=wlcsim_d%PHIA(I)+wlcsim_d%PHIB(I)
            EChi=VV*(wlcsim_p%CHI/wlcsim_p%beadVolume)*PHIPoly*(1.0_dp-PHIPoly)
            ECouple=VV*wlcsim_p%HP1_Bind*(wlcsim_d%PHIA(I))**2
            if(PHIPoly.GT.1.0_dp) then
            EKap=VV*(wlcsim_p%KAP/wlcsim_p%beadVolume)*(PHIPoly-1.0_dp)**2
            else
            cycle
            EKap=0.0_dp
            endif
            write(*,"(4f8.4,3f8.1)") wlcsim_d%PHIA(I), wlcsim_d%PHIB(I), &
                                wlcsim_d%PHIA(I)+wlcsim_d%PHIB(I),wlcsim_d%Vol(I),&
                                EKap,EChi,ECouple
        enddo
        print*,"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    end subroutine
    subroutine wlcsim_params_printWindowStats(wlcsim_p, wlcsim_d)
    ! For realtime feedback on adaptation
        IMPLICIT NONE
        type(wlcsim_params), intent(in) :: wlcsim_p
        type(wlcsim_data), intent(in) :: wlcsim_d
        integer I ! counter
        I=0
        print*, "Succes | MCAMP | WindoW| type "
        do I=1,wlcsim_p%movetypes
            if (wlcsim_p%MOVEON(i).eq.1) then
                write(*,"(f8.5,2f8.2,1I8)") wlcsim_d%phit(i), wlcsim_d%MCAMP(i),  wlcsim_d%WindoW(i), i
            endif
        enddo
        return
    end subroutine
    subroutine wlcsim_params_LoadField(wlcsim_p,wlcsim_d,fileName)
        IMPLICIT NONE
        type(wlcsim_params), intent(in) :: wlcsim_p
        type(wlcsim_data), intent(inout) :: wlcsim_d
        integer I
        character(MAXFILENAMELEN) fileName ! file name to load from
        open (unit = 1, file = fileName, status = 'OLD')
        do I=1,wlcsim_p%NBIN
            read(1,*) wlcsim_d%PHIH(I)
        enddo
        return
    end subroutine
    subroutine wlcsim_params_MakeField(wlcsim_p,wlcsim_d)
        IMPLICIT NONE
        type(wlcsim_params), intent(in) :: wlcsim_p
        type(wlcsim_data), intent(inout) :: wlcsim_d
        integer indBIN  ! index of bin
        integer IX,IY,IZ ! bin corrdinates

        do IX=1,wlcsim_p%NBINX(1)
            do IY=1,wlcsim_p%NBINX(2)
                do IZ=1,wlcsim_p%NBINX(3)
                    indBIN=IX+&
                        (IY-1)*wlcsim_p%NBINX(1)+&
                        (IZ-1)*wlcsim_p%NBINX(1)*wlcsim_p%NBINX(2)
                    wlcsim_d%PHIH(indBIN)=dsin(wlcsim_p%k_field*wlcsim_p%dbin*dble(IX))
                enddo
            enddo
        enddo
        return
    end subroutine
    subroutine wlcsim_params_loadAB(wlcsim_p,wlcsim_d,fileName)
    ! Loads AB for file...has not been tested
        IMPLICIT NONE
        type(wlcsim_params), intent(in) :: wlcsim_p
        type(wlcsim_data), intent(inout) :: wlcsim_d
        character(MAXFILENAMELEN), intent(in) :: fileName ! file name to load from
        integer IB, I, J ! counters
        open (unit = 1, file = fileName, status = 'OLD')
        IB=1
        do I=1,wlcsim_p%NP
        do J=1,wlcsim_p%NB
            read(1,"(I2)") wlcsim_d%AB(IB)
            IB=IB+1
            enddo
        enddo
        close(1)
    end subroutine
    subroutine wlcsim_params_saveR(wlcsim_p,wlcsim_d,fileName,repeatingBC)
    ! Writes R and AB to file for analysis
    ! Rx  Ry  Rz AB
        IMPLICIT NONE
        integer, intent(in) :: repeatingBC  ! 1 for reapeating boundary conditions
        integer I,J,IB  ! counters
        type(wlcsim_params), intent(in) :: wlcsim_p
        type(wlcsim_data), intent(in) :: wlcsim_d
        character(MAXFILENAMELEN), intent(in) :: fileName
        character(MAXFILENAMELEN) fullName
        fullName=  trim(fileName) // trim(wlcsim_p%repSuffix)
        fullName=trim(fullName)
        open (unit = 1, file = fullName, status = 'NEW')
        IB=1
        if (repeatingBC.eq.1) then
            do I=1,wlcsim_p%NP
                do J=1,wlcsim_p%NB
                        WRITE(1,"(3f10.3,I2)") &
                            wlcsim_d%R(IB,1)-0.*nint(wlcsim_d%R(IB,1)/wlcsim_p%lbox(1)-0.5_dp)*wlcsim_p%lbox(1), &
                            wlcsim_d%R(IB,2)-0.*nint(wlcsim_d%R(IB,2)/wlcsim_p%lbox(2)-0.5_dp)*wlcsim_p%lbox(2), &
                            wlcsim_d%R(IB,3)-0.*nint(wlcsim_d%R(IB,3)/wlcsim_p%lbox(3)-0.5_dp)*wlcsim_p%lbox(3), &
                            wlcsim_d%AB(IB)
                    IB=IB+1
                enddo
            enddo
            print*, "Error in wlcsim_params_saveR"
            print*, "Are you sure you want reapeating BC"
            stop 1
        else
            do I=1,wlcsim_p%NP
                do J=1,wlcsim_p%NB
                    if (wlcsim_p%simtype.eq.0) then
                        WRITE(1,"(3f10.3,I2)") wlcsim_d%R(IB,1),wlcsim_d%R(IB,2),wlcsim_d%R(IB,3),wlcsim_d%AB(IB)
                    else
                        WRITE(1,"(3f10.3,I2)") wlcsim_d%R(IB,1),wlcsim_d%R(IB,2),wlcsim_d%R(IB,3),wlcsim_d%AB(IB), wlcsim_d%METH(IB)
                    endif
                    IB=IB+1
                enddo
            enddo
        endif
        close(1)
    end subroutine
    subroutine wlcsim_params_savePHI(wlcsim_p,wlcsim_d,fileName)
    ! Saves PHIA and PHIB to file for analysis
        IMPLICIT NONE
        integer I  ! counters
        type(wlcsim_params), intent(in) :: wlcsim_p
        type(wlcsim_data), intent(in) :: wlcsim_d
        character(MAXFILENAMELEN), intent(in) :: fileName
        character(MAXFILENAMELEN) fullName
        fullName=  trim(fileName) // trim(wlcsim_p%repSuffix)
        open (unit = 1, file = fullName, status = 'NEW')
        do I=1,wlcsim_p%NBIN
            WRITE(1,"(2f7.2)") wlcsim_d%PHIA(I),wlcsim_d%PHIB(I)
        enddo
        close(1)
    end subroutine
    subroutine wlcsim_params_saveU(wlcsim_p,wlcsim_d,fileName)
    ! Saves U to ASCII file for analisys
        IMPLICIT NONE
        integer I,J,IB  ! counters
        type(wlcsim_params), intent(in) :: wlcsim_p
        type(wlcsim_data), intent(in) :: wlcsim_d
        character(MAXFILENAMELEN), intent(in) :: fileName
        character(MAXFILENAMELEN) fullName
        fullName=  trim(fileName) // trim(wlcsim_p%repSuffix)
        open (unit = 1, file = fullName, status = 'NEW')
        IB=1
        do I=1,wlcsim_p%NP
            do J=1,wlcsim_p%NB
                WRITE(1,"(3f8.3,2I2)") wlcsim_d%U(IB,1),wlcsim_d%U(IB,2),wlcsim_d%U(IB,3)
                IB=IB+1
            enddo
        enddo
        close(1)
    end subroutine
    subroutine wlcsim_params_saveparameters(wlcsim_p,fileName)
    ! Write a number of parameters ASCII variables to file for reccords
        IMPLICIT NONE
        type(wlcsim_params), intent(in) :: wlcsim_p
        character(MAXFILENAMELEN), intent(in) :: fileName
        character(MAXFILENAMELEN) fullName
        fullName=  trim(fileName) // trim(wlcsim_p%repSuffix)
        open (unit =1, file = fullName, status = 'NEW')
            WRITE(1,"(I8)") wlcsim_p%NT ! 1 Number of beads in simulation
            WRITE(1,"(I8)") wlcsim_p%nMpP  ! 2 Number of monomers in a polymer
            WRITE(1,"(I8)") wlcsim_p%NB ! 3 Number of beads in a polymer
            WRITE(1,"(I8)") wlcsim_p%NP ! 4 Number of polymers in simulation
            WRITE(1,"(I8)") wlcsim_p%NT ! 5 Number of beads in simulation
            WRITE(1,"(I8)") wlcsim_p%nBpM  ! 6 Number of beads per monomer

            WRITE(1,"(f10.5)") wlcsim_p%L0    ! Equilibrium segment length
            WRITE(1,"(f10.5)") wlcsim_p%CHI  ! 8  initail CHI parameter value
            WRITE(1,"(f10.5)") wlcsim_p%Fpoly ! Fraction polymer
            WRITE(1,"(f10.5)") wlcsim_p%lbox(1)  ! 10 Lenth of box
            WRITE(1,"(f10.5)") wlcsim_p%EU    ! Energy unmethalated
            WRITE(1,"(f10.5)") wlcsim_p%EM    ! 12 Energy methalated
            WRITE(1,"(f10.5)") wlcsim_p%HP1_Bind ! Energy of HP1 binding
            WRITE(1,"(f10.5)") (wlcsim_p%L0/wlcsim_p%EPS) ! 14 Khun lenth
            WRITE(1,"(A)") "-999"  ! for historic reasons
            WRITE(1,"(f10.5)") wlcsim_p%F_METH  ! methalation fraction
            WRITE(1,"(f10.5)") wlcsim_p%LAM_METH  ! methalation lambda
        close(1)
    end subroutine
    subroutine wlcsim_params_appendEnergyData(wlcsim_p, wlcsim_d, fileName)
    ! print Energy data
        IMPLICIT NONE
        type(wlcsim_params), intent(in) :: wlcsim_p
        type(wlcsim_data), intent(in) :: wlcsim_d
        LOGICAL isfile
        character(MAXFILENAMELEN), intent(in) :: fileName
        character(MAXFILENAMELEN) fullName
        fullName=  trim(fileName) // trim(wlcsim_p%repSuffix)
        inquire(file = fullName, exist=isfile)
        if (isfile) then
            open (unit = 1, file = fullName, status ='OLD', POSITION="append")
        else
            open (unit = 1, file = fullName, status = 'new')
            WRITE(1,*) "ind | id |",&
                       " ebend  | eparll | EShear | ECoupl | E Kap  | E Chi  |",&
                       " EField | ebind  |  x_Mu  | Couple |  Chi   |  mu    |",&
                       "  Kap   | Field  |"
        endif
        WRITE(1,"(2I5, 9f9.1,5f9.4)") wlcsim_d%ind, wlcsim_d%id, &
            wlcsim_d%EELAS(1), wlcsim_d%EELAS(2), wlcsim_d%EELAS(3), wlcsim_d%ECouple, &
            wlcsim_d%EKap, wlcsim_d%ECHI, wlcsim_d%EField, wlcsim_d%ebind, wlcsim_d%x_Mu, &
            wlcsim_p%HP1_Bind*wlcsim_p%Couple_on, wlcsim_p%CHI*wlcsim_p%CHI_ON, wlcsim_p%mu, wlcsim_p%KAP*wlcsim_p%KAP_ON,&
            wlcsim_p%hA
        close(1)
    end subroutine
    subroutine wlcsim_params_appendAdaptData(wlcsim_p, wlcsim_d, fileName)
    ! Appends wlcsim_p move adaptation data to the file
        IMPLICIT NONE
        type(wlcsim_params), intent(in) :: wlcsim_p
        type(wlcsim_data), intent(in) :: wlcsim_d
        LOGICAL isfile
        character(MAXFILENAMELEN), intent(in) :: fileName
        character(MAXFILENAMELEN) fullName
        fullName=  trim(fileName) // trim(wlcsim_p%repSuffix)
        inquire(file = fullName, exist=isfile)
        if (isfile) then
            open (unit = 1, file = fullName, status ='OLD', POSITION="append")
        else
            open (unit = 1, file = fullName, status = 'new')
            WRITE(1,*) "ind| id|",&
                       " WIN 1 | AMP 1 | SUC 1 | WIN 2 | AMP 2 | SUC 2 |",&
                       " WIN 3 | AMP 3 | SUC 3 | ON  4 | AMP 4 | SUC 4 |",&
                       " ON  5 | AMP 5 | SUC 5 | ON  6 | AMP 6 | SUC 6 |",&
                       " ON  7 | SUC 7 | ON  8 | SUC 8 |", &
                       " ON  9 | SUC 9 | ON 10 | SUC 10|"
        endif
        WRITE(1,"(2I4,26f8.3)") wlcsim_d%ind,wlcsim_d%id,&
            real(wlcsim_d%WindoW(1)),wlcsim_d%MCAMP(1),wlcsim_d%PHIT(1), &
            real(wlcsim_d%WindoW(2)),wlcsim_d%MCAMP(2),wlcsim_d%PHIT(2), &
            real(wlcsim_d%WindoW(3)),wlcsim_d%MCAMP(3),wlcsim_d%PHIT(3), &
            real(wlcsim_p%MOVEON(4)),wlcsim_d%MCAMP(4),wlcsim_d%PHIT(4), &
            real(wlcsim_p%MOVEON(5)),wlcsim_d%MCAMP(5),wlcsim_d%PHIT(5), &
            real(wlcsim_p%MOVEON(6)),wlcsim_d%MCAMP(6),wlcsim_d%PHIT(6), &
            real(wlcsim_p%MOVEON(7)),wlcsim_d%PHIT(7), &
            real(wlcsim_p%MOVEON(8)),wlcsim_d%PHIT(8), &
            real(wlcsim_p%MOVEON(9)),wlcsim_d%PHIT(9), &
            real(wlcsim_p%MOVEON(10)),wlcsim_d%PHIT(10)
        close(1)
    end subroutine
    subroutine wlcsim_params_writebinary(wlcsim_p,wlcsim_d,baceName)
    !    This function writes the contence of the structures wlcsim_p and wlcsim_d
    !  to a binary file.  if you add more variables to wlcsim_d you need to
    !  a seperate write command for them as it is not possible to write
    !  a structure with allocatables to a binar file.
    !    The contence are stored in
    !     baceName//'R'
    !     baceName//'U'
    !     etc.
        IMPLICIT NONE
        integer sizeOftype         ! for binary saving
        type(wlcsim_params), intent(in) :: wlcsim_p             ! to be save or filled
        type(wlcsim_data), intent(in) :: wlcsim_d             ! to be save or filled
        CHARACTER(LEN=16), intent(in) :: baceName ! for example 'record/'
        CHARACTER(LEN=16) fileName ! fileName
        CHARACTER(LEN=16) sufix    ! end of file name
        LOGICAL exists    ! does file already exist?

        !  ------parameters -----

        sizeOftype=int(SIZEOF(wlcsim_p))
        sufix='parameters'
        fileName=trim(baceName) // trim(sufix)
        inquire(file=fileName,exist=exists)
        if(exists) then
            open(unit=1,file=fileName, status='old', &
                form='unformatted',access='direct',recl=sizeOftype)
        else
            open(unit=1,file=fileName, status='new', &
                form='unformatted',access='direct',recl=sizeOftype)
        endif
        write(1,rec=1) wlcsim_p
        close(1)

        ! -------- R --------

        sizeOftype=int(SIZEOF(wlcsim_d%R))
        sufix='R'
        fileName=trim(baceName) // trim(sufix)
        inquire(file=fileName,exist=exists)
        if(exists) then
            open(unit=1,file=fileName, status='old', &
                form='unformatted',access='direct',recl=sizeOftype)
        else
            open(unit=1,file=fileName, status='new', &
                form='unformatted',access='direct',recl=sizeOftype)
        endif
        write(1,rec=1) wlcsim_d%R
        close(1)

        ! -------- U --------

        sizeOftype=int(SIZEOF(wlcsim_d%U))
        sufix='U'
        fileName=trim(baceName) // trim(sufix)
        inquire(file=fileName,exist=exists)
        if(exists) then
            open(unit=1,file=fileName, status='old', &
                form='unformatted',access='direct',recl=sizeOftype)
        else
            open(unit=1,file=fileName, status='new', &
                form='unformatted',access='direct',recl=sizeOftype)
        endif
        write(1,rec=1) wlcsim_d%U
        close(1)

        ! -------- AB --------

        sizeOftype=int(SIZEOF(wlcsim_d%AB))
        sufix='AB'
        fileName=trim(baceName) // trim(sufix)
        inquire(file=fileName,exist=exists)
        if(exists) then
            open(unit=1,file=fileName, status='old', &
                form='unformatted',access='direct',recl=sizeOftype)
        else
            open(unit=1,file=fileName, status='new', &
                form='unformatted',access='direct',recl=sizeOftype)
        endif
        write(1,rec=1) wlcsim_d%AB
        close(1)

        ! -------- Vol --------

        sizeOftype=int(SIZEOF(wlcsim_d%Vol))
        sufix='Vol'
        fileName=trim(baceName) // trim(sufix)
        inquire(file=fileName,exist=exists)
        if(exists) then
            open(unit=1,file=fileName, status='old', &
                form='unformatted',access='direct',recl=sizeOftype)
        else
            open(unit=1,file=fileName, status='new', &
                form='unformatted',access='direct',recl=sizeOftype)
        endif
        write(1,rec=1) wlcsim_d%Vol
        close(1)
    end subroutine

    subroutine wlcsim_params_readBindary(wlcsim_p,wlcsim_d,baceName)
    ! This function reads what wlcsim_params_writebinary writes and
    ! stores it to wlcsim_p and wlcsim_d.  Be sure to allocate wlcsim_d before
    ! calling this command.
        IMPLICIT NONE
        integer sizeOftype         ! for binary saving
        type(wlcsim_params) wlcsim_p             ! to be save or filled
        type(wlcsim_data) wlcsim_d             ! to be save or filled
        CHARACTER(LEN=16) baceName ! for example 'record/'
        CHARACTER(LEN=16) fileName ! fileName
        CHARACTER(LEN=16) sufix    ! end of file name
        LOGICAL exists    ! does file already exist?

        !  ------parameters -----

        sizeOftype=int(SIZEOF(wlcsim_p))
        sufix='parameters'
        fileName=trim(baceName) // trim(sufix)
        inquire(file=fileName,exist=exists)
        if(exists) then
            open(unit=1,file=fileName, status='old', &
                form='unformatted',access='direct',recl=sizeOftype)
        else
            print*, 'Error in wlcsim_params_readBinary. file ',fileName,'does not exist'
            stop 1
        endif
        read(1,rec=1) wlcsim_p
        close(1)

        ! -------- R --------

        sizeOftype=int(SIZEOF(wlcsim_d%R))
        sufix='R'
        fileName=trim(baceName) // trim(sufix)
        inquire(file=fileName,exist=exists)
        if(exists) then
            open(unit=1,file=fileName, status='old', &
                form='unformatted',access='direct',recl=sizeOftype)
        else
            print*, 'Error in wlcsim_params_readBinary. file ',fileName,'does not exist'
            stop 1
        endif
        read(1,rec=1) wlcsim_d%R
        close(1)

        ! -------- U --------

        sizeOftype=int(SIZEOF(wlcsim_d%U))
        sufix='U'
        fileName=trim(baceName) // trim(sufix)
        inquire(file=fileName,exist=exists)
        if(exists) then
            open(unit=1,file=fileName, status='old', &
                form='unformatted',access='direct',recl=sizeOftype)
        else
            print*, 'Error in wlcsim_params_readBinary. file ',fileName,'does not exist'
            stop 1
        endif
        read(1,rec=1) wlcsim_d%U
        close(1)

        ! -------- AB --------

        sizeOftype=int(SIZEOF(wlcsim_d%AB))
        sufix='AB'
        fileName=trim(baceName) // trim(sufix)
        inquire(file=fileName,exist=exists)
        if(exists) then
            open(unit=1,file=fileName, status='old', &
                form='unformatted',access='direct',recl=sizeOftype)
        else
            print*, 'Error in wlcsim_params_readBinary. file ',fileName,'does not exist'
            stop 1
        endif
        read(1,rec=1) wlcsim_d%AB
        close(1)

        ! -------- Vol --------

        sizeOftype=int(SIZEOF(wlcsim_d%Vol))
        sufix='Vol'
        fileName=trim(baceName) // trim(sufix)
        inquire(file=fileName,exist=exists)
        if(exists) then
            open(unit=1,file=fileName, status='old', &
                form='unformatted',access='direct',recl=sizeOftype)
        else
            print*, 'Error in wlcsim_params_readBinary. file ',fileName,'does not exist'
            stop 1
        endif
        read(1,rec=1) wlcsim_d%Vol
        close(1)
    end subroutine
end module params
