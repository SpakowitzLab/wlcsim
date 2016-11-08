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

    IMPLICIT NONE

    private
    public :: nMoveTypes, dp, pi, wlcsim_params

    ! hardcoded params. will need to change if certain parts of code change
    integer, parameter :: nMoveTypes = 10 ! number of MC move types

    ! precision of simulations
    integer, parameter :: dp = real64 ! preferred over SELECTED_real_Kind(15,307)

    ! universal constants
    real(dp) :: pi = 4 * atan(1.0_dp) ! fully accurate, adaptive precision
    real(dp) :: nan = IEEE_VALUE(IEEE_QUIET_NAN)

    ! for all parameters that cannot change during individual simulations
    ! these are documented more thoroughly where they are read in (see the
    ! subroutine get_input_from_file), in the docs (TODO), and often the default
    ! values will help with understanding how the variable is used.
    type wlcsim_params
    !   Simulation parameters
        integer nT                ! Total number of beads  NT=NP*N*G
        integer nB                ! Number of beads in a polymer NB=N*G
        integer nMpP              ! Number of monomers (NOT BEADS!) in a polymer
        integer nBpM              ! number beads per monomer
        integer nP                ! Number of polymers
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

    !   for passing 1st order phase transition in (quinn/shifan's) random copolymer MC sims
        real(dp) kField   ! wave vector of applied sinusoidal field (used in PT to step around 1st order phase transition)
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
        integer confinetype       ! type of Boundary Conditions
        integer settype           ! initial condition type
        logical field_interactions ! field-based self interactions on
        logical intrapolymer_stick_crossing_enforced ! field-based self interactions on
        logical FRMCHEM           ! read initial chemical sequence from file
        logical FRMMETH           ! read initial methylation from file
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
        Character*16 repSufix    ! prefix for writing files
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

contains

subroutine set_param_defaults(wlc_p)
    implicit none
    type(wlcsim_params), intent(out) :: wlc_p
    ! file IO
    wlc_p%FRMfile=.FALSE.      ! don't load initial bead positions from file
    wlc_p%FRMMETH=.FALSE.      ! don't load initial "methyl" status from file
    wlc_p%FRMFIELD=.false.     ! don't load initial field values from file
    wlc_p%saveU=.true.         ! do save orientation vectors (makes restart of ssWLC possible)
    wlc_p%savePhi=.FALSE.      ! don't save A/B density per bin (not needed for restart)
    wlc_p%FRMCHEM=.FALSE.      ! don't load initial a/b states from file
    wlc_p%restart=.FALSE.      ! don't restart from previously saved simulation

    ! geometry options
    wlc_p%NP  =1               ! one polymer
    wlc_p%nB  =200             ! 200 beads per polymer
    wlc_p%nBpM = 10
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
    wlc_p%min_accept=0.05 ! if a move succeeds < 5% of the time, start using it only every reduce_move cycles

    ! timing options
    wlc_p%NStep=400000  ! number of simulation steps to take
    wlc_p%NNoInt=100    ! number of simulation steps before turning on interactions in Quinn's MC scheduler
    wlc_p%indMAX=200    ! 2000 steps per save point
    wlc_p%reduce_move=10 ! use moves that fall below the min_accept threshold only once every 10 times they would otherwise be used
    wlc_p%useSchedule=.False. ! use Quinn's scheduler to modify MC params halfway through the simulation
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
    character(fileNameLen), intent(in) :: infile
    integer pf = infileUnit
    open(unit=PF,file=fileName,status='OLD')

    ! read in the keywords one line at a time
    do
       CALL readLINE(PF,fileend,NITEMS)
       if (fileend.and.nitems.eq.0) EXIT

       ! skip empty lines
       if (NITEMS.EQ.0) cycle

       ! read in the keyword for this line
       CALL readA(WORD,CASESET=1)

       ! Skip any empty lines or any comment lines
       if (WORD(1:1).EQ.'#') cycle

       SELECT CASE(WORD) ! pick which keyword
       CASE('SETtype')
           Call readI(wlc_p%settype)
           ! settype      |  Discription
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
           Call reado(wlc_p%FRMCHEM) ! Initial chemical sequence from file
       CASE('FRMfile')
           call reado(wlc_p%FRMfile) ! read configuration from file
       CASE('FRMMETH')
           Call reado(wlc_p%FRMMETH) ! read methalation from file
       CASE('PTON')
           CALL reado(wlc_p%PTON) ! parallel Tempering on
       CASE('SAVE_U')
           Call reado(wlc_p%saveU)  ! save u vectors to file (every savepoint)
       CASE('SAVE_PHI')
           Call reado(wlc_p%savePhi) ! save Phi vectors to file (every savepoint)
       CASE('L0')
           Call readF(wlc_p%L0)  ! Equilibrium segment length
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
       CASE('G')
           Call readI(wlc_p%G) ! Beads per monomer
       CASE('N')
           CALL readI(wlc_p%N) ! Number of monomers in a polymer
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
           Call readF(wlc_p%V) ! Bead volume
       CASE('FA')
           Call readF(wlc_p%FA) ! Fraction of A beads (fraction bound)
       CASE('LAM')
           Call readF(wlc_p%LAM) ! Chemical correlation parameter
       CASE('EPS')
           Call readF(wlc_p%EPS) ! Elasticity l0/(2lp)
       CASE('CHI')
           Call readF(wlc_p%CHI) ! CHI parameter (definition depends on  hamiltoniaon
       CASE('H_A')
           Call readF(wlc_p%h_A) ! strength of externally applied field
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
           print*, "Error in MCvar_setparams.  Unidentified keyword:", &
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
        if (wlc_p%settype.eq.4) then
            print*, "You shouldn't put a shpere in and unequal box!"
            stop 1
        endif
    endif
    if (wlc_p%NBINX(1)*wlc_p%NBINX(2)*wlc_p%NBINX(3).ne.wlc_p%NBIN) then
        print*, "error in MCsim. Wrong number of bins"
        stop 1
    endif
    if (wlc_p%NT.ne.wlc_p%N*wlc_p%NP*wlc_p%G) then
        print*, "error in MCsim.  NT=",wlc_p%NT," N=",wlc_p%N," NP=",wlc_p%NP," G=",wlc_p%G
        stop 1
    endif
    if (wlc_p%NB.ne.wlc_p%N*wlc_p%G) then
        print*, "error in MCsim.  NB=",wlc_p%NB," N=",wlc_p%N," G=",wlc_p%G
        stop 1
    endif
    if (wlc_p%NNoInt.gt.wlc_p%indStartRepAdapt) then
        print*, "error in MCsim. don't run adapt without int on"
        stop 1
    endif
    if (wlc_p%NNoInt.gt.wlc_p%N_CHI_ON) then
        print*, "error in MCsim. Can't have chi without int on"
        stop 1
    endif
    if (wlc_p%NNoInt.gt.wlc_p%N_KAP_ON) then
        print*, "error in MCsim. Can't have kap without int on"
        stop 1
    endif
end subroutine


subroutine get_input_from_file(infile, wlc_p)
! Based on Elena's readkeys subroutine
    IMPLICIT NONE
    type(wlcsim_params), intent(out) :: wlc_p
    character(1024), intent(in) :: infile  ! file with parameters
    integer :: PF   ! input file unit
    LOGICAL :: fileend = .FALSE. ! done reading file?
    CHARACTER*100 :: WORD ! parameter name currently being read in
    integer :: NITEMS ! number of items on the line in the parameter file

    call set_param_defaults(wlc_p)

    call tweak_mc_defaults(wlc_p)

    call read_from_file(infile, wlc_p)

    ! get derived parameters that aren't directly input from file
    call getpara(wlc_p)

    call idiot_checks(wlc_p)

end subroutine


subroutine initialize_wlcsim_data(wlc_d, wlc_p)
    implicit none
    type(wlcsim_params), intent(in)  :: wlc_p
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


subroutine MCvar_printDescription(mc)
    IMPLICIT NONE
    type(MCvar), intent(in) :: mc
    print*, "---------------System Description---------------"
    print*, "Bead variables:"
    print*, " Total number of beads, NT=", mc%NT
    print*, " Number of beads in a polymer, NB=", mc%NB
    print*, " Number of monomers in a polymer, N=", mc%N
    print*, " Number of polymers, NP=",mc%NP
    print*, " Number of beads in a monomer, G=", mc%G
    print*, " fraction Methalated", mc%F_METH
    print*, " LAM_METH", mc%LAM_METH
    print*, "Length and volume Variables:"
    print*, " persistance length =",(mc%L0/(2.0_dp*mc%EPS))
    print*, " lbox=", mc%lbox(1), mc%lbox(2), mc%lbox(3)
    print*, " Number of bins in x direction", &
             mc%NBINX(1), mc%NBINX(2),mc%NBINX(3)
    print*, " Number of bins", mc%NBIN
    print*, " spatial descritation dbin=",mc%dbin
    print*, " L0=", mc%L0
    print*, " volume fraction polymer =", mc%Fpoly
    print*, " bead volume V=", mc%V
    print*, "Energy Variables"
    print*, " elasticity EPS =", mc%EPS
    print*, " solvent-polymer CHI =",mc%CHI
    print*, " compression cof, KAP =", mc%KAP
    print*, " field strength, h_A =", mc%h_A

    print*, " -energy of binding unmethalated ", mc%EU," more positive for favorable binding"
    print*, " -energy of binding methalated",mc%EM
    print*, " HP1_Binding energy parameter", mc%HP1_Bind
    print*, " chemical potential of HP1", mc%mu
    print*, "Other:"
    print*, " confinetype:",mc%confinetype
    print*, " settype:",mc%settype
    print*, "---------------------------------------------"

end subroutine
subroutine MCvar_allocate(mc,md)
    IMPLICIT NONE
    type(MCvar), intent(in) :: mc
    type(MCData), intent(out) :: md
    integer NT  ! total number of beads
    integer NBIN ! total number of bins
    NT=mc%NT
    NBIN=mc%NBIN

    if ((NT.GT.200000).OR.(NT.lt.1)) then
        print*, "Tried to allocate ", NT," beads in MCvar_allocate"
        stop 1
    endif
    if ((NBIN.GT.20000).or.(NBIN.lt.1)) then
        print*, "Tried to allocate ",NBIN," bins in MCvar_allocate"
        stop 1
    endif
    ALLOCATE(md%R(NT,3))
    ALLOCATE(md%U(NT,3))
    Allocate(md%RP(NT,3))
    Allocate(md%UP(NT,3))
    ALLOCATE(md%AB(NT))   !Chemical identity aka binding state
    ALLOCATE(md%ABP(NT))   !Chemical identity aka binding state
    ALLOCATE(md%METH(NT)) !Underlying methalation profile
    ALLOCATE(md%PHIA(NBIN))
    ALLOCATE(md%PHIB(NBIN))
    ALLOCATE(md%DPHIA(NBIN))
    ALLOCATE(md%DPHIB(NBIN))
    ALLOCATE(md%Vol(NBIN))
    Allocate(md%indPHI(NBIN))
    Allocate(md%PhiH(NBIN))

end subroutine

subroutine tweak_mc_defaults(mc)
    IMPLICIT NONE
    type(wlcsim_params), intent(inout) :: mc
    integer MCtype ! type of move
!   Edit the following to optimize MC performance
    !  Monte-Carlo simulation parameters
    mc%MCAMP(1)=0.5_dp*PI
    mc%MCAMP(2)=0.3_dp*mc%L0
    mc%MCAMP(3)=0.5_dp*PI
    mc%MCAMP(4)=0.5_dp*PI
    mc%MCAMP(5)=0.5_dp*PI
    mc%MCAMP(6)=5.0_dp*mc%L0
    mc%MCAMP(7)=NAND
    mc%MCAMP(8)=NAND
    mc%MCAMP(9)=NAND
    mc%MCAMP(10)=NAND
    !switches to turn on various types of moves
    mc%MOVEON(1)=1  ! crank-shaft move
    mc%MOVEON(2)=1  ! slide move
    mc%MOVEON(3)=1  ! pivot move
    mc%MOVEON(4)=1  ! rotate move
    mc%MOVEON(5)=0  ! full chain rotation
    mc%MOVEON(6)=0  ! full chain slide
    mc%MOVEON(7)=1  ! Change in Binding state
    mc%MOVEON(8)=0  ! Chain flip
    mc%MOVEON(9)=0  ! Chain exchange
    mc%MOVEON(10)=0 ! Reptation

    !     Initial segment window for MC moves
    mc%WindoW(1)=15.0_dp ! used to be N*G
    mc%WindoW(2)=15.0_dp ! used to be N*G
    mc%WindoW(3)=15.0_dp ! used to be N*G
    mc%WindoW(4)=1.0_dp
    mc%WindoW(5)=dble(mc%N*mc%G)
    mc%WindoW(6)=dble(mc%N*mc%G)
    mc%WindoW(7)=15.0_dp ! used to be N*G
    mc%WindoW(8)=dble(mc%N*mc%G)
    mc%WindoW(9)=dble(mc%N*mc%G)
    mc%WindoW(9)=1.0_dp

    !    Maximum window size (large windows are expensive)
    mc%MAXWindoW(1)=dble(min(150,mc%NB))
    mc%MAXWindoW(2)=dble(min(150,mc%NB))
    mc%MAXWindoW(3)=dble(min(150,mc%NB))
    mc%MAXWindoW(4)=NAND
    mc%MAXWindoW(5)=NAND
    mc%MAXWindoW(6)=NAND
    mc%MAXWindoW(7)=dble(min(4,mc%NB))
    mc%MAXWindoW(8)=NAND
    mc%MAXWindoW(9)=NAND
    mc%MAXWindoW(9)=NAND ! need to chaige code to allow >1

    mc%MINWindoW(1)=dble(min(4,mc%NB))
    mc%MINWindoW(2)=dble(min(4,mc%NB))
    mc%MINWindoW(3)=dble(min(4,mc%NB))
    mc%MINWindoW(4)=NAND
    mc%MINWindoW(5)=NAND
    mc%MINWindoW(6)=NAND
    mc%MINWindoW(7)=dble(min(4,mc%NB))
    mc%MINWindoW(8)=NAND
    mc%MINWindoW(9)=NAND
    mc%MINWindoW(10)=NAND
    do MCtype=1,mc%movetypes
        mc%winTarget(MCtype)=8.0_dp
    enddo

    mc%MINAMP(1)=0.1_dp*PI
    mc%MINAMP(2)=0.2_dp*mc%L0
    mc%MINAMP(3)=0.2_dp*PI
    mc%MINAMP(4)=0.2_dp*PI
    mc%MINAMP(5)=0.05_dp*PI
    mc%MINAMP(6)=0.2_dp*mc%L0
    mc%MINAMP(7)=NAND
    mc%MINAMP(8)=NAND
    mc%MINAMP(9)=NAND
    mc%MINAMP(10)=NAND

    mc%MAXAMP(1)=1.0_dp*PI
    mc%MAXAMP(2)=1.0_dp*mc%L0
    mc%MAXAMP(3)=1.0_dp*PI
    mc%MAXAMP(4)=1.0_dp*PI
    mc%MAXAMP(5)=1.0_dp*PI
    mc%MAXAMP(6)=0.1*mc%lbox(1)
    mc%MAXAMP(7)=NAND
    mc%MAXAMP(8)=NAND
    mc%MAXAMP(9)=NAND
    mc%MAXAMP(9)=NAND

    do MCtype=1,mc%movetypes
        mc%NADAPT(MCtype)=1000 ! adapt after at most 1000 steps
        mc%PDESIRE(MCtype)=0.5_dp ! Target
        mc%SUCCESS(MCtype)=0
        mc%PHIT(MCtype)=0.0_dp
    enddo
end subroutine

subroutine MCvar_recenter(mc,md)
!  Prevents drift in periodic BC
    IMPLICIT NONE
    type(MCvar), intent(in) :: mc
    type(MCData), intent(inout) :: md
    integer IB, I, J   ! Couners
    real(dp) R0(3)  ! Offset to move by
    IB=1
    do I=1,mc%NP
       R0(1)=nint(md%R(IB,1)/mc%lbox(1)-0.5_dp)*mc%lbox(1)
       R0(2)=nint(md%R(IB,2)/mc%lbox(2)-0.5_dp)*mc%lbox(2)
       R0(3)=nint(md%R(IB,3)/mc%lbox(3)-0.5_dp)*mc%lbox(3)
       if (abs(R0(1)*R0(2)*R0(3)) .gt. 0.0001_dp) then
           do J=1,mc%NB
              md%R(IB,1)=md%R(IB,1)-R0(1)
              md%R(IB,2)=md%R(IB,2)-R0(2)
              md%R(IB,3)=md%R(IB,3)-R0(3)
              IB=IB+1
           enddo
      endif
    enddo
end subroutine
subroutine MCvar_printEnergies(mc)
! For realtime feedback on MC simulation
    IMPLICIT NONE
    type(MCvar), intent(in) :: mc
    print*, "ECouple:", mc%ECouple
    print*, "Bending energy", mc%EELAS(1)
    print*, "Par compression energy", mc%EELAS(2)
    print*, "Shear energy", mc%EELAS(3)
    print*, "ECHI", mc%ECHI
    print*, "EField", mc%EField
    print*, "EKAP", mc%EKAP
    print*, "ebind", mc%ebind
end subroutine
subroutine MCvar_printPhi(mc,md)
! prints densities for trouble shooting
    IMPLICIT NONE
    type(MCvar), intent(in) :: mc
    type(MCData), intent(in) :: md
    integer I
    real(dp) EKap, ECouple, EChi,VV, PHIPOly
    print*,"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    print*, " PHIA  | PHIB  | PPoly |  Vol  | EKap  | EChi  |ECouple|"
    do I=1,mc%NBIN
        VV=md%Vol(I)
        if (VV.le.0.1_dp) cycle
        PHIPOLY=md%PHIA(I)+md%PHIB(I)
        EChi=VV*(mc%CHI/mc%V)*PHIPoly*(1.0_dp-PHIPoly)
        ECouple=VV*mc%HP1_Bind*(md%PHIA(I))**2
        if(PHIPoly.GT.1.0_dp) then
           EKap=VV*(mc%KAP/mc%V)*(PHIPoly-1.0_dp)**2
        else
           cycle
           EKap=0.0_dp
        endif
        write(*,"(4f8.4,3f8.1)"), md%PHIA(I), md%PHIB(I), &
                            md%PHIA(I)+md%PHIB(I),md%Vol(I),&
                            EKap,EChi,ECouple
    enddo
    print*,"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
end subroutine
subroutine MCvar_printWindowStats(mc)
! For realtime feedback on adaptation
    IMPLICIT NONE
    type(MCvar), intent(in) :: mc
    integer I ! counter
    I=0
    print*, "Succes | MCAMP | WindoW| type "
    do I=1,mc%movetypes
        if (mc%MOVEON(i).eq.1) then
            write(*,"(f8.5,2f8.2,1I8)"), mc%phit(i), mc%MCAMP(i),  mc%WindoW(i), i
        endif
    enddo
    return
end subroutine
subroutine MCvar_LoadField(mc,md,fileName)
    IMPLICIT NONE
    type(MCvar), intent(in) :: mc
    type(MCData), intent(inout) :: md
    integer I
    character*16 fileName ! file name to load from
    open (unit = 1, file = fileName, status = 'OLD')
    do I=1,mc%NBIN
        read(1,*) md%PHIH(I)
    enddo
    return
end subroutine
subroutine MCvar_MakeField(mc,md)
    IMPLICIT NONE
    type(MCvar), intent(in) :: mc
    type(MCData), intent(inout) :: md
    integer indBIN  ! index of bin
    integer IX,IY,IZ ! bin corrdinates

    do IX=1,mc%NBINX(1)
        do IY=1,mc%NBINX(2)
            do IZ=1,mc%NBINX(3)
                indBIN=IX+&
                       (IY-1)*mc%NBINX(1)+&
                       (IZ-1)*mc%NBINX(1)*mc%NBINX(2)
                md%PHIH(indBIN)=dsin(mc%k_field*mc%dbin*dble(IX))
            enddo
        enddo
    enddo
    return
end subroutine
subroutine MCvar_loadAB(mc,md,fileName)
! Loads AB for file...has not been tested
    IMPLICIT NONE
    type(MCvar), intent(in) :: mc
    type(MCData), intent(inout) :: md
    character*16, intent(in) :: fileName ! file name to load from
    integer IB, I, J ! counters
    open (unit = 1, file = fileName, status = 'OLD')
    IB=1
    do I=1,mc%NP
       do J=1,mc%NB
          read(1,"(I2)") md%AB(IB)
          IB=IB+1
          enddo
    enddo
    close(1)
end subroutine
subroutine MCvar_saveR(mc,md,fileName,repeatingBC)
! Writes R and AB to file for analysis
! Rx  Ry  Rz AB
    IMPLICIT NONE
    integer, intent(in) :: repeatingBC  ! 1 for reapeating boundary conditions
    integer I,J,IB  ! counters
    type(MCvar), intent(in) :: mc
    type(MCData), intent(in) :: md
    character*16, intent(in) :: fileName
    character*32 fullName
    fullName=  trim(fileName) // trim(mc%repSufix)
    fullName=trim(fullName)
    open (unit = 1, file = fullName, status = 'NEW')
    IB=1
    if (repeatingBC.eq.1) then
        do I=1,mc%NP
            do J=1,mc%NB
                    WRITE(1,"(3f10.3,I2)") , &
                          md%R(IB,1)-0.*nint(md%R(IB,1)/mc%lbox(1)-0.5_dp)*mc%lbox(1), &
                          md%R(IB,2)-0.*nint(md%R(IB,2)/mc%lbox(2)-0.5_dp)*mc%lbox(2), &
                          md%R(IB,3)-0.*nint(md%R(IB,3)/mc%lbox(3)-0.5_dp)*mc%lbox(3), &
                          md%AB(IB)
                IB=IB+1
            enddo
        enddo
        print*, "Error in MCvar_saveR"
        print*, "Are you sure you want reapeating BC"
        stop 1
    else
        do I=1,mc%NP
            do J=1,mc%NB
                if (mc%simtype.eq.0) then
                    WRITE(1,"(3f10.3,I2)") md%R(IB,1),md%R(IB,2),md%R(IB,3),md%AB(IB)
                else
                    WRITE(1,"(3f10.3,I2)") md%R(IB,1),md%R(IB,2),md%R(IB,3),md%AB(IB), md%METH(IB)
                endif
                IB=IB+1
            enddo
        enddo
    endif
    close(1)
end subroutine
subroutine MCVar_savePHI(mc,md,fileName)
! Saves PHIA and PHIB to file for analysis
    IMPLICIT NONE
    integer I  ! counters
    type(MCvar), intent(in) :: mc
    type(MCData), intent(in) :: md
    character*16, intent(in) :: fileName
    character*32 fullName
    fullName=  trim(fileName) // trim(mc%repSufix)
    open (unit = 1, file = fullName, status = 'NEW')
    do I=1,mc%NBIN
        WRITE(1,"(2f7.2)") md%PHIA(I),md%PHIB(I)
    enddo
    close(1)
end subroutine
subroutine MCvar_saveU(mc,md,fileName)
! Saves U to ASCII file for analisys
    IMPLICIT NONE
    integer I,J,IB  ! counters
    type(MCvar), intent(in) :: mc
    type(MCData), intent(in) :: md
    character*16, intent(in) :: fileName
    character*32 fullName
    fullName=  trim(fileName) // trim(mc%repSufix)
    open (unit = 1, file = fullName, status = 'NEW')
    IB=1
    do I=1,mc%NP
        do J=1,mc%NB
            WRITE(1,"(3f8.3,2I2)") md%U(IB,1),md%U(IB,2),md%U(IB,3)
            IB=IB+1
        enddo
    enddo
    close(1)
end subroutine
subroutine MCvar_saveparameters(mc,fileName)
! Write a number of parameters ASCII variables to file for reccords
    IMPLICIT NONE
    type(MCvar), intent(in) :: mc
    character*16, intent(in) :: fileName
    character*32 fullName
    fullName=  trim(fileName) // trim(mc%repSufix)
    open (unit =1, file = fullName, status = 'NEW')
        WRITE(1,"(I8)") mc%NT ! 1 Number of beads in simulation
        WRITE(1,"(I8)") mc%N  ! 2 Number of monomers in a polymer
        WRITE(1,"(I8)") mc%NB ! 3 Number of beads in a polymer
        WRITE(1,"(I8)") mc%NP ! 4 Number of polymers in simulation
        WRITE(1,"(I8)") mc%NT ! 5 Number of beads in simulation
        WRITE(1,"(I8)") mc%G  ! 6 Number of beads per monomer

        WRITE(1,"(f10.5)") mc%L0    ! Equilibrium segment length
        WRITE(1,"(f10.5)") mc%CHI  ! 8  initail CHI parameter value
        WRITE(1,"(f10.5)") mc%Fpoly ! Fraction polymer
        WRITE(1,"(f10.5)") mc%lbox(1)  ! 10 Lenth of box
        WRITE(1,"(f10.5)") mc%EU    ! Energy unmethalated
        WRITE(1,"(f10.5)") mc%EM    ! 12 Energy methalated
        WRITE(1,"(f10.5)") mc%HP1_Bind ! Energy of HP1 binding
        WRITE(1,"(f10.5)") (mc%L0/mc%EPS) ! 14 Khun lenth
        WRITE(1,"(A)") "-999"  ! for historic reasons
        WRITE(1,"(f10.5)") mc%F_METH  ! methalation fraction
        WRITE(1,"(f10.5)") mc%LAM_METH  ! methalation lambda
    close(1)
end subroutine
subroutine MCvar_appendEnergyData(mc,fileName)
! print Energy data
    IMPLICIT NONE
    type(MCvar), intent(in) :: mc
    LOGICAL isfile
    character*16, intent(in) :: fileName
    character*32 fullName
    fullName=  trim(fileName) // trim(mc%repSufix)
    inquire(file = fullName, exist=isfile)
    if (isfile) then
        open (unit = 1, file = fullName, status ='OLD', POSITION="append")
    else
        open (unit = 1, file = fullName, status = 'new')
        WRITE(1,*), "ind | id |",&
                    " ebend  | eparll | EShear | ECoupl | E Kap  | E Chi  |",&
                    " EField | ebind  |   M    | Couple |  Chi   |  mu    |",&
                    "  Kap   | Field  |"
    endif
    WRITE(1,"(2I5, 9f9.1,5f9.4)") mc%ind, mc%id, &
           mc%EELAS(1), mc%EELAS(2), mc%EELAS(3), mc%ECouple, &
           mc%EKap, mc%ECHI, mc%EField, mc%ebind, mc%M, &
           mc%HP1_Bind*mc%Couple_on, mc%CHI*mc%CHI_ON, mc%mu, mc%KAP*mc%KAP_ON,&
           mc%h_A
    close(1)
end subroutine
subroutine MCvar_appendAdaptData(mc,fileName)
! Appends MC move adaptation data to the file
    IMPLICIT NONE
    type(MCvar), intent(in) :: mc
    LOGICAL isfile
    character*16, intent(in) :: fileName
    character*32 fullName
    fullName=  trim(fileName) // trim(mc%repSufix)
    inquire(file = fullName, exist=isfile)
    if (isfile) then
        open (unit = 1, file = fullName, status ='OLD', POSITION="append")
    else
        open (unit = 1, file = fullName, status = 'new')
        WRITE(1,*), "ind| id|",&
                    " WIN 1 | AMP 1 | SUC 1 | WIN 2 | AMP 2 | SUC 2 |",&
                    " WIN 3 | AMP 3 | SUC 3 | ON  4 | AMP 4 | SUC 4 |",&
                    " ON  5 | AMP 5 | SUC 5 | ON  6 | AMP 6 | SUC 6 |",&
                    " ON  7 | SUC 7 | ON  8 | SUC 8 |", &
                    " ON  9 | SUC 9 | ON 10 | SUC 10|"
    endif
    WRITE(1,"(2I4,26f8.3)") mc%ind,mc%id,&
          real(mc%WindoW(1)),mc%MCAMP(1),mc%PHIT(1), &
          real(mc%WindoW(2)),mc%MCAMP(2),mc%PHIT(2), &
          real(mc%WindoW(3)),mc%MCAMP(3),mc%PHIT(3), &
          real(mc%MOVEON(4)),mc%MCAMP(4),mc%PHIT(4), &
          real(mc%MOVEON(5)),mc%MCAMP(5),mc%PHIT(5), &
          real(mc%MOVEON(6)),mc%MCAMP(6),mc%PHIT(6), &
          real(mc%MOVEON(7)),mc%PHIT(7), &
          real(mc%MOVEON(8)),mc%PHIT(8), &
          real(mc%MOVEON(9)),mc%PHIT(9), &
          real(mc%MOVEON(10)),mc%PHIT(10)
    close(1)
end subroutine
subroutine MCvar_writebinary(mc,md,baceName)
!    This function writes the contence of the structures mc and md
!  to a binary file.  if you add more variables to md you need to
!  a seperate write command for them as it is not possible to write
!  a structure with allocatables to a binar file.
!    The contence are stored in
!     baceName//'R'
!     baceName//'U'
!     etc.
    IMPLICIT NONE
    integer sizeOftype         ! for binary saving
    type(MCvar), intent(in) :: mc             ! to be save or filled
    type(MCData), intent(in) :: md             ! to be save or filled
    CHARACTER(LEN=16), intent(in) :: baceName ! for example 'record/'
    CHARACTER(LEN=16) fileName ! fileName
    CHARACTER(LEN=16) sufix    ! end of file name
    LOGICAL exists    ! does file already exist?

    !  ------parameters -----

    sizeOftype=int(SIZEOF(mc))
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
    write(1,rec=1) mc
    close(1)

    ! -------- R --------

    sizeOftype=int(SIZEOF(md%R))
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
    write(1,rec=1) md%R
    close(1)

    ! -------- U --------

    sizeOftype=int(SIZEOF(md%U))
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
    write(1,rec=1) md%U
    close(1)

    ! -------- AB --------

    sizeOftype=int(SIZEOF(md%AB))
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
    write(1,rec=1) md%AB
    close(1)

    ! -------- Vol --------

    sizeOftype=int(SIZEOF(md%Vol))
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
    write(1,rec=1) md%Vol
    close(1)
end subroutine

subroutine MCvar_readBindary(mc,md,baceName)
! This function reads what MCvar_writebinary writes and
! stores it to mc and md.  Be sure to allocate md before
! calling this command.
    IMPLICIT NONE
    integer sizeOftype         ! for binary saving
    type(MCvar) mc             ! to be save or filled
    type(MCData) md             ! to be save or filled
    CHARACTER(LEN=16) baceName ! for example 'record/'
    CHARACTER(LEN=16) fileName ! fileName
    CHARACTER(LEN=16) sufix    ! end of file name
    LOGICAL exists    ! does file already exist?

    !  ------parameters -----

    sizeOftype=int(SIZEOF(mc))
    sufix='parameters'
    fileName=trim(baceName) // trim(sufix)
    inquire(file=fileName,exist=exists)
    if(exists) then
        open(unit=1,file=fileName, status='old', &
             form='unformatted',access='direct',recl=sizeOftype)
    else
        print*, 'Error in MCvar_readBinary. file ',fileName,'does not exist'
        stop 1
    endif
    read(1,rec=1) mc
    close(1)

    ! -------- R --------

    sizeOftype=int(SIZEOF(md%R))
    sufix='R'
    fileName=trim(baceName) // trim(sufix)
    inquire(file=fileName,exist=exists)
    if(exists) then
        open(unit=1,file=fileName, status='old', &
             form='unformatted',access='direct',recl=sizeOftype)
    else
        print*, 'Error in MCvar_readBinary. file ',fileName,'does not exist'
        stop 1
    endif
    read(1,rec=1) md%R
    close(1)

    ! -------- U --------

    sizeOftype=int(SIZEOF(md%U))
    sufix='U'
    fileName=trim(baceName) // trim(sufix)
    inquire(file=fileName,exist=exists)
    if(exists) then
        open(unit=1,file=fileName, status='old', &
             form='unformatted',access='direct',recl=sizeOftype)
    else
        print*, 'Error in MCvar_readBinary. file ',fileName,'does not exist'
        stop 1
    endif
    read(1,rec=1) md%U
    close(1)

    ! -------- AB --------

    sizeOftype=int(SIZEOF(md%AB))
    sufix='AB'
    fileName=trim(baceName) // trim(sufix)
    inquire(file=fileName,exist=exists)
    if(exists) then
        open(unit=1,file=fileName, status='old', &
             form='unformatted',access='direct',recl=sizeOftype)
    else
        print*, 'Error in MCvar_readBinary. file ',fileName,'does not exist'
        stop 1
    endif
    read(1,rec=1) md%AB
    close(1)

    ! -------- Vol --------

    sizeOftype=int(SIZEOF(md%Vol))
    sufix='Vol'
    fileName=trim(baceName) // trim(sufix)
    inquire(file=fileName,exist=exists)
    if(exists) then
        open(unit=1,file=fileName, status='old', &
             form='unformatted',access='direct',recl=sizeOftype)
    else
        print*, 'Error in MCvar_readBinary. file ',fileName,'does not exist'
        stop 1
    endif
    read(1,rec=1) md%Vol
    close(1)
end subroutine
end module simMod
