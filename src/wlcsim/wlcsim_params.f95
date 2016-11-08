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

    IMPLICIT NONE

    private
    public :: Nmovetypes, dp, pi, wlcsim_params !todo, rename below

    ! hardcoded params. will need to change if certain parts of code change
    integer, parameter :: Nmovetypes = 10 ! number of MC move types

    ! precision of simulations
    integer, parameter :: dp = real64 ! preferred over SELECTED_real_Kind(15,307)

    ! universal constants
    real(dp) :: pi = 4 * atan(1.0_dp) ! fully accurate, adaptive precision

    ! for all parameters that cannot change during individual simulations
    type wlcsim_params
    !     Simulation parameters
        integer NT                ! Total number of beads  NT=NP*N*G
        integer NB                ! Number of beads in a polymer NB=N*G
        integer N                 ! Number of monomers in a polymer
        integer G                 ! Beads per monomer
        integer NP                ! Number of polymers
        real(dp) lbox(3)  ! Box length (approximate)
        real(dp) dbin      ! Discretization size (approximate)
        real(dp) L0       ! Equilibrium segment length (same as gam)
        real(dp) V        ! Bead volume
        real(dp) FA       ! Fraction of A beads
        real(dp) LAM      ! Chemical correlation parameter
        real(dp) EPS      ! Elasticity l0/(2lp)
        real(dp) CHI      ! Chi parameter value (solvent-polymer)
        real(dp) KAP      ! Incompressibility parameter
        real(dp) h_A      ! fild strength
        real(dp) k_field  ! wave vector of template field
        real(dp) Fpoly    ! Fraction Polymer
        real(dp) EU       ! Energy of binding for unmethalated
        real(dp) EM       ! Energy of binding for methalated
        real(dp) HP1_Bind ! Energy of binding of HP1 to eachother
        real(dp) F_METH   ! Fraction methalated is using option 2
        real(dp) LAM_METH ! eigenvalue of methalation setup
        real(dp) mu       ! chemical potential of HP1
        integer NBIN     ! Number of bins
        integer NBINX(3) ! Number of bin on an edge
        real(dp) eb     ! Energy of bending
        real(dp) epar   ! energy "parallel" i.e. psring energy
        real(dp) eperp  ! energy "perpendicular" i.e. shear energy
        real(dp) gam    ! average equilibrium interbead spacing
        real(dp) eta    ! bend-shear coupling parameter
        real(dp) xir    ! drag per unit persistence length
        real(dp) xiu    ! rotational drag
        real(dp) rend   ! initial end-to-end distance (if applicable in initialization)
        real(dp) del    ! L0/(lp)


    !   Monte Carlo Variables (for adaptation)
        integer movetypes
        real(dp) PDesire(Nmovetypes) ! desired hit rate
        real(dp) MAXWindoW(Nmovetypes)         ! Max Size of window for bead selection
        real(dp) MINWindoW(Nmovetypes)         ! Min Size of window for bead selection
        real(dp) MinAMP(Nmovetypes) ! minium amplitude
        real(dp) MaxAMP(Nmovetypes) ! maximum amplitude
        integer MOVEON(Nmovetypes)         ! Is the move active
        real(dp) winTarget(Nmovetypes) ! target for ratio of window to anmplitude
        integer NADAPT(Nmovetypes) ! Nunber of steps between adapt
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
        real(dp), allocatable, dimension(:):: PHIH ! Applied Field
        real(dp), allocatable, dimension(:):: Vol ! Volume fraction of A
        integer, allocatable, dimension(:):: AB            ! Chemical identity of beads
        integer, allocatable, dimension(:):: ABP           ! Test Chemical identity of beads
        integer, allocatable, dimension(:):: METH          ! Methalation state of beads
        real(dp), allocatable, dimension(:):: DPHIA    ! Change in phi A
        real(dp), allocatable, dimension(:):: DPHIB    ! Change in phi A
        integer, allocatable, dimension(:) :: indPHI   ! indices of the phi

    !   Monte Carlo Variables (for adaptation)
        real(dp) MCAMP(Nmovetypes) ! Amplitude of random change
        real(dp) WindoW(Nmovetypes)         ! Size of window for bead selection
        integer SUCCESS(Nmovetypes)        ! Number of successes
        real(dp) PHit(Nmovetypes) ! hit rate

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

subroutine MCvar_setparams(mc,fileName)
! Based on Elena's readkeys subroutine
    use INPUTparaMS, only : readLINE, readA, readF, readI, reado
    IMPLICIT NONE
    real(dp), parameter :: PI=3.141592654_dp
    type(MCvar), intent(out) :: mc
    character*16, intent(in) :: fileName  ! file with parameters
    integer :: PF   ! input file unit
    LOGICAL :: fileend=.FALSE. ! done reading file?
    CHARACTER*100 :: WORD ! keyword
    integer :: NITEMS ! number of items on the line in the parameter file
    real(dp) NAN_dp
    NAN_dp=0;NAN_dp=NAN_dp/NAN_dp

    ! ----------------------------------------------------------
    !
    !  set Default values
    !
    ! ----------------------------------------------------------

    ! file IO
    mc%FRMfile=.FALSE.
    mc%FRMMETH=.FALSE.
    mc%FRMFIELD=.false.
    mc%saveU=.FALSE.
    mc%savePhi=.FALSE.
    mc%FRMCHEM=.FALSE.
    mc%restart=.FALSE.

    ! geometry options
    mc%NP  =1
    mc%N   =2000
    mc%G   =1
    mc%NB=mc%G*mc%N
    mc%lbox(1)=25.0_dp
    mc%lbox(2)=25.0_dp
    mc%lbox(3)=25.0_dp
    mc%dbin =1.0_dp
    mc%L0  =1.25_dp
    mc%V   =0.1_dp
    mc%FA  =0.5_dp
    mc%LAM =0.0_dp
    mc%F_METH=0.5_dp
    mc%LAM_METH=0.9_dp
    mc%Fpoly=0.025_dp
    mc%k_field=1.5708_dp !0.3145_dp

    ! energy parameters
    mc%EPS =0.3_dp
    mc%CHI =0.0_dp
    mc%h_A =0.0_dp
    mc%KAP =10.0_dp
    mc%EU  =-1.52_dp
    mc%EM  =0.01_dp
    mc%mu  =0.0_dp
    mc%HP1_Bind=0.0_dp !-28.0_dp

    ! options
    mc%movetypes=10
    mc%settype = 4 ! 4 for shereical
    mc%confinetype = 3 ! 3 for sherical
    mc%simtype=1
    mc%min_accept=0.05

    ! timing options
    mc%NStep=400000
    mc%NNoInt=100
    mc%indMAX=180
    mc%reduce_move=10
    mc%useSchedule=.False.
    mc%KAP_ON=1.0_dp
    mc%CHI_ON=1.0_dp
    mc%Couple_ON=1.0_dp
    mc%N_KAP_ON=1
    mc%N_CHI_ON=1
    mc%recenter_on=.TRUE.
    mc%INITIAL_MAX_S=0.1

    ! replica options
    mc%PTON=.TRUE.
    mc%NPT=100
    mc%NRepAdapt=1000
    mc%lowerRepExe=0.09
    mc%upperRepExe=0.18
    mc%lowerCofRail=0.005
    mc%upperCofRail=0.1
    mc%indStartRepAdapt=10
    mc%indendRepAdapt=20
    mc%repAnnealSpeed=0.01
    mc%replicaBounds=.TRUE.
    mc%PT_chi =.False.
    mc%PT_h =.False.
    mc%PT_kap =.False.
    mc%PT_mu =.False.
    mc%PT_couple =.False.



    call MCvar_defaultAmp(mc)


    ! -----------------------
    !
    !  read from file
    !
    !-------------------------
    PF=55
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
           Call readI(mc%settype)
           ! settype      |  Discription
           ! _____________|_________________________________
           !    1         |   straight line in y direction with random starting
           !    2         |   rerandomize when reaching boundary, slit in z dir
           !    3         |   rerandomize when reaching boundary, cube boundary
           !    4         |   rerandomize when reaching boundary, shpere
       CASE('CONFINEtype')
           Call readI(mc%confinetype)
           ! confinetype  |  Discription
           ! _____________|_________________________________
           !    0         |  No confinement, periodic cube
           !    1         |  Between two plates in Z direction at 0 and lbox
           !    2         |  Cube of size lbox**3,  range: 0-lbox
           !    3         |  Circle of radius lbox, centered at lbox/2
           !    4         |  Periodic, unequal dimensions
       CASE('RECENTER_ON')
           Call reado(mc%recenter_on) ! recenter in periodic boundary
       CASE('SIMtype')
           Call readI(mc%simtype)
           ! simtype      | Discription
           !______________|_________________________________
           !    0         | Melt density fluctuates around fixed mean
           !    1         | Solution (For DNA)
       CASE('FRMCHEM')
           Call reado(mc%FRMCHEM) ! Initial chemical sequence from file
       CASE('FRMfile')
           call reado(mc%FRMfile) ! read configuration from file
       CASE('FRMMETH')
           Call reado(mc%FRMMETH) ! read methalation from file
       CASE('PTON')
           CALL reado(mc%PTON) ! parallel Tempering on
       CASE('SAVE_U')
           Call reado(mc%saveU)  ! save u vectors to file (every savepoint)
       CASE('SAVE_PHI')
           Call reado(mc%savePhi) ! save Phi vectors to file (every savepoint)
       CASE('L0')
           Call readF(mc%L0)  ! Equilibrium segment length
       CASE('dbin')
           Call readF(mc%dbin) ! spaitial descretation length, not tested
       CASE('lbox')
           Call readF(mc%lbox(1)) ! side length of box
           mc%lbox(2)=mc%lbox(1)
           mc%lbox(3)=mc%lbox(1)
       CASE('lboxX')
           Call readF(mc%lbox(1)) ! side length of box in x direction
       CASE('lboxY')
           Call readF(mc%lbox(2)) ! side length of box in y direction
       CASE('lboxZ')
           Call readF(mc%lbox(3)) ! side length of box in z direction
       CASE('NP')
           CALL readI(mc%NP)  ! Number of polymers
       CASE('G')
           Call readI(mc%G) ! Beads per monomer
       CASE('N')
           CALL readI(mc%N) ! Number of monomers in a polymer
       CASE('NNOINT')
           Call readI(mc%NNoInt) ! save points before turning on interaction
       CASE('N_KAP_ON')
           call readI(mc%N_KAP_ON) ! when to turn compression energy on
       CASE('N_CHI_ON')
           call readI(mc%N_CHI_ON) ! when to turn CHI energy on
       CASE('indMAX')
           Call readI(mc%indMAX) ! total number of save points
       CASE('NSTEP')
           Call readI(mc%NStep) ! steps per save point
       CASE('NPT')
           Call readI(mc%NPT) ! number of steps between parallel tempering
       CASE('FPOLY')
           Call readF(mc%Fpoly) ! Fraction Polymer
       CASE('V')
           Call readF(mc%V) ! Bead volume
       CASE('FA')
           Call readF(mc%FA) ! Fraction of A beads (fraction bound)
       CASE('LAM')
           Call readF(mc%LAM) ! Chemical correlation parameter
       CASE('EPS')
           Call readF(mc%EPS) ! Elasticity l0/(2lp)
       CASE('CHI')
           Call readF(mc%CHI) ! CHI parameter (definition depends on  hamiltoniaon
       CASE('H_A')
           Call readF(mc%h_A) ! strength of externally applied field
       CASE('KAP')
           Call readF(mc%KAP) !  Incompressibility parameter
       CASE('EU')
           Call readF(mc%EU) ! Energy of binding for unmethalated
       CASE('EM')
           Call readF(mc%EM) ! Energy of binding for methalated
       CASE('MU')
           Call readF(mc%MU) ! chemical potential of HP1
       CASE('HP1_Bind')
           Call readF(mc%HP1_Bind) ! Energy of binding of HP1 to eachother
       CASE('F_METH')
           Call readF(mc%F_METH) ! Fraction methalated
       CASE('LAM_METH')
           Call readF(mc%LAM_METH) ! eigenvalue of methalation setup
       CASE('CRANK_SHAFT_ON')
           Call readI(mc%MOVEON(1)) ! is Crank shaft move on 1/0
       CASE('SLIDE_ON')
           Call readI(mc%MOVEON(2)) ! is Slide move on 1/0
       CASE('PIVOT_ON')
           Call readI(mc%MOVEON(3)) ! is Pivot move on 1/0
       CASE('ROTATE_ON')
           Call readI(mc%MOVEON(4)) ! is single bead rotate on 1/0
       CASE('FULL_CHAIN_ROTATION_ON')
           Call readI(mc%MOVEON(5)) ! is full chain rotate on 1/0
       CASE('FULL_CHAIN_SLIDE_ON')
           Call readI(mc%MOVEON(6)) ! is full chain slide on 1/0
       CASE('Bind_MOVE_ON')
           Call readI(mc%MOVEON(7)) ! is bind/unbind move on 1/0
       CASE('CHAIN_FLIP_MOVE_ON')
           Call readI(mc%MOVEON(8)) ! is flip move move on 1/0
       CASE('CHAIN_SWAP_MOVE_ON')
           Call readI(mc%MOVEON(9)) ! is chain swap move on 1/0
       CASE('REPTATION_MOVE_ON')
           Call readI(mc%MOVEON(10)) ! is reptation move on 1/0
       CASE('MIN_CRANK_SHAFT_WIN')
           Call readF(mc%MINWindoW(1)) ! min mean window size
       CASE('MIN_SLIDE_WIN')
           Call readF(mc%MINWindoW(2))
       CASE('MIN_PIVOT_WIN')
           Call readF(mc%MINWindoW(3))
       CASE('MIN_Bind_WIN')
           Call readF(mc%MINWindoW(7))
       CASE('REDUCE_MOVE')
           Call readI(mc%reduce_move) !  only exicute unlikely movetypes every ____ cycles
       CASE('MIN_ACCEPT')
           Call readF(mc%MIN_ACCEPT) ! below which moves are turned off
       CASE('CRANK_SHAFT_TARGET')
           Call readF(mc%winTarget(1)) ! target window size for crank shaft move
       CASE('SLIDE_TARGET')
           Call readF(mc%winTarget(2)) ! target window size for slide move
       CASE('PIVOT_TARGET')
           Call readF(mc%winTarget(3)) ! target window size for Pivot move
       CASE('STRENGTH_SCHEDULE')
           Call reado(mc%useSchedule) ! use scheduled ramp in interaction strength(s)
       CASE('N_REP_ADAPT')
           Call readI(mc%NRepAdapt)  ! number of exchange attemts between adapt
       CASE('LOWER_REP_EXE')
           Call readF(mc%lowerRepExe) ! when to decrease cof spacing
       CASE('UPPER_REP_EXE')
           Call readF(mc%upperRepExe) ! when to increase cof spacing
       CASE('LOWER_COF_RAIL')
           Call readF(mc%lowerCofRail) ! minumum acceptable Cof
       CASE('UPPER_COF_RAIL')
           Call readF(mc%upperCofRail) ! maximum acceptable Cof
       CASE('ind_START_REP_ADAPT')
           Call readI(mc%indStartRepAdapt) ! ind to start rep. cof. adaptiation on
       CASE('ind_end_REP_ADAPT')
           Call readI(mc%indendRepAdapt) ! turn off rep adapt
       CASE('REP_ANNEAL_SPEED')
           Call readF(mc%repAnnealSpeed)  ! max change in cof. every adjust
       CASE('FRMFIELD')
           Call reado(mc%FRMFIELD)  ! read field from file
       CASE('K_FIELD')
           Call readF(mc%k_field)  ! wave mode for default field
       CASE('REPLICA_BOUNDS')
           Call reado(mc%replicaBounds) ! insure that 0 < s < 1
       CASE('INITIAL_MAX_S')
           call readF(mc%INITIAL_MAX_S) ! inital s of rep with highest s
       CASE('PT_CHI')
           call reado(mc%PT_chi) ! parallel temper chi
       CASE('PT_H')
           call reado(mc%PT_h) ! parallel temper h
       CASE('PT_KAP')
           call reado(mc%PT_kap) ! parallel temper kap
       CASE('PT_MU')
           call reado(mc%PT_mu) ! parallel temper mu
       CASE('PT_COUPLE')
           call reado(mc%PT_couple) ! parallel temper HP1_bind
       CASE('RESTART')
           call reado(mc%restart) ! Restart from parallel tempering
       CASE DEFAULT
           print*, "Error in MCvar_setparams.  Unidentified keyword:", &
                   TRIM(WORD)
           stop 1
       endSELECT
    enddo
    close(PF)

    if ((mc%NBINX(1)-mc%NBINX(2).ne.0).or. &
        (mc%NBINX(1)-mc%NBINX(3).ne.0)) then
        if (mc%simtype.eq.1) then
            print*, "Solution not tested with non-cube box, more coding needed"
            stop 1
        endif
        if (mc%confinetype.ne.4) then
            print*, "Unequal boundaries require confinetype=4"
            stop 1
        endif
        if (mc%settype.eq.4) then
            print*, "You shouldn't put a shpere in and unequal box!"
            stop 1
        endif
    endif
    ! --------------------
    !
    ! Derived Variables, Reconcile inputs
    !
    ! --------------------
    if (mc%simtype.eq.1) then
        mc%NT=mc%N*mc%NP*mc%G
        mc%NB=mc%N*mc%G
        if (mc%confinetype.eq.3) then
            mc%lbox(1)=(mc%V*mc%NT*6/(mc%Fpoly*PI))**(1.0_dp/3.0_dp)
            mc%lbox(2)=mc%lbox(1)
            mc%lbox(3)=mc%lbox(1)
        else
            mc%lbox(1)=(mc%V*mc%NT/mc%Fpoly)**(1.0_dp/3.0_dp)
            mc%lbox(2)=mc%lbox(1)
            mc%lbox(3)=mc%lbox(1)
        endif
        mc%NBINX(1)=nint(mc%lbox(1)/mc%dbin)
        mc%NBINX(2)=nint(mc%lbox(2)/mc%dbin)
        mc%NBINX(3)=nint(mc%lbox(3)/mc%dbin)
        mc%NBIN=mc%NBINX(1)*mc%NBINX(2)*mc%NBINX(3)
        mc%lbox(1) = mc%NBINX(1)*mc%dbin! used to be: dbin=lbox/NBINX
        mc%lbox(2) = mc%NBINX(2)*mc%dbin! used to be: dbin=lbox/NBINX
        mc%lbox(3) = mc%NBINX(3)*mc%dbin! used to be: dbin=lbox/NBINX
    elseif (mc%simtype.eq.0) then
        if (mc%confinetype.eq.0) then
            mc%NP=nint(mc%lbox(1)*mc%lbox(2)*mc%lbox(3)/(mc%N*mc%G*mc%V))
            mc%lbox=(mc%V*mc%N*mc%G*mc%NP)**(1.0_dp/3.0_dp)
            mc%NBINX(1)=nint(mc%lbox(1)/mc%dbin)
            mc%NBINX(2)=nint(mc%lbox(2)/mc%dbin)
            mc%NBINX(3)=nint(mc%lbox(3)/mc%dbin)
            mc%dbin=mc%lbox(1)/mc%NBINX(1)
        elseif(mc%confinetype.eq.4) then
            mc%dbin=mc%lbox(1)/nint(mc%lbox(1)/mc%dbin)
            mc%NBINX(1)=nint(mc%lbox(1)/mc%dbin)
            mc%NBINX(2)=nint(mc%lbox(2)/mc%dbin)
            mc%NBINX(3)=nint(mc%lbox(3)/mc%dbin)
            mc%lbox(2)=mc%dbin*mc%NBINX(2)
            mc%lbox(3)=mc%dbin*mc%NBINX(3)
            mc%NP=nint(mc%lbox(1)*mc%lbox(2)*mc%lbox(3)/(mc%N*mc%G*mc%V))
            print*, "Density =", &
                  mc%N*mc%G*mc%V*mc%NP/(mc%lbox(1)*mc%lbox(2)*mc%lbox(3))
        endif
        mc%NB=mc%N*mc%G
        mc%NBIN=mc%NBINX(1)*mc%NBINX(2)*mc%NBINX(3)
        mc%NT=mc%N*mc%NP*mc%G
        mc%WindoW(5)=mc%NB
        mc%WindoW(6)=mc%NB
        mc%WindoW(8)=mc%NB
        mc%WindoW(9)=mc%NB
        mc%MCAMP(2)=0.3_dp*mc%L0
        mc%MCAMP(6)=5.0_dp*mc%L0
        mc%MINAMP(2)=0.2_dp*mc%L0
        mc%MINAMP(6)=0.2_dp*mc%L0
        mc%MAXWindoW(1)=min(150,mc%NB)
        mc%MAXWindoW(2)=min(150,mc%NB)
        mc%MAXWindoW(3)=min(150,mc%NB)
        mc%MAXWindoW(7)=min(150,mc%NB)
        mc%MAXAMP(2)=1.0_dp*mc%L0
        mc%MAXAMP(6)=0.1*mc%lbox(1)
    else
       print*, "Error in simMod: symtype",mc%simtype," not found"
    endif
    call getpara(mc%para,mc%EPS,mc%L0,NAN_dp)

    ! -----------------------
    !
    ! Set initial values
    !
    ! ----------------------
    mc%EElas(1)=0.0_dp
    mc%EElas(2)=0.0_dp
    mc%EElas(3)=0.0_dp
    mc%ECouple=0.0_dp
    mc%ebind=0.0_dp
    mc%EKap=0.0_dp
    mc%ECHI=0.0_dp
    mc%EField=0.0_dp
    mc%x_mu=0.0_dp
    mc%x_Field=0.0_dp
    mc%x_couple=0.0_dp
    mc%x_Kap=0.0_dp
    mc%x_Chi=0.0_dp

    !-----------------------------
    !
    !  Idiot checks
    !
    !-----------------------------
    if (mc%NBINX(1)*mc%NBINX(2)*mc%NBINX(3).ne.mc%NBIN) then
        print*, "error in MCsim. Wrong number of bins"
        stop 1
    endif
    if (mc%NT.ne.mc%N*mc%NP*mc%G) then
        print*, "error in MCsim.  NT=",mc%NT," N=",mc%N," NP=",mc%NP," G=",mc%G
        stop 1
    endif
    if (mc%NB.ne.mc%N*mc%G) then
        print*, "error in MCsim.  NB=",mc%NB," N=",mc%N," G=",mc%G
        stop 1
    endif
    if (mc%NNoInt.gt.mc%indStartRepAdapt) then
        print*, "error in MCsim. don't run adapt without int on"
        stop 1
    endif
    if (mc%NNoInt.gt.mc%N_CHI_ON) then
        print*, "error in MCsim. Can't have chi without int on"
        stop 1
    endif
    if (mc%NNoInt.gt.mc%N_KAP_ON) then
        print*, "error in MCsim. Can't have kap without int on"
        stop 1
    endif

end subroutine
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
subroutine MCvar_defaultAmp(mc)
    IMPLICIT NONE
    real(dp), parameter :: PI=3.141592654_dp
    type(MCvar), intent(inout) :: mc
    integer MCtype ! type of move
    integer NANI ! NaN
    real(dp) NAND !NAND
    NANI=0;NANI=NANI/NANI
    NAND=0.0_dp;NAND=NAND/NAND
!   ~~~~~~~~~~~~~~~~~~~
!    Edit the following variables for better performance
!   ~~~~~~~~~~~~~~~~~~~
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
