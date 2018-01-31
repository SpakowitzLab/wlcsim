#include "../defines.inc"
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
    use binning, only: constructBin, binType, addBead

    implicit none

    public

    !!!     hardcoded params. will need to change if certain parts of code change
    ! number of wlc_p move types
    integer, parameter :: nMoveTypes = 11
    integer, parameter :: nDim = 3

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
    real(dp), parameter :: HUGE_ENERGY = 9990000.0_dp

    character(len = 20), parameter, dimension(nMoveTypes) :: moveNames(nMoveTypes) = &
        (/'crank-shaft         ','slide               ',&
          'pivot               ','rotate              ',&
          'fullChainRotation   ','fullChianSlide      ',&
          'chem-identity       ','end-end filp        ',&
          'chain swap          ','reptation           ',&
          'superReptation      '/)

    !!!     universal constants
    ! fully accurate, adaptive precision
    real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
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
    real(dp) :: max_wlc_l0 = 0.01_dp
    real(dp) :: max_sswlc_delta = 10.0_dp
    integer, parameter :: INT_MIN = -HUGE(nMoveTypes)

    ! for all parameters that cannot change during individual simulations
    ! these are documented more thoroughly where they are read in (see the
    ! subroutine get_input_from_file), in the docs (TOdo), and often the default
    ! values will help with understanding how the variable is used.
    !
    ! many of these variables are used only in certain kinds of simulations
    type wlcsim_params
    !   Simulation parameters
        integer simType           ! whether to use WLC, ssWLC, or Gaussian Chain
        integer nT                ! Total number of beads  NT = nBpM*nMpP*np
        real(dp) dt ! sets time scale of simulation
        real(dp) l0       ! Path length between beads. (meaning unknown for gaussian chain?)
        real(dp) gam    ! average equilibrium interbead spacing
        real(dp) eta    ! bend-shear coupling parameter
        real(dp) xir    ! drag per unit persistence length
        real(dp) sigma  ! variance of interbead position distribution of appropriately renormalized gaussian chain
        real(dp) xiu    ! rotational drag
        real(dp) eps      ! number of kuhn lengths between beads
        real(dp) del      ! number of persistence lengths between beads
        real(dp) chi      ! Chi parameter value (solvent-polymer) (Flory-Huggins separation constant (how much A/B's hate each))
        real(dp) chi_l2   ! maier saupe parameter (possibly multiplied by 4pi or something like that)
        real(dp) kap      ! Incompressibility parameter of the melt
        real(dp) lhc    !TOdo something to do with intrapolymer interaction strength
        real(dp) vhc    !TOdo something to do with intrapolymer interaction strength, fill in defaults, etc
        real(dp) eb     ! effective bending energy for ssWLC
        real(dp) eperp  ! effective shearing energy for ssWLC
        real(dp) epar   ! effective stretch energy for ssWLC

    !   for passing 1st order phase transition in (quinn/shifan's) random copolymer wlc_p sims
        real(dp) hA       ! strength of applied sinusoidal field (used in PT to step around 1st order phase transition)
        real(dp) rend   ! initial end-to-end distance (if applicable in initialization)

    !   for simulating chromatin methylation
        real(dp) HP1_Bind ! Energy of binding of HP1 to eachother
        real(dp) mu       ! chemical potential of HP1

    !   boundary/box things
        integer NBin     ! Number of bins
        integer NBinX(nDim) ! Number of MC bins on an edge

    !   Monte Carlo Variables (for adaptation)
        real(dp) PDesire(nMoveTypes) ! desired hit rate
        real(dp) MAXWindoW(nMoveTypes)         ! Max Size of window for bead selection
        real(dp) MinWindoW(nMoveTypes)         ! Min Size of window for bead selection
        real(dp) MinAMP(nMoveTypes) ! minium amplitude
        real(dp) MaxAMP(nMoveTypes) ! maximum amplitude
        integer MOVEON(nMoveTypes)         ! Is the move active
        real(dp) winTarget(nMoveTypes) ! target for ratio of window to anmplitude
        integer NADAPT(nMoveTypes) ! Nunber of steps between adapt

    !   Timing variables
        integer movesPerStep(nMoveTypes) ! how many times each move is done per step

    !   Switches
        integer LK                ! Linking number
        real(dp) KAP_ON     ! fraction of KAP energy contributing to "calculated" energy
        real(dp) CHI_ON     ! fraction of CHI energy contributing to "calculated" energy
        real(dp) Couple_ON  ! fraction of Coupling energy contributing to "calculated" energy
        logical field_int_on ! include field interactions (e.g. A/B interactions) uses many of the same functions as the chemical identity/"meth"ylation code, but energies are calcualted via a field-based approach
        logical chi_l2_on

    !   parallel Tempering parameters

    !   Replica Dynamic Cof choice

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
        real(dp), allocatable, dimension(:,:):: PHI_l2 ! l=2 oreientational field
        real(dp), allocatable, dimension(:):: PHIH ! Quinn's sinusoidal field for passing 1st order phase transitions
        real(dp), allocatable, dimension(:):: Vol  ! Volume fraction of A
        integer, allocatable, dimension(:):: AB    ! Chemical identity of beads
        integer, allocatable, dimension(:):: ABP   ! Test Chemical identity of beads
        integer, allocatable, dimension(:):: METH  ! Methalation state of beads
        real(dp), allocatable, dimension(:):: DPHIA    ! Change in phi A
        real(dp), allocatable, dimension(:):: DPHIB    ! Change in phi A
        real(dp), allocatable, dimension(:,:):: DPHI_l2 ! change in l=2 oreientational field
        integer, allocatable, dimension(:) :: indPHI   ! indices of the phi
        ! simulation times at which (i,j)th bead pair first collided
        real(dp), allocatable, dimension(:,:) :: coltimes
        real(dp) :: wr

        type(binType) bin ! Structure for keeping track of neighbors

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
        integer ATTEMPTS(nMoveTypes)        ! Number of successes
        integer successTOTAL(nMoveTypes)               !Total number of successes
        real(dp) PHit(nMoveTypes) ! hit rate

    !   Energys
        real(dp) eElas(4) ! Elastic energy
        real(dp) eChi     ! CHI energy
        real(dp) eKap     ! KAP energy
        real(dp) eCouple  ! Coupling
        real(dp) eBind    ! binding energy
        real(dp) eMu      ! Chemical potential energy
        real(dp) eField   ! Field energy
        real(dp) eSelf    ! repulsive lennard jones on closest approach self-interaction energy (polymer on polymer)
        real(dp) eMaierSaupe ! Maier Saupe energy

    !   Congigate Energy variables (needed to avoid NaN when cof-> 0 in rep exchange)
        real(dp) x_Chi,   dx_Chi
        real(dp) x_Couple,dx_Couple
        real(dp) x_Kap,   dx_Kap
        real(dp) x_Field, dx_Field
        real(dp) x_Mu,    dx_Mu
        real(dp) x_maierSaupe, dx_maierSaupe ! Maier Saupe energy / chi_l2

    !   Move Variables
        real(dp) DEELAS(4) ! Change in bending energy
        real(dp) DECouple ! Coupling energy
        real(dp) DEChi    ! chi interaction energy
        real(dp) DEKap    ! compression energy
        real(dp) Debind   ! Change in binding energy
        real(dp) DEMu   ! Change in binding energy
        real(dp) DEField  ! Change in field energy
        real(dp) DESelf   ! change in self interaction energy
        real(dp) ECon     ! Confinement Energy
        real(dp) deMaierSaupe ! change in Maier Saupe energy
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
        integer ind_exchange            ! number of exchange moves since last save point
        integer time_ind                ! current time point
        real(dp) time
    end type


contains

    subroutine set_param_defaults(wlc_p)
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

        wlc_p%L0 = WLC_P__L/real(WLC_P__NB-1.0_dp) ! -1.0 because one fewer segments then beads
        wlc_p%EPS=wlc_p%L0/(2.0_dp*WLC_P__LP)

        ! parallel temper variables
        wlc_p%CHI      = WLC_P__CHI
        wlc_p%MU       = WLC_p__MU
        wlc_p%HA       = WLC_P__HA
        wlc_p%HP1_BIND = WLC_P__HP1_BIND
        wlc_p%KAP      = WLC_P__KAP
        wlc_p%CHI_L2   = WLC_P__CHI_L2

        wlc_p%lhc = NAN ! I have no idea what this does
        wlc_p%vhc = NAN ! I have no idea what this does
        wlc_p%couple_on = 1.0 ! on by default
        wlc_p%kap_on = 1.0 ! on by default
        wlc_p%chi_on = 1.0 ! on by default
        wlc_p%chi_l2_on = .TRUE. ! on by default
        wlc_p%field_int_on = WLC_P__FIELD_INT_ON ! on by default

        wlc_p%NBINX(1) = WLC_P__NBINX_X
        wlc_p%NBINX(2) = WLC_P__NBINX_Y
        wlc_p%NBINX(3) = WLC_P__NBINX_Z
        wlc_p%PDESIRE(1) = WLC_P__PDESIRE_CRANK_SHAFT
        wlc_p%PDESIRE(2) = WLC_P__PDESIRE_SLIDE_MOVE
        wlc_p%PDESIRE(3) = WLC_P__PDESIRE_PIVOT_MOVE
        wlc_p%PDESIRE(4) = WLC_P__PDESIRE_ROTATE_MOVE
        wlc_p%PDESIRE(5) = WLC_P__PDESIRE_FULL_CHAIN_ROTATION
        wlc_p%PDESIRE(6) = WLC_P__PDESIRE_FULL_CHAIN_SLIDE
        wlc_p%PDESIRE(7) = WLC_P__PDESIRE_CHANGE_BINDING_STATE
        wlc_p%PDESIRE(8) = WLC_P__PDESIRE_CHAIN_FLIP
        wlc_p%PDESIRE(9) = WLC_P__PDESIRE_CHAIN_EXCHANGE
        wlc_p%PDESIRE(10) = WLC_P__PDESIRE_REPTATION
        wlc_p%PDESIRE(11) = WLC_P__PDESIRE_SUPER_REPTATION
        wlc_p%MAXWINDOW(1) = WLC_P__MAXWINDOW_CRANK_SHAFT
        wlc_p%MAXWINDOW(2) = WLC_P__MAXWINDOW_SLIDE_MOVE
        wlc_p%MAXWINDOW(3) = WLC_P__MAXWINDOW_PIVOT_MOVE
        wlc_p%MAXWINDOW(4) = WLC_P__MAXWINDOW_ROTATE_MOVE
        wlc_p%MAXWINDOW(5) = WLC_P__MAXWINDOW_FULL_CHAIN_ROTATION
        wlc_p%MAXWINDOW(6) = WLC_P__MAXWINDOW_FULL_CHAIN_SLIDE
        wlc_p%MAXWINDOW(7) = WLC_P__MAXWINDOW_CHANGE_BINDING_STATE
        wlc_p%MAXWINDOW(8) = WLC_P__MAXWINDOW_CHAIN_FLIP
        wlc_p%MAXWINDOW(9) = WLC_P__MAXWINDOW_CHAIN_EXCHANGE
        wlc_p%MAXWINDOW(10) = WLC_P__MAXWINDOW_REPTATION
        wlc_p%MAXWINDOW(11) = WLC_P__MAXWINDOW_SUPER_REPTATION
        wlc_p%MINWINDOW(1) = WLC_P__MINWINDOW_CRANK_SHAFT
        wlc_p%MINWINDOW(2) = WLC_P__MINWINDOW_SLIDE_MOVE
        wlc_p%MINWINDOW(3) = WLC_P__MINWINDOW_PIVOT_MOVE
        wlc_p%MINWINDOW(4) = WLC_P__MINWINDOW_ROTATE_MOVE
        wlc_p%MINWINDOW(5) = WLC_P__MINWINDOW_FULL_CHAIN_ROTATION
        wlc_p%MINWINDOW(6) = WLC_P__MINWINDOW_FULL_CHAIN_SLIDE
        wlc_p%MINWINDOW(7) = WLC_P__MINWINDOW_CHANGE_BINDING_STATE
        wlc_p%MINWINDOW(8) = WLC_P__MINWINDOW_CHAIN_FLIP
        wlc_p%MINWINDOW(9) = WLC_P__MINWINDOW_CHAIN_EXCHANGE
        wlc_p%MINWINDOW(10) = WLC_P__MINWINDOW_REPTATION
        wlc_p%MINWINDOW(11) = WLC_P__MINWINDOW_SUPER_REPTATION
        wlc_p%MINAMP(1) = WLC_P__MINAMP_CRANK_SHAFT
        wlc_p%MINAMP(2) = WLC_P__MINAMP_SLIDE_MOVE
        wlc_p%MINAMP(3) = WLC_P__MINAMP_PIVOT_MOVE
        wlc_p%MINAMP(4) = WLC_P__MINAMP_ROTATE_MOVE
        wlc_p%MINAMP(5) = WLC_P__MINAMP_FULL_CHAIN_ROTATION
        wlc_p%MINAMP(6) = WLC_P__MINAMP_FULL_CHAIN_SLIDE
        wlc_p%MINAMP(7) = WLC_P__MINAMP_CHANGE_BINDING_STATE
        wlc_p%MINAMP(8) = WLC_P__MINAMP_CHAIN_FLIP
        wlc_p%MINAMP(9) = WLC_P__MINAMP_CHAIN_EXCHANGE
        wlc_p%MINAMP(10) = WLC_P__MINAMP_REPTATION
        wlc_p%MINAMP(11) = WLC_P__MINAMP_SUPER_REPTATION
        wlc_p%MAXAMP(1) = WLC_P__MAXAMP_CRANK_SHAFT
        wlc_p%MAXAMP(2) = WLC_P__MAXAMP_SLIDE_MOVE
        wlc_p%MAXAMP(3) = WLC_P__MAXAMP_PIVOT_MOVE
        wlc_p%MAXAMP(4) = WLC_P__MAXAMP_ROTATE_MOVE
        wlc_p%MAXAMP(5) = WLC_P__MAXAMP_FULL_CHAIN_ROTATION
        wlc_p%MAXAMP(6) = WLC_P__MAXAMP_FULL_CHAIN_SLIDE
        wlc_p%MAXAMP(7) = WLC_P__MAXAMP_CHANGE_BINDING_STATE
        wlc_p%MAXAMP(8) = WLC_P__MAXAMP_CHAIN_FLIP
        wlc_p%MAXAMP(9) = WLC_P__MAXAMP_CHAIN_EXCHANGE
        wlc_p%MAXAMP(10) = WLC_P__MAXAMP_REPTATION
        wlc_p%MAXAMP(11) = WLC_P__MAXAMP_SUPER_REPTATION
        wlc_p%MOVEON(1) = WLC_P__MOVEON_CRANK_SHAFT
        wlc_p%MOVEON(2) = WLC_P__MOVEON_SLIDE_MOVE
        wlc_p%MOVEON(3) = WLC_P__MOVEON_PIVOT_MOVE
        wlc_p%MOVEON(4) = WLC_P__MOVEON_ROTATE_MOVE
        wlc_p%MOVEON(5) = WLC_P__MOVEON_FULL_CHAIN_ROTATION
        wlc_p%MOVEON(6) = WLC_P__MOVEON_FULL_CHAIN_SLIDE
        wlc_p%MOVEON(7) = WLC_P__MOVEON_CHANGE_BINDING_STATE
        wlc_p%MOVEON(8) = WLC_P__MOVEON_CHAIN_FLIP
        wlc_p%MOVEON(9) = WLC_P__MOVEON_CHAIN_EXCHANGE
        wlc_p%MOVEON(10) = WLC_P__MOVEON_REPTATION
        wlc_p%MOVEON(11) = WLC_P__MOVEON_SUPER_REPTATION
        wlc_p%WINTARGET(1) = WLC_P__WINTARGET_CRANK_SHAFT
        wlc_p%WINTARGET(2) = WLC_P__WINTARGET_SLIDE_MOVE
        wlc_p%WINTARGET(3) = WLC_P__WINTARGET_PIVOT_MOVE
        wlc_p%WINTARGET(4) = WLC_P__WINTARGET_ROTATE_MOVE
        wlc_p%WINTARGET(5) = WLC_P__WINTARGET_FULL_CHAIN_ROTATION
        wlc_p%WINTARGET(6) = WLC_P__WINTARGET_FULL_CHAIN_SLIDE
        wlc_p%WINTARGET(7) = WLC_P__WINTARGET_CHANGE_BINDING_STATE
        wlc_p%WINTARGET(8) = WLC_P__WINTARGET_CHAIN_FLIP
        wlc_p%WINTARGET(9) = WLC_P__WINTARGET_CHAIN_EXCHANGE
        wlc_p%WINTARGET(10) = WLC_P__WINTARGET_REPTATION
        wlc_p%WINTARGET(11) = WLC_P__WINTARGET_SUPER_REPTATION
        wlc_p%NADAPT(1) = WLC_P__NADAPT_CRANK_SHAFT
        wlc_p%NADAPT(2) = WLC_P__NADAPT_SLIDE_MOVE
        wlc_p%NADAPT(3) = WLC_P__NADAPT_PIVOT_MOVE
        wlc_p%NADAPT(4) = WLC_P__NADAPT_ROTATE_MOVE
        wlc_p%NADAPT(5) = WLC_P__NADAPT_FULL_CHAIN_ROTATION
        wlc_p%NADAPT(6) = WLC_P__NADAPT_FULL_CHAIN_SLIDE
        wlc_p%NADAPT(7) = WLC_P__NADAPT_CHANGE_BINDING_STATE
        wlc_p%NADAPT(8) = WLC_P__NADAPT_CHAIN_FLIP
        wlc_p%NADAPT(9) = WLC_P__NADAPT_CHAIN_EXCHANGE
        wlc_p%NADAPT(10) = WLC_P__NADAPT_REPTATION
        wlc_p%NADAPT(11) = WLC_P__NADAPT_SUPER_REPTATION
        wlc_p%MOVESPERSTEP(1) = WLC_P__MOVESPERSTEP_CRANK_SHAFT
        wlc_p%MOVESPERSTEP(2) = WLC_P__MOVESPERSTEP_SLIDE_MOVE
        wlc_p%MOVESPERSTEP(3) = WLC_P__MOVESPERSTEP_PIVOT_MOVE
        wlc_p%MOVESPERSTEP(4) = WLC_P__MOVESPERSTEP_ROTATE_MOVE
        wlc_p%MOVESPERSTEP(5) = WLC_P__MOVESPERSTEP_FULL_CHAIN_ROTATION
        wlc_p%MOVESPERSTEP(6) = WLC_P__MOVESPERSTEP_FULL_CHAIN_SLIDE
        wlc_p%MOVESPERSTEP(7) = WLC_P__MOVESPERSTEP_CHANGE_BINDING_STATE
        wlc_p%MOVESPERSTEP(8) = WLC_P__MOVESPERSTEP_CHAIN_FLIP
        wlc_p%MOVESPERSTEP(9) = WLC_P__MOVESPERSTEP_CHAIN_EXCHANGE
        wlc_p%MOVESPERSTEP(10) = WLC_P__MOVESPERSTEP_REPTATION
        wlc_p%MOVESPERSTEP(11) = WLC_P__MOVESPERSTEP_SUPER_REPTATION

    end subroutine set_param_defaults





    subroutine idiot_checks(wlc_p, wlc_d)
#if MPI_VERSION
        use mpi
#endif
        implicit none
        type(wlcsim_params), intent(inout) :: wlc_p
        type(wlcsim_data), intent(inout) :: wlc_d
        logical err

        if (WLC_P__ASYMMETRICALTERNATINGCHEM .and. WLC_P__CHANGINGCHEMICALIDENTITY) then
            print*, "Asymmetric AlternatingChem and changing Chemical Identity is not avaiable."
            stop
        endif
        if (WLC_P__RING) then
            if (WLC_P__NP .gt. 1) then
                print*, "As of the writing of this error message"
                print*, "MC_eelals and possible energy_elas are"
                print*, "not capable of more than one rings"
                stop
            endif
            if (WLC_P__INITCONDTYPE == 'randomWalkWithBoundary') then
                print*, "initCondType = 7 doesn't know how to make a ring."
                stop
            endif
        endif

        if (WLC_P__LBOX_X .ne. WLC_P__LBOX_X) then
            print*, "No box size set.  If you need a box please specify it."
            call stop_if_err(WLC_P__INITCONDTYPE /= 'randomWalkWithBoundary', &
                'Only one initial polymer config supported if you''re not '//&
                'using LBOX to define a MC simulation box.')
        else
            if ((wlc_p%NBIN > 8000000).or.(wlc_p%NBIN.lt.1)) then
                print*, "ERROR: Requested ", wlc_p%NBIN," bins."
                print*, "You probably don't want this."
                print*, "Comment me out if you do."
                stop 1
            endif
        endif

        if (wlc_p%FIELD_INT_ON .and. (WLC_P__LBOX_X .ne. WLC_P__LBOX_Y .or. WLC_P__LBOX_Y .ne. WLC_P__LBOX_Z)) then
            call stop_if_err(.True., 'Bin-based fields not tested with non-cube boundary box size.')
        endif

        call stop_if_err(WLC_P__COLLISIONDETECTIONTYPE == 2, &
            'KD-tree based collision detection not yet implemented.')

        call stop_if_err(wlc_p%REND > WLC_P__L, &
            "Requesting initial end-to-end distance larger than polymer length.")

        if (WLC_P__CODENAME == 'quinn') then
           if ((wlc_p%NBINX(1)-wlc_p%NBINX(2).ne.0).or. &
                (wlc_p%NBINX(1)-wlc_p%NBINX(3).ne.0)) then
              err = WLC_P__CONFINETYPE.ne.'periodicUnequal'
              call stop_if_err(err, "Unequal boundaries require confinetype = periodicUnequal")
              err = WLC_P__INITCONDTYPE.eq.'randomLineSphereBoundary'
              call stop_if_err(err, "You shouldn't put a sphere in and unequal box!")
           endif

           err = wlc_p%NBINX(1)*wlc_p%NBINX(2)*wlc_p%NBINX(3).ne.wlc_p%NBIN
           call stop_if_err(err, "error in mcsim. Wrong number of bins")

           !TOdo: replace with semantic descriptions of error encountered, instead
           ! of simply outputting the input that the user put in
           if (wlc_p%NT.ne.WLC_P__NMPP*WLC_P__NP*WLC_P__NBPM) then
              print*, "error in mcsim.  NT = ",wlc_p%NT," nMpP = ",WLC_P__NMPP," NP = ",WLC_P__NP," nBpM = ",WLC_P__NBPM
              stop 1
           endif

           if (WLC_P__NB.ne.WLC_P__NMPP*WLC_P__NBPM) then
              print*, "error in mcsim.  NB = ",WLC_P__NB," nMpP = ",WLC_P__NMPP," nBpM = ",WLC_P__NBPM
              stop 1
           endif

           err = WLC_P__NNOINT.gt.WLC_P__INDSTARTREPADAPT
           call stop_if_err(err, "error in mcsim. don't run adapt without int on")

           err = WLC_P__NNOINT.gt.WLC_P__N_CHI_ON
           call stop_if_err(err, "error in mcsim. Can't have chi without int on")

           err = WLC_P__NNOINT.gt.WLC_P__N_CHI_L2_ON
           call stop_if_err(err, "error in mcsim. Can't have chi_l2 without int on")

           err = WLC_P__NNOINT.gt.WLC_P__N_KAP_ON
           call stop_if_err(err, "error in mcsim. Can't have kap without int on")

           err = (wlc_p%MOVEON(7) /= 0 .and. (.not. WLC_P__VARIABLE_CHEM_STATE))
           call stop_if_err(err,"You need bindon if you have bindmove on")

        endif


#if MPI_VERSION
    if (WLC_P__PT_TWIST) then
        if (.NOT.WLC_P__TWIST) then
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


    subroutine set_parameters(wlc_d, wlc_p)
        ! Based on Elena's readkeys subroutine
        implicit none
        type(wlcsim_params), intent(out) :: wlc_p
        type(wlcsim_data), intent(out) :: wlc_d

        ! baseline defaults
        call set_param_defaults(wlc_p)

        ! advanced defaults that require some inputs to specify
        call tweak_param_defaults(wlc_p, wlc_d)

        ! get derived parameters that aren't directly input from file
        call get_renormalized_chain_params(wlc_p)

        !If parallel tempering is on, read the Lks
        if (WLC_P__PT_TWIST) then
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
#if MPI_VERSION
        integer ( kind = 4 ) dest   !destination id for messages
        integer ( kind = 4 ) source  !source id for messages
        integer ( kind = 4 ) status(MPI_status_SIZE) ! MPI stuff
        integer ( kind = 4 ) error  ! error id for MIP functions
#endif
        character(MAXFILENAMELEN) iostr  ! string for file name
        real(dp) setBinSize(3)
        real(dp) setMinXYZ(3) ! location of corner of bin
        integer setBinShape(3)! Specify first level of binning
        nt = wlc_p%NT
        nbin = wlc_p%NBIN

#if MPI_VERSION
        call init_MPI(wlc_d)
#endif
        allocate(wlc_d%R(3,NT))
        allocate(wlc_d%U(3,NT))
        if (WLC_P__CODENAME /= 'bruno' .OR. WLC_P__NINITMCSTEPS /= 0) then
            allocate(wlc_d%RP(3,NT))
            allocate(wlc_d%UP(3,NT))
            wlc_d%RP=nan  ! To prevent accidental use
            wlc_d%UP=nan
        endif
        !TOdo these should in principle be inside the following if statement,
        !but it's not clear if that's possible without adding a bunch of dirty
        !if statements deep inside mc_move. which is fine, but I would want to
        !check with quinn *exactly* in which cases they're needed if i do that
        if (wlc_p%FIELD_INT_ON) then
            allocate(wlc_d%AB(NT))   !Chemical identity aka binding state
            if (WLC_P__CHANGINGCHEMICALIDENTITY) then
                allocate(wlc_d%ABP(NT))   !Chemical identity aka binding state
                wlc_d%ABP = INT_MIN
            endif
            if (wlc_p%CHI_L2_ON) then
                allocate(wlc_d%PHI_l2(-2:2,NBin))
                allocate(wlc_d%dPHI_l2(-2:2,NBin))
            endif
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
                wlc_d%indphi(I) = INT_MIN
            enddo
        endif
        if (WLC_P__VARIABLE_CHEM_STATE) then
            allocate(wlc_d%METH(NT)) !Underlying methalation profile
        endif
        !Allocate vector of writhe and elastic energies for replicas
        if (WLC_P__PT_TWIST) then
            allocate(wlc_d%Wrs(wlc_d%nLKs))
            allocate(wlc_d%eelasREPLICAS(wlc_d%nLKs,4))

        endif
        if (WLC_P__RING) then !TOdo this should be if ("knot")
            wlc_d%NCross = 0
            wlc_d%CrossSize = min(10000,WLC_P__NB)**2
            if (wlc_d%CrossSize == 10000) then
                print*, "Two many beas for ring calculation"
                stop 1
            endif
            ! In include
            allocate(wlc_d%Cross(wlc_d%CrossSize,6))
            allocate(wlc_d%CrossP(wlc_d%CrossSize,6))
        endif
        !If parallel tempering is on, initialize the nodeNumbers

        if (WLC_P__PT_TWIST) then

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

        if (WLC_P__COLLISIONDETECTIONTYPE /= 0) then
            allocate(wlc_d%coltimes(NT,NT))
            wlc_d%coltimes = -1.0_dp
        else
            allocate(wlc_d%coltimes(1,1))
            wlc_d%coltimes = nan
        endif

#if MPI_VERSION
        ! -----------------------------------------------
        !
        !   Generate thread safe random number seeds
        !
        !--------------------------------------------
        call MPI_COMM_RANK(MPI_COMM_WORLD,wlc_d%id,error)
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
        call initcond(wlc_d%R, wlc_d%U, wlc_p%NT, WLC_P__NB, &
            WLC_P__NP, WLC_P__FRMFILE, pack_as_para(wlc_p), &
            wlc_d%rand_stat, wlc_p)

        if (wlc_p%FIELD_INT_ON) then
            ! initialize a/b sequence
            if (WLC_P__CHEM_STATE_FROM_FILE) then
                iostr='input/ab'
                call MCvar_loadAB(wlc_d,iostr)
            elseif (WLC_P__ASYMMETRICALTERNATINGCHEM) then
                call alternChem(wlc_d%AB, wlc_p%NT, WLC_P__NMPP, WLC_P__NBPM, WLC_P__NP, WLC_P__FA, wlc_d%rand_stat)
            else
                call initchem(wlc_d%AB, wlc_p%NT, WLC_P__NMPP, WLC_P__NBPM, WLC_P__NP, WLC_P__FA, WLC_P__LAM, wlc_d%rand_stat)
            endif

            ! initialize methalation sequence
            if (WLC_P__VARIABLE_CHEM_STATE) then
                if (WLC_P__CHEM_SEQ_FROM_FILE) then
                    print*, "Loding input meth seq..."
                    iostr='input/meth'
                    call wlcsim_params_loadMeth(wlc_d,iostr)
                else
                    call initchem(wlc_d%meth, wlc_p%NT, WLC_P__NMPP, WLC_P__NBPM, WLC_P__NP, WLC_P__F_METH, WLC_P__LAM_METH, wlc_d%rand_stat)
                endif
            endif

            ! calculate volumes of bins
            if (WLC_P__CONFINETYPE.eq.'sphere') then
                call MC_calcVolume(wlc_p%NBINX, WLC_P__DBIN, &
                                WLC_P__LBOX_X, wlc_d%Vol, wlc_d%rand_stat)
            else
                do I = 1,NBin
                    wlc_d%Vol(I) = WLC_P__DBIN**3
                enddo
            endif
        endif

        ! -------------------------------------------
        !
        !  Set up binning proceedure for keeping track of neighbors
        !
        ! ------------------------------------------
        if (WLC_P__NEIGHBOR_BINS) then
            !  Set up binning object
            setBinSize = [WLC_P__LBOX_X, WLC_P__LBOX_Y, WLC_P__LBOX_Z] ! size of bin
            setMinXYZ = [0.0,0.0,0.0]  ! location of corner of bin
            setBinShape = [10,10,10]   ! Specify first level of binning
            call constructBin(wlc_d%bin,setBinShape,setMinXYZ,setBinSize)
            do i=1,NT
                call addBead(wlc_d%bin,wlc_d%R,NT,i)
            enddo
        endif

        ! initialize energies
        !-------------------------------------
        !
        !  Set all energies to zero in case they aren't set later
        !
        !--------------------------------------
        wlc_d%eElas       = 0.0_dp ! Elastic force
        wlc_d%eChi        = 0.0_dp ! CHI energy
        wlc_d%eKap        = 0.0_dp ! KAP energy
        wlc_d%eCouple     = 0.0_dp ! Coupling
        wlc_d%eBind       = 0.0_dp ! binding energy
        wlc_d%eMu         = 0.0_dp ! chemical potential energy
        wlc_d%eField      = 0.0_dp ! Field energy
        wlc_d%eSelf       = 0.0_dp ! repulsive lennard jones on closest approach self-interaction energy (polymer on polymer)
        wlc_d%eMaierSaupe = 0.0_dp ! Maier Saupe energy
        wlc_d%DEELAS      = 0.0_dp ! Change in bending energy
        wlc_d%DECouple    = 0.0_dp ! Coupling energy
        wlc_d%DEChi       = 0.0_dp ! chi interaction energy
        wlc_d%DEKap       = 0.0_dp ! compression energy
        wlc_d%Debind      = 0.0_dp ! Change in binding energy
        wlc_d%DeMu        = 0.0_dp ! Change in chemcial potential energy
        wlc_d%DEField     = 0.0_dp ! Change in field energy
        wlc_d%DESelf      = 0.0_dp ! change in self interaction energy
        wlc_d%ECon        = 0.0_dp ! Confinement Energy
        wlc_d%deMaierSaupe= 0.0_dp ! change in Maier Saupe energy
        wlc_d%NPHI = 0  ! NUMBER o phi values that change, i.e. number of bins that were affected

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
        para(8) = WLC_P__LBOX_X
        para(9) = wlc_p%LHC
        para(10) = wlc_p%VHC
    end function pack_as_para


    subroutine printDescription(wlc_p)
        implicit none
        type(wlcsim_params), intent(in) :: wlc_p
        print*, "---------------System Description---------------"
        print*, " type of simulation, codeName", WLC_P__CODENAME
        print*, " WLC, DSSWLC, GC, simType", wlc_p%SIMTYPE
        print*, "Bead variables:"
        print*, " Total number of beads, NT = ", wlc_p%NT
        print*, " Number of beads in a polymer, NB = ", WLC_P__NB
        print*, " Number of monomers in a polymer, nMpP = ", WLC_P__NMPP
        print*, " Number of polymers, NP = ",WLC_P__NP
        print*, " Number of beads in a monomer, nBpM = ", WLC_P__NBPM
        print*, " fraction Methalated", WLC_P__F_METH
        print*, " LAM_METH", WLC_P__LAM_METH
        print*, " "
        print*, "Length and volume Variables:"
        print*, " persistance length =",(wlc_p%L0/(2.0_dp*wlc_p%EPS))
        print*, " length of each polymer in simulation, l = ",WLC_P__L
        print*, " twist persistence length, lt", WLC_P__LT
        print*, " lbox = ", WLC_P__LBOX_X, WLC_P__LBOX_Y, WLC_P__LBOX_Z
        print*, " Number of bins in x direction", &
                   wlc_p%NBINX(1), wlc_p%NBINX(2),wlc_p%NBINX(3)
        print*, " Number of bins", wlc_p%NBIN
        print*, " spatial descritation dbin = ",WLC_P__DBIN
        print*, " L0 = ", wlc_p%L0
        print*, " bead volume V = ", WLC_P__BEADVOLUME
        print*, " number of kuhn lengths between beads, eps ", wlc_p%EPS
        print*, " "
        print*, "Energy Variables"
        print*, " elasticity EPS =", wlc_p%EPS
        print*, " solvent-polymer CHI =",wlc_p%CHI
        print*, " compression cof, KAP =", wlc_p%KAP
        print*, " field strength, hA =", wlc_p%HA
        print*, " -energy of binding unmethalated ", WLC_P__EU," more positive for favorable binding"
        print*, " -energy of binding methalated",WLC_P__EM
        print*, " HP1_Binding energy parameter", wlc_p%HP1_BIND
        print*, " chemical potential of HP1, mu", wlc_p%MU
        print*, " bend-shear coupling parameter, eta ", wlc_p%ETA
        print*, " "
        print*, "Time Variables"
        print*, " stepsPerExchange", WLC_P__STEPSPEREXCHANGE
        print*, " nReplicaExchangePerSavePoint", WLC_P__NREPLICAEXCHANGEPERSAVEPOINT
        print*, " numSavePoints", WLC_P__NUMSAVEPOINTS
        print*, " stepsPerSave", WLC_P__STEPSPERSAVE
        print*, " "
        print*, "Switches:"
        print*, " confinetype:",WLC_P__CONFINETYPE
        print*, " initCondType:",WLC_P__INITCONDTYPE
        print*, " ring:", WLC_P__RING
        print*, " twist:", WLC_P__TWIST
        print*, " "
        print*, "---------------------------------------------"

    end subroutine

    subroutine tweak_param_defaults(wlc_p, wlc_d)
        implicit none
        type(wlcsim_params), intent(inout) :: wlc_p
        type(wlcsim_data), intent(inout) :: wlc_d

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

        ! if we're not using field interactions
        ! energies, then this should never be on
        if ((.not. wlc_p%FIELD_INT_ON) .and. wlc_p%MOVEON(9)/=0) then
            wlc_p%MOVEON(7) = 0  ! Change in Binding state
            print*, "turning off movetype 7, binding, becuase unneeded"
        endif
        wlc_p%MOVEON(8) = 0  ! Chain flip ! TOdo not working
        ! if number of polymers is 1, this should never be on
        if (WLC_P__NP < 2 .and. wlc_p%MOVEON(9)/=0) then
            wlc_p%MOVEON(9) = 0
            print*, "Turning off movetype 9, chain exchange, because <2 polymers"
        endif

        if (isnan(wlc_p%MINWINDOW(1))) wlc_p%MINWINDOW(1) = dble(min(10,WLC_P__NB))
        if (isnan(wlc_p%MINWINDOW(2))) wlc_p%MINWINDOW(2) = dble(min(10,WLC_P__NB))
        if (isnan(wlc_p%MINWINDOW(3))) wlc_p%MINWINDOW(3) = dble(min(10,WLC_P__NB))
        if (isnan(wlc_p%MINWINDOW(7))) wlc_p%MINWINDOW(7) = dble(min(10,WLC_P__NB))

        ! Solution
        !WLC_P__LBOX_X = wlc_p%NBINX(1)*WLC_P__DBIN
        !WLC_P__LBOX_Y = wlc_p%NBINX(2)*WLC_P__DBIN
        !WLC_P__LBOX_Z = wlc_p%NBINX(3)*WLC_P__DBIN
        wlc_p%NBIN = wlc_p%NBINX(1)*wlc_p%NBINX(2)*wlc_p%NBINX(3)

        if (WLC_P__CODENAME == 'brad') then
            ! initialize windows to number of beads
            wlc_p%MAXWINDOW = WLC_P__NB         ! Max Size of window for bead selection
            wlc_p% MinWindoW  = 1.0_dp         ! Min Size of window for bead selection

            ! Window amplitudes
            wlc_p%MINAMP = 0.0_dp ! minium amplitude
            wlc_p%MINAMP(1) = 0.07_dp*pi
            wlc_p%MINAMP(2) = 0.01_dp*WLC_P__L/WLC_P__NB
            wlc_p%MAXAMP = 2.0_dp*pi
            wlc_p%MAXAMP(2) = WLC_P__LBOX_X
            wlc_p%MAXAMP(6) = WLC_P__LBOX_X

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
        if (WLC_P__RING) then
            wlc_p%MOVEON(3) = 0
        endif

    end subroutine

    subroutine wlcsim_params_recenter(wlc_d)
    !  Prevents drift in periodic BC
        implicit none
        type(wlcsim_data), intent(inout) :: wlc_d
        integer IB, I, J   ! Couners
        real(dp) R0(3)  ! Offset to move by
        do I = 1,WLC_P__NP
            IB=WLC_P__NB * (I-1) + 1
            R0(1) = floor(wlc_d%R(1,IB)/WLC_P__LBOX_X)*WLC_P__LBOX_X
            R0(2) = floor(wlc_d%R(2,IB)/WLC_P__LBOX_Y)*WLC_P__LBOX_Y
            R0(3) = floor(wlc_d%R(3,IB)/WLC_P__LBOX_Z)*WLC_P__LBOX_Z
            if (abs(R0(1)*R0(2)*R0(3)) .gt. 0.0001_dp) then
                do J = 1,WLC_P__NB
                    wlc_d%R(1,IB) = wlc_d%R(1,IB)-R0(1)
                    wlc_d%R(2,IB) = wlc_d%R(2,IB)-R0(2)
                    wlc_d%R(3,IB) = wlc_d%R(3,IB)-R0(3)
                    IB = IB + 1
                enddo
            endif
        enddo
    end subroutine

    subroutine printSimInfo(i, wlc_d)
    ! print out current simulation metainformation
        implicit none
        type(wlcsim_data), intent(in) :: wlc_d
        integer, intent(in) :: i
        print*, 'Current time ', wlc_d%time
        print*, 'Time point ', wlc_d%time_ind, ' out of ', WLC_P__STEPSPERSAVE*WLC_P__NUMSAVEPOINTS
        print*, 'Save point ', i, ' out of ', WLC_P__NUMSAVEPOINTS
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
        print*, "ECHI l=2", wlc_d%eMaierSaupe
        print*, "EField", wlc_d%EField
        print*, "EKAP", wlc_d%EKAP
        print*, "ebind", wlc_d%ebind
        print*, "eMu", wlc_d%eMu
    end subroutine

    subroutine calcTotalPolymerVolume(wlc_p,wlc_d,totalVpoly)
        implicit none
        type(wlcsim_params), intent(in) :: wlc_p
        type(wlcsim_data), intent(in) :: wlc_d
        integer i
        real(dp), intent(out) :: totalVpoly
        real(dp) VV
        totalVpoly=0.0_dp
        do I = 1,wlc_p%NBIN
            VV = wlc_d%Vol(I)
            !if (VV.le.0.1_dp) cycle
            totalVpoly = totalVpoly + VV*(wlc_d%PHIA(I) + wlc_d%PHIB(I))
        enddo
        print*, "Total volume of polymer from density", totalVpoly,&
                " and from beads ",wlc_p%NT*WLC_P__BEADVOLUME

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
        do I = 1,wlc_p%NBIN
            VV = wlc_d%Vol(I)
            if (VV.le.0.1_dp) cycle
            PHIPOLY = wlc_d%PHIA(I) + wlc_d%PHIB(I)
            EChi = VV*(wlc_p%CHI/WLC_P__BEADVOLUME)*PHIPoly*(1.0_dp-PHIPoly)
            ECouple = VV*wlc_p%HP1_BIND*(wlc_d%PHIA(I))**2
            if(PHIPoly > 1.0_dp) then
            EKap = VV*(wlc_p%KAP/WLC_P__BEADVOLUME)*(PHIPoly-1.0_dp)**2
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
        do I = 1,nMovetypes
            if (wlc_p%MOVEON(i).eq.1) then
                write(*,"(f8.5,2f8.2,1A1,1A20)") wlc_d%phit(i), wlc_d%MCAMP(i),&
                    wlc_d%WindoW(i),' ', moveNames(i)
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
        do I = 1,wlc_p%NBIN
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

        do IX = 1,wlc_p%NBINX(1)
            do IY = 1,wlc_p%NBINX(2)
                do IZ = 1,wlc_p%NBINX(3)
                    indBin = IX + &
                        (IY-1)*wlc_p%NBINX(1) + &
                        (IZ-1)*wlc_p%NBINX(1)*wlc_p%NBINX(2)
                    wlc_d%PHIH(indBin) = dsin(WLC_P__K_FIELD*WLC_P__DBIN*dble(IX))
                enddo
            enddo
        enddo
        return
    end subroutine

    subroutine wlcsim_params_loadAB(wlc_d,fileName)
    ! Loads AB for file...has not been tested
        implicit none
        type(wlcsim_data), intent(inout) :: wlc_d
        character(MAXFILENAMELEN), intent(in) :: fileName ! file name to load from
        integer IB, I, J ! counters
        open (unit = inFileUnit, file = fileName, status = 'OLD')
        IB = 1
        do I = 1,WLC_P__NP
            do J = 1,WLC_P__NB
                read(inFileUnit,"(I2)") wlc_d%AB(IB)
                IB = IB + 1
            enddo
        enddo
        close(inFileUnit)
    end subroutine
    subroutine wlcsim_params_loadMeth(wlc_d,fileName)
    ! Loads Methalation for file...has not been tested
        implicit none
        type(wlcsim_data), intent(inout) :: wlc_d
        character(MAXFILENAMELEN), intent(in) :: fileName ! file name to load from
        integer IB, I, J ! counters
        open (unit = inFileUnit, file = fileName, status = 'OLD')
        IB = 1
        do I = 1,WLC_P__NP
            do J = 1,WLC_P__NB
                read(inFileUnit,"(I2)") wlc_d%meth(IB)
                IB = IB + 1
            enddo
        enddo
        close(inFileUnit)
    end subroutine

    subroutine wlcsim_params_saveR(wlc_d,fileName,repeatingBC)
    ! Writes R and AB to file for analysis
    ! Rx  Ry  Rz AB
        implicit none
        logical, intent(in) :: repeatingBC  ! 1 for reapeating boundary conditions
        integer I,J,IB  ! counters
        type(wlcsim_data), intent(in) :: wlc_d
        character(MAXFILENAMELEN), intent(in) :: fileName
        character(MAXFILENAMELEN) fullName
        LOGICAL isfile
        fullName = trim(fileName) // trim(wlc_d%repSuffix)
        fullName = trim(fullName)
        inquire(file = fullName, exist = isfile)
        if (isfile) then
            open (unit = outFileUnit, file = fullName, status ='OLD', POSITION = "append")
        else
            open (unit = outFileUnit, file = fullName, status = 'NEW')
        endif
        IB = 1
        if (repeatingBC) then
           do I = 1,WLC_P__NP
              do J = 1,WLC_P__NB
                 if (WLC_P__SAVEAB) then
                    write(outFileUnit,"(3f10.3,I2)") &
                         wlc_d%R(1,IB)-floor(wlc_d%R(1,IB)/WLC_P__LBOX_X)*WLC_P__LBOX_X, &
                         wlc_d%R(2,IB)-floor(wlc_d%R(2,IB)/WLC_P__LBOX_Y)*WLC_P__LBOX_Y, &
                         wlc_d%R(3,IB)-floor(wlc_d%R(3,IB)/WLC_P__LBOX_Z)*WLC_P__LBOX_Z, &
                         wlc_d%AB(IB)
                 else
                    write(outFileUnit,"(3f10.3)") &
                         wlc_d%R(1,IB)-floor(wlc_d%R(1,IB)/WLC_P__LBOX_X)*WLC_P__LBOX_X, &
                         wlc_d%R(2,IB)-floor(wlc_d%R(2,IB)/WLC_P__LBOX_Y)*WLC_P__LBOX_Y, &
                         wlc_d%R(3,IB)-floor(wlc_d%R(3,IB)/WLC_P__LBOX_Z)*WLC_P__LBOX_Z
                 endif
                 IB = IB + 1
              enddo
           enddo
           print*, "Error in wlcsim_params_saveR"
           print*, "Are you sure you want repeating BC?"
           print*, "Quinn put this in ages ago but never implemented it...."
           stop 1
        else
           do I = 1,WLC_P__NP
              do J = 1,WLC_P__NB
                  if (WLC_P__SAVEAB) then
                     if (WLC_P__CHANGINGCHEMICALIDENTITY) then
                         write(outFileUnit,"(3f10.3,2I3)") &
                                wlc_d%R(1,IB),wlc_d%R(2,IB),wlc_d%R(3,IB),wlc_d%AB(IB),wlc_d%METH(IB)
                     else
                         write(outFileUnit,"(3f10.3,I2)") &
                                wlc_d%R(1,IB),wlc_d%R(2,IB),wlc_d%R(3,IB),wlc_d%AB(IB)
                     endif
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
        if (wlc_p%CHI_L2_ON) then
            do I = 1,wlc_p%NBIN
                write(outFileUnit,"(7f7.2)") wlc_d%PHIA(I),wlc_d%PHIB(I),wlc_d%PHI_l2(:,I)
            enddo
        else
            do I = 1,wlc_p%NBIN
                write(outFileUnit,"(2f7.2)") wlc_d%PHIA(I),wlc_d%PHIB(I)
            enddo
        endif

        close(outFileUnit)
    end subroutine

    subroutine wlcsim_params_saveU(wlc_d,fileName,stat)
    ! Saves U to ASCII file for analisys
        implicit none
        integer I,J,IB  ! counters
        type(wlcsim_data), intent(in) :: wlc_d
        character(MAXFILENAMELEN), intent(in) :: fileName
        character(MAXFILENAMELEN) fullName
        character(len = *), intent(in) :: stat
        fullName=  trim(fileName) // trim(wlc_d%repSuffix)
        open (unit = outFileUnit, file = fullName, status = stat)
        IB = 1
        do I = 1,WLC_P__NP
            do J = 1,WLC_P__NB
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
            write(outFileUnit,"(I8)") WLC_P__NMPP  ! 2 Number of monomers in a polymer
            write(outFileUnit,"(I8)") WLC_P__NB ! 3 Number of beads in a polymer
            write(outFileUnit,"(I8)") WLC_P__NP ! 4 Number of polymers in simulation
            write(outFileUnit,"(I8)") wlc_p%NT ! 5 Number of beads in simulation
            write(outFileUnit,"(I8)") WLC_P__NBPM  ! 6 Number of beads per monomer

            write(outFileUnit,"(f10.5)") wlc_p%L0    ! Equilibrium segment length
            write(outFileUnit,"(f10.5)") wlc_p%CHI  ! 8  initail CHI parameter value
            write(outFileUnit,"(f10.5)") WLC_P__LBOX_X  ! 10 Lenth of box
            write(outFileUnit,"(f10.5)") WLC_P__EU    ! Energy unmethalated
            write(outFileUnit,"(f10.5)") WLC_P__EM    ! 12 Energy methalated
            write(outFileUnit,"(f10.5)") wlc_p%HP1_BIND ! Energy of HP1 binding
            write(outFileUnit,"(f10.5)") (wlc_p%L0/wlc_p%EPS) ! 14 Khun lenth
            write(outFileUnit,"(A)") "-999"  ! for historic reasons
            write(outFileUnit,"(f10.5)") WLC_P__F_METH  ! methalation fraction
            write(outFileUnit,"(f10.5)") WLC_P__LAM_METH  ! methalation lambda
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
                       "  ebend    |  eparll   |  EShear   |  ECoupl   |  E Kap    |  E Chi    |",&
                       "  EField   |  ebind    |   x_Mu    |  Couple   |   Chi     |   mu      |",&
                       "   Kap     |  Field    |   x_MS    |  chi_l2   |  E_Mu     |"
        endif
        write(outFileUnit,"(2I5, 9f12.1,5f12.5,f12.1,f12.5,f12.2)") save_ind, wlc_d%id, &
            wlc_d%EELAS(1), wlc_d%EELAS(2), wlc_d%EELAS(3), wlc_d%ECouple, &
            wlc_d%EKap, wlc_d%ECHI, wlc_d%EField, wlc_d%ebind, wlc_d%x_Mu, &
            wlc_p%HP1_BIND*wlc_p%COUPLE_ON, wlc_p%CHI*wlc_p%CHI_ON, wlc_p%MU, wlc_p%KAP*wlc_p%KAP_ON,&
            wlc_p%HA, wlc_d%x_maierSaupe, wlc_p%CHI_L2,wlc_d%EMu
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

        if (WLC_P__SAVEPHI) then
            write(filename,num2strFormatString) save_ind
            filename = trim(adjustL(outfile_base)) // 'phi' // trim(adjustL(filename))
            call wlcsim_params_savePHI(wlc_p,wlc_d,filename)
        endif

        if (WLC_P__SAVER) then
            if (WLC_P__APPEND_R) then
                filename = trim(adjustL(outfile_base)) // 'rAll'
            else
                write(filename,num2strFormatString) save_ind
                filename = trim(adjustL(outfile_base)) // 'r' // trim(adjustL(filename))
            endif
            call wlcsim_params_saveR(wlc_d,filename,.false.)
        endif

        if (WLC_P__SAVEU) then
            write(filename,num2strFormatString) save_ind
            filename = trim(adjustL(outfile_base)) // 'u' // trim(adjustL(filename))
            call wlcsim_params_saveU(wlc_d,filename,stat)
        endif

        if (WLC_P__COLLISIONDETECTIONTYPE /= 0) then
            filename = trim(adjustL(outfile_base)) // 'coltimes'
            open(unit = outFileUnit, file = filename, status = 'REPLACE')
            do ind = 1,wlc_p%NT
                write(outFileUnit,*) (wlc_d%coltimes(ind,j), j = 1,wlc_p%NT)
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

    subroutine get_renormalized_chain_params(wlc_p)
    !     Setup the parameters for the simulation
    !
    !     1. Determine the simulation type
    !     2. Evaluate the polymer elastic parameters
    !     3. Determine the parameters for Brownian dynamics simulation
        implicit none

        integer i,ind
        real(dp) m

        type(wlcsim_params), intent(inout) :: wlc_p
        REAL(dp) :: pvec(679, 8) ! array holding dssWLC params calculated by Elena

        !Calculate total number of beads
        wlc_p%NT = WLC_P__NB*WLC_P__NP

        if (WLC_P__NB == 1.0d0) then
            ! since we use "DEL" as an intermediate, we need at least two beads
            PRinT*, 'Some intermediate calculations used require at least two beads, 1 requested.'
            STOP 1
        endif

        ! calculate metrics that don't change between WLC, ssWLC, GC
        if (WLC_P__RING) then
            wlc_p%DEL = WLC_P__L/WLC_P__LP/(WLC_P__NB)
        else
            wlc_p%DEL = WLC_P__L/WLC_P__LP/(WLC_P__NB-1.0_dp)
        ENDif
        ! std dev of interbead distribution of nearest possible GC, used to initialize sometimes
        wlc_p%SIGMA = sqrt(2.0_dp*WLC_P__LP*WLC_P__L/3.0_dp)/real(WLC_P__NB - 1)

    !     Load the tabulated parameters

        open (UNIT = 5,FILE = 'input/dssWLCparams',STATUS = 'OLD')
        do I = 1,679
            READ(5,*) PVEC(I,1),PVEC(I,2),PVEC(I,3),PVEC(I,4),PVEC(I,5),PVEC(I,6),PVEC(I,7),PVEC(I,8)
        ENDdo
        CLOSE(5)


    !     Setup the parameters for WLC simulation

        ! if del < 0.01
        if (wlc_p%DEL < PVEC(1,1)) then
            PRinT*, 'It has never been known if the WLC code actually works.'
            PRinT*, 'An entire summer student (Luis Nieves) was thrown at this'
            PRinT*, 'problem and it is still not solved.'
            stop 1
            wlc_p%EB = WLC_P__LP/wlc_p%DEL
            wlc_p%GAM = wlc_p%DEL
            wlc_p%XIR = WLC_P__L/WLC_P__LP/WLC_P__NB
            wlc_p%SIMTYPE = 1

    !    Setup the parameters for GC simulation

        ! if del > 10
        elseif (wlc_p%DEL > PVEC(679,1)) then
            wlc_p%EPAR = 1.5/wlc_p%DEL
            wlc_p%GAM = 0.0_dp
            wlc_p%SIMTYPE = 3
            wlc_p%XIR = WLC_P__L/WLC_P__NB/WLC_P__LP

    !    Setup the parameters for ssWLC simulation
        ! if 0.01 <= del <= 10
        else !  if (DEL >= PVEC(1,1).AND.DEL <= PVEC(679,1)) then
            wlc_p%SIMTYPE = 2

        ! find(del < pvec, 1, 'first')
        inD = 1
        do while (wlc_p%DEL > PVEC(inD,1))
            inD = inD + 1
        enddo

        !     Perform linear interpolations
        I = 2
        M = (PVEC(inD,I)-PVEC(inD-1,I))/(PVEC(inD,1)-PVEC(inD-1,1))
        wlc_p%EB = M*(wlc_p%DEL-PVEC(inD,1)) + PVEC(inD,I)

        I = 3
        M = (PVEC(inD,I)-PVEC(inD-1,I))/(PVEC(inD,1)-PVEC(inD-1,1))
        wlc_p%GAM = M*(wlc_p%DEL-PVEC(inD,1)) + PVEC(inD,I)

        I = 4
        M = (PVEC(inD,I)-PVEC(inD-1,I))/(PVEC(inD,1)-PVEC(inD-1,1))
        wlc_p%EPAR = M*(wlc_p%DEL-PVEC(inD,1)) + PVEC(inD,I)

        I = 5
        M = (PVEC(inD,I)-PVEC(inD-1,I))/(PVEC(inD,1)-PVEC(inD-1,1))
        wlc_p%EPERP = M*(wlc_p%DEL-PVEC(inD,1)) + PVEC(inD,I)

        I = 6
        M = (PVEC(inD,I)-PVEC(inD-1,I))/(PVEC(inD,1)-PVEC(inD-1,1))
        wlc_p%ETA = M*(wlc_p%DEL-PVEC(inD,1)) + PVEC(inD,I)

        I = 7
        M = (PVEC(inD,I)-PVEC(inD-1,I))/(PVEC(inD,1)-PVEC(inD-1,1))
        wlc_p%XIU = M*(wlc_p%DEL-PVEC(inD,1)) + PVEC(inD,I)

        ! The values read in from file are all non-dimentionalized by the
        ! persistance length.  We now re-dimentionalize them.
        ! We also divied by DEL which is also re-dimentionalized.

        wlc_p%EB = WLC_P__LP*wlc_p%EB/(wlc_p%DEL*WLC_P__LP)
        wlc_p%EPAR = wlc_p%EPAR/(wlc_p%DEL*WLC_P__LP*WLC_P__LP)
        wlc_p%EPERP = wlc_p%EPERP/(wlc_p%DEL*WLC_P__LP*WLC_P__LP)
        wlc_p%GAM = wlc_p%DEL*WLC_P__LP*wlc_p%GAM
        wlc_p%ETA = wlc_p%ETA/WLC_P__LP
        wlc_p%XIU = wlc_p%XIU*WLC_P__L/WLC_P__NB/WLC_P__LP
        wlc_p%XIR = WLC_P__L/WLC_P__LP/WLC_P__NB
        wlc_p%DT = 0.5*wlc_p%XIU/(wlc_p%EPERP*wlc_p%GAM**2.)

        ! wlc_p%L0 = wlc_p%GAM  ! not sure why this was included
        endif

        return
    end subroutine get_renormalized_chain_params

end module params
