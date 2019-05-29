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
    use precision, only: dp, eps, epsapprox, pi, nan
    use inputparams, only: MAXPARAMLEN
    use binning, only: constructBin, binType, addBead
    use precalc_spider, only: spider, load_precalc_spiders, get_highestNumberOfLegs

    implicit none

    public

    !!!     hardcoded params. will need to change if certain parts of code change
    ! number of wlc_p move types
    integer, parameter :: nMoveTypes = 12
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
          'superReptation      ','spider              '/)

    !!!     universal constants
    ! fully accurate, adaptive precision
    ! ! won't get optimized away by compiler, see e.g.
    ! ! https://software.intel.com/en-us/forums/intel-visual-fortran-compiler-for-windows/topic/294680
    ! real(dp) :: nan = transfer((/ Z'00000000', Z'7FF80000' /),1.0_dp)
    ! the following would be preferred, but generated compilation errors...
    ! real(dp) :: one = 1.0_dp
    ! real(dp) :: inf = -log(one - one)
    ! the following would be preferred, but generated compilation errors...
    !real(dp) :: inf = ieee_value(inf, ieee_positive_inf)
    real(dp) :: max_wlc_l0 = 0.01_dp
    real(dp) :: maxWlcDelta = 10.0_dp
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
        real(dp) dt ! sets time scale of simulation
        real(dp) l0       ! Path length between beads. (meaning unknown for gaussian chain?)
        real(dp) gam    ! average equilibrium interbead spacing
        real(dp) eta    ! bend-shear coupling parameter
        real(dp) xir    ! drag per unit persistence length
        real(dp) sigma  ! variance of interbead position distribution of appropriately renormalized gaussian chain
        real(dp) xiu    ! rotational drag
        real(dp) eps      ! number of kuhn lengths between beads
        real(dp) del      ! number of persistence lengths between beads
        real(dp) lhc    !TOdo something to do with intrapolymer interaction strength
        real(dp) vhc    !TOdo something to do with intrapolymer interaction strength, fill in defaults, etc
        real(dp) eb     ! effective bending energy for ssWLC
        real(dp) eperp  ! effective shearing energy for ssWLC
        real(dp) epar   ! effective stretch energy for ssWLC
        real(dp) etwist

    !   boundary/box things
        integer NBin     ! Number of bins

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
        logical field_int_on_currently ! include field interactions (e.g. A/B interactions) uses many of the same functions as the chemical identity/"meth"ylation code, but energies are calcualted via a field-based approach
        logical chi_l2_on

    end type

    real(dp), allocatable, dimension(:,:):: wlc_R   ! Conformation of polymer chains
    real(dp), allocatable, dimension(:,:):: wlc_R_period   ! Conformation of polymer chains subracted to first period with lower corner at the origin
    real(dp), allocatable, dimension(:,:):: wlc_U   ! Conformation of polymer chains
    real(dp), allocatable, dimension(:,:):: wlc_V   ! Conformation of polymer chains
    real(dp), allocatable, dimension(:,:):: wlc_RP !Test Bead positions - only valid from IT1 to IT2
    real(dp), allocatable, dimension(:,:):: wlc_UP !Test target vectors - only valid from IT1 to IT2
    real(dp), allocatable, dimension(:,:):: wlc_VP !Test target vectors - only valid from IT1 to IT2
    integer, allocatable, dimension(:):: wlc_ExplicitBindingPair ! List of other points bound to this one
    integer, allocatable, dimension(:):: wlc_network_start_index ! Index in of first in wlc_network_start_index
    integer, allocatable, dimension(:):: wlc_other_beads ! Other beads attached to beads see wlc_network_start_index
    logical, allocatable, dimension(:):: wlc_external_bind_points ! Random points attached to boundary
    real(dp), allocatable, dimension(:):: wlc_PHIA ! Volume fraction of A
    real(dp), allocatable, dimension(:):: wlc_PHIB ! Volume fraction of B
    real(dp), allocatable, dimension(:,:):: wlc_PHI_l2 ! l=2 oreientational field
    real(dp), allocatable, dimension(:):: wlc_PHIH ! Applied field in hamiltonian
    real(dp), allocatable, dimension(:,:):: wlc_PHIH_L2 ! Applied field in hamiltonian
    real(dp), allocatable, dimension(:):: wlc_Vol  ! Volume fraction of A
    integer, allocatable, dimension(:):: wlc_AB    ! Chemical identity of beads
    integer, allocatable, dimension(:):: wlc_ABP   ! Test Chemical identity of beads
    integer, allocatable, dimension(:):: wlc_METH  ! Methalation state of beads
    real(dp), allocatable, dimension(:):: wlc_DPHIA    ! Change in phi A
    real(dp), allocatable, dimension(:):: wlc_DPHIB    ! Change in phi A
    real(dp), allocatable, dimension(:,:):: wlc_DPHI_l2 ! change in l=2 oreientational field
    integer, allocatable, dimension(:) :: wlc_indPHI   ! indices of the phi
    integer, allocatable, dimension(:):: wlc_ind_in_list ! index in indPhi
    ! simulation times at which (i,j)th bead pair first collided
    real(dp), allocatable, dimension(:,:) :: wlc_coltimes
    real(dp) :: wlc_wr
    type(spider), allocatable, dimension(:) :: wlc_spiders ! spiders based on polymer network
    integer wlc_numberOfSpiders
    integer wlc_spider_id
    real(dp) wlc_spider_dr(3)

    type(binType) wlc_bin ! Structure for keeping track of neighbors

    !   Twist variables
    real(dp), ALLOCATABLE :: wlc_CROSS(:,:)   !Matrix of information for crossings in a 2-D projection of the polymer
    real(dp), ALLOCATABLE :: wlc_CROSSP(:,:)  !Matrix of crossings for the trial configuration
    integer wlc_NCROSS
    integer wlc_NCROSSP
    integer wlc_CrossSize


    !   Monte Carlo Variables (for adaptation)
    real(dp) wlc_MCAMP(nMoveTypes) ! Amplitude of random change
    real(dp) wlc_WindoW(nMoveTypes)         ! Size of window for bead selection
    integer wlc_SUCCESS(nMoveTypes)        ! Number of successes
    integer wlc_ATTEMPTS(nMoveTypes)        ! Number of successes
    integer wlc_successTOTAL(nMoveTypes)               !Total number of successes
    real(dp) wlc_PHit(nMoveTypes) ! hit rate

    integer wlc_NPHI  ! NUMBER o phi values that change, i.e. number of bins that were affected
    integer, allocatable, dimension(:) :: wlc_bendPoints ! index of left end of bends presint in chain
    integer wlc_nBend ! the number of points bent
    integer wlc_maxNBend ! the number of points bent
    integer, allocatable, dimension(:) :: wlc_pointsMoved  ! Indicies of points moved
    integer wlc_nPointsMoved

    !   Parallel tempering variables
    integer wlc_numProcesses !number of MPI processes running
    integer wlc_rep  ! which replica am I
    integer wlc_id   ! which thread am I
    integer wlc_error  ! MPI error
    integer, allocatable, dimension(:) ::  wlc_LKs    !Vector of linking numbers for replicas
    integer wlc_nLKs      !Number of linking number replicas to parallel temper over
    real(dp), allocatable, dimension(:) :: wlc_Wrs !Vector of writhe for each replica
    real(dp), allocatable, dimension(:,:) :: wlc_eelasREPLICAS !elastic energies of replicas
    integer wlc_replicaSTART !index for replica to start with for exchange loop
    integer wlc_replicaEND   !index for replica to end at for exchange loop
    integer, allocatable, dimension(:) :: wlc_nTRIALup !number of times this replica has attempted to swap with replica above
    integer, allocatable, dimension(:) :: wlc_nTRIALdown !number of times this replica has attempted to swap with replica below
    integer, allocatable, dimension(:) :: wlc_nSWAPup !number of times this replica has swapped with replica above
    integer, allocatable, dimension(:) :: wlc_nSWAPdown !number of times this replica has swapped with replica below
    integer, allocatable, dimension(:) :: wlc_nodeNUMBER !vector of replicas indices for nodes
    character(MAXFILENAMELEN) wlc_repSuffix    ! prefix for writing files


    !   random number generator state
    type(random_stat) wlc_rand_stat
    integer wlc_rand_seed

    !   indices
    integer wlc_mc_ind                  ! current save point index for mc
    integer wlc_ind_exchange            ! number of exchange moves since last save point
    integer wlc_time_ind                ! current time point
    real(dp) wlc_time

    !   nucleosomes
    integer, allocatable, dimension(:) :: wlc_basepairs
    integer, allocatable, dimension(:) :: wlc_nucleosomeWrap



contains

    subroutine set_param_defaults(wlc_p)
        use energies, only: set_up_energyOf
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


        wlc_p%EPS=WLC_P__L0/(2.0_dp*WLC_P__LP)

        call set_up_energyOf()

        wlc_p%lhc = NAN ! I have no idea what this does
        wlc_p%vhc = NAN ! I have no idea what this does
        wlc_p%field_int_on_currently = WLC_P__FIELD_INT_ON ! on by default

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
        wlc_p%PDESIRE(12) = WLC_P__PDESIRE_SPIDER
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
        wlc_p%MAXWINDOW(12) = NAN ! max window spider
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
        wlc_p%MINWINDOW(12) = NAN ! min window spider
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
        wlc_p%MINAMP(12) = WLC_P__MINAMP_SPIDER
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
        wlc_p%MAXAMP(12) = WLC_P__MAXAMP_SPIDER
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
        wlc_p%MOVEON(12) = WLC_P__MOVEON_SPIDER
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
        wlc_p%WINTARGET(12) = NAN
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
        wlc_p%NADAPT(12) = WLC_P__NADAPT_SPIDER
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
        wlc_p%MOVESPERSTEP(12) = WLC_P__MOVESPERSTEP_SPIDER

    end subroutine set_param_defaults





    subroutine idiot_checks(wlc_p)
#if MPI_VERSION
        use mpi
#endif
        implicit none
        type(wlcsim_params), intent(inout) :: wlc_p
        logical err

        if (WLC_P__NEIGHBOR_BINS .and. (WLC_P__CONFINETYPE .ne. 'excludedShpereInPeriodic')&
            .and. (WLC_P__CONFINETYPE .ne. 'sphere')) then
            print*, "The code is untested for Neighbor bins and other confinetypes"
            print*, "No confinement (e.g. infinite volume) should be OK.  As should a fixed confinement"
            print*, "However, if you want a different periodic confiment you should add it to places where R_period is used"
            stop
        endif
        if ((.not. WLC_P__FIELD_INT_ON) .and. (WLC_P__SAVEAB)) then
            print*, "SAVE_AB = True is currently incompatable with Field_int_on=False"
            print*, "because AB is only allocated if field_int_on=True"
        endif
        if (WLC_P__ASYMMETRICALTERNATINGCHEM .and. WLC_P__CHANGINGCHEMICALIDENTITY) then
            print*, "Asymmetric AlternatingChem and changing Chemical Identity is not avaiable."
            stop
        endif
        if ((.not. WLC_P__EXPLICIT_BINDING) .and. WLC_P__INITCONDTYPE ==  "multiRing") then
            print*, "multiRing initial condition not possible if explicit binding is off"
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

        err = WLC_P__EXPLICIT_BINDING .and. WLC_P__MOVEON_CHAIN_EXCHANGE == 1
        call stop_if_err(err, "Explicit binding not set up for exchange move")
        
        err = WLC_P__NETWORK .and. .not. WLC_P__EXPLICIT_BINDING
        call stop_if_err(err, "Network requeires explicit binding")

        err = WLC_P__APPLY_EXTERNAL_FIELD .and. WLC_P__MOVEON_CHAIN_EXCHANGE == 1
        call stop_if_err(err, "External field not set up for exchange move")

        err = WLC_P__MOVEON_REPTATION ==1 .and. WLC_P__LOCAL_TWIST
        call stop_if_err(err, "Reptation move energy calc not set up for twist.")

        err = WLC_P__FRACTIONAL_BIN .and. (WLC_P__CONFINETYPE .ne. 'sphere')
        call stop_if_err(err, "Fractional bin only implimented for sphere")

        err = WLC_P__ENSEMBLE_METH .and. WLC_P__PTON
        call stop_if_err(err,"Parallel tmpering isn't valid for differet Meth profiles")

        call stop_if_err(WLC_P__COLLISIONDETECTIONTYPE == 2, &
            'KD-tree based collision detection not yet implemented.')

        err = (WLC_P__BOUNDARY_TYPE == 'SolidEdgeBin') .and. &
              (WLC_P__FIELD_INT_ON) .and. &
              (WLC_P__CONFINETYPE == 'sphere' .or. WLC_P__CONFINETYPE == 'ecoli')
        call stop_if_err(err,"I don't know how to do SolidEdgeBin for curved boundary")

        if (WLC_P__CODENAME == 'quinn') then
           if ((WLC_P__NBIN_X-WLC_P__NBIN_Y.ne.0).or. &
                (WLC_P__NBIN_X-WLC_P__NBIN_Z.ne.0)) then
              err = (WLC_P__CONFINETYPE.eq.'sphere')
              call stop_if_err(err, "Don't use unequal boundaries for phsere")
              err = WLC_P__INITCONDTYPE.eq.'randomLineSphereBoundary'
              call stop_if_err(err, "You shouldn't put a sphere in and unequal box!")
           endif

           err = WLC_P__NBIN_X*WLC_P__NBIN_Y*WLC_P__NBIN_Z.ne.wlc_p%NBIN
           call stop_if_err(err, "error in mcsim. Wrong number of bins")

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

        err = (WLC_P__CODENAME == 'brad' .and. WLC_P__NT > WLC_P__NB)
        call stop_if_err(err,"Brad's code assumes one chain. Need to update all alexander ...")

        err = (WLC_P__POLY_DISP_TYPE .ne. "None" .and. WLC_P__RING)
        call stop_if_err(err,"writhe and possibly other functions not set up for polydispersity")

        err = (WLC_P__INTERP_BEAD_LENNARD_JONES .and. WLC_P__POLY_DISP_TYPE .ne. "None")
        call stop_if_err(err,"INTERP_BEAD_LENNARD_JONES not setup for polydispersity")

#if MPI_VERSION
    if (WLC_P__PT_TWIST) then
        if (.NOT.WLC_P__TWIST) then
            print *, 'parallel tempering on twist, but twist off'
            stop
        endif
        if (wlc_nLKs + 1.ne.wlc_numProcesses) then
            print *, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            print *, 'number of threads not equal to number of replicas!'
            print *, 'exiting...'
            print *, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            stop
        endif
  endif
#endif
    end subroutine


    subroutine set_parameters( wlc_p)
        use nucleosome, only: setup_nucleosome_constants
        ! Based on Elena's readkeys subroutine
        implicit none
        type(wlcsim_params), intent(out) :: wlc_p

        ! baseline defaults
        call set_param_defaults(wlc_p)

        ! advanced defaults that require some inputs to specify
        call tweak_param_defaults(wlc_p)

        ! get derived parameters that aren't directly input from file

        if (WLC_P__ELASTICITY_TYPE == "constant") then
            call get_renormalized_chain_params(wlc_p)
        elseif (WLC_P__ELASTICITY_TYPE == "nucleosomes") then
            call get_renormalized_chain_params(wlc_p) ! only so that there are constants for initization
            call setup_nucleosome_constants()
        else
            print*, "Elasticity Type ",WLC_P__ELASTICITY_TYPE," not recognized"
            stop
        endif

        !If parallel tempering is on, read the Lks
        if (WLC_P__PT_TWIST) then
            call get_LKs_from_file()
        endif

        call printDescription(wlc_p)
        call idiot_checks(wlc_p)

    end subroutine



    subroutine initialize_wlcsim_data( wlc_p)
        use nucleosome, only: loadNucleosomePositions
        use polydispersity, only: max_chain_length, setup_polydispersity
        use energies, only: set_all_energy_to_zero
#if MPI_VERSION
        use mpi
#endif
        implicit none
        type(wlcsim_params), intent(in)    :: wlc_p
        character(8) datedum  ! trash
        character(10) timedum ! trash
        character(5) zonedum  ! trash
        integer seedvalues(8) ! clock readings
        integer NBin ! total number of bins
        integer i, ii
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
        integer len_file
        nbin = wlc_p%NBIN

#if MPI_VERSION
        call init_MPI()
#endif

        call setup_polydispersity()
        allocate(wlc_R(3,WLC_P__NT))
        if (WLC_P__NEIGHBOR_BINS .and. ((WLC_P__CONFINETYPE == 'excludedShpereInPeriodic') .or. WLC_P__CONFINETYPE == 'none')) then
            allocate(wlc_R_period(3,WLC_P__NT))
        endif
        allocate(wlc_U(3,WLC_P__NT))
        if (WLC_P__LOCAL_TWIST) then
            allocate(wlc_V(3,WLC_P__NT))
        endif
        if (WLC_P__CODENAME /= 'bruno' .OR. WLC_P__NINITMCSTEPS /= 0) then
            allocate(wlc_RP(3,WLC_P__NT))
            allocate(wlc_UP(3,WLC_P__NT))
            if (WLC_P__LOCAL_TWIST) then
                allocate(wlc_VP(3,WLC_P__NT))
            endif
            wlc_RP=nan  ! To prevent accidental use
            wlc_UP=nan
            if (WLC_P__LOCAL_TWIST) then
                wlc_VP=nan
            endif
        endif
        !TOdo these should in principle be inside the following if statement,
        !but it's not clear if that's possible without adding a bunch of dirty
        !if statements deep inside mc_move. which is fine, but I would want to
        !check with quinn *exactly* in which cases they're needed if i do that
        if (WLC_P__FIELD_INT_ON) then
            allocate(wlc_AB(WLC_P__NT))   !Chemical identity aka binding state
            if (WLC_P__CHANGINGCHEMICALIDENTITY) then
                allocate(wlc_ABP(WLC_P__NT))   !Chemical identity aka binding state
                wlc_ABP = INT_MIN
            endif
            if (WLC_P__CHI_L2_ABLE) then
                allocate(wlc_PHI_l2(-2:2,NBin))
                allocate(wlc_dPHI_l2(-2:2,NBin))
            endif
            if (WLC_P__FIELDINTERACTIONTYPE=='AppliedAligningFieldMelt') then
                allocate(wlc_PHIH_l2(-2:2,NBin))
                call load_l2_field(wlc_p)
            endif
            allocate(wlc_PHIA(NBin))
            allocate(wlc_PHIB(NBin))
            allocate(wlc_DPHIA(NBin))
            allocate(wlc_DPHIB(NBin))
            allocate(wlc_indPHI(NBin))
            allocate(wlc_ind_in_list(NBin))
            wlc_ind_in_list = -1  ! -1 stands for not in list
            if (WLC_P__FIELDINTERACTIONTYPE == "ABmelt" .or.\
                WLC_P__FIELDINTERACTIONTYPE == "ABsoluution") then
                allocate(wlc_PhiH(NBin))
            endif
            if (WLC_P__FRACTIONAL_BIN) then
                allocate(wlc_Vol(NBin))
            endif
            do I = 1,NBin
                wlc_PHIA(I) = 0.0_dp
                wlc_PHIB(I) = 0.0_dp
                wlc_indphi(I) = INT_MIN
            enddo
        endif
        if (WLC_P__ELASTICITY_TYPE == "nucleosomes") then
            allocate(wlc_basepairs(WLC_P__NT))
            allocate(wlc_nucleosomeWrap(WLC_P__NT))
        endif
        if (WLC_P__EXPLICIT_BINDING) then
            if (WLC_P__NETWORK) then
                allocate(wlc_network_start_index(WLC_P__NT+1))
                open(unit = 5,file = "input/network_start_index",status = 'OLD')
                do I = 1, WLC_P__NT+1
                    READ(5,*) wlc_network_start_index(I)
                enddo
                close(5)

                open(unit = 5,file = "input/other_beads",status = 'OLD')
                READ(5,*) len_file ! First line is number of lines to follow
                allocate(wlc_other_beads(len_file))
                do I = 1, len_file
                    READ(5,*) wlc_other_beads(I)
                enddo
                close(5)
            else
                allocate(wlc_ExplicitBindingPair(WLC_P__NT))
                if (WLC_P__ENSEMBLE_BIND .and. wlc_id>0) then
                    write(iostr,"(I4)") wlc_id
                    iostr = adjustL(iostr)
                    iostr = trim(iostr)
                    iostr = "input/L393216nloops50000_"//trim(iostr)//".txt"
                    iostr = trim(iostr)
                    print*, "reading ",iostr
                    open(unit = 5,file = iostr,status = 'OLD')
                else
                    open(unit = 5,file = "input/bindpairs",status = 'OLD')
                endif
                do I = 1,WLC_P__NT
                    Read(5,*) wlc_ExplicitBindingPair(I)
                    if (wlc_ExplicitBindingPair(I) .gt. WLC_P__NT) then
                        print*, "Loop to nonexistant bead"
                        stop 1
                    endif
                enddo
                close(5)
                print*, "Read explidit binding"
                print*, wlc_ExplicitBindingPair(1:10)
                print*, "..."
            endif
        endif

        if (WLC_P__EXTERNAL_FIELD_TYPE == 'Random_to_cube_side' .and. WLC_P__APPLY_EXTERNAL_FIELD) then
            allocate(wlc_external_bind_points(WLC_P__NT))
        endif

        allocate(wlc_pointsMoved(WLC_P__NT))
        if (WLC_P__MOVEON_SPIDER .ne. 0) then
            iostr='input/spiders'
            call load_precalc_spiders(iostr,wlc_spiders,wlc_numberOfSpiders)
            wlc_maxNBend = 2000 + 4*get_highestNumberOfLegs(wlc_spiders,wlc_numberOfSpiders)
        else
            wlc_maxNBend = 2000
        endif
        allocate(wlc_bendPoints(wlc_maxNBend))
        wlc_nBend=0

        if (WLC_P__VARIABLE_CHEM_STATE) then
            allocate(wlc_METH(WLC_P__NT)) !Underlying methalation profile
        endif
        !Allocate vector of writhe and elastic energies for replicas
        if (WLC_P__PT_TWIST) then
            allocate(wlc_Wrs(wlc_nLKs))
            allocate(wlc_eelasREPLICAS(wlc_nLKs,4))

        endif
        if (WLC_P__RING) then !TOdo this should be if ("knot")
            wlc_NCross = 0
            wlc_CrossSize = min(10000,max_chain_length())**2
            if (wlc_CrossSize == 10000) then
                print*, "Two many beas for ring calculation"
                stop 1
            endif
            ! In include
            allocate(wlc_Cross(wlc_CrossSize,6))
            allocate(wlc_CrossP(wlc_CrossSize,6))
        endif
        !If parallel tempering is on, initialize the nodeNumbers

        if (WLC_P__PT_TWIST) then

            !Allocate node numbers
            allocate(wlc_nodeNUMBER(wlc_nLKs))
            do i = 1,wlc_nLKs
                wlc_nodeNUMBER(i) = i
            enddo

            !Initially, replica start and replica end are the first and second to last replicas for even
            !nLKs and the first and second to last for odd nLKs
            if (mod(wlc_nLKs,2).eq.0) then
                wlc_replicaSTART = 1
                wlc_replicaEND = wlc_nLKs - 1
            else
                wlc_replicaSTART = 1
                wlc_replicaEND = wlc_nLKs - 2
            endif

            !Allocate the number of replica exchange trials and successes and initialize to zero
            allocate(wlc_nSWAPup(wlc_nLKs))
            allocate(wlc_nSWAPdown(wlc_nLKs))
            allocate(wlc_nTRIALup(wlc_nLKs))
            allocate(wlc_nTRIALdown(wlc_nLKs))

            wlc_nSWAPup = 0
            wlc_nSWAPdown = 0
            wlc_nTRIALup = 0
            wlc_nTRIALdown = 0

        endif

        if (WLC_P__COLLISIONDETECTIONTYPE /= 0) then
            allocate(wlc_coltimes(WLC_P__NT,WLC_P__NT))
            wlc_coltimes = -1.0_dp
        else
            allocate(wlc_coltimes(1,1))
            wlc_coltimes = nan
        endif

#if MPI_VERSION
        ! -----------------------------------------------
        !
        !   Generate thread safe random number seeds
        !
        !--------------------------------------------
        call MPI_COMM_RANK(MPI_COMM_WORLD,wlc_id,error)
        if (wlc_id .eq. 0) then ! head node
            if (.false.) then ! set spedific seed
                Irand = 7171
            else ! seed from clock
                call date_and_time(datedum,timedum,zonedum,seedvalues)
                Irand = int(-seedvalues(5)*1E7-seedvalues(6)*1E5 &
                          -seedvalues(7)*1E3-seedvalues(8))
                Irand = mod(Irand,10000)
                print*, "MPI Random Intiger seed:",Irand
            endif
            call random_setseed(Irand*(wlc_id + 1),wlc_rand_stat) ! random seed for head node
            do dest = 1,wlc_numProcesses-1 ! send out the others
                call MPI_Send (Irand,1, MPI_integer, dest,   0, &
                                MPI_COMM_WORLD,error )
            enddo
        else ! worker node
            source = 0
            call MPI_Recv ( Irand, 1, MPI_integer, source, 0, &
                            MPI_COMM_WORLD, status, error )
            call random_setseed(Irand*(wlc_id + 1),wlc_rand_stat)
            !if (wlc_restart) then
            !    call pt_restart(wlc_p)
            !endif
        endif
#else
        if (.false.) then ! if you wanted to set specific seed
            wlc_rand_seed = 7171
        else ! seed from clock
            call date_and_time(datedum,timedum,zonedum,seedvalues)
            ! funny business
            wlc_rand_seed = int(-seedvalues(5)*1E7-seedvalues(6)*1E5 &
                      -seedvalues(7)*1E3-seedvalues(8))
            wlc_rand_seed = mod(wlc_rand_seed,10000)
            print*, "Random Intiger seed:",wlc_rand_seed
        endif

        call random_setseed(wlc_rand_seed, wlc_rand_stat)
#endif
        if (WLC_P__ELASTICITY_TYPE=='nucleosomes') then
            call loadNucleosomePositions(wlc_nucleosomeWrap,wlc_basepairs)
        endif

        call initcond(wlc_R, wlc_U, WLC_P__NT, &
            WLC_P__NP, WLC_P__FRMFILE, wlc_rand_stat,wlc_p)

        if (WLC_P__EXTERNAL_FIELD_TYPE == 'Random_to_cube_side' .and. WLC_P__APPLY_EXTERNAL_FIELD) then
            call set_external_bindpoints(wlc_rand_stat)
        endif

        if (WLC_P__NEIGHBOR_BINS .and. ((WLC_P__CONFINETYPE == 'excludedShpereInPeriodic') .or. WLC_P__CONFINETYPE == 'none')) then
            do ii=1,WLC_P__NT
                wlc_R_period(1,ii)=modulo(wlc_R(1,ii),WLC_P__LBOX_X)
                wlc_R_period(2,ii)=modulo(wlc_R(2,ii),WLC_P__LBOX_Y)
                wlc_R_period(3,ii)=modulo(wlc_R(3,ii),WLC_P__LBOX_Z)
            enddo
        endif
        if (WLC_P__FIELD_INT_ON) then
            ! initialize a/b sequence
            if (WLC_P__CHEM_STATE_FROM_FILE) then
                iostr='input/ab'
                call wlcsim_params_loadAB(iostr)
            else
                call init_chemical_state(wlc_AB,WLC_P__LAM,WLC_P__FA, WLC_P__ASYMMETRICALTERNATINGCHEM)
            endif

            ! initialize methalation sequence
            if (WLC_P__VARIABLE_CHEM_STATE) then
                if (WLC_P__CHEM_SEQ_FROM_FILE) then
                    print*, "Loding input meth seq..."
                    if (WLC_P__ENSEMBLE_METH .and. wlc_id>0) then
                        write(iostr,"(I4)") wlc_id
                        iostr = adjustL(iostr)
                        iostr = trim(iostr)
                        iostr = "input/meth_"//trim(iostr)
                        iostr = trim(iostr)
                        print*, "reading ",iostr
                        call wlcsim_params_loadMeth(iostr)
                    else
                        iostr='input/meth'
                        call wlcsim_params_loadMeth(iostr)
                    endif
                else
                    call init_chemical_state(wlc_meth,WLC_P__LAM_METH,WLC_P__F_METH,.False.)
                endif
            endif

            ! calculate volumes of bins
            if (WLC_P__CONFINETYPE.eq.'sphere' .and. WLC_P__FRACTIONAL_BIN) then
                call MC_calcVolume(WLC_P__DBIN, &
                                WLC_P__LBOX_X, wlc_Vol, wlc_rand_stat)
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
            setMinXYZ = [0.0_dp,0.0_dp,0.0_dp]  ! location of corner of bin
            setBinShape = [10,10,10]   ! Specify first level of binning
            call constructBin(wlc_bin,setBinShape,setMinXYZ,setBinSize)
            do i=1,WLC_P__NT
                if (WLC_P__NEIGHBOR_BINS .and.&
                    ((WLC_P__CONFINETYPE == 'excludedShpereInPeriodic')&
                    .or. WLC_P__CONFINETYPE == 'none')) then
                    call addBead(wlc_bin,wlc_R_period,WLC_P__NT,i)
                elseif (WLC_P__CONFINETYPE == 'sphere') then
                    call addBead(wlc_bin,wlc_R,WLC_P__NT,i)
                else
                    print*, "Not an option yet.  See params."
                    stop 1
                endif
            enddo
        endif

        ! initialize energies
        !-------------------------------------
        !
        !  Set all energies to zero in case they aren't set later
        !
        !--------------------------------------
        call set_all_energy_to_zero()
        wlc_NPHI = 0  ! NUMBER o phi values that change, i.e. number of bins that were affected

        wlc_time = 0.0_dp
        wlc_time_ind = 0
        wlc_mc_ind = 0

    end subroutine initialize_wlcsim_data

    function pack_as_para(wlc_p) result(para)
        implicit none
        type(wlcsim_params), intent(in) :: wlc_p
        real(dp) para(10)
        para(1) = wlc_p%eb
        para(2) = wlc_p%epar
        para(3) = wlc_p%eperp
        para(4) = wlc_p%GAM
        para(5) = wlc_p%ETA
        para(6) = wlc_p%XIR
        para(7) = wlc_p%XIU
        para(8) = WLC_P__LBOX_X
        para(9) = wlc_p%LHC
        para(10) = wlc_p%VHC
    end function pack_as_para


    subroutine printDescription(wlc_p)
        use energies
        implicit none
        type(wlcsim_params), intent(in) :: wlc_p
        print*, "---------------System Description---------------"
        print*, " type of simulation, codeName", WLC_P__CODENAME
        print*, " WLC, DSSWLC, GC, simType", wlc_p%SIMTYPE
        print*, "Bead variables:"
        print*, " Total number of beads, NT = ", WLC_P__NT
        print*, " Number of beads in a polymer, NB = ", WLC_P__NB
        print*, " Number of monomers in a polymer, nMpP = ", WLC_P__NMPP
        print*, " Number of polymers, NP = ",WLC_P__NP
        print*, " Number of beads in a monomer, nBpM = ", WLC_P__NBPM
        print*, " fraction Methalated", WLC_P__F_METH
        print*, " LAM_METH", WLC_P__LAM_METH
        print*, " "
        print*, "Length and volume Variables:"
        print*, " persistance length =",(WLC_P__L0/(2.0_dp*wlc_p%EPS))
        print*, " length of each polymer in simulation, l = ",WLC_P__L
        print*, " twist persistence length, lt", WLC_P__LT
        print*, " lbox = ", WLC_P__LBOX_X, WLC_P__LBOX_Y, WLC_P__LBOX_Z
        print*, " Number of bins in x direction", &
                   WLC_P__NBIN_X, WLC_P__NBIN_Y,WLC_P__NBIN_Z
        print*, " Number of bins", wlc_p%NBIN
        print*, " spatial descritation dbin = ",WLC_P__DBIN
        print*, " L0 = ", WLC_P__L0
        print*, " GAM = ", wlc_p%GAM
        print*, " bead volume V = ", WLC_P__BEADVOLUME
        print*, " number of kuhn lengths between beads, eps ", wlc_p%EPS
        print*, " "
        print*, "Energy Variables"
        print*, " elasticity EPS =", wlc_p%EPS
        print*, " solvent-polymer CHI =",energyOf(chi_)%cof
        print*, " compression cof, KAP =", energyOf(kap_)%cof
        print*, " field strength, hA =", energyOf(field_)%cof
        print*, " field strength, AEF =", energyOf(external_)%cof
        print*, " two body potential strength, A2B =", energyOf(twoBody_)%cof
        print*, " -energy of binding unmethalated ", WLC_P__EU," more positive for favorable binding"
        print*, " -energy of binding methalated",WLC_P__EM
        print*, " HP1_Binding energy parameter", energyOf(couple_)%cof
        print*, " chemical potential of HP1, mu", energyOf(mu_)%cof
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

    subroutine tweak_param_defaults(wlc_p)
        use polydispersity, only: max_chain_length
        implicit none
        integer ii
        type(wlcsim_params), intent(inout) :: wlc_p

        !  Edit the following to optimize wlc_p performance
        !  Monte-Carlo simulation parameters
        wlc_MCAMP(1) = 0.5_dp*PI
        wlc_MCAMP(2) = 0.3_dp*WLC_P__L0
        wlc_MCAMP(3) = 0.5_dp*PI
        wlc_MCAMP(4) = 0.5_dp*PI
        wlc_MCAMP(5) = 0.5_dp*PI
        wlc_MCAMP(6) = 5.0_dp*WLC_P__L0
        wlc_MCAMP(7) = nan
        wlc_MCAMP(8) = nan
        wlc_MCAMP(9) = nan
        wlc_MCAMP(10) = nan

        ! if we're not using field interactions
        ! energies, then this should never be on
        if ((.not. WLC_P__FIELD_INT_ON) .and. wlc_p%MOVEON(9)/=0) then
            wlc_p%MOVEON(7) = 0  ! Change in Binding state
            print*, "turning off movetype 7, binding, becuase unneeded"
        endif
        wlc_p%MOVEON(8) = 0  ! Chain flip ! TOdo not working
        ! if number of polymers is 1, this should never be on
        if (WLC_P__NP < 2 .and. wlc_p%MOVEON(9)/=0) then
            wlc_p%MOVEON(9) = 0
            print*, "Turning off movetype 9, chain exchange, because <2 polymers"
        endif

        if (isnan(wlc_p%MINWINDOW(1))) wlc_p%MINWINDOW(1) = dble(min(10,max_chain_length()))
        if (isnan(wlc_p%MINWINDOW(2))) wlc_p%MINWINDOW(2) = dble(min(10,max_chain_length()))
        if (isnan(wlc_p%MINWINDOW(3))) wlc_p%MINWINDOW(3) = dble(min(10,max_chain_length()))
        if (isnan(wlc_p%MINWINDOW(7))) wlc_p%MINWINDOW(7) = dble(min(10,max_chain_length()))

        ! Solution
        !WLC_P__LBOX_X = WLC_P__NBIN_X*WLC_P__DBIN
        !WLC_P__LBOX_Y = WLC_P__NBIN_Y*WLC_P__DBIN
        !WLC_P__LBOX_Z = WLC_P__NBIN_Z*WLC_P__DBIN
        wlc_p%NBIN = WLC_P__NBIN_X*WLC_P__NBIN_Y*WLC_P__NBIN_Z

        do ii = 1,nMovetypes
            wlc_window(ii)=wlc_p%MINWINDOW(ii)
        enddo

        if (WLC_P__CODENAME == 'brad') then
            ! initialize windows to number of beads
            wlc_p%MAXWINDOW = real(max_chain_length(),dp)! Max Size of window for bead selection
            wlc_p% MinWindoW  = 1.0_dp         ! Min Size of window for bead selection

            ! Window amplitudes
            wlc_p%MINAMP = 0.0_dp ! minium amplitude
            wlc_p%MINAMP(1) = 0.07_dp*pi
            wlc_p%MINAMP(2) = 0.01_dp*WLC_P__L0
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
            wlc_p%MOVEON(11) = 0 ! SuperReptation
            wlc_p%MOVEON(12) = 0 ! Spider

        endif

        ! If ring is on, turn off the pivot move
        if (WLC_P__RING) then
            wlc_p%MOVEON(3) = 0
        endif

    end subroutine

    subroutine wlcsim_params_recenter()
        use polydispersity, only: first_bead_of_chain, length_of_chain
    !  Prevents drift in periodic BC
        implicit none
        integer IB, I,ii, J   ! Couners
        real(dp) R0(3)  ! Offset to move by
        do I = 1,WLC_P__NP
            IB=first_bead_of_chain(I)
            R0(1) = wlc_R(1,IB) - MODULO(wlc_R(1,IB),WLC_P__LBOX_X)
            R0(2) = wlc_R(2,IB) - MODULO(wlc_R(2,IB),WLC_P__LBOX_Y)
            R0(3) = wlc_R(3,IB) - MODULO(wlc_R(3,IB),WLC_P__LBOX_Z)
            if ( abs(R0(1))+abs(R0(2))+abs(R0(3)) .gt. eps) then
                do J = 1,length_of_chain(I)
                    wlc_R(:,IB) = wlc_R(:,IB)-R0(:)
                    IB = IB + 1
                enddo
            endif
        enddo
    end subroutine

    subroutine printSimInfo(i)
    ! print out current simulation metainformation
        implicit none
        integer, intent(in) :: i
        print*, 'Current time ', wlc_time
        print*, 'Time point ', wlc_time_ind, ' out of ', WLC_P__STEPSPERSAVE*WLC_P__NUMSAVEPOINTS
        print*, 'Save point ', i, ' out of ', WLC_P__NUMSAVEPOINTS
    end subroutine

    subroutine printEnergies()
    ! For realtime feedback on wlc_p simulation
        use energies, only: energyOf, NUMBER_OF_ENERGY_TYPES
        implicit none
        integer ii
        do ii = 1,NUMBER_OF_ENERGY_TYPES
            print*, "Energy of ",energyOf(ii)%name_str,"=",energyOf(ii)%E
        enddo
    end subroutine

    subroutine printEnergyChanges()
        use energies, only: energyOf, NUMBER_OF_ENERGY_TYPES
        implicit none
        integer ii
        do ii = 1,NUMBER_OF_ENERGY_TYPES
            print*, "Change in energy of ",energyOf(ii)%name_str,"=",energyOf(ii)%dE
        enddo
    end subroutine

    subroutine calcTotalPolymerVolume(wlc_p,totalVpoly)
        implicit none
        type(wlcsim_params), intent(in) :: wlc_p
        integer i
        real(dp), intent(out) :: totalVpoly
        real(dp) VV
        totalVpoly=0.0_dp
        VV=WLC_P__DBIN**3
        do I = 1,wlc_p%NBIN
            if (WLC_P__FRACTIONAL_BIN) VV = wlc_Vol(I)
            !if (VV.le.0.1_dp) cycle
            totalVpoly = totalVpoly + VV*(wlc_PHIA(I) + wlc_PHIB(I))
        enddo
        print*, "Total volume of polymer from density", totalVpoly,&
                " and from beads ",WLC_P__NT*WLC_P__BEADVOLUME

    end subroutine

    subroutine wlcsim_params_printPhi(wlc_p)
        use energies, only: energyOf, chi_, kap_, couple_
    ! prints densities for trouble shooting
        implicit none
        type(wlcsim_params), intent(in) :: wlc_p
        integer I
        real(dp) EKap, ECouple, EChi,VV, PHIPOly
        print*,"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        print*, " PHIA  | PHIB  | PPoly |  Vol  | EKap  | EChi  |ECouple|"
        VV=WLC_P__DBIN**3
        do I = 1,wlc_p%NBIN
            if (WLC_P__FRACTIONAL_BIN) VV = wlc_Vol(I)
            if (VV.le.0.1_dp) cycle
            PHIPOLY = wlc_PHIA(I) + wlc_PHIB(I)
            EChi = VV*(energyOf(chi_)%cof/WLC_P__BEADVOLUME)*PHIPoly*(1.0_dp-PHIPoly)
            ECouple = VV*energyOf(couple_)%cof*(wlc_PHIA(I))**2
            if(PHIPoly > 1.0_dp) then
            EKap = VV*(energyOf(kap_)%cof/WLC_P__BEADVOLUME)*(PHIPoly-1.0_dp)**2
            else
            cycle
            EKap = 0.0_dp
            endif
            write(*,"(4f8.4,3f8.1)") wlc_PHIA(I), wlc_PHIB(I), &
                                wlc_PHIA(I) + wlc_PHIB(I),VV,&
                                EKap,EChi,ECouple
        enddo
        print*,"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    end subroutine

    subroutine printWindowStats(wlc_p)
    ! For realtime feedback on adaptation
        implicit none
        type(wlcsim_params), intent(in) :: wlc_p
        integer I ! counter
        I = 0
        print*, "Succes | MCAMP | WindoW| type "
        do I = 1,nMovetypes
            if (wlc_p%MOVEON(i).eq.1) then
                write(*,"(f8.5,2f8.2,1A1,1A20)") wlc_phit(i), wlc_MCAMP(i),&
                    wlc_WindoW(i),' ', moveNames(i)
            endif
        enddo
        return
    end subroutine

    subroutine wlcsim_params_LoadField(wlc_p,fileName)
        implicit none
        type(wlcsim_params), intent(in) :: wlc_p
        integer I
        character(MAXFILENAMELEN) fileName ! file name to load from
        open (unit = inFileUnit, file = fileName, status = 'OLD')
        do I = 1,wlc_p%NBIN
            read(inFileUnit,*) wlc_PHIH(I)
        enddo
        return
    end subroutine

    subroutine load_l2_field(wlc_p)
        implicit none
        type(wlcsim_params), intent(in) :: wlc_p
        integer I
        open (unit = inFileUnit, file = 'input/field_l2', status = 'OLD')
        do I = 1,wlc_p%NBIN
            read(inFileUnit,*) wlc_PHIH_l2(:,I)
        enddo
        close(inFileUnit)
        return
    end subroutine

    subroutine wlcsim_params_MakeField(wlc_p)
        implicit none
        type(wlcsim_params), intent(in) :: wlc_p
        integer indBin  ! index of bin
        integer IX,IY,IZ ! bin corrdinates

        do IX = 1,WLC_P__NBIN_X
            do IY = 1,WLC_P__NBIN_Y
                do IZ = 1,WLC_P__NBIN_Z
                    indBin = IX + &
                        (IY-1)*WLC_P__NBIN_X + &
                        (IZ-1)*WLC_P__NBIN_X*WLC_P__NBIN_Y
                    wlc_PHIH(indBin) = dsin(WLC_P__K_FIELD*WLC_P__DBIN*dble(IX))
                enddo
            enddo
        enddo
        return
    end subroutine

    subroutine wlcsim_params_loadAB(fileName)
    ! Loads AB for file...has not been tested
        use polydispersity, only: length_of_chain
        implicit none
        character(MAXFILENAMELEN), intent(in) :: fileName ! file name to load from
        integer IB, I, J ! counters
        open (unit = inFileUnit, file = fileName, status = 'OLD')
        IB = 1
        do I = 1,WLC_P__NP
            do J = 1,length_of_chain(I)
                read(inFileUnit,"(I2)") wlc_AB(IB)
                IB = IB + 1
            enddo
        enddo
        close(inFileUnit)
    end subroutine
    subroutine wlcsim_params_loadMeth(fileName)
    ! Loads Methalation for file...has not been tested
        use polydispersity, only: length_of_chain
        implicit none
        character(MAXFILENAMELEN), intent(in) :: fileName ! file name to load from
        integer IB, I, J ! counters
        open (unit = inFileUnit, file = fileName, status = 'OLD')
        IB = 1
        do I = 1,WLC_P__NP
            do J = 1,length_of_chain(I)
                read(inFileUnit,"(I2)") wlc_meth(IB)
                IB = IB + 1
            enddo
        enddo
        close(inFileUnit)
    end subroutine

    subroutine wlcsim_params_saveR(fileName,repeatingBC)
    ! Writes R and AB to file for analysis
    ! Rx  Ry  Rz AB
        use polydispersity, only: length_of_chain
        implicit none
        logical, intent(in) :: repeatingBC  ! 1 for reapeating boundary conditions
        integer I,J,IB  ! counters
        character(MAXFILENAMELEN), intent(in) :: fileName
        character(MAXFILENAMELEN) fullName
        LOGICAL isfile
        fullName = trim(fileName) // trim(wlc_repSuffix)
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
              do J = 1,length_of_chain(I)
                 if (WLC_P__SAVEAB) then
                    write(outFileUnit,"(3f10.3,I2)") &
                          modulo(wlc_R(1,IB),WLC_P__LBOX_X)&
                         ,modulo(wlc_R(2,IB),WLC_P__LBOX_Y)&
                         ,modulo(wlc_R(3,IB),WLC_P__LBOX_Z)&
                         ,wlc_AB(IB)
                 else
                    write(outFileUnit,"(3f10.3)") &
                          modulo(wlc_R(1,IB),WLC_P__LBOX_X)&
                         ,modulo(wlc_R(2,IB),WLC_P__LBOX_Y)&
                         ,modulo(wlc_R(3,IB),WLC_P__LBOX_Z)
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
              do J = 1,length_of_chain(I)
                  if (WLC_P__SAVEAB) then
                     if (WLC_P__CHANGINGCHEMICALIDENTITY) then
                         call print_11char_vec(outFileUnit, wlc_R(:,IB), .False.)
                         write(outFileUnit,"(2I3)") wlc_AB(IB), wlc_METH(IB)
                     else
                         call print_11char_vec(outFileUnit, wlc_R(:,IB), .False.)
                         write(outFileUnit,"(I2)") wlc_AB(IB)
                     endif
                  else
                      call print_11char_vec(outFileUnit, wlc_R(:,IB), .True.)
                  endif
                  IB = IB + 1
              enddo
           enddo
        endif
        close(outFileUnit)

      end subroutine wlcsim_params_saveR

    subroutine wlcsim_params_savePHI(wlc_p,fileName)
        use energies, only: energyOf, maierSaupe_
    ! Saves PHIA and PHIB to file for analysis
        implicit none
        integer I  ! counters
        type(wlcsim_params), intent(in) :: wlc_p
        character(MAXFILENAMELEN), intent(in) :: fileName
        character(MAXFILENAMELEN) fullName
        fullName=  trim(fileName) // trim(wlc_repSuffix)
        open (unit = outFileUnit, file = fullName, status = 'NEW')
        if (energyOf(maierSaupe_)%isOn) then
            do I = 1,wlc_p%NBIN
                write(outFileUnit,"(7f7.2)") wlc_PHIA(I),wlc_PHIB(I),wlc_PHI_l2(:,I)
            enddo
        else
            do I = 1,wlc_p%NBIN
                write(outFileUnit,"(2f7.2)") wlc_PHIA(I),wlc_PHIB(I)
            enddo
        endif

        close(outFileUnit)
    end subroutine

    subroutine wlcsim_params_saveU(fileName,stat)
    ! Saves U to ASCII file for analisys
        use polydispersity, only: length_of_chain
        implicit none
        integer I,J,IB  ! counters
        character(MAXFILENAMELEN), intent(in) :: fileName
        character(MAXFILENAMELEN) fullName
        character(len = *), intent(in) :: stat
        fullName=  trim(fileName) // trim(wlc_repSuffix)
        open (unit = outFileUnit, file = fullName, status = stat)
        IB = 1
        do I = 1,WLC_P__NP
            do J = 1,length_of_chain(I)
                if (WLC_P__LOCAL_TWIST) then
                    call print_11char_vec(outFileUnit, wlc_U(:,IB), .False.)
                    call print_11char_vec(outFileUnit, wlc_V(:,IB), .True.)
                else
                    call print_11char_vec(outFileUnit, wlc_U(:,IB), .True.)
                endif
                IB = IB + 1
            enddo
        enddo
        close(outFileUnit)
    end subroutine

    subroutine save_parameters(wlc_p,fileName)
        use energies, only: energyOf, chi_, couple_
        ! Write a number of parameters ASCII variables to file for reccords
        implicit none
        type(wlcsim_params), intent(in) :: wlc_p
        character(len=*), intent(in) :: fileName
        open (unit =outFileUnit, file = fileName, status = 'NEW')
            write(outFileUnit,"(I8)") WLC_P__NT ! 1 Number of beads in simulation
            write(outFileUnit,"(I8)") WLC_P__NMPP  ! 2 Number of monomers in a polymer
            write(outFileUnit,"(I8)") WLC_P__NB ! 3 Number of beads in a polymer
            write(outFileUnit,"(I8)") WLC_P__NP ! 4 Number of polymers in simulation
            write(outFileUnit,"(I8)") WLC_P__NT ! 5 Number of beads in simulation
            write(outFileUnit,"(I8)") WLC_P__NBPM  ! 6 Number of beads per monomer

            write(outFileUnit,"(f10.5)") WLC_P__L0    ! Equilibrium segment length
            write(outFileUnit,"(f10.5)") energyOf(chi_)%cof  ! 8  initail CHI parameter value
            write(outFileUnit,"(f10.5)") WLC_P__LBOX_X  ! 10 Lenth of box
            write(outFileUnit,"(f10.5)") WLC_P__EU    ! Energy unmethalated
            write(outFileUnit,"(f10.5)") WLC_P__EM    ! 12 Energy methalated
            write(outFileUnit,"(f10.5)") energyOf(couple_)%cof ! Energy of HP1 binding
            write(outFileUnit,"(f10.5)") (WLC_P__L0/wlc_p%EPS) ! 14 Khun lenth
            write(outFileUnit,"(A)") "-999"  ! for historic reasons
            write(outFileUnit,"(f10.5)") WLC_P__F_METH  ! methalation fraction
            write(outFileUnit,"(f10.5)") WLC_P__LAM_METH  ! methalation lambda
        close(outFileUnit)
    end subroutine

    subroutine print_11char_vec(outFileUnit,vec, nextLine)
        implicit none
        integer, intent(in) :: outFileUnit
        real(dp), intent(in), dimension(3) :: vec
        logical nextLine
        integer ii

        do ii = 1,3
            call print_11char_float(outFileUnit,vec(ii))
        enddo
        if (nextLine) then
            write(outFileUnit, "(a)") ""
        endif

    end subroutine
    subroutine print_11char_float(outFileUnit,x)
        implicit none
        integer, intent(in) :: outFileUnit
        real(dp), intent(in) :: x
        if (x > 999999999.9 .or. x < -99999999.9) then
            write(outFileUnit, "(E11.3)") x
        elseif (x > 9999999.99 .or. x < -999999.99) then
            write(outFileUnit, "(F11.1)", advance="no") x
        elseif (x > 999999.999 .or. x < -99999.999) then
            write(outFileUnit, "(F11.2)", advance="no") x
        elseif (x > 99999.9999 .or. x < -9999.9999) then
            write(outFileUnit, "(F11.3)", advance="no") x
        elseif (x > 9999.99999 .or. x < -999.99999) then
            write(outFileUnit, "(F11.4)", advance="no") x
        elseif (x > 999.999999 .or. x < -99.999999) then
            write(outFileUnit, "(F11.5)", advance="no") x
        elseif (x > 99.9999999 .or. x < -9.9999999) then
            write(outFileUnit, "(F11.6)", advance="no") x
        else
            write(outFileUnit, "(F11.7)", advance="no") x
        endif
    end subroutine

    subroutine wlcsim_params_appendEnergyData(save_ind, fileName)
        use energies, only: energyOf, NUMBER_OF_ENERGY_TYPES, kap_
    ! print Energy data
        implicit none
        integer, intent(in) :: save_ind
        character(MAXFILENAMELEN), intent(in) :: fileName
        LOGICAL isfile
        character(MAXFILENAMELEN) fullName
        integer ii
        fullName=  trim(fileName) // trim(wlc_repSuffix)
        inquire(file = fullName, exist = isfile)
        if (isfile) then
            open (unit = outFileUnit, file = fullName, status ='OLD', POSITION = "append")
        else
            open (unit = outFileUnit, file = fullName, status = 'new')
            write(outFileUnit,"(10A)",advance="no") "ind | id |"
            do ii = 1, NUMBER_OF_ENERGY_TYPES
                write(outFileUnit,"(12A)",advance="no")  " E-",energyOf(ii)%name_str, " "
                write(outFileUnit,"(12A)",advance="no")  " x-",energyOf(ii)%name_str, " "
                write(outFileUnit,"(12A)",advance="no")  " c-",energyOf(ii)%name_str, " "
            enddo
            write(outFileUnit,*) " "
        endif
        write(outFileUnit,"(2I5)",advance="no") save_ind, wlc_id
        do ii = 1, NUMBER_OF_ENERGY_TYPES
            call print_11char_float(outFileUnit, energyOf(ii)%E)
            write(outFileUnit,"(A)",advance="no") " "
            call print_11char_float(outFileUnit, energyOf(ii)%x)
            write(outFileUnit,"(A)",advance="no") " "
            call print_11char_float(outFileUnit, energyOf(ii)%cof)
            write(outFileUnit,"(A)",advance="no") " "
        enddo
        write(outFileUnit,*) " "
        close(outFileUnit)
    end subroutine

    subroutine wlcsim_params_appendAdaptData(save_ind, wlc_p, fileName)
    ! Appends wlc_p move adaptation data to the file
        implicit none
        type(wlcsim_params), intent(in) :: wlc_p
        integer, intent(in) :: save_ind
        LOGICAL isfile
        character(MAXFILENAMELEN), intent(in) :: fileName
        character(MAXFILENAMELEN) fullName
        fullName=  trim(fileName) // trim(wlc_repSuffix)
        inquire(file = fullName, exist = isfile)
        if (isfile) then
            open (unit = outFileUnit, file = fullName, status ='OLD', POSITION = "append")
        else
            open (unit = outFileUnit, file = fullName, status = 'new')
            write(outFileUnit,*) "        ",&
                       " ----- Crank ------    |  -----  Slide -----   |",&
                       " ----- Pivot ------    |  -----  rotate ----   |",&
                       " full Chain Rotation   |  full Chain Slide     |",&
                       "    chem-move  |  end-end flip |",&
                       "  chian swap   |  reptation    |super reptation|",&
                       "    spider     |"
            write(outFileUnit,*) "ind| id|",&
                       " Win 1 | AMP 1 | SUC 1 | Win 2 | AMP 2 | SUC 2 |",&
                       " Win 3 | AMP 3 | SUC 3 | ON  4 | AMP 4 | SUC 4 |",&
                       " ON  5 | AMP 5 | SUC 5 | ON  6 | AMP 6 | SUC 6 |",&
                       " ON  7 | SUC 7 | ON  8 | SUC 8 |", &
                       " ON  9 | SUC 9 | ON 10 | SUC 10| ON 11 | SUC11 |",&
                       " AMP12 | SUC12 |"
        endif
        write(outFileUnit,"(2I4,30f8.3)") save_ind,wlc_id,&
            real(wlc_WindoW(1)),wlc_MCAMP(1),wlc_PHIT(1), &
            real(wlc_WindoW(2)),wlc_MCAMP(2),wlc_PHIT(2), &
            real(wlc_WindoW(3)),wlc_MCAMP(3),wlc_PHIT(3), &
            real(wlc_p%MOVEON(4)),wlc_MCAMP(4),wlc_PHIT(4), &
            real(wlc_p%MOVEON(5)),wlc_MCAMP(5),wlc_PHIT(5), &
            real(wlc_p%MOVEON(6)),wlc_MCAMP(6),wlc_PHIT(6), &
            real(wlc_p%MOVEON(7)),wlc_PHIT(7), &
            real(wlc_p%MOVEON(8)),wlc_PHIT(8), &
            real(wlc_p%MOVEON(9)),wlc_PHIT(9), &
            real(wlc_p%MOVEON(10)),wlc_PHIT(10), &
            real(wlc_p%MOVEON(11)),wlc_PHIT(11), &
            wlc_MCAMP(12),wlc_PHIT(12)
        close(outFileUnit)
    end subroutine
    subroutine wlcsim_params_writebinary(wlc_p,baseName)
    !    This function writes the contence of the structures wlc_p and
    !  to a binary file.  if you add more variables to  you need to
    !  a seperate write command for them as it is not possible to write
    !  a structure with allocatables to a binar file.
    !    The contence are stored in
    !     baseName//'R'
    !     baseName//'U'
    !     etc.
        implicit none
        integer sizeOftype         ! for binary saving
        type(wlcsim_params), intent(in) :: wlc_p             ! to be save or filled
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

        sizeOftype = int(SIZEOF(wlc_R))
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
        write(outFileUnit,rec = 1) wlc_R
        close(outFileUnit)

        ! -------- U --------

        sizeOftype = int(SIZEOF(wlc_U))
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
        write(outFileUnit,rec = 1) wlc_U
        close(outFileUnit)

        ! -------- AB --------

        sizeOftype = int(SIZEOF(wlc_AB))
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
        write(outFileUnit,rec = 1) wlc_AB
        close(outFileUnit)

        ! -------- Vol --------

        sizeOftype = int(SIZEOF(wlc_Vol))
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
        write(outFileUnit,rec = 1) wlc_Vol
        close(outFileUnit)
    end subroutine

    subroutine wlcsim_params_readBinary(wlc_p, baseName)
    ! This function reads what wlcsim_params_writebinary writes and
    ! stores it to wlc_p and .  Be sure to allocate  before
    ! calling this command.
        implicit none
        integer sizeOftype         ! for binary saving
        type(wlcsim_params) wlc_p             ! to be save or filled
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

        sizeOftype = int(SIZEOF(wlc_R))
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
        read(inFileUnit,rec = 1) wlc_R
        close(inFileUnit)

        ! -------- U --------

        sizeOftype = int(SIZEOF(wlc_U))
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
        read(inFileUnit,rec = 1) wlc_U
        close(inFileUnit)

        ! -------- AB --------

        sizeOftype = int(SIZEOF(wlc_AB))
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
        read(inFileUnit,rec = 1) wlc_AB
        close(inFileUnit)

        ! -------- Vol --------

        sizeOftype = int(SIZEOF(wlc_Vol))
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
        read(inFileUnit,rec = 1) wlc_Vol
        close(inFileUnit)
    end subroutine

    subroutine save_simulation_state(save_ind, wlc_p, outfile_base, stat)
        implicit none
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
        call wlcsim_params_appendEnergyData(save_ind, filename)

        !part 2.5 - adaptations
        filename = trim(adjustL(outfile_base)) // 'adaptations'
        call wlcsim_params_appendAdaptData(save_ind, wlc_p, filename)

        if (WLC_P__SAVEPHI) then
            write(filename,num2strFormatString) save_ind
            filename = trim(adjustL(outfile_base)) // 'phi' // trim(adjustL(filename))
            call wlcsim_params_savePHI(wlc_p,filename)
        endif

        if (WLC_P__SAVER) then
            if (WLC_P__APPEND_R) then
                filename = trim(adjustL(outfile_base)) // 'rAll'
            else
                write(filename,num2strFormatString) save_ind
                filename = trim(adjustL(outfile_base)) // 'r' // trim(adjustL(filename))
            endif
            call wlcsim_params_saveR(filename,.false.)
        endif

        if (WLC_P__SAVEU) then
            write(filename,num2strFormatString) save_ind
            filename = trim(adjustL(outfile_base)) // 'u' // trim(adjustL(filename))
            call wlcsim_params_saveU(filename,stat)
        endif

        if (WLC_P__COLLISIONDETECTIONTYPE /= 0) then
            filename = trim(adjustL(outfile_base)) // 'coltimes'
            open(unit = outFileUnit, file = filename, status = 'REPLACE')
            do ind = 1,WLC_P__NT
                write(outFileUnit,*) (wlc_coltimes(ind,j), j = 1,WLC_P__NT)
            enddo
            close(outFileUnit)
        endif
    end subroutine save_simulation_state

    !Get Lks for parallel tempering from file
    subroutine get_LKs_from_file()
    implicit none
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

    wlc_nLKs = nLKs
    allocate(wlc_LKs(nLks))

    open(unit = 1, file = 'input/LKs')
    do i = 1, nLks
        read(unit = 1,fmt = *) wlc_Lks(i)
    enddo
    close(unit = 1)
    end subroutine get_LKs_from_file

    subroutine get_renormalized_chain_params(wlc_p)
    use MC_wlc, only: calc_elastic_constants
    !     Setup the parameters for the simulation
    !
    !     1. Determine the simulation type
    !     2. Evaluate the polymer elastic parameters
    !     3. Determine the parameters for Brownian dynamics simulation
    implicit none
    type(wlcsim_params), intent(inout) :: wlc_p

    ! calculate metrics that don't change between WLC, ssWLC, GC
    if (WLC_P__ELASTICITY_TYPE == "constant") then
        wlc_p%DEL = WLC_P__L0/WLC_P__LP

        call calc_elastic_constants(wlc_p%DEL,WLC_P__LP,WLC_P__LT,&
                                    wlc_p%eb,wlc_p%epar, &
                                    wlc_p%GAM,wlc_p%XIR,wlc_p%eperp,wlc_p%ETA, &
                                    wlc_p%XIU,wlc_p%DT, &
                                    wlc_p%SIGMA,wlc_p%etwist,wlc_p%simtype)
    elseif (WLC_P__ELASTICITY_TYPE == "nucleosomes") then
        wlc_p%DEL = 0.0 ! not used
        wlc_p%GAM = 0.0 ! not used
        wlc_p%XIR = 0.0 ! not used
        wlc_p%ETA = 0.0 ! not used
        wlc_p%DT = 0.0 ! not used
        wlc_p%SIGMA = 0.0 ! not used
        wlc_p%eb = 0.0 ! not used
        wlc_p%epar = 0.0 ! not used
        wlc_p%eperp = 0.0 ! not used
        wlc_p%etwist = 0.0 ! not used
        wlc_p%XIU = 0.0 ! not used
        wlc_p%DT = 0.0 ! not used
    endif
    end subroutine get_renormalized_chain_params
end module params
