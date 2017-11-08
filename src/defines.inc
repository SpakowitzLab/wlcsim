! VARIABLE COMMENT: ! which simulation code to run
! DEFAULTS COMMENT: ! This input is based on Quinn's code
! TYPE: CHARACTER(MAXPARAMLEN)
#define WLC_P__CODENAME "quinn"
! VARIABLE COMMENT: ! Number of beads in a polymer NB = nMpP*nBpM
! DEFAULTS COMMENT: ! From Chr 16 at 200bp per nucleosome
! TYPE: INTEGER
#define WLC_P__NB 393216
! VARIABLE COMMENT: ! Number of monomers (NOT BEADS!) in a polymer
! DEFAULTS COMMENT: None
! TYPE: INTEGER
#define WLC_P__NMPP WLC_P__NB/WLC_P__NBPM
! VARIABLE COMMENT: ! number beads per monomer
! DEFAULTS COMMENT: Each bead is a nucleosome
! TYPE: INTEGER
#define WLC_P__NBPM 1
! VARIABLE COMMENT: ! Number of polymers
! DEFAULTS COMMENT: ! sigle chromosome
! TYPE: INTEGER
#define WLC_P__NP 1
! VARIABLE COMMENT: ! length of each polymer in simulation
! DEFAULTS COMMENT: ! chr 16 / Delta
! TYPE: REAL(DP)
#define WLC_P__L 231716.571_dp
! VARIABLE COMMENT: ! persistence length
! DEFAULTS COMMENT: ! 50nm / Delta
! TYPE: REAL(DP)
#define WLC_P__LP 1.786_DP
! VARIABLE COMMENT: ! twist persistence length
! DEFAULTS COMMENT: ! twist persistence length equals persistence length by default
! TYPE: REAL(DP)
#define WLC_P__LT 1
! VARIABLE COMMENT: ! Bead volume
! DEFAULTS COMMENT: ! From nucleosome volume estimate of 520 nm^3
! TYPE: REAL(DP)
#define WLC_P__BEADVOLUME 0.0237_DP
! VARIABLE COMMENT: ! Fraction of A beads
! DEFAULTS COMMENT: ! Initial fraction methelated
! TYPE: REAL(DP)
#define WLC_P__FA 0.0_DP
! VARIABLE COMMENT: ! Chemical correlation parameter (eigenvalue of transition matrix that generates A/B's)
! DEFAULTS COMMENT: ! perfectly random sequence  (see generating_sequences.rst for details)
! TYPE: REAL(DP)
#define WLC_P__LAM 0.0_DP
! VARIABLE COMMENT: ! radius triggering collisions to be recorded in "coltimes"
! DEFAULTS COMMENT: ! consistant with bead volume
! TYPE: REAL(DP)
#define WLC_P__COLLISIONRADIUS 0.178_DP
! VARIABLE COMMENT: ! wave vector of applied sinusoidal field (used in PT to step around 1st order phase transition)
! DEFAULTS COMMENT: ! No field
! TYPE: REAL(DP)
#define WLC_P__K_FIELD 0.0_DP
! VARIABLE COMMENT: ! Energy of HP1 binding for unmethalated sites
! DEFAULTS COMMENT: ! From Pete's paper
! TYPE: REAL(DP)
#define WLC_P__EU -1.52_DP
! VARIABLE COMMENT: ! Energy of HP1 binding for methalated sites
! DEFAULTS COMMENT: ! by default, no hp1 binding energy included
! TYPE: REAL(DP)
#define WLC_P__EM 0.01_DP
! VARIABLE COMMENT: ! Fraction methalated is using option 2 ( for initialization )
! DEFAULTS COMMENT: ! not used
! TYPE: REAL(DP)
#define WLC_P__F_METH 0.0_DP
! VARIABLE COMMENT: ! eigenvalue of transition matrix generating initial methalation
! DEFAULTS COMMENT: ! not used
! TYPE: REAL(DP)
#define WLC_P__LAM_METH 0.0_DP
! VARIABLE COMMENT: ! Number of MC bins on an edge
! DEFAULTS COMMENT: ! set number density of nucleus
! TYPE: INTEGER
#define WLC_P__NBINX_X 64
! VARIABLE COMMENT: ! Number of MC bins on an edge
! DEFAULTS COMMENT: ! same
! TYPE: INTEGER
#define WLC_P__NBINX_Y 64
! VARIABLE COMMENT: ! Number of MC bins on an edge
! DEFAULTS COMMENT: ! same
! TYPE: INTEGER
#define WLC_P__NBINX_Z 64
! VARIABLE COMMENT: ! Number of collision-detection bins on each edge
! DEFAULTS COMMENT: ! not used.  equivalent to collisionDetectionType = 1
! TYPE: INTEGER
#define WLC_P__NCOLBIN 1
! VARIABLE COMMENT: ! monte carlo field bin total box length (approximate)
! DEFAULTS COMMENT: ! box size/confinment, *MUST* be set by user
! TYPE: REAL(DP)
#define WLC_P__LBOX_X 64.0_dp
! VARIABLE COMMENT: ! monte carlo field bin total box length (approximate)
! DEFAULTS COMMENT: ! box size/confinment, *MUST* be set by user
! TYPE: REAL(DP)
#define WLC_P__LBOX_Y 64.0_dp
! VARIABLE COMMENT: ! monte carlo field bin total box length (approximate)
! DEFAULTS COMMENT: ! box size/confinment, *MUST* be set by user
! TYPE: REAL(DP)
#define WLC_P__LBOX_Z 64.0_dp
! VARIABLE COMMENT: ! monte carlo field bin discretization size (approximate)
! DEFAULTS COMMENT: ! should always be 1 according to quinn
! TYPE: REAL(DP)
#define WLC_P__DBIN 1
! VARIABLE COMMENT: None
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)
#define WLC_P__CONFINEMENTPARAMETER_1 NAN
! VARIABLE COMMENT: None
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)
#define WLC_P__CONFINEMENTPARAMETER_2 NAN
! VARIABLE COMMENT: ! desired hit rate
! DEFAULTS COMMENT: ! Target
! TYPE: REAL(DP)
#define WLC_P__PDESIRE_CRANK_SHAFT 0.5_DP
! VARIABLE COMMENT: ! desired hit rate
! DEFAULTS COMMENT: ! Target
! TYPE: REAL(DP)
#define WLC_P__PDESIRE_SLIDE_MOVE 0.5_DP
! VARIABLE COMMENT: ! desired hit rate
! DEFAULTS COMMENT: ! Target
! TYPE: REAL(DP)
#define WLC_P__PDESIRE_PIVOT_MOVE 0.5_DP
! VARIABLE COMMENT: ! desired hit rate
! DEFAULTS COMMENT: ! Target
! TYPE: REAL(DP)
#define WLC_P__PDESIRE_ROTATE_MOVE 0.5_DP
! VARIABLE COMMENT: ! desired hit rate
! DEFAULTS COMMENT: ! Target
! TYPE: REAL(DP)
#define WLC_P__PDESIRE_FULL_CHAIN_ROTATION 0.5_DP
! VARIABLE COMMENT: ! desired hit rate
! DEFAULTS COMMENT: ! Target
! TYPE: REAL(DP)
#define WLC_P__PDESIRE_FULL_CHAIN_SLIDE 0.5_DP
! VARIABLE COMMENT: ! desired hit rate
! DEFAULTS COMMENT: ! Target
! TYPE: REAL(DP)
#define WLC_P__PDESIRE_CHANGE_BINDING_STATE 0.5_DP
! VARIABLE COMMENT: ! desired hit rate
! DEFAULTS COMMENT: ! Target
! TYPE: REAL(DP)
#define WLC_P__PDESIRE_CHAIN_FLIP 0.5_DP
! VARIABLE COMMENT: ! desired hit rate
! DEFAULTS COMMENT: ! Target
! TYPE: REAL(DP)
#define WLC_P__PDESIRE_CHAIN_EXCHANGE 0.5_DP
! VARIABLE COMMENT: ! desired hit rate
! DEFAULTS COMMENT: ! Target
! TYPE: REAL(DP)
#define WLC_P__PDESIRE_REPTATION 0.5_DP
! VARIABLE COMMENT: ! desired hit rate
! DEFAULTS COMMENT: ! Target
! TYPE: REAL(DP)
#define WLC_P__PDESIRE_SUPER_REPTATION 0.5_DP
! VARIABLE COMMENT: ! Max Size of window for bead selection
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)
#define WLC_P__MAXWINDOW_CRANK_SHAFT DBLE(MIN(150,WLC_P__NB))
! VARIABLE COMMENT: ! Max Size of window for bead selection
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)
#define WLC_P__MAXWINDOW_SLIDE_MOVE DBLE(MIN(150,WLC_P__NB))
! VARIABLE COMMENT: ! Max Size of window for bead selection
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)
#define WLC_P__MAXWINDOW_PIVOT_MOVE DBLE(MIN(150,INT(WLC_P__NB/2.0)))
! VARIABLE COMMENT: ! Max Size of window for bead selection
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)
#define WLC_P__MAXWINDOW_ROTATE_MOVE NAN
! VARIABLE COMMENT: ! Max Size of window for bead selection
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)
#define WLC_P__MAXWINDOW_FULL_CHAIN_ROTATION NAN
! VARIABLE COMMENT: ! Max Size of window for bead selection
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)
#define WLC_P__MAXWINDOW_FULL_CHAIN_SLIDE NAN
! VARIABLE COMMENT: ! Max Size of window for bead selection
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)
#define WLC_P__MAXWINDOW_CHANGE_BINDING_STATE DBLE(MIN(4,WLC_P__NB))
! VARIABLE COMMENT: ! Max Size of window for bead selection
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)
#define WLC_P__MAXWINDOW_CHAIN_FLIP NAN
! VARIABLE COMMENT: ! Max Size of window for bead selection
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)
#define WLC_P__MAXWINDOW_CHAIN_EXCHANGE NAN
! VARIABLE COMMENT: ! Max Size of window for bead selection
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)
#define WLC_P__MAXWINDOW_REPTATION NAN
! VARIABLE COMMENT: ! Max Size of window for bead selection
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)
#define WLC_P__MAXWINDOW_SUPER_REPTATION NAN
! VARIABLE COMMENT: ! Min Size of window for bead selection
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)
#define WLC_P__MINWINDOW_CRANK_SHAFT 30
! VARIABLE COMMENT: ! Min Size of window for bead selection
! DEFAULTS COMMENT: ! If this is too small equllibration will be very slow
! TYPE: REAL(DP)
#define WLC_P__MINWINDOW_SLIDE_MOVE 10
! VARIABLE COMMENT: ! Min Size of window for bead selection
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)
#define WLC_P__MINWINDOW_PIVOT_MOVE 50
! VARIABLE COMMENT: ! Min Size of window for bead selection
! DEFAULTS COMMENT: ! If this is too small equllibration will be very slow
! TYPE: REAL(DP)
#define WLC_P__MINWINDOW_ROTATE_MOVE NAN
! VARIABLE COMMENT: ! Min Size of window for bead selection
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)
#define WLC_P__MINWINDOW_FULL_CHAIN_ROTATION NAN
! VARIABLE COMMENT: ! Min Size of window for bead selection
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)
#define WLC_P__MINWINDOW_FULL_CHAIN_SLIDE NAN
! VARIABLE COMMENT: ! Min Size of window for bead selection
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)
#define WLC_P__MINWINDOW_CHANGE_BINDING_STATE 1
! VARIABLE COMMENT: ! Min Size of window for bead selection
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)
#define WLC_P__MINWINDOW_CHAIN_FLIP NAN
! VARIABLE COMMENT: ! Min Size of window for bead selection
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)
#define WLC_P__MINWINDOW_CHAIN_EXCHANGE NAN
! VARIABLE COMMENT: ! Min Size of window for bead selection
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)
#define WLC_P__MINWINDOW_REPTATION NAN
! VARIABLE COMMENT: ! Min Size of window for bead selection
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)
#define WLC_P__MINWINDOW_SUPER_REPTATION NAN
! VARIABLE COMMENT: ! minium amplitude
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)
#define WLC_P__MINAMP_CRANK_SHAFT 0.1_DP*PI
! VARIABLE COMMENT: ! minium amplitude
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)
#define WLC_P__MINAMP_SLIDE_MOVE 0.2_DP*wlc_p%L0
! VARIABLE COMMENT: ! minium amplitude
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)
#define WLC_P__MINAMP_PIVOT_MOVE 0.2_DP*PI
! VARIABLE COMMENT: ! minium amplitude
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)
#define WLC_P__MINAMP_ROTATE_MOVE 0.2_DP*PI
! VARIABLE COMMENT: ! minium amplitude
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)
#define WLC_P__MINAMP_FULL_CHAIN_ROTATION 0.05_DP*PI
! VARIABLE COMMENT: ! minium amplitude
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)
#define WLC_P__MINAMP_FULL_CHAIN_SLIDE 0.2_DP*wlc_p%L0
! VARIABLE COMMENT: ! minium amplitude
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)
#define WLC_P__MINAMP_CHANGE_BINDING_STATE NAN
! VARIABLE COMMENT: ! minium amplitude
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)
#define WLC_P__MINAMP_CHAIN_FLIP NAN
! VARIABLE COMMENT: ! minium amplitude
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)
#define WLC_P__MINAMP_CHAIN_EXCHANGE NAN
! VARIABLE COMMENT: ! minium amplitude
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)
#define WLC_P__MINAMP_REPTATION NAN
! VARIABLE COMMENT: ! minium amplitude
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)
#define WLC_P__MINAMP_SUPER_REPTATION NAN
! VARIABLE COMMENT: ! maximum amplitude
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)
#define WLC_P__MAXAMP_CRANK_SHAFT 1.0_DP*PI
! VARIABLE COMMENT: ! maximum amplitude
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)
#define WLC_P__MAXAMP_SLIDE_MOVE 1.0_DP*wlc_p%L0
! VARIABLE COMMENT: ! maximum amplitude
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)
#define WLC_P__MAXAMP_PIVOT_MOVE 1.0_DP*PI
! VARIABLE COMMENT: ! maximum amplitude
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)
#define WLC_P__MAXAMP_ROTATE_MOVE 1.0_DP*PI
! VARIABLE COMMENT: ! maximum amplitude
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)
#define WLC_P__MAXAMP_FULL_CHAIN_ROTATION 1.0_DP*PI
! VARIABLE COMMENT: ! maximum amplitude
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)
#define WLC_P__MAXAMP_FULL_CHAIN_SLIDE 0.1*wlc_p%LBOX(1)
! VARIABLE COMMENT: ! maximum amplitude
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)
#define WLC_P__MAXAMP_CHANGE_BINDING_STATE NAN
! VARIABLE COMMENT: ! maximum amplitude
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)
#define WLC_P__MAXAMP_CHAIN_FLIP NAN
! VARIABLE COMMENT: ! maximum amplitude
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)
#define WLC_P__MAXAMP_CHAIN_EXCHANGE NAN
! VARIABLE COMMENT: ! maximum amplitude
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)
#define WLC_P__MAXAMP_REPTATION NAN
! VARIABLE COMMENT: ! maximum amplitude
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)
#define WLC_P__MAXAMP_SUPER_REPTATION NAN
! VARIABLE COMMENT: ! Is the move active
! DEFAULTS COMMENT: ! crank-shaft move
! TYPE: INTEGER
#define WLC_P__MOVEON_CRANK_SHAFT 1
! VARIABLE COMMENT: ! Is the move active
! DEFAULTS COMMENT: ! slide move
! TYPE: INTEGER
#define WLC_P__MOVEON_SLIDE_MOVE 1
! VARIABLE COMMENT: ! pivot move
! DEFAULTS COMMENT: ! need to be able to move end of chain
! TYPE: INTEGER
#define WLC_P__MOVEON_PIVOT_MOVE 1
! VARIABLE COMMENT: ! Is the move active
! DEFAULTS COMMENT: ! rotate move
! TYPE: INTEGER
#define WLC_P__MOVEON_ROTATE_MOVE 1
! VARIABLE COMMENT: ! full chain rotation
! DEFAULTS COMMENT: ! DNA too long to bother
! TYPE: INTEGER
#define WLC_P__MOVEON_FULL_CHAIN_ROTATION 0
! VARIABLE COMMENT: ! full chain slide
! DEFAULTS COMMENT: ! DNA too long to bother
! TYPE: INTEGER
#define WLC_P__MOVEON_FULL_CHAIN_SLIDE 0
! VARIABLE COMMENT: ! Changes A/B type
! DEFAULTS COMMENT: ! HP1 binding state
! TYPE: INTEGER
#define WLC_P__MOVEON_CHANGE_BINDING_STATE 1
! VARIABLE COMMENT: ! Is the move active
! DEFAULTS COMMENT: ! Chain flip ! TOdo not working
! TYPE: INTEGER
#define WLC_P__MOVEON_CHAIN_FLIP 0
! VARIABLE COMMENT: ! Chain exchange
! DEFAULTS COMMENT: ! DNA too long / single chr
! TYPE: INTEGER
#define WLC_P__MOVEON_CHAIN_EXCHANGE 0
! VARIABLE COMMENT: ! Reptation
! DEFAULTS COMMENT: ! Chain too long to bother
! TYPE: INTEGER
#define WLC_P__MOVEON_REPTATION 0
! VARIABLE COMMENT: ! Is the move active
! DEFAULTS COMMENT: ! Super Reptation
! TYPE: INTEGER
#define WLC_P__MOVEON_SUPER_REPTATION 0
! VARIABLE COMMENT: ! target for ratio of window to anmplitude
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)
#define WLC_P__WINTARGET_CRANK_SHAFT 8.0_DP
! VARIABLE COMMENT: ! target for ratio of window to anmplitude
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)
#define WLC_P__WINTARGET_SLIDE_MOVE 8.0_DP
! VARIABLE COMMENT: ! target for ratio of window to anmplitude
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)
#define WLC_P__WINTARGET_PIVOT_MOVE 8.0_DP
! VARIABLE COMMENT: ! target for ratio of window to anmplitude
! DEFAULTS COMMENT: ! window of 1
! TYPE: REAL(DP)
#define WLC_P__WINTARGET_ROTATE_MOVE NAN
! VARIABLE COMMENT: ! target for ratio of window to anmplitude
! DEFAULTS COMMENT: ! No window
! TYPE: REAL(DP)
#define WLC_P__WINTARGET_FULL_CHAIN_ROTATION NAN
! VARIABLE COMMENT: ! target for ratio of window to anmplitude
! DEFAULTS COMMENT: ! No window
! TYPE: REAL(DP)
#define WLC_P__WINTARGET_FULL_CHAIN_SLIDE NAN
! VARIABLE COMMENT: ! target for ratio of window to anmplitude
! DEFAULTS COMMENT: ! no amplitude
! TYPE: REAL(DP)
#define WLC_P__WINTARGET_CHANGE_BINDING_STATE NAN
! VARIABLE COMMENT: ! target for ratio of window to anmplitude
! DEFAULTS COMMENT: ! no window or amplitude
! TYPE: REAL(DP)
#define WLC_P__WINTARGET_CHAIN_FLIP NAN
! VARIABLE COMMENT: ! target for ratio of window to anmplitude
! DEFAULTS COMMENT: ! No window or amplitude
! TYPE: REAL(DP)
#define WLC_P__WINTARGET_CHAIN_EXCHANGE NAN
! VARIABLE COMMENT: ! target for ratio of window to anmplitude
! DEFAULTS COMMENT: ! window or amplitiuge
! TYPE: REAL(DP)
#define WLC_P__WINTARGET_REPTATION NAN
! VARIABLE COMMENT: ! target for ratio of window to anmplitude
! DEFAULTS COMMENT: ! no window
! TYPE: REAL(DP)
#define WLC_P__WINTARGET_SUPER_REPTATION NAN
! VARIABLE COMMENT: ! Nunber of steps between adapt
! DEFAULTS COMMENT: ! adapt after at most 1000 steps
! TYPE: INTEGER
#define WLC_P__NADAPT_CRANK_SHAFT 1000
! VARIABLE COMMENT: ! Nunber of steps between adapt
! DEFAULTS COMMENT: ! adapt after at most 1000 steps
! TYPE: INTEGER
#define WLC_P__NADAPT_SLIDE_MOVE 1000
! VARIABLE COMMENT: ! Nunber of steps between adapt
! DEFAULTS COMMENT: ! adapt after at most 1000 steps
! TYPE: INTEGER
#define WLC_P__NADAPT_PIVOT_MOVE 1000
! VARIABLE COMMENT: ! Nunber of steps between adapt
! DEFAULTS COMMENT: ! adapt after at most 1000 steps
! TYPE: INTEGER
#define WLC_P__NADAPT_ROTATE_MOVE 1000
! VARIABLE COMMENT: ! Nunber of steps between adapt
! DEFAULTS COMMENT: ! adapt after at most 1000 steps
! TYPE: INTEGER
#define WLC_P__NADAPT_FULL_CHAIN_ROTATION 1000
! VARIABLE COMMENT: ! Nunber of steps between adapt
! DEFAULTS COMMENT: ! adapt after at most 1000 steps
! TYPE: INTEGER
#define WLC_P__NADAPT_FULL_CHAIN_SLIDE 1000
! VARIABLE COMMENT: ! Nunber of steps between adapt
! DEFAULTS COMMENT: ! adapt after at most 1000 steps
! TYPE: INTEGER
#define WLC_P__NADAPT_CHANGE_BINDING_STATE 1000
! VARIABLE COMMENT: ! Nunber of steps between adapt
! DEFAULTS COMMENT: ! adapt after at most 1000 steps
! TYPE: INTEGER
#define WLC_P__NADAPT_CHAIN_FLIP 1000
! VARIABLE COMMENT: ! Nunber of steps between adapt
! DEFAULTS COMMENT: ! adapt after at most 1000 steps
! TYPE: INTEGER
#define WLC_P__NADAPT_CHAIN_EXCHANGE 1000
! VARIABLE COMMENT: ! Nunber of steps between adapt
! DEFAULTS COMMENT: ! adapt after at most 1000 steps
! TYPE: INTEGER
#define WLC_P__NADAPT_REPTATION 1000
! VARIABLE COMMENT: ! Nunber of steps between adapt
! DEFAULTS COMMENT: ! adapt after at most 1000 steps
! TYPE: INTEGER
#define WLC_P__NADAPT_SUPER_REPTATION 1000
! VARIABLE COMMENT: ! threshold for deciding to usually not use a move
! DEFAULTS COMMENT: ! if a move succeeds < 5% of the time, start using it only every reduce_move cycles
! TYPE: REAL(DP)
#define WLC_P__MIN_ACCEPT 0.05_DP
! VARIABLE COMMENT: ! whether or not to stop usuing a move when it goes below min_accept success
! DEFAULTS COMMENT: ! use moves that fall below the min_accept threshold only once every ~ times they would otherwise be used, set to one for no effect
! TYPE: INTEGER
#define WLC_P__REDUCE_MOVE 1
! VARIABLE COMMENT: ! distributionof segment size in crankshaft move (unif = 0, exp = 1)
! DEFAULTS COMMENT: ! exponential fragment sizes mix better
! TYPE: INTEGER
#define WLC_P__WINTYPE 1
! VARIABLE COMMENT: ! number of steps between parallel tempering
! DEFAULTS COMMENT: ! 100 steps between parallel tempering is pretty frequent
! TYPE: INTEGER
#define WLC_P__STEPSPEREXCHANGE 100
! VARIABLE COMMENT: ! read teh variable name
! DEFAULTS COMMENT: ! make this large <e.g. 2000>
! TYPE: INTEGER
#define WLC_P__NREPLICAEXCHANGEPERSAVEPOINT 10000
! VARIABLE COMMENT: ! total number of save points
! DEFAULTS COMMENT: ! 200 total save points, i.e. 2000 steps per save point
! TYPE: INTEGER
#define WLC_P__NUMSAVEPOINTS 1000
! VARIABLE COMMENT: ! number of simulation steps to take, BD
! DEFAULTS COMMENT: ! This is MC simulation not BD
! TYPE: INTEGER
#define WLC_P__STEPSPERSAVE 0
! VARIABLE COMMENT: ! save points before turning on NNoInt
! DEFAULTS COMMENT: ! number of simulation steps before turning on interactions in Quinn's wlc_p scheduler
! TYPE: INTEGER
#define WLC_P__NNOINT 5
! VARIABLE COMMENT: ! when to turn KAP energy on
! DEFAULTS COMMENT: ! turn on compression energy immediately
! TYPE: INTEGER
#define WLC_P__N_KAP_ON 5
! VARIABLE COMMENT: ! when to turn CHI energy on
! DEFAULTS COMMENT: ! turn on chi energy immediately
! TYPE: INTEGER
#define WLC_P__N_CHI_ON 7
! VARIABLE COMMENT: ! when to turn CHI energy on
! DEFAULTS COMMENT: ! turn on chi energy immediately
! TYPE: INTEGER
#define WLC_P__N_CHI_L2_ON 7
! VARIABLE COMMENT: ! number of mc steps before starting BD
! DEFAULTS COMMENT: ! number of initilizing mc steps. 1000s x num polymers is good
! TYPE: INTEGER
#define WLC_P__NINITMCSTEPS 4000
! VARIABLE COMMENT: ! how many times each move is done per step
! DEFAULTS COMMENT: ! crank shaft <default 30>
! TYPE: INTEGER
#define WLC_P__MOVESPERSTEP_CRANK_SHAFT 30
! VARIABLE COMMENT: ! how many times each move is done per step
! DEFAULTS COMMENT: ! slide <default 60>
! TYPE: INTEGER
#define WLC_P__MOVESPERSTEP_SLIDE_MOVE 60
! VARIABLE COMMENT: ! how many times each move is done per step
! DEFAULTS COMMENT: ! pivot <default 30>
! TYPE: INTEGER
#define WLC_P__MOVESPERSTEP_PIVOT_MOVE 1
! VARIABLE COMMENT: ! how many times each move is done per step
! DEFAULTS COMMENT: ! rotate <default 60>
! TYPE: INTEGER
#define WLC_P__MOVESPERSTEP_ROTATE_MOVE 60
! VARIABLE COMMENT: ! how many times each move is done per step
! DEFAULTS COMMENT: ! full chain rotation <default 5>
! TYPE: INTEGER
#define WLC_P__MOVESPERSTEP_FULL_CHAIN_ROTATION 5
! VARIABLE COMMENT: ! how many times each move is done per step
! DEFAULTS COMMENT: ! full chain slide <default 6>
! TYPE: INTEGER
#define WLC_P__MOVESPERSTEP_FULL_CHAIN_SLIDE 6
! VARIABLE COMMENT: ! how many times each move is done per step
! DEFAULTS COMMENT: ! bind < default 10>
! TYPE: INTEGER
#define WLC_P__MOVESPERSTEP_CHANGE_BINDING_STATE 10
! VARIABLE COMMENT: ! how many times each move is done per step
! DEFAULTS COMMENT: ! filp 
! TYPE: INTEGER
#define WLC_P__MOVESPERSTEP_CHAIN_FLIP 0
! VARIABLE COMMENT: ! how many times each move is done per step
! DEFAULTS COMMENT: ! full chain swap < default 5>
! TYPE: INTEGER
#define WLC_P__MOVESPERSTEP_CHAIN_EXCHANGE 5
! VARIABLE COMMENT: ! how many times each move is done per step
! DEFAULTS COMMENT: ! reptation <default 10>
! TYPE: INTEGER
#define WLC_P__MOVESPERSTEP_REPTATION 10
! VARIABLE COMMENT: ! how many times each move is done per step
! DEFAULTS COMMENT: ! super-reptation <default 57>
! TYPE: INTEGER
#define WLC_P__MOVESPERSTEP_SUPER_REPTATION 57
! VARIABLE COMMENT: ! whether the polymer is a ring
! DEFAULTS COMMENT: ! not a ring by default
! TYPE: LOGICAL
#define WLC_P__RING .FALSE.
! VARIABLE COMMENT: ! whether to include twist (wlc_p only for now)
! DEFAULTS COMMENT: ! don't include twist by default
! TYPE: LOGICAL
#define WLC_P__TWIST .FALSE.
! VARIABLE COMMENT: ! type of Boundary Conditions
! DEFAULTS COMMENT: ! spherical nucleus
! TYPE: CHARACTER(MAXPARAMLEN)
#define WLC_P__CONFINETYPE 'sphere'
! VARIABLE COMMENT: ! initial condition type
! DEFAULTS COMMENT: ! because it is easy
! TYPE: CHARACTER(MAXPARAMLEN)
#define WLC_P__INITCONDTYPE 'randomLineSphereBoundary'
! VARIABLE COMMENT: ! field-based self interactions on
! DEFAULTS COMMENT: ! Crossings allowed
! TYPE: LOGICAL
#define WLC_P__INTRAPOLYMER_STICK_CROSSING_ENFORCED .FALSE.
! VARIABLE COMMENT: ! read initial chemical sequence from file
! DEFAULTS COMMENT: ! don't read binding state from file
! TYPE: LOGICAL
#define WLC_P__CHEM_STATE_FROM_FILE .FALSE.
! VARIABLE COMMENT: ! read initial chemical/methylation state from file
! DEFAULTS COMMENT: ! read in methelation state
! TYPE: LOGICAL
#define WLC_P__CHEM_SEQ_FROM_FILE .TRUE.
! VARIABLE COMMENT: ! read initial condition R from file
! DEFAULTS COMMENT: ! don't load initial bead positions from file
! TYPE: LOGICAL
#define WLC_P__FRMFILE .FALSE.
! VARIABLE COMMENT: ! read initial field from file
! DEFAULTS COMMENT: ! don't load initial field values from file
! TYPE: LOGICAL
#define WLC_P__FRMFIELD .FALSE.
! VARIABLE COMMENT: ! save first passage time vectors to file
! DEFAULTS COMMENT: ! don't track first passage time collisions between beads
! TYPE: INTEGER
#define WLC_P__COLLISIONDETECTIONTYPE 0
! VARIABLE COMMENT: ! stop sim with coltimes is full
! DEFAULTS COMMENT: ! stop sim when coltimes is full
! TYPE: LOGICAL
#define WLC_P__EXITWHENCOLLIDED .FALSE.
! VARIABLE COMMENT: ! save R vectors to file
! DEFAULTS COMMENT: ! don't save orientation vectors (saves disk space)
! TYPE: LOGICAL
#define WLC_P__SAVER .FALSE.
! VARIABLE COMMENT: ! save U vectors to file
! DEFAULTS COMMENT: ! don't save orientation vectors (makes restart of ssWLC possible)
! TYPE: LOGICAL
#define WLC_P__SAVEU .FALSE.
! VARIABLE COMMENT: ! save AB (chemical identity) to file
! DEFAULTS COMMENT: ! save AB sequence
! TYPE: LOGICAL
#define WLC_P__SAVEAB .TRUE.
! VARIABLE COMMENT: ! save Phi vectors to file
! DEFAULTS COMMENT: ! don't save A/B density per bin (not needed for restart)
! TYPE: LOGICAL
#define WLC_P__SAVEPHI .FALSE.
! VARIABLE COMMENT: ! whether we are restarting from a previous sim or not
! DEFAULTS COMMENT: ! don't restart from previously saved simulation
! TYPE: LOGICAL
#define WLC_P__RESTART .FALSE.
! VARIABLE COMMENT: ! whether to have inter bead lennard jones energies
! DEFAULTS COMMENT: ! no intrapolymer interactions by default
! TYPE: LOGICAL
#define WLC_P__INTERP_BEAD_LENNARD_JONES .FALSE.
! VARIABLE COMMENT: ! Allows chemical A/B sequence to change
! DEFAULTS COMMENT: ! Meth state needs to change with time
! TYPE: LOGICAL
#define WLC_P__CHANGINGCHEMICALIDENTITY .TRUE.
! VARIABLE COMMENT: ! initialize asymmetric alturnating chain
! DEFAULTS COMMENT: ! Nope
! TYPE: LOGICAL
#define WLC_P__ASYMMETRICALTERNATINGCHEM .FALSE.
! VARIABLE COMMENT: ! Choose the Hamiltonian you would like
! DEFAULTS COMMENT: ! See MC_Hamiltonian
! TYPE: CHARACTER(MAXPARAMLEN)
#define WLC_P__FIELDINTERACTIONTYPE 'chromatin'
! VARIABLE COMMENT: ! whether or not to parallel temper
! DEFAULTS COMMENT: ! use parallel if applicable
! TYPE: LOGICAL
#define WLC_P__PTON .FALSE.
! VARIABLE COMMENT: None
! DEFAULTS COMMENT: ! don't parallel temper linking number (twist) by default
! TYPE: LOGICAL
#define WLC_P__PT_TWIST .FALSE.
! VARIABLE COMMENT: None
! DEFAULTS COMMENT: ! don't parallel temper chi by default
! TYPE: LOGICAL
#define WLC_P__PT_CHI .FALSE.
! VARIABLE COMMENT: None
! DEFAULTS COMMENT: ! don't parallel temper h by default
! TYPE: LOGICAL
#define WLC_P__PT_H .FALSE.
! VARIABLE COMMENT: None
! DEFAULTS COMMENT: ! don't parallel temper kap by default
! TYPE: LOGICAL
#define WLC_P__PT_KAP .FALSE.
! VARIABLE COMMENT: None
! DEFAULTS COMMENT: ! don't parallel temper mu by default
! TYPE: LOGICAL
#define WLC_P__PT_MU .FALSE.
! VARIABLE COMMENT: None
! DEFAULTS COMMENT: ! don't parallel temper HP1 binding by default
! TYPE: LOGICAL
#define WLC_P__PT_COUPLE .TRUE.
! VARIABLE COMMENT: None
! DEFAULTS COMMENT: ! don't parallel temper maier Saupe by default
! TYPE: LOGICAL
#define WLC_P__PT_MAIERSAUPE .FALSE.
! VARIABLE COMMENT: ! number of exchange attemts between adapt
! DEFAULTS COMMENT: ! 1000 exchange attempts between adaptations
! TYPE: INTEGER
#define WLC_P__NREPADAPT 1000
! VARIABLE COMMENT: ! when to decrese cof spacing
! DEFAULTS COMMENT: ! TOdo: enter justification for these defaults, if any.
! TYPE: REAL(DP)
#define WLC_P__LOWERREPEXE 0.09_DP
! VARIABLE COMMENT: ! when to increase cof spacing
! DEFAULTS COMMENT: ! TOdo: fine if the only justification is "these just work"
! TYPE: REAL(DP)
#define WLC_P__UPPERREPEXE 0.18_DP
! VARIABLE COMMENT: ! minumum acceptable Cof
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)
#define WLC_P__LOWERCOFRAIL 0.005_DP
! VARIABLE COMMENT: ! maximum acceptable Cof
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)
#define WLC_P__UPPERCOFRAIL 0.1_DP
! VARIABLE COMMENT: Allow time for the chains to equilibrate without interaction
! DEFAULTS COMMENT: None
! TYPE: INTEGER
#define WLC_P__INDSTARTREPADAPT 100000000
! VARIABLE COMMENT: Stop changing values to allow for data collection
! DEFAULTS COMMENT: None
! TYPE: INTEGER
#define WLC_P__INDENDREPADAPT 100000000
! VARIABLE COMMENT: ! for annealing
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)
#define WLC_P__REPANNEALSPEED 0.00_DP
! VARIABLE COMMENT: None
! DEFAULTS COMMENT: None
! TYPE: LOGICAL
#define WLC_P__REPLICABOUNDS .TRUE.
! VARIABLE COMMENT: ! The upper end of the initial values you which to PT
! DEFAULTS COMMENT: 
! TYPE: REAL(DP)
#define WLC_P__INITIAL_MAX_S 5.0_DP
! VARIABLE COMMENT: ! Does model contain changing chemical state
! DEFAULTS COMMENT: ! Yes, binding state changes
! TYPE: LOGICAL
#define WLD_P__VARIABLE_CHEM_STATE .TRUE.
