! ----------------------------
!
!   General Switches
!
! ---------------------------

#define WLC_P__CODENAME "bruno_bd"
! VARIABLE COMMENT: ! which simulation code to run
! DEFAULTS COMMENT: ! not bruno, brad, or quinn, so will error unless specified elsewehre
! TYPE: CHARACTER(MAXPARAMLEN)

#define WLC_P__RING .FALSE.
! VARIABLE COMMENT: ! whether the polymer is a ring
! DEFAULTS COMMENT: ! not a ring by default
! TYPE: LOGICAL

#define WLC_P__TWIST .FALSE.
! VARIABLE COMMENT: ! whether to include twist (wlc_p only for now)
! DEFAULTS COMMENT: ! don't include twist by default
! TYPE: LOGICAL

#define WLC_P__CONFINETYPE 'NONE'
! VARIABLE COMMENT: ! type of Boundary Conditions
! DEFAULTS COMMENT: ! NONE for no confinement
! TYPE: CHARACTER(MAXPARAMLEN)

#define WLC_P__INITCONDTYPE 'randomWalkWithBoundary'
! VARIABLE COMMENT: ! initial condition type
! DEFAULTS COMMENT: ! r(0)=0, r(i+1) = r(i)+rnorm(); respect boundary condition
! TYPE: CHARACTER(MAXPARAMLEN)

#define WLC_P__INTRAPOLYMER_STICK_CROSSING_ENFORCED .FALSE.
! VARIABLE COMMENT: ! field-based self interactions on
! DEFAULTS COMMENT: ! not yet implemented, so better be false
! TYPE: LOGICAL

#define WLC_P__CHEM_STATE_FROM_FILE .FALSE.
! VARIABLE COMMENT: ! read initial chemical sequence from file
! DEFAULTS COMMENT: ! don't load initial a/b states from file
! TYPE: LOGICAL

#define WLC_P__CHEM_SEQ_FROM_FILE .FALSE.
! VARIABLE COMMENT: ! read initial chemical/methylation state from file
! DEFAULTS COMMENT: ! don't load from file
! TYPE: LOGICAL

#define WLC_P__FRMFILE .FALSE.
! VARIABLE COMMENT: ! read initial condition R from file
! DEFAULTS COMMENT: ! don't load initial bead positions from file
! TYPE: LOGICAL

#define WLC_P__FRMFIELD .FALSE.
! VARIABLE COMMENT: ! read initial field from file
! DEFAULTS COMMENT: ! don't load initial field values from file
! TYPE: LOGICAL

#define WLC_P__COLLISIONDETECTIONTYPE 0
! VARIABLE COMMENT: ! save first passage time vectors to file
! DEFAULTS COMMENT: ! don't track first passage time collisions between beads
! TYPE: INTEGER

#define WLC_P__RESTART .FALSE.
! VARIABLE COMMENT: ! whether we are restarting from a previous sim or not
! DEFAULTS COMMENT: ! don't restart from previously saved simulation
! TYPE: LOGICAL

#define WLC_P__INTERP_BEAD_LENNARD_JONES .FALSE.
! VARIABLE COMMENT: ! whether to have inter bead lennard jones energies
! DEFAULTS COMMENT: ! no intrapolymer interactions by default
! TYPE: LOGICAL

#define WLC_P__CHANGINGCHEMICALIDENTITY .FALSE.
! VARIABLE COMMENT: ! whether many types of bead exist
! DEFAULTS COMMENT: ! Only one type of bead
! TYPE: LOGICAL

#define WLC_P__VARIABLE_CHEM_STATE WLC_P__CHANGINGCHEMICALIDENTITY
! VARIABLE COMMENT: !TODO figure out what subtle difference if any r.t.  WLC_P__CHANGINGCHEMICALIDENTITY
! DEFAULTS COMMENT: !TODO
! TYPE: LOGICAL

#define WLC_P__FIELD_INT_ON .FALSE.
! VARIABLE COMMENT: whether or not to execute mc_int
! DEFAULTS COMMENT: ! no interactions by default
! TYPE: LOGICAL

#define WLC_P__FIELDINTERACTIONTYPE 'NONE'
! VARIABLE COMMENT: None
! DEFAULTS COMMENT: ! See MC_Hamiltonian
! TYPE: CHARACTER(MAXPARAMLEN)

#define WLC_P__RECENTER_ON .FALSE.
! VARIABLE COMMENT: ! Wheter to move drifted chains back in periodic BC
! DEFAULTS COMMENT: ! Not in periodic BC
! TYPE: LOGICAL

#define WLC_P__CYLINDRICAL_CHAIN_EXCLUSION .FALSE.
! VARIABLE COMMENT: ! Chains cannont occupy the same place at same time
! DEFAULTS COMMENT: ! We care about collisions
! TYPE: LOGICAL

#define WLC_P__NEIGHBOR_BINS .FALSE.
! VARIABLE COMMENT: ! Keep track of neighbors by recursive binning
! DEFAULTS COMMENT: ! Set off because expensive
! TYPE: LOGICAL

#define WLC_P__APPLY_EXTERNAL_FIELD .FALSE.
! VARIABLE COMMENT: ! Externally applied field (e.g. radial binding)
! DEFAULTS COMMENT: !
! TYPE: LOGICAL

! ----------------------------
!
!   Output
!
! ---------------------------

#define WLC_P__EXITWHENCOLLIDED .FALSE.
! VARIABLE COMMENT: ! stop sim after every n^2 loops have formed at least once
! DEFAULTS COMMENT: ! don't stop sim when coltimes is full by default
! TYPE: LOGICAL

#define WLC_P__SAVER .TRUE.
! VARIABLE COMMENT: ! save R vectors to file
! DEFAULTS COMMENT: ! do save orientation vectors (makes restart of ssWLC possible)
! TYPE: LOGICAL

#define WLC_P__SAVEU .TRUE.
! VARIABLE COMMENT: ! save U vectors to file
! DEFAULTS COMMENT: ! do save orientation vectors (makes restart of ssWLC possible)
! TYPE: LOGICAL

#define WLC_P__SAVEAB .FALSE.
! VARIABLE COMMENT: ! save AB (chemical identity) to file
! DEFAULTS COMMENT: ! dont' save AB by default, almost nobody uses this
! TYPE: LOGICAL

#define WLC_P__SAVEPHI .FALSE.
! VARIABLE COMMENT: ! save Phi vectors to file
! DEFAULTS COMMENT: ! don't save A/B density per bin (not needed for restart)
! TYPE: LOGICAL

#define WLC_P__APPEND_R .FALSE.
! VARIABLE COMMENT: ! Output all r(t) to same file
! DEFAULTS COMMENT: ! Use old output format
! TYPE: LOGICAL

! ----------------------------
!
!   Number / Chemical identity variables
!
! ---------------------------

#define WLC_P__NB 201
! VARIABLE COMMENT: ! Number of beads in a polymer NB = nMpP*nBpM
! DEFAULTS COMMENT: ! 201 beads per polymer
! TYPE: INTEGER

#define WLC_P__NMPP WLC_P__NB/WLC_P__NBPM
! VARIABLE COMMENT: ! Number of monomers (NOT BEADS!) in a polymer
! DEFAULTS COMMENT: None
! TYPE: INTEGER

#define WLC_P__NBPM 10
! VARIABLE COMMENT: ! number beads per monomer
! DEFAULTS COMMENT: None
! TYPE: INTEGER

#define WLC_P__NP 1
! VARIABLE COMMENT: ! Number of polymers
! DEFAULTS COMMENT: ! one polymer
! TYPE: INTEGER

#define WLC_P__FA 0.5_DP
! VARIABLE COMMENT: ! Fraction of A beads
! DEFAULTS COMMENT: ! half A, half B by default
! TYPE: REAL(DP)

#define WLC_P__LAM 0.0_DP
! VARIABLE COMMENT: ! Chemical correlation parameter (eigenvalue of transition matrix that generates A/B's)
! DEFAULTS COMMENT: ! perfectly random sequence  (see generating_sequences.rst for details)
! TYPE: REAL(DP)

#define WLC_P__F_METH 0.5_DP
! VARIABLE COMMENT: ! Fraction methalated is using option 2 ( for initialization )
! DEFAULTS COMMENT: ! half beads methylated by default
! TYPE: REAL(DP)

#define WLC_P__LAM_METH 0.9_DP
! VARIABLE COMMENT: ! eigenvalue of transition matrix generating initial methalation
! DEFAULTS COMMENT: ! highly alternating sequence by default
! TYPE: REAL(DP)

#define WLC_P__ASYMMETRICALTERNATINGCHEM .FALSE.
! VARIABLE COMMENT: None
! DEFAULTS COMMENT: None
! TYPE: LOGICAL

#define WLC_P__TWO_TAIL .FALSE.
! VARIABLE COMMENT: ! Chromatin with two histon tails per bead
! DEFAULTS COMMENT: ! Specific to Chromatin
! TYPE: LOGICAL


! ----------------------------
!
!   Length variables
!
! ---------------------------

#define WLC_P__L 199
! VARIABLE COMMENT: ! length of each polymer in simulation
! DEFAULTS COMMENT: ! set this yourself!
! TYPE: REAL(DP)

#define WLC_P__LP 1
! VARIABLE COMMENT: ! persistence length
! DEFAULTS COMMENT: ! units of lp by default
! TYPE: REAL(DP)

#define WLC_P__LT 1
! VARIABLE COMMENT: ! twist persistence length
! DEFAULTS COMMENT: ! twist persistence length equals persistence length by default
! TYPE: REAL(DP)

#define WLC_P__BEADVOLUME 0.1_DP
! VARIABLE COMMENT: ! Bead volume
! DEFAULTS COMMENT: ! much smaller than space between beads
! TYPE: REAL(DP)

#define WLC_P__COLLISIONRADIUS 0
! VARIABLE COMMENT: ! radius triggering collisions to be recorded in "coltimes"
! DEFAULTS COMMENT: ! never collide except on floating-point coincidence
! TYPE: REAL(DP)

#define WLC_P__K_FIELD 0.0_DP
! VARIABLE COMMENT: ! wave vector of applied sinusoidal field (used in PT to step around 1st order phase transition)
! DEFAULTS COMMENT: ! some previous values: !1.5708_dp !0.3145_dp
! TYPE: REAL(DP)

#define WLC_P__NCOLBIN 1
! VARIABLE COMMENT: ! Number of collision-detection bins on each edge
! DEFAULTS COMMENT: ! equivalent to collisionDetectionType = 1
! TYPE: INTEGER

#define WLC_P__NBINX 20
#define WLC_P__NBINX_X WLC_P__NBINX
! VARIABLE COMMENT: ! Number of MC bins on an edge
! DEFAULTS COMMENT: ! set this yourself!
! TYPE: INTEGER

#define WLC_P__NBINX_Y WLC_P__NBINX
! VARIABLE COMMENT: ! Number of MC bins on an edge
! DEFAULTS COMMENT: ! set this yourself!
! TYPE: INTEGER

#define WLC_P__NBINX_Z WLC_P__NBINX
! VARIABLE COMMENT: ! Number of MC bins on an edge
! DEFAULTS COMMENT: ! set this yourself!
! TYPE: INTEGER

#define WLC_P__LBOX 20.0_dp
#define WLC_P__LBOX_X WLC_P__LBOX
! VARIABLE COMMENT: ! monte carlo field bin total box length (approximate)
! DEFAULTS COMMENT: ! box size/conËfinment, *MUST* be set by user
! TYPE: REAL(DP)

#define WLC_P__LBOX_Y WLC_P__LBOX
! VARIABLE COMMENT: ! monte carlo field bin total box length (approximate)
! DEFAULTS COMMENT: ! box size/confinment, *MUST* be set by user
! TYPE: REAL(DP)

#define WLC_P__LBOX_Z WLC_P__LBOX
! VARIABLE COMMENT: ! monte carlo field bin total box length (approximate)
! DEFAULTS COMMENT: ! box size/confinment, *MUST* be set by user
! TYPE: REAL(DP)

#define WLC_P__DBIN 1.0_dp
! VARIABLE COMMENT: ! monte carlo field bin discretization size (approximate)
! DEFAULTS COMMENT: ! should always be 1.0_dp according to quinn
! TYPE: REAL(DP)

#define WLC_P__CONFINEMENT_CUBE_LENGTH NAN
! VARIABLE COMMENT: ! Size of confinement. Set = LBOX - 2*DBIN if bin interaction on
! DEFAULTS COMMENT: ! Must set if in use
! TYPE: REAL(DP)

#define WLC_P__CONFINEMENT_SLIT_WIDTH NAN
! VARIABLE COMMENT: !  Size of confinement. Set = LBOX -2*DBIN if bin interation is on
! DEFAULTS COMMENT: ! Must set if in use
! TYPE: REAL(DP)

#define WLC_P__CONFINEMENT_SPHERE_DIAMETER 10.0_dp
! VARIABLE COMMENT: ! Size of confinement. Set = LBOX - 2*DBIN if bin interaction on
! DEFAULTS COMMENT: ! Must set if in use
! TYPE: REAL(DP)

#define WLC_P__CONFINEMENT_ECOLI_LENGTH 0.0_dp
! VARIABLE COMMENT: !
! DEFAULTS COMMENT: ! zero for sphere
! TYPE: REAL(DP)

#define WLC_P__CONFINEMENT_ECOLI_DIAMETER NAN
! VARIABLE COMMENT: ! Size of confinement. Set = LBOX - 2*DBIN if bin interaction on
! DEFAULTS COMMENT: ! Must set if in use
! TYPE: REAL(DP)

#define WLC_P__BINDING_R 0.0_DP
! VARIABLE COMMENT: ! Binding Radius
! DEFAULTS COMMENT: ! No binding
! TYPE: REAL(DP)

#define WLC_P__CHAIN_D 0.0_DP
! VARIABLE COMMENT: Diameter of chain for exclusion
! DEFAULTS COMMENT: Set to zero for no exclusion
! TYPE: REAL(DP)


! ----------------------------
!
!   Energy variables
!
! ---------------------------

#define WLC_P__EU 0.0_DP
! VARIABLE COMMENT: ! Energy of HP1 binding for unmethalated sites
! DEFAULTS COMMENT: ! a function of coarse graining. This should be set by hand if needed.
! TYPE: REAL(DP)

#define WLC_P__EM 0.0_DP
! VARIABLE COMMENT: ! Energy of HP1 binding for methalated sites
! DEFAULTS COMMENT: ! by default, no hp1 binding energy included
! TYPE: REAL(DP)

#define WLC_P__CHI 0.0_DP
! VARIABLE COMMENT: ! Chi paremter
! DEFAULTS COMMENT: ! Not in simple Hamiltonian
! TYPE: REAL(DP)

#define WLC_p__MU 0.0_DP
! VARIABLE COMMENT: ! HP1 chemical potential
! DEFAULTS COMMENT: ! Not in simple Hamiltonian
! TYPE: REAL(DP)

#define WLC_P__HA 0.0_DP
! VARIABLE COMMENT: ! Applied Field
! DEFAULTS COMMENT: ! no applied field
! TYPE: REAL(DP)

#define WLC_P__COUPLING_ENERGY 0.0_dp
! VARIABLE COMMENT: ! Coupling energy between beads (negative = attractive)
! DEFAULTS COMMENT: ! no coupling energy
! TYPE: REAL(DP)

#define WLC_P__HP1_BIND 0.0_DP
! VARIABLE COMMENT: ! Strength of Coupling
! DEFAULTS COMMENT: ! Assumed unwanted
! TYPE: REAL(DP)

#define WLC_P__KAP 0.0_DP
! VARIABLE COMMENT: ! Degree if incompressibility
! DEFAULTS COMMENT: ! Assumed unwanted
! TYPE: REAL(DP)

#define WLC_P__CHI_L2 0.0_DP
! VARIABLE COMMENT: ! prop. to . Maier Saupe parameter
! DEFAULTS COMMENT: ! zero for simplitity
! TYPE: REAL(DP)



! ----------------------------
!
!   Time / steps variables
!
! ---------------------------

#define WLC_P__STEPSPEREXCHANGE 100
! VARIABLE COMMENT: ! number of steps between parallel tempering
! DEFAULTS COMMENT: ! 100 steps between parallel tempering is pretty frequent
! TYPE: INTEGER

#define WLC_P__NREPLICAEXCHANGEPERSAVEPOINT 1000
! VARIABLE COMMENT: ! read teh variable name
! DEFAULTS COMMENT: ! make this large
! TYPE: INTEGER

#define WLC_P__NUMSAVEPOINTS 200
! VARIABLE COMMENT: ! total number of save points
! DEFAULTS COMMENT: ! 200 total save points, i.e. 2000 steps per save point
! TYPE: INTEGER

#define WLC_P__STEPSPERSAVE 2000
! VARIABLE COMMENT: ! steps per save point
! DEFAULTS COMMENT: ! number of simulation steps to take
! TYPE: INTEGER

#define WLC_P__NNOINT 0
! VARIABLE COMMENT: ! save points before turning on NNoInt
! DEFAULTS COMMENT: ! number of simulation steps before turning on interactions in Quinn's wlc_p scheduler
! TYPE: INTEGER

#define WLC_P__N_KAP_ON 1
! VARIABLE COMMENT: ! when to turn KAP energy on
! DEFAULTS COMMENT: ! turn on compression energy immediately
! TYPE: INTEGER

#define WLC_P__N_CHI_ON 1
! VARIABLE COMMENT: ! when to turn CHI energy on
! DEFAULTS COMMENT: ! turn on chi energy immediately
! TYPE: INTEGER

#define WLC_P__N_CHI_L2_ON 1
! VARIABLE COMMENT: ! when to turn CHI energy on
! DEFAULTS COMMENT: ! turn on chi energy immediately
! TYPE: INTEGER

#define WLC_P__NINITMCSTEPS 4000
! VARIABLE COMMENT: ! number of mc steps before starting BD
! DEFAULTS COMMENT: ! number of initilizing mc steps. 1000s x num polymers is good
! TYPE: INTEGER


#define WLC_P__MOVESPERSTEP_CRANK_SHAFT 11
! VARIABLE COMMENT: ! how many times each move is done per step
! DEFAULTS COMMENT: ! crank shaft
! TYPE: INTEGER

#define WLC_P__MOVESPERSTEP_SLIDE_MOVE 25
! VARIABLE COMMENT: ! how many times each move is done per step
! DEFAULTS COMMENT: ! slide
! TYPE: INTEGER

#define WLC_P__MOVESPERSTEP_PIVOT_MOVE 19
! VARIABLE COMMENT: ! how many times each move is done per step
! DEFAULTS COMMENT: ! pivot
! TYPE: INTEGER

#define WLC_P__MOVESPERSTEP_ROTATE_MOVE 45
! VARIABLE COMMENT: ! how many times each move is done per step
! DEFAULTS COMMENT: ! rotate
! TYPE: INTEGER

#define WLC_P__MOVESPERSTEP_FULL_CHAIN_ROTATION 5
! VARIABLE COMMENT: ! how many times each move is done per step
! DEFAULTS COMMENT: ! full chain rotation
! TYPE: INTEGER

#define WLC_P__MOVESPERSTEP_FULL_CHAIN_SLIDE 6
! VARIABLE COMMENT: ! how many times each move is done per step
! DEFAULTS COMMENT: ! full chain slide
! TYPE: INTEGER

#define WLC_P__MOVESPERSTEP_CHANGE_BINDING_STATE 111
! VARIABLE COMMENT: ! how many times each move is done per step
! DEFAULTS COMMENT: ! bind
! TYPE: INTEGER

#define WLC_P__MOVESPERSTEP_CHAIN_FLIP 0
! VARIABLE COMMENT: ! how many times each move is done per step
! DEFAULTS COMMENT: ! filp
! TYPE: INTEGER

#define WLC_P__MOVESPERSTEP_CHAIN_EXCHANGE 5
! VARIABLE COMMENT: ! how many times each move is done per step
! DEFAULTS COMMENT: ! full chain swap
! TYPE: INTEGER

#define WLC_P__MOVESPERSTEP_REPTATION 10
! VARIABLE COMMENT: ! how many times each move is done per step
! DEFAULTS COMMENT: ! reptation
! TYPE: INTEGER

#define WLC_P__MOVESPERSTEP_SUPER_REPTATION 57
! VARIABLE COMMENT: ! how many times each move is done per step
! DEFAULTS COMMENT: ! super-reptation
! TYPE: INTEGER

#define WLC_P__PROBSINGLESWAP 0.5_dp
! VARIABLE COMMENT: ! What fraction of the time to try to change only one tail
! DEFAULTS COMMENT: ! 0< x  <1
! TYPE: LOGICAL



! --------------------------
!
!   Parallel Tempering settings
!
! --------------------------

#define WLC_P__PTON .FALSE.
! VARIABLE COMMENT: ! whether or not to parallel temper
! DEFAULTS COMMENT: ! use parallel if applicable
! TYPE: LOGICAL

#define WLC_P__PT_TWIST .FALSE.
! VARIABLE COMMENT: None
! DEFAULTS COMMENT: ! don't parallel temper linking number (twist) by default
! TYPE: LOGICAL

#define WLC_P__PT_CHI .FALSE.
! VARIABLE COMMENT: None
! DEFAULTS COMMENT: ! don't parallel temper chi by default
! TYPE: LOGICAL

#define WLC_P__PT_H .FALSE.
! VARIABLE COMMENT: ! Parallel temper applied field
! DEFAULTS COMMENT: ! bindng energy is external appled field
! TYPE: LOGICAL

#define WLC_P__PT_KAP .FALSE.
! VARIABLE COMMENT: None
! DEFAULTS COMMENT: ! don't parallel temper kap by default
! TYPE: LOGICAL

#define WLC_P__PT_MU .FALSE.
! VARIABLE COMMENT: None
! DEFAULTS COMMENT: ! don't parallel temper mu by default
! TYPE: LOGICAL

#define WLC_P__PT_COUPLE .FALSE.
! VARIABLE COMMENT: None
! DEFAULTS COMMENT: ! don't parallel temper HP1 binding by default
! TYPE: LOGICAL

#define WLC_P__PT_MAIERSAUPE .FALSE.
! VARIABLE COMMENT: None
! DEFAULTS COMMENT: ! don't parallel temper maier Saupe by default
! TYPE: LOGICAL

#define WLC_P__NREPADAPT 1000
! VARIABLE COMMENT: ! number of exchange attemts between adapt
! DEFAULTS COMMENT: ! 1000 exchange attempts between adaptations
! TYPE: INTEGER

#define WLC_P__LOWERREPEXE 0.09_DP
! VARIABLE COMMENT: ! when to decrese cof spacing
! DEFAULTS COMMENT: ! TOdo: enter justification for these defaults, if any.
! TYPE: REAL(DP)

#define WLC_P__UPPERREPEXE 0.18_DP
! VARIABLE COMMENT: ! when to increase cof spacing
! DEFAULTS COMMENT: ! TOdo: fine if the only justification is "these just work"
! TYPE: REAL(DP)

#define WLC_P__LOWERCOFRAIL 0.005_DP
! VARIABLE COMMENT: ! minumum acceptable Cof
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__UPPERCOFRAIL 0.1_DP
! VARIABLE COMMENT: ! maximum acceptable Cof
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__INDSTARTREPADAPT 10
! VARIABLE COMMENT: None
! DEFAULTS COMMENT: None
! TYPE: INTEGER

#define WLC_P__INDENDREPADAPT 20
! VARIABLE COMMENT: None
! DEFAULTS COMMENT: None
! TYPE: INTEGER

#define WLC_P__REPANNEALSPEED 0.01_DP
! VARIABLE COMMENT: ! for annealing
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__REPLICABOUNDS .TRUE.
! VARIABLE COMMENT: None
! DEFAULTS COMMENT: None
! TYPE: LOGICAL

#define WLC_P__INITIAL_MAX_S NAN
! VARIABLE COMMENT: ! Parameterized location along parallel temper path
! DEFAULTS COMMENT: ! You should set this variable if you intend to use it
! TYPE: REAL(DP)

! --------------------------
!
!   Available MC Moves
!
! --------------------------

#define WLC_P__MOVEON_CRANK_SHAFT 1
! VARIABLE COMMENT: ! Is the move active
! DEFAULTS COMMENT: ! crank-shaft move
! TYPE: INTEGER

#define WLC_P__MOVEON_SLIDE_MOVE 1
! VARIABLE COMMENT: ! Is the move active
! DEFAULTS COMMENT: ! slide move
! TYPE: INTEGER

#define WLC_P__MOVEON_PIVOT_MOVE 1
! VARIABLE COMMENT: ! Is the move active
! DEFAULTS COMMENT: ! pivot move
! TYPE: INTEGER

#define WLC_P__MOVEON_ROTATE_MOVE 1
! VARIABLE COMMENT: ! Is the move active
! DEFAULTS COMMENT: ! rotate move
! TYPE: INTEGER

#define WLC_P__MOVEON_FULL_CHAIN_ROTATION 1
! VARIABLE COMMENT: ! Is the move active
! DEFAULTS COMMENT: ! full chain rotation
! TYPE: INTEGER

#define WLC_P__MOVEON_FULL_CHAIN_SLIDE 1
! VARIABLE COMMENT: ! Is the move active
! DEFAULTS COMMENT: ! full chain slide
! TYPE: INTEGER

#define WLC_P__MOVEON_CHANGE_BINDING_STATE 1
! VARIABLE COMMENT: ! Is the move active
! DEFAULTS COMMENT: ! Change in Binding state
! TYPE: INTEGER

#define WLC_P__MOVEON_CHAIN_FLIP 0
! VARIABLE COMMENT: ! Is the move active
! DEFAULTS COMMENT: ! Chain flip ! TOdo not working
! TYPE: INTEGER

#define WLC_P__MOVEON_CHAIN_EXCHANGE 1
! VARIABLE COMMENT: ! Is the move active
! DEFAULTS COMMENT: ! Chain exchange
! TYPE: INTEGER

#define WLC_P__MOVEON_REPTATION 1
! VARIABLE COMMENT: ! Is the move active
! DEFAULTS COMMENT: ! Reptation
! TYPE: INTEGER

#define WLC_P__MOVEON_SUPER_REPTATION 0
! VARIABLE COMMENT: ! Is the move active
! DEFAULTS COMMENT: ! Super Reptation
! TYPE: INTEGER



! ----------------------------
!
!   MC Optimization
!
! ---------------------------

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
! VARIABLE COMMENT: ! desired hit rate
! DEFAULTS COMMENT: ! Target
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
! VARIABLE COMMENT: ! Max Size of window for bead selection
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__MINWINDOW_CRANK_SHAFT NAN
! VARIABLE COMMENT: ! Min Size of window for bead selection
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__MINWINDOW_SLIDE_MOVE NAN
! VARIABLE COMMENT: ! Min Size of window for bead selection
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__MINWINDOW_PIVOT_MOVE NAN
! VARIABLE COMMENT: ! Min Size of window for bead selection
! DEFAULTS COMMENT: None
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

#define WLC_P__MINWINDOW_CHANGE_BINDING_STATE NAN
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
! VARIABLE COMMENT: ! Min Size of window for bead selection
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
! VARIABLE COMMENT: ! minium amplitude
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

#define WLC_P__MAXAMP_FULL_CHAIN_SLIDE 0.1*WLC_P__L
! VARIABLE COMMENT: ! maximum amplitude
! DEFAULTS COMMENT: ! a fraction of the change length
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
! VARIABLE COMMENT: ! maximum amplitude
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
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__WINTARGET_ROTATE_MOVE 8.0_DP
! VARIABLE COMMENT: ! target for ratio of window to anmplitude
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__WINTARGET_FULL_CHAIN_ROTATION 8.0_DP
! VARIABLE COMMENT: ! target for ratio of window to anmplitude
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__WINTARGET_FULL_CHAIN_SLIDE 8.0_DP
! VARIABLE COMMENT: ! target for ratio of window to anmplitude
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__WINTARGET_CHANGE_BINDING_STATE 8.0_DP
! VARIABLE COMMENT: ! target for ratio of window to anmplitude
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__WINTARGET_CHAIN_FLIP 8.0_DP
! VARIABLE COMMENT: ! target for ratio of window to anmplitude
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__WINTARGET_CHAIN_EXCHANGE 8.0_DP
! VARIABLE COMMENT: ! target for ratio of window to anmplitude
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__WINTARGET_REPTATION 8.0_DP
! VARIABLE COMMENT: ! target for ratio of window to anmplitude
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__WINTARGET_SUPER_REPTATION 8.0_DP
! VARIABLE COMMENT: ! target for ratio of window to anmplitude
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

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
! VARIABLE COMMENT: ! Nunber of steps between adapt
! DEFAULTS COMMENT: ! adapt after at most 1000 steps
! TYPE: INTEGER

#define WLC_P__MIN_ACCEPT 0.05_DP
! VARIABLE COMMENT: ! threshold for deciding to usually not use a move
! DEFAULTS COMMENT: ! if a move succeeds < 5% of the time, start using it only every reduce_move cycles
! TYPE: REAL(DP)

#define WLC_P__REDUCE_MOVE 1
! VARIABLE COMMENT: ! whether or not to stop usuing a move when it goes below min_accept success
! DEFAULTS COMMENT: ! use moves that fall below the min_accept threshold only once every ~ times they would otherwise be used, set to one for no effect
! TYPE: INTEGER

#define WLC_P__WINTYPE 1
! VARIABLE COMMENT: ! distributionof segment size in crankshaft move (unif = 0, exp = 1)
! DEFAULTS COMMENT: ! exponential fragment sizes mix better
! TYPE: INTEGER



