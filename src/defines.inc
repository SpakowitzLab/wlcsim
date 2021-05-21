!-------------------------------------------
! This defines.inc file contains default values, it will not produce a meaning simulation.
! Values are replaced in codebase at compile time.
! See inputs/example_defines.inc for examples.
! ------------------------------------------


! ----------------------------
!
!   General Switches
!
! ---------------------------

#define WLC_P__CODENAME ""
! VARIABLE COMMENT: ! Which simulation code to run (from when we combined our code).  For example "quinn" will run the wlcsim_quinn simulation which uses MPI for parallel tempering.  Options include "quinn", "bruno", and "brad".
! DEFAULTS COMMENT: ! Left blank to require an input.
! TYPE: CHARACTER(MAXPARAMLEN)

#define WLC_P__WARNING_LEVEL 1
! VARIABLE COMMENT: ! Level of coded warnings to use.  0 = no warnings, 1 = moderate performance impact, 2 = considerable slowdown.
! DEFAULTS COMMENT: ! Warning level 1 is recommended.
! TYPE: INTEGER

#define WLC_P__RING .FALSE.
! VARIABLE COMMENT: ! Whether the polymer is a ring (i.e. ends connected).
! DEFAULTS COMMENT: ! Not a ring by default.
! TYPE: LOGICAL

#define WLC_P__LOCAL_TWIST .FALSE.
! VARIABLE COMMENT: ! Whether to include local twist energy, that is, calculated from the angle between sequential bead orientations.
! DEFAULTS COMMENT: ! Twist not included by default.
! TYPE: LOGICAL

#define WLC_P__TWIST .FALSE.
! VARIABLE COMMENT: ! Whether to include global twist energy via Twist = Link - Writhe.
! DEFAULTS COMMENT: ! Twist not included by default.
! TYPE: LOGICAL

#define WLC_P__CONFINETYPE 'NONE'
! VARIABLE COMMENT: ! Boundary condition: "platesInZperiodicXY", "cube", "sphere", "excludedShpereInPeriodic", "ecoli", "None". "None" produces periodic boundaries if interactions are included. 
! DEFAULTS COMMENT: ! No confinement by default.
! TYPE: CHARACTER(MAXPARAMLEN)

#define WLC_P__BOUNDARY_TYPE "ExtendBinsPast"
! VARIABLE COMMENT: ! How to handle boundaries for interpolated binned interaction choose from: "SolidEdgeBin", "ExtendBinsPast".  "ExtendBinsPast" treats all bins the same and requires bins to extend at least one bin beyond the confinement.  "SolidEdgeBin" adjusts boundary bins so that all density is maintained within the confinement.  "SolidEdgeBin" doesn’t work for curved boundaries.
! DEFAULTS COMMENT: ! "ExtendBinsPast" is default because it is compatible with more boundary types.
! TYPE: CHARACTER(MAXPARAMLEN)

#define WLC_P__FRACTIONAL_BIN .FALSE.
! VARIABLE COMMENT: ! Whether to use fractional bin volumes in field treatment.  Only an issue for curved (e.g. spherical) confinement.  If used, the density for parietal bins is calculated by dividing by the partial volume.
! DEFAULTS COMMENT: ! If in doubt, make it false.
! TYPE: LOGICAL

#define WLC_P__INITCONDTYPE 'randomWalkWithBoundary'
! VARIABLE COMMENT: ! Initial condition for bead positions.  Choices include: "lineInYFromOrigin", "lineInY", "randomLineSlitInZBoundary", "randomLineCubeBoundary", "randomLineOutsideOfSphere", "randomLineSphereBoundary", "randomlyDistributeBeadsInSphere", "ring", "multiRing"
! DEFAULTS COMMENT: ! 'randomWalkWithBoundary' works for many boundary conditions.
! TYPE: CHARACTER(MAXPARAMLEN)

#define WLC_P__INTRAPOLYMER_STICK_CROSSING_ENFORCED .FALSE.
! VARIABLE COMMENT: ! Included to prevent polymers from ever crossing each other and themselves.  The code to implement this is not fully implemented.
! DEFAULTS COMMENT: ! NOT IMPLEMENTED
! TYPE: LOGICAL

#define WLC_P__CHEM_STATE_FROM_FILE .FALSE.
! VARIABLE COMMENT: ! Whether to read the initial chemical sequence from a file.  In the chromatin simulations this corresponds to the HP1 binding profile and is therefore FALSE.  For block-copolymers this corresponds to A/B block bead pattern.  Reads from "input/ab".
! DEFAULTS COMMENT: ! Don't load initial a/b states from file
! TYPE: 

#define WLC_P__CHEM_SEQ_FROM_FILE .FALSE.
! VARIABLE COMMENT: ! Read the fixed methylation profile (or equivalent for other problems) from a file.  Read from "input/meth" or "input/meth_..." if WLC_P__ENSEMBLE_METH.
! DEFAULTS COMMENT: ! False by default.
! TYPE: LOGICAL

#define WLC_P__FRMFILE .FALSE.
! VARIABLE COMMENT: ! Read initial bead positions from a file.
! DEFAULTS COMMENT: ! don't load initial bead positions from file
! TYPE: LOGICAL

#define WLC_P__FRMFIELD .FALSE.
! VARIABLE COMMENT: ! Read initial field from file.  Useful for spatially dependent fields.  Note: this is not the only way to put in a potential.  Also see WLC_P__APPLY_EXTERNAL_FIELD.
! DEFAULTS COMMENT: ! Don't load initial field values from file.
! TYPE: LOGICAL

#define WLC_P__COLLISIONDETECTIONTYPE 0
! VARIABLE COMMENT: ! Save first passage time vectors to file.  Used for brownian dynamics simulation.
! DEFAULTS COMMENT: ! Don't track first passage time collisions between beads
! TYPE: INTEGER

#define WLC_P__RESTART .FALSE.
! VARIABLE COMMENT: ! Whether we are restarting from a previously saved simulation or not.  Code very out of date, don’t use.
! DEFAULTS COMMENT: ! Don't restart from previously saved simulation
! TYPE: LOGICAL

#define WLC_P__INTERP_BEAD_LENNARD_JONES .FALSE.
! VARIABLE COMMENT: ! Whether to have inter bead lennard jones energies.  That is lenerard jones repulsion from the point of closest approach between two sticks.  Code may be out of date.
! DEFAULTS COMMENT: ! No interactions by default.
! TYPE: LOGICAL

#define WLC_P__CHANGINGCHEMICALIDENTITY .FALSE.
! VARIABLE COMMENT: ! Whether the bead types change during the simulation.  For example, by HP1 binding and unbinding.
! DEFAULTS COMMENT: ! Only one type of bead by default.
! TYPE: LOGICAL

#define WLC_P__VARIABLE_CHEM_STATE WLC_P__CHANGINGCHEMICALIDENTITY
! VARIABLE COMMENT: ! Sames as WLC_P__CHANGINGCHEMICALIDENTITY.
! DEFAULTS COMMENT: ! Sames as WLC_P__CHANGINGCHEMICALIDENTITY.
! TYPE: LOGICAL

#define WLC_P__FIELD_INT_ON .FALSE.
! VARIABLE COMMENT: Whether or not to execute mc_int, i.e. whether to include the field theoretic (bined) interaction.
! DEFAULTS COMMENT: ! No interactions by default.
! TYPE: LOGICAL

#define WLC_P__FIELDINTERACTIONTYPE 'NONE'
! VARIABLE COMMENT: ! Choose the interbead interaction type you would like, e.g. "AppliedAligningField", "AppliedAligningFieldMelt",  "MaierSaupe", "ABmelt", "ABsolution", "chromatin", "chromatin2".  Interactions defined in src/mc/MC_Hamiltonian.f03.
! DEFAULTS COMMENT: ! No interactions by default.
! TYPE: CHARACTER(MAXPARAMLEN)

#define WLC_P__RECENTER_ON .FALSE.
! VARIABLE COMMENT: ! Whether to move drifted chains back in periodic BC to prevent rounding error.  Use only for periodic boundary conditions or only one chain with no boundaries or field interaction.
! DEFAULTS COMMENT: ! Often is not needed.  You will know if you need it from the output.
! TYPE: LOGICAL

#define WLC_P__CYLINDRICAL_CHAIN_EXCLUSION .FALSE.
! VARIABLE COMMENT: ! Explicitly simulaties chains as cylinders that can’t occupy the same space at the same time.
! DEFAULTS COMMENT: ! No interactions by default.
! TYPE: LOGICAL

#define WLC_P__GJK_STERICS .FALSE.
! VARIABLE COMMENT: ! GJK algorithm to model DNA + nucleosomes as convex hulls
! DEFAULTS COMMENT: ! Don't check for sterics via GJK by default.
! TYPE: LOGICAL

#define WLC_P__GJK_POLYGON 16
! VARIABLE COMMENT: ! sides of polygon to model with GJK (i.e 2*sides)
! DEFAULTS COMMENT: ! octagons seem natural choice for nucleosomes
! TYPE: INTEGER

#define WLC_P__NEIGHBOR_BINS .FALSE.
! VARIABLE COMMENT: ! Keep track of neighbors by recursive binning for pairwise interactions.  Not needed if interactions are accomplished via a field.
! DEFAULTS COMMENT: ! Set off because it is expensive.
! TYPE: LOGICAL

#define WLC_P__APPLY_EXTERNAL_FIELD .FALSE.
! VARIABLE COMMENT: ! Externally applied field (e.g. radial binding)
! DEFAULTS COMMENT: ! None applied
! TYPE: LOGICAL

#define WLC_P__EXTERNAL_FIELD_TYPE "none"
! VARIABLE COMMENT: ! Explicitly defined potential type to test various conditions.  Options: "nonSpecificToCubeSide",  "nonSpecificToSphere", "toExcludedSphereInPeriodic", "NSTCS_plus_binding".  Used, for example, for peripheral attraction.
! DEFAULTS COMMENT: ! No external field by default.
! TYPE: CHARACTER(MAXPARAMLEN)

#define WLC_P__APPLY_2body_potential .FALSE.
! VARIABLE COMMENT: ! Explicit pairwise interaction potential between beads.
! DEFAULTS COMMENT: ! No interactions by default.
! TYPE: LOGICAL

#define WLC_P__EXPLICIT_BINDING .FALSE.
! VARIABLE COMMENT: ! Include specific interactions between specified beads
! DEFAULTS COMMENT: ! Not included by default
! TYPE: LOGICAL

#define WLC_P__NETWORK .FALSE.
! VARIABLE COMMENT: ! Pin pairs of beads together (to enforce the network).  For example, used for loop extrusion factors.
! DEFAULTS COMMENT: ! Not included by default.
! TYPE: LOGICAL

#define WLC_P__NO_LEF_CROSSING .FALSE.
! VARIABLE COMMENT: ! No explicit binding (Loop Extrusion Factors) pairs can cross each other, e.g. 1-to-10 and 5-15.  This parameter only matters when loop extrusion factor positions are simulated rather than loaded explicitly.
! DEFAULTS COMMENT: ! LEF’s not included so it doesn’t matter.
! TYPE: LOGICAL

#define WLC_P__ELASTICITY_TYPE "constant"
! VARIABLE COMMENT: ! whether to allow differential elasticity along polymer, also allows for kinks: "constant", "nucleosome"
! DEFAULTS COMMENT: ! by default, 'constant' bead spacing
! TYPE: CHARACTER(MAXPARAMLEN)

#define WLC_P__INCLUDE_NUC_TRANS .TRUE.
! VARIABLE COMMENT: ! Whether to include the translation from entry to exit position of the nucleosome.
! DEFAULTS COMMENT: ! By default include it.
! TYPE: LOGICAL

#define WLC_P__INCLUDE_DISCRETIZE_LINKER .FALSE.
! VARIABLE COMMENT: ! Whether to discretize linker DNA and nucleosomes
! DEFAULTS COMMENT: ! turn on for risca, i.e. meso-scale sims
! TYPE: LOGICAL

#define WLC_P__CHI_L2_ABLE .FALSE.
! VARIABLE COMMENT: ! Include quadrupole interactions.  That is, whether the l=2 component of the field is ever calculated.
! DEFAULTS COMMENT: ! Not included by default.
! TYPE: LOGICAL

#define WLC_P__UMBRELLA .FALSE.
! VARIABLE COMMENT: ! Perfore umbrella sampling.
! DEFAULTS COMMENT: ! No umbrella sampling by default.
! TYPE: LOGICAL

#define WLC_P__NO_SELF_CROSSING .FALSE.
! VARIABLE COMMENT: ! Whether to prevent self-crossing of ring polymer by conserving the linking number.
! DEFAULTS COMMENT: ! Allow self-crossing by default.
! TYPE: LOGICAL

! ----------------------------
!
!   Output
!
! ---------------------------

#define WLC_P__EXITWHENCOLLIDED .FALSE.
! VARIABLE COMMENT: ! Stop sim after every n^2 loops have formed at least once (for brownian dynamics simulations).
! DEFAULTS COMMENT: ! Don't stop sim when coltimes is full by default.
! TYPE: LOGICAL

#define WLC_P__SAVER .TRUE.
! VARIABLE COMMENT: ! Save position vectors to file "data/r...".
! DEFAULTS COMMENT: ! Save orientation vectors by default.
! TYPE: LOGICAL

#define WLC_P__SAVEU .TRUE.
! VARIABLE COMMENT: ! Save orientation vectors to file "data/u..."
! DEFAULTS COMMENT: ! Save orientation vectors.  This is needed if you want to start another simulation based on the output of this simulation.
! TYPE: LOGICAL

#define WLC_P__SAVEAB .FALSE.
! VARIABLE COMMENT: ! Save AB (chemical identity/HP1 binding) to file.
! DEFAULTS COMMENT: ! Dont' save AB by default.
! TYPE: LOGICAL

#define WLC_P__SAVEPHI .FALSE.
! VARIABLE COMMENT: ! Save densities, Phi, values to file.
! DEFAULTS COMMENT: ! Don't save A/B density per bin (you can always be calculated from the positions anyway).
! TYPE: LOGICAL

#define WLC_P__APPEND_R .FALSE.
! VARIABLE COMMENT: ! Output positions from multiple savepoints to the same file.
! DEFAULTS COMMENT: ! Use the separate file per savepoint format.
! TYPE: LOGICAL

#define WLC_P__SAVE_LK .FALSE.
! VARIABLE COMMENT: ! Save linking number. Only for ring polymer.
! DEFAULTS COMMENT: ! Does not save linking number.
! TYPE: LOGICAL

! ----------------------------
!
!   Number / Chemical identity variables
!
! ---------------------------

#define WLC_P__NB 201
! VARIABLE COMMENT: ! Number of beads in a polymer.  Note that NB = nMpP*nBpM.
! DEFAULTS COMMENT: ! 201 beads per polymer
! TYPE: INTEGER

#define WLC_P__NMPP WLC_P__NB/WLC_P__NBPM
! VARIABLE COMMENT: ! Number of monomers (NOT BEADS!) in a polymer.  Monomers only matter if chemical identities are specific.  Typically useful for block copolymers where each block is considered a monomer.
! DEFAULTS COMMENT: ! Set by default.
! TYPE: INTEGER

#define WLC_P__NBPM 1
! VARIABLE COMMENT: ! Number beads per monomer.  Set to 1 in must circumstances.
! DEFAULTS COMMENT: ! Assume each bead is a monomer.
! TYPE: INTEGER

#define WLC_P__NP 1
! VARIABLE COMMENT: ! Number of polymers  in simulation.
! DEFAULTS COMMENT: ! Only one polymer.
! TYPE: INTEGER

#define WLC_P__NT WLC_P__NB*WLC_P__NP
! VARIABLE COMMENT: ! Total Number of beads.  Set to WLC_P__NB*WLC_P__NP for the same number of beads in each monomer.
! DEFAULTS COMMENT: ! Assuming same number of beads on each polymer.
! TYPE: INTEGER

#define WLC_P__NBPL 0
! VARIABLE COMMENT: ! Number of beads per linker segement, this determines discretization (only relevant for risca chromatin sims)
! DEFAULTS COMMENT: ! 0 for simulations outside risca lab / chromatin with discretized linkers
! TYPE: INTEGER

#define WLC_P__LINKING_NUMBER 0
! VARIABLE COMMENT: ! Linking number of ring for global twist energy calculation.
! DEFAULTS COMMENT: ! Simple loop.
! TYPE: REAL(DP)

#define WLC_P__POLY_DISP_TYPE "None"
! VARIABLE COMMENT: ! Whether to include chains of different lengths.  "None" for all the same length, "FromFile" for different lengths.
! DEFAULTS COMMENT: ! None by default.
! TYPE: CHARACTER(MAXPARAMLEN)

#define WLC_P__FA 0.5_DP
! VARIABLE COMMENT: ! Fraction of A beads.  For random copolymers.
! DEFAULTS COMMENT: ! Half A, half B by default.
! TYPE: REAL(DP)

#define WLC_P__LAM 0.0_DP
! VARIABLE COMMENT: ! Chemical correlation parameter (eigenvalue of transition matrix that generates A/B's).  Used for random copolymers (if not input from file).
! DEFAULTS COMMENT: ! Perfectly random sequence  (see generating_sequences.rst for details)
! TYPE: REAL(DP)

#define WLC_P__F_METH 0.5_DP
! VARIABLE COMMENT: ! Fraction methylated for random copolymers (if not input from file).
! DEFAULTS COMMENT: ! Half beads methylated by default.
! TYPE: REAL(DP)

#define WLC_P__LAM_METH 0.9_DP
! VARIABLE COMMENT: ! Eigenvalue of transition matrix generating initial methylation sequence (if not input from file).
! DEFAULTS COMMENT: ! Highly alternating sequence by default.
! TYPE: REAL(DP)

#define WLC_P__ASYMMETRICALTERNATINGCHEM .FALSE.
! VARIABLE COMMENT: ! initialize asymmetric alturnating chain
! DEFAULTS COMMENT: None
! TYPE: LOGICAL

#define WLC_P__TWO_TAIL .FALSE.
! VARIABLE COMMENT: ! Chromatin with two histone tails per bead.
! DEFAULTS COMMENT: ! Specific to chromatin simulation, therefore turned off.
! TYPE: LOGICAL

#define WLC_P__N_SPHERES_TO_SIDE 1
! VARIABLE COMMENT: ! Number of spheres along a slide.  This is used for simulations of polymers excluded from multiple spheres.
! DEFAULTS COMMENT: ! Set to one. Won’t be used unless you use corresponding B.C.
! TYPE: INTEGER

#define WLC_P__MAX_BACEPAIRS_PER_BEAD 350
! VARIABLE COMMENT: ! Maximum base pairs per bead.  For modeling DNA and using WLC_P__ELASTICITY_TYPE="nucleosome", i.e. with nucleosome genometry included.
! DEFAULTS COMMENT: ! 203 linker + 147 wrap
! TYPE: INTEGER

#define WLC_P__N_EXTERNAL_BIND_POINTS 10
! VARIABLE COMMENT: ! Number of points bound to boundary (originally for Quinn's lamin project).
! DEFAULTS COMMENT: ! 0 for off
! TYPE: INTEGER

#define WLC_P__TWIST_DENSITY 0.0_DP
! VARIABLE COMMENT: ! Baseline twist density.  E.g. natural twist of DNA if not already subtracted out. 
! DEFAULTS COMMENT: ! baseline twist density 0.0
! TYPE: REAL(DP)

#define WLC_P__NUCLEOSOME_WRAPPING 147
! VARIABLE COMMENT: ! How many basepairs wrap a nucleosome (this then defines the entry-exit geometry of the nucleosome)
! DEFAULTS COMMENT: ! 147 bp wrapping
! TYPE: REAL(DP)

! ----------------------------
!
!   Length variables
!   All lengths in simulations units unless otherwise state.
!
! ---------------------------

#define WLC_P__LENGTH_PER_BP 0.332_DP
! VARIABLE COMMENT: ! Length along DNA per bace pair.  Currently not in use.
! DEFAULTS COMMENT: ! known value in nm
! TYPE: REAL(DP)

#define WLC_P__LL 38_DP
! VARIABLE COMMENT: ! linker length in bp between adjacent nucleosomes.  For modeling DNA and using WLC_P__ELASTICITY_TYPE="nucleosome" i.e. with nucleosome genometry included.
! DEFAULTS COMMENT: ! 38bp is like mouse-ish
! TYPE: REAL(DP)

#define WLC_P__LINKER_DISCRETIZATION 0
! VARIABLE COMMENT: ! Discretization of computational beads (only relevant for risca chromatin sims)
! DEFAULTS COMMENT: ! 0 for simulations outside risca lab / chromatin with discretized linkers
! TYPE: REAL(DP)

#define WLC_P__NUM_NUCLEOSOMES 0
! VARIABLE COMMENT: ! how many nucleosomes per chain (only relevant for risca chromatin sims)
! DEFAULTS COMMENT: ! 0 for simulations outside risca lab / chromatin with discretized linkers
! TYPE: INTEGER

#define WLC_P__LINKER_TYPE 'phased'
! VARIABLE COMMENT: ! linker length geometry type for initialization (i.e. phased array or sampled from a file)
! DEFAULTS COMMENT: ! default to phased, i.e. set up a chain of constant linker length
! TYPE: REAL(DP)

#define WLC_P__L WLC_P__LENGTH_PER_BP*WLC_P__LL*(WLC_P__NB-1)
! VARIABLE COMMENT: ! Length of each polymer in simulation.
! DEFAULTS COMMENT: ! By default, constant linker lengths of LL assumed.
! TYPE: REAL(DP)

#define WLC_P__L0 (WLC_P__L/(WLC_P__NB-1.0_DP))
! VARIABLE COMMENT: ! Bead spacing, Remove the "-1.0_dp" if ring !TODO make #if
! DEFAULTS COMMENT: ! Derive from L and NB.
! TYPE: REAL(DP)

#define WLC_P__LP 50.0_DP
! VARIABLE COMMENT: ! Persistence length of polymer.
! DEFAULTS COMMENT: ! in nm
! TYPE: REAL(DP)

#define WLC_P__LT 100.0_DP
! VARIABLE COMMENT: ! twist persistence length
! DEFAULTS COMMENT: ! approximent twist persistance length of DNA in nm
! TYPE: REAL(DP)

#define WLC_P__BEADVOLUME 0.1_DP
! VARIABLE COMMENT: ! Bead volume (for binned density interaction).  In cubic simulation units.
! DEFAULTS COMMENT: ! much smaller than space between beads
! TYPE: REAL(DP)

#define WLC_P__COLLISIONRADIUS 0
! VARIABLE COMMENT: ! Radius triggering collisions to be recorded in "coltimes" (brownian simulation only)
! DEFAULTS COMMENT: ! Never collide except on floating-point coincidence.
! TYPE: REAL(DP)

#define WLC_P__K_FIELD 0.0_DP
! VARIABLE COMMENT: ! Wave vector of applied sinusoidal field (used in PT to step around 1st order phase transition).  This has rarely been used.
! DEFAULTS COMMENT: ! some previous values: !1.5708_dp !0.3145_dp
! TYPE: REAL(DP)

#define WLC_P__NCOLBIN 1
! VARIABLE COMMENT: ! Number of collision-detection bins on each edge (internal setting for binned based pairwise interactions).
! DEFAULTS COMMENT: ! Equivalent to collisionDetectionType = 1
! TYPE: INTEGER

#define WLC_P__NBIN_X 20
! VARIABLE COMMENT: ! Number of field theoretic density calculation bins in the X direction.
! DEFAULTS COMMENT: ! set this yourself!
! TYPE: INTEGER

#define WLC_P__NBIN_Y WLC_P__NBIN_X
! VARIABLE COMMENT: ! Number of field theoretic density calculation bins in the Y direction.
! DEFAULTS COMMENT: ! assume cube
! TYPE: INTEGER

#define WLC_P__NBIN_Z WLC_P__NBIN_X
! VARIABLE COMMENT: ! Number of field theoretic density calculation bins in the Z direction.
! DEFAULTS COMMENT: ! assume cube
! TYPE: INTEGER

#define WLC_P__LBOX_X 20.0_dp
! VARIABLE COMMENT: ! Field theoretic box length in X.  Usually you should make this NBIN_X*DBIN.
! DEFAULTS COMMENT: ! Box size/confinement, *MUST* be set by user.
! TYPE: REAL(DP)

#define WLC_P__LBOX_Y WLC_P__LBOX_X
! VARIABLE COMMENT: ! Field theoretic box length in Y.  Usually you should make this NBIN_Y*DBIN.
! DEFAULTS COMMENT: ! assume cube
! TYPE: REAL(DP)

#define WLC_P__LBOX_Z WLC_P__LBOX_X
! VARIABLE COMMENT: ! Field theoretic box length in Z.  Usually you should make this NBIN_Z*DBIN.
! DEFAULTS COMMENT: ! assume cube
! TYPE: REAL(DP)

#define WLC_P__DBIN 1.0_dp
! VARIABLE COMMENT: ! Monte carlo field bin discretization size.  Be sure you check the code if you set to anything other than 1.0.
! DEFAULTS COMMENT: ! Should always be 1.0_dp according to Quinn.
! TYPE: REAL(DP)

#define WLC_P__CONFINEMENT_CUBE_LENGTH WLC_P__LBOX_Z
! VARIABLE COMMENT: ! Size of confinement box. Often set to: LBOX.  Set to: LBOX - 2*DBIN if WLC_P__BOUNDARY_TYPE == "ExtendbinsPast".
! DEFAULTS COMMENT: ! Set to match LBOX_Z
! TYPE: REAL(DP)

#define WLC_P__CONFINEMENT_SLIT_WIDTH NAN
! VARIABLE COMMENT: !  Size of confinement. Set = LBOX -2*DBIN if bin interation is on
! DEFAULTS COMMENT: ! Must set if in use
! TYPE: REAL(DP)

#define WLC_P__CONFINEMENT_SPHERE_DIAMETER WLC_P__LBOX_X - 2*WLC_P__DBIN
! VARIABLE COMMENT: ! Size of confinement. Often set to: LBOX.  Set to: LBOX - 2*DBIN if WLC_P__BOUNDARY_TYPE == "ExtendbinsPast".
! DEFAULTS COMMENT: ! Must set if in use
! TYPE: REAL(DP)

#define WLC_P__CONFINEMENT_ECOLI_LENGTH 0.0_dp
! VARIABLE COMMENT: ! Length of the cylindrical part of ecoli shape.  Set to zero to get a sphere.
! DEFAULTS COMMENT: ! zero for sphere
! TYPE: REAL(DP)

#define WLC_P__CONFINEMENT_ECOLI_DIAMETER 20.0_dp
! VARIABLE COMMENT: ! Diameter of spherical caps for ecoli shape.  This shape is a cylinder with spherical caps.
! DEFAULTS COMMENT: ! Must set if in use.
! TYPE: REAL(DP)

#define WLC_P__BINDING_R 0.0_DP
! VARIABLE COMMENT: ! Binding Radius of extenral field.  For example, range of boundary potential for peripheral heterochromatin paper.
! DEFAULTS COMMENT: ! Not included.
! TYPE: REAL(DP)

#define WLC_P__HP1_BINDING_R 0.0_DP
! VARIABLE COMMENT: ! Binding Radius of HP1. I.e. 3-5nm / DBIN in nm.
! DEFAULTS COMMENT: ! Not included
! TYPE: REAL(DP)

#define WLC_P__CHAIN_D 0.0_DP
! VARIABLE COMMENT: Diameter of chain for cylindrical chain exclusion.
! DEFAULTS COMMENT: Set to zero for no exclusion.
! TYPE: REAL(DP)

#define WLC_P__DNA_RADIUS 1.0_DP
! VARIABLE COMMENT: Radius of DNA (GJK)
! DEFAULTS COMMENT: 1 nm
! TYPE: REAL(DP)

#define WLC_P__NUCLEOSOME_RADIUS 5.2_DP
! VARIABLE COMMENT: Radius of nuclelsome (GJK)
! DEFAULTS COMMENT: 5.2 nm
! TYPE: REAL(DP)

#define WLC_P__NUCLEOSOME_HEIGHT 5.5_DP
! VARIABLE COMMENT: Height of nucleosome (GJK)
! DEFAULTS COMMENT: 5.5. nm
! TYPE: REAL(DP)

#define WLC_P__GJK_RADIUS SQRT(WLC_P__NUCLEOSOME_RADIUS**2 + (WLC_P__NUCLEOSOME_HEIGHT/2)**2)
! VARIABLE COMMENT: ! search radius for GJK alorithm checking
! DEFAULTS COMMENT: ! set to hypotenuse of nuc rad + height
! TYPE: REAL(DP)

! ----------------------------
!
!   Energy variables
!   In kT’s unless otherwise stated.
!
! ---------------------------

#define WLC_P__EU (-1.52_DP)
! VARIABLE COMMENT: ! Energy of HP1 binding for unmethylated sites.
! DEFAULTS COMMENT: ! From Pete's paper.
! TYPE: REAL(DP)

#define WLC_P__EM 0.01_DP
! VARIABLE COMMENT: ! Energy of HP1 binding for methalated sites
! DEFAULTS COMMENT: ! by default, no hp1 binding energy included
! TYPE: REAL(DP)

#define WLC_P__CHI 0.0_DP
! VARIABLE COMMENT: ! Flory-Huggins chi parameter for field theoretic interaction.
! DEFAULTS COMMENT: ! Not in simple Hamiltonian
! TYPE: REAL(DP)

#define WLC_P__MU 0.0_DP
! VARIABLE COMMENT: ! HP1 chemical potential (relative to EM and EU).
! DEFAULTS COMMENT: ! Not in simple Hamiltonian
! TYPE: REAL(DP)

#define WLC_P__HA 0.0_DP
! VARIABLE COMMENT: ! Applied field strength.
! DEFAULTS COMMENT: ! no applied field
! TYPE: REAL(DP)

#define WLC_P__COUPLING_ENERGY 0.0_dp
! VARIABLE COMMENT: ! Coupling energy between beads (negative = attractive).  For example, for HP1 oligomerization.
! DEFAULTS COMMENT: ! no coupling energy
! TYPE: REAL(DP)

#define WLC_P__HP1_BIND 0.0_DP
! VARIABLE COMMENT: ! Strength of coupling from HP1 binding.  Related to coupling energy.  Two parameters needed to subtract self interaction.
! DEFAULTS COMMENT: ! Assumed unwanted
! TYPE: REAL(DP)

#define WLC_P__KAP 0.0_DP
! VARIABLE COMMENT: ! Degree of incompressibility.
! DEFAULTS COMMENT: ! Assumed unwanted.
! TYPE: REAL(DP)

#define WLC_P__CHI_L2 0.0_DP
! VARIABLE COMMENT: ! Proportional to Maier Saupe parameter for liquid crystal polymer solutions.
! DEFAULTS COMMENT: ! zero for simplicity
! TYPE: REAL(DP)

#define WLC_P__EXPLICIT_BIND_ENERGY 0.0_DP
! VARIABLE COMMENT: ! Include specific interactions between specified beads.
! DEFAULTS COMMENT: ! Not included by default
! TYPE: REAL(DP)

#define WLC_P__AmplitudeExternalField 0.0_DP
! VARIABLE COMMENT: ! Amplitude of external fields if used.  For example, non-specific attraction to the nuclear boundary in the corresponding paper.
! DEFAULTS COMMENT: ! Not included by default
! TYPE: REAL(DP)

#define WLC_P__Amplitude2beadPotential 0.0_DP
! VARIABLE COMMENT: ! Amplitude of pairwise inter-bead potential.
! DEFAULTS COMMENT: !
! TYPE: REAL(DP)

#define WLC_P__EXCLUDE_NEIGHBORS_IN_POTENTIAL 2
! VARIABLE COMMENT: ! Number of neighbors on either side to exclude in pairwise interaction.
! DEFAULTS COMMENT: !
! TYPE: INTERGER

#define WLC_P__MAX_RXN_COORD 1.0_DP
! VARIABLE COMMENT: ! Maximum of umbrella sampled reaction coordinate.
! DEFAULTS COMMENT: !
! TYPE: REAL(DP)

#define WLC_P__MIN_RXN_COORD 0.0_DP
! VARIABLE COMMENT: ! Minimum of umbrella sampled reaction coordinate.
! DEFAULTS COMMENT: !
! TYPE: REAL(DP)

#define WLC_P__INTERNUCLEOSOME_ON .FALSE.
! VARIABLE COMMENT: ! whether or not to turn internucleosome attraction term on
! DEFAULTS COMMENT: ! toggle on or off
! TYPE: REAL(DP)

#define WLC_P__INTERNUCLEOSOME_ENERGY 4.0_DP 
! VARIABLE COMMENT: ! nucleosome-nucleosome attraction scale-term for energy
! DEFAULTS COMMENT: ! inspired from de pablos group, look at nucleosome.f03 for more details
! TYPE: REAL(DP)

! ----------------------------
!
!   Time / steps variables
!   For example, the total number of crankshaft moves is: 
!   NUMSAVEPOINTS*NREPLICAEXCHANGEPERSAVEPOINT*STEPSPEREXCHANGE*MOVESPERSTEP_CRANK_SHAFT
!
! ---------------------------

#define WLC_P__STEPSPEREXCHANGE 1000
! VARIABLE COMMENT: ! The number of steps between parallel tempering exchanges.
! DEFAULTS COMMENT: ! 100 steps between parallel tempering is pretty frequent
! TYPE: INTEGER

#define WLC_P__NREPLICAEXCHANGEPERSAVEPOINT 100
! VARIABLE COMMENT: !  Number of parallel tempering exchanges per save point.
! DEFAULTS COMMENT: ! make this large
! TYPE: INTEGER

#define WLC_P__NUMSAVEPOINTS 200
! VARIABLE COMMENT: ! Total number of save points before simulation end.
! DEFAULTS COMMENT: ! 200 total save points, i.e. 2000 steps per save point
! TYPE: INTEGER

#define WLC_P__STEPSPERSAVE 2000
! VARIABLE COMMENT: ! Steps per save point in browning dynamics code.
! DEFAULTS COMMENT: ! number of simulation steps to take
! TYPE: INTEGER

#define WLC_P__NNOINT 0
! VARIABLE COMMENT: ! Number of save points before turning on interactions.  This allows this simulation to equilibrate the chain energy prior to turn on interbead interactions.  To save time, densities will not be calculated during this period.
! DEFAULTS COMMENT: ! Zero by default.
! TYPE: INTEGER

#define WLC_P__N_KAP_ON 0
! VARIABLE COMMENT: ! Number of save points before turning on the kappa density interaction.  This allows this simulation to equilibrate the chain energy prior to turn on interbead interactions.
! DEFAULTS COMMENT: ! turn on compression energy immediately
! TYPE: INTEGER

#define WLC_P__N_CHI_ON 0
! VARIABLE COMMENT: ! Number of save points before turning on the Flory-Huggins interaction.  This allows this simulation to equilibrate the chain energy prior to turn on interbead interactions.
! DEFAULTS COMMENT: ! Turn on chi energy immediately
! TYPE: INTEGER

#define WLC_P__N_CHI_L2_ON 0
! VARIABLE COMMENT: ! Number of save points before turning on the quadrupole liquid crystal interaction.  This allows this simulation to equilibrate the chain energy prior to turn on interbead interactions.
! DEFAULTS COMMENT: ! Turn on chi energy immediately
! TYPE: INTEGER

#define WLC_P__N_EXTERNAL_ON 0
! VARIABLE COMMENT: ! Number of save points before turning on the external field.
! DEFAULTS COMMENT: ! turn on immediately
! TYPE: INTEGER

#define WLC_P__N_FIELD_ON 0
! VARIABLE COMMENT: ! Number of save points before turning on the applied density external field. (Rarely used)
! DEFAULTS COMMENT: ! turn on immediately
! TYPE: INTEGER

#define WLC_P__N_UMBRELLA_ON 0
! VARIABLE COMMENT: ! Number of save points before starting umbrella energy.
! DEFAULTS COMMENT: ! turn on immediately
! TYPE: INTEGER

#define WLC_P__NINITMCSTEPS 4000
! VARIABLE COMMENT: ! Number of mc steps before starting  brownian dynamics.
! DEFAULTS COMMENT: ! number of initializing mc steps. 1000s x num polymers is good
! TYPE: INTEGER

#define WLC_P__MOVESPERSTEP_CRANK_SHAFT 11
! VARIABLE COMMENT: ! How many times each move is done per step.
! DEFAULTS COMMENT: ! From experience.
! TYPE: INTEGER

#define WLC_P__MOVESPERSTEP_SLIDE_MOVE 25
! VARIABLE COMMENT: ! How many times each move is done per step.
! DEFAULTS COMMENT: ! From experience.
! TYPE: INTEGER

#define WLC_P__MOVESPERSTEP_PIVOT_MOVE 19
! VARIABLE COMMENT: ! How many times each move is done per step.
! DEFAULTS COMMENT: ! From experience.
! TYPE: INTEGER

#define WLC_P__MOVESPERSTEP_ROTATE_MOVE 45
! VARIABLE COMMENT: ! How many times each move is done per step.
! DEFAULTS COMMENT: ! From experience.
! TYPE: INTEGER

#define WLC_P__MOVESPERSTEP_FULL_CHAIN_ROTATION 5
! VARIABLE COMMENT: ! How many times each move is done per step.
! DEFAULTS COMMENT: ! From experience.
! TYPE: INTEGER

#define WLC_P__MOVESPERSTEP_FULL_CHAIN_SLIDE 6
! VARIABLE COMMENT: ! how many times each move is done per step
! DEFAULTS COMMENT: ! full chain slide 
! TYPE: INTEGER

#define WLC_P__MOVESPERSTEP_CHANGE_BINDING_STATE 111
! VARIABLE COMMENT: ! How many times each move is done per step.
! DEFAULTS COMMENT: ! From experience.
! TYPE: INTEGER

#define WLC_P__MOVESPERSTEP_CHAIN_FLIP 0
! VARIABLE COMMENT: ! How many times each move is done per step.
! DEFAULTS COMMENT: ! Not in use.
! TYPE: INTEGER

#define WLC_P__MOVESPERSTEP_CHAIN_EXCHANGE 5
! VARIABLE COMMENT: ! How many times each move is done per step.
! DEFAULTS COMMENT: ! From experience.
! TYPE: INTEGER

#define WLC_P__MOVESPERSTEP_REPTATION 10
! VARIABLE COMMENT: ! How many times each move is done per step.
! DEFAULTS COMMENT: ! reptation
! TYPE: INTEGER

#define WLC_P__MOVESPERSTEP_SUPER_REPTATION 57
! VARIABLE COMMENT: ! How many times each move is done per step.
! DEFAULTS COMMENT: ! super-reptation <default 57>
! TYPE: INTEGER

#define WLC_P__MOVESPERSTEP_SPIDER 3
! VARIABLE COMMENT: ! How many times each move is done per step.
! DEFAULTS COMMENT: ! From experience.
! TYPE: INTEGER

#define WLC_P__MOVESPERSTEP_NUCLEOSOME_SLIDE 60
! VARIABLE COMMENT: ! how many times each move is done per step
! DEFAULTS COMMENT: ! nucleosome slide <default 60> (just copying bead slide)
! TYPE: INTEGER

#define WLC_P__MOVESPERSTEP_NUCLEOSOME_BREATHE 10
! VARIABLE COMMENT: ! how many times each move is done per step
! DEFAULTS COMMENT: ! nucleosome breathe <default ?> 
! TYPE: INTEGER

#define WLC_P__PROBSINGLESWAP 0.5_dp
! VARIABLE COMMENT: ! For TWO_TAIL simulations, what fraction of the BINDING_STATE moves to change the binding state of only one tail.
! DEFAULTS COMMENT: ! 0< x  <1
! TYPE: LOGICAL

! --------------------------
!
!   Parallel Tempering settings (via MPI)
!
! --------------------------

#define WLC_P__PTON .FALSE.
! VARIABLE COMMENT: ! whether or not to parallel temper
! DEFAULTS COMMENT: ! use parallel if applicable
! TYPE: LOGICAL

#define WLC_P__ENSEMBLE_BIND .FALSE.
! VARIABLE COMMENT: ! set but not parrallel temered
! DEFAULTS COMMENT: !
! TYPE: LOGICAL

#define WLC_P__ENSEMBLE_METH .FALSE.
! VARIABLE COMMENT: ! set but not parrallel temered
! DEFAULTS COMMENT: !
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
! VARIABLE COMMENT: ! Parallel temper applied field strength.
! DEFAULTS COMMENT: ! No parallel temper.
! TYPE: LOGICAL

#define WLC_P__PT_KAP .FALSE.
! VARIABLE COMMENT: ! Parallel temper kappa (incompressibility) energy.
! DEFAULTS COMMENT: ! don't parallel temper kap by default
! TYPE: LOGICAL
#define WLC_P__PT_MU .FALSE.
! VARIABLE COMMENT: ! Parallel temper chemical potential of HP1.
! DEFAULTS COMMENT: ! don't parallel temper mu by default
! TYPE: LOGICAL

#define WLC_P__PT_COUPLE .FALSE.
! VARIABLE COMMENT:  ! Parallel temper HP1 coupling energy.
! DEFAULTS COMMENT: ! don't parallel temper HP1 binding by default
! TYPE: LOGICAL
#define WLC_P__PT_MAIERSAUPE .FALSE.
! VARIABLE COMMENT: ! Parallel temper Maier Saupe liquid crystal interaction.
! DEFAULTS COMMENT: ! don't parallel temper maier Saupe by default
! TYPE: LOGICAL

#define WLC_P__PT_AEF .FALSE.
! VARIABLE COMMENT: ! Parallel temper over external applied field strength.
! DEFAULTS COMMENT: ! don't parallel temper h by default
! TYPE: LOGICAL
#define WLC_P__PT_A2B .FALSE.
! VARIABLE COMMENT: ! Parallel temper two body interaction strength.
! DEFAULTS COMMENT: ! don't parallel temper.
! TYPE: LOGICAL

#define WLC_P__PT_STERICS .FALSE.
! VARIABLE COMMENT: ! parallel temper steric collision penalty
! DEFAULTS COMMENT: !
! TYPE: LOGICAL

#define WLC_P__PT_INTERNUCLEOSOME .FALSE.
! VARIABLE COMMENT: ! parallel temper internucleosome attraction term
! DEFAULTS COMMENT: !
! TYPE: LOGICAL

#define WLC_P__NREPADAPT 1000
! VARIABLE COMMENT: ! Number of exchange attempts between adapting the spacing of the replicas (only for adaptive parallel tempering spacing).
! DEFAULTS COMMENT: ! 1000 exchange attempts between adaptations
! TYPE: INTEGER

#define WLC_P__LOWERREPEXE 0.09_DP
! VARIABLE COMMENT: ! Minimum acceptable exchange success rate (only for adaptive parallel tempering spacing).
! DEFAULTS COMMENT: ! Exchange still happens enough to obtain mixing.
! TYPE: REAL(DP)

#define WLC_P__UPPERREPEXE 0.18_DP
! VARIABLE COMMENT: ! Maximum acceptable exchange success rate (only for adaptive parallel tempering spacing).
! DEFAULTS COMMENT: ! Appreciable value but not wasting too much processor power.
! TYPE: REAL(DP)

#define WLC_P__LOWERCOFRAIL 0.005_DP
! VARIABLE COMMENT: ! Minimum acceptable coefficient spacing of replicas (only for adaptive parallel tempering spacing).
! DEFAULTS COMMENT: ! From experience.
! TYPE: REAL(DP)

#define WLC_P__UPPERCOFRAIL 0.1_DP
! VARIABLE COMMENT: ! Maximum acceptable coefficient spacing of replicas (only for adaptive parallel tempering spacing).
! DEFAULTS COMMENT: ! From experience.
! TYPE: REAL(DP)

#define WLC_P__INDSTARTREPADAPT -1
! VARIABLE COMMENT: ! Time point to start replica adapting (only for adaptive parallel tempering spacing).
! DEFAULTS COMMENT: ! Never adapt
! TYPE: INTEGER

#define WLC_P__INDENDREPADAPT -1
! VARIABLE COMMENT: ! Time point to end replica adapting (only for adaptive parallel tempering spacing).
! DEFAULTS COMMENT: ! Never adapt
! TYPE: INTEGER

#define WLC_P__REPANNEALSPEED 0.01_DP
! VARIABLE COMMENT: ! Maximum allowable change on parallel tempering coefficient at any one time (only for adaptive parallel tempering spacing).  This allows for annealing.
! DEFAULTS COMMENT: ! Small
! TYPE: REAL(DP)

#define WLC_P__REPLICABOUNDS .FALSE.
! VARIABLE COMMENT: ! Keep coefficient between 0 and 1.
! DEFAULTS COMMENT: ! Not in use.
! TYPE: LOGICAL

#define WLC_P__INITIAL_MAX_S NAN
! VARIABLE COMMENT: ! The s variable parameterized location along parallel temper path.  This sets it’s initial range.
! DEFAULTS COMMENT: ! You should set this variable if you intend to use it
! TYPE: REAL(DP)

! --------------------------
!
!   Available MC Moves
!
! --------------------------

#define WLC_P__MOVEON_CRANK_SHAFT 1
! VARIABLE COMMENT: ! Is the move active.
! DEFAULTS COMMENT: ! crank-shaft move
! TYPE: INTEGER

#define WLC_P__MOVEON_SLIDE_MOVE 1
! VARIABLE COMMENT: ! Is the move active.
! DEFAULTS COMMENT: ! slide move
! TYPE: INTEGER

#define WLC_P__MOVEON_PIVOT_MOVE 1
! VARIABLE COMMENT: ! Is the move active.
! DEFAULTS COMMENT: ! pivot move
! TYPE: INTEGER

#define WLC_P__MOVEON_ROTATE_MOVE 1
! VARIABLE COMMENT: ! Is the move active.
! DEFAULTS COMMENT: ! rotate move
! TYPE: INTEGER

#define WLC_P__MOVEON_FULL_CHAIN_ROTATION 1
! VARIABLE COMMENT: ! Is the move active.
! DEFAULTS COMMENT: ! full chain rotation
! TYPE: INTEGER

#define WLC_P__MOVEON_FULL_CHAIN_SLIDE 1
! VARIABLE COMMENT: ! Is the move active.
! DEFAULTS COMMENT: ! full chain slide
! TYPE: INTEGER

#define WLC_P__MOVEON_CHANGE_BINDING_STATE 0
! VARIABLE COMMENT: ! Is the move active.
! DEFAULTS COMMENT: ! Change in Binding state
! TYPE: INTEGER

#define WLC_P__MOVEON_CHAIN_FLIP 0
! VARIABLE COMMENT: ! Is the move active.
! DEFAULTS COMMENT: ! Chain flip ! TOdo not working
! TYPE: INTEGER

#define WLC_P__MOVEON_CHAIN_EXCHANGE 0
! VARIABLE COMMENT: ! Is the move active.
! DEFAULTS COMMENT: ! Chain exchange
! TYPE: INTEGER

#define WLC_P__MOVEON_REPTATION 1
! VARIABLE COMMENT: ! Is the move active.
! DEFAULTS COMMENT: ! Reptation
! TYPE: INTEGER

#define WLC_P__MOVEON_SUPER_REPTATION 0
! VARIABLE COMMENT: ! Is the move active.
! DEFAULTS COMMENT: ! Super Reptation
! TYPE: INTEGER

#define WLC_P__MOVEON_SPIDER 0
! VARIABLE COMMENT: ! Is the move active.
! DEFAULTS COMMENT: ! Super Reptation
! TYPE: INTEGER

#define WLC_P__MOVEON_NUCLEOSOME_SLIDE 0
! VARIABLE COMMENT: ! Is the move active
! DEFAULTS COMMENT: ! nucleosome slide
! TYPE: INTEGER

#define WLC_P__MOVEON_NUCLEOSOME_BREATHE 0
! VARIABLE COMMENT: ! Is the move active
! DEFAULTS COMMENT: ! nucleosome angle
! TYPE: INTEGER

! ----------------------------
!
!   MC Optimization
!
! ---------------------------

#define WLC_P__PDESIRE_CRANK_SHAFT 0.5_DP
! VARIABLE COMMENT: ! Desired Monte-Carlo acceptance probability for move.
! DEFAULTS COMMENT: ! Target
! TYPE: REAL(DP)

#define WLC_P__PDESIRE_SLIDE_MOVE 0.5_DP
! VARIABLE COMMENT: ! Desired Monte-Carlo acceptance probability for move.
! DEFAULTS COMMENT: ! Target
! TYPE: REAL(DP)

#define WLC_P__PDESIRE_PIVOT_MOVE 0.5_DP
! VARIABLE COMMENT: ! Desired Monte-Carlo acceptance probability for move.
! DEFAULTS COMMENT: ! Target
! TYPE: REAL(DP)

#define WLC_P__PDESIRE_ROTATE_MOVE 0.5_DP
! VARIABLE COMMENT: ! Desired Monte-Carlo acceptance probability for move.
! DEFAULTS COMMENT: ! Target
! TYPE: REAL(DP)

#define WLC_P__PDESIRE_FULL_CHAIN_ROTATION 0.5_DP
! VARIABLE COMMENT: ! Desired Monte-Carlo acceptance probability for move.
! DEFAULTS COMMENT: ! Target
! TYPE: REAL(DP)

#define WLC_P__PDESIRE_FULL_CHAIN_SLIDE 0.5_DP
! VARIABLE COMMENT: ! Desired Monte-Carlo acceptance probability for move.
! DEFAULTS COMMENT: ! Target
! TYPE: REAL(DP)

#define WLC_P__PDESIRE_CHANGE_BINDING_STATE 0.5_DP
! VARIABLE COMMENT: ! Desired Monte-Carlo acceptance probability for move.
! DEFAULTS COMMENT: ! Target
! TYPE: REAL(DP)

#define WLC_P__PDESIRE_CHAIN_FLIP 0.5_DP
! VARIABLE COMMENT: ! Desired Monte-Carlo acceptance probability for move.
! DEFAULTS COMMENT: ! Target
! TYPE: REAL(DP)

#define WLC_P__PDESIRE_CHAIN_EXCHANGE 0.5_DP
! VARIABLE COMMENT: ! Desired Monte-Carlo acceptance probability for move.
! DEFAULTS COMMENT: ! Target
! TYPE: REAL(DP)

#define WLC_P__PDESIRE_REPTATION 0.5_DP
! VARIABLE COMMENT: ! Desired Monte-Carlo acceptance probability for move.
! DEFAULTS COMMENT: ! Target
! TYPE: REAL(DP)

#define WLC_P__PDESIRE_SUPER_REPTATION 0.5_DP
! VARIABLE COMMENT: ! Desired Monte-Carlo acceptance probability for move.
! DEFAULTS COMMENT: ! Target
! TYPE: REAL(DP)

#define WLC_P__PDESIRE_SPIDER 0.25_DP
! VARIABLE COMMENT: ! Desired Monte-Carlo acceptance probability for move.
! DEFAULTS COMMENT: ! Target
! TYPE: REAL(DP)

#define WLC_P__PDESIRE_NUCLEOSOME_SLIDE 0.5_DP
! VARIABLE COMMENT: ! desired hit rate
! DEFAULTS COMMENT: ! Target

#define WLC_P__PDESIRE_NUCLEOSOME_BREATHE 0.5_DP
! VARIABLE COMMENT: ! desired hit rate
! DEFAULTS COMMENT: ! Target

#define WLC_P__MAXWINDOW_CRANK_SHAFT DBLE(MIN(150,WLC_P__NB))
! VARIABLE COMMENT: ! Maximum mean window size (i.e. how many beads are moved) for this move.
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__MAXWINDOW_SLIDE_MOVE DBLE(MIN(150,WLC_P__NB))
! VARIABLE COMMENT: ! Maximum mean window size (i.e. how many beads are moved) for this move.
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__MAXWINDOW_PIVOT_MOVE DBLE(MIN(150,INT(WLC_P__NB/2.0)))
! VARIABLE COMMENT: ! Max Size of mean window size (i.e. how many beads are moved) for this move.  
! DEFAULTS COMMENT: ! From expereince.
! TYPE: REAL(DP)

#define WLC_P__MAXWINDOW_ROTATE_MOVE NAN
! VARIABLE COMMENT: ! Maximum mean window size (i.e. how many beads are moved) for this move.
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__MAXWINDOW_FULL_CHAIN_ROTATION NAN
! VARIABLE COMMENT: ! Maximum mean window size (i.e. how many beads are moved) for this move.
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__MAXWINDOW_FULL_CHAIN_SLIDE NAN
! VARIABLE COMMENT: ! Maximum mean window size (i.e. how many beads are moved) for this move.
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__MAXWINDOW_CHANGE_BINDING_STATE DBLE(MIN(4,WLC_P__NB))
! VARIABLE COMMENT: ! Maximum mean window size (i.e. how many beads are moved) for this move.  
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__MAXWINDOW_CHAIN_FLIP NAN
! VARIABLE COMMENT: ! Maximum mean window size (i.e. how many beads are moved) for this move.
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__MAXWINDOW_CHAIN_EXCHANGE NAN
! VARIABLE COMMENT: ! Maximum mean window size (i.e. how many beads are moved) for this move.
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__MAXWINDOW_REPTATION NAN
! VARIABLE COMMENT: ! Maximum mean window size (i.e. how many beads are moved) for this move.
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__MAXWINDOW_SUPER_REPTATION NAN
! VARIABLE COMMENT: ! Maximum mean window size (i.e. how many beads are moved) for this move.
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__MINWINDOW_CRANK_SHAFT 30_DP
! VARIABLE COMMENT: ! Minimum mean window size (i.e. how many beads are moved) for this move. This setting is importnat! Depends on problem.
! DEFAULTS COMMENT: ! From Expereince.
! TYPE: REAL(DP)

#define WLC_P__MINWINDOW_SLIDE_MOVE 10_DP
! VARIABLE COMMENT: ! Minimum mean window size (i.e. how many beads are moved) for this move.
! DEFAULTS COMMENT: ! If this is too small equllibration will be very slow.
! TYPE: REAL(DP)

#define WLC_P__MINWINDOW_PIVOT_MOVE 50_DP
! VARIABLE COMMENT: ! Minimum mean window size (i.e. how many beads are moved) for this move.
! DEFAULTS COMMENT: ! From experenice.
! TYPE: REAL(DP)

#define WLC_P__MINWINDOW_ROTATE_MOVE NAN
! VARIABLE COMMENT: ! Minimum mean window size (i.e. how many beads are moved) for this move.
! DEFAULTS COMMENT: ! If this is too small equllibration will be very slow.
! TYPE: REAL(DP)

#define WLC_P__MINWINDOW_FULL_CHAIN_ROTATION NAN
! VARIABLE COMMENT: ! Minimum mean window size (i.e. how many beads are moved) for this move.
! DEFAULTS COMMENT: ! Not in use.
! TYPE: REAL(DP)

#define WLC_P__MINWINDOW_FULL_CHAIN_SLIDE NAN
! VARIABLE COMMENT: ! Minimum mean window size (i.e. how many beads are moved) for this move.
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__MINWINDOW_CHANGE_BINDING_STATE NAN
! VARIABLE COMMENT: ! Minimum mean window size (i.e. how many beads are moved) for this move.
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__MINWINDOW_CHAIN_FLIP NAN
! VARIABLE COMMENT: ! Minimum mean window size (i.e. how many beads are moved) for this move.
! DEFAULTS COMMENT: ! Not in use.
! TYPE: REAL(DP)

#define WLC_P__MINWINDOW_CHAIN_EXCHANGE NAN
! VARIABLE COMMENT: ! Minimum mean window size (i.e. how many beads are moved) for this move.
! DEFAULTS COMMENT: ! Not in use.
! TYPE: REAL(DP)

#define WLC_P__MINWINDOW_REPTATION NAN
! VARIABLE COMMENT: ! Minimum mean window size (i.e. how many beads are moved) for this move.
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__MINWINDOW_SUPER_REPTATION NAN
! VARIABLE COMMENT: ! Minimum acceptable amplitude range for this move.
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__MINAMP_CRANK_SHAFT 0.1_DP*PI
! VARIABLE COMMENT: ! Minimum acceptable amplitude range for this move.
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__MINAMP_SLIDE_MOVE 0.2_DP*WLC_P__L0
! VARIABLE COMMENT: ! Minimum acceptable amplitude range for this move.
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__MINAMP_PIVOT_MOVE 0.2_DP*PI
! VARIABLE COMMENT: ! Minimum acceptable amplitude range for this move.
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__MINAMP_ROTATE_MOVE 0.2_DP*PI
! VARIABLE COMMENT: ! Minimum acceptable amplitude range for this move.
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__MINAMP_FULL_CHAIN_ROTATION 0.05_DP*PI
! VARIABLE COMMENT: ! minium amplitude
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__MINAMP_FULL_CHAIN_SLIDE 0.2_DP*WLC_P__L0
! VARIABLE COMMENT: ! Minimum acceptable amplitude range for this move.
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__MINAMP_CHANGE_BINDING_STATE NAN
! VARIABLE COMMENT: ! Minimum acceptable amplitude range for this move.
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__MINAMP_CHAIN_FLIP NAN
! VARIABLE COMMENT: ! Minimum acceptable amplitude range for this move.
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__MINAMP_CHAIN_EXCHANGE NAN
! VARIABLE COMMENT: ! Minimum acceptable amplitude range for this move.
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__MINAMP_REPTATION NAN
! VARIABLE COMMENT: ! Minimum acceptable amplitude range for this move.
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__MINAMP_SUPER_REPTATION NAN
! VARIABLE COMMENT: ! Minimum acceptable amplitude range for this move.
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__MINAMP_SPIDER 1.0_dp*WLC_P__L0
! VARIABLE COMMENT: ! Minimum acceptable amplitude range for this move.
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__MINAMP_NUCLEOSOME_SLIDE NAN
! VARIABLE COMMENT: ! minium amplitude
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__MINAMP_NUCLEOSOME_BREATHE NAN
! VARIABLE COMMENT: ! minium amplitude
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__MAXAMP_CRANK_SHAFT 1.0_DP*PI
! VARIABLE COMMENT: ! Maximum acceptable amplitude range for this move.
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__MAXAMP_SLIDE_MOVE 1.0_DP*WLC_P__L0
! VARIABLE COMMENT: ! Maximum acceptable amplitude range for this move.
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__MAXAMP_PIVOT_MOVE 1.0_DP*PI
! VARIABLE COMMENT: ! Maximum acceptable amplitude range for this move.
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__MAXAMP_ROTATE_MOVE 1.0_DP*PI
! VARIABLE COMMENT: ! Maximum acceptable amplitude range for this move.
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__MAXAMP_FULL_CHAIN_ROTATION 1.0_DP*PI
! VARIABLE COMMENT: ! Maximum acceptable amplitude range for this move.
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__MAXAMP_FULL_CHAIN_SLIDE 0.1*WLC_P__L
! VARIABLE COMMENT: ! Maximum acceptable amplitude range for this move.
! DEFAULTS COMMENT: ! a fraction of the change length
! TYPE: REAL(DP)

#define WLC_P__MAXAMP_CHANGE_BINDING_STATE NAN
! VARIABLE COMMENT: ! Maximum acceptable amplitude range for this move.
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__MAXAMP_CHAIN_FLIP NAN
! VARIABLE COMMENT: ! Maximum acceptable amplitude range for this move.
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__MAXAMP_CHAIN_EXCHANGE NAN
! VARIABLE COMMENT: ! Maximum acceptable amplitude range for this move.
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__MAXAMP_REPTATION NAN
! VARIABLE COMMENT: ! Maximum acceptable amplitude range for this move.
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__MAXAMP_SUPER_REPTATION NAN
! VARIABLE COMMENT: ! Maximum acceptable amplitude range for this move.
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__MAXAMP_SPIDER 5.0_dp*WLC_P__L0
! VARIABLE COMMENT: ! Maximum acceptable amplitude range for this move.
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__MAXAMP_NUCLEOSOME_SLIDE NAN
! VARIABLE COMMENT: ! maximum amplitude
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__MAXAMP_NUCLEOSOME_BREATHE NAN
! VARIABLE COMMENT: ! maximum amplitude
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__WINTARGET_CRANK_SHAFT 8.0_DP
! VARIABLE COMMENT: ! Target for ratio of average window to amplitude range for this move.
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__WINTARGET_SLIDE_MOVE 8.0_DP
! VARIABLE COMMENT: ! Target for ratio of average window to amplitude range for this move.
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__WINTARGET_PIVOT_MOVE 8.0_DP
! VARIABLE COMMENT: ! Target for ratio of average window to amplitude range for this move.
! DEFAULTS COMMENT: None
! TYPE: REAL(DP)

#define WLC_P__WINTARGET_ROTATE_MOVE NAN
! VARIABLE COMMENT: ! Target for ratio of average window to amplitude range for this move.
! DEFAULTS COMMENT: ! Not in use.
! TYPE: REAL(DP)

#define WLC_P__WINTARGET_FULL_CHAIN_ROTATION NAN
! VARIABLE COMMENT: ! Target for ratio of average window to amplitude range for this move.
! DEFAULTS COMMENT: ! No window.
! TYPE: REAL(DP)

#define WLC_P__WINTARGET_FULL_CHAIN_SLIDE NAN
! VARIABLE COMMENT: ! Target for ratio of average window to amplitude range for this move.
! DEFAULTS COMMENT: ! No window.
! TYPE: REAL(DP)

#define WLC_P__WINTARGET_CHANGE_BINDING_STATE NAN
! VARIABLE COMMENT: ! Target for ratio of average window to amplitude range for this move.
! DEFAULTS COMMENT: ! No amplitude.
! TYPE: REAL(DP)

#define WLC_P__WINTARGET_CHAIN_FLIP NAN
! VARIABLE COMMENT: ! Target for ratio of average window to amplitude range for this move.
! DEFAULTS COMMENT: ! Not in use.
! TYPE: REAL(DP)

#define WLC_P__WINTARGET_CHAIN_EXCHANGE NAN
! VARIABLE COMMENT: ! Target for ratio of average window to amplitude range for this move.
! DEFAULTS COMMENT: ! No window or amplitude.
! TYPE: REAL(DP)

#define WLC_P__WINTARGET_REPTATION NAN
! VARIABLE COMMENT: ! Target for ratio of average window to amplitude range for this move.
! DEFAULTS COMMENT: ! No window or amplitude.
! TYPE: REAL(DP)

#define WLC_P__WINTARGET_SUPER_REPTATION NAN
! VARIABLE COMMENT: ! Target for ratio of average window to amplitude range for this move.
! DEFAULTS COMMENT: ! No window.
! TYPE: REAL(DP)

#define WLC_P__NADAPT_CRANK_SHAFT 1000
! VARIABLE COMMENT: ! Number of steps (not moves) between adaptations for this move.
! DEFAULTS COMMENT: ! adapt after at most 1000 steps
! TYPE: INTEGER

#define WLC_P__NADAPT_SLIDE_MOVE 1000
! VARIABLE COMMENT: ! Number of steps (not moves) between adaptations for this move.
! DEFAULTS COMMENT: ! adapt after at most 1000 steps
! TYPE: INTEGER

#define WLC_P__NADAPT_PIVOT_MOVE 1000
! VARIABLE COMMENT: ! Number of steps (not moves) between adaptations for this move.
! DEFAULTS COMMENT: ! adapt after at most 1000 steps
! TYPE: INTEGER

#define WLC_P__NADAPT_ROTATE_MOVE 1000
! VARIABLE COMMENT: ! Number of steps (not moves) between adaptations for this move.
! DEFAULTS COMMENT: ! adapt after at most 1000 steps
! TYPE: INTEGER

#define WLC_P__NADAPT_FULL_CHAIN_ROTATION 1000
! VARIABLE COMMENT: ! Number of steps (not moves) between adaptations for this move.
! DEFAULTS COMMENT: ! adapt after at most 1000 steps
! TYPE: INTEGER

#define WLC_P__NADAPT_FULL_CHAIN_SLIDE 1000
! VARIABLE COMMENT: ! Number of steps (not moves) between adaptations for this move.
! DEFAULTS COMMENT: ! adapt after at most 1000 steps
! TYPE: INTEGER

#define WLC_P__NADAPT_CHANGE_BINDING_STATE 1000
! VARIABLE COMMENT: ! Number of steps (not moves) between adaptations for this move.
! DEFAULTS COMMENT: ! adapt after at most 1000 steps
! TYPE: INTEGER

#define WLC_P__NADAPT_CHAIN_FLIP 1000
! VARIABLE COMMENT: ! Number of steps (not moves) between adaptations for this move.
! DEFAULTS COMMENT: ! adapt after at most 1000 steps
! TYPE: INTEGER

#define WLC_P__NADAPT_CHAIN_EXCHANGE 1000
! VARIABLE COMMENT: ! Number of steps (not moves) between adaptations for this move.
! DEFAULTS COMMENT: ! adapt after at most 1000 steps
! TYPE: INTEGER

#define WLC_P__NADAPT_REPTATION 1000
! VARIABLE COMMENT: ! Number of steps (not moves) between adaptations for this move.t
! DEFAULTS COMMENT: ! adapt after at most 1000 steps
! TYPE: INTEGER

#define WLC_P__NADAPT_SUPER_REPTATION 1000
! VARIABLE COMMENT: ! Number of steps (not moves) between adaptations for this move.
! DEFAULTS COMMENT: ! adapt after at most 1000 steps
! TYPE: INTEGER

#define WLC_P__NADAPT_SPIDER 1000
! VARIABLE COMMENT: ! Number of steps (not moves) between adaptations for this move.
! DEFAULTS COMMENT: ! adapt after at most 1000 steps
! TYPE: INTEGER

#define WLC_P__NADAPT_NUCLEOSOME_SLIDE 1000
! VARIABLE COMMENT: ! Nunber of steps between adapt
! DEFAULTS COMMENT: ! adapt after at most 1000 steps
! TYPE: INTEGER

#define WLC_P__NADAPT_NUCLEOSOME_BREATHE 1000
! VARIABLE COMMENT: ! Nunber of steps between adapt
! DEFAULTS COMMENT: ! adapt after at most 1000 steps
! TYPE: INTEGER

#define WLC_P__MIN_ACCEPT 0.05_DP
! VARIABLE COMMENT: ! Threshold for deciding to usually not use a move (if applicable).
! DEFAULTS COMMENT: ! if a move succeeds < 5% of the time, start using it only every reduce_move cycles
! TYPE: REAL(DP)

#define WLC_P__REDUCE_MOVE 1
! VARIABLE COMMENT: ! Factor by which to reduce the usage of a move if it goes below min_accept success rate.  Set to 1 for no reduction.
! DEFAULTS COMMENT: ! No effect.
! TYPE: INTEGER

#define WLC_P__WINTYPE 1
! VARIABLE COMMENT: ! Distribution type of segment size in crankshaft move (uniform = 0, exponential = 1).  Exponential recommended.
! DEFAULTS COMMENT: ! exponential fragment sizes mix better
! TYPE: INTEGER

#define WLC_P__PROB_BIND_RESPECTING_MOVE 0.5_DP
! VARIABLE COMMENT: ! Probability that we make sure the move doesn’t break loop (for network simulations).  Set to 1.0_DP of bonds are enforced strictly.
! DEFAULTS COMMENT: ! about half, only used if explicit binding turned on
! TYPE: REAL(DP)

#define WLC_P__SPIDER_LEG_LENGTH 5
! VARIABLE COMMENT: ! Number of beads from hip-knee (and knee to toe) for spider move.
! DEFAULTS COMMENT: ! several but not too many
! TYPE: Integer

#define WLC_P__MAX_SPIDER_LEG_LENGTH 25
! VARIABLE COMMENT: ! Max number of beads from knee-toe for spider move.
! DEFAULTS COMMENT:
! TYPE: Integer

#define WLC_P__SPIDER_WINDOW 30.0_DP
! VARIABLE COMMENT: ! Mean of segment size distribution for spider move.
! DEFAULTS COMMENT:
! TYPE: REAL(DP)

#define WLC_P__PROBABILITY_PRECALC_SPIDER 0.0_DP
! VARIABLE COMMENT: ! Probability that we use a precalculated set of beads for a spider move.  Set to >1.0 for all precalculated.
! DEFAULTS COMMENT: !
! TYPE: REAL(DP)

#define WLC_P__NO_CROSSING_CUTOFF 0.5_DP
! VARIABLE COMMENT: ! Deviation from initial linking number allowed when WLC_P__NO_SELF_CROSSING is true.
! DEFAULTS COMMENT: ! allow linking number to change by 0.5 for each move 
! TYPE: REAL(DP)
