module KEYS
  ! keyword parameters that are globally used in many different places in the code
  implicit none

   ! -------- Other ---------------
  CHARACTER*100 :: ACTION
  integer :: RNGSEED
  LOGICAL :: VERBOSE

  ! ----------------------
  ! chain geometry and energetics
  ! ---------------------
  LOGICAL :: SHEARABLE, STRETCHABLE, COUPLED, FinITEXT, useSTERICS
  real(dp) :: LP, LS, EC, EPERP, EPAR, GAM, FinITSHEAR, CONSTMOD,STERRAD
  integer :: STARTNPT, MAXNPT,STERSKIP
  real(dp) :: STERMOD
  integer :: NEDGESEG
  real(dp) :: EDGELS, EDGELP, EDGEGAM, EDGEEPAR, EDGEEPERP,EDGEEC


  ! --------- Output / input-------------
  CHARACTER*100 :: outFILE, SNAPSHOTFILE, RESTARTFILE, STRESSFILE
  LOGICAL :: DUMPSNAPSHOTS, RESTART, APPENDSNAPSHOTS, STARTEQUIL
  integer :: SNAPSHOTEVERY, SKIPREAD, EQUILSAMPLETYPE
  LOGICAL :: EQUILBEADROD, PARAMFROMSNAPSHOT
  real(dp) :: STARTEQUILLP

  ! Monte Carlo
  integer :: MCPRinTFREQ, MCTOTSTEPS, MCinITSTEPS, MCSTATSTEPS,MCoutPUTFREQ, ADJUSTEVERY
  real(dp) :: FACCTARGET, FACCTOL, ADJUSTSCL, inITRANGE(4)
  LOGICAL :: doREDISC, outPUTBEADWEIGHT
  real(dp) :: MinSEGLEN, MAXSEGLEN
  integer :: inTUWEIGHTNPT, inTRWEIGHTNPT

  ! brownian dynamics
  real(dp) :: FRICTR, FRICTU, DELTSCL
  LOGICAL :: FRICTPERLEN
  integer :: BDSTEPS, RUNGEKUTTA
  LOGICAL :: LOGRTERM,FIXBEAD1,FIXBEADMID
  real(dp) :: BDPRinTEVERY
  LOGICAL :: BDPRinTLOG, GAUSSIANCHAin, useBDENERGY, doBROWN, usePSEUdoforCE
  real(dp) :: BRCRELAX
  LOGICAL :: TRACKLOOPinG
  CHARACTER*100 :: LOOPFILE
  real(dp) :: LOOPRAD

  ! obstacles
  real(dp) :: RAdoB, MOdoB, FRICTOB

  real(dp) :: mu

  ! ------------ groups of chains ------------
  integer :: NforCE, NCONNECT, NCHAin
  integer, PARAMETER :: MAXNforCE = 100
  integer, PARAMETER :: MAXNCONNECT = 1000
  real(dp) :: forCE(MAXNforCE,3)
  integer :: forCEBEAD(MAXNforCE,2)
  integer :: CONNECTIONS(MAXNCONNECT,4)
  LOGICAL :: CONNECTPOS, CONNECTUVEC, SQUARELATTICE, FIXCONNECT
  real(dp) :: CONPOSMOD, CONUVECMOD
  integer :: TRACKDIST(4) ! track distance between which two points (bead&chain)?

  LOGICAL :: DIAMONDLATTICE ! set up a diamon lattice
  integer :: NDIAMOND(2) ! number of diamonds across and down
  real(dp) :: WIDTHDIAMOND ! Width of each diamond initially
  integer :: LENDIAMOND ! number of chain segments along diamond side

  ! for each fixed bead list: bead, chain, fix position, fix orientation
  ! (0 for not, >0 for yes)
  integer, PARAMETER :: MAXFIXBEAD = 100
  integer :: NFIXBEAD
  integer :: FIXBEAD(MAXFIXBEAD,4)
  ! fix all bead positions for top/bottom and or side boundaries
  ! 1) top/bottom fixed 2) sides fixed
  ! 3) fix positions 4) fix orientations
  LOGICAL :: FIXBOUNDARY(4)
  ! force a shear deformation
  LOGICAL :: SETSHEAR
  real(dp) :: SHEARGAMMA
  ! start with collapsed structure
  LOGICAL :: STARTCOLLAPSE

  LOGICAL :: doLOCALMOVES
END module KEYS
