!   --------------------------------------------------------------
!
!    This module if defines two strucutres.  MCvar is for simulation
!    parameters of a fixed size.  MCdata is for allocateable variables.
!    These strucutres contain most of the variables needed for a MC
!    simulation of a wlc polymer with beads interacting via a field.
!
!        By Quinn MacPherson ~ Spring 2016
!
!   --------------------------------------------------------------
Module simMod
    use setPrecision
    IMPLICIT NONE
    INTEGER, Parameter :: NmoveTypes = 10 ! ******* YOU MAY NEED TO CHAGE THIS ***

  Type MCvar   ! Structure for simulation variables of known size
!     Simulation parameters
    INTEGER NT                ! Total number of beads  NT=NP*N*G
    INTEGER NB                ! Number of beads in a polymer NB=N*G
    INTEGER N                 ! Number of monomers in a polymer
    INTEGER G                 ! Beads per monomer
    INTEGER NP                ! Number of polymers
    DOUBLE PRECISION LBOX(3)  ! Box length (approximate)
    DOUBLE PRECISION DEL      ! Discretization size (approximate)
    DOUBLE PRECISION L0       ! Equilibrium segment length
    DOUBLE PRECISION V        ! Bead volume
    DOUBLE PRECISION FA       ! Fraction of A beads
    DOUBLE PRECISION LAM      ! Chemical correlation parameter
    DOUBLE PRECISION EPS      ! Elasticity l0/(2lp)
    DOUBLE PRECISION CHI      ! Chi parameter value (solvent-polymer)        
    DOUBLE PRECISION KAP      ! Incompressibility parameter
    DOUBLE PRECISION h_A      ! fild strength
    Double Precision k_field  ! wave vector of template field
    DOUBLE PRECISION Fpoly    ! Fraction Polymer
    DOUBLE PRECISION EU       ! Energy of binding for unmethalated
    DOUBLE PRECISION EM       ! Energy of binding for methalated
    DOUBLE PRECISION HP1_Bind ! Energy of binding of HP1 to eachother
    DOUBLE PRECISION F_METH   ! Fraction methalated is using option 2
    DOUBLE PRECISION LAM_METH ! eigenvalue of methalation setup
    DOUBLE PRECISION mu       ! chemical potential of HP1
    INTEGER NBIN     ! Number of bins
    INTEGER NBINX(3) ! Number of bin on an edge
    DOUBLE PRECISION PARA(10) ! Parameters for sswlc
        ! EB, EPAR, EPERP, GAM, ETA, ...


!   Monte Carlo Variables (for adaptation)
    INTEGER moveTypes
    DOUBLE PRECISION MCAMP(NmoveTypes) ! Amplitude of random change
    DOUBLE PRECISION MinAMP(NmoveTypes) ! minium amplitude
    DOUBLE PRECISION MaxAMP(NmoveTypes) ! maximum amplitude
    INTEGER MOVEON(NmoveTypes)         ! Is the move active
    DOUBLE PRECISION WINDOW(NmoveTypes)         ! Size of window for bead selection
    Double precision MAXWINDOW(NmoveTypes)         ! Max Size of window for bead selection
    double precision MINWINDOW(NmoveTypes)         ! Min Size of window for bead selection
    DOUBLE PRECISION winTarget(NmoveTypes) ! target for ratio of window to anmplitude
    INTEGER SUCCESS(NmoveTypes)        ! Number of successes
    DOUBLE PRECISION PDesire(NmoveTypes) ! desired hit rate     
    DOUBLE PRECISION PHit(NmoveTypes) ! hit rate 
    ! see Timing variables for NADAPT(NmoveTypes)

!   Energys
    !DOUBLE PRECISION Eint     ! running Eint
    DOUBLE PRECISION EELAS(3) ! Elastic force
    DOUBLE PRECISION ECHI     ! CHI energy
    DOUBLE PRECISION EKAP     ! KAP energy
    DOUBLE PRECISION ECouple  ! Coupling 
    DOUBLE PRECISION EBind    ! binding energy
    DOUBLE PRECISION EField   ! Field energy

!   Congigate Energy variables (needed to avoid NaN when cof-> 0 in rep exchange)
    DOUBLE PRECISION x_Chi,   dx_Chi
    DOUBLE PRECISION x_Couple,dx_Couple
    DOUBLE PRECISION x_Kap,   dx_Kap
    DOUBLE PRECISION x_Field, dx_Field
    DOUBLE PRECISION x_Mu,    dx_Mu

!   Move Variables 
    DOUBLE PRECISION DEELAS(3)   ! Change in bending energy
!    DOUBLE PRECISION DEINT    ! Change in self energy
    DOUBLE PRECISION DECouple ! Coupling energy
    DOUBLE PRECISION DEChi    ! chi interaction energy
    DOUBLE PRECISION DEKap    ! compression energy
    DOUBLE PRECISION DEBind   ! Change in binding energy
    Double Precision DEField  ! Change in field energy
    DOUBLE PRECISION ECon     ! Confinement Energy
    INTEGER NPHI  ! NUMBER o phi values that change

!   Timing variables
    INTEGER NADAPT(NmoveTypes) ! Nunber of steps between adapt
    integer NPT                ! number of steps between parallel tempering
    integer INDMAX             ! total number of save points
    integer IND                ! save point
    integer NSTEP              ! steps per save point
    integer NNoInt             ! save points before turning on NNoInt
    integer N_KAP_ON           ! when to turn KAP energy on
    integer N_CHI_ON           ! when to turn CHI energy on
    
!   Switches
    INTEGER confineType       ! type of Boundary Conditions
    INTEGER setType           ! initial condition type
!    INTEGER INTON             ! interaction on
    logical FRMCHEM           ! Initial chemical sequence
    logical FRMMETH           ! Read methalation from file
    logical FRMFILE           ! Read Initial condition R
    logical FRMField          ! Read field from file
    logical saveU             ! save U vectors to file
    logical savePhi           ! save Phi vectors to file
    integer simType           ! Melt vs. Solution, Choose hamiltonian
    logical recenter_on       ! recenter in periodic boundary
    integer winType           ! how to choose random section of polymer to move
    integer reduce_move       ! only exicute unlikely movetypes every ____ cycles
    double precision min_accept
    logical UseSchedule       ! use scheduled change in interactions strength(s)
    double precision KAP_ON   
    double precision CHI_ON
    double precision Couple_ON
    logical restart

!   Parallel Tempering variables
    Character*16 repSufix    ! prefix for writing files
    integer rep  ! which replica am I
    integer id   ! which thread am I
    integer (kind = 4) error  ! MPI error
    double precision M ! M=\sum_i \sigma_i   like ising magnitization
    logical PTON
    ! see Timing variables for NPT\
    logical PT_chi
    logical PT_h
    logical PT_kap
    logical PT_mu
    logical PT_couple

!   Replica Dynamic Cof choice 
    integer NRepAdapt ! number of exchange attemts between adapt
    double precision lowerRepExe ! when to decrese cof spacing
    double precision upperRepExe ! when to increase cof spacing
    double precision lowerCofRail ! minumum acceptable Cof
    double precision upperCofRail ! maximum acceptable Cof
    integer indStartRepAdapt
    integer indEndRepAdapt
    double precision repAnnealSpeed  ! for annealing
    logical replicaBounds
    double precision INITIAL_MAX_S

  end Type

  Type MCData  ! for Allocateable variables
!   Configuration Data
    REAL(dp), ALLOCATABLE, DIMENSION(:,:):: R   ! Conformation of polymer chains
    REAL(dp), ALLOCATABLE, DIMENSION(:,:):: U   ! Conformation of polymer chains 
    REAL(dp), ALLOCATABLE, DIMENSION(:,:):: RP !Test Bead positions - only valid from IT1 to IT2
    REAL(dp), ALLOCATABLE, DIMENSION(:,:):: UP !Test target vectors - only valid from IT1 to IT2
    REAL(dp), ALLOCATABLE, DIMENSION(:):: PHIA ! Volume fraction of A
    REAL(dp), ALLOCATABLE, DIMENSION(:):: PHIB ! Volume fraction of B
    REAL(dp), ALLOCATABLE, DIMENSION(:):: PHIH ! Applied Field
    REAL(dp), ALLOCATABLE, DIMENSION(:):: Vol ! Volume fraction of A
    INTEGER, ALLOCATABLE, DIMENSION(:):: AB            ! Chemical identity of beads
    INTEGER, ALLOCATABLE, DIMENSION(:):: ABP           ! Test Chemical identity of beads
    INTEGER, ALLOCATABLE, DIMENSION(:):: METH          ! Methalation state of beads
    REAl(dp), Allocatable, Dimension(:):: DPHIA    ! Change in phi A
    REAl(dp), Allocatable, Dimension(:):: DPHIB    ! Change in phi A
    INTEGER, Allocatable, Dimension(:) :: INDPHI   ! Indices of the phi
  end TYPE
contains
Subroutine MCvar_setParams(mc,fileName)
! Based on Elena's readkeys subroutine
    USE INPUTPARAMS, ONLY : READLINE, READA, READF, READI, READO
    IMPLICIT NONE
    DOUBLE PRECISION, Parameter :: PI=3.141592654_dp
    TYPE(MCvar), intent(out) :: mc
    character*16, intent(in) :: fileName  ! file with parameters
    INTEGER :: PF   ! input file unit
    LOGICAL :: FILEEND=.FALSE. ! done reading file?
    CHARACTER*100 :: WORD ! keyword
    INTEGER :: NITEMS ! number of items on the line in the parameter file
    DOUBLE PRECISION NAN_dp
    NAN_dp=0;NAN_dp=NAN_dp/NAN_dp

    ! ----------------------------------------------------------
    !
    !  set Default values
    !
    ! ----------------------------------------------------------

    ! file IO
    mc%FRMFILE=.FALSE.
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
    mc%LBOX(1)=25.0_dp
    mc%LBOX(2)=25.0_dp
    mc%LBOX(3)=25.0_dp
    mc%DEL =1.0_dp
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
    mc%moveTypes=10
    mc%setType = 4 ! 4 for shereical
    mc%confineType = 3 ! 3 for sherical
    mc%simType=1
    mc%winType=1
    mc%min_accept=0.05

    ! timing options    
    mc%NStep=400000
    mc%NNoInt=100
    mc%INDMAX=180   
    mc%reduce_move=10
    mc%UseSchedule=.False.
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
    mc%indEndRepAdapt=20
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
    !  Read from file
    !
    !------------------------- 
    PF=55
    OPEN(UNIT=PF,FILE=fileName,STATUS='OLD') 

    ! read in the keywords one line at a time
    DO 
       CALL READLINE(PF,FILEEND,NITEMS)
       IF (FILEEND.and.nitems.eq.0) EXIT

       ! skip empty lines
       IF (NITEMS.EQ.0) CYCLE

       ! Read in the keyword for this line
       CALL READA(WORD,CASESET=1)

       ! Skip any empty lines or any comment lines
       IF (WORD(1:1).EQ.'#') CYCLE

       SELECT CASE(WORD) ! pick which keyword
       CASE('SETTYPE')
           Call READI(mc%setType)
           ! setType      |  Discription
           ! _____________|_________________________________
           !    1         |   straight line in y direction with random starting
           !    2         |   rerandomize when reaching boundary, slit in z dir
           !    3         |   rerandomize when reaching boundary, cube boundary
           !    4         |   rerandomize when reaching boundary, shpere
       CASE('CONFINETYPE')
           Call READI(mc%confineType)
           ! confineType  |  Discription
           ! _____________|_________________________________
           !    0         |  No confinement, periodic cube
           !    1         |  Between two plates in Z direction at 0 and LBox
           !    2         |  Cube of size LBox**3,  range: 0-LBox
           !    3         |  Circle of radius LBox, centered at LBox/2
           !    4         |  Periodic, unequal dimensions
       CASE('RECENTER_ON')
           Call READO(mc%recenter_on) ! recenter in periodic boundary
       CASE('SIMTYPE')
           Call READI(mc%simType) 
           ! simType      | Discription
           !______________|_________________________________
           !    0         | Melt density fluctuates around fixed mean
           !    1         | Solution (For DNA)
       CASE('FRMCHEM')
           Call READO(mc%FRMCHEM) ! Initial chemical sequence from file
       CASE('FRMFILE')
           call READO(mc%FRMFILE) ! Read configuration from file
       CASE('FRMMETH')
           Call READO(mc%FRMMETH) ! Read methalation from file
       CASE('PTON')
           CALL READO(mc%PTON) ! Parallel Tempering on
       CASE('SAVE_U')
           Call READO(mc%saveU)  ! save u vectors to file (every savepoint)
       CASE('SAVE_PHI')
           Call READO(mc%savePhi) ! save Phi vectors to file (every savepoint)
       CASE('L0')
           Call READF(mc%L0)  ! Equilibrium segment length
       CASE('DEL')
           Call READF(mc%DEL) ! spaitial descretation length, not tested
       CASE('LBOX')
           Call READF(mc%LBox(1)) ! side length of box
           mc%LBox(2)=mc%LBox(1)
           mc%LBox(3)=mc%LBox(1)
       CASE('LBOXX')
           Call READF(mc%LBox(1)) ! side length of box in x direction
       CASE('LBOXY')
           Call READF(mc%LBox(2)) ! side length of box in y direction
       CASE('LBOXZ')
           Call READF(mc%LBox(3)) ! side length of box in z direction
       CASE('NP')
           CALL READI(mc%NP)  ! Number of polymers
       CASE('G')
           Call READI(mc%G) ! Beads per monomer
       CASE('N')
           CALL READI(mc%N) ! Number of monomers in a polymer
       CASE('NNOINT')
           Call READI(mc%NNoInt) ! save points before turning on interaction
       CASE('N_KAP_ON')
           call readI(mc%N_KAP_ON) ! when to turn compression energy on
       CASE('N_CHI_ON')
           call readI(mc%N_CHI_ON) ! when to turn CHI energy on
       CASE('INDMAX')
           Call READI(mc%INDMAX) ! total number of save points
       CASE('NSTEP')
           Call READI(mc%NStep) ! steps per save point
       CASE('NPT')
           Call READI(mc%NPT) ! number of steps between parallel tempering
       CASE('FPOLY')
           Call READF(mc%Fpoly) ! Fraction Polymer
       CASE('V')
           Call READF(mc%V) ! Bead volume
       CASE('FA')
           Call READF(mc%FA) ! Fraction of A beads (fraction bound)
       CASE('LAM')
           Call READF(mc%LAM) ! Chemical correlation parameter
       CASE('EPS')
           Call READF(mc%EPS) ! Elasticity l0/(2lp) 
       CASE('CHI')
           Call READF(mc%CHI) ! CHI parameter (definition depends on  hamiltoniaon
       CASE('H_A')
           Call READF(mc%h_A) ! strength of externally applied field
       CASE('KAP')
           Call READF(mc%KAP) !  Incompressibility parameter 
       CASE('EU')
           Call READF(mc%EU) ! Energy of binding for unmethalated
       CASE('EM')
           Call READF(mc%EM) ! Energy of binding for methalated
       CASE('MU')
           Call READF(mc%MU) ! chemical potential of HP1
       CASE('HP1_BIND')
           Call READF(mc%HP1_BIND) ! Energy of binding of HP1 to eachother
       CASE('F_METH')
           Call READF(mc%F_METH) ! Fraction methalated
       CASE('LAM_METH')
           Call READF(mc%LAM_METH) ! eigenvalue of methalation setup
       CASE('CRANK_SHAFT_ON')
           Call READI(mc%MOVEON(1)) ! is Crank shaft move on 1/0
       CASE('SLIDE_ON')
           Call READI(mc%MOVEON(2)) ! is Slide move on 1/0
       CASE('PIVOT_ON')
           Call READI(mc%MOVEON(3)) ! is Pivot move on 1/0
       CASE('ROTATE_ON')
           Call READI(mc%MOVEON(4)) ! is single bead rotate on 1/0
       CASE('FULL_CHAIN_ROTATION_ON')
           Call READI(mc%MOVEON(5)) ! is full chain rotate on 1/0
       CASE('FULL_CHAIN_SLIDE_ON')
           Call READI(mc%MOVEON(6)) ! is full chain slide on 1/0
       CASE('BIND_MOVE_ON')
           Call READI(mc%MOVEON(7)) ! is bind/unbind move on 1/0
       CASE('CHAIN_FLIP_MOVE_ON')
           Call READI(mc%MOVEON(8)) ! is flip move move on 1/0
       CASE('CHAIN_SWAP_MOVE_ON')
           Call READI(mc%MOVEON(9)) ! is chain swap move on 1/0
       CASE('REPTATION_MOVE_ON')
           Call READI(mc%MOVEON(10)) ! is reptation move on 1/0
       CASE('WINTYPE')
           Call READI(mc%winType)  ! 0 for uniform, 1 for exponential
       CASE('MIN_CRANK_SHAFT_WIN') 
           Call READF(mc%MINWINDOW(1)) ! min mean window size
       CASE('MIN_SLIDE_WIN') 
           Call READF(mc%MINWINDOW(2)) 
       CASE('MIN_PIVOT_WIN') 
           Call READF(mc%MINWINDOW(3))
       CASE('MIN_BIND_WIN') 
           Call READF(mc%MINWINDOW(7))
       CASE('REDUCE_MOVE')
           Call READI(mc%reduce_move) !  only exicute unlikely movetypes every ____ cycles
       CASE('MIN_ACCEPT')
           Call READF(mc%MIN_ACCEPT) ! below which moves are turned off
       CASE('CRANK_SHAFT_TARGET')
           Call READF(mc%winTarget(1)) ! target window size for crank shaft move
       CASE('SLIDE_TARGET')
           Call READF(mc%winTarget(2)) ! target window size for slide move
       CASE('PIVOT_TARGET')
           Call READF(mc%winTarget(3)) ! target window size for Pivot move
       CASE('STRENGTH_SCHEDULE')
           Call READO(mc%UseSchedule) ! use scheduled ramp in interaction strength(s)
       CASE('N_REP_ADAPT')
           Call READI(mc%NRepAdapt)  ! number of exchange attemts between adapt 
       CASE('LOWER_REP_EXE')
           Call READF(mc%lowerRepExe) ! when to decrease cof spacing
       CASE('UPPER_REP_EXE')
           Call READF(mc%upperRepExe) ! when to increase cof spacing
       CASE('LOWER_COF_RAIL')
           Call READF(mc%lowerCofRail) ! minumum acceptable Cof
       CASE('UPPER_COF_RAIL')
           Call READF(mc%upperCofRail) ! maximum acceptable Cof
       CASE('IND_START_REP_ADAPT')
           Call READI(mc%indStartRepAdapt) ! ind to start rep. cof. adaptiation on
       CASE('IND_END_REP_ADAPT')
           Call READI(mc%indEndRepAdapt) ! turn off rep adapt
       CASE('REP_ANNEAL_SPEED')
           Call READF(mc%repAnnealSpeed)  ! max change in cof. every adjust
       CASE('FRMFIELD')
           Call READO(mc%FRMFIELD)  ! read field from file
       CASE('K_FIELD')
           Call READF(mc%k_field)  ! wave mode for default field
       CASE('REPLICA_BOUNDS')
           Call READO(mc%replicaBounds) ! insure that 0 < s < 1
       CASE('INITIAL_MAX_S')
           call READF(mc%INITIAL_MAX_S) ! inital s of rep with highest s
       CASE('PT_CHI')
           call READO(mc%PT_chi) ! parallel temper chi
       CASE('PT_H')
           call READO(mc%PT_h) ! parallel temper h
       CASE('PT_KAP')
           call READO(mc%PT_kap) ! parallel temper kap
       CASE('PT_MU')
           call READO(mc%PT_mu) ! parallel temper mu 
       CASE('PT_COUPLE')
           call READO(mc%PT_couple) ! parallel temper HP1_bind
       CASE('RESTART')
           call READO(mc%restart) ! Restart from parallel tempering
       CASE DEFAULT
           print*, "Error in MCvar_setParams.  Unidentified keyword:", &
                   TRIM(WORD)
           stop 1
       ENDSELECT
    ENDDO
    close(PF)

    if ((mc%NBINX(1)-mc%NBINX(2).ne.0).or. &
        (mc%NBINX(1)-mc%NBINX(3).ne.0)) then
        if (mc%simType.eq.1) then
            print*, "Solution not tested with non-cube box, more coding needed"
            stop 1
        endif
        if (mc%confineType.ne.4) then
            print*, "Unequal boundaries require confineType=4"
            stop 1
        endif    
        if (mc%setType.eq.4) then
            print*, "You shouldn't put a shpere in and unequal box!"
            stop 1
        endif    
    endif
    ! --------------------
    !
    ! Derived Variables, Reconcile inputs
    !
    ! --------------------
    if (mc%simType.eq.1) then
        mc%NT=mc%N*mc%NP*mc%G
        mc%NB=mc%N*mc%G
        if (mc%confineType.eq.3) then
            mc%LBOX(1)=(mc%V*mc%NT*6/(mc%Fpoly*PI))**(1.0_dp/3.0_dp)
            mc%LBOX(2)=mc%LBOX(1)
            mc%LBOX(3)=mc%LBOX(1)
        else
            mc%LBOX(1)=(mc%V*mc%NT/mc%Fpoly)**(1.0_dp/3.0_dp)
            mc%LBOX(2)=mc%LBOX(1)
            mc%LBOX(3)=mc%LBOX(1)
        endif
        mc%NBINX(1)=nint(mc%LBOX(1)/mc%DEL)
        mc%NBINX(2)=nint(mc%LBOX(2)/mc%DEL)
        mc%NBINX(3)=nint(mc%LBOX(3)/mc%DEL)
        mc%NBIN=mc%NBINX(1)*mc%NBINX(2)*mc%NBINX(3)
        mc%LBOX(1) = mc%NBINX(1)*mc%DEL! used to be: DEL=LBOX/NBINX
        mc%LBOX(2) = mc%NBINX(2)*mc%DEL! used to be: DEL=LBOX/NBINX
        mc%LBOX(3) = mc%NBINX(3)*mc%DEL! used to be: DEL=LBOX/NBINX
    elseif (mc%simType.eq.0) then
        if (mc%confineType.eq.0) then
            mc%NP=nint(mc%LBOX(1)*mc%LBOX(2)*mc%LBOX(3)/(mc%N*mc%G*mc%V))
            mc%LBOX=(mc%V*mc%N*mc%G*mc%NP)**(1.0_dp/3.0_dp)
            mc%NBINX(1)=nint(mc%LBOX(1)/mc%DEL)
            mc%NBINX(2)=nint(mc%LBOX(2)/mc%DEL)
            mc%NBINX(3)=nint(mc%LBOX(3)/mc%DEL)
            mc%DEL=mc%LBOX(1)/mc%NBINX(1)
        elseif(mc%confineType.eq.4) then
            mc%DEL=mc%LBOX(1)/nint(mc%LBOX(1)/mc%DEL)
            mc%NBINX(1)=nint(mc%LBOX(1)/mc%DEL)
            mc%NBINX(2)=nint(mc%LBOX(2)/mc%DEL)
            mc%NBINX(3)=nint(mc%LBOX(3)/mc%DEL)
            mc%LBOX(2)=mc%DEL*mc%NBINX(2)
            mc%LBOX(3)=mc%DEL*mc%NBINX(3)
            mc%NP=nint(mc%LBOX(1)*mc%LBOX(2)*mc%LBOX(3)/(mc%N*mc%G*mc%V))
            print*, "Density =", &
                  mc%N*mc%G*mc%V*mc%NP/(mc%LBOX(1)*mc%LBOX(2)*mc%LBOX(3))
        endif
        mc%NB=mc%N*mc%G
        mc%NBIN=mc%NBINX(1)*mc%NBINX(2)*mc%NBINX(3)
        mc%NT=mc%N*mc%NP*mc%G
        mc%WINDOW(5)=mc%NB
        mc%WINDOW(6)=mc%NB
        mc%WINDOW(8)=mc%NB
        mc%WINDOW(9)=mc%NB
        mc%MCAMP(2)=0.3_dp*mc%L0
        mc%MCAMP(6)=5.0_dp*mc%L0
        mc%MINAMP(2)=0.2_dp*mc%L0
        mc%MINAMP(6)=0.2_dp*mc%L0
        mc%MAXWINDOW(1)=min(150,mc%NB)
        mc%MAXWINDOW(2)=min(150,mc%NB)
        mc%MAXWINDOW(3)=min(150,mc%NB)
        mc%MAXWINDOW(7)=min(150,mc%NB)
        mc%MAXAMP(2)=1.0_dp*mc%L0
        mc%MAXAMP(6)=0.1*mc%LBOX(1)
    else
       print*, "Error in simMod: symType",mc%simType," not found"
    endif 
    call getpara(mc%PARA,mc%EPS,mc%L0,NAN_dp)
   
    ! -----------------------
    !
    ! Set initial values
    !
    ! ----------------------
    mc%EElas(1)=0.0_dp
    mc%EElas(2)=0.0_dp
    mc%EElas(3)=0.0_dp
    mc%ECouple=0.0_dp
    mc%EBind=0.0_dp
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
        print*, "error in MCsim. Don't run adapt without int on"
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
    
end Subroutine
Subroutine MCvar_printDescription(mc)
    IMPLICIT NONE
    TYPE(MCvar), intent(in) :: mc
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
    print*, " LBOX=", mc%LBOX(1), mc%LBOX(2), mc%LBOX(3)
    print*, " Number of bins in x direction", &
             mc%NBINX(1), mc%NBINX(2),mc%NBINX(3)
    print*, " Number of bins", mc%NBIN
    print*, " spatial descritation DEL=",mc%DEL
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
    print*, " confineType:",mc%confineType
    print*, " setType:",mc%setType
    print*, "---------------------------------------------"
    
end Subroutine
Subroutine MCvar_allocate(mc,md)
    IMPLICIT NONE
    TYPE(MCvar), intent(in) :: mc
    TYPE(MCData), intent(out) :: md
    INTEGER NT  ! total number of beads
    INTEGER NBIN ! total number of bins
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
    Allocate(md%INDPHI(NBIN))
    Allocate(md%PhiH(NBIN))

end subroutine
Subroutine MCvar_defaultAmp(mc)
    IMPLICIT NONE
    DOUBLE PRECISION, Parameter :: PI=3.141592654_dp
    TYPE(MCvar), intent(inout) :: mc
    INTEGER MCTYPE ! Type of move
    INTEGER NANI ! NaN
    DOUBLE PRECISION NAND !NAND
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
    mc%WINDOW(1)=15.0_dp ! used to be N*G
    mc%WINDOW(2)=15.0_dp ! used to be N*G
    mc%WINDOW(3)=15.0_dp ! used to be N*G
    mc%WINDOW(4)=1.0_dp
    mc%WINDOW(5)=dble(mc%N*mc%G)
    mc%WINDOW(6)=dble(mc%N*mc%G)
    mc%WINDOW(7)=15.0_dp ! used to be N*G
    mc%WINDOW(8)=dble(mc%N*mc%G)
    mc%WINDOW(9)=dble(mc%N*mc%G)
    mc%WINDOW(9)=1.0_dp

    !    Maximum window size (large windows are expensive)
    mc%MAXWINDOW(1)=dble(min(150,mc%NB))
    mc%MAXWINDOW(2)=dble(min(150,mc%NB))
    mc%MAXWINDOW(3)=dble(min(150,mc%NB))
    mc%MAXWINDOW(4)=NAND 
    mc%MAXWINDOW(5)=NAND 
    mc%MAXWINDOW(6)=NAND
    mc%MAXWINDOW(7)=dble(min(4,mc%NB))
    mc%MAXWINDOW(8)=NAND
    mc%MAXWINDOW(9)=NAND
    mc%MAXWINDOW(9)=NAND ! need to chaige code to allow >1

    mc%MINWINDOW(1)=dble(min(4,mc%NB))
    mc%MINWINDOW(2)=dble(min(4,mc%NB))
    mc%MINWINDOW(3)=dble(min(4,mc%NB))
    mc%MINWINDOW(4)=NAND 
    mc%MINWINDOW(5)=NAND 
    mc%MINWINDOW(6)=NAND
    mc%MINWINDOW(7)=dble(min(4,mc%NB))
    mc%MINWINDOW(8)=NAND
    mc%MINWINDOW(9)=NAND
    mc%MINWINDOW(10)=NAND
    Do MCTYPE=1,mc%moveTypes
        mc%winTarget(MCTYPE)=8.0_dp
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
    mc%MAXAMP(6)=0.1*mc%LBOX(1)
    mc%MAXAMP(7)=NAND
    mc%MAXAMP(8)=NAND
    mc%MAXAMP(9)=NAND
    mc%MAXAMP(9)=NAND
     
    DO MCTYPE=1,mc%moveTypes
        mc%NADAPT(MCTYPE)=1000 ! adapt after at most 1000 steps
        mc%PDESIRE(MCTYPE)=0.5_dp ! Target
        mc%SUCCESS(MCTYPE)=0
        mc%PHIT(MCTYPE)=0.0_dp
    ENDDO
end subroutine
Subroutine MCvar_recenter(mc,md)
!  Prevents drift in periodic BC
    IMPLICIT NONE
    TYPE(MCvar), intent(in) :: mc
    TYPE(MCData), intent(inout) :: md
    INTEGER IB, I, J   ! Couners
    DOUBLE PRECISION R0(3)  ! Offset to move by
    IB=1
    DO I=1,mc%NP
       R0(1)=nint(md%R(IB,1)/mc%LBOX(1)-0.5_dp)*mc%LBOX(1)
       R0(2)=nint(md%R(IB,2)/mc%LBOX(2)-0.5_dp)*mc%LBOX(2)
       R0(3)=nint(md%R(IB,3)/mc%LBOX(3)-0.5_dp)*mc%LBOX(3)
       if (abs(R0(1)*R0(2)*R0(3)) .gt. 0.0001_dp) then
           DO J=1,mc%NB
              md%R(IB,1)=md%R(IB,1)-R0(1)
              md%R(IB,2)=md%R(IB,2)-R0(2)
              md%R(IB,3)=md%R(IB,3)-R0(3)
              IB=IB+1
           Enddo
      endif
    enddo
end Subroutine
Subroutine MCvar_printEnergies(mc)
! For realtime feedback on MC simulation
    IMPLICIT NONE
    TYPE(MCvar), intent(in) :: mc
    print*, "ECouple:", mc%ECouple
    print*, "Bending energy", mc%EELAS(1)
    print*, "Par compression energy", mc%EELAS(2)
    print*, "Shear energy", mc%EELAS(3)
    print*, "ECHI", mc%ECHI
    print*, "EField", mc%EField
    print*, "EKAP", mc%EKAP
    print*, "EBind", mc%EBind
end subroutine
Subroutine MCvar_printPhi(mc,md)
! Prints densities for trouble shooting
    IMPLICIT NONE
    TYPE(MCvar), intent(in) :: mc
    TYPE(MCData), intent(in) :: md
    Integer I
    DOUBLE PRECISION EKap, ECouple, EChi,VV, PHIPOly
    print*,"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    print*, " PHIA  | PHIB  | PPoly |  Vol  | EKap  | EChi  |ECouple|" 
    Do I=1,mc%NBIN
        VV=md%Vol(I)
        if (VV.le.0.1_dp) CYCLE
        PHIPOLY=md%PHIA(I)+md%PHIB(I)
        EChi=VV*(mc%CHI/mc%V)*PHIPoly*(1.0_dp-PHIPoly)
        ECouple=VV*mc%HP1_Bind*(md%PHIA(I))**2
        if(PHIPoly.GT.1.0_dp) then
           EKap=VV*(mc%KAP/mc%V)*(PHIPoly-1.0_dp)**2
        else
           CYCLE
           EKap=0.0_dp
        endif
        write(*,"(4f8.4,3f8.1)"), md%PHIA(I), md%PHIB(I), & 
                            md%PHIA(I)+md%PHIB(I),md%Vol(I),&
                            EKap,EChi,ECouple
    enddo
    print*,"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
end subroutine
Subroutine MCvar_printWindowStats(mc)
! For realtime feedback on adaptation
    IMPLICIT NONE
    TYPE(MCvar), intent(in) :: mc
    INTEGER I ! counter
    I=0
    print*, "Succes | MCAMP | WINDOW| Type "
    Do I=1,mc%moveTypes
        if (mc%MOVEON(i).eq.1) then
            write(*,"(f8.5,2f8.2,1I8)"), mc%phit(i), mc%MCAMP(i),  mc%WINDOW(i), i
        endif
    enddo
    return
end subroutine
Subroutine MCvar_LoadField(mc,md,fileName)
    IMPLICIT NONE
    TYPE(MCvar), intent(in) :: mc
    TYPE(MCData), intent(inout) :: md
    Integer I
    character*16 fileName ! file name to load from
    OPEN (UNIT = 1, FILE = fileName, STATUS = 'OLD')      
    Do I=1,mc%NBIN
        READ(1,*) md%PHIH(I)
    enddo
    return
end subroutine
subroutine MCvar_MakeField(mc,md)
    IMPLICIT NONE
    TYPE(MCvar), intent(in) :: mc
    TYPE(MCData), intent(inout) :: md
    integer INDBIN  ! index of bin
    integer IX,IY,IZ ! bin corrdinates

    do IX=1,mc%NBINX(1)
        do IY=1,mc%NBINX(2)
            do IZ=1,mc%NBINX(3)
                INDBIN=IX+&
                       (IY-1)*mc%NBINX(1)+&
                       (IZ-1)*mc%NBINX(1)*mc%NBINX(2)
                md%PHIH(INDBIN)=dsin(mc%k_field*mc%DEL*dble(IX))
            enddo
        enddo
    enddo
    return
end subroutine
Subroutine MCvar_loadAB(mc,md,fileName)
! Loads AB for file...has not been tested
    IMPLICIT NONE
    TYPE(MCvar), intent(in) :: mc
    TYPE(MCData), intent(inout) :: md
    character*16, intent(in) :: fileName ! file name to load from
    INTEGER IB, I, J ! counters
    OPEN (UNIT = 1, FILE = fileName, STATUS = 'OLD')      
    IB=1
    DO I=1,mc%NP
       DO J=1,mc%NB
          READ(1,"(I2)") md%AB(IB)
          IB=IB+1
          enddo
    enddo 
    CLOSE(1)
end subroutine
Subroutine MCvar_saveR(mc,md,fileName,repeatingBC)
! Writes R and AB to file for analysis
! Rx  Ry  Rz AB
    IMPLICIT NONE
    INTEGER, intent(in) :: repeatingBC  ! 1 for reapeating boundary conditions
    INTEGER I,J,IB  ! counters
    TYPE(MCvar), intent(in) :: mc
    TYPE(MCData), intent(in) :: md
    character*16, intent(in) :: fileName
    character*32 fullName
    fullName=  trim(fileName) // trim(mc%repSufix) 
    fullName=trim(fullName)
    OPEN (UNIT = 1, FILE = fullName, STATUS = 'NEW')
    IB=1
    if (repeatingBC.eq.1) then 
        Do I=1,mc%NP
            Do J=1,mc%NB
                    WRITE(1,"(3f10.3,I2)") , &
                          md%R(IB,1)-0.*nint(md%R(IB,1)/mc%LBOX(1)-0.5_dp)*mc%LBOX(1), &
                          md%R(IB,2)-0.*nint(md%R(IB,2)/mc%LBOX(2)-0.5_dp)*mc%LBOX(2), &
                          md%R(IB,3)-0.*nint(md%R(IB,3)/mc%LBOX(3)-0.5_dp)*mc%LBOX(3), & 
                          md%AB(IB)
                IB=IB+1
            enddo
        enddo
        print*, "Error in MCvar_saveR"
        print*, "Are you sure you want reapeating BC"
        stop 1
    else
        Do I=1,mc%NP
            Do J=1,mc%NB
                if (mc%simtype.eq.0) then
                    WRITE(1,"(3f10.3,I2)") md%R(IB,1),md%R(IB,2),md%R(IB,3),md%AB(IB)
                else
                    WRITE(1,"(3f10.3,I2)") md%R(IB,1),md%R(IB,2),md%R(IB,3),md%AB(IB), md%METH(IB)
                endif
                IB=IB+1
            enddo
        enddo
    endif
    Close(1)
end subroutine
Subroutine MCVar_savePHI(mc,md,fileName)
! Saves PHIA and PHIB to file for analysis
    IMPLICIT NONE
    INTEGER I  ! counters
    TYPE(MCvar), intent(in) :: mc
    TYPE(MCData), intent(in) :: md
    character*16, intent(in) :: fileName 
    character*32 fullName
    fullName=  trim(fileName) // trim(mc%repSufix) 
    OPEN (UNIT = 1, FILE = fullName, STATUS = 'NEW')
    Do I=1,mc%NBIN
        WRITE(1,"(2f7.2)") md%PHIA(I),md%PHIB(I)
    enddo
    Close(1)
end subroutine
Subroutine MCvar_saveU(mc,md,fileName)
! Saves U to ASCII file for analisys
    IMPLICIT NONE
    INTEGER I,J,IB  ! counters
    TYPE(MCvar), intent(in) :: mc
    TYPE(MCData), intent(in) :: md
    character*16, intent(in) :: fileName
    character*32 fullName
    fullName=  trim(fileName) // trim(mc%repSufix) 
    OPEN (UNIT = 1, FILE = fullName, STATUS = 'NEW')
    IB=1
    Do I=1,mc%NP
        Do J=1,mc%NB
            WRITE(1,"(3f8.3,2I2)") md%U(IB,1),md%U(IB,2),md%U(IB,3)
            IB=IB+1
        enddo
    enddo
    Close(1)
end subroutine
Subroutine MCvar_saveParameters(mc,fileName)
! Write a number of parameters ASCII variables to file for reccords
    IMPLICIT NONE
    TYPE(MCvar), intent(in) :: mc
    character*16, intent(in) :: fileName
    character*32 fullName
    fullName=  trim(fileName) // trim(mc%repSufix) 
    OPEN (UNIT =1, FILE = fullName, STATUS = 'NEW')
        WRITE(1,"(I8)") mc%NT ! 1 Number of beads in simulation
        WRITE(1,"(I8)") mc%N  ! 2 Number of monomers in a polymer
        WRITE(1,"(I8)") mc%NB ! 3 Number of beads in a polymer
        WRITE(1,"(I8)") mc%NP ! 4 Number of polymers in simulation
        WRITE(1,"(I8)") mc%NT ! 5 Number of beads in simulation
        WRITE(1,"(I8)") mc%G  ! 6 Number of beads per monomer

        WRITE(1,"(f10.5)") mc%L0    ! Equilibrium segment length 
        WRITE(1,"(f10.5)") mc%CHI  ! 8  initail CHI parameter value 
        WRITE(1,"(f10.5)") mc%Fpoly ! Fraction polymer
        WRITE(1,"(f10.5)") mc%LBOX(1)  ! 10 Lenth of box
        WRITE(1,"(f10.5)") mc%EU    ! Energy unmethalated       
        WRITE(1,"(f10.5)") mc%EM    ! 12 Energy methalated
        WRITE(1,"(f10.5)") mc%HP1_Bind ! Energy of HP1 binding
        WRITE(1,"(f10.5)") (mc%L0/mc%EPS) ! 14 Khun lenth
        WRITE(1,"(A)") "-999"  ! for historic reasons
        WRITE(1,"(f10.5)") mc%F_METH  ! methalation fraction
        WRITE(1,"(f10.5)") mc%LAM_METH  ! methalation lambda
    CLOSE(1)
end subroutine
Subroutine MCvar_appendEnergyData(mc,fileName)
! Print Energy data
    IMPLICIT NONE
    TYPE(MCvar), intent(in) :: mc
    LOGICAL isfile
    character*16, intent(in) :: fileName
    character*32 fullName
    fullName=  trim(fileName) // trim(mc%repSufix) 
    inquire(file = fullName, exist=isfile)
    if (isfile) then
        OPEN (UNIT = 1, FILE = fullName, STATUS ='OLD', POSITION="append")
    else 
        OPEN (UNIT = 1, FILE = fullName, STATUS = 'new')
        WRITE(1,*), "IND | id |",&
                    " EBend  | EParll | EShear | ECoupl | E Kap  | E Chi  |",&
                    " EField | EBind  |   M    | Couple |  Chi   |  mu    |",&
                    "  Kap   | Field  |"
    endif
    WRITE(1,"(2I5, 9f9.1,5f9.4)") mc%IND, mc%id, &
           mc%EELAS(1), mc%EELAS(2), mc%EELAS(3), mc%ECouple, &
           mc%EKap, mc%ECHI, mc%EField, mc%EBind, mc%M, &
           mc%HP1_Bind*mc%Couple_on, mc%CHI*mc%CHI_ON, mc%mu, mc%KAP*mc%KAP_ON,&
           mc%h_A
    Close(1)
end subroutine
Subroutine MCvar_appendAdaptData(mc,fileName)
! Appends MC move adaptation data to the file  
    IMPLICIT NONE
    TYPE(MCvar), intent(in) :: mc
    LOGICAL isfile
    character*16, intent(in) :: fileName
    character*32 fullName
    fullName=  trim(fileName) // trim(mc%repSufix) 
    inquire(file = fullName, exist=isfile)
    if (isfile) then
        OPEN (UNIT = 1, FILE = fullName, STATUS ='OLD', POSITION="append")
    else 
        OPEN (UNIT = 1, FILE = fullName, STATUS = 'new')
        WRITE(1,*), "IND| id|",&
                    " WIN 1 | AMP 1 | SUC 1 | WIN 2 | AMP 2 | SUC 2 |",&
                    " WIN 3 | AMP 3 | SUC 3 | ON  4 | AMP 4 | SUC 4 |",&
                    " ON  5 | AMP 5 | SUC 5 | ON  6 | AMP 6 | SUC 6 |",&
                    " ON  7 | SUC 7 | ON  8 | SUC 8 |", &
                    " ON  9 | SUC 9 | ON 10 | SUC 10|"
    endif
    WRITE(1,"(2I4,26f8.3)") mc%IND,mc%id,& 
          REAL(mc%WINDOW(1)),mc%MCAMP(1),mc%PHIT(1), &
          REAL(mc%WINDOW(2)),mc%MCAMP(2),mc%PHIT(2), &
          REAL(mc%WINDOW(3)),mc%MCAMP(3),mc%PHIT(3), &
          REAL(mc%MOVEON(4)),mc%MCAMP(4),mc%PHIT(4), &
          REAL(mc%MOVEON(5)),mc%MCAMP(5),mc%PHIT(5), &
          REAL(mc%MOVEON(6)),mc%MCAMP(6),mc%PHIT(6), &
          REAL(mc%MOVEON(7)),mc%PHIT(7), &
          REAL(mc%MOVEON(8)),mc%PHIT(8), &
          REAL(mc%MOVEON(9)),mc%PHIT(9), &
          REAL(mc%MOVEON(10)),mc%PHIT(10)
    Close(1)
end subroutine
Subroutine MCvar_writeBinary(mc,md,baceName)
!    This function writes the contence of the structures mc and md
!  to a binary file.  If you add more variables to md you need to 
!  a seperate write command for them as it is not possible to write
!  a structure with allocatables to a binar file.
!    The contence are stored in 
!     baceName//'R'
!     baceName//'U'
!     etc.
    IMPLICIT NONE
    INTEGER sizeOfType         ! for binary saving
    TYPE(MCvar), intent(in) :: mc             ! to be save or filled
    TYPE(MCData), intent(in) :: md             ! to be save or filled
    CHARACTER(LEN=16), intent(in) :: baceName ! for example 'record/'
    CHARACTER(LEN=16) fileName ! fileName
    CHARACTER(LEN=16) sufix    ! end of file name
    LOGICAL exists    ! Does file already exist?

    !  ------parameters -----

    sizeOfType=int(SIZEOF(mc))
    sufix='parameters'
    fileName=trim(baceName) // trim(sufix)
    inquire(file=fileName,exist=exists)
    if(exists) then
        open(unit=1,file=fileName, status='old', &
             form='unformatted',access='direct',recl=sizeOfType)
    else
        open(unit=1,file=fileName, status='new', &
             form='unformatted',access='direct',recl=sizeOfType)
    endif
    write(1,rec=1) mc
    close(1)    

    ! -------- R --------

    sizeOfType=int(SIZEOF(md%R))
    sufix='R'
    fileName=trim(baceName) // trim(sufix)
    inquire(file=fileName,exist=exists)
    if(exists) then
        open(unit=1,file=fileName, status='old', &
             form='unformatted',access='direct',recl=sizeOfType)
    else
        open(unit=1,file=fileName, status='new', &
             form='unformatted',access='direct',recl=sizeOfType)
    endif
    write(1,rec=1) md%R
    close(1)

    ! -------- U --------

    sizeOfType=int(SIZEOF(md%U))
    sufix='U'
    fileName=trim(baceName) // trim(sufix)
    inquire(file=fileName,exist=exists)
    if(exists) then
        open(unit=1,file=fileName, status='old', &
             form='unformatted',access='direct',recl=sizeOfType)
    else
        open(unit=1,file=fileName, status='new', &
             form='unformatted',access='direct',recl=sizeOfType)
    endif
    write(1,rec=1) md%U
    close(1)

    ! -------- AB --------

    sizeOfType=int(SIZEOF(md%AB))
    sufix='AB'
    fileName=trim(baceName) // trim(sufix)
    inquire(file=fileName,exist=exists)
    if(exists) then
        open(unit=1,file=fileName, status='old', &
             form='unformatted',access='direct',recl=sizeOfType)
    else
        open(unit=1,file=fileName, status='new', &
             form='unformatted',access='direct',recl=sizeOfType)
    endif
    write(1,rec=1) md%AB
    close(1)

    ! -------- Vol --------

    sizeOfType=int(SIZEOF(md%Vol))
    sufix='Vol'
    fileName=trim(baceName) // trim(sufix)
    inquire(file=fileName,exist=exists)
    if(exists) then
        open(unit=1,file=fileName, status='old', &
             form='unformatted',access='direct',recl=sizeOfType)
    else
        open(unit=1,file=fileName, status='new', &
             form='unformatted',access='direct',recl=sizeOfType)
    endif
    write(1,rec=1) md%Vol
    close(1)
end Subroutine

Subroutine MCvar_readBindary(mc,md,baceName)
! This function reads what MCvar_writeBinary writes and 
! stores it to mc and md.  Be sure to allocate md before 
! calling this command.
    IMPLICIT NONE
    INTEGER sizeOfType         ! for binary saving
    TYPE(MCvar) mc             ! to be save or filled
    TYPE(MCData) md             ! to be save or filled
    CHARACTER(LEN=16) baceName ! for example 'record/'
    CHARACTER(LEN=16) fileName ! fileName
    CHARACTER(LEN=16) sufix    ! end of file name
    LOGICAL exists    ! Does file already exist?

    !  ------parameters -----

    sizeOfType=int(SIZEOF(mc))
    sufix='parameters'
    fileName=trim(baceName) // trim(sufix)
    inquire(file=fileName,exist=exists)
    if(exists) then
        open(unit=1,file=fileName, status='old', &
             form='unformatted',access='direct',recl=sizeOfType)
    else
        print*, 'Error in MCvar_readBinary. File ',fileName,'does not exist'
        stop 1
    endif
    read(1,rec=1) mc
    close(1)    

    ! -------- R --------

    sizeOfType=int(SIZEOF(md%R))
    sufix='R'
    fileName=trim(baceName) // trim(sufix)
    inquire(file=fileName,exist=exists)
    if(exists) then
        open(unit=1,file=fileName, status='old', &
             form='unformatted',access='direct',recl=sizeOfType)
    else
        print*, 'Error in MCvar_readBinary. File ',fileName,'does not exist'
        stop 1
    endif
    read(1,rec=1) md%R
    close(1)

    ! -------- U --------

    sizeOfType=int(SIZEOF(md%U))
    sufix='U'
    fileName=trim(baceName) // trim(sufix)
    inquire(file=fileName,exist=exists)
    if(exists) then
        open(unit=1,file=fileName, status='old', &
             form='unformatted',access='direct',recl=sizeOfType)
    else
        print*, 'Error in MCvar_readBinary. File ',fileName,'does not exist'
        stop 1
    endif
    read(1,rec=1) md%U
    close(1)

    ! -------- AB --------

    sizeOfType=int(SIZEOF(md%AB))
    sufix='AB'
    fileName=trim(baceName) // trim(sufix)
    inquire(file=fileName,exist=exists)
    if(exists) then
        open(unit=1,file=fileName, status='old', &
             form='unformatted',access='direct',recl=sizeOfType)
    else
        print*, 'Error in MCvar_readBinary. File ',fileName,'does not exist'
        stop 1
    endif
    read(1,rec=1) md%AB
    close(1)

    ! -------- Vol --------

    sizeOfType=int(SIZEOF(md%Vol))
    sufix='Vol'
    fileName=trim(baceName) // trim(sufix)
    inquire(file=fileName,exist=exists)
    if(exists) then
        open(unit=1,file=fileName, status='old', &
             form='unformatted',access='direct',recl=sizeOfType)
    else
        print*, 'Error in MCvar_readBinary. File ',fileName,'does not exist'
        stop 1
    endif
    read(1,rec=1) md%Vol
    close(1)
end Subroutine
end module simMod
