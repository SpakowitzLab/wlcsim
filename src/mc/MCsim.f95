!---------------------------------------------------------------*
!
!
!
!     This subroutine performs a Monte Carlo simulation on the
!     polymer chain.
!
!    Quinn Made Changes to this file starting on 12/15/15
!

SUBROUTINE MCsim(mc,md,NSTEP)

    !use mt19937, only : grnd, sgrnd, rnorm, mt, mti
    use mersenne_twister
    use params

    IMPLICIT NONE

    INTEGER, intent(in) :: NSTEP             ! Number of MC steps
    LOGICAL :: INTON             ! Include polymer interactions

!   Variables for the simulation

    INTEGER ISTEP             ! Current MC step index
    DOUBLE PRECISION PROB     ! Calculated test prob
    DOUBLE PRECISION TEST     ! Random test variable
    INTEGER IP                ! Test polymer
    INTEGER IB1               ! Test bead position 1
    INTEGER IT1               ! Index of test bead 1
    INTEGER IB2               ! Test bead position 2
    INTEGER IT2               ! Index of test bead 2
    INTEGER IT3, IT4          ! second polymer for polymer swap
    logical forward           ! direction of reptation move

    INTEGER I,J

    INTEGER MCTYPE                    ! Type of MC move

    DOUBLE PRECISION EB,EPAR,EPERP
    DOUBLE PRECISION GAM,ETA
    DOUBLE PRECISION XIR,XIU
    DOUBLE PRECISION LHC      ! Length of HC int
    DOUBLE PRECISION VHC      ! HC strength
    DOUBLE PRECISION phiTot  ! for testing

    DOUBLE PRECISION ENERGY
    logical isfile
! Things for random number generator
    real urnd(1) ! single random number
    type(random_stat) :: rand_stat
!   Load the input parameters
    Type(wlcsim_params), intent(inout) :: mc      ! system varibles
    Type(wlcsim_data), intent(inout) :: md     ! system allocated data
!     Alexander Polynomial Variables
    DOUBLE PRECISION, ALLOCATABLE :: CROSS(:,:)   !Matrix of information for crossings in a 2-D projection of the polymer
    DOUBLE PRECISION, ALLOCATABLE :: CROSSP(:,:)  !Matrix of crossings for the trial configuration
    INTEGER NCROSS
    INTEGER NCROSSP
    INTEGER CrossSize
    INTEGER DELTA             !Alexander polynomial evaluated at t=-1; used for knot checking
    INTEGER DELTAP            !Alexandper polynomial of trial configuration

    inton = mc%inton
    rand_stat = md%rand_stat


    EB=   mc%PARA(1)
    EPAR= mc%PARA(2)
    EPERP=mc%PARA(3)
    GAM=  mc%PARA(4)
    ETA=  mc%PARA(5)
    XIR=  mc%PARA(6)
    XIU=  mc%PARA(7)
    LHC=  mc%PARA(9)
    VHC=  mc%PARA(10)
! -------------------------------------
!
!   initialize densities and energies
!
! -------------------------------------
    ! --- Binding Energy ---
    md%ABP=0 ! set entire array to zero
    !  Notide that ABP and AB are intensionally swapped below
    IT1=1; IT2=mc%NT
    call MC_bind(mc%NT,mc%G,IT1,IT2,md%ABP,md%AB,md%METH, &
                 mc%EU,mc%EM,mc%DEBind,mc%mu,mc%dx_mu)

    inquire(file = "data/error", exist=isfile)
    if (isfile) then
        OPEN (UNIT = 3, FILE = "data/error", STATUS ='OLD', POSITION="append")
    else
        OPEN (UNIT = 3, FILE = "data/error", STATUS = 'new')
    endif

    if(abs(mc%EBind-mc%DEBind).gt.0.00001) then
        print*, "Warning. Integrated binding enrgy:", &
                mc%EBind," while absolute binding energy:", &
                mc%DEBind
        write(3,*), "Warning. Integrated binding enrgy:", &
                mc%EBind," while absolute binding energy:", &
                mc%DEBind
    endif
    mc%EBind=mc%DEBind
    mc%x_mu=mc%dx_mu


    ! --- Elastic Energy ---
    call energy_elas(mc%DEELAS,md%R,md%U,mc%NT,mc%NB,mc%NP,mc%Para)
    if(abs((mc%EElas(1)+  mc%EElas(2)+ mc%EElas(3))-&
           (mc%DEElas(1)+mc%DEElas(2)+mc%DEElas(3))).gt.0.0001) then
        print*, "Warning. Integrated elastic enrgy:", &
                (mc%EElas(1)+mc%EElas(2)+mc%EElas(3)),&
                " while absolute elastic energy:", &
                (mc%DEElas(1)+mc%DEElas(2)+mc%DEElas(3))
        write(3,*), "Warning. Integrated elastic enrgy:", &
                (mc%EElas(1)+mc%EElas(2)+mc%EElas(3)),&
                " while absolute elastic energy:", &
                (mc%DEElas(1)+mc%DEElas(2)+mc%DEElas(3))
    endif
    mc%EElas=mc%DEElas ! copy array

    ! --- Interaction Energy ---
    if (INTON) then
        ! initialize phi
        IT1=1
        IT2=mc%NT ! need to set up all beads
        do I=1,mc%NBIN
             md%PHIA(I)=0.0_dp
             md%PHIB(I)=0.0_dp
        enddo
        call MC_int(mc,md,IT1,IT2,.True.)
        do I=1,mc%NBIN
            phiTot=phiTot+(md%PHIA(I)+md%PHIB(I))*md%Vol(I)
        enddo
        ! test to see if sum of changes are same as calculating from scratch
        print*, "phiTot", phiTot," NT:",mc%NT
        if(abs(mc%EChi-mc%DEChi).gt. 0.0001_dp) then
             print*, "Warning. Intigrated chi energy:", &
                     mc%EChi,"  while absolute chi energy:", &
                     mc%DEChi
             write(3,*), "Warning. Intigrated chi energy:", &
                     mc%EChi,"  while absolute chi energy:", &
                     mc%DEChi
        endif
        mc%EChi=mc%DEChi
        mc%x_chi=mc%dx_chi
        if(abs(mc%ECouple-mc%DECouple).gt. 0.0001_dp) then
             print*, "Warning. Intigrated couple energy:", &
                     mc%ECouple,"  while absolute couple energy:", &
                     mc%DECouple
             write(3,*), "Warning. Intigrated couple energy:", &
                     mc%ECouple,"  while absolute couple energy:", &
                     mc%DECouple
        endif
        mc%ECouple=mc%DECouple
        mc%x_Couple=mc%dx_couple
        if(abs(mc%EKap-mc%DEKap).gt. 0.0001_dp) then
             print*, "Warning. Intigrated Kap energy:", &
                     mc%EKap,"  while absolute Kap energy:", &
                     mc%DEKap
             write(3,*), "Warning. Intigrated Kap energy:", &
                     mc%EKap,"  while absolute Kap energy:", &
                     mc%DEKap
        endif
        mc%EKap=mc%DEKap
        mc%x_Kap=mc%dx_Kap

        if(abs(mc%EField-mc%DEField).gt.0.00001) then
            print*, "Warning. Integrated field enrgy:", &
                    mc%EField," while absolute field energy:", &
                    mc%DEField
            write(3,*), "Warning. Integrated field enrgy:", &
                    mc%EField," while absolute field energy:", &
                    mc%DEField
        endif
        mc%EField=mc%DEField
        mc%x_Field=mc%dx_Field

        ! check for NaN
        do I=1,mc%NBIN
            if (abs(md%Vol(I)).lt.0.00001) Cycle
            if (isnan(md%PHIA(I))) then
                write(*,"(A,I5,A)"), "PHIA(",I,")=NaN"
                write(*,"(A,I5,A,f8.4)"), "Vol(",I,")=",md%Vol(I)
                stop 1
            endif
            if (isnan(md%PHIB(I))) then
                write(*,"(A,I5,A)"), "PHIB(",I,")=NaN"
                write(*,"(A,I5,A,f8.4)"), "Vol(",I,")=",md%Vol(I)
                stop 1
            endif
            if (isnan(md%Vol(I))) then
                write(*,"(A,I5,A)"), "Vol(",I,")=NaN"
                stop 1
            endif
        enddo

    else
        do I=1,mc%NBIN
             md%PHIA(I)=0.0_dp
             md%PHIB(I)=0.0_dp
        enddo
    endif
    close (3)
  IF (mc%RING .EQ. 1) then
     ! --- Initial Writhe
     call WRITHE(md%R,mc%N,md%Wr)

     !     Initialize the Cross matrix

     CrossSize=N**2
     ALLOCATE(Cross(CrossSize,6))
     ALLOCATE(CrossP(CrossSize,6))


     !     Get initial value of Alexander polynomial and Cross matrix
     NCross=0
     CALL ALEXANDERP(md%R,mc%N,DELTA,Cross,CrossSize,NCross)
     !     Begin Monte Carlo simulation
  ENDIF
! -------------------------------------
!
!   Begin Monte Carlo simulation
!
! -------------------------------------
    ISTEP=1
    DO WHILE (ISTEP.LE.NSTEP)

       DO MCTYPE=1,mc%moveTypes

          if (mc%MOVEON(MCTYPE).EQ.0) cycle

          ! Turn down poor moves
          if ((mc%PHit(MCTYPE).lt.mc%MIN_ACCEPT).and. &
              (mod(ISTEP,mc%reduce_move).ne.0).and. &
              ((MCTYPE.eq.5).or.(MCTYPE.eq.6))) then
              CYCLE
          endif

          call MC_move(md%R,md%U,md%RP,md%UP,mc%NT,mc%NB,mc%NP, &
                       IP,IB1,IB2,IT1,IT2,MCTYPE, &
                       mc%MCAMP,mc%WINDOW,md%AB,md%ABP,mc%G,&
                       rand_stat, mc%winType,IT3,IT4,forward)

!   Calculate the change in compression and bending energy
          if ((MCTYPE.NE.5) .and. &
              (MCTYPE.NE.6) .and. &
              (MCTYPE.NE.7) .and. &
              (MCTYPE.NE.8) .and. &
              (MCTYPE.NE.9) .and. &
              (MCTYPE.NE.10) )then
              call MC_eelas(mc%DEELAS,md%R,md%U,md%RP,md%UP,&
                            mc%NT,mc%NB,IB1,IB2, &
                            IT1,IT2,EB,EPAR,EPERP,GAM,ETA)
          else
              mc%DEELAS(1)=0.0
              mc%DEELAS(2)=0.0
              mc%DEELAS(3)=0.0
          endif
          if (MCTYPE.eq.8) then
              print*, "Flop move not working!  Chain energy isn't symmetric"
              stop 1
          endif
!   Calculate the change in the binding energy
          if (MCTYPE.EQ.7) then
              !print*, 'MCsim says EM:',EM,'EU',EU
              call MC_bind(mc%NT,mc%G,IT1,IT2,md%AB,md%ABP,md%METH,mc%EU,mc%EM, &
                           mc%DEBind,mc%mu,mc%dx_mu)
          else
              mc%DEBind=0.0
          endif
          if (INTERP_BEAD_LENNARD_JONES.EQ.1)
              !call MC_self(DESELF,md%R,md%U,md%RP,md%UP,mc%NT,mc%NB,mc%NP,IP,IB1,IB2,IT1,IT2,LHC,VHC,LBOX,GAM)
           IF (MCTYPE.EQ.1) THEN
              CALL DE_SELF_CRANK(DESELF,md%R,md%RP,mc%NT,mc%N,mc%NP,mc%PARA,mc%RING,IB1,IB2)

           ELSEIF (MCTYPE.EQ.2) THEN
              CALL ENERGY_SELF_SLIDE(ESELF,md%R,mc%NT,mc%N,mc%NP,mc%PARA,mc%RING,IB1,IB2)
              CALL ENERGY_SELF_SLIDE(ESELFP,md%R,mc%NT,mc%N,mc%NP,mc%PARA,mc%RING,IB1,IB2)

              DESELF=ESELFP-ESELF
           ELSEIF (MCTYPE.EQ.3) THEN
               CALL DE_SELF_CRANK(DESELF,md%R,md%RP,mc%NT,mc%N,mc%NP,mc%PARA,mc%RING,IB1,IB2)
           ELSEIF (MCTYPE.EQ.10) THEN
               PRINT *, 'Nobody has used this branch before. write a DE_SELF_CRANK to calculate change in self-interaction energy from this move, sorry!'
               STOP 1
           ELSE
              DESELF=0.
           ENDIF
          endif

!   Calculate the change in the self-interaction energy (actually all
!   interation energy, not just self?)
          if (FIELD_INTERACTIONS.EQ.1) then
             if (MCTYPE.EQ.9) then
                 !skip if doesn't do anything
                 if (abs(mc%CHI_ON).lt.0.00001) CYCLE
                 call MC_int_swap(mc,md,IT1,IT2,IT3,IT4)
                 if (abs(mc%DEKap).gt.0.0001) then
                     print*, "Error in MCsim.  Kappa energy shouldn't change on move 9"
                     print*, "DEKap", mc%DEKap
                     stop 1
                 endif
             elseif (MCTYPE.EQ.10) then
                 call MC_int_rep(mc,md,IT1,IT2,forward)
             else
                 call MC_int(mc,md,IT1,IT2,.false.)
             endif
          else
              mc%DEKap=0.0_dp
              mc%DECouple=0.0_dp
              mc%DEChi=0.0_dp
              mc%DEField=0.0_dp
          endif
          if ((MCTYPE.eq.8).and.(mc%DEKap.gt.0.00001)) then
              print*, "Error in MCsim. Kappa energy shouldn't change on move 8"
          endif

!   Calculate the change in confinement energy
          if ((MCTYPE.NE.7).and. &
              (MCTYPE.NE.8).and. &
              (MCTYPE.NE.9)) then
              call MC_confine(mc%confineType, mc%LBox, md%RP, mc%NT, &
                              IT1,IT2,mc%ECon)
          else
              mc%ECon=0.0_dp;
          endif

!   Change the position if appropriate
          ENERGY=mc%DEELAS(1)+mc%DEELAS(2)+mc%DEELAS(3) &
                 +mc%DEKap+mc%DECouple+mc%DEChi+mc%DEBind+mc%ECon+mc%DEField
          PROB=exp(-ENERGY)
          call random_number(urnd,rand_stat)
          TEST=urnd(1)
          if (TEST.LE.PROB) then
             if(MCTYPE.EQ.7) then
                 DO I=IT1,IT2
                      md%AB(I)=md%ABP(I)
                 ENDDO
             else
                 DO I=IT1,IT2
                     md%R(I,1)=md%RP(I,1)
                     md%R(I,2)=md%RP(I,2)
                     md%R(I,3)=md%RP(I,3)
                     md%U(I,1)=md%UP(I,1)
                     md%U(I,2)=md%UP(I,2)
                     md%U(I,3)=md%UP(I,3)
                 enddo
                 if (MCTYPE.EQ.9) then
                     DO I=IT3,IT4
                         md%R(I,1)=md%RP(I,1)
                         md%R(I,2)=md%RP(I,2)
                         md%R(I,3)=md%RP(I,3)
                         md%U(I,1)=md%UP(I,1)
                         md%U(I,2)=md%UP(I,2)
                         md%U(I,3)=md%UP(I,3)
                     enddo
                 endif
             endif
             if (mc%ECon.gt.0.0_dp) then
                 print*, "MCTYPE", MCType
                 call wlcsim_params_printEnergies(mc)
                 print*, "error in MCsim, out of bounds "
                 stop 1
             endif
             mc%EBind=mc%EBind+mc%DEBind
             mc%x_mu=mc%x_mu+mc%dx_mu
             mc%EELAS(1)=mc%EELAS(1)+mc%DEELAS(1)
             mc%EELAS(2)=mc%EELAS(2)+mc%DEELAS(2)
             mc%EELAS(3)=mc%EELAS(3)+mc%DEELAS(3)
             if (INTON) then
                DO I=1,mc%NPHI
                   J=md%INDPHI(I)
                   md%PHIA(J)=md%PHIA(J)+md%DPHIA(I)
                   md%PHIB(J)=md%PHIB(J)+md%DPHIB(I)
                   if ((md%PHIA(J).lt.-0.000001_dp) .or. (md%PHIB(J).lt.-0.00001_dp)) then
                       print*, "Error in MCsim. Negitive phi"
                       stop 1
                   endif
                enddo
                mc%ECouple=mc%ECouple+mc%DECouple
                mc%EKap=mc%EKap+mc%DEKap
                mc%EChi=mc%EChi+mc%DEChi
                mc%EField=mc%EField+mc%DEField
                mc%x_Couple=mc%x_couple+mc%dx_couple
                mc%x_kap=mc%x_Kap+mc%dx_kap
                mc%x_chi=mc%x_chi+mc%dx_chi
                mc%x_field=mc%x_field+mc%dx_field

             endif
             mc%SUCCESS(MCTYPE)=mc%SUCCESS(MCTYPE)+1
          endif
!   Adapt the amplitude of step every NADAPT steps

          !amplitude and window adaptations
          if (mod(ISTEP,mc%NADAPT(MCTYPE)).EQ.0) then  ! Addapt ever NADAPT moves
             call wlcsim_params_adapt(mc,MCTYPE)

             ! move each chain back if drifted though repeated BC
             if (mc%recenter_on) then
                 call wlcsim_params_recenter(mc,md)  ! You don't need to do this if there is confinement
            endif
          endif

       enddo ! End of movetype loop

       !  -----  Parallel tempering ----
       IF (mod(ISTEP,mc%NPT).eq.0) THEN
          call replicaExchange(mc)
       ENDIF

       ! seps in this subroutine
       ISTEP=ISTEP+1
    enddo ! end of ISTEP loop

    RETURN
END

!-------------------------------------------------------------*
