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
    integer dib               ! number of beads moved in a move
    logical forward           ! direction of reptation move

    INTEGER I,J

    INTEGER MCTYPE                    ! Type of MC move

    DOUBLE PRECISION EB,EPAR,EPERP,ESELFP
    DOUBLE PRECISION GAM,ETA
    DOUBLE PRECISION XIR,XIU
    DOUBLE PRECISION LHC      ! Length of HC int
    DOUBLE PRECISION VHC      ! HC strength
!    DOUBLE PRECISION phiTot  ! for testing
    real(dp) wrp ! proposed writhe

    DOUBLE PRECISION ENERGY
! Things for random number generator
    real urnd(1) ! single random number
    type(random_stat) :: rand_stat
!   Load the input parameters
    Type(wlcsim_params), intent(inout) :: mc      ! system varibles
    Type(wlcsim_data), intent(inout) :: md     ! system allocated data
    INTEGER DELTA             !Alexander polynomial evaluated at t=-1; used for knot checking
    real(dp) para(10)

    rand_stat = md%rand_stat


    para = pack_as_para(mc)
    EB=   PARA(1)
    EPAR= PARA(2)
    EPERP=PARA(3)
    GAM=  PARA(4)
    ETA=  PARA(5)
    XIR=  PARA(6)
    XIU=  PARA(7)
    LHC=  PARA(9)
    VHC=  PARA(10)
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
          if ((md%PHit(MCTYPE).lt.mc%MIN_ACCEPT).and. &
              (mod(ISTEP,mc%reduce_move).ne.0).and. &
              ((MCTYPE.eq.5).or.(MCTYPE.eq.6))) then
              CYCLE
          endif

          call MC_move(md%R,md%U,md%RP,md%UP,mc%NT,mc%NB,mc%NP, &
                       IP,IB1,IB2,IT1,IT2,MCTYPE, &
                       md%MCAMP,md%WINDOW,md%AB,md%ABP,mc%nBpM,&
                       rand_stat, mc%winType,IT3,IT4,forward,dib,mc%ring, &
                       mc%INTERP_BEAD_LENNARD_JONES)

        IF (mc%RING) THEN
           md%CrossP=md%Cross
           md%NCrossP=md%NCross
           IF (MCTYPE.EQ.1) THEN
              CALL alexanderp_crank(md%RP,mc%NB,DELTA,md%CrossP,md%CrossSize,md%NCrossP,IT1,IT2,DIB)
           ELSEIF (MCTYPE.EQ.2) THEN
              IF (DIB.NE.mc%NB) THEN
                 CALL alexanderp_slide(md%RP,mc%NB,DELTA,md%CrossP,md%CrossSize,md%NCrossP,IT1,IT2,DIB)
              ENDIF
           ELSE
              CALL ALEXANDERP(md%RP,mc%NB,DELTA,md%CrossP,md%CrossSize,md%NCrossP)
           ENDIF
           IF (DELTA.NE.1) THEN
              md%eKnot = inf
           else
               md%eKnot = 0.0_dp
           ENDIF
        ENDIF


!   Calculate the change in compression and bending energy
          if ((MCTYPE.NE.5) .and. &
              (MCTYPE.NE.6) .and. &
              (MCTYPE.NE.7) .and. &
              (MCTYPE.NE.8) .and. &
              (MCTYPE.NE.9) .and. &
              (MCTYPE.NE.10) )then
              call MC_eelas(md%DEElas,md%R,md%U,md%RP,md%UP,&
                            mc%NT,mc%NB,IB1,IB2, &
                            IT1,IT2,EB,EPAR,EPERP,GAM,ETA, &
                            mc%ring,mc%twist,mc%lk,mc%lt,mc%l, &
                            mctype,md%wr,wrp,mc%simtype)
          else
              md%DEElas(1)=0.0
              md%DEElas(2)=0.0
              md%DEElas(3)=0.0
          endif
          if (MCTYPE.eq.8) then
              print*, "Flop move not working!  Chain energy isn't symmetric"
              stop 1
          endif
!   Calculate the change in the binding energy
          if (MCTYPE.EQ.7) then
              !print*, 'MCsim says EM:',EM,'EU',EU
              call MC_bind(mc%NT,mc%NB,IT1,IT2,md%AB,md%ABP,md%METH,mc%EU,mc%EM, &
                          md%DEBind,mc%mu,md%dx_mu)
          else
              md%DEBind=0.0
          endif
          if (mc%INTERP_BEAD_LENNARD_JONES) then
              !call MC_self(DESELF,md%R,md%U,md%RP,md%UP,mc%NT,mc%NB,mc%NP,IP,IB1,IB2,IT1,IT2,LHC,VHC,LBOX,GAM)
              IF (MCTYPE.EQ.1) THEN
                  CALL DE_SELF_CRANK(md%DESELF,md%R,md%RP,mc%NT,mc%NB,mc%NP,pack_as_para(mc),mc%RING,IB1,IB2)

              ELSEIF (MCTYPE.EQ.2) THEN
                  CALL ENERGY_SELF_SLIDE(md%ESELF,md%R,mc%NT,mc%NB,mc%NP,pack_as_para(mc),mc%RING,IB1,IB2)
                  CALL ENERGY_SELF_SLIDE(ESELFP,md%R,mc%NT,mc%NB,mc%NP,pack_as_para(mc),mc%RING,IB1,IB2)

                  md%DESELF=ESELFP-md%ESELF
              ELSEIF (MCTYPE.EQ.3) THEN
                  CALL DE_SELF_CRANK(md%DESELF,md%R,md%RP,mc%NT,mc%NB,mc%NP,pack_as_para(mc),mc%RING,IB1,IB2)
              ELSEIF (MCTYPE.EQ.10) THEN
                  PRINT *, 'Nobody has used this branch before. write a DE_SELF_CRANK '
                  PRINT *, 'to calculate change in self-interaction energy from this move, sorry!'
                  STOP 1
              ELSE
                  md%DESELF=0.
              ENDIF
          endif

!   Calculate the change in the self-interaction energy (actually all
!   interation energy, not just self?)
          if (mc%FIELD_INT_ON) then
             if (MCTYPE.EQ.9) then
                 !skip if doesn't do anything
                 if (abs(mc%CHI_ON).lt.0.00001) CYCLE
                 call MC_int_swap(mc,md,IT1,IT2,IT3,IT4)
                 if (abs(md%DEKap).gt.0.0001) then
                     print*, "Error in MCsim.  Kappa energy shouldn't change on move 9"
                     print*, "DEKap", md%DEKap
                     stop 1
                 endif
             elseif (MCTYPE.EQ.10) then
                 call MC_int_rep(mc,md,IT1,IT2,forward)
             else
                 call MC_int_update(mc,md,IT1,IT2,.false.)
             endif
          else
              md%DEKap=0.0_dp
              md%DECouple=0.0_dp
              md%DEChi=0.0_dp
              md%DEField=0.0_dp
          endif
          if ((MCTYPE.eq.8).and.(md%DEKap.gt.0.00001)) then
              print*, "Error in MCsim. Kappa energy shouldn't change on move 8"
          endif

!   Calculate the change in confinement energy
          if ((MCTYPE.NE.7).and. &
              (MCTYPE.NE.8).and. &
              (MCTYPE.NE.9)) then
              call MC_confine(mc%confineType, mc%LBox, md%RP, mc%NT, &
                              IT1,IT2,md%ECon)
          else
              md%ECon=0.0_dp;
          endif

!   Change the position if appropriate
          ENERGY=md%DEElas(1)+md%DEElas(2)+md%DEElas(3) &
                 +md%DEKap+md%DECouple+md%DEChi+md%DEBind+md%ECon+md%DEField &
                 +md%eKnot
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
             if (md%ECon.gt.0.0_dp) then
                 print*, "MCTYPE", MCType
                 call printEnergies(md)
                 print*, "error in MCsim, out of bounds "
                 stop 1
             endif
             md%EBind=md%EBind+md%DEBind
             md%x_mu=md%x_mu+md%dx_mu
             md%EElas(1)=md%EElas(1)+md%DEElas(1)
             md%EElas(2)=md%EElas(2)+md%DEElas(2)
             md%EElas(3)=md%EElas(3)+md%DEElas(3)
             if (mc%FIELD_INT_ON) then
                DO I=1,md%NPHI
                   J=md%INDPHI(I)
                   md%PHIA(J)=md%PHIA(J)+md%DPHIA(I)
                   md%PHIB(J)=md%PHIB(J)+md%DPHIB(I)
                   if ((md%PHIA(J).lt.-0.000001_dp) .or. (md%PHIB(J).lt.-0.00001_dp)) then
                       print*, "Error in MCsim. Negitive phi"
                       stop 1
                   endif
                enddo
                md%ECouple=md%ECouple+md%DECouple
                md%EKap=md%EKap+md%DEKap
                md%EChi=md%EChi+md%DEChi
                md%EField=md%EField+md%DEField
                md%x_Couple=md%x_couple+md%dx_couple
                md%x_kap=md%x_Kap+md%dx_kap
                md%x_chi=md%x_chi+md%dx_chi
                md%x_field=md%x_field+md%dx_field

             endif
             md%WR=WRP
             md%NCross=md%NCrossP
             md%Cross=md%CrossP
             md%SUCCESS(MCTYPE)=md%SUCCESS(MCTYPE)+1
          endif
!   Adapt the amplitude of step every NADAPT steps

          !amplitude and window adaptations
          if (mod(ISTEP,mc%NADAPT(MCTYPE)).EQ.0) then  ! Addapt ever NADAPT moves
             call mc_adapt(mc,MCTYPE)

             ! move each chain back if drifted though repeated BC
             if (mc%recenter_on) then
                 call wlcsim_params_recenter(mc,md)  ! You don't need to do this if there is confinement
            endif
          endif

       enddo ! End of movetype loop

      ! Parallel tempereing used to happen here but now has been moved to wlcsim_name
      ! !  -----  Parallel tempering ----
      ! IF (mod(ISTEP,mc%NPT).eq.0) THEN
      !    call replicaExchange(mc)
      ! ENDIF

       ! seps in this subroutine
       ISTEP=ISTEP+1
    enddo ! end of ISTEP loop

    RETURN
END

!-------------------------------------------------------------*
