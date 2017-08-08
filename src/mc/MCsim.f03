!---------------------------------------------------------------*
!
!
!
!     This subroutine performs a Monte Carlo simulation on the
!     polymer chain.
!
!    Quinn Made Changes to this file starting on 12/15/15
!

subroutine MCsim(wlc_p,wlc_d,NSTEP)

    !use mt19937, only : grnd, sgrnd, rnorm, mt, mti
    use mersenne_twister
    use params

    implicit none

    integer, intent(in) :: NSTEP             ! Number of MC steps

!   Variables for the simulation

    integer ISTEP             ! Current MC step index
    real(dp) PROB     ! Calculated test prob
    real(dp) TEST     ! Random test variable
    integer IP                ! Test polymer
    integer IB1               ! Test bead position 1
    integer IT1               ! Index of test bead 1
    integer IB2               ! Test bead position 2
    integer IT2               ! Index of test bead 2
    integer IT3, IT4          ! second polymer for polymer swap
    integer dib               ! number of beads moved in a move
    logical forward           ! direction of reptation move

    integer I,J

    integer MCTYPE                    ! Type of MC move

    real(dp) EB,EPAR,EPERP,ESELFP
    real(dp) GAM,ETA
    real(dp) XIR,XIU
    real(dp) LHC      ! Length of HC int
    real(dp) VHC      ! HC strength
!    real(dp) phiTot  ! for testing
    real(dp) wrp ! proposed writhe

    real(dp) ENERGY
! Things for random number generator
    real urnd(1) ! single random number
    type(random_stat) :: rand_stat
!   Load the input parameters
    Type(wlcsim_params), intent(inout) :: wlc_p      ! system varibles
    Type(wlcsim_data), intent(inout) :: wlc_d     ! system allocated data
    integer DELTA             !Alexander polynomial evaluated at t = -1; used for knot checking
    real(dp) para(10)
    integer m_plus3

    rand_stat = wlc_d%rand_stat


    para = pack_as_para(wlc_p)
    EB=   PARA(1)
    EPAR= PARA(2)
    EPERP = PARA(3)
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
    ISTEP = 1

    do while (ISTEP <= NSTEP)

       do MCTYPE = 1,wlc_p%moveTypes
          if (wlc_p%MOVEON(MCTYPE) == 0) cycle

          ! Turn down poor moves
          if ((wlc_d%PHit(MCTYPE).lt.wlc_p%Min_ACCEPT).and. &
              (mod(ISTEP,wlc_p%reduce_move).ne.0).and. &
              ((MCTYPE.eq.5).or.(MCTYPE.eq.6))) then
              CYCLE
          endif
          call MC_move(wlc_d%R,wlc_d%U,wlc_d%RP,wlc_d%UP,wlc_p%NT,wlc_p%NB,wlc_p%NP, &
                       IP,IB1,IB2,IT1,IT2,MCTYPE, &
                       wlc_d%MCAMP,wlc_d%WindoW,wlc_p%nBpM,&
                       rand_stat, wlc_p%winType,IT3,IT4,forward,dib,wlc_p%ring, &
                       wlc_p%inTERP_BEAD_LENNARD_JONES,wlc_d)

        if (wlc_p%RinG) then
           wlc_d%CrossP = wlc_d%Cross
           wlc_d%NCrossP = wlc_d%NCross
           if (MCTYPE == 1) then
              CALL alexanderp_crank(wlc_d%RP,wlc_p%NB,DELTA,wlc_d%CrossP,wlc_d%CrossSize,wlc_d%NCrossP,IT1,IT2,DIB)
           elseif (MCTYPE == 2) then
              if (DIB /= wlc_p%NB) then
                 CALL alexanderp_slide(wlc_d%RP,wlc_p%NB,DELTA,wlc_d%CrossP,wlc_d%CrossSize,wlc_d%NCrossP,IT1,IT2,DIB)
              ENDif
           else
              CALL ALEXANDERP(wlc_d%RP,wlc_p%NB,DELTA,wlc_d%CrossP,wlc_d%CrossSize,wlc_d%NCrossP)
           ENDif
           if (DELTA /= 1) then
              wlc_d%eKnot = inf
           else
               wlc_d%eKnot = 0.0_dp
           ENDif
        else
            wlc_d%eKnot=0.0_dp
        ENDif


!   Calculate the change in compression and bending energy
          if ((MCTYPE /= 5) .and. &
              (MCTYPE /= 6) .and. &
              (MCTYPE /= 7) .and. &
              (MCTYPE /= 8) .and. &
              (MCTYPE /= 9) .and. &
              (MCTYPE /= 10) )then
              call MC_eelas(wlc_d%DEElas,wlc_d%R,wlc_d%U,wlc_d%RP,wlc_d%UP,&
                            wlc_p%NT,wlc_p%NB,IB1,IB2, &
                            IT1,IT2,EB,EPAR,EPERP,GAM,ETA, &
                            wlc_p%ring,wlc_p%twist,wlc_p%lk,wlc_p%lt,wlc_p%l, &
                            mctype,wlc_d%wr,wrp,wlc_p%simType)
          else
              wlc_d%DEElas(1) = 0.0
              wlc_d%DEElas(2) = 0.0
              wlc_d%DEElas(3) = 0.0
          endif
          if (MCTYPE.eq.8) then
              print*, "Flop move not working!  Chain energy isn't symmetric"
              stop 1
          endif
!   Calculate the change in the binding energy
          if (MCTYPE == 7) then
              !print*, 'MCsim says EM:',EM,'EU',EU
              call MC_bind(wlc_p%NT,wlc_p%NB,IT1,IT2,wlc_d%AB,wlc_d%ABP,wlc_d%METH,wlc_p%EU,wlc_p%EM, &
                          wlc_d%DEBind,wlc_p%mu,wlc_d%dx_mu)
          else
              wlc_d%DEBind = 0.0
          endif
          if (wlc_p%inTERP_BEAD_LENNARD_JONES) then
              !call MC_self(DESELF,wlc_d%R,wlc_d%U,wlc_d%RP,wlc_d%UP,wlc_p%NT,wlc_p%NB,wlc_p%NP,IP,IB1,IB2,IT1,IT2,LHC,VHC,LBOX,GAM)
              if (MCTYPE == 1) then
                  CALL DE_SELF_CRANK(wlc_d%DESELF,wlc_d%R,wlc_d%RP,wlc_p%NT,wlc_p%NB,wlc_p%NP, &
                      pack_as_para(wlc_p),wlc_p%RinG,IB1,IB2)

              elseif (MCTYPE == 2) then
                  CALL ENERGY_SELF_SLIDE(wlc_d%ESELF,wlc_d%R,wlc_p%NT,wlc_p%NB,wlc_p%NP, &
                      pack_as_para(wlc_p),wlc_p%RinG,IB1,IB2)
                  CALL ENERGY_SELF_SLIDE(ESELFP,wlc_d%R,wlc_p%NT,wlc_p%NB,wlc_p%NP, &
                      pack_as_para(wlc_p),wlc_p%RinG,IB1,IB2)

                  wlc_d%DESELF = ESELFP-wlc_d%ESELF
              elseif (MCTYPE == 3) then
                  CALL DE_SELF_CRANK(wlc_d%DESELF,wlc_d%R,wlc_d%RP,wlc_p%NT,wlc_p%NB,wlc_p%NP,&
                      pack_as_para(wlc_p),wlc_p%RinG,IB1,IB2)
              elseif (MCTYPE == 10) then
                  PRinT *, 'Nobody has used this branch before. write a DE_SELF_CRANK '
                  PRinT *, 'to calculate change in self-interaction energy from this move, sorry!'
                  STOP 1
              else
                  wlc_d%DESELF = 0.0
              ENDif
          else
              wlc_d%DESELF=0.0
          endif

!   Calculate the change in the self-interaction energy (actually all
!   interation energy, not just self?)
          if (wlc_p%FIELD_inT_ON) then
             if (MCTYPE == 9) then !swap move
                 !skip if doesn't do anything
                 if (abs(wlc_p%CHI_ON).lt.0.00001) CYCLE
                 call MC_int_swap(wlc_p,wlc_d,IT1,IT2,IT3,IT4)
             elseif (MCTYPE == 7) then
                 call MC_int_chem(wlc_p,wlc_d,IT1,IT2)
             elseif (MCTYPE == 10) then ! reptation move
                 call MC_int_rep(wlc_p,wlc_d,IT1,IT2,forward)
             else ! motion of chain
                 call MC_int_update(wlc_p,wlc_d,IT1,IT2,.false.)
             endif
          else
              wlc_d%DEKap = 0.0_dp
              wlc_d%DECouple = 0.0_dp
              wlc_d%DEChi = 0.0_dp
              wlc_d%DEField = 0.0_dp
              wlc_d%deMaierSaupe = 0.0_dp
          endif
          if ((MCTYPE.eq.8).and.(wlc_d%DEKap.gt.0.00001)) then
              print*, "Error in MCsim. Kappa energy shouldn't change on move 8"
          endif

!   Calculate the change in confinement energy
          if ((MCTYPE /= 7).and. &
              (MCTYPE /= 8).and. &
              (MCTYPE /= 9)) then
              call MC_confine(wlc_p%confineType, wlc_p%LBox, wlc_d%RP, wlc_p%NT, &
                              IT1,IT2,wlc_d%ECon)
          else
              wlc_d%ECon = 0.0_dp;
          endif

!   Change the position if appropriate
          ENERGY = wlc_d%DEElas(1) + wlc_d%DEElas(2) + wlc_d%DEElas(3) &
                 +wlc_d%DEKap + wlc_d%DECouple + wlc_d%DEChi + wlc_d%DEBind + wlc_d%ECon + wlc_d%DEField &
                 +wlc_d%eKnot + wlc_d%deMaierSaupe
          PROB = exp(-ENERGY)
          call random_number(urnd,rand_stat)
          TEST = urnd(1)
          if (TEST <= PROB) then

             if(MCTYPE == 7) then
                 if (.not.wlc_p%ChangingChemicalIdentity) then
                     call stop_if_err(1, "Tried to change chemical Identity when you can't")
                 endif
                 do I = IT1,IT2
                      wlc_d%AB(I) = wlc_d%ABP(I)
                 ENDdo
             else
                 do I = IT1,IT2
                     wlc_d%R(1,I) = wlc_d%RP(1,I)
                     wlc_d%R(2,I) = wlc_d%RP(2,I)
                     wlc_d%R(3,I) = wlc_d%RP(3,I)
                     wlc_d%U(1,I) = wlc_d%UP(1,I)
                     wlc_d%U(2,I) = wlc_d%UP(2,I)
                     wlc_d%U(3,I) = wlc_d%UP(3,I)
                 enddo
                 if (MCTYPE == 9) then
                     do I = IT3,IT4
                         wlc_d%R(1,I) = wlc_d%RP(1,I)
                         wlc_d%R(2,I) = wlc_d%RP(2,I)
                         wlc_d%R(3,I) = wlc_d%RP(3,I)
                         wlc_d%U(1,I) = wlc_d%UP(1,I)
                         wlc_d%U(2,I) = wlc_d%UP(2,I)
                         wlc_d%U(3,I) = wlc_d%UP(3,I)
                     enddo
                 endif
             endif
             if (wlc_d%ECon.gt.0.0_dp) then
                 print*, "MCTYPE", MCType
                 call printEnergies(wlc_d)
                 print*, "error in MCsim, out of bounds "
                 stop 1
             endif
             wlc_d%EBind = wlc_d%EBind + wlc_d%DEBind
             wlc_d%x_mu = wlc_d%x_mu + wlc_d%dx_mu
             wlc_d%EElas(1) = wlc_d%EElas(1) + wlc_d%DEElas(1)
             wlc_d%EElas(2) = wlc_d%EElas(2) + wlc_d%DEElas(2)
             wlc_d%EElas(3) = wlc_d%EElas(3) + wlc_d%DEElas(3)
             if (wlc_p%FIELD_inT_ON) then
                do I = 1,wlc_d%NPHI
                   J = wlc_d%inDPHI(I)
                   if (wlc_p%chi_l2_on) then
                       do m_plus3 = 1,5
                           wlc_d%PHI_l2(m_plus3,J) =  wlc_d%PHI_l2(m_plus3,J) + wlc_d%DPHI_l2(m_plus3,J)
                       enddo
                   endif
                   wlc_d%PHIA(J) = wlc_d%PHIA(J) + wlc_d%DPHIA(I)
                   wlc_d%PHIB(J) = wlc_d%PHIB(J) + wlc_d%DPHIB(I)
                   if ((wlc_d%PHIA(J).lt.-0.0001_dp) .or. (wlc_d%PHIB(J).lt.-0.00001_dp)) then
                       print*, "IT1-4",IT1,IT2,IT3,IT4
                       print*, "Vol", wlc_d%Vol(I)
                       print*, "MCTYPE", MCTYPE
                       print*, "DPHIA ",wlc_d%DPHIA(I)," DPHIB",wlc_d%DPHIB(I)
                       print*, "PHIA(J) ", wlc_d%PHIA(J), " PHIB(J) ", wlc_d%PHIB(J)
                       print*, "I", I,"J",J
                       print*, "Error in MCsim. Negative phi"
                       stop 1
                   endif
                enddo
                wlc_d%ECouple = wlc_d%ECouple + wlc_d%DECouple
                wlc_d%EKap = wlc_d%EKap + wlc_d%DEKap
                wlc_d%EChi = wlc_d%EChi + wlc_d%DEChi
                wlc_d%EField = wlc_d%EField + wlc_d%DEField
                wlc_d%EmaierSaupe = wlc_d%EmaierSaupe + wlc_d%demaierSaupe

                wlc_d%x_Couple = wlc_d%x_couple + wlc_d%dx_couple
                wlc_d%x_kap = wlc_d%x_Kap + wlc_d%dx_kap
                wlc_d%x_chi = wlc_d%x_chi + wlc_d%dx_chi
                wlc_d%x_field = wlc_d%x_field + wlc_d%dx_field
                wlc_d%x_maierSaupe = wlc_d%x_maierSaupe + wlc_d%dx_maierSaupe

             endif
             if (wlc_p%ring) then
                wlc_d%WR = WRP
                wlc_d%NCross = wlc_d%NCrossP
                wlc_d%Cross = wlc_d%CrossP
            endif
             wlc_d%SUCCESS(MCTYPE) = wlc_d%SUCCESS(MCTYPE) + 1
          endif
          wlc_d%ATTEMPTS(MCTYPE) = wlc_d%ATTEMPTS(MCTYPE) + 1
!   Adapt the amplitude of step every NADAPT steps

          !amplitude and window adaptations
          if (mod(ISTEP+wlc_d%ind_exchange*NSTEP,wlc_p%NADAPT(MCTYPE)) == 0) then  ! Addapt ever NADAPT moves
             call mc_adapt(wlc_p,wlc_d,MCTYPE)

             ! move each chain back if drifted though repeated BC
             if (wlc_p%recenter_on) then
                 call wlcsim_params_recenter(wlc_p,wlc_d)  ! You don't need to do this if there is confinement
            endif
          endif

       enddo ! End of movetype loop

      ! Parallel tempereing used to happen here but now has been moved to wlcsim_name
      ! !  -----  Parallel tempering ----
      ! if (mod(ISTEP,wlc_p%NPT).eq.0) then
      !    call replicaExchange(wlc_p)
      ! ENDif

       ! seps in this subroutine
       ISTEP = ISTEP + 1
    enddo ! end of ISTEP loop

    RETURN
END

!-------------------------------------------------------------*
