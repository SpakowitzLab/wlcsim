#include "../defines.inc"
!---------------------------------------------------------------*
!
!
!
!     This subroutine performs a Monte Carlo simulation on the
!     polymer chain.
!
!    Quinn Made Changes to this file starting on 12/15/15
!

subroutine MCsim(wlc_p)
! values from wlcsim_data
use params, only: wlc_PHit, wlc_CrossP, wlc_dx_Externalfield, wlc_ABP, wlc_WR&
    , wlc_EExternalField, wlc_EChi, wlc_DEExplicitBinding, wlc_deMu, wlc_x_Externalfield, wlc_DPHI_l2&
    , wlc_DEKap, wlc_AB, wlc_DEExternalField, wlc_NCross, wlc_x_Kap, wlc_x_chi&
    , wlc_dx_chi, wlc_EKap, wlc_DEBind, wlc_ESELF, wlc_ind_exchange, wlc_inDPHI&
    , wlc_dx_field, wlc_PHIA, wlc_x_kap, wlc_x_couple, wlc_rand_stat, wlc_Cross&
    , wlc_demu, wlc_U, wlc_Vol, wlc_PHI_l2, wlc_DECouple, wlc_NPHI&
    , wlc_x_Couple, wlc_DPHIB, wlc_x_field, wlc_ECon, wlc_DESELF, wlc_ATTEMPTS&
    , wlc_UP, wlc_dx_maierSaupe, wlc_x_mu, wlc_dx_couple, wlc_EmaierSaupe, wlc_CrossSize&
    , wlc_ECouple, wlc_demaierSaupe, wlc_EBind, wlc_DEMu, wlc_NCrossP, wlc_DEField&
    , wlc_DEElas, wlc_deMaierSaupe, wlc_R, wlc_EField, wlc_SUCCESS, wlc_dx_kap&
    , wlc_RP, wlc_METH, wlc_spiders, wlc_EElas, wlc_DPHIA, wlc_x_maierSaupe&
    , wlc_DEChi, wlc_DEELAS, wlc_EMu, wlc_eExplicitBinding, wlc_x_ExternalField, wlc_dx_mu&
    , wlc_PHIB

    !use mt19937, only : grnd, sgrnd, rnorm, mt, mti
    use mersenne_twister
    use params
    use binning, only: addBead, removeBead
    use updateRU, only: updateR

    implicit none

    !integer, intent(in) :: NSTEP             ! Number of MC steps

!   Variables for the simulation

    integer ISTEP             ! Current MC step index
    real(dp) PROB     ! Calculated test prob
    real(dp) TEST     ! Random test variable
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
    real(dp) urnd(1) ! single random number
!   Load the input parameters
    Type(wlcsim_params), intent(inout) :: wlc_p      ! system varibles
    integer DELTA             !Alexander polynomial evaluated at t = -1; used for knot checking
    real(dp) para(10)
    integer m_index  ! m is the m from spherical harmonics (z component)
    integer sweepIndex
    logical in_confinement
    logical collide
    logical success
    integer section_n, spider_id
    logical positions_have_changed

    !TODO: unpack parameters in MC_elas
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

    do while (ISTEP <= WLC_P__STEPSPEREXCHANGE)
       do MCTYPE = 1,nMoveTypes
       if (wlc_p%MOVEON(MCTYPE) == 0) cycle
       do sweepIndex = 1,wlc_p%MOVESPERSTEP(MCTYPE)
          wlc_ECon = 0.0_dp
          wlc_DEBind = 0.0_dp
          wlc_DEMu = 0.0_dp
          wlc_dx_mu = 0.0_dp
          wlc_DESELF=0.0_dp
          wlc_DEKap = 0.0_dp
          wlc_DECouple = 0.0_dp
          wlc_DEChi = 0.0_dp
          wlc_DEField = 0.0_dp
          wlc_DEExternalField = 0.0_dp
          wlc_deMaierSaupe = 0.0_dp
          wlc_DEElas=0.0_dp
          wlc_DEExplicitBinding = 0.0_dp

          ! Turn down poor moves
          if ((wlc_PHit(MCTYPE).lt.WLC_P__MIN_ACCEPT).and. &
              (mod(ISTEP,WLC_P__REDUCE_MOVE).ne.0).and. &
              ((MCTYPE.eq.5).or.(MCTYPE.eq.6))) then
              CYCLE
          endif
          call MC_move(IB1,IB2,IT1,IT2,IT3,IT4,&
                       MCTYPE,forward,wlc_rand_stat,dib,spider_id,success)
          if (.not. success) then
              wlc_ATTEMPTS(MCTYPE) = wlc_ATTEMPTS(MCTYPE) + 1
              cycle
          endif
          if ((MCTYPE == 4) .or. (MCTYPE == 7) .or. (MCTYPE == 8) ) then
              positions_have_changed = .False.
          else
              positions_have_changed = .True.
          endif

!   Calculate the change in confinement energy
          if (positions_have_changed .and. &
              (MCTYPE /= 9).and. &
              (MCTYPE /= 12)) then
              !call MC_confine(wlc_RP, WLC_P__NT,IT1,IT2,wlc_ECon)
              ! Completely skip move if outside confinement
              if (.not. in_confinement(wlc_RP, WLC_P__NT, IT1, IT2)) then
                  wlc_ATTEMPTS(MCTYPE) = wlc_ATTEMPTS(MCTYPE) + 1
                  cycle
              endif
          elseif (MCTYPE == 12) then
              do section_n = 1, wlc_spiders(spider_id)%nSections
                  IT1 = wlc_spiders(spider_id)%moved_sections(1,section_n)
                  IT2 = wlc_spiders(spider_id)%moved_sections(2,section_n)
                  if (.not. in_confinement(wlc_RP, WLC_P__NT, IT1, IT2)) then
                      wlc_ATTEMPTS(MCTYPE) = wlc_ATTEMPTS(MCTYPE) + 1
                      success = .False.
                      exit
                  endif
              enddo
              if (.not. success) cycle
          endif

          if(WLC_P__CYLINDRICAL_CHAIN_EXCLUSION) then
              call MC_cylinder(wlc_p,collide,IB1,IB2,IT1,IT2,MCTYPE,forward)
              if (collide) then
                  wlc_ATTEMPTS(MCTYPE) = wlc_ATTEMPTS(MCTYPE) + 1
                  cycle
              endif
          endif

          if (WLC_P__RING) then
              wlc_CrossP = wlc_Cross
              wlc_NCrossP = wlc_NCross
              if (MCTYPE == 1) then
                 CALL alexanderp_crank(wlc_p,wlc_RP,DELTA,wlc_CrossP,wlc_CrossSize,wlc_NCrossP,IT1,IT2,DIB)
              elseif (MCTYPE == 2) then
                 if (DIB /= WLC_P__NB) then
                    CALL alexanderp_slide(wlc_p,wlc_RP,DELTA,wlc_CrossP,wlc_CrossSize,wlc_NCrossP,IT1,IT2,DIB)
                 ENDif
              else
                 CALL ALEXANDERP(wlc_RP,WLC_P__NB,DELTA,wlc_CrossP,wlc_CrossSize,wlc_NCrossP)
              ENDif
              if (DELTA /= 1) then
                 wlc_ATTEMPTS(MCTYPE) = wlc_ATTEMPTS(MCTYPE) + 1
                 cycle
              ENDif
          ENDif


!   Calculate the change in compression and bending energy
          if (MCTYPE<5) then
              call MC_eelas(wlc_p,IB1,IB2,IT1,IT2,EB,EPAR,EPERP,GAM,ETA,MCTYPE,WRP)
          elseif (MCTYPE==12) then
              call MC_eelas_spider(wlc_p,wlc_DEELAS,spider_id,&
                                   EB,EPAR,EPERP,GAM,ETA)
          endif


          if (MCTYPE.eq.8) then
              print*, "Flop move not working!  Chain energy isn't symmetric"
              stop 1
          endif
!   Calculate the change in the binding energy
          if (MCTYPE == 7 .or. MCTYPE == 11) then
              !print*, 'MCsim says EM:',EM,'EU',EU
              call MC_bind(wlc_p,IT1,IT2,wlc_AB,wlc_ABP,wlc_METH,&
                           wlc_DEBind,wlc_dx_mu,wlc_demu)
          endif
          if (WLC_P__INTERP_BEAD_LENNARD_JONES) then
              !call MC_self(DESELF,wlc_R,wlc_U,wlc_RP,wlc_UP,WLC_P__NT,WLC_P__NB,WLC_P__NP,IP,IB1,IB2,IT1,IT2,LHC,VHC,LBOX,GAM)
              if (MCTYPE == 1) then
                  CALL DE_SELF_CRANK(wlc_DESELF,wlc_R,wlc_RP,WLC_P__NT,WLC_P__NB,WLC_P__NP, &
                      para,WLC_P__RING,IB1,IB2)

              elseif (MCTYPE == 2) then
                  CALL ENERGY_SELF_SLIDE(wlc_ESELF,wlc_R,WLC_P__NT,WLC_P__NB,WLC_P__NP, &
                      para,WLC_P__RING,IB1,IB2)
                  CALL ENERGY_SELF_SLIDE(ESELFP,wlc_R,WLC_P__NT,WLC_P__NB,WLC_P__NP, &
                      para,WLC_P__RING,IB1,IB2)

                  wlc_DESELF = ESELFP-wlc_ESELF
              elseif (MCTYPE == 3) then
                  CALL DE_SELF_CRANK(wlc_DESELF,wlc_R,wlc_RP,WLC_P__NT,WLC_P__NB,WLC_P__NP,&
                      para,WLC_P__RING,IB1,IB2)
              elseif (MCTYPE == 10) then
                  PRinT *, 'Nobody has used this branch before. write a DE_SELF_CRANK '
                  PRinT *, 'to calculate change in self-interaction energy from this move, sorry!'
                  STOP 1
              ENDif
          endif

!   Calculate the change in the self-interaction energy (actually all
!   interation energy, not just self?)
          if (wlc_p%field_int_on_currently .and. WLC_P__FIELD_INT_ON) then
             if (MCTYPE == 9) then !swap move
                 !skip if doesn't do anything
                 if (abs(wlc_p%CHI_ON).lt.0.00001_dp) CYCLE
                 call MC_int_swap(wlc_p,IT1,IT2,IT3,IT4)
             elseif (MCTYPE == 7) then
                 call MC_int_chem(wlc_p,IT1,IT2)
             elseif (MCTYPE == 10) then ! reptation move
                 call MC_int_rep(wlc_p,IT1,IT2,forward)
             elseif (MCTYPE == 11) then ! super reptation move
                 call MC_int_super_rep(wlc_p,IT1,IT2,forward)
             elseif (MCTYPE == 12) then
                 call MC_int_update_spider(wlc_p,spider_id)
             else ! motion of chain
                 call MC_int_update(wlc_p,IT1,IT2)
             endif
          endif

          if (WLC_P__APPLY_EXTERNAL_FIELD .and. positions_have_changed) then
              if (MCTYPE == 12) then
                  call MC_external_field_spider(wlc_p,spider_id)
              else
                  wlc_dx_Externalfield = 0.0_dp
                  call MC_external_field(wlc_p,IT1,IT2)
              endif
          endif

          if (WLC_P__EXPLICIT_BINDING .and. positions_have_changed) then
              if (MCTYPE == 12) then
                  call MC_excplicit_binding_spider(wlc_p,spider_id)
              else
                  call MC_explicit_binding(IT1,IT2,MCTYPE)
              endif
          endif

!   Change the position if appropriate
          ENERGY = wlc_DEElas(1) + wlc_DEElas(2) + wlc_DEElas(3) &
                 + wlc_DEKap + wlc_DECouple + wlc_DEChi + wlc_DEBind &
                 + wlc_deMu &
                 + wlc_ECon + wlc_DEField &
                 + wlc_DEExternalField &
                 + wlc_deMaierSaupe &
                 + wlc_DEExplicitBinding
          PROB = exp(-ENERGY)
          call random_number(urnd,wlc_rand_stat)
          TEST = urnd(1)
          if (TEST <= PROB) then
             if(MCTYPE == 7 .or. MCTYPE == 11) then
                 if (.not.WLC_P__CHANGINGCHEMICALIDENTITY) then
                     call stop_if_err(1, "Tried to change chemical Identity when you can't")
                 endif
                 do I = IT1,IT2
                      wlc_AB(I) = wlc_ABP(I)
                 ENDdo
             endif
             if(MCTYPE /= 7 .and. MCTYPE /= 12) then
                 do I = IT1,IT2
                     call updateR(I)
                 enddo
                 if (MCTYPE == 9) then
                     do I = IT3,IT4
                         call updateR(I)
                     enddo
                 endif
             elseif(MCTYPE == 12) then
                 do section_n = 1, wlc_spiders(spider_id)%nSections
                     IT1 = wlc_spiders(spider_id)%moved_sections(1,section_n)
                     IT2 = wlc_spiders(spider_id)%moved_sections(2,section_n)
                     do I = IT1,IT2
                         call updateR(I)
                     enddo
                 enddo
             endif
             if (wlc_ECon.gt.0.0_dp) then
                 print*, "MCTYPE", MCType
                 call printEnergies()
                 print*, "error in MCsim, out of bounds "
                 stop 1
             endif
             wlc_EBind = wlc_EBind + wlc_DEBind
             wlc_EMu = wlc_EMu + wlc_DEMu
             wlc_x_mu = wlc_x_mu + wlc_dx_mu
             wlc_eExplicitBinding = wlc_eExplicitBinding + wlc_DEExplicitBinding
             wlc_EElas = wlc_EElas + wlc_DEElas
             if ((MCTYPE .ne. 4) .and. (MCTYPE .ne. 7) .and. &
                 (MCTYPE .ne. 8) .and. (MCTYPE .ne. 9) .and. &
                 WLC_P__APPLY_EXTERNAL_FIELD) then
                 wlc_EExternalField = wlc_EExternalField + wlc_DEExternalField
                 wlc_x_ExternalField = wlc_x_Externalfield + wlc_dx_Externalfield
             endif
             if (wlc_p%field_int_on_currently .and. WLC_P__FIELD_INT_ON) then
                do I = 1,wlc_NPHI
                   J = wlc_inDPHI(I)
                   if (wlc_p%CHI_L2_ON) then
                       do m_index = -2,2
                           wlc_PHI_l2(m_index,J) =  wlc_PHI_l2(m_index,J) + wlc_DPHI_l2(m_index,I)
                       enddo
                   endif
                   wlc_PHIA(J) = wlc_PHIA(J) + wlc_DPHIA(I)
                   wlc_PHIB(J) = wlc_PHIB(J) + wlc_DPHIB(I)

                   if ((wlc_PHIA(J).lt.-0.0001_dp) .or. &
                       (wlc_PHIB(J).lt.-0.00001_dp .and. (.not. WLC_P__TWO_TAIL))) then
                       print*, "IT1-4",IT1,IT2,IT3,IT4
                       if(WLC_P__FRACTIONAL_BIN) print*, "Vol", wlc_Vol(I)
                       print*, "MCTYPE", MCTYPE
                       print*, "DPHIA ",wlc_DPHIA(I)," DPHIB",wlc_DPHIB(I)
                       print*, "PHIA(J) ", wlc_PHIA(J), " PHIB(J) ", wlc_PHIB(J)
                       print*, "I", I,"J",J
                       print*, "Error in MCsim. Negative phi"
                       stop 1
                   endif
                enddo
                wlc_ECouple = wlc_ECouple + wlc_DECouple
                wlc_EKap = wlc_EKap + wlc_DEKap
                wlc_EChi = wlc_EChi + wlc_DEChi
                wlc_EField = wlc_EField + wlc_DEField
                wlc_EmaierSaupe = wlc_EmaierSaupe + wlc_demaierSaupe

                wlc_x_Couple = wlc_x_couple + wlc_dx_couple
                wlc_x_kap = wlc_x_Kap + wlc_dx_kap
                wlc_x_chi = wlc_x_chi + wlc_dx_chi
                wlc_x_field = wlc_x_field + wlc_dx_field
                wlc_x_maierSaupe = wlc_x_maierSaupe + wlc_dx_maierSaupe

             endif
             if (WLC_P__RING) then
                wlc_WR = WRP
                wlc_NCross = wlc_NCrossP
                wlc_Cross = wlc_CrossP
            endif
             wlc_SUCCESS(MCTYPE) = wlc_SUCCESS(MCTYPE) + 1
          endif
          wlc_ATTEMPTS(MCTYPE) = wlc_ATTEMPTS(MCTYPE) + 1

          !  vvvvvvvvvv Beginning of hold check vvvvvvvvvvvvvvvvvv
          do I = IT1,IT2
              wlc_RP(1,I) = nan
              wlc_RP(2,I) = nan
              wlc_RP(3,I) = nan
              wlc_UP(1,I) = nan
              wlc_UP(2,I) = nan
              wlc_UP(3,I) = nan
              if (WLC_P__VARIABLE_CHEM_STATE) then
                  wlc_ABP(I) = INT_MIN
              endif
          enddo
          if (MCTYPE == 9) then
              do I = IT3,IT4
                  wlc_RP(1,I) = nan
                  wlc_RP(2,I) = nan
                  wlc_RP(3,I) = nan
                  wlc_UP(1,I) = nan
                  wlc_UP(2,I) = nan
                  wlc_UP(3,I) = nan
              enddo
          endif
          !^^^^^^^^^^^End of check ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!   Adapt the amplitude of step every NADAPT steps
       enddo ! End of sweepIndex loop
          !amplitude and window adaptations
          if (mod(ISTEP+wlc_ind_exchange*WLC_P__STEPSPEREXCHANGE,wlc_p%NADAPT(MCTYPE)) == 0) then  ! Addapt ever NADAPT moves
             call mc_adapt(wlc_p,MCTYPE)

             ! move each chain back if drifted though repeated BC
             if (WLC_P__RECENTER_ON) then
                 call wlcsim_params_recenter()  ! You don't need to do this if there is confinement
            endif
          endif
       enddo ! End of movetype loop

      ! Parallel tempereing used to happen here but now has been moved to wlcsim_name
      ! !  -----  Parallel tempering ----
      ! if (mod(ISTEP,WLC_P__NPT).eq.0) then
      !    call replicaExchange(wlc_p)
      ! ENDif

       ! seps in this subroutine
       ISTEP = ISTEP + 1
    enddo ! end of ISTEP loop

    RETURN
END

!-------------------------------------------------------------*
