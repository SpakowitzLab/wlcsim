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
use params, only: wlc_PHit, wlc_CrossP, wlc_ABP &
    , wlc_DPHI_l2, wlc_AB, wlc_NCross &
    , wlc_ind_exchange, wlc_inDPHI, wlc_rand_stat, wlc_Cross&
    , wlc_Vol, wlc_PHI_l2, wlc_NPHI, wlc_DPHIB, wlc_ATTEMPTS&
    , wlc_UP, wlc_CrossSize, wlc_NCrossP, wlc_R, wlc_SUCCESS &
    , wlc_RP, wlc_METH, wlc_DPHIA, wlc_PHIB, printEnergies&
    , wlcsim_params, wlc_PHIA, int_min, NAN, wlc_nBend, wlc_nPointsMoved&
    , pack_as_para, nMoveTypes, wlc_pointsMoved, wlc_bendPoints&
    , wlcsim_params_recenter, wlc_Lk0, wlc_Lk, wlc_Tw, wlc_Wr, wlc_basepairs, wlc_nucleosomeWrap
    use energies
    use umbrella, only: umbrella_energy

    !use mt19937, only : grnd, sgrnd, rnorm, mt, mti
    use mersenne_twister
    use updateRU, only: updateR
    use polydispersity, only: length_of_chain, chain_ID, leftmost_from
    use linkingNumber, only: getDelTw_Wr_Lk

    implicit none
    interface
        pure function list_confinement()
            logical list_confinement
        end function
    end interface

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

    integer I,J, I_left

    integer MCTYPE                    ! Type of MC move

    real(dp) EB,EPAR,EPERP,ESELFP
    real(dp) GAM,ETA
    real(dp) XIR,XIU
    real(dp) LHC      ! Length of HC int
    real(dp) VHC      ! HC strength
!    real(dp) phiTot  ! for testing

    real(dp) ENERGY
! Things for random number generator
    real(dp) urnd(1) ! single random number
!   Load the input parameters
    Type(wlcsim_params), intent(inout) :: wlc_p      ! system varibles
    integer DELTA             !Alexander polynomial evaluated at t = -1; used for knot checking
    real(dp) para(10)
    integer m_index  ! m is the m from spherical harmonics (z component)
    integer sweepIndex
    logical success
    logical wlc_AlexanderP
    integer collisions

    real(dp) delTw      ! change in twist
    real(dp) delWr      ! change in writhe
    real(dp) delLk      ! change in linking number
    real(dp) TwP        ! twist of the proposed configuration
    real(dp) WrP        ! writhe of the proposed configuration
    real(dp) LkP        ! linking number of the proposed configuration

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
          call set_all_dEnergy_to_zero()
          wlc_nPointsMoved = 0
          wlc_nBend = 0

          ! ------------------------------------------
          !
          !  Each energy function adds to applicable energyOf(*_)%dx.
          !  Many functions assume enertyOf(*_)%dx starts at zero.
          !
          !-------------------------------------------

          ! Turn down poor moves
          if ((wlc_PHit(MCTYPE).lt.WLC_P__MIN_ACCEPT).and. &
              (mod(ISTEP,WLC_P__REDUCE_MOVE).ne.0).and. &
              ((MCTYPE.eq.5).or.(MCTYPE.eq.6))) then
              goto 10 ! skip move, return RP to nan
          endif
          call MC_move(IB1,IB2,IT1,IT2,IT3,IT4,&
                       MCTYPE,forward,wlc_rand_stat,dib,success, 1.0*ISTEP/WLC_P__STEPSPEREXCHANGE)
          if (.not. success) then
              wlc_ATTEMPTS(MCTYPE) = wlc_ATTEMPTS(MCTYPE) + 1
              goto 10 ! skip move, return RP to nan
          endif

!   Calculate the change in confinement energy
          if ((MCTYPE /= 4) .and. wlc_nPointsMoved > 0 .and. (MCTYPE /= 7)) then
              if (.not. list_confinement()) then
                  wlc_ATTEMPTS(MCTYPE) = wlc_ATTEMPTS(MCTYPE) + 1
                  success = .False.
                  goto 10 ! skip move, return RP to nan
              endif
          endif


! sterics check here !
          if(WLC_P__GJK_STERICS) then
            call MC_sterics(collisions,IB1,IB2,IT1,IT2,MCTYPE,forward)
            ! ascribe collision penalty
            energyOf(sterics_)%dx = collisions ! 0 
            !if (collisions > 0) then 
            !    wlc_ATTEMPTS(MCTYPE) = wlc_ATTEMPTS(MCTYPE) + 1
            !    goto 10 ! skip move, return RP to nan
            !endif
          endif
    
          call check_RP_for_NAN(success,MCTYPE)
          if (.not. success) then
              wlc_ATTEMPTS(MCTYPE) = wlc_ATTEMPTS(MCTYPE) + 1
              goto 10 ! skip move, return RP to nan
          endif

	  wlc_AlexanderP = .FALSE.

          if (WLC_P__RING .and. WLC_P__TWIST) then
	    if(wlc_AlexanderP) then !unsure if correct
              wlc_CrossP = wlc_Cross
              wlc_NCrossP = wlc_NCross
              if(MCTYPE == 1) then !was MCTYPE == 1
                  CALL alexanderp_crank(wlc_p,wlc_RP,DELTA,wlc_CrossP,wlc_CrossSize,wlc_NCrossP,IT1,IT2,DIB)
              elseif (MCTYPE == 2) then
                  if (DIB /= length_of_chain(chain_ID(IT1))) then
                     CALL alexanderp_slide(wlc_p,wlc_RP,DELTA,wlc_CrossP,wlc_CrossSize,wlc_NCrossP,IT1,IT2,DIB)
                  ENDif
              else
                  CALL ALEXANDERP(wlc_RP,WLC_P__NB,DELTA,wlc_CrossP,wlc_CrossSize,wlc_NCrossP)
              ENDif
            ENDif
            if (DELTA /= 1) then
                wlc_ATTEMPTS(MCTYPE) = wlc_ATTEMPTS(MCTYPE) + 1
                goto 10 ! skip move, return RP to nan
            ENDif
          ENDif

!   Calculate the change in compression and bending energy
          if (wlc_nBend>0) then
              call MC_eelas(wlc_p)
              if (WLC_P__RING.AND.WLC_P__TWIST) then
                  print*, "Change this to new global twist energy!!!"
                  stop
                  call MC_global_twist(IT1,IT2,MCTYPE)
              endif
          endif

          if (MCTYPE.eq.8) then
              print*, "Flop move not working!  Chain energy isn't symmetric"
              stop 1
          endif
!   Calculate the change in the binding energy
          if (WLC_P__CHANGINGCHEMICALIDENTITY .and. MCTYPE == 7 .or. MCTYPE == 11) then
              !print*, 'MCsim says EM:',EM,'EU',EU
              call MC_bind(IT1,IT2,wlc_AB,wlc_ABP,wlc_METH)
          endif
          if (WLC_P__INTERP_BEAD_LENNARD_JONES) then
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

              !call MC_self(DESELF,wlc_R,wlc_U,wlc_RP,wlc_UP,WLC_P__NT,WLC_P__NB,WLC_P__NP,IP,IB1,IB2,IT1,IT2,LHC,VHC,LBOX,GAM)
              if (MCTYPE == 1) then
                  CALL DE_SELF_CRANK(energyOf(self_)%dx,wlc_R,wlc_RP,WLC_P__NT,WLC_P__NB,WLC_P__NP, &
                      para,WLC_P__RING,IB1,IB2)

              elseif (MCTYPE == 2) then
                  CALL ENERGY_SELF_SLIDE(energyOf(self_)%x,wlc_R,WLC_P__NT,WLC_P__NB,WLC_P__NP, &
                      para,WLC_P__RING,IB1,IB2)
                  CALL ENERGY_SELF_SLIDE(ESELFP,wlc_R,WLC_P__NT,WLC_P__NB,WLC_P__NP, &
                      para,WLC_P__RING,IB1,IB2)

                  energyOf(self_)%dx = ESELFP-energyOf(self_)%x
              elseif (MCTYPE == 3) then
                  CALL DE_SELF_CRANK(energyOf(self_)%dx,wlc_R,wlc_RP,WLC_P__NT,WLC_P__NB,WLC_P__NP,&
                                     para,WLC_P__RING,IB1,IB2)
              endif
          endif

!   Calculate the change in the self-interaction energy (actually all
!   interation energy, not just self?)
          if (wlc_p%field_int_on_currently .and. WLC_P__FIELD_INT_ON) then
             if (MCTYPE == 7) then !
                 call MC_int_chem(wlc_p,IT1,IT2)
             elseif (MCTYPE == 10) then ! reptation move
                 call MC_int_rep(wlc_p,IT1,IT2,forward)
             elseif (MCTYPE == 11) then ! super reptation move
                 call MC_int_super_rep(wlc_p,IT1,IT2,forward)
             else ! motion of chain
                 call MC_int_update(wlc_p)
             endif
          endif

          if (WLC_P__APPLY_EXTERNAL_FIELD .and. energyOf(external_)%isOn &
              .and. wlc_nPointsMoved>0 .and. MCTYPE .ne. 4 .and. (MCTYPE /= 7)) then
              call MC_external_field()
          endif

          if (WLC_P__APPLY_2body_potential .and. wlc_nPointsMoved>0) then
              call MC_2bead_potential(MCTYPE)
          endif

          if (WLC_P__EXPLICIT_BINDING .and. wlc_nPointsMoved>0 .and. MCTYPE .ne. 4 .and. (MCTYPE /= 7)) then
              call MC_explicit_binding()
          endif

          if (WLC_P__UMBRELLA .and. energyOf(umbrella_)%isOn &
              .and. wlc_nPointsMoved>0 .and. MCTYPE .ne. 4 .and. (MCTYPE /= 7)) then
              call umbrella_energy()
          endif

          ! When self-crossing is not allowed, linking number must be conserved.
          if (WLC_P__NO_SELF_CROSSING) then
              call getDelTw_Wr_Lk(IB1, IB2, MCTYPE, delTw, delWr, delLk)
              TwP = wlc_Tw + delTw
              WrP = wlc_Wr + delWr
              LkP = wlc_Lk + delLk
              if (abs(LkP - wlc_Lk0) > WLC_P__NO_CROSSING_CUTOFF) then
                  wlc_ATTEMPTS(MCTYPE) = wlc_ATTEMPTS(MCTYPE) + 1
                  goto 10 ! skip move, return RP to nan
              endif
          endif

!   Change the position if appropriate
          call apply_energy_isOn()
          call calc_all_dE_from_dx()
          call sum_all_dEnergies(ENERGY)
          !call MC_save_energy_data(MCTYPE)
          PROB = exp(-ENERGY)
          call random_number(urnd,wlc_rand_stat)
          TEST = urnd(1)
          if (TEST <= PROB) then
             call accept_all_energies()
             if(MCTYPE == 7 .or. MCTYPE == 11) then
                 if (.not.WLC_P__CHANGINGCHEMICALIDENTITY) then
                     call stop_if_err(1, "Tried to change chemical Identity when you can't")
                 endif
                 do I = IT1,IT2
                      wlc_AB(I) = wlc_ABP(I)
                 ENDdo
             endif
             if(MCTYPE /= 7) then
                 do I = 1,wlc_nPointsMoved
                     J = wlc_pointsMoved(I)
                     call updateR(J)
                 enddo
             endif
             if (energyOf(confine_)%dE.gt.0.0_dp) then
                 print*, "MCTYPE", MCType
                 call printEnergies()
                 print*, "error in MCsim, out of bounds "
                 stop 1
             endif
             if (wlc_p%field_int_on_currently .and. WLC_P__FIELD_INT_ON) then
                do I = 1,wlc_NPHI
                   J = wlc_inDPHI(I)
                   if (WLC_P__CHI_L2_ABLE .and. energyOf(maierSaupe_)%isOn) then
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
             endif
             if (WLC_P__RING) then
                wlc_NCross = wlc_NCrossP
                wlc_Cross = wlc_CrossP
             endif
             if (WLC_P__NO_SELF_CROSSING) then
                wlc_Lk = LkP
                wlc_Tw = TwP
                wlc_Wr = WrP
             endif
             wlc_SUCCESS(MCTYPE) = wlc_SUCCESS(MCTYPE) + 1
          endif
          wlc_ATTEMPTS(MCTYPE) = wlc_ATTEMPTS(MCTYPE) + 1

10        continue

          !  vvvvvvvvvv Beginning of hold check vvvvvvvvvvvvvvvvvv
          if (isnan(ENERGY)) then
              print*, "Energy = NAN"
          endif
          !  It is now assumed that RP=R or nan after this.  Do not remove this loop.
          do J = 1,wlc_nPointsMoved
              I = wlc_pointsMoved(J)
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

          if (.not. WLC_P__RING) then
             do J = 1,wlc_nBend
                I = wlc_bendPoints(J)
                wlc_RP(:,I:I+1) = nan
                wlc_UP(:,I:I+1) = nan
                if (WLC_P__VARIABLE_CHEM_STATE) then
                    wlc_ABP(I:I+1) = INT_MIN
                endif
             enddo
          else
             do J = 1,wlc_nBend
                I = wlc_bendPoints(J)
                if (I == length_of_chain(chain_ID(IT1))) then
                    I_left=leftmost_from(IT1)
                    wlc_RP(:, I) = nan
                    wlc_RP(:, I_left) = nan
                    wlc_UP(:, I) = nan
                    wlc_UP(:, I_left) = nan
                else
                    wlc_RP(:,I:I+1) = nan
                    wlc_UP(:,I:I+1) = nan
                endif
                if (WLC_P__VARIABLE_CHEM_STATE) then
                    wlc_ABP(I:I+1) = INT_MIN
                endif
             enddo
          endif

          if (.False.) then
              do I = 1,WLC_P__NT
                  if (.not. isnan(wlc_RP(1,I))) then
                      print*, "should be NAN at", I
                      stop
                  endif
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
