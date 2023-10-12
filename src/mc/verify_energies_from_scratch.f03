#include "../defines.inc"
! -------------------------------------------------------------------
!
!  Calculate Binding energies and x values from scratch
!  Puts output in energyOf(*_)%dE and energyOf(*_)%dx
!  Sets non-relivant dE's and dx's to zero
! -------------------------------------------------------------------
subroutine calculate_energies_from_scratch(wlc_p)
   use params, only: NAN, wlc_METH, wlc_Cross, wlc_AB &
                     , wlc_NCross, wlc_PHIB, wlc_PHIA, wlc_CrossSize, wlc_ABP, WLC_GJK &
                     , wlc_R, wlc_ind_in_list, dp, wlc_bin, wlc_R_GJK, wlc_U, wlc_V
   use params, only: wlcsim_params, wlc_nucleosomeWrap, wlc_basepairs
   use umbrella, only: umbrella_energy_from_scratch
   use linkingNumber, only: link_twist_writhe_from_scratch
   use nucleosome, only: internucleosome_energy
   use iso_fortran_env
   use energies
   use binning, only: find_neighbors
   implicit none
   integer IT1, IT2, I, j, k
   real(dp) phiTot
   type(wlcsim_params), intent(in) :: wlc_p
   integer Delta !transh
   real(dp) EELAS(4) ! Elastic force
   real(dp) distances(WLC_P__NT) ! Returned distances
   integer neighbors(WLC_P__NT) ! ID of neighboring beads
   integer nn ! number of neighbors
   integer collisions
   integer ignore(WLC_P__NT)
   real(dp), parameter :: cutoff = 25.0 ! nm (cutoff for internuc attraction)
   real(dp) delInt

   call set_all_dEnergy_to_zero()

   !----------------------------
   !  For each energy type calculate x as dx from scratch
   !  or set E to 0.
   !----------------------------

   if (WLC_P__VARIABLE_CHEM_STATE .and. WLC_P__CHANGINGCHEMICALIDENTITY) then
      wlc_ABP = 0 ! set entire array to zero
      !  Notide that ABP and AB are intensionally swapped below
      IT1 = 1; IT2 = WLC_P__NT
      call mc_bind(IT1, IT2, wlc_ABP, wlc_AB, wlc_METH)
   endif

   call energy_elas(EELAS, wlc_p)
   energyOf(bend_)%dx = EELAS(1)
   energyOf(stretch_)%dx = EELAS(2)
   energyOf(shear_)%dx = EELAS(3)
   energyOf(twist_)%dx = EELAS(4)

   ! ---- External Field Energy ---
   if (WLC_P__APPLY_EXTERNAL_FIELD) then
      call mc_external_field_from_scratch()
   endif
   ! ---- 2Body potentials Field Energy ---
   if (WLC_P__APPLY_2body_potential) then
      call mc_2bead_potential_from_scratch()
   endif
   ! --- Interaction Energy ---
   if (wlc_p%field_int_on_currently) then
      do I = 1, wlc_p%NBIN
         if (wlc_ind_in_list(I) .ne. -1) then
            ! Quinn put this check in to make sure that wlc_ind_in_list
            ! is reset to -1.  The program only resects values that have
            ! been changed in the move so if there is some problem and
            ! it isn't reset and incorrect answer could be produced!
            print *, "wlc_ind_in_list(", I, ") should have be reset to -1"
            print *, "instead it was ", wlc_ind_in_list(I)
            stop
         endif
      enddo
      ! initialize phi
      call mc_int_intitialize(wlc_p)
      phiTot = 0.0_dp
      do I = 1, wlc_p%NBIN
         phiTot = phiTot + (wlc_PHIA(I) + wlc_PHIB(I))*(WLC_P__DBIN**3)
      enddo
      print *, "N-Tot", phiTot*(WLC_P__DBIN**3)/WLC_P__BEADVOLUME, " NT:", WLC_P__NT
   endif

   if (WLC_P__EXPLICIT_BINDING) then
      call mc_explicit_binding_from_scratch()
   endif
   if (WLC_P__RING) then
      if (WLC_P__TWIST) then
         ! --- Initial Writhe
         call WRITHE(wlc_R, WLC_P__NB, energyOf(global_twistLiner_)%x)

         !     Get initial value of Alexander polynomial and Cross matrix
         CALL ALEXANDERP(wlc_R, WLC_P__NB, DELTA, wlc_Cross, wlc_CrossSize, wlc_NCross)
         !     Begin Monte Carlo simulation
         print *, "Inside calculate_energies_from_scratch"
         print *, "Did I do the correct thing with Delta, NCross, Wr, ...?"
         print *, "Add the correct checks to verify_energies_from_scratch"
!           stop 1
      endif
   ENDif

   if (WLC_P__UMBRELLA) then
      call umbrella_energy_from_scratch()
   endif

   ! ToDo: Put from scratch calculate of self_ and confine_ energy here
   if (WLC_P__NO_SELF_CROSSING) then

      call link_twist_writhe_from_scratch()
   endif

   if (WLC_P__GJK_STERICS) then
      collisions = 0
      ignore = NAN
      ! check for neighbors on new beads
      do i = 1, WLC_P__NT
         ignore(i) = i
         if (WLC_P__NEIGHBOR_BINS) then
            nn = 0
            call find_neighbors(wlc_bin, wlc_R_GJK(:, i), 2*WLC_P__GJK_RADIUS, wlc_R_GJK, WLC_P__NT, &
                                WLC_P__NT, neighbors, distances, nn)
         else
            call findNeighbors(wlc_R_GJK(:, i), 2*WLC_P__GJK_RADIUS, wlc_R_GJK, WLC_P__NT, &
                               WLC_P__NT, neighbors, distances, nn)
         endif
         ! check for collisions
         call sterics_check(collisions, wlc_R, wlc_U, wlc_V, wlc_GJK, wlc_basepairs, wlc_nucleosomeWrap, ignore, &
                            i, nn, neighbors(1:nn), distances(1:nn), .true.)
      enddo
      ! ascribe collision penalty
      energyOf(sterics_)%dx = collisions
   endif

   if (WLC_P__INTERNUCLEOSOME_ON) then
      delInt = 0
      ignore = NAN
      k = 1
      do i = 1, WLC_P__NT
         if (wlc_nucleosomeWrap(i) == 0) cycle
         ignore(k) = i
         k = k + 1
         if (WLC_P__NEIGHBOR_BINS) then
            nn = 0
            call find_neighbors(wlc_bin, wlc_R_GJK(:, i), cutoff, wlc_R_GJK, WLC_P__NT, &
                                WLC_P__NT, neighbors, distances, nn)
         else
            call findNeighbors(wlc_R_GJK(:, i), cutoff, wlc_R_GJK, WLC_P__NT, &
                               WLC_P__NT, neighbors, distances, nn)
         endif
         do j = 1, nn
            if (wlc_nucleosomeWrap(neighbors(j)) == 0 .or. ANY(ignore == neighbors(j))) cycle
             delInt = delInt + internucleosome_energy(wlc_R(:, i), wlc_R(:, neighbors(j)), &
                                                     wlc_U(:, i), wlc_U(:, neighbors(j)), &
                                                     wlc_V(:, i), wlc_V(:, neighbors(j)))
         enddo
      enddo
      energyOf(internucleosome_)%dx = delInt
   endif

   call apply_energy_isOn()
   call calc_all_dE_from_dx()

end subroutine

subroutine initialize_energies_from_scratch(wlc_p)
   use params, only: wlcsim_params
   use energies, only: energyOf, NUMBER_OF_ENERGY_TYPES
   implicit none
   type(wlcsim_params), intent(in) :: wlc_p
   integer ii
   ! identical to verify_energies_from_scratch, but instead of checkign if they
   ! match previous values, the values are simply updated into energyOf(*)%E
   call calculate_energies_from_scratch(wlc_p)
   do ii = 1, NUMBER_OF_ENERGY_TYPES
      energyOf(ii)%E = energyOf(ii)%dE
   enddo
end subroutine

subroutine verify_energies_from_scratch(wlc_p)
! values from wlcsim_data
   use params, only: wlcsim_params, epsapprox, ERROR_UNIT, wlc_mc_ind, &
                     wlc_Lk, wlc_Tw, wlc_Wr, wlc_LkScratch, wlc_TwScratch, wlc_WrScratch
   use energies, only: energyOf, NUMBER_OF_ENERGY_TYPES
   implicit none
   type(wlcsim_params), intent(in) :: wlc_p
   integer ii
! -------------------------------------
!
!   recalculate all energies from scratch, check them against the values they've
!   taken after being updated through various monte carlo moves
!
! -------------------------------------
   ! to save RAM...
   ! after this call, the DE's will hold the true values of the energies
   ! currently. we can then compare these to the E's
   call calculate_energies_from_scratch(wlc_p)

   do ii = 1, NUMBER_OF_ENERGY_TYPES
      if (abs(energyOf(ii)%E - energyOf(ii)%dE) > epsapprox) then
         write (ERROR_UNIT, *) "Warning. Integerated ", &
            energyOf(ii)%name_str, " energy:", &
            energyOf(ii)%E, " while absolute ", &
            energyOf(ii)%name_str, " energy:", &
            energyOf(ii)%dE, " save points mc_ind = ", wlc_mc_ind
      endif
      energyOf(ii)%E = energyOf(ii)%dE
      energyOf(ii)%x = energyOf(ii)%dx
   enddo

   if (WLC_P__NO_SELF_CROSSING) then
      if (abs(wlc_Lk - wlc_LkScratch) > epsapprox) then
         write (ERROR_UNIT, *) "Warning. Integrated linking number:", &
            wlc_Lk, "while absolute linking number:", &
            wlc_LkScratch
      endif
      if (abs(wlc_Tw - wlc_TwScratch) > epsapprox) then
         write (ERROR_UNIT, *) "Warning. Integrated twist:", &
            wlc_Tw, "while absolute twist:", &
            wlc_TwScratch
      endif
      if (abs(wlc_Wr - wlc_WrScratch) > epsapprox) then
         write (ERROR_UNIT, *) "Warning. Integrated writhe:", &
            wlc_Wr, "while absolute writhe:", &
            wlc_WrScratch
      endif
      wlc_Lk = wlc_LkScratch
      wlc_Tw = wlc_TwScratch
      wlc_Wr = wlc_WrScratch
   endif
end subroutine
