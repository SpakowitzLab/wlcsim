#include "../defines.inc"
! -------------------------------------------------------------------
!
!  Calculate Binding energies and x values from scratch
!  Puts output in energyOf(*_)%dE and energyOf(*_)%dx
!  Sets non-relivant dE's and dx's to zero
! -------------------------------------------------------------------
subroutine CalculateEnergiesFromScratch(wlc_p)
use params, only: wlc_METH, wlc_Cross, wlc_Wr, wlc_AB&
    , wlc_NCross, wlc_PHIB, wlc_PHIA, wlc_CrossSize, wlc_ABP&
    , wlc_R, wlc_ind_in_list, dp
use params, only: wlcsim_params
    use iso_fortran_env
    use energies
    implicit none
    integer IT1, IT2, I
    real(dp) phiTot
    type(wlcsim_params), intent(in) :: wlc_p
    integer Delta !transh
    real(dp) EELAS(4) ! Elastic force


    call set_all_dEnergy_to_zero()

    !----------------------------
    !  For each energy type calculate x as dx from scratch
    !  or set E to 0.
    !----------------------------

    if (WLC_P__VARIABLE_CHEM_STATE.and.WLC_P__CHANGINGCHEMICALIDENTITY) then
        wlc_ABP = 0 ! set entire array to zero
        !  Notide that ABP and AB are intensionally swapped below
        IT1 = 1; IT2 = WLC_P__NT
        call MC_bind(IT1,IT2,wlc_ABP,wlc_AB,wlc_METH)
    endif

    call energy_elas(EELAS,wlc_p)
    energyOf(bend_)%dx = EELAS(1)
    energyOf(stretch_)%dx = EELAS(2)
    energyOf(shear_)%dx = EELAS(3)
    energyOf(twist_)%dx = EELAS(4)

    ! ---- External Field Energy ---
    if(WLC_P__APPLY_EXTERNAL_FIELD) then
        call MC_external_field_from_scratch()
    endif
    ! ---- 2Body potentials Field Energy ---
    if(WLC_P__APPLY_2body_potential) then
        call MC_2bead_potential_from_scratch()
    endif
    ! --- Interaction Energy ---
    if (wlc_p%field_int_on_currently) then
        do I = 1,wlc_p%NBIN
            if (wlc_ind_in_list(I) .ne. -1) then
                ! Quinn put this check in to make sure that wlc_ind_in_list
                ! is reset to -1.  The program only resects values that have
                ! been changed in the move so if there is some problem and
                ! it isn't reset and incorrect answer could be produced!
                print*, "wlc_ind_in_list(",I,") should have be reset to -1"
                print*, "instead it was ",wlc_ind_in_list(I)
                stop
            endif
        enddo
        ! initialize phi
        call MC_int_initialize(wlc_p)
        phiTot=0.0_dp
        do I = 1,wlc_p%NBIN
            phiTot = phiTot + (wlc_PHIA(I) + wlc_PHIB(I))*(WLC_P__DBIN**3)
        enddo
        print*, "N-Tot", phiTot*(WLC_P__DBIN**3)/WLC_P__BEADVOLUME," NT:",WLC_P__NT
    endif

    if (WLC_P__EXPLICIT_BINDING) then
        call MC_explicit_binding_from_scratch()
    endif
    if (WLC_P__RING) then
        ! --- Initial Writhe
        call WRITHE(wlc_R,WLC_P__NB,wlc_Wr)

        !     Get initial value of Alexander polynomial and Cross matrix
        CALL ALEXANDERP(wlc_R,WLC_P__NB,DELTA,wlc_Cross,wlc_CrossSize,wlc_NCross)
        !     Begin Monte Carlo simulation

        print*, "Inside CalculateEnergiesFromScratch"
        print*, "Did I do the correct thing with Delta, NCross, Wr, ...?"
        print*, "Add the correct checks to VerifyEnergiesFromScratch"
        stop 1
    ENDif

    ! ToDo: Put from scratch calculate of self_ and confine_ energy here

    call calc_all_dE_from_dx()
end subroutine

subroutine InitializeEnergiesForVerifier(wlc_p)
    use params, only: wlcsim_params
    use energies, only: energyOf, NUMBER_OF_ENERGY_TYPES
    implicit none
    type(wlcsim_params), intent(in) :: wlc_p
    integer ii
    ! identical to VerifyEnergiesFromScratch, but instead of checkign if they
    ! match previous values, the values are simply updated into energyOf(*)%E
    call CalculateEnergiesFromScratch(wlc_p)
    do ii = 1,NUMBER_OF_ENERGY_TYPES
        energyOf(ii)%E = energyOf(ii)%dE
    enddo
end subroutine

subroutine VerifyEnergiesFromScratch(wlc_p)
! values from wlcsim_data
    use params, only : wlcsim_params,  epsapprox, ERROR_UNIT, wlc_mc_ind
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
   call CalculateEnergiesFromScratch(wlc_p)

   do ii = 1,NUMBER_OF_ENERGY_TYPES
       if(abs(energyOf(ii)%E-energyOf(ii)%dE) > epsapprox) then
           write(ERROR_UNIT,*) "Warning. Integerated ",&
               energyOf(ii)%name_str," energy:", &
               energyOf(ii)%E," while absolute ",&
               energyOf(ii)%name_str," energy:", &
               energyOf(ii)%dE," save points mc_ind = ",wlc_mc_ind
       endif
       energyOf(ii)%E = energyOf(ii)%dE
       energyOf(ii)%x = energyOf(ii)%dx
   enddo
end subroutine
