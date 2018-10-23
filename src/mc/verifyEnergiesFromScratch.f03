#include "../defines.inc"
! -------------------------------------------------------------------
!
!  Calculate Binding energies and x values from scratch
!  Puts output in wlc_DEElas, wlc_DE_... and wlc_dx_...
! -------------------------------------------------------------------
subroutine CalculateEnergiesFromScratch(wlc_p)
! values from wlcsim_data
use params, only: wlc_METH, wlc_Cross, wlc_Wr, wlc_AB, wlc_dx_mu&
    , wlc_NCross, wlc_PHIB, wlc_PHIA, wlc_DEELAS, wlc_CrossSize, wlc_ABP&
    , wlc_demu, wlc_DEBind, wlc_R, wlc_ind_in_list, dp
use params, only: wlcsim_params
    use iso_fortran_env
    implicit none
    integer IT1, IT2, I
    real(dp) phiTot
    type(wlcsim_params), intent(in) :: wlc_p
    integer Delta !transh

    if (WLC_P__VARIABLE_CHEM_STATE.and.WLC_P__CHANGINGCHEMICALIDENTITY) then
        wlc_ABP = 0 ! set entire array to zero
        !  Notide that ABP and AB are intensionally swapped below
        IT1 = 1; IT2 = WLC_P__NT
        call MC_bind(wlc_p,IT1,IT2,wlc_ABP,wlc_AB,wlc_METH, &
                     wlc_DEBind,wlc_dx_mu,wlc_demu)
    endif

    call energy_elas(wlc_DEELAS,wlc_p)

    ! ---- External Field Energy ---
    if(WLC_P__APPLY_EXTERNAL_FIELD) then
        call MC_external_field_from_scratch(wlc_p)
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

end subroutine

subroutine InitializeEnergiesForVerifier(wlc_p)
! values from wlcsim_data
use params, only: wlc_x_Kap, wlc_dx_mu, wlc_x_Field, wlc_dx_Kap, wlc_EChi&
    , wlc_EExternalField, wlc_ECouple, wlc_DEExplicitBinding, wlc_EBind, wlc_DEField, wlc_eExplicitBinding&
    , wlc_EKap, wlc_DEBind, wlc_x_chi, wlc_x_mu, wlc_dx_ExternalField, wlc_dx_couple&
    , wlc_dx_Field, wlc_DECouple, wlc_x_ExternalField, wlc_DEExternalField, wlc_EField, wlc_DEElas&
    , wlc_x_Couple, wlc_DEKap, wlc_EElas, wlc_DEChi, wlc_dx_chi
    use params
    implicit none
    type(wlcsim_params), intent(in) :: wlc_p
    ! identical to VerifyEnergiesFromScratch, but instead of checkign if they
    ! match previous values, the values are simply updated
    call CalculateEnergiesFromScratch(wlc_p)
    wlc_EBind = wlc_DEBind
    wlc_x_mu = wlc_dx_mu
    wlc_EElas = wlc_DEElas ! copy array
    wlc_eExplicitBinding = wlc_DEExplicitBinding
    ! ---- External Field Energy ---
    if(WLC_P__APPLY_EXTERNAL_FIELD) then
        wlc_EExternalField = wlc_DEExternalField
        wlc_x_ExternalField = wlc_dx_ExternalField
    endif
    ! --- Interaction Energy ---
    if (wlc_p%field_int_on_currently) then
        wlc_EChi = wlc_DEChi
        wlc_x_chi = wlc_dx_chi
        wlc_ECouple = wlc_DECouple
        wlc_x_Couple = wlc_dx_couple
        wlc_EKap = wlc_DEKap
        wlc_x_Kap = wlc_dx_Kap
        wlc_EField = wlc_DEField
        wlc_x_Field = wlc_dx_Field
    endif
end subroutine

subroutine VerifyEnergiesFromScratch(wlc_p)
! values from wlcsim_data
use params, only: wlc_x_Kap, wlc_dx_mu, wlc_x_Field, wlc_dx_Kap, wlc_x_maierSaupe&
    , wlc_EChi, wlc_DEMu, wlc_EExternalField, wlc_ECouple, wlc_DEExplicitBinding, wlc_EBind&
    , wlc_DEField, wlc_eExplicitBinding, wlc_DEBind, wlc_EKap, wlc_eExternalField, wlc_x_chi&
    , wlc_EMaiersaupe, wlc_x_mu, wlc_dx_ExternalField, wlc_dx_couple, wlc_EMu, wlc_dx_Field&
    , wlc_DECouple, wlc_mc_ind, wlc_x_ExternalField, wlc_DEExternalField, wlc_EField, wlc_dx_maierSaupe&
    , wlc_deMaierSaupe, wlc_DEElas, wlc_x_Couple, wlc_DEKap, wlc_EElas, wlc_DEChi&
    , wlc_dx_chi, dp
    use params, only : wlcsim_params,  epsapprox, ERROR_UNIT
    implicit none
    type(wlcsim_params), intent(in) :: wlc_p
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

    ! --- Binding Energy ---
    if(abs(wlc_EBind-wlc_DEBind) > epsapprox) then
        write(ERROR_UNIT,*) "Warning. Integrated binding enrgy:", &
                wlc_EBind," while absolute binding energy:", &
                wlc_DEBind," save point mc_ind = ",wlc_mc_ind
    endif
    wlc_EBind = wlc_DEBind
    if(abs(wlc_EMu-wlc_DEMu) > epsapprox) then
        write(ERROR_UNIT,*) "Warning. Integrated chemical potential enrgy:", &
                wlc_EMu," while absolute chemical potential energy:", &
                wlc_DEMu," save point mc_ind = ",wlc_mc_ind
    endif
    wlc_EMu = wlc_DEMu
    wlc_x_mu = wlc_dx_mu

    ! --- Elastic Energy ---
    if(abs((wlc_EElas(1) +  wlc_EElas(2) + wlc_EElas(3))-&
           (wlc_DEElas(1) + wlc_DEElas(2) + wlc_DEElas(3))) .gt. 0.0001_dp) then
        write(ERROR_UNIT,*) "Warning. Integrated elastic enrgy:", &
                (wlc_EElas(1) + wlc_EElas(2) + wlc_EElas(3)),&
                " while absolute elastic energy:", &
                (wlc_DEElas(1) + wlc_DEElas(2) + wlc_DEElas(3))
        write(ERROR_UNIT,*) " save point mc_ind = ",wlc_mc_ind
    endif
    wlc_EElas = wlc_DEElas ! copy array

    ! --- Explicit Binding energy ---
    if(abs(wlc_eExplicitBinding - wlc_DEExplicitBinding) .gt. 0.0001_dp) then
        write(ERROR_UNIT,*) "Warning. Explicit Binding enrgy:", &
                wlc_eExplicitBinding, &
                " while absolute explicit binding energy:", &
                wlc_DEExplicitBinding
        write(ERROR_UNIT,*) " save point mc_ind = ",wlc_mc_ind
    endif
    wlc_eExplicitBinding = wlc_DEExplicitBinding


    ! ---- External Field Energy ---
    if(WLC_P__APPLY_EXTERNAL_FIELD) then
        if(abs(wlc_eExternalField-wlc_DEExternalField) > epsapprox) then
            write(ERROR_UNIT,*) "Warning. Integrated external field enrgy:", &
                    wlc_EExternalField," while absolute external field energy:", &
                    wlc_DEExternalField," save point mc_ind = ",wlc_mc_ind
        endif
        wlc_EExternalField = wlc_DEExternalField
        wlc_x_ExternalField = wlc_dx_ExternalField

    endif

    ! --- Interaction Energy ---
    if (wlc_p%field_int_on_currently) then
        ! test to see if sum of changes are same as calculating from scratch
        if(abs(wlc_EChi-wlc_DEChi) > epsapprox) then
             write(ERROR_UNIT,*) "Warning. Intigrated chi energy:", &
                     wlc_EChi,"  while absolute chi energy:", &
                     wlc_DEChi," save point mc_ind = ",wlc_mc_ind
        endif
        wlc_EChi = wlc_DEChi
        wlc_x_chi = wlc_dx_chi
        if(abs(wlc_ECouple-wlc_DECouple) > epsapprox) then
             write(ERROR_UNIT,*) "Warning. Intigrated couple energy:", &
                     wlc_ECouple,"  while absolute couple energy:", &
                     wlc_DECouple," save point mc_ind = ",wlc_mc_ind
        endif
        wlc_ECouple = wlc_DECouple
        wlc_x_Couple = wlc_dx_couple
        if(abs(wlc_EKap-wlc_DEKap) > epsapprox) then
             write(ERROR_UNIT,*) "Warning. Intigrated Kap energy:", &
                     wlc_EKap,"  while absolute Kap energy:", &
                     wlc_DEKap," save point mc_ind = ",wlc_mc_ind
        endif
        wlc_EKap = wlc_DEKap
        wlc_x_Kap = wlc_dx_Kap

        if(abs(wlc_EField-wlc_DEField) > epsapprox) then
            write(ERROR_UNIT,*) "Warning. Integrated field enrgy:", &
                    wlc_EField," while absolute field energy:", &
                    wlc_DEField," save point mc_ind = ",wlc_mc_ind
        endif
        wlc_EField = wlc_DEField
        wlc_x_Field = wlc_dx_Field

        if(wlc_p%CHI_L2_ON) then
            if(abs(wlc_EMaiersaupe-wlc_deMaierSaupe) > epsapprox) then
                write(ERROR_UNIT,*) "Warning. Integerated Maier Saupe energy:", &
                    wlc_EMaiersaupe," while absolute Maier Saupe energy:", &
                    wlc_deMaierSaupe," save points mc_ind = ",wlc_mc_ind
            endif
        endif
        wlc_EMaiersaupe = wlc_deMaierSaupe
        wlc_x_maierSaupe = wlc_dx_maierSaupe
    endif
end subroutine
