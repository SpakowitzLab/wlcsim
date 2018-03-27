#include "../defines.inc"
! -------------------------------------------------------------------
!
!  Calculate Binding energies and x values from scratch
!  Puts output in wlc_d%DEElas, wlc_d%DE_... and wlc_d%dx_...
! -------------------------------------------------------------------
subroutine CalculateEnergiesFromScratch(wlc_p, wlc_d)
    use params
    use iso_fortran_env
    implicit none
    integer IT1, IT2, I
    real(dp) phiTot
    type(wlcsim_params), intent(in) :: wlc_p
    type(wlcsim_data), intent(inout) :: wlc_d
    integer Delta !transh

    if (WLC_P__VARIABLE_CHEM_STATE.and.WLC_P__CHANGINGCHEMICALIDENTITY) then
        wlc_d%ABP = 0 ! set entire array to zero
        !  Notide that ABP and AB are intensionally swapped below
        IT1 = 1; IT2 = wlc_p%NT
        call MC_bind(wlc_p,IT1,IT2,wlc_d%ABP,wlc_d%AB,wlc_d%METH, &
                     wlc_d%DEBind,wlc_d%dx_mu,wlc_d%demu)
    endif

    call energy_elas(wlc_d%DEELAS,wlc_d%R,wlc_d%U,wlc_p%NT,WLC_P__NB,WLC_P__NP,pack_as_para(wlc_p),&
                     WLC_P__RING,WLC_P__TWIST,wlc_p%LK,WLC_P__LT,WLC_P__L)

    ! --- Interaction Energy ---
    if (wlc_p%field_int_on_currently) then
        ! initialize phi
        call MC_int_initialize(wlc_p, wlc_d)
        phiTot=0.0_dp
        do I = 1,wlc_p%NBIN
            phiTot = phiTot + (wlc_d%PHIA(I) + wlc_d%PHIB(I))*(WLC_P__DBIN**3)
        enddo
        print*, "N-Tot", phiTot*(WLC_P__DBIN**3)/WLC_P__BEADVOLUME," NT:",wlc_p%NT
    endif

    if (WLC_P__EXPLICIT_BINDING) then
        call MC_explicit_binding_from_scratch(wlc_p,wlc_d)
    endif
  if (WLC_P__RING) then
     ! --- Initial Writhe
     call WRITHE(wlc_d%R,WLC_P__NB,wlc_d%Wr)

     !     Get initial value of Alexander polynomial and Cross matrix
     CALL ALEXANDERP(wlc_d%R,WLC_P__NB,DELTA,wlc_d%Cross,wlc_d%CrossSize,wlc_d%NCross)
     !     Begin Monte Carlo simulation

     print*, "Inside CalculateEnergiesFromScratch"
     print*, "Did I do the correct thing with Delta, NCross, Wr, ...?"
     print*, "Add the correct checks to VerifyEnergiesFromScratch"
     stop 1
  ENDif

end subroutine

subroutine InitializeEnergiesForVerifier(wlc_p, wlc_d)
    use params
    implicit none
    type(wlcsim_params), intent(in) :: wlc_p
    type(wlcsim_data), intent(inout) :: wlc_d
    ! identical to VerifyEnergiesFromScratch, but instead of checkign if they
    ! match previous values, the values are simply updated
    call CalculateEnergiesFromScratch(wlc_p, wlc_d)
    wlc_d%EBind = wlc_d%DEBind
    wlc_d%x_mu = wlc_d%dx_mu
    wlc_d%EElas = wlc_d%DEElas ! copy array
    wlc_d%eExplicitBinding = wlc_d%DEExplicitBinding
    ! --- Interaction Energy ---
    if (wlc_p%field_int_on_currently) then
        wlc_d%EChi = wlc_d%DEChi
        wlc_d%x_chi = wlc_d%dx_chi
        wlc_d%ECouple = wlc_d%DECouple
        wlc_d%x_Couple = wlc_d%dx_couple
        wlc_d%EKap = wlc_d%DEKap
        wlc_d%x_Kap = wlc_d%dx_Kap
        wlc_d%EField = wlc_d%DEField
        wlc_d%x_Field = wlc_d%dx_Field
    endif
end subroutine

subroutine VerifyEnergiesFromScratch(wlc_p, wlc_d)
    use params, only : wlcsim_params, wlcsim_data, eps, ERROR_UNIT
    implicit none
    type(wlcsim_params), intent(in) :: wlc_p
    type(wlcsim_data), intent(inout) :: wlc_d
! -------------------------------------
!
!   recalculate all energies from scratch, check them against the values they've
!   taken after being updated through various monte carlo moves
!
! -------------------------------------
    ! to save RAM...
    ! after this call, the DE's will hold the true values of the energies
    ! currently. we can then compare these to the E's
   call CalculateEnergiesFromScratch(wlc_p, wlc_d)

    ! --- Binding Energy ---
    if(abs(wlc_d%EBind-wlc_d%DEBind) > eps) then
        write(ERROR_UNIT,*) "Warning. Integrated binding enrgy:", &
                wlc_d%EBind," while absolute binding energy:", &
                wlc_d%DEBind," save point mc_ind = ",wlc_d%mc_ind
    endif
    wlc_d%EBind = wlc_d%DEBind
    if(abs(wlc_d%EMu-wlc_d%DEMu) > eps) then
        write(ERROR_UNIT,*) "Warning. Integrated chemical potential enrgy:", &
                wlc_d%EMu," while absolute chemical potential energy:", &
                wlc_d%DEMu," save point mc_ind = ",wlc_d%mc_ind
    endif
    wlc_d%EMu = wlc_d%DEMu
    wlc_d%x_mu = wlc_d%dx_mu

    ! --- Elastic Energy ---
    if(abs((wlc_d%EElas(1) +  wlc_d%EElas(2) + wlc_d%EElas(3))-&
           (wlc_d%DEElas(1) + wlc_d%DEElas(2) + wlc_d%DEElas(3))).gt.0.0001) then
        write(ERROR_UNIT,*) "Warning. Integrated elastic enrgy:", &
                (wlc_d%EElas(1) + wlc_d%EElas(2) + wlc_d%EElas(3)),&
                " while absolute elastic energy:", &
                (wlc_d%DEElas(1) + wlc_d%DEElas(2) + wlc_d%DEElas(3))
        write(ERROR_UNIT,*) " save point mc_ind = ",wlc_d%mc_ind
    endif
    wlc_d%EElas = wlc_d%DEElas ! copy array

    ! --- Explicit Binding energy ---
    if(abs(wlc_d%eExplicitBinding - wlc_d%DEExplicitBinding).gt.0.0001) then
        write(ERROR_UNIT,*) "Warning. Explicit Binding enrgy:", &
                wlc_d%eExplicitBinding, &
                " while absolute explicit binding energy:", &
                wlc_d%DEExplicitBinding
        write(ERROR_UNIT,*) " save point mc_ind = ",wlc_d%mc_ind
    endif
    wlc_d%eExplicitBinding = wlc_d%DEExplicitBinding




    ! --- Interaction Energy ---
    if (wlc_p%field_int_on_currently) then
        ! test to see if sum of changes are same as calculating from scratch
        if(abs(wlc_d%EChi-wlc_d%DEChi) > eps) then
             write(ERROR_UNIT,*) "Warning. Intigrated chi energy:", &
                     wlc_d%EChi,"  while absolute chi energy:", &
                     wlc_d%DEChi," save point mc_ind = ",wlc_d%mc_ind
        endif
        wlc_d%EChi = wlc_d%DEChi
        wlc_d%x_chi = wlc_d%dx_chi
        if(abs(wlc_d%ECouple-wlc_d%DECouple) > eps) then
             write(ERROR_UNIT,*) "Warning. Intigrated couple energy:", &
                     wlc_d%ECouple,"  while absolute couple energy:", &
                     wlc_d%DECouple," save point mc_ind = ",wlc_d%mc_ind
        endif
        wlc_d%ECouple = wlc_d%DECouple
        wlc_d%x_Couple = wlc_d%dx_couple
        if(abs(wlc_d%EKap-wlc_d%DEKap) > eps) then
             write(ERROR_UNIT,*) "Warning. Intigrated Kap energy:", &
                     wlc_d%EKap,"  while absolute Kap energy:", &
                     wlc_d%DEKap," save point mc_ind = ",wlc_d%mc_ind
        endif
        wlc_d%EKap = wlc_d%DEKap
        wlc_d%x_Kap = wlc_d%dx_Kap

        if(abs(wlc_d%EField-wlc_d%DEField) > eps) then
            write(ERROR_UNIT,*) "Warning. Integrated field enrgy:", &
                    wlc_d%EField," while absolute field energy:", &
                    wlc_d%DEField," save point mc_ind = ",wlc_d%mc_ind
        endif
        wlc_d%EField = wlc_d%DEField
        wlc_d%x_Field = wlc_d%dx_Field

        if(wlc_p%CHI_L2_ON) then
            if(abs(wlc_d%EMaiersaupe-wlc_d%deMaierSaupe) > eps) then
                write(ERROR_UNIT,*) "Warning. Integerated Maier Saupe energy:", &
                    wlc_d%EMaiersaupe," while absolute Maier Saupe energy:", &
                    wlc_d%deMaierSaupe," save points mc_ind = ",wlc_d%mc_ind
            endif
        endif
        wlc_d%EMaiersaupe = wlc_d%deMaierSaupe
        wlc_d%x_maierSaupe = wlc_d%dx_maierSaupe
    endif
end subroutine
