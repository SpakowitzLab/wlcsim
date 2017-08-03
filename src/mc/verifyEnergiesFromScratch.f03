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

    if (wlc_p%bind_on.and.wlc_p%ChangingChemicalIdentity) then
        wlc_d%ABP = 0 ! set entire array to zero
        !  Notide that ABP and AB are intensionally swapped below
        IT1 = 1; IT2 = wlc_p%NT
        call MC_bind(wlc_p%NT,wlc_p%nBpM,IT1,IT2,wlc_d%ABP,wlc_d%AB,wlc_d%METH, &
                     wlc_p%EU,wlc_p%EM,wlc_d%DEBind,wlc_p%mu,wlc_d%dx_mu)
    endif

    call energy_elas(wlc_d%DEELAS,wlc_d%R,wlc_d%U,wlc_p%NT,wlc_p%NB,wlc_p%NP,pack_as_para(wlc_p),&
                     wlc_p%Ring,wlc_p%TWIST,wlc_p%lk,wlc_p%lt,wlc_p%L)

    ! --- Interaction Energy ---
    if (wlc_p%field_int_on) then
        ! initialize phi
        call MC_int_initialize(wlc_p, wlc_d)
        do I = 1,wlc_p%NBin
            phiTot = phiTot + (wlc_d%PHIA(I) + wlc_d%PHIB(I))*wlc_d%Vol(I)
        enddo
        print*, "phiTot", phiTot," NT:",wlc_p%NT
    endif
  if (wlc_p%RinG) then
     ! --- Initial Writhe
     call WRITHE(wlc_d%R,wlc_p%NB,wlc_d%Wr)

     !     Get initial value of Alexander polynomial and Cross matrix
     CALL ALEXANDERP(wlc_d%R,wlc_p%NB,DELTA,wlc_d%Cross,wlc_d%CrossSize,wlc_d%NCross)
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
    ! --- Interaction Energy ---
    if (wlc_p%field_int_on) then
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
    use params
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
    if(abs(wlc_d%EBind-wlc_d%DEBind).gt.0.00001) then
        write(ERROR_UNIT,*) "Warning. Integrated binding enrgy:", &
                wlc_d%EBind," while absolute binding energy:", &
                wlc_d%DEBind," save point mc_ind = ",wlc_d%mc_ind
    endif
    wlc_d%EBind = wlc_d%DEBind
    wlc_d%x_mu = wlc_d%dx_mu

    ! --- Elastic Energy ---
    if(abs((wlc_d%EElas(1)+  wlc_d%EElas(2)+ wlc_d%EElas(3))-&
           (wlc_d%DEElas(1) + wlc_d%DEElas(2) + wlc_d%DEElas(3))).gt.0.0001) then
        write(ERROR_UNIT,*) "Warning. Integrated elastic enrgy:", &
                (wlc_d%EElas(1) + wlc_d%EElas(2) + wlc_d%EElas(3)),&
                " while absolute elastic energy:", &
                (wlc_d%DEElas(1) + wlc_d%DEElas(2) + wlc_d%DEElas(3))
        write(ERROR_UNIT,*) " save point mc_ind = ",wlc_d%mc_ind
    endif
    wlc_d%EElas = wlc_d%DEElas ! copy array

    ! --- Interaction Energy ---
    if (wlc_p%field_int_on) then
        ! test to see if sum of changes are same as calculating from scratch
        if(abs(wlc_d%EChi-wlc_d%DEChi).gt. 0.0001_dp) then
             write(ERROR_UNIT,*) "Warning. Intigrated chi energy:", &
                     wlc_d%EChi,"  while absolute chi energy:", &
                     wlc_d%DEChi," save point mc_ind = ",wlc_d%mc_ind
        endif
        wlc_d%EChi = wlc_d%DEChi
        wlc_d%x_chi = wlc_d%dx_chi
        if(abs(wlc_d%ECouple-wlc_d%DECouple).gt. 0.0001_dp) then
             write(ERROR_UNIT,*) "Warning. Intigrated couple energy:", &
                     wlc_d%ECouple,"  while absolute couple energy:", &
                     wlc_d%DECouple," save point mc_ind = ",wlc_d%mc_ind
        endif
        wlc_d%ECouple = wlc_d%DECouple
        wlc_d%x_Couple = wlc_d%dx_couple
        if(abs(wlc_d%EKap-wlc_d%DEKap).gt. 0.0001_dp) then
             write(ERROR_UNIT,*) "Warning. Intigrated Kap energy:", &
                     wlc_d%EKap,"  while absolute Kap energy:", &
                     wlc_d%DEKap," save point mc_ind = ",wlc_d%mc_ind
        endif
        wlc_d%EKap = wlc_d%DEKap
        wlc_d%x_Kap = wlc_d%dx_Kap

        if(abs(wlc_d%EField-wlc_d%DEField).gt.0.00001) then
            write(ERROR_UNIT,*) "Warning. Integrated field enrgy:", &
                    wlc_d%EField," while absolute field energy:", &
                    wlc_d%DEField," save point mc_ind = ",wlc_d%mc_ind
        endif
        wlc_d%EField = wlc_d%DEField
        wlc_d%x_Field = wlc_d%dx_Field
    endif
end subroutine
