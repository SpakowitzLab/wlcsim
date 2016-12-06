! -------------------------------------------------------------------
!
!  Calculate Binding energies and x values from scratch
!  Puts output in md%DEElas, md%DE_... and md%dx_...
! -------------------------------------------------------------------
subroutine CalculateEnergiesFromScratch(mc, md)
    use params
    use iso_fortran_env
    implicit none
    integer IT1, IT2, I
    real(dp) phiTot
    type(wlcsim_params), intent(in) :: mc
    type(wlcsim_data), intent(out) :: md
    integer Delta !transh

    if (mc%bind_on) then
        md%ABP=0 ! set entire array to zero
        !  Notide that ABP and AB are intensionally swapped below
        IT1=1; IT2=mc%NT
        call MC_bind(mc%NT,mc%nBpM,IT1,IT2,md%ABP,md%AB,md%METH, &
                     mc%EU,mc%EM,md%DEBind,mc%mu,md%dx_mu)
    endif

    call energy_elas(md%DEELAS,md%R,md%U,mc%NT,mc%NB,mc%NP,pack_as_para(mc))

    ! --- Interaction Energy ---
    if (mc%field_int_on) then
        ! initialize phi
        call MC_int_initialize(mc,md)
        do I=1,mc%NBIN
            phiTot=phiTot+(md%PHIA(I)+md%PHIB(I))*md%Vol(I)
        enddo
        print*, "phiTot", phiTot," NT:",mc%NT
    endif
  IF (mc%RING) then
     ! --- Initial Writhe
     call WRITHE(md%R,mc%NB,md%Wr)

     !     Get initial value of Alexander polynomial and Cross matrix
     CALL ALEXANDERP(md%R,mc%NB,DELTA,md%Cross,md%CrossSize,md%NCross)
     !     Begin Monte Carlo simulation

     print*, "Inside CalculateEnergiesFromScratch"
     print*, "Did I do the correct thing with Delta, NCross, Wr, ...?"
     print*, "Add the correct checks to VerifyEnergiesFromScratch"
     stop 1
  ENDIF

end subroutine
subroutine VerifyEnegiesFromScratch(mc, md)
    use params
    implicit none
    type(wlcsim_params), intent(in) :: mc
    type(wlcsim_data), intent(inout) :: md
! -------------------------------------
!
!   recalculate all energies from scratch, check them against the values they've
!   taken after being updated through various monte carlo moves
!
! -------------------------------------
    ! to save RAM...
    ! after this call, the DE's will hold the true values of the energies
    ! currently. we can then compare these to the E's
   call CalculateEnergiesFromScratch(mc, md)

    ! --- Binding Energy ---
    if(abs(md%EBind-md%DEBind).gt.0.00001) then
        write(ERROR_UNIT,*) "Warning. Integrated binding enrgy:", &
                md%EBind," while absolute binding energy:", &
                md%DEBind," save point mc_ind=",md%mc_ind
    endif
    md%EBind=md%DEBind
    md%x_mu=md%dx_mu

    ! --- Elastic Energy ---
    if(abs((md%EElas(1)+  md%EElas(2)+ md%EElas(3))-&
           (md%DEElas(1)+md%DEElas(2)+md%DEElas(3))).gt.0.0001) then
        write(ERROR_UNIT,*) "Warning. Integrated elastic enrgy:", &
                (md%EElas(1)+md%EElas(2)+md%EElas(3)),&
                " while absolute elastic energy:", &
                (md%DEElas(1)+md%DEElas(2)+md%DEElas(3))
        write(ERROR_UNIT,*) " save point mc_ind=",md%mc_ind
    endif
    md%EElas=md%DEElas ! copy array

    ! --- Interaction Energy ---
    if (mc%field_int_on) then
        ! test to see if sum of changes are same as calculating from scratch
        if(abs(md%EChi-md%DEChi).gt. 0.0001_dp) then
             write(ERROR_UNIT,*) "Warning. Intigrated chi energy:", &
                     md%EChi,"  while absolute chi energy:", &
                     md%DEChi," save point mc_ind=",md%mc_ind
        endif
        md%EChi=md%DEChi
        md%x_chi=md%dx_chi
        if(abs(md%ECouple-md%DECouple).gt. 0.0001_dp) then
             write(ERROR_UNIT,*) "Warning. Intigrated couple energy:", &
                     md%ECouple,"  while absolute couple energy:", &
                     md%DECouple," save point mc_ind=",md%mc_ind
        endif
        md%ECouple=md%DECouple
        md%x_Couple=md%dx_couple
        if(abs(md%EKap-md%DEKap).gt. 0.0001_dp) then
             write(ERROR_UNIT,*) "Warning. Intigrated Kap energy:", &
                     md%EKap,"  while absolute Kap energy:", &
                     md%DEKap," save point mc_ind=",md%mc_ind
        endif
        md%EKap=md%DEKap
        md%x_Kap=md%dx_Kap

        if(abs(md%EField-md%DEField).gt.0.00001) then
            write(ERROR_UNIT,*) "Warning. Integrated field enrgy:", &
                    md%EField," while absolute field energy:", &
                    md%DEField," save point mc_ind=",md%mc_ind
        endif
        md%EField=md%DEField
        md%x_Field=md%dx_Field
    endif
end subroutine
