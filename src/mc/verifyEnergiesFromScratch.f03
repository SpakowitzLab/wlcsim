subroutine CalculateEnergiesFromScratch(EBind, EElas, EChi, ECouple, &
                                        EField, EKap, ESelf, EKnot, &
                                        x_mu, x_chi, x_field, x_couple, &
                                        x_kap, phiTot, mc, md)
    use params
    use iso_fortran_env

    real(dp), intent(out) :: x_mu
    real(dp), intent(out) :: x_chi
    real(dp), intent(out) :: x_field
    real(dp), intent(out) :: x_couple
    real(dp), intent(out) :: x_kap
    real(dp), intent(out) :: EBind
    real(dp), intent(out) :: EElas(N_ELAS_ENERGIES)
    real(dp), intent(out) :: EChi
    real(dp), intent(out) :: ECouple
    real(dp), intent(out) :: EField
    real(dp), intent(out) :: EKap
    real(dp), intent(out) :: ESelf
    real(dp), intent(out) :: EKnot
    real(dp), intent(out) :: phiTot ! also check total volume is consistent while we're at it
    type(wlcsim_params), intent(in) :: mc
    type(wlcsim_data), intent(out) :: md

    md%ABP=0 ! set entire array to zero
    !  Notide that ABP and AB are intensionally swapped below
    IT1=1; IT2=mc%NT
    call MC_bind(mc%NT,mc%nBpM,IT1,IT2,md%ABP,md%AB,md%METH, &
                 mc%EU,mc%EM,EBind,mc%mu,x_mu)



    call energy_elas(ELAS,md%R,md%U,mc%NT,mc%NB,mc%NP,pack_as_para(mc))

    ! --- Interaction Energy ---
    if (mc%field_int_on) then
        ! initialize phi
        IT1=1
        IT2=mc%NT ! need to set up all beads
        do I=1,mc%NBIN
             md%PHIA(I)=0.0_dp
             md%PHIB(I)=0.0_dp
        enddo
        call MC_int(mc,md,IT1,IT2,.True.)
        x_chi = md%dx_chi
        do I=1,mc%NBIN
            phiTot=phiTot+(md%PHIA(I)+md%PHIB(I))*md%Vol(I)
        enddo
    endif
  IF (mc%RING) then
     ! --- Initial Writhe
     call WRITHE(md%R,mc%NB,md%Wr)

     !     Get initial value of Alexander polynomial and Cross matrix
     CALL ALEXANDERP(md%R,mc%NB,DELTA,md%Cross,md%CrossSize,md%NCross)
     !     Begin Monte Carlo simulation
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
    CalculateEnergiesFromScratch(md%DEBind, md%DEElas, md%DEChi, md%DECouple, &
        md%DEField, md%DEKap, md%DESelf, md%DEKnot, phitot mc, md)

    ! --- Binding Energy ---
    if(abs(md%EBind-md%DEBind).gt.0.00001) then
        write(ERROR_UNIT,*) "Warning. Integrated binding enrgy:", &
                md%EBind," while absolute binding energy:", &
                md%DEBind
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
    endif
    md%EElas=md%DEElas ! copy array

    ! --- Interaction Energy ---
    if (mc%field_int_on) then
        ! test to see if sum of changes are same as calculating from scratch
        print*, "phiTot", phiTot," NT:",mc%NT
        if(abs(md%EChi-md%DEChi).gt. 0.0001_dp) then
             write(ERROR_UNIT,*) "Warning. Intigrated chi energy:", &
                     md%EChi,"  while absolute chi energy:", &
                     md%DEChi
        endif
        md%EChi=md%DEChi
        md%x_chi=md%dx_chi
        if(abs(md%ECouple-md%DECouple).gt. 0.0001_dp) then
             write(ERROR_UNIT,*) "Warning. Intigrated couple energy:", &
                     md%ECouple,"  while absolute couple energy:", &
                     md%DECouple
        endif
        md%ECouple=md%DECouple
        md%x_Couple=md%dx_couple
        if(abs(md%EKap-md%DEKap).gt. 0.0001_dp) then
             write(ERROR_UNIT,*) "Warning. Intigrated Kap energy:", &
                     md%EKap,"  while absolute Kap energy:", &
                     md%DEKap
        endif
        md%EKap=md%DEKap
        md%x_Kap=md%dx_Kap

        if(abs(md%EField-md%DEField).gt.0.00001) then
            write(ERROR_UNIT,*) "Warning. Integrated field enrgy:", &
                    md%EField," while absolute field energy:", &
                    md%DEField
        endif
        md%EField=md%DEField
        md%x_Field=md%dx_Field
    endif
end subroutine
